"""Build a VCC submission from the predict-no-change baseline: for each
perturbation in each context, resample 400 control cells from that context and
relabel them with the perturbation.

Written in per-context chunks and concatenated on disk; a full submission is
~1.8e9 stored entries, ~14 GB as a single in-memory CSR.

Usage:
    # smoke test, 5 perturbations
    python baseline_submit.py --controls-dir controls --n-perts 5 \
        --out prediction_smoke.h5ad

    # full run
    python baseline_submit.py --controls-dir controls --out prediction.h5ad

Validate before uploading:
    vcc prep prediction.h5ad -g controls/gene_names.csv \
        --perts controls/pert_counts.csv -o prediction.vcc --dry-run
"""

import argparse
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    from anndata.experimental import concat_on_disk
except ImportError:  # newer anndata promoted it
    from anndata import concat_on_disk

CELLS_PER_PERT = 400
STORED_ENTRY_CAP = 4_750_000_000
MAX_COUNTS_PER_CELL = 1_000_000
CONTEXTS = ("A", "B", "C")


def load_gene_names(path: Path) -> list[str]:
    """Read the ordered gene symbol list, tolerating an optional header row."""
    first = pd.read_csv(path, header=None, nrows=1).iloc[0, 0]
    header = 0 if str(first).strip().lower() in {"gene", "gene_name", "gene_names"} else None
    return pd.read_csv(path, header=header).iloc[:, 0].astype(str).tolist()


def load_controls(path: Path, gene_names: list[str]) -> sp.csr_matrix:
    """Load one context's control counts as float32 CSR, columns in
    gene_names order (reindexed if the file's var order differs)."""
    adata = ad.read_h5ad(path)
    if list(adata.var_names) != gene_names:
        missing = set(gene_names) - set(adata.var_names)
        if missing:
            raise ValueError(
                f"{path.name} is missing {len(missing)} genes required by "
                f"gene_names.csv, e.g. {sorted(missing)[:5]}"
            )
        adata = adata[:, gene_names]
    X = adata.X
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    return X.astype(np.float32)


def apply_multipliers(
    block: sp.csr_matrix, mult: np.ndarray, rng: np.random.Generator
) -> sp.csr_matrix:
    """Scale each stored count by its gene's multiplier, then stochastically
    round (floor + Bernoulli on the fraction) so counts stay integral with
    the expected value exact. Only nonzero entries are touched: a gene at
    zero stays zero even if its multiplier says up."""
    vals = block.data * mult[block.indices]
    floor = np.floor(vals)
    block.data = (floor + (rng.random(vals.size) < (vals - floor))).astype(np.float32)
    block.eliminate_zeros()
    return block


def load_multipliers(
    path: Path, gene_names: list[str], target_residual: float
) -> dict[str, np.ndarray]:
    """Load per-perturbation LFC signatures (perts x genes AnnData from
    compute_signatures.py) as fold-change multiplier vectors over the full
    gene panel. The target gene's own multiplier is forced to
    target_residual, matching the challenge's stated CRISPRi efficacy."""
    sig = ad.read_h5ad(path)
    col_of = {g: i for i, g in enumerate(gene_names)}
    sig_cols = np.array([col_of[g] for g in sig.var_names])
    lfc = np.asarray(sig.X)

    mults = {}
    for i, pert in enumerate(sig.obs_names):
        m = np.ones(len(gene_names), dtype=np.float32)
        m[sig_cols] = 2.0 ** lfc[i]
        if pert in col_of:
            m[col_of[pert]] = target_residual
        mults[str(pert)] = m
    return mults


def knock_down_gene(
    block: sp.csr_matrix, col: int, residual: float, rng: np.random.Generator
) -> sp.csr_matrix:
    """Binomially thin one gene's counts: each molecule survives with
    probability `residual`. Keeps counts integral and zeros sparse."""
    block = block.tocsc()
    start, end = block.indptr[col], block.indptr[col + 1]
    vals = block.data[start:end]
    block.data[start:end] = rng.binomial(vals.astype(np.int64), residual).astype(np.float32)
    block.eliminate_zeros()
    return block.tocsr()


def build_chunk(
    ctrl: sp.csr_matrix,
    perts: list[str],
    context: str,
    gene_names: list[str],
    rng: np.random.Generator,
    knockdown_residual: float | None = None,
    multipliers: dict[str, np.ndarray] | None = None,
) -> ad.AnnData:
    """Predictions for one block of perturbations in one context.

    Base behavior (predict-no-change): for each perturbation, sample
    CELLS_PER_PERT control cells with replacement and relabel them. With
    knockdown_residual set, additionally thin the target gene's own counts
    to that fraction (target-kd). With multipliers set, scale every gene by
    that perturbation's fold-change vector instead (global-mean); perts
    without a signature fall back to target-kd. Returns an AnnData with obs
    columns `target_gene` and `context`. Any real model replaces this
    function and keeps the same signature.
    """
    n_ctrl = ctrl.shape[0]
    col_of = {g: i for i, g in enumerate(gene_names)}
    blocks, target_gene, obs_names = [], [], []

    for pert in perts:
        idx = rng.choice(n_ctrl, size=CELLS_PER_PERT, replace=True)
        block = ctrl[idx]
        if multipliers is not None and pert in multipliers:
            block = apply_multipliers(block, multipliers[pert], rng)
        elif knockdown_residual is not None and pert in col_of:
            block = knock_down_gene(block, col_of[pert], knockdown_residual, rng)
        blocks.append(block)
        target_gene.extend([pert] * CELLS_PER_PERT)
        obs_names.extend(f"{context}_{pert}_{i:04d}" for i in range(CELLS_PER_PERT))

    X = sp.vstack(blocks, format="csr")
    obs = pd.DataFrame(
        {"target_gene": target_gene, "context": [context] * X.shape[0]},
        index=pd.Index(obs_names, name=None),
    )
    var = pd.DataFrame(index=pd.Index(gene_names, name=None))
    return ad.AnnData(X=X, obs=obs, var=var)


def verify(path: Path, perts: list[str], gene_names: list[str]) -> None:
    """Re-open the written file in backed mode and check it against the
    submission rules: row count, gene order, no non-targeting rows, 400 cells
    per (context, perturbation), stored-entry cap, and non-negative integer
    counts on a sample slice. Advisory only; `vcc prep` is the authority."""
    print(f"\n=== verifying {path.name} ===")
    adata = ad.read_h5ad(path, backed="r")

    n_expected = len(perts) * CELLS_PER_PERT * len(CONTEXTS)
    print(f"  shape                 {adata.shape[0]:,} x {adata.shape[1]:,}")
    print(f"  expected rows         {n_expected:,}  "
          f"{'OK' if adata.shape[0] == n_expected else 'MISMATCH'}")
    print(f"  gene order matches    {list(adata.var_names) == gene_names}")

    obs = adata.obs
    print(f"  obs columns           {list(obs.columns)}")
    print(f"  contexts              {sorted(obs['context'].unique())}")
    print(f"  n perturbations       {obs['target_gene'].nunique()}")

    has_ntc = obs["target_gene"].astype(str).str.lower().isin(
        {"non-targeting", "non_targeting", "ntc", "control"}
    ).any()
    print(f"  contains NTC rows     {has_ntc}  {'REJECTED' if has_ntc else 'OK'}")

    per = obs.groupby(["context", "target_gene"], observed=True).size()
    bad = per[per != CELLS_PER_PERT]
    print(f"  cells per (ctx,pert)  "
          f"{'all 400, OK' if bad.empty else f'{len(bad)} groups wrong'}")
    if not bad.empty:
        print(f"    e.g. {bad.head().to_dict()}")

    # Backed sparse matrices don't always expose .nnz; read it from the file.
    X = adata.X
    if hasattr(X, "nnz"):
        nnz = X.nnz
    else:
        with h5py.File(path, "r") as f:
            nnz = int(f["X/indptr"][-1]) if "X/indptr" in f else None
    if nnz is not None:
        print(f"  stored entries        {nnz:,}  cap {STORED_ENTRY_CAP:,}  "
              f"{'OK' if nnz <= STORED_ENTRY_CAP else 'OVER CAP'}")
        print(f"  per cell              {nnz / adata.shape[0]:,.0f}")

    # Spot-check the count requirements on a slice rather than the whole matrix.
    sample = adata[:2000].to_memory().X
    sample = sp.csr_matrix(sample) if not sp.issparse(sample) else sample
    data = sample.data
    print(f"  sample non-negative   {bool(np.all(data >= 0))}")
    print(f"  sample integral       {bool(np.all(data == np.floor(data)))}")
    print(f"  sample finite         {bool(np.all(np.isfinite(data)))}")
    totals = np.asarray(sample.sum(axis=1)).ravel()
    print(f"  sample max cell total {totals.max():,.0f}  "
          f"cap {MAX_COUNTS_PER_CELL:,}  "
          f"{'OK' if totals.max() <= MAX_COUNTS_PER_CELL else 'OVER CAP'}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls-dir", type=Path, default=Path("controls"))
    ap.add_argument("--out", type=Path, default=Path("prediction.h5ad"))
    ap.add_argument("--chunk-perts", type=int, default=25,
                    help="perturbations per on-disk chunk; lower it if you run out of RAM")
    ap.add_argument("--n-perts", type=int, default=None,
                    help="use only the first N perturbations, for smoke tests")
    ap.add_argument("--model", choices=["no-change", "target-kd", "global-mean"],
                    default="no-change",
                    help="no-change: resampled controls as-is; "
                         "target-kd: also thin the target gene's counts; "
                         "global-mean: apply per-pert LFC signatures "
                         "(context-blind), needs --signatures")
    ap.add_argument("--signatures", type=Path, default=None,
                    help="perts x genes LFC AnnData from compute_signatures.py")
    ap.add_argument("--residual", type=float, default=0.15,
                    help="surviving fraction of target-gene counts "
                         "(target-kd, and forced on-target value for global-mean)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    residual = args.residual if args.model in ("target-kd", "global-mean") else None

    d = args.controls_dir
    gene_names = load_gene_names(d / "gene_names.csv")
    perts = pd.read_csv(d / "pert_counts.csv")["target_gene"].astype(str).tolist()
    if args.n_perts:
        perts = perts[: args.n_perts]
        print(f"SMOKE TEST: {len(perts)} perturbations only, not a valid submission")

    multipliers = None
    if args.model == "global-mean":
        if args.signatures is None:
            ap.error("--model global-mean requires --signatures")
        multipliers = load_multipliers(args.signatures, gene_names, args.residual)
        covered = sum(p in multipliers for p in perts)
        print(f"signatures cover {covered}/{len(perts)} perturbations "
              f"(rest fall back to target-kd)")

    print(f"model: {args.model}"
          + (f" (residual {residual})" if residual is not None else ""))
    print(f"{len(gene_names):,} genes, {len(perts):,} perturbations, "
          f"{len(CONTEXTS)} contexts")
    print(f"target rows: {len(perts) * CELLS_PER_PERT * len(CONTEXTS):,}")

    rng = np.random.default_rng(args.seed)
    tmpdir = Path(tempfile.mkdtemp(prefix="vcc_chunks_"))
    chunk_paths = []

    try:
        for ctx in CONTEXTS:
            ctrl = load_controls(d / f"context_{ctx}.h5ad", gene_names)
            print(f"\ncontext {ctx}: {ctrl.shape[0]:,} control cells, "
                  f"{ctrl.nnz / ctrl.shape[0]:,.0f} nnz/cell")

            for start in range(0, len(perts), args.chunk_perts):
                block = perts[start : start + args.chunk_perts]
                chunk = build_chunk(ctrl, block, ctx, gene_names, rng,
                                    knockdown_residual=residual,
                                    multipliers=multipliers)
                p = tmpdir / f"chunk_{ctx}_{start:05d}.h5ad"
                chunk.write_h5ad(p, compression="gzip")
                chunk_paths.append(p)
                print(f"  wrote {p.name}  {chunk.shape[0]:,} cells  "
                      f"{chunk.X.nnz:,} nnz")
                del chunk

            del ctrl

        print(f"\nconcatenating {len(chunk_paths)} chunks on disk -> {args.out}")
        if args.out.exists():
            args.out.unlink()
        concat_on_disk([str(p) for p in chunk_paths], str(args.out), axis=0, join="inner")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    verify(args.out, perts, gene_names)

    print("\nnext:")
    print(f"  vcc prep {args.out} -g {d / 'gene_names.csv'} "
          f"--perts {d / 'pert_counts.csv'} -o prediction.vcc --dry-run")


if __name__ == "__main__":
    main()
