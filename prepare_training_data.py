"""Harmonize the five training datasets into cell-load's expected layout.

Output: one h5ad per cell type under --out-dir, grouped by dataset:

    replogle/K562.h5ad, replogle/RPE1.h5ad
    nadig/HepG2.h5ad,  nadig/Jurkat.h5ad
    h1/H1.h5ad

Each file: X = raw counts (CSR float32) on the shared gene vocabulary (the
intersection of gene symbols across all five datasets), obs columns
`gene` (perturbation, 'non-targeting' for controls), `gem_group` (batch),
`cell_type`. Written in row chunks and concatenated on disk to bound memory.

Usage:
    python prepare_training_data.py --data-dir C:/Users/aaron/vcc_data \
        --out-dir C:/Users/aaron/vcc_data/state_train [--only K562]
"""

import argparse
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    from anndata.experimental import concat_on_disk
except ImportError:
    from anndata import concat_on_disk

# cell-load reads obs/_index with raw h5py and only understands the classic
# string-array encoding; anndata's default nullable-string-array breaks it.
ad.settings.allow_write_nullable_strings = False

CHUNK_ROWS = 50_000

# (dataset, cell_type, relative path, symbol source, pert col, batch col)
SPECS = [
    ("replogle", "K562", "replogle/K562_essential_raw_singlecell.h5ad", "gene_name", "gene", "gem_group"),
    ("replogle", "RPE1", "replogle/rpe1_raw_singlecell.h5ad", "gene_name", "gene", "gem_group"),
    ("nadig", "HepG2", "nadig/hepg2_raw_singlecell.h5ad", "gene_name", "gene", "gem_group"),
    ("nadig", "Jurkat", "nadig/jurkat_raw_singlecell.h5ad", "gene_name", "gene", "gem_group"),
    ("h1", "H1", "h1_2025_adata_Training.h5ad", None, "target_gene", "batch"),
]


def symbols_of(a: ad.AnnData, sym_col: str | None) -> pd.Series:
    return a.var[sym_col].astype(str) if sym_col else a.var_names.to_series().astype(str)


def shared_gene_vocabulary(data_dir: Path) -> list[str]:
    sets = []
    for _, _, rel, sym_col, _, _ in SPECS:
        a = ad.read_h5ad(data_dir / rel, backed="r")
        sets.append(set(symbols_of(a, sym_col)))
    return sorted(set.intersection(*sets))


def convert(data_dir: Path, out: Path, spec: tuple, genes: list[str]) -> None:
    _, cell_type, rel, sym_col, pert_col, batch_col = spec
    a = ad.read_h5ad(data_dir / rel, backed="r")

    # first column position for each symbol (duplicate symbols are rare and
    # ambiguous; keep the first)
    pos: dict[str, int] = {}
    for i, s in enumerate(symbols_of(a, sym_col)):
        pos.setdefault(s, i)
    cols = np.array([pos[g] for g in genes])
    var = pd.DataFrame(index=pd.Index(genes, name=None))

    n = a.shape[0]
    tmpdir = Path(tempfile.mkdtemp(prefix=f"state_prep_{cell_type}_"))
    chunks = []
    try:
        for start in range(0, n, CHUNK_ROWS):
            sub = a[start:min(start + CHUNK_ROWS, n)].to_memory()
            X = sub.X if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
            X = X.tocsr()[:, cols].astype(np.float32)
            obs = pd.DataFrame({
                "gene": pd.Categorical(sub.obs[pert_col].astype(str).values),
                "gem_group": pd.Categorical(sub.obs[batch_col].astype(str).values),
                "cell_type": pd.Categorical([cell_type] * X.shape[0]),
            }, index=pd.Index([f"{cell_type}_{start + i}" for i in range(X.shape[0])]))
            p = tmpdir / f"chunk_{start:07d}.h5ad"
            ad.AnnData(X=X, obs=obs, var=var).write_h5ad(p, compression="gzip")
            chunks.append(p)
            print(f"  {cell_type}: {min(start + CHUNK_ROWS, n):,}/{n:,} cells", flush=True)

        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists():
            out.unlink()
        concat_on_disk([str(p) for p in chunks], str(out), axis=0, join="inner")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    b = ad.read_h5ad(out, backed="r")
    nt = int((b.obs["gene"] == "non-targeting").sum())
    print(f"  wrote {out}  {b.shape[0]:,} x {b.shape[1]:,}  "
          f"({b.obs['gene'].nunique():,} perts, {nt:,} control cells)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=Path("C:/Users/aaron/vcc_data"))
    ap.add_argument("--out-dir", type=Path, default=Path("C:/Users/aaron/vcc_data/state_train"))
    ap.add_argument("--only", type=str, default=None,
                    help="convert just one cell type (smoke test)")
    args = ap.parse_args()

    genes = shared_gene_vocabulary(args.data_dir)
    print(f"shared gene vocabulary: {len(genes):,} symbols")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.Series(genes).to_csv(args.out_dir / "gene_vocab.csv", index=False, header=False)

    for spec in SPECS:
        dataset, cell_type = spec[0], spec[1]
        if args.only and cell_type != args.only:
            continue
        print(f"converting {dataset}/{cell_type}")
        convert(args.data_dir, args.out_dir / dataset / f"{cell_type}.h5ad", spec, genes)


if __name__ == "__main__":
    main()
