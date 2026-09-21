"""Convert Replogle K562 genome-wide single-cell data into the pretrained
ST-HVG-Replogle representation, for fine-tuning with challenge-perturbation
coverage.

Representation (reverse-engineered from Arc's published eval outputs and
verified against raw data): per cell, normalize total counts over the file's
full gene panel to 10,000, log1p, then subset to the checkpoint's 2,000 HVGs
in its exact order. Stored as X (the model's input and output are both this
2,000-dim space). obs: `gene` (perturbation), `gem_group`, `cell_line`.

Usage:
    python prepare_finetune_data.py \
        --gwps C:/Users/aaron/vcc_data/replogle/K562_gwps_raw_singlecell.h5ad \
        --var-dims C:/Users/aaron/vcc_data/pretrained/st_hvg_replogle_fewshot_jurkat/var_dims.pkl \
        --out C:/Users/aaron/vcc_data/state_finetune/replogle_gwps/K562gw.h5ad
"""

import argparse
import pickle
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

ad.settings.allow_write_nullable_strings = False  # classic encoding for cell-load

CHUNK_ROWS = 100_000
TARGET_SUM = 1e4


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwps", type=Path, required=True)
    ap.add_argument("--var-dims", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    with open(args.var_dims, "rb") as f:
        hvgs = [str(g) for g in pickle.load(f)["gene_names"]]

    a = ad.read_h5ad(args.gwps, backed="r")
    sym = a.var["gene_name"].astype(str)
    pos: dict[str, int] = {}
    for i, s in enumerate(sym):
        pos.setdefault(s, i)
    missing = [g for g in hvgs if g not in pos]
    if missing:
        raise SystemExit(f"{len(missing)} HVGs missing from gwps: {missing[:5]}")
    cols = np.array([pos[g] for g in hvgs])
    var = pd.DataFrame(index=pd.Index(hvgs, name=None))

    n = a.shape[0]
    print(f"{n:,} cells; normalizing over {a.shape[1]:,} genes -> {len(hvgs)} HVGs")
    tmpdir = Path(tempfile.mkdtemp(prefix="finetune_prep_"))
    chunks = []
    try:
        for start in range(0, n, CHUNK_ROWS):
            sub = a[start:min(start + CHUNK_ROWS, n)].to_memory()
            X = sub.X if sp.issparse(sub.X) else sp.csr_matrix(np.asarray(sub.X))
            totals = np.asarray(X.sum(axis=1)).ravel()
            totals[totals == 0] = 1.0
            H = X.tocsr()[:, cols].astype(np.float32)
            H = sp.csr_matrix(H.multiply(TARGET_SUM / totals[:, None]))
            H.data = np.log1p(H.data)
            obs = pd.DataFrame({
                "gene": pd.Categorical(sub.obs["gene"].astype(str).values),
                "gem_group": pd.Categorical(sub.obs["gem_group"].astype(str).values),
                "cell_line": pd.Categorical(["K562gw"] * H.shape[0]),
            }, index=pd.Index([f"K562gw_{start + i}" for i in range(H.shape[0])]))
            p = tmpdir / f"chunk_{start:08d}.h5ad"
            ad.AnnData(X=H, obs=obs, var=var).write_h5ad(p)
            chunks.append(p)
            print(f"  {min(start + CHUNK_ROWS, n):,}/{n:,}", flush=True)

        args.out.parent.mkdir(parents=True, exist_ok=True)
        if args.out.exists():
            args.out.unlink()
        concat_on_disk([str(p) for p in chunks], str(args.out), axis=0, join="inner")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    b = ad.read_h5ad(args.out, backed="r")
    nt = int((b.obs["gene"] == "non-targeting").sum())
    print(f"wrote {args.out}: {b.shape[0]:,} x {b.shape[1]:,}, "
          f"{b.obs['gene'].nunique():,} perts, {nt:,} control cells")


if __name__ == "__main__":
    main()
