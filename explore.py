"""Inspect the VCC control bundle: shapes, count depth, sparsity, gene order,
and the projected stored-entry count of a full submission vs the 4.75e9 cap.

Usage:
    python explore.py --controls-dir path/to/controls
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_GENES_EXPECTED = 18_533
N_CTRL_CELLS_EXPECTED = 18_400
CELLS_PER_PERT = 400
STORED_ENTRY_CAP = 4_750_000_000


def load_gene_names(path: Path) -> list:
    """Read the ordered gene symbol list, tolerating an optional header row."""
    first = pd.read_csv(path, header=None, nrows=1).iloc[0, 0]
    header = 0 if str(first).strip().lower() in {"gene", "gene_name", "gene_names"} else None
    return pd.read_csv(path, header=header).iloc[:, 0].astype(str).tolist()


def summarize_context(path: Path) -> dict:
    """Print a QC summary of one context's control file.

    Covers shape, obs columns, count dtype/integrality, UMIs and genes
    detected per cell, and overall density. Returns the per-context stats
    (median/mean genes detected, gene index) that main() needs for the gene
    order check and the projected submission size.
    """
    adata = ad.read_h5ad(path)
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr()

    umis = np.asarray(X.sum(axis=1)).ravel()
    genes_detected = np.diff(X.indptr)

    print(f"\n=== {path.name} ===")
    print(f"  shape                {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print(f"  obs columns          {list(adata.obs.columns)}")
    print(f"  dtype of X           {X.dtype}")
    print(f"  integer counts       {bool(np.all(X.data == np.floor(X.data)))}")
    print(f"  UMIs per cell        median {np.median(umis):,.0f}  "
          f"p5 {np.percentile(umis, 5):,.0f}  p95 {np.percentile(umis, 95):,.0f}")
    print(f"  genes detected/cell  median {np.median(genes_detected):,.0f}  "
          f"p95 {np.percentile(genes_detected, 95):,.0f}  max {genes_detected.max():,}")
    print(f"  overall density       {X.nnz / (X.shape[0] * X.shape[1]):.4f}")

    if "target_gene" in adata.obs:
        vals = adata.obs["target_gene"].unique()
        print(f"  target_gene values   {list(vals)[:5]}{' ...' if len(vals) > 5 else ''}")
    if "ntc_id" in adata.obs:
        counts = adata.obs["ntc_id"].value_counts()
        print(f"  ntc guides           {len(counts)} guides, "
              f"{counts.min()}-{counts.max()} cells each")

    if adata.shape[0] != N_CTRL_CELLS_EXPECTED:
        print(f"  WARNING expected {N_CTRL_CELLS_EXPECTED:,} control cells")
    if adata.shape[1] != N_GENES_EXPECTED:
        print(f"  WARNING expected {N_GENES_EXPECTED:,} genes")

    return {
        "name": path.stem,
        "median_genes_detected": float(np.median(genes_detected)),
        "mean_genes_detected": float(genes_detected.mean()),
        "gene_index": adata.var_names,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls-dir", type=Path, default=Path("controls"))
    args = ap.parse_args()

    d = args.controls_dir

    gene_names = load_gene_names(d / "gene_names.csv")
    perts = pd.read_csv(d / "pert_counts.csv")["target_gene"].astype(str).tolist()

    print(f"gene_names.csv       {len(gene_names):,} genes")
    print(f"pert_counts.csv      {len(perts):,} perturbations")
    print(f"first 5 genes        {gene_names[:5]}")
    print(f"first 5 perts        {perts[:5]}")
    print(f"duplicate perts      {len(perts) - len(set(perts))}")

    overlap = len(set(perts) & set(gene_names))
    print(f"perts also in gene set {overlap}/{len(perts)}")

    summaries = []
    for ctx in ("A", "B", "C"):
        p = d / f"context_{ctx}.h5ad"
        if p.exists():
            summaries.append(summarize_context(p))
        else:
            print(f"\nmissing {p}")

    # Gene order must match gene_names.csv exactly in the submission.
    for s in summaries:
        same = list(s["gene_index"]) == list(gene_names)
        print(f"\n{s['name']} var order matches gene_names.csv: {same}")
        if not same:
            print("  reindex to gene_names.csv order before submitting")

    if summaries:
        mean_detected = float(np.mean([s["mean_genes_detected"] for s in summaries]))
        n_rows = len(perts) * CELLS_PER_PERT * len(summaries)
        projected = n_rows * mean_detected
        print("\n=== projected submission size ===")
        print(f"  rows                 {n_rows:,}")
        print(f"  mean nnz per cell    {mean_detected:,.0f}")
        print(f"  projected stored     {projected:,.0f}")
        print(f"  cap                  {STORED_ENTRY_CAP:,}")
        print(f"  headroom             {STORED_ENTRY_CAP / projected:.2f}x"
              if projected else "")
        print(f"  CSR memory if held   ~{projected * 8 / 1e9:.1f} GB "
              "(float32 data + int32 indices)")
        print("\n  A dense prediction would store "
              f"{n_rows * N_GENES_EXPECTED:,} entries, "
              f"{n_rows * N_GENES_EXPECTED / STORED_ENTRY_CAP:.1f}x over the cap.")


if __name__ == "__main__":
    main()
