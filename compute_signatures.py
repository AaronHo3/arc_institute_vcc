"""Compute per-perturbation effect signatures from Replogle 2022 K562
genome-wide pseudobulk, restricted to the challenge gene panel.

For each challenge perturbation covered by the screen, the signature is the
log2 fold change of every measured gene vs the mean of the non-targeting
control pseudobulks. Output is a small AnnData (perturbations x genes,
gene symbols as var_names) consumed by the global-mean baseline.

Usage:
    python compute_signatures.py \
        --bulk C:/Users/aaron/vcc_data/replogle/K562_gwps_raw_bulk.h5ad \
        --controls-dir controls \
        --out C:/Users/aaron/vcc_data/signatures_k562gwps.h5ad
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Pseudocount for the fold-change ratio. X is mean UMI per cell (typical
# range 0.01-10), so 0.1 damps ratios for barely-expressed genes without
# swamping real signal.
PSEUDOCOUNT = 0.1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", type=Path, required=True)
    ap.add_argument("--controls-dir", type=Path, default=Path("controls"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    perts = pd.read_csv(args.controls_dir / "pert_counts.csv")["target_gene"].astype(str).unique()

    a = ad.read_h5ad(args.bulk)
    X = np.asarray(a.X)

    # obs_names look like '10023_ZC3H18_P1P2_ENSG00000158545'
    symbols = np.array([n.split("_")[1] for n in a.obs_names])

    ctrl = X[a.obs["core_control"].values].mean(axis=0)
    print(f"{a.shape[0]:,} pseudobulks, {int(a.obs['core_control'].sum())} controls, "
          f"ctrl mean expr median {np.median(ctrl):.3f}")

    # One row per covered challenge perturbation (average if multiple).
    rows, kept = [], []
    for p in perts:
        mask = symbols == p
        if mask.any():
            rows.append(X[mask].mean(axis=0))
            kept.append(p)
    print(f"covered perturbations: {len(kept)}/{len(perts)}")

    lfc = np.log2((np.vstack(rows) + PSEUDOCOUNT) / (ctrl + PSEUDOCOUNT))

    # Reindex columns to gene symbols; drop genes without a symbol match in
    # the challenge panel. Duplicated symbols keep the higher-expressed row.
    gene_syms = a.var["gene_name"].astype(str).values
    challenge_genes = set(
        pd.read_csv(args.controls_dir / "gene_names.csv").iloc[:, 0].astype(str)
    )
    order = np.argsort(-ctrl)  # prefer higher-expressed duplicate
    seen, keep_idx = set(), []
    for i in order:
        s = gene_syms[i]
        if s in challenge_genes and s not in seen:
            seen.add(s)
            keep_idx.append(i)
    keep_idx = np.array(sorted(keep_idx))
    print(f"genes mapped to challenge panel: {len(keep_idx):,}/{a.shape[1]:,}")

    sig = ad.AnnData(
        X=lfc[:, keep_idx].astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(kept, name="target_gene")),
        var=pd.DataFrame(index=pd.Index(gene_syms[keep_idx], name="gene")),
    )
    sig.write_h5ad(args.out)
    print(f"wrote {args.out}")

    # Sanity: each perturbation's own gene should be strongly down (CRISPRi).
    own = []
    var_pos = {g: i for i, g in enumerate(sig.var_names)}
    for i, p in enumerate(kept):
        if p in var_pos:
            own.append(sig.X[i, var_pos[p]])
    own = np.array(own)
    print(f"\ntarget gene's own LFC: median {np.median(own):.2f}, "
          f"p90 {np.percentile(own, 90):.2f} "
          f"({(own < -1).mean():.0%} below -1, i.e. <50% residual)")
    big = (np.abs(sig.X) > 1).sum(axis=1)
    print(f"genes with |LFC|>1 per pert: median {int(np.median(big))}, "
          f"p90 {int(np.percentile(big, 90))}, max {big.max()}")


if __name__ == "__main__":
    main()
