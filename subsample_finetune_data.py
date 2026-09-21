"""Subsample the preprocessed gwps fine-tune data so it fits the OS page
cache (random-read training on the full 2M-cell file is disk-bound at
~19 s/step; a ~4 GB subsample is RAM-speed after one pass).

Per perturbation keeps up to --per-pert cells (up to --challenge-per-pert
for perturbations on the challenge panel, which the fine-tune must learn
best), and up to --controls non-targeting cells. Seeded.

Usage:
    python subsample_finetune_data.py \
        --src C:/Users/aaron/vcc_data/state_finetune/replogle_gwps/K562gw.h5ad \
        --out C:/Users/aaron/vcc_data/state_finetune_small/replogle_gwps/K562gw.h5ad
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ad.settings.allow_write_nullable_strings = False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--per-pert", type=int, default=48)
    ap.add_argument("--challenge-per-pert", type=int, default=192)
    ap.add_argument("--controls", type=int, default=25_000)
    ap.add_argument("--pert-counts", type=Path,
                    default=Path("C:/Users/aaron/vcc_data/pert_counts.csv"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    challenge = set(pd.read_csv(args.pert_counts)["target_gene"].astype(str))
    a = ad.read_h5ad(args.src, backed="r")
    genes = a.obs["gene"].astype(str).values
    rng = np.random.default_rng(args.seed)

    keep: list[np.ndarray] = []
    for pert, idx in pd.Series(range(len(genes))).groupby(genes).groups.items():
        idx = np.asarray(idx)
        cap = (args.controls if pert == "non-targeting"
               else args.challenge_per_pert if pert in challenge
               else args.per_pert)
        if len(idx) > cap:
            idx = rng.choice(idx, size=cap, replace=False)
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    print(f"keeping {len(keep):,} of {len(genes):,} cells")

    sub = a[keep].to_memory()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sub.write_h5ad(args.out)
    b = ad.read_h5ad(args.out, backed="r")
    print(f"wrote {args.out}: {b.shape[0]:,} x {b.shape[1]:,}, "
          f"{b.obs['gene'].nunique():,} perts, "
          f"{int((b.obs['gene'] == 'non-targeting').sum()):,} controls")


if __name__ == "__main__":
    main()
