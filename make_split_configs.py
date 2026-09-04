"""Generate the paired cell-load TOML configs for the core experiment.

Both configs train on the same five cell types and are evaluated on the SAME
held-out (Jurkat, perturbation) pairs; they differ in exactly one respect:

  fewshot.toml  — Jurkat is in training, but the held-out perturbations are
                  excluded from it (model has seen Jurkat perturbed).
  zeroshot.toml — Jurkat is excluded from training entirely (model has never
                  seen Jurkat perturbed).

The gap between the two models' scores on the shared held-out set is the
zero-shot transfer penalty.

Usage:
    python make_split_configs.py [--n-test 200] [--n-val 50] [--seed 0]
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np

TRAIN_DIR = "C:/Users/aaron/vcc_data/state_train"


def toml_list(genes: list[str]) -> str:
    inner = ", ".join(f'"{g}"' for g in genes)
    return f"[{inner}]"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-test", type=int, default=200)
    ap.add_argument("--n-val", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", type=Path, default=Path("configs"))
    args = ap.parse_args()

    a = ad.read_h5ad(r"C:\Users\aaron\vcc_data\nadig\jurkat_raw_singlecell.h5ad", backed="r")
    perts = sorted(set(a.obs["gene"].astype(str)) - {"non-targeting"})
    rng = np.random.default_rng(args.seed)
    held = rng.choice(perts, size=args.n_test + args.n_val, replace=False)
    test, val = sorted(held[: args.n_test]), sorted(held[args.n_test:])
    print(f"Jurkat panel {len(perts)} perts -> {len(test)} test + {len(val)} val held out (seed {args.seed})")

    common = f"""[datasets]
replogle = "{TRAIN_DIR}/replogle/"
nadig = "{TRAIN_DIR}/nadig/"
h1 = "{TRAIN_DIR}/h1/"

[training]
replogle = "train"
nadig = "train"
h1 = "train"
"""

    fewshot = common + f"""
[fewshot."nadig.Jurkat"]
val = {toml_list(val)}
test = {toml_list(test)}
"""

    zeroshot = common + """
[zeroshot]
"nadig.Jurkat" = "test"
"""

    args.out_dir.mkdir(exist_ok=True)
    (args.out_dir / "fewshot.toml").write_text(fewshot)
    (args.out_dir / "zeroshot.toml").write_text(zeroshot)
    (args.out_dir / "heldout_perts.txt").write_text("\n".join(test) + "\n")
    print(f"wrote {args.out_dir}/fewshot.toml, zeroshot.toml, heldout_perts.txt")


if __name__ == "__main__":
    main()
