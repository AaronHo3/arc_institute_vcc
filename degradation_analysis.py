"""Fewshot vs zeroshot degradation, from Arc's published ST-HVG-Replogle runs.

For each held-out cell line, Arc trained STATE twice (fewshot: line in
training with a subset of perturbations held out; zeroshot: line excluded
entirely) and published per-perturbation evaluation results. This pairs the
two runs on the perturbations they both evaluated and reports, per metric,
how much performance degrades when the cell type is unseen.

Usage:
    python degradation_analysis.py \
        --evals-dir C:/Users/aaron/vcc_data/st_hvg_replogle_evals
"""

import argparse
from pathlib import Path

import pandas as pd

LINES = ["k562", "rpe1", "hepg2", "jurkat"]
METRICS = [
    "pearson_delta",
    "de_direction_match",
    "de_spearman_sig",
    "overlap_at_100",
    "de_sig_genes_recall",
    "discrimination_score_l1",
    "mae",
]


EFFECT_BINS = [0, 10, 100, 1000, float("inf")]
EFFECT_LABELS = ["<10", "10-100", "100-1000", ">1000"]


def stratify(evals_dir: Path) -> None:
    """Degradation by true effect size (n significant DE genes in the real
    data, from the zeroshot run's DE so binning is identical across runs)."""
    frames = []
    for line in LINES:
        fs = pd.read_csv(evals_dir / "fewshot" / line / f"{line}_results.csv")
        zs = pd.read_csv(evals_dir / "zeroshot" / line / f"{line}_results.csv")
        key = fs.columns[0]
        shared = sorted(set(fs[key]) & set(zs[key]))
        fs = fs.set_index(key).loc[shared]
        zs = zs.set_index(key).loc[shared]
        d = pd.DataFrame({
            "line": line,
            "effect": pd.cut(zs["de_nsig_counts_real"], EFFECT_BINS,
                             labels=EFFECT_LABELS, right=False),
            "fs_pearson_delta": fs["pearson_delta"],
            "zs_pearson_delta": zs["pearson_delta"],
            "fs_direction": fs["de_direction_match"],
            "zs_direction": zs["de_direction_match"],
            "fs_overlap100": fs["overlap_at_100"],
            "zs_overlap100": zs["overlap_at_100"],
        })
        frames.append(d)
    d = pd.concat(frames)

    print("=== stratified by true effect size (n significant DE genes) ===")
    print("perturbations per bin:",
          d.groupby("effect", observed=True).size().to_dict(), "\n")
    g = d.groupby("effect", observed=True).mean(numeric_only=True)
    for name, fs_c, zs_c in [("pearson_delta", "fs_pearson_delta", "zs_pearson_delta"),
                             ("direction_match", "fs_direction", "zs_direction"),
                             ("overlap@100", "fs_overlap100", "zs_overlap100")]:
        out = pd.DataFrame({
            "fewshot": g[fs_c], "zeroshot": g[zs_c],
            "degradation": g[fs_c] - g[zs_c],
            "zs_retained_frac": g[zs_c] / g[fs_c],
        })
        print(f"--- {name} ---")
        print(out.round(3).to_string(), "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evals-dir", type=Path,
                    default=Path("C:/Users/aaron/vcc_data/st_hvg_replogle_evals"))
    ap.add_argument("--stratify", action="store_true",
                    help="also break degradation down by true effect size")
    args = ap.parse_args()

    if args.stratify:
        stratify(args.evals_dir)
        return

    rows = []
    for line in LINES:
        fs = pd.read_csv(args.evals_dir / "fewshot" / line / f"{line}_results.csv")
        zs = pd.read_csv(args.evals_dir / "zeroshot" / line / f"{line}_results.csv")
        key = fs.columns[0]
        shared = sorted(set(fs[key]) & set(zs[key]))
        fs = fs.set_index(key).loc[shared]
        zs = zs.set_index(key).loc[shared]
        for m in METRICS:
            if m not in fs.columns:
                continue
            rows.append({
                "line": line, "metric": m, "n_perts": len(shared),
                "fewshot": fs[m].mean(), "zeroshot": zs[m].mean(),
                "degradation": fs[m].mean() - zs[m].mean(),
            })

    df = pd.DataFrame(rows)
    print(f"perturbations shared per line: "
          f"{df.groupby('line')['n_perts'].first().to_dict()}\n")
    for m in METRICS:
        sub = df[df.metric == m]
        if sub.empty:
            continue
        print(f"=== {m} ===")
        print(sub.pivot_table(index="line", values=["fewshot", "zeroshot", "degradation"])
                 .reindex(LINES).round(3).to_string())
        print()


if __name__ == "__main__":
    main()
