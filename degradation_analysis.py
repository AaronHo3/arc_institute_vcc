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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evals-dir", type=Path,
                    default=Path("C:/Users/aaron/vcc_data/st_hvg_replogle_evals"))
    args = ap.parse_args()

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
