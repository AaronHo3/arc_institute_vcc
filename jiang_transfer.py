"""Cross-cell-line conservation of perturbation effects in Jiang 2025.

The Jiang DE results (Zenodo 14518762) give per-cell-line log2FC for 218
pathway regulators across six cancer lines under cytokine stimulation. For
each (regulator, line-pair), correlates the log2FC vectors between the two
lines over responsive genes, under two filters that bracket the selection
effect:

  any-responsive:  genes with max |log2FC| > 0.5 in ANY line (deflates r via
                   line-specific responders)
  both-responsive: genes with |log2FC| > 0.5 in BOTH lines of the pair
                   (conditions on shared response)

Caveats: cells were cytokine-stimulated (a different regime from resting
CRISPRi screens), values come through Mixscale weighted DE, and regulators
are signaling genes — the gene class most expected to be context-dependent.

Usage:
    python jiang_transfer.py --zip C:/Users/aaron/vcc_data/jiang/DE_results_all_pathway.zip
"""

import argparse
import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
LFC_THRESHOLD = 0.5
MIN_GENES = 20


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", type=Path,
                    default=Path("C:/Users/aaron/vcc_data/jiang/DE_results_all_pathway.zip"))
    args = ap.parse_args()

    z = zipfile.ZipFile(args.zip)
    files = [n for n in z.namelist()
             if re.match(r"DE_results_all_pathway/Parse_\w+/[A-Za-z0-9-]+_\w+_pathway_DE_results\.txt$", n)]
    cols = [f"log2FC_{l}" for l in LINES]

    any_rs, both_rs = [], []
    for n in files:
        df = pd.read_csv(io.BytesIO(z.read(n)), sep=r"\s+")
        if not all(c in df.columns for c in cols):
            continue
        lfc = df[cols].to_numpy(dtype=float)
        lfc = lfc[np.isfinite(lfc).all(axis=1)]

        m = np.abs(lfc).max(axis=1) > LFC_THRESHOLD
        if m.sum() >= MIN_GENES:
            c = np.corrcoef(lfc[m].T)
            iu = np.triu_indices(len(LINES), k=1)
            any_rs.extend(c[iu])

        for i in range(len(LINES)):
            for j in range(i + 1, len(LINES)):
                m2 = (np.abs(lfc[:, i]) > LFC_THRESHOLD) & (np.abs(lfc[:, j]) > LFC_THRESHOLD)
                if m2.sum() >= MIN_GENES:
                    both_rs.append(np.corrcoef(lfc[m2, i], lfc[m2, j])[0, 1])

    for name, rs in [("any-responsive", np.array(any_rs)),
                     ("both-responsive", np.array(both_rs))]:
        print(f"=== {name} filter ({len(rs)} (pert, line-pair) combos) ===")
        print(f"  median r {np.median(rs):.2f} | "
              f"IQR {np.percentile(rs, 25):.2f}..{np.percentile(rs, 75):.2f} | "
              f"r>0.7: {(rs > 0.7).mean():.0%} | r<0.3: {(rs < 0.3).mean():.0%}")


if __name__ == "__main__":
    main()
