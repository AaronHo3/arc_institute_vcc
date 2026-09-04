"""Which known cell lines do the anonymous contexts A/B/C resemble?

Computes a control-cell pseudobulk (mean expression per gene) for each
challenge context and for each public reference line, normalizes to log1p
CPM on the shared gene set, and reports pairwise correlations. Contexts and
references come from different platforms (10x Flex vs 10x v3 vs others), so
only the ranking within a context is meaningful, not absolute values.

Usage:
    python identify_contexts.py --controls-dir controls --data-dir C:/Users/aaron/vcc_data
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr


def log_cpm(mean_counts: np.ndarray) -> np.ndarray:
    return np.log1p(mean_counts / mean_counts.sum() * 1e6)


def context_pseudobulk(path: Path) -> pd.Series:
    a = ad.read_h5ad(path)
    m = np.asarray(a.X.mean(axis=0)).ravel()
    return pd.Series(m, index=a.var_names.astype(str))


def replogle_bulk_pseudobulk(path: Path) -> pd.Series:
    """Mean of the non-targeting control pseudobulks, indexed by gene symbol."""
    a = ad.read_h5ad(path)
    m = np.asarray(a.X[a.obs["core_control"].values].mean(axis=0)).ravel()
    s = pd.Series(m, index=a.var["gene_name"].astype(str))
    return s.groupby(level=0).mean()  # collapse duplicate symbols


def nadig_pseudobulk(path: Path) -> pd.Series:
    """Mean over non-targeting cells, from the raw single-cell file."""
    a = ad.read_h5ad(path, backed="r")
    mask = (a.obs["gene"].astype(str) == "non-targeting").values
    sub = a[mask].to_memory()
    X = sub.X if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    m = np.asarray(X.mean(axis=0)).ravel()
    if "gene_name" in sub.var.columns:
        idx = sub.var["gene_name"].astype(str)
    else:
        idx = sub.var_names.astype(str)
    return pd.Series(m, index=idx).groupby(level=0).mean()


def depmap_matches(contexts: dict, depmap_dir: Path, n_hvg: int = 2000) -> None:
    """Rank all DepMap lines by similarity to each context's pseudobulk.

    DepMap values are log2(TPM+1) bulk RNA-seq; contexts are 10x Flex
    single-cell counts. TPM is transcript-length normalized and counts are
    not, so Spearman (rank) correlation on the top-variance genes is used;
    absolute values still aren't comparable across chemistries — only the
    ranking matters.
    """
    expr = pd.read_csv(depmap_dir / "expression_tpm_logp1.csv", index_col=0)
    expr.columns = [c.split(" ")[0] for c in expr.columns]
    expr = expr.loc[:, ~expr.columns.duplicated()]
    model = pd.read_csv(depmap_dir / "Model.csv").set_index("ModelID")
    print(f"\nDepMap: {expr.shape[0]:,} lines x {expr.shape[1]:,} genes")

    shared = sorted(set(contexts["A"].index) & set(expr.columns))
    expr = expr[shared]
    hvg = expr.var(axis=0).nlargest(n_hvg).index
    E = expr[hvg].to_numpy()
    from scipy.stats import rankdata
    ER = rankdata(E, axis=1)
    ER = (ER - ER.mean(axis=1, keepdims=True)) / ER.std(axis=1, keepdims=True)
    print(f"shared genes {len(shared):,}; matching on top {n_hvg} variable")

    for cname, c in contexts.items():
        cv = log_cpm(c.reindex(shared).to_numpy())
        cv = pd.Series(cv, index=shared).reindex(hvg).to_numpy()
        cr = rankdata(cv)
        cr = (cr - cr.mean()) / cr.std()
        r = ER @ cr / len(cr)
        order = np.argsort(-r)
        print(f"\ncontext {cname} — top matches:")
        for k in order[:8]:
            mid = expr.index[k]
            name = model.loc[mid, "StrippedCellLineName"] if mid in model.index else "?"
            lineage = model.loc[mid, "OncotreeLineage"] if mid in model.index else "?"
            print(f"  {r[k]:.3f}  {name:<14} {lineage}")
        for known in ("K562", "RPE1SS48", "HEPG2", "JURKAT"):
            hits = model.index[model["StrippedCellLineName"].str.upper() == known]
            hits = [m for m in hits if m in expr.index]
            if hits:
                k = expr.index.get_loc(hits[0])
                rank = int((r > r[k]).sum()) + 1
                print(f"  [{known}: rank {rank}/{len(r)}, r={r[k]:.3f}]")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls-dir", type=Path, default=Path("controls"))
    ap.add_argument("--data-dir", type=Path, default=Path("C:/Users/aaron/vcc_data"))
    ap.add_argument("--depmap", action="store_true",
                    help="also match contexts against all DepMap lines")
    args = ap.parse_args()

    if args.depmap:
        contexts = {c: context_pseudobulk(args.controls_dir / f"context_{c}.h5ad")
                    for c in ("A", "B", "C")}
        depmap_matches(contexts, args.data_dir / "depmap")
        return

    contexts = {c: context_pseudobulk(args.controls_dir / f"context_{c}.h5ad")
                for c in ("A", "B", "C")}

    ref_sources = {
        "K562": lambda: replogle_bulk_pseudobulk(args.data_dir / "replogle" / "K562_gwps_raw_bulk.h5ad"),
        "RPE1": lambda: replogle_bulk_pseudobulk(args.data_dir / "replogle" / "rpe1_raw_bulk.h5ad"),
        "HepG2": lambda: nadig_pseudobulk(args.data_dir / "nadig" / "hepg2_raw_singlecell.h5ad"),
        "Jurkat": lambda: nadig_pseudobulk(args.data_dir / "nadig" / "jurkat_raw_singlecell.h5ad"),
    }
    refs = {}
    for name, load in ref_sources.items():
        try:
            refs[name] = load()
            print(f"loaded {name}: {len(refs[name]):,} genes")
        except Exception as e:
            print(f"skipping {name}: {e}")

    shared = set(contexts["A"].index)
    for r in refs.values():
        shared &= set(r.index)
    shared = sorted(shared)
    print(f"\nshared genes across all sources: {len(shared):,}")

    rows = []
    for cname, c in contexts.items():
        cv = log_cpm(c.reindex(shared).to_numpy())
        for rname, r in refs.items():
            rv = log_cpm(r.reindex(shared).to_numpy())
            rows.append({
                "context": cname, "reference": rname,
                "pearson": np.corrcoef(cv, rv)[0, 1],
                "spearman": spearmanr(cv, rv).statistic,
            })
    df = pd.DataFrame(rows)

    for metric in ("pearson", "spearman"):
        print(f"\n{metric} (rows: contexts, cols: references):")
        print(df.pivot(index="context", columns="reference", values=metric)
                .round(3).to_string())
    print("\nbest match per context (spearman):")
    for cname in contexts:
        sub = df[df.context == cname].sort_values("spearman", ascending=False)
        top, second = sub.iloc[0], sub.iloc[1]
        print(f"  {cname}: {top.reference} ({top.spearman:.3f}), "
              f"runner-up {second.reference} ({second.spearman:.3f})")


if __name__ == "__main__":
    main()
