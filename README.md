# Arc Virtual Cell Challenge 2026

Predicting single-cell expression responses to CRISPRi knockdowns in cell
contexts with no perturbation training data. Beyond the leaderboard, the
question I'm interested in: how much do predictions degrade when a model has
never seen the target cell type perturbed, and do the challenge metrics
detect that failure?

## Setup

Requires [uv](https://docs.astral.sh/uv/).

```bash
uv sync
uv run vcc login --token-stdin   # token from the challenge site's Credentials page
```

## Data

Keep data outside the repo (it's large and gitignored anyway):

```bash
uv run vcc datasets list
uv run vcc datasets download <validation-bundle-id> -o /path/to/vcc_data
```

Expected layout:

```
vcc_data/
  context_A.h5ad    # unperturbed control cells, contexts anonymized
  context_B.h5ad
  context_C.h5ad
  gene_names.csv    # 18,533 genes, defines submission column order
  pert_counts.csv   # 300 target genes per context
```

## Workflow

```bash
# 1. sanity-check the bundle and the sparsity budget
uv run python explore.py --controls-dir /path/to/vcc_data

# 2. smoke-test the submission pipeline (5 perturbations, not a valid submission)
uv run python baseline_submit.py --controls-dir /path/to/vcc_data --n-perts 5 --out prediction_smoke.h5ad
uv run vcc prep prediction_smoke.h5ad -g /path/to/vcc_data/gene_names.csv \
    --perts /path/to/vcc_data/pert_counts.csv -o smoke.vcc --dry-run

# 3. full baseline
uv run python baseline_submit.py --controls-dir /path/to/vcc_data --out prediction.h5ad
uv run vcc prep prediction.h5ad -g /path/to/vcc_data/gene_names.csv \
    --perts /path/to/vcc_data/pert_counts.csv -o prediction.vcc
uv run vcc submit prediction.vcc
```

The dry-run in step 2 will complain about the perturbation count (expected);
what matters is that gene order, dtypes, context labels, and the raw-counts
check pass.

## Models

Every model implements the `build_chunk` interface from `baseline_submit.py`,
so the submission harness never changes.

- [x] **Predict no change** — resample control cells, relabel. The floor every
  model has to beat.
- [ ] **Global mean effect** — average log-fold-change per target gene learned
  from public Perturb-seq data, pooled across cell types, applied to each
  context's controls. Context-blind by design.
- [ ] **Nearest-context transfer** — match each anonymized context to a public
  cell line by pseudobulk similarity, then transfer that line's measured
  effects.
- [ ] **STATE**, trained with matched zero-shot (held-out cell type) and
  few-shot (held-out perturbations) splits to measure the transfer gap.

## Training data

| Source | Contexts | Size |
|---|---|---|
| VCC 2025 (H1 hESC) | 1 | via Arc Virtual Cell Atlas |
| Replogle 2022 | K562, RPE1 | 9.9 + 8.1 GB |
| Nadig 2025 (GSE264667) | HepG2, Jurkat | 5.2 + 8.7 GB |
| Jiang 2025 (Zenodo 14518762) | 6 cancer lines | — |

## Status

Environment and submission tooling verified; scripts not yet run against the
real validation bundle. Constraints they enforce: exactly 400 cells per
perturbation per context, raw integer counts, sparse storage under the
4.75e9 stored-entry cap, no control cells in the submission.
