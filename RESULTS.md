# Submission log

Validation partition (contexts A/B/C, panel vcc2026-val-1). Scores are the
server's baseline-scaled values: 0 = the server's internal reference baseline,
positive = better than it. Overall is the unweighted mean of the six metrics.

| Date | Model | Entry | Rank | Overall | PDS | MSE | NMAE | FID | Reach | JAC |
|---|---|---|---|---|---|---|---|---|---|---|
| 2026-09-03 | predict-no-change | YBSJ2Ua87fYVq24VhXHb | 568 | -0.3039 | -0.0088 | 0.0000 | -0.0053 | -1.7158 | -0.0111 | -0.0826 |
| 2026-09-04 | target-knockdown | 8SnxzMGp74VQRoU7cJzk | 552 | -0.2987 | +0.0215 | 0.0000 | -0.0064 | -1.7161 | -0.0087 | -0.0826 |
| 2026-09-04 | global-mean-effect | DWEB4FyVrq7t7NLJh6ko | 265 | +0.0575 | +0.3990 | 0.0000 | +0.0596 | -0.1927 | +0.0980 | -0.0191 |

## Notes

**2026-09-03, predict-no-change.** Resampled non-targeting control cells per
context, relabeled per perturbation. No model. Reads on the six metrics:

- Sits at ~0 (the server's reference) on five of six metrics, which suggests
  the server's own scaling baseline is itself a no-change-style predictor;
  MSE = 0 exactly.
- DE direction fidelity (fid = -1.72) accounts for ~94% of the overall
  penalty. Predicting "no change" is punished far harder for asserting no
  direction on truly-moved genes than for anything else it gets wrong.
- Perturbation discrimination is at chance, as expected: all 300 predictions
  per context are draws from the same pool.
- Any future model must clear these numbers per metric, not just overall.

**2026-09-04, target-knockdown.** Same as predict-no-change, plus the target
gene's own counts binomially thinned to 15% in every cell (CRISPRi efficacy).
Differs from the null by exactly one directionally-correct gene per
perturbation. Reads:

- fid unchanged (-1.7161 vs -1.7158) and jac unchanged to 5 decimals. The DE
  metrics appear to exclude the target gene itself from evaluation (standard
  practice, since the target's drop is trivially known) — to be confirmed
  against cell-eval source. Consequence: the only route to improving fid is
  predicting downstream gene directions, i.e. actual perturbation knowledge.
- PDS is the one metric that rewarded it: -0.009 -> +0.022. Thinning one gene
  makes each prediction identifiable among the 300.
- nmae/reach shifts (~±0.003) are likely run-to-run resampling noise (the RNG
  stream differs from the null run), giving a rough noise floor for reading
  future score differences.

**2026-09-04, global-mean-effect.** Replogle K562 genome-wide LFC signatures
applied multiplicatively to resampled controls, identically in all three
contexts (context-blind); 272/300 perts covered, rest fell back to
target-knockdown. Rank 552 -> 265. Reads:

- fid recovered from -1.716 to -0.193: ~89% of the null's direction-fidelity
  gap closed by borrowing one cell line's measured effects with zero context
  modeling. Cross-context conservation of perturbation direction is
  substantial in these contexts.
- pds jumped to +0.40 — borrowed signatures make predictions strongly
  identifiable.
- jac still negative (-0.019): predicting *which genes reach significance*
  remains below reference even when directions are right. Signatures are
  sparse (median 2 genes |LFC|>1), so predicted DE sets are thin.
- mse = 0 for the third consecutive submission despite three very different
  prediction sets. The metric appears saturated/clamped at the reference
  value — a metric-behavior question to investigate (cell-eval source).
- Aggregate crossed zero: a lookup table with no training now beats the
  server reference. The bar for STATE is now this, not the null.

## STATE from-scratch training: negative result (2026-09-06)

Trained STATE (132M params) from scratch on the harmonized 5-context corpus,
paired fewshot/zeroshot Jurkat splits, 16,000 steps x batch 8 (~2-4% of the
400k-step default; chosen to fit one RTX 4070 Ti overnight). Both models
evaluated on the shared 200 held-out (Jurkat, pert) pairs, and the fewshot
model additionally on 100 pairs it trained on:

| eval set | pearson_delta | DE direction | pred sig genes | discrim (chance=0.5) |
|---|---|---|---|---|
| fewshot, held-out 200 | 0.006 | 0.51 | 0.005 (real: 218) | 0.54 |
| zeroshot, held-out 200 | 0.003 | 0.50 | 0.5 (real: 218) | 0.52 |
| fewshot, trained-on 100 | 0.006 | 0.50 | 0.0 (real: 211) | 0.53 |

The model fails identically on data it trained on, so this is a training
failure (insufficient steps and/or missing the recommended preprocessing
recipe), not a transfer measurement. No fewshot-vs-zeroshot conclusion can
be drawn from these runs. Next: fine-tune from Arc's pretrained STATE
weights (init_from) instead of training from scratch on consumer hardware.

## Zeroshot-vs-fewshot degradation (2026-09-06, from Arc's published runs)

Arc's ST-HVG-Replogle release (HuggingFace) contains their preprint's paired
fewshot/zeroshot STATE runs for all four essential-panel lines, with
per-perturbation evaluation CSVs. `degradation_analysis.py` pairs each
line's two runs on their shared ~1,000 perturbations. Mean degradation
(fewshot minus zeroshot), by metric, across held-out lines:

| metric | fewshot | zeroshot | degradation |
|---|---|---|---|
| pearson_delta | 0.36-0.49 | 0.26-0.32 | **-0.10 to -0.17** |
| DE direction match | 0.75-0.76 | 0.67-0.70 | -0.06 to -0.09 |
| DE overlap@100 | 0.18-0.23 | 0.08-0.11 | **roughly halved** |
| DE sig-gene recall | 0.26-0.28 | 0.14-0.25 | -0.01 to -0.14 |
| discrimination (L1) | 0.73-0.78 | 0.58-0.69 | -0.07 to -0.17 |
| MAE (raw error) | 0.05-0.06 | 0.06-0.08 | ~0 (slightly worse) |

Reads:

- Zero-shot transfer degrades every effect-level metric consistently across
  all four lines, but is far from destroyed: direction match stays at
  0.67-0.70 vs 0.5 chance. Consistent with our leaderboard finding that
  cross-context effect conservation is substantial.
- **The metric-sensitivity split is the thesis result**: raw-expression
  error (MAE) barely moves under the same regime change that halves DE
  overlap. A leaderboard weighted toward expression-error metrics would not
  see zero-shot failure; DE-level metrics see it clearly. Mirrors our
  submission-side finding that mse pinned at 0 across radically different
  predictions.
- Attribution: models and eval outputs are Arc's (ST-HVG-Replogle); the
  pairing analysis and framing are ours. Our own from-scratch runs (above)
  could not measure this due to compute-floor training failure.

### Stratified by true effect size (pooled over 4 lines)

`degradation_analysis.py --stratify`, binning by n significant DE genes in
the real data:

| effect bin | n perts | pearson_delta fs/zs | direction fs/zs | overlap@100 fs/zs |
|---|---|---|---|---|
| <10 sig genes | 203 | 0.02 / 0.02 | 0.53 / 0.52 | 0.01 / 0.01 |
| 10-100 | 1,773 | 0.19 / 0.13 | 0.67 / 0.62 | 0.10 / 0.07 |
| 100-1000 | 2,050 | 0.64 / 0.44 | 0.86 / 0.75 | 0.31 / 0.15 |
| >1000 | 32 | 0.82 / 0.60 | 0.87 / 0.75 | 0.47 / 0.23 |

- **Below ~10 significant genes, nothing is predictable in either regime**
  (delta correlation ~0.02, direction at chance) — small-effect failure is
  fundamental, not a transfer problem.
- Predictive utility starts around >100 significant genes; that is also
  where the zero-shot cost concentrates in absolute terms (pearson_delta
  -0.20, overlap@100 halved).
- Boundary quantified: the "reliable" regime is large-effect perturbations
  (~1/3 of the panel), and zero-shot transfer retains ~68-74% of
  delta-correlation and ~87% of direction accuracy there, but only ~50% of
  DE-set overlap.

## Pipeline verification and fine-tune diagnosis (2026-09-07)

Ran Arc's untouched ST-HVG-Replogle fewshot/jurkat checkpoint through our
full local pipeline (our Nadig conversion, our reverse-engineered
normalization, our cell-load setup, our eval) on 100 of their held-out
Jurkat perturbations: **pearson_delta 0.377 vs their published 0.381**
(direction 0.785 vs 0.762, overlap@100 0.168 vs 0.183). The local pipeline
is verified end to end; all earlier "dead model" results were real model
failures, not evaluation artifacts.

Fine-tune v1/v2 failure diagnosed (both floored even on trained-on perts):

- Exposure starvation: 9,867 freshly-initialized one-hot perturbation
  embeddings x 6k steps = ~19 examples/pert (Arc's runs: ~3,000). The
  model regressed to its residual predict-the-controls path.
- Holdout design flaw: with one-hot pert encoding, a held-out pert's
  embedding is untrained noise — heldout evals were floored by
  construction unless the pert was seen in another training context.

Fix (fine-tune v3, running): restrict to the 272 challenge perturbations +
controls (66k cells), all perts in training (~2,300 exposures each at 20k
steps); validate via reproduction of known K562 effects; context transfer
is then tested on the leaderboard.

Also: the 50-pert challenge-gene eval revealed the challenge panel's
covered genes are mostly small-effect in K562 (mean 6 significant DE genes
vs 200+ for essential-panel perts) — the regime our stratification shows is
hard for any model. Tempers expected leaderboard gains for all approaches.

## Public-data census (2026-09-04)

Coverage of the 300 validation perturbations by public CRISPRi datasets:

| Source | Contexts | Coverage |
|---|---|---|
| Replogle 2022 K562 genome-wide | K562 | **272/300** |
| VCC 2025 training | H1 hESC | 13/300 (all within the 272) |
| Jiang 2025 | 6 cancer lines | 9/300 (all within the 272) |
| Replogle 2022 essential panels | K562, RPE1 | 0/300 |
| Nadig 2025 | HepG2, Jurkat | 0/300 |

28 perturbations have no measured public effect in any checked source. The
panel appears designed so that exactly one dataset (the genome-wide K562
screen) provides broad coverage; per-context lookup ("nearest-context
transfer") is not buildable for this panel.

Context identities (identify_contexts.py --depmap, Spearman on top-2000
variable genes vs DepMap 24Q4): A = Jurkat (r=0.881, rank 1/1673),
B = HeLa (0.875), C = CAL-33 (0.856); A independently confirmed against
Nadig's Jurkat CRISPRi controls. No context resembles K562, so the
global-mean fid recovery reflects genuinely cross-lineage effect
conservation.
