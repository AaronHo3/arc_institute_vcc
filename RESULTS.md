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
