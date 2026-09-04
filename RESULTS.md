# Submission log

Validation partition (contexts A/B/C, panel vcc2026-val-1). Scores are the
server's baseline-scaled values: 0 = the server's internal reference baseline,
positive = better than it. Overall is the unweighted mean of the six metrics.

| Date | Model | Entry | Rank | Overall | PDS | MSE | NMAE | FID | Reach | JAC |
|---|---|---|---|---|---|---|---|---|---|---|
| 2026-09-03 | predict-no-change | YBSJ2Ua87fYVq24VhXHb | 568 | -0.3039 | -0.0088 | 0.0000 | -0.0053 | -1.7158 | -0.0111 | -0.0826 |
| 2026-09-04 | target-knockdown | 8SnxzMGp74VQRoU7cJzk | 552 | -0.2987 | +0.0215 | 0.0000 | -0.0064 | -1.7161 | -0.0087 | -0.0826 |

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
