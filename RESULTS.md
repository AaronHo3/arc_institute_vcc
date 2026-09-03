# Submission log

Validation partition (contexts A/B/C, panel vcc2026-val-1). Scores are the
server's baseline-scaled values: 0 = the server's internal reference baseline,
positive = better than it. Overall is the unweighted mean of the six metrics.

| Date | Model | Entry | Rank | Overall | PDS | MSE | NMAE | FID | Reach | JAC |
|---|---|---|---|---|---|---|---|---|---|---|
| 2026-09-03 | predict-no-change | YBSJ2Ua87fYVq24VhXHb | 568 | -0.3039 | -0.0088 | 0.0000 | -0.0053 | -1.7158 | -0.0111 | -0.0826 |

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
