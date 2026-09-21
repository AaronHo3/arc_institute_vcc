# When do perturbation-response predictions transfer to a new cell context — and do the metrics tell you when they've failed?

> **STATUS: skeleton.** Structure, claims, numbers, and figure references are
> final; prose is to be written. Each section
> lists its claims with evidence pointers. Nothing goes in the final text
> that isn't backed by a pointer here.

## Abstract (5-6 sentences, write last)

- The question; the VCC 2026 setting (zero-shot contexts).
- Method: baseline ladder on the live leaderboard + paired fewshot/zeroshot
  analysis of Arc's published STATE runs + effect-size stratification.
- Finding 1: transfer is real but partial, and concentrated in
  large-effect perturbations.
- Finding 2: expression-error metrics are structurally blind to the
  failure that DE-level metrics see clearly.
- One honest limitation sentence (observational reuse of Arc's runs; own
  from-scratch training hit a compute floor).

## 1. Setting and question

- VCC 2026: predict post-CRISPRi expression distributions for 300 genes x 3
  anonymized contexts, no perturbation data from target contexts.
- Why metric validity under context shift is the open problem (cite Arc's
  own framing: "model accuracy, metric design and biological
  generalization"; six-metric aggregate).
- Analogy to external validation of predictive models (held-out site ->
  held-out cell type). One paragraph, no overclaiming.

## 2. The baseline ladder (leaderboard evidence)

Table from RESULTS.md submission log: null / target-kd / global-mean, with
per-metric server-scaled scores and ranks.

Claims:
- C2.1 The server's reference baseline is itself a no-change predictor
  (MSE pinned at 0 for three radically different submissions).
- C2.2 The DE metrics exclude the target gene (fid/jac unchanged to ~5
  decimals between null and target-kd) — measured from outside via a
  two-submission experiment.
- C2.3 A context-blind K562 lookup table (global-mean) recovered 89% of the
  null's fid deficit and reached rank 265, mid-pack among trained models.
- C2.4 Metric disagreement is observable between our own submissions (PDS
  rewarded target-kd; five metrics ignored it).

## 3. What the challenge forecloses (census + context identity)

- Census table: K562-gwps 272/300; H1 13; Jiang 9; all essential-panel
  datasets 0. 28 perts have no public measurement anywhere checked.
- Context de-anonymization: A=Jurkat (r=0.881 vs 1,673 DepMap lines,
  independently confirmed vs Nadig Jurkat controls), B=HeLa (0.875),
  C=CAL-33 (0.856). Method: Spearman on top-2000-variance genes,
  pseudobulk vs bulk.
- Implication: per-context lookup is impossible by design; the global-mean
  result is genuinely cross-lineage (no context resembles K562).

## 4. The core result: zero-shot degradation (Figures 1-2)

Source: Arc's published ST-HVG-Replogle paired runs (their models, our
pairing analysis; ~1,000 shared perts per line, 4 lines).

- Fig 1: per-metric fewshot vs zeroshot dumbbells.
- C4.1 Every effect-level metric degrades consistently across all four
  held-out lines (pearson_delta -0.10..-0.17; overlap@100 halved).
- C4.2 Zero-shot is far from destroyed: direction 0.67-0.70 vs 0.5 chance.
- C4.3 MAE barely moves (0.05->0.07) under the same regime change — the
  raw-error family cannot see zero-shot failure. Ties to C2.1.
- Fig 2: stratified by true effect size.
- C4.4 Below ~10 significant DE genes nothing is predictable in either
  regime (delta-r 0.02, direction at chance; n=203).
- C4.5 Utility begins >~100 sig genes; the zero-shot cost concentrates
  there (retention: 87% direction, ~70% delta-r, ~50% overlap).
- C4.6 Within the reliable regime, the size of the zero-shot penalty
  depends on which metric you ask (13% vs 30% vs 50%).

## 5. Corroboration: cross-line effect conservation without models

- jiang_transfer.py on Jiang 2025 (6 lines, 218 stimulated regulators):
  median pairwise r 0.06 (any-responsive) / 0.22 (both-responsive); 19%
  of pairs r>0.7, 57% r<0.3; ~90% of (pert, pair) combos lack even 20
  shared responsive genes.
- Read: conservation is bimodal and gene-class dependent; consistent with
  C4.4/C4.5 (transfer lives in large, core-program effects).
- Caveats box: stimulated conditions, Mixscale weighting, signaling-gene
  panel.

## 6. Negative result: the compute floor

- From-scratch STATE (132M params) at 16k steps x batch 8 (2-4% of default
  budget, one RTX 4070 Ti) fails to fit its own training data
  (delta-r 0.006 on trained-on perts). Table from RESULTS.md.
- Read honestly: consumer-hardware from-scratch training cannot support
  transfer conclusions; pretrained initialization is a practical necessity.
- This is why Section 4 uses Arc's published runs.

## 7. Synthesis: what the metrics can't see

- Three independent evidence lines (leaderboard MSE, Arc-run MAE, metric
  spread within reliable regime) -> expression-error metrics are
  insensitive to exactly the failure mode the challenge is about.
- What this recommends for metric design (modest, concrete: report
  DE-level metrics stratified by effect size; treat raw-error metrics as
  format checks, not accuracy signals).

## 8. Limitations (bullet-for-bullet, no softening)

- Section 4 reuses Arc's published runs: not an independent replication;
  training choices were theirs.
- One model family (STATE HVG variant); one holdout design (essential
  panel).
- Context identities for B/C are high-confidence but not certain (could be
  close sibling lines).
- Jiang regime differences (see Section 5 caveats).
- Leaderboard findings use validation contexts A-C; final-test
  generalization unknown until Oct 22.

## 9. Reproducibility

- Repo map: explore.py, baseline_submit.py (3 models), compute_signatures.py,
  identify_contexts.py, jiang_transfer.py, degradation_analysis.py,
  make_figures.py; configs/; RESULTS.md as the running lab log.
- One command per figure; pinned uv environment; fixed seeds.
- Data provenance table: challenge bundle, Replogle (Figshare 20029387),
  Nadig (GSE264667), Jiang (Zenodo 14518762), DepMap 24Q4, ST-HVG-Replogle
  (HuggingFace, CC BY-NC-SA — attribution).

## Acknowledgments / attribution

- Arc Institute for STATE, cell-eval/cell-load, and published checkpoints
  with evaluation outputs; challenge organizers; dataset authors (Replogle,
  Nadig, Jiang, DepMap).
