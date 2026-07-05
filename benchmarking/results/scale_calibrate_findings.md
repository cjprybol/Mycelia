# Stage 0 k-mer classification — scale grid + threshold calibration

Two follow-ups to the ground-truth comparison matrix
(`benchmarking/stage0_classification_matrix.jl`):

1. **Scale grid** — run the matrix across seeds {42, 43, 44} × error rates
   {0.005, 0.01, 0.02}, all at `MYCELIA_S0_REFLEN=10000` (cov 30, readlen 100,
   k 21). Nine runs; per-strategy mean ± sd of balanced accuracy and AUC.
2. **Threshold calibration** — for one base config (reflen 10000, seed 42,
   err 0.01), sweep each arm's continuous score over its observed range to find
   the balanced-accuracy-maximizing threshold, and measure the gain over each
   arm's default cutoff (`benchmarking/stage0_threshold_calibration.jl`).

## 1. Scale grid (n = 9: 3 seeds × 3 error rates, reflen 10000)

| strategy               | balanced accuracy | AUC             |
| ---------------------- | ----------------- | --------------- |
| FixedCoverageThreshold | 0.932 ± 0.019     | 0.946 ± 0.025   |
| QualityThreshold       | 0.500 ± 0.000     | 0.946 ± 0.024   |
| BayesianMixture        | **0.945 ± 0.024** | 0.946 ± 0.024   |
| EffectiveCoverage      | 0.939 ± 0.020     | 0.945 ± 0.025   |
| LogisticFusion         | 0.760 ± 0.019     | 0.946 ± 0.024   |

**AUC is effectively identical across all five arms (~0.946)** — the scores rank
genomic vs error k-mers equally well regardless of arm. The spread in *balanced
accuracy* is therefore entirely a function of where each arm's DEFAULT threshold
falls, not of discriminating power.

### Balanced accuracy by error rate (mean over 3 seeds)

| strategy               | err 0.005 | err 0.01 | err 0.02 |
| ---------------------- | --------- | -------- | -------- |
| FixedCoverageThreshold | 0.908     | 0.938    | 0.950    |
| QualityThreshold       | 0.500     | 0.500    | 0.500    |
| BayesianMixture        | 0.915     | 0.950    | 0.969    |
| EffectiveCoverage      | 0.913     | 0.946    | 0.957    |
| LogisticFusion         | 0.743     | 0.763    | 0.774    |

- **Best default arm is BayesianMixture at every error rate** (0.915 / 0.950 /
  0.969). Its margin over the coverage-only baseline grows with error rate
  (+0.007 at 0.005 → +0.019 at 0.02): as more error k-mers appear, fusing the
  quality signal separates them from genomic k-mers more effectively than a
  fixed coverage cutoff.
- **All arms improve as error rate rises** (except QualityThreshold, pinned at
  0.5). Higher error rate spreads the genomic and error coverage distributions
  further apart (error k-mers cluster near coverage 1), so any coverage-aware
  rule gains specificity.
- **QualityThreshold's default Q20 cutoff is degenerate** — mean combined-Phred
  exceeds 20 for essentially every k-mer (quality is accumulated across
  observations), so specificity is 0.0 and balanced accuracy is exactly 0.5 in
  all nine runs. It is a "call everything solid" classifier at its default.
- **LogisticFusion's 0.5 posterior cutoff is off-optimum under class imbalance**
  (error k-mers outnumber genomic ~6–10:1). It trades recall for specificity at
  the wrong point (sensitivity ~0.58, specificity ~0.94), landing at ~0.76.

## 2. Threshold calibration (reflen 10000, seed 42, err 0.01)

Sweeping each arm's score over its observed range and picking the
balanced-accuracy-maximizing threshold:

| strategy               | default bal_acc | calibrated bal_acc | gain      | best threshold |
| ---------------------- | --------------- | ------------------ | --------- | -------------- |
| FixedCoverageThreshold | 0.938           | 0.950              | +0.012    | 2.5 (cov)      |
| QualityThreshold       | 0.500           | 0.950              | **+0.450** | 106.19 (Phred) |
| BayesianMixture        | 0.950           | 0.950              | +0.000    | 0.9965 (post.) |
| EffectiveCoverage      | 0.945           | 0.950              | +0.005    | 0.0 (post.)    |
| LogisticFusion         | 0.763           | 0.950              | **+0.187** | 0.081 (post.)  |

**Every arm converges to ~0.950 balanced accuracy once its threshold is
calibrated** — exactly the shared AUC ceiling. This is the headline: the arms do
not differ in discriminating power; they differ only in default-threshold
calibration.

- **QualityThreshold recovers the most (+0.45)**: from a degenerate 0.5 to 0.95.
  Its optimal cutoff is ~106 Phred (far above the default 20), because
  accumulated combined-Phred for genuine multi-observation genomic k-mers is very
  high; the cutoff that separates them from single-observation error k-mers sits
  well above any fixed per-base quality intuition.
- **LogisticFusion recovers +0.187**: its optimal posterior threshold is ~0.08,
  not 0.5 — the class imbalance pushes the balanced-accuracy-optimal cutoff far
  below the naive midpoint.
- **BayesianMixture is already near-calibrated** (+0.000): its default 0.5
  posterior cutoff essentially coincides with the balanced-accuracy optimum, the
  reason it tops the default-threshold grid.

## Takeaway

The scale grid and the calibration agree: at reflen 10000 the five arms share a
~0.95 AUC, so the classification problem at this scale is "solved" at the ranking
level. The only lever that matters for the operating point is threshold
selection. Two of the five arms (QualityThreshold, LogisticFusion) leave large
balanced-accuracy gains on the table with their defaults (+0.45, +0.19); a
calibration step that picks the balanced-accuracy-optimal cutoff per arm erases
the differences between them. For a default that needs no calibration,
BayesianMixture is the best-performing arm out of the box and its edge over
coverage-only widens with error rate.

### Reproduce

```bash
# Grid (one cell shown; loop seeds {42,43,44} × err {0.005,0.01,0.02}):
MYCELIA_S0_REFLEN=10000 MYCELIA_S0_SEED=42 MYCELIA_S0_ERRORRATE=0.01 \
  LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
  benchmarking/stage0_classification_matrix.jl

# Calibration:
MYCELIA_S0_REFLEN=10000 MYCELIA_S0_SEED=42 MYCELIA_S0_ERRORRATE=0.01 \
  LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
  benchmarking/stage0_threshold_calibration.jl
```
