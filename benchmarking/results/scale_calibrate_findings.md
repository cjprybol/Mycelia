# Stage 0 k-mer classification — scale grid + threshold calibration

Two follow-ups to the ground-truth comparison matrix
(`benchmarking/stage0_classification_matrix.jl`), regenerated on the corrected,
depth-independent quality signal (per-position mean Phred, `memory_profile=:full`)
and the full seven-arm set:

1. **Scale grid** — run the matrix across seeds {42, 43, 44} × error rates
   {0.005, 0.01, 0.02}, all at `MYCELIA_S0_REFLEN=10000` (cov 30, readlen 100,
   k 21); per-strategy mean ± sd of balanced accuracy and AUC.
2. **Threshold calibration** — for one base config (reflen 10000, seed 42,
   err 0.01), sweep each arm's continuous score over its observed range to find
   the balanced-accuracy-maximizing threshold, and measure the gain over each
   arm's default cutoff (`benchmarking/stage0_threshold_calibration.jl`).

## 1. Scale grid (3 seeds × 3 error rates, reflen 10000)

| strategy               | balanced accuracy | AUC           |
| ---------------------- | ----------------- | ------------- |
| FixedCoverageThreshold | 0.984 ± 0.005     | 1.000 ± 0.000 |
| QualityThreshold       | 0.500 ± 0.000     | 0.809 ± 0.014 |
| MixtureModel           | **0.998 ± 0.001** | 1.000 ± 0.000 |
| BloomFilter            | 0.884 ± 0.017     | 1.000 ± 0.000 |
| BayesianMixture        | 0.997 ± 0.001     | 1.000 ± 0.000 |
| EffectiveCoverage      | 0.992 ± 0.005     | 1.000 ± 0.000 |
| LogisticFusion         | 0.992 ± 0.002     | 1.000 ± 0.000 |

- **Coverage-based ranking is saturated (AUC ≈ 1.0)** for every arm except the
  quality-only threshold, which — now that quality is the depth-independent
  per-base mean rather than an accumulated coverage proxy — is honestly weaker
  (AUC 0.809): per-base quality alone separates genomic from error k-mers far
  less cleanly than coverage.
- **Best default arms are MixtureModel (0.998) and BayesianMixture (0.997).**
  The auto-threshold mixture-model (cutoff at the k-mer-spectrum valley) and the
  quality-fusion Bayesian mixture both essentially solve the problem at their
  default operating point.
- **BloomFilter (0.884) trades accuracy for memory**: at a fixed filter size its
  approximate membership admits false positives (specificity ~0.75), a cost that
  grows with the solid-set cardinality — the expected behavior of the
  memory-efficient tier.
- **QualityThreshold's default Q20 cutoff is degenerate** (specificity 0.0,
  balanced accuracy exactly 0.5): with real per-base quality, essentially every
  k-mer's mean Phred still clears 20, so it calls everything solid at its default.

## 2. Threshold calibration (reflen 10000, seed 42, err 0.01)

Sweeping each arm's score over its observed range and picking the
balanced-accuracy-maximizing threshold:

| strategy               | default bal_acc | calibrated bal_acc | gain      | best threshold |
| ---------------------- | --------------- | ------------------ | --------- | -------------- |
| FixedCoverageThreshold | 0.987           | 0.999              | +0.013    | 2.5 (cov)      |
| QualityThreshold       | 0.500           | 0.840              | **+0.340** | 34.67 (Phred)  |
| MixtureModel           | 0.999           | 0.999              | +0.000    | 0.194 (post.)  |
| BloomFilter            | 0.893           | 0.999              | +0.106    | 2.5 (cov)      |
| BayesianMixture        | 0.998           | 1.000              | +0.002    | 0.0 (post.)    |
| EffectiveCoverage      | 0.996           | 0.999              | +0.004    | 0.0 (post.)    |
| LogisticFusion         | 0.993           | 0.999              | +0.007    | 0.139 (post.)  |

- **QualityThreshold recovers the most (+0.34)** and its optimal cutoff is **34.67
  Phred — a physically realistic per-base value** in the simulated [2, 40] range.
  (Under the old confounded accumulated-Phred signal this optimum was ~106,
  impossible for a per-base mean; the realistic 34.67 is direct confirmation the
  quality signal is now depth-independent.) Even calibrated it tops out at 0.84,
  its AUC ceiling — per-base quality alone is a weaker discriminator.
- **BloomFilter recovers +0.106** once its cutoff is calibrated — the false
  positives are a threshold artifact, not a ranking failure (AUC 1.0).
- **The coverage-based and fusion arms are already near-calibrated** (gains
  ≤ +0.013): their default operating points essentially coincide with the
  balanced-accuracy optimum, which is why they top the default-threshold grid.

## Takeaway

At reflen 10000 the coverage-based arms share a ~1.0 AUC — the classification
problem is solved at the ranking level for any arm that uses coverage. Quality
alone is honestly weaker (AUC ~0.81). The one lever that matters for the operating
point is threshold selection: calibrating each arm to its balanced-accuracy-optimal
cutoff lifts the weaker defaults (QualityThreshold +0.34, BloomFilter +0.11) to the
near-perfect ceiling. For a default that needs no calibration, MixtureModel (auto
valley threshold) and BayesianMixture are the best-performing arms out of the box.

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
