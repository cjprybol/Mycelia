# opt1: serial-vs-parallel corrector scaling (phix, 20x, batch_size=50)

Records the measured wall-clock effect of the opt1 `Threads.@threads` per-read
decode change on the `:scalable` corrector, and confirms byte-identical accuracy
across thread counts. This is a smoke-scale sanity measurement, not a scaling
study — see caveats below before drawing any conclusion beyond "identical
output, directional timing."

## Setup

- Benchmark: `benchmarking/rhizomorph_correction_accuracy_benchmark.jl`
- Genome: phiX174 (`NC_001422`, 5386 bp)
- Coverage: 20x target (719 simulated reads, 150bp, effective 20.02x)
- Error rate: 0.01 substitution-only (single cell — not the default 4-rate
  sweep)
- `MYCELIA_RCA_BATCH_SIZE=50` (forced multi-batch: ~15 batches of ~48 reads
  each, vs the 10000 default which would put all reads in one batch)
- `MYCELIA_RCA_SEED=42`, `MYCELIA_RCA_SCALE_FLOOR=1` (bypass the default scale
  guard at this genome size)
- Machine: Apple M5 Max, 18 logical cores (macOS 25.5.0 / Darwin arm64)
- Julia 1.10.10, `--project=.`, `LD_LIBRARY_PATH=''`
- 3 repetitions per thread count (`-t 1` and `-t 8`), run back-to-back on an
  otherwise-idle laptop, no other heavy processes observed running concurrently

## Accuracy identity (load-bearing result)

All 6 runs (3x `-t 1`, 3x `-t 8`) produced **byte-identical** accuracy metrics:

- `recall = 0.9276773296244785`
- `precision = 1.0`
- `over_correction_rate = 0.0`
- `correction_rate = 0.012369031061659713`

No divergence across thread counts on any of the 6 runs. This confirms the opt1
parallel decode path is deterministic and produces the same corrected output
regardless of thread count, at least at this scale and seed.

## Timing (`corrector_runtime_s`, directional — not a scaling study)

| Rep      | `-t 1` (s)  | `-t 8` (s)  |
| -------- | ----------- | ----------- |
| 1        | 17.76       | 17.961      |
| 2        | 15.283      | 16.84       |
| 3        | 17.092      | 15.8        |
| **mean** | **16.71**   | **16.87**   |
| min–max  | 15.28–17.76 | 15.80–17.96 |

Mean speedup ratio (`-t 1` mean / `-t 8` mean) = 16.71 / 16.87 ≈ **0.99x** —
i.e. no measurable speedup at this scale, and the two distributions overlap
entirely (15.28–17.96s spans both thread counts). `-t 8` is not reliably faster
than `-t 1` here; the run-to-run spread (~2.5s, ~15% of the mean) is comparable
to or larger than any thread-count effect.

`null_runtime_s` (Control B, read-scramble arm, same corrector code path) showed
more spread and no consistent pattern either (`-t 1`: 8.907, 6.684, 5.381s;
`-t 8`: 5.273, 3.869, 4.057s) — directionally lower under `-t 8` across all 3
pairs, but 3 reps on a shared machine is not enough to call that a confirmed
effect, and it is not the primary metric this task tracks.

## Honest interpretation

- **This is not evidence that opt1 is broken or ineffective.** The task brief
  anticipated this outcome: "phix is small → limited parallel headroom." At 719
  reads / ~15 batches of ~48 reads each, per-batch Viterbi decode work is small
  relative to thread-spawn/sync overhead and to the serial portions of the
  pipeline (graph construction, I/O, reference download/parsing), so
  `Threads.@threads` over per-read decode doesn't have enough parallel work per
  batch to show a clear win at this scale.
- **The value of opt1 demonstrated here is the accuracy-identity result, not a
  timing win.** Parallel and serial decode produce byte-identical
  recall/precision/correction_rate across 6 independent runs — the correctness
  bar opt1 needed to clear.
- **Do not extrapolate a multiplier from this run.** A single-genome,
  single-scale A/B on a shared laptop with ~15% run-to-run timing noise is not
  sufficient to report a headline speedup number. Larger genomes / higher
  coverage (more reads per batch, more total decode work) are expected to show
  more parallel headroom, since the fixed per-run overhead (reference
  acquisition, graph build) amortizes over more per-read decode work — but that
  expectation is untested here and should be measured separately, not asserted
  from this result.
- **Noise caveat:** timing was collected on a general-purpose laptop with no
  isolation from background processes; 3 reps per thread count is a minimal
  sample for characterizing variance, not a rigorous benchmark.
