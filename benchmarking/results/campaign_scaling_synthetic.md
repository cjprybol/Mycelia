# Cross-campaign corrector wall-clock benchmark (opt5 + opt1 + opt2)

Measures the CUMULATIVE wall-clock effect of the three merged corrector
performance PRs — #433 (opt5, between-batch GC opt-in), #434 (opt1, thread-safe
parallel decode), #435 (opt2, typed transitions + fused edge-scan) — by
comparing the pre-campaign baseline against the fully-merged campaign checkout
on a decode-dominated fixture, and confirms byte-identical corrected output.

**Fixture note:** this run uses a **synthetic random genome, NOT E. coli** —
the "E. coli-scale" wall-clock target motivating the campaign remains
unmeasured. See "Caveats" below: this run does not reach, and must not be
extrapolated to, E. coli @30x.

## Setup

- **BASELINE checkout:** `5b7b495e` (parent of PR #433), scratch worktree at
  `/tmp/precampaign-bench`.
- **CAMPAIGN checkout:** `218bfa37` (opt5+opt1+opt2 merged), worktree at
  `/Users/cameronprybol/workspace/Mycelia/.worktrees/campaign-bench`.
- **Corrector entry point:** `Mycelia.mycelia_iterative_assemble`, called
  identically on both checkouts with:
  `max_k=21, n_k_rungs=3, max_iterations_per_k=2, graph_mode=:canonical,
  skip_solid=true, cheap_correct=true, hard_window=true, soft_em=true,
  verbose=false, enable_checkpointing=false`. `gc_between_batches` and
  `MYCELIA_RCA_BATCH_SIZE` (opt5-only knobs) were **never** passed;
  `batch_size` was left at its shared default (10000), so all runs place every
  read of this fixture in a single decode batch.
  - **Consequence — opt5 is _not exercised_ in this configuration.** opt5
    gates the between-batch `GC.gc()`, which fires only when
    `batch_end < total_reads` (i.e. across ≥2 batches). With every read in a
    single batch, that GC fires in **neither** checkout (verified: baseline
    `iterative-assembly.jl:4941`, campaign `:5044`), so opt5 contributes
    nothing measurable here. The serial gain below is therefore **opt2
    alone**, and the total is **opt1+opt2**; opt5's between-batch-GC benefit
    requires a multi-batch run (`batch_size` < read count) and remains
    untested at this fixture scale.
- **`enable_parallel` was passed explicitly** (`true` for the campaign-parallel
  arm, `false` elsewhere) — see "Deviation from the original plan" below for
  why.
- **Fixture:** self-contained, Mycelia-independent generator
  (`benchmarking/campaign_scaling_fixture_gen.jl`, Base + Random + FASTX only,
  sequential deterministic read IDs, seed=42) — a 2000 bp synthetic random DNA
  genome, 150 bp single-end reads, 400x target coverage → **5334 reads**,
  1% per-base substitution error rate. Generated once and reused byte-identically
  across every run (`benchmarking/results/` does not carry the fixture itself;
  it is reproducible from the seed and the committed generator script).
- **Timing script:** `benchmarking/campaign_scaling_corrector_bench.jl`. Does a
  warm-up call on an unrelated tiny (300 bp) fixture first (compilation only,
  reported separately, NOT included in the measured runtime), then times
  **only** the real `mycelia_iterative_assemble` call:
  `t0=time(); result=mycelia_iterative_assemble(...); runtime_s=time()-t0`.
  `GC.gc()` is forced immediately before `t0`. The reported `runtime_s` therefore
  excludes Julia package precompilation/JIT warm-up but is otherwise an
  unmodified wall-clock corrector call (I/O, graph construction, decode,
  reassembly all included).
- **Machine:** Apple M5 Max, 18 logical cores, macOS 26.5.2, Julia 1.10.10,
  `LD_LIBRARY_PATH=''`. **This machine was shared with 13 concurrent logged-in
  users during the run** (see Caveats) — not an isolated benchmark box.
- **Reps:** 2 per configuration, run strictly sequentially (never overlapping
  in wall-clock time) to avoid cross-run CPU contention from the benchmark
  itself.

## Deviation from the original plan (documented, not silent)

The task brief assumed `mycelia_iterative_assemble` itself defaults
`enable_parallel` based on thread count post-opt1 ("campaign defaults
`enable_parallel=nthreads>1`"). That is **not where the default lives**: reading
`src/iterative-assembly.jl` in both checkouts, `mycelia_iterative_assemble`'s
own keyword default is `enable_parallel::Bool = false` in **both** the baseline
and the campaign checkout, byte-identically. The nthreads-aware auto-default
(`enable_parallel = nthreads>1` for the `:scalable` strategy) lives one layer up,
in `_corrector_strategy_knobs(:scalable)` inside `src/rhizomorph/assembly.jl`
(added by opt1, PR #434) — a function used by `assemble_genome`, not by
`mycelia_iterative_assemble` directly.

Consequence: calling `mycelia_iterative_assemble` with no `enable_parallel`
argument runs **serially regardless of `--threads`** on both checkouts, which
would make a pure default-passthrough comparison measure nothing about opt1.
To still exercise opt1's thread-safe decode path through the stable,
version-independent `mycelia_iterative_assemble` signature (present unchanged,
same keyword, same semantics, in both checkouts — this is not one of the
opt5-only knobs the brief excluded), `enable_parallel` is passed **explicitly**:
`true` only for the "campaign parallel" arm, `false` (baseline) elsewhere. This
is a keyword that predates the whole campaign; only the code path that
auto-sets its default changed.

One more relevant fact found in the same read: pre-opt1, baseline's decode
function additionally special-cased soft-EM —
`use_parallel = enable_parallel && soft_weights === nothing`, i.e. **baseline
silently reverts to serial (with a `@warn`) whenever `soft_em=true`**, which
this fixture's call always sets. Post-opt1, campaign's decode function is
`use_parallel = enable_parallel` unconditionally — soft-EM no longer disables
parallelism. This is exactly opt1's fix and is why "campaign parallel" is a
meaningful arm at all; passing `enable_parallel=true` to the *baseline* would
not have produced a real parallel run (soft_em=true forces it back to serial
there), so that combination was correctly not run.

## Results (raw)

| Config                         | `--threads` | `enable_parallel` | rep | runtime_s   | final FASTQ sha256 (first 16 hex) |
| ------------------------------- | ----------- | ------------------ | --- | ----------- | ---------------------------------- |
| baseline (`5b7b495e`) serial     | 1           | false               | 1   | 106.559     | `e67cc706da51a59f`                 |
| baseline (`5b7b495e`) serial     | 1           | false               | 2   | 105.874     | `e67cc706da51a59f`                 |
| campaign (`218bfa37`) parallel   | 8           | true                | 1   | 42.887      | `e67cc706da51a59f`                 |
| campaign (`218bfa37`) parallel   | 8           | true                | 2   | 64.427      | `e67cc706da51a59f`                 |
| campaign (`218bfa37`) serial     | 1           | false               | 1   | 90.317      | `e67cc706da51a59f`                 |
| campaign (`218bfa37`) serial     | 1           | false               | 2   | 101.378     | `e67cc706da51a59f`                 |

All 6 runs produced **the same 64-hex-character SHA256**
(`e67cc706da51a59f3f317352cd2e966ee5acdb6d2334eec4d69f531f0adff062`) of the
final corrected FASTQ (5334 reads each). SHA256 equality is full-file
byte-identity, not a sampled check.

> The campaign-parallel rep swing (42.9 → 64.4s) and campaign-serial swing
> (90.3 → 101.4s) above are contention artifacts, not corrector variance — see
> the "Shared machine, heavy concurrent load" caveat below.

Per-config summary (mean, min–max, and spread across 2 reps):

| Config                        | mean runtime_s | min–max         | spread as % of mean |
| ------------------------------ | --------------- | ---------------- | -------------------- |
| baseline (`5b7b495e`) `-t1`     | 106.22           | 105.87–106.56    | 0.6%                 |
| campaign (`218bfa37`) `-t1`     | 95.85            | 90.32–101.38     | 11.5%                |
| campaign (`218bfa37`) `-t8`     | 53.66            | 42.89–64.43      | 40.1%                |

## Speedup ratios

Using **per-rep-mean** runtimes:

| Comparison                                            | Ratio                        | Interpretation                          |
| ------------------------------------------------------ | ----------------------------- | ---------------------------------------- |
| Serial gain: baseline `-t1` / campaign `-t1`            | 106.22 / 95.85 = **1.11x**    | opt2 serial-path contribution (opt5 not exercised — single batch) |
| Parallel gain: campaign `-t1` / campaign `-t8`          | 95.85 / 53.66 = **1.79x**     | opt1 thread-scaling contribution         |
| **Total campaign speedup:** baseline `-t1` / campaign `-t8` | 106.22 / 53.66 = **1.98x**  | opt1+opt2 (opt5 not exercised at this scale) |

Using the **least-contended (min) rep per config** — an alternative estimate
under the assumption that contention noise, not real variance, dominates the
spread (see Caveats):

| Comparison                                  | Ratio                         |
| --------------------------------------------- | ------------------------------ |
| Serial gain (min reps)                        | 105.87 / 90.32 = **1.17x**     |
| Parallel gain (min reps)                      | 90.32 / 42.89 = **2.11x**      |
| **Total campaign speedup (min reps)**         | 105.87 / 42.89 = **2.47x**     |

Both the mean-based (1.98x) and min-based (2.47x) totals are reported because
the two arms have very different noise sensitivity (see below) — **neither
number should be quoted alone as "the" campaign speedup**; the honest range on
this fixture, this machine, this run is **roughly 2.0x–2.5x total, with the
serial contribution (opt2 alone; opt5 not exercised here) around 1.1x–1.2x and
the parallel contribution (opt1, 8 vs 1 thread) around 1.8x–2.1x**. Note the
min-selection is asymmetric: baseline's spread is 0.6% (min≈mean) while the
campaign arms swing 11.5%/40.1%, so "min" corrects the campaign arms far more
than the baseline, inflating the min-based ratio partly by selection rather
than purely by noise removal — and with only 2 reps, "min" is itself a weak
noise-floor estimator.

## Byte-identity (6/6 runs identical, this fixture scale)

All 6 runs (2 baseline-serial, 2 campaign-parallel, 2 campaign-serial) produced
the exact same corrected-FASTQ SHA256. The campaign's byte-identity claim holds
at this fixture's scale (5334 reads, 2 kb genome, k=21/3 rungs/2 iterations,
:canonical mode, soft-EM on).

## Caveats (read before quoting any number above)

- **This is NOT E. coli, and NOT @30x coverage.** The fixture is a 2 kb
  synthetic random-DNA genome at 400x coverage / 5334 150 bp reads — several
  orders of magnitude smaller than an E. coli genome (~4.6 Mb) at 30x. Do not
  extrapolate these ratios to that scale; superlinear terms in the `:scalable`
  corrector's decode path (documented separately in
  `benchmarking/shortread_scalable_scaling_profile.jl`, which found short-read
  `:doublestrand`/windowed-decode runtime scaling from 18s@1kb to 1346s@5kb,
  alpha≈2.7) mean genome-length scaling is emphatically NOT linear for this
  corrector family, and thread-count scaling at a larger genome is untested
  here.
- **Shared machine, heavy concurrent load — this is the dominant noise
  source.** `uptime` during the run showed **13 concurrent logged-in users**
  and a load average that spiked to **47.6** (baseline: ~16–17 outside the
  spike) on an 18-logical-core machine — i.e. the box was oversubscribed by
  ~2.6x at the worst-observed moment, which landed during campaign-parallel
  rep 2 (64.43s, 50% slower than rep 1's 42.89s) and bled into campaign-serial
  rep 2 (101.38s vs rep 1's 90.32s). Baseline's two reps were tighter (0.6%
  spread) only because they happened to run before the load spike, not because
  baseline is inherently less noisy.
  - **The parallel arm is structurally more exposed to this kind of noise**
    than the serial arm: an 8-thread job competing for 8 of 18 logical cores
    against other users' processes will lose wall-clock time to preemption in
    a way a 1-thread job mostly won't. This is the most likely explanation for
    campaign-`t8`'s 40% rep spread vs campaign-`t1`'s 11.5% and baseline's
    0.6% — not a property of opt1 itself.
- **runtime_s excludes Julia package compilation** (a separate warm-up call on
  an unrelated tiny fixture is timed and reported per-run but subtracted out by
  construction — only the real corrector call on the real fixture is inside
  `t0`/`runtime_s`). It does NOT exclude method-specialization compilation
  triggered only by code paths the warm-up doesn't hit (e.g. `enable_parallel`
  branches only exercised when `enable_parallel=true`); the warm-up call is run
  with the same `enable_parallel` value as the timed call specifically to
  minimize this, but a small residual first-call JIT cost inside `runtime_s`
  cannot be fully ruled out.
- **2 reps per configuration** is a minimal sample for characterizing variance,
  chosen to fit comfortably inside the time budget (each run: 43–107s, total
  battery ~8.5 minutes of corrector wall time). Given the observed spread
  (up to 40% on the `-t8` arm), a 2-rep mean should be read as a rough
  central estimate, not a precise measurement — this is why both the mean- and
  min-based ratios are reported above rather than a single headline number.
- **No accuracy/quality claim is made here.** Both checkouts logged a large
  volume of `path_to_sequence: no (k-1) overlap ... falling back to
  stored/forward orientation` warnings on this synthetic-random-genome fixture
  (hundreds of thousands of lines per run) — expected on a repeat-free random
  genome under `:canonical` mode and not investigated further, since this
  benchmark only measures wall-clock and byte-identity, not corrector quality.
  These warnings are identical in kind and (as far as observed) volume across
  both checkouts, so they should not differentially bias the timing comparison,
  but they were not independently verified to be volume-matched line-for-line.
- **Single machine, no repeated-measures statistics.** No confidence interval,
  no outlier-rejection procedure beyond reporting min/max/mean; this is a
  directional cross-campaign measurement, not a publication-grade benchmark.

## Honest bottom line

On a small (2 kb / 5334-read), decode-dominated, `:canonical`-mode fixture with
soft-EM enabled, the merged campaign checkout is **directionally and
substantially faster** than the pre-campaign baseline (roughly **2.0x–2.5x**
total, split into a ~1.1–1.2x serial contribution from opt2 and a ~1.8–2.1x
parallel contribution from opt1's 1→8-thread decode; opt5's between-batch GC is
not exercised in this single-batch config), and **produces byte-identical
corrected output** across all 6 runs at this scale. The measurement is noisy — dominated by a
heavily shared 13-user machine during part of the run, most visibly on the
8-thread arm — and does not reach, and must not be read as evidence about,
E. coli-scale genomes or realistic 30x whole-genome coverage.

## Reproduce

```bash
# Fixture (deterministic; ~1s). FX_COVERAGE=400 is REQUIRED to get 5334 reads
# — the generator default (50) yields only 667.
FX_OUT=/tmp/campaign_fixture.fastq FX_GENOME_LEN=2000 FX_READLEN=150 \
  FX_COVERAGE=400 FX_ERR=0.01 FX_SEED=42 \
  LD_LIBRARY_PATH='' julia --project=. benchmarking/campaign_scaling_fixture_gen.jl

# Timed run (example: campaign checkout, 8 threads, parallel):
BENCH_FASTQ=/tmp/campaign_fixture.fastq BENCH_OUTDIR=/tmp/campaign_out \
  BENCH_ENABLE_PARALLEL=true \
  LD_LIBRARY_PATH='' julia --project=. --threads=8 \
  benchmarking/campaign_scaling_corrector_bench.jl
```
