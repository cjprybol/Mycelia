# td-4mbg evidence — read-level indel predicate, before/after

Committed snapshot of the measurement that argues for gating complete-span
windowing on `_read_runs_indel_decode`. It exists because the first version of
this evidence came from an uncommitted scratch script: nothing in the repo
reproduced the numbers, and a reviewer could not tell which read generator had
been used — which is decisive, because the two committed candidates differ in
whether the reads contain indels at all.

## How to regenerate

```bash
# against an unfixed checkout
MYCELIA_PROJECT=/path/to/master LD_LIBRARY_PATH='' \
  julia --project=/path/to/master benchmarking/indel_reroute_evidence.jl \
    --seeds 42,43,44 --label master --out evidence_master.csv

# against this branch
MYCELIA_PROJECT=/path/to/branch LD_LIBRARY_PATH='' \
  julia --project=/path/to/branch benchmarking/indel_reroute_evidence.jl \
    --seeds 42,43,44 --label fix --out evidence_fix.csv
```

`master` rows were produced at `4e068933`; `fix` rows at `2e618d94`.

## What was measured

Three cells × three seeds × two arms, per code version. Within a cell+seed both
arms receive **byte-identical** reads and differ only in `sequencing_tech`:

- `illumina` — substitution-only profile. Builds no window-source map, so the
  change cannot reach it. This is the **control**.
- `nanopore` — indel profile. Under `:frontier_budgeted` its scheduler emits the
  window-source map that the change reroutes.

Reads come from `Mycelia.observe(...; tech = :nanopore)` — 40% mismatch /
30% insertion / 30% deletion — so they carry real indels and the indel arm has
something to recover.

## Control check

**0 of 9 `illumina` rows differ between `master` and `fix`.** The change
provably does not touch the substitution-only arm. Without this the nanopore
comparison would not be interpretable.

## Result — largest-contig ratio, nanopore arm

| cell | seed | illumina | master | fix | Δ | divergences | anchor rej |
| ---- | ---- | -------- | ------ | --- | - | ----------- | ---------- |
| 2000/1200/8x | 42 | 0.830 | 0.452 | 0.700 | +0.247 | 6 → 6 | 5 → 5 |
| 2000/1200/8x | 43 | 0.825 | 0.952 | 0.952 | +0.000 | 8 → 8 | 5 → 5 |
| 2000/1200/8x | 44 | 0.718 | 0.830 | 0.830 | +0.000 | 25 → 22 | 9 → 4 |
| | mean | | 0.744 | 0.827 | **+0.083** | sd 0.261 → 0.126 | |
| 2000/1200/20x | 42 | 0.904 | 0.680 | 0.964 | +0.284 | 53 → 34 | 23 → 13 |
| 2000/1200/20x | 43 | 0.906 | 0.476 | 0.694 | +0.218 | 21 → 6 | 46 → 2 |
| 2000/1200/20x | 44 | 0.971 | 0.982 | 0.725 | **−0.257** | 40 → 26 | 38 → 8 |
| | mean | | 0.713 | 0.794 | **+0.082** | sd 0.254 → 0.148 | |
| 6000/1200/20x | 42 | 0.909 | 0.442 | 0.465 | +0.022 | 164 → 108 | 126 → 22 |
| 6000/1200/20x | 43 | 0.930 | 0.208 | 0.839 | +0.631 | 177 → 98 | 128 → 13 |
| 6000/1200/20x | 44 | 0.907 | 0.258 | 0.898 | +0.641 | 108 → 44 | 108 → 21 |
| | mean | | 0.303 | 0.734 | **+0.431** | sd 0.124 → 0.235 | |

## What this supports, and what it does not

**Supports.** At 6000/1200/20x — the cell the epic cares about — the mean
largest-contig ratio moves 0.303 → 0.734 against an illumina control of ~0.915,
and all three seeds improve. Divergences and anchor rejections fall alongside it.

**Does not support.**

- **A regression exists.** 2000/1200/20x seed 44 goes 0.982 → 0.725. It is the
  only negative delta in nine and is unexplained. It is reported, not averaged
  away.
- **Magnitude is not characterized.** n=3. The spread at 6000/1200/20x actually
  *widens* (sd 0.124 → 0.235), so the fixed arm is more variable, not less.
- **Single-seed readings are unreliable here, in both directions.** Seed 42 at
  6000/1200/20x gives +0.022 and seed 43 gives +0.631 — from the same cell.
  Largest-contig is a discontinuous order statistic; one correction can join or
  split a contig. Do not quote one seed as an effect size.
- **Two of three 8x seeds are exact no-ops.** The change acts only where the
  scheduler emits demoted windows over spans exceeding `max_window`, which tracks
  coverage and genome size.
- **Contiguity does not answer under-correction.** The truncated partition
  decodes fewer bases, so it could leave genuinely correctable bases uncorrected
  without contiguity showing it. That needs per-base identity against the
  reference and is tracked separately.
