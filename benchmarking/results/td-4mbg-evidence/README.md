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

Provenance, stated per artifact because the runs happened at different commits:

| artifact | master rows | fix rows |
| -------- | ----------- | -------- |
| `evidence_master.csv` / `evidence_fix.csv` | `4e068933` | `d817276b` |
| `undercorrection_identity.csv` | `4e068933` | `2e618d94` |

The contiguity run was launched before the review-fix commit landed, so its `fix`
rows reflect `d817276b`. The review fixes that followed (`2e618d94`) derive
`_INDEL_ADMITTED_WINDOW_SOURCES` and switch the window-level kernel test from
`!= :substitution` to ADMITTED membership — semantically identical while the two
tuples are exact complements, which a committed test now pins. The
byte-identity oracle WAS re-measured after that change and the run is committed
here as `oracle_after_review_fixes.txt` (18/18 gates, sha256 `d36e3b6a…d311`) --
not inferred from this sentence. So the contiguity numbers
are expected to carry over, but they were not re-measured at `2e618d94` and this
table should not be read as if they were.

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

## Under-correction — the trade-off, measured

Contiguity is a downstream proxy and cannot tell whether truncating the partition
leaves genuinely correctable bases uncorrected. `undercorrection_probe.jl`
measures it directly: per-read local identity (Smith-Waterman, `EDNAFULL`, both
orientations, denominator `matches + mismatches + insertions + deletions`) of the
**corrected** reads against the known reference, at 6000/1200/20x, 100 reads per
seed, byte-identical input per seed.

| seed | raw | master corrected | fix corrected | fix − master | master gain | fix gain |
| ---- | --- | ---------------- | ------------- | ------------ | ----------- | -------- |
| 42 | 0.98887 | 0.99095 | 0.99099 | +0.00004 | +0.00208 | +0.00212 |
| 43 | 0.98906 | 0.99306 | 0.99272 | −0.00034 | +0.00400 | +0.00366 |
| 44 | 0.98954 | 0.99479 | 0.99387 | −0.00092 | +0.00525 | +0.00433 |
| **mean** | | 0.99293 | 0.99253 | **−0.00041** | +0.00378 | **+0.00337** |

**The fix does under-correct.** Its corrected reads are lower-identity than
master's on 2 of 3 seeds (the third, +0.00004, is noise), and it captures **89%**
of master's per-read correction gain. It still improves on raw in 3 of 3, so it
is not destroying correction — it captures slightly less of what is available,
consistent with decoding fewer hard-span bases.

**So the trade is: ~11% less per-read correction gain, for 0.303 → 0.734 assembly
contiguity at 6 kb.** That direction is coherent rather than paradoxical —
over-correction can manufacture wrong k-mers that fragment the graph even while
nominally correcting more bases, so fewer better-targeted corrections can
assemble far better while scoring slightly worse per read.

Whether the trade is worth taking depends on what the corrector feeds. For an
assembly-directed pipeline (which is what td-jt7r is) contiguity dominates and
the trade is favourable. For a consumer of corrected reads directly, the
0.041-percentage-point identity loss — about half a base per 1200 bp read —
matters more. n=3, one cell; direction is consistent, magnitude is not
characterized.
