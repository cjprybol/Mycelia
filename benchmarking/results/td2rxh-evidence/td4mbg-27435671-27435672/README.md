# Demoted-read reroute (PR #487): BASE vs FIX evidence

Pre-registered paired measurement of PR #487 ("stop forcing fully demoted reads into indel windows") against its merge base, on the td2rxh grid.

**Verdict: VOID.** A pre-registered control failed (nanopore rows differ between BASE and FIX in the 2000 bp cells, where every read is long and the PR's routing change should not apply).
The outcome statistic is therefore not a verdict; it is reported below descriptively only.

## What was run

| Run | Commit | Lawrencium job | Notes |
| --- | --- | --- | --- |
| Oracle (fixed toy, both commits) | FIX then BASE | 27435670 | one task, both commits on one node |
| BASE | `ca9eb736636101e6a4b3e6adf87687574608466b` (PR #487 merge base) | 27435672_[1-16] | |
| FIX | `dae40feaf6973c087ce956d08a25a7f7a46d5bfd` (PR #487 head) | 27435671_[1-16] | seed 11 failed before any arm ran; retried as 27437121_11 |

- Both checkouts used one byte-identical `Manifest.toml` (sha256 `85b75e3eb5c8…`, copied here). Both run labels carry the suffix `manifest85b75e3eb5c8`.
- Grid: `run_td2rxh_evidence_lrc.sbatch` with `--cells td2rxh`: 6000 bp genome × read length {1000, 2000} × 20x × error {0.01, 0.02, 0.05}, seeds 1-16, arms `:illumina` and `:nanopore` on identical reads.
- FIX seed 11 (27435671_11) failed on node n0018 while loading packages ("Unable to find compatible target in cached code image … runtime-disabled features"), with 0 rows written. It was rerun at the same commit, Manifest, and seed, excluding that node, as 27437121_11. Its logs are kept in `logs/`.
- All 32 seed tasks used for the analysis exited 0:0 with 12 data rows each (`sacct.txt`, `shard_line_counts.txt`).

## Decision rule (fixed before the runs)

- **G0a**: illumina rows identical BASE vs FIX in all 6 cells.
- **G0b**: nanopore rows identical BASE vs FIX in the three 2000 bp cells.
- **G0c**: the fixed-toy oracle passes every gate at both commits, with equal illumina sha.
- "Identical" means exact CSV field equality on every column except `label`, `wall_seconds`, `wall_seconds_per_read`.
- Any G0 failure voids the comparison.
- **G2** (outcome): at 1000 bp / err 0.01, per seed `d = nanopore − illumina` best-contig reference coverage in each run, `x = |d_BASE| − |d_FIX|`; exact one-sided sign-flip Wilcoxon signed-rank on the 16 `x`, alpha 0.05, H1: x > 0.
- The fail-open counters are descriptive only. The new reroute counter is not exported by the harness, so it was not captured.

The analysis code (`benchmarking/td4mbg_base_fix_analysis.jl`, with tests) was committed before any shard was analysed.

## Results

| Control | Result |
| --- | --- |
| G0a illumina identical, 6 cells | PASS |
| G0b nanopore identical, 2000 bp cells | **FAIL**: 61 fields over 20 of 48 (cell, seed) pairs; e0.01 5 seeds, e0.02 8, e0.05 7 |
| G0c oracle | PASS: 19/19 rows at both commits; illumina sha `d36e3b6a…d311` at both, equal to the historical value |

Most G0b differences are small (`n_contigs`, `total_contig_bases`).
The largest is err 0.02 seed 2: reassembly k 19 vs 17, coverage 0.670167 vs 0.644833.
Full list: `analysis/td4mbg_g0_mismatches.csv` (columns `a` = BASE, `b` = FIX).

Descriptive only (not a verdict, because G0 failed):

| Cell | n | one-sided p | median x | seeds with x > 0 | median d_BASE | median d_FIX |
| --- | --- | --- | --- | --- | --- | --- |
| 1000 bp / err 0.01 | 16 | 0.9998 | −0.1912 | 2 | −0.4261 | −0.5700 |
| 1000 bp / err 0.02 | 16 | 0.9533 | −0.0327 | 5 | −0.1085 | −0.2208 |
| 1000 bp / err 0.05 | 16 (1 zero) | 0.8738 | −0.0095 | 5 | 0.0203 | 0.0243 |

At err 0.01, FIX moved the nanopore-minus-illumina delta away from 0 in 14 of 16 seeds.
Nanopore wall time per read was lower under FIX in every cell (`analysis/td4mbg_descriptive_counters.csv`).

Context (enters no gate): an earlier run on master (`../27118659/`) vs BASE.
Illumina is identical in 16/16 seeds in every cell; nanopore is identical in 7-13 of 16 seeds per cell.

## Interpretation and caveats

- Reading the routing code, PR #487 changes the route only for reads with at most 1024 k-mer observations. The 2000 bp reads are full-length fragments with at least 1970 observations at the k values used. So the G0b failure is not explained by the routing change as written.
- The leading alternative is run-to-run nondeterminism in the nanopore arm. The context comparison is consistent with that, but it is not established: that comparison also spans a code and dependency change.
- A BASE-vs-BASE replicate of seeds 2 and 10 at `ca9eb736` (same Manifest) was submitted to separate the two explanations. The reading of each outcome was fixed before it ran:
  - If BASE does not reproduce itself, exact-identity G0b is unachievable and the protocol needs a determinism-controlled redesign.
  - If it does, the PR changes long-read results through a path the code reading missed.
- Until then, these data support neither merging nor rejecting PR #487 on measured grounds.
- Single genome size (6 kb), synthetic reads, one coverage; nothing here generalizes beyond the grid.

## Files

- `shards/{base,fix}/seed*.csv`: per-seed harness output (12 rows each).
- `oracle/{base,fix}/`: fixed-toy acceptance checks, arm summary, run manifest, rung telemetry.
- `analysis/`: `report.md`, G0/G2/descriptive/context CSVs, `td4mbg_g2_paired.{png,svg}`.
- `logs/`: SLURM stdout/stderr for every task, including the failed FIX seed 11 attempt and its retry.
- `provenance/`: the oracle sbatch and the Lawrencium setup script and log.
- `Manifest.toml`, `sacct.txt`, `shard_line_counts.txt`.

Reproduce the analysis from the repository root on a compute node:

```bash
julia +lts --startup-file=no --project=. benchmarking/td4mbg_base_fix_analysis.jl \
  --base-dir benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/shards/base \
  --fix-dir benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/shards/fix \
  --base-oracle benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/oracle/base/fixed_toy_acceptance_checks.csv \
  --fix-oracle benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/oracle/fix/fixed_toy_acceptance_checks.csv \
  --base-head ca9eb7366361 --fix-head dae40feaf697 --manifest 85b75e3eb5c8 \
  --context-dir benchmarking/results/td2rxh-evidence/27118659/shards \
  --output-dir /tmp/td4mbg-analysis
```
