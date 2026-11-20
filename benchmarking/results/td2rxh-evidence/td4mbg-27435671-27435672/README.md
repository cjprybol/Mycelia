# Demoted-read reroute (PR #487): BASE vs FIX evidence

Pre-registered paired measurement of PR #487 ("stop forcing fully demoted reads into indel windows") against its merge base, on the td2rxh grid.

**Verdict: VOID.** A pre-registered control failed: nanopore rows differ between BASE and FIX in the 2000 bp cells, where every read is long and the pre-registration assumed the PR's routing change would not apply.
The outcome statistic is therefore not a verdict; it is reported below descriptively only.

## What was run

| Run | Commit | Lawrencium job | Notes |
| --- | --- | --- | --- |
| Oracle (fixed toy, both commits) | FIX then BASE | 27435670 | one task, both commits on one node |
| BASE | `ca9eb736636101e6a4b3e6adf87687574608466b` (PR #487 merge base) | 27435672_[1-16] | |
| FIX | `dae40feaf6973c087ce956d08a25a7f7a46d5bfd` (PR #487 head) | 27435671_[1-16] | seed 11 failed before any arm ran; retried as 27437121_11 |

- Both checkouts used one byte-identical `Manifest.toml` (sha256 `85b75e3eb5c8…`, copied here). Both run labels carry the suffix `manifest85b75e3eb5c8`.
- Grid: `run_td2rxh_evidence_lrc.sbatch` with `--cells td2rxh`: 6000 bp genome × read length {1000, 2000} × 20x × error {0.01, 0.02, 0.05}, seeds 1-16, arms `:illumina` and `:nanopore` on identical reads.
- FIX seed 11 (27435671_11) failed on node n0018 while loading packages ("Unable to find compatible target in cached code image … runtime-disabled features"), with 0 rows written. It was rerun at the same commit, Manifest, and seed as 27437121_11, submitted from the FIX checkout's `benchmarking/` with `sbatch --array=11 --exclude=n0018.lr6 run_td2rxh_evidence_lrc.sbatch`; it ran on n0341. Both attempts' logs are kept in `logs/`.
- All 32 seed tasks used for the analysis exited 0:0 with 12 data rows each (`sacct.txt`, `shard_line_counts.txt`).

## Decision rule

The pre-registration and its amendments are kept in the project's tracking notes, not in this repository; what follows is a summary.
Amendments 1-4 (the gates below, the BASE-vs-FIX design, the shared Manifest) were recorded on 2026-10-28, before any run.
Amendment 5 was recorded on 2026-10-30 after the runs were launched and before any shard was analysed. It pinned the "identical" column rule and the one-sided test details, and recorded that the reroute counter is not exported by the harness.

- **G0a**: illumina rows identical BASE vs FIX in all 6 cells.
- **G0b**: nanopore rows identical BASE vs FIX in the three 2000 bp cells.
- **G0c**: the fixed-toy oracle passes at both commits, with the illumina sha equal between them.
  - Amendment 4 says the oracle "passes 19/19" and that equality with the historical prewiring sha is reported descriptively: a mismatch there points to dependency drift and does not void the comparison. The one `advisory` row is the host-load-dependent wall-clock check.
  - Amendment 5 says every `severity=gate` row passes, and that the illumina sha in `illumina_byte_identical_to_prewiring_oracle` is equal between the commits. That check is itself a `gate` row, so amendment 5's "every gate row" wording would gate the historical comparison. That conflicts with amendment 4.
  - The implementation follows amendment 4: the advisory row and the historical-sha gate row are reported, not gated, and cross-commit equality of the illumina sha is still required. In this run all 19 rows passed at both commits, so neither reading changes the result.
- "Identical" means exact CSV field equality on every column except `label`, `wall_seconds`, `wall_seconds_per_read`.
- Any G0 failure voids the comparison.
- **G2** (outcome): at 1000 bp / err 0.01, per seed `d = nanopore − illumina` best-contig reference coverage in each run, `x = |d_BASE| − |d_FIX|`; exact one-sided sign-flip Wilcoxon signed-rank on the 16 `x`, alpha 0.05, H1: x > 0.
- The fail-open counters are descriptive only. The new reroute counter is not exported by the harness, so it was not captured.

The analysis code (`benchmarking/td4mbg_base_fix_analysis.jl`, with tests) was first committed at `94c0cbb9` (2026-10-30 10:28 PDT), before any array task had ended: FIX seed 11's failed attempt ended first, at 10:36 PDT; the first task to complete did so at 10:53 PDT and the last at 11:58 PDT (`sacct.txt`).
The oracle job had already finished (09:52 PDT), and its checks (G0c inputs) had been read when the code was committed.
After code review the script was revised, again with the decision rule unchanged. The revision adds input validation, oracle provenance and shape checks, and the read-integrity checks of the earlier td2rxh analysis, and fixes report formatting and the figure.
Re-running the script on these inputs reproduces the original verdict and every G0 and G2 number. The files in `analysis/` are the outputs of the script as run for this evidence. A later edit of the script's wall-time caveat no longer asserts node placement, which this analysis does not record; that edit does not rewrite `analysis/report.md`, and the verdict line there is unchanged.

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

Median nanopore wall time per read was lower under FIX in every cell, including the 2000 bp cells where the code path should be unchanged (`analysis/td4mbg_descriptive_counters.csv`).
BASE and FIX tasks ran on different nodes, so this difference is confounded with node placement and is not attributable to the PR.

Context (enters no gate): an earlier run on master (`../27118659/`) vs BASE.
Illumina is identical in 16/16 seeds in every cell; nanopore is identical in 7-13 of 16 seeds per cell.

## Interpretation and caveats

- **The routing change, as written, does not explain G0b.**
  - In PR #487, a read is "demoted" only in a pass that has a window-source map, and such a map exists only when bounded windowed decoding is active, which is the condition under which long reads are windowed anyway.
  - So the new route differs from the old one only for reads with at most 1024 k-mer observations.
  - The 2000 bp reads are full-length simulated fragments. Nanopore errors change their length by a few tens of bases, so at the k values used (at most 31) they have about 1970 observations, far above 1024.
- **Three explanations remain, and they are not yet separated:**
  - run-to-run nondeterminism in the nanopore arm;
  - hardware-dependent numerics: the lr6 partition mixes CPU feature sets (the seed 11 failure on n0018 shows this), and no BASE/FIX pair of a seed ran on the same node;
  - a path in the PR that this code reading missed.
  The context comparison is consistent with either of the first two, but it also spans a code and dependency change.
- A BASE-vs-BASE replicate of seeds 2 and 10 at `ca9eb736` (same Manifest, `--exclude=n0018.lr6`, job 27438474) was run to start separating them. It ran on n0230 (seed 2) and n0211 (seed 10); the original BASE tasks for those seeds ran on n0216 and n0375. The reading was fixed before it ran, as amendment 6 records it (diagnostic; it enters no gate). Compare the replicate with the BASE shards under the G0 field rule:
  - If the replicate's nanopore rows differ from BASE in the 2000 bp cells: the nanopore arm is nondeterministic, exact-identity G0b is unachievable as designed, the VOID stands, and the measurement needs a redesigned protocol (for example a determinism-controlled run or a replicate-based noise floor) before any merge decision rests on it.
  - If those nanopore rows are identical: PR #487 changes long-read results through a path the code reading missed, amendment 1's premise is falsified, and the PR needs investigation before merge.
  - Illumina rows must also replicate; if they do not, the whole grid is suspect.

## Replicate result: BASE does not reproduce itself

Both replicate tasks completed 0:0 with 12 rows each (`replicate-27438474/`).
Compared with the original BASE rows under the G0 field rule (`replicate-27438474/comparison.txt`, produced by `compare_replicate.py`):

| Arm | Rows compared | Rows differing |
| --- | --- | --- |
| illumina | 12 | 0 |
| nanopore | 12 | 3, all seed 10 |

The three differing nanopore rows:

- 1000 bp / err 0.01: `total_contig_bases` 26690 → 26691.
- 1000 bp / err 0.05: `substitution_length_divergences` 587 → 585, `n_contigs` 2837 → 2836, `total_contig_bases` 66562 → 66552.
- 2000 bp / err 0.01: `indel_engaged` 28 → 29, `n_contigs` 470 → 471, `total_contig_bases` 34871 → 34904.

This is the first pre-registered outcome: the nanopore arm does not reproduce at the same commit and Manifest on this setup, so exact-identity G0b is unachievable as designed and the VOID stands.
Because the replicate ran on different nodes, it does not separate run-to-run nondeterminism from node-dependent numerics.

Descriptive and not pre-registered, from two seeds only:

- Seed 2, including the 2000 bp / err 0.02 cell where BASE and FIX differed most (reassembly k 19 vs 17), reproduced exactly.
- In all 12 replicated nanopore rows, best-contig reference coverage was unchanged; the differences were at most one contig and 33 bases, plus changes of 1-2 in two corrector counters (`substitution_length_divergences`, `indel_engaged`).
- At 2000 bp / err 0.01 seed 10, `total_contig_bases` differs by 2 between BASE and FIX (34871 vs 34873) but by 33 between BASE and its replicate (34871 vs 34904): there, run-to-run noise exceeds the BASE-vs-FIX difference.
- That is small next to the BASE-vs-FIX coverage changes at 1000 bp / err 0.01 (median −0.19), but two seeds cannot bound the run-to-run variance.
- Until then, these data support neither merging nor rejecting PR #487 on measured grounds.
- Single genome size (6 kb), synthetic reads, one coverage; nothing here generalizes beyond the grid.

## Files

- `shards/{base,fix}/seed*.csv`: per-seed harness output (12 rows each).
- `oracle/{base,fix}/`: fixed-toy acceptance checks, arm summary, run manifest, rung telemetry.
- `analysis/`: `report.md` and the G0/G2/descriptive/context CSVs. The figure `td4mbg_g2_paired.{png,svg}` is not committed (the repository ignores plots under `td2rxh-evidence/*/analysis/`); the reproduce command below regenerates it.
- `logs/`: SLURM stdout/stderr for every task, including the failed FIX seed 11 attempt and its retry.
- `provenance/`: the oracle sbatch and the Lawrencium setup script and log.
- `replicate-27438474/`: the BASE replicate's shards, logs, `sacct.txt`, the comparison script and its output.
- `Manifest.toml`, `sacct.txt`, `shard_line_counts.txt`.

`Manifest.toml` and the SLURM logs match `.gitignore` patterns and were force-added (`git add -f`).
A later cleanup that re-applies `.gitignore` would drop them, so keep them when tidying this directory.

Reproduce the replicate comparison from the repository root (no Julia needed):

```bash
python3 benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/replicate-27438474/compare_replicate.py \
  benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/shards/base \
  benchmarking/results/td2rxh-evidence/td4mbg-27435671-27435672/replicate-27438474/shards \
  2 10
```

Its output should equal `replicate-27438474/comparison.txt`.

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
