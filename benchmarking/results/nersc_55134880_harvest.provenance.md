# NERSC benchmark harvest — job 55134880 and the runs beside it

Harvested 2026-08-05 from Perlmutter, resolving the open question of whether
`rhizomorph_correction_validation_sweep_20260711_125013.csv` was the harvest of
job 55134880 or a distinct run.

## Answer

They are unrelated. Job 55134880 is `mycelia_talk_bench`; the 2026-07-11 CSV is
a separate, later correction-validation sweep with a different schema
(`reference,genome_len,error_rate,regime,readlen,tech,…`). The two share no
artifacts, and the CSV does not exist anywhere under the Mycelia tree on
Perlmutter.

The strongest evidence is the job's own stdout, which ends with a
self-inventory (`--- result JSONs ---`) naming exactly the four 2026-06-27
files below.

## Job facts

Confirmed via `sacct` on 2026-08-05:

| Field | 55134880 | 55170557 |
| --- | --- | --- |
| Job name | `mycelia_talk_bench` | `mycelia-pr331` |
| State | COMPLETED (0:0) | COMPLETED (0:0) |
| Start | 2026-06-27T08:52:17 | 2026-06-28T01:58:24 |
| End | 2026-06-27T09:02:08 | 2026-06-28T02:08:30 |
| Elapsed | 00:09:51 | 00:10:06 |
| Node | nid006360 | nid004408 |

Source directory on Perlmutter (also the job's `WorkDir`):
`/global/cfs/cdirs/m4269/cjprybol/Mycelia/benchmarking/`

## What was harvested

**All twenty result JSONs** of these two suites present on CFS — five runs,
none of which had ever been tracked in this repo:

| Run | Files | Note |
| --- | --- | --- |
| 2026-06-19 12:01 / 12:09 | 4 | |
| 2026-06-19 12:33 / 12:44 | 4 | |
| 2026-06-20 07:43 / 07:51 | 4 | |
| 2026-06-22 10:41 / 10:52 | 4 | |
| 2026-06-27 08:55 / 09:02 | 4 | job 55134880 |

Verified byte-identical to the Perlmutter originals (per-file `md5`, aggregate
`086dbdd0a485f2d940ad48dce82274bc`), and confirmed absent from all prior git
history via `git log --all --diff-filter=A`.

An earlier draft of this note harvested only the 2026-06-27 run and listed the
unharvested remainder as "2026-06-22, 2026-06-20." That list was **wrong**: it
was derived from the job stdout's `--- result JSONs ---` block, which truncates
after ten entries and therefore could not show the two 2026-06-19 runs. Four
runs were unharvested, not two — so the accompanying sentence claiming their
absence was "a recorded decision rather than a silent drop" was false for
exactly the runs it failed to name. Harvesting all five removes the arbitrary
boundary rather than restating it. Provenance notes are only worth having if
they are checked against the directory, not against another artifact's summary.

## Measurement conditions

Measured on AMD EPYC 7763 (64-core), Julia 1.10.10, 540 GB total memory.
K-mer analysis covered k ∈ {11, 15, 19, 21} across genome sizes
{10000, 50000, 100000}, total sequence length 800000. Data processing covered
5 files / 4.24 MB at an estimated 32.58 MB/s, peak 177.99 MB.

## How to read these numbers

They are **smoke-test-grade, not benchmark-grade**, and should not be cited as
performance results without re-measurement:

1. **Timings are dominated by JIT warm-up, not by algorithmic cost.** In
   `kmer_analysis_benchmark`, `kmer_counting_k11_size50000` reports
   `min 11.4 ms`, `median 11.7 ms`, `mean 183 ms`, `max 527 ms`,
   `std 298 ms` — a standard deviation ~25× the median. The same shape recurs
   across every `kmer_counting_*` and `count_records_file1`. Sample counts are
   small.
2. **MPI was not available.** `libmpi_gnu_123.so` failed to load
   (`libfabric.so.1: cannot open shared object file`). These are single-node,
   non-MPI measurements.
3. **The memory-fit estimator was emitting nonsense.** The run logged requests
   of 32.000 TiB and 32768.000 PiB against 478 GiB free, from
   `src/utility-functions.jl:740` (still the `@warn "Matrix may not fit in
   available memory!"` line in the current tree). The benchmark completed
   regardless, but the estimator looks defective and is tracked separately.
4. **`k_mer_throughput_per_second` covers only k ∈ {11, 15, 19}.** k=21 is in
   `k_values_tested` but absent from the throughput block, so that series is
   partial.

They are committed because they are hardware- and date-specific observations
that cannot be reconstructed from the repo — not because they are strong
measurements. `benchmarking/results/.gitignore` excludes
`*_perkmer_*.jsonl` and `stage0_calibration_*.json`; the repo's root
`.gitignore` additionally excludes several bulky `benchmarking/results/`
subtrees. None of those patterns match these files, confirmed via
`git check-ignore`.
