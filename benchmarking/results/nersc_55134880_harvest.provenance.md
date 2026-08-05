# NERSC job 55134880 — harvest and reconciliation

Harvested 2026-08-05 from Perlmutter, resolving the open question of whether
`rhizomorph_correction_validation_sweep_20260711_125013.csv` was the harvest of
job 55134880 or a distinct run.

## Answer

They are unrelated. Job 55134880 is `mycelia_talk_bench`; the 2026-07-11 CSV is
a separate, later correction-validation sweep. The two share no artifacts.

## Job facts

Confirmed via `sacct` on 2026-08-05:

| Field | 55134880 | 55170557 |
| --- | --- | --- |
| Job name | `mycelia_talk_bench` | `mycelia-pr331` |
| State | COMPLETED (0:0) | COMPLETED (0:0) |
| Start | 2026-06-27T08:52:17 | 2026-06-28T01:58:24 |
| End | 2026-06-27T09:02:08 | 2026-06-28T02:08:30 |
| Elapsed | 00:09:51 | 00:10:06 |
| Node | nid006360 | — |

Source path on Perlmutter:
`/global/cfs/cdirs/m4269/cjprybol/Mycelia/benchmarking/`

## Artifacts harvested

Four result JSONs, all previously absent from this repo:

- `kmer_analysis_benchmark_2026-06-27_08-55-36.json`
- `kmer_analysis_comprehensive_2026-06-27_08-55-43.json`
- `data_processing_benchmark_2026-06-27_09-02-05.json`
- `data_processing_enhanced_2026-06-27_09-02-05.json`

Measured on AMD EPYC 7763 (64-core), Julia 1.10.10, 540 GB total memory.
K-mer analysis covered k ∈ {11, 15, 19, 21} across genome sizes
{10000, 50000, 100000}, total sequence length 800000. Data processing covered
5 files / 4.24 MB at an estimated 32.58 MB/s, peak 177.99 MB.

These are hardware- and date-specific measurements, not deterministically
regenerable output, so they are committed rather than treated as cache.

## Caveats carried forward

The job's stderr recorded two conditions that do not invalidate the timings
above but do bound their interpretation:

1. **MPI was not available.** `libmpi_gnu_123.so` failed to load
   (`libfabric.so.1: cannot open shared object file`). These numbers are
   single-node, non-MPI measurements.
2. **Memory-fit estimates were nonsensical** — the run logged requests of
   32.000 TiB and 32768.000 PiB against 478 GiB free, from
   `src/utility-functions.jl:740`. The benchmark completed regardless, but the
   estimator itself looks defective and is tracked separately.

## Not harvested

Older runs of the same suites remain on Perlmutter and are unharvested
(2026-06-22, 2026-06-20). They are noted here so their absence is a recorded
decision rather than a silent drop.
