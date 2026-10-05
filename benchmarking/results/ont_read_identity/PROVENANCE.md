# Provenance of the committed read-identity tables

The three TSVs in this directory were produced by one SLURM job on Lawrencium on
2026-10-05, using `benchmarking/run_ont_read_identity_lrc.sbatch` with the
default parameters of `benchmarking/ont_read_identity.jl` (Lambda NC_001416,
simulated ONT reads at 30x, seed 42).

Verbatim lines from that job's log (`td-p9sh-26687382.out`, not committed):

```text
Job: 26687382  Node: n0348.lr6  CPUs: 4  Start: Mon Oct  5 07:46:54 PDT 2026
Repo: /global/home/users/cjprybol/workspace/Mycelia/.worktrees/td-p9sh  HEAD: 236bffae57a365cc9daf97f93c631f2f14620ecb
JULIA: /global/home/users/cjprybol/.juliaup/bin/julia +lts  (julia version 1.10.10)
  version:                Badread v0.4.1
  default --identity:     95,99,2.5
  default --error_model:  nanopore2023
  default --length:       15000,13000
>>> OK  ont_read_identity.jl
=== complete: Mon Oct  5 07:49:24 PDT 2026 ===
```

`sacct -j 26687382` reported `COMPLETED`, exit code `0:0`, elapsed `00:02:32`.

## Comparison with the previous (Lovelace, Badread 0.4.2) tables

Checked column by column against the tables these replaced:

- `per_read_identity.tsv`: all 139 rows and all 8 previously present columns are
  identical. The only new column is `read_length`.
- `read_identity_summary.tsv`: every quantile, `n`, `n_reads_total`,
  `n_unmapped`, `coverage` and `seed` are identical. The recorded
  `badread_version` reads v0.4.1 instead of v0.4.2. New columns: `read_source`,
  `reads_file`, `total_read_bases`, `scored_read_bases`,
  `aligned_base_fraction`.
- `kmer_survival_ladder.tsv`: `p_error_free_kmer` and `error_free_kmer_coverage`
  are identical. New columns: `aligned_base_fraction`,
  `p_error_free_kmer_all_reads`, `error_free_kmer_coverage_all_reads`.

So the simulated read set is unchanged, and the new values are the read-length
accounting and the quantities derived from it.
