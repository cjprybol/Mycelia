# Frozen: pre-k-stamping B8 run

This directory holds the B8 Viterbi accuracy run as it existed before `k` became
a per-run parameter. **No code path writes here any more**, and no invocation of
`benchmarking/viterbi_accuracy_benchmark.jl` can reproduce it.

| | |
| --- | --- |
| `k` | **11** (not recorded in the tables; recoverable only as `editable_positions / (observation_count - 2)`) |
| `benchmark_run_id` | `b8_viterbi_accuracy_local_20260625` — no `_k` suffix |
| `benchmark_git_commit` | `806f0748` |
| `benchmark_k` column | absent |
| `metadata.k` | absent |

## Why it cannot be reproduced

The generator now stamps `k` into the run id, the metadata, a per-row
`benchmark_k` column, and the output directory. A run at k=11 therefore emits
run id `b8_viterbi_accuracy_local_20260625_k11` into
`benchmarking/results/viterbi_accuracy_b8_k11/`. Nothing produces the unsuffixed
run id above.

The directory is kept rather than regenerated because it is the artifact several
external references cite by path and by run id; freezing keeps every one of them
resolving. It is superseded for all new work.

## The gap this leaves

The `benchmark_k` column exists because provenance sidecars do not survive a
hand-copy across a repo boundary while a column does. That reasoning does not
reach these tables: they predate the column, and since no command regenerates
them, they can never acquire it in place. Any consumer that needs `k` alongside
these rows must either read it from this file or re-harvest from a k-stamped
run.

Regenerating into `viterbi_accuracy_b8_k11/` and re-pointing the consumers is
the real fix. It is a benchmark execution, so it routes to SLURM rather than a
laptop, and it is deliberately not bundled into the review-fix change that
created this note.
