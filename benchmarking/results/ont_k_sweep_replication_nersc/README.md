# ONT k-sweep — independent cross-host replication (NERSC)

Date: 2026-08-08.
Replicates `../ont_k_sweep/` on a different host with a different toolchain.

**This is a SEPARATE table. Do not merge it into `../ont_k_sweep/`.**
Its entire value is that it was produced independently; pooling the two
destroys the comparison and leaves one table of unclear provenance.

## Why this exists

Two open questions motivated it, and it answers both.

1. **No committed invocation reproduced the 240-cell grid** (`td-7cpm`). The
   driver's defaults yield 96 cells, not 240 (`td-4blm`), so the committed
   deliverable could not be regenerated from anything in the repo.
2. **Whether the measurements are host- and toolchain-stable at all.** The
   project had already discarded a set of laptop-computed cells on the grounds
   that a different Badread and a macOS QUAST made them incomparable. That
   reasoning was never tested against a full grid.

## The invocation

Reconstructed from the committed TSV, then confirmed by re-running it.
Values are **space-separated**; commas collapse the shard fan-out (`td-rrw6`).
Both arms run sequentially into one shared output directory.

```bash
OUTPUT_DIR=<dir> KS="11 13 15 17 19 21 31" COVERAGES="10 30 50 100" \
  TECHNOLOGIES="illumina ont" ORGANISMS="Lambda" \
  ./benchmarking/run_ont_k_sweep_shards.sh          # 168 cells

OUTPUT_DIR=<dir> KS="11 15 21 31" COVERAGES="10 30 50" \
  TECHNOLOGIES="illumina ont" ORGANISMS="T4" \
  ./benchmarking/run_ont_k_sweep_shards.sh          #  72 cells
```

Lambda carries a denser k ladder (7 values, adding 13/17/19) and T4 a narrower
coverage set (no 100x). Neither is expressible through the driver's defaults,
which is why nothing in the repo reproduced the grid.

## Provenance

The two rows that differ — Badread and Julia — are the whole point of the
comparison; everything else is held constant.

| component | this replication                    | committed table (`../ont_k_sweep/`) |
| --------- | ----------------------------------- | ----------------------------------- |
| host      | NERSC Perlmutter, Linux x86_64      | Lovelace, 112-core Linux x86_64     |
| Badread   | 0.4.1 (differs)                     | 0.4.2 (differs)                     |
| QUAST     | 5.3.0                               | 5.3.0                               |
| minimap2  | 2.31-r1302                          | —                                   |
| Julia     | 1.10.10 (differs)                   | 1.10.11 (differs)                   |
| commit    | `3479d1fbb5b2`                      | —                                   |
| SLURM job | `56466504`, COMPLETED 0:0, 2 h 11 m | —                                   |

`3479d1fbb5b2` differs from this branch's head only by a merge of master and an
artifact-harvest commit; `ont_k_sweep.jl`, `run_ont_k_sweep_shards.sh` and
`src/simulation.jl` are byte-identical across that range, so the measurement
path is the same code.

## Result

The grid matches: 240 cells, 240 complete checkpoints, Lambda 168 / T4 72,
illumina 120 / ont 120 — and 81-line summary, 46-line verdict stats, the same
as the committed files.

Compared cell-by-cell on the `(organism, technology, k, coverage, seed, arm)`
key, the key sets are identical and **20 of 21 columns are byte-identical
across all 240 cells**: `n_reads`, `n_contigs`, `asm_contigs_ge_min`,
`asm_max_contig`, `asm_total_bp`, `quast_contigs`, `total_length`, `N50`,
`largest_contig`, `NGA50`, `nga50_status`, `NA50`, `largest_alignment`,
`genome_fraction`, `duplication_ratio`, `misassemblies`, `unaligned_contigs`,
`unaligned_length`, `outcome`, `status`.

The one column that differs is `wall_seconds`, which is host-dependent by
definition and differs in all 240 rows.

## What this establishes, and what it does not

**Establishes.** The sweep is deterministic and reproducible across hosts. A
Badread **minor-version** difference (0.4.1 vs 0.4.2) and a Julia **patch**
difference (1.10.10 vs 1.10.11) produced zero measurable change across 240
cells. `../ont_k_sweep/README.md` already asserted Badread 0.4.1/0.4.2
output-equivalence, but rested it on two small files being unchanged across one
host switch; this is the same claim tested across the whole grid on a third
host.

**Does not establish.** This says nothing about the **macOS QUAST** half of the
rationale for discarding the earlier laptop-computed cells — that build was
forced to a single thread by a known Python 3.8+ bug, a variable untouched
here. It is therefore **not** grounds to reinstate those cells. It also does
not establish toolchain-independence in general: two minor/patch version
deltas held, which says nothing about a major version, a different aligner, or
a different simulator.

## Reproducing

Both arms into a fresh output directory, on a SLURM host — never the laptop
(see `protocols/compute-routing-lbnl-benchmarks.md`). Two prerequisites that
are not obvious and cost a day when missed (`td-j8bh`):

- **Instantiate the project against the host's depot first.** The dispatch
  recipe goes straight from `git worktree add` to `sbatch`; a branch carrying
  deps the depot lacks fails seconds in, writing its real error to
  `<output-dir>/shard-logs/prewarm.log` while the job's `.err` stays empty.
- **Do not `module load julia/1.10.2` on Lawrencium.** It is the only julia
  module there and it cannot load this Manifest — `import Mycelia` dies with
  `KeyError: "Graphs"` in `insert_extension_triggers`. Use juliaup's 1.10.10.

One driver defect to know when sharding: the serial pre-warm forwards `KS`,
`COVERAGES` and `ORGANISMS` from the environment but **not** `TECHNOLOGIES`, so
a run that narrows technologies gets extra default-technology cells written
into the same tree. A full run passing both technologies is unaffected, because
the pre-warm's cells are then a subset and dedupe via checkpoint resume.

Raw per-cell checkpoints live on NERSC scratch
(`$SCRATCH/ont_k_sweep_nersc_repro_20260807`), which auto-purges; only these
three aggregate files are durable.
