# ONT k-sweep — independent cross-host replication (NERSC)

Date: 2026-08-08.
Replicates `../ont_k_sweep/` on a different host with a different toolchain.

**This is a SEPARATE table. Do not merge it into `../ont_k_sweep/`.**
Its entire value is that it was produced independently; pooling the two
destroys the comparison and leaves one table of unclear provenance.

## Why this exists

Two open questions motivated it, and it answers both.

1. **No committed invocation reproduced the 240-cell grid** (`td-7cpm`). At the
   time, the driver's defaults yielded 96 cells, not 240 (`td-4blm`), so the
   committed deliverable could not be regenerated from anything in the repo.
   `../ont_k_sweep/README.md` now describes reproducing it as separate Lambda
   and T4 runs into one output directory.
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
which is why nothing in the repo reproduced the grid at the time.

## Provenance

Apart from the host, the two rows that differ — Badread and Julia — are the
whole point of the comparison. QUAST matches; versions not recorded for the
committed table ("—") are unknown, not shown to be held constant.

| component | this replication                    | committed table (`../ont_k_sweep/`) |
| --------- | ----------------------------------- | ----------------------------------- |
| host      | NERSC Perlmutter, Linux x86_64      | Lovelace, 112-core Linux x86_64     |
| Badread   | 0.4.1 (differs)                     | 0.4.2 (differs)                     |
| QUAST     | 5.3.0                               | 5.3.0                               |
| minimap2  | 2.31-r1302                          | —                                   |
| Julia     | 1.10.10 (differs)                   | 1.10.11 (differs)                   |
| commit    | `3479d1fbb5b2`                      | —                                   |
| SLURM job | `56466504`, COMPLETED 0:0, 2 h 11 m | —                                   |

`3479d1fbb5b2` is the code the replication ran. It is not on master; it is
kept on the branch `cp/ont-sweep-nersc-replication` so this provenance stays
resolvable. Between it and that branch's head there are only merges of master
and the commit adding these results, and `ont_k_sweep.jl`,
`run_ont_k_sweep_shards.sh` and `src/simulation.jl` are byte-identical across
that range, so the measurement path is the same code.

## Result

The grid matches: 240 cells, 240 complete checkpoints, Lambda 168 / T4 72,
illumina 120 / ont 120. The 46-line verdict stats are byte-identical to the
committed file; the 81-line summary matches it in every column except
`median_wall_seconds`.

Compared cell-by-cell on the
`(organism, technology, k, coverage, seed, decoder_arm)` key, the key sets are
identical and **21 of the 22 non-key columns are byte-identical across all 240
cells**: `accession`, `n_reads`, `n_contigs`, `asm_contigs_ge_min`,
`asm_max_contig`, `asm_total_bp`, `quast_contigs`, `total_length`, `N50`,
`largest_contig`, `NGA50`, `nga50_status`, `NA50`, `largest_alignment`,
`genome_fraction`, `duplication_ratio`, `misassemblies`, `unaligned_contigs`,
`unaligned_length`, `outcome`, `status`.

The one column that differs is `wall_seconds`, which is host-dependent by
definition and differs in all 240 rows.

## What this establishes, and what it does not

**Establishes.** Across these two hosts the sweep is deterministic and
reproducible. Two **patch-level** differences — Badread 0.4.1 vs 0.4.2 and
Julia 1.10.10 vs 1.10.11 — produced zero measurable change across 240 cells.
`../ont_k_sweep/README.md` asserts Badread 0.4.1/0.4.2 output-equivalence
from smaller checks; this tests the same claim across the whole grid on
another host.

**Does not establish.** This says nothing about the **macOS QUAST** half of the
rationale for discarding the earlier laptop-computed cells — that build was
forced to a single thread by a known Python 3.8+ bug, a variable untouched
here. It is therefore **not** grounds to reinstate those cells. It also does
not establish toolchain-independence in general: two patch-level version
deltas held, which says nothing about a minor or major version, a different
aligner, or a different simulator.

## Reproducing

Both arms into a fresh output directory, on a SLURM host — never the laptop.
Two prerequisites that are not obvious and cost a day when missed (`td-j8bh`):

- **Instantiate the project against the host's depot first.** The dispatch
  recipe goes straight from `git worktree add` to `sbatch`; a branch carrying
  deps the depot lacks fails seconds in, writing its real error to
  `<output-dir>/shard-logs/prewarm.log` while the job's `.err` stays empty.
- **Do not `module load julia/1.10.2` on Lawrencium.** It is the only julia
  module there and it cannot load this Manifest — `import Mycelia` dies with
  `KeyError: "Graphs"` in `insert_extension_triggers`. Use juliaup's 1.10.10.

One driver behaviour to know when sharding: the serial pre-warm runs a fixed
cell (`--technologies illumina,ont --ks 31 --coverages 10 --seeds 42`) and
forwards only `ORGANISMS` from the environment. That is true at
`3479d1fbb5b2` and on master. It writes up to two cells per organism, one per
technology, and the final `--aggregate-only` pass reads the whole `cells/`
tree. So a run whose grid excludes either technology, k = 31, 10x or seed 42
also reports those pre-warm cells. Both runs above include all four, so their
pre-warm cells are a subset and dedupe via checkpoint resume.

Raw per-cell checkpoints live on NERSC scratch
(`$SCRATCH/ont_k_sweep_nersc_repro_20260807`), which auto-purges; only these
three aggregate files are durable.
