# gqc v0.1 — reproducible HPC runs

How to run `Mycelia.run_gqc` on Lawrencium, NERSC (Perlmutter), and Lovelace
on the **post-v0.1 stacked fixes** branch (based on `cursor/gqc-v0.1` @
`f2ece978`, plus GraphML escaping, collision-free seeds, unified node TSV).
Pin to the consolidated PR head, not a floating tip.

Companion API/spec: [`docs/gqc-v0.1-spec.md`](gqc-v0.1-spec.md).

## Julia version and pin

| Item | Value |
|------|-------|
| `Project.toml` `[compat] julia` | `"1.10"` |
| `Manifest.toml` `julia_version` | `1.10.12` (when a Manifest is present) |
| CI matrix | `julia-version: ['lts']` (Julia 1.10 LTS line) |

Use **Julia 1.10.x** (prefer 1.10.12 to match the Manifest). Instantiate **at the
pinned Mycelia commit** for the run (consolidated stack head, or `f2ece978` only
if you intentionally accept legacy seeds):

```bash
git clone https://github.com/cjprybol/Mycelia.git
cd Mycelia
git checkout <consolidated-stack-sha>   # or f2ece978 for #472-only
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

Keep `JULIA_PROJECT` / `--project=` pointed at that checkout for the job.

## Headless plots

`run_gqc(...; make_plots=true)` uses Plots.jl. On batch nodes set:

```bash
export GKSwstype=100
```

(or `ENV["GKSwstype"]="100"` inside the Julia `-e` string). Without this, plot
creation can fail or hang on HPC.

## Nulls and references (required before trusting q-values)

Defaults (`src/gqc.jl`):

- `null_source=:auto` — if any reference records load, nulls use **reference**;
  else segments (k-mer-node graphs then set `null_applicable=false`)
- Reference dir: `<graph_dir>/reference/`
- Extensions: `.fa`, `.fas`, `.fna`, `.ffn`, `.fasta`, optional `.gz`
- Recommended for HVP-style k-mer graphs: put genomes under `reference/*.fasta`
  and keep `null_source=:auto` (or `:reference`)

**Before trusting any dinucleotide / BH q-value**, open the packet
`manifest.json` and confirm:

| Field | Required for trusted nulls |
|-------|----------------------------|
| `params.resolved_null_source` | `"reference"` |
| `params.null_applicable` | `true` |
| `params.null_source` | `"auto"` or `"reference"` (what you requested) |
| `params.seed_scheme` | `"gqc-null-v1-sha256"` on this stack (see below) |

Summary rows also carry `null_source` / `null_applicable`. If
`resolved_null_source` is `"segments"` or `null_applicable` is `false`, do
**not** treat p/q as paper-style reference nulls.

Cross-genus FDR is **not** computed inside `run_gqc` (within-graph BH only;
column `dinucleotide_bh_q_within_graph`). Use a post-pass
(`Mycelia.gqc_cross_bh` or your own BH across packets).

### Seed scheme (differs from `f2ece978`)

This stack uses **`gqc-null-v1-sha256`**: one `Xoshiro` per
`(shuffle_seed, salt, genus_index, null_kind, replicate)` via
`Mycelia.gqc_derive_null_seed`. Null numeric outputs **differ** from packets
produced at `f2ece978` (legacy `shuffle_seed + 1000*salt + genus_index`). Reruns
under the same scheme remain bit-identical. Confirm `params.seed_scheme` in the
manifest before comparing to older HVP tables.

## Threads and memory

`run_gqc` does **not** call `Threads.@threads` for null shuffles; each replicate
gets its own seeded RNG and runs sequentially. Extra CPUs help Julia/BLAS and
Plots more than the null loop itself.

Practical starting points (tune with `seff` / `sacct`):

| Site | Suggested starting point |
|------|---------------------------|
| Lawrencium | `cpus_per_task=8`–`16` (`LRC_THREADS` default in `lawrencium_sbatch` is **16**), `mem_gb=64` |
| NERSC shared | `nersc_sbatch_shared` with modest `cpus_per_task` and `mem_gb ≈ 2×cpus` (helper default) |
| NERSC regular | Avoid login-node compute; use `nersc_sbatch` with `constraint="cpu"` |
| Lovelace | Dedicated node via `lovelace_run` (no SLURM); size to graph-set RAM |

Large graph-sets: prefer scratch for `output_root`; use `overwrite=:skip` to
resume after partial failures (`gqc_graph_set_errors.jsonl`).

## Recommended `run_gqc` call

```julia
import Mycelia
ENV["GKSwstype"] = "100"

Mycelia.run_gqc(
    "/path/to/graph_set_or_single_graph",
    "/scratch/gqc_out";
    ks=[5, 11, 21],
    bh_ks=[11],
    n_shuffles=100,              # GQC_N_SHUFFLES default
    shuffle_seed=20261021,       # GQC_SHUFFLE_SEED default
    null_source=:auto,
    strand_mode=:canonical,      # prefer for DNA k-mer-node HVP graphs (I6)
    overwrite=:skip,
    make_plots=true,
)
```

Single-graph folder: `graph.gfa` + optional `metadata.json` + optional
`reference/`. Graph-set: immediate children each containing `graph.gfa`.

## Lawrencium (SLURM)

Helper: `Mycelia.lawrencium_sbatch` (`src/slurm-sbatch.jl`). Defaults include
`partition="lr6"`, `qos="lr_normal"`, `cpus_per_task=LRC_THREADS` (16),
`mem_gb=64`. Requires account (`account=` or `ENV["LRC_ACCOUNT"]`) and mail
user (`mail_user=`, `ENV["SLURM_MAIL_USER"]`, or `git config user.email`).

```julia
import Mycelia
ENV["GKSwstype"] = "100"

mycelia = "/path/to/Mycelia"   # pinned consolidated-stack SHA
cmd = """
julia --project=$(mycelia) -e '
ENV["GKSwstype"]="100"
import Mycelia
Mycelia.run_gqc("/path/to/graph_set", "/scratch/gqc_out";
    ks=[5, 11, 21], bh_ks=[11], n_shuffles=100,
    null_source=:auto, strand_mode=:canonical,
    overwrite=:skip, make_plots=true)
'
"""

Mycelia.lawrencium_sbatch(;
    job_name="gqc-v0.1",
    cmd=cmd,
    time="0-04:00:00",
    cpus_per_task=8,
    mem_gb=64,
    dry_run=true,   # set false to submit
)
```

## NERSC Perlmutter (SLURM)

Helpers: `Mycelia.nersc_sbatch` / `Mycelia.nersc_sbatch_shared`. **Account
required** (`account=` or `ENV["NERSC_ACCOUNT"]`). Do not run heavy gqc work on
login nodes (cgroup limits). Prefer `constraint="cpu"` for this CPU packet.

```julia
import Mycelia
ENV["GKSwstype"] = "100"

mycelia = "/path/to/Mycelia"   # pinned consolidated-stack SHA
cmd = """
julia --project=$(mycelia) -e '
ENV["GKSwstype"]="100"
import Mycelia
Mycelia.run_gqc("/path/to/graph_set", "/pscratch/gqc_out";
    ks=[5, 11, 21], bh_ks=[11], n_shuffles=100,
    null_source=:auto, strand_mode=:canonical,
    overwrite=:skip, make_plots=true)
'
"""

Mycelia.nersc_sbatch_shared(;
    job_name="gqc-v0.1",
    mail_user="you@example.com",
    cmd=cmd,
    time="0-04:00:00",
    cpus_per_task=8,
    mem_gb=32,
    constraint="cpu",
    dry_run=true,
)
```

## Lovelace (no SLURM)

Helper: `Mycelia.lovelace_run` — runs `bash -c cmd` on the dedicated host and
writes timestamped stdout/stderr under `~/workspace/slurmlogs` by default.

```julia
import Mycelia
ENV["GKSwstype"] = "100"

mycelia = "/path/to/Mycelia"   # pinned consolidated-stack SHA
cmd = """
julia --project=$(mycelia) -e '
ENV["GKSwstype"]="100"
import Mycelia
Mycelia.run_gqc("/path/to/graph_set", "/data/gqc_out";
    ks=[5, 11, 21], bh_ks=[11], n_shuffles=100,
    null_source=:auto, strand_mode=:canonical,
    overwrite=:skip, make_plots=true)
'
"""

Mycelia.lovelace_run(; job_name="gqc-v0.1", cmd=cmd)
```

## Post-run checklist

1. Packet(s) under `output_root/<graph_id_slug>/`
2. `manifest.json`: `resolved_null_source == "reference"`, `null_applicable == true`,
   `seed_scheme == "gqc-null-v1-sha256"`
3. `summary.csv` / `summary.json` for metrics; `nulls/` for replicate tables
4. Node TSVs (`nodes/<slug>.tsv` and `formats/nodes.tsv`) share one header
   (`GQC_NODE_TSV_COLUMNS`); GraphML IDs/attrs are XML-escaped
5. Graph-set: inspect `output_root/gqc_graph_set_errors.jsonl` if present
6. Re-run with the same commit + `(shuffle_seed, salt, genus_index)` for
   bit-identical nulls under `gqc-null-v1-sha256`

## See also

- [`docs/gqc-v0.1-spec.md`](gqc-v0.1-spec.md) — input layout, null_source, I6
- [`docs/src/hpc-job-templates.md`](src/hpc-job-templates.md) — JobSpec / site notes
