# HPC Job Templates

Mycelia now supports a single canonical `JobSpec` that renders to multiple execution backends:

- SLURM batch scripts (`render_sbatch`)
- SLURM interactive allocation commands (`render_salloc`)
- Docker scripts (`render_docker_run`)
- Cloud Build YAML (`render_cloudbuild`)

## Why This Exists

The goal is one validated job definition that can target:

- NERSC Perlmutter
- LBL Lawrencium
- Stanford SCG
- local Docker
- Google Cloud Build

## Quick Example

```julia
job = Mycelia.JobSpec(
    job_name="example-job",
    cmd="julia pipeline.jl",
    site=:scg,
    partition="nih_s10",
    account="PI_SUNetID",
    time_limit="24:00:00",
    nodes=1,
    cpus_per_task=12,
    mem_gb=96
)

Mycelia.submit(job; dry_run=true)
```

## NERSC Notes

### Login Node Limits

NERSC login nodes are constrained by cgroups; do not run heavy compute directly on login nodes.

- Memory limit per user: 56 GB
- Aggregate CPU throttling: 25%

Use SLURM allocations for substantial compute.

### Constraint + GPU Request Rules

- CPU nodes: `constraint="cpu"`
- GPU nodes: `constraint="gpu"`
- Large-memory GPU nodes: `constraint="gpu&hbm80g"`
- Standard-memory GPU nodes can use `constraint="gpu&hbm40g"`

If you request a GPU constraint, you must also request GPUs (`gpus_per_node` or `gpus_per_task`).

### Interactive Jobs

Use `qos="interactive"` and render with `render_salloc`.

```julia
cmd = Mycelia.render_salloc(job)
```

For `gpu&hbm80g`, command rendering quotes the constraint correctly (`-C "gpu&hbm80g"`).

### Monitoring

Use these commands to observe jobs and resource usage:

- `squeue -u $USER`
- `scontrol show job <jobid>`
- `sacct -j <jobid> --format=JobID,JobName,AllocCPUS,Elapsed,State,ExitCode,MaxRSS`
- `sstat -j <jobid>.batch --format=JobID,MaxRSS,MaxVMSize,AvgCPU`
- `seff <jobid>` (if available)
- `nvidia-smi`
- `nvidia-smi topo -m`

`Mycelia.summarize_job(jobid)` provides a best-effort aggregated view.

## Lawrencium Notes

Lawrencium requires both partition and account. QoS must match your allowed associations.

Helpers:

```julia
Mycelia.list_lawrencium_associations()
Mycelia.list_lawrencium_qos_limits()
```

Association/QoS parsing is best-effort and produces warnings instead of hard failures when lookup is unavailable.

## SCG Notes

Supported partitions:

- `batch`
- `interactive`
- `nih_s10`

Rules encoded by validator:

- `nodes=1` expected for SCG
- `batch` and `nih_s10` require `account` (`PI_SUNetID` format)
- `time_limit` required

## Container Portability

### Docker

```julia
docker_script = Mycelia.render_docker_run(job)
```

### Cloud Build

```julia
yaml = Mycelia.render_cloudbuild(job)
```

Cloud Build rendering supports either:

- build + push from local `Dockerfile`, or
- prebuilt container image (`container_image`)

Cloud Build is CPU-oriented unless additional GPU infrastructure is configured separately.

### HPC Containers

When `container_image` is supplied, SLURM renderers can execute via site-appropriate engines
(defaults vary by site, override via `container_engine`).
