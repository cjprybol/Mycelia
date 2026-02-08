# SLURM Template Migration

This document maps legacy Mycelia SLURM usage onto the new `JobSpec` + template renderer.

## What Changed

- Added a canonical `JobSpec` abstraction for site-aware job definitions.
- Added file-backed templates under `templates/` for SLURM, Docker, and Cloud Build.
- Added strict validation (NERSC/Lawrencium/SCG policy checks) before submission.
- Existing public APIs (`lawrencium_sbatch`, `scg_sbatch`, `nersc_sbatch`, `submit_job`) still exist and now call the new backend.

## Old Usage -> New Usage

### SCG batch or NIHS10

Old:

```julia
Mycelia.scg_sbatch(
    job_name="my-job",
    mail_user="user@example.org",
    account="PI_SUNetID",
    partition="nih_s10",
    cmd="julia run.jl"
)
```

New canonical form:

```julia
job = Mycelia.JobSpec(
    job_name="my-job",
    cmd="julia run.jl",
    site=:scg,
    partition="nih_s10",
    account="PI_SUNetID",
    time_limit="24:00:00",
    nodes=1,
    cpus_per_task=8,
    mem_gb=64
)

Mycelia.submit(job; dry_run=true)   # render only
Mycelia.submit(job; dry_run=false)  # submit
```

### Lawrencium

```julia
job = Mycelia.JobSpec(
    job_name="lawrencium-index",
    cmd="julia build_index.jl",
    site=:lawrencium,
    partition="lr6",
    qos="lr_normal",
    account="pc_example",
    time_limit="12:00:00",
    cpus_per_task=16,
    mem_gb=64
)

Mycelia.submit(job; dry_run=true)
```

### NERSC GPU (1 task / 1 GPU)

```julia
job = Mycelia.JobSpec(
    job_name="pm-gpu",
    cmd="python train.py",
    site=:nersc,
    account="m1234_g",
    qos="regular",
    constraint="gpu",
    gpus_per_node=1,
    ntasks=1,
    cpus_per_task=16,
    mem_gb=64,
    time_limit="04:00:00"
)

Mycelia.submit(job; dry_run=true)
```

### Docker local run

```julia
job = Mycelia.JobSpec(
    job_name="local-dev",
    cmd="python script.py --input /data/input.tsv",
    site=:local,
    time_limit="01:00:00",
    container_image="ghcr.io/example/my-image:latest",
    container_mounts=["$(pwd()):/workspace"],
    container_workdir="/workspace"
)

docker_script = Mycelia.render_docker_run(job)
```

### Cloud Build render

```julia
job = Mycelia.JobSpec(
    job_name="ci-job",
    cmd="pytest -q",
    site=:cloudbuild,
    time_limit="01:00:00",
    container_image="gcr.io/my-project/my-image:sha256-deadbeef"
)

yaml = Mycelia.render_cloudbuild(job)
```

## Serialization

```julia
json_text = Mycelia.job_spec_to_json(job)
Mycelia.write_job_spec_json(job, "jobs/my-job.json")
loaded = Mycelia.job_spec_from_json("jobs/my-job.json")
```

## Dry Run and Auditing

- `submit(...; dry_run=true)` prints the rendered artifact and exact submit command.
- `write_sbatch(job)` writes scripts for review/commit before submission.

## Add a New Queue Template

1. Add a template file under `templates/slurm/<site>/`.
2. Extend template selection logic in `src/slurm-templates.jl` (`_selected_template` and site-specific helpers).
3. Add or update policy rules in `validate`.
4. Add fixture + test coverage under `test/`.
