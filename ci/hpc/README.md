# HPC CI Runner

This directory contains the Mycelia HPC CI driver script for running full tests,
extended tutorials/benchmarks, and Codecov uploads on a compute node.

## Prerequisites
- Julia + git on the HPC node.
- Codecov token available as `CODECOV_TOKEN`.
- Optional: Codecov CLI (`codecov`) on PATH. If missing, the script falls back
  to Coverage.jl for upload.

## Basic Usage

1) Export your Codecov token (or load from your scheduler environment):

```bash
export CODECOV_TOKEN="your-token-here"
```

2) Submit the job to your scheduler (example for SLURM):

```bash
sbatch --time=04:00:00 --mem=32G --cpus-per-task=8 \
  --wrap="bash ci/hpc/run_hpc_ci.sh"
```

Adjust time/memory as needed for your cluster.

## Interactive SLURM Session

If you are already on a compute node via `salloc`/`srun`, run the driver directly:

```bash
export CODECOV_TOKEN="your-token-here"
bash ci/hpc/run_hpc_ci.sh
```

Use `--tests-only` or `--benchmarks-only` if you are limited on walltime or memory.

## Lawrencium sbatch Helper

This repo provides a Julia helper that calls the Mycelia Lawrencium wrapper:

```bash
export HPC_MAIL_USER="you@example.org"
export HPC_ACCOUNT="your_project"
export CODECOV_TOKEN="your-token-here"
julia --project=. ci/hpc/submit_lawrencium_hpc_ci.jl
```

Optional environment knobs (defaults shown):
- `HPC_JOB_NAME` = `mycelia-hpc-ci`
- `HPC_PARTITION` = `lr4`
- `HPC_QOS` = `lr_normal`
- `HPC_TIME` = `4-00:00:00`
- `HPC_CPUS_PER_TASK` = `8`
- `HPC_MEM_GB` = `32`
- `HPC_LOGDIR` = `~/workspace/slurmlogs`
- `HPC_CMD` = `bash <repo>/ci/hpc/run_hpc_ci.sh`

## Optional Environment Knobs
- `HPC_CODECOV_FLAG`: Override the Codecov flag (default: `hpc-extended`).
- `HPC_RESULTS_DIR`: Override output directory (default: `hpc-ci`).
- `HPC_CI_COMMIT`: Commit SHA to report in artifacts/Codecov.
- `HPC_CI_BRANCH`: Branch name to report in artifacts/Codecov.
- `HPC_CI_CHECKOUT=1`: Fetch/checkout `HPC_CI_COMMIT`/`HPC_CI_BRANCH` first.

## Outputs
- `hpc-ci/coverage/lcov.info`: LCOV coverage report.
- `hpc-ci/coverage/summary.json`: Coverage summary.
- `hpc-ci/logs/*.log`: Test/tutorial/benchmark logs.
- `hpc-ci/hpc-results.json`: JSON summary for downstream ingestion.
