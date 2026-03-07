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

## Publish Results Branch

After `run_hpc_ci.sh` writes `hpc-ci/hpc-results.json`, publish a lightweight
results branch with badge endpoints:

```bash
bash ci/hpc/publish_hpc_results.sh
```

This updates the `hpc-results` branch with:
- `latest-hpc-results.json`: raw Phase 1/2 summary emitted by `run_hpc_ci.sh`
- `latest-tests.json`: Shields endpoint JSON for the HPC test badge
- `latest-benchmarks.json`: Shields endpoint JSON for the HPC benchmark badge
- `latest-meta.json`: commit/timestamp/cluster metadata
- `<commit>/...`: archived copies for each published HPC run

Run `bash ci/hpc/publish_hpc_results.sh --no-push` to stage the branch locally
for inspection before pushing.

To make this part of the same HPC job, export `HPC_RESULTS_PUBLISH=1` before
running `run_hpc_ci.sh`. The driver will publish the `hpc-results` branch after
writing `hpc-ci/hpc-results.json`.

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
- `HPC_RESULTS_PUBLISH=1`: Automatically run `publish_hpc_results.sh` after summary generation.
- `HPC_RESULTS_BRANCH`: Branch used by `publish_hpc_results.sh` (default: `hpc-results`).
- `HPC_RESULTS_REMOTE`: Remote used by `publish_hpc_results.sh` (default: `origin`).
- `HPC_RESULTS_REPO_URL`: Explicit repo URL for cloning/pushing the results branch.
- `HPC_RESULTS_REPOSITORY_SLUG`: Explicit `owner/repo` for badge URLs if remote parsing is insufficient.

## Outputs
- `hpc-ci/coverage/lcov.info`: LCOV coverage report.
- `hpc-ci/coverage/summary.json`: Coverage summary.
- `hpc-ci/logs/*.log`: Test/tutorial/benchmark logs.
- `hpc-ci/hpc-results.json`: JSON summary for downstream ingestion.
- `hpc-results/latest-*.json`: published badge endpoints and latest run metadata.
