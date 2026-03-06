#!/bin/bash
# Default SBATCH headers target Lawrencium `lr6`.
# Override on the `sbatch` CLI for other sites.
#SBATCH --job-name=mycelia_rhizo_suite
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=32
#SBATCH --output=rhizo_suite_%j.out
#SBATCH --error=rhizo_suite_%j.err

set -euo pipefail

detect_hpc_site() {
    if [[ -n "${MYCELIA_BENCHMARK_HPC_SITE:-}" ]]; then
        printf '%s\n' "${MYCELIA_BENCHMARK_HPC_SITE}"
        return
    fi

    case "${SLURM_JOB_PARTITION:-}" in
        lr*|es1)
            printf '%s\n' "lawrencium"
            ;;
        nih_*|normal|owners|serc)
            printf '%s\n' "scg"
            ;;
        *)
            printf '%s\n' "lawrencium"
            ;;
    esac
}

configure_site_environment() {
    local site="$1"

    module purge || true
    if [[ "${site}" == "lawrencium" ]]; then
        module load julia/1.10.10 || module load julia || module load julia/1.9.0 || true
    else
        module load julia/1.9.0 || module load julia || module load julia/1.10.10 || true
    fi
    module load bioconda || true
    source activate mycelia-bench || true
}

PROFILE="${1:-${BENCHMARK_PROFILE:-medium}}"
if [[ "${PROFILE}" != "medium" && "${PROFILE}" != "large" ]]; then
    echo "Profile must be 'medium' or 'large'. Got: ${PROFILE}" >&2
    exit 1
fi

HPC_SITE="$(detect_hpc_site)"

echo "=== Mycelia Rhizomorph Suite Launcher ==="
echo "HPC site: ${HPC_SITE}"
echo "Profile: ${PROFILE}"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Partition: ${SLURM_JOB_PARTITION:-unset}"
echo "Start time: $(date)"

configure_site_environment "${HPC_SITE}"

export JULIA_PROJECT=.
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
fi

mkdir -p results

if [[ "${PROFILE}" == "medium" ]]; then
    : "${BENCHMARK_SCALE:=medium}"
    : "${MYCELIA_BENCHMARK_MATRIX_STRATEGY:=tiered}"
    : "${MYCELIA_BENCHMARK_MAX_READS_PER_CONDITION:=50000}"
    : "${MYCELIA_BENCHMARK_MAX_RECORDS_PER_GENOME:=15000}"
    : "${MYCELIA_BENCHMARK_MAX_TOTAL_RECORDS:=120000}"
    : "${RHIZO_ISOLATE_TIMEOUT:=10h}"
    : "${RHIZO_METAGENOME_TIMEOUT:=12h}"
else
    : "${BENCHMARK_SCALE:=large}"
    : "${MYCELIA_BENCHMARK_MATRIX_STRATEGY:=tiered}"
    : "${MYCELIA_BENCHMARK_MAX_READS_PER_CONDITION:=100000}"
    : "${MYCELIA_BENCHMARK_MAX_RECORDS_PER_GENOME:=30000}"
    : "${MYCELIA_BENCHMARK_MAX_TOTAL_RECORDS:=250000}"
    : "${RHIZO_ISOLATE_TIMEOUT:=14h}"
    : "${RHIZO_METAGENOME_TIMEOUT:=16h}"
fi

export BENCHMARK_SCALE
export MYCELIA_BENCHMARK_MATRIX_STRATEGY
export MYCELIA_BENCHMARK_MAX_READS_PER_CONDITION
export MYCELIA_BENCHMARK_MAX_RECORDS_PER_GENOME
export MYCELIA_BENCHMARK_MAX_TOTAL_RECORDS

echo "JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-unset}"
echo "BENCHMARK_SCALE=${BENCHMARK_SCALE}"
echo "MYCELIA_BENCHMARK_MATRIX_STRATEGY=${MYCELIA_BENCHMARK_MATRIX_STRATEGY}"
echo "MYCELIA_BENCHMARK_MAX_READS_PER_CONDITION=${MYCELIA_BENCHMARK_MAX_READS_PER_CONDITION}"
echo "MYCELIA_BENCHMARK_MAX_RECORDS_PER_GENOME=${MYCELIA_BENCHMARK_MAX_RECORDS_PER_GENOME}"
echo "MYCELIA_BENCHMARK_MAX_TOTAL_RECORDS=${MYCELIA_BENCHMARK_MAX_TOTAL_RECORDS}"

echo "=== Running Isolate Suite (09) ==="
timeout "${RHIZO_ISOLATE_TIMEOUT}" julia --project=. benchmarking/09_rhizomorph_reconstruction_suite.jl \
    2>&1 | tee "results/09_rhizomorph_reconstruction_suite_${PROFILE}_$(date +%Y%m%d_%H%M%S).log"

echo "=== Running Metagenome Suite (10) ==="
timeout "${RHIZO_METAGENOME_TIMEOUT}" julia --project=. benchmarking/10_rhizomorph_metagenome_suite.jl \
    2>&1 | tee "results/10_rhizomorph_metagenome_suite_${PROFILE}_$(date +%Y%m%d_%H%M%S).log"

echo "End time: $(date)"
