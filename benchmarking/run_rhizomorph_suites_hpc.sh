#!/bin/bash
#SBATCH --job-name=mycelia_rhizo_suite
#SBATCH --partition=compute
#SBATCH --time=36:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --output=rhizo_suite_%j.out
#SBATCH --error=rhizo_suite_%j.err

set -euo pipefail

PROFILE="${1:-${BENCHMARK_PROFILE:-medium}}"
if [[ "${PROFILE}" != "medium" && "${PROFILE}" != "large" ]]; then
    echo "Profile must be 'medium' or 'large'. Got: ${PROFILE}" >&2
    exit 1
fi

echo "=== Mycelia Rhizomorph Suite Launcher ==="
echo "Profile: ${PROFILE}"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"

module purge || true
module load julia/1.9.0 || true
module load bioconda || true
source activate mycelia-bench || true

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
