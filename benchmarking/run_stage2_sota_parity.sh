#!/usr/bin/env bash
# Run the Stage-2 corrector-validation benchmark on the REPEAT-BEARING E. coli
# K-12 fixture, emitting the NGA50 + misassembly-proxy contiguity metrics added
# for the SOTA-parity comparison (td-3e54).
#
# Repeat-free phiX/lambda make "contiguity parity" trivial; E. coli K-12's rRNA
# operons are where short-read assemblers fragment, so the arm-to-arm NGA50 /
# misassembly spread there is what makes parity a meaningful claim. Heavy
# (4.64 Mb @ 50x + assembly + dnadiff), so this is intended for an HPC node, not
# the CPU-only laptop.
#
# This is a thin, cluster-agnostic runner. The SLURM header below is a TEMPLATE:
# fill in --account / --partition / resource limits for the target cluster
# (Lawrencium or NERSC) before `sbatch`-ing it. Verify partition names with
# `sinfo` on the cluster first — do NOT assume. Run interactively with:
#   bash benchmarking/run_stage2_sota_parity.sh
#
# --- SLURM template (edit before use; values below are PLACEHOLDERS) ----------
# #SBATCH --job-name=rhizo-stage2-parity
# #SBATCH --account=<FILL_IN>
# #SBATCH --partition=<FILL_IN>
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=16
# #SBATCH --mem=64G
# #SBATCH --time=08:00:00
# #SBATCH --output=rhizo-stage2-parity-%j.log
# -----------------------------------------------------------------------------

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: run_stage2_sota_parity.sh [--targets LIST] [--coverage N] [--seed N] [--smoke]

Runs benchmarking/real_data_corrector_validation.jl against a repeat-bearing
fixture and writes a committed results bundle under benchmarking/results/.

Options (all have env-var equivalents consumed by the harness):
  --targets LIST   comma list of REGISTRY names (default: ecoli_k12)
  --coverage N     fold coverage for ART        (default: 50)
  --seed N         ART rndSeed                   (default: 42)
  --smoke          quick plumbing check (phix174 only, coverage 30)
  --no-sota        disable the standalone SPAdes SOTA arm (on by default here)
  --help           show this help

The harness self-verifies each reference length (fail-closed) and asserts every
arm produced complete dnadiff metrics before writing results. NO parity claim is
asserted by running this; it produces the evidence the claim would rest on.
USAGE
}

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGETS="ecoli_k12"
COVERAGE="50"
SEED="42"
SMOKE="false"
SOTA="true"   # this runner IS the SOTA-parity run; the independent SPAdes arm is on

while [[ $# -gt 0 ]]; do
    case "$1" in
        --targets)  TARGETS="$2"; shift 2 ;;
        --coverage) COVERAGE="$2"; shift 2 ;;
        --seed)     SEED="$2"; shift 2 ;;
        --smoke)    SMOKE="true"; shift ;;
        --no-sota)  SOTA="false"; shift ;;
        --help)     usage; exit 0 ;;
        *)          echo "unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

echo "Repo         : ${REPO_ROOT}"
echo "Targets      : ${TARGETS}"
echo "Coverage     : ${COVERAGE}x  seed=${SEED}  smoke=${SMOKE}  sota-arm=${SOTA}"
echo "Julia        : $(command -v julia)"

# LD_LIBRARY_PATH="" is required for Julia on HPC (see global CLAUDE.md).
cd "${REPO_ROOT}"
LD_LIBRARY_PATH="" \
MYCELIA_RDV_TARGETS="${TARGETS}" \
MYCELIA_RDV_COVERAGE="${COVERAGE}" \
MYCELIA_RDV_SEED="${SEED}" \
MYCELIA_RDV_SMOKE="${SMOKE}" \
MYCELIA_RDV_SOTA="${SOTA}" \
    julia --project="${REPO_ROOT}" \
        benchmarking/real_data_corrector_validation.jl
