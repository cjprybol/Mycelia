#!/usr/bin/env bash
#
# hpc-setup.sh — Julia environment preflight for HPC benchmark runs.
#
# Resolves, instantiates, and (by default) precompiles the Mycelia project so
# that the first `import Mycelia` inside a batch job does not fail. This guards
# against the failure mode where a cluster checkout has a stale, untracked
# Manifest.toml (e.g. missing a dependency such as HDF5 that was added to
# Project.toml after the manifest was last built) — `Pkg.instantiate()` alone
# refuses to proceed, so `Pkg.resolve()` must run first to reconcile the
# manifest with Project.toml.
#
# Two modes:
#
#   (default)  resolve -> instantiate -> precompile -> import sanity check.
#              Run ONCE on a login node before submitting a benchmark job.
#              Precompilation is CPU-heavy and can take 10-30 min the first
#              time; do NOT bury it inside the (time-limited) batch window. Safe
#              on clusters whose login nodes permit long CPU jobs (e.g.
#              NERSC/Perlmutter).
#
#   --submit-only <sbatch-file>
#              resolve -> instantiate ONLY (cheap, login-node-safe), then
#              `sbatch <sbatch-file>` and print the job id — SKIPPING the
#              login-node precompile. Use this on clusters whose login nodes run
#              a CPU watchdog that SIGKILLs long-running precompiles (and the SSH
#              session carrying them), e.g. Lawrencium. The compute node then
#              precompiles on demand within the job's walltime, where no such
#              killer exists and the shared depot caches the result. The sbatch
#              file path is resolved relative to the current directory, and the
#              job is submitted from CWD.
#
# Prerequisite: a `julia` (1.10.x) must already be on PATH — load the cluster's
# module first, e.g. `module load julia/1.10.10` (NERSC) or
# `module load julia/1.10.2-11.4` (Lawrencium). On HPC, LD_LIBRARY_PATH is
# cleared here to avoid system libstdc++ conflicts.
#
# Usage:
#   benchmarking/hpc-setup.sh                                # full preflight
#   benchmarking/hpc-setup.sh --submit-only run_x.sbatch    # resolve+instantiate, then submit
#   benchmarking/hpc-setup.sh --help
#
set -euo pipefail

usage() {
    # Print the header comment block as help, skipping the shebang line.
    grep '^#' "$0" | grep -v '^#!' | sed 's/^# \{0,1\}//'
}

# --- argument parsing ---------------------------------------------------------
submit_only=""
case "${1:-}" in
    --help | -h)
        usage
        exit 0
        ;;
    --submit-only)
        submit_only="${2:-}"
        if [[ -z "${submit_only}" ]]; then
            echo "ERROR: --submit-only requires an sbatch file argument." >&2
            echo "       e.g. 'benchmarking/hpc-setup.sh --submit-only run_x.sbatch'." >&2
            exit 1
        fi
        ;;
    "")
        : # default: full preflight (resolve + instantiate + precompile)
        ;;
    *)
        echo "ERROR: unknown argument '${1}'. See --help." >&2
        exit 1
        ;;
esac

# Project root is the parent of this script's benchmarking/ directory.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"

if ! command -v julia >/dev/null 2>&1; then
    echo "ERROR: julia not found on PATH. Load the cluster's julia module first," >&2
    echo "       e.g. 'module load julia/1.10.10'." >&2
    exit 1
fi

echo "=== Mycelia HPC env preflight ==="
echo "julia:       $(command -v julia) ($(LD_LIBRARY_PATH="" julia --version))"
echo "project:     ${project_dir}"
echo "depot:       ${JULIA_DEPOT_PATH:-<julia default>}"
echo "mode:        ${submit_only:+submit-only ($submit_only)}${submit_only:-full preflight}"
echo "start:       $(date)"

# resolve (repair stale/missing manifest) -> instantiate. Cheap and
# login-node-safe in BOTH modes — this is never the step that trips a
# login-node CPU watchdog.
LD_LIBRARY_PATH="" julia --project="${project_dir}" -e '
    import Pkg
    Pkg.resolve()
    Pkg.instantiate()
'

# --- submit-only: skip login-node precompile, submit the job ------------------
if [[ -n "${submit_only}" ]]; then
    if [[ ! -f "${submit_only}" ]]; then
        echo "ERROR: sbatch file '${submit_only}' not found in $(pwd)." >&2
        exit 1
    fi
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "ERROR: sbatch not found on PATH (is this a SLURM login node?)." >&2
        exit 1
    fi
    echo "--- submit-only: skipping login-node precompile (compute node will"
    echo "    precompile on demand within walltime) ---"
    echo "--- submitting ${submit_only} ---"
    sbatch "${submit_only}"
    echo "=== submit-only complete: $(date) ==="
    exit 0
fi

# --- default: precompile (CPU-heavy) + import sanity --------------------------
echo "--- precompile ---"
LD_LIBRARY_PATH="" julia --project="${project_dir}" -e 'import Pkg; Pkg.precompile()'

echo "--- sanity: import Mycelia ---"
LD_LIBRARY_PATH="" julia --project="${project_dir}" -e 'import Mycelia; println("MYCELIA_IMPORT_OK")'

echo "=== preflight complete: $(date) ==="
