#!/usr/bin/env bash
#
# hpc-setup.sh — one-shot Julia environment preflight for HPC benchmark runs.
#
# Resolves, instantiates, and precompiles the Mycelia project so that the first
# `import Mycelia` inside a batch job does not fail. This guards against the
# failure mode where a cluster checkout has a stale, untracked Manifest.toml
# (e.g. missing a dependency such as HDF5 that was added to Project.toml after
# the manifest was last built) — `Pkg.instantiate()` alone refuses to proceed,
# so `Pkg.resolve()` must run first to reconcile the manifest with Project.toml.
#
# Run this ONCE on a login node before submitting a benchmark sbatch job.
# Precompilation is CPU-only and allowed on login nodes; it can take 10-30 min
# the first time, so do NOT bury it inside the (time-limited) batch window.
#
# Prerequisite: a `julia` (1.10.x) must already be on PATH — load the cluster's
# module first, e.g. `module load julia/1.10.10` (NERSC) or
# `module load julia/1.10.2-11.4` (Lawrencium). On HPC, LD_LIBRARY_PATH is
# cleared here to avoid system libstdc++ conflicts.
#
# Usage:
#   benchmarking/hpc-setup.sh            # set up the parent Mycelia project
#   benchmarking/hpc-setup.sh --help
#
set -euo pipefail

usage() {
    grep '^#' "$0" | sed 's/^# \{0,1\}//'
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    usage
    exit 0
fi

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
echo "start:       $(date)"

# resolve (repair stale/missing manifest) -> instantiate -> precompile.
LD_LIBRARY_PATH="" julia --project="${project_dir}" -e '
    import Pkg
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()
'

echo "--- sanity: import Mycelia ---"
LD_LIBRARY_PATH="" julia --project="${project_dir}" -e 'import Mycelia; println("MYCELIA_IMPORT_OK")'

echo "=== preflight complete: $(date) ==="
