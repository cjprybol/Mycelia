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
# Prerequisite: a runnable `julia` (1.10.x). On NERSC, load the cluster
# module first: `module load julia/1.10.10`. On Lawrencium, do NOT
# `module load julia/1.10.2-11.4` — it downgrades to a Julia that cannot load
# current master's Manifest (extension-trigger KeyError; td-j8bh). This
# script prefers juliaup's `+lts` channel (already on the host at
# ~/.juliaup/bin, resolving to 1.10.10 at last check) when present, invoked
# directly rather than via PATH + `command -v` — the latter would silently
# follow the host's mutable `juliaup default` instead of a pinned channel
# (td-j8bh round 2). Falls back to whatever `julia` the NERSC module put on
# PATH when juliaup is absent. On HPC, LD_LIBRARY_PATH is cleared here to
# avoid system libstdc++ conflicts.
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

# Prefer juliaup's pinned lts channel when present (Lawrencium): resolving
# via PATH + `command -v` would silently follow the host's mutable
# `juliaup default` rather than the version the four Lawrencium sbatch
# wrappers pin (td-j8bh round 2) -- this preflight prepares the SAME shared
# depot those jobs consume, so a skewed preflight version reintroduces the
# same Manifest-incompatibility class one script over. On NERSC there is no
# juliaup install; fall back to whatever `julia` the cluster's module put on
# PATH, unchanged from before.
if [[ -x "${HOME}/.juliaup/bin/julia" ]]; then
    JULIA="${HOME}/.juliaup/bin/julia"
    JULIA_ARGS=(+lts)
else
    JULIA="$(command -v julia || true)"
    JULIA_ARGS=()
fi
if [[ -z "${JULIA}" ]] || ! LD_LIBRARY_PATH="" "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
    echo "ERROR: julia not runnable. On NERSC, load the cluster's julia module first," >&2
    echo "       e.g. 'module load julia/1.10.10'. On Lawrencium, ensure juliaup's lts" >&2
    echo "       channel is installed (~/.juliaup/bin/julia +lts)." >&2
    exit 1
fi

if [[ -n "${submit_only}" ]]; then
    mode_desc="submit-only (${submit_only})"
else
    mode_desc="full preflight"
fi

echo "=== Mycelia HPC env preflight ==="
echo "julia:       ${JULIA} ${JULIA_ARGS[*]}  ($(LD_LIBRARY_PATH="" "${JULIA}" "${JULIA_ARGS[@]}" --version 2>&1 | head -1))"
echo "project:     ${project_dir}"
echo "depot:       ${JULIA_DEPOT_PATH:-<julia default>}"
echo "mode:        ${mode_desc}"
echo "start:       $(date)"

# resolve (repair stale/missing manifest) -> instantiate. Login-node-safe in
# BOTH modes ONLY with auto-precompile disabled: Pkg.instantiate() triggers
# precompilation by default, which is the step the Lawrencium login-node CPU
# watchdog SIGKILLs, not resolve/instantiate themselves (td-j8bh round 2).
JULIA_PKG_PRECOMPILE_AUTO=0 LD_LIBRARY_PATH="" "${JULIA}" "${JULIA_ARGS[@]}" --project="${project_dir}" -e '
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
LD_LIBRARY_PATH="" "${JULIA}" "${JULIA_ARGS[@]}" --project="${project_dir}" -e 'import Pkg; Pkg.precompile()'

echo "--- sanity: import Mycelia ---"
LD_LIBRARY_PATH="" "${JULIA}" "${JULIA_ARGS[@]}" --project="${project_dir}" -e 'import Mycelia; println("MYCELIA_IMPORT_OK")'

echo "=== preflight complete: $(date) ==="
