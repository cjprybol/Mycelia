#!/bin/bash
#SBATCH --job-name=ecoli_qualmer_memgate
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
# lr_normal MaxWall is 3 days. The previous 12h wall was known-insufficient: job
# 24024728 held a flat ~34.9 GB peak but was CANCELLED at 12:00:26 by the wall
# before the corrector finished, so no verdict was printed.
#SBATCH --time=3-00:00:00
#SBATCH --output=ecoli_qualmer_memgate_%j.out
#SBATCH --error=ecoli_qualmer_memgate_%j.err

# td-n8ax E. coli qualmer-path MEMORY gate on lr6 (~95GB node — the node the
# un-fixed runs OOM'd on at 74-200GB). GATE: run COMPLETES and peak RSS <= 32GB
# (aspiration 16GB) with the COMPOSED config (:ultralight_quality aggregate storage
# + qualmer_prefilter_min_count=2). Driver: benchmarking/ecoli_qualmer_memory_gate.jl.
#
# Usage (from anywhere inside the checkout to measure; it is the default run dir):
#   cd ~/workspace/Mycelia && sbatch benchmarking/run_ecoli_qualmer_memory_gate_lrc.sbatch
# Overrides (sbatch exports the submitting environment):
#   QUALMER_GATE_RUN_DIR=<checkout>   checkout whose code runs (default: the git
#                                     toplevel of the submit directory)
#   QUALMER_GATE_MANIFEST=resolve     (default) Pkg.resolve() + instantiate against
#                                     <checkout>/Project.toml inside the job
#   QUALMER_GATE_MANIFEST=existing    use <checkout>/Manifest.toml as-is (instantiate only)
#   QUALMER_GATE_MANIFEST=<path>      copy that Manifest.toml in, then instantiate only
#   QUALMER_GATE_TMP_BASE=<dir>       node-local parent for TMPDIR (default
#                                     ${SLURM_TMPDIR:-/tmp})
# Driver knobs (MYCELIA_ECOLI_COVERAGE, MYCELIA_MEMORY_CEILING_GB, ...) pass through.
# Example:
#   QUALMER_GATE_RUN_DIR=$HOME/workspace/Mycelia QUALMER_GATE_MANIFEST=resolve \
#     sbatch benchmarking/run_ecoli_qualmer_memory_gate_lrc.sbatch
#
# Why resolve by default: the previous version copied the MAIN checkout's
# Manifest.toml into a hardcoded worktree. That manifest goes stale whenever
# Project.toml gains or moves a dependency, and precompile then fails. Resolving
# against the checkout's own Project.toml keeps the two in step. Resolve and
# instantiate run with auto-precompile off; the `import Mycelia` warmup below
# precompiles on the compute node (the Lawrencium login-node watchdog kills long
# precompiles; see benchmarking/hpc-setup.sh).
#
# Mycelia bootstraps ART + QUAST via Bioconda into the group depot on first use.
# `-e` is deliberately omitted: every step that matters is guarded explicitly.
set -uo pipefail
echo "=== E. coli qualmer memory gate (Lawrencium lr6) ==="
echo "Job ${SLURM_JOB_ID}  Node ${SLURM_NODELIST}  Start $(date)"

# julia/1.10.2-11.4 DOWNGRADES to a Julia that cannot load current master's
# Manifest (extension-trigger KeyError; td-j8bh). Invoke the juliaup launcher
# DIRECTLY with `+lts` rather than PATH + `command -v` -- the latter silently
# falls through to whatever julia is already on the inherited PATH if juliaup
# is missing on a node, and even when it finds the shim, the shim resolves to
# the host's mutable `juliaup default` -- not necessarily lts. `+lts` is
# resolved by juliaup itself, pinned regardless of the default (td-j8bh review).
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable (juliaup lts channel missing/broken on this host)" >&2
  exit 1
fi
echo "JULIA: ${JULIA} ${JULIA_ARGS[*]}  ($("${JULIA}" "${JULIA_ARGS[@]}" --version 2>&1 | head -1))"

# Group depot (precompiled deps + Bioconda envs install here).
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
export LD_LIBRARY_PATH=""
# Julia's thread count is DERIVED from the allocation, never hardcoded. The old
# `${SLURM_CPUS_PER_TASK:-16}` default was a silent oversubscription hazard: with
# no --cpus-per-task directive Slurm allocates 1 CPU per task, the variable is
# unset, and Julia would still spawn 16 threads on 1 core — enough contention to
# turn the completion gate into a false timeout. --cpus-per-task=16 above is the
# single source of truth for both Slurm and Julia; if it is ever missing the job
# fails loudly here rather than silently oversubscribing.
if [[ -z "${SLURM_CPUS_PER_TASK:-}" ]]; then
  echo "ERROR: SLURM_CPUS_PER_TASK is unset — the --cpus-per-task directive is missing." >&2
  echo "       Refusing to guess a thread count (oversubscription risks a false timeout)." >&2
  exit 1
fi
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}  JULIA_NUM_THREADS=${JULIA_NUM_THREADS}"

# --- Run directory: the checkout whose code is measured ----------------------
# Default is the checkout this was submitted from, never a hardcoded worktree
# (the old `cd .worktrees/td-n8ax-gate` pinned every run to a stale commit).
if [[ -n "${QUALMER_GATE_RUN_DIR:-}" ]]; then
  run_dir="${QUALMER_GATE_RUN_DIR}"
elif ! run_dir="$(git -C "${SLURM_SUBMIT_DIR}" rev-parse --show-toplevel)" || [[ -z "${run_dir}" ]]; then
  echo "ERROR: ${SLURM_SUBMIT_DIR} is not inside a git checkout; submit from the checkout or set QUALMER_GATE_RUN_DIR." >&2
  exit 1
fi
if [[ ! -f "${run_dir}/benchmarking/ecoli_qualmer_memory_gate.jl" ]] || [[ ! -f "${run_dir}/Project.toml" ]]; then
  echo "ERROR: ${run_dir} is not a Mycelia checkout with benchmarking/ecoli_qualmer_memory_gate.jl." >&2
  exit 1
fi
cd "${run_dir}/benchmarking" || exit 1
# Provenance for the log (informational): what code actually ran.
echo "Run dir: ${run_dir}  HEAD: $(git -C "${run_dir}" rev-parse HEAD 2>&1)"
echo "Uncommitted tracked changes (empty if clean):"
git -C "${run_dir}" status --short --untracked-files=no 2>&1 || true

# --- Node-local TMPDIR --------------------------------------------------------
# The driver's mktempdir() and Mycelia's iterative corrector both follow TMPDIR.
# Lawrencium's default TMPDIR resolves onto Lustre scratch, where the corrector's
# renameat2(RENAME_NOREPLACE) cleanup returns EINVAL (td-2o7h; smoke jobs
# 27116288, 27118360), and its artifact-identity check rejects a symlinked temp
# root. Same construction as run_td2rxh_evidence_lrc.sbatch: a fresh private
# directory on node-local disk, resolved to its physical path, with a local
# filesystem required up front.
tmp_base="${QUALMER_GATE_TMP_BASE:-${SLURM_TMPDIR:-/tmp}}"
# mktemp -d creates a fresh 0700 directory and refuses an existing path, so a
# requeued job never inherits stale contents.
if ! tmp_root="$(mktemp -d "${tmp_base}/ecoli_qualmer_memgate_${SLURM_JOB_ID}.XXXXXX")"; then
  echo "ERROR: cannot create a TMPDIR under ${tmp_base}." >&2
  exit 1
fi
# Always removed on exit (clean run, failure, or walltime/scancel SIGTERM): the
# E. coli reads and graph scratch are large and node-local disk is shared with
# later jobs. Results are written to benchmarking/results/, not here. A SIGKILL
# after KillWait skips the trap and leaves the directory on the node.
trap 'rm -rf "${tmp_root}"' EXIT
if ! TMPDIR="$(cd "${tmp_root}" && pwd -P)" || [[ -z "${TMPDIR}" ]]; then
  echo "ERROR: cannot resolve TMPDIR ${tmp_root}." >&2
  exit 1
fi
export TMPDIR
if ! tmp_fs="$(stat -f -c %T "${TMPDIR}")"; then
  echo "ERROR: cannot determine the filesystem type of ${TMPDIR}." >&2
  exit 1
fi
# Allowlist, not denylist: an unrecognised type must not slip through to fail mid-run.
case "${tmp_fs}" in
  ext2/ext3|xfs|tmpfs|btrfs|overlayfs) ;;
  *)
    echo "ERROR: TMPDIR ${TMPDIR} is on ${tmp_fs}; Mycelia needs a local filesystem (ext/xfs/tmpfs/btrfs/overlayfs)." >&2
    exit 1 ;;
esac
echo "TMPDIR=${TMPDIR} (${tmp_fs}, $(df -h "${TMPDIR}" | awk 'NR==2 {print $4 " free"}'))"
# Confirm Julia itself resolves tempdir() to the node-local directory, not a
# value baked in elsewhere (startup file, inherited JULIA_* settings).
if ! julia_tmp="$("${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no -e 'print(tempdir())')" ||
   [[ "${julia_tmp%/}" != "${TMPDIR%/}" ]]; then
  echo "ERROR: Julia tempdir() is '${julia_tmp:-}', expected ${TMPDIR}." >&2
  exit 1
fi
echo "Julia tempdir(): ${julia_tmp}"

# --- Manifest: resolve against this checkout's Project.toml by default --------
manifest_mode="${QUALMER_GATE_MANIFEST:-resolve}"
case "${manifest_mode}" in
  resolve)
    pkg_expr='import Pkg; Pkg.resolve(); Pkg.instantiate()' ;;
  existing)
    if [[ ! -f ../Manifest.toml ]]; then
      echo "ERROR: QUALMER_GATE_MANIFEST=existing but ${run_dir}/Manifest.toml does not exist." >&2
      exit 1
    fi
    pkg_expr='import Pkg; Pkg.instantiate()' ;;
  *)
    if [[ ! -f "${manifest_mode}" ]]; then
      echo "ERROR: QUALMER_GATE_MANIFEST='${manifest_mode}' is not resolve, existing, or an existing file." >&2
      exit 1
    fi
    if ! cp -f "${manifest_mode}" ../Manifest.toml; then
      echo "ERROR: cannot copy ${manifest_mode} to ${run_dir}/Manifest.toml." >&2
      exit 1
    fi
    echo ">>> copied Manifest.toml from ${manifest_mode}"
    pkg_expr='import Pkg; Pkg.instantiate()' ;;
esac
echo ">>> Pkg (${manifest_mode}) against ${run_dir}/Project.toml $(date)"
if ! JULIA_PKG_PRECOMPILE_AUTO=0 "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project=.. -e "${pkg_expr}"; then
  echo ">>> Pkg ${manifest_mode} FAILED" >&2
  exit 1
fi

echo ">>> precompile warmup (import Mycelia) $(date)"
if ! "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project=.. -e 'import Mycelia; println("Mycelia loaded")'; then
  echo ">>> precompile/load FAILED" >&2
  exit 1
fi

echo ">>> RUN ecoli_qualmer_memory_gate.jl $(date)"
"${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project=.. ecoli_qualmer_memory_gate.jl
status=$?

echo "=== gate exit=${status} complete $(date) ==="
echo "Slurm-side peak (after the job ends): sacct -j ${SLURM_JOB_ID} --format=JobID,State,Elapsed,MaxRSS"
# exit 0 = PASS (completes, <=32GB); 1 = memory FAIL (OOM/over-ceiling); 2 = code error.
exit "${status}"
