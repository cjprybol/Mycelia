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
# Relative to the submit directory; with the documented usage (submit from
# benchmarking/) the logs land in benchmarking/, not the repo root.
#SBATCH --output=ecoli_qualmer_memgate_%j.out
#SBATCH --error=ecoli_qualmer_memgate_%j.err

# td-n8ax E. coli qualmer-path MEMORY gate on lr6 (~95GB node — the node the
# un-fixed runs OOM'd on at 74-200GB). GATE: run COMPLETES and peak RSS <= 32GB
# (aspiration 16GB) with the COMPOSED config (:ultralight_quality aggregate storage
# + qualmer_prefilter_min_count=2). Driver: benchmarking/ecoli_qualmer_memory_gate.jl.
#
# Usage (submit from the checkout's benchmarking/ directory; that checkout is the
# default run dir and the Slurm logs land beside this script):
#   cd ~/workspace/Mycelia/benchmarking && sbatch run_ecoli_qualmer_memory_gate_lrc.sbatch
# Overrides (sbatch exports the submitting environment):
#   QUALMER_GATE_RUN_DIR=<checkout>   checkout whose code runs (default: the git
#                                     toplevel of the submit directory)
#   QUALMER_GATE_MANIFEST=resolve     (default) resolve the checkout's dependencies
#                                     from scratch into a private environment
#   QUALMER_GATE_MANIFEST=existing    pin every version in <checkout>/Manifest.toml
#   QUALMER_GATE_MANIFEST=<path>      pin every version in that Manifest.toml
#                                     (relative paths are taken from the submit dir)
#   QUALMER_GATE_TMP_BASE=<dir>       node-local parent for TMPDIR (default
#                                     ${SLURM_TMPDIR:-/tmp})
# Driver knobs (MYCELIA_ECOLI_COVERAGE, MYCELIA_MEMORY_CEILING_GB, ...) pass through.
#
# PRIVATE ENVIRONMENT. The job never writes the checkout's Manifest.toml: other
# jobs read it concurrently (run_td2rxh_evidence_lrc.sbatch hashes it into its run
# label). Instead it builds ${TMPDIR}/env, whose Project.toml carries the
# checkout's [deps] and [compat] plus Mycelia developed from the checkout by path,
# so the code that runs is the checkout's src/. The previous version copied the
# main checkout's Manifest.toml into a hardcoded worktree; that manifest went stale
# when Project.toml gained XML and moved UMAP 0.1 -> 0.3, and precompile failed.
# CAUTION: resolving from scratch once hit a Pkg extension bug (KeyError "Graphs",
# job 24000519, July 2026). That has not been shown fixed; if it recurs, rerun with
# QUALMER_GATE_MANIFEST=existing or a known-good Manifest path.
# The sha256 of the Manifest actually used and the checkout HEAD are printed and
# written into the driver's results CSV (git_head, manifest_sha256,
# manifest_src_sha256); the env Project/Manifest are saved to results/<jobid>/.
#
# A dirty tracked tree is reported in the log but does NOT stop the run: this is a
# memory gate on whatever is checked out, not a provenance-locked evidence run.
#
# Pkg steps run with auto-precompile off; the `import Mycelia` warmup precompiles
# on the compute node (the Lawrencium login-node watchdog kills long precompiles;
# see benchmarking/hpc-setup.sh). Mycelia bootstraps ART + QUAST via Bioconda into
# the group depot on first use.
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
# Environment hygiene, as in run_td2rxh_evidence_lrc.sbatch: sbatch exports the
# submit shell, so a stray MYCELIA_PROJECT / JULIA_PROJECT could select another
# environment, and Julia's default load path would stack the global environment.
unset MYCELIA_PROJECT JULIA_PROJECT
export JULIA_LOAD_PATH="@:@stdlib"
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
# Made absolute and physical before any cd; a relative override is taken from the
# job's starting directory, which is the submit directory.
if [[ -n "${QUALMER_GATE_RUN_DIR:-}" ]]; then
  if ! run_dir="$(cd "${QUALMER_GATE_RUN_DIR}" && pwd -P)" || [[ -z "${run_dir}" ]]; then
    echo "ERROR: cannot resolve QUALMER_GATE_RUN_DIR='${QUALMER_GATE_RUN_DIR}'." >&2
    exit 1
  fi
else
  if ! top="$(git -C "${SLURM_SUBMIT_DIR}" rev-parse --show-toplevel)" || [[ -z "${top}" ]] ||
     ! run_dir="$(cd "${top}" && pwd -P)" || [[ -z "${run_dir}" ]]; then
    echo "ERROR: ${SLURM_SUBMIT_DIR} is not inside a git checkout; submit from the checkout or set QUALMER_GATE_RUN_DIR." >&2
    exit 1
  fi
fi
if [[ ! -f "${run_dir}/benchmarking/ecoli_qualmer_memory_gate.jl" ]] || [[ ! -f "${run_dir}/Project.toml" ]]; then
  echo "ERROR: ${run_dir} is not a Mycelia checkout with benchmarking/ecoli_qualmer_memory_gate.jl." >&2
  exit 1
fi
if ! head="$(git -C "${run_dir}" rev-parse HEAD)" || [[ -z "${head}" ]]; then
  echo "ERROR: cannot resolve HEAD of ${run_dir}." >&2
  exit 1
fi
echo "Run dir: ${run_dir}  HEAD: ${head}"
echo "Uncommitted tracked changes (empty if clean; reported, not refused):"
git -C "${run_dir}" status --short --untracked-files=no 2>&1 || true

# --- Manifest source, made absolute before any cd ---------------------------
manifest_mode="${QUALMER_GATE_MANIFEST:-resolve}"
manifest_src=""
case "${manifest_mode}" in
  resolve) ;;
  existing) manifest_src="${run_dir}/Manifest.toml" ;;
  *)
    manifest_path="${manifest_mode}"
    [[ "${manifest_path}" == /* ]] || manifest_path="${SLURM_SUBMIT_DIR}/${manifest_path}"
    if [[ ! -f "${manifest_path}" ]] ||
       ! manifest_dir="$(cd "$(dirname "${manifest_path}")" && pwd -P)"; then
      echo "ERROR: QUALMER_GATE_MANIFEST='${manifest_mode}' is not resolve, existing, or an existing file." >&2
      exit 1
    fi
    manifest_src="${manifest_dir}/$(basename "${manifest_path}")"
    if [[ "${manifest_src}" == "${run_dir}/Manifest.toml" ]]; then
      manifest_mode=existing
    fi ;;
esac
if [[ -n "${manifest_src}" ]] && [[ ! -f "${manifest_src}" ]]; then
  echo "ERROR: QUALMER_GATE_MANIFEST=${manifest_mode} but ${manifest_src} does not exist." >&2
  exit 1
fi

cd "${run_dir}/benchmarking" || exit 1

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
# On every exit the private env's Project.toml + Manifest.toml (the only record
# of the resolved dependency versions; most Mycelia deps carry no [compat]) and
# the driver's contigs FASTA and QUAST report (if any) are copied to
# benchmarking/results/<jobid>/. TMPDIR is then REMOVED unless keep_tmp=1, which
# is set only for a setup failure (Pkg step, manifest hash, warmup) or a driver
# CODE error (exit 2) -- the cases whose scratch is evidence. A gate PASS (0) or
# memory FAIL (1) is an expected result and its multi-GB scratch is removed, as
# it is on an early guard exit or a walltime/scancel SIGTERM. A SIGKILL after
# KillWait skips the trap.
status=0
keep_tmp=0
env_dir=""
save_dir="${run_dir}/benchmarking/results/${SLURM_JOB_ID}"
# shellcheck disable=SC2329  # invoked from the EXIT trap below
save_artifacts() {
  local d f
  mkdir -p "${save_dir}" || { echo "WARN: cannot create ${save_dir}" >&2; return; }
  if [[ -n "${env_dir}" ]]; then
    for f in Project.toml Manifest.toml; do
      [[ -f "${env_dir}/${f}" ]] && cp -p "${env_dir}/${f}" "${save_dir}/env_${f}"
    done
  fi
  for d in "${tmp_root}"/ecoli_qualmer_memory_gate_*/; do
    [[ -d "${d}" ]] || continue
    cp -p "${d}"*_contigs.fasta "${save_dir}/" 2>/dev/null || true
    [[ -d "${d}quast" ]] && cp -pR "${d}quast" "${save_dir}/" 2>/dev/null
  done
  echo "saved run artifacts to ${save_dir}"
}
trap 'save_artifacts; if [[ "${keep_tmp}" == "1" ]]; then echo "kept TMPDIR ${tmp_root} on $(hostname) for debugging" >&2; else rm -rf "${tmp_root}"; fi' EXIT
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

# --- Private environment ------------------------------------------------------
env_dir="${TMPDIR}/env"
mkdir -p "${env_dir}" || exit 1
if [[ -n "${manifest_src}" ]] && ! cp -f "${manifest_src}" "${env_dir}/Manifest.toml"; then
  echo "ERROR: cannot copy ${manifest_src} into ${env_dir}." >&2
  exit 1
fi
# Env Project.toml = the checkout's [deps] (so the driver's direct imports
# resolve) + [compat] for those deps, then Mycelia developed from the checkout by
# path. With a pinned manifest, PRESERVE_ALL keeps every version it records.
setup_jl="${TMPDIR}/gate_env_setup.jl"
cat > "${setup_jl}" <<'JL'
import Pkg, TOML
run_dir, env_dir, mode = ARGS[1], ARGS[2], ARGS[3]
proj = TOML.parsefile(joinpath(run_dir, "Project.toml"))
deps = get(proj, "deps", Dict{String, Any}())
compat = Dict{String, Any}(k => v for (k, v) in get(proj, "compat", Dict{String, Any}())
                           if k == "julia" || haskey(deps, k))
open(io -> TOML.print(io, Dict{String, Any}("deps" => deps, "compat" => compat)),
     joinpath(env_dir, "Project.toml"), "w")
Pkg.activate(env_dir)
if mode == "resolve"
    Pkg.develop(path = run_dir)
else
    Pkg.develop(path = run_dir; preserve = Pkg.PRESERVE_ALL)
end
Pkg.instantiate()
JL
# Two hashes, not interchangeable: manifest_src_sha256 is the Manifest as supplied
# (comparable with a checkout's Manifest.toml; empty in resolve mode), while
# manifest_sha256 is the env Manifest that actually ran, which includes the
# Mycelia develop-path entry and so never matches a checkout's file.
manifest_src_sha=""
if [[ -n "${manifest_src}" ]]; then
  if ! manifest_src_sha="$(sha256sum "${manifest_src}" | cut -d' ' -f1)" ||
     [[ ! "${manifest_src_sha}" =~ ^[0-9a-f]{64}$ ]]; then
    echo "ERROR: could not hash ${manifest_src} (got '${manifest_src_sha:-}')." >&2
    exit 1
  fi
fi
echo ">>> Pkg (${manifest_mode}${manifest_src:+ from ${manifest_src}}) into ${env_dir} $(date)"
if ! JULIA_PKG_PRECOMPILE_AUTO=0 "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no \
    "${setup_jl}" "${run_dir}" "${env_dir}" "${manifest_mode}"; then
  echo ">>> Pkg ${manifest_mode} FAILED" >&2
  status=1; keep_tmp=1
  exit 1
fi
if ! manifest_sha="$(sha256sum "${env_dir}/Manifest.toml" | cut -d' ' -f1)" ||
   [[ ! "${manifest_sha}" =~ ^[0-9a-f]{64}$ ]]; then
  echo "ERROR: could not hash ${env_dir}/Manifest.toml (got '${manifest_sha:-}')." >&2
  status=1; keep_tmp=1
  exit 1
fi
# QUALMER_GATE_HEAD also tells the driver it is wrapped, so it leaves its workdir
# for the EXIT trap to save instead of deleting it at process exit.
export QUALMER_GATE_HEAD="${head}"
export QUALMER_GATE_MANIFEST_SHA="${manifest_sha}"
export QUALMER_GATE_MANIFEST_SRC_SHA="${manifest_src_sha}"
echo "Provenance: HEAD=${QUALMER_GATE_HEAD}  manifest_sha256=${QUALMER_GATE_MANIFEST_SHA}  manifest_src_sha256=${QUALMER_GATE_MANIFEST_SRC_SHA:-<resolve mode>}"

echo ">>> precompile warmup (import Mycelia) $(date)"
if ! "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project="${env_dir}" \
    -e 'import Mycelia; println("Mycelia loaded from ", pathof(Mycelia))'; then
  echo ">>> precompile/load FAILED" >&2
  status=1; keep_tmp=1
  exit 1
fi

echo ">>> RUN ecoli_qualmer_memory_gate.jl $(date)"
"${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project="${env_dir}" ecoli_qualmer_memory_gate.jl
status=$?
# Exit 2 is a driver CODE error: keep its scratch. 0 and 1 are gate results.
[[ "${status}" -eq 2 ]] && keep_tmp=1

echo "=== gate exit=${status} complete $(date) ==="
echo "Slurm-side peak (after the job ends): sacct -j ${SLURM_JOB_ID} --format=JobID,State,Elapsed,MaxRSS"
# exit 0 = PASS (completes, <=32GB); 1 = memory FAIL (OOM/over-ceiling); 2 = code error.
exit "${status}"
