#!/bin/bash
#SBATCH --job-name=mycelia_td2rxh
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --array=1-16
#SBATCH --time=1-00:00:00
#SBATCH --output=td2rxh_lrc_%A_%a.out
#SBATCH --error=td2rxh_lrc_%A_%a.err

# td-2rxh pre-registered evidence run on Lawrencium: does the merged sparse-rung
# indel staging (#420) plus the #450 reroute fix make the :nanopore corrector
# profile beat :illumina at 1000bp+? Runs
# indel_reroute_evidence.jl --cells td2rxh (6 kb x {1000,2000} bp x 20x x
# err {0.01,0.02,0.05}). One array task = one seed (1..16, per pre-registration
# amendment 1) over all six cells. Both arms of every pair run in the same
# process on an exclusive node; the harness runs GC before each timed arm and
# alternates arm order by seed parity.
#
# Usage (from benchmarking/ of the checkout to measure):
#   sbatch run_td2rxh_evidence_lrc.sbatch                    # seeds 1-16
#   sbatch --array=1 --export=ALL,TD2RXH_CELLS=smoke run_td2rxh_evidence_lrc.sbatch  # smoke
# Then analyse ONLY after every task exited 0 (see `sacct -j <ARRAY_JOB_ID>`):
#   julia +lts --project=.. td2rxh_paired_analysis.jl \
#     --csv "$(ls results/td2rxh/<ARRAY_JOB_ID>/seed*.csv | paste -sd, -)" \
#     --output-dir results/td2rxh/<ARRAY_JOB_ID>/analysis
# The analysis itself refuses mixed run labels and lists every missing pair.
#
# Prereq (once on a login node, same checkout): benchmarking/hpc-setup.sh.
# `set -e` is deliberately off so the row-count guard below can report; every
# command whose failure matters is checked explicitly instead.
set -uo pipefail
echo "=== td-2rxh evidence (Lawrencium) ==="
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Pinned juliaup lts launcher; see run_rgv_sweep_lrc.sbatch (td-j8bh) for why a
# PATH lookup or `module load julia` is not used.
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable on this host." >&2
  exit 1
fi

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
export LD_LIBRARY_PATH=""

# Measure the checkout this was SUBMITTED from, never a hardcoded path: the
# provenance below must describe the code that actually ran.
cd "${SLURM_SUBMIT_DIR}" || exit 1
if [[ ! -f indel_reroute_evidence.jl ]]; then
  echo "ERROR: submit from a checkout's benchmarking/ directory (got ${SLURM_SUBMIT_DIR})." >&2
  exit 1
fi
# Each git call is checked: a failing git (not a repo, safe.directory refusal)
# prints nothing, and an unchecked empty status would read as "clean".
if ! commit="$(git -C .. rev-parse HEAD)" || [[ -z "${commit}" ]]; then
  echo "ERROR: cannot resolve HEAD of $(cd .. && pwd)." >&2
  exit 1
fi
if ! status="$(git -C .. status --porcelain --untracked-files=no)"; then
  echo "ERROR: git status failed in $(cd .. && pwd)." >&2
  exit 1
fi
if [[ -n "${status}" ]]; then
  echo "ERROR: tracked files are modified; results would not be attributable to ${commit}:" >&2
  echo "${status}" >&2
  exit 1
fi
# Manifest.toml is untracked, so dependency versions are not covered by HEAD.
# Fold its hash into the run label so shards resolved against different
# dependency sets cannot be pooled (the analysis rejects mixed labels).
if [[ ! -f ../Manifest.toml ]]; then
  echo "ERROR: ../Manifest.toml missing; run benchmarking/hpc-setup.sh first." >&2
  exit 1
fi
if ! manifest_sha="$(sha256sum ../Manifest.toml | cut -c1-12)" ||
   [[ ! "${manifest_sha}" =~ ^[0-9a-f]{12}$ ]]; then
  echo "ERROR: could not hash ../Manifest.toml (got '${manifest_sha:-}')." >&2
  exit 1
fi
label="td2rxh@${commit:0:12}+manifest${manifest_sha}"
echo "Checkout: $(cd .. && pwd)  label=${label}"

# Mycelia's Stage-1 corrector needs a temp root that (a) does not resolve through
# a symlink (_owned_assembly_artifact_identity) and (b) supports
# renameat2(RENAME_NOREPLACE) for its fail-closed cleanup. Lawrencium's default
# TMPDIR fails both: it sits under the symlinked $HOME/workspace on Lustre scratch,
# where RENAME_NOREPLACE returns EINVAL (smoke jobs 27116288, 27118360). So every
# task gets a private directory on node-local disk, resolved to its physical path,
# and a network/parallel filesystem is refused up front instead of failing later.
tmp_base="${TD2RXH_TMP_BASE:-${SLURM_TMPDIR:-/tmp}}"
# mktemp -d creates a fresh 0700 directory and refuses an existing path, so a
# requeued task never inherits stale contents and a pre-planted symlink at a
# predictable name cannot redirect the cleanup below.
if ! tmp_root="$(mktemp -d "${tmp_base}/td2rxh_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")"; then
  echo "ERROR: cannot create a TMPDIR under ${tmp_base}." >&2
  exit 1
fi
# Cleanup is armed immediately. The directory is KEPT only when julia or the row
# check failed (rc != 0 below), because Mycelia's temp artifacts are the evidence
# for why. It is REMOVED on a clean run, on an early guard exit, and on a signal
# (walltime/scancel SIGTERM reaches the trap with rc still 0); a SIGKILL after
# KillWait skips the trap and leaves it on the node.
rc=0
trap 'if [[ "${rc}" == "0" ]]; then rm -rf "${tmp_root}"; else echo "kept TMPDIR ${tmp_root} for debugging" >&2; fi' EXIT
if ! TMPDIR="$(cd "${tmp_root}" && pwd -P)" || [[ -z "${TMPDIR}" ]]; then
  echo "ERROR: cannot resolve TMPDIR ${tmp_root}." >&2
  exit 1
fi
export TMPDIR
if ! tmp_fs="$(stat -f -c %T "${TMPDIR}")"; then
  echo "ERROR: cannot determine the filesystem type of ${TMPDIR}." >&2
  exit 1
fi
# Allowlist, not denylist: an unrecognised type (BeeGFS prints "UNKNOWN (0x...)")
# must not slip through to fail mid-run.
case "${tmp_fs}" in
  ext2/ext3|xfs|tmpfs|btrfs|overlayfs) ;;
  *)
    echo "ERROR: TMPDIR ${TMPDIR} is on ${tmp_fs}; Mycelia needs a local filesystem (ext/xfs/tmpfs/btrfs/overlayfs)." >&2
    exit 1 ;;
esac
echo "TMPDIR=${TMPDIR} (${tmp_fs})"

cells="${TD2RXH_CELLS:-td2rxh}"
case "${cells}" in
  td2rxh) expected=12 ;;   # 6 cells x 1 seed x 2 arms
  smoke)  expected=4 ;;    # 2 cells x 1 seed x 2 arms
  all)    expected=6 ;;    # 3 cells x 1 seed x 2 arms
  *) echo "ERROR: unknown TD2RXH_CELLS=${cells}" >&2; exit 1 ;;
esac
seed="${SLURM_ARRAY_TASK_ID}"
outdir="results/td2rxh/${SLURM_ARRAY_JOB_ID}"
mkdir -p "${outdir}" || exit 1
out="${outdir}/seed${seed}.csv"
# A requeued task must not inherit the previous attempt's rows.
rm -f "${out}"

echo ">>> RUN cells=${cells} seed=${seed} -> ${out}  $(date)"
if "${JULIA}" "${JULIA_ARGS[@]}" --project=.. indel_reroute_evidence.jl \
    --cells "${cells}" --seeds "${seed}" --out "${out}" --label "${label}"; then
  julia_rc=0
else
  julia_rc=$?
fi

# A green exit with too few rows means some arms threw (each failure is a @warn
# + skip in the harness). Treat that as a failure, not a result.
rows=0
[[ -f "${out}" ]] && rows=$(( $(wc -l < "${out}") - 1 ))
echo "julia_rc=${julia_rc} rows=${rows} expected=${expected}"
rc=0
if [[ "${julia_rc}" -ne 0 ]]; then
  echo ">>> FAIL: julia exited ${julia_rc} (see .err)" >&2
  rc=1
elif [[ "${rows}" -lt "${expected}" ]]; then
  echo ">>> FAIL: ${out} has ${rows} of ${expected} arm rows; some arms threw (see .err)" >&2
  rc=1
fi
echo "=== done rc=${rc} $(date) ==="
exit "${rc}"
