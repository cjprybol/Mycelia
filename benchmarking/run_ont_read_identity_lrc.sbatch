#!/bin/bash
#SBATCH --job-name=ont_read_identity
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=ont_read_identity_lrc_%j.out
#SBATCH --error=ont_read_identity_lrc_%j.err

# ONT read-identity measurement (td-4e19d.29) — Lawrencium.
#
# Regenerates the three TSVs under results/ont_read_identity/ with the script's
# DEFAULT parameters (Lambda NC_001416, simulated ONT reads at 30x, seed 42),
# which are the parameters of the committed run. Do not pass --coverage/--seed
# here unless the intent is a DIFFERENT measurement.
#
# Output goes to a fresh per-job directory, never to the repo's results/
# directory: the script writes fixed filenames, so pointing two runs at one
# directory would let the later silently replace the earlier. Copy the three
# TSVs into the repo deliberately once they have been checked.
#
# Submit from the repo root (so --project=. resolves):
#   sbatch --chdir="$WT" benchmarking/run_ont_read_identity_lrc.sbatch
# Override the output root:
#   ONT_IDENTITY_OUT_ROOT=/some/dir sbatch --chdir="$WT" benchmarking/run_ont_read_identity_lrc.sbatch
# Prereq: deps instantiated once on a login node (JULIA_PKG_PRECOMPILE_AUTO=0).
# Self-test (stubbed julia, no cluster needed):
#   /bin/bash benchmarking/run_ont_read_identity_lrc.selftest.sh
#
# A green exit means: julia exited 0, all three TSVs exist and are non-empty,
# the summary records a non-empty Badread version and no "unknown" field, and
# ${OUT_DIR}/PROVENANCE records the git commit, tree state and Julia version
# with no "unknown" value. Julia exiting 0 alone is not treated as success.
#
# `-e` is deliberately omitted: every step that matters is guarded explicitly.
# Setup steps before PROVENANCE exists exit 1 directly; after it exists, every
# exit goes through `finish` so the file always ends with exit_status. A new
# step added below must be guarded the same way.
set -uo pipefail
echo "=== ONT read-identity measurement (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Every provenance value is either measured or the literal "unknown", which
# fails the job at the end. A failed git call must never read as a clean tree.
if ! GIT_HEAD="$(git rev-parse HEAD 2>/dev/null)" || [[ -z "${GIT_HEAD}" ]]; then
  GIT_HEAD=unknown
fi
if GIT_STATUS="$(git status --porcelain --untracked-files=no 2>/dev/null)"; then
  if [[ -z "${GIT_STATUS}" ]]; then
    GIT_DIRTY=0
  else
    GIT_DIRTY="$(printf '%s\n' "${GIT_STATUS}" | wc -l | tr -d ' ')"
  fi
else
  GIT_DIRTY=unknown
fi
echo "Repo: $(pwd)  HEAD: ${GIT_HEAD}  tracked files modified: ${GIT_DIRTY}"

# juliaup lts CHANNEL: this pins the channel, not a fixed Julia version, so the
# resolved version is echoed below and written to PROVENANCE. See
# run_track_a_baseline_lrc.sbatch for why the julia module and a bare PATH
# lookup are both avoided (td-j8bh).
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable (juliaup lts channel missing/broken on this host)." >&2
  exit 1
fi
# stdout only, first line only: juliaup can print notices on stderr, and a
# second line would break PROVENANCE's one-key-per-line format.
JULIA_VERSION="$("${JULIA}" "${JULIA_ARGS[@]}" --version 2>/dev/null | head -n 1)"
[[ -n "${JULIA_VERSION}" ]] || JULIA_VERSION=unknown
echo "JULIA: ${JULIA} ${JULIA_ARGS[*]}  (${JULIA_VERSION})"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

if [[ ! -f benchmarking/ont_read_identity.jl ]]; then
  echo "ERROR: run from the repo root (sbatch --chdir=<repo>); not found: $(pwd)/benchmarking/ont_read_identity.jl" >&2
  exit 1
fi

# Manifest.toml is not tracked, so git_head alone does not pin the resolved
# dependency versions. "absent" is a measured state, not an unknown one.
if [[ -f Manifest.toml ]]; then
  MANIFEST_SHA="$( (sha256sum Manifest.toml 2>/dev/null || shasum -a 256 Manifest.toml 2>/dev/null) | cut -d ' ' -f 1)"
  [[ -n "${MANIFEST_SHA}" ]] || MANIFEST_SHA=unknown
else
  MANIFEST_SHA=absent
fi

OUT_ROOT="${ONT_IDENTITY_OUT_ROOT:-/global/scratch/users/${USER}/ont_read_identity}"
OUT_DIR="${OUT_ROOT}/${SLURM_JOB_ID}"
PROVENANCE="${OUT_DIR}/PROVENANCE"
# The script reuses any non-empty reads file and SAM it finds in OUT_DIR, so a
# directory left by an interrupted attempt would be scored as if it were
# complete. Refuse a non-empty OUT_DIR rather than resume into it.
if [[ -d "${OUT_DIR}" ]] && [[ -n "$(ls -A "${OUT_DIR}" 2>/dev/null)" ]]; then
  echo "ERROR: ${OUT_DIR} already exists and is not empty; refusing to reuse a partial run" >&2
  exit 1
fi
mkdir -p "${OUT_DIR}" || { echo "ERROR: cannot create ${OUT_DIR}" >&2; exit 1; }
echo "Output: ${OUT_DIR}"
{
  echo "slurm_job_id=${SLURM_JOB_ID}"
  echo "node=${SLURM_NODELIST}"
  echo "repo=$(pwd)"
  echo "git_head=${GIT_HEAD}"
  echo "git_tracked_files_modified=${GIT_DIRTY}"
  echo "manifest_sha256=${MANIFEST_SHA}"
  echo "julia=${JULIA_VERSION}"
  echo "started=$(date)"
} > "${PROVENANCE}" || { echo "ERROR: cannot write ${PROVENANCE}" >&2; exit 1; }

# From here on PROVENANCE exists; every exit records its status there. If that
# final append itself fails, the job fails rather than leaving a file that
# reads as still running.
finish() {
  local code="$1"
  if ! echo "exit_status=${code}" >> "${PROVENANCE}"; then
    echo "ERROR: cannot append exit_status to ${PROVENANCE}" >&2
    code=1
  fi
  exit "${code}"
}

# The script records `badread --version` BEFORE it simulates, and simulation is
# what installs the env. On a cold depot that probe would record "unknown" in
# read_identity_summary.tsv, so the envs are provisioned first.
echo ">>> PREFLIGHT conda envs  $(date)"
if ! "${JULIA}" "${JULIA_ARGS[@]}" --project=. -e \
  'import Mycelia; Mycelia.add_bioconda_env("badread"); Mycelia.add_bioconda_env("minimap2")'; then
  echo ">>> FAIL preflight (conda env provisioning)" >&2
  finish 1
fi

echo ">>> RUN ont_read_identity.jl  $(date)"
if "${JULIA}" "${JULIA_ARGS[@]}" --project=. benchmarking/ont_read_identity.jl --output-dir "${OUT_DIR}"; then
  echo ">>> OK  ont_read_identity.jl"
  status=0
else
  echo ">>> FAIL ont_read_identity.jl"
  status=1
fi

# julia exiting 0 does not prove the outputs exist: kmer_survival_ladder.tsv is
# written only when at least one read aligned, and the Badread version/default
# probes fall back to "unknown" (or an empty string) with at most a warning.
if [[ "${status}" -eq 0 ]]; then
  for f in per_read_identity.tsv read_identity_summary.tsv kmer_survival_ladder.tsv; do
    if [[ ! -s "${OUT_DIR}/${f}" ]]; then
      echo ">>> FAIL missing or empty: ${OUT_DIR}/${f}" >&2
      status=1
    fi
  done
  SUMMARY="${OUT_DIR}/read_identity_summary.tsv"
  if [[ -s "${SUMMARY}" ]]; then
    if grep -q -w unknown "${SUMMARY}"; then
      echo ">>> FAIL read_identity_summary.tsv records 'unknown' provenance" >&2
      status=1
    fi
    # Every data row must carry a non-empty badread_version.
    if ! awk -F '\t' '
      NR == 1 { for (i = 1; i <= NF; i++) if ($i == "badread_version") col = i; next }
      col == 0 || $col == "" { bad = 1 }
      END { exit (col == 0 || bad || NR < 2) }' "${SUMMARY}"; then
      echo ">>> FAIL read_identity_summary.tsv has a missing or empty badread_version" >&2
      status=1
    fi
  fi
fi
if grep -q '=unknown$' "${PROVENANCE}"; then
  echo ">>> FAIL PROVENANCE records an unknown value:" >&2
  grep '=unknown$' "${PROVENANCE}" >&2
  status=1
fi

echo ""
echo "=== complete: $(date) ==="
echo "--- result artifacts ---"
ls -l "${OUT_DIR}"/*.tsv || true

finish "${status}"
