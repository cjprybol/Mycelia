#!/bin/bash
#SBATCH --job-name=rdv_ecoli_peak_rss
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr8
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
# lr_normal MaxWall is 3 days. Earlier E. coli runs did not finish a single
# native arm inside 110 min on lr_bigmem (td-ck03), so a short wall would only
# report a timeout.
#SBATCH --time=3-00:00:00
#SBATCH --no-requeue
#SBATCH --output=rdv_ecoli_peak_rss_%j.out
#SBATCH --error=rdv_ecoli_peak_rss_%j.err

# td-ck03 per-arm peak-RSS measurement: E. coli K-12 MG1655 (NC_000913.3) at 50x
# through benchmarking/real_data_corrector_validation.jl, on an lr8 node
# (128 cores, ~773 GB) so the measurement is not truncated by the 95 GB lr6
# ceiling that OOM-killed job 23984771.
#
# Four arms, each in its OWN julia process so each gets its own peak:
#   sota_spades  SPAdes on the raw reads (external baseline)
#   hybrid_olc   corrector=:iterative, strategy=:scalable, layout=:olc (MEGAHIT)
#   scalable     corrector=:iterative, strategy=:scalable
#   naive        corrector=:none
# They run in that order, cheapest first, so a naive arm that exhausts the wall
# does not cost the other three their measurement. A failing arm does not stop
# the later ones; the job exits non-zero if any arm failed.
#
# Peak RSS per arm comes from /usr/bin/time -v ("Maximum resident set size") of
# that arm's julia process. getrusage reports the largest single process in the
# tree, NOT a sum: for sota_spades it is max(julia, spades), not julia + spades.
# sacct MaxRSS for the batch step (printed at the end) is the job-wide peak; it
# is polled, so it can miss a short spike.
#
# DRIVER REQUIREMENT. master's real_data_corrector_validation.jl has neither the
# ecoli_k12 target nor the sota_spades arm nor MYCELIA_RDV_ARMS arm selection.
# They exist on branch cp/stage2-sota-parity. The preflight below refuses a
# driver that lacks any of them rather than letting it skip the unknown target
# and fail after the fact. Point RDV_RUN_DIR at a checkout that has them.
#
# Usage:
#   RDV_RUN_DIR=<checkout with the 4-arm driver> \
#     sbatch benchmarking/run_rdv_ecoli_peak_rss_lrc.sbatch
# Overrides (sbatch exports the submitting environment):
#   RDV_RUN_DIR=<checkout>      checkout whose driver runs (default: the git
#                               toplevel of the submit directory)
#   RDV_ARMS=a,b                subset/order of arms (default: the four above);
#                               one arm per job lets arms run in parallel
#   RDV_MANIFEST=resolve        (default) Pkg.resolve() + instantiate against
#                               <checkout>/Project.toml inside the job
#   RDV_MANIFEST=existing       use <checkout>/Manifest.toml as-is (instantiate
#                               only) -- use this for concurrent jobs on one
#                               checkout after resolving it once
#   RDV_MANIFEST=<path>         copy that Manifest.toml in, then instantiate only
#   RDV_TMP_BASE=<dir>          node-local parent for TMPDIR (default ${SLURM_TMPDIR:-/tmp})
#   RDV_OUT_ROOT=<dir>          where per-arm logs + the peak-RSS TSV go
#                               (default /global/scratch/users/${USER}/rdv_ecoli_peak_rss)
#   MYCELIA_RDV_COVERAGE / _K / _SEED pass through (defaults 50 / 21 / 42).
#
# `-e` is deliberately omitted: every step that matters is guarded explicitly.
set -uo pipefail
echo "=== td-ck03 E. coli 50x per-arm peak RSS (Lawrencium lr8) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK:-unset}  Start: $(date)"

# Pinned juliaup lts launcher; see run_rgv_sweep_lrc.sbatch (td-j8bh) for why a
# PATH lookup or `module load julia` is not used.
JULIA="${HOME}/.juliaup/bin/julia"
JULIA_ARGS=(+lts)
if [[ ! -x "${JULIA}" ]] || ! "${JULIA}" "${JULIA_ARGS[@]}" --version >/dev/null 2>&1; then
  echo "ERROR: '${JULIA} ${JULIA_ARGS[*]}' not runnable (juliaup lts channel missing/broken on this host)." >&2
  exit 1
fi
echo "JULIA: ${JULIA} ${JULIA_ARGS[*]}  ($("${JULIA}" "${JULIA_ARGS[@]}" --version 2>&1 | head -1))"
TIME_BIN=/usr/bin/time
if [[ ! -x "${TIME_BIN}" ]]; then
  echo "ERROR: ${TIME_BIN} not found; it is the per-arm peak-RSS instrument." >&2
  exit 1
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""
# Thread count comes from the allocation; refuse to guess (see
# run_ecoli_qualmer_memory_gate_lrc.sbatch). Thread count may itself move peak
# memory (td-ck03 notes, lr6 vs lr8), so it is logged with every arm.
if [[ -z "${SLURM_CPUS_PER_TASK:-}" ]]; then
  echo "ERROR: SLURM_CPUS_PER_TASK is unset — the --cpus-per-task directive is missing." >&2
  exit 1
fi
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

# --- Run directory ------------------------------------------------------------
if [[ -n "${RDV_RUN_DIR:-}" ]]; then
  run_dir="${RDV_RUN_DIR}"
elif ! run_dir="$(git -C "${SLURM_SUBMIT_DIR}" rev-parse --show-toplevel)" || [[ -z "${run_dir}" ]]; then
  echo "ERROR: ${SLURM_SUBMIT_DIR} is not inside a git checkout; submit from the checkout or set RDV_RUN_DIR." >&2
  exit 1
fi
driver="${run_dir}/benchmarking/real_data_corrector_validation.jl"
if [[ ! -f "${driver}" ]] || [[ ! -f "${run_dir}/Project.toml" ]]; then
  echo "ERROR: ${run_dir} is not a Mycelia checkout with benchmarking/real_data_corrector_validation.jl." >&2
  exit 1
fi
# Capability preflight (see DRIVER REQUIREMENT above). A text check, so it can
# only refuse a driver that obviously lacks a feature; it cannot prove one works.
missing=()
grep -q '"ecoli_k12"' "${driver}" || missing+=("ecoli_k12 target")
grep -q 'MYCELIA_RDV_ARMS' "${driver}" || missing+=("MYCELIA_RDV_ARMS arm selection")
grep -q ':sota_spades' "${driver}" || missing+=("sota_spades arm")
if [[ "${#missing[@]}" -gt 0 ]]; then
  echo "ERROR: ${driver} lacks: ${missing[*]}." >&2
  echo "       Point RDV_RUN_DIR at a checkout whose driver has them (branch cp/stage2-sota-parity)." >&2
  exit 1
fi
cd "${run_dir}" || exit 1
echo "Run dir: ${run_dir}  HEAD: $(git rev-parse HEAD 2>&1)"
echo "Uncommitted tracked changes (empty if clean):"
git status --short --untracked-files=no 2>&1 || true

arms_csv="${RDV_ARMS:-sota_spades,hybrid_olc,scalable,naive}"
IFS=',' read -r -a arms <<< "${arms_csv}"
if [[ "${#arms[@]}" -eq 0 ]]; then
  echo "ERROR: RDV_ARMS is empty." >&2
  exit 1
fi
for arm in "${arms[@]}"; do
  case "${arm}" in
    naive|scalable|hybrid_olc|sota_spades) ;;
    *) echo "ERROR: unknown arm '${arm}' in RDV_ARMS (valid: naive, scalable, hybrid_olc, sota_spades)." >&2; exit 1 ;;
  esac
done

export MYCELIA_RDV_TARGETS=ecoli_k12
export MYCELIA_RDV_SMOKE=false
export MYCELIA_RDV_COVERAGE="${MYCELIA_RDV_COVERAGE:-50}"
export MYCELIA_RDV_K="${MYCELIA_RDV_K:-21}"
export MYCELIA_RDV_SEED="${MYCELIA_RDV_SEED:-42}"

out_dir="${RDV_OUT_ROOT:-/global/scratch/users/${USER}/rdv_ecoli_peak_rss}/${SLURM_JOB_ID}"
if [[ -d "${out_dir}" ]] && [[ -n "$(ls -A "${out_dir}" 2>/dev/null)" ]]; then
  echo "ERROR: ${out_dir} already exists and is not empty; refusing to mix runs." >&2
  exit 1
fi
mkdir -p "${out_dir}" || { echo "ERROR: cannot create ${out_dir}" >&2; exit 1; }
echo "Output: ${out_dir}"

# --- Node-local TMPDIR (same construction as run_td2rxh_evidence_lrc.sbatch) --
# The driver's mktempdir() (reads, per-arm assemblies, SPAdes scratch) and the
# iterative corrector follow TMPDIR. Lawrencium's default TMPDIR is on Lustre,
# where the corrector's renameat2(RENAME_NOREPLACE) cleanup returns EINVAL
# (td-2o7h), and a symlinked temp root is rejected. Always removed on exit; a
# SIGKILL after KillWait skips the trap and leaves it on the node.
tmp_base="${RDV_TMP_BASE:-${SLURM_TMPDIR:-/tmp}}"
if ! tmp_root="$(mktemp -d "${tmp_base}/rdv_ecoli_${SLURM_JOB_ID}.XXXXXX")"; then
  echo "ERROR: cannot create a TMPDIR under ${tmp_base}." >&2
  exit 1
fi
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
case "${tmp_fs}" in
  ext2/ext3|xfs|tmpfs|btrfs|overlayfs) ;;
  *)
    echo "ERROR: TMPDIR ${TMPDIR} is on ${tmp_fs}; Mycelia needs a local filesystem (ext/xfs/tmpfs/btrfs/overlayfs)." >&2
    exit 1 ;;
esac
echo "TMPDIR=${TMPDIR} (${tmp_fs}, $(df -h "${TMPDIR}" | awk 'NR==2 {print $4 " free"}'))"
if ! julia_tmp="$("${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no -e 'print(tempdir())')" ||
   [[ "${julia_tmp%/}" != "${TMPDIR%/}" ]]; then
  echo "ERROR: Julia tempdir() is '${julia_tmp:-}', expected ${TMPDIR}." >&2
  exit 1
fi
echo "Julia tempdir(): ${julia_tmp}"

# --- Manifest -------------------------------------------------------------------
manifest_mode="${RDV_MANIFEST:-resolve}"
case "${manifest_mode}" in
  resolve)
    pkg_expr='import Pkg; Pkg.resolve(); Pkg.instantiate()' ;;
  existing)
    if [[ ! -f Manifest.toml ]]; then
      echo "ERROR: RDV_MANIFEST=existing but ${run_dir}/Manifest.toml does not exist." >&2
      exit 1
    fi
    pkg_expr='import Pkg; Pkg.instantiate()' ;;
  *)
    if [[ ! -f "${manifest_mode}" ]]; then
      echo "ERROR: RDV_MANIFEST='${manifest_mode}' is not resolve, existing, or an existing file." >&2
      exit 1
    fi
    if ! cp -f "${manifest_mode}" Manifest.toml; then
      echo "ERROR: cannot copy ${manifest_mode} to ${run_dir}/Manifest.toml." >&2
      exit 1
    fi
    pkg_expr='import Pkg; Pkg.instantiate()' ;;
esac
echo ">>> Pkg (${manifest_mode}) against ${run_dir}/Project.toml $(date)"
if ! JULIA_PKG_PRECOMPILE_AUTO=0 "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project=. -e "${pkg_expr}"; then
  echo ">>> Pkg ${manifest_mode} FAILED" >&2
  exit 1
fi
# Precompile outside the timed arms so the first arm's peak is not a compile peak.
echo ">>> precompile warmup (import Mycelia) $(date)"
if ! "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project=. -e 'import Mycelia; println("Mycelia loaded")'; then
  echo ">>> precompile/load FAILED" >&2
  exit 1
fi

# --- Arms -------------------------------------------------------------------------
summary="${out_dir}/peak_rss.tsv"
printf 'arm\texit\tpeak_rss_kb\tpeak_rss_gb\telapsed\tthreads\tcoverage\tseed\tnode\n' > "${summary}"
rc=0
for arm in "${arms[@]}"; do
  log="${out_dir}/${arm}.log"
  timefile="${out_dir}/${arm}.time"
  echo ">>> ARM ${arm} start $(date)"
  MYCELIA_RDV_ARMS="${arm}" "${TIME_BIN}" -v -o "${timefile}" \
    "${JULIA}" "${JULIA_ARGS[@]}" --startup-file=no --project=. \
    benchmarking/real_data_corrector_validation.jl > "${log}" 2>&1
  arm_rc=$?
  peak_kb="$(awk -F': ' '/Maximum resident set size/ {print $2}' "${timefile}" 2>/dev/null)"
  elapsed="$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "${timefile}" 2>/dev/null)"
  if [[ ! "${peak_kb}" =~ ^[0-9]+$ ]]; then
    echo ">>> ARM ${arm}: no peak RSS recorded in ${timefile}" >&2
    peak_kb="NA"; peak_gb="NA"; rc=1
  else
    peak_gb="$(awk -v kb="${peak_kb}" 'BEGIN { printf "%.1f", kb / 1048576 }')"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${arm}" "${arm_rc}" "${peak_kb}" "${peak_gb}" \
    "${elapsed:-NA}" "${JULIA_NUM_THREADS}" "${MYCELIA_RDV_COVERAGE}" "${MYCELIA_RDV_SEED}" \
    "${SLURM_NODELIST}" >> "${summary}"
  echo "PEAK_RSS arm=${arm} exit=${arm_rc} peak_rss_gb=${peak_gb} elapsed=${elapsed:-NA}"
  if [[ "${arm_rc}" -ne 0 ]]; then
    echo ">>> ARM ${arm} FAILED (exit ${arm_rc}); tail of ${log}:" >&2
    tail -n 30 "${log}" >&2
    rc=1
  fi
done

echo ""
echo "=== complete rc=${rc} $(date) ==="
cat "${summary}"
echo "Job-wide Slurm peak (after the job ends): sacct -j ${SLURM_JOB_ID} --format=JobID,State,Elapsed,MaxRSS"
exit "${rc}"
