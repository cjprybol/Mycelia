#!/bin/bash
#SBATCH --job-name=mycelia_rgv_sweep
#SBATCH --account=m4269
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=rgv_sweep_nersc_%j.out
#SBATCH --error=rgv_sweep_nersc_%j.err

# Rhizomorph correction-core real-genome validation sweep (td-wlfq / td-4e19d) —
# NERSC/Perlmutter variant of run_rgv_sweep_lrc.sbatch. Same sweep, NERSC-site
# directives (account m4269, --constraint=cpu, --qos=regular, julia/1.10.10,
# ~/.julia depot). See the Lawrencium version for full rationale.
#
# Run from benchmarking/:  sbatch run_rgv_sweep_nersc.sbatch
# Prereq (ONCE on a login node, same checkout):
#   module load julia/1.10.10
#   benchmarking/hpc-setup.sh
# NOTE: keep `set -uo pipefail` WITHOUT `-e`. The green-but-empty guard below
# relies on awk/grep non-zero exit codes NOT aborting the script; adding `set -e`
# would abort before the guard can set rc and echo its failure reason.
set -uo pipefail
echo "=== Mycelia Rhizomorph correction-validation sweep (NERSC) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

# Do NOT `module purge`: it strips Cray PE defaults and breaks NERSC's juliaup
# modulefile. Restore PE defaults defensively, then load julia.
source /opt/cray/pe/cpe/25.09/restore_lmod_system_defaults.sh 2>/dev/null || true
module load julia/1.10.10

JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]]; then
  JULIA=/global/common/software/nersc9/julia/1.10.10/bin/julia
fi
if ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable at '${JULIA}' — module load failed or path stale." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${HOME}/.julia}"
export LD_LIBRARY_PATH=""

# === Sweep configuration (real Lambda-phage validation scale) =============
export MYCELIA_RGV_SMOKE="${MYCELIA_RGV_SMOKE:-false}"
export MYCELIA_RGV_ACCESSION="${MYCELIA_RGV_ACCESSION:-NC_001416}"
export MYCELIA_RGV_ERR="${MYCELIA_RGV_ERR:-0.01,0.05,0.10}"
export MYCELIA_RGV_READLEN="${MYCELIA_RGV_READLEN:-150,5000}"
export MYCELIA_RGV_COVERAGE="${MYCELIA_RGV_COVERAGE:-30}"
export MYCELIA_RGV_K="${MYCELIA_RGV_K:-21}"
export MYCELIA_RGV_SEED="${MYCELIA_RGV_SEED:-42}"
export MYCELIA_RUN_EXTERNAL="${MYCELIA_RUN_EXTERNAL:-true}"

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

echo ""
echo ">>> RUN rhizomorph_correction_validation_sweep.jl  $(date)"
echo "    ERR=${MYCELIA_RGV_ERR}  READLEN=${MYCELIA_RGV_READLEN}  COVERAGE=${MYCELIA_RGV_COVERAGE}x  K=${MYCELIA_RGV_K}  QUAST=${MYCELIA_RUN_EXTERNAL}"

# Capture julia stdout to a runlog while still streaming it to the job .out (tee),
# so the guard can bind to THIS run's CSV by the path the .jl prints
# ("CSV written: <path>") — not `ls -t` which can grab a concurrent/stale CSV.
runlog="rgv_sweep_${SLURM_JOB_ID}.runlog"
if LD_LIBRARY_PATH="" "${JULIA}" --project=.. rhizomorph_correction_validation_sweep.jl | tee "${runlog}"; then
  echo ">>> OK  sweep complete  $(date)"
  rc=0
else
  echo ">>> FAIL sweep  $(date)" >&2
  rc=1
fi

echo ""
# Post-check: a green exit is a silent failure UNLESS at least one
# (error_rate, readlen) cell has BOTH the naive AND the iterative arm passing —
# ok=true, mode=VERDICT, and n_contigs>0. The .jl records a failed/empty arm as
# ok=false (or n_contigs=0) and still writes the CSV + exits 0; and the sweep's
# whole purpose IS the naive-vs-iterative comparison, so one lone naive arm is
# NOT a result. Fail loud so sacct shows FAILED rather than green-but-empty.
if [[ "${rc}" -eq 0 ]]; then
  # Bind to THIS job's CSV via the path the .jl printed; fail loud if absent.
  csv_line="$(grep '^CSV written: ' "${runlog}" 2>/dev/null | tail -1)"
  csv="${csv_line#CSV written: }"
  if [[ -z "${csv}" ]]; then
    echo ">>> FAIL: sweep exited 0 but printed no 'CSV written:' path" >&2; rc=1
  elif [[ ! -f "${csv}" ]]; then
    echo ">>> FAIL: printed CSV path does not exist: ${csv}" >&2; rc=1
  # awk guard assumes comma-free CSV fields: CSV.jl quotes any embedded comma,
  # which would shift these -F, column positions. Safe for the current schema
  # (only numeric + bare-word columns: reference/regime/tech/arm/mode).
  elif ! awk -F, '
      NR==1 {
        for (i=1;i<=NF;i++) {
          if ($i=="arm")        a=i
          if ($i=="error_rate") er=i
          if ($i=="readlen")    rl=i
          if ($i=="ok")         ok=i
          if ($i=="mode")       m=i
          if ($i=="n_contigs")  nc=i
        }
        if (!a || !er || !rl || !ok || !m || !nc) exit 2
        next
      }
      ($ok=="true") && ($m=="VERDICT") && (($nc+0)>0) {
        key = $er SUBSEP $rl
        if ($a=="naive")          naive[key]=1
        else if ($a=="iterative") iter[key]=1
      }
      END { for (key in naive) if (key in iter) found=1; exit found ? 0 : 1 }
    ' "${csv}"; then
    echo ">>> FAIL: ${csv} has no cell with BOTH naive and iterative ok=true VERDICT arms (n_contigs>0) — green-but-empty" >&2; rc=1
  else
    echo ">>> OK: verified >=1 cell with BOTH naive+iterative ok=true VERDICT n_contigs>0 arms in ${csv}"
  fi
fi

echo "=== complete: $(date) ==="
if [[ "${rc}" -eq 0 ]]; then
  ls -t results/rhizomorph_correction_validation_sweep_*.csv 2>/dev/null | head || true
else
  echo "(sweep failed — not listing stale CSVs)" >&2
fi
exit "${rc}"
