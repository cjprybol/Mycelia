#!/bin/bash
#SBATCH --job-name=mycelia_rgv_sweep
#SBATCH --account=pc_mfnanofabio
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
# 2 days. The 4h that was here is not merely tight, it is known-insufficient: the
# 2026-07-25 submissions had to be overridden with --time=2-00:00:00 because 4h
# was too short for a SINGLE seed. MYCELIA_RGV_SEED now defaults to three seeds,
# so the sweep iterates the full (err x readlen x arm) matrix three times. A
# timeout kills the job mid-sweep and the green-but-empty guard below never runs,
# so the failure would present as a missing CSV rather than as a clear timeout.
# lr_normal MaxWall is 3 days, so 2 days leaves headroom without hitting the cap.
#SBATCH --time=2-00:00:00
#SBATCH --output=rgv_sweep_lrc_%j.out
#SBATCH --error=rgv_sweep_lrc_%j.err

# Rhizomorph correction-core real-genome validation sweep (td-wlfq / td-4e19d) —
# Lawrencium. Runs benchmarking/rhizomorph_correction_validation_sweep.jl on the
# REAL Lambda-phage reference (NC_001416, ~48.5 kb) across the
# (error_rate x read-regime) x {naive, iterative} matrix at validation-scale
# coverage, so the scale guard emits a VERDICT (not SMOKE-ONLY).
#
# This answers the manuscript's central preliminary question: does iterative
# correction improve final ASSEMBLY quality over naive across error regimes, and
# in which regime does correction deliver value?
#
# Run from benchmarking/:  sbatch run_rgv_sweep_lrc.sbatch
# Prereq (ONCE on a login node, in the same checkout):
#   module load julia/1.10.2-11.4
#   benchmarking/hpc-setup.sh          # resolve + instantiate + precompile
# NOTE: keep `set -uo pipefail` WITHOUT `-e`. The green-but-empty guard below
# relies on awk/grep non-zero exit codes NOT aborting the script; adding `set -e`
# would abort before the guard can set rc and echo its failure reason.
set -uo pipefail
echo "=== Mycelia Rhizomorph correction-validation sweep (Lawrencium) ==="
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURM_NODELIST}  CPUs: ${SLURM_CPUS_PER_TASK}  Start: $(date)"

module load julia/1.10.2-11.4

# Resolve the julia binary explicitly and fail LOUDLY if the module did not put
# julia on PATH — otherwise an empty ${JULIA} would run the sweep args as a
# command, fail, and still exit green with no results.
JULIA="$(command -v julia || true)"
if [[ -z "${JULIA}" ]] || ! "${JULIA}" --version >/dev/null 2>&1; then
  echo "ERROR: julia not runnable after 'module load julia/1.10.2-11.4'." >&2
  exit 1
fi
echo "JULIA: ${JULIA}  ($("${JULIA}" --version 2>&1 | head -1))"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Reuse the group depot where hpc-setup.sh precompiled the deps (SLURM does not
# reliably inherit interactive-shell exports). `:-` keeps any pre-exported value.
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot}"
# HPC Julia rule: clear LD_LIBRARY_PATH to avoid system libstdc++ conflicts.
export LD_LIBRARY_PATH=""

# === Sweep configuration (real Lambda-phage validation scale) =============
# Override any of these on the sbatch line with `--export=ALL,VAR=val,...`.
export MYCELIA_RGV_SMOKE="${MYCELIA_RGV_SMOKE:-false}"          # FULL run: download NC_001416
export MYCELIA_RGV_ACCESSION="${MYCELIA_RGV_ACCESSION:-NC_001416}"
export MYCELIA_RGV_ERR="${MYCELIA_RGV_ERR:-0.01,0.05,0.10}"     # low / mid / high error regimes
export MYCELIA_RGV_READLEN="${MYCELIA_RGV_READLEN:-150,5000}"   # short-Illumina / long-Nanopore
export MYCELIA_RGV_COVERAGE="${MYCELIA_RGV_COVERAGE:-30}"       # 48.5kb x 30x ~= 1.46 Mbase > scale floor
export MYCELIA_RGV_K="${MYCELIA_RGV_K:-21}"
# Three pre-registered replicate seeds. The sweep loops over this list, so one
# submission yields a table the paired analysis can pair over; a single seed
# would need three manual submissions plus a three-way --csv merge.
export MYCELIA_RGV_SEED="${MYCELIA_RGV_SEED:-42,123,456}"
# Enable QUAST (Mycelia auto-provisions it via Conda.jl); harness degrades
# gracefully to internal metrics if QUAST/Conda is unavailable.
export MYCELIA_RUN_EXTERNAL="${MYCELIA_RUN_EXTERNAL:-true}"

cd "${HOME}/workspace/Mycelia/benchmarking" || exit 1

echo ""
echo ">>> RUN rhizomorph_correction_validation_sweep.jl  $(date)"
echo "    ERR=${MYCELIA_RGV_ERR}  READLEN=${MYCELIA_RGV_READLEN}  COVERAGE=${MYCELIA_RGV_COVERAGE}x  K=${MYCELIA_RGV_K}  SEED=${MYCELIA_RGV_SEED}  QUAST=${MYCELIA_RUN_EXTERNAL}"

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
# (error_rate, readlen, seed) cell has BOTH the naive AND the iterative arm
# passing — ok=true, mode=VERDICT, and n_contigs>0. The .jl records a
# failed/empty arm as ok=false (or n_contigs=0) and still writes the CSV + exits
# 0; and the sweep's whole purpose IS the naive-vs-iterative comparison, so one
# lone naive arm is NOT a result. Fail loud so sacct shows FAILED rather than
# green-but-empty.
#
# The cell key includes `seed`. MYCELIA_RGV_SEED now defaults to THREE seeds, so
# an (error_rate, readlen)-only key can be satisfied by a naive arm from seed 42
# paired with an iterative arm from seed 123 — green, yet no single replicate
# ever produced a comparable pair. The seed column is OPTIONAL (its absence is
# NOT an exit-2 header error) so this guard still validates CSVs written before
# the seed column existed; those degrade correctly to the old single-seed key.
#
# The guard ALSO asserts, for every ok=true `iterative` row, that
# `corrector_sequencing_tech` (read back from the assembler's own stamp) equals
# the `tech` the reads were simulated with. That column exists to make a re-break
# of the corrector wiring visible, but a column nothing reads is provenance, not
# a detector — this post-check is its automated consumer. Both columns are
# OPTIONAL: `corrector_sequencing_tech` is absent from pre-wiring CSVs, and that
# absence is precisely what identifies them, so a missing column skips the check
# rather than failing it.
#
# The awk status is CAPTURED into `guard_rc`, never consumed by `! awk ...`:
# `!` is true for ANY non-zero status, so exit 1 (green-but-empty), exit 2 (bad
# header) and exit 3 (tech mismatch) all took the same branch and printed the
# same green-but-empty message. sacct showing FAILED is only half the point —
# the job .out has to say WHY, and the awk-level exit codes have no other
# consumer.
if [[ "${rc}" -eq 0 ]]; then
  # Bind to THIS job's CSV via the path the .jl printed; fail loud if absent.
  csv_line="$(grep '^CSV written: ' "${runlog}" 2>/dev/null | tail -1)"
  csv="${csv_line#CSV written: }"
  if [[ -z "${csv}" ]]; then
    echo ">>> FAIL: sweep exited 0 but printed no 'CSV written:' path" >&2; rc=1
  elif [[ ! -f "${csv}" ]]; then
    echo ">>> FAIL: printed CSV path does not exist: ${csv}" >&2; rc=1
  else
    # awk guard assumes comma-free CSV fields: CSV.jl quotes any embedded comma,
    # which would shift these -F, column positions. Safe for the current schema
    # (only numeric + bare-word columns: reference/regime/tech/arm/mode).
    guard_rc=0
    awk -F, '
      NR==1 {
        for (i=1;i<=NF;i++) {
          if ($i=="arm")        a=i
          if ($i=="error_rate") er=i
          if ($i=="readlen")    rl=i
          if ($i=="ok")         ok=i
          if ($i=="mode")       m=i
          if ($i=="n_contigs")  nc=i
          if ($i=="seed")       sd=i
          if ($i=="tech")       t=i
          if ($i=="corrector_sequencing_tech") ct=i
        }
        # `hdr_bad` is re-checked in END: awk runs END even after `exit`, and an
        # END-block `exit` OVERRIDES an earlier one — so a bare `exit 2` here was
        # unreachable and a malformed header got reported with the exit-1
        # "green-but-empty" reason. The caller branches on the captured status,
        # so the distinction now reaches the operator instead of dying at the
        # awk layer.
        if (!a || !er || !rl || !ok || !m || !nc) { hdr_bad=1; exit 2 }
        next
      }
      ($ok=="true") && ($m=="VERDICT") && (($nc+0)>0) {
        key = $er SUBSEP $rl SUBSEP (sd ? $sd : "-")
        if ($a=="naive")          naive[key]=1
        else if ($a=="iterative") iter[key]=1
      }
      t && ct && ($ok=="true") && ($a=="iterative") {
        checked++
        if ($ct != $t) {
          tech_bad=1
          if (!tech_shown) {
            print "GUARD-DETAIL: line " NR ": iterative row simulated tech=" $t " but corrector ran corrector_sequencing_tech=" $ct
            tech_shown=1
          }
        }
      }
      END {
        if (hdr_bad)  exit 2
        # Say out loud whether the corrector-profile check actually RAN. A
        # summary that silently covered zero rows is the same overclaim this
        # post-check exists to prevent, one level up.
        if (!ct)      print "GUARD-NOTE: no corrector_sequencing_tech column (pre-wiring CSV) — corrector-profile consistency NOT checked"
        else if (!t)  print "GUARD-NOTE: no tech column — corrector-profile consistency NOT checked"
        else          print "GUARD-NOTE: corrector-profile consistency checked on " checked+0 " ok=true iterative row(s)"
        if (tech_bad) exit 3
        for (key in naive) if (key in iter) found=1
        exit found ? 0 : 1
      }
    ' "${csv}" || guard_rc=$?
    if [[ "${guard_rc}" -eq 2 ]]; then
      echo ">>> FAIL: ${csv} is missing a required column (need arm/error_rate/readlen/ok/mode/n_contigs) — SCHEMA error, the table was never checked; this is NOT green-but-empty" >&2; rc=1
    elif [[ "${guard_rc}" -eq 3 ]]; then
      echo ">>> FAIL: ${csv} has an ok=true iterative row whose corrector_sequencing_tech disagrees with the simulated tech — the CORRECTOR ran a different error profile than the reads (see GUARD-DETAIL above); results are not comparable" >&2; rc=1
    elif [[ "${guard_rc}" -ne 0 ]]; then
      echo ">>> FAIL: ${csv} has no (error_rate, readlen, seed) cell with BOTH naive and iterative ok=true VERDICT arms (n_contigs>0) — green-but-empty" >&2; rc=1
    else
      echo ">>> OK: verified >=1 (error_rate, readlen, seed) cell with BOTH naive+iterative ok=true VERDICT n_contigs>0 arms in ${csv} — see the GUARD-NOTE line above for how many iterative rows the corrector-profile check actually covered"
    fi
  fi
fi

echo "=== complete: $(date) ==="
if [[ "${rc}" -eq 0 ]]; then
  echo "--- result CSVs ---"
  ls -t results/rhizomorph_correction_validation_sweep_*.csv 2>/dev/null | head || true
else
  echo "(sweep failed — not listing stale CSVs)" >&2
fi

exit "${rc}"
