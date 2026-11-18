# Node-local temp-root cleanup for run_d5_grid_nersc.sbatch (td-yddpr.5).
# Sourced, not executed. GNU coreutils (long flags). Tested by
# test/4_assembly/d5_ont_adaptive_grid_run_test.jl on Linux.

# d5_log_memory_peak
# Log memory at the cgroup level that enforces the job's --mem. Under SLURM's
# cgroup v2 plugin this shell sits in a leaf (job_<id>/step_batch/user/task_0)
# whose own memory.max is normally "max"; the limit is set further up, on
# step_batch or job_<id>. So walk up from the leaf to the first level whose
# memory.max is a number, but never above the job_<id> level: higher levels
# (slurmstepd.scope, system.slice) hold every job on the node, and on a shared
# node their peak is not this job's. Log that level's memory.max, memory.peak
# and oom_kill count (memory.events counts the level and everything below it).
# memory.peak includes tmpfs (the node-local temp root) and page cache, which
# peak_rss_bytes does not, so it is an upper bound to size --mem from.
# D5_CGROUP_ROOT and D5_PROC_CGROUP exist for tests.
d5_log_memory_peak() {
  local base="${D5_CGROUP_ROOT:-/sys/fs/cgroup}" rel dir max peak oom
  local by_hand="sacct MaxRSS excludes tmpfs, so add the temp root's peak size by hand before sizing --mem"
  rel="$(cut --delimiter=: --fields=3 "${D5_PROC_CGROUP:-/proc/self/cgroup}" 2>/dev/null | head --lines=1)"
  dir="${base}${rel}"
  dir="${dir%/}"
  if [[ -z "${rel}" || ! -d "${dir}" ]]; then
    echo "cgroup memory.peak unavailable (no cgroup v2 directory for this shell under ${base}); ${by_hand}" >&2
    return 0
  fi
  if [[ "${rel}" != */job_* ]]; then
    echo "cgroup memory.peak unavailable (${rel} is not inside a SLURM job_<id> cgroup); ${by_hand}" >&2
    return 0
  fi
  while :; do
    max=""
    if [[ -e "${dir}/memory.max" ]]; then
      if ! max="$(cat -- "${dir}/memory.max" 2>&1)" || [[ -z "${max}" ]]; then
        echo "cgroup memory.max at ${dir} unreadable (${max:-empty}); not reading higher levels; ${by_hand}" >&2
        return 0
      fi
      [[ "${max}" != "max" ]] && break
      max=""
    fi
    [[ "${dir##*/}" == job_* || "${dir}" != "${base}"/* ]] && break
    dir="$(dirname -- "${dir}")"
  done
  if [[ -z "${max}" ]]; then
    peak="$(cat -- "${base}${rel}/memory.peak" 2>/dev/null)"
    peak="${peak:+${peak} bytes}"
    echo "cgroup: no numeric memory.max from ${base}${rel} up to its job level; leaf memory.peak: ${peak:-unreadable}; ${by_hand}" >&2
    return 0
  fi
  oom="$(grep --max-count=1 '^oom_kill ' "${dir}/memory.events" 2>/dev/null | cut --delimiter=' ' --fields=2)"
  if peak="$(cat -- "${dir}/memory.peak" 2>/dev/null)" && [[ -n "${peak}" ]]; then
    echo "cgroup memory at the limiting level ${dir}: memory.max ${max} bytes; memory.peak (RSS + tmpfs + page cache) ${peak} bytes; oom_kill ${oom:-unreadable}" >&2
  else
    echo "cgroup memory at the limiting level ${dir}: memory.max ${max} bytes; memory.peak unreadable (memory.current now $(cat -- "${dir}/memory.current" 2>/dev/null || echo unreadable), not a peak); oom_kill ${oom:-unreadable}; ${by_hand}" >&2
  fi
}

# d5_temp_root_cleanup RC TEMP_ROOT EVIDENCE_DIR
# Log the cgroup memory peak. When RC is non-zero and the temp root is not
# empty, copy it to EVIDENCE_DIR: the corrector leaves fail-closed evidence
# under TMPDIR, cell error rows cite those paths, and the node is unreachable
# after the job. The copy goes to EVIDENCE_DIR.partial and is renamed only when
# complete, so a directory without .partial is whole. Then log the temp root's
# size and remove it.
# Further TERM/INT/HUP are ignored here: a second scancel during the copy would
# otherwise exit mid-cleanup, leaving a whole copy named .partial and the temp
# root on tmpfs. SLURM's SIGKILL after KillWait (30 s on Perlmutter) still
# bounds the cleanup; the "copying" line marks where a killed copy stopped.
d5_temp_root_cleanup() {
  trap '' TERM INT HUP
  local rc="$1" root="$2" keep="$3" entries resolved size total
  d5_log_memory_peak
  if [[ "${rc}" -ne 0 ]]; then
    if ! entries="$(ls --almost-all -- "${root}" 2>&1)"; then
      echo "WARN: temp root ${root} unreadable, failure evidence not checked: ${entries}" >&2
    elif [[ -n "${entries}" ]]; then
      resolved="$(realpath -- "${root}" 2>&1)" || resolved="unresolved: ${resolved}"
      if [[ -e "${keep}" || -e "${keep}.partial" ]]; then
        echo "WARN: ${keep} already exists; evidence in ${root} not kept" >&2
      elif ! mkdir --parents -- "${keep}.partial"; then
        echo "WARN: cannot create ${keep}.partial; evidence in ${root} not kept" >&2
      else
        echo "copying ${root} -> ${keep}.partial (rc=${rc}, started $(date --iso-8601=seconds))" >&2
        if cp --archive -- "${root}/." "${keep}.partial/" &&
          printf '%s\n' "${root}" >"${keep}.partial/ORIGINAL_TEMP_ROOT" &&
          printf '%s\n' "${resolved}" >"${keep}.partial/RESOLVED_TEMP_ROOT" &&
          mv --no-target-directory -- "${keep}.partial" "${keep}"; then
          echo "temp root preserved at ${keep} (rc=${rc})" >&2
        else
          echo "WARN: could not preserve ${root}; the partial copy is at ${keep}.partial" >&2
        fi
      fi
    fi
  fi
  if size="$(du --summarize --bytes -- "${root}" 2>&1)"; then
    echo "temp root ${root} at exit (rc=${rc}): $(printf '%s\n' "${size}" | cut --fields=1) bytes (end state, not the peak)" >&2
  else
    # A partial failure (an unreadable subdirectory) still prints a total that
    # undercounts; a total failure prints only errors.
    total="$(printf '%s\n' "${size}" | grep --extended-regexp $'^[0-9]+\t' | tail --lines=1 | cut --fields=1)"
    if [[ -n "${total}" ]]; then
      echo "temp root ${root} at exit (rc=${rc}): at least ${total} bytes (du incomplete: $(printf '%s' "${size}" | grep --invert-match --extended-regexp $'^[0-9]+\t' | tr '\n' ' '))" >&2
    else
      echo "temp root ${root} at exit (rc=${rc}): size unknown, du failed: $(printf '%s' "${size}" | tr '\n' ' ')" >&2
    fi
  fi
  rm --recursive --force -- "${root}" ||
    echo "WARN: temp root ${root} not fully removed; what is left stays on node-local tmpfs" >&2
}

# d5_install_temp_root_traps TEMP_ROOT EVIDENCE_DIR
# Clean up on every exit. A signal (SLURM time limit, scancel, preemption) is
# turned into a non-zero exit first: without that, the EXIT trap would see the
# status of the last completed command (0) and delete the evidence uncopied.
# Bash runs the TERM trap only once the foreground child (julia) exits, and the
# copy then has what is left of SLURM's KillWait before SIGKILL; a larger temp
# root may be cut off, which the "copying" line without a "preserved" line shows.
# Returns non-zero if any trap fails to install.
d5_install_temp_root_traps() {
  local root keep
  root="$(printf '%q' "$1")"
  keep="$(printf '%q' "$2")"
  # shellcheck disable=SC2064  # expand the paths now; $? when the trap fires
  trap "d5_temp_root_cleanup \"\$?\" ${root} ${keep}" EXIT &&
    trap 'exit 143' TERM &&
    trap 'exit 130' INT &&
    trap 'exit 129' HUP
}
