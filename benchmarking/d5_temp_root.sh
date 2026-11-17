# Node-local temp-root cleanup for run_d5_grid_nersc.sbatch (td-yddpr.5).
# Sourced, not executed. GNU coreutils (long flags). Tested by
# test/4_assembly/d5_ont_adaptive_grid_run_test.jl on Linux.

# d5_log_memory_peak
# Log the job cgroup's peak memory. On cgroup v2 memory.peak includes tmpfs
# (the node-local temp root) and page cache, which peak_rss_bytes does not, so
# it is the upper bound to size --mem from.
d5_log_memory_peak() {
  local rel cg
  rel="$(cut --delimiter=: --fields=3 /proc/self/cgroup 2>/dev/null | head --lines=1)"
  cg="/sys/fs/cgroup${rel}"
  if [[ -n "${rel}" && -r "${cg}/memory.peak" ]]; then
    echo "cgroup memory.peak (RSS + tmpfs + page cache): $(cat "${cg}/memory.peak") bytes (${cg})" >&2
  else
    echo "cgroup memory.peak unavailable (looked in ${cg}); size --mem from sacct and the temp root by hand" >&2
  fi
}

# d5_temp_root_cleanup RC TEMP_ROOT EVIDENCE_DIR
# Log the temp root's size and the cgroup memory peak. When RC is non-zero and
# the temp root is not empty, copy it to EVIDENCE_DIR: the corrector leaves
# fail-closed evidence under TMPDIR, cell error rows cite those paths, and the
# node is unreachable after the job. The copy goes to EVIDENCE_DIR.partial and
# is renamed only when complete, so a directory without .partial is whole.
# Always removes TEMP_ROOT.
d5_temp_root_cleanup() {
  local rc="$1" root="$2" keep="$3" entries
  echo "temp root ${root} at exit (rc=${rc}): $(du --summarize --bytes -- "${root}" 2>&1 | cut --fields=1) bytes (end state, not the peak)" >&2
  d5_log_memory_peak
  if [[ "${rc}" -ne 0 ]]; then
    if ! entries="$(ls --almost-all -- "${root}" 2>&1)"; then
      echo "WARN: temp root ${root} unreadable, failure evidence not checked: ${entries}" >&2
    elif [[ -n "${entries}" ]]; then
      if [[ -e "${keep}" || -e "${keep}.partial" ]]; then
        echo "WARN: ${keep} already exists; evidence in ${root} not kept" >&2
      elif mkdir --parents -- "${keep}.partial" &&
        cp --archive -- "${root}/." "${keep}.partial/" &&
        printf '%s\n' "${root}" >"${keep}.partial/ORIGINAL_TEMP_ROOT" &&
        mv --no-target-directory -- "${keep}.partial" "${keep}"; then
        echo "temp root preserved at ${keep} (rc=${rc})" >&2
      else
        echo "WARN: could not preserve ${root}; any partial copy is at ${keep}.partial" >&2
      fi
    fi
  fi
  rm --recursive --force -- "${root}"
}

# d5_install_temp_root_traps TEMP_ROOT EVIDENCE_DIR
# Clean up on every exit. A signal (SLURM time limit, scancel, preemption) is
# turned into a non-zero exit first: without that, the EXIT trap would see the
# status of the last completed command (0) and delete the evidence uncopied.
# The copy then runs inside SLURM's KillWait window before SIGKILL.
d5_install_temp_root_traps() {
  local root keep
  root="$(printf '%q' "$1")"
  keep="$(printf '%q' "$2")"
  # shellcheck disable=SC2064  # expand the paths now; $? when the trap fires
  trap "d5_temp_root_cleanup \"\$?\" ${root} ${keep}" EXIT
  trap 'exit 143' TERM
  trap 'exit 130' INT
  trap 'exit 129' HUP
}
