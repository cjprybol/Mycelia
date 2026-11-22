#!/usr/bin/env bash
#
# check-sbatch-pattern.sh — lint SLURM sbatch scripts for footguns that static
# `bash -n` cannot catch.
#
# Ported from panlingual-metagraphs scripts/check-sbatch-pattern.sh (GitLab
# nanofabio/panlingual-metagraphs MR !19) for Mycelia bead td-hhkd. The checks
# are the classes that slipped through review on the durable-HPC-pattern work
# (td-2rfi / td-b05w):
#
#   1. `julia -u`                   — Julia has NO unbuffered flag; it errors at
#                                     startup. Use a `script -e` pseudo-TTY wrap.
#   2. JULIA_NUM_THREADS=auto       — oversubscribes a shared node and segfaults
#                                     (exit 139). Pin to SLURM_CPUS_PER_TASK /
#                                     SLURM_CPUS_ON_NODE.
#   3. $VAR / ${VAR} in a #SBATCH   — SLURM does NOT expand shell vars in #SBATCH
#      output/error path              directives, so logs scatter to a literal
#                                     `$HOME/` (or `$SCRATCH/`...) dir. Use a
#                                     relative path.
#   4. NERSC-ism in a Lawrencium    — a Lawrencium script that names Perlmutter
#      script                         hardware (milan) or the NERSC CFS depot path
#                                     is a copy-paste from the NERSC variant.
#
# A script is Lawrencium when its name ends in `_lrc.sbatch`, its path contains
# `/lawrencium/`, or it carries a `# Site profile: lawrencium` line (the line
# Mycelia's renderer emits).
#
# Full-line `#` comments are exempt. Corrected scripts legitimately mention
# these strings in explanatory comments (for example `# "julia -u" errors at
# startup`, or `# as on Perlmutter` as contrast documentation). `#SBATCH`
# directives are not comments: check 3 matches the directive value, and check 4
# still flags NERSC hardware or the CFS depot path on a directive. A Perlmutter
# or `/global/cfs` path in code is still a hit.
#
# This script only reads files. It does not call sbatch or any scheduler.
# It is intentionally not a GitHub Actions job: `.github/workflows/ci.yml` is
# owned by other pull requests. Invoke it from the command line or from
# `test/8_tool_integration/sbatch_pattern_lint.jl`.
#
# Usage:
#   scripts/check-sbatch-pattern.sh [PATH...]   lint *.sbatch (default: templates/slurm)
#   scripts/check-sbatch-pattern.sh --selftest  run built-in known-bad/good cases
#   scripts/check-sbatch-pattern.sh --help       show this header
#
# Exit: 0 clean, 1 violation found, 2 usage/environment error.

set -euo pipefail

# --- regex patterns (single source of truth; reused by lint + --selftest) ------
# `-u` as a standalone flag token anywhere after `julia ` on the same command,
# stopping at any `#` so a trailing comment is never matched.
RX_JULIA_U='julia[[:space:]]([^#]*[[:space:]])?-u([[:space:]]|$)'
RX_THREADS_AUTO='JULIA_NUM_THREADS=auto'
# Any unexpanded $VAR or ${VAR} in the #SBATCH output/error VALUE (before
# whitespace/inline-comment): $HOME, ${HOME}, $SCRATCH, ... are all the same bug.
RX_SBATCH_VAR='^#SBATCH --(output|error)=[^[:space:]#]*[$][{]?[A-Za-z_]'
RX_LRC_NERSC='[Mm]ilan|[Pp]erlmutter|/global/cfs'
RX_SITE_LAWRENCIUM='^# Site profile:[[:space:]]*lawrencium[[:space:]]*$'

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEFAULT_TARGET="$REPO_ROOT/templates/slurm"

usage() { sed -n '/^# Usage:/,/^# Exit:/p' "$0"; }

# Fail loud if the environment's grep cannot do what the checks require, rather
# than letting an unsupported flag error to /dev/null and read as "no hits"
# (a silent pass). Targets BusyBox + GNU grep alike (ERE, -n, -H).
preflight() {
    if ! printf 'aXb\n' | grep -HnE 'a[[:space:]]?X' >/dev/null 2>&1; then
        printf 'error: grep -HnE (POSIX ERE) is unavailable in this environment.\n' >&2
        exit 2
    fi
}

# drop `grep -Hn` hits whose content is a full-line comment (so a string
# mentioned in a `#` comment is not flagged as code).
strip_comment_hits() { grep -vE '^[^:]+:[0-9]+:[[:space:]]*#' || true; }

# Same as strip_comment_hits, but keep `#SBATCH` directives. POSIX ERE has no
# lookahead, so this matches a full-line `#` that is not the `#SBATCH` token.
strip_prose_comment_hits() {
    grep -vE '^[^:]+:[0-9]+:[[:space:]]*#($|[^S]|S($|[^B])|SB($|[^A])|SBA($|[^T])|SBAT($|[^C])|SBATC($|[^H])|SBATCH($|[^[:space:]]))' || true
}

is_lawrencium_script() {
    local f="$1"
    case "$f" in
        *_lrc.sbatch) return 0 ;;
    esac
    case "/$f/" in
        */lawrencium/*) return 0 ;;
    esac
    grep -qE "$RX_SITE_LAWRENCIUM" "$f"
}

# lint_paths <path...> -> 0 clean / 1 violations / 2 usage. Prints a report.
# Directories are searched recursively for *.sbatch. Uses full paths so
# `grep -Hn` always labels with the filename (no mapfile — bash 3.2 ok).
lint_paths() {
    local -a wrappers=() lrc=()
    local path f fail=0

    if [[ $# -eq 0 ]]; then
        printf 'error: no paths to lint\n' >&2
        return 2
    fi

    for path in "$@"; do
        if [[ -d "$path" ]]; then
            while IFS= read -r -d '' f; do
                wrappers+=("$f")
            done < <(find "$path" -type f -name '*.sbatch' -print0)
        elif [[ -f "$path" ]]; then
            case "$path" in
                *.sbatch) wrappers+=("$path") ;;
                *)
                    printf 'error: not an sbatch file: %s\n' "$path" >&2
                    return 2
                    ;;
            esac
        else
            printf 'error: not found: %s\n' "$path" >&2
            return 2
        fi
    done

    if [[ ${#wrappers[@]} -eq 0 ]]; then
        printf 'error: no *.sbatch files found\n' >&2
        return 2
    fi

    for f in "${wrappers[@]}"; do
        if is_lawrencium_script "$f"; then
            lrc+=("$f")
        fi
    done

    _report() {  # _report <label> <grep-output>
        if [[ -n "$2" ]]; then
            printf '\n[FAIL] %s\n' "$1"
            printf '%s\n' "$2" | sed 's/^/    /'
            fail=1
        fi
    }

    _report "invalid 'julia -u' flag — Julia has no unbuffered flag; use a 'script -e' PTY-wrap" \
        "$(grep -HnE -- "$RX_JULIA_U" "${wrappers[@]}" 2>/dev/null | strip_comment_hits)"
    _report "JULIA_NUM_THREADS=auto — oversubscribes shared nodes (exit-139 segfault); pin to SLURM_CPUS_*" \
        "$(grep -HnE -- "$RX_THREADS_AUTO" "${wrappers[@]}" 2>/dev/null | strip_comment_hits)"
    _report "unexpanded shell var (\$VAR/\${VAR}) in a #SBATCH path — SLURM does not expand it; use a relative path" \
        "$(grep -HnE -- "$RX_SBATCH_VAR" "${wrappers[@]}" 2>/dev/null)"
    if [[ ${#lrc[@]} -gt 0 ]]; then
        _report "NERSC-ism in a Lawrencium script — 'milan'/'Perlmutter'/CFS depot are NERSC-only" \
            "$(grep -HnE -- "$RX_LRC_NERSC" "${lrc[@]}" 2>/dev/null | strip_prose_comment_hits)"
    fi

    if [[ $fail -eq 0 ]]; then
        printf 'OK: %d sbatch scripts clean (no julia -u, no =auto, no unexpanded var in #SBATCH, no Lawrencium NERSC-isms).\n' \
            "${#wrappers[@]}"
    fi
    return "$fail"
}

lint_dir() { lint_paths "$1"; }

# --- --selftest: round-trip known-bad/good fixtures through the real checks ----
# The fixtures inside selftest deliberately embed literal '$HOME', '${HOME}',
# '$SCRATCH' etc. as the bad-input strings the lint must catch — they must NOT
# expand, so SC2016 is a false positive throughout this function.
# shellcheck disable=SC2016
selftest() {
    local tmp; tmp="$(mktemp -d)"
    local rc=0

    _mk() {  # _mk <dir> <julia-line> <output-line>  -> minimal wrapper
        mkdir -p "$tmp/$1"
        {
            printf '#!/bin/bash\n#SBATCH --job-name=t\n%s\nmkdir -p slurmlogs\n%s\n' "$3" "$2"
        } > "$tmp/$1/run_nb01.sbatch"
    }
    _expect() {  # _expect <name> <dir> <want-rc> [want-substr]
        local out got
        out="$(lint_dir "$tmp/$2" 2>&1)" && got=0 || got=$?
        if [[ "$got" -ne "$3" ]] || { [[ -n "${4:-}" ]] && ! printf '%s' "$out" | grep -qF "$4"; }; then
            printf '  SELFTEST FAIL: %s (rc=%s want=%s)\n' "$1" "$got" "$3"; rc=1
            printf '%s\n' "$out" | sed 's/^/    | /'
        else
            printf '  ok: %s\n' "$1"
        fi
    }

    local good_out='#SBATCH --output=slurmlogs/%x-%j.out'
    local good_jl='script -e -q -c "julia x.jl" /dev/null'

    _mk clean   "$good_jl"          "$good_out";  _expect "clean wrapper passes"          clean   0 "clean"
    _mk u       'julia -u x.jl'     "$good_out";  _expect "julia -u caught"               u       1 "julia -u"
    _mk u2      'julia -O3 -u x.jl' "$good_out";  _expect "julia -O3 -u (non-adjacent)"   u2      1 "julia -u"
    _mk au      "$good_jl"          "$good_out";  printf 'export JULIA_NUM_THREADS=auto\n' >>"$tmp/au/run_nb01.sbatch"
    _expect "JULIA_NUM_THREADS=auto caught" au 1 "=auto"
    _mk home    "$good_jl" '#SBATCH --output=$HOME/workspace/slurmlogs/%x.out'
    _expect "\$HOME caught"          home    1 "unexpanded shell var"
    _mk braced  "$good_jl" '#SBATCH --output=${HOME}/logs/%x.out'
    _expect "\${HOME} braced caught" braced  1 "unexpanded shell var"
    _mk scratch "$good_jl" '#SBATCH --output=$SCRATCH/logs/%x.out'
    _expect "\$SCRATCH caught"       scratch 1 "unexpanded shell var"
    _mk cmt     "$good_jl"          "$good_out"
    printf '# note: "julia -u" errors at startup\n# $HOME is NOT expanded in #SBATCH\n' >>"$tmp/cmt/run_nb01.sbatch"
    _expect "comment mentions are exempt" cmt 0 "clean"
    mkdir -p "$tmp/lrc"
    printf '#!/bin/bash\n#SBATCH --output=slurmlogs/%%x.out\n# as on Perlmutter\n%s\n' "$good_jl" \
        > "$tmp/lrc/run_nb01_lrc.sbatch"
    _expect "comment Perlmutter is exempt" lrc 0 "clean"

    # Mycelia discovery: templates live under lawrencium/, and render_sbatch
    # stamps `# Site profile: lawrencium` even when the filename is not *_lrc.
    # `#SBATCH` is still a directive, so a milan constraint is a real hit.
    mkdir -p "$tmp/lawpath/lawrencium"
    printf '#!/bin/bash\n#SBATCH --output=slurmlogs/%%x.out\n#SBATCH --constraint=milan\n%s\n' "$good_jl" \
        > "$tmp/lawpath/lawrencium/cpu.sbatch"
    _expect "#SBATCH milan caught" lawpath 1 "NERSC-ism"
    mkdir -p "$tmp/profile"
    printf '#!/bin/bash\n# Site profile: lawrencium\n#SBATCH --output=slurmlogs/%%x.out\nexport JULIA_DEPOT_PATH=/global/cfs/cdirs/m4269/Perlmutter/.julia_depot\n%s\n' "$good_jl" \
        > "$tmp/profile/rendered.sbatch"
    _expect "Perlmutter path in code caught" profile 1 "NERSC-ism"
    mkdir -p "$tmp/profileok"
    printf '#!/bin/bash\n# Site profile: lawrencium\n#SBATCH --output=slurmlogs/%%x.out\n%s\n' "$good_jl" \
        > "$tmp/profileok/rendered.sbatch"
    _expect "site profile lawrencium clean" profileok 0 "clean"

    rm -rf "$tmp"
    if [[ $rc -eq 0 ]]; then printf 'SELFTEST: all cases pass.\n'; fi
    return "$rc"
}

# --- entrypoint ----------------------------------------------------------------
case "${1:-}" in
    -h|--help)  usage; exit 0 ;;
    --selftest) preflight; selftest; exit $? ;;
    *)
        preflight
        if [[ $# -eq 0 ]]; then
            lint_paths "$DEFAULT_TARGET"
        else
            lint_paths "$@"
        fi
        exit $?
        ;;
esac
