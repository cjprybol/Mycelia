#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: reset_lustre_tmp.sh [--pool POOL] [--base-dir DIR] [--tmp-path PATH] [--dry-run]

Repoint a Lustre scratch tmp path to a new directory with explicit striping.
Defaults:
  POOL=ddn_hdd
  BASE_DIR=/global/scratch/users/$USER/workspace
  TMP_PATH=$BASE_DIR/tmp

Environment overrides:
  POOL, BASE_DIR, TMP_PATH, STRIPE_COUNT, STRIPE_SIZE, DRY_RUN
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

POOL="${POOL:-ddn_hdd}"
BASE_DIR="${BASE_DIR:-/global/scratch/users/${USER}/workspace}"
TMP_PATH="${TMP_PATH:-${BASE_DIR}/tmp}"
STRIPE_COUNT="${STRIPE_COUNT:-1}"
STRIPE_SIZE="${STRIPE_SIZE:-1M}"
DRY_RUN="${DRY_RUN:-false}"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
NEW_DIR="${BASE_DIR}/tmp_${POOL}_${TIMESTAMP}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pool)
            POOL="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --tmp-path)
            TMP_PATH="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="true"
            shift
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if ! command -v lfs >/dev/null 2>&1; then
    echo "ERROR: lfs not found in PATH." >&2
    exit 1
fi

if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: base dir does not exist: $BASE_DIR" >&2
    exit 1
fi

run_cmd() {
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "+ $*"
    else
        "$@"
    fi
}

run_cmd mkdir -p "$NEW_DIR"
run_cmd lfs setstripe -p "$POOL" -c "$STRIPE_COUNT" -S "$STRIPE_SIZE" "$NEW_DIR"
run_cmd lfs getstripe -d "$NEW_DIR"

if [[ -e "$TMP_PATH" || -L "$TMP_PATH" ]]; then
    BACKUP="${BASE_DIR}/tmp_old_${TIMESTAMP}"
    run_cmd mv "$TMP_PATH" "$BACKUP"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would move existing tmp to $BACKUP"
    else
        echo "Moved existing tmp to $BACKUP"
    fi
fi

run_cmd ln -s "$NEW_DIR" "$TMP_PATH"
if [[ "$DRY_RUN" == "true" ]]; then
    echo "Would update TMPDIR: $TMP_PATH -> $NEW_DIR"
    echo "Would export TMPDIR=\"$TMP_PATH\""
else
    echo "TMPDIR updated: $TMP_PATH -> $NEW_DIR"
    echo "export TMPDIR=\"$TMP_PATH\""
fi
