#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: publish_hpc_results.sh [--input <path>] [--branch <name>] [--remote <name>] [--no-push]

Publish hpc-ci/hpc-results.json to a lightweight results branch with Shields badge endpoints.

Environment:
  HPC_RESULTS_DIR                Defaults to hpc-ci
  HPC_RESULTS_BRANCH             Defaults to hpc-results
  HPC_RESULTS_REMOTE             Defaults to origin
  HPC_RESULTS_REPO_URL           Optional explicit push/pull URL
  HPC_RESULTS_REPOSITORY_SLUG    Optional owner/repo override for badge URLs
EOF
}

INPUT_JSON=""
RESULTS_BRANCH="${HPC_RESULTS_BRANCH:-hpc-results}"
REMOTE_NAME="${HPC_RESULTS_REMOTE:-origin}"
PUSH_RESULTS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT_JSON="$2"
      shift 2
      ;;
    --branch)
      RESULTS_BRANCH="$2"
      shift 2
      ;;
    --remote)
      REMOTE_NAME="$2"
      shift 2
      ;;
    --no-push)
      PUSH_RESULTS=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

ROOT_DIR="$(git rev-parse --show-toplevel)"
INPUT_JSON="${INPUT_JSON:-${HPC_RESULTS_DIR:-${ROOT_DIR}/hpc-ci}/hpc-results.json}"
INPUT_JSON="$(cd "$(dirname "${INPUT_JSON}")" && pwd)/$(basename "${INPUT_JSON}")"

if [[ ! -f "${INPUT_JSON}" ]]; then
  echo "Input JSON not found: ${INPUT_JSON}" >&2
  exit 1
fi

REPO_URL="${HPC_RESULTS_REPO_URL:-$(git config --get "remote.${REMOTE_NAME}.url")}"
if [[ -z "${REPO_URL}" ]]; then
  echo "Could not determine remote URL for ${REMOTE_NAME}" >&2
  exit 1
fi

repository_slug="${HPC_RESULTS_REPOSITORY_SLUG:-}"
if [[ -z "${repository_slug}" ]]; then
  if [[ "${REPO_URL}" =~ github\.com[:/]([^/]+/[^/.]+)(\.git)?$ ]]; then
    repository_slug="${BASH_REMATCH[1]}"
  else
    repository_slug=""
  fi
fi

scratch_dir="$(mktemp -d)"
cleanup() {
  rm -rf "${scratch_dir}"
}
trap cleanup EXIT

branch_exists=0
if git ls-remote --exit-code --heads "${REPO_URL}" "${RESULTS_BRANCH}" >/dev/null 2>&1; then
  branch_exists=1
fi

if [[ "${branch_exists}" -eq 1 ]]; then
  git clone --quiet --single-branch --branch "${RESULTS_BRANCH}" "${REPO_URL}" "${scratch_dir}/repo"
else
  git clone --quiet --no-checkout "${REPO_URL}" "${scratch_dir}/repo"
  (
    cd "${scratch_dir}/repo"
    git checkout --orphan "${RESULTS_BRANCH}"
    find . -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
  )
fi

export HPC_RESULTS_BRANCH="${RESULTS_BRANCH}"
export HPC_RESULTS_REPOSITORY_SLUG="${repository_slug}"

julia --project="${ROOT_DIR}" -e 'import Pkg; Pkg.instantiate()'

julia --project="${ROOT_DIR}" \
  "${ROOT_DIR}/ci/hpc/render_hpc_results_site.jl" \
  "${INPUT_JSON}" \
  "${scratch_dir}/repo"

commit_sha="$(
  julia --project="${ROOT_DIR}" -e 'import JSON; println(string(JSON.parsefile(ARGS[1])["commit"]))' "${INPUT_JSON}"
)"

(
  cd "${scratch_dir}/repo"

  git add README.md .nojekyll latest-hpc-results.json latest-tests.json latest-benchmarks.json latest-meta.json "${commit_sha}"

  if git diff --cached --quiet; then
    echo "No hpc-results changes to publish."
    exit 0
  fi

  if [[ -z "$(git config user.name)" ]]; then
    git config user.name "Mycelia HPC CI"
  fi
  if [[ -z "$(git config user.email)" ]]; then
    git config user.email "actions@users.noreply.github.com"
  fi

  git commit -m "update: publish HPC CI results for ${commit_sha}"

  if [[ "${PUSH_RESULTS}" -eq 1 ]]; then
    git push origin "${RESULTS_BRANCH}"
  else
    echo "Prepared ${RESULTS_BRANCH} branch in ${scratch_dir}/repo"
  fi
)
