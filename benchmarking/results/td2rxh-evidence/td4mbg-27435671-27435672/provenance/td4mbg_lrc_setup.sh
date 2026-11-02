#!/usr/bin/env bash
# LRC setup for PR #487 measurement (amendment 4): BASE ca9eb736 + FIX dae40fea, one shared Manifest.
set -euo pipefail
R=/global/home/groups-sw/pc_mfnanofabio/cjprybol/Mycelia
FIX="${R}/.worktrees/td-4mbg-fix"
BASE="${R}/.worktrees/td-4mbg-base"
export JULIA_DEPOT_PATH=/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot
export LD_LIBRARY_PATH=""
export TMPDIR=/global/scratch/users/cjprybol/tmp_td4mbg
mkdir -p "${TMPDIR}"

git -C "${R}" fetch origin
git -C "${R}" cat-file -e 'dae40feaf6973c087ce956d08a25a7f7a46d5bfd^{commit}'
git -C "${R}" cat-file -e 'ca9eb736^{commit}'
[[ -e "${FIX}" ]] || git -C "${R}" worktree add --detach "${FIX}" dae40feaf6973c087ce956d08a25a7f7a46d5bfd
[[ -e "${BASE}" ]] || git -C "${R}" worktree add --detach "${BASE}" ca9eb736
echo "FIX  HEAD $(git -C "${FIX}" rev-parse HEAD)"
echo "BASE HEAD $(git -C "${BASE}" rev-parse HEAD)"

[[ -f "${FIX}/Manifest.toml" ]] || cp "${R}/Manifest.toml" "${FIX}/Manifest.toml"
JULIA_PKG_PRECOMPILE_AUTO=0 "${HOME}/.juliaup/bin/julia" +lts --startup-file=no --project="${FIX}" \
  -e 'import Pkg; Pkg.update("UMAP"; preserve = Pkg.PRESERVE_ALL); Pkg.resolve(); Pkg.instantiate(); println("INSTANTIATE_OK")'
cp "${FIX}/Manifest.toml" "${BASE}/Manifest.toml"
sha256sum "${FIX}/Manifest.toml" "${BASE}/Manifest.toml"
echo "FIX dirty:  [$(git -C "${FIX}" status --porcelain --untracked-files=no)]"
echo "BASE dirty: [$(git -C "${BASE}" status --porcelain --untracked-files=no)]"
echo "SETUP_DONE"
