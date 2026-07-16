# td-jt7r.2 evidence generated from `20400a02027daa3f8b8070ae8a25002ecb12a778`

This generation-specific snapshot preserves twelve source-generation/capture
artifacts produced while the tracked worktree was clean at source commit
`20400a02027daa3f8b8070ae8a25002ecb12a778`. The subsequent artifact-only
commit that adds this directory does not change the benchmark implementation.
All three embedded manifests therefore retain source
`code_sha=20400a02027daa3f8b8070ae8a25002ecb12a778` instead of claiming the later snapshot commit.

The source includes all pre-merge review fixes through `20400a02`, including
initial deletion-band pruning in both the score-free topology probe and the
scored pair-HMM, publication of valid unchanged-window soft-EM E-step
responsibilities, and checkpoint count/rate consistency validation. This
snapshot is regenerated from that reviewed code and supersedes the immediately
prior `d3a4ed00` measurements rather than extrapolating from them or the earlier
`8adbb91e` generation.

## Branching/frontier runtime calibration

- The runtime matrix measures warmed pair-HMM runtime against topology/frontier
  work at 250, 500, and 1,200 bp. It never computes or uses correction accuracy
  (`accuracy_metric_used=false`).
- The 64-vertex, maximally branching 500 bp control is rejected with
  `probe_reason=work_limit`.
- The 10,001-vertex linear 500 bp control is admitted with
  `probe_reason=complete`; all six warmup/measurement rows are complete,
  full, trace-valid, and nontruncated.
- Its warmed p95 is `328.719` ms, or `0.8796`x the same-run
  2,001-vertex control, below the 1,000 ms and 3x gates.
- The figure is preserved as both SVG and PNG. Overlapping point labels remain
  a presentation-only limitation; the CSV tables are the authoritative values.

## Fixed toy and oracle matrix

- Fixture: 2000 bp reference,
  1200 bp source reads,
  8x target coverage, 0.05
  nanopore error, observed read lengths `1190-1214` bp.
- Nanopore completed `25` indel decodes above the initial
  `k=3` rung and recorded requested/attempted/completed/truncated/engaged
  totals of `207/25/25/0/14`.
- Nanopore finished in `79.749` s and achieved
  identity `0.1`, versus `0.0655` for the
  identical-read Illumina arm. All 19 acceptance checks pass.
- Explicit Illumina and the default substitution oracle are byte-identical to
  each other and to the pre-wiring SHA-256
  `d36e3b6a10685346aa7b0238b48b4ab7fcefbed88f82cad7d959b0a831cdd311`.
- A separate tiny-fixture oracle matrix uses a 60 bp reference, 40 bp reads,
  30x coverage, 5% error, and seeds 1-5. It reports mean nanopore identity
  `0.9467759563` versus Illumina
  `0.9334426230`, with nanopore
  winning or tying `4/5` seeds. Accuracy is a
  validation outcome only and is never an input to the runtime classifier.
- `oracle_matrix_capture.jl` is the exact capture-runner input preserved
  by the later artifact commit. Prior snapshot copies are evidence rather than
  benchmark implementation. The wrapper included the committed oracle test
  unchanged while running against `20400a02` and does not relabel the
  measured implementation.

### Reproducing the oracle capture

The preserved runner's relative paths are valid only from its original ignored
canonical directory, `benchmarking/results/td-jt7r-2-oracle-matrix/`. **Do not
run the archived copy in place:** it would try to replace the checked-in CSVs
before its relocated path lookup fails. To reproduce the exact source-bound
capture without changing this snapshot:

```bash
set -euo pipefail
REPO=$(git rev-parse --show-toplevel)
TMP_PARENT=$(mktemp -d "${TMPDIR:-/tmp}/mycelia-td-jt7r-2.XXXXXX")
SOURCE_WT="$TMP_PARENT/source"
test -f "$REPO/Manifest.toml"
printf '%s  %s\n' \
  bd814ac9d374e22d43fc391f243e1527041728f5cacf2ca7b071cc25b52e09e3 \
  "$REPO/Manifest.toml" | shasum -a 256 -c -
git -C "$REPO" worktree add "$SOURCE_WT" \
  20400a02027daa3f8b8070ae8a25002ecb12a778
cp "$REPO/Manifest.toml" "$SOURCE_WT/Manifest.toml"
mkdir -p "$SOURCE_WT/benchmarking/results/td-jt7r-2-oracle-matrix"
git -C "$REPO" show 438aeba3e59e469c7e1d8b43cd8ca494099bfa37:benchmarking/results/td_jt7r_2_evidence_20400a02/oracle-matrix/oracle_matrix_capture.jl \
  > "$SOURCE_WT/benchmarking/results/td-jt7r-2-oracle-matrix/oracle_matrix_capture.jl"
cd "$SOURCE_WT"
LD_LIBRARY_PATH="" julia --project=. \
  benchmarking/results/td-jt7r-2-oracle-matrix/oracle_matrix_capture.jl
printf 'Reproduction worktree retained at %s\n' "$SOURCE_WT"
# After inspecting the ignored outputs:
# cd "$REPO"
# git worktree remove --force "$SOURCE_WT"
# rmdir "$TMP_PARENT"
```

This restores the hash-bound dependency manifest and runner to their original
locations and writes only to the ignored canonical output directory in the
clean `20400a02` worktree. The rerun manifest will still describe its actual
host and Julia runtime. Exact reproduction therefore requires access to the
ignored Manifest with the recorded SHA-256 above.

## Claim boundary

This is **interim engineering validation**, not the Nature Methods four-tier
benchmark or a manuscript H5 result. “1,000 bp+” refers to input read length,
not assembled output span. The toy assemblies remain fragmented: nanopore has
`372` contigs with largest span
`200` bp; Illumina has `526`
contigs with largest span `131` bp. Do not generalize
these deterministic toy results to production assembly accuracy.

`artifact-index.json` records SHA-256 digests, byte sizes, CSV schemas and row
counts, and the exact validation values. It omits its own digest by construction.
