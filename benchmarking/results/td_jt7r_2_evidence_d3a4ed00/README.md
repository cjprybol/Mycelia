# td-jt7r.2 evidence generated from `d3a4ed001b364e8e0dc50c5161d94837724ed972`

This generation-specific snapshot preserves twelve source-generation/capture
artifacts produced while the tracked worktree was clean at source commit
`d3a4ed001b364e8e0dc50c5161d94837724ed972`. The subsequent artifact-only
commit that adds this directory does not change the benchmark implementation.
All three embedded manifests therefore retain source
`code_sha=d3a4ed001b364e8e0dc50c5161d94837724ed972` instead of claiming the later snapshot commit.

The source includes all pre-merge review fixes through `d3a4ed00`, including
initial deletion-band pruning in both the score-free topology probe and the
scored pair-HMM. This snapshot is regenerated from that reviewed code rather
than extrapolated from the earlier `8adbb91e` measurements.

## Branching/frontier runtime calibration

- The runtime matrix measures warmed pair-HMM runtime against topology/frontier
  work at 250, 500, and 1,200 bp. It never computes or uses correction accuracy
  (`accuracy_metric_used=false`).
- The 64-vertex, maximally branching 500 bp control is rejected with
  `probe_reason=work_limit`.
- The 10,001-vertex linear 500 bp control is admitted with
  `probe_reason=complete`; all six warmup/measurement rows are complete,
  full, trace-valid, and nontruncated.
- Its warmed p95 is `352.941` ms, or `0.9819`x the same-run
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
- Nanopore finished in `69.733` s and achieved
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
- `oracle_matrix_capture.jl` is a capture-only wrapper stored by the later
  artifact commit. It included the committed oracle test unchanged while
  running against `d3a4ed00`; it was not itself a tracked source file at
  that commit and does not relabel the measured implementation.

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
