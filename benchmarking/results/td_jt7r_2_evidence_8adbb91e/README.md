# td-jt7r.2 evidence generated from `8adbb91e1e17ecfec3391248c7ab5dcd0dc1dd42`

This generation-specific snapshot preserves the nine byte-identical artifacts
produced from a clean tracked worktree at source commit
`8adbb91e1e17ecfec3391248c7ab5dcd0dc1dd42`. The subsequent artifact-only
commit that adds this directory does not change the benchmark implementation.
Both embedded manifests therefore retain source
`code_sha=8adbb91e1e17ecfec3391248c7ab5dcd0dc1dd42` instead of claiming the later
snapshot commit.

## Branching/frontier runtime calibration

- The runtime matrix measures warmed pair-HMM runtime against topology/frontier
  work at 250, 500, and 1,200 bp. It never computes or uses correction accuracy
  (`accuracy_metric_used=false`).
- The 64-vertex, maximally branching 500 bp control is rejected with
  `probe_reason=work_limit`.
- The 10,001-vertex linear 500 bp control is admitted with
  `probe_reason=complete`; all six warmup/measurement rows are complete,
  full, trace-valid, and nontruncated.
- Its warmed p95 is `492.580` ms, or `1.0047`x the same-run 2,001-vertex
  control, below the 1,000 ms and 3x gates.
- The figure is preserved as both SVG and PNG. Overlapping point labels remain
  a presentation-only limitation; the CSV tables are the authoritative values.

## Fixed toy and oracle matrix

- Fixture: 2,000 bp reference, 1,200 bp source reads, 8x target coverage,
  5% nanopore error, observed read lengths `1190-1214` bp.
- Nanopore completed `26` indel decodes above the initial
  `k=3` rung and recorded requested/attempted/completed/truncated/engaged
  totals of `207/26/26/0/15`.
- Nanopore finished in `100.056` s and achieved identity `0.1`, versus
  `0.0655` for the identical-read Illumina arm. All 19 acceptance checks pass.
- Explicit Illumina and the default substitution oracle are byte-identical to
  each other and to the pre-wiring SHA-256
  `d36e3b6a10685346aa7b0238b48b4ab7fcefbed88f82cad7d959b0a831cdd311`.

## Claim boundary

This is **interim engineering validation**, not the Nature Methods four-tier
benchmark or a manuscript H5 result. “1,000 bp+” refers to input read length,
not assembled output span. The toy assemblies remain fragmented: nanopore has
`369` contigs with largest span `200` bp; Illumina has `526` contigs with
largest span `131` bp. Do not generalize these deterministic toy results to
production assembly accuracy.

`artifact-index.json` records SHA-256 digests, byte sizes, CSV schemas and row
counts, and the exact validation values. It omits its own digest by construction.
