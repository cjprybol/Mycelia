# opt2 hot-loop allocation reduction

Measurement of `Mycelia._viterbi_correct_observation` allocations before
(pre-opt2, `origin/master`) and after (opt2, this branch) the typed-transition
refactor. opt2 replaced `Vector{Any}` + `Dict{Symbol,Any}` transition
representations in the hot Viterbi decode loop with a concrete-eltype
`_Transition{L,E}` struct built in one fused pass
(`_valid_transitions_and_total`), and hoisted per-depth diagnostics counters
from boxed values to typed `Int` locals.

## Result

| | Bytes allocated | Commit |
| --- | --- | --- |
| Before (`origin/master`) | 44,528 | `0804d0ecc` (perf: parallelize `:scalable` corrector decode...) |
| After (this branch, opt2) | 32,816 | `3eecafce` (fix: type-stabilize hoisted per-depth diagnostics locals with `::Int`) |

**Reduction: 11,712 bytes (26.30%)** for a single call to
`_viterbi_correct_observation` on the fixture below.

Both measurements are exactly reproducible (7 repeated post-warmup calls on
each branch returned identical byte counts every time — 32,816 bytes on
`opt2-alloc`, 44,528 bytes on `origin/master`; zero variance).

## Output identity confirmation

Ran the identical fixture and call on both branches and compared every field
of the result. All of the following matched byte-for-byte / bit-for-bit
between before and after:

- `path.steps[].vertex_label` — `[ATG, TGC, GCG, CGT, GTA, TAA]`
- `path.steps[].strand` — all `Forward`
- `path.steps[].probability` — `[1.0, 1.0, 0.6666666666666666, 1.0, 1.0, 1.0]`
- `path.steps[].cumulative_probability` —
  `[1.0, 1.0, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666]`
- `path.total_probability` — `0.6666666666666666`
- `result.score` — `-0.5863711534711905`
- `result.diagnostics` — all 21 key/value pairs identical, including
  `successor_bounded => 1` (confirms the `_top_b_transitions` bounded-branch
  code path fired on this fixture, with the same count on both branches)

**The allocation reduction did not change decode output on this fixture.** This
also matches
the pinned expected values in
`test/4_assembly/viterbi_decode_golden_test.jl` (captured pre-refactor at
commit `047142aa8`), which continues to pass on this branch.

## Fixture

Same fixture as the opt2 golden-master characterization test
(`test/4_assembly/viterbi_decode_golden_test.jl`):

- k=3 kmer graph built from 3 reads: `ATGCGTAA` (x2), `ATGCTTAA` (x1) — reads
  share prefix `ATG,TGC`, branch at `TGC` into two out-edges (`GCG`
  coverage=2, `GCT` coverage=1), then reconverge at `TAA`.
- `mode = :singlestrand`, alphabet resolved to `:DNA`.
- Observation: 6-kmer path `[ATG, TGC, GCG, CGT, GTA, TAA]`.
- `Mycelia.ViterbiCorrectionConfig(max_successors_per_state = 1)` — forces
  `_top_b_transitions` to bound the 2-way branch to one successor per depth,
  so the fixture exercises both opt2-touched call sites: the fused
  `_valid_transitions_and_total` scan and the `_top_b_transitions`
  field-access sort.
- Graph converted via `Mycelia.Rhizomorph.weighted_graph_from_rhizomorph`
  before decode, matching production call shape.
- Measurement protocol: build the fixture once, run
  `_viterbi_correct_observation` once to force JIT compilation (warm-up, not
  measured), then measure `@allocated` on the same fixture object for the
  actual call. Repeated 7x per branch to confirm determinism.

## Machine

- Darwin 25.5.0 (arm64), Apple M5 Max
- Julia 1.10.10
- Project instantiated independently in both worktrees
  (`.worktrees/opt2-alloc` for after, a scratch worktree of `origin/master`
  for before) from the same `Project.toml` dependency versions available at
  measurement time.

## Caveats (read before generalizing this number)

- **Micro-measurement, single fixture.** This is one small graph (3 reads,
  k=3, 6-step observation, single branch point) run through one function.
  It is not a benchmark of end-to-end corrector runtime, and it is not
  necessarily representative of allocation behavior on larger graphs, longer
  reads, wider branching factor, or `max_successors_per_state > 1` (which
  exercises more of the beam-retention code paths that this fixture, by
  design of `max_successors_per_state = 1`, does not stress).
- **Directional, not a throughput claim.** A 26% reduction in bytes allocated
  per call does not imply a 26% wall-clock speedup — allocation counts and
  GC-driven wall time are correlated but not linearly so, and this report
  makes no wall-clock claim.
- **`@allocated` measures only the marked call, post-warmup.** Compilation
  allocations are excluded by construction (the warm-up call absorbs them).
  Allocations inside called library code (e.g. `MetaGraphsNext`,
  `Rhizomorph` graph accessors invoked from within the decode) are included
  if triggered by this call, but this report does not attribute the 11,712
  byte reduction to any specific sub-expression beyond what opt2's commit
  history documents (typed `_Transition{L,E}` construction, hoisted `Int`
  diagnostics locals).
- **No statistical variance reported** because none was observed across 7
  repeated post-warmup calls per branch (identical byte count every time on
  each branch) — this is expected for a purely deterministic, single-threaded
  decode over a fixed fixture, not evidence about variance on other
  workloads.
