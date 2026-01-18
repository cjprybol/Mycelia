# Development Code Triage and Coverage Plan

This document captures the current triage of `src/development` and the plan to
make this code production-ready with full test coverage. Per the current
coverage roadmap, we will leave this code in place for now and return after
core coverage (steps 2 and 3) and external wrapper coverage (step 4) are done.

## Goals
- Migrate all functional and complete code into `src/` and add tests.
- Make a concrete plan to uncomment and finish the commented-out code, then test.
- Refactor and generalize any weak implementations, then test.

## Current Triage Snapshot
### Active code (not included in `src/Mycelia.jl`, no tests)
- `src/development/genomic-graph-algorithms.jl`
  - Active implementations for genomic Dijkstra, bidirectional search, rerouting.
  - Likely belongs under Rhizomorph algorithms or a new graph algorithms module.
  - TODO: isolate reusable pieces (path cost, reconstruction) and add unit tests.
- `src/development/intelligent-assembly.jl`
  - Active implementation with duplicated helpers (sparsity, prime k selection).
  - Error-correction hooks are placeholders; needs integration with Viterbi/error-correction code.
  - TODO: refactor to remove duplication, wire to production error correction, then test.
- `src/development/cross-validation.jl`
  - Active pipeline with heavy IO and assembly dependency chain.
  - TODO: extract reusable partitioning + metrics into `src/`, keep pipeline under
    benchmarking or HPC-only workflow tests.

### Placeholder/commented-out code (non-executable)
Fully commented files (no executable lines detected):
- `src/development/annealing-error-correction.jl`
- `src/development/circular-genome-reconstruction.jl`
- `src/development/neo4jl.jl`
- `src/development/pangenome-core-genome.jl`
- `src/development/performance-benchmarks.jl`
- `src/development/protein-clustering-pangenome.jl`
- `src/development/reinforcement-learning-comparison.jl`
- `src/development/reinforcement-learning-mcts.jl`
- `src/development/reinforcement-learning-pomdp.jl`
- `src/development/reinforcement-learning-rl-jl.jl`
- `src/development/reinforcement-learning.jl`
- `src/development/strain-resolved-types.jl`
- `src/development/therapeutic-optimization.jl`

### Inventory (quick status pass)
- Active (callable but not included in Mycelia): `cross-validation.jl`, `genomic-graph-algorithms.jl`, `intelligent-assembly.jl`
- Commented/prototype only: `annealing-error-correction.jl`, `circular-genome-reconstruction.jl`, `neo4jl.jl`,
  `pangenome-core-genome.jl`, `performance-benchmarks.jl`, `protein-clustering-pangenome.jl`,
  `reinforcement-learning*.jl`, `strain-resolved-types.jl`, `therapeutic-optimization.jl`

## Planned Work (deferred until after steps 2–4)
### Phase 1: Promotion + Tests
- Migrate `genomic-graph-algorithms.jl` into `src/` (likely `sequence-graphs-next`
  or a new `graph-algorithms.jl` file) and add tests under `test/4_assembly`.
- Split `cross-validation.jl`:
  - Keep k-fold partitioning and metrics in `src/` with unit tests.
  - Keep full pipeline under `benchmarking/` or a workflow module with opt-in
    tests for HPC (`MYCELIA_RUN_ALL=true`).
- Refactor `intelligent-assembly.jl` to remove placeholder hooks and integrate
  production error-correction or Viterbi logic; add deterministic tests.

### Phase 2: Uncomment + Implement
- For each fully commented file, decide one of:
  - Convert to active code with real implementations.
  - Move to planning docs and remove from `src/` if it is not ready.
- When activating:
  - Add docstrings and small deterministic tests with `StableRNGs`.
  - Avoid external tools in core tests; use HPC-gated integration tests instead.

### Phase 3: Refactor + Generalize
- Replace ad-hoc or hard-coded implementations with reusable utilities.
- Add tests for edge cases and error handling.
- Ensure any external IO or tool integration is isolated and testable.

## Sequencing and Gate
1) Complete core coverage work (zero-coverage modules + largest gaps).
2) Complete external wrapper coverage on HPC.
3) Then start migrating `src/development` code, re-running coverage at each phase.

## Triage TODOs (to schedule after coverage pass)
- [ ] Decide destination module(s) for `genomic-graph-algorithms.jl` and migrate active functions.
- [ ] Refactor `intelligent-assembly.jl` to reuse production helpers and remove placeholders; add unit tests.
- [ ] Split `cross-validation.jl` into reusable library helpers vs. HPC pipeline; add tests for helpers.
- [ ] For each commented prototype, decide: implement or archive to planning docs; document the decision.
- [ ] If implementing commented prototypes, add docstrings + deterministic tests (StableRNGs).
