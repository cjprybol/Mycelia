# Topological Data Analysis (TDA) Integration Plan

This document turns “TDA-native assembly” into concrete Mycelia work items, scoped to the current architecture (Rhizomorph graphs, quality-aware fixed-length graphs, staged tests, and the existing self-optimizing assembly vision).

## Goals

- Make topology a first-class signal for assembly graph cleaning and parameter selection.
- Provide a stable, backend-agnostic API for TDA features (so we can start with simple graph invariants and later add persistent homology via a Julia TDA library).
- Integrate TDA outputs into assembly QC/optimization without forcing Ripserer (or any specific backend) on all users.

## Design Principles

1. **Thin adapter layer**
   - Mycelia APIs should return Mycelia-owned types (`TDAConfig`, `TDAMetrics`, `TDARunSummary`), not backend-specific diagram types.
2. **Backend optionality**
   - Start with graph-theoretic Betti approximations (connected components + cycle rank) that require no new dependencies.
   - Add persistent homology later as an optional backend (planned: Ripserer.jl).
3. **Filtration-first**
   - TDA should operate on *filtrations* we already use heuristically (coverage/quality/confidence thresholds), producing curves and stability summaries.
4. **Assembly-native hooks**
   - TDA is computed on the same assembly graphs we already build (k-mer/qualmer/string/FASTA/FASTQ) and is optionally included in QC and optimization results.

## Proposed API Surface (initial)

Implemented as a new core file: `src/tda.jl`.

### Types

- `TDAConfig`
  - `max_dim` (start with 1)
  - `thresholds` (coverage/quality/confidence thresholds)
  - `max_points` (cap for expensive backends; unused in invariant-only backend)
  - `backend` (`:graph_betti` initially; later `:ripserer`)
- `TDAMetrics`
  - `thresholds`, `betti0`, `betti1`
  - persistence summaries (placeholder until a PH backend lands)
- `TDARunSummary`
  - `config`, `graph_stats`, `metrics`

### Functions

- `tda_betti_numbers(g)`
  - Returns `(betti0, betti1)` for the underlying undirected graph (Betti₀ = #components; Betti₁ = cycle rank).
- `tda_betti_curves(g; thresholds, vertex_weights)`
  - Computes Betti curves across a vertex-weight filtration (e.g. coverage thresholding).
- `tda_on_graph(g, cfg; vertex_weights)`
  - Convenience wrapper returning a `TDARunSummary`.
- `tda_graph_score(metrics; α, β)`
  - Scalar “simplicity” objective used by graph cleaning / parameter selection.

## Phase Plan

### Phase 0 — Optional persistent homology backend

**Outcome:** Add an optional PH backend (planned: Ripserer.jl) behind the same Mycelia API.

- [ ] Add Ripserer as an optional dependency (prefer Julia extensions / weakdeps).
- [ ] Implement backend dispatch in `tda_on_graph` / `tda_*` functions.
- [ ] Add toy PH tests (circle vs noise; cycle graph vs path graph).

**Acceptance:** Core tests pass without Ripserer; PH tests run only when the backend is available.

### Phase 1 — Graph-level TDA for assembly graphs

**Outcome:** For any assembly graph snapshot, compute topology summaries across thresholds.

- [ ] Implement coverage/quality/confidence vertex-weight extraction adapters for:
  - qualmer graphs (old pipeline: `coverage`, `joint_probability`)
  - Rhizomorph graphs (new pipeline: coverage computed from evidence)
- [ ] Define standard “filtration recipes” per graph type (what weight drives thresholding).

**Acceptance:** Tutorials can print a `TDARunSummary` alongside existing QC metrics on small datasets.

### Phase 2 — TDA-guided graph cleaning

**Outcome:** Replace ad-hoc threshold loops with a topology objective.

- [ ] Implement `choose_threshold_via_tda(g; cfg, vertex_weights)` returning best threshold + summary.
- [ ] Replace heuristic “increase coverage threshold until topology changes” loops with the above.
- [ ] Add synthetic regression tests: noisy spur removal without fragmenting the main contig.

**Acceptance:** Existing assembly/QC metrics do not regress on synthetic fixtures; graph cleaning behavior is deterministic.

### Phase 3 — Topology-aware self-optimizing assembly

**Outcome:** Use TDA metrics in the scalar objective for parameter search/learning.

- [ ] Extend assembly result logging to record `TDARunSummary` (stored in `Mycelia.Rhizomorph.AssemblyResult.assembly_stats` initially).
- [ ] Define an `assembly_objective(result; weights)` that mixes classical QC metrics and TDA simplicity.
- [ ] Add a simple black-box search driver (random / Latin hypercube) over `Mycelia.Rhizomorph.AssemblyConfig`.

**Acceptance:** `optimize_assembly` returns the best config on small simulated datasets and logs topology metrics.

### Phase 4 — Local TDA for bubble/repeat/strain resolution

**Outcome:** Use topology in *local subgraphs* to guide bubble/repeat resolution.

- [ ] Bubble extraction API returns subgraphs suitable for per-bubble scoring.
- [ ] Local TDA scoring combined with read support / quality to pick resolutions or report multi-allelic outputs.

**Acceptance:** Synthetic bubble fixtures show improved path choice vs coverage-only heuristics.

### Phase 5 — Docs and tutorials

**Outcome:** Make TDA discoverable and reproducible.

- [ ] Add a tutorial: “TDA for Qualmer Graphs”
- [ ] Add a narrative doc page: “Topological Assembly Optimization in Mycelia”
- [ ] Document all new TDA APIs and their intended usage

## Implementation Status (as of this doc)

- `src/tda.jl`: implemented graph Betti invariants and filtration curves (graph-only backend).
- `test/4_assembly/tda_metrics_test.jl`: tests exist for cycle vs path and filtration sanity.
