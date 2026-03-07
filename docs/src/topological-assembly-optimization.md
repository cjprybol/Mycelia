# Topological Assembly Optimization in Mycelia

Mycelia's TDA layer treats topology as a quantitative assembly signal rather
than a purely descriptive graph property. The current implementation is backend
agnostic and starts with graph Betti invariants, which makes the workflow cheap
to run on small and medium assembly graphs without introducing a persistent
homology dependency.

## Core idea

For an assembly graph, Mycelia measures:

- `betti0`: the number of connected components, which is a proxy for assembly fragmentation
- `betti1`: the cycle rank, which is a proxy for unresolved repeats, bubbles, or circular structure

These values are evaluated across threshold filtrations driven by coverage,
confidence, or quality weights assigned to vertices. In practice this means you
can ask not only "what is the topology of the current graph?" but also "how does
the topology change as weakly supported vertices are removed?"

## Current API

The graph-only backend lives in `src/tda.jl` and exposes:

- `Mycelia.TDAConfig`
- `Mycelia.tda_betti_numbers`
- `Mycelia.tda_betti_curves`
- `Mycelia.tda_on_graph`
- `Mycelia.tda_graph_score`

The returned `Mycelia.TDARunSummary` packages the configuration, lightweight
graph statistics, and threshold-wise Betti curves.

## How to use it

The intended workflow is:

1. Build an assembly graph with an existing Mycelia or Rhizomorph API.
2. Extract a vertex-aligned support signal such as coverage, evidence count, or quality.
3. Evaluate `Mycelia.tda_on_graph(graph, config; vertex_weights=weights)`.
4. Compare candidate thresholds or graph-cleaning outputs with `Mycelia.tda_graph_score`.

Lower scores are better. A useful threshold often removes cycles without
increasing the number of connected components.

## Tutorial

For a concrete walkthrough on synthetic assembly-like graphs, see
[Tutorial 21: Topological Assembly Optimization](generated/tutorials/21_tda_topological_assembly_optimization.md).

## Roadmap

The long-term TDA plan is documented in `planning-docs/TDA_INTEGRATION_PLAN.md`.
The next major step is an optional persistent homology backend layered behind
the same API, so current graph-Betti workflows remain stable while richer
topological summaries become available.
