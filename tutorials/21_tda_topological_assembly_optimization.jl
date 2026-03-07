# # Tutorial 21: Topological Assembly Optimization
#
# This tutorial introduces Mycelia's graph-based TDA utilities on small
# assembly-like graphs. The current backend is intentionally lightweight:
# instead of full persistent homology, it tracks graph-theoretic Betti curves
# across threshold filtrations.
#
# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/21_tda_topological_assembly_optimization.jl", "tutorials/notebooks", execute=false)'
# ```
#
# ## Learning objectives
#
# By the end of this tutorial, you will be able to:
# 1. Interpret `Mycelia.tda_betti_numbers` on assembly-like graphs
# 2. Build a threshold filtration with `Mycelia.TDAConfig`
# 3. Compare candidate graph states with `Mycelia.tda_on_graph`
# 4. Use `Mycelia.tda_graph_score` as a simple topology-aware objective

import Graphs
import Mycelia
import Printf

println("=== Topological Assembly Optimization Tutorial ===")

# ## Step 1: Build synthetic assembly-like graphs
#
# We will compare three common topologies:
# - a clean contig path
# - a bubble graph with one extra cycle
# - a fragmented graph with two disconnected components

function build_path_graph()
    graph = Graphs.SimpleGraph(6)
    for (src, dst) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6))
        Graphs.add_edge!(graph, src, dst)
    end
    return graph
end

function build_bubble_graph()
    graph = Graphs.SimpleGraph(6)
    for (src, dst) in ((1, 2), (2, 3), (3, 5), (2, 4), (4, 5), (5, 6))
        Graphs.add_edge!(graph, src, dst)
    end
    return graph
end

function build_fragmented_graph()
    graph = Graphs.SimpleGraph(6)
    for (src, dst) in ((1, 2), (2, 3), (4, 5), (5, 6))
        Graphs.add_edge!(graph, src, dst)
    end
    return graph
end

graphs = [
    ("clean_path", build_path_graph()),
    ("bubble", build_bubble_graph()),
    ("fragmented", build_fragmented_graph()),
]

println("\nStep 1: Base graph summaries")
for (name, graph) in graphs
    betti0, betti1 = Mycelia.tda_betti_numbers(graph)
    stats = Mycelia.tda_graph_stats(graph)
    println("  $(name): nv=$(stats.nv), ne=$(stats.ne), betti0=$(betti0), betti1=$(betti1)")
end

# ## Step 2: Interpret Betti numbers
#
# `betti0` counts connected components. In assembly terms, that is a proxy for
# fragmentation. `betti1` counts independent cycles, which can indicate repeat
# structures, unresolved bubbles, or circular genomes depending on context.

println("\nStep 2: Interpretation")
println("  clean_path  -> one component, no cycles")
println("  bubble      -> one component, one cycle")
println("  fragmented  -> two components, no cycles")

# ## Step 3: Add a vertex-weight filtration
#
# TDA in Mycelia is filtration-first. Here we simulate per-vertex support
# values, such as coverage or confidence scores, and ask how topology changes as
# we raise the threshold.

thresholds = [0.0, 5.0, 10.0]
config = Mycelia.TDAConfig(thresholds = thresholds)

vertex_weights = Dict(
    "clean_path" => [12.0, 12.0, 11.0, 11.0, 10.0, 10.0],
    "bubble" => [12.0, 12.0, 3.0, 11.0, 11.0, 10.0],
    "fragmented" => [12.0, 11.0, 10.0, 4.0, 3.0, 2.0],
)

function print_summary(name, summary)
    println("  $(name):")
    println("    thresholds = $(summary.metrics.thresholds)")
    println("    betti0     = $(summary.metrics.betti0)")
    println("    betti1     = $(summary.metrics.betti1)")
    println("    score      = $(Mycelia.tda_graph_score(summary.metrics))")
end

println("\nStep 3: Filtration summaries")
summaries = Dict{String, Mycelia.TDARunSummary}()
for (name, graph) in graphs
    summary = Mycelia.tda_on_graph(graph, config; vertex_weights = vertex_weights[name])
    summaries[name] = summary
    print_summary(name, summary)
end

# ## Step 4: See what the thresholds are doing
#
# The bubble graph carries one weakly supported branch. At threshold `5.0`, that
# branch disappears and the cycle vanishes. The fragmented graph behaves
# differently: increasing the threshold removes one entire component, so the
# graph becomes smaller rather than cleaner.

println("\nStep 4: Threshold-by-threshold interpretation")
for threshold_index in eachindex(thresholds)
    threshold = thresholds[threshold_index]
    println("  threshold >= $(threshold)")
    for name in ("clean_path", "bubble", "fragmented")
        summary = summaries[name]
        println(
            "    $(name): betti0=$(summary.metrics.betti0[threshold_index]), " *
            "betti1=$(summary.metrics.betti1[threshold_index])"
        )
    end
end

# ## Step 5: Rank candidate graph states
#
# `Mycelia.tda_graph_score` is a simple scalar objective:
# - higher `betti1` means more cyclic complexity
# - higher `betti0` means more fragmentation
#
# Lower scores are preferred when choosing between graph-cleaning thresholds or
# candidate assembly states.

candidate_scores = [
    (name, Mycelia.tda_graph_score(summary.metrics))
    for (name, summary) in summaries
]
sort!(candidate_scores; by = last)

println("\nStep 5: Topology-aware ranking")
for (rank, (name, score)) in enumerate(candidate_scores)
    Printf.@printf("  %d. %-11s score=%.1f\n", rank, name, score)
end

# ## Step 6: Practical guidance
#
# In real assembly workflows, the same pattern applies to graphs built from
# k-mers, qualmers, or Rhizomorph evidence graphs:
# 1. extract a support signal per vertex
# 2. evaluate Betti curves across thresholds
# 3. choose thresholds that reduce spurious cycles without over-fragmenting the graph

println("\nStep 6: Practical guidance")
println("  1. Use coverage, confidence, or quality as vertex weights.")
println("  2. Look for thresholds where betti1 drops while betti0 stays stable.")
println("  3. Treat persistent cycles as candidates for repeats or unresolved bubbles.")
println("  4. Treat rising betti0 as a warning that graph cleaning is too aggressive.")

# ## Next steps
#
# - Replace these toy graphs with assembly graphs derived from FASTQ or qualmer data.
# - Record `TDARunSummary` alongside assembly QC metrics during optimization.
# - Extend to a persistent homology backend when a backend such as Ripserer is enabled.
