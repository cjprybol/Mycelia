"""
Topological Data Analysis (TDA) utilities.

This file defines lightweight, backend-agnostic types and functions for computing
topology summaries on Mycelia graphs.

Initial implementation focuses on graph-theoretic Betti invariants:
- Betti₀: number of connected components (on the underlying undirected graph)
- Betti₁: cycle rank (dimension of the cycle space) on the underlying undirected graph

Persistent homology backends (e.g. Ripserer.jl) are planned behind the same API.
See `planning-docs/TDA_INTEGRATION_PLAN.md`.
"""

struct TDAConfig
    max_dim::Int
    thresholds::Vector{Float64}
    max_points::Int
    backend::Symbol

    function TDAConfig(;
        max_dim::Int = 1,
        thresholds::AbstractVector{<:Real} = Float64[],
        max_points::Int = 10_000,
        backend::Symbol = :graph_betti,
    )
        max_dim < 0 && throw(ArgumentError("max_dim must be ≥ 0, got $(max_dim)"))
        max_points < 1 && throw(ArgumentError("max_points must be ≥ 1, got $(max_points)"))

        threshold_vec = sort!(collect(Float64.(thresholds)))
        return new(max_dim, threshold_vec, max_points, backend)
    end
end

struct TDAPersistenceSummary
    longest_persistence_0::Float64
    longest_persistence_1::Float64
    total_persistence_0::Float64
    total_persistence_1::Float64
end

struct TDAMetrics
    thresholds::Vector{Float64}
    betti0::Vector{Int}
    betti1::Vector{Int}
    persistence::Union{Nothing, TDAPersistenceSummary}

    function TDAMetrics(
        thresholds::Vector{Float64},
        betti0::Vector{Int},
        betti1::Vector{Int};
        persistence::Union{Nothing, TDAPersistenceSummary} = nothing,
    )
        length(thresholds) == length(betti0) || throw(ArgumentError("thresholds and betti0 must have the same length"))
        length(thresholds) == length(betti1) || throw(ArgumentError("thresholds and betti1 must have the same length"))
        return new(thresholds, betti0, betti1, persistence)
    end
end

struct TDARunSummary
    config::TDAConfig
    graph_stats::NamedTuple
    metrics::TDAMetrics
end

function _underlying_simple_graph(g::Graphs.AbstractGraph)
    undirected = Graphs.SimpleGraph(Graphs.nv(g))
    for edge in Graphs.edges(g)
        Graphs.add_edge!(undirected, Graphs.src(edge), Graphs.dst(edge))
    end
    return undirected
end

"""
    tda_betti_numbers(g::Graphs.AbstractGraph) -> (betti0::Int, betti1::Int)

Compute graph-theoretic Betti numbers on the underlying undirected simple graph:
- Betti₀ = number of connected components
- Betti₁ = cycle rank (|E| - |V| + Betti₀)
"""
function tda_betti_numbers(g::Graphs.AbstractGraph)
    if Graphs.nv(g) == 0
        return (0, 0)
    end

    undirected = _underlying_simple_graph(g)
    components = Graphs.connected_components(undirected)
    betti0 = length(components)
    betti1 = Graphs.ne(undirected) - Graphs.nv(undirected) + betti0
    return (betti0, max(betti1, 0))
end

"""
    tda_graph_stats(g::Graphs.AbstractGraph) -> NamedTuple

Return lightweight graph statistics useful for logging.
"""
function tda_graph_stats(g::Graphs.AbstractGraph)
    vertex_count = Graphs.nv(g)
    edge_count = Graphs.ne(g)
    directed = Graphs.is_directed(g)
    density =
        if vertex_count <= 1
            0.0
        elseif directed
            edge_count / (vertex_count * (vertex_count - 1))
        else
            edge_count / (vertex_count * (vertex_count - 1) / 2)
        end

    return (
        nv = vertex_count,
        ne = edge_count,
        directed = directed,
        density = density,
    )
end

"""
    tda_betti_curves(g::Graphs.AbstractGraph; thresholds, vertex_weights) -> TDAMetrics

Compute Betti₀ and Betti₁ across a vertex-weight filtration, where for each threshold `t`
we take the induced subgraph on vertices with `vertex_weights[v] >= t`.

This provides an initial “topology signal” without requiring persistent homology.
"""
function tda_betti_curves(
    g::Graphs.AbstractGraph;
    thresholds::AbstractVector{<:Real},
    vertex_weights::AbstractVector{<:Real},
)
    Graphs.nv(g) == length(vertex_weights) || throw(ArgumentError("vertex_weights must have length nv(g)"))

    threshold_vec = sort!(collect(Float64.(thresholds)))
    isempty(threshold_vec) && (threshold_vec = Float64[-Inf])

    betti0 = Vector{Int}(undef, length(threshold_vec))
    betti1 = Vector{Int}(undef, length(threshold_vec))

    for (i, thr) in enumerate(threshold_vec)
        keep_vertices = findall(w -> Float64(w) >= thr, vertex_weights)
        if isempty(keep_vertices)
            betti0[i] = 0
            betti1[i] = 0
            continue
        end

        subgraph, _ = Graphs.induced_subgraph(g, keep_vertices)
        b0, b1 = tda_betti_numbers(subgraph)
        betti0[i] = b0
        betti1[i] = b1
    end

    return TDAMetrics(threshold_vec, betti0, betti1)
end

"""
    tda_on_graph(g::Graphs.AbstractGraph, cfg::TDAConfig; vertex_weights) -> TDARunSummary

Compute a standardized TDA summary for a graph.

`vertex_weights` should typically be a coverage/quality/confidence proxy aligned with
vertex indices. If omitted, all vertices are given weight 1.0.
"""
function tda_on_graph(
    g::Graphs.AbstractGraph,
    cfg::TDAConfig;
    vertex_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    cfg.backend == :graph_betti || throw(ArgumentError("Unsupported TDA backend $(cfg.backend); currently supported: :graph_betti"))
    cfg.max_dim <= 1 || throw(ArgumentError("max_dim=$(cfg.max_dim) is not supported by backend :graph_betti (max_dim ≤ 1)"))

    weights =
        if isnothing(vertex_weights)
            ones(Float64, Graphs.nv(g))
        else
            vertex_weights
        end

    metrics = tda_betti_curves(g; thresholds = cfg.thresholds, vertex_weights = weights)
    return TDARunSummary(cfg, tda_graph_stats(g), metrics)
end

"""
    tda_graph_score(metrics::TDAMetrics; α=1.0, β=1.0) -> Float64

Scalar “simplicity” score for topology-driven graph cleaning / parameter selection.

Lower is better. This initial implementation uses only Betti curves:
- `cycle_penalty` uses the maximum Betti₁ across thresholds (proxy for persistent cycles)
- `fragmentation_penalty` uses the maximum Betti₀ across thresholds

When a persistent homology backend is added, this can be extended to use true
total persistence in H₁.
"""
function tda_graph_score(metrics::TDAMetrics; α::Float64 = 1.0, β::Float64 = 1.0)
    cycle_penalty = α * (isempty(metrics.betti1) ? 0.0 : Float64(maximum(metrics.betti1)))
    fragmentation_penalty = β * (isempty(metrics.betti0) ? 0.0 : Float64(maximum(metrics.betti0)))
    return cycle_penalty + fragmentation_penalty
end
