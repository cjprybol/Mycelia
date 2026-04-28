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
            backend::Symbol = :graph_betti
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
            persistence::Union{Nothing, TDAPersistenceSummary} = nothing
    )
        length(thresholds) == length(betti0) ||
            throw(ArgumentError("thresholds and betti0 must have the same length"))
        length(thresholds) == length(betti1) ||
            throw(ArgumentError("thresholds and betti1 must have the same length"))
        return new(thresholds, betti0, betti1, persistence)
    end
end

struct TDARunSummary
    config::TDAConfig
    graph_stats::NamedTuple
    metrics::TDAMetrics
end

function _tda_graph(g)
    graph = g isa MetaGraphsNext.MetaGraph ? g.graph : g
    graph isa Graphs.AbstractGraph ||
        throw(ArgumentError(
            "TDA metrics require a Graphs.AbstractGraph or MetaGraphsNext.MetaGraph"))
    return graph
end

function _underlying_simple_graph(g)
    graph = _tda_graph(g)
    undirected = Graphs.SimpleGraph(Graphs.nv(graph))
    for edge in Graphs.edges(graph)
        Graphs.add_edge!(undirected, Graphs.src(edge), Graphs.dst(edge))
    end
    return undirected
end

function _tda_vertex_weights(g, vertex_weights)
    graph = _tda_graph(g)
    vertex_count = Graphs.nv(graph)

    if isnothing(vertex_weights)
        return ones(Float64, vertex_count)
    elseif vertex_weights isa AbstractVector
        length(vertex_weights) == vertex_count ||
            throw(ArgumentError("vertex_weights must have length nv(g)"))
        return collect(Float64.(vertex_weights))
    elseif vertex_weights isa AbstractDict
        weights = Vector{Float64}(undef, vertex_count)
        if g isa MetaGraphsNext.MetaGraph
            for label in MetaGraphsNext.labels(g)
                haskey(vertex_weights, label) ||
                    throw(ArgumentError("vertex_weights is missing label $(repr(label))"))
                weights[MetaGraphsNext.code_for(g, label)] = Float64(vertex_weights[label])
            end
        else
            for vertex in Graphs.vertices(graph)
                haskey(vertex_weights, vertex) ||
                    throw(ArgumentError("vertex_weights is missing vertex $(vertex)"))
                weights[vertex] = Float64(vertex_weights[vertex])
            end
        end
        return weights
    else
        throw(ArgumentError("vertex_weights must be nothing, an AbstractVector, or an AbstractDict"))
    end
end

function _tda_weight_stats(weights::AbstractVector{<:Real})
    if isempty(weights)
        return (weight_min = NaN, weight_max = NaN, weight_mean = NaN)
    end

    float_weights = Float64.(weights)
    return (
        weight_min = minimum(float_weights),
        weight_max = maximum(float_weights),
        weight_mean = Statistics.mean(float_weights)
    )
end

"""
    tda_betti_numbers(g) -> (betti0::Int, betti1::Int)

Compute graph-theoretic Betti numbers on the underlying undirected simple graph:
- Betti₀ = number of connected components
- Betti₁ = cycle rank (|E| - |V| + Betti₀)

Accepts `Graphs.AbstractGraph` and Rhizomorph `MetaGraphsNext.MetaGraph` inputs.

# Arguments
- `g`: `Graphs.AbstractGraph` or `MetaGraphsNext.MetaGraph` to summarize.

# Returns
- A tuple `(betti0, betti1)` with integer component and cycle-rank counts.

# Example
```julia
graph = Graphs.SimpleGraph(3)
Graphs.add_edge!(graph, 1, 2)
betti0, betti1 = Mycelia.tda_betti_numbers(graph)
```
"""
function tda_betti_numbers(g)
    graph = _tda_graph(g)
    if Graphs.nv(graph) == 0
        return (0, 0)
    end

    undirected = _underlying_simple_graph(graph)
    components = Graphs.connected_components(undirected)
    betti0 = length(components)
    betti1 = Graphs.ne(undirected) - Graphs.nv(undirected) + betti0
    return (betti0, max(betti1, 0))
end

"""
    tda_graph_stats(g) -> NamedTuple

Return lightweight graph statistics useful for logging.

# Arguments
- `g`: `Graphs.AbstractGraph` or `MetaGraphsNext.MetaGraph` to summarize.

# Returns
- A named tuple with `nv`, `ne`, `directed`, and `density`.

# Example
```julia
stats = Mycelia.tda_graph_stats(Graphs.SimpleGraph(4))
```
"""
function tda_graph_stats(g)
    graph = _tda_graph(g)
    vertex_count = Graphs.nv(graph)
    edge_count = Graphs.ne(graph)
    directed = Graphs.is_directed(graph)
    density = if vertex_count <= 1
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
        density = density
    )
end

"""
    tda_betti_curves(g; thresholds, vertex_weights) -> TDAMetrics

Compute Betti₀ and Betti₁ across a vertex-weight filtration, where for each threshold `t`
we take the induced subgraph on vertices with `vertex_weights[v] >= t`.

This provides an initial “topology signal” without requiring persistent homology.

# Arguments
- `g`: `Graphs.AbstractGraph` or `MetaGraphsNext.MetaGraph` to filter.
- `thresholds`: Real-valued filtration thresholds.
- `vertex_weights`: Vector aligned with graph vertex codes, integer-keyed dictionary
  for plain graphs, or label-keyed dictionary for `MetaGraphsNext.MetaGraph`.

# Returns
- `TDAMetrics` containing sorted thresholds and Betti curves.

# Example
```julia
graph = Graphs.SimpleGraph(3)
metrics = Mycelia.tda_betti_curves(
    graph;
    thresholds = [0.0, 1.0],
    vertex_weights = ones(Float64, Graphs.nv(graph))
)
```
"""
function tda_betti_curves(
        g;
        thresholds::AbstractVector{<:Real},
        vertex_weights
)
    graph = _tda_graph(g)
    weights = _tda_vertex_weights(g, vertex_weights)

    threshold_vec = sort!(collect(Float64.(thresholds)))
    isempty(threshold_vec) && (threshold_vec = Float64[-Inf])

    betti0 = Vector{Int}(undef, length(threshold_vec))
    betti1 = Vector{Int}(undef, length(threshold_vec))

    for (i, thr) in enumerate(threshold_vec)
        keep_vertices = findall(w -> Float64(w) >= thr, weights)
        if isempty(keep_vertices)
            betti0[i] = 0
            betti1[i] = 0
            continue
        end

        subgraph, _ = Graphs.induced_subgraph(graph, keep_vertices)
        b0, b1 = tda_betti_numbers(subgraph)
        betti0[i] = b0
        betti1[i] = b1
    end

    return TDAMetrics(threshold_vec, betti0, betti1)
end

"""
    tda_on_graph(g, cfg::TDAConfig; vertex_weights) -> TDARunSummary

Compute a standardized TDA summary for a graph.

`vertex_weights` should typically be a coverage/quality/confidence proxy. For
`Graphs.AbstractGraph`, pass either a vector aligned with vertex indices or a
dictionary keyed by integer vertex id. For Rhizomorph `MetaGraphsNext.MetaGraph`,
pass either a vector aligned with internal vertex codes or a dictionary keyed by
vertex label. If omitted, all vertices are given weight 1.0.

# Arguments
- `g`: `Graphs.AbstractGraph` or `MetaGraphsNext.MetaGraph` to summarize.
- `cfg`: `TDAConfig` controlling backend, maximum dimension, and thresholds.
- `vertex_weights = nothing`: Optional weights for the vertex filtration.

# Returns
- `TDARunSummary` with configuration, graph statistics, and `TDAMetrics`.

# Example
```julia
graph = Graphs.SimpleGraph(3)
cfg = Mycelia.TDAConfig(thresholds = [0.0])
summary = Mycelia.tda_on_graph(graph, cfg)
```
"""
function tda_on_graph(
        g,
        cfg::TDAConfig;
        vertex_weights = nothing
)
    cfg.backend == :graph_betti ||
        throw(ArgumentError("Unsupported TDA backend $(cfg.backend); currently supported: :graph_betti"))
    cfg.max_dim <= 1 ||
        throw(ArgumentError("max_dim=$(cfg.max_dim) is not supported by backend :graph_betti (max_dim ≤ 1)"))

    weights = _tda_vertex_weights(g, vertex_weights)

    metrics = tda_betti_curves(g; thresholds = cfg.thresholds, vertex_weights = weights)
    return TDARunSummary(cfg, tda_graph_stats(g), metrics)
end

"""
    tda_metric_rows(summary::TDARunSummary; kwargs...) -> Vector{NamedTuple}

Convert a `TDARunSummary` into stable, table-ready rows. One row is emitted per
filtration threshold, with graph statistics, Betti values, summary score, and a
provenance tuple that records the TDA configuration and filtration choices.

# Arguments
- `summary`: `TDARunSummary` returned by `tda_on_graph`.
- `graph_id = "graph"`: Identifier copied into every row.
- `filtration = :vertex_weight_threshold`: Filtration recipe name.
- `weight_name = :uniform`: Name of the vertex weight signal.
- `weight_stats = (weight_min = NaN, weight_max = NaN, weight_mean = NaN)`: Weight summary.
- `provenance = (; )`: Extra user provenance merged into the generated provenance tuple.

# Returns
- `Vector{NamedTuple}` with one table-ready row per threshold.

# Example
```julia
graph = Graphs.SimpleGraph(3)
summary = Mycelia.tda_on_graph(graph, Mycelia.TDAConfig(thresholds = [0.0]))
rows = Mycelia.tda_metric_rows(summary; graph_id = "assembly_graph")
```
"""
function tda_metric_rows(
        summary::TDARunSummary;
        graph_id = "graph",
        filtration::Symbol = :vertex_weight_threshold,
        weight_name::Symbol = :uniform,
        weight_stats::NamedTuple = (weight_min = NaN, weight_max = NaN, weight_mean = NaN),
        provenance::NamedTuple = (; )
)
    metrics = summary.metrics
    config = summary.config
    graph_stats = summary.graph_stats
    score = tda_graph_score(metrics)
    generated_provenance = (
        backend = config.backend,
        max_dim = config.max_dim,
        max_points = config.max_points,
        thresholds = copy(config.thresholds),
        filtration = filtration,
        weight_name = weight_name
    )
    row_provenance = merge(provenance, generated_provenance)

    rows = Vector{NamedTuple}(undef, length(metrics.thresholds))
    for i in eachindex(metrics.thresholds)
        rows[i] = (
            graph_id = graph_id,
            threshold_index = i,
            threshold = metrics.thresholds[i],
            betti0 = metrics.betti0[i],
            betti1 = metrics.betti1[i],
            nv = graph_stats.nv,
            ne = graph_stats.ne,
            directed = graph_stats.directed,
            density = graph_stats.density,
            score = score,
            backend = config.backend,
            max_dim = config.max_dim,
            max_points = config.max_points,
            filtration = filtration,
            weight_name = weight_name,
            weight_min = weight_stats.weight_min,
            weight_max = weight_stats.weight_max,
            weight_mean = weight_stats.weight_mean,
            provenance = row_provenance
        )
    end

    return rows
end

"""
    tda_metric_table(summary::TDARunSummary; kwargs...) -> DataFrames.DataFrame

Return `tda_metric_rows(summary; kwargs...)` as a `DataFrames.DataFrame`.

# Arguments
- `summary`: `TDARunSummary` returned by `tda_on_graph`.
- `kwargs...`: Keyword arguments forwarded to `tda_metric_rows`.

# Returns
- `DataFrames.DataFrame` with one row per filtration threshold.

# Example
```julia
graph = Graphs.SimpleGraph(3)
summary = Mycelia.tda_on_graph(graph, Mycelia.TDAConfig(thresholds = [0.0]))
table = Mycelia.tda_metric_table(summary)
```
"""
function tda_metric_table(summary::TDARunSummary; kwargs...)
    return DataFrames.DataFrame(tda_metric_rows(summary; kwargs...))
end

"""
    extract_tda_metrics(g, cfg::TDAConfig = TDAConfig(); kwargs...) -> DataFrames.DataFrame

Compute graph Betti/TDA metrics and return one table-ready row per filtration
threshold. This is the convenience entry point for downstream assembly-quality
correlation analyses.

# Arguments
- `g`: `Graphs.AbstractGraph` or `MetaGraphsNext.MetaGraph` to summarize.
- `cfg = TDAConfig()`: TDA configuration.
- `vertex_weights = nothing`: Uniform weights, a vector, or a graph-label dictionary.
- `graph_id = "graph"`: Assembly graph identifier copied into every row.
- `filtration = :vertex_weight_threshold`: Filtration recipe name recorded in provenance.
- `weight_name`: Defaults to `:uniform` when `vertex_weights` is omitted, otherwise
  `:vertex_weight`.
- `provenance = (; )`: Extra run metadata merged into each row's provenance tuple.

# Returns
- `DataFrames.DataFrame` with stable TDA metric rows and provenance.

# Example
```julia
graph = Graphs.SimpleGraph(3)
table = Mycelia.extract_tda_metrics(
    graph,
    Mycelia.TDAConfig(thresholds = [0.0]);
    graph_id = "assembly_graph"
)
```
"""
function extract_tda_metrics(
        g,
        cfg::TDAConfig = TDAConfig();
        vertex_weights = nothing,
        graph_id = "graph",
        filtration::Symbol = :vertex_weight_threshold,
        weight_name::Symbol = isnothing(vertex_weights) ? :uniform : :vertex_weight,
        provenance::NamedTuple = (; )
)
    weights = _tda_vertex_weights(g, vertex_weights)
    summary = tda_on_graph(g, cfg; vertex_weights = weights)
    return tda_metric_table(
        summary;
        graph_id = graph_id,
        filtration = filtration,
        weight_name = weight_name,
        weight_stats = _tda_weight_stats(weights),
        provenance = provenance
    )
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
    fragmentation_penalty = β * (isempty(metrics.betti0) ? 0.0 :
                             Float64(maximum(metrics.betti0)))
    return cycle_penalty + fragmentation_penalty
end
