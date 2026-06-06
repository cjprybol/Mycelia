# Path Finding and Sequence Reconstruction Algorithms
#
# This file contains algorithms for finding paths through assembly graphs
# and reconstructing sequences from those paths.
#
# Key algorithms:
# - Eulerian path finding (Hierholzer's algorithm)
# - Sequence reconstruction from paths
# - Type-stable path traversal
# - K-shortest paths (Yen's algorithm)
#
# Based on legacy sequence-graph utilities now ported to Rhizomorph

# ============================================================================
# Path Data Structures
# ============================================================================

"""
    WalkStep{T}

Represents a single step in a graph walk.

# Fields
- `vertex_label::T`: Vertex visited at this step
- `strand::StrandOrientation`: Strand orientation used for traversal
- `probability::Float64`: Transition probability for this step
- `cumulative_probability::Float64`: Product of step probabilities up to this step

# Example
```julia
step = Mycelia.Rhizomorph.WalkStep("AAA", Mycelia.Rhizomorph.Forward, 0.5, 0.5)
```
"""
struct WalkStep{T}
    vertex_label::T
    strand::StrandOrientation
    probability::Float64
    cumulative_probability::Float64
end

"""
    GraphPath{T}

Represents a complete walk through a Rhizomorph graph.

# Fields
- `steps::Vector{WalkStep{T}}`: Ordered walk steps
- `total_probability::Float64`: Cumulative probability of the full walk

# Notes
`GraphPath` stores traversal metadata only. Use
`Mycelia.Rhizomorph.path_to_sequence` to reconstruct the assembled sequence.

# Example
```julia
steps = [Mycelia.Rhizomorph.WalkStep("AAA", Mycelia.Rhizomorph.Forward, 1.0, 1.0)]
path = Mycelia.Rhizomorph.GraphPath(steps)
```
"""
struct GraphPath{T}
    steps::Vector{WalkStep{T}}
    total_probability::Float64

    function GraphPath{T}(steps::Vector{WalkStep{T}}) where {T}
        total_prob = isempty(steps) ? 0.0 : last(steps).cumulative_probability
        new{T}(steps, total_prob)
    end
end

"""
    GraphPath(steps::Vector{WalkStep{T}}) where T

Construct a `GraphPath` and infer the vertex label type from `steps`.

# Arguments
- `steps::Vector{WalkStep{T}}`: Ordered traversal steps

# Returns
- `GraphPath{T}`: Path with `total_probability` derived from the final step

# Example
```julia
steps = [Mycelia.Rhizomorph.WalkStep("AAA", Mycelia.Rhizomorph.Forward, 1.0, 1.0)]
path = Mycelia.Rhizomorph.GraphPath(steps)
```
"""
function GraphPath(steps::Vector{WalkStep{T}}) where {T}
    return GraphPath{T}(steps)
end

"""
    LocalPathEnumerationProvenance{T}

Records how a local path enumeration was bounded and whether it was truncated.

# Fields
- `entry_vertex::T`: Local region entry vertex
- `exit_vertex::T`: Local region exit vertex
- `method::Symbol`: Enumeration method used
- `max_paths::Int`: Maximum alternatives requested
- `max_depth::Int`: Maximum edge depth explored per path
- `max_expansions::Int`: Maximum partial-path expansions allowed
- `expansions::Int`: Number of partial paths expanded
- `truncated::Bool`: Whether a guard stopped enumeration before exhaustion
- `reason::Symbol`: Completion reason (`:complete`, `:max_paths`,
  `:max_depth`, or `:max_expansions`)

# Example
```julia
result = Mycelia.Rhizomorph.enumerate_local_paths(weighted, source, target)
result.provenance.truncated
```
"""
struct LocalPathEnumerationProvenance{T}
    entry_vertex::T
    exit_vertex::T
    method::Symbol
    max_paths::Int
    max_depth::Int
    max_expansions::Int
    expansions::Int
    truncated::Bool
    reason::Symbol
end

"""
    RankedPathAlternative{T}

A ranked local path alternative with score and provenance metadata.

# Fields
- `rank::Int`: One-based rank after deterministic sorting
- `path::GraphPath{T}`: Enumerated graph path
- `score::Float64`: Path score; currently the path total probability
- `provenance::LocalPathEnumerationProvenance{T}`: Shared enumeration metadata

# Example
```julia
alternative = first(Mycelia.Rhizomorph.enumerate_local_paths(weighted, "A", "D").alternatives)
alternative.score == alternative.path.total_probability
```
"""
struct RankedPathAlternative{T}
    rank::Int
    path::GraphPath{T}
    score::Float64
    provenance::LocalPathEnumerationProvenance{T}
end

"""
    LocalPathEnumerationResult{T}

Result container for bounded local path enumeration.

# Fields
- `alternatives::Vector{RankedPathAlternative{T}}`: Ranked alternatives
- `provenance::LocalPathEnumerationProvenance{T}`: Search bounds and status

# Example
```julia
result = Mycelia.Rhizomorph.enumerate_local_paths(weighted, source, target; max_paths=5)
paths = [alternative.path for alternative in result.alternatives]
```
"""
struct LocalPathEnumerationResult{T}
    alternatives::Vector{RankedPathAlternative{T}}
    provenance::LocalPathEnumerationProvenance{T}
end

"""
    StrandWeightedEdgeData

Edge payload for probabilistic path algorithms with explicit strand tracking.

# Fields
- `weight::Float64`: Edge transition weight
- `src_strand::StrandOrientation`: Strand at the source vertex
- `dst_strand::StrandOrientation`: Strand at the destination vertex

# Example
```julia
edge = Mycelia.Rhizomorph.StrandWeightedEdgeData(
    3.0,
    Mycelia.Rhizomorph.Forward,
    Mycelia.Rhizomorph.Forward,
)
```
"""
struct StrandWeightedEdgeData
    weight::Float64
    src_strand::StrandOrientation
    dst_strand::StrandOrientation
end

"""
    edge_data_weight(edge_data)

Extract the numeric weight used for path scoring and random-walk sampling.

# Arguments
- `edge_data`: Rhizomorph edge payload or `StrandWeightedEdgeData`

# Returns
- `Float64`: Edge weight for traversal algorithms

# Example
```julia
edge = Mycelia.Rhizomorph.StrandWeightedEdgeData(
    2.5,
    Mycelia.Rhizomorph.Forward,
    Mycelia.Rhizomorph.Forward,
)
weight = Mycelia.Rhizomorph.edge_data_weight(edge)
```
"""
function edge_data_weight(edge_data::StrandWeightedEdgeData)
    return edge_data.weight
end

function edge_data_weight(edge_data)
    return compute_edge_weight(edge_data)
end

"""
    edge_quality_weight(edge_data)

Compute an edge weight from quality evidence when available.
Falls back to evidence counts for non-quality edges.

# Arguments
- `edge_data`: Rhizomorph edge payload with quality or evidence annotations

# Returns
- `Float64`: Quality-derived edge weight, or evidence count fallback

# Example
```julia
weight = Mycelia.Rhizomorph.edge_quality_weight(edge_data)
```
"""
function edge_quality_weight(edge_data)
    # TYPE-CHECK-AUDIT: super-union guard — must cover all reduced edge types
    if edge_data isa AllReducedEdgeData
        return Float64(count_evidence(edge_data))
    end
    entries = collect_evidence_entries(edge_data.evidence)
    if isempty(entries)
        return 0.0
    end

    total = 0.0
    quality_found = false

    for entry in entries
        if entry isa EdgeQualityEvidenceEntry
            quality_found = true
            combined = vcat(
                decode_quality_scores(entry.from_quality),
                decode_quality_scores(entry.to_quality)
            )
            if !isempty(combined)
                total += Statistics.mean(combined)
            end
        end
    end

    if quality_found && total > 0
        return total
    end

    return Float64(count_evidence(edge_data))
end

"""
    weighted_graph_from_rhizomorph(source_graph; default_weight=1e-10, edge_weight=count_evidence)

Convert a Rhizomorph evidence graph into a weighted graph suitable for
probabilistic path algorithms. Edge weights are derived from evidence counts
by default, and strand orientation is inferred from the first evidence entry
when present. Undirected graphs are expanded into bidirectional edges.

# Arguments
- `source_graph::MetaGraphsNext.MetaGraph`: Evidence-backed graph to convert
- `default_weight::Float64=1e-10`: Fallback weight for zero-support edges
- `edge_weight::Function=count_evidence`: Callback used to score each edge

# Returns
- `MetaGraphsNext.MetaGraph`: Directed graph with `StrandWeightedEdgeData` edges

# Example
```julia
weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(raw_graph)
```
"""
function weighted_graph_from_rhizomorph(
        source_graph::MetaGraphsNext.MetaGraph;
        default_weight::Float64 = 1e-10,
        edge_weight::Function = count_evidence
)
    labels = collect(MetaGraphsNext.labels(source_graph))
    label_type = isempty(labels) ? String : typeof(first(labels))

    weighted = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = label_type,
        vertex_data_type = Any,
        edge_data_type = StrandWeightedEdgeData,
        weight_function = edge_data_weight,
        default_weight = 0.0
    )

    for label in labels
        weighted[label] = nothing
    end

    is_directed = Graphs.is_directed(source_graph.graph)

    for (src, dst) in MetaGraphsNext.edge_labels(source_graph)
        edge_data = source_graph[src, dst]
        weight_value = Float64(edge_weight(edge_data))
        weight = weight_value > 0 ? weight_value : default_weight
        strand = if edge_data isa AllReducedEdgeData  # TYPE-CHECK-AUDIT: super-union guard
            Forward  # Reduced types don't track strand
        else
            first_evidence_strand(edge_data.evidence; default = Forward)
        end
        weighted[src, dst] = StrandWeightedEdgeData(weight, strand, strand)
        if !is_directed
            weighted[dst, src] = StrandWeightedEdgeData(weight, strand, strand)
        end
    end

    return weighted
end

"""
    probabilistic_walk_next(graph::MetaGraphsNext.MetaGraph, start_vertex, max_steps; seed=nothing)

Generate a weighted random walk from `start_vertex`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted Rhizomorph graph
- `start_vertex`: Vertex label to start from
- `max_steps::Int`: Maximum number of edge traversals
- `seed::Union{Nothing, Int}=nothing`: Optional RNG seed for reproducibility

# Returns
- `GraphPath`: Walk annotated with per-step and cumulative probabilities

# Example
```julia
weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(raw_graph)
path = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, first(MetaGraphsNext.labels(weighted)), 25)
```
"""
function probabilistic_walk_next(
        graph::MetaGraphsNext.MetaGraph,
        start_vertex::T,
        max_steps::Int;
        seed::Union{Nothing, Int} = nothing
) where {T}
    if seed !== nothing
        Mycelia.Random.seed!(seed)
    end

    if !(start_vertex in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Start vertex $start_vertex not found in graph"))
    end

    steps = WalkStep{T}[]
    current_vertex = start_vertex
    current_strand = Forward
    cumulative_prob = 1.0

    push!(steps, WalkStep(current_vertex, current_strand, 1.0, cumulative_prob))

    for _ in 1:max_steps
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)

        if isempty(valid_transitions)
            break
        end

        transition_probs = _calculate_transition_probabilities(valid_transitions)
        next_transition = _sample_transition(valid_transitions, transition_probs)

        step_prob = next_transition[:probability]
        cumulative_prob *= step_prob

        current_vertex = next_transition[:target_vertex]
        current_strand = next_transition[:target_strand]

        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end

    return GraphPath(steps)
end

"""
    maximum_weight_walk_next(graph::MetaGraphsNext.MetaGraph, start_vertex, max_steps; weight_function=edge_data_weight)

Generate a greedy walk that chooses the highest-weight outgoing transition at
each step.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted Rhizomorph graph
- `start_vertex`: Vertex label to start from
- `max_steps::Int`: Maximum number of edge traversals
- `weight_function::Function=edge_data_weight`: Scoring callback for transitions

# Returns
- `GraphPath`: Greedy walk through the graph

# Example
```julia
weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(raw_graph)
path = Mycelia.Rhizomorph.maximum_weight_walk_next(weighted, first(MetaGraphsNext.labels(weighted)), 25)
```
"""
function maximum_weight_walk_next(
        graph::MetaGraphsNext.MetaGraph,
        start_vertex::T,
        max_steps::Int;
        weight_function::Function = edge_data_weight
) where {T}
    if !(start_vertex in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Start vertex $start_vertex not found in graph"))
    end

    steps = WalkStep{T}[]
    current_vertex = start_vertex
    current_strand = Forward
    cumulative_prob = 1.0

    push!(steps, WalkStep(current_vertex, current_strand, 1.0, cumulative_prob))

    for _ in 1:max_steps
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)

        if isempty(valid_transitions)
            break
        end

        best_transition = nothing
        max_weight = -Inf

        for transition in valid_transitions
            weight = weight_function(transition[:edge_data])
            if weight > max_weight
                max_weight = weight
                best_transition = transition
            end
        end

        if best_transition === nothing
            break
        end

        step_prob = best_transition[:probability]
        cumulative_prob *= step_prob

        current_vertex = best_transition[:target_vertex]
        current_strand = best_transition[:target_strand]

        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end

    return GraphPath(steps)
end

function _normalize_strand(strand)
    if strand == Forward || strand == Reverse
        return strand
    end
    if isdefined(Mycelia, :Forward) && strand == Mycelia.Forward
        return Forward
    end
    if isdefined(Mycelia, :Reverse) && strand == Mycelia.Reverse
        return Reverse
    end
    return Forward
end

function _get_valid_transitions(graph, vertex_label, strand)
    transitions = []

    for edge_labels in MetaGraphsNext.edge_labels(graph)
        if length(edge_labels) == 2 && edge_labels[1] == vertex_label
            target_vertex = edge_labels[2]
            edge_data = graph[edge_labels...]

            edge_src_strand = _normalize_strand(edge_data.src_strand)
            if edge_src_strand == strand
                probability = edge_data.weight > 0 ? edge_data.weight : 1e-10

                push!(transitions,
                    Dict(
                        :target_vertex => target_vertex,
                        :target_strand => _normalize_strand(edge_data.dst_strand),
                        :probability => probability,
                        :edge_data => edge_data
                    ))
            end
        end
    end

    return transitions
end

function _calculate_transition_probabilities(transitions)
    if isempty(transitions)
        return Float64[]
    end

    weights = [t[:probability] for t in transitions]
    total_weight = sum(weights)

    if total_weight == 0
        return fill(1.0 / length(transitions), length(transitions))
    end

    return weights ./ total_weight
end

function _sample_transition(transitions, probabilities)
    if isempty(transitions)
        return nothing
    end

    if length(transitions) == 1
        return first(transitions)
    end

    r = Mycelia.Random.rand()
    cumulative = 0.0

    for (i, prob) in enumerate(probabilities)
        cumulative += prob
        if r <= cumulative
            return transitions[i]
        end
    end

    return last(transitions)
end

"""
    shortest_probability_path_next(graph::MetaGraphsNext.MetaGraph, source, target)

Find the highest-probability path between `source` and `target` by minimizing
negative log transition probabilities.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted Rhizomorph graph
- `source`: Source vertex label
- `target`: Target vertex label

# Returns
- `Union{Nothing, GraphPath}`: Most probable path, or `nothing` if unreachable

# Example
```julia
weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(raw_graph)
labels = collect(MetaGraphsNext.labels(weighted))
path = Mycelia.Rhizomorph.shortest_probability_path_next(weighted, labels[1], labels[end])
```
"""
function shortest_probability_path_next(
        graph::MetaGraphsNext.MetaGraph,
        source::T,
        target::T
) where {T}
    if !(source in MetaGraphsNext.labels(graph)) ||
       !(target in MetaGraphsNext.labels(graph))
        return nothing
    end

    distances = Dict{Tuple{T, StrandOrientation}, Float64}()
    predecessors = Dict{
        Tuple{T, StrandOrientation}, Union{Nothing, Tuple{T, StrandOrientation}}}()
    visited = Set{Tuple{T, StrandOrientation}}()
    pq = DataStructures.PriorityQueue{Tuple{T, StrandOrientation}, Float64}()

    for strand in (Forward, Reverse)
        start_state = (source, strand)
        distances[start_state] = 0.0
        predecessors[start_state] = nothing
        DataStructures.enqueue!(pq, start_state, 0.0)
    end

    while !isempty(pq)
        current_state = DataStructures.dequeue!(pq)
        if current_state in visited
            continue
        end

        push!(visited, current_state)
        current_vertex, current_strand = current_state

        if current_vertex == target
            return _reconstruct_shortest_path(
                predecessors, distances, (source, Forward), current_state, graph)
        end

        transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        for transition in transitions
            next_vertex = transition[:target_vertex]
            next_strand = transition[:target_strand]
            next_state = (next_vertex, next_strand)

            edge_prob = transition[:probability]
            distance = -log(edge_prob)
            new_distance = distances[current_state] + distance

            if !haskey(distances, next_state) || new_distance < distances[next_state]
                distances[next_state] = new_distance
                predecessors[next_state] = current_state
                DataStructures.enqueue!(pq, next_state, new_distance)
            end
        end
    end

    return nothing
end

function _reconstruct_shortest_path(predecessors, distances, start_state, end_state, graph)
    path_states = []
    current_state = end_state

    while current_state !== nothing
        pushfirst!(path_states, current_state)
        current_state = predecessors[current_state]
    end

    if isempty(path_states)
        return GraphPath(WalkStep{Any}[])
    end
    vertex_type = typeof(path_states[1][1])
    steps = WalkStep{vertex_type}[]

    for (i, (vertex, strand)) in enumerate(path_states)
        step_prob = if i == 1
            1.0
        else
            prev_state = path_states[i - 1]
            step_distance = distances[(vertex, strand)] - distances[prev_state]
            exp(-step_distance)
        end

        cumulative_prob = exp(-distances[(vertex, strand)])
        push!(steps, WalkStep(vertex, strand, step_prob, cumulative_prob))
    end

    return GraphPath(steps)
end

# ============================================================================
# K-Shortest Paths (Yen's Algorithm)
# ============================================================================

# Minimum weight floor to prevent -log(0) = Inf in Dijkstra distance calculations.
# This is a numerical stability guard, not a real probability.
const _KSP_MIN_WEIGHT = 1e-10

"""
    _total_outgoing_weight(graph, vertex)

Sum the original outgoing edge weights from `vertex`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted Rhizomorph graph
- `vertex`: Vertex label whose outgoing mass should be measured

# Returns
- `Float64`: Total outgoing edge weight from `vertex`, floored at `_KSP_MIN_WEIGHT`

# Example
```julia
total = Mycelia.Rhizomorph._total_outgoing_weight(weighted, source)
```
"""
function _total_outgoing_weight(graph::MetaGraphsNext.MetaGraph, vertex)
    total = 0.0
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        if src == vertex
            total += max(edge_data_weight(graph[src, dst]), _KSP_MIN_WEIGHT)
        end
    end
    return max(total, _KSP_MIN_WEIGHT)
end

function _total_outgoing_weight(
        graph::MetaGraphsNext.MetaGraph,
        vertex,
        strand::StrandOrientation
)
    total = 0.0
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        if src != vertex
            continue
        end
        edge_data = graph[src, dst]
        edge_src_strand = _normalize_strand(edge_data.src_strand)
        if edge_src_strand != strand
            continue
        end
        total += max(edge_data_weight(edge_data), _KSP_MIN_WEIGHT)
    end
    return max(total, _KSP_MIN_WEIGHT)
end

"""
    _total_outgoing_weight_excluding(graph, vertex, excluded_edges)

Sum the remaining outgoing edge weights from `vertex` after suppressing
`excluded_edges`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted Rhizomorph graph
- `vertex`: Vertex label whose outgoing mass should be measured
- `excluded_edges::Set{Tuple{T, T}}`: Edge label pairs `(src, dst)` to omit

# Returns
- `Float64`: Remaining outgoing edge weight from `vertex`, floored at `_KSP_MIN_WEIGHT`

# Example
```julia
excluded = Set([(source, blocked_neighbor)])
total = Mycelia.Rhizomorph._total_outgoing_weight_excluding(weighted, source, excluded)
```
"""
function _total_outgoing_weight_excluding(
        graph::MetaGraphsNext.MetaGraph,
        vertex,
        excluded_edges::Set{Tuple{T, T}}
) where {T}
    total = 0.0
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        if src == vertex && !((src, dst) in excluded_edges)
            w = edge_data_weight(graph[src, dst])
            total += max(w, _KSP_MIN_WEIGHT)
        end
    end
    return max(total, _KSP_MIN_WEIGHT)
end

"""
    _build_graph_path_from_vertices(graph, vertices)

Construct a `GraphPath` from an ordered list of vertex labels or `(vertex, strand)`
states, computing transition probabilities from the original graph weights.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted Rhizomorph graph
- `vertices`: Ordered path as vertex labels or `(vertex, strand)` states

# Returns
- `GraphPath`: Path annotated with per-step and cumulative probabilities

# Example
```julia
path = Mycelia.Rhizomorph._build_graph_path_from_vertices(weighted, [source, target])
```
"""
function _build_graph_path_from_vertices(
        graph::MetaGraphsNext.MetaGraph,
        vertices::Vector{T}
) where {T}
    path_states = Tuple{T, StrandOrientation}[(vertex, Forward) for vertex in vertices]
    return _build_graph_path_from_vertices(graph, path_states)
end

function _build_graph_path_from_vertices(
        graph::MetaGraphsNext.MetaGraph,
        path_states::Vector{Tuple{T, StrandOrientation}}
) where {T}
    if isempty(path_states)
        return GraphPath(WalkStep{T}[])
    end

    steps = WalkStep{T}[]
    cumulative_prob = 1.0

    for (i, (vertex, strand)) in enumerate(path_states)
        if i == 1
            push!(steps, WalkStep(vertex, strand, 1.0, 1.0))
        else
            prev_vertex, prev_strand = path_states[i - 1]
            if !haskey(graph, prev_vertex, vertex)
                throw(ArgumentError(
                    "Cannot build path: no edge from $(prev_vertex) to $(vertex) in graph"))
            end
            edge_data = graph[prev_vertex, vertex]
            edge_src_strand = _normalize_strand(edge_data.src_strand)
            edge_dst_strand = _normalize_strand(edge_data.dst_strand)
            if edge_src_strand != prev_strand || edge_dst_strand != strand
                throw(ArgumentError(
                    "Cannot build path: edge $(prev_vertex) => $(vertex) uses " *
                    "strand $(edge_src_strand) => $(edge_dst_strand), not " *
                    "$(prev_strand) => $(strand)"))
            end
            total_out = _total_outgoing_weight(graph, prev_vertex, prev_strand)
            edge_w = max(edge_data_weight(edge_data), _KSP_MIN_WEIGHT)
            step_prob = edge_w / total_out
            cumulative_prob *= step_prob
            push!(steps, WalkStep(vertex, strand, step_prob, cumulative_prob))
        end
    end

    return GraphPath(steps)
end

"""
    _shortest_path_excluding(graph, source, target, excluded_vertices, excluded_edges)

Dijkstra shortest path that skips excluded vertices and edges. Edge weights
are `-log(probability)` where probability = edge_weight / total_outgoing_weight.
Returns a `GraphPath` or `nothing` if no path exists.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted graph from `weighted_graph_from_rhizomorph`
- `source`: Source vertex label
- `target`: Target vertex label
- `excluded_vertices::Set`: Vertex labels to skip during traversal
- `excluded_edges::Set{Tuple}`: Edge label pairs `(src, dst)` to skip

# Returns
- `Union{Nothing, GraphPath}`: Most probable constrained path, or `nothing` if unreachable

# Example
```julia
path = Mycelia.Rhizomorph._shortest_path_excluding(
    weighted,
    source,
    target,
    Set([blocked_vertex]),
    Set([(source, blocked_neighbor)]),
)
```
"""
function _shortest_path_excluding(
        graph::MetaGraphsNext.MetaGraph,
        source::T,
        target::T,
        excluded_vertices::Set{T},
        excluded_edges::Set{Tuple{T, T}};
        initial_states::Union{Nothing, Vector{Tuple{T, StrandOrientation}}} = nothing
) where {T}
    if source in excluded_vertices || target in excluded_vertices
        return nothing
    end
    if !(source in MetaGraphsNext.labels(graph)) ||
       !(target in MetaGraphsNext.labels(graph))
        return nothing
    end

    start_states = if isnothing(initial_states)
        [(source, Forward), (source, Reverse)]
    else
        initial_states
    end

    distances = Dict{Tuple{T, StrandOrientation}, Float64}()
    predecessors = Dict{
        Tuple{T, StrandOrientation}, Union{Nothing, Tuple{T, StrandOrientation}}}()
    visited = Set{Tuple{T, StrandOrientation}}()
    pq = DataStructures.PriorityQueue{Tuple{T, StrandOrientation}, Float64}()

    for start_state in start_states
        distances[start_state] = 0.0
        predecessors[start_state] = nothing
        DataStructures.enqueue!(pq, start_state, 0.0)
    end

    while !isempty(pq)
        current_state = DataStructures.dequeue!(pq)

        if current_state in visited
            continue
        end
        push!(visited, current_state)

        current, current_strand = current_state

        if current == target
            path_states = Tuple{T, StrandOrientation}[]
            node = current_state
            while node !== nothing
                pushfirst!(path_states, node)
                node = predecessors[node]
            end
            return _build_graph_path_from_vertices(graph, path_states)
        end

        total_out = _total_outgoing_weight(graph, current, current_strand)

        for transition in _get_valid_transitions(graph, current, current_strand)
            src = current
            dst = transition[:target_vertex]
            if dst in excluded_vertices || (src, dst) in excluded_edges
                continue
            end

            next_state = (dst, transition[:target_strand])
            edge_w = max(edge_data_weight(transition[:edge_data]), _KSP_MIN_WEIGHT)
            transition_prob = edge_w / total_out
            if transition_prob <= 0.0
                continue
            end
            distance = -log(transition_prob)
            new_distance = distances[current_state] + distance

            if !haskey(distances, next_state) || new_distance < distances[next_state]
                distances[next_state] = new_distance
                predecessors[next_state] = current_state
                if haskey(pq, next_state)
                    pq[next_state] = new_distance
                else
                    DataStructures.enqueue!(pq, next_state, new_distance)
                end
            end
        end
    end

    return nothing
end

"""
    _paths_equal(a, b)

Check if two `GraphPath`s traverse the same sequence of vertex/strand states.

# Arguments
- `a::GraphPath`: First path to compare
- `b::GraphPath`: Second path to compare

# Returns
- `Bool`: `true` when both paths visit the same vertices with the same strands

# Example
```julia
same_path = Mycelia.Rhizomorph._paths_equal(path_a, path_b)
```
"""
function _paths_equal(a::GraphPath, b::GraphPath)
    if length(a.steps) != length(b.steps)
        return false
    end
    for (sa, sb) in zip(a.steps, b.steps)
        if sa.vertex_label != sb.vertex_label || sa.strand != sb.strand
            return false
        end
    end
    return true
end

"""
    k_shortest_paths(graph, source, target, k)

Find up to `k` most probable loopless paths between `source` and `target`
using Yen's algorithm. Internally, edge transition probabilities are converted
to distances via `-log(p)`, so the shortest distance path corresponds to the
most probable path. Returns paths sorted by **descending** total_probability
(most probable first).

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Weighted graph from `weighted_graph_from_rhizomorph`
- `source`: Source vertex label
- `target`: Target vertex label
- `k::Int`: Maximum number of paths to return (must be non-negative)

# Returns
- `Vector{GraphPath}`: Up to `k` paths, sorted by descending probability

# Throws
- `ArgumentError`: If `source` or `target` is not in the graph, or `k < 0`

# Notes
- When `source == target`, returns a single trivial single-vertex path with
  probability 1.0. Cycle detection (finding non-trivial paths from a vertex
  back to itself) is not currently supported.

# Example
```julia
graph = Mycelia.Rhizomorph.build_kmer_graph(records, k; mode = :singlestrand)
weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
paths = Mycelia.Rhizomorph.k_shortest_paths(weighted, source, target, 5)
for (i, path) in enumerate(paths)
    println("Path \$i: probability = \$(path.total_probability)")
end
```
"""
function k_shortest_paths(
        graph::MetaGraphsNext.MetaGraph,
        source::T,
        target::T,
        k::Int
) where {T}
    if !(source in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Source vertex $source not found in graph"))
    end
    if !(target in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Target vertex $target not found in graph"))
    end
    if k < 0
        throw(ArgumentError("k must be non-negative, got $k"))
    end
    if k == 0
        return GraphPath{T}[]
    end

    # Handle trivial source == target case
    if source == target
        step = WalkStep(source, Forward, 1.0, 1.0)
        return [GraphPath(WalkStep{T}[step])]
    end

    # A holds the confirmed k-shortest paths
    A = GraphPath{T}[]

    # Find the shortest path (A_1)
    empty_verts = Set{T}()
    empty_edges = Set{Tuple{T, T}}()
    a1 = _shortest_path_excluding(graph, source, target, empty_verts, empty_edges)
    if a1 === nothing
        return A
    end
    push!(A, a1)

    # B is the candidate list
    B = GraphPath{T}[]

    for ki in 2:k
        prev_path = A[ki - 1]
        prev_states = [(s.vertex_label, s.strand) for s in prev_path.steps]
        prev_vertices = [state[1] for state in prev_states]

        for i in 1:(length(prev_vertices) - 1)
            spur_node, spur_strand = prev_states[i]
            root_path_vertices = prev_vertices[1:i]
            root_path_states = prev_states[1:i]

            # Exclude edges that share the same root path prefix in already-found paths
            excluded_edges = Set{Tuple{T, T}}()
            for path in A
                path_verts = [s.vertex_label for s in path.steps]
                if length(path_verts) >= i + 1 && path_verts[1:i] == root_path_vertices
                    push!(excluded_edges, (path_verts[i], path_verts[i + 1]))
                end
            end

            # Exclude root path vertices except the spur node
            excluded_vertices = Set{T}()
            for j in 1:(i - 1)
                push!(excluded_vertices, root_path_vertices[j])
            end

            spur_path = _shortest_path_excluding(
                graph,
                spur_node,
                target,
                excluded_vertices,
                excluded_edges;
                initial_states = [(spur_node, spur_strand)]
            )

            if spur_path !== nothing
                spur_states = [(s.vertex_label, s.strand) for s in spur_path.steps]
                # Combine root + spur (root_path_vertices[1:i-1] + spur_path)
                combined_states = vcat(root_path_states[1:(i - 1)], spur_states)
                # Note: Path probabilities are computed against the full graph (all edges),
                # not the constrained graph used during spur path search. This is intentional:
                # the search uses exclusions to find alternative routes, but the final
                # probabilities reflect each path's likelihood in the original graph.
                candidate = _build_graph_path_from_vertices(graph, combined_states)

                # Only add if not a duplicate
                is_duplicate = false
                for existing in B
                    if _paths_equal(existing, candidate)
                        is_duplicate = true
                        break
                    end
                end
                for existing in A
                    if _paths_equal(existing, candidate)
                        is_duplicate = true
                        break
                    end
                end
                if !is_duplicate
                    push!(B, candidate)
                end
            end
        end

        if isempty(B)
            break
        end

        # Pick the candidate with highest probability (shortest -log distance)
        sort!(B; by = p -> p.total_probability, rev = true)
        push!(A, popfirst!(B))
    end

    # Sort results by descending probability
    sort!(A; by = p -> p.total_probability, rev = true)
    return A
end

# ============================================================================
# Bounded Local Path Enumeration
# ============================================================================

function _path_state_sort_key(states::Vector{Tuple{T, StrandOrientation}}) where {T}
    return join((string(vertex) * ":" * string(strand) for (vertex, strand) in states), ">")
end

function _graph_path_sort_key(path::GraphPath)
    return join(
        (string(step.vertex_label) * ":" * string(step.strand) for step in path.steps),
        ">"
    )
end

function _is_strand_weighted_graph(graph::MetaGraphsNext.MetaGraph)
    for edge_label in MetaGraphsNext.edge_labels(graph)
        return graph[edge_label...] isa StrandWeightedEdgeData
    end
    return true
end

function _path_enumeration_graph(
        graph::MetaGraphsNext.MetaGraph;
        default_weight::Float64,
        edge_weight::Function
)
    if _is_strand_weighted_graph(graph)
        return graph
    end
    return weighted_graph_from_rhizomorph(
        graph; default_weight = default_weight, edge_weight = edge_weight)
end

function _validate_local_path_enumeration_args(
        graph::MetaGraphsNext.MetaGraph,
        source::T,
        target::T,
        max_paths::Int,
        max_depth::Int,
        max_expansions::Int,
        min_probability::Float64
) where {T}
    labels = collect(MetaGraphsNext.labels(graph))
    if !(source in labels)
        throw(ArgumentError("Source vertex $source not found in graph"))
    end
    if !(target in labels)
        throw(ArgumentError("Target vertex $target not found in graph"))
    end
    if max_paths < 0
        throw(ArgumentError("max_paths must be non-negative, got $max_paths"))
    end
    if max_depth < 0
        throw(ArgumentError("max_depth must be non-negative, got $max_depth"))
    end
    if max_expansions < 0
        throw(ArgumentError("max_expansions must be non-negative, got $max_expansions"))
    end
    if min_probability < 0.0 || min_probability > 1.0
        throw(ArgumentError("min_probability must be in [0, 1], got $min_probability"))
    end
    return nothing
end

function _rank_local_path_alternatives(
        paths::Vector{GraphPath{T}},
        provenance::LocalPathEnumerationProvenance{T}
) where {T}
    sort!(paths; by = path -> (-path.total_probability, _graph_path_sort_key(path)))

    alternatives = RankedPathAlternative{T}[]
    for (rank, path) in enumerate(paths)
        push!(alternatives, RankedPathAlternative(rank, path, path.total_probability, provenance))
    end

    return alternatives
end

function _local_path_result(
        paths::Vector{GraphPath{T}},
        entry_vertex::T,
        exit_vertex::T,
        method::Symbol,
        max_paths::Int,
        max_depth::Int,
        max_expansions::Int,
        expansions::Int,
        truncated::Bool,
        reason::Symbol
) where {T}
    provenance = LocalPathEnumerationProvenance(
        entry_vertex,
        exit_vertex,
        method,
        max_paths,
        max_depth,
        max_expansions,
        expansions,
        truncated,
        reason
    )
    alternatives = _rank_local_path_alternatives(paths, provenance)
    return LocalPathEnumerationResult(alternatives, provenance)
end

function _enumerate_local_paths_impl(
        graph::MetaGraphsNext.MetaGraph,
        source::T,
        target::T,
        method::Symbol;
        max_paths::Int = 10,
        max_depth::Int = 100,
        max_expansions::Int = 10_000,
        min_probability::Float64 = 0.0,
        default_weight::Float64 = 1e-10,
        edge_weight::Function = count_evidence
) where {T}
    _validate_local_path_enumeration_args(
        graph, source, target, max_paths, max_depth, max_expansions, min_probability)

    search_graph = _path_enumeration_graph(
        graph; default_weight = default_weight, edge_weight = edge_weight)

    if max_paths == 0
        return _local_path_result(
            GraphPath{T}[],
            source,
            target,
            method,
            max_paths,
            max_depth,
            max_expansions,
            0,
            true,
            :max_paths
        )
    end

    if source == target
        path = _build_graph_path_from_vertices(search_graph, [(source, Forward)])
        return _local_path_result(
            GraphPath{T}[path],
            source,
            target,
            method,
            max_paths,
            max_depth,
            max_expansions,
            0,
            false,
            :complete
        )
    end

    paths = GraphPath{T}[]
    seen_paths = Set{String}()
    candidate_states = Dict{Int, Vector{Tuple{T, StrandOrientation}}}()
    candidate_distances = Dict{Int, Float64}()
    candidate_probabilities = Dict{Int, Float64}()
    queue = DataStructures.PriorityQueue{Int, Tuple{Float64, String}}()
    next_candidate_id = 0

    for strand in (Forward, Reverse)
        states = Tuple{T, StrandOrientation}[(source, strand)]
        next_candidate_id += 1
        candidate_states[next_candidate_id] = states
        candidate_distances[next_candidate_id] = 0.0
        candidate_probabilities[next_candidate_id] = 1.0
        DataStructures.enqueue!(queue, next_candidate_id, (0.0, _path_state_sort_key(states)))
    end

    expansions = 0
    truncated = false
    reason = :complete

    while !isempty(queue)
        candidate_id = DataStructures.dequeue!(queue)
        states = candidate_states[candidate_id]
        distance = candidate_distances[candidate_id]
        probability = candidate_probabilities[candidate_id]
        current_vertex, current_strand = last(states)

        delete!(candidate_states, candidate_id)
        delete!(candidate_distances, candidate_id)
        delete!(candidate_probabilities, candidate_id)

        if current_vertex == target
            path_key = _path_state_sort_key(states)
            if !(path_key in seen_paths) && probability >= min_probability
                path = _build_graph_path_from_vertices(search_graph, states)
                push!(paths, path)
                push!(seen_paths, path_key)
            end
            if length(paths) >= max_paths
                if !isempty(queue)
                    truncated = true
                    reason = :max_paths
                end
                break
            end
            continue
        end

        if length(states) - 1 >= max_depth
            if !isempty(_get_valid_transitions(search_graph, current_vertex, current_strand))
                truncated = true
                reason = :max_depth
            end
            continue
        end

        if expansions >= max_expansions
            truncated = true
            reason = :max_expansions
            break
        end
        expansions += 1

        transitions = _get_valid_transitions(search_graph, current_vertex, current_strand)
        sort!(transitions; by = transition -> (
            string(transition[:target_vertex]),
            string(transition[:target_strand])
        ))

        total_out = _total_outgoing_weight(search_graph, current_vertex, current_strand)

        for transition in transitions
            next_vertex = transition[:target_vertex]::T
            next_strand = transition[:target_strand]::StrandOrientation
            if any(state -> state[1] == next_vertex, states)
                continue
            end

            edge_w = max(edge_data_weight(transition[:edge_data]), _KSP_MIN_WEIGHT)
            step_probability = edge_w / total_out
            if step_probability <= 0.0
                continue
            end

            next_probability = probability * step_probability
            if next_probability < min_probability
                continue
            end

            next_states = copy(states)
            push!(next_states, (next_vertex, next_strand))
            next_distance = distance - log(step_probability)
            next_candidate_id += 1
            candidate_states[next_candidate_id] = next_states
            candidate_distances[next_candidate_id] = next_distance
            candidate_probabilities[next_candidate_id] = next_probability
            DataStructures.enqueue!(
                queue,
                next_candidate_id,
                (next_distance, _path_state_sort_key(next_states))
            )
        end
    end

    return _local_path_result(
        paths,
        source,
        target,
        method,
        max_paths,
        max_depth,
        max_expansions,
        expansions,
        truncated,
        reason
    )
end

"""
    enumerate_local_paths(graph, source, target; max_paths=10, max_depth=100, max_expansions=10000, min_probability=0.0)

Enumerate ranked loopless alternatives between two local graph vertices with
explicit bounds to prevent global path explosion.

The search uses a best-first queue ordered by negative log transition
probability, then deterministically breaks equal-score ties by the vertex/strand
path key. Source graphs with evidence-backed Rhizomorph edge data are converted
with `weighted_graph_from_rhizomorph`; already weighted graphs with
`StrandWeightedEdgeData` are used directly.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Rhizomorph or weighted graph
- `source`: Entry vertex label
- `target`: Exit vertex label

# Keywords
- `max_paths::Int=10`: Maximum number of alternatives to return
- `max_depth::Int=100`: Maximum number of edges per path
- `max_expansions::Int=10000`: Maximum partial paths to expand
- `min_probability::Float64=0.0`: Prune paths below this cumulative probability
- `default_weight::Float64=1e-10`: Fallback edge weight during conversion
- `edge_weight::Function=count_evidence`: Edge scoring callback during conversion

# Returns
- `LocalPathEnumerationResult`: Ranked alternatives plus provenance metadata

# Notes
- Enumerated paths are vertex-simple. Non-trivial cycles from a vertex back to
  itself are intentionally not enumerated.
- Highly repetitive regions should be bounded with `max_paths`, `max_depth`, and
  `max_expansions`. If a guard stops the search, `result.provenance.truncated`
  is `true` and `result.provenance.reason` identifies the guard.

# Example
```julia
result = Mycelia.Rhizomorph.enumerate_local_paths(
    weighted_graph,
    "A",
    "D";
    max_paths = 5,
    max_depth = 20,
)
scores = [alternative.score for alternative in result.alternatives]
```
"""
function enumerate_local_paths(
        graph::MetaGraphsNext.MetaGraph,
        source::T,
        target::T;
        max_paths::Int = 10,
        max_depth::Int = 100,
        max_expansions::Int = 10_000,
        min_probability::Float64 = 0.0,
        default_weight::Float64 = 1e-10,
        edge_weight::Function = count_evidence
) where {T}
    return _enumerate_local_paths_impl(
        graph,
        source,
        target,
        :bounded_best_first;
        max_paths = max_paths,
        max_depth = max_depth,
        max_expansions = max_expansions,
        min_probability = min_probability,
        default_weight = default_weight,
        edge_weight = edge_weight
    )
end

# ============================================================================
# Eulerian Path Finding
# ============================================================================

"""
    find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, <:Any, <:Any})

Efficient Eulerian path finder that works directly with the underlying Graphs.jl structure.
No temporary graph creation - works with vertex indices and translates back to labels.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: The graph to find paths in

# Returns
- `Vector{Vector{T}}`: Vector containing the Eulerian path (or empty if no path exists)

# Details
Uses Hierholzer's algorithm to find an Eulerian path. An Eulerian path visits every edge
exactly once. The algorithm checks that:
- At most one vertex has (out-degree - in-degree) = 1 (start vertex)
- At most one vertex has (out-degree - in-degree) = -1 (end vertex)
- All other vertices have equal in-degree and out-degree

# Example
```julia
paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
```
"""
function find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{
        <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any}) where {T}
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return Vector{Vector{T}}()
    end

    adj_list = Dict{T, Vector{T}}()
    in_degrees = Dict{T, Int}()
    out_degrees = Dict{T, Int}()

    for label in labels
        adj_list[label] = T[]
        in_degrees[label] = 0
        out_degrees[label] = 0
    end

    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        # Treat each edge once for Eulerian feasibility; evidence counts are coverage, not topology.
        out_degrees[src] += 1
        in_degrees[dst] += 1
        push!(adj_list[src], dst)
    end

    # Check Eulerian path conditions
    start_vertices = T[]
    end_vertices = T[]

    for label in labels
        diff = out_degrees[label] - in_degrees[label]
        if diff == 1
            push!(start_vertices, label)
        elseif diff == -1
            push!(end_vertices, label)
        elseif diff != 0
            return Vector{Vector{T}}()
        end
    end

    if length(start_vertices) > 1 || length(end_vertices) > 1
        return Vector{Vector{T}}()
    end
    if length(start_vertices) != length(end_vertices)
        return Vector{Vector{T}}()
    end

    start_vertex = if !isempty(start_vertices)
        start_vertices[1]
    else
        first(labels)
    end

    path = _find_eulerian_path_hierholzer_labels(adj_list, start_vertex)
    if isempty(path)
        return Vector{Vector{T}}()
    end

    return [path]
end

"""
Hierholzer's algorithm implementation working directly on Graphs.jl structure.

# Arguments
- `graph::Graphs.AbstractGraph`: The underlying graph
- `start_vertex::Int`: Starting vertex index

# Returns
- `Vector{Int}`: Path as vertex indices (empty if no Eulerian path exists)

# Algorithm
Hierholzer's algorithm for finding Eulerian paths:
1. Start from the designated starting vertex
2. Follow edges, marking them as used
3. When stuck (no unused outgoing edges), backtrack
4. Insert subcircuits when encountering vertices with unused edges
5. Continue until all edges are used
"""
function _find_eulerian_path_hierholzer(graph::Graphs.AbstractGraph, start_vertex::Int)
    # Create mutable adjacency lists (using Sets for easier removal)
    adj_list = Dict{Int, Set{Int}}()

    for v in Graphs.vertices(graph)
        adj_list[v] = Set(Graphs.outneighbors(graph, v))
    end

    # Stack-based Hierholzer's algorithm
    path = Int[]
    stack = [start_vertex]

    while !isempty(stack)
        v = stack[end]

        if !isempty(adj_list[v])
            # Take an edge (any edge) from v
            next_v = pop!(adj_list[v])
            push!(stack, next_v)
        else
            # No more edges from v - add to path
            push!(path, pop!(stack))
        end
    end

    # Path is built in reverse order
    reverse!(path)

    # Verify all edges were used
    total_edges_used = Graphs.ne(graph) -
                       sum(length(neighbors) for neighbors in values(adj_list))

    if total_edges_used != Graphs.ne(graph)
        return Int[]  # Failed to find Eulerian path
    end

    return path
end

function _find_eulerian_path_hierholzer_labels(adj_list::Dict{T, Vector{T}}, start_vertex::T) where {T}
    path = T[]
    stack = T[start_vertex]

    while !isempty(stack)
        v = stack[end]
        if !isempty(adj_list[v])
            next_v = pop!(adj_list[v])
            push!(stack, next_v)
        else
            push!(path, pop!(stack))
        end
    end

    reverse!(path)

    if any(values(adj_list)) do neighbors
        !isempty(neighbors)
    end
        return T[]
    end

    return path
end

"""
    _find_eulerian_path_dfs(graph, start_vertex) -> Vector{Int}

Alternative Eulerian path finder using depth-first search with backtracking.

This is a simpler alternative to Hierholzer's algorithm that may handle
certain graph structures more robustly. Uses recursive DFS with edge removal.

# Algorithm
1. Start from the given vertex
2. Recursively follow edges, removing them as we go
3. Add vertices to path in post-order (after exploring all edges)
4. Reverse the final path
"""
function _find_eulerian_path_dfs(graph::Graphs.AbstractGraph, start_vertex::Int)
    # Create mutable adjacency lists
    adj_list = Dict{Int, Vector{Int}}()

    for v in Graphs.vertices(graph)
        adj_list[v] = collect(Graphs.outneighbors(graph, v))
    end

    # DFS with edge removal
    path = Int[]

    function dfs(v::Int)
        while !isempty(adj_list[v])
            next_v = pop!(adj_list[v])
            dfs(next_v)
        end
        push!(path, v)
    end

    dfs(start_vertex)
    reverse!(path)

    # Verify all edges used
    total_edges_remaining = sum(length(neighbors) for neighbors in values(adj_list))

    if total_edges_remaining > 0
        return Int[]  # Failed - not all edges used
    end

    return path
end

# ============================================================================
# Sequence Reconstruction from Paths
# ============================================================================

"""
    path_to_sequence(path::Vector{T}, graph::MetaGraphsNext.MetaGraph) where T

Convert a vector of k-mer labels (from find_eulerian_paths_next) to a biological sequence.

This is an overload that accepts raw k-mer vectors directly, which is the output format
of find_eulerian_paths_next. For k-mer graphs, reconstructs by overlapping k-mers.

# Arguments
- `path::Vector{T}`: Path as vector of k-mer labels
- `graph::MetaGraphsNext.MetaGraph`: Assembly graph

# Returns
Appropriate sequence type based on the k-mer type.

# Example
```julia
paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
sequence = Mycelia.Rhizomorph.path_to_sequence(first(paths), graph)
```
"""
function path_to_sequence(path::Vector{T}, graph::MetaGraphsNext.MetaGraph) where {T}
    if isempty(path)
        return _get_empty_sequence_for_label_type(T)
    end

    # Determine the sequence type from the first k-mer
    SequenceType = _sequence_type_from_kmer_type(T)

    # For k-mer paths: first k-mer in full, then add last base of each subsequent k-mer
    if SequenceType <: BioSequences.LongSequence
        # Start with the first k-mer
        result = SequenceType(string(path[1]))

        # Add last base of each subsequent k-mer
        for i in 2:length(path)
            kmer_str = string(path[i])
            last_base = kmer_str[end:end]
            result = result * SequenceType(last_base)
        end

        return result
    else
        # String type: concatenate with overlaps
        result = string(path[1])
        for i in 2:length(path)
            result *= string(path[i])[end:end]
        end
        return result
    end
end

"""
    path_to_sequence(path::GraphPath{T}, graph::MetaGraphsNext.MetaGraph) where T

Convert a graph path to a biological sequence.

# Arguments
- `path::GraphPath{T}`: Path through the assembly graph
- `graph::MetaGraphsNext.MetaGraph`: Assembly graph

# Returns
Appropriate sequence type based on the graph content (maintains type stability).
For k-mer graphs, reconstructs by overlapping k-mers. For variable-length graphs,
concatenates sequences.

# Details
- For k-mer graphs: Takes first k-mer in full, then adds last base of each subsequent k-mer
- For BioSequence graphs: Handles overlap based on edge data
- Maintains type stability - returns same sequence type as graph contains
- Handles strand orientation correctly (forward/reverse complement)

# Example
```julia
weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(raw_graph)
start_vertex = first(MetaGraphsNext.labels(weighted))
path = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start_vertex, 25)
sequence = Mycelia.Rhizomorph.path_to_sequence(path, weighted)
```
"""
function path_to_sequence(path::GraphPath{T}, graph::MetaGraphsNext.MetaGraph) where {T}
    if isempty(path.steps)
        # Return appropriate empty sequence type
        return _get_empty_sequence_for_label_type(T)
    end

    # Determine the sequence type from the first vertex
    first_step = path.steps[1]
    SequenceType = _get_sequence_type_for_path(first_step, graph)

    # Collect sequence parts maintaining proper types
    sequence_parts = Vector{SequenceType}()

    for (i, step) in enumerate(path.steps)
        if i == 1
            # First k-mer: add the full sequence
            seq_part = _extract_sequence_from_step(step, graph, SequenceType)
            push!(sequence_parts, seq_part)
        else
            # Subsequent k-mers: add only the last symbol (overlap by k-1)
            seq_part = _extract_sequence_from_step(step, graph, SequenceType)
            if length(seq_part) > 0
                last_symbol = seq_part[end:end]
                push!(sequence_parts, last_symbol)
            end
        end
    end

    # Concatenate all parts
    if isempty(sequence_parts)
        return _get_empty_sequence_for_label_type(T)
    end

    result = sequence_parts[1]
    for part in sequence_parts[2:end]
        result = result * part
    end

    return result
end

"""
Extract sequence from a single step, maintaining type stability.
"""
function _extract_sequence_from_step(step::WalkStep{T}, graph, SequenceType) where {T}
    vertex_data = graph[step.vertex_label]

    # Handle different vertex data types
    if hasfield(typeof(vertex_data), :sequence)
        # BioSequence graph: vertex data contains the sequence directly
        return vertex_data.sequence
    elseif hasfield(typeof(vertex_data), :Kmer)
        # K-mer graph: extract k-mer and convert to appropriate sequence type
        kmer = vertex_data.Kmer
        return _convert_kmer_to_sequence(kmer, SequenceType, step.strand)
    elseif hasfield(typeof(vertex_data), :canonical_kmer)
        # Legacy K-mer graph: extract k-mer and convert to appropriate sequence type
        kmer = vertex_data.canonical_kmer
        return _convert_kmer_to_sequence(kmer, SequenceType, step.strand)
    else
        # Direct k-mer labels: vertex label IS the k-mer
        return _convert_kmer_to_sequence(step.vertex_label, SequenceType, step.strand)
    end
end

"""
Convert k-mer to appropriate sequence type, handling strand orientation.
"""
function _convert_kmer_to_sequence(kmer, SequenceType, strand::StrandOrientation)
    # Convert k-mer to sequence
    if SequenceType <: BioSequences.LongSequence
        seq = SequenceType(string(kmer))
        # Handle reverse strand for biological sequences
        if strand == Reverse
            try
                seq = BioSequences.reverse_complement(seq)
            catch
                # If reverse complement fails (e.g., for AA), keep original
            end
        end
        return seq
    else
        # For string types or others, convert to string
        kmer_str = string(kmer)
        return strand == Forward ? kmer_str : kmer_str  # No rev comp for strings
    end
end

"""
Determine the appropriate sequence type for reconstruction from the first step.
"""
function _get_sequence_type_for_path(step::WalkStep{T}, graph) where {T}
    vertex_data = graph[step.vertex_label]

    if hasfield(typeof(vertex_data), :sequence)
        # BioSequence graph: use the sequence type from vertex data
        return typeof(vertex_data.sequence)
    elseif hasfield(typeof(vertex_data), :Kmer)
        # K-mer graph: determine sequence type from k-mer
        kmer = vertex_data.Kmer
        return _sequence_type_from_kmer_type(typeof(kmer))
    elseif hasfield(typeof(vertex_data), :canonical_kmer)
        # Legacy K-mer graph: determine sequence type from k-mer
        kmer = vertex_data.canonical_kmer
        return _sequence_type_from_kmer_type(typeof(kmer))
    else
        # Direct k-mer labels
        return _sequence_type_from_kmer_type(T)
    end
end

"""
Map k-mer type to corresponding sequence type.
"""
function _sequence_type_from_kmer_type(kmer_type::Type)
    if kmer_type <: AbstractString
        return String
    end

    # For Kmers.jl types
    if kmer_type <: Kmers.DNAKmer
        return BioSequences.LongDNA{4}
    elseif kmer_type <: Kmers.RNAKmer
        return BioSequences.LongRNA{4}
    elseif kmer_type <: Kmers.AAKmer
        return BioSequences.LongAA
    end

    # Fallback: try to infer from type parameters
    if length(kmer_type.parameters) >= 1
        alphabet_type = kmer_type.parameters[1]
        if alphabet_type <: BioSequences.DNAAlphabet
            return BioSequences.LongDNA{4}
        elseif alphabet_type <: BioSequences.RNAAlphabet
            return BioSequences.LongRNA{4}
        elseif alphabet_type <: BioSequences.AminoAcidAlphabet
            return BioSequences.LongAA
        end
    end

    # Default fallback
    return String
end

"""
Get empty sequence for a given label type.
"""
function _get_empty_sequence_for_label_type(::Type{T}) where {T}
    if T <: Kmers.DNAKmer
        return BioSequences.LongDNA{4}()
    elseif T <: Kmers.RNAKmer
        return BioSequences.LongRNA{4}()
    elseif T <: Kmers.AAKmer
        return BioSequences.LongAA()
    elseif T <: BioSequences.LongDNA
        return BioSequences.LongDNA{4}()
    elseif T <: BioSequences.LongRNA
        return BioSequences.LongRNA{4}()
    elseif T <: BioSequences.LongAA
        return BioSequences.LongAA()
    else
        return ""
    end
end
