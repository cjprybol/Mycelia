# Batch sequence generation via weighted random walks.
#
# Builds on existing probabilistic_walk_next and path_to_sequence from
# path-finding.jl to provide a high-level generation API with temperature
# control, start vertex selection strategies, and sequence likelihood scoring.

# ============================================================================
# Batch Sequence Generation
# ============================================================================

"""
    generate_sequences(graph, n; walk_length=100, seed=nothing, temperature=1.0, start_vertices=nothing)

Generate `n` sequences via weighted random walks on a Rhizomorph graph.

The graph should be a weighted graph (output of `weighted_graph_from_rhizomorph`).
If a non-weighted Rhizomorph evidence graph is passed, it will be converted automatically.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to walk on
- `n::Int`: Number of sequences to generate
- `walk_length::Int=100`: Maximum steps per walk
- `seed::Union{Nothing, Int}=nothing`: Random seed for reproducibility
- `temperature::Float64=1.0`: Sampling temperature (>1 = more random, <1 = more greedy)
- `start_vertices::Union{Nothing, Vector}=nothing`: Specific start vertices (cycles through if fewer than n)

# Returns
- `Vector{NamedTuple{(:sequence, :path, :walk_probability, :length)}}` where:
  - `sequence`: The generated sequence (String or BioSequence)
  - `path`: The GraphPath object
  - `walk_probability`: Total probability of the walk
  - `length`: Length of the generated sequence
"""
function generate_sequences(
        graph::MetaGraphsNext.MetaGraph,
        n::Int;
        walk_length::Int = 100,
        seed::Union{Nothing, Int} = nothing,
        temperature::Float64 = 1.0,
        start_vertices::Union{Nothing, Vector} = nothing
)
    if seed !== nothing
        Mycelia.Random.seed!(seed)
    end

    working_graph = _ensure_weighted_graph(graph)

    vertices = if start_vertices === nothing
        select_start_vertices(working_graph, n)
    else
        start_vertices
    end

    results = NamedTuple{(:sequence, :path, :walk_probability, :length)}[]

    for i in 1:n
        start_vertex = vertices[mod1(i, length(vertices))]

        path = if temperature == 1.0
            probabilistic_walk_next(working_graph, start_vertex, walk_length)
        else
            _tempered_walk(working_graph, start_vertex, walk_length, temperature)
        end

        try
            seq = path_to_sequence(path, working_graph)
            push!(results,
                (
                    sequence = seq,
                    path = path,
                    walk_probability = path.total_probability,
                    length = length(seq)
                ))
        catch e
            @warn "Failed to convert path to sequence at iteration $i" exception = e
            continue
        end
    end

    return results
end

"""
    _ensure_weighted_graph(graph::MetaGraphsNext.MetaGraph)

Check whether `graph` already has `StrandWeightedEdgeData` edges. If not,
convert it via `weighted_graph_from_rhizomorph`.
"""
function _ensure_weighted_graph(graph::MetaGraphsNext.MetaGraph)
    edge_iter = MetaGraphsNext.edge_labels(graph)
    first_edge = iterate(edge_iter)

    if first_edge === nothing
        return graph
    end

    (src, dst), _ = first_edge
    edge_data = graph[src, dst]

    if edge_data isa StrandWeightedEdgeData
        return graph
    end

    return weighted_graph_from_rhizomorph(graph)
end

# ============================================================================
# Temperature-Scaled Random Walk
# ============================================================================

"""
    _tempered_walk(graph, start_vertex, max_steps, temperature)

Temperature-scaled random walk. Reuses internal helpers from path-finding.jl.

Temperature > 1.0 increases randomness; < 1.0 makes walks more greedy.
Temperature approaching 0 approximates the maximum-weight walk.
"""
function _tempered_walk(
        graph::MetaGraphsNext.MetaGraph,
        start_vertex::T,
        max_steps::Int,
        temperature::Float64
) where {T}
    if temperature <= 0.0
        throw(ArgumentError("Temperature must be positive, got $temperature"))
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

        base_probs = _calculate_transition_probabilities(valid_transitions)
        tempered = base_probs .^ (1.0 / temperature)
        tempered_sum = sum(tempered)

        if tempered_sum > 0.0
            tempered ./= tempered_sum
        else
            tempered = fill(1.0 / length(valid_transitions), length(valid_transitions))
        end

        next_transition = _sample_transition(valid_transitions, tempered)

        if next_transition === nothing
            break
        end

        step_prob = next_transition[:probability]
        cumulative_prob *= step_prob

        current_vertex = next_transition[:target_vertex]
        current_strand = next_transition[:target_strand]

        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end

    return GraphPath(steps)
end

# ============================================================================
# Start Vertex Selection
# ============================================================================

"""
    select_start_vertices(graph, n; strategy=:degree_weighted)

Select `n` start vertices for generation.

# Strategies
- `:degree_weighted`: Probability proportional to out-degree (default)
- `:uniform`: Uniform random selection with replacement
- `:largest_component`: Restrict to largest connected component, then degree-weighted
"""
function select_start_vertices(
        graph::MetaGraphsNext.MetaGraph,
        n::Int;
        strategy::Symbol = :degree_weighted
)
    labels = collect(MetaGraphsNext.labels(graph))

    if isempty(labels)
        return eltype(labels)[]
    end

    if strategy == :uniform
        return [labels[Mycelia.Random.rand(1:length(labels))] for _ in 1:n]
    elseif strategy == :degree_weighted
        return _select_degree_weighted(graph, labels, n)
    elseif strategy == :largest_component
        return _select_largest_component(graph, labels, n)
    else
        throw(ArgumentError("Unknown strategy: $strategy"))
    end
end

function _select_degree_weighted(graph::MetaGraphsNext.MetaGraph, labels::Vector, n::Int)
    degrees = Float64[]
    for label in labels
        transitions = _get_valid_transitions(graph, label, Forward)
        push!(degrees, Float64(length(transitions)) + 1.0)
    end

    return _weighted_sample_with_replacement(labels, degrees, n)
end

function _select_largest_component(graph::MetaGraphsNext.MetaGraph, labels::Vector, n::Int)
    components = Graphs.connected_components(graph.graph)

    if isempty(components)
        return eltype(labels)[]
    end

    largest_component = components[argmax(length.(components))]
    component_labels = [MetaGraphsNext.label_for(graph, code) for code in largest_component]

    return _select_degree_weighted(graph, component_labels, n)
end

function _weighted_sample_with_replacement(items::Vector, weights::Vector{Float64}, n::Int)
    total_weight = sum(weights)
    if total_weight == 0.0
        return [items[Mycelia.Random.rand(1:length(items))] for _ in 1:n]
    end

    cumulative = cumsum(weights ./ total_weight)

    selected = Vector{eltype(items)}(undef, n)
    for i in 1:n
        r = Mycelia.Random.rand()
        idx = searchsortedfirst(cumulative, r)
        idx = clamp(idx, 1, length(items))
        selected[i] = items[idx]
    end

    return selected
end

# ============================================================================
# Sequence Likelihood Scoring
# ============================================================================

"""
    compute_sequence_likelihood(sequence::AbstractString, graph; k::Int)

Compute log-likelihood of a sequence under the graph's transition model.

Decomposes the sequence into k-mers and sums log₂ transition probabilities.
Returns `-Inf` if any k-mer transition is not found in the graph.
Returns `0.0` for sequences shorter than or equal to `k`.
"""
function compute_sequence_likelihood(
        sequence::AbstractString,
        graph::MetaGraphsNext.MetaGraph;
        k::Int
)
    seq_len = length(sequence)
    if seq_len <= k
        return 0.0
    end

    working_graph = _ensure_weighted_graph(graph)
    kmers = [sequence[i:(i + k - 1)] for i in 1:(seq_len - k + 1)]

    log_likelihood = 0.0

    for i in 1:(length(kmers) - 1)
        src_kmer = kmers[i]
        dst_kmer = kmers[i + 1]

        if !haskey(working_graph, src_kmer, dst_kmer)
            return -Inf
        end

        transition_prob = _compute_edge_transition_probability(
            working_graph, src_kmer, dst_kmer)

        if transition_prob <= 0.0
            return -Inf
        end

        log_likelihood += log2(transition_prob)
    end

    return log_likelihood
end

function _compute_edge_transition_probability(
        graph::MetaGraphsNext.MetaGraph, src, dst)
    if !haskey(graph, src, dst)
        return 0.0
    end

    edge_data = graph[src, dst]
    edge_weight = edge_data isa StrandWeightedEdgeData ? edge_data.weight : 1.0

    total_outgoing_weight = 0.0
    for (edge_src, edge_dst) in MetaGraphsNext.edge_labels(graph)
        if edge_src == src
            out_data = graph[edge_src, edge_dst]
            w = out_data isa StrandWeightedEdgeData ? out_data.weight : 1.0
            total_outgoing_weight += w
        end
    end

    if total_outgoing_weight == 0.0
        return 0.0
    end

    return edge_weight / total_outgoing_weight
end
