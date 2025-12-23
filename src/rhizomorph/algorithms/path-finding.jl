# Path Finding and Sequence Reconstruction Algorithms
#
# This file contains algorithms for finding paths through assembly graphs
# and reconstructing sequences from those paths.
#
# Key algorithms:
# - Eulerian path finding (Hierholzer's algorithm)
# - Sequence reconstruction from paths
# - Type-stable path traversal
#
# Based on functions from src/sequence-graphs-next.jl

# ============================================================================
# Path Data Structures
# ============================================================================

"""
Represents a step in a probabilistic walk through the graph.
Type parameter T allows for different vertex label types (k-mers, strings, biosequences, etc.)
"""
struct WalkStep{T}
    vertex_label::T
    strand::StrandOrientation
    probability::Float64
    cumulative_probability::Float64
end

"""
Represents a complete path through the graph.
Type parameter T allows for different vertex label types (k-mers, strings, biosequences, etc.)
No pre-computed sequence - use path_to_sequence() to reconstruct sequences from paths.
"""
struct GraphPath{T}
    steps::Vector{WalkStep{T}}
    total_probability::Float64

    function GraphPath{T}(steps::Vector{WalkStep{T}}) where T
        total_prob = isempty(steps) ? 0.0 : last(steps).cumulative_probability
        new{T}(steps, total_prob)
    end
end

# Convenience constructor that infers the type parameter
function GraphPath(steps::Vector{WalkStep{T}}) where T
    return GraphPath{T}(steps)
end

"""
Edge data for probabilistic path algorithms with strand awareness.
"""
struct StrandWeightedEdgeData
    weight::Float64
    src_strand::StrandOrientation
    dst_strand::StrandOrientation
end

"""
    edge_quality_weight(edge_data)

Compute an edge weight from quality evidence when available.
Falls back to evidence counts for non-quality edges.
"""
function edge_quality_weight(edge_data)
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
"""
function weighted_graph_from_rhizomorph(
    source_graph::MetaGraphsNext.MetaGraph;
    default_weight::Float64=1e-10,
    edge_weight::Function=count_evidence,
)
    labels = collect(MetaGraphsNext.labels(source_graph))
    label_type = isempty(labels) ? String : typeof(first(labels))

    weighted = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=label_type,
        vertex_data_type=Any,
        edge_data_type=StrandWeightedEdgeData,
        weight_function=Mycelia.edge_data_weight,
        default_weight=0.0,
    )

    for label in labels
        weighted[label] = nothing
    end

    is_directed = Graphs.is_directed(source_graph.graph)

    for (src, dst) in MetaGraphsNext.edge_labels(source_graph)
        edge_data = source_graph[src, dst]
        weight_value = Float64(edge_weight(edge_data))
        weight = weight_value > 0 ? weight_value : default_weight
        strand = first_evidence_strand(edge_data.evidence; default=Forward)
        weighted[src, dst] = StrandWeightedEdgeData(weight, strand, strand)
        if !is_directed
            weighted[dst, src] = StrandWeightedEdgeData(weight, strand, strand)
        end
    end

    return weighted
end

function probabilistic_walk_next(
    graph::MetaGraphsNext.MetaGraph,
    start_vertex::T,
    max_steps::Int;
    seed::Union{Nothing, Int}=nothing,
) where T
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

function maximum_weight_walk_next(
    graph::MetaGraphsNext.MetaGraph,
    start_vertex::T,
    max_steps::Int;
    weight_function::Function = Mycelia.edge_data_weight,
) where T
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

                push!(transitions, Dict(
                    :target_vertex => target_vertex,
                    :target_strand => _normalize_strand(edge_data.dst_strand),
                    :probability => probability,
                    :edge_data => edge_data,
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

function shortest_probability_path_next(
    graph::MetaGraphsNext.MetaGraph,
    source::T,
    target::T,
) where T
    if !(source in MetaGraphsNext.labels(graph)) || !(target in MetaGraphsNext.labels(graph))
        return nothing
    end

    distances = Dict{Tuple{T, StrandOrientation}, Float64}()
    predecessors = Dict{Tuple{T, StrandOrientation}, Union{Nothing, Tuple{T, StrandOrientation}}}()
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
            return _reconstruct_shortest_path(predecessors, distances, (source, Forward), current_state, graph)
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
"""
function find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any}) where T
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
        edge_data = graph[src, dst]
        multiplicity = max(1, count_evidence(edge_data))
        out_degrees[src] += multiplicity
        in_degrees[dst] += multiplicity
        for _ in 1:multiplicity
            push!(adj_list[src], dst)
        end
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
    total_edges_used = Graphs.ne(graph) - sum(length(neighbors) for neighbors in values(adj_list))

    if total_edges_used != Graphs.ne(graph)
        return Int[]  # Failed to find Eulerian path
    end

    return path
end

function _find_eulerian_path_hierholzer_labels(adj_list::Dict{T, Vector{T}}, start_vertex::T) where T
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
"""
function path_to_sequence(path::Vector{T}, graph::MetaGraphsNext.MetaGraph) where T
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
"""
function path_to_sequence(path::GraphPath{T}, graph::MetaGraphsNext.MetaGraph) where T
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
function _extract_sequence_from_step(step::WalkStep{T}, graph, SequenceType) where T
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
function _get_sequence_type_for_path(step::WalkStep{T}, graph) where T
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
function _get_empty_sequence_for_label_type(::Type{T}) where T
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
