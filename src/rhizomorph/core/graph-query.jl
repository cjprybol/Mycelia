# Graph Query and Traversal Functions
#
# Functions for querying graph properties and traversing paths through
# strand-specific k-mer and qualmer graphs.

# ============================================================================
# Basic Graph Properties
# ============================================================================

"""
    vertex_count(graph)

Return the number of vertices in the graph.
"""
function vertex_count(graph)
    return length(MetaGraphsNext.labels(graph))
end

"""
    edge_count(graph)

Return the number of edges in the graph.
"""
function edge_count(graph)
    return MetaGraphsNext.ne(graph)
end

"""
    has_vertex(graph, kmer)

Check if a vertex exists in the graph for the given k-mer.
"""
function has_vertex(graph, kmer)
    return haskey(graph, kmer)
end

"""
    has_edge(graph, src, dst)

Check if an edge exists between two k-mers.
"""
function has_edge(graph, src, dst)
    return haskey(graph, src, dst)
end

# ============================================================================
# Vertex and Edge Data Access
# ============================================================================

"""
    get_vertex_data(graph, kmer)

Get the vertex data for a k-mer. Returns `nothing` if vertex doesn't exist.
"""
function get_vertex_data(graph, kmer)
    if !haskey(graph, kmer)
        return nothing
    end
    return graph[kmer]
end

"""
    get_edge_data(graph, src, dst)

Get the edge data between two k-mers. Returns `nothing` if edge doesn't exist.
"""
function get_edge_data(graph, src, dst)
    if !haskey(graph, src, dst)
        return nothing
    end
    return graph[src, dst]
end

"""
    get_vertex_observation_count(graph, kmer)

Get the total number of observations for a vertex across all datasets.
"""
function get_vertex_observation_count(graph, kmer)
    vertex_data = get_vertex_data(graph, kmer)
    if isnothing(vertex_data)
        return 0
    end
    return count_total_observations(vertex_data)
end

"""
    get_edge_observation_count(graph, src, dst)

Get the total number of observations for an edge across all datasets.
"""
function get_edge_observation_count(graph, src, dst)
    edge_data = get_edge_data(graph, src, dst)
    if isnothing(edge_data)
        return 0
    end
    return count_total_observations(edge_data)
end

# ============================================================================
# Graph Traversal - Neighbors
# ============================================================================

"""
    get_outgoing_neighbors(graph, kmer)

Get all k-mers that are direct successors of the given k-mer.
"""
function get_outgoing_neighbors(graph, kmer)
    if !haskey(graph, kmer)
        return []
    end
    return collect(MetaGraphsNext.outneighbor_labels(graph, kmer))
end

"""
    get_incoming_neighbors(graph, kmer)

Get all k-mers that are direct predecessors of the given k-mer.
"""
function get_incoming_neighbors(graph, kmer)
    if !haskey(graph, kmer)
        return []
    end
    return collect(MetaGraphsNext.inneighbor_labels(graph, kmer))
end

"""
    get_outgoing_degree(graph, kmer)

Get the number of outgoing edges from a k-mer.
"""
function get_outgoing_degree(graph, kmer)
    return length(get_outgoing_neighbors(graph, kmer))
end

"""
    get_incoming_degree(graph, kmer)

Get the number of incoming edges to a k-mer.
"""
function get_incoming_degree(graph, kmer)
    return length(get_incoming_neighbors(graph, kmer))
end

# ============================================================================
# Path Classification
# ============================================================================

"""
    is_linear_path(graph, kmer)

Check if a k-mer is part of a linear path (exactly one incoming and one outgoing edge).
"""
function is_linear_path(graph, kmer)
    return get_incoming_degree(graph, kmer) == 1 && get_outgoing_degree(graph, kmer) == 1
end

"""
    is_source_vertex(graph, kmer)

Check if a k-mer is a source vertex (no incoming edges).
"""
function is_source_vertex(graph, kmer)
    return get_incoming_degree(graph, kmer) == 0
end

"""
    is_sink_vertex(graph, kmer)

Check if a k-mer is a sink vertex (no outgoing edges).
"""
function is_sink_vertex(graph, kmer)
    return get_outgoing_degree(graph, kmer) == 0
end

"""
    get_all_sources(graph)

Get all source vertices in the graph (vertices with no incoming edges).
"""
function get_all_sources(graph)
    sources = []
    for kmer in MetaGraphsNext.labels(graph)
        if is_source_vertex(graph, kmer)
            push!(sources, kmer)
        end
    end
    return sources
end

"""
    get_all_sinks(graph)

Get all sink vertices in the graph (vertices with no outgoing edges).
"""
function get_all_sinks(graph)
    sinks = []
    for kmer in MetaGraphsNext.labels(graph)
        if is_sink_vertex(graph, kmer)
            push!(sinks, kmer)
        end
    end
    return sinks
end

# ============================================================================
# Source and Sink Sequence Helpers
# ============================================================================

"""
    find_source_sequences(graph::MetaGraphsNext.MetaGraph)

Return labels for vertices with no incoming edges (assembly start points).
"""
function find_source_sequences(graph::MetaGraphsNext.MetaGraph)
    return get_all_sources(graph)
end

"""
    find_sink_sequences(graph::MetaGraphsNext.MetaGraph)

Return labels for vertices with no outgoing edges (assembly end points).
"""
function find_sink_sequences(graph::MetaGraphsNext.MetaGraph)
    return get_all_sinks(graph)
end

# """
#     find_source_sequences(graph)

# Find sequences that have no incoming edges (potential assembly starting points).

# # Returns
# - `Vector`: Sequences with no predecessors

# # Examples
# ```julia
# sources = find_source_sequences(graph)
# println("Found \$(length(sources)) source sequences")
# ```
# """
# function find_source_sequences(graph::MetaGraphsNext.MetaGraph)
#     sources = []

#     for seq_label in MetaGraphsNext.labels(graph)
#         # Check if this vertex has any incoming edges
#         has_incoming = false
#         for (src, dst) in MetaGraphsNext.edge_labels(graph)
#             if dst == seq_label
#                 has_incoming = true
#                 break
#             end
#         end

#         if !has_incoming
#             push!(sources, seq_label)
#         end
#     end

#     return sources
# end

# """
#     find_sink_sequences(graph)

# Find sequences that have no outgoing edges (potential assembly ending points).

# # Returns
# - `Vector`: Sequences with no successors

# # Examples
# ```julia
# sinks = find_sink_sequences(graph)
# println("Found \$(length(sinks)) sink sequences")
# ```
# """
# function find_sink_sequences(graph::MetaGraphsNext.MetaGraph)
#     sinks = []

#     for seq_label in MetaGraphsNext.labels(graph)
#         # Check if this vertex has any outgoing edges
#         has_outgoing = false
#         for (src, dst) in MetaGraphsNext.edge_labels(graph)
#             if src == seq_label
#                 has_outgoing = true
#                 break
#             end
#         end

#         if !has_outgoing
#             push!(sinks, seq_label)
#         end
#     end

#     return sinks
# end

# ============================================================================
# Unitig Extension
# ============================================================================

"""
    extend_unitig_forward(graph, start_kmer)

Extend a unitig forward from the starting k-mer, following unique linear paths.
Stops when reaching a branch point, sink, or cycle.

Returns a vector of k-mers representing the path.
"""
function extend_unitig_forward(graph, start_kmer)
    if !haskey(graph, start_kmer)
        return []
    end

    path = [start_kmer]
    visited = Set([start_kmer])
    current = start_kmer

    while true
        neighbors = get_outgoing_neighbors(graph, current)

        # Stop at sink or branch point
        if length(neighbors) != 1
            break
        end

        next_kmer = neighbors[1]

        # Stop at cycle
        if next_kmer in visited
            break
        end

        # Stop if next k-mer has multiple incoming edges (join point)
        if get_incoming_degree(graph, next_kmer) != 1
            break
        end

        push!(path, next_kmer)
        push!(visited, next_kmer)
        current = next_kmer
    end

    return path
end

"""
    extend_unitig_backward(graph, start_kmer)

Extend a unitig backward from the starting k-mer, following unique linear paths.
Stops when reaching a branch point, source, or cycle.

Returns a vector of k-mers representing the path (in reverse order from source).
"""
function extend_unitig_backward(graph, start_kmer)
    if !haskey(graph, start_kmer)
        return []
    end

    path = [start_kmer]
    visited = Set([start_kmer])
    current = start_kmer

    while true
        neighbors = get_incoming_neighbors(graph, current)

        # Stop at source or join point
        if length(neighbors) != 1
            break
        end

        prev_kmer = neighbors[1]

        # Stop at cycle
        if prev_kmer in visited
            break
        end

        # Stop if previous k-mer has multiple outgoing edges (branch point)
        if get_outgoing_degree(graph, prev_kmer) != 1
            break
        end

        push!(path, prev_kmer)
        push!(visited, prev_kmer)
        current = prev_kmer
    end

    return path
end

"""
    get_maximal_unitig(graph, kmer)

Get the maximal unitig containing the given k-mer.
Extends both forward and backward to find the full linear path.

Returns a vector of k-mers representing the complete unitig.
"""
function get_maximal_unitig(graph, kmer)
    if !haskey(graph, kmer)
        return []
    end

    # Extend backward (returns path in reverse from kmer to source)
    backward_path = extend_unitig_backward(graph, kmer)

    # Extend forward (returns path from kmer to sink)
    forward_path = extend_unitig_forward(graph, kmer)

    # Combine: reverse(backward) + forward[2:end] (excluding duplicate kmer)
    # backward_path[1] == forward_path[1] == kmer
    if length(backward_path) > 1
        # Reverse backward path (excluding the starting kmer which is already in forward)
        unitig = reverse(backward_path[2:end])
        append!(unitig, forward_path)
    else
        unitig = forward_path
    end

    return unitig
end

"""
    assemble_path_sequence(path::Vector{<:Kmers.Kmer})

Reconstruct the sequence represented by a path of k-mers.
K-mers must overlap by k-1 bases.
"""
function assemble_path_sequence(path::Vector{T}) where T <: Kmers.Kmer
    if isempty(path)
        error("Cannot assemble empty path")
    end

    if length(path) == 1
        # Single k-mer - return as sequence
        kmer = path[1]
        if T <: Kmers.DNAKmer
            return BioSequences.LongDNA{4}(kmer)
        elseif T <: Kmers.RNAKmer
            return BioSequences.LongRNA{4}(kmer)
        elseif T <: Kmers.AAKmer
            return BioSequences.LongAA(kmer)
        else
            error("Unsupported k-mer type: $(T)")
        end
    end

    # Start with first k-mer
    first_kmer = path[1]
    if T <: Kmers.DNAKmer
        sequence = BioSequences.LongDNA{4}(first_kmer)
    elseif T <: Kmers.RNAKmer
        sequence = BioSequences.LongRNA{4}(first_kmer)
    elseif T <: Kmers.AAKmer
        sequence = BioSequences.LongAA(first_kmer)
    else
        error("Unsupported k-mer type: $(T)")
    end

    # Add last base of each subsequent k-mer
    for i in 2:length(path)
        kmer = path[i]
        last_base = kmer[Kmers.ksize(T)]  # Get last base
        push!(sequence, last_base)
    end

    return sequence
end

# ============================================================================
# Graph Analysis
# ============================================================================

"""
    get_all_vertices(graph)

Get all vertex labels (k-mers) in the graph.
"""
function get_all_vertices(graph)
    return collect(MetaGraphsNext.labels(graph))
end

"""
    filter_vertices_by_observation_count(graph, min_count::Int)

Get all vertices with at least `min_count` observations.
"""
function filter_vertices_by_observation_count(graph, min_count::Int)
    filtered = []
    for kmer in MetaGraphsNext.labels(graph)
        if get_vertex_observation_count(graph, kmer) >= min_count
            push!(filtered, kmer)
        end
    end
    return filtered
end

"""
    get_graph_statistics(graph)

Compute basic statistics about the graph structure.

Returns a dictionary with:
- `:vertex_count`: Number of vertices
- `:edge_count`: Number of edges
- `:source_count`: Number of source vertices
- `:sink_count`: Number of sink vertices
- `:linear_count`: Number of linear path vertices
- `:branch_count`: Number of branch vertices (out-degree > 1)
- `:join_count`: Number of join vertices (in-degree > 1)
"""
function get_graph_statistics(graph)
    stats = Dict{Symbol, Int}()

    stats[:vertex_count] = vertex_count(graph)
    stats[:edge_count] = edge_count(graph)
    stats[:source_count] = 0
    stats[:sink_count] = 0
    stats[:linear_count] = 0
    stats[:branch_count] = 0
    stats[:join_count] = 0

    for kmer in MetaGraphsNext.labels(graph)
        in_deg = get_incoming_degree(graph, kmer)
        out_deg = get_outgoing_degree(graph, kmer)

        if in_deg == 0
            stats[:source_count] += 1
        end
        if out_deg == 0
            stats[:sink_count] += 1
        end
        if in_deg == 1 && out_deg == 1
            stats[:linear_count] += 1
        end
        if out_deg > 1
            stats[:branch_count] += 1
        end
        if in_deg > 1
            stats[:join_count] += 1
        end
    end

    return stats
end
