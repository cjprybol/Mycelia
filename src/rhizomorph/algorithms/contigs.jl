# Contig Extraction Algorithms
#
# Ports contig finding utilities from legacy graph code and adapts them
# to Rhizomorph's type-safe vertex data (k-mers, BioSequences, and strings).

"""
    ContigPath{T, S}

Represents a linear path through the graph forming a contig.
"""
struct ContigPath{T, S}
    vertices::Vector{T}
    sequence::S
    coverage_profile::Vector{Float64}
    length::Int
    n50_contribution::Float64

    function ContigPath(vertices::Vector{T}, sequence::S, coverage::Vector{Float64}) where {
            T, S}
        new{T, S}(vertices, sequence, coverage, length(sequence), 0.0)
    end
end

"""
    find_contigs_next(graph::MetaGraphsNext.MetaGraph; min_contig_length::Int=500)

Extract linear contigs from the assembly graph. Returns contigs whose assembled
sequence length is at least `min_contig_length`.
"""
function find_contigs_next(graph::MetaGraphsNext.MetaGraph; min_contig_length::Int = 500)
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return ContigPath{Any, Any}[]
    end

    contigs = Vector{ContigPath}()
    visited = Set{eltype(labels)}()

    for start_vertex in labels
        if start_vertex in visited
            continue
        end

        path = find_linear_path(graph, start_vertex, visited)

        if !isempty(path)
            sequence = generate_contig_sequence(graph, path)
            coverage_profile = generate_coverage_profile(graph, path)

            if length(sequence) >= min_contig_length
                push!(contigs, ContigPath(path, sequence, coverage_profile))
            end

            for vertex in path
                push!(visited, vertex)
            end
        end
    end

    return sort_contigs_by_length(contigs)
end

"""
Find a linear path through the graph starting from a vertex.
"""
function find_linear_path(graph::MetaGraphsNext.MetaGraph, start_vertex, visited::Set)
    if start_vertex in visited
        return Vector{typeof(start_vertex)}()
    end

    path = [start_vertex]
    current = start_vertex

    while true
        out_neighbors = Rhizomorph.get_outgoing_neighbors(graph, current)
        valid_neighbors = [n for n in out_neighbors if !(n in visited) && !(n in path)]

        next_vertex = _select_linear_neighbor(graph, current, valid_neighbors; direction = :out)
        if next_vertex === nothing
            break
        end

        if !_edge_is_dominant_among_incoming(graph, current, next_vertex)
            break
        end

        push!(path, next_vertex)
        current = next_vertex
    end

    current = start_vertex
    backward_path = Vector{typeof(start_vertex)}()

    while true
        in_neighbors = Rhizomorph.get_incoming_neighbors(graph, current)
        valid_neighbors = [n
                           for n in in_neighbors
                           if !(n in visited) && !(n in path) && !(n in backward_path)]

        prev_vertex = _select_linear_neighbor(graph, current, valid_neighbors; direction = :in)
        if prev_vertex === nothing
            break
        end

        if !_edge_is_dominant_among_outgoing(graph, prev_vertex, current)
            break
        end

        pushfirst!(backward_path, prev_vertex)
        current = prev_vertex
    end

    return [backward_path; path]
end

const _DOMINANT_EDGE_RATIO = 2.0

function _edge_support(graph, src, dst)
    if !haskey(graph, src, dst)
        return 0
    end
    edge_data = graph[src, dst]
    if hasfield(typeof(edge_data), :evidence)
        return Rhizomorph.count_evidence(edge_data)
    end
    if hasfield(typeof(edge_data), :coverage)
        return length(edge_data.coverage)
    end
    return 1
end

function _select_linear_neighbor(graph, current, neighbors; direction::Symbol)
    if isempty(neighbors)
        return nothing
    end
    if length(neighbors) == 1
        return neighbors[1]
    end
    return _dominant_neighbor(graph, current, neighbors; direction = direction)
end

function _dominant_neighbor(graph, current, neighbors; direction::Symbol)
    supports = Vector{Tuple{eltype(neighbors), Int}}(undef, length(neighbors))
    for (idx, neighbor) in enumerate(neighbors)
        support = direction == :out ? _edge_support(graph, current, neighbor) :
                  _edge_support(graph, neighbor, current)
        supports[idx] = (neighbor, support)
    end
    sort!(supports, by = last, rev = true)

    best_neighbor, best_support = supports[1]
    second_support = length(supports) > 1 ? supports[2][2] : 0

    if best_support >= _DOMINANT_EDGE_RATIO * max(second_support, 1)
        return best_neighbor
    end

    return nothing
end

function _edge_is_dominant_among_incoming(graph, src, dst)
    incoming = Rhizomorph.get_incoming_neighbors(graph, dst)
    if length(incoming) <= 1
        return true
    end
    dominant = _dominant_neighbor(graph, dst, incoming; direction = :in)
    return dominant == src
end

function _edge_is_dominant_among_outgoing(graph, src, dst)
    outgoing = Rhizomorph.get_outgoing_neighbors(graph, src)
    if length(outgoing) <= 1
        return true
    end
    dominant = _dominant_neighbor(graph, src, outgoing; direction = :out)
    return dominant == dst
end

"""
Generate an assembled sequence for a contig path.
"""
function generate_contig_sequence(
        graph::MetaGraphsNext.MetaGraph{
            <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any},
        path::Vector{T}) where {T}
    if isempty(path)
        return ""
    end

    first_vertex_data = graph[path[1]]

    if hasfield(typeof(first_vertex_data), :Kmer)
        # K-mer graph: reconstruct using overlap logic without string conversion
        return Rhizomorph.assemble_path_sequence(path)
    elseif hasfield(typeof(first_vertex_data), :sequence)
        # Variable-length BioSequence graph
        assembled = first_vertex_data.sequence
        SequenceType = typeof(assembled)

        for i in 2:length(path)
            src = path[i - 1]
            dst = path[i]
            edge_overlap = _edge_overlap_length(graph, src, dst)
            next_sequence = graph[dst].sequence

            overlap = min(edge_overlap, length(next_sequence))
            append_part = overlap < length(next_sequence) ?
                          next_sequence[(overlap + 1):end] : SequenceType()
            assembled = assembled * append_part
        end

        return assembled
    else
        # String graphs or other label-only paths
        sequence = string(path[1])

        for i in 2:length(path)
            sequence *= string(path[i])[end:end]
        end

        return sequence
    end
end

"""
Generate coverage profile for a contig path based on evidence counts.
"""
function generate_coverage_profile(graph::MetaGraphsNext.MetaGraph, path::Vector)
    coverage = Float64[]

    for vertex in path
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            vertex_coverage = Float64(Rhizomorph.count_evidence(vertex_data))
            push!(coverage, vertex_coverage)
        else
            push!(coverage, 0.0)
        end
    end

    return coverage
end

"""
Sort contigs by sequence length (descending).
"""
function sort_contigs_by_length(contigs::Vector{ContigPath})
    return sort(contigs, by = c -> c.length, rev = true)
end

"""
Retrieve overlap length from edge data when available.
"""
function _edge_overlap_length(graph::MetaGraphsNext.MetaGraph, src, dst)
    if haskey(graph, src, dst)
        edge_data = graph[src, dst]
        if hasfield(typeof(edge_data), :overlap_length)
            return edge_data.overlap_length
        end
    end
    return 0
end
