# Repeat Resolution Algorithms
#
# Ports repeat detection utilities from legacy graph code into the
# Rhizomorph namespace while adopting the new evidence-centric vertex types.

"""
    RepeatRegion{T}

Represents a repetitive region in the assembly graph.

# Fields
- `repeat_vertices::Vector{T}`: Vertices that belong to the repeat
- `incoming_edges::Vector{Tuple{T, T}}`: Edges entering the repeat region
- `outgoing_edges::Vector{Tuple{T, T}}`: Edges leaving the repeat region
- `copy_number_estimate::Float64`: Estimated repeat copy number
- `repeat_type::Symbol`: Repeat classification (`:tandem`, `:interspersed`, `:palindromic`)
- `confidence::Float64`: Confidence score in `[0, 1]`
"""
struct RepeatRegion{T}
    repeat_vertices::Vector{T}
    incoming_edges::Vector{Tuple{T, T}}
    outgoing_edges::Vector{Tuple{T, T}}
    copy_number_estimate::Float64
    repeat_type::Symbol
    confidence::Float64

    function RepeatRegion(vertices::Vector{T},
            incoming::Vector{Tuple{T, T}},
            outgoing::Vector{Tuple{T, T}},
            copy_num::Float64,
            rep_type::Symbol,
            conf::Float64) where {T}
        @assert rep_type in (:tandem, :interspersed, :palindromic) "Invalid repeat type"
        @assert 0.0 <= conf <= 1.0 "Confidence must be in [0,1]"
        new{T}(vertices, incoming, outgoing, copy_num, rep_type, conf)
    end
end

"""
    resolve_repeats_next(graph::MetaGraphsNext.MetaGraph; min_repeat_length::Int=10)

Identify and characterize repetitive regions in the assembly graph.
Returns a vector of `RepeatRegion` objects specific to the graph's vertex label type.
"""
function resolve_repeats_next(graph::MetaGraphsNext.MetaGraph; min_repeat_length::Int = 10)
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return RepeatRegion{Any}[]
    end

    vertex_type = eltype(labels)
    repeats = RepeatRegion{vertex_type}[]

    in_degrees, out_degrees = calculate_degrees(graph)
    repeat_candidates = find_repeat_candidates(in_degrees, out_degrees)

    for candidate_vertex in repeat_candidates
        repeat_region = analyze_repeat_region(graph, candidate_vertex, min_repeat_length)
        if repeat_region !== nothing
            push!(repeats, repeat_region)
        end
    end

    return merge_overlapping_repeats(repeats)
end

"""
Calculate in-degrees and out-degrees for all vertices.
"""
function calculate_degrees(graph::MetaGraphsNext.MetaGraph{
        <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any}) where {T}
    vertices = collect(MetaGraphsNext.labels(graph))
    in_degrees = Dict{T, Int}()
    out_degrees = Dict{T, Int}()

    for vertex in vertices
        in_degrees[vertex] = 0
        out_degrees[vertex] = 0
    end

    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        out_degrees[src] = get(out_degrees, src, 0) + 1
        in_degrees[dst] = get(in_degrees, dst, 0) + 1
    end

    return in_degrees, out_degrees
end

"""
Find vertices that could be part of repeat regions.
"""
function find_repeat_candidates(in_degrees::Dict{T, Int}, out_degrees::Dict{
        T, Int}) where {T}
    candidates = T[]

    for vertex in keys(in_degrees)
        total_degree = in_degrees[vertex] + out_degrees[vertex]
        if total_degree > 4  # Threshold for repeat consideration
            push!(candidates, vertex)
        end
    end

    return candidates
end

"""
Analyze a potential repeat region starting from a vertex.
"""
function analyze_repeat_region(
        graph::MetaGraphsNext.MetaGraph{
            <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any},
        start_vertex::T,
        min_length::Int) where {T}
    local_vertices = get_local_subgraph(graph, start_vertex, min_length)

    if length(local_vertices) < min_length
        return nothing
    end

    incoming_edges = Tuple{T, T}[]
    outgoing_edges = Tuple{T, T}[]

    for vertex in local_vertices
        in_neighbors = Rhizomorph.get_incoming_neighbors(graph, vertex)
        out_neighbors = Rhizomorph.get_outgoing_neighbors(graph, vertex)

        for neighbor in in_neighbors
            if !(neighbor in local_vertices)
                push!(incoming_edges, (neighbor, vertex))
            end
        end

        for neighbor in out_neighbors
            if !(neighbor in local_vertices)
                push!(outgoing_edges, (vertex, neighbor))
            end
        end
    end

    copy_number = estimate_copy_number(graph, local_vertices)
    repeat_type = classify_repeat_type(graph, local_vertices, incoming_edges, outgoing_edges)
    confidence = calculate_repeat_confidence(graph, local_vertices, copy_number)

    return RepeatRegion(local_vertices, incoming_edges, outgoing_edges,
        copy_number, repeat_type, confidence)
end

"""
Get local subgraph around a vertex within the provided radius.
"""
function get_local_subgraph(
        graph::MetaGraphsNext.MetaGraph{
            <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any},
        center_vertex::T,
        radius::Int) where {T}
    visited = Set{T}()
    queue = DataStructures.Queue{Tuple{T, Int}}()
    DataStructures.enqueue!(queue, (center_vertex, 0))

    while !isempty(queue)
        vertex, distance = DataStructures.dequeue!(queue)

        if vertex in visited || distance > radius
            continue
        end

        push!(visited, vertex)

        for neighbor in Rhizomorph.get_outgoing_neighbors(graph, vertex)
            if !(neighbor in visited)
                DataStructures.enqueue!(queue, (neighbor, distance + 1))
            end
        end

        for neighbor in Rhizomorph.get_incoming_neighbors(graph, vertex)
            if !(neighbor in visited)
                DataStructures.enqueue!(queue, (neighbor, distance + 1))
            end
        end
    end

    return collect(visited)
end

"""
Estimate copy number for repeat region based on evidence counts.
"""
function estimate_copy_number(
        graph::MetaGraphsNext.MetaGraph{
            <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any},
        vertices::Vector{T}) where {T}
    if isempty(vertices)
        return 1.0
    end

    total_evidence = 0.0
    for vertex in vertices
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            total_evidence += Rhizomorph.count_evidence(vertex_data)
        end
    end

    avg_support = total_evidence / length(vertices)
    expected_single_copy = 10.0  # Placeholder until dataset-derived estimates are available

    return max(1.0, avg_support / expected_single_copy)
end

"""
Classify the type of repeat based on connectivity.
"""
function classify_repeat_type(
        graph::MetaGraphsNext.MetaGraph{
            <:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any},
        vertices::Vector{T},
        incoming_edges::Vector{Tuple{T, T}},
        outgoing_edges::Vector{Tuple{T, T}}) where {T}
    n_incoming = length(incoming_edges)
    n_outgoing = length(outgoing_edges)

    if n_incoming <= 2 && n_outgoing <= 2
        return :tandem
    elseif n_incoming > 2 || n_outgoing > 2
        return :interspersed
    else
        return :palindromic
    end
end

"""
Calculate confidence in repeat identification.
"""
function calculate_repeat_confidence(graph::MetaGraphsNext.MetaGraph,
        vertices::Vector,
        copy_number::Float64)
    size_score = min(1.0, length(vertices) / 20.0)
    copy_score = min(1.0, (copy_number - 1.0) / 5.0)

    return (size_score + copy_score) / 2.0
end

"""
Merge overlapping repeat regions.
"""
function merge_overlapping_repeats(repeats::Vector{RepeatRegion{T}}) where {T}
    if isempty(repeats)
        return repeats
    end

    merged = RepeatRegion{T}[]
    used = falses(length(repeats))

    for i in 1:length(repeats)
        if used[i]
            continue
        end

        current_repeat = repeats[i]
        overlapping_indices = [i]

        for j in (i + 1):length(repeats)
            if !used[j] && regions_overlap(current_repeat, repeats[j])
                push!(overlapping_indices, j)
                used[j] = true
            end
        end

        if length(overlapping_indices) > 1
            merged_repeat = merge_repeat_regions([repeats[idx]
                                                  for idx in overlapping_indices])
            push!(merged, merged_repeat)
        else
            push!(merged, current_repeat)
        end

        used[i] = true
    end

    return merged
end

"""
Check if two repeat regions overlap.
"""
function regions_overlap(r1::RepeatRegion{T}, r2::RepeatRegion{T}) where {T}
    return !isempty(intersect(Set(r1.repeat_vertices), Set(r2.repeat_vertices)))
end

"""
Merge multiple repeat regions into one.
"""
function merge_repeat_regions(regions::Vector{RepeatRegion{T}}) where {T}
    all_vertices = Vector{T}()
    all_incoming = Tuple{T, T}[]
    all_outgoing = Tuple{T, T}[]

    for region in regions
        append!(all_vertices, region.repeat_vertices)
        append!(all_incoming, region.incoming_edges)
        append!(all_outgoing, region.outgoing_edges)
    end

    unique_vertices = unique(all_vertices)
    unique_incoming = unique(all_incoming)
    unique_outgoing = unique(all_outgoing)

    avg_copy_number = Statistics.mean(r.copy_number_estimate for r in regions)
    avg_confidence = Statistics.mean(r.confidence for r in regions)
    repeat_type = _mode_symbol([r.repeat_type for r in regions])

    return RepeatRegion(unique_vertices, unique_incoming, unique_outgoing,
        avg_copy_number, repeat_type, avg_confidence)
end

"""
Compute the most common symbol in a collection (ties pick the first max).
"""
function _mode_symbol(symbols::Vector{Symbol})
    counts = Dict{Symbol, Int}()
    for sym in symbols
        counts[sym] = get(counts, sym, 0) + 1
    end
    max_sym = first(keys(counts))
    max_count = counts[max_sym]
    for (sym, count) in counts
        if count > max_count
            max_sym = sym
            max_count = count
        end
    end
    return max_sym
end
