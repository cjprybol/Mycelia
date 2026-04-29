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

# Example
```julia
repeats = Mycelia.Rhizomorph.resolve_repeats_next(graph; min_repeat_length=10)
```
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

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to analyze
- `min_repeat_length::Int=10`: Minimum local subgraph size to retain as a repeat

# Returns
- `Vector{RepeatRegion}`: Merged repeat calls for the input graph

# Example
```julia
repeats = Mycelia.Rhizomorph.resolve_repeats_next(graph; min_repeat_length=8)
```
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
    calculate_degrees(graph::MetaGraphsNext.MetaGraph)

Compute in-degree and out-degree tables for each vertex label in `graph`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to summarize

# Returns
- `Tuple{Dict, Dict}`: `(in_degrees, out_degrees)` keyed by vertex label

# Example
```julia
in_degrees, out_degrees = Mycelia.Rhizomorph.calculate_degrees(graph)
```
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
    find_repeat_candidates(in_degrees::Dict{T, Int}, out_degrees::Dict{T, Int}) where T

Select vertices whose local degree pattern suggests repeat structure.

# Arguments
- `in_degrees::Dict{T, Int}`: In-degree per vertex
- `out_degrees::Dict{T, Int}`: Out-degree per vertex

# Returns
- `Vector{T}`: Candidate repeat vertices

# Example
```julia
candidates = Mycelia.Rhizomorph.find_repeat_candidates(in_degrees, out_degrees)
```
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
    analyze_repeat_region(graph::MetaGraphsNext.MetaGraph, start_vertex, min_length)

Analyze the local neighborhood around `start_vertex` and, when supported,
construct a `RepeatRegion`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to analyze
- `start_vertex`: Candidate repeat seed
- `min_length::Int`: Minimum region size required for a call

# Returns
- `Union{Nothing, RepeatRegion}`: Characterized repeat region or `nothing`

# Example
```julia
region = Mycelia.Rhizomorph.analyze_repeat_region(graph, candidate, 10)
```
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
    get_local_subgraph(graph::MetaGraphsNext.MetaGraph, center_vertex, radius)

Collect vertices reachable within `radius` steps of `center_vertex`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to traverse
- `center_vertex`: Seed vertex label
- `radius::Int`: Search radius in graph steps

# Returns
- `Vector`: Vertex labels in the local neighborhood

# Example
```julia
local_vertices = Mycelia.Rhizomorph.get_local_subgraph(graph, candidate, 5)
```
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
    estimate_copy_number(graph::MetaGraphsNext.MetaGraph, vertices::Vector)

Estimate repeat copy number from the aggregate support of `vertices`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Source graph
- `vertices::Vector`: Vertices assigned to the repeat region

# Returns
- `Float64`: Approximate copy-number estimate bounded below by `1.0`

# Example
```julia
copy_number = Mycelia.Rhizomorph.estimate_copy_number(graph, local_vertices)
```
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
    classify_repeat_type(graph::MetaGraphsNext.MetaGraph, vertices, incoming_edges, outgoing_edges)

Classify a repeat as tandem, interspersed, or palindromic from its boundary
connectivity.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Source graph
- `vertices::Vector`: Repeat vertices
- `incoming_edges::Vector{Tuple}`: Edges entering the repeat
- `outgoing_edges::Vector{Tuple}`: Edges leaving the repeat

# Returns
- `Symbol`: One of `:tandem`, `:interspersed`, or `:palindromic`

# Example
```julia
repeat_type = Mycelia.Rhizomorph.classify_repeat_type(graph, vertices, incoming_edges, outgoing_edges)
```
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
    calculate_repeat_confidence(graph::MetaGraphsNext.MetaGraph, vertices::Vector, copy_number::Float64)

Score confidence for a repeat call from region size and copy-number estimate.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Source graph
- `vertices::Vector`: Repeat vertices
- `copy_number::Float64`: Estimated copy number

# Returns
- `Float64`: Confidence score in `[0, 1]`

# Example
```julia
confidence = Mycelia.Rhizomorph.calculate_repeat_confidence(graph, vertices, copy_number)
```
"""
function calculate_repeat_confidence(graph::MetaGraphsNext.MetaGraph,
        vertices::Vector,
        copy_number::Float64)
    size_score = min(1.0, length(vertices) / 20.0)
    copy_score = min(1.0, (copy_number - 1.0) / 5.0)

    return (size_score + copy_score) / 2.0
end

"""
    merge_overlapping_repeats(repeats::Vector{RepeatRegion{T}}) where T

Merge repeat calls that share one or more vertices.

# Arguments
- `repeats::Vector{RepeatRegion{T}}`: Repeat calls to merge

# Returns
- `Vector{RepeatRegion{T}}`: Repeat regions after overlap consolidation

# Example
```julia
merged = Mycelia.Rhizomorph.merge_overlapping_repeats(repeats)
```
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
    regions_overlap(r1::RepeatRegion{T}, r2::RepeatRegion{T}) where T

Test whether two repeat calls share any vertex labels.

# Arguments
- `r1::RepeatRegion{T}`: First repeat region
- `r2::RepeatRegion{T}`: Second repeat region

# Returns
- `Bool`: `true` when the regions overlap

# Example
```julia
overlap = Mycelia.Rhizomorph.regions_overlap(repeats[1], repeats[2])
```
"""
function regions_overlap(r1::RepeatRegion{T}, r2::RepeatRegion{T}) where {T}
    return !isempty(intersect(Set(r1.repeat_vertices), Set(r2.repeat_vertices)))
end

"""
    merge_repeat_regions(regions::Vector{RepeatRegion{T}}) where T

Merge multiple repeat calls into a single summary region.

# Arguments
- `regions::Vector{RepeatRegion{T}}`: Repeat regions to combine

# Returns
- `RepeatRegion{T}`: Union of vertices, boundary edges, and averaged scores

# Example
```julia
merged = Mycelia.Rhizomorph.merge_repeat_regions(repeats[1:2])
```
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
