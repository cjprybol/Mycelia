"""
Advanced Graph Algorithms for Next-Generation Assembly

This module provides sophisticated graph algorithms for assembly graph analysis,
including Eulerian path finding, bubble detection, repeat resolution, and
graph simplification operations.

All algorithms are designed to work with strand-aware MetaGraphsNext graphs
and maintain biological correctness constraints.
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import DataStructures
import Statistics

export find_eulerian_paths_next, detect_bubbles_next, resolve_repeats_next
export simplify_graph_next, find_contigs_next, scaffold_contigs_next
export BubbleStructure, RepeatRegion, ContigPath, ScaffoldResult

"""
    BubbleStructure

Represents a bubble (alternative path) in the assembly graph.
"""
struct BubbleStructure
    entry_vertex::String
    exit_vertex::String
    path1::Vector{String}
    path2::Vector{String}
    path1_support::Int
    path2_support::Int
    complexity_score::Float64
    
    function BubbleStructure(entry::String, exit::String, 
                           p1::Vector{String}, p2::Vector{String},
                           s1::Int, s2::Int, complexity::Float64)
        new(entry, exit, p1, p2, s1, s2, complexity)
    end
end

"""
    RepeatRegion

Represents a repetitive region in the assembly graph.
"""
struct RepeatRegion
    repeat_vertices::Vector{String}
    incoming_edges::Vector{Tuple{String, String}}
    outgoing_edges::Vector{Tuple{String, String}}
    copy_number_estimate::Float64
    repeat_type::Symbol  # :tandem, :interspersed, :palindromic
    confidence::Float64
    
    function RepeatRegion(vertices::Vector{String}, 
                         incoming::Vector{Tuple{String, String}},
                         outgoing::Vector{Tuple{String, String}},
                         copy_num::Float64, rep_type::Symbol, conf::Float64)
        @assert rep_type in [:tandem, :interspersed, :palindromic] "Invalid repeat type"
        @assert 0.0 <= conf <= 1.0 "Confidence must be in [0,1]"
        new(vertices, incoming, outgoing, copy_num, rep_type, conf)
    end
end

"""
    ContigPath

Represents a linear path through the graph forming a contig.
"""
struct ContigPath
    vertices::Vector{String}
    sequence::String
    coverage_profile::Vector{Float64}
    length::Int
    n50_contribution::Float64
    
    function ContigPath(vertices::Vector{String}, sequence::String,
                       coverage::Vector{Float64})
        new(vertices, sequence, coverage, length(sequence), 0.0)
    end
end

"""
    ScaffoldResult

Results from scaffolding analysis.
"""
struct ScaffoldResult
    scaffolds::Vector{Vector{ContigPath}}
    gap_estimates::Vector{Tuple{Int, Int, Float64}}  # (min_gap, max_gap, confidence)
    scaffold_n50::Float64
    total_length::Int
    
    function ScaffoldResult(scaffolds::Vector{Vector{ContigPath}},
                          gaps::Vector{Tuple{Int, Int, Float64}})
        total_len = sum(sum(contig.length for contig in scaffold) for scaffold in scaffolds)
        # Calculate N50 (simplified)
        lengths = [sum(contig.length for contig in scaffold) for scaffold in scaffolds]
        sort!(lengths, rev=true)
        cumsum_lengths = cumsum(lengths)
        n50_idx = findfirst(x -> x >= total_len / 2, cumsum_lengths)
        n50 = n50_idx !== nothing ? lengths[n50_idx] : 0.0
        
        new(scaffolds, gaps, n50, total_len)
    end
end

"""
    find_eulerian_paths_next(graph::MetaGraph) -> Vector{Vector{String}}

Find Eulerian paths in the assembly graph. An Eulerian path visits every edge exactly once.
"""
function find_eulerian_paths_next(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData})
    if isempty(MetaGraphsNext.labels(graph))
        return Vector{String}[]
    end
    
    # Check if Eulerian path exists
    in_degrees, out_degrees = calculate_degrees(graph)
    eulerian_info = check_eulerian_conditions(in_degrees, out_degrees)
    
    if !eulerian_info.has_path
        println("Graph does not have Eulerian path")
        return Vector{String}[]
    end
    
    # Find starting vertices
    start_vertices = find_eulerian_start_vertices(in_degrees, out_degrees, eulerian_info)
    
    paths = Vector{String}[]
    for start_vertex in start_vertices
        path = find_eulerian_path_from_vertex(graph, start_vertex, in_degrees, out_degrees)
        if !isempty(path)
            push!(paths, path)
        end
    end
    
    return paths
end

"""
Calculate in-degrees and out-degrees for all vertices.
"""
function calculate_degrees(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData})
    vertices = collect(MetaGraphsNext.labels(graph))
    in_degrees = Dict{String, Int}()
    out_degrees = Dict{String, Int}()
    
    # Initialize
    for vertex in vertices
        in_degrees[vertex] = 0
        out_degrees[vertex] = 0
    end
    
    # Count degrees
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        out_degrees[src] += 1
        in_degrees[dst] += 1
    end
    
    return in_degrees, out_degrees
end

"""
Check conditions for Eulerian path existence.
"""
function check_eulerian_conditions(in_degrees::Dict{String, Int}, out_degrees::Dict{String, Int})
    odd_vertices = String[]
    start_vertices = String[]
    end_vertices = String[]
    
    for vertex in keys(in_degrees)
        in_deg = in_degrees[vertex]
        out_deg = out_degrees[vertex]
        
        if in_deg != out_deg
            push!(odd_vertices, vertex)
            if out_deg > in_deg
                push!(start_vertices, vertex)
            elseif in_deg > out_deg
                push!(end_vertices, vertex)
            end
        end
    end
    
    # Eulerian path conditions:
    # - All vertices have equal in/out degree (Eulerian cycle), OR
    # - Exactly one vertex has out_degree = in_degree + 1 (start)
    # - Exactly one vertex has in_degree = out_degree + 1 (end)
    # - All other vertices have equal in/out degree
    
    has_path = length(odd_vertices) == 0 || 
               (length(start_vertices) == 1 && length(end_vertices) == 1)
    
    return (has_path=has_path, start_vertices=start_vertices, end_vertices=end_vertices)
end

"""
Find valid starting vertices for Eulerian paths.
"""
function find_eulerian_start_vertices(in_degrees::Dict{String, Int}, 
                                     out_degrees::Dict{String, Int},
                                     eulerian_info)
    if !isempty(eulerian_info.start_vertices)
        return eulerian_info.start_vertices
    else
        # If it's an Eulerian cycle, can start from any vertex with edges
        return [v for v in keys(out_degrees) if out_degrees[v] > 0]
    end
end

"""
Find Eulerian path starting from a specific vertex using Hierholzer's algorithm.
"""
function find_eulerian_path_from_vertex(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                       start_vertex::String,
                                       in_degrees::Dict{String, Int},
                                       out_degrees::Dict{String, Int})
    # Create adjacency list with edge tracking
    adj_list = Dict{String, Vector{Tuple{String, Bool}}}()  # (neighbor, used)
    
    for vertex in keys(out_degrees)
        adj_list[vertex] = Tuple{String, Bool}[]
    end
    
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        push!(adj_list[src], (dst, false))
    end
    
    path = String[]
    stack = [start_vertex]
    current_vertex = start_vertex
    
    while !isempty(stack) || any(any(!used for (_, used) in neighbors) for neighbors in values(adj_list))
        if any(!used for (_, used) in adj_list[current_vertex])
            push!(stack, current_vertex)
            
            # Find first unused edge
            for (i, (neighbor, used)) in enumerate(adj_list[current_vertex])
                if !used
                    adj_list[current_vertex][i] = (neighbor, true)
                    current_vertex = neighbor
                    break
                end
            end
        else
            push!(path, current_vertex)
            if !isempty(stack)
                current_vertex = pop!(stack)
            end
        end
    end
    
    reverse!(path)
    return path
end

"""
    detect_bubbles_next(graph::MetaGraph, min_bubble_length::Int=2, max_bubble_length::Int=100) -> Vector{BubbleStructure}

Detect bubble structures (alternative paths) in the assembly graph.
"""
function detect_bubbles_next(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData};
                           min_bubble_length::Int=2,
                           max_bubble_length::Int=100)
    bubbles = BubbleStructure[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    for entry_vertex in vertices
        # Find potential bubble entry points (vertices with out-degree > 1)
        out_neighbors = get_out_neighbors(graph, entry_vertex)
        
        if length(out_neighbors) >= 2
            # Look for bubbles starting from this vertex
            bubble_candidates = find_bubble_paths(graph, entry_vertex, out_neighbors, 
                                                min_bubble_length, max_bubble_length)
            
            for bubble in bubble_candidates
                if is_valid_bubble(graph, bubble)
                    push!(bubbles, bubble)
                end
            end
        end
    end
    
    return remove_duplicate_bubbles(bubbles)
end

"""
Get outgoing neighbors of a vertex.
"""
function get_out_neighbors(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData}, vertex::String)
    neighbors = String[]
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        if src == vertex
            push!(neighbors, dst)
        end
    end
    return neighbors
end

"""
Get incoming neighbors of a vertex.
"""
function get_in_neighbors(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData}, vertex::String)
    neighbors = String[]
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        if dst == vertex
            push!(neighbors, src)
        end
    end
    return neighbors
end

"""
Find potential bubble paths from an entry vertex.
"""
function find_bubble_paths(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                          entry_vertex::String, 
                          out_neighbors::Vector{String},
                          min_length::Int, max_length::Int)
    bubbles = BubbleStructure[]
    
    # Try all pairs of outgoing paths
    for i in 1:length(out_neighbors)
        for j in (i+1):length(out_neighbors)
            path1_start = out_neighbors[i]
            path2_start = out_neighbors[j]
            
            # Find paths from each starting point
            path1 = find_limited_path(graph, path1_start, max_length)
            path2 = find_limited_path(graph, path2_start, max_length)
            
            # Check if paths reconverge
            convergence_point = find_path_convergence(path1, path2)
            
            if convergence_point !== nothing && 
               length(path1) >= min_length && length(path2) >= min_length
                
                # Extract paths up to convergence
                conv_idx1 = findfirst(v -> v == convergence_point, path1)
                conv_idx2 = findfirst(v -> v == convergence_point, path2)
                
                if conv_idx1 !== nothing && conv_idx2 !== nothing
                    bubble_path1 = path1[1:conv_idx1]
                    bubble_path2 = path2[1:conv_idx2]
                    
                    # Calculate support (simplified - could use actual coverage)
                    support1 = calculate_path_support(graph, bubble_path1)
                    support2 = calculate_path_support(graph, bubble_path2)
                    
                    complexity = calculate_bubble_complexity(bubble_path1, bubble_path2)
                    
                    bubble = BubbleStructure(entry_vertex, convergence_point,
                                           bubble_path1, bubble_path2,
                                           support1, support2, complexity)
                    push!(bubbles, bubble)
                end
            end
        end
    end
    
    return bubbles
end

"""
Find a limited-length path from a starting vertex.
"""
function find_limited_path(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                          start_vertex::String, max_length::Int)
    path = [start_vertex]
    current = start_vertex
    
    for _ in 1:max_length
        neighbors = get_out_neighbors(graph, current)
        if length(neighbors) == 1
            next_vertex = neighbors[1]
            if next_vertex in path  # Avoid cycles
                break
            end
            push!(path, next_vertex)
            current = next_vertex
        else
            break  # Multiple or no neighbors
        end
    end
    
    return path
end

"""
Find where two paths converge.
"""
function find_path_convergence(path1::Vector{String}, path2::Vector{String})
    # Find first common vertex (excluding starting vertices)
    for vertex1 in path1[2:end]
        if vertex1 in path2[2:end]
            return vertex1
        end
    end
    return nothing
end

"""
Calculate support for a path (simplified version).
"""
function calculate_path_support(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                               path::Vector{String})
    if isempty(path)
        return 0
    end
    
    # Use vertex data coverage if available
    total_coverage = 0.0
    for vertex in path
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            # Use metadata if available, otherwise default to 1
            coverage = length(vertex_data.source_positions)
            total_coverage += coverage
        end
    end
    
    return round(Int, total_coverage / length(path))
end

"""
Calculate complexity score for a bubble.
"""
function calculate_bubble_complexity(path1::Vector{String}, path2::Vector{String})
    # Simple complexity metric based on path length difference and sequence similarity
    length_diff = abs(length(path1) - length(path2))
    avg_length = (length(path1) + length(path2)) / 2
    
    length_score = length_diff / max(avg_length, 1)
    
    # Could add sequence similarity scoring here
    similarity_score = 0.5  # Placeholder
    
    return length_score + (1 - similarity_score)
end

"""
Check if a bubble structure is valid.
"""
function is_valid_bubble(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                        bubble::BubbleStructure)
    # Check that paths are distinct
    if bubble.path1 == bubble.path2
        return false
    end
    
    # Check that entry and exit vertices exist
    if !(haskey(graph, bubble.entry_vertex) && haskey(graph, bubble.exit_vertex))
        return false
    end
    
    # Check that paths are valid in graph
    is_valid_path(graph, [bubble.entry_vertex; bubble.path1]) &&
    is_valid_path(graph, [bubble.entry_vertex; bubble.path2])
end

"""
Check if a path is valid in the graph.
"""
function is_valid_path(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                      path::Vector{String})
    if length(path) < 2
        return true
    end
    
    for i in 1:(length(path)-1)
        if !haskey(graph, path[i], path[i+1])
            return false
        end
    end
    return true
end

"""
Remove duplicate bubbles.
"""
function remove_duplicate_bubbles(bubbles::Vector{BubbleStructure})
    unique_bubbles = BubbleStructure[]
    
    for bubble in bubbles
        is_duplicate = false
        for existing in unique_bubbles
            if are_equivalent_bubbles(bubble, existing)
                is_duplicate = true
                break
            end
        end
        
        if !is_duplicate
            push!(unique_bubbles, bubble)
        end
    end
    
    return unique_bubbles
end

"""
Check if two bubbles are equivalent.
"""
function are_equivalent_bubbles(b1::BubbleStructure, b2::BubbleStructure)
    return (b1.entry_vertex == b2.entry_vertex && 
            b1.exit_vertex == b2.exit_vertex &&
            ((b1.path1 == b2.path1 && b1.path2 == b2.path2) ||
             (b1.path1 == b2.path2 && b1.path2 == b2.path1)))
end

"""
    resolve_repeats_next(graph::MetaGraph, min_repeat_length::Int=10) -> Vector{RepeatRegion}

Identify and characterize repetitive regions in the assembly graph.
"""
function resolve_repeats_next(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData};
                            min_repeat_length::Int=10)
    repeats = RepeatRegion[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    # Find vertices with high in-degree or out-degree (potential repeat boundaries)
    in_degrees, out_degrees = calculate_degrees(graph)
    
    # Find potential repeat vertices (high connectivity)
    repeat_candidates = find_repeat_candidates(in_degrees, out_degrees)
    
    # Analyze each candidate region
    for candidate_vertex in repeat_candidates
        repeat_region = analyze_repeat_region(graph, candidate_vertex, min_repeat_length)
        if repeat_region !== nothing
            push!(repeats, repeat_region)
        end
    end
    
    # Merge overlapping repeat regions
    merged_repeats = merge_overlapping_repeats(repeats)
    
    return merged_repeats
end

"""
Find vertices that could be part of repeat regions.
"""
function find_repeat_candidates(in_degrees::Dict{String, Int}, out_degrees::Dict{String, Int})
    candidates = String[]
    
    for vertex in keys(in_degrees)
        # High connectivity suggests repeat
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
function analyze_repeat_region(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                              start_vertex::String, min_length::Int)
    # Get local subgraph around the vertex
    local_vertices = get_local_subgraph(graph, start_vertex, min_length)
    
    if length(local_vertices) < min_length
        return nothing
    end
    
    # Analyze connectivity patterns
    incoming_edges = Tuple{String, String}[]
    outgoing_edges = Tuple{String, String}[]
    
    for vertex in local_vertices
        in_neighbors = get_in_neighbors(graph, vertex)
        out_neighbors = get_out_neighbors(graph, vertex)
        
        # Count external connections
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
    
    # Estimate copy number based on coverage
    copy_number = estimate_copy_number(graph, local_vertices)
    
    # Classify repeat type
    repeat_type = classify_repeat_type(graph, local_vertices, incoming_edges, outgoing_edges)
    
    # Calculate confidence
    confidence = calculate_repeat_confidence(graph, local_vertices, copy_number)
    
    return RepeatRegion(local_vertices, incoming_edges, outgoing_edges,
                       copy_number, repeat_type, confidence)
end

"""
Get local subgraph around a vertex.
"""
function get_local_subgraph(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                           center_vertex::String, radius::Int)
    visited = Set{String}()
    queue = DataStructures.Queue{Tuple{String, Int}}()
    DataStructures.enqueue!(queue, (center_vertex, 0))
    
    while !isempty(queue)
        vertex, distance = DataStructures.dequeue!(queue)
        
        if vertex in visited || distance > radius
            continue
        end
        
        push!(visited, vertex)
        
        # Add neighbors
        for neighbor in get_out_neighbors(graph, vertex)
            if !(neighbor in visited)
                DataStructures.enqueue!(queue, (neighbor, distance + 1))
            end
        end
        
        for neighbor in get_in_neighbors(graph, vertex)
            if !(neighbor in visited)
                DataStructures.enqueue!(queue, (neighbor, distance + 1))
            end
        end
    end
    
    return collect(visited)
end

"""
Estimate copy number for repeat region.
"""
function estimate_copy_number(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                             vertices::Vector{String})
    if isempty(vertices)
        return 1.0
    end
    
    # Use coverage information if available
    total_coverage = 0.0
    for vertex in vertices
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            coverage = length(vertex_data.source_positions)
            total_coverage += coverage
        end
    end
    
    avg_coverage = total_coverage / length(vertices)
    
    # Estimate based on coverage relative to expected single-copy coverage
    expected_single_copy = 10.0  # Could be estimated from graph statistics
    copy_number = max(1.0, avg_coverage / expected_single_copy)
    
    return copy_number
end

"""
Classify the type of repeat.
"""
function classify_repeat_type(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                             vertices::Vector{String},
                             incoming_edges::Vector{Tuple{String, String}},
                             outgoing_edges::Vector{Tuple{String, String}})
    # Simple classification based on connectivity pattern
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
function calculate_repeat_confidence(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                   vertices::Vector{String}, copy_number::Float64)
    # Higher copy number and larger region size increase confidence
    size_score = min(1.0, length(vertices) / 20.0)
    copy_score = min(1.0, (copy_number - 1.0) / 5.0)
    
    return (size_score + copy_score) / 2.0
end

"""
Merge overlapping repeat regions.
"""
function merge_overlapping_repeats(repeats::Vector{RepeatRegion})
    if isempty(repeats)
        return repeats
    end
    
    merged = RepeatRegion[]
    used = falses(length(repeats))
    
    for i in 1:length(repeats)
        if used[i]
            continue
        end
        
        current_repeat = repeats[i]
        overlapping_indices = [i]
        
        # Find overlapping repeats
        for j in (i+1):length(repeats)
            if !used[j] && regions_overlap(current_repeat, repeats[j])
                push!(overlapping_indices, j)
                used[j] = true
            end
        end
        
        # Merge overlapping regions
        if length(overlapping_indices) > 1
            merged_repeat = merge_repeat_regions([repeats[idx] for idx in overlapping_indices])
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
function regions_overlap(r1::RepeatRegion, r2::RepeatRegion)
    return !isempty(intersect(Set(r1.repeat_vertices), Set(r2.repeat_vertices)))
end

"""
Merge multiple repeat regions into one.
"""
function merge_repeat_regions(regions::Vector{RepeatRegion})
    all_vertices = String[]
    all_incoming = Tuple{String, String}[]
    all_outgoing = Tuple{String, String}[]
    
    for region in regions
        append!(all_vertices, region.repeat_vertices)
        append!(all_incoming, region.incoming_edges)
        append!(all_outgoing, region.outgoing_edges)
    end
    
    # Remove duplicates
    unique_vertices = unique(all_vertices)
    unique_incoming = unique(all_incoming)
    unique_outgoing = unique(all_outgoing)
    
    # Average copy number and confidence
    avg_copy_number = Statistics.mean(r.copy_number_estimate for r in regions)
    avg_confidence = Statistics.mean(r.confidence for r in regions)
    
    # Use most common repeat type
    repeat_types = [r.repeat_type for r in regions]
    repeat_type = mode(repeat_types)
    
    return RepeatRegion(unique_vertices, unique_incoming, unique_outgoing,
                       avg_copy_number, repeat_type, avg_confidence)
end

"""
    find_contigs_next(graph::MetaGraph, min_contig_length::Int=500) -> Vector{ContigPath}

Extract linear contigs from the assembly graph.
"""
function find_contigs_next(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData};
                          min_contig_length::Int=500)
    contigs = ContigPath[]
    visited = Set{String}()
    vertices = collect(MetaGraphsNext.labels(graph))
    
    for start_vertex in vertices
        if start_vertex in visited
            continue
        end
        
        # Find linear path starting from this vertex
        path = find_linear_path(graph, start_vertex, visited)
        
        if length(path) >= 2  # At least 2 k-mers for a meaningful contig
            # Generate sequence and coverage profile
            sequence = generate_contig_sequence(graph, path)
            coverage_profile = generate_coverage_profile(graph, path)
            
            if length(sequence) >= min_contig_length
                contig = ContigPath(path, sequence, coverage_profile)
                push!(contigs, contig)
            end
            
            # Mark vertices as visited
            for vertex in path
                push!(visited, vertex)
            end
        end
    end
    
    return sort_contigs_by_length(contigs)
end

"""
Find a linear path through the graph.
"""
function find_linear_path(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                         start_vertex::String, visited::Set{String})
    if start_vertex in visited
        return String[]
    end
    
    path = [start_vertex]
    current = start_vertex
    
    # Extend forward
    while true
        out_neighbors = get_out_neighbors(graph, current)
        valid_neighbors = [n for n in out_neighbors if !(n in visited) && !(n in path)]
        
        if length(valid_neighbors) == 1
            next_vertex = valid_neighbors[1]
            # Check if next vertex has only one incoming edge (to current)
            in_neighbors = get_in_neighbors(graph, next_vertex)
            if length(in_neighbors) == 1 && in_neighbors[1] == current
                push!(path, next_vertex)
                current = next_vertex
            else
                break
            end
        else
            break
        end
    end
    
    # Try to extend backward from start
    current = start_vertex
    backward_path = String[]
    
    while true
        in_neighbors = get_in_neighbors(graph, current)
        valid_neighbors = [n for n in in_neighbors if !(n in visited) && !(n in path) && !(n in backward_path)]
        
        if length(valid_neighbors) == 1
            prev_vertex = valid_neighbors[1]
            # Check if prev vertex has only one outgoing edge (to current)
            out_neighbors = get_out_neighbors(graph, prev_vertex)
            if length(out_neighbors) == 1 && out_neighbors[1] == current
                pushfirst!(backward_path, prev_vertex)
                current = prev_vertex
            else
                break
            end
        else
            break
        end
    end
    
    return [backward_path; path]
end

"""
Generate sequence for a contig path.
"""
function generate_contig_sequence(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                 path::Vector{String})
    if isempty(path)
        return ""
    end
    
    # Start with first k-mer
    sequence = path[1]
    
    # Add last character of each subsequent k-mer
    for i in 2:length(path)
        kmer = path[i]
        sequence *= kmer[end]  # Add last character
    end
    
    return sequence
end

"""
Generate coverage profile for a contig path.
"""
function generate_coverage_profile(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                  path::Vector{String})
    coverage = Float64[]
    
    for vertex in path
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            vertex_coverage = Float64(length(vertex_data.source_positions))
            push!(coverage, vertex_coverage)
        else
            push!(coverage, 0.0)
        end
    end
    
    return coverage
end

"""
Sort contigs by length (descending).
"""
function sort_contigs_by_length(contigs::Vector{ContigPath})
    return sort(contigs, by=c -> c.length, rev=true)
end

"""
    simplify_graph_next(graph::MetaGraph, bubbles::Vector{BubbleStructure}) -> MetaGraph

Simplify the graph by resolving bubbles and removing low-confidence paths.
"""
function simplify_graph_next(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                           bubbles::Vector{BubbleStructure})
    # Create a copy of the graph
    simplified_graph = deepcopy(graph)
    
    # Process bubbles by confidence
    sorted_bubbles = sort(bubbles, by=b -> b.complexity_score)
    
    for bubble in sorted_bubbles
        # Choose which path to keep based on support
        if bubble.path1_support > bubble.path2_support
            remove_path_from_graph!(simplified_graph, bubble.path2, bubble.entry_vertex, bubble.exit_vertex)
        elseif bubble.path2_support > bubble.path1_support
            remove_path_from_graph!(simplified_graph, bubble.path1, bubble.entry_vertex, bubble.exit_vertex)
        else
            # Equal support - could use other criteria or keep both
            continue
        end
    end
    
    # Remove isolated vertices
    remove_isolated_vertices!(simplified_graph)
    
    return simplified_graph
end

"""
Remove a path from the graph.
"""
function remove_path_from_graph!(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                path::Vector{String}, entry_vertex::String, exit_vertex::String)
    # Remove edges in the path
    prev_vertex = entry_vertex
    for vertex in path
        if haskey(graph, prev_vertex, vertex)
            delete!(graph, prev_vertex, vertex)
        end
        prev_vertex = vertex
    end
    
    # Remove final edge to exit
    if !isempty(path) && haskey(graph, path[end], exit_vertex)
        delete!(graph, path[end], exit_vertex)
    end
    
    # Remove vertices that are now isolated
    for vertex in path
        if is_isolated_vertex(graph, vertex)
            delete!(graph, vertex)
        end
    end
end

"""
Check if a vertex is isolated (no edges).
"""
function is_isolated_vertex(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData}, vertex::String)
    return isempty(get_in_neighbors(graph, vertex)) && isempty(get_out_neighbors(graph, vertex))
end

"""
Remove all isolated vertices from the graph.
"""
function remove_isolated_vertices!(graph::MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData})
    vertices_to_remove = String[]
    
    for vertex in MetaGraphsNext.labels(graph)
        if is_isolated_vertex(graph, vertex)
            push!(vertices_to_remove, vertex)
        end
    end
    
    for vertex in vertices_to_remove
        delete!(graph, vertex)
    end
end
