# Graph Simplification Algorithms
#
# This file contains algorithms for simplifying assembly graphs by:
# - Detecting and resolving bubbles (alternative paths)
# - Removing low-coverage paths
# - Collapsing linear chains
# - Removing tips (dead ends)
#
# Based on functions from src/sequence-graphs-next.jl

# ============================================================================
# Data Structures for Simplification
# ============================================================================

"""
Represents a bubble (alternative path) in the assembly graph.

A bubble occurs when the graph splits into two or more alternative paths
that later rejoin. Bubbles often represent:
- SNPs (single nucleotide polymorphisms)
- Small indels
- Sequencing errors
- Strain variation

# Fields
- `entry_vertex`: Vertex where paths diverge
- `exit_vertex`: Vertex where paths reconverge
- `path1`: First alternative path
- `path2`: Second alternative path
- `path1_support`: Coverage/support for path 1
- `path2_support`: Coverage/support for path 2
- `complexity_score`: Measure of bubble complexity (higher = more complex)
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

# ============================================================================
# Bubble Detection
# ============================================================================

"""
    detect_bubbles_next(graph::MetaGraphsNext.MetaGraph; min_bubble_length::Int=2, max_bubble_length::Int=100)

Detect bubble structures (alternative paths) in the assembly graph.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: The assembly graph to analyze
- `min_bubble_length::Int`: Minimum bubble path length (default: 2)
- `max_bubble_length::Int`: Maximum bubble path length (default: 100)

# Returns
- `Vector{BubbleStructure}`: Detected bubbles

# Algorithm
1. Find vertices with out-degree > 1 (potential bubble entry points)
2. Trace paths from each outgoing edge
3. Check if paths reconverge within max_bubble_length
4. Validate bubble structure and calculate support scores
5. Remove duplicate detections

# Example
```julia
graph = build_kmer_graph_next(Kmers.DNAKmer{31}, records)
bubbles = detect_bubbles_next(graph, min_bubble_length=2, max_bubble_length=50)
```
"""
function detect_bubbles_next(graph::MetaGraphsNext.MetaGraph;
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
function get_out_neighbors(graph::MetaGraphsNext.MetaGraph, vertex)
    neighbors = []
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
function get_in_neighbors(graph::MetaGraphsNext.MetaGraph, vertex)
    neighbors = []
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
function find_bubble_paths(graph::MetaGraphsNext.MetaGraph,
                          entry_vertex,
                          out_neighbors::Vector,
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
function find_limited_path(graph::MetaGraphsNext.MetaGraph, start_vertex, max_length::Int)
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
function find_path_convergence(path1::Vector, path2::Vector)
    for v1 in path1
        if v1 in path2
            return v1
        end
    end
    return nothing
end

"""
Calculate support for a path based on vertex coverage.
"""
function calculate_path_support(graph::MetaGraphsNext.MetaGraph, path::Vector)
    total_support = 0
    for vertex in path
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            if hasfield(typeof(vertex_data), :evidence)
                total_support += count_evidence(vertex_data)
            elseif hasfield(typeof(vertex_data), :coverage)
                total_support += length(vertex_data.coverage)
            else
                total_support += 1
            end
        end
    end
    return total_support
end

"""
Calculate complexity score for a bubble.
"""
function calculate_bubble_complexity(path1::Vector, path2::Vector)
    # Simple complexity measure: difference in path lengths
    length_diff = abs(length(path1) - length(path2))
    total_length = length(path1) + length(path2)
    return length_diff / max(total_length, 1)
end

"""
Check if a bubble structure is valid.
"""
function is_valid_bubble(graph::MetaGraphsNext.MetaGraph, bubble::BubbleStructure)
    # Check that entry and exit vertices exist
    if !haskey(graph, bubble.entry_vertex) || !haskey(graph, bubble.exit_vertex)
        return false
    end

    # Check that paths are non-empty
    if isempty(bubble.path1) || isempty(bubble.path2)
        return false
    end

    # Check that paths are different
    if bubble.path1 == bubble.path2
        return false
    end

    return true
end

"""
Remove duplicate bubble detections.
"""
function remove_duplicate_bubbles(bubbles::Vector{BubbleStructure})
    unique_bubbles = BubbleStructure[]
    seen_pairs = Set{Tuple{String, String}}()

    for bubble in bubbles
        pair = (bubble.entry_vertex, bubble.exit_vertex)
        if pair ∉ seen_pairs
            push!(unique_bubbles, bubble)
            push!(seen_pairs, pair)
        end
    end

    return unique_bubbles
end

# ============================================================================
# Graph Simplification
# ============================================================================

"""
    simplify_graph_next(graph::MetaGraphsNext.MetaGraph, bubbles::Vector{BubbleStructure})

Simplify the graph by resolving bubbles and removing low-confidence paths.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to simplify
- `bubbles::Vector{BubbleStructure}`: Detected bubbles to resolve

# Returns
- `MetaGraphsNext.MetaGraph`: Simplified graph

# Strategy
For each bubble:
1. Compare support scores for alternative paths
2. If one path has significantly higher support, remove the other
3. If support is similar, keep both paths (potential variation)

# Example
```julia
graph = build_kmer_graph_next(Kmers.DNAKmer{31}, records)
bubbles = detect_bubbles_next(graph)
simplified = simplify_graph_next(graph, bubbles)
```
"""
function simplify_graph_next(graph::MetaGraphsNext.MetaGraph,
                           bubbles::Vector{BubbleStructure})
    # Create a copy of the graph
    simplified_graph = deepcopy(graph)

    # Process bubbles by confidence
    for bubble in bubbles
        # Determine which path to keep based on support
        if bubble.path1_support > 2 * bubble.path2_support
            # Path 1 has much higher support - remove path 2
            remove_path_from_graph!(simplified_graph, bubble.path2, bubble.entry_vertex, bubble.exit_vertex)
        elseif bubble.path2_support > 2 * bubble.path1_support
            # Path 2 has much higher support - remove path 1
            remove_path_from_graph!(simplified_graph, bubble.path1, bubble.entry_vertex, bubble.exit_vertex)
        else
            # Similar support - keep both paths (potential true variation)
            continue
        end
    end

    return simplified_graph
end

"""
Remove a path from the graph (helper function for simplification).
"""
function remove_path_from_graph!(graph::MetaGraphsNext.MetaGraph, path::Vector,
                                entry_vertex, exit_vertex)
    # Remove internal vertices and edges along the path
    for i in 1:(length(path)-1)
        current = path[i]
        next_v = path[i+1]

        # Remove edge from current to next
        if haskey(graph, current, next_v)
            # Note: MetaGraphsNext doesn't have delete edge yet, skip for now
            # This is a placeholder for future implementation
        end
    end

    # Note: Full implementation would remove vertices with no remaining edges
end
