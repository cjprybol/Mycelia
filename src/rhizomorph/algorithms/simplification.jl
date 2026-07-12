# Graph Simplification Algorithms
#
# This file contains algorithms for simplifying assembly graphs by:
# - Detecting and resolving bubbles (alternative paths)
# - Removing low-coverage paths
# - Collapsing linear chains
# - Removing tips (dead ends)
#
# Based on legacy graph utilities now ported to Rhizomorph

# ============================================================================
# Data Structures for Simplification
# ============================================================================

"""
    BubbleStructure{T}

Represents a bubble (alternative path) in the assembly graph.

A bubble occurs when the graph splits into two or more alternative paths
that later rejoin. Bubbles often represent:
- SNPs (single nucleotide polymorphisms)
- Small indels
- Sequencing errors
- Strain variation

# Type Parameter
- `T`: The vertex label type (e.g., `String`, `Kmers.Kmer`, etc.)

# Fields
- `entry_vertex`: Vertex where paths diverge
- `exit_vertex`: Vertex where paths reconverge
- `path1`: First alternative path
- `path2`: Second alternative path
- `path1_support`: Coverage/support for path 1
- `path2_support`: Coverage/support for path 2
- `complexity_score`: Measure of bubble complexity (higher = more complex)

# Example
```julia
bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph)
```
"""
struct BubbleStructure{T}
    entry_vertex::T
    exit_vertex::T
    path1::Vector{T}
    path2::Vector{T}
    path1_support::Int
    path2_support::Int
    complexity_score::Float64

    function BubbleStructure(entry::T, exit::T,
            p1::Vector{T}, p2::Vector{T},
            s1::Int, s2::Int, complexity::Float64) where {T}
        new{T}(entry, exit, p1, p2, s1, s2, complexity)
    end
end

# ============================================================================
# Bubble Detection
# ============================================================================

"""
    _build_out_adjacency_index(graph::MetaGraphsNext.MetaGraph, ::Type{T}) where {T}

Build a vertex -> outgoing-neighbor-labels index in a single `O(E)` pass over the
graph's edge labels.

This replaces the per-vertex `get_out_neighbors` scan (each of which is itself
`O(E)`) inside bubble detection. Building the index once up front turns the
detection pass from `O(V*E)` into `O(V+E)`. For any given source vertex the
neighbor order matches what `get_out_neighbors` would return — both consume
`MetaGraphsNext.edge_labels(graph)` in the same order — so downstream bubble
results are identical to the pre-index scan.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to index
- `T`: Vertex label type

# Returns
- `Dict{T, Vector{T}}`: Map from each source vertex to its outgoing neighbors
"""
function _build_out_adjacency_index(graph::MetaGraphsNext.MetaGraph, ::Type{T}) where {T}
    out_index = Dict{T, Vector{T}}()
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        push!(get!(() -> T[], out_index, src), dst)
    end
    return out_index
end

"""
    _indexed_out_neighbors(out_index::AbstractDict{T, Vector{T}}, vertex) where {T}

Look up the outgoing neighbors of `vertex` in a prebuilt out-adjacency index in
`O(1)`. Returns an empty vector for vertices with no outgoing edges. Mirrors the
result of `get_out_neighbors(graph, vertex)` without rescanning the edge set.
"""
function _indexed_out_neighbors(out_index::AbstractDict{T, Vector{T}}, vertex) where {T}
    return haskey(out_index, vertex) ? out_index[vertex] : T[]
end

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
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31; dataset_id="dataset_01", mode=:singlestrand)
bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=2, max_bubble_length=50)
```
"""
function detect_bubbles_next(graph::MetaGraphsNext.MetaGraph;
        min_bubble_length::Int = 2,
        max_bubble_length::Int = 100)
    vertices = collect(MetaGraphsNext.labels(graph))
    T = eltype(vertices)
    bubbles = BubbleStructure{T}[]

    # Build the out-adjacency index ONCE in O(E), then reuse it for every
    # entry-point out-degree check and every forward path trace. This is the
    # fix for the former O(V*E) scaling: each vertex previously triggered a
    # full O(E) get_out_neighbors edge scan.
    out_index = _build_out_adjacency_index(graph, T)

    for entry_vertex in vertices
        # Find potential bubble entry points (vertices with out-degree > 1)
        out_neighbors = _indexed_out_neighbors(out_index, entry_vertex)

        if length(out_neighbors) >= 2
            # Look for bubbles starting from this vertex
            bubble_candidates = find_bubble_paths(graph, entry_vertex, out_neighbors,
                min_bubble_length, max_bubble_length; out_index = out_index)

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
    get_out_neighbors(graph::MetaGraphsNext.MetaGraph, vertex)

Return outgoing neighbor labels for `vertex`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to query
- `vertex`: Vertex label

# Returns
- `Vector`: Outgoing neighbor labels

# Example
```julia
neighbors = Mycelia.Rhizomorph.get_out_neighbors(graph, vertex)
```
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
    get_in_neighbors(graph::MetaGraphsNext.MetaGraph, vertex)

Return incoming neighbor labels for `vertex`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to query
- `vertex`: Vertex label

# Returns
- `Vector`: Incoming neighbor labels

# Example
```julia
neighbors = Mycelia.Rhizomorph.get_in_neighbors(graph, vertex)
```
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
    find_bubble_paths(graph::MetaGraphsNext.MetaGraph, entry_vertex, out_neighbors, min_length, max_length)

Enumerate candidate bubbles that originate at `entry_vertex`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to analyze
- `entry_vertex`: Candidate bubble entry point
- `out_neighbors::Vector`: Outgoing neighbors from `entry_vertex`
- `min_length::Int`: Minimum path length required for a bubble
- `max_length::Int`: Maximum path length to trace before giving up

# Returns
- `Vector{BubbleStructure}`: Candidate bubble structures

# Example
```julia
candidates = Mycelia.Rhizomorph.find_bubble_paths(graph, vertex, neighbors, 2, 50)
```
"""
function find_bubble_paths(graph::MetaGraphsNext.MetaGraph,
        entry_vertex::T,
        out_neighbors::Vector,
        min_length::Int, max_length::Int;
        out_index::AbstractDict = _build_out_adjacency_index(graph, T)) where {T}
    bubbles = BubbleStructure{T}[]

    # Try all pairs of outgoing paths
    for i in 1:length(out_neighbors)
        for j in (i + 1):length(out_neighbors)
            path1_start = out_neighbors[i]
            path2_start = out_neighbors[j]

            # Find paths from each starting point (reuse the shared O(E) index)
            path1 = find_limited_path(graph, path1_start, max_length; out_index = out_index)
            path2 = find_limited_path(graph, path2_start, max_length; out_index = out_index)

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
    find_limited_path(graph::MetaGraphsNext.MetaGraph, start_vertex, max_length::Int)

Trace a simple path forward from `start_vertex` while the out-degree stays `1`.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to traverse
- `start_vertex`: Seed vertex label
- `max_length::Int`: Maximum number of extension steps

# Returns
- `Vector`: Limited path beginning at `start_vertex`

# Example
```julia
path = Mycelia.Rhizomorph.find_limited_path(graph, vertex, 25)
```
"""
function find_limited_path(graph::MetaGraphsNext.MetaGraph, start_vertex, max_length::Int;
        out_index::AbstractDict = _build_out_adjacency_index(
            graph, eltype(collect(MetaGraphsNext.labels(graph)))))
    path = [start_vertex]
    current = start_vertex

    for _ in 1:max_length
        neighbors = _indexed_out_neighbors(out_index, current)
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
    find_path_convergence(path1::Vector, path2::Vector)

Return the first shared vertex between `path1` and `path2`.

# Arguments
- `path1::Vector`: First candidate path
- `path2::Vector`: Second candidate path

# Returns
- Shared vertex label or `nothing` if the paths never reconverge

# Example
```julia
vertex = Mycelia.Rhizomorph.find_path_convergence(path1, path2)
```
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
    calculate_path_support(graph::MetaGraphsNext.MetaGraph, path::Vector)

Estimate support for `path` by summing per-vertex evidence.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph containing evidence-bearing vertices
- `path::Vector`: Vertex labels along the path

# Returns
- `Int`: Total support score

# Example
```julia
support = Mycelia.Rhizomorph.calculate_path_support(graph, path)
```
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
    calculate_bubble_complexity(path1::Vector, path2::Vector)

Compute a simple length-based complexity score for a bubble.

# Arguments
- `path1::Vector`: First bubble branch
- `path2::Vector`: Second bubble branch

# Returns
- `Float64`: Normalized path-length difference in `[0, 1]`

# Example
```julia
complexity = Mycelia.Rhizomorph.calculate_bubble_complexity(path1, path2)
```
"""
function calculate_bubble_complexity(path1::Vector, path2::Vector)
    # Simple complexity measure: difference in path lengths
    length_diff = abs(length(path1) - length(path2))
    total_length = length(path1) + length(path2)
    return length_diff / max(total_length, 1)
end

"""
    is_valid_bubble(graph::MetaGraphsNext.MetaGraph, bubble::BubbleStructure{T}) where T

Validate that a bubble references existing vertices and distinct branches.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph that produced `bubble`
- `bubble::BubbleStructure{T}`: Candidate bubble

# Returns
- `Bool`: `true` if the bubble is structurally valid

# Example
```julia
valid = Mycelia.Rhizomorph.is_valid_bubble(graph, bubble)
```
"""
function is_valid_bubble(graph::MetaGraphsNext.MetaGraph, bubble::BubbleStructure{T}) where {T}
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
    remove_duplicate_bubbles(bubbles::Vector{BubbleStructure{T}}) where T

Drop duplicate bubble calls with the same entry and exit vertices.

# Arguments
- `bubbles::Vector{BubbleStructure{T}}`: Candidate bubble calls

# Returns
- `Vector{BubbleStructure{T}}`: Deduplicated bubble collection

# Example
```julia
unique_bubbles = Mycelia.Rhizomorph.remove_duplicate_bubbles(bubbles)
```
"""
function remove_duplicate_bubbles(bubbles::Vector{BubbleStructure{T}}) where {T}
    unique_bubbles = BubbleStructure{T}[]
    seen_pairs = Set{Tuple{T, T}}()

    for bubble in bubbles
        pair = (bubble.entry_vertex, bubble.exit_vertex)
        if pair ∉ seen_pairs
            push!(unique_bubbles, bubble)
            push!(seen_pairs, pair)
        end
    end

    return unique_bubbles
end

"""
    enumerate_superbubble_paths(graph, bubble; max_paths=10, max_depth=100, max_expansions=10000)

Enumerate ranked alternatives across a detected superbubble-like local region.

This wrapper uses `bubble.entry_vertex` and `bubble.exit_vertex` as local
boundaries, then delegates to `enumerate_local_paths`. It can return more than
the two branch paths stored in `BubbleStructure`, which makes it useful for
multi-branch superbubbles and repeat-adjacent local regions.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Rhizomorph or weighted graph
- `bubble::BubbleStructure`: Detected local bubble boundary

# Keywords
- `max_paths::Int=10`: Maximum number of alternatives to return
- `max_depth::Int=100`: Maximum number of edges per path
- `max_expansions::Int=10000`: Maximum partial paths to expand
- `min_probability::Float64=0.0`: Prune paths below this cumulative probability
- `default_weight::Float64=1e-10`: Fallback edge weight during conversion
- `edge_weight::Function=count_evidence`: Edge scoring callback during conversion

# Returns
- `LocalPathEnumerationResult`: Ranked alternatives plus provenance metadata

# Example
```julia
bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph)
result = Mycelia.Rhizomorph.enumerate_superbubble_paths(graph, first(bubbles))
```
"""
function enumerate_superbubble_paths(
        graph::MetaGraphsNext.MetaGraph,
        bubble::BubbleStructure{T};
        max_paths::Int = 10,
        max_depth::Int = 100,
        max_expansions::Int = 10_000,
        min_probability::Float64 = 0.0,
        default_weight::Float64 = 1e-10,
        edge_weight::Function = count_evidence
) where {T}
    return _enumerate_local_paths_impl(
        graph,
        bubble.entry_vertex,
        bubble.exit_vertex,
        :superbubble;
        max_paths = max_paths,
        max_depth = max_depth,
        max_expansions = max_expansions,
        min_probability = min_probability,
        default_weight = default_weight,
        edge_weight = edge_weight
    )
end

"""
    enumerate_superbubble_paths(graph; min_bubble_length=2, max_bubble_length=100, max_bubbles=typemax(Int), ...)

Detect local bubble boundaries and enumerate ranked alternatives for each one.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Rhizomorph or weighted graph

# Keywords
- `min_bubble_length::Int=2`: Minimum branch length passed to bubble detection
- `max_bubble_length::Int=100`: Maximum branch trace length passed to detection
- `max_bubbles::Int=typemax(Int)`: Maximum detected regions to enumerate
- `max_paths::Int=10`: Maximum alternatives per region
- `max_depth::Int=100`: Maximum edge depth per local enumeration
- `max_expansions::Int=10000`: Maximum partial-path expansions per region
- `min_probability::Float64=0.0`: Prune paths below this cumulative probability
- `default_weight::Float64=1e-10`: Fallback edge weight during conversion
- `edge_weight::Function=count_evidence`: Edge scoring callback during conversion

# Returns
- `Vector{LocalPathEnumerationResult}`: One enumeration result per detected region

# Notes
Use `max_bubbles`, `max_paths`, `max_depth`, and `max_expansions` together on
highly repetitive graphs. Each result records whether a guard truncated local
search in its provenance.
"""
function enumerate_superbubble_paths(
        graph::MetaGraphsNext.MetaGraph;
        min_bubble_length::Int = 2,
        max_bubble_length::Int = 100,
        max_bubbles::Int = typemax(Int),
        max_paths::Int = 10,
        max_depth::Int = 100,
        max_expansions::Int = 10_000,
        min_probability::Float64 = 0.0,
        default_weight::Float64 = 1e-10,
        edge_weight::Function = count_evidence
)
    if max_bubbles < 0
        throw(ArgumentError("max_bubbles must be non-negative, got $max_bubbles"))
    end

    labels = collect(MetaGraphsNext.labels(graph))
    results = LocalPathEnumerationResult{eltype(labels)}[]
    if max_bubbles == 0 || isempty(labels)
        return results
    end

    bubbles = detect_bubbles_next(
        graph; min_bubble_length = min_bubble_length, max_bubble_length = max_bubble_length)

    for bubble in bubbles
        if length(results) >= max_bubbles
            break
        end
        push!(
            results,
            enumerate_superbubble_paths(
                graph,
                bubble;
                max_paths = max_paths,
                max_depth = max_depth,
                max_expansions = max_expansions,
                min_probability = min_probability,
                default_weight = default_weight,
                edge_weight = edge_weight
            )
        )
    end

    return results
end

# ============================================================================
# Tip removal (dead-end pruning)
# ============================================================================

"""
    remove_tips!(graph::MetaGraphsNext.MetaGraph; min_support::Int=1)

Remove dead-end vertices (tips) whose evidence support is ≤ `min_support`.
Tips are vertices with indegree == 0 or outdegree == 0. Returns the modified graph.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to prune in place
- `min_support::Int=1`: Maximum evidence count allowed for removable tips

# Returns
- `MetaGraphsNext.MetaGraph`: The modified input graph

# Example
```julia
Mycelia.Rhizomorph.remove_tips!(graph; min_support=2)
```
"""
function remove_tips!(graph::MetaGraphsNext.MetaGraph; min_support::Int = 1)
    if Graphs.nv(graph.graph) == 0
        return graph
    end

    labels = collect(MetaGraphsNext.labels(graph))
    to_delete = Set{eltype(labels)}()

    for label in labels
        idx = MetaGraphsNext.code_for(graph, label)
        indeg = Graphs.indegree(graph.graph, idx)
        outdeg = Graphs.outdegree(graph.graph, idx)
        if indeg == 0 || outdeg == 0
            support = 0
            data = graph[label]
            if hasproperty(data, :evidence)
                support = Rhizomorph.count_evidence(data)
            end
            if support <= min_support
                push!(to_delete, label)
            end
        end
    end

    for label in to_delete
        idx = MetaGraphsNext.code_for(graph, label)
        MetaGraphsNext.rem_vertex!(graph, idx)
    end

    return graph
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
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31; dataset_id="dataset_01", mode=:singlestrand)
bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph)
simplified = Mycelia.Rhizomorph.simplify_graph_next(graph, bubbles)
```
"""
function simplify_graph_next(graph::MetaGraphsNext.MetaGraph,
        bubbles::Vector{BubbleStructure{T}}) where {T}
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

    remove_isolated_vertices!(simplified_graph)

    return simplified_graph
end

"""
    remove_path_from_graph!(graph::MetaGraphsNext.MetaGraph, path::Vector, entry_vertex, exit_vertex)

Delete a bubble branch from `graph` and remove any vertices that become isolated.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to mutate
- `path::Vector`: Interior branch vertices to remove
- `entry_vertex`: Bubble entry vertex
- `exit_vertex`: Bubble exit vertex

# Returns
- `Nothing`: Mutates `graph` in place

# Example
```julia
Mycelia.Rhizomorph.remove_path_from_graph!(graph, bubble.path2, bubble.entry_vertex, bubble.exit_vertex)
```
"""
function remove_path_from_graph!(graph::MetaGraphsNext.MetaGraph, path::Vector,
        entry_vertex, exit_vertex)
    prev_vertex = entry_vertex

    for vertex in path
        if haskey(graph, prev_vertex, vertex)
            delete!(graph, prev_vertex, vertex)
        end
        prev_vertex = vertex
    end

    if !isempty(path) && haskey(graph, path[end], exit_vertex)
        delete!(graph, path[end], exit_vertex)
    end

    for vertex in path
        if is_isolated_vertex(graph, vertex)
            delete!(graph, vertex)
        end
    end
end

"""
    is_isolated_vertex(graph::MetaGraphsNext.MetaGraph, vertex)

Check whether `vertex` has neither incoming nor outgoing edges.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to query
- `vertex`: Vertex label

# Returns
- `Bool`: `true` if `vertex` is isolated

# Example
```julia
isolated = Mycelia.Rhizomorph.is_isolated_vertex(graph, vertex)
```
"""
function is_isolated_vertex(graph::MetaGraphsNext.MetaGraph, vertex)
    return isempty(get_in_neighbors(graph, vertex)) &&
           isempty(get_out_neighbors(graph, vertex))
end

"""
    remove_isolated_vertices!(graph::MetaGraphsNext.MetaGraph)

Remove every isolated vertex from `graph` in place.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to prune

# Returns
- `MetaGraphsNext.MetaGraph`: The modified input graph

# Example
```julia
Mycelia.Rhizomorph.remove_isolated_vertices!(graph)
```
"""
function remove_isolated_vertices!(graph::MetaGraphsNext.MetaGraph)
    vertices_to_remove = Vector{eltype(MetaGraphsNext.labels(graph))}()

    for vertex in MetaGraphsNext.labels(graph)
        if is_isolated_vertex(graph, vertex)
            push!(vertices_to_remove, vertex)
        end
    end

    for vertex in vertices_to_remove
        delete!(graph, vertex)
    end

    return graph
end

# ============================================================================
# Linear chain collapsing
# ============================================================================

"""
    collapse_linear_chains!(graph::MetaGraphsNext.MetaGraph)

Collapse maximal linear unitigs (vertices with indegree=1 and outdegree=1) into
single variable-length vertices while preserving aggregated evidence.
Currently operates on graphs whose vertex data include a `sequence` field
(`BioSequenceVertexData`/`QualityBioSequenceVertexData`, their
lightweight/ultralight variants, and `StringVertexData` plus the lightweight
and ultralight string variants).

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Graph to simplify in place

# Returns
- `MetaGraphsNext.MetaGraph`: The modified input graph

# Example
```julia
Mycelia.Rhizomorph.collapse_linear_chains!(graph)
```
"""
function collapse_linear_chains!(graph::MetaGraphsNext.MetaGraph)
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return graph
    end

    visited = Set{eltype(labels)}()

    for label in labels
        if label in visited
            continue
        end

        path = Rhizomorph.get_maximal_unitig(graph, label)
        for v in path
            push!(visited, v)
        end

        if length(path) < 2
            continue
        end

        # Cyclic unitigs can revisit the same label (e.g. 2-cycles in string graphs).
        # Those are not linear chains and should be left untouched.
        if length(unique(path)) < length(path)
            continue
        end

        first_vertex_data = graph[path[1]]
        if hasfield(typeof(first_vertex_data), :Kmer)
            # K-mer graphs cannot be collapsed in-place: the assembled sequence
            # is longer than k, so the label type changes from Kmer to BioSequence,
            # violating MetaGraphsNext's parametric label_type constraint.
            # Use convert_fixed_to_variable() first, then collapse.
            continue
        end

        new_sequence, offsets = _assemble_linear_chain_sequence(graph, path)
        collapsed_vertex = _build_collapsed_vertex(first_vertex_data, new_sequence)
        _merge_path_evidence!(collapsed_vertex, graph, path, offsets)

        # Capture external edges before removal
        edge_copies = Dict{Tuple{eltype(path), eltype(path)}, Any}()
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            in_src = haskey(offsets, src)
            in_dst = haskey(offsets, dst)
            if xor(in_src, in_dst)
                edge_copies[(src, dst)] = graph[src, dst]
            end
        end

        # Remove internal vertices (and their edges)
        # Must use code_for() since rem_vertex! takes integer codes, not labels.
        # Re-lookup each time because removal reassigns vertex codes.
        for v in path
            idx = MetaGraphsNext.code_for(graph, v)
            MetaGraphsNext.rem_vertex!(graph, idx)
        end

        new_label = new_sequence
        graph[new_label] = collapsed_vertex

        # Reconnect external edges with position-adjusted evidence
        for ((src, dst), edge_data) in edge_copies
            if haskey(offsets, src) && !haskey(offsets, dst)
                shifted = _shift_edge_data(edge_data, offsets[src]; shift_from = true)
                graph[new_label, dst] = shifted
            elseif !haskey(offsets, src) && haskey(offsets, dst)
                shifted = _shift_edge_data(edge_data, offsets[dst]; shift_to = true)
                graph[src, new_label] = shifted
            end
        end
    end

    return graph
end

"""
    _get_vertex_content(vertex_data)

Extract the core content (sequence or string) from any vertex data type.
Works for BioSequence types (`.sequence` field) and String types (`.string_value` field).
"""
function _get_vertex_content(vertex_data)
    if hasfield(typeof(vertex_data), :sequence)
        return vertex_data.sequence
    elseif hasfield(typeof(vertex_data), :string_value)
        return vertex_data.string_value
    else
        error("Cannot extract content from vertex data type: $(typeof(vertex_data))")
    end
end

function _assemble_linear_chain_sequence(graph, path::Vector)
    first_data = graph[path[1]]
    sequence = _get_vertex_content(first_data)
    offsets = Dict{eltype(path), Int}()
    offsets[path[1]] = 0

    for i in 2:length(path)
        src = path[i - 1]
        dst = path[i]
        overlap = _edge_overlap_length(graph, src, dst)
        append_sequence = _get_vertex_content(graph[dst])

        offset = offsets[src] + length(_get_vertex_content(graph[src])) - overlap
        offsets[dst] = offset

        if overlap < length(append_sequence)
            sequence = sequence * _suffix_after_overlap(append_sequence, overlap)
        end
    end

    return sequence, offsets
end

function _suffix_after_overlap(sequence::AbstractString, overlap::Int)
    overlap <= 0 && return sequence

    current_index = firstindex(sequence)
    for _ in 1:overlap
        current_index > lastindex(sequence) && return ""
        current_index = nextind(sequence, current_index)
    end

    return current_index > lastindex(sequence) ? "" : sequence[current_index:end]
end

function _suffix_after_overlap(sequence, overlap::Int)
    return overlap < length(sequence) ? sequence[(overlap + 1):end] : typeof(sequence)()
end

function _build_collapsed_vertex(first_vertex_data, sequence)
    # TYPE-CHECK-AUDIT: factory dispatch — add branch for each new vertex type
    # BioSequence vertex types
    if first_vertex_data isa BioSequenceVertexData
        return BioSequenceVertexData(sequence)
    elseif first_vertex_data isa QualityBioSequenceVertexData
        return QualityBioSequenceVertexData(sequence)
    elseif first_vertex_data isa LightweightBioSequenceVertexData
        return LightweightBioSequenceVertexData(sequence)
    elseif first_vertex_data isa UltralightBioSequenceVertexData
        return UltralightBioSequenceVertexData(sequence)
    elseif first_vertex_data isa UltralightQualityBioSequenceVertexData
        return UltralightQualityBioSequenceVertexData(sequence)
    elseif first_vertex_data isa LightweightQualityBioSequenceVertexData
        return LightweightQualityBioSequenceVertexData(sequence)
        # String vertex types
    elseif first_vertex_data isa StringVertexData
        return StringVertexData(sequence)
    elseif first_vertex_data isa LightweightStringVertexData
        return LightweightStringVertexData(sequence)
    elseif first_vertex_data isa UltralightStringVertexData
        return UltralightStringVertexData(sequence)
    else
        error("Unsupported vertex data type for collapsing: $(typeof(first_vertex_data))")
    end
end

function _merge_path_evidence!(target_vertex, graph, path, offsets)
    for vertex in path
        vertex_data = graph[vertex]
        offset = offsets[vertex]

        # TYPE-CHECK-AUDIT: super-union guard — must cover all reduced vertex types
        if vertex_data isa AllReducedVertexData
            target_vertex.total_count += vertex_data.total_count
            for (ds_id, ds_count) in vertex_data.dataset_counts
                target_vertex.dataset_counts[ds_id] = get(
                    target_vertex.dataset_counts, ds_id, 0) + ds_count
            end
            if hasfield(typeof(vertex_data), :dataset_observations)
                for (ds_id, obs_set) in vertex_data.dataset_observations
                    if !haskey(target_vertex.dataset_observations, ds_id)
                        target_vertex.dataset_observations[ds_id] = Set{String}()
                    end
                    union!(target_vertex.dataset_observations[ds_id], obs_set)
                end
            end
            if hasfield(typeof(vertex_data), :joint_quality) &&
               hasfield(typeof(target_vertex), :joint_quality)
                if isempty(target_vertex.joint_quality)
                    append!(target_vertex.joint_quality, vertex_data.joint_quality)
                elseif !isempty(vertex_data.joint_quality)
                    for i in eachindex(vertex_data.joint_quality)
                        target_vertex.joint_quality[i] = UInt8(clamp(
                            Int(target_vertex.joint_quality[i]) +
                            Int(vertex_data.joint_quality[i]), 0, 255))
                    end
                end
            end
        else
            for (dataset_id, dataset_evidence) in vertex_data.evidence
                for (obs_id, evidence_set) in dataset_evidence
                    for entry in evidence_set
                        shifted_entry = _shift_evidence_entry(entry, offset)
                        add_evidence!(target_vertex, dataset_id, obs_id, shifted_entry)
                    end
                end
            end
        end
    end
end

function _shift_evidence_entry(entry::EvidenceEntry, offset::Int)
    return EvidenceEntry(entry.position + offset, entry.strand)
end

function _shift_evidence_entry(entry::QualityEvidenceEntry, offset::Int)
    return QualityEvidenceEntry(entry.position + offset, entry.strand, entry.quality_scores)
end

function _shift_edge_data(edge_data, offset::Int; shift_from::Bool = false, shift_to::Bool = false)
    # TYPE-CHECK-AUDIT: super-union guard — must cover all reduced edge types
    if edge_data isa AllReducedEdgeData
        new_edge = if hasfield(typeof(edge_data), :overlap_length)
            typeof(edge_data)(edge_data.overlap_length)
        else
            typeof(edge_data)()
        end
        new_edge.total_count = edge_data.total_count
        for (ds_id, ds_count) in edge_data.dataset_counts
            new_edge.dataset_counts[ds_id] = ds_count
        end
        if hasfield(typeof(edge_data), :dataset_observations)
            for (ds_id, obs_set) in edge_data.dataset_observations
                new_edge.dataset_observations[ds_id] = copy(obs_set)
            end
        end
        if hasfield(typeof(edge_data), :from_joint_quality)
            append!(new_edge.from_joint_quality, edge_data.from_joint_quality)
            append!(new_edge.to_joint_quality, edge_data.to_joint_quality)
        end
        return new_edge
    end

    if hasfield(typeof(edge_data), :overlap_length)
        new_edge = typeof(edge_data)(edge_data.overlap_length)
    else
        new_edge = typeof(edge_data)()
    end

    for (dataset_id, dataset_evidence) in edge_data.evidence
        for (obs_id, evidence_set) in dataset_evidence
            for entry in evidence_set
                shifted_entry = _shift_edge_entry(entry, offset; shift_from = shift_from, shift_to = shift_to)
                add_evidence!(new_edge, dataset_id, obs_id, shifted_entry)
            end
        end
    end

    return new_edge
end

function _shift_edge_entry(entry::EdgeEvidenceEntry, offset::Int; shift_from::Bool, shift_to::Bool)
    from_pos = shift_from ? entry.from_position + offset : entry.from_position
    to_pos = shift_to ? entry.to_position + offset : entry.to_position
    return EdgeEvidenceEntry(from_pos, to_pos, entry.strand)
end

function _shift_edge_entry(entry::EdgeQualityEvidenceEntry, offset::Int; shift_from::Bool, shift_to::Bool)
    from_pos = shift_from ? entry.from_position + offset : entry.from_position
    to_pos = shift_to ? entry.to_position + offset : entry.to_position
    return EdgeQualityEvidenceEntry(from_pos, to_pos, entry.strand,
        entry.from_quality_scores,
        entry.to_quality_scores)
end

# ============================================================================
# Linear-time corrector-graph defragmentation cleanup (td-969e)
# ============================================================================
#
# The scalable corrector's final graph retains error-induced BRANCH POINTS
# (dead-end tips + low-coverage bubbles). find_contigs_next breaks a unitig at
# every branch vertex, so each residual branch point becomes a contig boundary
# (~6 contigs for a 1kb genome after RC-dedup). This module removes ONLY the
# unambiguous errors before contig extraction, in O(V+E), while never collapsing
# a data-supported variant (the td-h6w9 variation-preservation invariant).
#
# Two safe passes:
#   1. clip_error_tips!        — coverage-1 dead-end branches dangling off a
#                                junction (never a stand-alone run / genome
#                                terminus, never a rejoining variant).
#   2. collapse_error_bubbles! — a bubble branch collapsed ONLY when its mean
#                                coverage is ~1 (error) AND the sibling is well
#                                supported. Both-supported bubbles (balanced or
#                                skewed) retain BOTH branches.

"""
    _build_in_adjacency_index(graph::MetaGraphsNext.MetaGraph, ::Type{T}) where {T}

Build a vertex -> incoming-neighbor-labels index in a single `O(E)` pass over the
graph's edge labels. The incoming-edge counterpart to
[`_build_out_adjacency_index`](@ref); together they give `O(1)` in/out-degree and
neighbor lookups so tip/bubble cleanup stays `O(V+E)`.
"""
function _build_in_adjacency_index(graph::MetaGraphsNext.MetaGraph, ::Type{T}) where {T}
    in_index = Dict{T, Vector{T}}()
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        push!(get!(() -> T[], in_index, dst), src)
    end
    return in_index
end

"""
    _vertex_support(graph::MetaGraphsNext.MetaGraph, vertex) -> Int

Coverage support (number of evidence entries) on `vertex`, or `0` when the vertex
data carries no `evidence` field. This is the per-vertex coverage used to
separate coverage-1 errors from data-supported sequence.
"""
function _vertex_support(graph::MetaGraphsNext.MetaGraph, vertex)
    haskey(graph, vertex) || return 0
    data = graph[vertex]
    return hasproperty(data, :evidence) ? count_evidence(data) : 0
end

"""
    _collect_dangling_tip(deadend, in_index, out_index, indeg, outdeg, max_len, ::Type{T}; backward) where {T}

Walk inward from a dead-end vertex, collecting the maximal linear chain until it
either reaches a branch JUNCTION (the vertex the tip dangles off, whose relevant
degree is > 1) or terminates without one.

`backward=true` walks a SINK tip (dead end with no out-edges) back toward its
in-neighbors; `backward=false` walks a SOURCE tip (dead end with no in-edges)
forward. Returns `(vertices, reached_junction)`. The junction vertex itself is
NOT included in `vertices`. A chain that reaches the other terminus (a stand-
alone linear run / genome end) without a junction returns `reached_junction=false`
and is therefore NEVER clipped — the variation/terminus-preservation guard.
"""
function _collect_dangling_tip(deadend, in_index, out_index, indeg, outdeg,
        max_len::Int, ::Type{T}; backward::Bool) where {T}
    tip = T[deadend]
    cur = deadend
    reached_junction = false
    step_index = backward ? in_index : out_index
    # The junction is detected on the vertex's degree on the SIDE we walk toward:
    # a sink tip dangles off a vertex with out-degree > 1; a source tip off a
    # vertex with in-degree > 1.
    branch_degree = backward ? outdeg : indeg
    while length(tip) <= max_len
        steps = get(step_index, cur, T[])
        length(steps) == 1 || break            # 0 neighbors => other terminus, no junction
        nxt = steps[1]
        if branch_degree(nxt) > 1
            reached_junction = true             # dangling tip off a real branch junction
            break
        end
        # Continue only through purely linear interior (indeg==1 && outdeg==1).
        if indeg(nxt) == 1 && outdeg(nxt) == 1
            push!(tip, nxt)
            cur = nxt
        else
            break                               # convergence/other structure: not a clean tip
        end
    end
    return (vertices = tip, reached_junction = reached_junction)
end

"""
    clip_error_tips!(graph; max_tip_length, max_tip_support=1, max_rounds=10) -> (; removed, rounds)

Remove short, low-coverage dead-end branches ("tips") that dangle off a branch
vertex. Runs in `O(V+E)` per round using prebuilt in/out adjacency indices.

A tip is clipped only when BOTH unambiguous-error conditions hold:
  * every vertex on the tip has coverage `<= max_tip_support` (coverage-1 error
    support), and
  * the tip has at most `max_tip_length` vertices.

A maximal dead-end chain that reaches the OTHER dead end (a stand-alone linear
run or genome terminus) WITHOUT passing a junction is NEVER clipped. Real variants
rejoin (bubbles) and are never dead ends, so tip clipping structurally cannot
collapse a real variant — this is the safe subset of the cleanup.

Iterates up to `max_rounds` times because removing a tip can turn a former
junction into a linear vertex and expose a nested tip one layer up.

# Returns
- `(; removed::Int, rounds::Int)`: total vertices removed and rounds executed.
"""
function clip_error_tips!(graph::MetaGraphsNext.MetaGraph;
        max_tip_length::Int,
        max_tip_support::Real = 1,
        max_rounds::Int = 10)
    total_removed = 0
    rounds = 0
    for _ in 1:max_rounds
        labels = collect(MetaGraphsNext.labels(graph))
        isempty(labels) && break
        T = eltype(labels)
        out_index = _build_out_adjacency_index(graph, T)
        in_index = _build_in_adjacency_index(graph, T)
        outdeg = v -> length(get(out_index, v, T[]))
        indeg = v -> length(get(in_index, v, T[]))

        to_remove = Set{T}()
        visited = Set{T}()

        for v in labels
            (v in visited || v in to_remove) && continue
            od = outdeg(v)
            id = indeg(v)
            tip = if od == 0 && id >= 1
                _collect_dangling_tip(v, in_index, out_index, indeg, outdeg,
                    max_tip_length, T; backward = true)
            elseif id == 0 && od >= 1
                _collect_dangling_tip(v, in_index, out_index, indeg, outdeg,
                    max_tip_length, T; backward = false)
            else
                nothing
            end
            tip === nothing && continue
            for u in tip.vertices
                push!(visited, u)
            end
            if tip.reached_junction &&
               all(u -> _vertex_support(graph, u) <= max_tip_support, tip.vertices)
                for u in tip.vertices
                    push!(to_remove, u)
                end
            end
        end

        isempty(to_remove) && break
        for label in to_remove
            haskey(graph, label) && delete!(graph, label)
        end
        total_removed += length(to_remove)
        rounds += 1
    end
    return (; removed = total_removed, rounds = rounds)
end

"""
    _branch_mean_coverage(graph, path::Vector) -> Union{Float64, Nothing}

Mean per-vertex coverage across a bubble branch, excluding the shared
convergence vertex (`path[end]`). Returns `nothing` for an empty branch. Used to
separate a coverage-1 error branch from a data-supported (real) allele.
"""
function _branch_mean_coverage(graph::MetaGraphsNext.MetaGraph, path::Vector)
    isempty(path) && return nothing
    interior = length(path) > 1 ? @view(path[1:(end - 1)]) : path
    total = 0.0
    n = 0
    for v in interior
        haskey(graph, v) || continue
        total += _vertex_support(graph, v)
        n += 1
    end
    return n == 0 ? nothing : total / n
end

"""
    collapse_error_bubbles!(graph; max_error_support=2, min_real_support=3, ...) -> (; collapsed)

Collapse a bubble branch ONLY when it is an unambiguous error: its mean per-vertex
coverage `<= max_error_support` (≈ coverage-1) AND the sibling branch's mean
per-vertex coverage `>= min_real_support`. When both branches are data-supported —
a balanced (e.g. 15x/15x) or skewed-but-supported (e.g. 20x/10x) bubble — BOTH
branches are retained (neither branch coverage is `<= max_error_support`), so real
variation is preserved (the td-h6w9 invariant).

Bubble detection is the `O(V+E)` `detect_bubbles_next`; each error branch is
removed with [`remove_path_from_graph!`](@ref).

# Returns
- `(; collapsed::Int)`: number of error branches removed.
"""
function collapse_error_bubbles!(graph::MetaGraphsNext.MetaGraph;
        max_error_support::Real = 2,
        min_real_support::Real = 3,
        min_bubble_length::Int = 2,
        max_bubble_length::Int = 100)
    isempty(collect(MetaGraphsNext.labels(graph))) && return (; collapsed = 0)
    bubbles = detect_bubbles_next(graph;
        min_bubble_length = min_bubble_length, max_bubble_length = max_bubble_length)
    collapsed = 0
    for bubble in bubbles
        # A prior collapse may have already removed this bubble's boundary.
        (haskey(graph, bubble.entry_vertex) && haskey(graph, bubble.exit_vertex)) || continue
        m1 = _branch_mean_coverage(graph, bubble.path1)
        m2 = _branch_mean_coverage(graph, bubble.path2)
        (m1 === nothing || m2 === nothing) && continue
        if m1 <= max_error_support && m2 >= min_real_support && m2 > m1
            remove_path_from_graph!(graph, bubble.path1, bubble.entry_vertex, bubble.exit_vertex)
            collapsed += 1
        elseif m2 <= max_error_support && m1 >= min_real_support && m1 > m2
            remove_path_from_graph!(graph, bubble.path2, bubble.entry_vertex, bubble.exit_vertex)
            collapsed += 1
        end
    end
    return (; collapsed = collapsed)
end

# ----------------------------------------------------------------------------
# Disconnected error-component pruning (td-byva)
# ----------------------------------------------------------------------------
#
# clip_error_tips!/collapse_error_bubbles! only touch structure WITHIN a
# component (dangling tips off a junction, error branches of a bubble). After
# RC-dedup + those passes the main genome collapses to a single full-length
# contig, but the graph still carries standalone DISCONNECTED components — small
# coverage-1/2 islands of sequencing-error k-mers that never joined the main
# path. Each is its own contig at extraction, so a 1kb genome can extract ~5
# contigs (1 real + ~4 error islands) and a 48kb genome ~98 (1 real + ~97).
#
# A standalone linear run is deliberately protected by the terminus guard in
# clip_error_tips! (it reaches its other dead end without a junction), so tip
# clipping cannot remove these islands — hence this dedicated pass. It removes a
# whole component ONLY when it is provably error under ALL of:
#   * SEPARATE component (never the main/largest component),
#   * SMALL (span <= max_component_length),
#   * uniformly LOW coverage (MAX per-vertex coverage <= max_component_support,
#     the ~1-2 error floor), AND
#   * REAL-SEQUENCE GUARD: shares no canonical k-mer with any well-supported
#     (coverage >= min_real_support) vertex. A standalone island that COULD be
#     genuine low-coverage minor sequence — because it re-uses real k-mer content
#     — is RETAINED. This biases to UNDER-pruning (td-h6w9): leave a debris
#     contig rather than risk removing a real low-coverage variant.

"""
    _label_dna_string(label) -> String

Best-effort DNA string for a vertex label. K-mer / BioSequence labels stringify
to their sequence; an already-string label is returned as-is.
"""
function _label_dna_string(label)
    label isa AbstractString && return label
    return String(label)
end

"""
    _accumulate_canonical_kmers!(kset::Set{String}, label, k::Int)

Add the canonical `k`-mer string keys of a vertex `label`'s sequence to `kset`.
Canonicalization is the lexicographic min of a k-mer and its reverse complement
(strand-agnostic key), so a run and its reverse-complement copy hash together —
this is what lets the real-sequence guard recognise reverse-strand overlap with a
supported contig. For a fixed-length k-mer graph each label contributes exactly
one key.
"""
function _accumulate_canonical_kmers!(kset::Set{String}, label, k::Int)
    s = _label_dna_string(label)
    n = length(s)
    n < k && return kset
    for i in 1:(n - k + 1)
        sub = s[i:(i + k - 1)]
        rc = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(sub)))
        push!(kset, min(sub, rc))
    end
    return kset
end

"""
    _connected_components_labels(labels, out_index, in_index, ::Type{T}) where {T}

Weakly-connected components of the graph, as label vectors, via BFS over the
UNDIRECTED closure of the prebuilt in/out adjacency indices. `O(V+E)`: every
vertex is enqueued once and every edge is followed a constant number of times.
"""
function _connected_components_labels(labels::AbstractVector{T},
        out_index::AbstractDict, in_index::AbstractDict, ::Type{T}) where {T}
    components = Vector{Vector{T}}()
    seen = Set{T}()
    queue = T[]
    for start in labels
        start in seen && continue
        push!(seen, start)
        empty!(queue)
        push!(queue, start)
        comp = T[]
        head = 1
        while head <= length(queue)
            v = queue[head]
            head += 1
            push!(comp, v)
            for nbr in get(out_index, v, T[])
                if !(nbr in seen)
                    push!(seen, nbr)
                    push!(queue, nbr)
                end
            end
            for nbr in get(in_index, v, T[])
                if !(nbr in seen)
                    push!(seen, nbr)
                    push!(queue, nbr)
                end
            end
        end
        push!(components, comp)
    end
    return components
end

"""
    prune_disconnected_error_components!(graph; k, max_component_length=2000,
        max_component_support=2, min_real_support=3) -> (; removed, components_pruned)

Remove standalone DISCONNECTED components that are provably sequencing-error
debris, in `O(V+E)`. A component is pruned ONLY when it is NOT the main (largest)
component AND its span `<= max_component_length` AND its MAXIMUM per-vertex
coverage `<= max_component_support` (uniformly at the ~1-2 error floor) AND it
shares no canonical k-mer with any well-supported (`coverage >= min_real_support`)
vertex (the real-sequence guard). Any component failing even one condition —
including a data-supported genome that happens to be disconnected — is RETAINED.

The design deliberately UNDER-prunes (td-h6w9): the coverage gate keys on the
component MAXIMUM (not mean), so a single supported k-mer anywhere in the
component vetoes removal, and the canonical-k-mer guard retains anything that
re-uses real sequence content. Span is estimated as `n_vertices + k - 1` (the
length of a linear k-mer chain), an over-estimate for branchy components, which
only makes the size cap MORE conservative.

The two SAFETY proofs are the coverage floor (`<= max_component_support`) and the
real-sequence guard (no shared canonical k-mer with a `>= min_real_support`
vertex). The size cap is a separate, OPTIONAL conservatism knob whose only job is
to spare a large but shallow island that might be a genuinely under-sequenced
real replicon; pass `max_component_length = nothing` to disable it and prune
large provably-error islands too (td-byva). Disabling the size cap never weakens
the two safety proofs — a well-supported or real-content-sharing component is
retained regardless of span.

# Empirical scope note (td-byva, 10-48 kb toy genomes at 20x illumina)
This pass removes only SEPARATE components; on the corrector's final k=21 graph at
these scales the corrected read graph is typically a SINGLE connected component
(0-2 disconnected islands), and the residual debris-contig count is dominated by
WITHIN-main-component branch fragmentation (find_contigs breaks a unitig at every
branch vertex), which is out of this pass's scope. When disconnected islands do
occur they are usually at real-sequence coverage (correctly retained by the
coverage floor). Extending this pass therefore cannot, by construction, reduce
within-component fragmentation; the `max_component_length = nothing` opt-in only
widens reach over genuine large disconnected error islands.

# Returns
- `(; removed::Int, components_pruned::Int)`: vertices removed and components
  pruned.
"""
function prune_disconnected_error_components!(graph::MetaGraphsNext.MetaGraph;
        k::Int,
        max_component_length::Union{Int, Nothing} = 2000,
        max_component_support::Real = 2,
        min_real_support::Real = 3)
    labels = collect(MetaGraphsNext.labels(graph))
    isempty(labels) && return (; removed = 0, components_pruned = 0)
    T = eltype(labels)
    out_index = _build_out_adjacency_index(graph, T)
    in_index = _build_in_adjacency_index(graph, T)

    components = _connected_components_labels(labels, out_index, in_index, T)
    # A single component means nothing is disconnected: this pass only removes
    # SEPARATE components, so structure within the one component is out of scope.
    length(components) <= 1 && return (; removed = 0, components_pruned = 0)

    # The largest component is the main assembly and is never a prune candidate.
    main_idx = argmax(map(length, components))

    # Canonical k-mers of every well-supported vertex, built once in O(V). A
    # candidate that touches any of these is retained (real-sequence guard).
    # Candidate vertices (max coverage <= max_component_support < min_real_support)
    # are never themselves in this set, so it cannot self-veto a genuine error.
    well_supported_kmers = Set{String}()
    for v in labels
        if _vertex_support(graph, v) >= min_real_support
            _accumulate_canonical_kmers!(well_supported_kmers, v, k)
        end
    end

    removed = 0
    components_pruned = 0
    for (ci, comp) in enumerate(components)
        ci == main_idx && continue
        # Size gate: an OPTIONAL conservatism knob (not one of the two safety
        # proofs). Its sole job is to spare a LARGE low-coverage island that might
        # be a genuinely under-sequenced real replicon (a real plasmid/segment can
        # be big yet shallow). `max_component_length === nothing` disables it, for
        # callers that know their input carries no low-coverage real replicon
        # (e.g. single high-coverage isolate) and want large error islands pruned.
        # The coverage floor and the real-sequence guard below remain binding
        # either way, so disabling the span gate never removes real sequence that
        # is well-supported OR shares canonical k-mer content with supported
        # sequence — it only widens the reach over provably-error debris.
        if max_component_length !== nothing
            span = length(comp) + k - 1
            span <= max_component_length || continue
        end
        all(v -> _vertex_support(graph, v) <= max_component_support, comp) || continue

        comp_kmers = Set{String}()
        for v in comp
            _accumulate_canonical_kmers!(comp_kmers, v, k)
        end
        any(km -> km in well_supported_kmers, comp_kmers) && continue

        for v in comp
            haskey(graph, v) && delete!(graph, v)
        end
        removed += length(comp)
        components_pruned += 1
    end
    return (; removed = removed, components_pruned = components_pruned)
end

"""
    clean_corrector_graph!(graph; k=21, clip_tips=true, collapse_bubbles=true, ...) -> Dict{String,Any}

Linear-time defragmentation of the corrector's final graph, run BEFORE contig
extraction (td-969e). Combines coverage-1 dead-end tip clipping, guarded
low-coverage bubble collapse, and disconnected error-component pruning (td-byva):
first tips/bubbles are alternated to a fixpoint (re-clipping any tips exposed by a
collapse), then whole standalone error-debris components are pruned. Every removal
targets an unambiguous error; a data-supported variant is never collapsed (the
td-h6w9 variation-preservation invariant is enforced by
[`clip_error_tips!`](@ref), [`collapse_error_bubbles!`](@ref), and
[`prune_disconnected_error_components!`](@ref)).

# Keywords
- `k::Int=21`: k-mer size; the tip-length ceiling is `tip_length_multiple * k`.
- `clip_tips::Bool=true`: run coverage-1 dead-end tip clipping.
- `collapse_bubbles::Bool=true`: run the guarded low-coverage bubble collapse.
- `prune_components::Bool=true`: run disconnected error-component pruning.
- `max_tip_support::Real=1`: max coverage for a removable tip vertex.
- `tip_length_multiple::Int=3`: tip-length ceiling in units of `k`.
- `max_error_support::Real=2`: max mean coverage for a collapsible bubble branch.
- `min_real_support::Real=3`: min mean coverage the surviving sibling must have.
- `max_bubble_length::Int=100`: max branch trace length for bubble detection.
- `max_rounds::Int=10`: max tip-clipping rounds per invocation.
- `max_component_length::Union{Int,Nothing}=2000`: max span of a prunable
  disconnected component, or `nothing` to disable the span gate and prune large
  provably-error islands too (the coverage floor + real-sequence guard stay
  binding). Default retains large shallow islands as possible under-sequenced
  real replicons.
- `max_component_support::Real=2`: max per-vertex coverage in a prunable component.

# Returns
- `Dict{String,Any}`: cleanup telemetry (vertices before/after, tips removed,
  bubbles collapsed, components pruned, component vertices removed).
"""
function clean_corrector_graph!(graph::MetaGraphsNext.MetaGraph;
        k::Int = 21,
        clip_tips::Bool = true,
        collapse_bubbles::Bool = true,
        prune_components::Bool = true,
        max_tip_support::Real = 2,
        tip_length_multiple::Int = 3,
        max_error_support::Real = 2,
        min_real_support::Real = 3,
        max_bubble_length::Int = 100,
        max_rounds::Int = 10,
        max_outer_rounds::Int = 5,
        max_component_length::Union{Int, Nothing} = 2000,
        max_component_support::Real = 2)
    n_before = length(collect(MetaGraphsNext.labels(graph)))
    tip_max_length = max(1, tip_length_multiple * k)
    tips_removed = 0
    bubbles_collapsed = 0
    components_pruned = 0
    component_vertices_removed = 0

    # Alternate tip clipping and guarded bubble collapse to a fixpoint: clipping a
    # tip can expose a bubble whose error branch is now reachable, and collapsing a
    # bubble can turn a former junction into a linear vertex that exposes a nested
    # tip. Each pass only removes unambiguous errors, so iterating to convergence
    # cannot cross into real-variant territory (the td-h6w9 invariant holds
    # per-pass). Bounded by max_outer_rounds.
    for _ in 1:max_outer_rounds
        round_tips = 0
        round_bubbles = 0
        if clip_tips
            r = clip_error_tips!(graph; max_tip_length = tip_max_length,
                max_tip_support = max_tip_support, max_rounds = max_rounds)
            round_tips += r.removed
        end
        if collapse_bubbles
            r = collapse_error_bubbles!(graph; max_error_support = max_error_support,
                min_real_support = min_real_support, max_bubble_length = max_bubble_length)
            round_bubbles += r.collapsed
        end
        tips_removed += round_tips
        bubbles_collapsed += round_bubbles
        (round_tips == 0 && round_bubbles == 0) && break
    end

    # Prune standalone error-debris components AFTER the tip/bubble fixpoint.
    # Tip clipping and bubble collapse only reshape structure WITHIN a component
    # and never split one component into two, so the disconnected-component set is
    # stable and a single pruning pass suffices. The main (largest) component is
    # always retained; only small, uniformly-low-coverage islands that share no
    # real k-mer content are removed (td-byva / td-h6w9).
    if prune_components
        r = prune_disconnected_error_components!(graph; k = k,
            max_component_length = max_component_length,
            max_component_support = max_component_support,
            min_real_support = min_real_support)
        components_pruned += r.components_pruned
        component_vertices_removed += r.removed
    end

    return Dict{String, Any}(
        "graph_cleanup_vertices_before" => n_before,
        "graph_cleanup_vertices_after" => length(collect(MetaGraphsNext.labels(graph))),
        "graph_cleanup_tips_removed" => tips_removed,
        "graph_cleanup_bubbles_collapsed" => bubbles_collapsed,
        "graph_cleanup_components_pruned" => components_pruned,
        "graph_cleanup_component_vertices_removed" => component_vertices_removed,
    )
end
