"""
Graph Cleanup and Statistical Tip Clipping

This module implements advanced graph cleanup algorithms based on Mycelia-Dev
research, particularly the statistical tip clipping algorithm from the
2021-02-02 notebook that uses coverage-based statistics for error removal.

Key innovation: Remove low-evidence tips using statistical thresholds relative
to the connected component's coverage distribution, rather than simple cutoffs.
"""

# =============================================================================
# Core Statistical Tip Clipping Algorithm
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Perform statistical tip clipping on assembly graphs using coverage-based thresholds.

Based on the novel approach from Mycelia-Dev 2021-02-02 notebook. Removes graph
tips (degree-1 nodes) that have significantly lower coverage than the rest of 
their connected component, using statistical analysis rather than hard cutoffs.

# Arguments
- `graph`: Assembly graph with nodes and coverage information
- `min_coverage_threshold::Int`: Hard minimum coverage threshold (default: 1)
- `std_dev_multiplier::Float64`: Standard deviation multiplier for statistical threshold (default: 3.0)
- `preserve_high_quality::Bool`: Preserve tips with high quality scores (default: true)

# Returns
- Modified graph with low-evidence tips removed
- `Dict`: Statistics about removed tips and cleanup effectiveness

# Algorithm
1. Identify all connected components in the graph
2. For each connected component:
   - Calculate coverage statistics (median, standard deviation)
   - Identify degree-1 nodes (tips)
   - Remove tips that are either:
     - Single-coverage (1x coverage) OR
     - More than `std_dev_multiplier` standard deviations below median coverage
3. Preserve high-quality tips if specified

# Statistical Foundation
Uses the principle that erroneous tips will have coverage distributions significantly
different from the true backbone of the assembly graph.

# Example
```julia
# Standard tip clipping with 3σ threshold
cleaned_graph, stats = statistical_tip_clipping(graph)

# More aggressive cleanup with 2σ threshold  
cleaned_graph, stats = statistical_tip_clipping(graph, 1, 2.0)

println("Removed \$(stats[:tips_removed]) tips from \$(stats[:connected_components]) components")
```

# References
Based on 2021-02-02-tip-clipping.ipynb from Mycelia-Dev research demonstrating
coverage-based statistical tip removal for assembly graph cleanup.
"""
function statistical_tip_clipping(graph; 
                                 min_coverage_threshold::Int = 1,
                                 std_dev_multiplier::Float64 = 3.0,
                                 preserve_high_quality::Bool = true)
    
    # Initialize statistics tracking
    stats = Dict{Symbol, Any}(
        :tips_removed => 0,
        :connected_components => 0,
        :single_coverage_removed => 0,
        :statistical_threshold_removed => 0,
        :high_quality_preserved => 0
    )
    
    # Find all connected components
    connected_components = find_connected_components(graph)
    stats[:connected_components] = length(connected_components)
    
    # Process each connected component independently
    for component in connected_components
        component_stats = process_connected_component_tips(
            graph, component, min_coverage_threshold, std_dev_multiplier, preserve_high_quality
        )
        
        # Aggregate statistics
        stats[:tips_removed] += component_stats[:tips_removed]
        stats[:single_coverage_removed] += component_stats[:single_coverage_removed]
        stats[:statistical_threshold_removed] += component_stats[:statistical_threshold_removed]
        stats[:high_quality_preserved] += component_stats[:high_quality_preserved]
    end
    
    return graph, stats
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Process tip clipping for a single connected component.

# Arguments
- `graph`: Assembly graph
- `component::Vector{Int}`: Node indices in the connected component
- `min_coverage_threshold::Int`: Hard minimum coverage
- `std_dev_multiplier::Float64`: Statistical threshold multiplier
- `preserve_high_quality::Bool`: Whether to preserve high-quality tips

# Returns
- `Dict`: Statistics for this component's tip clipping
"""
function process_connected_component_tips(graph, component::Vector{Int},
                                        min_coverage_threshold::Int,
                                        std_dev_multiplier::Float64,
                                        preserve_high_quality::Bool)
    
    stats = Dict{Symbol, Int}(
        :tips_removed => 0,
        :single_coverage_removed => 0,
        :statistical_threshold_removed => 0,
        :high_quality_preserved => 0
    )
    
    # Extract coverage values for this connected component
    coverage_values = get_coverage_values(graph, component)
    
    if length(coverage_values) < 2
        # Skip components that are too small for statistical analysis
        return stats
    end
    
    # Calculate coverage statistics
    median_coverage = Statistics.median(coverage_values)
    coverage_std = Statistics.std(coverage_values)
    statistical_threshold = median_coverage - std_dev_multiplier * coverage_std
    
    # Identify tips (degree-1 nodes) in this component
    tips = identify_tips_in_component(graph, component)
    
    # Process each tip
    for tip_node in tips
        tip_coverage = get_node_coverage(graph, tip_node)
        tip_degree = get_node_degree(graph, tip_node)
        
        # Confirm this is actually a tip (degree = 1)
        if tip_degree != 1
            continue
        end
        
        # Determine if tip should be removed
        should_remove = false
        removal_reason = :none
        
        # Hard minimum coverage threshold
        if tip_coverage <= min_coverage_threshold
            should_remove = true
            removal_reason = :single_coverage
        
        # Statistical threshold based on component distribution
        elseif tip_coverage < statistical_threshold
            should_remove = true
            removal_reason = :statistical_threshold
        end
        
        # Check for high-quality preservation
        if should_remove && preserve_high_quality && is_high_quality_tip(graph, tip_node)
            should_remove = false
            stats[:high_quality_preserved] += 1
        end
        
        # Remove the tip if criteria are met
        if should_remove
            remove_graph_node!(graph, tip_node)
            stats[:tips_removed] += 1
            
            if removal_reason == :single_coverage
                stats[:single_coverage_removed] += 1
            elseif removal_reason == :statistical_threshold
                stats[:statistical_threshold_removed] += 1
            end
        end
    end
    
    return stats
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Check if a tip should be preserved due to high quality indicators.

# Arguments
- `graph`: Assembly graph
- `tip_node::Int`: Node index of the tip to evaluate

# Returns
- `Bool`: true if tip should be preserved due to high quality

# Quality indicators considered:
- High-quality k-mer sequences (low ambiguity)
- Strong connecting edge weight
- Consistent with long-read evidence (if available)
"""
function is_high_quality_tip(graph, tip_node::Int)::Bool
    # Get node properties
    node_data = get_node_data(graph, tip_node)
    
    # Check for high-quality sequence (low N content, reasonable complexity)
    if haskey(node_data, :sequence)
        sequence_quality = assess_sequence_quality(node_data.sequence)
        if sequence_quality > 0.9  # High quality threshold
            return true
        end
    end
    
    # Check for strong connecting edge
    if haskey(node_data, :edge_weights)
        max_edge_weight = maximum(node_data.edge_weights)
        if max_edge_weight > 10.0  # Strong connection threshold
            return true
        end
    end
    
    # Check for long-read support (if available)
    if haskey(node_data, :long_read_support) && node_data.long_read_support > 5
        return true
    end
    
    return false
end

# =============================================================================
# Graph Structure Analysis Functions
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Find all connected components in the assembly graph.

# Arguments
- `graph`: Assembly graph structure

# Returns
- `Vector{Vector{Int}}`: List of connected components, each containing node indices

# Example
```julia
components = find_connected_components(graph)
println("Found \$(length(components)) connected components")
```
"""
function find_connected_components(graph)::Vector{Vector{Int}}
    # This would use the graph structure to identify connected components
    # Implementation depends on the specific graph representation
    
    # Placeholder implementation - would use actual graph traversal
    components = Vector{Vector{Int}}()
    
    # Example: if graph has an adjacency list representation
    visited = Set{Int}()
    
    for node in get_all_nodes(graph)
        if node ∉ visited
            component = depth_first_search_component(graph, node, visited)
            push!(components, component)
        end
    end
    
    return components
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Identify tip nodes (degree-1) within a connected component.

# Arguments
- `graph`: Assembly graph
- `component::Vector{Int}`: Node indices in the component

# Returns
- `Vector{Int}`: Node indices of tips in this component
"""
function identify_tips_in_component(graph, component::Vector{Int})::Vector{Int}
    tips = Int[]
    
    for node in component
        if get_node_degree(graph, node) == 1
            push!(tips, node)
        end
    end
    
    return tips
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Perform depth-first search to find all nodes in a connected component.

# Arguments
- `graph`: Assembly graph
- `start_node::Int`: Starting node for DFS
- `visited::Set{Int}`: Set of already visited nodes (modified in place)

# Returns  
- `Vector{Int}`: All nodes in the connected component containing start_node
"""
function depth_first_search_component(graph, start_node::Int, visited::Set{Int})::Vector{Int}
    component = Int[]
    stack = Int[start_node]
    
    while !isempty(stack)
        current_node = pop!(stack)
        
        if current_node ∉ visited
            push!(visited, current_node)
            push!(component, current_node)
            
            # Add all neighbors to stack
            for neighbor in get_node_neighbors(graph, current_node)
                if neighbor ∉ visited
                    push!(stack, neighbor)
                end
            end
        end
    end
    
    return component
end

# =============================================================================
# Bubble Removal and Complex Structure Cleanup
# =============================================================================

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Remove simple bubbles from the assembly graph.

Bubbles are alternative paths between the same start and end nodes, often
caused by sequencing errors or minor variants. This function identifies and
removes the lower-coverage path in simple bubble structures.

# Arguments
- `graph`: Assembly graph
- `max_bubble_length::Int`: Maximum length of bubbles to consider (default: 10)
- `coverage_ratio_threshold::Float64`: Minimum ratio between paths to remove bubble (default: 0.3)

# Returns
- Modified graph with bubbles removed
- `Dict`: Statistics about bubble removal

# Algorithm
1. Identify potential bubble start nodes (out-degree > 1)
2. Trace alternative paths of limited length
3. Check if paths reconverge at a common node
4. Compare coverage between alternative paths
5. Remove lower-coverage path if ratio exceeds threshold

# Example
```julia
cleaned_graph, bubble_stats = remove_simple_bubbles(graph, 15, 0.25)
println("Removed \$(bubble_stats[:bubbles_removed]) bubbles")
```
"""
function remove_simple_bubbles(graph; 
                              max_bubble_length::Int = 10,
                              coverage_ratio_threshold::Float64 = 0.3)
    
    stats = Dict{Symbol, Int}(
        :bubbles_removed => 0,
        :nodes_removed => 0,
        :bubbles_examined => 0
    )
    
    # Find all potential bubble start points (nodes with out-degree > 1)
    bubble_starts = find_bubble_start_candidates(graph)
    
    for start_node in bubble_starts
        bubble_result = identify_and_remove_bubble(
            graph, start_node, max_bubble_length, coverage_ratio_threshold
        )
        
        stats[:bubbles_examined] += 1
        if bubble_result[:removed]
            stats[:bubbles_removed] += 1
            stats[:nodes_removed] += bubble_result[:nodes_removed]
        end
    end
    
    return graph, stats
end

"""
    $(DocStringExtensions.TYPEDSIGNATURES)

Identify and potentially remove a bubble starting from a given node.

# Arguments
- `graph`: Assembly graph
- `start_node::Int`: Potential bubble start node
- `max_length::Int`: Maximum bubble length to consider
- `coverage_threshold::Float64`: Coverage ratio threshold for removal

# Returns
- `Dict`: Information about bubble detection and removal
"""
function identify_and_remove_bubble(graph, start_node::Int, max_length::Int, 
                                   coverage_threshold::Float64)::Dict{Symbol, Any}
    
    result = Dict{Symbol, Any}(
        :removed => false,
        :nodes_removed => 0,
        :bubble_detected => false,
        :coverage_ratio => 0.0
    )
    
    # Get alternative paths from start node
    out_neighbors = get_node_neighbors(graph, start_node, :out)
    
    if length(out_neighbors) != 2
        # Not a simple bubble candidate
        return result
    end
    
    # Trace both paths up to max_length
    path1 = trace_path(graph, out_neighbors[1], max_length)
    path2 = trace_path(graph, out_neighbors[2], max_length)
    
    # Check if paths reconverge
    if !paths_reconverge(path1, path2)
        return result
    end
    
    result[:bubble_detected] = true
    
    # Calculate coverage for both paths
    coverage1 = calculate_path_coverage(graph, path1)
    coverage2 = calculate_path_coverage(graph, path2)
    
    # Determine which path to remove
    if coverage1 > coverage2
        ratio = coverage2 / coverage1
        path_to_remove = path2
    else
        ratio = coverage1 / coverage2
        path_to_remove = path1
    end
    
    result[:coverage_ratio] = ratio
    
    # Remove lower-coverage path if ratio is below threshold
    if ratio < coverage_threshold
        for node in path_to_remove[2:end-1]  # Exclude start and end nodes
            remove_graph_node!(graph, node)
            result[:nodes_removed] += 1
        end
        result[:removed] = true
    end
    
    return result
end

# =============================================================================
# Helper Functions (Graph Interface Abstraction)
# =============================================================================

"""
Helper functions that provide an abstraction layer for different graph implementations.
These would need to be implemented based on the actual graph data structure used.
"""

function get_coverage_values(graph, component::Vector{Int})::Vector{Float64}
    # Implementation depends on graph structure
    return [get_node_coverage(graph, node) for node in component]
end

function get_node_coverage(graph, node::Int)::Float64
    # Implementation depends on graph structure
    node_data = get_node_data(graph, node)
    return haskey(node_data, :coverage) ? node_data.coverage : 1.0
end

function get_node_degree(graph, node::Int)::Int
    # Implementation depends on graph structure
    neighbors = get_node_neighbors(graph, node)
    return length(neighbors)
end

function get_node_data(graph, node::Int)
    # Implementation depends on graph structure
    # Returns node metadata/properties
    return Dict{Symbol, Any}()
end

function get_node_neighbors(graph, node::Int, direction::Symbol = :all)::Vector{Int}
    # Implementation depends on graph structure
    # direction can be :all, :in, :out for directed graphs
    return Int[]
end

function get_all_nodes(graph)::Vector{Int}
    # Implementation depends on graph structure
    return Int[]
end

function remove_graph_node!(graph, node::Int)
    # Implementation depends on graph structure
    # Removes node and associated edges
    return nothing
end

function assess_sequence_quality(sequence::String)::Float64
    # Simple sequence quality assessment
    if isempty(sequence)
        return 0.0
    end
    
    # Count N's and ambiguous bases
    n_count = count(c -> c == 'N' || c == 'n', sequence)
    n_fraction = n_count / length(sequence)
    
    # Basic complexity assessment (avoid homopolymers)
    unique_chars = length(Set(uppercase(sequence)))
    complexity_score = unique_chars / 4.0  # Normalize by DNA alphabet size
    
    # Combine metrics
    quality = (1.0 - n_fraction) * complexity_score
    return min(1.0, quality)
end

function find_bubble_start_candidates(graph)::Vector{Int}
    # Find nodes with out-degree > 1 (potential bubble starts)
    candidates = Int[]
    
    for node in get_all_nodes(graph)
        out_neighbors = get_node_neighbors(graph, node, :out)
        if length(out_neighbors) > 1
            push!(candidates, node)
        end
    end
    
    return candidates
end

function trace_path(graph, start_node::Int, max_length::Int)::Vector{Int}
    # Trace a path from start_node up to max_length
    path = [start_node]
    current_node = start_node
    
    for _ in 1:max_length
        neighbors = get_node_neighbors(graph, current_node, :out)
        
        if length(neighbors) != 1
            # Path branches or ends
            break
        end
        
        next_node = neighbors[1]
        push!(path, next_node)
        current_node = next_node
    end
    
    return path
end

function paths_reconverge(path1::Vector{Int}, path2::Vector{Int})::Bool
    # Check if two paths end at the same node
    return !isempty(path1) && !isempty(path2) && last(path1) == last(path2)
end

function calculate_path_coverage(graph, path::Vector{Int})::Float64
    # Calculate average coverage along a path
    if isempty(path)
        return 0.0
    end
    
    total_coverage = sum(get_node_coverage(graph, node) for node in path)
    return total_coverage / length(path)
end