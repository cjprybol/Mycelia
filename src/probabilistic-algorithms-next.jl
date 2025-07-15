"""
Probabilistic algorithms for strand-aware k-mer graph traversal and assembly.

This module implements the core algorithms for Phase 2 of the assembly roadmap:
- Probabilistic walks with strand-consistent transitions
- Shortest probability paths (distance ∝ -log(probability))
- Maximum weight walks for high-confidence path finding
- Quality-aware scoring and path validation
"""

import MetaGraphsNext
import Graphs
import DataStructures
import Statistics
import Random
using DocStringExtensions

"""
Represents a step in a probabilistic walk through the graph.
"""
struct WalkStep
    vertex_label::String
    strand::StrandOrientation
    probability::Float64
    cumulative_probability::Float64
end

"""
Represents a complete path through the k-mer graph.
"""
struct GraphPath
    steps::Vector{WalkStep}
    total_probability::Float64
    sequence::String
    
    function GraphPath(steps::Vector{WalkStep})
        total_prob = isempty(steps) ? 0.0 : last(steps).cumulative_probability
        sequence = _reconstruct_sequence_from_path(steps)
        new(steps, total_prob, sequence)
    end
end

"""
$(TYPEDSIGNATURES)

Perform a probabilistic walk through the strand-aware k-mer graph.

This algorithm follows edges based on their probability weights, respecting strand
orientation constraints. The walk continues until max_steps is reached or no valid
transitions are available.

# Arguments
- `graph`: MetaGraphsNext k-mer graph with strand-aware edges
- `start_vertex`: Starting k-mer (vertex label)
- `max_steps`: Maximum number of steps to take
- `seed`: Random seed for reproducibility (optional)

# Returns
- `GraphPath`: Complete path with probability information

# Algorithm
1. Start at given vertex with forward strand orientation
2. At each step, calculate transition probabilities based on edge weights
3. Sample next vertex according to probabilities
4. Update cumulative probability and continue
5. Respect strand orientation constraints from edge metadata

# Example
```julia
graph = build_kmer_graph_next(DNAKmer{15}, observations)
path = probabilistic_walk_next(graph, "ATCGATCGATCGATC", 100)
println("Assembled sequence: \$(path.sequence)")
println("Path probability: \$(path.total_probability)")
```
"""
function probabilistic_walk_next(graph::MetaGraphsNext.MetaGraph, 
                                start_vertex::String, 
                                max_steps::Int;
                                seed::Union{Nothing, Int}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end
    
    if !(start_vertex in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Start vertex $start_vertex not found in graph"))
    end
    
    steps = Vector{WalkStep}()
    current_vertex = start_vertex
    current_strand = Forward  # Start with forward orientation
    cumulative_prob = 1.0
    
    # Add starting step
    push!(steps, WalkStep(current_vertex, current_strand, 1.0, cumulative_prob))
    
    for step in 1:max_steps
        # Get all valid outgoing edges from current vertex
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        
        if isempty(valid_transitions)
            # No valid transitions available
            break
        end
        
        # Calculate transition probabilities
        transition_probs = _calculate_transition_probabilities(valid_transitions)
        
        # Sample next transition
        next_transition = _sample_transition(valid_transitions, transition_probs)
        
        # Update path
        step_prob = next_transition[:probability]
        cumulative_prob *= step_prob
        
        current_vertex = next_transition[:target_vertex]
        current_strand = next_transition[:target_strand]
        
        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end
    
    return GraphPath(steps)
end

"""
$(TYPEDSIGNATURES)

Find the shortest path in probability space between two vertices.

Uses Dijkstra's algorithm where edge distances are -log(probability), so the
shortest path corresponds to the highest probability path.

# Arguments
- `graph`: MetaGraphsNext k-mer graph
- `source`: Source vertex label
- `target`: Target vertex label

# Returns
- `Union{GraphPath, Nothing}`: Shortest probability path, or nothing if no path exists

# Algorithm
1. Convert edge weights to -log(probability) distances
2. Run Dijkstra's algorithm with strand-aware edge traversal
3. Reconstruct path maintaining strand information
4. Convert back to probability space for final result
"""
function shortest_probability_path_next(graph::MetaGraphsNext.MetaGraph, 
                                       source::String, 
                                       target::String)
    if !(source in MetaGraphsNext.labels(graph)) || !(target in MetaGraphsNext.labels(graph))
        return nothing
    end
    
    # State: (vertex_label, strand_orientation)
    distances = Dict{Tuple{String, StrandOrientation}, Float64}()
    predecessors = Dict{Tuple{String, StrandOrientation}, Union{Nothing, Tuple{String, StrandOrientation}}}()
    visited = Set{Tuple{String, StrandOrientation}}()
    
    # Priority queue: (distance, vertex_label, strand)
    pq = DataStructures.PriorityQueue{Tuple{String, StrandOrientation}, Float64}()
    
    # Initialize
    start_state = (source, Forward)
    distances[start_state] = 0.0
    predecessors[start_state] = nothing
    pq[start_state] = 0.0
    
    while !isempty(pq)
        current_state = DataStructures.dequeue!(pq)
        current_vertex, current_strand = current_state
        
        if current_state in visited
            continue
        end
        push!(visited, current_state)
        
        # Check if we reached target
        if current_vertex == target
            # Reconstruct path
            return _reconstruct_shortest_path(predecessors, distances, start_state, current_state, graph)
        end
        
        # Explore neighbors
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        
        for transition in valid_transitions
            neighbor_vertex = transition[:target_vertex]
            neighbor_strand = transition[:target_strand]
            neighbor_state = (neighbor_vertex, neighbor_strand)
            
            if neighbor_state in visited
                continue
            end
            
            # Distance is -log(probability)
            edge_prob = transition[:probability]
            edge_distance = edge_prob > 0 ? -log(edge_prob) : Inf
            new_distance = distances[current_state] + edge_distance
            
            if !haskey(distances, neighbor_state) || new_distance < distances[neighbor_state]
                distances[neighbor_state] = new_distance
                predecessors[neighbor_state] = current_state
                pq[neighbor_state] = new_distance
            end
        end
    end
    
    return nothing  # No path found
end

"""
$(TYPEDSIGNATURES)

Perform a maximum weight walk prioritizing highest confidence edges.

This greedy algorithm always chooses the edge with the highest weight (coverage)
at each step, useful for finding high-confidence assembly paths.

# Arguments
- `graph`: MetaGraphsNext k-mer graph
- `start_vertex`: Starting vertex label
- `max_steps`: Maximum steps to take
- `weight_function`: Function to extract weight from edge data (default: uses edge.weight)

# Returns
- `GraphPath`: Path following maximum weight edges
"""
function maximum_weight_walk_next(graph::MetaGraphsNext.MetaGraph,
                                 start_vertex::String,
                                 max_steps::Int;
                                 weight_function::Function = edge_data -> edge_data.weight)
    if !(start_vertex in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Start vertex $start_vertex not found in graph"))
    end
    
    steps = Vector{WalkStep}()
    current_vertex = start_vertex
    current_strand = Forward
    cumulative_prob = 1.0
    
    # Add starting step
    push!(steps, WalkStep(current_vertex, current_strand, 1.0, cumulative_prob))
    
    for step in 1:max_steps
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        
        if isempty(valid_transitions)
            break
        end
        
        # Find transition with maximum weight
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
        
        # Follow best transition
        step_prob = best_transition[:probability]
        cumulative_prob *= step_prob
        
        current_vertex = best_transition[:target_vertex]
        current_strand = best_transition[:target_strand]
        
        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end
    
    return GraphPath(steps)
end

"""
Helper function to get valid transitions from a vertex with given strand orientation.
"""
function _get_valid_transitions(graph, vertex_label, strand)
    transitions = []
    
    # Get all outgoing edges from this vertex
    for edge_labels in MetaGraphsNext.edge_labels(graph)
        if length(edge_labels) == 2 && edge_labels[1] == vertex_label
            target_vertex = edge_labels[2]
            edge_data = graph[edge_labels...]
            
            # Check if this edge is valid for our current strand
            if edge_data.src_strand == strand
                probability = edge_data.weight > 0 ? edge_data.weight : 1e-10
                
                push!(transitions, Dict(
                    :target_vertex => target_vertex,
                    :target_strand => edge_data.dst_strand,
                    :probability => probability,
                    :edge_data => edge_data
                ))
            end
        end
    end
    
    return transitions
end

"""
Helper function to calculate normalized transition probabilities.
"""
function _calculate_transition_probabilities(transitions)
    if isempty(transitions)
        return Float64[]
    end
    
    weights = [t[:probability] for t in transitions]
    total_weight = sum(weights)
    
    if total_weight == 0
        # Equal probability for all transitions
        return fill(1.0 / length(transitions), length(transitions))
    else
        return weights ./ total_weight
    end
end

"""
Helper function to sample a transition based on probabilities.
"""
function _sample_transition(transitions, probabilities)
    if isempty(transitions)
        return nothing
    end
    
    if length(transitions) == 1
        return first(transitions)
    end
    
    # Sample according to probabilities
    r = rand()
    cumulative = 0.0
    
    for (i, prob) in enumerate(probabilities)
        cumulative += prob
        if r <= cumulative
            return transitions[i]
        end
    end
    
    # Fallback to last transition
    return last(transitions)
end

"""
Helper function to reconstruct sequence from a graph path.
"""
function _reconstruct_sequence_from_path(steps)
    if isempty(steps)
        return ""
    end
    
    # Start with first k-mer
    first_step = first(steps)
    first_vertex_data = first_step.vertex_label  # This should be the k-mer sequence
    
    if first_step.strand == Forward
        sequence = first_vertex_data
    else
        # Reverse complement for reverse strand
        sequence = string(BioSequences.reverse_complement(BioSequences.DNASequence(first_vertex_data)))
    end
    
    # Add subsequent nucleotides
    for i in 2:length(steps)
        step = steps[i]
        vertex_sequence = step.vertex_label
        
        # Get the sequence considering strand
        if step.strand == Forward
            kmer_seq = vertex_sequence
        else
            kmer_seq = string(BioSequences.reverse_complement(BioSequences.DNASequence(vertex_sequence)))
        end
        
        # Add the last nucleotide (assuming k-mer overlap)
        if length(kmer_seq) > 0
            sequence *= last(kmer_seq)
        end
    end
    
    return sequence
end

"""
Helper function to reconstruct shortest path from Dijkstra's algorithm.
"""
function _reconstruct_shortest_path(predecessors, distances, start_state, end_state, graph)
    path_states = []
    current_state = end_state
    
    while current_state !== nothing
        pushfirst!(path_states, current_state)
        current_state = predecessors[current_state]
    end
    
    # Convert to WalkStep format
    steps = Vector{WalkStep}()
    total_distance = distances[end_state]
    total_probability = exp(-total_distance)
    
    for (i, (vertex, strand)) in enumerate(path_states)
        step_prob = if i == 1
            1.0
        else
            # Calculate step probability from distance difference
            prev_state = path_states[i-1]
            step_distance = distances[(vertex, strand)] - distances[prev_state]
            exp(-step_distance)
        end
        
        cumulative_prob = exp(-distances[(vertex, strand)])
        push!(steps, WalkStep(vertex, strand, step_prob, cumulative_prob))
    end
    
    return GraphPath(steps)
end

# Export the main functions
export WalkStep, GraphPath
export probabilistic_walk_next, shortest_probability_path_next, maximum_weight_walk_next
