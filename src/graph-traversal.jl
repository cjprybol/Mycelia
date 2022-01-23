"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_edge_probabilities(graph, strand)
    kmers = graph.gprops[:kmers]
    outgoing_edge_probabilities = SparseArrays.spzeros(length(kmers), length(kmers))
    
    for (kmer_index, kmer) in enumerate(kmers)
        if !strand
            kmer = BioSequences.reverse_complement(kmer)
        end
        
        downstream_neighbor_indices = Int[]
        for neighbor in BioSequences.neighbors(kmer)
            index = get_kmer_index(kmers, BioSequences.canonical(neighbor))
            # kmer must be in our dataset and there must be a connecting edge
            if !isnothing(index) && Graphs.has_edge(graph, ordered_edge(kmer_index, index))
                push!(downstream_neighbor_indices, index)
            end
        end
        sort!(unique!(downstream_neighbor_indices))
        
        downstream_edge_weights = Int[
            length(get(graph.edge_evidence, ordered_edge(kmer_index, neighbor_index), EdgeEvidence[])) for neighbor_index in downstream_neighbor_indices
        ]
        
        non_zero_indices = downstream_edge_weights .> 0
        downstream_neighbor_indices = downstream_neighbor_indices[non_zero_indices]
        downstream_edge_weights = downstream_edge_weights[non_zero_indices]
        
        downstream_edge_likelihoods = downstream_edge_weights ./ sum(downstream_edge_weights)
        
        for (neighbor_index, likelihood) in zip(downstream_neighbor_indices, downstream_edge_likelihoods)
            outgoing_edge_probabilities[kmer_index, neighbor_index] = likelihood
        end
    end
    return outgoing_edge_probabilities
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_edge_probabilities(graph)
    outgoing_edge_probabilities = determine_edge_probabilities(graph, true)
    incoming_edge_probabilities = determine_edge_probabilities(graph, false)
    return outgoing_edge_probabilities, incoming_edge_probabilities
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

```jldoctest
julia> 1 + 1
2
```
"""
function oriented_unbranching_walk(kmer_graph, vertex, orientation)
    walk = []
    viable_neighbors = find_unbranched_neighbors(kmer_graph, vertex, orientation)
    while length(viable_neighbors) == 1
#         @show "found a viable neighbor!!"
        viable_neighbor = first(viable_neighbors)
        edge = Graphs.Edge(vertex, viable_neighbor)
        push!(walk, edge)
        vertex = edge.dst
        viable_neighbors = Set{Int}()
        destination_orientations = [o.destination_orientation for o in kmer_graph.eprops[edge][:orientations]]
        for destination_orientation in destination_orientations
            union!(viable_neighbors, find_unbranched_neighbors(kmer_graph, vertex, destination_orientation))
        end
    end
    return walk
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

```jldoctest
julia> 1 + 1
2
```
"""
function find_unbranched_neighbors(kmer_graph, vertex, orientation)
    downstream_vertices = find_downstream_vertices(kmer_graph, vertex, orientation)
#     backtrack_vertices
#     @show downstream_vertices
    if length(downstream_vertices) == 1
        downstream_vertex = first(downstream_vertices)
#         @show downstream_vertex
        edge = Graphs.Edge(vertex, downstream_vertex)
#         @show edge
        destination_orientations = [o.destination_orientation for o in kmer_graph.eprops[edge][:orientations]]
#         @show destination_orientations
        for destination_orientation in destination_orientations
            backtrack_vertices = find_downstream_vertices(kmer_graph, downstream_vertex, !destination_orientation)
#             @show backtrack_vertices
            # if the only backtrack is the vertex we're on, then we can simplify
            if backtrack_vertices == Set([vertex])
                return downstream_vertices
            end
        end
    end
    return Int[]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

```jldoctest
julia> 1 + 1
2
```
"""
function find_downstream_vertices(kmer_graph, vertex, orientation)
    viable_neighbors = Set{Int}()
    for neighbor in Graphs.neighbors(kmer_graph, vertex)
        not_same_vertex = vertex != neighbor
        candidate_edge = Graphs.Edge(vertex, neighbor)
        # palindromes can have multiple viable orientations
        # check each viable orientation individually
        edge_src_orientations = [e.source_orientation for e in kmer_graph.eprops[candidate_edge][:orientations]]
        for edge_src_orientation in edge_src_orientations
#             edge_src_orientation = kmer_graph.eprops[candidate_edge][:orientations].source_orientation
            viable_orientation = edge_src_orientation == orientation
            if not_same_vertex && viable_orientation
                push!(viable_neighbors, neighbor)
            end
        end
    end
    return viable_neighbors
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function reverse_oriented_path(oriented_path)
#     reversed_path = copy(oriented_path)
#     for (index, state) in enumerate(oriented_path)
#         reversed_path[index] = Mycelia.OrientedKmer(index = state.index, orientation = !state.orientation)
#     end
#     return reverse!(reversed_path)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function maximum_likelihood_walk(graph, connected_component)
#     max_count = maximum(graph.counts[connected_component])
#     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
#     initial_node_index = rand(max_count_indices)
#     initial_node = connected_component[initial_node_index]
#     outgoing_edge_probabilities, incoming_edge_probabilities = Mycelia.determine_edge_probabilities(graph)
#     forward_walk = maximum_likelihood_walk(graph, [Mycelia.OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
#     reverse_walk = maximum_likelihood_walk(graph, [Mycelia.OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
#     reversed_reverse_walk = reverse!(
#         [
#             Mycelia.OrientedKmer(index = oriented_kmer.index, orientation = oriented_kmer.orientation)
#             for oriented_kmer in reverse_walk[2:end]
#         ]
#         )
#     full_path = [reversed_reverse_walk..., forward_walk...]
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function maximum_likelihood_walk(graph, path::Vector{Mycelia.OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
#     done = false
#     while !done
#         maximum_path_likelihood = 0.0
#         maximum_likelihood_path = Vector{Mycelia.OrientedKmer}()
#         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
#             this_path = [last(path).index, neighbor]
#             this_oriented_path, this_path_likelihood = 
#                 Mycelia.assess_path(this_path,
#                     graph.kmers,
#                     graph.counts,
#                     last(path).orientation,
#                     outgoing_edge_probabilities,
#                     incoming_edge_probabilities)
#             if this_path_likelihood > maximum_path_likelihood
#                 maximum_path_likelihood = this_path_likelihood
#                 maximum_likelihood_path = this_oriented_path
#             end
#         end
#         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
#             done = true
#         else
#             append!(path, maximum_likelihood_path[2:end])
#         end
#     end
#     return path
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function take_a_walk(graph, connected_component)
#     max_count = maximum(graph.counts[connected_component])
#     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
#     initial_node_index = rand(max_count_indices)
#     initial_node = connected_component[initial_node_index]
#     outgoing_edge_probabilities, incoming_edge_probabilities = Mycelia.determine_edge_probabilities(graph)
    
#     # walk forwards from the initial starting node
#     forward_walk = take_a_walk(graph, [Mycelia.OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
#     # walk backwards from the initial starting node
#     reverse_walk = take_a_walk(graph, [Mycelia.OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
#     # we need to reverse everything to re-orient against the forward walk
#     reverse_walk = reverse_oriented_path(reverse_walk)
    
#     # also need to drop the last node, which is equivalent to the first node of the 
#     @assert last(reverse_walk) == first(forward_walk)
#     full_path = [reverse_walk[1:end-1]..., forward_walk...]
# #     @show full_path
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function take_a_walk(graph, path::Vector{Mycelia.OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
#     done = false
#     while !done
#         maximum_path_likelihood = 0.0
#         maximum_likelihood_path = Vector{Mycelia.OrientedKmer}()
#         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
#             this_path = [last(path).index, neighbor]
#             this_oriented_path, this_path_likelihood = 
#                 Mycelia.assess_path(this_path,
#                     graph.kmers,
#                     graph.counts,
#                     last(path).orientation,
#                     outgoing_edge_probabilities,
#                     incoming_edge_probabilities)
#             if this_path_likelihood > maximum_path_likelihood
#                 maximum_path_likelihood = this_path_likelihood
#                 maximum_likelihood_path = this_oriented_path
#             end
#         end
#         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
#             done = true
#         else
#             append!(path, maximum_likelihood_path[2:end])
#         end
#     end
#     return path
# end


# use default cost 1 for breadth first seach
# use default cost 0 for depth first search
function dijkstra_step!(graph, distances, arrival_paths, queue; default_cost = 1)
    current_kmer, cost = DataStructures.dequeue_pair!(queue)
    current_orientation = BioSequences.iscanonical(current_kmer)
    current_canonical_kmer = BioSequences.canonical(current_kmer)
    current_index = graph[current_canonical_kmer, :kmer]

    all_possible_neighbors = collect(BioSequences.neighbors(current_kmer))
    present_neighbors = 
        filter(neighbor -> 
            BioSequences.canonical(neighbor) in keys(MetaGraphs.get_prop(graph, :kmer_counts)), 
        all_possible_neighbors)

    true_neighbors =
        filter(neighbor -> 
            Graphs.has_edge(
                graph,
                graph[BioSequences.canonical(current_kmer), :kmer],
                graph[BioSequences.canonical(neighbor), :kmer]), 
            present_neighbors)
    
    if !isempty(true_neighbors)
        canonical_true_neighbors = BioSequences.canonical.(true_neighbors)
        true_neighbor_is_canonical = BioSequences.iscanonical.(true_neighbors)
        neighbor_indices = map(canonical_neighbor -> graph[canonical_neighbor, :kmer], canonical_true_neighbors)
        neighbor_weights = map(neighbor_i -> MetaGraphs.get_prop(graph, current_index, neighbor_i, :weight), neighbor_indices)


        # inverse probability!!
        neighbor_costs = 1 .- (neighbor_weights ./ sum(neighbor_weights))
        neighbor_costs .+= default_cost
        # afterwards, 100% likelihood edges cost 1
        # 50% likelihood costs 1.5
        # 0% likelihood is impossible and not included

        for (neighbor, cost) in zip(true_neighbors, neighbor_costs)    
            candidate_path = vcat(arrival_paths[current_kmer], current_kmer)
            candidate_path_cost = distances[current_kmer] + cost  
            if distances[neighbor] > candidate_path_cost
                arrival_paths[neighbor] = candidate_path
                DataStructures.enqueue!(queue, neighbor => cost)
                distances[neighbor] = candidate_path_cost
            end  
        end
    end
    return distances, arrival_paths, queue
end

function evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)

    lowest_cost = Inf
    optimal_path = Int[]

    for hit in hits
        reverse_complement_hit = BioSequences.reverse_complement(hit)
        total_cost = forward_distances[hit] + reverse_distances[reverse_complement_hit]
        if total_cost < lowest_cost    
            forward_path = vcat(forward_arrival_paths[hit], hit)   
            reverse_path = vcat(reverse_arrival_paths[reverse_complement_hit], reverse_complement_hit)
            reverse_path = reverse!(BioSequences.reverse_complement.(reverse_path))
            full_path = vcat(forward_path, reverse_path[2:end])

            lowest_cost = total_cost
            optimal_path = full_path
        end
    end

    return lowest_cost, optimal_path
end

function bidirectional_dijkstra(graph, a::T, b::T) where T <: BioSequences.AbstractMer
    forward_distances = DataStructures.DefaultDict{T, Float64}(Inf)
    forward_distances[a] = 0

    forward_arrival_paths = Dict{T, Vector{T}}()
    forward_arrival_paths[a] = Int[]

    forward_queue = DataStructures.PriorityQueue{T, Float64}()
    DataStructures.enqueue!(forward_queue, a, 0)

    reverse_b = BioSequences.reverse_complement(b)
    reverse_distances = DataStructures.DefaultDict{T, Float64}(Inf)
    reverse_distances[reverse_b] = 0

    reverse_arrival_paths = Dict{T, Vector{T}}()
    reverse_arrival_paths[reverse_b] = Int[]

    reverse_queue = DataStructures.PriorityQueue{T, Float64}()
    DataStructures.enqueue!(reverse_queue, reverse_b, 0)
    
    hits = Set{T}()
    
    while isempty(hits)
        dijkstra_step!(graph, forward_distances, forward_arrival_paths, forward_queue)
        dijkstra_step!(graph, reverse_distances, reverse_arrival_paths, reverse_queue)
        hits = intersect(
                    keys(forward_arrival_paths), 
                    BioSequences.reverse_complement.(keys(reverse_arrival_paths)))
    end
    lowest_cost, optimal_path = evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)
    return lowest_cost, optimal_path
end

# # simple point to point
# function dijkstra(graph, a::T, b::T) where {T <: BioSequences.AbstractMer}
#     distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     distances[a] = 0

#     arrival_paths = Dict{T, Vector{T}}()
#     arrival_paths[a] = Int[]

#     queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(queue, a, 0)

#     current_kmer, cost = a, Inf
#     while current_kmer != b
#         distances, arrival_paths, queue = dijkstra_step!(graph, distances, arrival_paths, queue)
#     end
#     shortest_path = vcat(arrival_paths[b], b)
#     return shortest_path, distances[b]    
# end

function dijkstra(graph, a::T, targets::Set{T}; search_strategy::Union{Symbol, Missing}=missing) where {T <: BioSequences.AbstractMer}
    if ismissing(search_strategy)
        # local breadth first search to find the single target
        # a-b shortest path
        if length(targets) == 1
            search_strategy = :BFS
        # a - no targets - go find the ends
        # a - series of any target
        else
            search_strategy = :DFS
        end
    end 
    if search_strategy == :BFS
        default_cost = 1
    elseif search_strategy == :DFS
        default_cost = 0
    else
        error("please specify either :BFS (breadth-first search) or :DFS (depth-first search)")
    end
    
    distances = DataStructures.DefaultDict{T, Float64}(Inf)
    distances[a] = 0

    arrival_paths = Dict{T, Vector{T}}()
    arrival_paths[a] = Int[]

    queue = DataStructures.PriorityQueue{T, Float64}()
    DataStructures.enqueue!(queue, a, 0)
    
    while !isempty(queue) && !(first(DataStructures.peek(queue)) in targets)
        distances, arrival_paths, queue = 
            Mycelia.dijkstra_step!(graph, distances, arrival_paths, queue, default_cost = default_cost)
    end
    if !isempty(queue)
        discovered_destination = DataStructures.dequeue!(queue)
        @assert discovered_destination in targets
        shortest_path = vcat(arrival_paths[discovered_destination], discovered_destination)
        return shortest_path, distances[discovered_destination]
    else
        longest_path = Int[]
        distance = Inf
        for (destination, arrival_path) in arrival_paths
            this_distance = distances[destination]
            this_path = vcat(arrival_path, destination)
            if (length(this_path) > length(longest_path)) && (this_distance <= distance)
                longest_path = this_path
                distance = this_distance
            end
            
        end
        return longest_path, distance
    end
end