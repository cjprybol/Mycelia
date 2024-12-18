# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_edge_probabilities(graph, strand)
#     kmers = graph.gprops[:kmers]
#     outgoing_edge_probabilities = SparseArrays.spzeros(length(kmers), length(kmers))
    
#     for (kmer_index, kmer) in enumerate(kmers)
#         if !strand
#             kmer = BioSequences.reverse_complement(kmer)
#         end
        
#         downstream_neighbor_indices = Int[]
#         for neighbor in BioSequences.neighbors(kmer)
#             index = get_kmer_index(kmers, BioSequences.canonical(neighbor))
#             # kmer must be in our dataset and there must be a connecting edge
#             if !isnothing(index) && Graphs.has_edge(graph, ordered_edge(kmer_index, index))
#                 push!(downstream_neighbor_indices, index)
#             end
#         end
#         sort!(unique!(downstream_neighbor_indices))
        
#         downstream_edge_weights = Int[
#             length(get(graph.edge_evidence, ordered_edge(kmer_index, neighbor_index), EdgeEvidence[])) for neighbor_index in downstream_neighbor_indices
#         ]
        
#         non_zero_indices = downstream_edge_weights .> 0
#         downstream_neighbor_indices = downstream_neighbor_indices[non_zero_indices]
#         downstream_edge_weights = downstream_edge_weights[non_zero_indices]
        
#         downstream_edge_likelihoods = downstream_edge_weights ./ sum(downstream_edge_weights)
        
#         for (neighbor_index, likelihood) in zip(downstream_neighbor_indices, downstream_edge_likelihoods)
#             outgoing_edge_probabilities[kmer_index, neighbor_index] = likelihood
#         end
#     end
#     return outgoing_edge_probabilities
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_edge_probabilities(graph)
#     outgoing_edge_probabilities = determine_edge_probabilities(graph, true)
#     incoming_edge_probabilities = determine_edge_probabilities(graph, false)
#     return outgoing_edge_probabilities, incoming_edge_probabilities
# end





# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function reverse_oriented_path(oriented_path)
# #     reversed_path = copy(oriented_path)
# #     for (index, state) in enumerate(oriented_path)
# #         reversed_path[index] = OrientedKmer(index = state.index, orientation = !state.orientation)
# #     end
# #     return reverse!(reversed_path)
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function maximum_likelihood_walk(graph, connected_component)
# #     max_count = maximum(graph.counts[connected_component])
# #     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
# #     initial_node_index = rand(max_count_indices)
# #     initial_node = connected_component[initial_node_index]
# #     outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
# #     forward_walk = maximum_likelihood_walk(graph, [OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
# #     reverse_walk = maximum_likelihood_walk(graph, [OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
# #     reversed_reverse_walk = reverse!(
# #         [
# #             OrientedKmer(index = oriented_kmer.index, orientation = oriented_kmer.orientation)
# #             for oriented_kmer in reverse_walk[2:end]
# #         ]
# #         )
# #     full_path = [reversed_reverse_walk..., forward_walk...]
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function maximum_likelihood_walk(graph, path::Vector{OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
# #     done = false
# #     while !done
# #         maximum_path_likelihood = 0.0
# #         maximum_likelihood_path = Vector{OrientedKmer}()
# #         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
# #             this_path = [last(path).index, neighbor]
# #             this_oriented_path, this_path_likelihood = 
# #                 assess_path(this_path,
# #                     graph.kmers,
# #                     graph.counts,
# #                     last(path).orientation,
# #                     outgoing_edge_probabilities,
# #                     incoming_edge_probabilities)
# #             if this_path_likelihood > maximum_path_likelihood
# #                 maximum_path_likelihood = this_path_likelihood
# #                 maximum_likelihood_path = this_oriented_path
# #             end
# #         end
# #         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
# #             done = true
# #         else
# #             append!(path, maximum_likelihood_path[2:end])
# #         end
# #     end
# #     return path
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function take_a_walk(graph, connected_component)
# #     max_count = maximum(graph.counts[connected_component])
# #     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
# #     initial_node_index = rand(max_count_indices)
# #     initial_node = connected_component[initial_node_index]
# #     outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
    
# #     # walk forwards from the initial starting node
# #     forward_walk = take_a_walk(graph, [OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
# #     # walk backwards from the initial starting node
# #     reverse_walk = take_a_walk(graph, [OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
# #     # we need to reverse everything to re-orient against the forward walk
# #     reverse_walk = reverse_oriented_path(reverse_walk)
    
# #     # also need to drop the last node, which is equivalent to the first node of the 
# #     @assert last(reverse_walk) == first(forward_walk)
# #     full_path = [reverse_walk[1:end-1]..., forward_walk...]
# # #     @show full_path
# # end


# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function take_a_walk(graph, path::Vector{OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
# #     done = false
# #     while !done
# #         maximum_path_likelihood = 0.0
# #         maximum_likelihood_path = Vector{OrientedKmer}()
# #         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
# #             this_path = [last(path).index, neighbor]
# #             this_oriented_path, this_path_likelihood = 
# #                 assess_path(this_path,
# #                     graph.kmers,
# #                     graph.counts,
# #                     last(path).orientation,
# #                     outgoing_edge_probabilities,
# #                     incoming_edge_probabilities)
# #             if this_path_likelihood > maximum_path_likelihood
# #                 maximum_path_likelihood = this_path_likelihood
# #                 maximum_likelihood_path = this_oriented_path
# #             end
# #         end
# #         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
# #             done = true
# #         else
# #             append!(path, maximum_likelihood_path[2:end])
# #         end
# #     end
# #     return path
# # end


# # use default cost 1 for breadth first seach
# # use default cost 0 for depth first search
# function dijkstra_step!(graph, distances, arrival_paths, queue; default_cost = 1)
#     current_kmer, cost = DataStructures.dequeue_pair!(queue)
#     current_orientation = BioSequences.iscanonical(current_kmer)
#     current_canonical_kmer = BioSequences.canonical(current_kmer)
#     current_index = graph[current_canonical_kmer, :kmer]

#     present_neighbors = Vector{typeof(current_kmer)}()
#     for neighbor in BioSequences.neighbors(current_kmer)
#         try
#             graph[BioSequences.canonical(neighbor), :kmer]
#             push!(present_neighbors, neighbor)
#         catch
#             continue
#         end
#     end
    
#     true_neighbors =
#         filter(neighbor -> 
#             Graphs.has_edge(
#                 graph,
#                 graph[BioSequences.canonical(current_kmer), :kmer],
#                 graph[BioSequences.canonical(neighbor), :kmer]), 
#             present_neighbors)
    
#     if !isempty(true_neighbors)
#         canonical_true_neighbors = BioSequences.canonical.(true_neighbors)
#         true_neighbor_is_canonical = BioSequences.iscanonical.(true_neighbors)
#         neighbor_indices = map(canonical_neighbor -> graph[canonical_neighbor, :kmer], canonical_true_neighbors)
#         neighbor_weights = map(neighbor_i -> MetaGraphs.get_prop(graph, current_index, neighbor_i, :weight), neighbor_indices)


#         # inverse probability!!
#         neighbor_costs = 1 .- (neighbor_weights ./ sum(neighbor_weights))
#         neighbor_costs .+= default_cost
#         # afterwards, 100% likelihood edges cost 1
#         # 50% likelihood costs 1.5
#         # 0% likelihood is impossible and not included

#         for (neighbor, cost) in zip(true_neighbors, neighbor_costs)    
#             candidate_path = vcat(arrival_paths[current_kmer], current_kmer)
#             candidate_path_cost = distances[current_kmer] + cost  
#             if distances[neighbor] > candidate_path_cost
#                 arrival_paths[neighbor] = candidate_path
#                 DataStructures.enqueue!(queue, neighbor => cost)
#                 distances[neighbor] = candidate_path_cost
#             end  
#         end
#     end
#     return distances, arrival_paths, queue
# end

# function evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)

#     lowest_cost = Inf
#     optimal_path = Int[]

#     for hit in hits
#         reverse_complement_hit = BioSequences.reverse_complement(hit)
#         total_cost = forward_distances[hit] + reverse_distances[reverse_complement_hit]
#         if total_cost < lowest_cost    
#             forward_path = vcat(forward_arrival_paths[hit], hit)   
#             reverse_path = vcat(reverse_arrival_paths[reverse_complement_hit], reverse_complement_hit)
#             reverse_path = reverse!(BioSequences.reverse_complement.(reverse_path))
#             full_path = vcat(forward_path, reverse_path[2:end])

#             lowest_cost = total_cost
#             optimal_path = full_path
#         end
#     end

#     return lowest_cost, optimal_path
# end

# function bidirectional_dijkstra(graph, a::T, b::T) where T <: Kmers.Kmer
#     forward_distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     forward_distances[a] = 0

#     forward_arrival_paths = Dict{T, Vector{T}}()
#     forward_arrival_paths[a] = Int[]

#     forward_queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(forward_queue, a, 0)

#     reverse_b = BioSequences.reverse_complement(b)
#     reverse_distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     reverse_distances[reverse_b] = 0

#     reverse_arrival_paths = Dict{T, Vector{T}}()
#     reverse_arrival_paths[reverse_b] = Int[]

#     reverse_queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(reverse_queue, reverse_b, 0)
    
#     hits = Set{T}()
    
#     while isempty(hits)
#         dijkstra_step!(graph, forward_distances, forward_arrival_paths, forward_queue)
#         dijkstra_step!(graph, reverse_distances, reverse_arrival_paths, reverse_queue)
#         hits = intersect(
#                     keys(forward_arrival_paths), 
#                     BioSequences.reverse_complement.(keys(reverse_arrival_paths)))
#     end
#     lowest_cost, optimal_path = evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)
#     return lowest_cost, optimal_path
# end

# # # simple point to point
# # function dijkstra(graph, a::T, b::T) where {T <: Kmers.Kmer}
# #     distances = DataStructures.DefaultDict{T, Float64}(Inf)
# #     distances[a] = 0

# #     arrival_paths = Dict{T, Vector{T}}()
# #     arrival_paths[a] = Int[]

# #     queue = DataStructures.PriorityQueue{T, Float64}()
# #     DataStructures.enqueue!(queue, a, 0)

# #     current_kmer, cost = a, Inf
# #     while current_kmer != b
# #         distances, arrival_paths, queue = dijkstra_step!(graph, distances, arrival_paths, queue)
# #     end
# #     shortest_path = vcat(arrival_paths[b], b)
# #     return shortest_path, distances[b]    
# # end

# function dijkstra(graph, a::T, targets::Set{T}; search_strategy::Union{Symbol, Missing}=missing) where {T <: Kmers.Kmer}
#     if ismissing(search_strategy)
#         # local breadth first search to find the single target
#         # a-b shortest path
#         if length(targets) == 1
#             search_strategy = :BFS
#         # a - no targets - go find the ends
#         # a - series of any target
#         else
#             search_strategy = :DFS
#         end
#     end 
#     if search_strategy == :BFS
#         default_cost = 1
#     elseif search_strategy == :DFS
#         default_cost = 0
#     else
#         error("please specify either :BFS (breadth-first search) or :DFS (depth-first search)")
#     end
    
#     distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     distances[a] = 0

#     arrival_paths = Dict{T, Vector{T}}()
#     arrival_paths[a] = Int[]

#     queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(queue, a, 0)
    
#     while !isempty(queue) && !(first(DataStructures.peek(queue)) in targets)
#         distances, arrival_paths, queue = 
#             dijkstra_step!(graph, distances, arrival_paths, queue, default_cost = default_cost)
#     end
#     if !isempty(queue)
#         discovered_destination = DataStructures.dequeue!(queue)
#         @assert discovered_destination in targets
#         shortest_path = vcat(arrival_paths[discovered_destination], discovered_destination)
#         return shortest_path, distances[discovered_destination]
#     else
#         longest_path = Int[]
#         distance = Inf
#         for (destination, arrival_path) in arrival_paths
#             this_distance = distances[destination]
#             this_path = vcat(arrival_path, destination)
#             if (length(this_path) > length(longest_path)) && (this_distance <= distance)
#                 longest_path = this_path
#                 distance = this_distance
#             end
            
#         end
#         return longest_path, distance
#     end
# end

# function update_remaining_targets(current_walk::AbstractVector{T}, remaining_targets::AbstractSet{T}) where T <: Kmers.Kmer
#     # assess whether targets have been hit in the canonical space
#     remaining_targets = setdiff(BioSequences.canonical.(remaining_targets), BioSequences.canonical.(current_walk))
#     # blow back out into forward and reverse_complement space
#     remaining_targets = Set{T}(vcat(remaining_targets, BioSequences.reverse_complement.(remaining_targets)))
#     return remaining_targets
# end

# function assess_downstream_weight(graph, kmer)
#     # here we look to see if walking forward or backward from the initial node gets us to heavier weight options
#     score = 0
#     for neighbor in BioSequences.neighbors(kmer)
#         try
#             score += MetaGraphs.get_prop(graph, graph[BioSequences.canonical(neighbor), :kmer], :count)
#         catch
#             continue
#         end
#     end
#     return score
# end

# # vertices should either be entire graph (by default) or a connected component
# # if people want to work on just the connected component, let them induce a subgraph
# function find_graph_core(graph; seed=rand(Int))
    
#     Random.seed!(seed)
    
#     T = typeof(MetaGraphs.get_prop(graph, 1, :kmer))
    
#     # targets = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
#     # sortperm
    
#     # targets = [
#     #         MetaGraphs.get_prop(graph, v_index, :count) for v_index in 
#     #         sortperm([MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)], rev=true)[1:Int(ceil(Graphs.nv(graph)/10))]]
    
#     targets = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph) if MetaGraphs.get_prop(graph, v, :count) > 1]
#             # sortperm([MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)], rev=true)[1:Int(ceil(Graphs.nv(graph)/10))]]
    
#     @show length(targets)
#     starting_kmer = first(targets)
#     max_degree = 0
#     for node in targets
#         node_degree = Graphs.degree(graph, graph[node, :kmer])
#         if node_degree > max_degree
#             max_degree = node_degree
#             starting_kmer = node
#         end
#     end
        
#     current_walk = [starting_kmer]
#     prior_walk_length = length(current_walk)
#     remaining_targets = update_remaining_targets(current_walk, Set(targets))
#     done = isempty(remaining_targets)
    
#     while !done
#         # here we look to see if walking forward or backward from the current ends gets us to heavier weight options
#         # we want to prioritize walks toward higher coverage nodes
#         forward_score = assess_downstream_weight(graph, last(current_walk))
#         reverse_score = assess_downstream_weight(graph, BioSequences.reverse_complement(first(current_walk)))
#         if reverse_score > forward_score
#             current_walk = reverse(BioSequences.reverse_complement.(current_walk))
#         end
        
#         forward_source = last(current_walk)
#         forward_walk, forward_distance = dijkstra(graph, forward_source, remaining_targets, search_strategy=:DFS)
#         current_walk = vcat(current_walk, forward_walk[2:end])
#         remaining_targets = update_remaining_targets(current_walk, remaining_targets)
#         if isempty(remaining_targets)
#             done = true
#         else
#             reverse_source = BioSequences.reverse_complement(first(current_walk))
#             reverse_walk, reverse_distance = dijkstra(graph, reverse_source, remaining_targets, search_strategy=:DFS)
#             current_walk = vcat(reverse(BioSequences.reverse_complement.(reverse_walk))[1:end-1], current_walk)
#             remaining_targets = update_remaining_targets(current_walk, remaining_targets)
#         end
#         @show length(current_walk)
#         failed_this_expansion = length(current_walk) == prior_walk_length
#         prior_walk_length = length(current_walk)
#         if isempty(remaining_targets)
#             done = true
#         elseif failed_this_expansion
#             done = true
#         end
#     end
#     return current_walk
# end