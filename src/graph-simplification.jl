# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simplify_kmer_graph(kmer_graph)
#     @info "simplifying kmer graph"
#     @info "resolving untigs..."
#     @time untigs = resolve_untigs(kmer_graph)
#     @info "determining untig orientations..."
#     oriented_untigs = determine_oriented_untigs(kmer_graph, untigs)
#     simplified_graph = MetaGraphs.MetaDiGraph(length(oriented_untigs))
#     MetaGraphs.set_prop!(simplified_graph, :k, kmer_graph.gprops[:k])
#     @info "initializing graph node metadata"
#     for (vertex, untig) in enumerate(oriented_untigs)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :sequence, untig.sequence)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :path, untig.path)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :orientations, untig.orientations)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :evidence, untig.evidence)
#     end
    
#     # determine oriented edges of simplified graph
#     simplified_untigs = []
#     @info "creating simplified untigs to resolve connections"
#     for vertex in Graphs.vertices(simplified_graph)
#         in_kmer = simplified_graph.vprops[vertex][:path][1] => simplified_graph.vprops[vertex][:orientations][1]
#         out_kmer = simplified_graph.vprops[vertex][:path][end] => simplified_graph.vprops[vertex][:orientations][end]
#     #     @show vertex, in_kmer, out_kmer
#         push!(simplified_untigs, in_kmer => out_kmer)
#     end

#     @info "resolving connections"
#     ProgressMeter.@showprogress for (ui, u) in enumerate(simplified_untigs)
#         for (vi, v) in enumerate(simplified_untigs)
#     #         + => +
#             source_kmer_index, source_orientation = last(u)
#             destination_kmer_index, destination_orientation = first(v)
#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! + +"

#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         + => -
#             source_kmer_index, source_orientation = last(u)
#             destination_kmer_index, destination_orientation = last(v)
#             destination_orientation = !destination_orientation

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! + -"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         - => +
#             source_kmer_index, source_orientation = first(u)
#             source_orientation = !source_orientation
#             destination_kmer_index, destination_orientation = first(v)

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! - +"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         - => -
#             source_kmer_index, source_orientation = first(u)
#             source_orientation = !source_orientation
#             destination_kmer_index, destination_orientation = last(v)
#             destination_orientation = !destination_orientation

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! - -"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#         end
#     end
#     return simplified_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simplify_kmer_graph(kmer_graph)
#     @info "simplifying kmer graph"
#     @info "resolving untigs..."
#     @time untigs = Mycelia.resolve_untigs(kmer_graph)

#     @info "determining untig orientations..."
#     @time oriented_untigs = Mycelia.determine_oriented_untigs(kmer_graph, untigs)

#     simplified_graph = MetaGraphs.MetaDiGraph(length(oriented_untigs))
#     MetaGraphs.set_prop!(simplified_graph, :k, kmer_graph.gprops[:k])
#     @info "initializing graph node metadata"
#     ProgressMeter.@showprogress for (vertex, untig) in enumerate(oriented_untigs)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :sequence, untig.sequence)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :path, untig.path)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :orientations, untig.orientations)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :weight, untig.weight)
#     end

#     # determine oriented edges of simplified graph
#     simplified_untigs = Vector{Pair{Pair{Int64,Bool},Pair{Int64,Bool}}}(undef, length(Graphs.vertices(simplified_graph)))
#     @info "creating simplified unitgs to help resolve connections"
#     # use a pre-allocated array here to speed up
#     ProgressMeter.@showprogress for vertex in Graphs.vertices(simplified_graph)
#         in_kmer = simplified_graph.vprops[vertex][:path][1] => simplified_graph.vprops[vertex][:orientations][1]
#         out_kmer = simplified_graph.vprops[vertex][:path][end] => simplified_graph.vprops[vertex][:orientations][end]
#     #     @show vertex, in_kmer, out_kmer
#         simplified_untigs[vertex] = in_kmer => out_kmer
#     #     push!(simplified_untigs, )
#     end

#     # make a dictionary mapping endcap to oriented_untig index

#     end_mer_map = Dict()
#     ProgressMeter.@showprogress for (i, oriented_untig) in enumerate(oriented_untigs)
#         end_mer_map[first(oriented_untig.path)] = i
#         end_mer_map[last(oriented_untig.path)] = i
#     end

#     ProgressMeter.@showprogress for (untig_index, oriented_untig) in enumerate(oriented_untigs)
#     #     @show untig_index
#         true_in_overlap = oriented_untig.sequence[1:simplified_graph.gprops[:k]-1]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[1])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[2])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_out_overlap = neighboring_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_true_out_overlap == true_in_overlap
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = true => true
#                 o = (source_orientation = true, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_out_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_false_out_overlap == true_in_overlap        
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = false => true
#                 o = (source_orientation = false, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         true_out_overlap = oriented_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[end])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[end-1])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_in_overlap = neighboring_untig.sequence[1:simplified_graph.gprops[:k]-1]
#             if true_out_overlap == neighbor_true_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = true => true
#                 o = (source_orientation = true, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_in_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[1:simplified_graph.gprops[:k]-1]
#             if true_out_overlap == neighbor_false_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = true => false
#                 o = (source_orientation = true, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         false_in_overlap = BioSequences.reverse_complement(oriented_untig.sequence)[1:simplified_graph.gprops[:k]-1]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[end])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[end-1])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_out_overlap = neighboring_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_true_out_overlap == false_in_overlap
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = true => false
#                 o = (source_orientation = true, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_out_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_false_out_overlap == false_in_overlap        
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = false => false
#                 o = (source_orientation = false, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         false_out_overlap = BioSequences.reverse_complement(oriented_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[1])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[2])
#         end

#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_in_overlap = neighboring_untig.sequence[1:simplified_graph.gprops[:k]-1]
#             if false_out_overlap == neighbor_true_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = false => true
#                 o = (source_orientation = false, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_in_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[1:simplified_graph.gprops[:k]-1]
#             if false_out_overlap == neighbor_false_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = false => false
#                 o = (source_orientation = false, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 Mycelia.set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end
#     end
#     return simplified_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
function determine_oriented_untigs(kmer_graph, untigs)
    oriented_untigs = []
    for path in untigs
        sequence = BioSequences.LongDNASeq(kmer_graph.vprops[first(path)][:kmer])
        if length(path) == 1
            orientations = [true]
        elseif length(path) > 1
            initial_edge = Graphs.Edge(path[1], path[2])
            initial_orientation = first(kmer_graph.eprops[initial_edge][:orientations]).source_orientation
            orientations = [initial_orientation]
            if !initial_orientation
                sequence = BioSequences.reverse_complement(sequence)
            end

            for (src, dst) in zip(path[1:end-1], path[2:end])
                edge = Graphs.Edge(src, dst)
                destination = BioSequences.LongDNASeq(kmer_graph.vprops[edge.dst][:kmer])
                destination_orientation = first(kmer_graph.eprops[edge][:orientations]).destination_orientation
                push!(orientations, destination_orientation)
                if !destination_orientation
                    destination = BioSequences.reverse_complement(destination)
                end
                sequence_suffix = sequence[end-length(destination)+2:end]
                destination_prefix = destination[1:end-1]
                @assert sequence_suffix == destination_prefix
                push!(sequence, destination[end])
            end
        end

        oriented_untig = 
        (
            sequence = BioSequences.canonical(sequence),
            path = BioSequences.iscanonical(sequence) ? path : reverse(path),
            orientations = BioSequences.iscanonical(sequence) ? orientations : reverse(.!orientations),
            weight = Statistics.median([kmer_graph.vprops[v][:weight] for v in path])
        )

        push!(oriented_untigs, oriented_untig)
    end
    return oriented_untigs
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
function resolve_untigs(kmer_graph)
    untigs = Vector{Int}[]
    visited = falses(Graphs.nv(kmer_graph))
    first_unvisited = findfirst(!, visited)
    while first_unvisited != nothing
        forward_walk = oriented_unbranching_walk(kmer_graph, first_unvisited, true)
        reverse_walk = oriented_unbranching_walk(kmer_graph, first_unvisited, false)
        inverted_reverse_walk = [Graphs.Edge(e.dst, e.src) for e in reverse(reverse_walk)]
        edges = vcat(inverted_reverse_walk, forward_walk)
        if isempty(edges)
            untig = [first_unvisited]
        else
            untig = vcat([first(edges).src], [edge.dst for edge in edges])
        end
        push!(untigs, untig)
        for vertex in untig
            visited[vertex] = true
        end
        first_unvisited = findfirst(!, visited)
    end
    return untigs
end

function apply_kmedoids_treshold(graph)
    kmer_counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]

    kmer_counts_histogram = sort(collect(StatsBase.countmap(values(kmer_counts))), by=x->x[1])

#     scale = 250
#     p = Mycelia.plot_kmer_frequency_spectra(values(kmer_counts), size=(2scale,scale), log_scale=log2, title="kmer frequencies")
#     display(p)

#     p = StatsPlots.scatter(log2.(first.(kmer_counts_histogram)))
#     display(p)

    kmer_depth_of_coverage_bins = log2.(first.(kmer_counts_histogram))

    distance_matrix = zeros((length(kmer_depth_of_coverage_bins), length(kmer_depth_of_coverage_bins)))
    for (row, depth_of_coverage_bin_1) in enumerate(kmer_depth_of_coverage_bins)
        for (col, depth_of_coverage_bin_2) in enumerate(kmer_depth_of_coverage_bins)
            distance = abs(depth_of_coverage_bin_1 - depth_of_coverage_bin_2)
            distance_matrix[row, col] = distance
        end
    end
    distance_matrix

    # max out k at the same max k we use for DNAMers
    max_k = min(length(kmer_depth_of_coverage_bins), 63)
    ks = Primes.primes(2, max_k)
    ys = map(k ->
                Statistics.mean(Statistics.mean(Clustering.silhouettes(Clustering.kmedoids(distance_matrix, k), distance_matrix)) for i in 1:100),
                ks)

    p = StatsPlots.plot(ks, ys, label="silhouette score", ylabel = "silhouette score", xlabel = "number of clusters")
    display(p)

    ymax, ymax_index = findmax(ys)
    optimal_k = ks[ymax_index]
    clusterings = [Clustering.kmedoids(distance_matrix, optimal_k) for i in 1:10]
    max_value, max_value_index = findmax(clustering -> Statistics.mean(Clustering.silhouettes(clustering, distance_matrix)), clusterings)
    optimal_clustering = clusterings[max_value_index]
    # optimal_clustering.assignments
    min_medoid_value, min_medoid_index = findmin(optimal_clustering.medoids)
    indices_to_include = map(assignment -> assignment .!= min_medoid_index, optimal_clustering.assignments)
    # kmer_depth_of_coverage_bins
    threshold = Int(ceil(2^maximum(kmer_depth_of_coverage_bins[.!indices_to_include]))) + 1

    scale = 250
    p = Mycelia.plot_kmer_frequency_spectra(values(kmer_counts), log_scale = log2, size=(2scale,scale), title="kmer frequencies")
    StatsPlots.vline!(p, log2.([threshold]))
    display(p)

    # find all vertices with count > threshold
    vertices_to_keep = [v for v in Graphs.vertices(graph) if (MetaGraphs.get_prop(graph, v, :count) > threshold)]
    # induce subgraph
    induced_subgraph, vertex_map = Graphs.induced_subgraph(graph, vertices_to_keep)

    # set kmer as indexing prop
    MetaGraphs.set_indexing_prop!(induced_subgraph, :kmer)
    return induced_subgraph
end