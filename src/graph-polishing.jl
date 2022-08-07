# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simple_polish_fastq(simple_kmer_graph, fastq_file; min_depth=3)
#     solid_vertices = filter(v -> simple_kmer_graph.vprops[v][:weight] >= min_depth, Graphs.vertices(simple_kmer_graph))
#     filtered_simple_kmer_graph, vertex_map = Graphs.induced_subgraph(simple_kmer_graph, solid_vertices)
# #     display(simple_kmer_graph)
# #     display(filtered_simple_kmer_graph)
#     kmers = sort!(graph_to_kmers(filtered_simple_kmer_graph))
    
#     old_kmers = graph_to_kmers(simple_kmer_graph)
    
#     k = filtered_simple_kmer_graph.gprops[:k]
    
#     polished_fastq_file = replace(fastq_file, ".fastq" => ".k$k.d$(min_depth).fastq")

#     transition_probabilities = initialize_transition_probabilities(filtered_simple_kmer_graph)
#     state_likelihoods = [Float64(filtered_simple_kmer_graph.vprops[v][:weight]) for v in Graphs.vertices(filtered_simple_kmer_graph)]
#     state_likelihoods ./= sum(state_likelihoods)
    
# #     @info "counting the number of records to establish runtime estimate"
#     number_of_records = 0
#     for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
#         number_of_records += 1
#     end
#     progress_bar = ProgressMeter.Progress(number_of_records, 1)

#     fastq_reader = FASTX.FASTQ.Reader(open(fastq_file))
#     fastq_writer = FASTX.FASTQ.Writer(open(polished_fastq_file, "w"))

#     for fastx_record in fastq_reader
#         ProgressMeter.next!(progress_bar)
#         bubble_start = 0
#         updated_path = Vector{Pair{Int, Bool}}()
#     #     @show FASTX.sequence(fastx_record)
#         for (i, kmer) in enumerate(BioSequences.each(BioSequences.BigDNAMer{k}, FASTX.sequence(fastx_record)))
#             canonical_kmer = min(kmer.fw, kmer.bw)
#             orientation = canonical_kmer == kmer.fw
#             kmer_index_range = searchsorted(kmers, canonical_kmer)
            
#             old_kmer_index = searchsortedfirst(old_kmers, canonical_kmer)
            
# #             FASTX.identifier(fastx_record) == "4" && @show kmer_index_range
# #             FASTX.identifier(fastx_record) == "4" && @show kmers[kmer_index_range]
# #             FASTX.identifier(fastx_record) == "4" && @show old_kmers[old_kmer_index]
#             kmer_is_solid = !isempty(kmer_index_range)
# #             FASTX.identifier(fastx_record) == "4" && @show canonical_kmer
# #             FASTX.identifier(fastx_record) == "4" &&  @show kmer_is_solid

#             if kmer_is_solid
#                 kmer_index = first(kmer_index_range)
# #                 FASTX.identifier(fastx_record) == "4" && @show filtered_simple_kmer_graph.vprops[kmer_index], simple_kmer_graph.vprops[old_kmer_index]
#             else
#                 kmer_index = 0
#             end

#             in_bubble = bubble_start > 0

#             if !kmer_is_solid
#                 if !in_bubble
# #                     FASTX.identifier(fastx_record) == "4" && @show "starting a bubble"
#                     bubble_start = i
#                 else
# #                     FASTX.identifier(fastx_record) == "4" && @show "continuing in a bubble"
#                 end
#             else
#                 if !in_bubble
# #                     FASTX.identifier(fastx_record) == "4" && @show "pushing solid kmer to updated path"
#                     push!(updated_path, kmer_index => orientation)
#                 else
#                     if bubble_start == 1
# #                         FASTX.identifier(fastx_record) == "4" && @show "ending an opening bubble"
#                         # we're in a bubble that began at the beginning of the read
#                         # we'll do nothing and just remove this
#                         # equivalent to tip clipping
# #                         FASTX.identifier(fastx_record) == "4" && @show "pushing solid kmer to updated path"
#                         push!(updated_path, kmer_index => orientation)
#                         bubble_start = 0
#                     else
# #                         FASTX.identifier(fastx_record) == "4" && @show "found end of an anchored bubble -- correcting"
#                         source_vertex, source_orientation = last(updated_path)
#                         destination_vertex, destination_orientation = kmer_index, orientation                

#                         shortest_paths = Graphs.yen_k_shortest_paths(
#                             filtered_simple_kmer_graph,
#                             source_vertex,
#                             destination_vertex,
#                             Graphs.weights(filtered_simple_kmer_graph),
#                             10).paths

#                         if isempty(shortest_paths)
#                             @show source_vertex, destination_vertex
#                             @warn "no valid alternate paths found for $(FASTX.identifier(fastx_record)), continuing"
#                             break
# #                             @show fastx_record
# #                             error("try increasing min_depth")
#                         end
#                         candidate_path_probabilities = ones(length(shortest_paths))
#                         oriented_candidate_paths = [
#                             [last(updated_path)] for i in 1:length(shortest_paths)
#                         ]

#                         for (i, candidate_path) in enumerate(shortest_paths)
#                             for dest_vertex in candidate_path[2:end]
#                                 source_vertex, source_orientation = last(oriented_candidate_paths[i])
#                                 candidate_path_probabilities[i] *= transition_probabilities[source_orientation][source_vertex, dest_vertex]
#                                 candidate_path_probabilities[i] *= state_likelihoods[dest_vertex]
#                                 if candidate_path_probabilities[i] > 0
#                                     edge = Graphs.Edge(source_vertex, dest_vertex)
#                                     destination_orientation = 
#                                     first(
#                                         filter(o -> o.source_orientation == source_orientation,
#                                             filtered_simple_kmer_graph.eprops[edge][:orientations])).destination_orientation
#                                     push!(oriented_candidate_paths[i], (dest_vertex => destination_orientation))
#                                 else
#                                     break # this path is no good, evaluate the next
#                                 end
#                             end
#                         end
#                         non_zero_indices = findall(p -> p > 0, candidate_path_probabilities)
#                         if isempty(non_zero_indices)
# #                             @show candidate_path_probabilities
#                             @warn "no resolution found for read $(FASTX.identifier(fastx_record))"
#                             break
# #                             error("no valid alternate path probabilities")
# #                             error("try increasing min_depth?")
#                         end

#                         candidate_path_probabilities = candidate_path_probabilities[non_zero_indices]
#                         oriented_candidate_paths = oriented_candidate_paths[non_zero_indices]

#                         # offset is for debugging
#                         # make sure that anchors on both sides are the same
#                         offset = 0
#                         observed_sequence = FASTX.sequence(fastx_record)[bubble_start+k-1-offset:i-1+offset]                    
#                         for (i, oriented_candidate_path) in enumerate(oriented_candidate_paths)
#                             candidate_sequence = oriented_path_to_sequence(
#                                 filtered_simple_kmer_graph, 
#                                 oriented_candidate_path)
#                             candidate_sequence = candidate_sequence[k+1-offset:end-k+offset]
#                             alignment_result = BioAlignments.pairalign(
#                                 BioAlignments.LevenshteinDistance(),
#                                 candidate_sequence,
#                                 observed_sequence)
# #                             @show alignment_result
#     #                         @show alignment_result
#                             average_error_rate = Statistics.mean(q_value_to_error_rate.(FASTX.quality(fastx_record)))
#                             for error in 1:alignment_result.value
#                                 candidate_path_probabilities[i] *= average_error_rate
#                             end
#                             for match in 1:BioAlignments.count_matches(alignment_result.aln)
#                                 candidate_path_probabilities[i] *= (1 - average_error_rate)
#                             end
#                         end

#                         chosen_replacement = StatsBase.sample(oriented_candidate_paths, StatsBase.weights(candidate_path_probabilities))

#                         for i in 2:length(chosen_replacement)
#                             oriented_state = chosen_replacement[i]
#                             push!(updated_path, oriented_state)
#                         end
#                         bubble_start = 0
#                     end
#                 end
#             end
#         end
#     #     @show updated_path
#         sequence = oriented_path_to_sequence(filtered_simple_kmer_graph, updated_path)
#         alignment_result = BioAlignments.pairalign(
#             BioAlignments.LevenshteinDistance(),
#             sequence,
#             FASTX.sequence(fastx_record))
# #         if alignment_result.value > 0
# # #             @show alignment_result
# #         end
#         quality = StatsBase.sample(FASTX.quality(fastx_record), length(sequence), replace=true, ordered=true)
#         description =  join(filter(!isempty, (FASTX.description(fastx_record), "k$k.d$(min_depth)")), '.')
#         identifier = FASTX.identifier(fastx_record)
#         new_record = FASTX.FASTQ.Record(identifier, description, sequence, quality)
#         write(fastq_writer, new_record)
#     end
#     close(fastq_reader)
#     close(fastq_writer)
#     return polished_fastq_file
# end



# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function assess_observations(graph::KmerGraph{KMER_TYPE}, observations, error_rate; verbose = isinteractive()) where {KMER_TYPE}
#     k = last(KMER_TYPE.parameters)
#     total_edits_accepted = 0
#     total_bases_evaluated = 0
#     reads_processed = 0
#     maximum_likelihood_observations = Vector{BioSequences.LongDNASeq}(undef, length(observations))
#     for (observation_index, observation) in enumerate(observations)
#         if length(observation) >= k
#             optimal_path, edit_distance, relative_likelihood = viterbi_maximum_likelihood_path(graph, observation, error_rate)
#             maximum_likelihood_observation = oriented_path_to_sequence(optimal_path, graph.kmers)
#             maximum_likelihood_observations[observation_index] = maximum_likelihood_observation
#             reads_processed += 1
#             total_bases_evaluated += length(observation)
#             total_edits_accepted += edit_distance
#         else
#             maximum_likelihood_observations[observation_index] = observation
#         end
#     end
#     inferred_error_rate = round(total_edits_accepted / total_bases_evaluated, digits = 3)
#     if verbose
#         display("reads_processed = $(reads_processed)")
#         display("total_edits_accepted = $(total_edits_accepted)")
#         display("inferred_error_rate = $(inferred_error_rate)")
#     end
#     if total_edits_accepted == 0
#         has_converged = true
#     else
#         has_converged = false
#     end
#     return maximum_likelihood_observations, has_converged
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function iterate_until_convergence(ks, observations, error_rate)
#     graph = nothing
#     for k in ks
#         graph = KmerGraph(BioSequences.DNAMer{k}, observations)
#         observations, has_converged = assess_observations(graph, observations, error_rate; verbose = verbose)
#     end
#     return graph, observations
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function clip_low_coverage_tips(graph, observations)
#     connected_components = Graphs.connected_components(graph.graph)
#     vertices_to_keep = Int[]
#     for connected_component in connected_components
        
#         component_coverage = graph.counts[connected_component]
#         median = Statistics.median(component_coverage)
#         standard_deviation = Statistics.std(component_coverage)
        
#         for vertex in connected_component
#             keep = true
#             if Graphs.degree(graph.graph, vertex) == 1
#                 this_coverage = graph.counts[vertex]
#                 is_low_coverage = (graph.counts[vertex] == 1) || 
#                                     (median-this_coverage) > (3*standard_deviation)
#                 if is_low_coverage
#                     keep = false
#                 end
#             end
#             if keep
#                 push!(vertices_to_keep, vertex)
#             end
#         end
#     end
    
#     KmerType = first(typeof(graph).parameters)
#     pruned_graph = KmerGraph(KmerType, observations, graph.kmers[vertices_to_keep], graph.counts[vertices_to_keep])
    
#     return pruned_graph
# end