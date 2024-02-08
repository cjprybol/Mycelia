# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function initialize_transition_probabilities(kmer_graph)
    
#     total_kmers = Graphs.nv(kmer_graph)
#     transition_likelihoods = Dict(
#         true => SparseArrays.spzeros(total_kmers, total_kmers),
#         false => SparseArrays.spzeros(total_kmers, total_kmers)
#     )

#     for edge in collect(Graphs.edges(kmer_graph))
# #         weight = length(kmer_graph.eprops[edge][:evidence])
#         weight = kmer_graph.eprops[edge][:weight]
#         for o in kmer_graph.eprops[edge][:orientations]
#             transition_likelihoods[o.source_orientation][edge.src, edge.dst] = weight
#         end
#     end

#     for source_orientation in (true, false)
#         for src in 1:total_kmers
#             transition_weights = transition_likelihoods[source_orientation][src, :]
#             total_weight = sum(transition_weights)
#             dsts, vals = SparseArrays.findnz(transition_weights)
#             for (dst, val) in zip(dsts, vals) 
#                 transition_likelihoods[source_orientation][src, dst] = val / total_weight
#             end
#             normalized_probability = sum(transition_likelihoods[source_orientation][src, :])
#             @assert isapprox(normalized_probability, 0) || isapprox(normalized_probability, 1)
#         end
#     end
#     return transition_likelihoods
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function set_initial_state_likelihoods!(
#         kmer_graph,
#         initial_state,
#         kmer_likelihoods,
#         error_rate,
#         state_likelihoods,
#         arrival_paths
#     )
#     for vertex in collect(Graphs.vertices(kmer_graph))
#         hidden_kmer = kmer_graph.vprops[vertex][:kmer]

#         fw_alignment = 
#             BioAlignments.pairalign(
#                 BioAlignments.LevenshteinDistance(), 
#                 initial_state.fw, 
#                 hidden_kmer)

#         fw_probability = kmer_likelihoods[vertex]

#         for match in 1:BioAlignments.count_matches(BioAlignments.alignment(fw_alignment))
#             fw_probability *= 1 - error_rate
#         end

#         for edit in 1:fw_alignment.value
#             fw_probability *= error_rate
#         end

#         bw_alignment = 
#             BioAlignments.pairalign(
#                 BioAlignments.LevenshteinDistance(),
#                 initial_state.bw,
#                 hidden_kmer)

#         bw_probability = kmer_likelihoods[vertex]

#         for match in 1:BioAlignments.count_matches(BioAlignments.alignment(bw_alignment))
#             bw_probability *= 1 - error_rate
#         end

#         for edit in 1:bw_alignment.value
#             bw_probability *= error_rate
#         end

#         if fw_probability > bw_probability
#             state_probability = fw_probability
#             state_orientation = true
#         elseif fw_probability < bw_probability
#             state_probability = bw_probability
#             state_orientation = false
#         else fw_probability == bw_probability
#             state_probability = fw_probability
#             state_orientation = missing
#         end
#         state_likelihoods[vertex, 1] = state_probability
#         arrival_paths[vertex, 1] = [vertex => state_orientation]
#     end
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_viterbi!(
#         current_state,
#         prior_state,
#         observed_nucleotide,
#         observed_quality_score,
#         observed_error_rate,
#         current_vertex,
#         prior_vertex,
#         state_likelihoods,
#         transition_likelihoods,
#         shortest_paths,
#         arrival_paths,
#         kmer_graph,
#         kmer_likelihoods
#         )
#     # if probability of prior state is lower than current probability, skip

# #     @show current_state
# #     @show prior_state
# #     @show current_vertex
# #     @show prior_vertex
    
    
#     current_state_likelihood = state_likelihoods[current_vertex, current_state]
#     prior_state_likelihood = state_likelihoods[prior_vertex, prior_state]

#     # if we already have a better possible path, skip calculating anything
#     if prior_state_likelihood < current_state_likelihood
# #         @show prior_state_likelihood < current_state_likelihood
#         return
#     end

#     # take shortest path and assume it's the maximum likelihood path
#     # this assumption seems fair because in an ideal situation
#     # we're just moving to an adjacent kmer
#     # and the shortest path and most likely path should be the same
#     shortest_path = shortest_paths[prior_vertex][current_vertex]
    
# #     no path & not considering insertion
#     if isempty(shortest_path) && (prior_vertex != current_vertex)
# #         @show "no path, skipping"
#         return
#     end
    
#     # if shortest path isn't viable, exit
#     if !isempty(shortest_path)
# #         @show "checking if path is viable"

#         terminal_orientation_prior_state = last(last(arrival_paths[prior_vertex, prior_state]))
# #         @show arrival_paths[prior_vertex, prior_state]
# #         @show "we were at vertex $(prior_vertex) in orientation $(terminal_orientation_prior_state)"
#         candidate_edge = Graphs.Edge(shortest_path[1], shortest_path[2])
                
#         if !ismissing(terminal_orientation_prior_state) && 
#             !any(o -> o.source_orientation == terminal_orientation_prior_state, kmer_graph.eprops[candidate_edge][:orientations])
            
# #             @show "no viable orientation matching edges detected between $(candidate_edge)"
# #             @show "full candidate path was $(shortest_path)"
# #             @show "orientation options were:"
# #             @show kmer_graph.eprops[candidate_edge][:orientations]
#             return
#         end
#     end
    
#     # zero step path - insertion in observed sequence relative to kmer graph
#     is_same_vertex = (current_vertex == prior_vertex)
#     has_edge = Graphs.has_edge(kmer_graph, Graphs.Edge(prior_vertex, current_vertex))
#     if is_same_vertex && has_edge
#         shortest_path = [prior_vertex, current_vertex]
#     end
    
#     if is_same_vertex
# #         @show "same vertex, considering insertion potential"
#         emission_likelihood = observed_error_rate
#         transition_likelihood = observed_error_rate
#         state_likelihood = kmer_likelihoods[current_vertex]
#         path_likelihood = prior_state_likelihood * emission_likelihood * transition_likelihood * state_likelihood
#         path = [last(arrival_paths[prior_vertex, prior_state])]

#         if current_state_likelihood > state_likelihoods[current_vertex, current_state]
# #             @show "selecting path"
# #             @show path
# #             @show path_likelihood
#             state_likelihoods[current_vertex, current_state] = path_likelihood
#             arrival_paths[current_vertex, current_state] = path
#         end
#     # one or more step path - match, mismatch, or deletion in observed sequence relative to kmer graph
#     elseif !isempty(shortest_path)
# #         @show "path is viable!"
# #         @show "considering shortest path: $(shortest_path)"

#         initial_path_state = last(arrival_paths[prior_vertex, prior_state])

#         path = Vector{typeof(initial_path_state)}(undef, length(shortest_path))
#         path[1] = initial_path_state

#         path_likelihood::Float64 = state_likelihoods[prior_vertex, prior_state]

#         for i in 2:length(shortest_path)

#             this_vertex = shortest_path[i]
#             prior_vertex, prior_orientation = path[i-1]
#             edge = Graphs.Edge(prior_vertex, this_vertex)

#             possible_edge_orientations::Set{NamedTuple{(:source_orientation, :destination_orientation), Tuple{Bool, Bool}}} = kmer_graph.eprops[edge][:orientations]
            
# #             @show possible_edge_orientations
            
#             if !ismissing(prior_orientation)
#                 possible_edge_orientations = filter(o -> o.source_orientation == prior_orientation, possible_edge_orientations)
#             end
            
# #             @show possible_edge_orientations
            
#             if isempty(possible_edge_orientations)
#                 path_likelihood *= 0.0
#                 path = Vector{eltype(path)}()
# #                 @show "no possible orientations, bailing early"
#                 break
#             end

# #             @show prior_orientation
#             if ismissing(prior_orientation)
#                 if transition_likelihoods[true][prior_vertex, this_vertex] > transition_likelihoods[false][prior_vertex, this_vertex]
#                     prior_orientation = true
#                     transition_likelihood = transition_likelihoods[true][prior_vertex, this_vertex]::Float64
#                 elseif transition_likelihoods[true][prior_vertex, this_vertex] < transition_likelihoods[false][prior_vertex, this_vertex]
#                     prior_orientation = false
#                     transition_likelihood = transition_likelihoods[false][prior_vertex, this_vertex]::Float64
#                 else transition_likelihoods[true][prior_vertex, this_vertex] == transition_likelihoods[false][prior_vertex, this_vertex]
#                     prior_orientation = missing
#                     transition_likelihood = transition_likelihoods[true][prior_vertex, this_vertex]::Float64
#                 end
#             else
#                 transition_likelihood = transition_likelihoods[prior_orientation][prior_vertex, this_vertex]::Float64
#             end
#             state_likelihood::Float64 = kmer_likelihoods[this_vertex]
#             path_likelihood *= transition_likelihood * state_likelihood
            
#             if length(possible_edge_orientations) == 1
#                 orientation = first(possible_edge_orientations).destination_orientation
#                 path[i] = this_vertex => orientation
#             else
#                 path[i] = this_vertex => missing
#             end
#         end

#         # see if new nucleotide is a match or mismatch to terminal kmer in path
#         if !isempty(path) && path_likelihood > 0
#             terminal_kmer_index, terminal_kmer_orientation = last(path)
#             terminal_kmer = BioSequences.LongDNASeq(kmer_graph.vprops[terminal_kmer_index][:kmer])::BioSequences.LongDNASeq
#             if ismissing(terminal_kmer_orientation)
#                 fw_is_match = observed_nucleotide == last(terminal_kmer)
#                 bw_is_match = observed_nucleotide == last(BioSequences.reverse_complement!(terminal_kmer))
#                 if fw_ismatch && !bw_is_match
#                     path[end] = terminal_kmer_index => true
#                     path_likelihood *= 1 - observed_error_rate
#                 elseif !fw_ismatch && bw_is_match
#                     path[end] = terminal_kmer_index => false
#                     path_likelihood *= 1 - observed_error_rate
#                 elseif fw_ismatch && bw_is_match
#                     path_likelihood *= 1 - observed_error_rate
#                 elseif !fw_ismatch && !bw_is_match
#                     path_likelihood *= observed_error_rate
#                 end
#             elseif terminal_kmer_orientation
#                 is_match = observed_nucleotide == last(terminal_kmer)
#                 if is_match
#                     path_likelihood *= 1 - observed_error_rate
#                 else
#                     path_likelihood *= observed_error_rate
#                 end
#             else
#                 terminal_kmer = BioSequences.reverse_complement!(terminal_kmer)
#                 is_match = observed_nucleotide == last(terminal_kmer)
#                 if is_match
#                     path_likelihood *= 1 - observed_error_rate
#                 else
#                     path_likelihood *= observed_error_rate
#                 end
#             end
#         end

#         if path_likelihood > state_likelihoods[current_vertex, current_state]
# #             @show "selecting path"
# #             @show path
# #             @show path_likelihood
#             state_likelihoods[current_vertex, current_state] = path_likelihood
#             arrival_paths[current_vertex, current_state] = path
#         end
#     end
#     return
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function polish_fastq(kmer_graph, fastq_file)

# #     @info "Assessing kmer likelihoods"
#     kmers = [kmer_graph.vprops[v][:kmer] for v in Graphs.vertices(kmer_graph)]
# #     kmer_counts = [length(kmer_graph.vprops[v][:evidence]) for v in Graphs.vertices(kmer_graph)]
#     kmer_counts = [kmer_graph.vprops[v][:weight] for v in Graphs.vertices(kmer_graph)]
#     kmer_likelihoods = kmer_counts ./ sum(kmer_counts)
#     k = kmer_graph.gprops[:k]
#     kmer_type = BioSequences.BigDNAMer{k}
#     total_kmers = length(kmers)
    
# #     @info "determining shortest paths between kmers"
#     shortest_paths = Graphs.enumerate_paths(Graphs.floyd_warshall_shortest_paths(kmer_graph));
    
#     @info "counting the number of records to establish runtime estimate"
#     number_of_records = 0
#     for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
#         number_of_records += 1
#     end
#     progress_bar = ProgressMeter.Progress(number_of_records, 1)
    
#     output_fastq_file = replace(fastq_file, ".fastq" => ".k$(kmer_graph.gprops[:k]).fastq")
#     fastq_writer = FASTX.FASTQ.Writer(open(output_fastq_file, "w"))
#     for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
#         ProgressMeter.next!(progress_bar)
        
# #         @info "Initializing matrices"
#         total_states = length(FASTX.sequence(fastq_record))-k+1
#         transition_likelihoods = initialize_transition_probabilities(kmer_graph)
#         state_likelihoods = zeros(total_kmers, total_states)
#         arrival_paths = fill(Pair{Int, Union{Bool, Missing}}[], total_kmers, total_states)

# #         @info "Determining Likelihoods of initial states"
#         initial_state = first(BioSequences.each(kmer_type, FASTX.sequence(fastq_record)))
#         current_state = 1
#         # note this is a place for potential improvement, use the q value at each base to guide probability rather than median
#         median_q_value = Statistics.median(Int.(FASTX.quality(fastq_record)[1:k]))
#         current_error_rate = q_value_to_error_rate(median_q_value)
#         # canonical_kmer = BioSequences.canonical(initial_state.fw)
#         set_initial_state_likelihoods!(
#                 kmer_graph,
#                 initial_state,
#                 kmer_likelihoods,
#                 current_error_rate,
#                 state_likelihoods,
#                 arrival_paths
#             )

# #         @info "Determining likelihood of downstream states"

# #         non_singleton_states = findall(kmer_counts .> 1)

#         for current_state in 2:total_states
#             prior_state = current_state - 1

#         #     observed_kmer = BioSequences.BigDNAMer{k}(FASTX.sequence(fastq_record)[current_state:current_state+k-1])

#         #     @assert observed_kmer == collect(BioSequences.each(kmer_type, FASTX.sequence(fastq_record)))[current_state].fw

#         #     canonical_kmer = BioSequences.canonical(observed_kmer)

#             observed_nucleotide = FASTX.sequence(fastq_record)[k-1+current_state]
#         #     observed_nucleotide = last(observed_kmer)
#             observed_quality_score = FASTX.quality(fastq_record)[k-1+current_state]
#             observed_error_rate = q_value_to_error_rate(observed_quality_score)

#             # we'll assess prior states in order of decreasing likelihood
#             # such that we maximize how frequently we are able to utilize the
#             # current_state_likelihood > candidate prior state
#             # break that won't bother evaluating lower likelihood possibilities
#             prior_states_in_decreasing_likelihood = sortperm(state_likelihoods[:, prior_state], rev=true)

#             # and skip all prior states with zero probability

#             for current_vertex in total_states
#                 for prior_vertex in prior_states_in_decreasing_likelihood
#                     if state_likelihoods[prior_vertex, prior_state] > 0
#                         run_viterbi!(
#                                 current_state,
#                                 prior_state,
#                                 observed_nucleotide,
#                                 observed_quality_score,
#                                 observed_error_rate,
#                                 current_vertex,
#                                 prior_vertex,
#                                 state_likelihoods,
#                                 transition_likelihoods,
#                                 shortest_paths,
#                                 arrival_paths,
#                                 kmer_graph,
#                                 kmer_likelihoods
#                                 )
#                     end
#                 end
#             end
#         end

# #         try
#         maximum_likelihood_path, maximum_likelihood_value = 
#             determine_maximum_likelihood_path(
#                 state_likelihoods,
#                 arrival_paths
#                 )
# #         catch
# #             return state_likelihoods, arrival_paths
# #         end

#         sequence = oriented_path_to_sequence(kmer_graph, maximum_likelihood_path)

# #         @info "comparing to original path"
#         original_sequence_likelihood = oriented_path_to_likelihood(kmer_graph, kmers, kmer_likelihoods, transition_likelihoods, fastq_record)
#         relative_likelihood = maximum_likelihood_value / original_sequence_likelihood
# #         relative_likelihood_formatted = NumericIO.formatted(relative_likelihood, ndigits=1, charset=:ASCII)
# #         println("relative likelihood of new path to old path is $(relative_likelihood_formatted)")

# #         @info "writing updated record"
#         identifier = FASTX.identifier(fastq_record) * "_k$(k)"
#         description = string(relative_likelihood)
#         # because the sequences won't always be the same length, we take an ordered sampling with replacement
#         # which introduces some random error but preserves overall patterns and areas of high/low accuracy
#         quality_scores = StatsBase.sample(FASTX.quality(fastq_record), length(sequence), ordered=true)

#         new_fastq_record = FASTX.FASTQ.Record(
#             identifier,
#             description,
#             sequence,
#             quality_scores
#         )
#         write(fastq_writer, new_fastq_record)
#     end
#     close(fastq_writer)
#     return output_fastq_file
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_maximum_likelihood_path(
#     state_likelihoods,
#     arrival_paths
#     )
#     maximum_likelihood_value = maximum(state_likelihoods[:, end])

#     maximum_likelihood_path_indices = findall(state_likelihoods[:, end] .== maximum_likelihood_value)

#     # if multiple paths are tied, randomly choose one
#     maximum_likelihood_path_index = rand(maximum_likelihood_path_indices)

#     maximum_likelihood_path = arrival_paths[maximum_likelihood_path_index, end]

#     for state_index in size(state_likelihoods, 2)-1:-1:1
#         next_kmer, next_orientation = first(maximum_likelihood_path)
#         maximum_likelihood_arrival_path = arrival_paths[next_kmer, state_index]
        
#         is_match = last(maximum_likelihood_arrival_path) == (next_kmer => next_orientation)
#         if !ismissing(is_match) && !is_match
#             error("breaking")
#         end
#         maximum_likelihood_path = vcat(maximum_likelihood_arrival_path[1:end-1], maximum_likelihood_path)
#     end
#     return maximum_likelihood_path, maximum_likelihood_value
# end


"""
document me
"""
function sequence_to_stranded_path(stranded_kmers, sequence)
    KMER_TYPE = typeof(first(stranded_kmers))
    path = Vector{Pair{Int, Bool}}()
    for (i, kmer) in Kmers.EveryKmer{KMER_TYPE}(BioSequences.LongDNA{4}(sequence))
        kmer_index = findfirst(stranded_kmer -> kmer == stranded_kmer, stranded_kmers)
        orientation = true
        push!(path, kmer_index => orientation)
    end
    return path
end

"""
document me
"""
function path_to_sequence(kmers, path)
    # @show path
    sequence = BioSequences.LongDNA{4}(kmers[first(first(path))])
    for i in 2:length(path)
        push!(sequence, kmers[first(path[i])][end])
    end
    return sequence
end

# """
#     edge_probability(stranded_kmer_graph, edge)
#
# Compute the probability of an edge relative to all edges of a given source vertex
# """
function edge_probability(stranded_kmer_graph, edge)
    neighbors = Graphs.outneighbors(stranded_kmer_graph, edge.src)
    neighbor_under_consideration = findfirst(neighbor -> neighbor == edge.dst, neighbors)
    edge_weights = [length(stranded_kmer_graph.eprops[Graphs.Edge(edge.src, neighbor)][:coverage]) for neighbor in neighbors]
    if sum(edge_weights) == 0
        p = 0.0
    else
        edge_probabilities = edge_weights ./ sum(edge_weights)
        p = edge_probabilities[neighbor_under_consideration]
    end
    return p
end

# """
#     viterbi_maximum_likelihood_traversal(stranded_kmer_graph; error_rate=1-1/k + 1))
#
#
#
# ```jldoctest
# import Primes
# julia> for prime in Primes.primes(3, 101)
#            println("$prime assumes accuracy of $(1-1/(prime + 1))")
#        end
# 3 assumes accuracy of 0.6666666666666667
# 5 assumes accuracy of 0.8
# 7 assumes accuracy of 0.8571428571428572
# 11 assumes accuracy of 0.9090909090909091
# 13 assumes accuracy of 0.9230769230769231
# 17 assumes accuracy of 0.9411764705882353
# 19 assumes accuracy of 0.9473684210526316
# 23 assumes accuracy of 0.9565217391304348
# 29 assumes accuracy of 0.9655172413793104
# 31 assumes accuracy of 0.967741935483871
# 37 assumes accuracy of 0.972972972972973
# 41 assumes accuracy of 0.975609756097561
# 43 assumes accuracy of 0.9767441860465116
# 47 assumes accuracy of 0.9787234042553191
# 53 assumes accuracy of 0.9811320754716981
# 59 assumes accuracy of 0.9830508474576272
# 61 assumes accuracy of 0.9836065573770492
# 67 assumes accuracy of 0.9850746268656716
# 71 assumes accuracy of 0.9859154929577465
# 73 assumes accuracy of 0.9863013698630136
# 79 assumes accuracy of 0.9873417721518988
# 83 assumes accuracy of 0.9879518072289156
# 89 assumes accuracy of 0.9887640449438202
# 97 assumes accuracy of 0.9896907216494846
# 101 assumes accuracy of 0.9900990099009901
# ```
# """
function viterbi_maximum_likelihood_traversals(stranded_kmer_graph;
                                               error_rate::Float64=1/(stranded_kmer_graph.gprops[:k] + 1),
                                               verbosity::String="dataset")
    @assert verbosity in ["debug", "reads", "dataset"]
    if error_rate >= .5
        error("Error rate >= 50%. Did you enter the accuracy by mistake?")
    end

    if verbosity in ["debug", "reads", "dataset"]
        println("computing kmer counts...")
    end
    stranded_kmer_counts = [length(stranded_kmer_graph.vprops[vertex][:coverage]) for vertex in Graphs.vertices(stranded_kmer_graph)]
    if verbosity in ["debug", "reads", "dataset"]
        println("computing kmer state likelihoods...")
    end
    stranded_kmer_likelihoods = stranded_kmer_counts ./ sum(stranded_kmer_counts)
    accuracy = 1 - error_rate

    if verbosity in ["debug"]
        println("STATE LIKELIHOODS:")
        println("\tkmer\tcount\tlikelihood")
        for vertex in Graphs.vertices(stranded_kmer_graph)
            kmer = stranded_kmer_graph.gprops[:stranded_kmers][vertex]
            count = stranded_kmer_counts[vertex]
            likelihood = stranded_kmer_likelihoods[vertex]
            println("\t$kmer\t$count\t$likelihood")
        end
    end
    if verbosity in ["debug", "reads", "dataset"]
        println("finding shortest paths between kmers...")
    end
    shortest_paths = Graphs.enumerate_paths(Graphs.floyd_warshall_shortest_paths(stranded_kmer_graph))
    K = stranded_kmer_graph.gprops[:K]
    for K1 in 1:K
        for K2 in 1:K
            if K1 != K2
                shortest_path = shortest_paths[K1][K2]
                path_likelihood = 1.0
                for ui in 1:length(shortest_path)-1
                    u = shortest_path[ui]
                    v = shortest_path[ui + 1]
                    # likelihood of the transition
                    path_likelihood *= edge_probability(stranded_kmer_graph, Graphs.Edge(u, v))
                end
                if path_likelihood == 0.0
                    shortest_paths[K1][K2] = Vector{Int}()
                end
            elseif K1 == K2
                # the shortest path from a kmer to itself is an insertion (no edge)
                # so need to manually check for self loops
                if Graphs.has_edge(stranded_kmer_graph, Graphs.Edge(K1, K2))
                    if edge_probability(stranded_kmer_graph, Graphs.Edge(K1, K2)) != 0.0
                        shortest_paths[K1][K2] = [K1, K2]
                    else
                        shortest_paths[K1][K2] = Vector{Int}()
                    end
                # otherwise, check to see if any outneighbors connect back to the kmer
                else
                    connected_outneighbors = filter(outneighbor -> Graphs.has_path(stranded_kmer_graph, outneighbor, K2), Graphs.outneighbors(stranded_kmer_graph, K1))
                    if !isempty(connected_outneighbors)
                        outneighbor_cycles = [[K1, shortest_paths[outneighbor][K2]...] for outneighbor in connected_outneighbors]
                        cycle_likelihoods = ones(length(outneighbor_cycles))
                        for (i, cycle) in enumerate(outneighbor_cycles)
                            for ui in 1:length(cycle)-1
                                u = cycle[ui]
                                v = cycle[ui + 1]
                                # likelihood of the transition
                                cycle_likelihoods[i] *= edge_probability(stranded_kmer_graph, Graphs.Edge(u, v))
                            end
                            # include likelihoods of states
                            for vertex in cycle[2:end-1]
                                cycle_likelihoods[i] *= stranded_kmer_likelihoods[vertex]
                            end
                        end
                        path_likelihood = maximum(cycle_likelihoods)
                        max_likelihood_cycle_indices = findall(cycle_likelihoods .== path_likelihood)
                        shortest_paths[K1][K2] = outneighbor_cycles[first(max_likelihood_cycle_indices)]
                    else
                        shortest_paths[K1][K2] = Vector{Int}()
                    end
                end
            end
            if length(shortest_paths[K1][K2]) == 1
                shortest_paths[K1][K2] = Vector{Int}()
            end
        end
    end

    if verbosity in ["debug"]
        for K1 in 1:K
            for K2 in 1:K
                println("\t$K1\t$K2\t$(shortest_paths[K1][K2])")
            end
        end
    end

    total_bases_observed = 0
    total_edits_accepted = 0

    corrected_observations = FASTX.FASTA.Record[]
    if verbosity in ["debug", "reads", "dataset"]
        println("finding viterbi maximum likelihood paths for observed sequences...")
    end
    # p = Progress(length(stranded_kmer_graph.gprops[:observed_paths]))
    for (observation_index, observed_path) in enumerate(stranded_kmer_graph.gprops[:observed_paths])
        if verbosity in ["debug", "reads"]
            println("\nevaluating sequence $observation_index of $(length(stranded_kmer_graph.gprops[:observed_paths]))")
        end
        # consider switching to log transform
        kmer_likelihoods = zeros(Graphs.nv(stranded_kmer_graph), length(observed_path))
        kmer_arrival_paths = Array{Vector{Int}}(undef, Graphs.nv(stranded_kmer_graph), length(observed_path))
        edit_distances = zeros(Int, Graphs.nv(stranded_kmer_graph), length(observed_path))
        # changed here!!
        observed_kmer_index, observed_kmer_orientation = observed_path[1]
        
        observed_kmer_sequence = stranded_kmer_graph.gprops[:stranded_kmers][observed_kmer_index]
        for hidden_kmer_index in Graphs.vertices(stranded_kmer_graph)
            hidden_kmer_sequence = stranded_kmer_graph.gprops[:stranded_kmers][hidden_kmer_index]
            alignment_result = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), observed_kmer_sequence, hidden_kmer_sequence)
            number_of_matches = BioAlignments.count_matches(BioAlignments.alignment(alignment_result))
            number_of_edits = stranded_kmer_graph.gprops[:k] - number_of_matches
            kmer_likelihoods[hidden_kmer_index, 1] = stranded_kmer_likelihoods[hidden_kmer_index]
            for match in 1:number_of_matches
                kmer_likelihoods[hidden_kmer_index, 1] *= accuracy
            end
            for edit in 1:number_of_edits
                kmer_likelihoods[hidden_kmer_index, 1] *= error_rate
            end
            kmer_arrival_paths[hidden_kmer_index, 1] = Vector{Int}()
            edit_distances[hidden_kmer_index, 1] = number_of_edits
        end
        kmer_likelihoods[:, 1] ./= sum(kmer_likelihoods[:, 1])
        # from here on, all probabilities are log transformed
        kmer_likelihoods[:, 1] .= log.(kmer_likelihoods[:, 1])
        if verbosity in ["debug"]
            println("\tconsidering path state 1")
            println("\t\tobserved kmer $observed_kmer_sequence")
            println("\t\tInitial state log likelihoods:")
            for line in split(repr(MIME("text/plain"), kmer_likelihoods[:, 1]), '\n')
                println("\t\t\t$line")
            end
        end
        for observed_path_index in 2:length(observed_path)
            # changed!!
            observed_kmer_index, observed_kmer_orientation = observed_path[observed_path_index]
            observed_base = stranded_kmer_graph.gprops[:stranded_kmers][observed_kmer_index][end]

            if verbosity in ["debug"]
                println("\tconsidering path state $observed_path_index")
                println("\t\tobserved base $observed_base")
            end

            MATCH = 1
            MISMATCH = 2
            DELETION = 3
            INSERTION = 4
            arrival_likelihoods = ones(K, 4)
            arrival_paths = fill(Vector{Int}(), K, 4)

            for K2 in 1:K
                kmer_base = stranded_kmer_graph.gprops[:stranded_kmers][K2][end]
                base_is_match = kmer_base == observed_base

                maximum_likelihood = log(0.0)
                maximum_likelihood_path = Vector{Int}()
                maximum_likelihood_edit_distance = 0

                for K1 in 1:K
                    shortest_path = shortest_paths[K1][K2]
                    if length(shortest_path) >= 2
                        edit_distance = Int(!base_is_match) + length(shortest_path) - 2
                        if edit_distance == 0
                            p = kmer_likelihoods[K1, observed_path_index-1] +
                                log(accuracy) + log(stranded_kmer_likelihoods[K2])
                        else
                            p = kmer_likelihoods[K1, observed_path_index-1] +
                                log(error_rate^edit_distance) + log(stranded_kmer_likelihoods[K2])
                        end
                        edit_distance += edit_distances[K1, observed_path_index-1]
                    else
                        p = log(0.0)
                    end
                    if K1 == K2 # consider insertion
                        # in theory, I don't think we should care if the base
                        # matches or not because it's an inserted & erroneous
                        # base, but in practice it's necessary to balance
                        # insertion probabilities with deletion probabilities
                        insertion_p = kmer_likelihoods[K1, observed_path_index-1] +
                                      log(error_rate^(1 + Int(!base_is_match))) + log(stranded_kmer_likelihoods[K2])
                        if insertion_p > p
                            p = insertion_p
                            edit_distance = edit_distances[K1, observed_path_index-1] + 1
                            shortest_path = [K2]
                        end
                    end
                    if p > maximum_likelihood
                        maximum_likelihood = p
                        maximum_likelihood_path = shortest_path
                        maximum_likelihood_edit_distance = edit_distance
                    end
                end
                kmer_likelihoods[K2, observed_path_index] = maximum_likelihood
                kmer_arrival_paths[K2, observed_path_index] = maximum_likelihood_path
                edit_distances[K2, observed_path_index] = maximum_likelihood_edit_distance
            end

            if verbosity in ["debug"]
                println("\t\tkmer log likelihoods")
                for line in split(repr(MIME("text/plain"), kmer_likelihoods), '\n')
                    println("\t\t\t$line")
                end
                println("\t\tarrival paths")
                for line in split(repr(MIME("text/plain"), kmer_arrival_paths), '\n')
                    println("\t\t\t$line")
                end
            end
        end

        if verbosity in ["debug"]
            println("\n\tInputs for viterbi maximum likelihood traversal evaluation:")
            println("\t\tkmer log likelihoods")
            for line in split(repr(MIME("text/plain"), kmer_likelihoods), '\n')
                println("\t\t\t$line")
            end
            println("\t\tkmer arrival paths")
            for line in split(repr(MIME("text/plain"), kmer_arrival_paths), '\n')
                println("\t\t\t$line")
            end
            println("\t\tedit distances")
            for line in split(repr(MIME("text/plain"), edit_distances), '\n')
                println("\t\t\t$line")
            end
        end

        ## backtrack
        maximum_likelihood_path_value = maximum(kmer_likelihoods[:, end])
        maximum_likelihood_path_indices = findall(kmer_likelihoods[:, end] .== maximum_likelihood_path_value)
        # if multiple paths are tied, randomly choose one
        maximum_likelihood_path_index = rand(maximum_likelihood_path_indices)
        maximum_likelihood_edit_distance = edit_distances[maximum_likelihood_path_index, end]

        if length(kmer_arrival_paths[maximum_likelihood_path_index, end]) > 0
            maximum_likelihood_path = last(kmer_arrival_paths[maximum_likelihood_path_index, end])
            for observed_path_index in length(observed_path):-1:1
                maximum_likelihood_arrival_path = kmer_arrival_paths[maximum_likelihood_path_index, observed_path_index]
                maximum_likelihood_path = vcat(maximum_likelihood_arrival_path[1:end-1], maximum_likelihood_path)
                maximum_likelihood_path_index = first(maximum_likelihood_path)
            end
        else
            maximum_likelihood_path = [maximum_likelihood_path_index]
        end
        observed_sequence = path_to_sequence(stranded_kmer_graph.gprops[:stranded_kmers], observed_path)
        maximum_likelihood_sequence = path_to_sequence(stranded_kmer_graph.gprops[:stranded_kmers], maximum_likelihood_path)
        if verbosity in ["debug", "reads"]
            println("\tobserved sequence                 $observed_sequence")
            println("\tmaximum likelihood sequence       $maximum_likelihood_sequence")
            println("\tmaximum likelihood edit distance  $maximum_likelihood_edit_distance")
        end
        total_bases_observed += length(observed_sequence)
        total_edits_accepted += maximum_likelihood_edit_distance
        id = stranded_kmer_graph.gprops[:observation_ids][observation_index]
        kmer_stamped_id = id * "_" * string(stranded_kmer_graph.gprops[:k])
        push!(corrected_observations, FASTX.FASTA.Record(kmer_stamped_id, maximum_likelihood_sequence))
        # progress meter
        # next!(p)
    end
    if verbosity in ["debug", "reads", "dataset"]
        println("\nDATASET STATISTICS:")
        println("\tassumed error rate    $(error_rate * 100)%")
        println("\ttotal bases observed  $total_bases_observed")
        println("\ttotal edits accepted  $total_edits_accepted")
        println("\tinferred error rate   $((total_edits_accepted/total_bases_observed) * 100)%")
    end
    return corrected_observations
end

"""
    build_stranded_kmer_graph(canonical_kmers, observations)

Create a weighted, strand-specific kmer (de bruijn) graph from a set of kmers
and a series of sequence observations in FASTA format.
"""
function build_stranded_kmer_graph(canonical_kmers, observations)
    if isempty(canonical_kmers)
        @error "isempty(canonical_kmers) = $(isempty(canonical_kmers))"
    elseif isempty(observations)
        @error "isempty(observations) = $(isempty(observations))"
    end
    stranded_kmers = sort!(vcat(canonical_kmers, [BioSequences.reverse_complement(kmer) for kmer in canonical_kmers]))
    stranded_kmer_to_reverse_complement_map = [
        findfirst(stranded_kmer -> BioSequences.reverse_complement(stranded_kmer) == kmer, stranded_kmers) for kmer in stranded_kmers
    ]
    stranded_kmer_graph = MetaGraphs.MetaDiGraph(length(stranded_kmers))
    stranded_kmer_graph.gprops[:stranded_kmers] = stranded_kmers
    stranded_kmer_graph.gprops[:reverse_complement_map] = stranded_kmer_to_reverse_complement_map
    stranded_kmer_graph.gprops[:k] = length(first(stranded_kmers))
    stranded_kmer_graph.gprops[:K] = length(stranded_kmers)
    stranded_kmer_graph.gprops[:observation_color_map] = Vector{Int}()
    stranded_kmer_graph.gprops[:observation_ids] = Vector{String}()
    stranded_kmer_graph.gprops[:observed_paths] = Vector{Vector{Pair{Int, Bool}}}()
    for vertex in 1:Graphs.nv(stranded_kmer_graph)
        stranded_kmer_graph.vprops[vertex] = Dict(:coverage => Vector{Pair{Int, Pair{Int, Bool}}}())
    end
    for (observation_index, observation) in enumerate(observations)
        observation_id = FASTX.FASTA.identifier(observation)
        observed_sequence = FASTX.FASTA.sequence(observation)
        if length(observed_sequence) < stranded_kmer_graph.gprops[:k]
            @error "skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        else
            observed_path = sequence_to_stranded_path(stranded_kmer_graph.gprops[:stranded_kmers], observed_sequence)
            i = 1
            ui, ui_orientation = observed_path[i]
            ui_coverage = (observation_index => (i => ui_orientation ))
            push!(stranded_kmer_graph.vprops[ui][:coverage], ui_coverage)
            for i in 2:length(observed_path)
                vi, vi_orientation = observed_path[i]
                vi_coverage = (observation_index => (i => vi_orientation))
                push!(stranded_kmer_graph.vprops[vi][:coverage], vi_coverage)
                edge_coverage = ui_coverage => vi_coverage
                if Graphs.has_edge(stranded_kmer_graph, ui, vi)
                    push!(stranded_kmer_graph.eprops[Graphs.Edge(ui, vi)][:coverage], edge_coverage)
                else
                    Graphs.add_edge!(stranded_kmer_graph, ui, vi, Dict(:coverage => [edge_coverage]))
                end
                # not sure this is necessary
#                 ui′ = stranded_kmer_graph.gprops[:reverse_complement_map][ui]
#                 vi′ = stranded_kmer_graph.gprops[:reverse_complement_map][vi]
#                 if !LightGraphs.has_edge(stranded_kmer_graph, vi′, ui′)
#                     LightGraphs.add_edge!(stranded_kmer_graph, vi′, ui′, Dict(:coverage => Vector{typeof(edge_coverage)}()))
#                 end
                ui, ui_orientation = vi, vi_orientation
                ui_coverage = vi_coverage
            end
            push!(stranded_kmer_graph.gprops[:observed_paths], observed_path)
            push!(stranded_kmer_graph.gprops[:observation_ids], observation_id)
            push!(stranded_kmer_graph.gprops[:observation_color_map], observation_index)
        end
    end
    return stranded_kmer_graph
end