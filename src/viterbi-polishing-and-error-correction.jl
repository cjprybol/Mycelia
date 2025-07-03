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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Finds maximum likelihood paths through a stranded k-mer graph using the Viterbi algorithm
to correct sequencing errors.

# Arguments
- `stranded_kmer_graph`: A directed graph where vertices represent k-mers and edges represent overlaps
- `error_rate::Float64`: Expected per-base error rate (default: 1/(k+1)). Must be < 0.5
- `verbosity::String`: Output detail level ("debug", "reads", or "dataset")

# Returns
Vector of FASTX.FASTA.Record containing error-corrected sequences

# Details
- Uses dynamic programming to find most likely path through k-mer graph
- Accounts for matches, mismatches, insertions and deletions
- State likelihoods based on k-mer coverage counts
- Transition probabilities derived from error rate
- Progress tracking based on verbosity level

# Notes
- Error rate should be probability of error (e.g. 0.01 for 1%), not accuracy
- Higher verbosity levels ("debug", "reads") provide detailed path finding information
- "dataset" verbosity shows only summary statistics
"""
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
$(DocStringExtensions.TYPEDSIGNATURES)

Process and error-correct a FASTQ sequence record using a kmer graph and path resampling.

# Arguments
- `record`: FASTQ record containing the sequence to process
- `kmer_graph`: MetaGraph containing the kmer network and associated properties
- `yen_k_shortest_paths_and_weights`: Cache of pre-computed k-shortest paths between nodes
- `yen_k`: Number of alternative paths to consider during resampling (default: 3)

# Description
Performs error correction by:
1. Trimming low-quality sequence ends
2. Identifying stretches requiring resampling between solid branching kmers
3. Selecting alternative paths through the kmer graph based on:
   - Path quality scores
   - Transition likelihoods
   - Path length similarity to original sequence

# Returns
- Modified FASTQ record with error-corrected sequence and updated quality scores
- Original record if no error correction was needed

# Required Graph Properties
The kmer_graph must contain the following properties:
- :ordered_kmers
- :likely_valid_kmer_indices  
- :kmer_indices
- :branching_nodes
- :assembly_k
- :transition_likelihoods
- :kmer_mean_quality
- :kmer_total_quality
"""
function process_fastq_record(;record, kmer_graph, yen_k_shortest_paths_and_weights, yen_k=3)
    ordered_kmers = MetaGraphs.get_prop(kmer_graph, :ordered_kmers)
    likely_valid_kmers = Set(ordered_kmers[MetaGraphs.get_prop(kmer_graph, :likely_valid_kmer_indices)])
    kmer_to_index_map = MetaGraphs.get_prop(kmer_graph, :kmer_indices)
    branching_nodes_set = MetaGraphs.get_prop(kmer_graph, :branching_nodes)
    assembly_k = MetaGraphs.get_prop(kmer_graph, :assembly_k)
    transition_likelihoods = MetaGraphs.get_prop(kmer_graph, :transition_likelihoods)
    kmer_mean_quality = MetaGraphs.get_prop(kmer_graph, :kmer_mean_quality)
    kmer_total_quality = MetaGraphs.get_prop(kmer_graph, :kmer_total_quality)
    
    new_record_identifier = FASTX.identifier(record) * ".k$(assembly_k)"
    record_sequence = BioSequences.LongDNA{4}(FASTX.sequence(record))

    kmer_type = Kmers.DNAKmer{assembly_k}
    record_kmers = last.(collect(Kmers.EveryKmer{kmer_type}(record_sequence)))
    record_quality_scores = collect(FASTX.quality_scores(record))
    record_kmer_quality_scores = [record_quality_scores[i:i+assembly_k-1] for i in 1:length(record_quality_scores)-assembly_k+1]
    
    record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
    record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
    record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    
    # trim beginning of fastq
    initial_solid_kmer = findfirst(record_kmer_solidity)
    if isnothing(initial_solid_kmer)
        return record
    elseif initial_solid_kmer > 1
        record_kmers = record_kmers[initial_solid_kmer:end]
        record_kmer_quality_scores = record_kmer_quality_scores[initial_solid_kmer:end]
        record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
        record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
        record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    end
    initial_solid_kmer = 1
    
    # trim end of fastq
    last_solid_kmer = findlast(record_kmer_solidity)
    if last_solid_kmer != length(record_kmer_solidity)
        record_kmers = record_kmers[1:last_solid_kmer]
        record_kmer_quality_scores = record_kmer_quality_scores[1:last_solid_kmer]
        record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
        record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
        record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    end
    
    # identify low quality runs and the solid branchpoints we will use for resampling
    solid_branching_kmer_indices = findall(record_solid_branching_kmers)
    resampling_stretches = find_resampling_stretches(;record_kmer_solidity, solid_branching_kmer_indices)

    # nothing to do
    if isempty(resampling_stretches)
        return record
    end
    trusted_range = 1:max(first(first(resampling_stretches))-1, 1)
    
    new_record_kmers = record_kmers[trusted_range]
    new_record_kmer_qualities = record_kmer_quality_scores[trusted_range]
    
    
    for (i, resampling_stretch) in enumerate(resampling_stretches)
        starting_solid_kmer = record_kmers[first(resampling_stretch)]
        ending_solid_kmer = record_kmers[last(resampling_stretch)]
        
        current_quality_scores = record_quality_scores[resampling_stretch]
        u = kmer_to_index_map[starting_solid_kmer]
        v = kmer_to_index_map[ending_solid_kmer]
        if !haskey(yen_k_shortest_paths_and_weights, u => v)
            yen_k_result = Graphs.yen_k_shortest_paths(kmer_graph, u, v, Graphs.weights(kmer_graph), yen_k)
            yen_k_shortest_paths_and_weights[u => v] = Vector{Pair{Vector{Int}, Float64}}()
            for path in yen_k_result.paths
                path_weight = Statistics.mean([kmer_total_quality[ordered_kmers[node]] for node in path])
                path_transition_likelihoods = 1.0
                for (a, b) in zip(path[1:end-1], path[2:end])
                    path_transition_likelihoods *= transition_likelihoods[a, b]
                end
                joint_weight = path_weight * path_transition_likelihoods
                push!(yen_k_shortest_paths_and_weights[u => v], path => joint_weight)
            end
        end
        yen_k_path_weights = yen_k_shortest_paths_and_weights[u => v]      
        if length(yen_k_path_weights) > 1
            current_distance = length(resampling_stretch)
            initial_weights = last.(yen_k_path_weights)
            path_lengths = length.(first.(yen_k_path_weights))
            deltas = map(l -> abs(l-current_distance), path_lengths)
            adjusted_weights = initial_weights .* map(d -> exp(-d * log(2)), deltas)
            # make it more severe?
            # adjusted_weights = adjusted_weights.^2
            
            # and a bonus for usually being correct
            
            selected_path_index = StatsBase.sample(StatsBase.weights(adjusted_weights))
            selected_path, selected_path_weights = yen_k_path_weights[selected_path_index]
            selected_path_kmers = [ordered_kmers[kmer_index] for kmer_index in selected_path]
            
            if last(new_record_kmers) == first(selected_path_kmers)
                selected_path_kmers = selected_path_kmers[2:end]
            end
            append!(new_record_kmers, selected_path_kmers)
            selected_kmer_qualities = [Int8.(min.(typemax(Int8), floor.(kmer_mean_quality[kmer]))) for kmer in selected_path_kmers]
            append!(new_record_kmer_qualities, selected_kmer_qualities)
        else
            selected_path_kmers = record_kmers[resampling_stretch]
            if last(new_record_kmers) == first(selected_path_kmers)
                selected_path_kmers = selected_path_kmers[2:end]
            end
            append!(new_record_kmers, selected_path_kmers)
            selected_kmer_qualities = [Int8.(min.(typemax(Int8), floor.(kmer_mean_quality[kmer]))) for kmer in selected_path_kmers]
            append!(new_record_kmer_qualities, selected_kmer_qualities)
        end
        if i < length(resampling_stretches) # append high quality gap
            next_solid_start = last(resampling_stretch)+1
            next_resampling_stretch = resampling_stretches[i+1]
            next_solid_stop = first(next_resampling_stretch)-1
            if !isempty(next_solid_start:next_solid_stop)
                selected_path_kmers = record_kmers[next_solid_start:next_solid_stop]
                append!(new_record_kmers, selected_path_kmers)
                selected_kmer_qualities = record_kmer_quality_scores[next_solid_start:next_solid_stop]
                append!(new_record_kmer_qualities, selected_kmer_qualities)
            end
        else # append remainder of sequence
            @assert i == length(resampling_stretches)
            next_solid_start = last(resampling_stretch)+1
            if next_solid_start < length(record_kmers)
                selected_path_kmers = record_kmers[next_solid_start:end]
                append!(new_record_kmers, selected_path_kmers)
                selected_kmer_qualities = record_kmer_quality_scores[next_solid_start:end]
                append!(new_record_kmer_qualities, selected_kmer_qualities)
            end
        end
    end
    
    for (a, b) in zip(new_record_kmers[1:end-1], new_record_kmers[2:end])
        @assert a != b
    end
    new_record_sequence = Mycelia.kmer_path_to_sequence(new_record_kmers)
    new_record_quality_scores = new_record_kmer_qualities[1]
    for new_record_kmer_quality in new_record_kmer_qualities[2:end]
        push!(new_record_quality_scores, last(new_record_kmer_quality))
    end
    new_record = fastq_record(identifier=new_record_identifier, sequence=new_record_sequence, quality_scores=new_record_quality_scores)
    return new_record
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Polish FASTQ reads using a k-mer graph-based approach to correct potential sequencing errors.

# Arguments
- `fastq::String`: Path to input FASTQ file
- `k::Int=1`: Initial k-mer size parameter. Final assembly k-mer size may differ.

# Process
1. Builds a directed k-mer graph from input reads
2. Processes each read through the graph to find optimal paths
3. Writes corrected reads to a new FASTQ file
4. Automatically compresses output with gzip

# Returns
Named tuple with:
- `fastq::String`: Path to output gzipped FASTQ file
- `k::Int`: Final assembly k-mer size used
"""
function polish_fastq(;fastq, k=1)
    kmer_graph = build_directed_kmer_graph(fastq=fastq, k=k)
    assembly_k = MetaGraphs.get_prop(kmer_graph, :assembly_k)
    @info "polishing with k = $(assembly_k)"
    revised_records = []
    yen_k_shortest_paths_and_weights = Dict{Pair{Int, Int}, Vector{Pair{Vector{Int}, Float64}}}()
    ProgressMeter.@showprogress for record in collect(Mycelia.open_fastx(fastq))
        revised_record = process_fastq_record(;record, kmer_graph, yen_k_shortest_paths_and_weights)
        push!(revised_records, revised_record)
    end
    
    fastq_out = replace(fastq, Mycelia.FASTQ_REGEX => ".k$(assembly_k).fq")
    open(fastq_out, "w") do io
        fastx_io = FASTX.FASTQ.Writer(io)
        for record in revised_records
            write(fastx_io, record)
        end
        close(fastx_io)
    end
    run(`gzip --force $(fastq_out)`)
    return (fastq = fastq_out * ".gz", k=assembly_k)
end

# selected after trialing previous and next ks and finding those to be too unstable
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Performs iterative error correction on FASTQ sequences using progressively larger k-mer sizes.

Starting with the default k-mer size, this function repeatedly applies polishing steps,
incrementing the k-mer size until either reaching max_k or encountering instability.

# Arguments
- `fastq`: Path to input FASTQ file or FastqRecord object
- `max_k`: Maximum k-mer size to attempt (default: 89)
- `plot`: Whether to generate diagnostic plots (default: false)

# Returns
Vector of polishing results, where each element contains:
- k: k-mer size used
- fastq: resulting polished sequences
"""
function iterative_polishing(fastq, max_k = 89, plot=false)
    # initial polishing
    polishing_results = [polish_fastq(fastq=fastq)]
    while (!ismissing(last(polishing_results).k)) && (last(polishing_results).k < max_k)
        next_k = first(filter(k -> k > last(polishing_results).k, Mycelia.ks()))
        # @show next_k
        push!(polishing_results, polish_fastq(fastq=last(polishing_results).fastq, k=next_k))
    end
    return polishing_results
end