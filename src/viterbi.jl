"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function initialize_transition_probabilities(kmer_graph)
    
    total_kmers = Graphs.nv(kmer_graph)
    transition_likelihoods = Dict(
        true => SparseArrays.spzeros(total_kmers, total_kmers),
        false => SparseArrays.spzeros(total_kmers, total_kmers)
    )

    for edge in collect(Graphs.edges(kmer_graph))
#         weight = length(kmer_graph.eprops[edge][:evidence])
        weight = kmer_graph.eprops[edge][:weight]
        for o in kmer_graph.eprops[edge][:orientations]
            transition_likelihoods[o.source_orientation][edge.src, edge.dst] = weight
        end
    end

    for source_orientation in (true, false)
        for src in 1:total_kmers
            transition_weights = transition_likelihoods[source_orientation][src, :]
            total_weight = sum(transition_weights)
            dsts, vals = SparseArrays.findnz(transition_weights)
            for (dst, val) in zip(dsts, vals) 
                transition_likelihoods[source_orientation][src, dst] = val / total_weight
            end
            normalized_probability = sum(transition_likelihoods[source_orientation][src, :])
            @assert isapprox(normalized_probability, 0) || isapprox(normalized_probability, 1)
        end
    end
    return transition_likelihoods
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function set_initial_state_likelihoods!(
        kmer_graph,
        initial_state,
        kmer_likelihoods,
        error_rate,
        state_likelihoods,
        arrival_paths
    )
    for vertex in collect(Graphs.vertices(kmer_graph))
        hidden_kmer = kmer_graph.vprops[vertex][:kmer]

        fw_alignment = 
            BioAlignments.pairalign(
                BioAlignments.LevenshteinDistance(), 
                initial_state.fw, 
                hidden_kmer)

        fw_probability = kmer_likelihoods[vertex]

        for match in 1:BioAlignments.count_matches(BioAlignments.alignment(fw_alignment))
            fw_probability *= 1 - error_rate
        end

        for edit in 1:fw_alignment.value
            fw_probability *= error_rate
        end

        bw_alignment = 
            BioAlignments.pairalign(
                BioAlignments.LevenshteinDistance(),
                initial_state.bw,
                hidden_kmer)

        bw_probability = kmer_likelihoods[vertex]

        for match in 1:BioAlignments.count_matches(BioAlignments.alignment(bw_alignment))
            bw_probability *= 1 - error_rate
        end

        for edit in 1:bw_alignment.value
            bw_probability *= error_rate
        end

        if fw_probability > bw_probability
            state_probability = fw_probability
            state_orientation = true
        elseif fw_probability < bw_probability
            state_probability = bw_probability
            state_orientation = false
        else fw_probability == bw_probability
            state_probability = fw_probability
            state_orientation = missing
        end
        state_likelihoods[vertex, 1] = state_probability
        arrival_paths[vertex, 1] = [vertex => state_orientation]
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function run_viterbi!(
        current_state,
        prior_state,
        observed_nucleotide,
        observed_quality_score,
        observed_error_rate,
        current_vertex,
        prior_vertex,
        state_likelihoods,
        transition_likelihoods,
        shortest_paths,
        arrival_paths,
        kmer_graph,
        kmer_likelihoods
        )
    # if probability of prior state is lower than current probability, skip

#     @show current_state
#     @show prior_state
#     @show current_vertex
#     @show prior_vertex
    
    
    current_state_likelihood = state_likelihoods[current_vertex, current_state]
    prior_state_likelihood = state_likelihoods[prior_vertex, prior_state]

    # if we already have a better possible path, skip calculating anything
    if prior_state_likelihood < current_state_likelihood
#         @show prior_state_likelihood < current_state_likelihood
        return
    end

    # take shortest path and assume it's the maximum likelihood path
    # this assumption seems fair because in an ideal situation
    # we're just moving to an adjacent kmer
    # and the shortest path and most likely path should be the same
    shortest_path = shortest_paths[prior_vertex][current_vertex]
    
#     no path & not considering insertion
    if isempty(shortest_path) && (prior_vertex != current_vertex)
#         @show "no path, skipping"
        return
    end
    
    # if shortest path isn't viable, exit
    if !isempty(shortest_path)
#         @show "checking if path is viable"

        terminal_orientation_prior_state = last(last(arrival_paths[prior_vertex, prior_state]))
#         @show arrival_paths[prior_vertex, prior_state]
#         @show "we were at vertex $(prior_vertex) in orientation $(terminal_orientation_prior_state)"
        candidate_edge = Graphs.Edge(shortest_path[1], shortest_path[2])
                
        if !ismissing(terminal_orientation_prior_state) && 
            !any(o -> o.source_orientation == terminal_orientation_prior_state, kmer_graph.eprops[candidate_edge][:orientations])
            
#             @show "no viable orientation matching edges detected between $(candidate_edge)"
#             @show "full candidate path was $(shortest_path)"
#             @show "orientation options were:"
#             @show kmer_graph.eprops[candidate_edge][:orientations]
            return
        end
    end
    
    # zero step path - insertion in observed sequence relative to kmer graph
    is_same_vertex = (current_vertex == prior_vertex)
    has_edge = Graphs.has_edge(kmer_graph, Graphs.Edge(prior_vertex, current_vertex))
    if is_same_vertex && has_edge
        shortest_path = [prior_vertex, current_vertex]
    end
    
    if is_same_vertex
#         @show "same vertex, considering insertion potential"
        emission_likelihood = observed_error_rate
        transition_likelihood = observed_error_rate
        state_likelihood = kmer_likelihoods[current_vertex]
        path_likelihood = prior_state_likelihood * emission_likelihood * transition_likelihood * state_likelihood
        path = [last(arrival_paths[prior_vertex, prior_state])]

        if current_state_likelihood > state_likelihoods[current_vertex, current_state]
#             @show "selecting path"
#             @show path
#             @show path_likelihood
            state_likelihoods[current_vertex, current_state] = path_likelihood
            arrival_paths[current_vertex, current_state] = path
        end
    # one or more step path - match, mismatch, or deletion in observed sequence relative to kmer graph
    elseif !isempty(shortest_path)
#         @show "path is viable!"
#         @show "considering shortest path: $(shortest_path)"

        initial_path_state = last(arrival_paths[prior_vertex, prior_state])

        path = Vector{typeof(initial_path_state)}(undef, length(shortest_path))
        path[1] = initial_path_state

        path_likelihood::Float64 = state_likelihoods[prior_vertex, prior_state]

        for i in 2:length(shortest_path)

            this_vertex = shortest_path[i]
            prior_vertex, prior_orientation = path[i-1]
            edge = Graphs.Edge(prior_vertex, this_vertex)

            possible_edge_orientations::Set{NamedTuple{(:source_orientation, :destination_orientation), Tuple{Bool, Bool}}} = kmer_graph.eprops[edge][:orientations]
            
#             @show possible_edge_orientations
            
            if !ismissing(prior_orientation)
                possible_edge_orientations = filter(o -> o.source_orientation == prior_orientation, possible_edge_orientations)
            end
            
#             @show possible_edge_orientations
            
            if isempty(possible_edge_orientations)
                path_likelihood *= 0.0
                path = Vector{eltype(path)}()
#                 @show "no possible orientations, bailing early"
                break
            end

#             @show prior_orientation
            if ismissing(prior_orientation)
                if transition_likelihoods[true][prior_vertex, this_vertex] > transition_likelihoods[false][prior_vertex, this_vertex]
                    prior_orientation = true
                    transition_likelihood = transition_likelihoods[true][prior_vertex, this_vertex]::Float64
                elseif transition_likelihoods[true][prior_vertex, this_vertex] < transition_likelihoods[false][prior_vertex, this_vertex]
                    prior_orientation = false
                    transition_likelihood = transition_likelihoods[false][prior_vertex, this_vertex]::Float64
                else transition_likelihoods[true][prior_vertex, this_vertex] == transition_likelihoods[false][prior_vertex, this_vertex]
                    prior_orientation = missing
                    transition_likelihood = transition_likelihoods[true][prior_vertex, this_vertex]::Float64
                end
            else
                transition_likelihood = transition_likelihoods[prior_orientation][prior_vertex, this_vertex]::Float64
            end
            state_likelihood::Float64 = kmer_likelihoods[this_vertex]
            path_likelihood *= transition_likelihood * state_likelihood
            
            if length(possible_edge_orientations) == 1
                orientation = first(possible_edge_orientations).destination_orientation
                path[i] = this_vertex => orientation
            else
                path[i] = this_vertex => missing
            end
        end

        # see if new nucleotide is a match or mismatch to terminal kmer in path
        if !isempty(path) && path_likelihood > 0
            terminal_kmer_index, terminal_kmer_orientation = last(path)
            terminal_kmer = BioSequences.LongDNASeq(kmer_graph.vprops[terminal_kmer_index][:kmer])::BioSequences.LongDNASeq
            if ismissing(terminal_kmer_orientation)
                fw_is_match = observed_nucleotide == last(terminal_kmer)
                bw_is_match = observed_nucleotide == last(BioSequences.reverse_complement!(terminal_kmer))
                if fw_ismatch && !bw_is_match
                    path[end] = terminal_kmer_index => true
                    path_likelihood *= 1 - observed_error_rate
                elseif !fw_ismatch && bw_is_match
                    path[end] = terminal_kmer_index => false
                    path_likelihood *= 1 - observed_error_rate
                elseif fw_ismatch && bw_is_match
                    path_likelihood *= 1 - observed_error_rate
                elseif !fw_ismatch && !bw_is_match
                    path_likelihood *= observed_error_rate
                end
            elseif terminal_kmer_orientation
                is_match = observed_nucleotide == last(terminal_kmer)
                if is_match
                    path_likelihood *= 1 - observed_error_rate
                else
                    path_likelihood *= observed_error_rate
                end
            else
                terminal_kmer = BioSequences.reverse_complement!(terminal_kmer)
                is_match = observed_nucleotide == last(terminal_kmer)
                if is_match
                    path_likelihood *= 1 - observed_error_rate
                else
                    path_likelihood *= observed_error_rate
                end
            end
        end

        if path_likelihood > state_likelihoods[current_vertex, current_state]
#             @show "selecting path"
#             @show path
#             @show path_likelihood
            state_likelihoods[current_vertex, current_state] = path_likelihood
            arrival_paths[current_vertex, current_state] = path
        end
    end
    return
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function polish_fastq(kmer_graph, fastq_file)

#     @info "Assessing kmer likelihoods"
    kmers = [kmer_graph.vprops[v][:kmer] for v in Graphs.vertices(kmer_graph)]
#     kmer_counts = [length(kmer_graph.vprops[v][:evidence]) for v in Graphs.vertices(kmer_graph)]
    kmer_counts = [kmer_graph.vprops[v][:weight] for v in Graphs.vertices(kmer_graph)]
    kmer_likelihoods = kmer_counts ./ sum(kmer_counts)
    k = kmer_graph.gprops[:k]
    kmer_type = BioSequences.BigDNAMer{k}
    total_kmers = length(kmers)
    
#     @info "determining shortest paths between kmers"
    shortest_paths = Graphs.enumerate_paths(Graphs.floyd_warshall_shortest_paths(kmer_graph));
    
    @info "counting the number of records to establish runtime estimate"
    number_of_records = 0
    for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
        number_of_records += 1
    end
    progress_bar = ProgressMeter.Progress(number_of_records, 1)
    
    output_fastq_file = replace(fastq_file, ".fastq" => ".k$(kmer_graph.gprops[:k]).fastq")
    fastq_writer = FASTX.FASTQ.Writer(open(output_fastq_file, "w"))
    for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
        ProgressMeter.next!(progress_bar)
        
#         @info "Initializing matrices"
        total_states = length(FASTX.sequence(fastq_record))-k+1
        transition_likelihoods = initialize_transition_probabilities(kmer_graph)
        state_likelihoods = zeros(total_kmers, total_states)
        arrival_paths = fill(Pair{Int, Union{Bool, Missing}}[], total_kmers, total_states)

#         @info "Determining Likelihoods of initial states"
        initial_state = first(BioSequences.each(kmer_type, FASTX.sequence(fastq_record)))
        current_state = 1
        # note this is a place for potential improvement, use the q value at each base to guide probability rather than median
        median_q_value = Statistics.median(Int.(FASTX.quality(fastq_record)[1:k]))
        current_error_rate = q_value_to_error_rate(median_q_value)
        # canonical_kmer = BioSequences.canonical(initial_state.fw)
        set_initial_state_likelihoods!(
                kmer_graph,
                initial_state,
                kmer_likelihoods,
                current_error_rate,
                state_likelihoods,
                arrival_paths
            )

#         @info "Determining likelihood of downstream states"

#         non_singleton_states = findall(kmer_counts .> 1)

        for current_state in 2:total_states
            prior_state = current_state - 1

        #     observed_kmer = BioSequences.BigDNAMer{k}(FASTX.sequence(fastq_record)[current_state:current_state+k-1])

        #     @assert observed_kmer == collect(BioSequences.each(kmer_type, FASTX.sequence(fastq_record)))[current_state].fw

        #     canonical_kmer = BioSequences.canonical(observed_kmer)

            observed_nucleotide = FASTX.sequence(fastq_record)[k-1+current_state]
        #     observed_nucleotide = last(observed_kmer)
            observed_quality_score = FASTX.quality(fastq_record)[k-1+current_state]
            observed_error_rate = q_value_to_error_rate(observed_quality_score)

            # we'll assess prior states in order of decreasing likelihood
            # such that we maximize how frequently we are able to utilize the
            # current_state_likelihood > candidate prior state
            # break that won't bother evaluating lower likelihood possibilities
            prior_states_in_decreasing_likelihood = sortperm(state_likelihoods[:, prior_state], rev=true)

            # and skip all prior states with zero probability

            for current_vertex in total_states
                for prior_vertex in prior_states_in_decreasing_likelihood
                    if state_likelihoods[prior_vertex, prior_state] > 0
                        run_viterbi!(
                                current_state,
                                prior_state,
                                observed_nucleotide,
                                observed_quality_score,
                                observed_error_rate,
                                current_vertex,
                                prior_vertex,
                                state_likelihoods,
                                transition_likelihoods,
                                shortest_paths,
                                arrival_paths,
                                kmer_graph,
                                kmer_likelihoods
                                )
                    end
                end
            end
        end

#         try
        maximum_likelihood_path, maximum_likelihood_value = 
            determine_maximum_likelihood_path(
                state_likelihoods,
                arrival_paths
                )
#         catch
#             return state_likelihoods, arrival_paths
#         end

        sequence = oriented_path_to_sequence(kmer_graph, maximum_likelihood_path)

#         @info "comparing to original path"
        original_sequence_likelihood = oriented_path_to_likelihood(kmer_graph, kmers, kmer_likelihoods, transition_likelihoods, fastq_record)
        relative_likelihood = maximum_likelihood_value / original_sequence_likelihood
#         relative_likelihood_formatted = NumericIO.formatted(relative_likelihood, ndigits=1, charset=:ASCII)
#         println("relative likelihood of new path to old path is $(relative_likelihood_formatted)")

#         @info "writing updated record"
        identifier = FASTX.identifier(fastq_record) * "_k$(k)"
        description = string(relative_likelihood)
        # because the sequences won't always be the same length, we take an ordered sampling with replacement
        # which introduces some random error but preserves overall patterns and areas of high/low accuracy
        quality_scores = StatsBase.sample(FASTX.quality(fastq_record), length(sequence), ordered=true)

        new_fastq_record = FASTX.FASTQ.Record(
            identifier,
            description,
            sequence,
            quality_scores
        )
        write(fastq_writer, new_fastq_record)
    end
    close(fastq_writer)
    return output_fastq_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_maximum_likelihood_path(
    state_likelihoods,
    arrival_paths
    )
    maximum_likelihood_value = maximum(state_likelihoods[:, end])

    maximum_likelihood_path_indices = findall(state_likelihoods[:, end] .== maximum_likelihood_value)

    # if multiple paths are tied, randomly choose one
    maximum_likelihood_path_index = rand(maximum_likelihood_path_indices)

    maximum_likelihood_path = arrival_paths[maximum_likelihood_path_index, end]

    for state_index in size(state_likelihoods, 2)-1:-1:1
        next_kmer, next_orientation = first(maximum_likelihood_path)
        maximum_likelihood_arrival_path = arrival_paths[next_kmer, state_index]
        
        is_match = last(maximum_likelihood_arrival_path) == (next_kmer => next_orientation)
        if !ismissing(is_match) && !is_match
            error("breaking")
        end
        maximum_likelihood_path = vcat(maximum_likelihood_arrival_path[1:end-1], maximum_likelihood_path)
    end
    return maximum_likelihood_path, maximum_likelihood_value
end