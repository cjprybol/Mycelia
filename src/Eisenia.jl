module Eisenia

import BioAlignments
import BioSequences
import BioSymbols
import DataStructures
import DocStringExtensions
import GraphRecipes
import LightGraphs
import Plots
import PrettyTables
import Random
import SparseArrays
import StatsBase
import StatsPlots

# preserve definitions between code jldoctest code blocks
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Preserving-Definitions-Between-Blocks
# use this to build up a story as we go, where outputs of earlier defined functions feed into
# tests for further downstream functions

# if things fall out of date but look correct, update them automatically
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Fixing-Outdated-Doctests

const DNA_ALPHABET = BioSymbols.ACGT
const RNA_ALPHABET = BioSymbols.ACGU
const AA_ALPHABET = filter(
    x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x) || BioSymbols.isterm(x)),
    BioSymbols.alphabet(BioSymbols.AminoAcid))


"""
$(DocStringExtensions.TYPEDEF)
$(DocStringExtensions.TYPEDFIELDS)

A short description of the Type

```jldoctest
julia> 1 + 1
2
```
"""
struct OrientedKmer
    index::Int
    orientation::Union{Missing, Bool}
    function OrientedKmer(;index, orientation)
        return new(index, orientation)
    end
end


"""
$(DocStringExtensions.TYPEDEF)
$(DocStringExtensions.TYPEDFIELDS)

A short description of the Type

```jldoctest
julia> 1 + 1
2
```
"""
struct EdgeEvidence
    observation_index::Int
    edge_index::Int
    function EdgeEvidence(;observation_index, edge_index)
        return new(observation_index, edge_index)
    end
end


"""
$(DocStringExtensions.TYPEDEF)
$(DocStringExtensions.TYPEDFIELDS)

A short description of the Type

```jldoctest
julia> 1 + 1
2
```
"""
struct KmerGraph{KmerType}
    graph::LightGraphs.SimpleGraphs.SimpleGraph{Int}
    edge_evidence::Dict{LightGraphs.SimpleGraphs.SimpleEdge{Int}, Vector{EdgeEvidence}}
    kmers::AbstractVector{KmerType}
    counts::AbstractVector{Int}
    function KmerGraph(;graph, edge_evidence, kmers::AbstractVector{KmerType}, counts) where {KmerType <: BioSequences.AbstractMer}
        new{KmerType}(graph, edge_evidence, kmers, counts)
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
function observe(sequence::BioSequences.LongSequence{T}; error_rate = 0.0) where T
    
    if T <: BioSequences.DNAAlphabet
        alphabet = DNA_ALPHABET
    elseif T <: BioSequences.RNAAlphabet
        alphabet = RNA_ALPHABET
    else
        @assert T <: BioSequences.AminoAcidAlphabet
        alphabet = AA_ALPHABET
    end
    
    new_seq = BioSequences.LongSequence{T}()
    for character in sequence
        if rand() > error_rate
            # match
            push!(new_seq, character)
        else
            error_type = rand(1:3)
            if error_type == 1
                # mismatch
                push!(new_seq, rand(setdiff(alphabet, character)))
            elseif error_type == 2
                # insertion
                push!(new_seq, rand(alphabet))
                push!(new_seq, character)
            else
                # deletion
                continue
            end
        end
    end
    if (T <: BioSequences.DNAAlphabet || T <: BioSequences.RNAAlphabet) && rand(Bool)
        BioSequences.reverse_complement!(new_seq)
    end
    return new_seq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function get_kmer_index(kmers, kmer)
    index_range = searchsorted(kmers, kmer)
    if !isempty(index_range)
        return first(index_range)
    else
        return nothing
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
function ordered_edge(a, b)
    if a <= b
        return LightGraphs.Edge(a, b)
    else
        return LightGraphs.Edge(b, a)
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
function KmerGraph(::Type{KMER_TYPE}, observations) where {KMER_TYPE <: BioSequences.AbstractMer{A, K}} where {A, K}
    
    if !isodd(K)
        error("Even kmers are not supported")
    end

    kmer_counts = count_kmers(KMER_TYPE, observations)
    kmers = collect(keys(kmer_counts))
    counts = collect(values(kmer_counts))

    # initalize graph
    graph = LightGraphs.SimpleGraph(length(kmers))

    # an individual piece of evidence takes the form of
    # (observation_index = observation #, edge_index = edge # starting from beginning of the observation)
    
    # evidence takes the form of Edge => [(evidence_1), (evidence_2), ..., (evidence_N)]    
    edge_evidence = Dict{LightGraphs.SimpleGraphs.SimpleEdge{Int}, Vector{EdgeEvidence}}()
    for (observation_index, observation) in enumerate(observations)
        for edge_index in 1:length(observation)-K
            a_to_b_connection = observation[edge_index:edge_index+K]
            a = BioSequences.canonical(KMER_TYPE(a_to_b_connection[1:end-1]))
            b = BioSequences.canonical(KMER_TYPE(a_to_b_connection[2:end]))
            a_index = get_kmer_index(kmers, a)
            b_index = get_kmer_index(kmers, b)
            edge = ordered_edge(a_index, b_index)
            LightGraphs.add_edge!(graph, edge)
            evidence = EdgeEvidence(;observation_index, edge_index)
            edge_evidence[edge] = push!(get(edge_evidence, edge, EdgeEvidence[]), evidence)
        end
    end
    return KmerGraph(;graph, edge_evidence, kmers, counts)
end


LightGraphs.has_edge(kmer_graph::KmerGraph, edge) = LightGraphs.has_edge(kmer_graph.graph, edge)
LightGraphs.has_path(kmer_graph::KmerGraph, u, v) = LightGraphs.has_path(kmer_graph.graph, u, v)

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

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_edge_probabilities(graph, strand)
    outgoing_edge_probabilities = SparseArrays.spzeros(length(graph.kmers), length(graph.kmers))
    
    for (kmer_index, kmer) in enumerate(graph.kmers)
        if !strand
            kmer = BioSequences.reverse_complement(kmer)
        end
        
        downstream_neighbor_indices = Int[]
        for neighbor in BioSequences.neighbors(kmer)
            index = get_kmer_index(graph.kmers, BioSequences.canonical(neighbor))
            # kmer must be in our dataset and there must be a connecting edge
            if !isnothing(index) && LightGraphs.has_edge(graph, ordered_edge(kmer_index, index))
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
function assess_alignment(a, b)
    pairwise_alignment = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), a, b)
    alignment_result = BioAlignments.alignment(pairwise_alignment)
    total_aligned_bases = BioAlignments.count_aligned(alignment_result)
    total_matches = Int(BioAlignments.count_matches(alignment_result))
    total_edits = Int(total_aligned_bases - total_matches)
    return (total_matches = total_matches, total_edits = total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_suffix_match(kmer_a, kmer_b, orientation)
    if !orientation
        kmer_b = BioSequences.reverse_complement(kmer_b)
    end
    la = length(kmer_a)
    lb = length(kmer_b)
    all_match = true
    for (ia, ib) in zip(2:la, 1:lb-1)
        this_match = kmer_a[ia] == kmer_b[ib]
        all_match &= this_match
    end
    return all_match
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_path_orientations(path, kmers, initial_orientation)
    
    orientations = Vector{Union{Bool, Missing}}(missing, length(path))
    orientations[1] = initial_orientation
    
    for (i, (a, b)) in enumerate(zip(path[1:end-1], path[2:end]))
        kmer_a = kmers[a]
        kmer_b = kmers[b]
        if !ismissing(orientations[i]) && !orientations[i]
            kmer_a = BioSequences.reverse_complement(kmer_a)
        end
        forward_match = assess_suffix_match(kmer_a, kmer_b, true)
        reverse_match = assess_suffix_match(kmer_a, kmer_b, false)
        if forward_match && reverse_match
            # ambiguous orientation
            this_orientation = missing
        elseif forward_match && !reverse_match
            this_orientation = true
        elseif !forward_match && reverse_match
            this_orientation = false
        else
#             error("neither orientation matches $kmer_a, $kmer_b")
            # I think this only applies when the path is possible but in the wrong orientation
            return nothing
        end
        orientations[i+1] = this_orientation
    end
    return orientations
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_path_likelihood(oriented_path, kmers, counts, outgoing_edge_probabilities, incoming_edge_probabilities)
    # we take it as a given that we are at the current node, so initialize with p=1 and then update with other node likelihoods
    likelihood = 1.0
    total_kmer_count = sum(counts)
    for node in oriented_path[2:end]
        likelihood *= counts[node.index] / total_kmer_count
    end
    for (a, b) in zip(oriented_path[1:end-1], oriented_path[2:end])
        if ismissing(a.orientation)
            # ambiguous orientation
            likelihood *= max(incoming_edge_probabilities[a.index, b.index], incoming_edge_probabilities[a.index, b.index])
        elseif a.orientation
            likelihood *= outgoing_edge_probabilities[a.index, b.index]
        elseif !a.orientation
            likelihood *= incoming_edge_probabilities[a.index, b.index]
        else
            error("unreachable condition")
        end
    end
    return likelihood
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function orient_path(path, orientations)
    oriented_path = OrientedKmer[OrientedKmer(index = i, orientation = o) for (i, o) in zip(path, orientations)] 
    return oriented_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_emission_match(a, b, orientation)
    if !orientation
        a = BioSequences.reverse_complement(a)
    end
    return a[end] == b[end]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_emission(current_orientation, current_kmer_index, observed_kmer, kmers)
    
    current_kmer_sequence = kmers[current_kmer_index]
    if observed_kmer.orientation
        observed_kmer_sequence = kmers[observed_kmer.index]
    else
        observed_kmer_sequence = BioSequences.reverse_complement(kmers[observed_kmer.index])
    end
    
    if !ismissing(current_orientation)
        emission_match = assess_emission_match(current_kmer_sequence, observed_kmer_sequence, current_orientation)
        evaluated_orientation = current_orientation
    else
        forward_emission_match = assess_emission_match(current_kmer_sequence, observed_kmer_sequence, true)
        reverse_emission_match = assess_emission_match(current_kmer_sequence, observed_kmer_sequence, false)
        emission_match = forward_emission_match || reverse_emission_match
        if forward_emission_match && !reverse_emission_match
            evaluated_orientation = true
        elseif !forward_emission_match && reverse_emission_match
            evaluated_orientation = false
        else
            evaluated_orientation = missing
        end
    end
    return (emission_match, evaluated_orientation)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_path(path,
    observed_kmer,
    kmers,
    counts,
    initial_orientation,
    outgoing_edge_probabilities,
    incoming_edge_probabilities,
    error_rate)
    
    orientations = assess_path_orientations(path, kmers, initial_orientation)
    if orientations == nothing
        path_likelihood = 0.0
        edit_distance = Inf
        oriented_path = OrientedKmer[]
    else
        emission_match, evaluated_orientation = assess_emission(last(orientations), last(path), observed_kmer, kmers)
        # assert an orientation if we found one
        if ismissing(last(orientations)) && !ismissing(evaluated_orientation)
            orientations[end] = evaluated_orientation
        end
        oriented_path = orient_path(path, orientations)
        path_likelihood = assess_path_likelihood(oriented_path, kmers, counts, outgoing_edge_probabilities, incoming_edge_probabilities)

        insertion_deletion_magnitue = abs(length(path) - 2)
        edit_distance = !emission_match + insertion_deletion_magnitue

        path_likelihood *= edit_distance == 0 ? (1.0 - error_rate) : (error_rate^edit_distance)
    end
    return (oriented_path = oriented_path, path_likelihood = path_likelihood, edit_distance = edit_distance)    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function find_outneighbors(orientation, kmer_index, outgoing_edge_probabilities, incoming_edge_probabilities)
    if ismissing(orientation)
        outneighbors = vcat(
            first(SparseArrays.findnz(outgoing_edge_probabilities[kmer_index, :])),
            first(SparseArrays.findnz(incoming_edge_probabilities[kmer_index, :]))
        )
    elseif orientation
        outneighbors = first(SparseArrays.findnz(outgoing_edge_probabilities[kmer_index, :]))
    else
        outneighbors = first(SparseArrays.findnz(incoming_edge_probabilities[kmer_index, :]))
    end
    return filter!(x -> x != kmer_index, unique!(outneighbors))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_insertion(previous_orientation, current_kmer_index, observed_kmer, kmers, counts, error_rate)    
    
    # I'm not sure if I should evaluate the emission or just accept it as an error
    # here we accept it as wrong without 
    transition_likelihood = emission_likelihood = error_rate
    edit_distance = 2
    oriented_path = [OrientedKmer(index = current_kmer_index, orientation = previous_orientation)]
    
#     # here we actually check
#     emission_match, evaluated_orientation = assess_emission(previous_orientation, current_kmer_index, observed_kmer, kmers)
#     oriented_path = [OrientedKmer(index = current_kmer_index, orientation = evaluated_orientation)]
#     transition_likelihood = error_rate
#     emission_likelihood = emission_match ? (1.0 - error_rate) : error_rate
#     transition_edit_distance = 1
#     edit_distance = !emission_match + 1
    
    # universal downstream
    state_likelihood = counts[current_kmer_index] / sum(counts)
    path_likelihood = state_likelihood * transition_likelihood * emission_likelihood
    
    return (oriented_path = oriented_path, path_likelihood = path_likelihood, edit_distance = edit_distance)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_alignment_accuracy(alignment_result)
    return alignment_result.total_matches / (alignment_result.total_matches + alignment_result.total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_optimal_alignment(kmer, observed_kmer)

    forward_alignment_result = assess_alignment(kmer, observed_kmer)
    forward_alignment_accuracy = assess_alignment_accuracy(forward_alignment_result)

    reverse_alignment_result = assess_alignment(kmer, BioSequences.reverse_complement(observed_kmer))
    reverse_alignment_accuracy = assess_alignment_accuracy(reverse_alignment_result)

    if forward_alignment_accuracy > reverse_alignment_accuracy
        alignment_result = forward_alignment_result
        orientation = true
    elseif forward_alignment_accuracy < reverse_alignment_accuracy
        alignment_result = reverse_alignment_result
        orientation = false
    elseif forward_alignment_accuracy == reverse_alignment_accuracy
        alignment_result, orientation = rand(((forward_alignment_result, missing), (reverse_alignment_result, missing)))
    end

    return (alignment_result, orientation)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function count_kmers(::Type{KMER_TYPE}, sequence::BioSequences.LongSequence) where KMER_TYPE
    canonical_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    canonical_kmer_iterator = (BioSequences.canonical(kmer.fw) for kmer in BioSequences.each(KMER_TYPE, sequence))
    for canonical_kmer in canonical_kmer_iterator
        canonical_kmer_counts[canonical_kmer] = get(canonical_kmer_counts, canonical_kmer, 0) + 1
    end
    return canonical_kmer_counts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function count_kmers(::Type{KMER_TYPE}, sequences) where KMER_TYPE
    joint_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    for sequence in sequences
        sequence_kmer_counts = count_kmers(KMER_TYPE, sequence)
        merge!(+, joint_kmer_counts, sequence_kmer_counts)
    end
    sort!(joint_kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function orient_oriented_kmer(kmers, kmer)
    oriented_kmer_sequence = kmers[kmer.index]
    if !ismissing(kmer.orientation) && !kmer.orientation
        oriented_kmer_sequence = BioSequences.reverse_complement(oriented_kmer_sequence)
    end
    return oriented_kmer_sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function oriented_path_to_sequence(oriented_path, kmers)
    initial_kmer_sequence = orient_oriented_kmer(kmers, oriented_path[1])
    
    sequence = BioSequences.LongDNASeq(initial_kmer_sequence)
    for i in 2:length(oriented_path)
        this_kmer_sequence = orient_oriented_kmer(kmers, oriented_path[i])
        push!(sequence, this_kmer_sequence[end])
    end
    return sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function sequence_to_oriented_path(sequence, kmers::Vector{T}) where {T <: BioSequences.AbstractMer{A, K}} where {A, K}
    observed_path = Vector{OrientedKmer}(undef, length(sequence)-K+1)
    for (i, kmer) in enumerate(BioSequences.each(T, sequence))
        canonical_kmer = BioSequences.canonical(kmer.fw)
        index = get_kmer_index(kmers, canonical_kmer)
        orientation = kmer.fw == canonical_kmer
        observed_path[i] = OrientedKmer(index = index, orientation = orientation)
    end
    return observed_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function find_optimal_path(observed_kmer,
    previous_kmer_index,
    previous_orientation,
    current_kmer_index,
    graph, 
    shortest_paths, 
    outgoing_edge_probabilities, 
    incoming_edge_probabilities,
    error_rate)

    path_likelihood = 0.0
    oriented_path = Vector{OrientedKmer}()
    edit_distance = 0

    if current_kmer_index == previous_kmer_index
        # could be a self loop, check this first
        if LightGraphs.has_edge(graph, LightGraphs.Edge(previous_kmer_index, current_kmer_index))
            this_path = [previous_kmer_index, current_kmer_index]
            this_oriented_path, this_likelihood, this_edit_distance = 
                assess_path(this_path,
                            observed_kmer,
                            graph.kmers,
                            graph.counts,
                            previous_orientation,
                            outgoing_edge_probabilities,
                            incoming_edge_probabilities,
                            error_rate)
            if this_likelihood > path_likelihood
                path_likelihood = this_likelihood
                oriented_path = this_oriented_path
                edit_distance = this_edit_distance
            end
        end
        
        # consider an insertion in observed sequence relative to the reference graph
        this_oriented_path, this_likelihood, this_edit_distance =
            assess_insertion(previous_orientation, current_kmer_index, observed_kmer, graph.kmers, graph.counts, error_rate)
        if this_likelihood > path_likelihood
            path_likelihood = this_likelihood
            oriented_path = this_oriented_path
            edit_distance = this_edit_distance
        end
        
        # consider a deletion in observed sequene relative to the reference graph
        # see if this has any neighbors that circle back, and evaluate the likelihood for each
        outneighbors = find_outneighbors(previous_orientation, previous_kmer_index, outgoing_edge_probabilities, incoming_edge_probabilities)
        for outneighbor in outneighbors
            if LightGraphs.has_path(graph, outneighbor, current_kmer_index)
                # manually build path
                this_path = [previous_kmer_index, shortest_paths[previous_kmer_index][current_kmer_index]...]
                this_oriented_path, this_likelihood, this_edit_distance = 
                    assess_path(this_path,
                                observed_kmer,
                                graph.kmers,
                                graph.counts,
                                previous_orientation,
                                outgoing_edge_probabilities,
                                incoming_edge_probabilities,
                                error_rate)
                if this_likelihood > path_likelihood
                    path_likelihood = this_likelihood
                    oriented_path = this_oriented_path
                    edit_distance = this_edit_distance
                end
            end
        end
    elseif LightGraphs.has_path(graph, previous_kmer_index, current_kmer_index)   
        this_path = shortest_paths[previous_kmer_index][current_kmer_index]
        this_oriented_path, this_likelihood, this_edit_distance = 
            assess_path(this_path,
                        observed_kmer,
                        graph.kmers,
                        graph.counts,
                        previous_orientation,
                        outgoing_edge_probabilities,
                        incoming_edge_probabilities,
                        error_rate)
        if this_likelihood > path_likelihood
            path_likelihood = this_likelihood
            oriented_path = this_oriented_path
            edit_distance = this_edit_distance
        end
    end
    
    return (oriented_path = oriented_path, path_likelihood = path_likelihood, edit_distance = edit_distance)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function backtrack_optimal_path(kmer_likelihoods, arrival_paths, edit_distances)
    
    maximum_likelihood = maximum(kmer_likelihoods[:, end])

    maximum_likelihood_path_index = Int(rand(findall(x -> x == maximum_likelihood, kmer_likelihoods[:, end])))::Int

    maximum_likelihood_edit_distance = edit_distances[maximum_likelihood_path_index, end]

    maximum_likelihood_path = arrival_paths[maximum_likelihood_path_index, end]::Vector{OrientedKmer}


    for observed_path_index in size(kmer_likelihoods, 2)-1:-1:1
        maximum_likelihood_arrival_path = arrival_paths[first(maximum_likelihood_path).index, observed_path_index][1:end-1]
        maximum_likelihood_path = vcat(maximum_likelihood_arrival_path, maximum_likelihood_path)
    end
    
    unlogged_maximum_likelihood = MathConstants.e^BigFloat(maximum_likelihood)
    unlogged_total_likelihood = sum(likelihood -> MathConstants.e^BigFloat(likelihood), kmer_likelihoods[:, end])
    relative_likelihood = Float64(unlogged_maximum_likelihood / unlogged_total_likelihood)
    
    return (
        maximum_likelihood_path,
        maximum_likelihood_edit_distance,
        relative_likelihood)
    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function initialize_viterbi(graph, observed_path, error_rate)

    edit_distances = Array{Union{Int, Missing}}(missing, LightGraphs.nv(graph.graph), length(observed_path))
    arrival_paths = Array{Union{Vector{OrientedKmer}, Missing}}(missing, LightGraphs.nv(graph.graph), length(observed_path))
    kmer_likelihoods = Array{Float64}(undef, LightGraphs.nv(graph.graph), length(observed_path)) .= -Inf
    kmer_likelihoods[:, 1] .= graph.counts ./ sum(graph.counts)

    observed_kmer_sequence = orient_oriented_kmer(graph.kmers, first(observed_path))
    for (kmer_index, kmer) in enumerate(graph.kmers)
        
        alignment_result, orientation = assess_optimal_alignment(kmer, observed_kmer_sequence)

        for match in 1:alignment_result.total_matches
            kmer_likelihoods[kmer_index, 1] *= (1.0 - error_rate)
        end
        for edit in 1:alignment_result.total_edits
            kmer_likelihoods[kmer_index, 1] *= error_rate
        end
        if kmer_likelihoods[kmer_index, 1] > 0.0
            arrival_paths[kmer_index, 1] = [OrientedKmer(index = kmer_index, orientation = orientation)]
            edit_distances[kmer_index, 1] = alignment_result.total_edits
        end
    end
    kmer_likelihoods[:, 1] ./= sum(kmer_likelihoods[:, 1])
    kmer_likelihoods[:, 1] .= log.(kmer_likelihoods[:, 1])
    
    return (kmer_likelihoods = kmer_likelihoods, arrival_paths = arrival_paths, edit_distances = edit_distances)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function viterbi_maximum_likelihood_path(graph, observation, error_rate; debug = false)

    observed_path = sequence_to_oriented_path(observation, graph.kmers)
    
    outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
    
    shortest_paths = LightGraphs.enumerate_paths(LightGraphs.floyd_warshall_shortest_paths(graph.graph))
    
    kmer_likelihoods, arrival_paths, edit_distances = initialize_viterbi(graph, observed_path, error_rate)

    if debug
        my_show(arrival_paths, graph.kmers, title = "Arrival Paths")
        my_show(kmer_likelihoods, graph.kmers, title="Kmer Likelihoods")
        my_show(edit_distances, graph.kmers, title="Edit Distances")
    end
    
    for current_observation_index in 2:length(observed_path)
        observed_kmer = observed_path[current_observation_index]
        if debug
            println("current_observation_index = $(current_observation_index)")
            println("observed_kmer = $(observed_kmer)")
        end
        for (current_kmer_index, current_kmer) in enumerate(graph.kmers)
            if debug
                println("\t\tcurrent_kmer_index = $(current_kmer_index)")
                println("\t\tcurrent_kmer = $(current_kmer)")
            end
            for (previous_kmer_index, previous_kmer) in enumerate(graph.kmers)
                if debug
                    println("\tprevious_kmer_index = $(previous_kmer_index)")
                    println("\tprevious_kmer = $(previous_kmer)")
                end

                previous_likelihood = kmer_likelihoods[previous_kmer_index, current_observation_index - 1]
                current_likelihood = kmer_likelihoods[current_kmer_index, current_observation_index]
                previous_arrival_path = arrival_paths[previous_kmer_index, current_observation_index - 1]
                if current_likelihood > previous_likelihood && !ismissing(previous_arrival_path)
                    # we've determined that there is no way to get to this kmer being evaluated
                    continue
                else
                    if ismissing(previous_arrival_path)
                        previous_orientation = missing
                    else
                        previous_orientation = last(arrival_paths[previous_kmer_index, current_observation_index - 1]).orientation
                    end
                    optimal_path_result =
                        find_optimal_path(observed_kmer,
                            previous_kmer_index,
                            previous_orientation,
                            current_kmer_index,
                            graph, 
                            shortest_paths,
                            outgoing_edge_probabilities, 
                            incoming_edge_probabilities,
                            error_rate)
                    if debug
                        for propertyname in propertynames(optimal_path_result)
                            println("\t\t\t$propertyname = $(getproperty(optimal_path_result, propertyname))")
                        end
                    end
                    this_likelihood = previous_likelihood + log(optimal_path_result.path_likelihood)

                    if this_likelihood > current_likelihood
                        kmer_likelihoods[current_kmer_index, current_observation_index] = this_likelihood
                        arrival_paths[current_kmer_index, current_observation_index] = optimal_path_result.oriented_path
                        previous_edit_distance = edit_distances[previous_kmer_index, current_observation_index - 1]
                        edit_distances[current_kmer_index, current_observation_index] = previous_edit_distance + optimal_path_result.edit_distance 
                    end
                end
            end
        end
    end
    if debug
        my_show(arrival_paths, graph.kmers, title = "Arrival Paths")
        my_show(kmer_likelihoods, graph.kmers, title="Kmer Likelihoods")
        my_show(edit_distances, graph.kmers, title="Edit Distances")
    end
    return backtrack_optimal_path(kmer_likelihoods, arrival_paths, edit_distances)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function plot_kmer_frequency_spectra(counts)
    kmer_counts_hist = StatsBase.countmap(c for c in counts)
    xs = collect(keys(kmer_counts_hist))
    ys = collect(values(kmer_counts_hist))

    StatsPlots.plot(
        xs,
        ys,
        xlims = (0, maximum(xs) + 1),
        ylims = (0, maximum(ys) + 1),
        seriestype = :scatter,
        legend = false,
        xlabel = "observed frequency",
        ylabel = "# of kmers"
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function my_show(x::Dict{LightGraphs.SimpleGraphs.SimpleEdge{Int},Array{NamedTuple{(:observation_index, :edge_index),Tuple{Int,Int}},1}})
    for (k, vs) in x
        println(k)
        for v in vs
            println("\t$v")
        end
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
function my_show(vector::AbstractVector{T}; kwargs...) where T <: OrientedKmer
    return PrettyTables.pretty_table(
    vector,
    tf = PrettyTables.tf_markdown;
    kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function my_show(array::AbstractMatrix, kmers; kwargs...)
    return PrettyTables.pretty_table(
    array,
    ["$i" for i in 1:size(array, 2)],
    row_names = kmers,
    tf = PrettyTables.tf_markdown;
    kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function my_show(arrival_paths::AbstractMatrix{T}, kmers; kwargs...) where {T <: Union{Missing, Vector{OrientedKmer}}}
    string_arrival_paths = Array{Union{String, Missing}, 2}(missing, size(arrival_paths)...)
    for index in eachindex(arrival_paths)
        cell = arrival_paths[index]
        if !ismissing(cell)
            cell_strings = []
            for node in cell
                node_string = "(" * join(["$(getfield(node, field))" for field in propertynames(node)], ", ") * ")"
                push!(cell_strings, node_string)
            end
            cell_string = '[' * join(cell_strings, ", ") * ']'
            string_arrival_paths[index] = cell_string
        end
    end
    return my_show(string_arrival_paths, kmers; kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function plot_graph(graph)
    graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))

    p = GraphRecipes.graphplot(
        graph.graph,
        names = 1:length(graph.kmers),
        node_weights = graph.counts,
        markersize = 0.2,
        hover=false,
        fontsize=12)
end

end # module