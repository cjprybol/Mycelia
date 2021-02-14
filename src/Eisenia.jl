module Eisenia

import BioAlignments
import BioSequences
import BioSymbols
import DataStructures
import Distributions
import DocStringExtensions
import GraphRecipes
import LightGraphs
import Plots
import PrettyTables
import Random
import SparseArrays
import Statistics
import StatsBase
import StatsPlots
import FASTX

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
$(DocStringExtensions.TYPEDSIGNATURES)

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

function KmerGraph(::Type{KMER_TYPE}, observations) where {KMER_TYPE <: BioSequences.AbstractMer{A, K}} where {A, K}
    
    if !isodd(K)
        error("Even kmers are not supported")
    end

    kmer_counts = count_kmers(KMER_TYPE, observations)
    kmers = collect(keys(kmer_counts))
    counts = collect(values(kmer_counts))

    return KmerGraph(KMER_TYPE, observations, kmers, counts)
end

function KmerGraph(::Type{KMER_TYPE}, observations, kmers, counts) where {KMER_TYPE <: BioSequences.AbstractMer{A, K}} where {A, K}
    
    if !isodd(K)
        error("Even kmers are not supported")
    end
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
            if (a_index != nothing) && (b_index != nothing)
                edge = ordered_edge(a_index, b_index)
                LightGraphs.add_edge!(graph, edge)
                evidence = EdgeEvidence(;observation_index, edge_index)
                edge_evidence[edge] = push!(get(edge_evidence, edge, EdgeEvidence[]), evidence)
            end
        end
    end
    return KmerGraph(;graph, edge_evidence, kmers, counts)
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
                total_insertions = 1 + rand(Distributions.Poisson(error_rate))
                for i in 1:total_insertions
                    push!(new_seq, rand(alphabet))
                end
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
#     total_kmer_count = sum(counts)
#     for node in oriented_path[2:end-1]
#         likelihood *= counts[node.index] / total_kmer_count
#     end
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
function assess_path(path,
    kmers,
    counts,
    initial_orientation,
    outgoing_edge_probabilities,
    incoming_edge_probabilities)
    
    orientations = assess_path_orientations(path, kmers, initial_orientation)
    if orientations == nothing
        path_likelihood = 0.0
        oriented_path = OrientedKmer[]
    else
        oriented_path = orient_path(path, orientations)
        path_likelihood = assess_path_likelihood(oriented_path, kmers, counts, outgoing_edge_probabilities, incoming_edge_probabilities)
    end
    return (oriented_path = oriented_path, path_likelihood = path_likelihood)    
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

function count_kmers(::Type{KMER_TYPE}, record::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    return count_kmers(KMER_TYPE, FASTX.sequence(record))    
end

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
    k = length(first(kmers))
    initial_kmer_sequence = orient_oriented_kmer(kmers, oriented_path[1])
    sequence = BioSequences.LongDNASeq(initial_kmer_sequence)
    for i in 2:length(oriented_path)
        this_kmer_sequence = orient_oriented_kmer(kmers, oriented_path[i])
        @assert sequence[end-k+2:end] == BioSequences.LongDNASeq(this_kmer_sequence)[1:end-1]     
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

    if current_kmer_index == previous_kmer_index
        # could be a self loop, check this first
        if LightGraphs.has_edge(graph, LightGraphs.Edge(previous_kmer_index, current_kmer_index))
            this_path = [previous_kmer_index, current_kmer_index]
            this_oriented_path, this_likelihood = 
                assess_path(this_path,
                    graph.kmers,
                    graph.counts,
                    previous_orientation,
                    outgoing_edge_probabilities,
                    incoming_edge_probabilities)
            if this_likelihood > path_likelihood
                path_likelihood = this_likelihood
                oriented_path = this_oriented_path
            end
        end
        
        # consider a deletion in observed sequene relative to the reference graph
        # see if this has any neighbors that circle back, and evaluate the likelihood for each
        outneighbors = find_outneighbors(previous_orientation, previous_kmer_index, outgoing_edge_probabilities, incoming_edge_probabilities)
        for outneighbor in outneighbors
            if LightGraphs.has_path(graph, outneighbor, current_kmer_index)
                # manually build path
                this_path = [previous_kmer_index, shortest_paths[outneighbor][current_kmer_index]...]
                this_oriented_path, this_likelihood = 
                    assess_path(this_path,
                        graph.kmers,
                        graph.counts,
                        previous_orientation,
                        outgoing_edge_probabilities,
                        incoming_edge_probabilities)
                if this_likelihood > path_likelihood
                    path_likelihood = this_likelihood
                    oriented_path = this_oriented_path
                end
            end
        end
    elseif LightGraphs.has_path(graph, previous_kmer_index, current_kmer_index)   
        this_path = shortest_paths[previous_kmer_index][current_kmer_index]
        this_oriented_path, this_likelihood = 
             assess_path(this_path,
                graph.kmers,
                graph.counts,
                previous_orientation,
                outgoing_edge_probabilities,
                incoming_edge_probabilities)
        if this_likelihood > path_likelihood
            path_likelihood = this_likelihood
            oriented_path = this_oriented_path
        end
    end
    
    return (oriented_path = oriented_path, path_likelihood = path_likelihood)
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

    # WHEN CODE STABILIZES, SWITCH BIGFLOAT STATE LIKELIHOODS BACK INTO LOG-ED FLOAT64!!!
    edit_distances = Array{Union{Int, Missing}}(missing, LightGraphs.nv(graph.graph), length(observed_path))
    arrival_paths = Array{Union{Vector{OrientedKmer}, Missing}}(missing, LightGraphs.nv(graph.graph), length(observed_path))
    state_likelihoods = zeros(BigFloat, LightGraphs.nv(graph.graph), length(observed_path))
    state_likelihoods[:, 1] .= graph.counts ./ sum(graph.counts)

    observed_kmer_sequence = orient_oriented_kmer(graph.kmers, first(observed_path))
    for (kmer_index, kmer) in enumerate(graph.kmers)
        
        alignment_result, orientation = assess_optimal_alignment(kmer, observed_kmer_sequence)

        for match in 1:alignment_result.total_matches
            state_likelihoods[kmer_index, 1] *= (1.0 - error_rate)
        end
        for edit in 1:alignment_result.total_edits
            state_likelihoods[kmer_index, 1] *= error_rate
        end
        if state_likelihoods[kmer_index, 1] > 0.0
            arrival_paths[kmer_index, 1] = [OrientedKmer(index = kmer_index, orientation = orientation)]
            edit_distances[kmer_index, 1] = alignment_result.total_edits
        end
    end
    state_likelihoods[:, 1] ./= sum(state_likelihoods[:, 1])
#     kmer_likelihoods[:, 1] .= log.(kmer_likelihoods[:, 1])
    
    return (state_likelihoods = state_likelihoods, arrival_paths = arrival_paths, edit_distances = edit_distances)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_emission(a, b, kmers)
    if a.orientation
        a_sequence = kmers[a.index]
    else
        a_sequence = BioSequences.reverse_complement(kmers[a.index])
    end
    if !ismissing(b.orientation)
        if b.orientation
            b_sequence = kmers[b.index]
        else
            b_sequence = BioSequences.reverse_complement(kmers[b.index])
        end
        is_match = a_sequence[end] == b_sequence[end]
        orientation = b.orientation
    else
        forward_match = a_sequence[end] == kmers[b.index][end]
        reverse_match = a_sequence[end] == BioSequences.reverse_complement(kmers[b.index])[end]
        is_match = forward_match || reverse_match
        if forward_match && !reverse_match
            orientation = true
        elseif !forward_match && reverse_match
            orientation = false
        else
            orientation = missing
        end
    end
    return is_match, orientation
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_optimal_path(
    graph,
    oriented_path,
    path_likelihood,
    observed_kmer,
    current_kmer_likelihood,
    prior_state_likelihood,
    previous_kmer_index,
    current_observation_index,
    prior_edit_distance,
    error_rate)
    
    path_likelihood *= current_kmer_likelihood
    path_likelihood *= prior_state_likelihood
    emission_match, orientation = assess_emission(observed_kmer, last(oriented_path), graph.kmers)
    # assert final orientation if it's missing
    if ismissing(last(oriented_path).orientation) && !ismissing(orientation)
        oriented_path[end] = OrientedKmer(index = last(oriented_path).index, orientation = orientation)
    end
    edit_distance = !emission_match + abs(length(oriented_path) - 2)
    if edit_distance == 0
        path_likelihood *= (1 - error_rate)
    else
        for i in 1:edit_distance
            path_likelihood *= error_rate
        end
    end
    edit_distance += prior_edit_distance
    return oriented_path, path_likelihood, edit_distance
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Given a graph with known edge probabilities, determine the median likelihood of non-zero 2-step transitions.
This value is then used to set the likelihood of 0-step transitions, which cannot be measured directly.
This enables an approximately balanced likelihood of insertions and deletions

It may be better to fit a probability distribution and allow sampling from that probability distribution
rather than taking the median non-zero value and fixing all 0-step transitions to that likelihood, but
that would require simulations to confirm

```jldoctest
julia> 1 + 1
2
```
"""
function determine_insertion_path_likelihoods(graph, shortest_paths, outgoing_edge_probabilities, incoming_edge_probabilities)
    # could switch this to use online stats to be more efficient
    likelihoods = Float64[]
    for i in 1:length(shortest_paths)
        for j in 1:length(shortest_paths[i])
            if length(shortest_paths[i][j]) == 3
                path = shortest_paths[i][j]
                for orientation in (true, false)
                    this_oriented_path, this_likelihood =
                        assess_path(path,
                            graph.kmers,
                            graph.counts,
                            orientation,
                            outgoing_edge_probabilities,
                            incoming_edge_probabilities)
                    push!(likelihoods, this_likelihood)
                end
            end
        end
    end
    if all(x -> x == 0.0, likelihoods)
        return 0.0
    else
        return Statistics.median(filter!(x -> x > 0.0, likelihoods))
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
function viterbi_maximum_likelihood_path(graph, observation, error_rate; debug = false)

    observed_path = sequence_to_oriented_path(observation, graph.kmers)
    
    outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
    
    shortest_paths = LightGraphs.enumerate_paths(LightGraphs.floyd_warshall_shortest_paths(graph.graph))
    
    insertion_path_likelihoods = determine_insertion_path_likelihoods(graph, shortest_paths, outgoing_edge_probabilities, incoming_edge_probabilities)
#    @show insertion_path_likelihoods
    
    state_likelihoods, arrival_paths, edit_distances = initialize_viterbi(graph, observed_path, error_rate)

    if debug
        my_show(arrival_paths, graph.kmers, title = "Arrival Paths")
        my_show(state_likelihoods, graph.kmers, title="State Likelihoods")
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
            
            current_kmer_likelihood = BigFloat(graph.counts[current_kmer_index] / sum(graph.counts))
            
            best_state_likelihood = BigFloat(0.0)
            best_arrival_path = Vector{Eisenia.OrientedKmer}()
            best_edit_distance = 0

            for (previous_kmer_index, previous_kmer) in enumerate(graph.kmers)
                
                if debug
                    println("\tprevious_kmer_index = $(previous_kmer_index)")
                    println("\tprevious_kmer = $(previous_kmer)")
                end
                
                if state_likelihoods[previous_kmer_index, current_observation_index - 1] == 0.0
                    continue
                    # unreachable condition
                else
                    previous_arrival_path = arrival_paths[previous_kmer_index, current_observation_index - 1]
                    previous_orientation = last(arrival_paths[previous_kmer_index, current_observation_index - 1]).orientation
                end

                oriented_path, path_likelihood =
                    Eisenia.find_optimal_path(observed_kmer,
                        previous_kmer_index,
                        previous_orientation,
                        current_kmer_index,
                        graph, 
                        shortest_paths,
                        outgoing_edge_probabilities, 
                        incoming_edge_probabilities,
                        error_rate)
                if !isempty(oriented_path)
                    prior_state_likelihood = state_likelihoods[previous_kmer_index, current_observation_index - 1]
                    prior_edit_distance = edit_distances[previous_kmer_index, current_observation_index-1]
                    oriented_path, path_likelihood, edit_distance = 
                        assess_optimal_path(
                            graph,
                            oriented_path,
                            path_likelihood,
                            observed_kmer,
                            current_kmer_likelihood,
                            prior_state_likelihood,
                            previous_kmer_index,
                            current_observation_index,
                            prior_edit_distance,
                            error_rate)
                    if path_likelihood > best_state_likelihood
                        best_state_likelihood = path_likelihood
                        best_arrival_path = oriented_path
                        best_edit_distance = edit_distance
                    end
                end
            end
            # consider insertion
            # not totally sure if this is correct
            insertion_arrival_path = arrival_paths[current_kmer_index, current_observation_index - 1]
            if !ismissing(insertion_arrival_path) && !isempty(insertion_arrival_path)
                prior_state_likelihood = state_likelihoods[current_kmer_index, current_observation_index - 1]
                prior_edit_distance = edit_distances[current_kmer_index, current_observation_index-1]
                oriented_path = [last(insertion_arrival_path)]
#                 path_likelihood = prior_state_likelihood * rand(insertion_path_likelihoods) * current_kmer_likelihood * error_rate
                path_likelihood = prior_state_likelihood * insertion_path_likelihoods * current_kmer_likelihood * error_rate
                edit_distance = 1 + prior_edit_distance

                if path_likelihood > best_state_likelihood
                    best_state_likelihood = path_likelihood
                    best_arrival_path = oriented_path
                    best_edit_distance = edit_distance
                end
            end
            
            state_likelihoods[current_kmer_index, current_observation_index] = best_state_likelihood
            arrival_paths[current_kmer_index, current_observation_index] = best_arrival_path
            edit_distances[current_kmer_index, current_observation_index] = best_edit_distance
        end
    end
    if debug
        my_show(arrival_paths, graph.kmers, title = "Arrival Paths")
        my_show(state_likelihoods, graph.kmers, title="State Likelihoods")
        my_show(edit_distances, graph.kmers, title="Edit Distances")
    end
    return backtrack_optimal_path(state_likelihoods, arrival_paths, edit_distances)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function plot_kmer_frequency_spectra(counts; log_scaled = true, kwargs...)
    kmer_counts_hist = StatsBase.countmap(c for c in counts)
    xs = collect(keys(kmer_counts_hist))
    ys = collect(values(kmer_counts_hist))
    if log_scaled
        xs = log.(xs)
        ys = log.(ys)
    end
    
    StatsPlots.plot(
        xs,
        ys,
        xlims = (0, maximum(xs) + max(1, ceil(0.1 * maximum(xs)))),
        ylims = (0, maximum(ys) + max(1, ceil(0.1 * maximum(ys)))),
        seriestype = :scatter,
        legend = false,
        xlabel = log_scaled ? "log(observed frequency)" : "observed frequency",
        ylabel = log_scaled ? "log(# of kmers)" : "observed frequency",
        ;kwargs...
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

    n = length(graph.kmers)
    p = GraphRecipes.graphplot(
        graph.graph,
        markersize = 1/log2(n),
        size = (100 * log(n), 66 * log(n)),
        node_weights = graph.counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function path_to_sequence(kmer_graph, path)
    sequence = BioSequences.LongDNASeq(oriented_kmer_to_sequence(kmer_graph, first(path)))
    for oriented_kmer in path[2:end]
        nucleotide = last(oriented_kmer_to_sequence(kmer_graph, oriented_kmer))
        push!(sequence, nucleotide)
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
function oriented_kmer_to_sequence(kmer_graph, oriented_kmer)
    kmer_sequence = kmer_graph.kmers[oriented_kmer.index]
    if !oriented_kmer.orientation
        kmer_sequence = BioSequences.reverse_complement(kmer_sequence)
    end
    return kmer_sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function maximum_likelihood_walk(graph, connected_component)
    max_count = maximum(graph.counts[connected_component])
    max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
    initial_node_index = rand(max_count_indices)
    initial_node = connected_component[initial_node_index]
    outgoing_edge_probabilities, incoming_edge_probabilities = Eisenia.determine_edge_probabilities(graph)
    forward_walk = maximum_likelihood_walk(graph, [Eisenia.OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
    reverse_walk = maximum_likelihood_walk(graph, [Eisenia.OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
    reversed_reverse_walk = reverse!(
        [
            Eisenia.OrientedKmer(index = oriented_kmer.index, orientation = oriented_kmer.orientation)
            for oriented_kmer in reverse_walk[2:end]
        ]
        )
    full_path = [reversed_reverse_walk..., forward_walk...]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function maximum_likelihood_walk(graph, path::Vector{Eisenia.OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
    done = false
    while !done
        maximum_path_likelihood = 0.0
        maximum_likelihood_path = Vector{Eisenia.OrientedKmer}()
        for neighbor in LightGraphs.neighbors(graph.graph, last(path).index)
            this_path = [last(path).index, neighbor]
            this_oriented_path, this_path_likelihood = 
                Eisenia.assess_path(this_path,
                    graph.kmers,
                    graph.counts,
                    last(path).orientation,
                    outgoing_edge_probabilities,
                    incoming_edge_probabilities)
            if this_path_likelihood > maximum_path_likelihood
                maximum_path_likelihood = this_path_likelihood
                maximum_likelihood_path = this_oriented_path
            end
        end
        if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
            done = true
        else
            append!(path, maximum_likelihood_path[2:end])
        end
    end
    return path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function my_plot(graph::Eisenia.KmerGraph)
    graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))
    filename = "/assets/images/$(graph_hash).svg"
    p = Eisenia.plot_graph(graph)
    Plots.savefig(p, dirname(pwd()) * filename)
    display(p)
    display("text/markdown", "![]($filename)")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_observations(graph::Eisenia.KmerGraph{KMER_TYPE}, observations, error_rate; verbose = isinteractive()) where {KMER_TYPE}
    k = last(KMER_TYPE.parameters)
    total_edits_accepted = 0
    total_bases_evaluated = 0
    reads_processed = 0
    maximum_likelihood_observations = Vector{BioSequences.LongDNASeq}(undef, length(observations))
    for (observation_index, observation) in enumerate(observations)
        if length(observation) >= k
            optimal_path, edit_distance, relative_likelihood = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
            maximum_likelihood_observation = Eisenia.oriented_path_to_sequence(optimal_path, graph.kmers)
            maximum_likelihood_observations[observation_index] = maximum_likelihood_observation
            reads_processed += 1
            total_bases_evaluated += length(observation)
            total_edits_accepted += edit_distance
        else
            maximum_likelihood_observations[observation_index] = observation
        end
    end
    inferred_error_rate = round(total_edits_accepted / total_bases_evaluated, digits = 3)
    if verbose
        display("reads_processed = $(reads_processed)")
        display("total_edits_accepted = $(total_edits_accepted)")
        display("inferred_error_rate = $(inferred_error_rate)")
    end
    if total_edits_accepted == 0
        has_converged = true
    else
        has_converged = false
    end
    return maximum_likelihood_observations, has_converged
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function iterate_until_convergence(ks, observations, error_rate; verbose = isinteractive())
    for k in ks
        graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, observations)
        if verbose
            display("k = $k")
            my_plot(graph)
        end
        observations, has_converged = assess_observations(graph, observations, error_rate; verbose = verbose)
    end
    graph = Eisenia.KmerGraph(BioSequences.DNAMer{last(ks)}, observations)
    if verbose
        display("final graph")
        my_plot(graph)
    end
    return graph, observations
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function clip_low_coverage_tips(graph, observations)
    connected_components = LightGraphs.connected_components(graph.graph)
    vertices_to_keep = Int[]
    for connected_component in connected_components
        
        component_coverage = graph.counts[connected_component]
        median = Statistics.median(component_coverage)
        standard_deviation = Statistics.std(component_coverage)
        
        for vertex in connected_component
            keep = true
            if LightGraphs.degree(graph.graph, vertex) == 1
                this_coverage = graph.counts[vertex]
                is_low_coverage = (graph.counts[vertex] == 1) || 
                                    (median-this_coverage) > (3*standard_deviation)
                if is_low_coverage
                    keep = false
                end
            end
            if keep
                push!(vertices_to_keep, vertex)
            end
        end
    end
    
    KmerType = first(typeof(graph).parameters)
    pruned_graph = Eisenia.KmerGraph(KmerType, observations, graph.kmers[vertices_to_keep], graph.counts[vertices_to_keep])
    
    return pruned_graph
end


function reverse_oriented_path(oriented_path)
    reversed_path = copy(oriented_path)
    for (index, state) in enumerate(oriented_path)
        reversed_path[index] = Eisenia.OrientedKmer(index = state.index, orientation = !state.orientation)
    end
    return reverse!(reversed_path)
end


function take_a_walk(graph, connected_component)
    max_count = maximum(graph.counts[connected_component])
    max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
    initial_node_index = rand(max_count_indices)
    initial_node = connected_component[initial_node_index]
    outgoing_edge_probabilities, incoming_edge_probabilities = Eisenia.determine_edge_probabilities(graph)
    
    # walk forwards from the initial starting node
    forward_walk = take_a_walk(graph, [Eisenia.OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
    # walk backwards from the initial starting node
    reverse_walk = take_a_walk(graph, [Eisenia.OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
    # we need to reverse everything to re-orient against the forward walk
    reverse_walk = reverse_oriented_path(reverse_walk)
    
    # also need to drop the last node, which is equivalent to the first node of the 
    @assert last(reverse_walk) == first(forward_walk)
    full_path = [reverse_walk[1:end-1]..., forward_walk...]
#     @show full_path
end

function take_a_walk(graph, path::Vector{Eisenia.OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
    done = false
    while !done
        maximum_path_likelihood = 0.0
        maximum_likelihood_path = Vector{Eisenia.OrientedKmer}()
        for neighbor in LightGraphs.neighbors(graph.graph, last(path).index)
            this_path = [last(path).index, neighbor]
            this_oriented_path, this_path_likelihood = 
                Eisenia.assess_path(this_path,
                    graph.kmers,
                    graph.counts,
                    last(path).orientation,
                    outgoing_edge_probabilities,
                    incoming_edge_probabilities)
            if this_path_likelihood > maximum_path_likelihood
                maximum_path_likelihood = this_path_likelihood
                maximum_likelihood_path = this_oriented_path
            end
        end
        if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
            done = true
        else
            append!(path, maximum_likelihood_path[2:end])
        end
    end
    return path
end

function is_equivalent(a, b)
    a == b || a == BioSequences.reverse_complement(b)
end

end # module
