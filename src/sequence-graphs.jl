"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build the legacy stranded DNA k-mer graph required by the canonical Viterbi
maximum-likelihood corrector.

This compatibility surface intentionally stays narrow: it supports the B0
DNA golden-master oracle and preserves the original `MetaGraphs.MetaDiGraph`
property contract (`:stranded_kmers`, `:observed_paths`, coverage props) that
`viterbi_maximum_likelihood_traversals` consumes. General graph/alphabet
interfaces belong to the follow-on B1+ beads.
"""
function build_stranded_kmer_graph(
        kmer_type::Type{KMER_TYPE},
        observations::AbstractVector{<:Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
)::MetaGraphs.MetaDiGraph where {KMER_TYPE}
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, observations)
    canonical_kmers = collect(keys(canonical_kmer_counts))
    stranded_kmers = sort!(vcat(
        canonical_kmers,
        [BioSequences.reverse_complement(kmer) for kmer in canonical_kmers]
    ))
    stranded_kmer_to_reverse_complement_map = [
        findfirst(candidate -> BioSequences.reverse_complement(candidate) == kmer, stranded_kmers)
        for kmer in stranded_kmers
    ]

    stranded_kmer_graph = MetaGraphs.MetaDiGraph(length(stranded_kmers))
    stranded_kmer_graph.gprops[:stranded_kmers] = stranded_kmers
    stranded_kmer_graph.gprops[:reverse_complement_map] = stranded_kmer_to_reverse_complement_map
    stranded_kmer_graph.gprops[:k] = length(first(stranded_kmers))
    stranded_kmer_graph.gprops[:K] = length(stranded_kmers)
    stranded_kmer_graph.gprops[:observation_color_map] = Vector{Int}()
    stranded_kmer_graph.gprops[:observation_ids] = Vector{String}()
    stranded_kmer_graph.gprops[:observed_paths] = Vector{Vector{Pair{Int, Bool}}}()

    for vertex in Graphs.vertices(stranded_kmer_graph)
        stranded_kmer_graph.vprops[vertex] = Dict(:coverage => Vector{Pair{Int, Pair{Int, Bool}}}())
    end

    for (observation_index, observation) in enumerate(observations)
        observation_id = String(FASTX.identifier(observation))
        observed_sequence = FASTX.sequence(observation)
        if length(observed_sequence) < stranded_kmer_graph.gprops[:k]
            @warn "skipping sequence shorter than k" observation_id length=length(observed_sequence)
            continue
        end

        observed_path = sequence_to_stranded_path(
            stranded_kmer_graph.gprops[:stranded_kmers],
            observed_sequence
        )
        isempty(observed_path) && continue

        ui, ui_orientation = observed_path[1]
        ui_coverage = observation_index => (1 => ui_orientation)
        push!(stranded_kmer_graph.vprops[ui][:coverage], ui_coverage)

        for path_index in 2:length(observed_path)
            vi, vi_orientation = observed_path[path_index]
            vi_coverage = observation_index => (path_index => vi_orientation)
            push!(stranded_kmer_graph.vprops[vi][:coverage], vi_coverage)
            edge_coverage = ui_coverage => vi_coverage
            if Graphs.has_edge(stranded_kmer_graph, ui, vi)
                push!(stranded_kmer_graph.eprops[Graphs.Edge(ui, vi)][:coverage], edge_coverage)
            else
                Graphs.add_edge!(stranded_kmer_graph, ui, vi, Dict(:coverage => [edge_coverage]))
            end
            ui = vi
            ui_coverage = vi_coverage
        end

        push!(stranded_kmer_graph.gprops[:observed_paths], observed_path)
        push!(stranded_kmer_graph.gprops[:observation_ids], observation_id)
        push!(stranded_kmer_graph.gprops[:observation_color_map], observation_index)
    end

    return stranded_kmer_graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a DNA sequence into a path through legacy stranded k-mer vertices.
"""
function sequence_to_stranded_path(
        stranded_kmers::AbstractVector{KMER_TYPE},
        sequence
)::Vector{Pair{Int, Bool}} where {KMER_TYPE}
    k = Kmers.ksize(KMER_TYPE)
    path = Vector{Pair{Int, Bool}}()
    for (kmer, _index) in Kmers.UnambiguousDNAMers{k}(BioSequences.LongDNA{4}(sequence))
        kmer_index = findfirst(stranded_kmer -> kmer == stranded_kmer, stranded_kmers)
        if kmer_index === nothing
            reverse_kmer = BioSequences.reverse_complement(kmer)
            kmer_index = findfirst(stranded_kmer -> reverse_kmer == stranded_kmer, stranded_kmers)
            orientation = false
        else
            orientation = true
        end
        kmer_index === nothing && continue
        push!(path, kmer_index => orientation)
    end
    return path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Reconstruct a DNA sequence from a legacy stranded k-mer path.
"""
function path_to_sequence(kmers::AbstractVector, path::AbstractVector)::BioSequences.LongDNA{4}
    sequence = BioSequences.LongDNA{4}(kmers[first(first(path))])
    for step in path[2:end]
        push!(sequence, kmers[first(step)][end])
    end
    return sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return the empirical transition probability for a legacy stranded graph edge.
"""
function edge_probability(stranded_kmer_graph::MetaGraphs.MetaDiGraph, edge::Graphs.Edge)::Float64
    neighbors = Graphs.outneighbors(stranded_kmer_graph, edge.src)
    neighbor_index = findfirst(neighbor -> neighbor == edge.dst, neighbors)
    neighbor_index === nothing && return 0.0

    edge_weights = [
        length(stranded_kmer_graph.eprops[Graphs.Edge(edge.src, neighbor)][:coverage])
        for neighbor in neighbors
    ]
    total_weight = sum(edge_weights)
    total_weight == 0 && return 0.0
    return edge_weights[neighbor_index] / total_weight
end
