"""
Qualmer graph utilities extracted from `sequence-graphs-next.jl`.

Provides quality-aware graph builders that operate directly on FASTQ
records along with helper functions for integrating qualmer paths.
"""

import MetaGraphsNext
import BioSequences
import FASTX
import Statistics

function build_qualmer_graph_next(
    kmer_type,
    fastq_records::Vector{FASTX.FASTQ.Record};
    graph_mode::GraphMode = DoubleStrand,
)
    kmer_counts = graph_mode == DoubleStrand ?
        Mycelia.count_canonical_qualmers(kmer_type, fastq_records) :
        Mycelia.count_qualmers(kmer_type, fastq_records)

    graph_kmers = collect(keys(kmer_counts))
    if isempty(graph_kmers)
        @warn "No k-mers found in FASTQ observations"
        return _create_empty_qualmer_graph(kmer_type)
    end

    actual_kmer_type = eltype(graph_kmers)
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = actual_kmer_type,
        vertex_data_type = QualmerVertexData{actual_kmer_type},
        edge_data_type = QualmerEdgeData,
        weight_function = edge_data -> edge_data.weight,
        default_weight = 0.0,
    )

    for kmer in graph_kmers
        graph[kmer] = QualmerVertexData(kmer)
    end

    for (obs_idx, observation) in enumerate(fastq_records)
        _add_qualmer_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
    end

    return graph
end

function build_qualmer_biosequence_graph_next(
    sequence_type::Type{<:BioSequences.BioSequence},
    fastq_records::Vector{FASTX.FASTQ.Record};
    graph_mode::GraphMode = DoubleStrand,
)
    if isempty(fastq_records)
        throw(ArgumentError("Cannot build graph from empty FASTQ records"))
    end

    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = sequence_type,
        vertex_data_type = QualityBioSequenceVertexData{sequence_type},
        edge_data_type = QualityBioSequenceEdgeData,
        weight_function = edge_data -> edge_data.weight,
        default_weight = 0.0,
    )

    for record in fastq_records
        sequence = FASTX.sequence(sequence_type, record)
        quality = FASTX.quality(record)

        canonical_sequence =
            if graph_mode == DoubleStrand && sequence isa Union{BioSequences.LongDNA, BioSequences.LongRNA}
                rc_sequence = BioSequences.reverse_complement(sequence)
                sequence <= rc_sequence ? sequence : rc_sequence
            else
                sequence
            end

        if !haskey(graph, canonical_sequence)
            graph[canonical_sequence] = QualityBioSequenceVertexData(canonical_sequence, [quality])
        else
            vertex_data = graph[canonical_sequence]
            push!(vertex_data.quality_scores, quality)
        end
    end

    return graph
end

function _create_empty_qualmer_graph(kmer_type)
    MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = kmer_type,
        vertex_data_type = QualmerVertexData{kmer_type},
        edge_data_type = QualmerEdgeData,
        weight_function = edge_data -> edge_data.weight,
        default_weight = 0.0,
    )
end

function _add_qualmer_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
    observed_sequence = FASTX.sequence(observation)
    quality_scores = FASTX.quality(observation)
    k = length(first(graph_kmers))

    if length(observed_sequence) < k
        observation_id = FASTX.identifier(observation)
        @warn "Skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        return
    end

    observed_path = _sequence_to_canonical_qualmer_path(graph_kmers, observed_sequence, quality_scores, graph_mode)
    isempty(observed_path) && return

    for i in 1:length(observed_path)
        canonical_kmer, strand, qual = observed_path[i]
        _add_qualmer_vertex_coverage!(graph, canonical_kmer, obs_idx, i, strand, qual)

        if i < length(observed_path)
            next_canonical_kmer, next_strand, next_qual = observed_path[i + 1]
            _add_qualmer_edge!(graph, canonical_kmer, next_canonical_kmer, strand, next_strand, qual, next_qual)
        end
    end
end

