"""
K-mer graph utilities extracted from `sequence-graphs-next.jl`.

Provides type-stable MetaGraphsNext k-mer graph construction along with
helper functions for strand-aware coverage tracking and GFA export.
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers

"""
Metadata for k-mer graph vertices with strand-aware coverage.
"""
struct KmerVertexData{KmerT}
    coverage::Vector{Tuple{Int, Int, StrandOrientation}}
    canonical_kmer::KmerT

    KmerVertexData(canonical_kmer::KmerT) where {KmerT} =
        new{KmerT}(Vector{Tuple{Int, Int, StrandOrientation}}(), canonical_kmer)
end

"""
Metadata for k-mer graph edges that track strand orientation.
"""
struct KmerEdgeData
    coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}
    weight::Float64
    src_strand::StrandOrientation
    dst_strand::StrandOrientation

    function KmerEdgeData(
        coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}},
        src_strand::StrandOrientation,
        dst_strand::StrandOrientation,
    )
        new(coverage, Float64(length(coverage)), src_strand, dst_strand)
    end

    KmerEdgeData(src_strand::StrandOrientation, dst_strand::StrandOrientation) =
        new(
            Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}(),
            0.0,
            src_strand,
            dst_strand,
        )
end

"""
Construct a strand-aware MetaGraphsNext k-mer graph.
"""
function build_kmer_graph_next(
    kmer_type,
    observations::AbstractVector{<:Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}};
    graph_mode::GraphMode = DoubleStrand,
)
    kmer_counts = graph_mode == DoubleStrand ?
        Mycelia.count_canonical_kmers(kmer_type, observations) :
        Mycelia.count_kmers(kmer_type, observations)

    graph_kmers = collect(keys(kmer_counts))
    if isempty(graph_kmers)
        @warn "No k-mers found in observations"
        return _create_empty_kmer_graph(kmer_type)
    end

    actual_kmer_type = eltype(graph_kmers)
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = actual_kmer_type,
        vertex_data_type = KmerVertexData{actual_kmer_type},
        edge_data_type = KmerEdgeData,
        weight_function = edge_data -> edge_data.weight,
        default_weight = 0.0,
    )

    for kmer in graph_kmers
        graph[kmer] = KmerVertexData(kmer)
    end

    for (obs_idx, observation) in enumerate(observations)
        _add_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
    end

    return graph
end

function _add_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
    observed_sequence = FASTX.sequence(observation)
    k = length(first(graph_kmers))

    if length(observed_sequence) < k
        observation_id = FASTX.identifier(observation)
        @warn "Skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        return
    end

    observed_path = _sequence_to_canonical_path(graph_kmers, observed_sequence, graph_mode)
    isempty(observed_path) && return

    first_canonical_kmer, first_strand = observed_path[1]
    _add_vertex_coverage!(graph, first_canonical_kmer, obs_idx, 1, first_strand)

    for i in 2:length(observed_path)
        curr_canonical_kmer, curr_strand = observed_path[i]
        prev_canonical_kmer, prev_strand = observed_path[i - 1]

        _add_vertex_coverage!(graph, curr_canonical_kmer, obs_idx, i, curr_strand)

        _add_strand_aware_edge!(
            graph,
            prev_canonical_kmer,
            curr_canonical_kmer,
            prev_strand,
            curr_strand,
            (obs_idx, i - 1, prev_strand),
            (obs_idx, i, curr_strand),
        )
    end
end

function _add_vertex_coverage!(graph, kmer, obs_idx, position, strand_orientation)
    vertex_data = graph[kmer]
    push!(vertex_data.coverage, (obs_idx, position, strand_orientation))
end

function _add_strand_aware_edge!(graph, src_kmer, dst_kmer, src_strand, dst_strand, src_coverage, dst_coverage)
    edge_exists = false
    existing_edge_data = nothing

    if haskey(graph, src_kmer, dst_kmer)
        existing_edge_data = graph[src_kmer, dst_kmer]
        if existing_edge_data.src_strand == src_strand && existing_edge_data.dst_strand == dst_strand
            edge_exists = true
        end
    end

    if !edge_exists
        graph[src_kmer, dst_kmer] = KmerEdgeData(src_strand, dst_strand)
        existing_edge_data = graph[src_kmer, dst_kmer]
    end

    push!(existing_edge_data.coverage, (src_coverage, dst_coverage))
    graph[src_kmer, dst_kmer] = KmerEdgeData(
        existing_edge_data.coverage,
        existing_edge_data.src_strand,
        existing_edge_data.dst_strand,
    )
end

function _create_empty_kmer_graph(kmer_type)
    MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = kmer_type,
        vertex_data_type = KmerVertexData{kmer_type},
        edge_data_type = KmerEdgeData,
        weight_function = edge_data -> edge_data.weight,
        default_weight = 0.0,
    )
end

function _sequence_to_canonical_path(graph_kmers, sequence, graph_mode)
    k = length(first(graph_kmers))
    path = Vector{Tuple{eltype(graph_kmers), StrandOrientation}}()
    graph_kmer_set = Set(graph_kmers)

    for i in 1:(length(sequence) - k + 1)
        subseq = sequence[i:i + k - 1]
        observed_kmer = typeof(first(graph_kmers))(subseq)

        if graph_mode == SingleStrand
            if observed_kmer in graph_kmer_set
                push!(path, (observed_kmer, Forward))
            else
                @warn "K-mer $observed_kmer not found in graph at position $i (SingleStrand mode)"
            end
        else
            if typeof(observed_kmer) <: Kmers.Kmer{<:BioSequences.NucleicAcidAlphabet}
                canonical_kmer = BioSequences.canonical(observed_kmer)
                strand_orientation = BioSequences.iscanonical(observed_kmer) ? Forward : Reverse
            else
                if observed_kmer in graph_kmer_set
                    push!(path, (observed_kmer, Forward))
                else
                    @warn "K-mer $observed_kmer not found in graph at position $i (amino acid in DoubleStrand mode)"
                end
                continue
            end

            if canonical_kmer in graph_kmer_set
                push!(path, (canonical_kmer, strand_orientation))
            else
                @warn "Canonical k-mer $canonical_kmer not found in graph at position $i (DoubleStrand mode)"
            end
        end
    end

    return path
end

function _is_valid_transition(src_kmer, dst_kmer, src_strand, dst_strand, k)
    src_seq = src_strand == Forward ? src_kmer : BioSequences.reverse_complement(src_kmer)
    dst_seq = dst_strand == Forward ? dst_kmer : BioSequences.reverse_complement(dst_kmer)

    src_suffix = src_seq[2:end]
    dst_prefix = dst_seq[1:end-1]

    return src_suffix == dst_prefix
end

"""
Export a MetaGraphsNext graph to GFA format.
"""
function write_gfa_next(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString)
    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0\tMY:Z:Mycelia-Next")

        vertex_id_map = Dict()
        for (i, label) in enumerate(MetaGraphsNext.labels(graph))
            vertex_id_map[label] = i
            vertex_data = graph[label]
            depth = length(vertex_data.coverage)

            sequence = if hasfield(typeof(vertex_data), :canonical_kmer)
                string(vertex_data.canonical_kmer)
            elseif hasfield(typeof(vertex_data), :canonical_qualmer)
                string(vertex_data.canonical_qualmer.kmer)
            elseif hasfield(typeof(vertex_data), :sequence)
                string(vertex_data.sequence)
            else
                error("Unknown vertex data type: $(typeof(vertex_data))")
            end

            println(io, "S\t$i\t$(sequence)\tDP:f:$depth")
        end

        if !isempty(MetaGraphsNext.labels(graph))
            first_label = first(MetaGraphsNext.labels(graph))
            first_vertex_data = graph[first_label]
            k = if hasfield(typeof(first_vertex_data), :canonical_kmer)
                length(first_vertex_data.canonical_kmer)
            elseif hasfield(typeof(first_vertex_data), :canonical_qualmer)
                length(first_vertex_data.canonical_qualmer.kmer)
            elseif hasfield(typeof(first_vertex_data), :sequence)
                if !isempty(first_vertex_data.constituent_kmers)
                    length(first_vertex_data.constituent_kmers[1])
                elseif hasfield(typeof(graph), :graph_data) && !isnothing(graph.graph_data) && haskey(graph.graph_data, :k)
                    graph.graph_data[:k]
                else
                    length(first_vertex_data.sequence)
                end
            else
                error("Unknown vertex data type: $(typeof(first_vertex_data))")
            end
            overlap = k - 1
        else
            overlap = 0
        end

        for edge_labels in MetaGraphsNext.edge_labels(graph)
            if length(edge_labels) == 2
                src_label, dst_label = edge_labels
                edge_data = graph[src_label, dst_label]
                src_id = vertex_id_map[src_label]
                dst_id = vertex_id_map[dst_label]

                src_orientation = edge_data.src_strand == Forward ? '+' : '-'
                dst_orientation = edge_data.dst_strand == Forward ? '+' : '-'

                println(io, "L\t$src_id\t$src_orientation\t$dst_id\t$dst_orientation\t$(overlap)M")
            end
        end
    end

    return outfile
end

