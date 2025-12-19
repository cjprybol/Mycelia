# Graph Type Conversion Helpers
#
# Conversions between fixed-length (k-mer/qualmer) graphs and variable-length
# OLC graphs, plus utilities for dropping quality scores.

"""
    convert_fixed_to_variable(graph::MetaGraphsNext.MetaGraph)

Convert a fixed-length k-mer or qualmer graph into a variable-length graph
using full BioSequence vertices and overlap lengths derived from k.
Preserves evidence (quality-aware graphs keep quality evidence).
"""
function convert_fixed_to_variable(graph::MetaGraphsNext.MetaGraph)
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return MetaGraphsNext.MetaGraph(Graphs.DiGraph(); label_type=String)
    end

    first_label = first(labels)
    sequence_type = _sequence_type_from_kmer(first_label)
    k = Kmers.ksize(typeof(first_label))

    first_vertex_data = graph[first_label]
    is_quality = first_vertex_data isa QualmerVertexData

    graph_type = typeof(graph.graph)
    new_graph = MetaGraphsNext.MetaGraph(
        graph_type();
        label_type=sequence_type,
        vertex_data_type=is_quality ? QualityBioSequenceVertexData{sequence_type} : BioSequenceVertexData{sequence_type},
        edge_data_type=is_quality ? QualityBioSequenceEdgeData : BioSequenceEdgeData
    )

    # Add vertices
    for label in labels
        sequence = sequence_type(string(label))
        if !haskey(new_graph, sequence)
            vertex_data = is_quality ? QualityBioSequenceVertexData(sequence) : BioSequenceVertexData(sequence)
            new_graph[sequence] = vertex_data
        end
        _copy_vertex_evidence!(new_graph[sequence], graph[label]; drop_quality=false)
    end

    # Add edges with overlap length = k - 1
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        new_src = sequence_type(string(src))
        new_dst = sequence_type(string(dst))
        edge_data = graph[src, dst]
        new_edge = is_quality ? QualityBioSequenceEdgeData(k - 1) : BioSequenceEdgeData(k - 1)
        _copy_edge_evidence!(new_edge, edge_data; drop_quality=false)
        new_graph[new_src, new_dst] = new_edge
    end

    return new_graph
end

"""
    convert_variable_to_fixed(graph::MetaGraphsNext.MetaGraph, ::Type{KmerType})

Convert a variable-length graph (BioSequence or quality-aware BioSequence)
into a fixed-length k-mer graph of the specified `KmerType`.
All sequences must match the requested k.
Quality-aware evidence is preserved unless `drop_quality=true`.
"""
function convert_variable_to_fixed(graph::MetaGraphsNext.MetaGraph, ::Type{KmerType}; drop_quality::Bool=false) where {KmerType<:Kmers.Kmer}
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return MetaGraphsNext.MetaGraph(Graphs.DiGraph(); label_type=KmerType)
    end

    k = Kmers.ksize(KmerType)

    # Validate lengths
    for seq in labels
        if length(seq) != k
            error("Cannot convert variable-length graph to $(KmerType): sequence length $(length(seq)) != k=$(k)")
        end
    end

    graph_type = typeof(graph.graph)
    source_is_quality = graph[first(labels)] isa QualityBioSequenceVertexData
    effective_drop_quality = drop_quality || source_is_quality

    new_graph = MetaGraphsNext.MetaGraph(
        graph_type();
        label_type=KmerType,
        vertex_data_type=KmerVertexData{KmerType},
        edge_data_type=KmerEdgeData
    )

    for seq in labels
        kmer_label = KmerType(string(seq))
        if !haskey(new_graph, kmer_label)
            new_graph[kmer_label] = KmerVertexData(kmer_label)
        end
        _copy_vertex_evidence!(new_graph[kmer_label], graph[seq]; drop_quality=effective_drop_quality)
    end

    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        new_src = KmerType(string(src))
        new_dst = KmerType(string(dst))
        edge_data = graph[src, dst]
        new_edge = KmerEdgeData()
        _copy_edge_evidence!(new_edge, edge_data; drop_quality=effective_drop_quality)
        new_graph[new_src, new_dst] = new_edge
    end

    return new_graph
end

"""
    drop_quality_scores(graph::MetaGraphsNext.MetaGraph)

Drop quality scores from qualmer or quality-aware BioSequence graphs while
preserving topology and positional evidence.
"""
function drop_quality_scores(graph::MetaGraphsNext.MetaGraph)
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return graph
    end

    first_vertex_data = graph[first(labels)]

    if first_vertex_data isa QualmerVertexData
        return _drop_quality_kmer_graph(graph)
    elseif first_vertex_data isa QualityBioSequenceVertexData
        return _drop_quality_biosequence_graph(graph)
    else
        return graph  # Already quality-unaware
    end
end

function _drop_quality_kmer_graph(graph::MetaGraphsNext.MetaGraph)
    labels = collect(MetaGraphsNext.labels(graph))
    first_label = first(labels)
    kmer_type = typeof(first_label)
    new_graph = MetaGraphsNext.MetaGraph(
        typeof(graph.graph)();
        label_type=kmer_type,
        vertex_data_type=KmerVertexData{kmer_type},
        edge_data_type=KmerEdgeData
    )

    for label in labels
        if !haskey(new_graph, label)
            new_graph[label] = KmerVertexData(label)
        end
        _copy_vertex_evidence!(new_graph[label], graph[label]; drop_quality=true)
    end

    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        new_edge = KmerEdgeData()
        _copy_edge_evidence!(new_edge, graph[src, dst]; drop_quality=true)
        new_graph[src, dst] = new_edge
    end

    return new_graph
end

function _drop_quality_biosequence_graph(graph::MetaGraphsNext.MetaGraph)
    labels = collect(MetaGraphsNext.labels(graph))
    first_label = first(labels)
    seq_type = typeof(first_label)
    new_graph = MetaGraphsNext.MetaGraph(
        typeof(graph.graph)();
        label_type=seq_type,
        vertex_data_type=BioSequenceVertexData{seq_type},
        edge_data_type=BioSequenceEdgeData
    )

    for label in labels
        if !haskey(new_graph, label)
            new_graph[label] = BioSequenceVertexData(label)
        end
        _copy_vertex_evidence!(new_graph[label], graph[label]; drop_quality=true)
    end

    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        edge_data = graph[src, dst]
        if hasfield(typeof(edge_data), :overlap_length)
            new_edge = BioSequenceEdgeData(edge_data.overlap_length)
        else
            new_edge = BioSequenceEdgeData(0)
        end
        _copy_edge_evidence!(new_edge, edge_data; drop_quality=true)
        new_graph[src, dst] = new_edge
    end

    return new_graph
end

# ============================================================================
# Evidence copying utilities (position-preserving, optional quality drop)
# ============================================================================

function _copy_vertex_evidence!(target, source; drop_quality::Bool)
    for (dataset_id, dataset_evidence) in source.evidence
        for (obs_id, evidence_set) in dataset_evidence
            for entry in evidence_set
                converted = drop_quality ? _strip_quality(entry) : entry
                add_evidence!(target, dataset_id, obs_id, converted)
            end
        end
    end
end

function _copy_edge_evidence!(target, source; drop_quality::Bool)
    for (dataset_id, dataset_evidence) in source.evidence
        for (obs_id, evidence_set) in dataset_evidence
            for entry in evidence_set
                converted = drop_quality ? _strip_quality(entry) : entry
                add_evidence!(target, dataset_id, obs_id, converted)
            end
        end
    end
end

function _strip_quality(entry::QualityEvidenceEntry)
    return EvidenceEntry(entry.position, entry.strand)
end

function _strip_quality(entry::EvidenceEntry)
    return entry
end

function _strip_quality(entry::EdgeQualityEvidenceEntry)
    return EdgeEvidenceEntry(entry.from_position, entry.to_position, entry.strand)
end

function _strip_quality(entry::EdgeEvidenceEntry)
    return entry
end

"""
Determine BioSequence type for a k-mer label.
"""
function _sequence_type_from_kmer(kmer)
    kmer_type = typeof(kmer)

    if kmer_type <: Kmers.DNAKmer
        return BioSequences.LongDNA{4}
    elseif kmer_type <: Kmers.RNAKmer
        return BioSequences.LongRNA{4}
    elseif kmer_type <: Kmers.AAKmer
        return BioSequences.LongAA
    else
        error("Unsupported k-mer type for conversion: $(kmer_type)")
    end
end
