"""
FASTQ Graph Implementation (Quality-Aware BioSequence Graphs)

This file implements variable-length quality-aware BioSequence graphs created by 
simplification/reduction of Qualmer graphs, following the 6-graph hierarchy specification.

Key features:
- Variable-length BioSequence vertices with per-base quality scores
- Quality retention throughout assembly process
- Convertible back to FASTQ records with full quality information
- Strand-aware edges with quality-weighted transitions
- Created by simplification of Qualmer graphs
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import DocStringExtensions
import Statistics

"""
Quality-aware BioSequence vertex data.
"""
struct QualityBioSequenceVertexData{SeqT<:BioSequences.BioSequence}
    sequence::SeqT                              # Variable-length BioSequence
    quality_scores::Vector{UInt8}               # Per-base PHRED quality scores
    coverage::Vector{Int}                       # Coverage from original qualmer observations
    constituent_qualmers::Vector{String}        # Original qualmers that were collapsed
    joint_probability::Float64                  # Joint probability of sequence correctness
    mean_quality::Float64                       # Mean quality across all bases
    length::Int                                 # Sequence length
    
    function QualityBioSequenceVertexData(sequence::SeqT, quality_scores::Vector{UInt8}, 
                                         coverage::Vector{Int}=Int[], 
                                         constituent_qualmers::Vector{String}=String[]) where {SeqT<:BioSequences.BioSequence}
        @assert length(sequence) == length(quality_scores) "Sequence and quality lengths must match"
        
        # Calculate joint probability from quality scores
        joint_prob = 1.0
        for q in quality_scores
            prob_correct = 1.0 - 10.0^(-q / 10.0)
            joint_prob *= prob_correct
        end
        
        mean_qual = isempty(quality_scores) ? 0.0 : Statistics.mean(quality_scores)
        
        new{SeqT}(sequence, quality_scores, coverage, constituent_qualmers, 
                 joint_prob, mean_qual, length(sequence))
    end
end

"""
Quality-aware BioSequence edge data.
"""
struct QualityBioSequenceEdgeData
    overlap_length::Int                         # Length of overlap between sequences
    src_strand::StrandOrientation               # Source sequence strand
    dst_strand::StrandOrientation               # Destination sequence strand
    weight::Float64                             # Edge weight from quality and coverage
    quality_weight::Float64                     # Quality-specific weight
    original_qualmer_count::Int                 # Number of qualmer transitions
    overlap_quality::Float64                    # Quality score of overlap region
    
    function QualityBioSequenceEdgeData(overlap_length::Int, src_strand::StrandOrientation, 
                                       dst_strand::StrandOrientation, weight::Float64=1.0, 
                                       quality_weight::Float64=1.0, original_qualmer_count::Int=1,
                                       overlap_quality::Float64=40.0)
        new(overlap_length, src_strand, dst_strand, weight, quality_weight, 
            original_qualmer_count, overlap_quality)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a quality-aware BioSequence graph directly from FASTQ records.

This creates a variable-length quality-aware BioSequence graph where vertices are 
the input sequences with their quality scores and edges represent quality-weighted 
relationships between sequences.

# Arguments
- `fastq_records`: Vector of FASTQ records with quality scores
- `graph_mode`: SingleStrand or DoubleStrand mode (default: DoubleStrand)
- `min_overlap`: Minimum overlap length for creating edges (default: 100)
- `min_quality`: Minimum mean quality to include sequence (default: 20)

# Returns
- MetaGraphsNext.MetaGraph with quality-aware BioSequence vertices

# Example
```julia
records = [FASTX.FASTQ.Record("seq1", "ATCGATCGATCG", "IIIIIIIIIIII")]
graph = build_quality_biosequence_graph(records)
```
"""
function build_quality_biosequence_graph(fastq_records::Vector{FASTX.FASTQ.Record}; 
                                        graph_mode::GraphMode=DoubleStrand,
                                        min_overlap::Int=100,
                                        min_quality::UInt8=UInt8(20))
    
    # Determine sequence type from first record
    if isempty(fastq_records)
        throw(ArgumentError("Cannot build graph from empty FASTQ records"))
    end
    
    first_seq = convert_sequence(FASTX.FASTQ.sequence(fastq_records[1]))
    biosequence_type = typeof(first_seq)
    
    # Create MetaGraphsNext graph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=biosequence_type,
        vertex_data_type=QualityBioSequenceVertexData{biosequence_type},
        edge_data_type=QualityBioSequenceEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices for each sequence (filtering by quality)
    for record in fastq_records
        sequence = convert_sequence(FASTX.FASTQ.sequence(record))
        quality_scores = collect(UInt8, FASTX.FASTQ.quality_scores(record))
        
        # Check minimum quality threshold
        mean_qual = Statistics.mean(quality_scores)
        if mean_qual >= min_quality
            # For double-strand mode, use canonical representation
            canonical_sequence, canonical_quality = if graph_mode == DoubleStrand && 
                                                       sequence isa Union{BioSequences.LongDNA, BioSequences.LongRNA}
                rc_sequence = BioSequences.reverse_complement(sequence)
                if sequence <= rc_sequence
                    (sequence, quality_scores)
                else
                    (rc_sequence, reverse(quality_scores))
                end
            else
                (sequence, quality_scores)
            end
            
            # Create vertex data
            vertex_data = QualityBioSequenceVertexData(canonical_sequence, canonical_quality, [1])
            graph[canonical_sequence] = vertex_data
        end
    end
    
    # Add edges based on quality-weighted overlaps
    _add_quality_biosequence_edges!(graph, fastq_records, graph_mode, min_overlap, min_quality)
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Qualmer graph to a quality-aware BioSequence graph by collapsing linear paths.

This is the primary method for creating quality-aware BioSequence graphs from Qualmer graphs,
following the 6-graph hierarchy where FASTQ graphs are simplifications of Qualmer graphs
with quality retention.

# Arguments
- `qualmer_graph`: MetaGraphsNext Qualmer graph to convert
- `min_path_length`: Minimum path length to keep (default: 2)

# Returns
- MetaGraphsNext.MetaGraph with quality-aware BioSequence vertices

# Example
```julia
# Start with qualmer graph
qualmer_graph = build_qualmer_graph(fastq_records)

# Convert to quality-aware BioSequence graph
quality_graph = qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
```
"""
function qualmer_graph_to_quality_biosequence_graph(qualmer_graph::MetaGraphsNext.MetaGraph; 
                                                   min_path_length::Int=2)
    
    # Determine sequence type from qualmer graph
    qualmer_labels = collect(MetaGraphsNext.labels(qualmer_graph))
    if isempty(qualmer_labels)
        throw(ArgumentError("Cannot convert empty qualmer graph"))
    end
    
    # Get first vertex to determine sequence type
    first_vertex_data = qualmer_graph[first(qualmer_labels)]
    first_qualmer = first_vertex_data.canonical_qualmer
    
    # Determine corresponding BioSequence type
    biosequence_type = if first_qualmer.kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_qualmer.kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else
        BioSequences.LongAA
    end
    
    # Create new quality-aware BioSequence graph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=biosequence_type,
        vertex_data_type=QualityBioSequenceVertexData{biosequence_type},
        edge_data_type=QualityBioSequenceEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Find linear paths in qualmer graph
    paths = _find_linear_qualmer_paths(qualmer_graph, min_path_length)
    
    # Convert each path to a quality-aware BioSequence
    for path in paths
        if length(path) >= min_path_length
            # Reconstruct sequence and quality from qualmer path
            sequence, quality_scores = _reconstruct_sequence_and_quality_from_qualmer_path(path, qualmer_graph, biosequence_type)
            
            # Calculate coverage from constituent qualmers
            coverage = _calculate_qualmer_path_coverage(path, qualmer_graph)
            
            # Create vertex data
            vertex_data = QualityBioSequenceVertexData(
                sequence, 
                quality_scores,
                coverage, 
                [string(vertex_label) for vertex_label in path]
            )
            
            graph[sequence] = vertex_data
        end
    end
    
    # Add edges between sequences based on qualmer transitions
    _add_quality_biosequence_edges_from_qualmer_paths!(graph, paths, qualmer_graph)
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert quality-aware BioSequence vertices back to FASTQ records.

This function demonstrates the key feature of FASTQ graphs - they maintain per-base 
quality information and can be converted back to FASTQ format.

# Arguments
- `graph`: Quality-aware BioSequence graph
- `output_file`: Path to output FASTQ file (optional)

# Returns
- Vector of FASTX.FASTQ.Record objects

# Example
```julia
# Convert graph back to FASTQ
fastq_records = quality_biosequence_graph_to_fastq(graph)

# Or write directly to file
quality_biosequence_graph_to_fastq(graph, "output.fastq")
```
"""
function quality_biosequence_graph_to_fastq(graph::MetaGraphsNext.MetaGraph, output_file::Union{Nothing,AbstractString}=nothing)
    
    fastq_records = FASTX.FASTQ.Record[]
    
    # Convert each vertex to FASTQ record
    for (i, vertex_label) in enumerate(MetaGraphsNext.labels(graph))
        vertex_data = graph[vertex_label]
        
        # Create identifier
        identifier = "contig_$i"
        
        # Get sequence and quality
        sequence = string(vertex_data.sequence)
        quality_chars = [Char(q + 33) for q in vertex_data.quality_scores]  # Convert to ASCII
        
        # Create FASTQ record
        record = FASTX.FASTQ.Record(identifier, sequence, String(quality_chars))
        push!(fastq_records, record)
    end
    
    # Write to file if specified
    if !isnothing(output_file)
        FASTX.FASTQ.open(output_file, "w") do writer
            for record in fastq_records
                write(writer, record)
            end
        end
    end
    
    return fastq_records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a quality-aware BioSequence graph to GFA format with quality information.

# Arguments
- `graph`: Quality-aware BioSequence graph
- `output_file`: Path to output GFA file

# Example
```julia
write_quality_biosequence_gfa(graph, "assembly_with_quality.gfa")
```
"""
function write_quality_biosequence_gfa(graph::MetaGraphsNext.MetaGraph, output_file::AbstractString)
    open(output_file, "w") do io
        # Write header
        println(io, "H\tVN:Z:1.0")
        
        # Write segments (vertices) with quality information
        sequence_to_id = Dict{Any, String}()
        for (i, vertex_label) in enumerate(MetaGraphsNext.labels(graph))
            vertex_data = graph[vertex_label]
            seg_id = "s$i"
            sequence_to_id[vertex_label] = seg_id
            
            # Include quality as optional field
            quality_str = String([Char(q + 33) for q in vertex_data.quality_scores])
            mean_qual = round(vertex_data.mean_quality, digits=2)
            joint_prob = round(vertex_data.joint_probability, digits=6)
            
            println(io, "S\t$seg_id\t$(string(vertex_label))\tQL:Z:$quality_str\tMQ:f:$mean_qual\tJP:f:$joint_prob")
        end
        
        # Write links (edges) with quality information
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            edge_data = graph[src, dst]
            src_id = sequence_to_id[src]
            dst_id = sequence_to_id[dst]
            
            # Convert strand orientations to GFA format
            src_orient = edge_data.src_strand == Forward ? "+" : "-"
            dst_orient = edge_data.dst_strand == Forward ? "+" : "-"
            
            # Use overlap length as CIGAR (simplified)
            overlap_cigar = "$(edge_data.overlap_length)M"
            
            # Include quality weight as optional field
            quality_weight = round(edge_data.quality_weight, digits=6)
            overlap_quality = round(edge_data.overlap_quality, digits=2)
            
            println(io, "L\t$src_id\t$src_orient\t$dst_id\t$dst_orient\t$overlap_cigar\tQW:f:$quality_weight\tOQ:f:$overlap_quality")
        end
    end
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
Add quality-weighted edges to BioSequence graph.
"""
function _add_quality_biosequence_edges!(graph::MetaGraphsNext.MetaGraph, 
                                        fastq_records::Vector{FASTX.FASTQ.Record},
                                        graph_mode::GraphMode, min_overlap::Int, min_quality::UInt8)
    
    vertices = collect(MetaGraphsNext.labels(graph))
    
    # Check all pairs of sequences for quality-weighted overlaps
    for i in 1:length(vertices)
        for j in (i+1):length(vertices)
            seq1 = vertices[i]
            seq2 = vertices[j]
            
            # Get vertex data for quality information
            vertex_data1 = graph[seq1]
            vertex_data2 = graph[seq2]
            
            # Check for overlap between sequences
            overlap_info = _find_quality_sequence_overlap(seq1, seq2, vertex_data1, vertex_data2, min_overlap)
            
            if !isnothing(overlap_info)
                overlap_length, src_strand, dst_strand, overlap_quality = overlap_info
                
                # Calculate quality weight
                quality_weight = vertex_data1.joint_probability * vertex_data2.joint_probability
                
                # Create edge data
                edge_data = QualityBioSequenceEdgeData(
                    overlap_length, src_strand, dst_strand, 
                    Float64(overlap_length), quality_weight, 1, overlap_quality
                )
                
                # Add edge
                graph[seq1, seq2] = edge_data
            end
        end
    end
end

"""
Find linear paths in a qualmer graph.
"""
function _find_linear_qualmer_paths(qualmer_graph::MetaGraphsNext.MetaGraph, min_path_length::Int)
    paths = Vector{Vector{Any}}()
    visited = Set{Any}()
    
    for start_vertex in MetaGraphsNext.labels(qualmer_graph)
        if start_vertex in visited
            continue
        end
        
        # Find linear path starting from this vertex
        path = [start_vertex]
        current = start_vertex
        
        # Extend path forward
        while true
            # Find outgoing edges sorted by quality weight
            outgoing = []
            for (src, dst) in MetaGraphsNext.edge_labels(qualmer_graph)
                if src == current && dst ∉ visited
                    edge_data = qualmer_graph[src, dst]
                    push!(outgoing, (dst, edge_data.quality_weight))
                end
            end
            
            # Sort by quality weight (descending)
            sort!(outgoing, by=x -> x[2], rev=true)
            
            # Only continue if exactly one high-quality outgoing edge
            if length(outgoing) == 1
                next_vertex, _ = outgoing[1]
                
                # Check if next vertex has exactly one high-quality incoming edge
                incoming = []
                for (src, dst) in MetaGraphsNext.edge_labels(qualmer_graph)
                    if dst == next_vertex
                        edge_data = qualmer_graph[src, dst]
                        push!(incoming, (src, edge_data.quality_weight))
                    end
                end
                
                if length(incoming) == 1
                    push!(path, next_vertex)
                    current = next_vertex
                else
                    break
                end
            else
                break
            end
        end
        
        # Mark all vertices in path as visited
        for vertex in path
            push!(visited, vertex)
        end
        
        # Add path if it meets minimum length
        if length(path) >= min_path_length
            push!(paths, path)
        end
    end
    
    return paths
end

"""
Reconstruct sequence and quality from qualmer path.
"""
function _reconstruct_sequence_and_quality_from_qualmer_path(qualmer_path::Vector, 
                                                           qualmer_graph::MetaGraphsNext.MetaGraph, 
                                                           biosequence_type::Type)
    if isempty(qualmer_path)
        return biosequence_type(), UInt8[]
    end
    
    # Get first qualmer data
    first_vertex_data = qualmer_graph[qualmer_path[1]]
    first_qualmer = first_vertex_data.canonical_qualmer
    
    # Start with first qualmer
    sequence_str = string(first_qualmer.kmer)
    quality_scores = collect(UInt8, first_qualmer.qualities)
    
    # Add one base + quality from each subsequent qualmer
    for i in 2:length(qualmer_path)
        vertex_data = qualmer_graph[qualmer_path[i]]
        qualmer = vertex_data.canonical_qualmer
        
        # Add last base and quality
        kmer_str = string(qualmer.kmer)
        if length(kmer_str) > 0
            sequence_str *= kmer_str[end]
            push!(quality_scores, qualmer.qualities[end])
        end
    end
    
    # Convert to BioSequence
    sequence = biosequence_type(sequence_str)
    
    return sequence, quality_scores
end

"""
Calculate coverage for qualmer path.
"""
function _calculate_qualmer_path_coverage(qualmer_path::Vector, qualmer_graph::MetaGraphsNext.MetaGraph)
    coverage = Int[]
    
    for vertex_label in qualmer_path
        vertex_data = qualmer_graph[vertex_label]
        push!(coverage, vertex_data.coverage)
    end
    
    return coverage
end

"""
Add edges between quality-aware BioSequences based on qualmer path relationships.
"""
function _add_quality_biosequence_edges_from_qualmer_paths!(graph::MetaGraphsNext.MetaGraph, 
                                                          paths::Vector{Vector{Any}}, 
                                                          qualmer_graph::MetaGraphsNext.MetaGraph)
    
    # Create mapping from qualmer to path
    qualmer_to_path = Dict{Any, Int}()
    for (path_idx, path) in enumerate(paths)
        for vertex_label in path
            qualmer_to_path[vertex_label] = path_idx
        end
    end
    
    # Check qualmer graph edges for inter-path connections
    for (src_vertex, dst_vertex) in MetaGraphsNext.edge_labels(qualmer_graph)
        src_path_idx = get(qualmer_to_path, src_vertex, 0)
        dst_path_idx = get(qualmer_to_path, dst_vertex, 0)
        
        # If vertices are in different paths, create edge between sequences
        if src_path_idx != dst_path_idx && src_path_idx > 0 && dst_path_idx > 0
            src_path = paths[src_path_idx]
            dst_path = paths[dst_path_idx]
            
            # Reconstruct sequences
            biosequence_type = typeof(first(MetaGraphsNext.labels(graph)))
            src_seq, src_qual = _reconstruct_sequence_and_quality_from_qualmer_path(src_path, qualmer_graph, biosequence_type)
            dst_seq, dst_qual = _reconstruct_sequence_and_quality_from_qualmer_path(dst_path, qualmer_graph, biosequence_type)
            
            # Get edge data from qualmer graph
            qualmer_edge_data = qualmer_graph[src_vertex, dst_vertex]
            
            # Create quality-aware BioSequence edge
            if haskey(graph, src_seq) && haskey(graph, dst_seq)
                # Get qualmer data for overlap quality
                src_vertex_data = qualmer_graph[src_vertex]
                dst_vertex_data = qualmer_graph[dst_vertex]
                
                # Calculate overlap quality (mean of connecting qualmers)
                overlap_quality = Statistics.mean([
                    Statistics.mean(src_vertex_data.canonical_qualmer.qualities),
                    Statistics.mean(dst_vertex_data.canonical_qualmer.qualities)
                ])
                
                # Use qualmer size as overlap length approximation
                k_size = length(src_vertex_data.canonical_qualmer.kmer)
                edge_data = QualityBioSequenceEdgeData(
                    k_size - 1,  # qualmer overlap
                    qualmer_edge_data.src_strand,
                    qualmer_edge_data.dst_strand,
                    qualmer_edge_data.weight,
                    qualmer_edge_data.quality_weight,
                    1,
                    overlap_quality
                )
                
                graph[src_seq, dst_seq] = edge_data
            end
        end
    end
end

"""
Find quality-weighted overlap between two sequences.
"""
function _find_quality_sequence_overlap(seq1, seq2, vertex_data1::QualityBioSequenceVertexData, 
                                      vertex_data2::QualityBioSequenceVertexData, min_overlap::Int)
    seq1_str = string(seq1)
    seq2_str = string(seq2)
    
    # Check suffix of seq1 against prefix of seq2
    for overlap_len in min_overlap:min(length(seq1_str), length(seq2_str))
        suffix = seq1_str[end-overlap_len+1:end]
        prefix = seq2_str[1:overlap_len]
        
        if suffix == prefix
            # Calculate overlap quality
            overlap_quality = Statistics.mean([
                Statistics.mean(vertex_data1.quality_scores[end-overlap_len+1:end]),
                Statistics.mean(vertex_data2.quality_scores[1:overlap_len])
            ])
            
            return (overlap_len, Forward, Forward, overlap_quality)
        end
    end
    
    # Check reverse complement overlaps if applicable
    if seq1 isa Union{BioSequences.LongDNA, BioSequences.LongRNA}
        rc_seq1_str = string(BioSequences.reverse_complement(seq1))
        rc_seq2_str = string(BioSequences.reverse_complement(seq2))
        
        # Check rc_seq1 suffix against seq2 prefix
        for overlap_len in min_overlap:min(length(rc_seq1_str), length(seq2_str))
            suffix = rc_seq1_str[end-overlap_len+1:end]
            prefix = seq2_str[1:overlap_len]
            
            if suffix == prefix
                # Calculate overlap quality (with reverse complement adjustment)
                overlap_quality = Statistics.mean([
                    Statistics.mean(reverse(vertex_data1.quality_scores)[end-overlap_len+1:end]),
                    Statistics.mean(vertex_data2.quality_scores[1:overlap_len])
                ])
                
                return (overlap_len, Reverse, Forward, overlap_quality)
            end
        end
        
        # Check seq1 suffix against rc_seq2 prefix
        for overlap_len in min_overlap:min(length(seq1_str), length(rc_seq2_str))
            suffix = seq1_str[end-overlap_len+1:end]
            prefix = rc_seq2_str[1:overlap_len]
            
            if suffix == prefix
                # Calculate overlap quality
                overlap_quality = Statistics.mean([
                    Statistics.mean(vertex_data1.quality_scores[end-overlap_len+1:end]),
                    Statistics.mean(reverse(vertex_data2.quality_scores)[1:overlap_len])
                ])
                
                return (overlap_len, Forward, Reverse, overlap_quality)
            end
        end
    end
    
    return nothing
end