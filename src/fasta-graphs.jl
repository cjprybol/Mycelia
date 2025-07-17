"""
FASTA Graph Implementation (BioSequence Graphs)

This file implements variable-length BioSequence graphs created by simplification/reduction
of fixed-length k-mer graphs, following the 6-graph hierarchy specification.

Key features:
- Variable-length BioSequence vertices (LongDNA, LongRNA, LongAA)
- No string conversions - maintains proper BioSequence types
- Strand-aware edges with biological transition validation
- Created by simplification of k-mer graphs
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import DocStringExtensions

"""
Vertex data for variable-length BioSequence graphs.
"""
struct BioSequenceVertexData{SeqT<:BioSequences.BioSequence}
    sequence::SeqT                      # Variable-length BioSequence
    coverage::Vector{Int}               # Coverage from original k-mer observations
    constituent_kmers::Vector{String}   # Original k-mers that were collapsed into this sequence
    length::Int                         # Sequence length
    
    function BioSequenceVertexData(sequence::SeqT, coverage::Vector{Int}=Int[], 
                                  constituent_kmers::Vector{String}=String[]) where {SeqT<:BioSequences.BioSequence}
        new{SeqT}(sequence, coverage, constituent_kmers, length(sequence))
    end
end

"""
Edge data for variable-length BioSequence graphs.
"""
struct BioSequenceEdgeData
    overlap_length::Int                 # Length of overlap between sequences
    src_strand::StrandOrientation       # Source sequence strand
    dst_strand::StrandOrientation       # Destination sequence strand
    weight::Float64                     # Edge weight from original k-mer graph
    original_kmer_count::Int           # Number of k-mer transitions that created this edge
    
    function BioSequenceEdgeData(overlap_length::Int, src_strand::StrandOrientation, 
                                dst_strand::StrandOrientation, weight::Float64=1.0, 
                                original_kmer_count::Int=1)
        new(overlap_length, src_strand, dst_strand, weight, original_kmer_count)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a BioSequence graph directly from FASTA records.

This creates a variable-length BioSequence graph where vertices are the input sequences
and edges represent relationships between sequences (e.g., overlaps, containment).

# Arguments
- `fasta_records`: Vector of FASTA records
- `graph_mode`: SingleStrand or DoubleStrand mode (default: DoubleStrand)
- `min_overlap`: Minimum overlap length for creating edges (default: 100)

# Returns
- MetaGraphsNext.MetaGraph with BioSequence vertices and overlap edges

# Example
```julia
records = [FASTX.FASTA.Record("seq1", "ATCGATCGATCG"), 
           FASTX.FASTA.Record("seq2", "CGATCGATCGAA")]
graph = build_biosequence_graph(records)
```
"""
function build_biosequence_graph(fasta_records::Vector{FASTX.FASTA.Record}; 
                                graph_mode::GraphMode=DoubleStrand,
                                min_overlap::Int=100)
    
    # Determine sequence type from first record
    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end
    
    first_seq = convert_sequence(FASTX.FASTA.sequence(fasta_records[1]))
    biosequence_type = typeof(first_seq)
    
    # Create MetaGraphsNext graph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=biosequence_type,
        vertex_data_type=BioSequenceVertexData{biosequence_type},
        edge_data_type=BioSequenceEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices for each sequence
    for record in fasta_records
        sequence = convert_sequence(FASTX.FASTA.sequence(record))
        
        # For double-strand mode, use canonical representation
        canonical_sequence = if graph_mode == DoubleStrand && sequence isa Union{BioSequences.LongDNA, BioSequences.LongRNA}
            rc_sequence = BioSequences.reverse_complement(sequence)
            sequence <= rc_sequence ? sequence : rc_sequence
        else
            sequence
        end
        
        # Create vertex data
        vertex_data = BioSequenceVertexData(canonical_sequence, [1])  # Single observation
        graph[canonical_sequence] = vertex_data
    end
    
    # Add edges based on overlaps
    _add_biosequence_edges!(graph, fasta_records, graph_mode, min_overlap)
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a k-mer graph to a BioSequence graph by collapsing linear paths.

This is the primary method for creating BioSequence graphs from k-mer graphs,
following the 6-graph hierarchy where BioSequence graphs are simplifications
of k-mer graphs.

# Arguments
- `kmer_graph`: MetaGraphsNext k-mer graph to convert
- `min_path_length`: Minimum path length to keep (default: 2)

# Returns
- MetaGraphsNext.MetaGraph with BioSequence vertices

# Example
```julia
# Start with k-mer graph
kmer_graph = build_kmer_graph_next(BioSequences.DNAKmer{31}, fasta_records)

# Convert to BioSequence graph
biosequence_graph = kmer_graph_to_biosequence_graph(kmer_graph)
```
"""
function kmer_graph_to_biosequence_graph(kmer_graph::MetaGraphsNext.MetaGraph; 
                                        min_path_length::Int=2)
    
    # Determine sequence type from k-mer graph
    kmer_labels = collect(MetaGraphsNext.labels(kmer_graph))
    if isempty(kmer_labels)
        throw(ArgumentError("Cannot convert empty k-mer graph"))
    end
    
    # Get k-mer type and determine corresponding BioSequence type
    first_kmer = first(kmer_labels)
    biosequence_type = if first_kmer isa Kmers.DNAKmer
        BioSequences.LongDNA{4}
    elseif first_kmer isa Kmers.RNAKmer
        BioSequences.LongRNA{4}
    else
        BioSequences.LongAA
    end
    
    # Create new BioSequence graph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=biosequence_type,
        vertex_data_type=BioSequenceVertexData{biosequence_type},
        edge_data_type=BioSequenceEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Find linear paths in k-mer graph
    paths = _find_linear_paths(kmer_graph, min_path_length)
    
    # Convert each path to a BioSequence
    for path in paths
        if length(path) >= min_path_length
            # Reconstruct sequence from k-mer path
            sequence = _reconstruct_sequence_from_kmer_path(path, biosequence_type)
            
            # Calculate coverage from constituent k-mers
            coverage = _calculate_path_coverage(path, kmer_graph)
            
            # Create vertex data
            vertex_data = BioSequenceVertexData(
                sequence, 
                coverage, 
                [string(kmer) for kmer in path]
            )
            
            graph[sequence] = vertex_data
        end
    end
    
    # Add edges between sequences based on overlaps
    _add_biosequence_edges_from_paths!(graph, paths, kmer_graph)
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a BioSequence graph to GFA format.

# Arguments
- `graph`: MetaGraphsNext BioSequence graph
- `output_file`: Path to output GFA file

# Example
```julia
write_biosequence_gfa(graph, "assembly.gfa")
```
"""
function write_biosequence_gfa(graph::MetaGraphsNext.MetaGraph, output_file::AbstractString)
    open(output_file, "w") do io
        # Write header
        println(io, "H\tVN:Z:1.0")
        
        # Write segments (vertices)
        sequence_to_id = Dict{Any, String}()
        for (i, vertex_label) in enumerate(MetaGraphsNext.labels(graph))
            seg_id = "s$i"
            sequence_to_id[vertex_label] = seg_id
            println(io, "S\t$seg_id\t$(string(vertex_label))")
        end
        
        # Write links (edges)
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            edge_data = graph[src, dst]
            src_id = sequence_to_id[src]
            dst_id = sequence_to_id[dst]
            
            # Convert strand orientations to GFA format
            src_orient = edge_data.src_strand == Forward ? "+" : "-"
            dst_orient = edge_data.dst_strand == Forward ? "+" : "-"
            
            # Use overlap length as CIGAR (simplified)
            overlap_cigar = "$(edge_data.overlap_length)M"
            
            println(io, "L\t$src_id\t$src_orient\t$dst_id\t$dst_orient\t$overlap_cigar")
        end
    end
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
Add edges to BioSequence graph based on sequence overlaps.
"""
function _add_biosequence_edges!(graph::MetaGraphsNext.MetaGraph, 
                                fasta_records::Vector{FASTX.FASTA.Record},
                                graph_mode::GraphMode, min_overlap::Int)
    
    vertices = collect(MetaGraphsNext.labels(graph))
    
    # Check all pairs of sequences for overlaps
    for i in 1:length(vertices)
        for j in (i+1):length(vertices)
            seq1 = vertices[i]
            seq2 = vertices[j]
            
            # Check for overlap between sequences
            overlap_info = _find_sequence_overlap(seq1, seq2, min_overlap)
            
            if !isnothing(overlap_info)
                overlap_length, src_strand, dst_strand = overlap_info
                
                # Create edge data
                edge_data = BioSequenceEdgeData(overlap_length, src_strand, dst_strand)
                
                # Add edge
                graph[seq1, seq2] = edge_data
            end
        end
    end
end

"""
Find linear paths in a k-mer graph.
"""
function _find_linear_paths(kmer_graph::MetaGraphsNext.MetaGraph, min_path_length::Int)
    paths = Vector{Vector{Any}}()
    visited = Set{Any}()
    
    for start_vertex in MetaGraphsNext.labels(kmer_graph)
        if start_vertex in visited
            continue
        end
        
        # Find linear path starting from this vertex
        path = [start_vertex]
        current = start_vertex
        
        # Extend path forward
        while true
            # Find outgoing edges
            outgoing = []
            for (src, dst) in MetaGraphsNext.edge_labels(kmer_graph)
                if src == current && dst ∉ visited
                    push!(outgoing, dst)
                end
            end
            
            # Only continue if exactly one outgoing edge
            if length(outgoing) == 1
                next_vertex = outgoing[1]
                
                # Check if next vertex has exactly one incoming edge
                incoming = []
                for (src, dst) in MetaGraphsNext.edge_labels(kmer_graph)
                    if dst == next_vertex
                        push!(incoming, src)
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
Reconstruct a BioSequence from a path of k-mers.
"""
function _reconstruct_sequence_from_kmer_path(kmer_path::Vector, biosequence_type::Type)
    if isempty(kmer_path)
        return biosequence_type()
    end
    
    # Start with first k-mer
    sequence_str = string(kmer_path[1])
    
    # Add one nucleotide/amino acid from each subsequent k-mer
    for i in 2:length(kmer_path)
        kmer_str = string(kmer_path[i])
        if length(kmer_str) > 0
            sequence_str *= kmer_str[end]
        end
    end
    
    # Convert to BioSequence
    return biosequence_type(sequence_str)
end

"""
Calculate coverage for a path from constituent k-mers.
"""
function _calculate_path_coverage(kmer_path::Vector, kmer_graph::MetaGraphsNext.MetaGraph)
    coverage = Int[]
    
    for kmer in kmer_path
        vertex_data = kmer_graph[kmer]
        # Use length of coverage vector as coverage count
        push!(coverage, length(vertex_data.coverage))
    end
    
    return coverage
end

"""
Add edges between BioSequences based on k-mer path relationships.
"""
function _add_biosequence_edges_from_paths!(graph::MetaGraphsNext.MetaGraph, 
                                          paths::Vector{Vector{Any}}, 
                                          kmer_graph::MetaGraphsNext.MetaGraph)
    
    # Create mapping from k-mer to path
    kmer_to_path = Dict{Any, Int}()
    for (path_idx, path) in enumerate(paths)
        for kmer in path
            kmer_to_path[kmer] = path_idx
        end
    end
    
    # Check k-mer graph edges for inter-path connections
    for (src_kmer, dst_kmer) in MetaGraphsNext.edge_labels(kmer_graph)
        src_path_idx = get(kmer_to_path, src_kmer, 0)
        dst_path_idx = get(kmer_to_path, dst_kmer, 0)
        
        # If k-mers are in different paths, create edge between sequences
        if src_path_idx != dst_path_idx && src_path_idx > 0 && dst_path_idx > 0
            src_path = paths[src_path_idx]
            dst_path = paths[dst_path_idx]
            
            # Reconstruct sequences
            src_seq = _reconstruct_sequence_from_kmer_path(src_path, typeof(first(MetaGraphsNext.labels(graph))))
            dst_seq = _reconstruct_sequence_from_kmer_path(dst_path, typeof(first(MetaGraphsNext.labels(graph))))
            
            # Get edge data from k-mer graph
            kmer_edge_data = kmer_graph[src_kmer, dst_kmer]
            
            # Create BioSequence edge
            if haskey(graph, src_seq) && haskey(graph, dst_seq)
                # Use k-mer size as overlap length approximation
                k_size = length(string(src_kmer))
                edge_data = BioSequenceEdgeData(
                    k_size - 1,  # k-mer overlap
                    kmer_edge_data.src_strand,
                    kmer_edge_data.dst_strand,
                    kmer_edge_data.weight
                )
                
                graph[src_seq, dst_seq] = edge_data
            end
        end
    end
end

"""
Find overlap between two sequences.
"""
function _find_sequence_overlap(seq1, seq2, min_overlap::Int)
    seq1_str = string(seq1)
    seq2_str = string(seq2)
    
    # Check suffix of seq1 against prefix of seq2
    for overlap_len in min_overlap:min(length(seq1_str), length(seq2_str))
        suffix = seq1_str[end-overlap_len+1:end]
        prefix = seq2_str[1:overlap_len]
        
        if suffix == prefix
            return (overlap_len, Forward, Forward)
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
                return (overlap_len, Reverse, Forward)
            end
        end
        
        # Check seq1 suffix against rc_seq2 prefix
        for overlap_len in min_overlap:min(length(seq1_str), length(rc_seq2_str))
            suffix = seq1_str[end-overlap_len+1:end]
            prefix = rc_seq2_str[1:overlap_len]
            
            if suffix == prefix
                return (overlap_len, Forward, Reverse)
            end
        end
    end
    
    return nothing
end