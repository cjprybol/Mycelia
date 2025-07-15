"""
GFA I/O functionality for MetaGraphsNext-based strand-aware k-mer graphs.

This module provides reading and writing of GFA (Graphical Fragment Assembly) format
files for the next-generation MetaGraphsNext implementation with strand-aware edges.
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
using DocStringExtensions

"""
$(TYPEDSIGNATURES)

Write a MetaGraphsNext k-mer graph to GFA (Graphical Fragment Assembly) format.

This function exports strand-aware k-mer graphs to standard GFA format, handling:
- Canonical k-mer vertices as segments (S lines)
- Strand-aware edges as links (L lines) with proper orientations
- Coverage information as depth annotations

# Arguments
- `graph`: MetaGraphsNext k-mer graph with strand-aware edges
- `outfile`: Path where the GFA file should be written

# Returns
- Path to the written GFA file

# GFA Format
The output follows GFA v1.0 specification:
- Header (H) line with version
- Segment (S) lines: vertex_id, canonical_k-mer_sequence, depth
- Link (L) lines: source_id, source_orientation, target_id, target_orientation, overlap

# Example
```julia
graph = build_kmer_graph_next(DNAKmer{3}, observations)
write_gfa_next(graph, "assembly.gfa")
```
"""
function write_gfa_next(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString)
    open(outfile, "w") do io
        # Write GFA header
        println(io, "H\tVN:Z:1.0\tMY:Z:Mycelia-Next")
        
        # Write segments (vertices) - canonical k-mers
        vertex_id_map = Dict{String, Int}()
        for (i, label) in enumerate(MetaGraphsNext.labels(graph))
            vertex_id_map[label] = i
            vertex_data = graph[label]
            
            # Calculate depth from coverage
            depth = length(vertex_data.coverage)
            
            # Write segment line: S<tab>id<tab>sequence<tab>optional_fields
            println(io, "S\t$i\t$(vertex_data.canonical_kmer)\tDP:f:$depth")
        end
        
        # Calculate k-mer size for overlap
        if !isempty(MetaGraphsNext.labels(graph))
            first_label = first(MetaGraphsNext.labels(graph))
            k = length(graph[first_label].canonical_kmer)
            overlap = k - 1
        else
            overlap = 0
        end
        
        # Write links (edges) - strand-aware transitions
        for edge_labels in MetaGraphsNext.edge_labels(graph)
            if length(edge_labels) == 2
                src_label, dst_label = edge_labels
                edge_data = graph[src_label, dst_label]
                
                src_id = vertex_id_map[src_label]
                dst_id = vertex_id_map[dst_label]
                
                # Convert strand orientations to GFA format
                src_orientation = edge_data.src_strand == Forward ? '+' : '-'
                dst_orientation = edge_data.dst_strand == Forward ? '+' : '-'
                
                # Write link line: L<tab>src<tab>src_orient<tab>dst<tab>dst_orient<tab>overlap
                println(io, "L\t$src_id\t$src_orientation\t$dst_id\t$dst_orientation\t$(overlap)M")
            end
        end
    end
    
    return outfile
end

"""
$(TYPEDSIGNATURES)

Read a GFA file and convert it to a MetaGraphsNext strand-aware k-mer graph.

This function parses GFA format files and creates a strand-aware k-mer graph compatible
with the next-generation implementation.

# Arguments
- `gfa_file`: Path to input GFA file
- `graph_mode`: GraphMode (SingleStrand or DoubleStrand, default: DoubleStrand)

# Returns
- MetaGraphsNext.MetaGraph with strand-aware edges

# GFA Format Support
Supports GFA v1.0 with:
- Header (H) lines (ignored)
- Segment (S) lines: parsed as canonical k-mer vertices
- Link (L) lines: parsed as strand-aware edges
- Path (P) lines: stored as metadata (future use)

# Example
```julia
graph = read_gfa_next("assembly.gfa")
# Or with specific mode
graph = read_gfa_next("assembly.gfa", SingleStrand)
```
"""
function read_gfa_next(gfa_file::AbstractString, graph_mode::GraphMode=DoubleStrand)
    # Parse GFA file content
    segments = Dict{String, String}()  # id -> sequence
    links = Vector{Tuple{String, Bool, String, Bool}}()  # (src_id, src_forward, dst_id, dst_forward)
    paths = Dict{String, Vector{String}}()  # path_name -> vertex_ids
    
    for line in eachline(gfa_file)
        fields = split(line, '\t')
        if isempty(fields)
            continue
        end
        
        line_type = first(fields)
        
        if line_type == "H"
            # Header line - skip for now
            continue
        elseif line_type == "S"
            # Segment line: S<tab>id<tab>sequence<tab>optional_fields
            if length(fields) >= 3
                seg_id = fields[2]
                sequence = fields[3]
                segments[seg_id] = sequence
            end
        elseif line_type == "L"
            # Link line: L<tab>src<tab>src_orient<tab>dst<tab>dst_orient<tab>overlap
            if length(fields) >= 6
                src_id = fields[2]
                src_orient = fields[3] == "+"
                dst_id = fields[4]
                dst_orient = fields[5] == "+"
                push!(links, (src_id, src_orient, dst_id, dst_orient))
            end
        elseif line_type == "P"
            # Path line: P<tab>path_name<tab>path<tab>overlaps
            if length(fields) >= 3
                path_name = fields[2]
                # Parse path string (removes +/- orientations for now)
                path_vertices = split(replace(fields[3], r"[+-]" => ""), ',')
                paths[path_name] = string.(path_vertices)
            end
        elseif line_type == "A"
            # Assembly info line (hifiasm) - skip for now
            continue
        else
            @warn "Unknown GFA line type: $line_type in line: $line"
        end
    end
    
    # Create MetaGraphsNext graph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=KmerVertexData,
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices (segments)
    for (seg_id, sequence) in segments
        # Determine canonical sequence based on graph mode
        canonical_seq = if graph_mode == DoubleStrand
            # Use canonical representation for double-strand mode
            dna_seq = DNASequence(sequence)
            rc_seq = reverse_complement(dna_seq)
            string(dna_seq <= rc_seq ? dna_seq : rc_seq)
        else
            # Use sequence as-is for single-strand mode
            sequence
        end
        
        # Create vertex with empty coverage (will be populated if we have observations)
        graph[seg_id] = KmerVertexData(canonical_seq)
    end
    
    # Add edges (links)
    for (src_id, src_forward, dst_id, dst_forward) in links
        if src_id in MetaGraphsNext.labels(graph) && dst_id in MetaGraphsNext.labels(graph)
            # Convert GFA orientations to StrandOrientation
            src_strand = src_forward ? Forward : Reverse
            dst_strand = dst_forward ? Forward : Reverse
            
            # Create strand-aware edge
            graph[src_id, dst_id] = KmerEdgeData(src_strand, dst_strand)
        else
            @warn "Link references unknown segment: $src_id -> $dst_id"
        end
    end
    
    # Store paths as graph metadata if needed (future enhancement)
    # MetaGraphsNext doesn't have global properties like MetaGraphs, so we'd need
    # a different approach for storing paths
    
    return graph
end

"""
$(TYPEDSIGNATURES)

Convert a legacy MetaGraphs GFA to next-generation MetaGraphsNext format.

This convenience function reads a GFA file using the legacy parser and converts
it to the new strand-aware format.

# Arguments
- `gfa_file`: Path to GFA file
- `graph_mode`: GraphMode for the output graph

# Returns
- MetaGraphsNext.MetaGraph in strand-aware format
"""
function convert_legacy_gfa_to_next(gfa_file::AbstractString, graph_mode::GraphMode=DoubleStrand)
    # Use legacy parser
    legacy_graph = parse_gfa(gfa_file)
    
    # Convert to next-generation format
    next_graph = legacy_to_next_graph(legacy_graph)
    
    return next_graph
end

# Export the functions
export write_gfa_next, read_gfa_next, convert_legacy_gfa_to_next
