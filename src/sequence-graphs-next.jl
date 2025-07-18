"""
MetaGraphsNext-based Sequence Graph Implementation

This file contains the next-generation implementation of sequence graphs using MetaGraphsNext.jl
with type-stable metadata structures, replacing the deprecated MetaGraphs.jl implementation.
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import DocStringExtensions

# Type-stable metadata structures for vertices and edges

"""
Strand orientation for k-mer observations and transitions.
- `Forward`: k-mer as observed (5' to 3')
- `Reverse`: reverse complement of k-mer (3' to 5')
"""
@enum StrandOrientation Forward=true Reverse=false

"""
Graph mode for handling strand information.
- `SingleStrand`: Sequences are single-stranded (RNA, amino acids, or directional DNA)
- `DoubleStrand`: Sequences are double-stranded DNA/RNA with canonical representation
"""
@enum GraphMode SingleStrand DoubleStrand

"""
Type-stable metadata for k-mer graph vertices.

Vertices always represent canonical k-mers for memory efficiency and cleaner graphs.
Strand information is tracked in the coverage data and edge transitions.

Fields:
- `coverage`: Vector of observation coverage data as (observation_id, position, strand_orientation) tuples
- `canonical_kmer`: The canonical k-mer (BioSequence type - NO string conversion)
"""
struct KmerVertexData{KmerT}
    coverage::Vector{Tuple{Int, Int, StrandOrientation}}  # (observation_id, position, strand_orientation)
    canonical_kmer::KmerT  # Actual k-mer type (DNAKmer, RNAKmer, AAKmer)
    
    # Constructor with default empty coverage
    KmerVertexData(canonical_kmer::KmerT) where {KmerT} = new{KmerT}(Vector{Tuple{Int, Int, StrandOrientation}}(), canonical_kmer)
end

"""
Type-stable metadata for k-mer graph edges.

Edges represent valid strand-aware transitions between canonical k-mers.
The transition is valid only if the strand orientations allow for proper overlap.

Fields:
- `coverage`: Vector of edge traversal observations with strand information
- `weight`: Edge weight/confidence score based on coverage
- `src_strand`: Required strand orientation of source k-mer for this transition
- `dst_strand`: Required strand orientation of destination k-mer for this transition
"""
struct KmerEdgeData
    coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}
    weight::Float64
    src_strand::StrandOrientation  # Required orientation of source k-mer
    dst_strand::StrandOrientation  # Required orientation of destination k-mer
    
    # Constructor with default weight calculation
    function KmerEdgeData(coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}},
                         src_strand::StrandOrientation, dst_strand::StrandOrientation)
        new(coverage, Float64(length(coverage)), src_strand, dst_strand)
    end
    
    # Default constructor
    KmerEdgeData(src_strand::StrandOrientation, dst_strand::StrandOrientation) = 
        new(Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}(), 
            0.0, src_strand, dst_strand)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a next-generation, type-stable k-mer graph using MetaGraphsNext.

This implementation uses canonical k-mers as vertices with strand-aware edges that respect
biological transition constraints for both single-strand and double-strand sequences.

# Arguments
- `kmer_type`: Type of k-mer (e.g., DNAKmer{K})
- `observations`: Vector of FASTA/FASTQ records
- `graph_mode`: `SingleStrand` for directional sequences, `DoubleStrand` for DNA (default)

# Returns
- `MetaGraphsNext.MetaGraph` with canonical vertices and strand-aware edges

# Details
- **Vertices**: Always canonical k-mers (lexicographically smaller of kmer/reverse_complement)
- **Edges**: Strand-aware transitions that respect biological constraints
- **SingleStrand mode**: Only forward-strand transitions allowed
- **DoubleStrand mode**: Both forward and reverse-complement transitions allowed
"""
function build_kmer_graph_next(kmer_type, observations::AbstractVector{<:Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}; 
                                graph_mode::GraphMode=DoubleStrand)
    
    # Count canonical k-mers from observations
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, observations)
    canonical_kmers = collect(keys(canonical_kmer_counts))
    
    if isempty(canonical_kmers)
        @warn "No k-mers found in observations"
        return _create_empty_kmer_graph(kmer_type)
    end
    
    # Infer the actual k-mer type with storage parameter from the data
    actual_kmer_type = eltype(canonical_kmers)
    
    # Create the MetaGraphsNext graph with type-stable metadata
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=actual_kmer_type,
        vertex_data_type=KmerVertexData{actual_kmer_type},
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices for each canonical k-mer (NO string conversion)
    for kmer in canonical_kmers
        graph[kmer] = KmerVertexData(kmer)
    end
    
    # Process observations to build strand-aware edges
    for (obs_idx, observation) in enumerate(observations)
        _add_observation_to_graph!(graph, observation, obs_idx, canonical_kmers, graph_mode)
    end
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add a sequence observation to an existing k-mer graph with strand-aware edge creation.

# Arguments
- `graph`: MetaGraphsNext k-mer graph with canonical vertices
- `observation`: FASTA/FASTQ record
- `obs_idx`: Observation index
- `canonical_kmers`: Vector of canonical k-mers in the graph
- `graph_mode`: SingleStrand or DoubleStrand mode
"""
function _add_observation_to_graph!(graph, observation, obs_idx, canonical_kmers, graph_mode)
    observed_sequence = FASTX.sequence(observation)
    k = length(first(canonical_kmers))
    
    if length(observed_sequence) < k
        observation_id = FASTX.identifier(observation)
        @warn "Skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        return
    end
    
    # Convert sequence to path through canonical k-mer graph with strand information
    observed_path = _sequence_to_canonical_path(canonical_kmers, observed_sequence, graph_mode)
    
    if isempty(observed_path)
        return
    end
    
    # Add coverage to first vertex
    first_canonical_kmer, first_strand = observed_path[1]
    _add_vertex_coverage!(graph, first_canonical_kmer, obs_idx, 1, first_strand)
    
    # Add strand-aware edges and vertex coverage for the rest of the path
    for i in 2:length(observed_path)
        curr_canonical_kmer, curr_strand = observed_path[i]
        prev_canonical_kmer, prev_strand = observed_path[i-1]
        
        # Add vertex coverage
        _add_vertex_coverage!(graph, curr_canonical_kmer, obs_idx, i, curr_strand)
        
        # Add or update strand-aware edge
        _add_strand_aware_edge!(graph, prev_canonical_kmer, curr_canonical_kmer, 
                               prev_strand, curr_strand,
                               (obs_idx, i-1, prev_strand), (obs_idx, i, curr_strand))
    end
end

"""
Helper function to add coverage data to a vertex.
"""
function _add_vertex_coverage!(graph, kmer, obs_idx, position, strand_orientation)
    vertex_data = graph[kmer]
    new_coverage = (obs_idx, position, strand_orientation)
    push!(vertex_data.coverage, new_coverage)
end

"""
Helper function to add strand-aware coverage data to an edge.

This function creates edges that respect strand orientation constraints.
Each edge represents a biologically valid transition between k-mers.
"""
function _add_strand_aware_edge!(graph, src_kmer, dst_kmer, src_strand, dst_strand, src_coverage, dst_coverage)
    # Create edge key that includes strand information
    edge_exists = false
    existing_edge_data = nothing
    
    # Check if an edge with this strand configuration already exists
    if haskey(graph, src_kmer, dst_kmer)
        existing_edge_data = graph[src_kmer, dst_kmer]
        # Check if this edge represents the same strand transition
        if existing_edge_data.src_strand == src_strand && existing_edge_data.dst_strand == dst_strand
            edge_exists = true
        end
    end
    
    if !edge_exists
        # Create new strand-aware edge
        graph[src_kmer, dst_kmer] = KmerEdgeData(src_strand, dst_strand)
        existing_edge_data = graph[src_kmer, dst_kmer]
    end
    
    # Add coverage to the edge
    new_coverage = (src_coverage, dst_coverage)
    push!(existing_edge_data.coverage, new_coverage)
    
    # Update weight based on coverage count
    graph[src_kmer, dst_kmer] = KmerEdgeData(existing_edge_data.coverage, 
                                            existing_edge_data.src_strand, 
                                            existing_edge_data.dst_strand)
end

"""
Helper function to create an empty k-mer graph.
"""
function _create_empty_kmer_graph(kmer_type)
    return MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=kmer_type,
        vertex_data_type=KmerVertexData{kmer_type},
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a sequence to a path through canonical k-mer space with strand awareness.

This is the key function that handles the distinction between single-strand and 
double-strand modes while maintaining canonical k-mer vertices.

# Arguments
- `canonical_kmers`: Vector of canonical k-mers available in the graph
- `sequence`: DNA/RNA sequence to convert
- `graph_mode`: SingleStrand or DoubleStrand mode

# Returns
- Vector of (canonical_kmer, strand_orientation) pairs representing the path

# Details
In DoubleStrand mode:
- Each observed k-mer is converted to its canonical form
- Strand orientation tracks whether the canonical form matches the observed k-mer
- Edges respect the biological constraint that transitions must maintain proper overlap

In SingleStrand mode:
- K-mers are used as-is (no reverse complement consideration)
- All strand orientations are Forward
- Suitable for RNA, amino acids, or directional DNA analysis
"""
function _sequence_to_canonical_path(canonical_kmers, sequence, graph_mode)
    k = length(first(canonical_kmers))
    path = Vector{Tuple{eltype(canonical_kmers), StrandOrientation}}()
    canonical_kmer_set = Set(canonical_kmers)
    
    # Extract k-mers from sequence
    for i in 1:(length(sequence) - k + 1)
        subseq = sequence[i:i+k-1]
        observed_kmer = typeof(first(canonical_kmers))(subseq)
        
        if graph_mode == SingleStrand
            # Single-strand mode: use k-mers as-is
            if observed_kmer in canonical_kmer_set
                push!(path, (observed_kmer, Forward))
            else
                @warn "K-mer $observed_kmer not found in graph at position $i (SingleStrand mode)"
            end
        else  # DoubleStrand mode
            # Double-strand mode: find canonical representation
            # Only valid for nucleic acids, not amino acids
            if typeof(observed_kmer) <: Kmers.Kmer{<:BioSequences.NucleicAcidAlphabet}
                rc_kmer = BioSequences.reverse_complement(observed_kmer)
            else
                # For amino acids, treat as single-strand even in DoubleStrand mode
                if observed_kmer in canonical_kmer_set
                    push!(path, (observed_kmer, Forward))
                    continue
                else
                    @warn "K-mer $observed_kmer not found in graph at position $i (amino acid in DoubleStrand mode)"
                    continue
                end
            end
            
            # Determine canonical k-mer (lexicographically smaller)
            if observed_kmer <= rc_kmer
                canonical_kmer = observed_kmer
                strand_orientation = Forward
            else
                canonical_kmer = rc_kmer
                strand_orientation = Reverse
            end
            
            if canonical_kmer in canonical_kmer_set
                push!(path, (canonical_kmer, strand_orientation))
            else
                @warn "Canonical k-mer $canonical_kmer not found in graph at position $i (DoubleStrand mode)"
            end
        end
    end
    
    return path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Validate that a transition between two k-mers with given strand orientations is biologically valid.

For a transition to be valid, the suffix of the source k-mer must match the prefix of the 
destination k-mer when accounting for strand orientations.

# Arguments
- `src_kmer`: Source canonical k-mer
- `dst_kmer`: Destination canonical k-mer  
- `src_strand`: Strand orientation of source k-mer
- `dst_strand`: Strand orientation of destination k-mer
- `k`: K-mer size

# Returns
- `Bool`: true if transition is biologically valid
"""
function _is_valid_transition(src_kmer, dst_kmer, src_strand, dst_strand, k)
    # Get the actual k-mer sequences considering strand orientation
    src_seq = src_strand == Forward ? src_kmer : BioSequences.reverse_complement(src_kmer)
    dst_seq = dst_strand == Forward ? dst_kmer : BioSequences.reverse_complement(dst_kmer)
    
    # Check if suffix of src matches prefix of dst (k-1 overlap)
    # For BioSequences, we can use slicing directly
    src_suffix = src_seq[2:end]  # Remove first nucleotide
    dst_prefix = dst_seq[1:end-1]  # Remove last nucleotide
    
    return src_suffix == dst_prefix
end

# Compatibility layer for migrating from legacy MetaGraphs implementation

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a legacy MetaGraphs-based k-mer graph to the next-generation MetaGraphsNext format.

This function provides a migration path from the deprecated MetaGraphs.jl implementation
to the new type-stable MetaGraphsNext.jl format.

# Arguments
- `legacy_graph`: MetaGraphs.MetaDiGraph from the old implementation

# Returns
- `MetaGraphsNext.MetaGraph` with equivalent structure and type-stable metadata
"""
function legacy_to_next_graph(legacy_graph, kmer_type=nothing)
    # Extract metadata from legacy graph
    stranded_kmers = legacy_graph.gprops[:stranded_kmers]
    k = legacy_graph.gprops[:k]
    
    # Determine k-mer type if not provided
    if kmer_type === nothing
        # Default to DNAKmer for backward compatibility
        kmer_type = Kmers.DNAKmer{k}
    end
    
    # Create next-generation graph
    next_graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=kmer_type,
        vertex_data_type=KmerVertexData{kmer_type},
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Migrate vertices to canonical representation
    for v in Graphs.vertices(legacy_graph)
        kmer_str = string(stranded_kmers[v])
        
        # Convert to canonical k-mer representation
        kmer_seq = kmer_type(kmer_str)
        
        # Only calculate reverse complement for nucleic acids
        if kmer_type <: Union{Kmers.DNAKmer, Kmers.RNAKmer}
            rc_kmer_seq = BioSequences.reverse_complement(kmer_seq)
            canonical_kmer = kmer_seq <= rc_kmer_seq ? kmer_seq : rc_kmer_seq
        else
            # For amino acids, no reverse complement
            canonical_kmer = kmer_seq
        end
        
        # Convert legacy coverage format to next-generation format with strand info
        legacy_coverage = get(legacy_graph.vprops[v], :coverage, [])
        next_coverage = Vector{Tuple{Int, Int, StrandOrientation}}()
        
        for (obs_idx, (pos, orientation)) in legacy_coverage
            strand_orientation = orientation ? Forward : Reverse
            push!(next_coverage, (obs_idx, pos, strand_orientation))
        end
        
        # Create or update canonical vertex
        if canonical_kmer in MetaGraphsNext.labels(next_graph)
            # Merge coverage with existing canonical vertex
            existing_vertex = next_graph[canonical_kmer]
            append!(existing_vertex.coverage, next_coverage)
        else
            next_graph[canonical_kmer] = KmerVertexData(canonical_kmer)
            existing_vertex = next_graph[canonical_kmer]
            append!(existing_vertex.coverage, next_coverage)
        end
    end
    
    # Migrate edges with strand awareness
    for edge in Graphs.edges(legacy_graph)
        src_kmer_str = string(stranded_kmers[Graphs.src(edge)])
        dst_kmer_str = string(stranded_kmers[Graphs.dst(edge)])
        
        # Convert to canonical representation and determine strand orientations
        src_kmer_seq = kmer_type(src_kmer_str)
        dst_kmer_seq = kmer_type(dst_kmer_str)
        
        # Handle canonical representation differently for nucleic acids vs amino acids
        if kmer_type <: Union{Kmers.DNAKmer, Kmers.RNAKmer}
            src_rc = BioSequences.reverse_complement(src_kmer_seq)
            dst_rc = BioSequences.reverse_complement(dst_kmer_seq)
            
            src_canonical = src_kmer_seq <= src_rc ? src_kmer_seq : src_rc
            dst_canonical = dst_kmer_seq <= dst_rc ? dst_kmer_seq : dst_rc
            
            src_strand = src_kmer_seq <= src_rc ? Forward : Reverse
            dst_strand = dst_kmer_seq <= dst_rc ? Forward : Reverse
        else
            # For amino acids, no reverse complement concept
            src_canonical = src_kmer_seq
            dst_canonical = dst_kmer_seq
            
            src_strand = Forward  # Always forward for amino acids
            dst_strand = Forward
        end
        
        # Convert legacy edge coverage to next-generation format
        legacy_edge_coverage = get(legacy_graph.eprops[edge], :coverage, [])
        next_edge_coverage = Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}()
        
        for (src_cov, dst_cov) in legacy_edge_coverage
            src_obs, (src_pos, src_orient) = src_cov
            dst_obs, (dst_pos, dst_orient) = dst_cov
            
            src_strand_cov = src_orient ? Forward : Reverse
            dst_strand_cov = dst_orient ? Forward : Reverse
            
            push!(next_edge_coverage, ((src_obs, src_pos, src_strand_cov), (dst_obs, dst_pos, dst_strand_cov)))
        end
        
        # Create strand-aware edge
        if !haskey(next_graph, src_canonical, dst_canonical)
            next_graph[src_canonical, dst_canonical] = KmerEdgeData(src_strand, dst_strand)
        end
        
        edge_data = next_graph[src_canonical, dst_canonical]
        append!(edge_data.coverage, next_edge_coverage)
        
        # Update with new coverage count
        next_graph[src_canonical, dst_canonical] = KmerEdgeData(edge_data.coverage, src_strand, dst_strand)
    end
    
    return next_graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Check if a graph is using the legacy MetaGraphs format.

# Arguments
- `graph`: Graph to check

# Returns
- `Bool`: true if legacy format, false if next-generation format
"""
function is_legacy_graph(graph)
    return graph isa MetaGraphs.MetaDiGraph || graph isa MetaGraphs.MetaGraph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Automatically convert legacy graphs to next-generation format if needed.

This is a convenience function that checks the graph type and converts if necessary.

# Arguments
- `graph`: Graph in either legacy or next-generation format

# Returns
- `MetaGraphsNext.MetaGraph`: Graph in next-generation format
"""
function ensure_next_graph(graph)
    if is_legacy_graph(graph)
        @info "Converting legacy MetaGraphs format to next-generation MetaGraphsNext format"
        return legacy_to_next_graph(graph)
    else
        return graph
    end
end

# Note: We don't export specific types/functions - use fully qualified names like Mycelia.build_kmer_graph_next

# GFA I/O functionality for MetaGraphsNext-based strand-aware k-mer graphs

"""
$(DocStringExtensions.TYPEDSIGNATURES)

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
        vertex_id_map = Dict()  # Generic dict to handle any k-mer type
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
$(DocStringExtensions.TYPEDSIGNATURES)

Read a GFA file and convert it to a MetaGraphsNext k-mer graph with fixed-length vertices.

This function parses GFA format files and creates a strand-aware k-mer graph compatible
with the next-generation implementation using fixed-length k-mer vertices.

# Arguments
- `gfa_file`: Path to input GFA file
- `kmer_type`: Type of k-mer to use (e.g., Kmers.DNAKmer{31})
- `graph_mode`: GraphMode (SingleStrand or DoubleStrand, default: DoubleStrand)

# Returns
- MetaGraphsNext.MetaGraph with k-mer vertices and strand-aware edges

# GFA Format Support
Supports GFA v1.0 with:
- Header (H) lines (ignored)
- Segment (S) lines: parsed as fixed-length k-mer vertices
- Link (L) lines: parsed as strand-aware edges
- Path (P) lines: stored as metadata (future use)

# Example
```julia
# Fixed-length k-mer graph
graph = read_gfa_next("assembly.gfa", Kmers.DNAKmer{31})
# Or with specific mode
graph = read_gfa_next("assembly.gfa", Kmers.DNAKmer{31}, SingleStrand)
```
"""
function read_gfa_next(gfa_file::AbstractString, kmer_type::Type, graph_mode::GraphMode=DoubleStrand)
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
        label_type=kmer_type,
        vertex_data_type=KmerVertexData{kmer_type},
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices (segments) - convert sequences to k-mer types
    for (seg_id, sequence) in segments
        # Convert sequence to k-mer type
        kmer_seq = kmer_type(sequence)
        
        # Determine canonical k-mer based on graph mode
        canonical_kmer = if graph_mode == DoubleStrand && kmer_type <: Union{Kmers.DNAKmer, Kmers.RNAKmer}
            # Use canonical representation for double-strand mode
            rc_kmer = BioSequences.reverse_complement(kmer_seq)
            kmer_seq <= rc_kmer ? kmer_seq : rc_kmer
        else
            # Use k-mer as-is for single-strand mode or amino acids
            kmer_seq
        end
        
        # Create vertex with empty coverage (will be populated if we have observations)
        graph[canonical_kmer] = KmerVertexData(canonical_kmer)
    end
    
    # Create mapping from segment IDs to k-mers
    id_to_kmer = Dict{String, kmer_type}()
    for (seg_id, sequence) in segments
        kmer_seq = kmer_type(sequence)
        canonical_kmer = if graph_mode == DoubleStrand && kmer_type <: Union{Kmers.DNAKmer, Kmers.RNAKmer}
            rc_kmer = BioSequences.reverse_complement(kmer_seq)
            kmer_seq <= rc_kmer ? kmer_seq : rc_kmer
        else
            kmer_seq
        end
        id_to_kmer[seg_id] = canonical_kmer
    end
    
    # Add edges (links)
    for (src_id, src_forward, dst_id, dst_forward) in links
        if haskey(id_to_kmer, src_id) && haskey(id_to_kmer, dst_id)
            src_kmer = id_to_kmer[src_id]
            dst_kmer = id_to_kmer[dst_id]
            
            # Convert GFA orientations to StrandOrientation
            src_strand = src_forward ? Forward : Reverse
            dst_strand = dst_forward ? Forward : Reverse
            
            # Create strand-aware edge
            graph[src_kmer, dst_kmer] = KmerEdgeData(src_strand, dst_strand)
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
$(DocStringExtensions.TYPEDSIGNATURES)

Read a GFA file and auto-detect whether to create a k-mer graph or BioSequence graph.

This function parses GFA format files and intelligently chooses between:
1. **Fixed-length k-mer graph** (if all segments have the same length)
2. **Variable-length BioSequence graph** (if segments have different lengths)

# Arguments
- `gfa_file`: Path to input GFA file
- `graph_mode`: GraphMode (SingleStrand or DoubleStrand, default: DoubleStrand)
- `force_biosequence_graph`: Force creation of variable-length BioSequence graph (default: false)

# Returns
- MetaGraphsNext.MetaGraph with either k-mer vertices or BioSequence vertices

# Auto-Detection Logic
- **Fixed-length detection**: If all segments are the same length k, creates `DNAKmer{k}`/`RNAKmer{k}`/`AAKmer{k}` graph
- **Variable-length fallback**: If segments have different lengths, creates `BioSequence` graph
- **Override**: Use `force_biosequence_graph=true` to force variable-length graph

# GFA Format Support
Supports GFA v1.0 with:
- Header (H) lines (ignored)
- Segment (S) lines: parsed as vertices (k-mer or BioSequence)
- Link (L) lines: parsed as strand-aware edges
- Path (P) lines: stored as metadata (future use)

# Examples
```julia
# Auto-detect graph type
graph = read_gfa_next("assembly.gfa")

# Force variable-length BioSequence graph
graph = read_gfa_next("assembly.gfa", force_biosequence_graph=true)

# SingleStrand mode with auto-detection
graph = read_gfa_next("assembly.gfa", SingleStrand)
```
"""
function read_gfa_next(gfa_file::AbstractString, graph_mode::GraphMode=DoubleStrand; 
                       force_biosequence_graph::Bool=false)
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
    
    # Check if all segments are the same length (fixed-length -> k-mer graph)
    segment_lengths = [length(seq) for seq in values(segments)]
    all_same_length = !isempty(segment_lengths) && all(l -> l == segment_lengths[1], segment_lengths)
    k_value = isempty(segment_lengths) ? 0 : segment_lengths[1]
    
    # Auto-detect graph type unless forced
    if !force_biosequence_graph && all_same_length && k_value > 0
        @info "Auto-detected fixed-length sequences (k=$k_value), creating k-mer graph"
        
        # Determine k-mer type from first segment
        first_seq = first(values(segments))
        if all(c -> c in "ACGT", uppercase(first_seq))
            kmer_type = Kmers.DNAKmer{k_value}
        elseif all(c -> c in "ACGU", uppercase(first_seq))
            kmer_type = Kmers.RNAKmer{k_value}
        else
            kmer_type = Kmers.AAKmer{k_value}
        end
        
        # Call the k-mer graph reader
        return read_gfa_next(gfa_file, kmer_type, graph_mode)
    end
    
    # Determine BioSequence type for variable-length graph
    biosequence_type = if isempty(segments)
        BioSequences.LongDNA{4}  # Default to DNA
    else
        first_seq = first(values(segments))
        if all(c -> c in "ACGT", uppercase(first_seq))
            BioSequences.LongDNA{4}
        elseif all(c -> c in "ACGU", uppercase(first_seq))
            BioSequences.LongRNA{4}
        else
            BioSequences.LongAA
        end
    end
    
    # Create MetaGraphsNext graph with BioSequence vertices
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=biosequence_type,
        vertex_data_type=KmerVertexData{biosequence_type},
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices (segments) - convert sequences to BioSequence types
    id_to_sequence = Dict{String, biosequence_type}()
    for (seg_id, sequence) in segments
        # Convert sequence to BioSequence type
        biosequence = biosequence_type(sequence)
        
        # Determine canonical representation based on graph mode
        canonical_seq = if graph_mode == DoubleStrand && biosequence_type <: Union{BioSequences.LongDNA, BioSequences.LongRNA}
            # Use canonical representation for double-strand mode
            rc_seq = BioSequences.reverse_complement(biosequence)
            biosequence <= rc_seq ? biosequence : rc_seq
        else
            # Use sequence as-is for single-strand mode or amino acids
            biosequence
        end
        
        # Create vertex with empty coverage (will be populated if we have observations)
        graph[canonical_seq] = KmerVertexData(canonical_seq)
        id_to_sequence[seg_id] = canonical_seq
    end
    
    # Add edges (links)
    for (src_id, src_forward, dst_id, dst_forward) in links
        if haskey(id_to_sequence, src_id) && haskey(id_to_sequence, dst_id)
            src_seq = id_to_sequence[src_id]
            dst_seq = id_to_sequence[dst_id]
            
            # Convert GFA orientations to StrandOrientation
            src_strand = src_forward ? Forward : Reverse
            dst_strand = dst_forward ? Forward : Reverse
            
            # Create strand-aware edge
            graph[src_seq, dst_seq] = KmerEdgeData(src_strand, dst_strand)
        else
            @warn "Link references unknown segment: $src_id -> $dst_id"
        end
    end
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

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
"""
Advanced Graph Algorithms for Next-Generation Assembly

This module provides sophisticated graph algorithms for assembly graph analysis,
including Eulerian path finding, bubble detection, repeat resolution, and
graph simplification operations.

All algorithms are designed to work with strand-aware MetaGraphsNext graphs
and maintain biological correctness constraints.
"""

import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import DataStructures
import Statistics

# Note: We don't export specific types/functions - use fully qualified names like Mycelia.find_eulerian_paths_next

"""
    BubbleStructure

Represents a bubble (alternative path) in the assembly graph.
"""
struct BubbleStructure
    entry_vertex::String
    exit_vertex::String
    path1::Vector{String}
    path2::Vector{String}
    path1_support::Int
    path2_support::Int
    complexity_score::Float64
    
    function BubbleStructure(entry::String, exit::String, 
                           p1::Vector{String}, p2::Vector{String},
                           s1::Int, s2::Int, complexity::Float64)
        new(entry, exit, p1, p2, s1, s2, complexity)
    end
end

"""
    RepeatRegion

Represents a repetitive region in the assembly graph.
"""
struct RepeatRegion
    repeat_vertices::Vector{String}
    incoming_edges::Vector{Tuple{String, String}}
    outgoing_edges::Vector{Tuple{String, String}}
    copy_number_estimate::Float64
    repeat_type::Symbol  # :tandem, :interspersed, :palindromic
    confidence::Float64
    
    function RepeatRegion(vertices::Vector{String}, 
                         incoming::Vector{Tuple{String, String}},
                         outgoing::Vector{Tuple{String, String}},
                         copy_num::Float64, rep_type::Symbol, conf::Float64)
        @assert rep_type in [:tandem, :interspersed, :palindromic] "Invalid repeat type"
        @assert 0.0 <= conf <= 1.0 "Confidence must be in [0,1]"
        new(vertices, incoming, outgoing, copy_num, rep_type, conf)
    end
end

"""
    ContigPath

Represents a linear path through the graph forming a contig.
"""
struct ContigPath
    vertices::Vector{String}
    sequence::String
    coverage_profile::Vector{Float64}
    length::Int
    n50_contribution::Float64
    
    function ContigPath(vertices::Vector{String}, sequence::String,
                       coverage::Vector{Float64})
        new(vertices, sequence, coverage, length(sequence), 0.0)
    end
end

"""
    ScaffoldResult

Results from scaffolding analysis.
"""
struct ScaffoldResult
    scaffolds::Vector{Vector{ContigPath}}
    gap_estimates::Vector{Tuple{Int, Int, Float64}}  # (min_gap, max_gap, confidence)
    scaffold_n50::Float64
    total_length::Int
    
    function ScaffoldResult(scaffolds::Vector{Vector{ContigPath}},
                          gaps::Vector{Tuple{Int, Int, Float64}})
        total_len = sum(sum(contig.length for contig in scaffold) for scaffold in scaffolds)
        # Calculate N50 (simplified)
        lengths = [sum(contig.length for contig in scaffold) for scaffold in scaffolds]
        sort!(lengths, rev=true)
        cumsum_lengths = cumsum(lengths)
        n50_idx = findfirst(x -> x >= total_len / 2, cumsum_lengths)
        n50 = n50_idx !== nothing ? lengths[n50_idx] : 0.0
        
        new(scaffolds, gaps, n50, total_len)
    end
end

"""
    find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph) -> Vector{Vector{String}}

Find Eulerian paths in the assembly graph. An Eulerian path visits every edge exactly once.
"""
function find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData})
    if isempty(MetaGraphsNext.labels(graph))
        return Vector{String}[]
    end
    
    # Check if Eulerian path exists
    in_degrees, out_degrees = calculate_degrees(graph)
    eulerian_info = check_eulerian_conditions(in_degrees, out_degrees)
    
    if !eulerian_info.has_path
        println("Graph does not have Eulerian path")
        return Vector{String}[]
    end
    
    # Find starting vertices
    start_vertices = find_eulerian_start_vertices(in_degrees, out_degrees, eulerian_info)
    
    paths = Vector{String}[]
    for start_vertex in start_vertices
        path = find_eulerian_path_from_vertex(graph, start_vertex, in_degrees, out_degrees)
        if !isempty(path)
            push!(paths, path)
        end
    end
    
    return paths
end

"""
Calculate in-degrees and out-degrees for all vertices.
"""
function calculate_degrees(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData})
    vertices = collect(MetaGraphsNext.labels(graph))
    in_degrees = Dict{String, Int}()
    out_degrees = Dict{String, Int}()
    
    # Initialize
    for vertex in vertices
        in_degrees[vertex] = 0
        out_degrees[vertex] = 0
    end
    
    # Count degrees
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        out_degrees[src] += 1
        in_degrees[dst] += 1
    end
    
    return in_degrees, out_degrees
end

"""
Check conditions for Eulerian path existence.
"""
function check_eulerian_conditions(in_degrees::Dict{String, Int}, out_degrees::Dict{String, Int})
    odd_vertices = String[]
    start_vertices = String[]
    end_vertices = String[]
    
    for vertex in keys(in_degrees)
        in_deg = in_degrees[vertex]
        out_deg = out_degrees[vertex]
        
        if in_deg != out_deg
            push!(odd_vertices, vertex)
            if out_deg > in_deg
                push!(start_vertices, vertex)
            elseif in_deg > out_deg
                push!(end_vertices, vertex)
            end
        end
    end
    
    # Eulerian path conditions:
    # - All vertices have equal in/out degree (Eulerian cycle), OR
    # - Exactly one vertex has out_degree = in_degree + 1 (start)
    # - Exactly one vertex has in_degree = out_degree + 1 (end)
    # - All other vertices have equal in/out degree
    
    has_path = length(odd_vertices) == 0 || 
               (length(start_vertices) == 1 && length(end_vertices) == 1)
    
    return (has_path=has_path, start_vertices=start_vertices, end_vertices=end_vertices)
end

"""
Find valid starting vertices for Eulerian paths.
"""
function find_eulerian_start_vertices(in_degrees::Dict{String, Int}, 
                                     out_degrees::Dict{String, Int},
                                     eulerian_info)
    if !isempty(eulerian_info.start_vertices)
        return eulerian_info.start_vertices
    else
        # If it's an Eulerian cycle, can start from any vertex with edges
        return [v for v in keys(out_degrees) if out_degrees[v] > 0]
    end
end

"""
Find Eulerian path starting from a specific vertex using Hierholzer's algorithm.
"""
function find_eulerian_path_from_vertex(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                       start_vertex::String,
                                       in_degrees::Dict{String, Int},
                                       out_degrees::Dict{String, Int})
    # Create adjacency list with edge tracking
    adj_list = Dict{String, Vector{Tuple{String, Bool}}}()  # (neighbor, used)
    
    for vertex in keys(out_degrees)
        adj_list[vertex] = Tuple{String, Bool}[]
    end
    
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        push!(adj_list[src], (dst, false))
    end
    
    path = String[]
    stack = [start_vertex]
    current_vertex = start_vertex
    
    while !isempty(stack) || any(any(!used for (_, used) in neighbors) for neighbors in values(adj_list))
        if any(!used for (_, used) in adj_list[current_vertex])
            push!(stack, current_vertex)
            
            # Find first unused edge
            for (i, (neighbor, used)) in enumerate(adj_list[current_vertex])
                if !used
                    adj_list[current_vertex][i] = (neighbor, true)
                    current_vertex = neighbor
                    break
                end
            end
        else
            push!(path, current_vertex)
            if !isempty(stack)
                current_vertex = pop!(stack)
            end
        end
    end
    
    reverse!(path)
    return path
end

"""
    detect_bubbles_next(graph::MetaGraphsNext.MetaGraph, min_bubble_length::Int=2, max_bubble_length::Int=100) -> Vector{BubbleStructure}

Detect bubble structures (alternative paths) in the assembly graph.
"""
function detect_bubbles_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData};
                           min_bubble_length::Int=2,
                           max_bubble_length::Int=100)
    bubbles = BubbleStructure[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    for entry_vertex in vertices
        # Find potential bubble entry points (vertices with out-degree > 1)
        out_neighbors = get_out_neighbors(graph, entry_vertex)
        
        if length(out_neighbors) >= 2
            # Look for bubbles starting from this vertex
            bubble_candidates = find_bubble_paths(graph, entry_vertex, out_neighbors, 
                                                min_bubble_length, max_bubble_length)
            
            for bubble in bubble_candidates
                if is_valid_bubble(graph, bubble)
                    push!(bubbles, bubble)
                end
            end
        end
    end
    
    return remove_duplicate_bubbles(bubbles)
end

"""
Get outgoing neighbors of a vertex.
"""
function get_out_neighbors(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData}, vertex::String)
    neighbors = String[]
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        if src == vertex
            push!(neighbors, dst)
        end
    end
    return neighbors
end

"""
Get incoming neighbors of a vertex.
"""
function get_in_neighbors(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData}, vertex::String)
    neighbors = String[]
    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        if dst == vertex
            push!(neighbors, src)
        end
    end
    return neighbors
end

"""
Find potential bubble paths from an entry vertex.
"""
function find_bubble_paths(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                          entry_vertex::String, 
                          out_neighbors::Vector{String},
                          min_length::Int, max_length::Int)
    bubbles = BubbleStructure[]
    
    # Try all pairs of outgoing paths
    for i in 1:length(out_neighbors)
        for j in (i+1):length(out_neighbors)
            path1_start = out_neighbors[i]
            path2_start = out_neighbors[j]
            
            # Find paths from each starting point
            path1 = find_limited_path(graph, path1_start, max_length)
            path2 = find_limited_path(graph, path2_start, max_length)
            
            # Check if paths reconverge
            convergence_point = find_path_convergence(path1, path2)
            
            if convergence_point !== nothing && 
               length(path1) >= min_length && length(path2) >= min_length
                
                # Extract paths up to convergence
                conv_idx1 = findfirst(v -> v == convergence_point, path1)
                conv_idx2 = findfirst(v -> v == convergence_point, path2)
                
                if conv_idx1 !== nothing && conv_idx2 !== nothing
                    bubble_path1 = path1[1:conv_idx1]
                    bubble_path2 = path2[1:conv_idx2]
                    
                    # Calculate support (simplified - could use actual coverage)
                    support1 = calculate_path_support(graph, bubble_path1)
                    support2 = calculate_path_support(graph, bubble_path2)
                    
                    complexity = calculate_bubble_complexity(bubble_path1, bubble_path2)
                    
                    bubble = BubbleStructure(entry_vertex, convergence_point,
                                           bubble_path1, bubble_path2,
                                           support1, support2, complexity)
                    push!(bubbles, bubble)
                end
            end
        end
    end
    
    return bubbles
end

"""
Find a limited-length path from a starting vertex.
"""
function find_limited_path(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                          start_vertex::String, max_length::Int)
    path = [start_vertex]
    current = start_vertex
    
    for _ in 1:max_length
        neighbors = get_out_neighbors(graph, current)
        if length(neighbors) == 1
            next_vertex = neighbors[1]
            if next_vertex in path  # Avoid cycles
                break
            end
            push!(path, next_vertex)
            current = next_vertex
        else
            break  # Multiple or no neighbors
        end
    end
    
    return path
end

"""
Find where two paths converge.
"""
function find_path_convergence(path1::Vector{String}, path2::Vector{String})
    # Find first common vertex (excluding starting vertices)
    for vertex1 in path1[2:end]
        if vertex1 in path2[2:end]
            return vertex1
        end
    end
    return nothing
end

"""
Calculate support for a path (simplified version).
"""
function calculate_path_support(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                               path::Vector{String})
    if isempty(path)
        return 0
    end
    
    # Use vertex data coverage if available
    total_coverage = 0.0
    for vertex in path
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            # Use metadata if available, otherwise default to 1
            coverage = length(vertex_data.coverage)
            total_coverage += coverage
        end
    end
    
    return round(Int, total_coverage / length(path))
end

"""
Calculate complexity score for a bubble.
"""
function calculate_bubble_complexity(path1::Vector{String}, path2::Vector{String})
    # Simple complexity metric based on path length difference and sequence similarity
    length_diff = abs(length(path1) - length(path2))
    avg_length = (length(path1) + length(path2)) / 2
    
    length_score = length_diff / max(avg_length, 1)
    
    # Could add sequence similarity scoring here
    similarity_score = 0.5  # Placeholder
    
    return length_score + (1 - similarity_score)
end

"""
Check if a bubble structure is valid.
"""
function is_valid_bubble(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                        bubble::BubbleStructure)
    # Check that paths are distinct
    if bubble.path1 == bubble.path2
        return false
    end
    
    # Check that entry and exit vertices exist
    if !(haskey(graph, bubble.entry_vertex) && haskey(graph, bubble.exit_vertex))
        return false
    end
    
    # Check that paths are valid in graph
    is_valid_path(graph, [bubble.entry_vertex; bubble.path1]) &&
    is_valid_path(graph, [bubble.entry_vertex; bubble.path2])
end

"""
Check if a path is valid in the graph.
"""
function is_valid_path(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                      path::Vector{String})
    if length(path) < 2
        return true
    end
    
    for i in 1:(length(path)-1)
        if !haskey(graph, path[i], path[i+1])
            return false
        end
    end
    return true
end

"""
Remove duplicate bubbles.
"""
function remove_duplicate_bubbles(bubbles::Vector{BubbleStructure})
    unique_bubbles = BubbleStructure[]
    
    for bubble in bubbles
        is_duplicate = false
        for existing in unique_bubbles
            if are_equivalent_bubbles(bubble, existing)
                is_duplicate = true
                break
            end
        end
        
        if !is_duplicate
            push!(unique_bubbles, bubble)
        end
    end
    
    return unique_bubbles
end

"""
Check if two bubbles are equivalent.
"""
function are_equivalent_bubbles(b1::BubbleStructure, b2::BubbleStructure)
    return (b1.entry_vertex == b2.entry_vertex && 
            b1.exit_vertex == b2.exit_vertex &&
            ((b1.path1 == b2.path1 && b1.path2 == b2.path2) ||
             (b1.path1 == b2.path2 && b1.path2 == b2.path1)))
end

"""
    resolve_repeats_next(graph::MetaGraphsNext.MetaGraph, min_repeat_length::Int=10) -> Vector{RepeatRegion}

Identify and characterize repetitive regions in the assembly graph.
"""
function resolve_repeats_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData};
                            min_repeat_length::Int=10)
    repeats = RepeatRegion[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    # Find vertices with high in-degree or out-degree (potential repeat boundaries)
    in_degrees, out_degrees = calculate_degrees(graph)
    
    # Find potential repeat vertices (high connectivity)
    repeat_candidates = find_repeat_candidates(in_degrees, out_degrees)
    
    # Analyze each candidate region
    for candidate_vertex in repeat_candidates
        repeat_region = analyze_repeat_region(graph, candidate_vertex, min_repeat_length)
        if repeat_region !== nothing
            push!(repeats, repeat_region)
        end
    end
    
    # Merge overlapping repeat regions
    merged_repeats = merge_overlapping_repeats(repeats)
    
    return merged_repeats
end

"""
Find vertices that could be part of repeat regions.
"""
function find_repeat_candidates(in_degrees::Dict{String, Int}, out_degrees::Dict{String, Int})
    candidates = String[]
    
    for vertex in keys(in_degrees)
        # High connectivity suggests repeat
        total_degree = in_degrees[vertex] + out_degrees[vertex]
        if total_degree > 4  # Threshold for repeat consideration
            push!(candidates, vertex)
        end
    end
    
    return candidates
end

"""
Analyze a potential repeat region starting from a vertex.
"""
function analyze_repeat_region(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                              start_vertex::String, min_length::Int)
    # Get local subgraph around the vertex
    local_vertices = get_local_subgraph(graph, start_vertex, min_length)
    
    if length(local_vertices) < min_length
        return nothing
    end
    
    # Analyze connectivity patterns
    incoming_edges = Tuple{String, String}[]
    outgoing_edges = Tuple{String, String}[]
    
    for vertex in local_vertices
        in_neighbors = get_in_neighbors(graph, vertex)
        out_neighbors = get_out_neighbors(graph, vertex)
        
        # Count external connections
        for neighbor in in_neighbors
            if !(neighbor in local_vertices)
                push!(incoming_edges, (neighbor, vertex))
            end
        end
        
        for neighbor in out_neighbors
            if !(neighbor in local_vertices)
                push!(outgoing_edges, (vertex, neighbor))
            end
        end
    end
    
    # Estimate copy number based on coverage
    copy_number = estimate_copy_number(graph, local_vertices)
    
    # Classify repeat type
    repeat_type = classify_repeat_type(graph, local_vertices, incoming_edges, outgoing_edges)
    
    # Calculate confidence
    confidence = calculate_repeat_confidence(graph, local_vertices, copy_number)
    
    return RepeatRegion(local_vertices, incoming_edges, outgoing_edges,
                       copy_number, repeat_type, confidence)
end

"""
Get local subgraph around a vertex.
"""
function get_local_subgraph(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                           center_vertex::String, radius::Int)
    visited = Set{String}()
    queue = DataStructures.Queue{Tuple{String, Int}}()
    DataStructures.enqueue!(queue, (center_vertex, 0))
    
    while !isempty(queue)
        vertex, distance = DataStructures.dequeue!(queue)
        
        if vertex in visited || distance > radius
            continue
        end
        
        push!(visited, vertex)
        
        # Add neighbors
        for neighbor in get_out_neighbors(graph, vertex)
            if !(neighbor in visited)
                DataStructures.enqueue!(queue, (neighbor, distance + 1))
            end
        end
        
        for neighbor in get_in_neighbors(graph, vertex)
            if !(neighbor in visited)
                DataStructures.enqueue!(queue, (neighbor, distance + 1))
            end
        end
    end
    
    return collect(visited)
end

"""
Estimate copy number for repeat region.
"""
function estimate_copy_number(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                             vertices::Vector{String})
    if isempty(vertices)
        return 1.0
    end
    
    # Use coverage information if available
    total_coverage = 0.0
    for vertex in vertices
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            coverage = length(vertex_data.coverage)
            total_coverage += coverage
        end
    end
    
    avg_coverage = total_coverage / length(vertices)
    
    # Estimate based on coverage relative to expected single-copy coverage
    expected_single_copy = 10.0  # Could be estimated from graph statistics
    copy_number = max(1.0, avg_coverage / expected_single_copy)
    
    return copy_number
end

"""
Classify the type of repeat.
"""
function classify_repeat_type(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                             vertices::Vector{String},
                             incoming_edges::Vector{Tuple{String, String}},
                             outgoing_edges::Vector{Tuple{String, String}})
    # Simple classification based on connectivity pattern
    n_incoming = length(incoming_edges)
    n_outgoing = length(outgoing_edges)
    
    if n_incoming <= 2 && n_outgoing <= 2
        return :tandem
    elseif n_incoming > 2 || n_outgoing > 2
        return :interspersed
    else
        return :palindromic
    end
end

"""
Calculate confidence in repeat identification.
"""
function calculate_repeat_confidence(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                   vertices::Vector{String}, copy_number::Float64)
    # Higher copy number and larger region size increase confidence
    size_score = min(1.0, length(vertices) / 20.0)
    copy_score = min(1.0, (copy_number - 1.0) / 5.0)
    
    return (size_score + copy_score) / 2.0
end

"""
Merge overlapping repeat regions.
"""
function merge_overlapping_repeats(repeats::Vector{RepeatRegion})
    if isempty(repeats)
        return repeats
    end
    
    merged = RepeatRegion[]
    used = falses(length(repeats))
    
    for i in 1:length(repeats)
        if used[i]
            continue
        end
        
        current_repeat = repeats[i]
        overlapping_indices = [i]
        
        # Find overlapping repeats
        for j in (i+1):length(repeats)
            if !used[j] && regions_overlap(current_repeat, repeats[j])
                push!(overlapping_indices, j)
                used[j] = true
            end
        end
        
        # Merge overlapping regions
        if length(overlapping_indices) > 1
            merged_repeat = merge_repeat_regions([repeats[idx] for idx in overlapping_indices])
            push!(merged, merged_repeat)
        else
            push!(merged, current_repeat)
        end
        
        used[i] = true
    end
    
    return merged
end

"""
Check if two repeat regions overlap.
"""
function regions_overlap(r1::RepeatRegion, r2::RepeatRegion)
    return !isempty(intersect(Set(r1.repeat_vertices), Set(r2.repeat_vertices)))
end

"""
Merge multiple repeat regions into one.
"""
function merge_repeat_regions(regions::Vector{RepeatRegion})
    all_vertices = String[]
    all_incoming = Tuple{String, String}[]
    all_outgoing = Tuple{String, String}[]
    
    for region in regions
        append!(all_vertices, region.repeat_vertices)
        append!(all_incoming, region.incoming_edges)
        append!(all_outgoing, region.outgoing_edges)
    end
    
    # Remove duplicates
    unique_vertices = unique(all_vertices)
    unique_incoming = unique(all_incoming)
    unique_outgoing = unique(all_outgoing)
    
    # Average copy number and confidence
    avg_copy_number = Statistics.mean(r.copy_number_estimate for r in regions)
    avg_confidence = Statistics.mean(r.confidence for r in regions)
    
    # Use most common repeat type
    repeat_types = [r.repeat_type for r in regions]
    repeat_type = mode(repeat_types)
    
    return RepeatRegion(unique_vertices, unique_incoming, unique_outgoing,
                       avg_copy_number, repeat_type, avg_confidence)
end

"""
    find_contigs_next(graph::MetaGraphsNext.MetaGraph, min_contig_length::Int=500) -> Vector{ContigPath}

Extract linear contigs from the assembly graph.
"""
function find_contigs_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData};
                          min_contig_length::Int=500)
    contigs = ContigPath[]
    visited = Set{String}()
    vertices = collect(MetaGraphsNext.labels(graph))
    
    for start_vertex in vertices
        if start_vertex in visited
            continue
        end
        
        # Find linear path starting from this vertex
        path = find_linear_path(graph, start_vertex, visited)
        
        if length(path) >= 2  # At least 2 k-mers for a meaningful contig
            # Generate sequence and coverage profile
            sequence = generate_contig_sequence(graph, path)
            coverage_profile = generate_coverage_profile(graph, path)
            
            if length(sequence) >= min_contig_length
                contig = ContigPath(path, sequence, coverage_profile)
                push!(contigs, contig)
            end
            
            # Mark vertices as visited
            for vertex in path
                push!(visited, vertex)
            end
        end
    end
    
    return sort_contigs_by_length(contigs)
end

"""
Find a linear path through the graph.
"""
function find_linear_path(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                         start_vertex::String, visited::Set{String})
    if start_vertex in visited
        return String[]
    end
    
    path = [start_vertex]
    current = start_vertex
    
    # Extend forward
    while true
        out_neighbors = get_out_neighbors(graph, current)
        valid_neighbors = [n for n in out_neighbors if !(n in visited) && !(n in path)]
        
        if length(valid_neighbors) == 1
            next_vertex = valid_neighbors[1]
            # Check if next vertex has only one incoming edge (to current)
            in_neighbors = get_in_neighbors(graph, next_vertex)
            if length(in_neighbors) == 1 && in_neighbors[1] == current
                push!(path, next_vertex)
                current = next_vertex
            else
                break
            end
        else
            break
        end
    end
    
    # Try to extend backward from start
    current = start_vertex
    backward_path = String[]
    
    while true
        in_neighbors = get_in_neighbors(graph, current)
        valid_neighbors = [n for n in in_neighbors if !(n in visited) && !(n in path) && !(n in backward_path)]
        
        if length(valid_neighbors) == 1
            prev_vertex = valid_neighbors[1]
            # Check if prev vertex has only one outgoing edge (to current)
            out_neighbors = get_out_neighbors(graph, prev_vertex)
            if length(out_neighbors) == 1 && out_neighbors[1] == current
                pushfirst!(backward_path, prev_vertex)
                current = prev_vertex
            else
                break
            end
        else
            break
        end
    end
    
    return [backward_path; path]
end

"""
Generate sequence for a contig path.
"""
function generate_contig_sequence(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                 path::Vector{String})
    if isempty(path)
        return ""
    end
    
    # Start with first k-mer
    sequence = path[1]
    
    # Add last character of each subsequent k-mer
    for i in 2:length(path)
        kmer = path[i]
        sequence *= kmer[end]  # Add last character
    end
    
    return sequence
end

"""
Generate coverage profile for a contig path.
"""
function generate_coverage_profile(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                  path::Vector{String})
    coverage = Float64[]
    
    for vertex in path
        if haskey(graph, vertex)
            vertex_data = graph[vertex]
            vertex_coverage = Float64(length(vertex_data.coverage))
            push!(coverage, vertex_coverage)
        else
            push!(coverage, 0.0)
        end
    end
    
    return coverage
end

"""
Sort contigs by length (descending).
"""
function sort_contigs_by_length(contigs::Vector{ContigPath})
    return sort(contigs, by=c -> c.length, rev=true)
end

"""
    simplify_graph_next(graph::MetaGraphsNext.MetaGraph, bubbles::Vector{BubbleStructure}) -> MetaGraphsNext.MetaGraph

Simplify the graph by resolving bubbles and removing low-confidence paths.
"""
function simplify_graph_next(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                           bubbles::Vector{BubbleStructure})
    # Create a copy of the graph
    simplified_graph = deepcopy(graph)
    
    # Process bubbles by confidence
    sorted_bubbles = sort(bubbles, by=b -> b.complexity_score)
    
    for bubble in sorted_bubbles
        # Choose which path to keep based on support
        if bubble.path1_support > bubble.path2_support
            remove_path_from_graph!(simplified_graph, bubble.path2, bubble.entry_vertex, bubble.exit_vertex)
        elseif bubble.path2_support > bubble.path1_support
            remove_path_from_graph!(simplified_graph, bubble.path1, bubble.entry_vertex, bubble.exit_vertex)
        else
            # Equal support - could use other criteria or keep both
            continue
        end
    end
    
    # Remove isolated vertices
    remove_isolated_vertices!(simplified_graph)
    
    return simplified_graph
end

"""
Remove a path from the graph.
"""
function remove_path_from_graph!(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData},
                                path::Vector{String}, entry_vertex::String, exit_vertex::String)
    # Remove edges in the path
    prev_vertex = entry_vertex
    for vertex in path
        if haskey(graph, prev_vertex, vertex)
            delete!(graph, prev_vertex, vertex)
        end
        prev_vertex = vertex
    end
    
    # Remove final edge to exit
    if !isempty(path) && haskey(graph, path[end], exit_vertex)
        delete!(graph, path[end], exit_vertex)
    end
    
    # Remove vertices that are now isolated
    for vertex in path
        if is_isolated_vertex(graph, vertex)
            delete!(graph, vertex)
        end
    end
end

"""
Check if a vertex is isolated (no edges).
"""
function is_isolated_vertex(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData}, vertex::String)
    return isempty(get_in_neighbors(graph, vertex)) && isempty(get_out_neighbors(graph, vertex))
end

"""
Remove all isolated vertices from the graph.
"""
function remove_isolated_vertices!(graph::MetaGraphsNext.MetaGraph{<:Integer, String, KmerVertexData, KmerEdgeData})
    vertices_to_remove = String[]
    
    for vertex in MetaGraphsNext.labels(graph)
        if is_isolated_vertex(graph, vertex)
            push!(vertices_to_remove, vertex)
        end
    end
    
    for vertex in vertices_to_remove
        delete!(graph, vertex)
    end
end
"""
Probabilistic algorithms for strand-aware k-mer graph traversal and assembly.

This module implements the core algorithms for Phase 2 of the assembly roadmap:
- Probabilistic walks with strand-consistent transitions
- Shortest probability paths (distance ∝ -log(probability))
- Maximum weight walks for high-confidence path finding
- Quality-aware scoring and path validation
"""

import MetaGraphsNext
import Graphs
import DataStructures
import Statistics
import Random
using DocStringExtensions

"""
Represents a step in a probabilistic walk through the graph.
"""
struct WalkStep
    vertex_label::String
    strand::StrandOrientation
    probability::Float64
    cumulative_probability::Float64
end

"""
Represents a complete path through the k-mer graph.
"""
struct GraphPath
    steps::Vector{WalkStep}
    total_probability::Float64
    sequence::String
    
    function GraphPath(steps::Vector{WalkStep})
        total_prob = isempty(steps) ? 0.0 : last(steps).cumulative_probability
        sequence = _reconstruct_sequence_from_path(steps)
        new(steps, total_prob, sequence)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform a probabilistic walk through the strand-aware k-mer graph.

This algorithm follows edges based on their probability weights, respecting strand
orientation constraints. The walk continues until max_steps is reached or no valid
transitions are available.

# Arguments
- `graph`: MetaGraphsNext k-mer graph with strand-aware edges
- `start_vertex`: Starting k-mer (vertex label)
- `max_steps`: Maximum number of steps to take
- `seed`: Random seed for reproducibility (optional)

# Returns
- `GraphPath`: Complete path with probability information

# Algorithm
1. Start at given vertex with forward strand orientation
2. At each step, calculate transition probabilities based on edge weights
3. Sample next vertex according to probabilities
4. Update cumulative probability and continue
5. Respect strand orientation constraints from edge metadata

# Example
```julia
graph = build_kmer_graph_next(DNAKmer{15}, observations)
path = probabilistic_walk_next(graph, "ATCGATCGATCGATC", 100)
println("Assembled sequence: \$(path.sequence)")
println("Path probability: \$(path.total_probability)")
```
"""
function probabilistic_walk_next(graph::MetaGraphsNext.MetaGraph, 
                                start_vertex::String, 
                                max_steps::Int;
                                seed::Union{Nothing, Int}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end
    
    if !(start_vertex in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Start vertex $start_vertex not found in graph"))
    end
    
    steps = Vector{WalkStep}()
    current_vertex = start_vertex
    current_strand = Forward  # Start with forward orientation
    cumulative_prob = 1.0
    
    # Add starting step
    push!(steps, WalkStep(current_vertex, current_strand, 1.0, cumulative_prob))
    
    for step in 1:max_steps
        # Get all valid outgoing edges from current vertex
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        
        if isempty(valid_transitions)
            # No valid transitions available
            break
        end
        
        # Calculate transition probabilities
        transition_probs = _calculate_transition_probabilities(valid_transitions)
        
        # Sample next transition
        next_transition = _sample_transition(valid_transitions, transition_probs)
        
        # Update path
        step_prob = next_transition[:probability]
        cumulative_prob *= step_prob
        
        current_vertex = next_transition[:target_vertex]
        current_strand = next_transition[:target_strand]
        
        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end
    
    return GraphPath(steps)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find the shortest path in probability space between two vertices.

Uses Dijkstra's algorithm where edge distances are -log(probability), so the
shortest path corresponds to the highest probability path.

# Arguments
- `graph`: MetaGraphsNext k-mer graph
- `source`: Source vertex label
- `target`: Target vertex label

# Returns
- `Union{GraphPath, Nothing}`: Shortest probability path, or nothing if no path exists

# Algorithm
1. Convert edge weights to -log(probability) distances
2. Run Dijkstra's algorithm with strand-aware edge traversal
3. Reconstruct path maintaining strand information
4. Convert back to probability space for final result
"""
function shortest_probability_path_next(graph::MetaGraphsNext.MetaGraph, 
                                       source::String, 
                                       target::String)
    if !(source in MetaGraphsNext.labels(graph)) || !(target in MetaGraphsNext.labels(graph))
        return nothing
    end
    
    # State: (vertex_label, strand_orientation)
    distances = Dict{Tuple{String, StrandOrientation}, Float64}()
    predecessors = Dict{Tuple{String, StrandOrientation}, Union{Nothing, Tuple{String, StrandOrientation}}}()
    visited = Set{Tuple{String, StrandOrientation}}()
    
    # Priority queue: (distance, vertex_label, strand)
    pq = DataStructures.PriorityQueue{Tuple{String, StrandOrientation}, Float64}()
    
    # Initialize
    start_state = (source, Forward)
    distances[start_state] = 0.0
    predecessors[start_state] = nothing
    pq[start_state] = 0.0
    
    while !isempty(pq)
        current_state = DataStructures.dequeue!(pq)
        current_vertex, current_strand = current_state
        
        if current_state in visited
            continue
        end
        push!(visited, current_state)
        
        # Check if we reached target
        if current_vertex == target
            # Reconstruct path
            return _reconstruct_shortest_path(predecessors, distances, start_state, current_state, graph)
        end
        
        # Explore neighbors
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        
        for transition in valid_transitions
            neighbor_vertex = transition[:target_vertex]
            neighbor_strand = transition[:target_strand]
            neighbor_state = (neighbor_vertex, neighbor_strand)
            
            if neighbor_state in visited
                continue
            end
            
            # Distance is -log(probability)
            edge_prob = transition[:probability]
            edge_distance = edge_prob > 0 ? -log(edge_prob) : Inf
            new_distance = distances[current_state] + edge_distance
            
            if !haskey(distances, neighbor_state) || new_distance < distances[neighbor_state]
                distances[neighbor_state] = new_distance
                predecessors[neighbor_state] = current_state
                pq[neighbor_state] = new_distance
            end
        end
    end
    
    return nothing  # No path found
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform a maximum weight walk prioritizing highest confidence edges.

This greedy algorithm always chooses the edge with the highest weight (coverage)
at each step, useful for finding high-confidence assembly paths.

# Arguments
- `graph`: MetaGraphsNext k-mer graph
- `start_vertex`: Starting vertex label
- `max_steps`: Maximum steps to take
- `weight_function`: Function to extract weight from edge data (default: uses edge.weight)

# Returns
- `GraphPath`: Path following maximum weight edges
"""
function maximum_weight_walk_next(graph::MetaGraphsNext.MetaGraph,
                                 start_vertex::String,
                                 max_steps::Int;
                                 weight_function::Function = edge_data -> edge_data.weight)
    if !(start_vertex in MetaGraphsNext.labels(graph))
        throw(ArgumentError("Start vertex $start_vertex not found in graph"))
    end
    
    steps = Vector{WalkStep}()
    current_vertex = start_vertex
    current_strand = Forward
    cumulative_prob = 1.0
    
    # Add starting step
    push!(steps, WalkStep(current_vertex, current_strand, 1.0, cumulative_prob))
    
    for step in 1:max_steps
        valid_transitions = _get_valid_transitions(graph, current_vertex, current_strand)
        
        if isempty(valid_transitions)
            break
        end
        
        # Find transition with maximum weight
        best_transition = nothing
        max_weight = -Inf
        
        for transition in valid_transitions
            weight = weight_function(transition[:edge_data])
            if weight > max_weight
                max_weight = weight
                best_transition = transition
            end
        end
        
        if best_transition === nothing
            break
        end
        
        # Follow best transition
        step_prob = best_transition[:probability]
        cumulative_prob *= step_prob
        
        current_vertex = best_transition[:target_vertex]
        current_strand = best_transition[:target_strand]
        
        push!(steps, WalkStep(current_vertex, current_strand, step_prob, cumulative_prob))
    end
    
    return GraphPath(steps)
end

"""
Helper function to get valid transitions from a vertex with given strand orientation.
"""
function _get_valid_transitions(graph, vertex_label, strand)
    transitions = []
    
    # Get all outgoing edges from this vertex
    for edge_labels in MetaGraphsNext.edge_labels(graph)
        if length(edge_labels) == 2 && edge_labels[1] == vertex_label
            target_vertex = edge_labels[2]
            edge_data = graph[edge_labels...]
            
            # Check if this edge is valid for our current strand
            if edge_data.src_strand == strand
                probability = edge_data.weight > 0 ? edge_data.weight : 1e-10
                
                push!(transitions, Dict(
                    :target_vertex => target_vertex,
                    :target_strand => edge_data.dst_strand,
                    :probability => probability,
                    :edge_data => edge_data
                ))
            end
        end
    end
    
    return transitions
end

"""
Helper function to calculate normalized transition probabilities.
"""
function _calculate_transition_probabilities(transitions)
    if isempty(transitions)
        return Float64[]
    end
    
    weights = [t[:probability] for t in transitions]
    total_weight = sum(weights)
    
    if total_weight == 0
        # Equal probability for all transitions
        return fill(1.0 / length(transitions), length(transitions))
    else
        return weights ./ total_weight
    end
end

"""
Helper function to sample a transition based on probabilities.
"""
function _sample_transition(transitions, probabilities)
    if isempty(transitions)
        return nothing
    end
    
    if length(transitions) == 1
        return first(transitions)
    end
    
    # Sample according to probabilities
    r = rand()
    cumulative = 0.0
    
    for (i, prob) in enumerate(probabilities)
        cumulative += prob
        if r <= cumulative
            return transitions[i]
        end
    end
    
    # Fallback to last transition
    return last(transitions)
end

"""
Helper function to reconstruct sequence from a graph path.
"""
function _reconstruct_sequence_from_path(steps)
    if isempty(steps)
        return ""
    end
    
    # Start with first k-mer
    first_step = first(steps)
    first_vertex_data = first_step.vertex_label  # This should be the k-mer sequence
    
    if first_step.strand == Forward
        sequence = first_vertex_data
    else
        # Reverse complement for reverse strand
        sequence = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(first_vertex_data)))
    end
    
    # Add subsequent nucleotides
    for i in 2:length(steps)
        step = steps[i]
        vertex_sequence = step.vertex_label
        
        # Get the sequence considering strand
        if step.strand == Forward
            kmer_seq = vertex_sequence
        else
            kmer_seq = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(vertex_sequence)))
        end
        
        # Add the last nucleotide (assuming k-mer overlap)
        if length(kmer_seq) > 0
            sequence *= last(kmer_seq)
        end
    end
    
    return sequence
end

"""
Helper function to reconstruct shortest path from Dijkstra's algorithm.
"""
function _reconstruct_shortest_path(predecessors, distances, start_state, end_state, graph)
    path_states = []
    current_state = end_state
    
    while current_state !== nothing
        pushfirst!(path_states, current_state)
        current_state = predecessors[current_state]
    end
    
    # Convert to WalkStep format
    steps = Vector{WalkStep}()
    total_distance = distances[end_state]
    total_probability = exp(-total_distance)
    
    for (i, (vertex, strand)) in enumerate(path_states)
        step_prob = if i == 1
            1.0
        else
            # Calculate step probability from distance difference
            prev_state = path_states[i-1]
            step_distance = distances[(vertex, strand)] - distances[prev_state]
            exp(-step_distance)
        end
        
        cumulative_prob = exp(-distances[(vertex, strand)])
        push!(steps, WalkStep(vertex, strand, step_prob, cumulative_prob))
    end
    
    return GraphPath(steps)
end

# Note: We don't export specific types/functions - use fully qualified names like Mycelia.probabilistic_walk_next
