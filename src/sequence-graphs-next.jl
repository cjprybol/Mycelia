"""
MetaGraphsNext-based Sequence Graph Implementation

This file contains the next-generation implementation of sequence graphs using MetaGraphsNext.jl
with type-stable metadata structures, replacing the deprecated MetaGraphs.jl implementation.
"""

# TODO: https://juliagraphs.org/Graphs.jl/stable/algorithms/traversals/#Graphs.eulerian

# Type-stable metadata structures for vertices and edges

# """
# Strand orientation for k-mer observations and transitions.
# - `Forward`: k-mer as observed (5' to 3')
# - `Reverse`: reverse complement of k-mer (3' to 5')
# """
# @enum StrandOrientation Forward=true Reverse=false

# """
# Graph mode for handling strand information.
# - `SingleStrand`: Sequences are single-stranded (RNA, amino acids, or directional DNA)
# - `DoubleStrand`: Sequences are double-stranded DNA/RNA with canonical representation
# """
# @enum GraphMode SingleStrand DoubleStrand

"""
Type-stable metadata for k-mer graph vertices.

Vertices always represent canonical k-mers for memory efficiency and cleaner graphs.
Strand information is tracked in the coverage data and edge transitions.

Fields:
- `coverage`: Vector of observation coverage data as (observation_id, position, strand_orientation) tuples
- `canonical_kmer`: The canonical k-mer (BioSequence type - NO string conversion)
"""
# >>> moved to src/kmer-graphs.jl (implementation commented below)
# struct KmerVertexData{KmerT}
#     coverage::Vector{Tuple{Int, Int, StrandOrientation}}  # (observation_id, position, strand_orientation)
#     canonical_kmer::KmerT  # Actual k-mer type (DNAKmer, RNAKmer, AAKmer)
#     
#     # Constructor with default empty coverage
#     KmerVertexData(canonical_kmer::KmerT) where {KmerT} = new{KmerT}(Vector{Tuple{Int, Int, StrandOrientation}}(), canonical_kmer)
# end

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
# struct KmerEdgeData
#     coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}
#     weight::Float64
#     src_strand::StrandOrientation  # Required orientation of source k-mer
#     dst_strand::StrandOrientation  # Required orientation of destination k-mer
#     
#     # Constructor with default weight calculation
#     function KmerEdgeData(coverage::Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}},
#                          src_strand::StrandOrientation, dst_strand::StrandOrientation)
#         new(coverage, Float64(length(coverage)), src_strand, dst_strand)
#     end
#     
#     # Default constructor
#     KmerEdgeData(src_strand::StrandOrientation, dst_strand::StrandOrientation) = 
#         new(Vector{Tuple{Tuple{Int, Int, StrandOrientation}, Tuple{Int, Int, StrandOrientation}}}(), 
#             0.0, src_strand, dst_strand)
# end

# Additional vertex data structures for string graphs only (others exist in their respective files)

# """
# Vertex data for string graphs.
# """
# struct StringVertexData
#     string_value::String
#     coverage::Vector{Int}
# 
#     function StringVertexData(string_value::String, coverage::Vector{Int}=Int[])
#         new(string_value, coverage)
#     end
# end
# 
# """
# Edge data for string graphs.
# """
# struct StringEdgeData
#     overlap_length::Int
#     weight::Float64
# 
#     function StringEdgeData(overlap_length::Int, weight::Float64=1.0)
#         new(overlap_length, weight)
#     end
# end

# """
# Vertex data for quality-aware string graphs.
# """
# struct QualityStringVertexData
#     string_value::String
#     quality_scores::Vector{Vector{Int}}
#     coverage::Vector{Tuple{Int, Int, StrandOrientation, Vector{Int}}}
# 
#     function QualityStringVertexData(
#         string_value::String,
#         quality_scores::Vector{Vector{Int}}=Vector{Int}[],
#         coverage::Vector{Tuple{Int, Int, StrandOrientation, Vector{Int}}}=Vector{Tuple{Int, Int, StrandOrientation, Vector{Int}}}()
#     )
#         new(string_value, quality_scores, coverage)
#     end
# end
# 
# """
# Edge data for quality-aware string graphs.
# """
# struct QualityStringEdgeData
#     overlap_length::Int
#     weight::Float64
#     quality_scores::Vector{Vector{Int}}
# 
#     function QualityStringEdgeData(overlap_length::Int, weight::Float64=1.0)
#         new(overlap_length, weight, Vector{Int}[])
#     end
# end
# >>> moved to src/kmer-graphs.jl
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
# function build_kmer_graph_next(kmer_type, observations::AbstractVector{<:Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}; 
#                                 graph_mode::GraphMode=DoubleStrand)
#     
#     # Count k-mers from observations (canonical for DoubleStrand, as-is for SingleStrand)
#     if graph_mode == DoubleStrand
#         kmer_counts = Mycelia.count_canonical_kmers(kmer_type, observations)
#     else  # SingleStrand mode
#         kmer_counts = Mycelia.count_kmers(kmer_type, observations)
#     end
#     graph_kmers = collect(keys(kmer_counts))
#     
#     if isempty(graph_kmers)
#         @warn "No k-mers found in observations"
#         return _create_empty_kmer_graph(kmer_type)
#     end
# 
#     # Infer the actual k-mer type with storage parameter from the data
#     actual_kmer_type = eltype(graph_kmers)
#     
#     # Create the MetaGraphsNext graph with type-stable metadata
#     graph = MetaGraphsNext.MetaGraph(
#         MetaGraphsNext.DiGraph(),
#        label_type=actual_kmer_type,
#         vertex_data_type=KmerVertexData{actual_kmer_type},
#         edge_data_type=KmerEdgeData,
#         weight_function=edge_data -> edge_data.weight,
#         default_weight=0.0
#     )
#     
#     # Add vertices for each k-mer (canonical in DoubleStrand, as-is in SingleStrand)
#     for kmer in graph_kmers
#         graph[kmer] = KmerVertexData(kmer)
#     end
# 
#     # Process observations to build strand-aware edges
#     for (obs_idx, observation) in enumerate(observations)
#         _add_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
#     end
#     
#     return graph
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add a sequence observation to an existing k-mer graph with strand-aware edge creation.

# Arguments
- `graph`: MetaGraphsNext k-mer graph with vertices
- `observation`: FASTA/FASTQ record
- `obs_idx`: Observation index
- `graph_kmers`: Vector of k-mers in the graph (canonical in DoubleStrand, as-is in SingleStrand)
- `graph_mode`: SingleStrand or DoubleStrand mode
"""
# function _add_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
#     observed_sequence = FASTX.sequence(observation)
#     k = length(first(graph_kmers))
#     
#     if length(observed_sequence) < k
#         observation_id = FASTX.identifier(observation)
#         @warn "Skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
#         return
#     end
#     
#     # Convert sequence to path through k-mer graph with strand information
#     observed_path = _sequence_to_canonical_path(graph_kmers, observed_sequence, graph_mode)
#     
#     if isempty(observed_path)
#         return
#     end
#     
#     # Add coverage to first vertex
#     first_canonical_kmer, first_strand = observed_path[1]
#     _add_vertex_coverage!(graph, first_canonical_kmer, obs_idx, 1, first_strand)
#     
#     # Add strand-aware edges and vertex coverage for the rest of the path
#     for i in 2:length(observed_path)
#         curr_canonical_kmer, curr_strand = observed_path[i]
#         prev_canonical_kmer, prev_strand = observed_path[i-1]
#         
#         # Add vertex coverage
#         _add_vertex_coverage!(graph, curr_canonical_kmer, obs_idx, i, curr_strand)
#         
#         # Add or update strand-aware edge
#         _add_strand_aware_edge!(graph, prev_canonical_kmer, curr_canonical_kmer, 
#                                prev_strand, curr_strand,
#                                (obs_idx, i-1, prev_strand), (obs_idx, i, curr_strand))
#     end
# end

"""
Helper function to add coverage data to a vertex.
"""
# function _add_vertex_coverage!(graph, kmer, obs_idx, position, strand_orientation)
#     vertex_data = graph[kmer]
#     new_coverage = (obs_idx, position, strand_orientation)
#     push!(vertex_data.coverage, new_coverage)
# end

"""
Helper function to add strand-aware coverage data to an edge.

This function creates edges that respect strand orientation constraints.
Each edge represents a biologically valid transition between k-mers.
"""
# function _add_strand_aware_edge!(graph, src_kmer, dst_kmer, src_strand, dst_strand, src_coverage, dst_coverage)
#     # Create edge key that includes strand information
#     edge_exists = false
#     existing_edge_data = nothing
#     
#     # Check if an edge with this strand configuration already exists
#     if haskey(graph, src_kmer, dst_kmer)
#         existing_edge_data = graph[src_kmer, dst_kmer]
#         # Check if this edge represents the same strand transition
#         if existing_edge_data.src_strand == src_strand && existing_edge_data.dst_strand == dst_strand
#             edge_exists = true
#         end
#     end
#     
#     if !edge_exists
#         # Create new strand-aware edge
#         graph[src_kmer, dst_kmer] = KmerEdgeData(src_strand, dst_strand)
#         existing_edge_data = graph[src_kmer, dst_kmer]
#     end
#     
#     # Add coverage to the edge
#     new_coverage = (src_coverage, dst_coverage)
#     push!(existing_edge_data.coverage, new_coverage)
#     
#     # Update weight based on coverage count
#     graph[src_kmer, dst_kmer] = KmerEdgeData(existing_edge_data.coverage, 
#                                             existing_edge_data.src_strand, 
#                                             existing_edge_data.dst_strand)
# end

"""
Helper function to create an empty k-mer graph.
"""
# function _create_empty_kmer_graph(kmer_type)
#     return MetaGraphsNext.MetaGraph(
#         MetaGraphsNext.DiGraph(),
#         label_type=kmer_type,
#         vertex_data_type=KmerVertexData{kmer_type},
#         edge_data_type=KmerEdgeData,
#         weight_function=edge_data -> edge_data.weight,
#         default_weight=0.0
#     )
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a sequence to a path through k-mer space with strand awareness.

This is the key function that handles the distinction between single-strand and
double-strand modes. In DoubleStrand mode, uses canonical k-mers; in SingleStrand mode, uses k-mers as-is.

# Arguments
- `graph_kmers`: Vector of k-mers available in the graph (canonical in DoubleStrand, as-is in SingleStrand)
- `sequence`: DNA/RNA sequence to convert
- `graph_mode`: SingleStrand or DoubleStrand mode

# Returns
- Vector of (graph_kmer, strand_orientation) pairs representing the path

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
# function _sequence_to_canonical_path(graph_kmers, sequence, graph_mode)
#     k = length(first(graph_kmers))
#     path = Vector{Tuple{eltype(graph_kmers), StrandOrientation}}()
#     graph_kmer_set = Set(graph_kmers)
#     
#     # Extract k-mers from sequence
#     for i in 1:(length(sequence) - k + 1)
#         subseq = sequence[i:i+k-1]
#         observed_kmer = typeof(first(graph_kmers))(subseq)
#
#         if graph_mode == SingleStrand
#             # Single-strand mode: use k-mers as-is
#             if observed_kmer in graph_kmer_set
#                 push!(path, (observed_kmer, Forward))
#             else
#                 @warn "K-mer $observed_kmer not found in graph at position $i (SingleStrand mode)"
#             end
#         else  # DoubleStrand mode
#             # Double-strand mode: find canonical representation
#             # Only valid for nucleic acids, not amino acids
#             if typeof(observed_kmer) <: Kmers.Kmer{<:BioSequences.NucleicAcidAlphabet}
#                 rc_kmer = BioSequences.reverse_complement(observed_kmer)
#             else
#                 # For amino acids, treat as single-strand even in DoubleStrand mode
#                 if observed_kmer in graph_kmer_set
#                     push!(path, (observed_kmer, Forward))
#                     continue
#                 else
#                     @warn "K-mer $observed_kmer not found in graph at position $i (amino acid in DoubleStrand mode)"
#                     continue
#                 end
#             end
#
#             # Determine canonical k-mer using BioSequences.canonical() to match graph construction
#             canonical_kmer = BioSequences.canonical(observed_kmer)
#             strand_orientation = BioSequences.iscanonical(observed_kmer) ? Forward : Reverse
#
#             if canonical_kmer in graph_kmer_set
#                 push!(path, (canonical_kmer, strand_orientation))
#             else
#                 @warn "Canonical k-mer $canonical_kmer not found in graph at position $i (DoubleStrand mode)"
#             end
#         end
#     end
#     
#     return path
# end

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
# function _is_valid_transition(src_kmer, dst_kmer, src_strand, dst_strand, k)
#     # Get the actual k-mer sequences considering strand orientation
#     src_seq = src_strand == Forward ? src_kmer : BioSequences.reverse_complement(src_kmer)
#     dst_seq = dst_strand == Forward ? dst_kmer : BioSequences.reverse_complement(dst_kmer)
#     
#     # Check if suffix of src matches prefix of dst (k-1 overlap)
#     # For BioSequences, we can use slicing directly
#     src_suffix = src_seq[2:end]  # Remove first nucleotide
#     dst_prefix = dst_seq[1:end-1]  # Remove last nucleotide
#     
#     return src_suffix == dst_prefix
# end

# Compatibility layer for migrating from legacy MetaGraphs implementation


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
# function write_gfa_next(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString)
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0\tMY:Z:Mycelia-Next")
#         vertex_id_map = Dict()
#         for (i, label) in enumerate(MetaGraphsNext.labels(graph))
#             vertex_id_map[label] = i
#             vertex_data = graph[label]
#             depth = length(vertex_data.coverage)
#             sequence = if hasfield(typeof(vertex_data), :canonical_kmer)
#                 string(vertex_data.canonical_kmer)
#             elseif hasfield(typeof(vertex_data), :canonical_qualmer)
#                 string(vertex_data.canonical_qualmer.kmer)
#             elseif hasfield(typeof(vertex_data), :sequence)
#                 string(vertex_data.sequence)
#             else
#                 error("Unknown vertex data type: $(typeof(vertex_data))")
#             end
#             println(io, "S\t$i\t$(sequence)\tDP:f:$depth")
#         end
# 
#         overlap = 0
#         if !isempty(MetaGraphsNext.labels(graph))
#             first_label = first(MetaGraphsNext.labels(graph))
#             first_vertex_data = graph[first_label]
#             if hasfield(typeof(first_vertex_data), :canonical_kmer)
#                 overlap = length(first_vertex_data.canonical_kmer) - 1
#             elseif hasfield(typeof(first_vertex_data), :canonical_qualmer)
#                 overlap = length(first_vertex_data.canonical_qualmer.kmer) - 1
#             elseif hasfield(typeof(first_vertex_data), :sequence) && !isempty(first_vertex_data.constituent_kmers)
#                 overlap = length(first_vertex_data.constituent_kmers[1]) - 1
#             end
#         end
# 
#         for edge_labels in MetaGraphsNext.edge_labels(graph)
#             if length(edge_labels) == 2
#                 src_label, dst_label = edge_labels
#                 edge_data = graph[src_label, dst_label]
#                 src_id = vertex_id_map[src_label]
#                 dst_id = vertex_id_map[dst_label]
#                 src_orientation = edge_data.src_strand == Forward ? '+' : '-'
#                 dst_orientation = edge_data.dst_strand == Forward ? '+' : '-'
#                 println(io, "L\t$src_id\t$src_orientation\t$dst_id\t$dst_orientation\t$(overlap)M")
#             end
#         end
#     end
# 
#     return outfile
# end

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
function read_gfa_next(
    gfa_file::AbstractString,
    kmer_type::Type,
    graph_mode::Union{GraphMode,Rhizomorph.GraphMode}=DoubleStrand,
)
    # Normalize legacy GraphMode to Rhizomorph.GraphMode
    rhizo_mode = graph_mode isa Rhizomorph.GraphMode ? graph_mode :
        Rhizomorph.GraphMode(Int(graph_mode))
    return Rhizomorph.read_gfa_next(gfa_file, kmer_type, rhizo_mode)
end

"""
function read_gfa_next(
    gfa_file::AbstractString,
    graph_mode::Rhizomorph.GraphMode;
    force_biosequence_graph::Bool=false,
)
    return Rhizomorph.read_gfa_next(
        gfa_file,
        graph_mode;
        force_biosequence_graph=force_biosequence_graph,
    )
end

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
    first_seq = first(values(segments))
    seq_type = detect_alphabet(first_seq)
    
    # Auto-detect graph type unless forced
    if !force_biosequence_graph && all_same_length && k_value > 0
        @info "Auto-detected fixed-length sequences (k=$k_value), creating k-mer graph"
        # Determine k-mer type from first segment
        if seq_type == :DNA
            sample_bio_seq = BioSequences.LongDNA{4}(first_seq)
            sample_kmer = Kmers.DNAKmer{k_value}(sample_bio_seq)
            kmer_type = typeof(sample_kmer)
        elseif seq_type == :RNA
            sample_bio_seq = BioSequences.LongRNA{4}(first_seq)
            sample_kmer = Kmers.RNAKmer{k_value}(sample_bio_seq)
            kmer_type = typeof(sample_kmer)
        else
            @assert seq_type == :AA "Unknown sequence type for k-mer graph"
            sample_bio_seq = BioSequences.LongAA(first_seq)
            sample_kmer = Kmers.AAKmer{k_value}(sample_bio_seq)
            kmer_type = typeof(sample_kmer)
        end
        
        # Call the k-mer graph reader
        return read_gfa_next(gfa_file, kmer_type, graph_mode)
    end
    
    # Determine BioSequence type for variable-length graph
    biosequence_type = if isempty(segments)
        BioSequences.LongDNA{4}  # Default to DNA
    else
        if seq_type == :DNA
            BioSequences.LongDNA{4}
        elseif seq_type == :RNA
            BioSequences.LongRNA{4}
        else
            @assert seq_type == :AA "Unknown sequence type for BioSequence graph"
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
        
        # Preserve sequence orientation as observed in the GFA.
        graph[biosequence] = KmerVertexData(biosequence)
        id_to_sequence[seg_id] = biosequence
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
Calculate in-degrees and out-degrees for all vertices.
"""
function calculate_degrees(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any}) where T
    vertices = collect(MetaGraphsNext.labels(graph))
    in_degrees = Dict{T, Int}()
    out_degrees = Dict{T, Int}()

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
function check_eulerian_conditions(in_degrees::Dict{T, Int}, out_degrees::Dict{T, Int}) where T
    odd_vertices = T[]
    start_vertices = T[]
    end_vertices = T[]

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
function find_eulerian_start_vertices(in_degrees::Dict{T, Int},
                                     out_degrees::Dict{T, Int},
                                     eulerian_info) where T
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
function find_eulerian_path_from_vertex(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any},
                                       start_vertex::T,
                                       in_degrees::Dict{T, Int},
                                       out_degrees::Dict{T, Int}) where T
    # Create adjacency list with edge tracking
    adj_list = Dict{T, Vector{Tuple{T, Bool}}}()  # (neighbor, used)

    for vertex in keys(out_degrees)
        adj_list[vertex] = Tuple{T, Bool}[]
    end

    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        push!(adj_list[src], (dst, false))
    end

    path = T[]
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
const RHIZOMORPH_VERTEX_TYPES = Union{
    Rhizomorph.KmerVertexData,
    Rhizomorph.QualmerVertexData,
    Rhizomorph.BioSequenceVertexData,
    Rhizomorph.QualityBioSequenceVertexData,
    Rhizomorph.StringVertexData,
    Rhizomorph.QualityStringVertexData,
}

const RHIZOMORPH_EDGE_TYPES = Union{
    Rhizomorph.KmerEdgeData,
    Rhizomorph.QualmerEdgeData,
    Rhizomorph.BioSequenceEdgeData,
    Rhizomorph.QualityBioSequenceEdgeData,
    Rhizomorph.StringEdgeData,
    Rhizomorph.QualityStringEdgeData,
}

# Delegate GFA export for Rhizomorph graphs to the evidence-aware implementation.
function write_gfa_next(
    graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, <:Any, V, E},
    outfile::AbstractString,
) where {V<:RHIZOMORPH_VERTEX_TYPES, E<:RHIZOMORPH_EDGE_TYPES}
    return Rhizomorph.write_gfa_next(graph, outfile)
end

# Fallback for Rhizomorph graphs that use evidence-based metadata rather than coverage vectors.
function detect_bubbles_next(
    graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, <:Any, V, E};
    min_bubble_length::Int=2,
    max_bubble_length::Int=100,
) where {V<:RHIZOMORPH_VERTEX_TYPES, E<:RHIZOMORPH_EDGE_TYPES}
    return Rhizomorph.detect_bubbles_next(
        graph;
        min_bubble_length=min_bubble_length,
        max_bubble_length=max_bubble_length,
    )
end

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
Type parameter T allows for different vertex label types (k-mers, strings, biosequences, etc.)
"""
struct WalkStep{T}
    vertex_label::T
    strand::StrandOrientation
    probability::Float64
    cumulative_probability::Float64
end

"""
Represents a complete path through the graph.
Type parameter T allows for different vertex label types (k-mers, strings, biosequences, etc.)
No pre-computed sequence - use path_to_sequence() to reconstruct sequences from paths.
"""
struct GraphPath{T}
    steps::Vector{WalkStep{T}}
    total_probability::Float64

    function GraphPath{T}(steps::Vector{WalkStep{T}}) where T
        total_prob = isempty(steps) ? 0.0 : last(steps).cumulative_probability
        new{T}(steps, total_prob)
    end
end

# Convenience constructor that infers the type parameter
function GraphPath(steps::Vector{WalkStep{T}}) where T
    return GraphPath{T}(steps)
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform probabilistic assembly with confidence intervals for each base.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Assembly graph with k-mer vertices
- `start_vertices::Vector{String}`: Starting vertices for assembly paths
- `confidence_threshold::Float64`: Minimum confidence threshold for base calls (default: 0.95)
- `max_paths::Int`: Maximum number of paths to consider (default: 1000)

# Returns
Named tuple containing:
- `assembly_sequences::Vector{BioSequences.LongDNA{4}}`: Assembled sequences
- `confidence_intervals::Vector{Vector{Float64}}`: Per-base confidence scores
- `alternative_bases::Vector{Dict{Int, Vector{BioSequences.DNA}}}`: Alternative bases at low-confidence positions

# Details
- Implements maximum likelihood path selection through assembly graphs
- Calculates per-base confidence intervals based on path probabilities
- Identifies ambiguous regions where multiple paths have similar likelihoods
- Supports probabilistic variant calling from graph ambiguities
"""
function probabilistic_assembly_with_confidence(graph::MetaGraphsNext.MetaGraph, 
                                               start_vertices::Vector{String}; 
                                               confidence_threshold::Float64=0.95, 
                                               max_paths::Int=1000)
    
    assembly_sequences = Vector{BioSequences.LongDNA{4}}()
    confidence_intervals = Vector{Vector{Float64}}()
    alternative_bases = Vector{Dict{Int, Vector{BioSequences.DNA}}}()
    
    for start_vertex in start_vertices
        @info "Processing assembly from vertex: $start_vertex"
        
        # Generate multiple probabilistic paths
        paths = Vector{GraphPath}()
        path_probabilities = Vector{Float64}()
        
        for i in 1:max_paths
            try
                path = probabilistic_walk_next(graph, start_vertex, 10000; seed=i)
                path_prob = path.steps[end].cumulative_probability
                
                if path_prob > 1e-10  # Filter out very low probability paths
                    push!(paths, path)
                    push!(path_probabilities, path_prob)
                end
            catch e
                # Handle dead ends or other path generation issues
                break
            end
        end
        
        if isempty(paths)
            @warn "No valid paths found from vertex $start_vertex"
            continue
        end
        
        # Normalize path probabilities
        path_probabilities = path_probabilities ./ sum(path_probabilities)
        
        # Convert paths to sequences and calculate confidence
        path_sequences = Vector{BioSequences.LongDNA{4}}()
        for path in paths
            seq = path_to_sequence(path, graph)
            push!(path_sequences, seq)
        end
        
        if isempty(path_sequences)
            continue
        end
        
        # Find consensus sequence and confidence intervals
        consensus_seq, confidence_scores, alt_bases = calculate_consensus_with_confidence(
            path_sequences, path_probabilities, confidence_threshold
        )
        
        push!(assembly_sequences, consensus_seq)
        push!(confidence_intervals, confidence_scores)
        push!(alternative_bases, alt_bases)
    end
    
    return (;assembly_sequences, confidence_intervals, alternative_bases)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate consensus sequence with per-base confidence intervals.

# Arguments
- `sequences::Vector{BioSequences.LongDNA{4}}`: Input sequences from different paths
- `weights::Vector{Float64}`: Probability weights for each sequence
- `confidence_threshold::Float64`: Minimum confidence threshold

# Returns
Tuple containing consensus sequence, confidence scores, and alternative bases
"""
function calculate_consensus_with_confidence(sequences::Vector{BioSequences.LongDNA{4}}, 
                                           weights::Vector{Float64}, 
                                           confidence_threshold::Float64)
    
    if isempty(sequences)
        return BioSequences.LongDNA{4}(), Float64[], Dict{Int, Vector{BioSequences.DNA}}()
    end
    
    # Find maximum sequence length
    max_length = maximum(length.(sequences))
    
    # Initialize base count matrices
    base_counts = Dict{BioSequences.DNA, Vector{Float64}}()
    for base in [BioSequences.DNA_A, BioSequences.DNA_C, BioSequences.DNA_G, BioSequences.DNA_T]
        base_counts[base] = zeros(Float64, max_length)
    end
    
    # Count weighted bases at each position
    for (seq_idx, seq) in enumerate(sequences)
        weight = weights[seq_idx]
        for (pos, base) in enumerate(seq)
            if base in keys(base_counts) && pos <= max_length
                base_counts[base][pos] += weight
            end
        end
    end
    
    # Build consensus sequence and calculate confidence
    consensus_bases = Vector{BioSequences.DNA}()
    confidence_scores = Vector{Float64}()
    alternative_bases = Dict{Int, Vector{BioSequences.DNA}}()
    
    for pos in 1:max_length
        # Find most frequent base at this position
        position_counts = Dict(base => counts[pos] for (base, counts) in base_counts)
        
        # Filter out zero counts
        position_counts = filter(x -> x.second > 0, position_counts)
        
        if isempty(position_counts)
            # No coverage at this position
            push!(consensus_bases, BioSequences.DNA_N)
            push!(confidence_scores, 0.0)
            continue
        end
        
        # Sort bases by frequency
        sorted_bases = sort(collect(position_counts), by=x->x[2], rev=true)
        consensus_base = sorted_bases[1][1]
        consensus_count = sorted_bases[1][2]
        
        # Calculate confidence as proportion of consensus base
        total_count = sum(values(position_counts))
        confidence = consensus_count / total_count
        
        push!(consensus_bases, consensus_base)
        push!(confidence_scores, confidence)
        
        # Record alternative bases for low-confidence positions
        if confidence < confidence_threshold && length(sorted_bases) > 1
            alt_bases = [base for (base, count) in sorted_bases[2:end] if count / total_count > 0.1]
            if !isempty(alt_bases)
                alternative_bases[pos] = alt_bases
            end
        end
    end
    
    consensus_sequence = BioSequences.LongDNA{4}(consensus_bases)
    return consensus_sequence, confidence_scores, alternative_bases
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a graph path to a DNA sequence.

# Arguments
- `path::GraphPath`: Path through the assembly graph
- `graph::MetaGraphsNext.MetaGraph`: Assembly graph

# Returns
Appropriate sequence type based on the graph content (maintains type stability)
"""
function path_to_sequence(path::GraphPath{T}, graph::MetaGraphsNext.MetaGraph) where T
    if isempty(path.steps)
        # Return appropriate empty sequence type
        return _get_empty_sequence_for_label_type(T)
    end

    # Determine the sequence type from the first vertex
    first_step = path.steps[1]
    SequenceType = _get_sequence_type_for_path(first_step, graph)

    # Collect sequence parts maintaining proper types
    sequence_parts = Vector{SequenceType}()

    for (i, step) in enumerate(path.steps)
        if i == 1
            # First k-mer: add the full sequence
            seq_part = _extract_sequence_from_step(step, graph, SequenceType)
            push!(sequence_parts, seq_part)
        else
            # Subsequent k-mers: add only the last symbol (overlap by k-1)
            seq_part = _extract_sequence_from_step(step, graph, SequenceType)
            if length(seq_part) > 0
                last_symbol = seq_part[end:end]
                push!(sequence_parts, last_symbol)
            end
        end
    end

    # Concatenate all parts
    if isempty(sequence_parts)
        return _get_empty_sequence_for_label_type(T)
    end

    result = sequence_parts[1]
    for part in sequence_parts[2:end]
        result = result * part
    end

    return result
end

"""
Extract sequence from a single step, maintaining type stability.
"""
function _extract_sequence_from_step(step::WalkStep{T}, graph, SequenceType) where T
    vertex_data = graph[step.vertex_label]

    # Handle different vertex data types
    if hasfield(typeof(vertex_data), :sequence)
        # BioSequence graph: vertex data contains the sequence directly
        return vertex_data.sequence
    elseif hasfield(typeof(vertex_data), :canonical_kmer)
        # K-mer graph: extract k-mer and convert to appropriate sequence type
        kmer = vertex_data.canonical_kmer
        return _convert_kmer_to_sequence(kmer, SequenceType, step.strand)
    else
        # Direct k-mer labels: vertex label IS the k-mer
        return _convert_kmer_to_sequence(step.vertex_label, SequenceType, step.strand)
    end
end

"""
Convert k-mer to appropriate sequence type, handling strand orientation.
"""
function _convert_kmer_to_sequence(kmer, SequenceType, strand::StrandOrientation)
    # Convert k-mer to sequence
    if SequenceType <: BioSequences.LongSequence
        seq = SequenceType(string(kmer))
        # Handle reverse strand for biological sequences
        if strand == Reverse
            try
                seq = BioSequences.reverse_complement(seq)
            catch
                # If reverse complement fails (e.g., for AA), keep original
            end
        end
        return seq
    else
        # For string types or others, convert to string
        kmer_str = string(kmer)
        return strand == Forward ? kmer_str : kmer_str  # No rev comp for strings
    end
end

"""
Determine the appropriate sequence type for reconstruction from the first step.
"""
function _get_sequence_type_for_path(step::WalkStep{T}, graph) where T
    vertex_data = graph[step.vertex_label]

    if hasfield(typeof(vertex_data), :sequence)
        # BioSequence graph: use the sequence type from vertex data
        return typeof(vertex_data.sequence)
    elseif hasfield(typeof(vertex_data), :canonical_kmer)
        # K-mer graph: determine sequence type from k-mer
        kmer = vertex_data.canonical_kmer
        return _sequence_type_from_kmer_type(typeof(kmer))
    else
        # Direct k-mer labels
        return _sequence_type_from_kmer_type(T)
    end
end

"""
Map k-mer type to corresponding sequence type.
"""
function _sequence_type_from_kmer_type(kmer_type::Type)
    if kmer_type <: AbstractString
        return String
    end

    # Extract alphabet from k-mer type parameters
    if length(kmer_type.parameters) >= 1
        alphabet_type = kmer_type.parameters[1]
        if alphabet_type <: BioSequences.DNAAlphabet
            return BioSequences.LongDNA{4}
        elseif alphabet_type <: BioSequences.RNAAlphabet
            return BioSequences.LongRNA{4}
        elseif alphabet_type <: BioSequences.AminoAcidAlphabet
            return BioSequences.LongAA
        end
    end

    # Default fallback
    return BioSequences.LongDNA{4}
end

"""
Get appropriate empty sequence for the given label type.
"""
function _get_empty_sequence_for_label_type(::Type{T}) where T
    if T <: AbstractString
        return ""
    else
        SequenceType = _sequence_type_from_kmer_type(T)
        return SequenceType()
    end
end

"""
    find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, <:Any, <:Any})

Efficient Eulerian path finder that works directly with the underlying Graphs.jl structure.
No temporary graph creation - works with vertex indices and translates back to labels.
"""
function find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, T, <:Any, <:Any, <:Any, <:Any, <:Any}) where T
    underlying_graph = graph.graph  # Direct access to underlying Graphs.jl structure

    if Graphs.nv(underlying_graph) == 0
        return Vector{Vector{T}}()
    end

    # Work with vertex indices directly for efficiency
    vertex_indices = Graphs.vertices(underlying_graph)

    # Calculate in-degrees and out-degrees using Graphs.jl functions
    in_degrees = Dict{Int, Int}()
    out_degrees = Dict{Int, Int}()

    for v in vertex_indices
        in_degrees[v] = Graphs.indegree(underlying_graph, v)
        out_degrees[v] = Graphs.outdegree(underlying_graph, v)
    end

    # Check Eulerian path conditions
    start_vertices = Int[]
    end_vertices = Int[]

    for v in vertex_indices
        diff = out_degrees[v] - in_degrees[v]
        if diff == 1
            push!(start_vertices, v)
        elseif diff == -1
            push!(end_vertices, v)
        elseif diff != 0
            # Graph doesn't have Eulerian path
            return Vector{Vector{T}}()
        end
    end

    # Validate Eulerian conditions
    if length(start_vertices) > 1 || length(end_vertices) > 1
        return Vector{Vector{T}}()
    end
    if length(start_vertices) != length(end_vertices)
        return Vector{Vector{T}}()
    end

    # Find starting vertex
    start_vertex = if !isempty(start_vertices)
        start_vertices[1]
    else
        # Eulerian cycle - start from any vertex
        first(vertex_indices)
    end

    # Hierholzer's algorithm on underlying graph
    path = _find_eulerian_path_hierholzer(underlying_graph, start_vertex)

    if isempty(path)
        return Vector{Vector{T}}()
    end

    # Convert path from indices to labels
    label_path = [MetaGraphsNext.label_for(graph, idx) for idx in path]

    return [label_path]
end

"""
Hierholzer's algorithm implementation working directly on Graphs.jl structure.
"""
function _find_eulerian_path_hierholzer(graph::Graphs.AbstractGraph, start_vertex::Int)
    # Create adjacency list with edge usage tracking
    adj_list = Dict{Int, Vector{Tuple{Int, Bool}}}()  # vertex -> [(neighbor, used)]

    for v in Graphs.vertices(graph)
        adj_list[v] = [(dst, false) for dst in Graphs.outneighbors(graph, v)]
    end

    circuit = Int[start_vertex]
    current_vertex = start_vertex

    while !isempty(circuit)
        if any(!used for (_, used) in adj_list[current_vertex])
            # Find unused edge
            for (i, (neighbor, used)) in enumerate(adj_list[current_vertex])
                if !used
                    # Mark edge as used
                    adj_list[current_vertex][i] = (neighbor, true)
                    push!(circuit, neighbor)
                    current_vertex = neighbor
                    break
                end
            end
        else
            # No unused edges - backtrack
            if length(circuit) == 1
                break
            end
            current_vertex = circuit[end-1]

            # Check if we can extend from any vertex in current path
            extended = false
            for (i, v) in enumerate(circuit[end:-1:1])
                if any(!used for (_, used) in adj_list[v])
                    # Found vertex with unused edges - create subcircuit
                    subcircuit_start = length(circuit) - i + 1
                    temp_circuit = Int[v]
                    temp_current = v

                    while any(!used for (_, used) in adj_list[temp_current])
                        for (j, (neighbor, used)) in enumerate(adj_list[temp_current])
                            if !used
                                adj_list[temp_current][j] = (neighbor, true)
                                push!(temp_circuit, neighbor)
                                temp_current = neighbor
                                break
                            end
                        end
                    end

                    # Insert subcircuit into main circuit
                    if length(temp_circuit) > 1
                        splice!(circuit, subcircuit_start:subcircuit_start, temp_circuit[2:end])
                        current_vertex = temp_circuit[end]
                        extended = true
                        break
                    end
                end
            end

            if !extended
                break
            end
        end
    end

    # Verify all edges were used
    total_edges_used = 0
    for v in keys(adj_list)
        total_edges_used += count(used for (_, used) in adj_list[v])
    end

    if total_edges_used != Graphs.ne(graph)
        return Int[]  # Failed to find Eulerian path
    end

    return circuit
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a BioSequence graph directly from FASTA records using MetaGraphsNext.

This creates variable-length BioSequence vertices without going through k-mer decomposition first.
"""
function build_biosequence_graph_next(sequence_type::Type{<:BioSequences.BioSequence},
                                     fasta_records::Vector{FASTX.FASTA.Record};
                                     graph_mode::GraphMode=DoubleStrand,
                                     min_overlap::Int=100)

    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end

    # Create MetaGraphsNext graph
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=sequence_type,
        vertex_data_type=BioSequenceVertexData{sequence_type},
        edge_data_type=BioSequenceEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )

    # Add vertices for each sequence
    for record in fasta_records
        sequence = FASTX.sequence(sequence_type, record)

        # For double-strand mode, use canonical representation
        canonical_sequence = if graph_mode == DoubleStrand && sequence isa Union{BioSequences.LongDNA, BioSequences.LongRNA}
            rc_sequence = BioSequences.reverse_complement(sequence)
            sequence <= rc_sequence ? sequence : rc_sequence
        else
            sequence
        end

        if !haskey(graph, canonical_sequence)
            graph[canonical_sequence] = BioSequenceVertexData(canonical_sequence)
        end
    end

    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a BioSequence graph for strings using MetaGraphsNext.
"""
function build_biosequence_graph_next(sequence_type::Type{String},
                                     fasta_records::Vector{FASTX.FASTA.Record};
                                     graph_mode::GraphMode=SingleStrand,
                                     min_overlap::Int=100)

    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end

    # Create MetaGraphsNext graph for strings
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=StringVertexData,
        edge_data_type=StringEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )

    # Add vertices for each sequence
    for record in fasta_records
        sequence_str = FASTX.sequence(String, record)

        if !haskey(graph, sequence_str)
            graph[sequence_str] = StringVertexData(sequence_str, Int[])
        end
    end

    return graph
end

## >>> moved to src/qualmer-graphs.jl
#= begin_qualmer_graph_next
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a qualmer graph from FASTQ records with quality scores.
"""
function build_qualmer_graph_next(kmer_type, fastq_records::Vector{FASTX.FASTQ.Record};
                                 graph_mode::GraphMode=DoubleStrand)
    # Count k-mers from FASTQ observations with quality integration
    if graph_mode == DoubleStrand
        kmer_counts = Mycelia.count_canonical_qualmers(kmer_type, fastq_records)
    else  # SingleStrand mode
        kmer_counts = Mycelia.count_qualmers(kmer_type, fastq_records)
    end
    graph_kmers = collect(keys(kmer_counts))

    if isempty(graph_kmers)
        @warn "No k-mers found in FASTQ observations"
        return _create_empty_qualmer_graph(kmer_type)
    end

    actual_kmer_type = eltype(graph_kmers)

    # Create the MetaGraphsNext graph with qualmer metadata
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=actual_kmer_type,
        vertex_data_type=QualmerVertexData{actual_kmer_type},
        edge_data_type=QualmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )

    # Add vertices for each k-mer with quality information
    for kmer in graph_kmers
        graph[kmer] = QualmerVertexData(kmer)
    end

    # Process FASTQ observations to build quality-aware edges
    for (obs_idx, observation) in enumerate(fastq_records)
        _add_qualmer_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
    end

    return graph
end
=#  # end begin_qualmer_graph_next

#= begin_qualmer_biosequence_graph_next
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a quality-aware BioSequence graph from FASTQ records.
"""
# function build_qualmer_biosequence_graph_next(sequence_type::Type{<:BioSequences.BioSequence},
#                                              fastq_records::Vector{FASTX.FASTQ.Record};
#                                              graph_mode::GraphMode=DoubleStrand)

    if isempty(fastq_records)
        throw(ArgumentError("Cannot build graph from empty FASTQ records"))
    end

    # Create MetaGraphsNext graph with quality-aware BioSequence data
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=sequence_type,
        vertex_data_type=QualityBioSequenceVertexData{sequence_type},
        edge_data_type=QualityBioSequenceEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )

    # Add vertices for each sequence with quality information
    for record in fastq_records
        sequence = FASTX.sequence(sequence_type, record)
        quality = FASTX.quality(record)

        # For double-strand mode, use canonical representation
        canonical_sequence = if graph_mode == DoubleStrand && sequence isa Union{BioSequences.LongDNA, BioSequences.LongRNA}
            rc_sequence = BioSequences.reverse_complement(sequence)
            sequence <= rc_sequence ? sequence : rc_sequence
        else
            sequence
        end

        if !haskey(graph, canonical_sequence)
            graph[canonical_sequence] = QualityBioSequenceVertexData(canonical_sequence, [quality])
        else
            # Add quality information to existing vertex
            vertex_data = graph[canonical_sequence]
            push!(vertex_data.quality_scores, quality)
        end
    end

    return graph
# end
=#

#= begin_string_graph_next
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a string n-gram graph using MetaGraphsNext.
"""
# function build_string_graph_next(fasta_records::Vector{FASTX.FASTA.Record}, ngram_length::Int;
#                                 graph_mode::GraphMode=SingleStrand)

    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end

    # Create MetaGraphsNext graph for string n-grams
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=StringVertexData,
        edge_data_type=StringEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )

    # Extract n-grams from each string sequence
    ngram_counts = Dict{String, Int}()
    for record in fasta_records
        sequence_str = FASTX.sequence(String, record)

        # Generate n-grams
        for i in 1:(length(sequence_str) - ngram_length + 1)
            ngram = sequence_str[i:i+ngram_length-1]
            ngram_counts[ngram] = get(ngram_counts, ngram, 0) + 1
        end
    end

    # Add vertices for each n-gram
    for ngram in keys(ngram_counts)
        graph[ngram] = StringVertexData(ngram, Int[])
    end

    # Add edges between consecutive n-grams
    for record in fasta_records
        sequence_str = FASTX.sequence(String, record)
        ngrams = [sequence_str[i:i+ngram_length-1] for i in 1:(length(sequence_str) - ngram_length + 1)]

        for i in 1:(length(ngrams) - 1)
            curr_ngram = ngrams[i]
            next_ngram = ngrams[i + 1]

            if haskey(graph, curr_ngram) && haskey(graph, next_ngram)
                # Create edge
                if !MetaGraphsNext.has_edge(graph, curr_ngram, next_ngram)
                    graph[curr_ngram, next_ngram] = StringEdgeData(ngram_length - 1, 1.0)
                end
            end
        end
    end

    return graph
# end
=#  # end begin_string_graph_next

#= begin_qualmer_string_graph_next
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a quality-aware string graph from FASTQ records.
"""
# function build_qualmer_string_graph_next(fastq_records::Vector{FASTX.FASTQ.Record}, ngram_length::Int;
#                                         graph_mode::GraphMode=SingleStrand)

    if isempty(fastq_records)
        throw(ArgumentError("Cannot build graph from empty FASTQ records"))
    end

    # Create MetaGraphsNext graph for quality-aware strings
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=QualityStringVertexData,
        edge_data_type=QualityStringEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )

    ngram_quality_data = Dict{String, Vector{Vector{Int}}}()
    ngram_coverage_data = Dict{String, Vector{Tuple{Int, Int, StrandOrientation, Vector{Int}}}}()

    for (obs_idx, record) in enumerate(fastq_records)
        sequence_str = FASTX.sequence(String, record)
        phred_scores = Int.(collect(FASTX.quality_scores(record)))
        observed_ngrams = ngrams(sequence_str, ngram_length)

        for (pos, ngram) in enumerate(observed_ngrams)
            quality_slice = phred_scores[pos:pos + ngram_length - 1]
            push!(get!(Vector{Vector{Int}}, ngram_quality_data, ngram), quality_slice)
            push!(
                get!(Vector{Tuple{Int, Int, StrandOrientation, Vector{Int}}}, ngram_coverage_data, ngram),
                (obs_idx, pos, Forward, quality_slice)
            )
        end
    end

    for (ngram, quality_vectors) in ngram_quality_data
        coverage = get(ngram_coverage_data, ngram, Tuple{Int, Int, StrandOrientation, Vector{Int}}[])
        graph[ngram] = QualityStringVertexData(ngram, quality_vectors, coverage)
    end

    return graph
# end
=#  # end begin_qualmer_string_graph_next

# Helper function to create empty qualmer graph
#= begin_create_empty_qualmer_graph
function _create_empty_qualmer_graph(kmer_type)
    return MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=kmer_type,
        vertex_data_type=QualmerVertexData{kmer_type},
        edge_data_type=QualmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
end
=#  # end begin_create_empty_qualmer_graph

# Helper function to add qualmer observation to graph
#= begin_add_qualmer_observation
function _add_qualmer_observation_to_graph!(graph, observation, obs_idx, graph_kmers, graph_mode)
    # Similar to _add_observation_to_graph! but with quality score integration
    observed_sequence = FASTX.sequence(observation)
    quality_scores = FASTX.quality(observation)
    k = length(first(graph_kmers))

    if length(observed_sequence) < k
        observation_id = FASTX.identifier(observation)
        @warn "Skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        return
    end

    # Convert sequence to path through k-mer graph with quality and strand information
    observed_path = _sequence_to_canonical_qualmer_path(graph_kmers, observed_sequence, quality_scores, graph_mode)

    if isempty(observed_path)
        return
    end

    # Add quality-aware coverage to vertices and edges
    for i in 1:length(observed_path)
        canonical_kmer, strand, qual = observed_path[i]
        _add_qualmer_vertex_coverage!(graph, canonical_kmer, obs_idx, i, strand, qual)

        if i < length(observed_path)
            next_canonical_kmer, next_strand, next_qual = observed_path[i + 1]
            _add_qualmer_edge!(graph, canonical_kmer, next_canonical_kmer, strand, next_strand, qual, next_qual)
        end
    end
end
=#  # end begin_add_qualmer_observation

# Note: We don't export specific types/functions - use fully qualified names like Mycelia.probabilistic_walk_next
