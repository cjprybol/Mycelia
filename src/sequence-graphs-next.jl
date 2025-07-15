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
- `canonical_kmer`: The canonical k-mer sequence (lexicographically smaller of kmer and reverse complement)
"""
struct KmerVertexData
    coverage::Vector{Tuple{Int, Int, StrandOrientation}}  # (observation_id, position, strand_orientation)
    canonical_kmer::String
    
    # Constructor with default empty coverage
    KmerVertexData(canonical_kmer::String) = new(Vector{Tuple{Int, Int, StrandOrientation}}(), canonical_kmer)
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
$(TYPEDSIGNATURES)

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
        return _create_empty_kmer_graph()
    end
    
    # Create the MetaGraphsNext graph with type-stable metadata
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=KmerVertexData,
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices for each canonical k-mer
    for kmer in canonical_kmers
        canonical_kmer_str = string(kmer)
        graph[canonical_kmer_str] = KmerVertexData(canonical_kmer_str)
    end
    
    # Process observations to build strand-aware edges
    for (obs_idx, observation) in enumerate(observations)
        _add_observation_to_graph!(graph, observation, obs_idx, canonical_kmers, graph_mode)
    end
    
    return graph
end
    
    # Count canonical k-mers from observations
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, observations)
    canonical_kmers = collect(keys(canonical_kmer_counts))
    
    if isempty(canonical_kmers)
        @warn "No k-mers found in observations"
        return _create_empty_kmer_graph()
    end
    
    # Create stranded k-mers (forward + reverse complement)
    stranded_kmers = if graph_type == :stranded
        sort!(vcat(canonical_kmers, [BioSequences.reverse_complement(kmer) for kmer in canonical_kmers]))
    else
        sort!(canonical_kmers)
    end
    
    # Create the MetaGraphsNext graph with type-stable metadata
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=KmerVertexData,
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Add vertices for each k-mer
    for (i, kmer) in enumerate(stranded_kmers)
        kmer_str = string(kmer)
        graph[kmer_str] = KmerVertexData(kmer_str)
    end
    
    # Process observations to build edges
    for (obs_idx, observation) in enumerate(observations)
        _add_observation_to_graph!(graph, observation, obs_idx, stranded_kmers, graph_type)
    end
    
    return graph
end

"""
$(TYPEDSIGNATURES)

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
    first_kmer_str = string(first_canonical_kmer)
    _add_vertex_coverage!(graph, first_kmer_str, obs_idx, 1, first_strand)
    
    # Add strand-aware edges and vertex coverage for the rest of the path
    for i in 2:length(observed_path)
        curr_canonical_kmer, curr_strand = observed_path[i]
        prev_canonical_kmer, prev_strand = observed_path[i-1]
        
        curr_kmer_str = string(curr_canonical_kmer)
        prev_kmer_str = string(prev_canonical_kmer)
        
        # Add vertex coverage
        _add_vertex_coverage!(graph, curr_kmer_str, obs_idx, i, curr_strand)
        
        # Add or update strand-aware edge
        _add_strand_aware_edge!(graph, prev_kmer_str, curr_kmer_str, 
                               prev_strand, curr_strand,
                               (obs_idx, i-1, prev_strand), (obs_idx, i, curr_strand))
    end
end

"""
Helper function to add coverage data to a vertex.
"""
function _add_vertex_coverage!(graph, kmer_str, obs_idx, position, strand_orientation)
    vertex_data = graph[kmer_str]
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
    if MetaGraphsNext.has_edge(graph, src_kmer, dst_kmer)
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
function _create_empty_kmer_graph()
    return MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=KmerVertexData,
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
end

"""
$(TYPEDSIGNATURES)

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
            rc_kmer = BioSequences.reverse_complement(observed_kmer)
            
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
$(TYPEDSIGNATURES)

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
    src_seq = src_strand == Forward ? string(src_kmer) : string(BioSequences.reverse_complement(src_kmer))
    dst_seq = dst_strand == Forward ? string(dst_kmer) : string(BioSequences.reverse_complement(dst_kmer))
    
    # Check if suffix of src matches prefix of dst
    src_suffix = src_seq[2:end]  # Remove first nucleotide
    dst_prefix = dst_seq[1:end-1]  # Remove last nucleotide
    
    return src_suffix == dst_prefix
end

# Compatibility layer for migrating from legacy MetaGraphs implementation

"""
$(TYPEDSIGNATURES)

Convert a legacy MetaGraphs-based k-mer graph to the next-generation MetaGraphsNext format.

This function provides a migration path from the deprecated MetaGraphs.jl implementation
to the new type-stable MetaGraphsNext.jl format.

# Arguments
- `legacy_graph`: MetaGraphs.MetaDiGraph from the old implementation

# Returns
- `MetaGraphsNext.MetaGraph` with equivalent structure and type-stable metadata
"""
function legacy_to_next_graph(legacy_graph)
    # Extract metadata from legacy graph
    stranded_kmers = legacy_graph.gprops[:stranded_kmers]
    k = legacy_graph.gprops[:k]
    
    # Create next-generation graph
    next_graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=KmerVertexData,
        edge_data_type=KmerEdgeData,
        weight_function=edge_data -> edge_data.weight,
        default_weight=0.0
    )
    
    # Migrate vertices to canonical representation
    for v in Graphs.vertices(legacy_graph)
        kmer_str = string(stranded_kmers[v])
        
        # Convert to canonical k-mer representation
        kmer_seq = BioSequences.DNAKmer{length(kmer_str)}(kmer_str)
        rc_kmer_seq = BioSequences.reverse_complement(kmer_seq)
        canonical_kmer = string(kmer_seq <= rc_kmer_seq ? kmer_seq : rc_kmer_seq)
        
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
        src_kmer_seq = BioSequences.DNAKmer{length(src_kmer_str)}(src_kmer_str)
        dst_kmer_seq = BioSequences.DNAKmer{length(dst_kmer_str)}(dst_kmer_str)
        
        src_rc = BioSequences.reverse_complement(src_kmer_seq)
        dst_rc = BioSequences.reverse_complement(dst_kmer_seq)
        
        src_canonical = string(src_kmer_seq <= src_rc ? src_kmer_seq : src_rc)
        dst_canonical = string(dst_kmer_seq <= dst_rc ? dst_kmer_seq : dst_rc)
        
        src_strand = src_kmer_seq <= src_rc ? Forward : Reverse
        dst_strand = dst_kmer_seq <= dst_rc ? Forward : Reverse
        
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
        if !MetaGraphsNext.has_edge(next_graph, src_canonical, dst_canonical)
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
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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

# Export the main functions and types
export build_kmer_graph_next, legacy_to_next_graph, is_legacy_graph, ensure_next_graph
export KmerVertexData, KmerEdgeData, StrandOrientation, GraphMode
export Forward, Reverse, SingleStrand, DoubleStrand
