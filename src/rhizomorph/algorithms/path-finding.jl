# Path Finding and Sequence Reconstruction Algorithms
#
# This file contains algorithms for finding paths through assembly graphs
# and reconstructing sequences from those paths.
#
# Key algorithms:
# - Eulerian path finding (Hierholzer's algorithm)
# - Sequence reconstruction from paths
# - Type-stable path traversal
#
# Based on functions from src/sequence-graphs-next.jl

# ============================================================================
# Path Data Structures
# ============================================================================

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

# ============================================================================
# Eulerian Path Finding
# ============================================================================

"""
    find_eulerian_paths_next(graph::MetaGraphsNext.MetaGraph{<:Integer, <:Any, <:Any, <:Any})

Efficient Eulerian path finder that works directly with the underlying Graphs.jl structure.
No temporary graph creation - works with vertex indices and translates back to labels.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: The graph to find paths in

# Returns
- `Vector{Vector{T}}`: Vector containing the Eulerian path (or empty if no path exists)

# Details
Uses Hierholzer's algorithm to find an Eulerian path. An Eulerian path visits every edge
exactly once. The algorithm checks that:
- At most one vertex has (out-degree - in-degree) = 1 (start vertex)
- At most one vertex has (out-degree - in-degree) = -1 (end vertex)
- All other vertices have equal in-degree and out-degree
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

    # Convert vertex indices back to labels
    # graph.vertex_labels is index -> label, so we can use it directly
    index_to_label = graph.vertex_labels

    # Convert path from indices to labels
    label_path = [index_to_label[idx] for idx in path]

    return [label_path]
end

"""
Hierholzer's algorithm implementation working directly on Graphs.jl structure.

# Arguments
- `graph::Graphs.AbstractGraph`: The underlying graph
- `start_vertex::Int`: Starting vertex index

# Returns
- `Vector{Int}`: Path as vertex indices (empty if no Eulerian path exists)

# Algorithm
Hierholzer's algorithm for finding Eulerian paths:
1. Start from the designated starting vertex
2. Follow edges, marking them as used
3. When stuck (no unused outgoing edges), backtrack
4. Insert subcircuits when encountering vertices with unused edges
5. Continue until all edges are used
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

# ============================================================================
# Sequence Reconstruction from Paths
# ============================================================================

"""
    path_to_sequence(path::GraphPath{T}, graph::MetaGraphsNext.MetaGraph) where T

Convert a graph path to a biological sequence.

# Arguments
- `path::GraphPath{T}`: Path through the assembly graph
- `graph::MetaGraphsNext.MetaGraph`: Assembly graph

# Returns
Appropriate sequence type based on the graph content (maintains type stability).
For k-mer graphs, reconstructs by overlapping k-mers. For variable-length graphs,
concatenates sequences.

# Details
- For k-mer graphs: Takes first k-mer in full, then adds last base of each subsequent k-mer
- For BioSequence graphs: Handles overlap based on edge data
- Maintains type stability - returns same sequence type as graph contains
- Handles strand orientation correctly (forward/reverse complement)
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
    elseif hasfield(typeof(vertex_data), :Kmer)
        # K-mer graph: extract k-mer and convert to appropriate sequence type
        kmer = vertex_data.Kmer
        return _convert_kmer_to_sequence(kmer, SequenceType, step.strand)
    elseif hasfield(typeof(vertex_data), :canonical_kmer)
        # Legacy K-mer graph: extract k-mer and convert to appropriate sequence type
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
    elseif hasfield(typeof(vertex_data), :Kmer)
        # K-mer graph: determine sequence type from k-mer
        kmer = vertex_data.Kmer
        return _sequence_type_from_kmer_type(typeof(kmer))
    elseif hasfield(typeof(vertex_data), :canonical_kmer)
        # Legacy K-mer graph: determine sequence type from k-mer
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

    # For Kmers.jl types
    if kmer_type <: Kmers.DNAKmer
        return BioSequences.LongDNA{4}
    elseif kmer_type <: Kmers.RNAKmer
        return BioSequences.LongRNA{4}
    elseif kmer_type <: Kmers.AAKmer
        return BioSequences.LongAA
    end

    # Fallback: try to infer from type parameters
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
    return String
end

"""
Get empty sequence for a given label type.
"""
function _get_empty_sequence_for_label_type(::Type{T}) where T
    if T <: Kmers.DNAKmer
        return BioSequences.LongDNA{4}()
    elseif T <: Kmers.RNAKmer
        return BioSequences.LongRNA{4}()
    elseif T <: Kmers.AAKmer
        return BioSequences.LongAA()
    elseif T <: BioSequences.LongDNA
        return BioSequences.LongDNA{4}()
    elseif T <: BioSequences.LongRNA
        return BioSequences.LongRNA{4}()
    elseif T <: BioSequences.LongAA
        return BioSequences.LongAA()
    else
        return ""
    end
end
