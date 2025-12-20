# Look back at this for additional inspiration
# https://github.com/fargolo/TextGraphs.jl

import MetaGraphsNext
import FASTX

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract all n-grams from a string.
Returns a vector containing all contiguous substrings of length `n` from the input string.

# Arguments
- `s`: Input string to extract n-grams from
- `n`: Length of each n-gram (must be positive)

# Returns
- `Vector{String}`: All n-grams of length `n` from the string

# Examples
```julia
ngrams("hello", 3)  # Returns ["hel", "ell", "llo"]
ngrams("abcd", 2)   # Returns ["ab", "bc", "cd"]
```
"""
function ngrams(s::AbstractString, n::Int)
    # Convert to character array to handle Unicode properly
    chars = collect(s)
    len = length(chars)
    count = max(len - n + 1, 0)
    result = Vector{String}(undef, count)
    for i in 1:count
        result[i] = String(chars[i:(i + n - 1)])
    end
    return result
end

# function ngrams(s::AbstractString, n::Int)
#     [s[i:i+n-1] for i in 1:length(s)-n+1]
# end

"""
Vertex metadata for string graphs (ported from `sequence-graphs-next.jl`).
"""
struct StringVertexData
    string_value::String
    coverage::Vector{Int}

    function StringVertexData(string_value::String, coverage::Vector{Int}=Int[])
        new(string_value, coverage)
    end
end

"""
Edge metadata for string graphs.
"""
struct StringEdgeData
    overlap_length::Int
    weight::Float64

    function StringEdgeData(overlap_length::Int, weight::Float64=1.0)
        new(overlap_length, weight)
    end
end

"""
Build a string n-gram graph with typed metadata.
"""
function build_string_graph_next(
    fasta_records::Vector{FASTX.FASTA.Record},
    ngram_length::Int;
    graph_mode::GraphMode = SingleStrand,
)
    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end

    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type = String,
        vertex_data_type = StringVertexData,
        edge_data_type = StringEdgeData,
        weight_function = edge_data_weight,
        default_weight = 0.0,
    )

    ngram_counts = Dict{String, Int}()
    for record in fasta_records
        sequence_str = FASTX.sequence(String, record)
        for i in 1:(length(sequence_str) - ngram_length + 1)
            ngram = sequence_str[i:i + ngram_length - 1]
            ngram_counts[ngram] = get(ngram_counts, ngram, 0) + 1
        end
    end

    for ngram in keys(ngram_counts)
        graph[ngram] = StringVertexData(ngram, Int[])
    end

    for record in fasta_records
        sequence_str = FASTX.sequence(String, record)
        ngrams_local = [sequence_str[i:i + ngram_length - 1] for i in 1:(length(sequence_str) - ngram_length + 1)]

        for i in 1:(length(ngrams_local) - 1)
            curr_ngram = ngrams_local[i]
            next_ngram = ngrams_local[i + 1]

            if haskey(graph, curr_ngram) && haskey(graph, next_ngram) && !MetaGraphsNext.has_edge(graph, curr_ngram, next_ngram)
                graph[curr_ngram, next_ngram] = StringEdgeData(ngram_length - 1, 1.0)
            end
        end
    end

    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a directed graph from a collection of string n-grams with edge and vertex weights.
Creates a MetaGraph where vertices represent unique n-grams and edges represent
transitions between consecutive n-grams across all input strings.

# Arguments
- `strings`: Collection of input strings to build graph from (Vector{String} or similar)
- `n`: N-gram size (length of each substring)

# Returns
- `MetaGraphsNext.MetaGraph`: Directed graph with:
  - Vertex labels: String n-grams
  - Vertex data: Occurrence counts across all strings
  - Edge data: Transition counts between consecutive n-grams
  - Weight function: Uses edge counts as weights

# Details
This function processes multiple strings simultaneously, building a unified graph
that captures n-gram patterns and transitions from the entire collection. This is
particularly useful for:
- Multi-sequence assembly where reads come from the same underlying sequence
- Comparative analysis of related sequences
- Building consensus graphs from multiple observations

# Examples
```julia
strings = ["ATCGATCG", "TCGATCGA", "CGATCGAT"]
graph = strings_to_ngram_graph(strings=strings, n=3)

# Single string (backward compatibility)
graph = strings_to_ngram_graph(strings=["ATCGATCG"], n=3)
```
"""
function strings_to_ngram_graph(; strings::AbstractVector{<:AbstractString}, n::Int)
    if isempty(strings)
        # Return empty graph
        return MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Int,
            edge_data_type=Int,
            weight_function=identity,
            default_weight=0
        )
    end

    # Collect all n-grams from all strings
    all_observed_ngrams = String[]
    for s in strings
        append!(all_observed_ngrams, ngrams(s, n))
    end

    # Create base graph
    mg = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=Int,
        edge_data_type=Int,
        weight_function=identity,
        default_weight=0
    )

    # Set vertex properties (n-gram labels and counts)
    ngram_counts = sort(StatsBase.countmap(all_observed_ngrams))
    for (ngram, count) in ngram_counts
        mg[ngram] = count
    end

    # Collect all edges from all strings
    all_edges = Tuple{String, String}[]
    for s in strings
        string_ngrams = ngrams(s, n)
        if length(string_ngrams) > 1
            string_edges = collect(zip(string_ngrams[1:end-1], string_ngrams[2:end]))
            append!(all_edges, string_edges)
        end
    end

    # Set edge properties (transition counts)
    edge_counts = StatsBase.countmap(all_edges)
    for (edge, count) in edge_counts
        src, dst = edge
        mg[src, dst] = count
    end

    return mg
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a directed graph from string n-grams with edge and vertex weights.
Creates a MetaGraph where vertices represent unique n-grams and edges represent
transitions between consecutive n-grams in the original string.

# Arguments
- `s`: Input string to build graph from
- `n`: N-gram size (length of each substring)

# Returns
- `MetaGraphsNext.MetaGraph`: Directed graph with:
  - Vertex labels: String n-grams
  - Vertex data: Integer counts of n-gram occurrences
  - Edge data: Integer counts of n-gram transitions

# Details
Backward compatibility function - now delegates to the multi-string version.
The graph construction process:
1. Extract all n-grams from the string
2. Count occurrences of each unique n-gram (stored as vertex data)
3. Count transitions between consecutive n-grams (stored as edge data)
4. Build a directed graph representing the n-gram sequence structure

# Examples
```julia
graph = string_to_ngram_graph(s="abcabc", n=2)
# Creates graph with vertices "ab", "bc", "ca" and weighted edges
```
"""
function string_to_ngram_graph(; s, n)
    return strings_to_ngram_graph(strings=[s], n=n)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a visual plot of an n-gram graph using GraphPlot.jl.
Displays the graph structure with vertex labels showing the n-grams.

# Arguments
- `g`: MetaGraph containing n-gram data (typically from `string_to_ngram_graph`)

# Returns
- Plot object showing the graph layout with labeled vertices

# Details
Uses a spring layout algorithm to position vertices and displays n-gram strings
as vertex labels for easy interpretation of the graph structure.

# Examples
```julia
graph = string_to_ngram_graph(s="hello", n=2)
plot_ngram_graph(graph)  # Shows visual representation
```
"""
function plot_ngram_graph(g)
    GraphPlot.gplot(
        g,
        nodelabel=collect(MetaGraphsNext.labels(g)),
        # edgelabel=edge_labels,
        layout=GraphPlot.spring_layout,
    )
end

"""
Find all connected components in the graph.
Returns a vector of vectors of vertex indices.

For directed graphs, uses weakly connected components.
For undirected graphs, uses standard connected components.
"""
function find_connected_components(graph::MetaGraphsNext.MetaGraph)
    underlying = graph.graph
    if Graphs.is_directed(underlying)
        return Graphs.weakly_connected_components(underlying)
    else
        return Graphs.connected_components(underlying)
    end
end

"""
Collapse linear paths (vertices with one incoming and one outgoing edge) into a
simpler graph where sequences are concatenated.
"""
function collapse_unbranching_paths(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    collapsed = MetaGraphsNext.copy(graph)
    vertices = _find_collapsible_vertices(collapsed)
    for v in vertices
        if MetaGraphsNext.has_vertex(collapsed, v)
            _collapse_vertex!(collapsed, v)
        end
    end
    return collapsed
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find vertices that can be collapsed in graph simplification.
Identifies vertices with exactly one incoming and one outgoing edge (directed graphs)
or vertices with degree 2 (undirected graphs) that represent linear path segments.

# Arguments
- `graph`: MetaGraph to analyze for collapsible vertices

# Returns
- `Vector{Int}`: Vertex indices that can be safely collapsed

# Details
For directed graphs: finds vertices with in-degree=1 and out-degree=1
For undirected graphs: finds vertices with degree=2
These represent linear segments that can be collapsed to simplify the graph structure.
"""
function _find_collapsible_vertices(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    verts = Int[]
    g = MetaGraphsNext.underlying_graph(graph)
    for v in MetaGraphsNext.vertices(graph)
        if _is_collapsible(g, v)
            push!(verts, v)
        end
    end
    return verts
end

_is_collapsible(g::Graphs.SimpleDiGraph, v) =
    Graphs.indegree(g, v) == 1 && Graphs.outdegree(g, v) == 1

_is_collapsible(g::Graphs.SimpleGraph, v) = Graphs.degree(g, v) == 2

function _collapse_vertex!(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph}, v)
    g = MetaGraphsNext.underlying_graph(graph)
    preds = collect(Graphs.inneighbors(g, v))
    succs = collect(Graphs.outneighbors(g, v))
    if length(preds) == 1 && length(succs) == 1
        pred, succ = preds[1], succs[1]
        _concatenate_sequences!(graph, pred, v)
        MetaGraphsNext.rem_vertex!(graph, v)
        if !haskey(graph, pred, succ)
            MetaGraphsNext.add_edge!(graph, pred, succ)
        end
    end
end

function _collapse_vertex!(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph}, v)
    g = MetaGraphsNext.underlying_graph(graph)
    neighbors = collect(Graphs.neighbors(g, v))
    if length(neighbors) == 2
        n1, n2 = neighbors
        _concatenate_sequences!(graph, n1, v)
        MetaGraphsNext.rem_vertex!(graph, v)
        if !haskey(graph, n1, n2)
            MetaGraphsNext.add_edge!(graph, n1, n2)
        end
    end
end

function _concatenate_sequences!(graph::MetaGraphsNext.MetaGraph, target, source)
    vdata = MetaGraphsNext.get_vertex_data(graph, source)
    tdata = MetaGraphsNext.get_vertex_data(graph, target)
    if haskey(tdata, :sequence) && haskey(vdata, :sequence)
        new_seq = tdata[:sequence] * vdata[:sequence]
        MetaGraphsNext.set_vertex_data!(graph, target, merge(tdata, Dict(:sequence => new_seq)))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble strings by performing graph walks from each connected component.

Traverses the string graph to reconstruct sequences by walking through connected
components and assembling overlapping n-grams back into full strings.

# Arguments
- `graph`: MetaGraphsNext.MetaGraph containing n-gram vertices and transition edges

# Returns
- `Vector{String}`: Assembled sequences, one for each connected component

# Examples
```julia
# Create a string graph from overlapping sequences
graph = build_string_ngram_metagraph("ATCGATCG", 3)
sequences = assemble_strings(graph)
```

This function is useful for reconstructing original sequences from n-gram decompositions
and for string-based genome assembly approaches.
"""

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble string sequences from an n-gram graph returning ALL contigs.
Returns all assembled sequences as a vector for comprehensive assembly results.

# Arguments
- `graph`: MetaGraph containing n-gram data (typically from `strings_to_ngram_graph`)

# Returns
- `Vector{String}`: All assembled sequences from the graph

# Details
This function returns ALL assembly results, not just the first one. This is the
proper behavior for bioinformatics assembly where multiple contigs are expected.
Use this instead of selecting just the first result for comprehensive assembly.

# Examples
```julia
strings = ["ATCGATCG", "TCGATCGA"]
graph = strings_to_ngram_graph(strings=strings, n=3)
all_contigs = assemble_string_graph(graph)  # Returns all contigs, not just first
```
"""
function assemble_string_graph(graph::MetaGraphsNext.MetaGraph)
    return assemble_strings(graph)  # Return ALL contigs, not just first
end

function assemble_strings(graph::MetaGraphsNext.MetaGraph)
    assemblies = String[]
    @show find_connected_components(graph)
    for comp in find_connected_components(graph)
        @show comp
        seqs = _assemble_component(graph, comp)
        @show seqs
        append!(assemblies, seqs)
    end
    return assemblies
end

function _assemble_component(graph::MetaGraphsNext.MetaGraph{Code, Graph}, comp::Vector{Int}) where {Code,Graph<:Graphs.SimpleDiGraph}
    # g = MetaGraphsNext.underlying_graph(graph)
    g = graph.graph
    @assert Graphs.is_directed(g) "Graph must be directed for this assembly method"
    assemblies = String[]
    sources = [v for v in comp if Graphs.indegree(g, v) == 0]
    isempty(sources) && (sources = [comp[1]])
    for s in sources
        visited = Set{Int}()
        parts = String[]
        function dfs(v)
            v in visited && return
            push!(visited, v)
            # data = MetaGraphsNext.get_vertex_data(graph, v)
            # haskey(data, :sequence) && push!(parts, data[:sequence])
            # push!()
            for n in Graphs.outneighbors(g, v)
                dfs(n)
            end
        end
        dfs(s)
        !isempty(parts) && push!(assemblies, join(parts, ""))
    end
    return assemblies
end

function _assemble_component(graph::MetaGraphsNext.MetaGraph{Code, Graph}, comp::Vector{Int}) where {Code,Graph<:Graphs.SimpleGraph}
    # g = MetaGraphsNext.underlying_graph(graph)
    g = graph.graph
    @assert !Graphs.is_directed(g) "Graph must be undirected for this assembly method"
    # TODO check this
    # @assert Graphs.is_simple(g) "Graph must be simple (no self-loops or multiple edges)"
    assemblies = String[]
    terminals = [v for v in comp if Graphs.degree(g, v) == 1]
    isempty(terminals) && (terminals = [comp[1]])
    for t in terminals
        visited = Set{Int}()
        parts = String[]
        function dfs(v)
            v in visited && return
            push!(visited, v)
            data = MetaGraphsNext.get_vertex_data(graph, v)
            haskey(data, :sequence) && push!(parts, data[:sequence])
            for n in Graphs.neighbors(g, v)
                dfs(n)
            end
        end
        dfs(t)
        !isempty(parts) && push!(assemblies, join(parts, ""))
    end
    return assemblies
end

# """
# Find all connected components in the graph.
# Returns a vector of vectors, where each inner vector contains the vertices
# belonging to one connected component.
# """
# function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
#     # Get the underlying graph from MetaGraphsNext
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
#     # Use Graphs.jl's connected_components function
#     components = Graphs.connected_components(underlying_graph)
    
#     return components
# end

# # Specialized version for directed graphs using weakly connected components
# function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph})
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
#     # For directed graphs, use weakly connected components
#     components = Graphs.weakly_connected_components(underlying_graph)
    
#     return components
# end

# # Specialized version for undirected graphs
# function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph})
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
#     # For undirected graphs, use standard connected components
#     components = Graphs.connected_components(underlying_graph)
    
#     return components
# end

# """
# Collapse linear paths (nodes with exactly one incoming and one outgoing edge)
# into single edges with concatenated sequence data.
# Returns a new collapsed graph.
# """
# function collapse_unbranching_paths(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
#     # Create a copy of the graph to modify
#     collapsed_graph = MetaGraphsNext.copy(graph)
    
#     # Find all vertices that can be collapsed
#     vertices_to_collapse = _find_collapsible_vertices(collapsed_graph)
    
#     # Process each collapsible vertex
#     for vertex in vertices_to_collapse
#         if MetaGraphsNext.has_vertex(collapsed_graph, vertex)  # Check if still exists
#             _collapse_vertex!(collapsed_graph, vertex)
#         end
#     end
    
#     return collapsed_graph
# end

# # Helper function to find collapsible vertices (works for both directed and undirected)
# function _find_collapsible_vertices(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
#     vertices_to_collapse = []
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
#     for vertex in MetaGraphsNext.vertices(graph)
#         if _is_collapsible(underlying_graph, vertex)
#             push!(vertices_to_collapse, vertex)
#         end
#     end
    
#     return vertices_to_collapse
# end

# # For directed graphs: check in-degree = 1 and out-degree = 1
# function _is_collapsible(graph::Graphs.SimpleDiGraph, vertex)
#     return Graphs.indegree(graph, vertex) == 1 && Graphs.outdegree(graph, vertex) == 1
# end

# # For undirected graphs: check degree = 2
# function _is_collapsible(graph::Graphs.SimpleGraph, vertex)
#     return Graphs.degree(graph, vertex) == 2
# end

# # Helper function to collapse a vertex (directed version)
# function _collapse_vertex!(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph}, vertex)
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
#     # Get predecessor and successor
#     predecessors = collect(Graphs.inneighbors(underlying_graph, vertex))
#     successors = collect(Graphs.outneighbors(underlying_graph, vertex))
    
#     if length(predecessors) == 1 && length(successors) == 1
#         pred = predecessors[1]
#         succ = successors[1]
        
#         # Concatenate sequence data
#         _concatenate_sequences!(graph, pred, vertex)
        
#         # Remove the collapsed vertex and add direct edge from pred to succ
#         MetaGraphsNext.rem_vertex!(graph, vertex)
#         if !MetaGraphsNext.has_edge(graph, pred, succ)
#             MetaGraphsNext.add_edge!(graph, pred, succ)
#         end
#     end
# end

# # Helper function to collapse a vertex (undirected version)
# function _collapse_vertex!(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph}, vertex)
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
#     # Get neighbors
#     neighbors = collect(Graphs.neighbors(underlying_graph, vertex))
    
#     if length(neighbors) == 2
#         neighbor1, neighbor2 = neighbors
        
#         # Concatenate sequence data with the first neighbor
#         _concatenate_sequences!(graph, neighbor1, vertex)
        
#         # Remove the collapsed vertex and add direct edge between neighbors
#         MetaGraphsNext.rem_vertex!(graph, vertex)
#         if !MetaGraphsNext.has_edge(graph, neighbor1, neighbor2)
#             MetaGraphsNext.add_edge!(graph, neighbor1, neighbor2)
#         end
#     end
# end

# # Helper function to concatenate sequences
# function _concatenate_sequences!(graph::MetaGraphsNext.MetaGraph, target_vertex, source_vertex)
#     vertex_data = MetaGraphsNext.get_vertex_data(graph, source_vertex)
#     target_data = MetaGraphsNext.get_vertex_data(graph, target_vertex)
    
#     # Concatenate sequences (assuming they have a :sequence field)
#     if haskey(target_data, :sequence) && haskey(vertex_data, :sequence)
#         new_sequence = target_data[:sequence] * vertex_data[:sequence]
#         MetaGraphsNext.set_vertex_data!(graph, target_vertex, 
#             merge(target_data, Dict(:sequence => new_sequence)))
#     end
# end

# """
# Perform graph walks to assemble strings from each connected component.
# Returns a vector of assembled sequences.
# """
# function assemble_strings(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
#     assembled_sequences = String[]
    
#     # Get connected components
#     components = find_connected_components(graph)
    
#     for component in components
#         sequences = _assemble_component(graph, component)
#         append!(assembled_sequences, sequences)
#     end
    
#     return assembled_sequences
# end

# # Directed graph assembly
# function _assemble_component(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph}, component::Vector{Int})
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
#     assembled_sequences = String[]
    
#     # Find source nodes (in-degree = 0) in this component
#     sources = []
#     for vertex in component
#         if Graphs.indegree(underlying_graph, vertex) == 0
#             push!(sources, vertex)
#         end
#     end
    
#     # If no sources, pick the first vertex as starting point
#     if isempty(sources)
#         sources = [component[1]]
#     end
    
#     # Perform depth-first traversal from each source
#     for source in sources
#         visited = Set{Int}()
#         sequence_parts = String[]
        
#         function dfs_walk(vertex)
#             if vertex in visited
#                 return
#             end
#             push!(visited, vertex)
            
#             # Get sequence data from vertex
#             vertex_data = MetaGraphsNext.get_vertex_data(graph, vertex)
#             if haskey(vertex_data, :sequence)
#                 push!(sequence_parts, vertex_data[:sequence])
#             end
            
#             # Continue to unvisited neighbors
#             for neighbor in Graphs.outneighbors(underlying_graph, vertex)
#                 if !(neighbor in visited)
#                     dfs_walk(neighbor)
#                 end
#             end
#         end
        
#         dfs_walk(source)
        
#         if !isempty(sequence_parts)
#             assembled_sequence = join(sequence_parts, "")
#             push!(assembled_sequences, assembled_sequence)
#         end
#     end
    
#     return assembled_sequences
# end

# # Undirected graph assembly
# function _assemble_component(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph}, component::Vector{Int})
#     underlying_graph = MetaGraphsNext.underlying_graph(graph)
#     assembled_sequences = String[]
    
#     # For undirected graphs, find terminal nodes (degree = 1) as starting points
#     terminals = []
#     for vertex in component
#         if Graphs.degree(underlying_graph, vertex) == 1
#             push!(terminals, vertex)
#         end
#     end
    
#     # If no terminals, pick the first vertex as starting point
#     if isempty(terminals)
#         terminals = [component[1]]
#     end
    
#     # Perform depth-first traversal from each terminal
#     for terminal in terminals
#         visited = Set{Int}()
#         sequence_parts = String[]
        
#         function dfs_walk(vertex)
#             if vertex in visited
#                 return
#             end
#             push!(visited, vertex)
            
#             # Get sequence data from vertex
#             vertex_data = MetaGraphsNext.get_vertex_data(graph, vertex)
#             if haskey(vertex_data, :sequence)
#                 push!(sequence_parts, vertex_data[:sequence])
#             end
            
#             # Continue to unvisited neighbors
#             for neighbor in Graphs.neighbors(underlying_graph, vertex)
#                 if !(neighbor in visited)
#                     dfs_walk(neighbor)
#                 end
#             end
#         end
        
#         dfs_walk(terminal)
        
#         if !isempty(sequence_parts)
#             assembled_sequence = join(sequence_parts, "")
#             push!(assembled_sequences, assembled_sequence)
#         end
#     end
    
#     return assembled_sequences
# end




"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a string n-gram graph from FASTA records using existing string graph infrastructure.

This is a wrapper around `strings_to_ngram_graph` that takes FASTA records as input,
compatible with the comprehensive graph construction API.
"""
function build_string_ngram_graph_next(fasta_records::Vector{FASTX.FASTA.Record}, ngram_length::Int;
                                      graph_mode::GraphMode=SingleStrand)

    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end

    # Extract strings from FASTA records
    strings = [FASTX.sequence(String, record) for record in fasta_records]

    # Use existing strings_to_ngram_graph function
    return strings_to_ngram_graph(strings=strings, n=ngram_length)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a string BioSequence-like graph from FASTA records.

This treats strings as sequence-like entities for graph construction,
creating a graph where each entire string sequence becomes a vertex.
"""
function build_string_biosequence_graph_next(fasta_records::Vector{FASTX.FASTA.Record};
                                            graph_mode::GraphMode=SingleStrand)

    if isempty(fasta_records)
        throw(ArgumentError("Cannot build graph from empty FASTA records"))
    end

    # Create MetaGraphsNext graph for string sequences
    graph = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=StringVertexData,
        edge_data_type=StringEdgeData,
        weight_function=edge_data_weight,
        default_weight=0.0
    )

    # Add vertices for each string sequence
    for record in fasta_records
        sequence_str = FASTX.sequence(String, record)

        if !haskey(graph, sequence_str)
            graph[sequence_str] = StringVertexData(sequence_str, Int[])
        end
    end

    return graph
end
