# Look back at this for additional inspiration
# https://github.com/fargolo/TextGraphs.jl

function ngrams(s::AbstractString, n::Int)
    len = lastindex(s)
    count = max(len - n + 1, 0)
    result = Vector{String}(undef, count)
    for i in 1:count
        result[i] = s[i:(i + n - 1)]
    end
    return result
end

# function ngrams(s::AbstractString, n::Int)
#     [s[i:i+n-1] for i in 1:length(s)-n+1]
# end

function string_to_ngram_graph(; s, n)
    observed_ngrams = ngrams(s, n)
    unique_ngrams = sort(collect(Set(observed_ngrams)))
    # ngram_to_vertex = Dict(ng => i for (i, ng) in enumerate(unique_ngrams))

    # Create base graph and metagraph
    # g = Graphs.SimpleDiGraph(length(unique_ngrams))
    mg = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=Int,
        edge_data_type=Int,
        weight_function=identity,
        default_weight=0
    )

    # Set node properties (label and count)
    ngram_counts = sort(StatsBase.countmap(observed_ngrams))
    for (ngram, count) in ngram_counts
        # display(count)
        mg[ngram] = count
    end

    ngram_edges = zip(observed_ngrams[1:end-1], observed_ngrams[2:end])
    edge_counts = StatsBase.countmap(ngram_edges)
    for (edge, count) in edge_counts
        src, dst = edge
        mg[src, dst] = count
    end

    return mg
end

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
Assemble strings by performing graph walks from each connected component.
"""
# function assemble_strings(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
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


