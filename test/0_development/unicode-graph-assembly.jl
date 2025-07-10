# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=test -e 'include("test/0_development/unicode-graph-assembly.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=test -e 'import Literate; Literate.notebook("test/0_development/unicode-graph-assembly.jl", "test/0_development", execute=false)'
# ````

import Pkg

Pkg.activate("..")

using Revise
using Test
import Mycelia
import Random
import Graphs
import Graphs
import GraphPlot
import StatsBase
import MetaGraphsNext

# to implement

"""
Find all connected components in the graph.
Returns a vector of vectors, where each inner vector contains the vertices
belonging to one connected component.
"""
function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    # Get the underlying graph from MetaGraphsNext
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
    # Use Graphs.jl's connected_components function
    components = Graphs.connected_components(underlying_graph)
    
    return components
end

# Specialized version for directed graphs using weakly connected components
function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph})
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
    # For directed graphs, use weakly connected components
    components = Graphs.weakly_connected_components(underlying_graph)
    
    return components
end

# Specialized version for undirected graphs
function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph})
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
    # For undirected graphs, use standard connected components
    components = Graphs.connected_components(underlying_graph)
    
    return components
end

"""
Collapse linear paths (nodes with exactly one incoming and one outgoing edge)
into single edges with concatenated sequence data.
Returns a new collapsed graph.
"""
function collapse_unbranching_paths(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    # Create a copy of the graph to modify
    collapsed_graph = MetaGraphsNext.copy(graph)
    
    # Find all vertices that can be collapsed
    vertices_to_collapse = _find_collapsible_vertices(collapsed_graph)
    
    # Process each collapsible vertex
    for vertex in vertices_to_collapse
        if MetaGraphsNext.has_vertex(collapsed_graph, vertex)  # Check if still exists
            _collapse_vertex!(collapsed_graph, vertex)
        end
    end
    
    return collapsed_graph
end

# Helper function to find collapsible vertices (works for both directed and undirected)
function _find_collapsible_vertices(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    vertices_to_collapse = []
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
    for vertex in MetaGraphsNext.vertices(graph)
        if _is_collapsible(underlying_graph, vertex)
            push!(vertices_to_collapse, vertex)
        end
    end
    
    return vertices_to_collapse
end

# For directed graphs: check in-degree = 1 and out-degree = 1
function _is_collapsible(graph::Graphs.SimpleDiGraph, vertex)
    return Graphs.indegree(graph, vertex) == 1 && Graphs.outdegree(graph, vertex) == 1
end

# For undirected graphs: check degree = 2
function _is_collapsible(graph::Graphs.SimpleGraph, vertex)
    return Graphs.degree(graph, vertex) == 2
end

# Helper function to collapse a vertex (directed version)
function _collapse_vertex!(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph}, vertex)
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
    # Get predecessor and successor
    predecessors = collect(Graphs.inneighbors(underlying_graph, vertex))
    successors = collect(Graphs.outneighbors(underlying_graph, vertex))
    
    if length(predecessors) == 1 && length(successors) == 1
        pred = predecessors[1]
        succ = successors[1]
        
        # Concatenate sequence data
        _concatenate_sequences!(graph, pred, vertex)
        
        # Remove the collapsed vertex and add direct edge from pred to succ
        MetaGraphsNext.rem_vertex!(graph, vertex)
        if !MetaGraphsNext.has_edge(graph, pred, succ)
            MetaGraphsNext.add_edge!(graph, pred, succ)
        end
    end
end

# Helper function to collapse a vertex (undirected version)
function _collapse_vertex!(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph}, vertex)
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    
    # Get neighbors
    neighbors = collect(Graphs.neighbors(underlying_graph, vertex))
    
    if length(neighbors) == 2
        neighbor1, neighbor2 = neighbors
        
        # Concatenate sequence data with the first neighbor
        _concatenate_sequences!(graph, neighbor1, vertex)
        
        # Remove the collapsed vertex and add direct edge between neighbors
        MetaGraphsNext.rem_vertex!(graph, vertex)
        if !MetaGraphsNext.has_edge(graph, neighbor1, neighbor2)
            MetaGraphsNext.add_edge!(graph, neighbor1, neighbor2)
        end
    end
end

# Helper function to concatenate sequences
function _concatenate_sequences!(graph::MetaGraphsNext.MetaGraph, target_vertex, source_vertex)
    vertex_data = MetaGraphsNext.get_vertex_data(graph, source_vertex)
    target_data = MetaGraphsNext.get_vertex_data(graph, target_vertex)
    
    # Concatenate sequences (assuming they have a :sequence field)
    if haskey(target_data, :sequence) && haskey(vertex_data, :sequence)
        new_sequence = target_data[:sequence] * vertex_data[:sequence]
        MetaGraphsNext.set_vertex_data!(graph, target_vertex, 
            merge(target_data, Dict(:sequence => new_sequence)))
    end
end

"""
Perform graph walks to assemble strings from each connected component.
Returns a vector of assembled sequences.
"""
function assemble_strings(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    assembled_sequences = String[]
    
    # Get connected components
    components = find_connected_components(graph)
    
    for component in components
        sequences = _assemble_component(graph, component)
        append!(assembled_sequences, sequences)
    end
    
    return assembled_sequences
end

# Directed graph assembly
function _assemble_component(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph}, component::Vector{Int})
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    assembled_sequences = String[]
    
    # Find source nodes (in-degree = 0) in this component
    sources = []
    for vertex in component
        if Graphs.indegree(underlying_graph, vertex) == 0
            push!(sources, vertex)
        end
    end
    
    # If no sources, pick the first vertex as starting point
    if isempty(sources)
        sources = [component[1]]
    end
    
    # Perform depth-first traversal from each source
    for source in sources
        visited = Set{Int}()
        sequence_parts = String[]
        
        function dfs_walk(vertex)
            if vertex in visited
                return
            end
            push!(visited, vertex)
            
            # Get sequence data from vertex
            vertex_data = MetaGraphsNext.get_vertex_data(graph, vertex)
            if haskey(vertex_data, :sequence)
                push!(sequence_parts, vertex_data[:sequence])
            end
            
            # Continue to unvisited neighbors
            for neighbor in Graphs.outneighbors(underlying_graph, vertex)
                if !(neighbor in visited)
                    dfs_walk(neighbor)
                end
            end
        end
        
        dfs_walk(source)
        
        if !isempty(sequence_parts)
            assembled_sequence = join(sequence_parts, "")
            push!(assembled_sequences, assembled_sequence)
        end
    end
    
    return assembled_sequences
end

# Undirected graph assembly
function _assemble_component(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph}, component::Vector{Int})
    underlying_graph = MetaGraphsNext.underlying_graph(graph)
    assembled_sequences = String[]
    
    # For undirected graphs, find terminal nodes (degree = 1) as starting points
    terminals = []
    for vertex in component
        if Graphs.degree(underlying_graph, vertex) == 1
            push!(terminals, vertex)
        end
    end
    
    # If no terminals, pick the first vertex as starting point
    if isempty(terminals)
        terminals = [component[1]]
    end
    
    # Perform depth-first traversal from each terminal
    for terminal in terminals
        visited = Set{Int}()
        sequence_parts = String[]
        
        function dfs_walk(vertex)
            if vertex in visited
                return
            end
            push!(visited, vertex)
            
            # Get sequence data from vertex
            vertex_data = MetaGraphsNext.get_vertex_data(graph, vertex)
            if haskey(vertex_data, :sequence)
                push!(sequence_parts, vertex_data[:sequence])
            end
            
            # Continue to unvisited neighbors
            for neighbor in Graphs.neighbors(underlying_graph, vertex)
                if !(neighbor in visited)
                    dfs_walk(neighbor)
                end
            end
        end
        
        dfs_walk(terminal)
        
        if !isempty(sequence_parts)
            assembled_sequence = join(sequence_parts, "")
            push!(assembled_sequences, assembled_sequence)
        end
    end
    
    return assembled_sequences
end

# create basic directed graph, find connected components, collapse unbranching paths, return minimal graph

# create basic undirectected graph, find connected components, collapse unbranching paths, return minimal graph

# add longer n's, find connected components, collapse unbranching paths, return minimal string graph
@testset "banana string graph" begin
    s = "banana"
    n = 2
    g = Mycelia.string_to_ngram_graph(;s,n)
    display(Mycelia.plot_ngram_graph(g))
    @test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    @test collect(MetaGraphsNext.labels(g)) == ["an", "ba", "na"]
    @test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
end

# add longer n's, find connected components, collapse unbranching paths, return minimal string graph

@testset "Mycelia string graph" begin
    s = "mycelia"
    n = 3
    g = Mycelia.string_to_ngram_graph(;s,n)
    display(Mycelia.plot_ngram_graph(g))
    @test collect(Graphs.weights(g)) == [0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 1; 1 0 0 0 0]
    @test collect(MetaGraphsNext.labels(g)) == ["cel", "eli", "lia", "myc", "yce"]
    @test collect(MetaGraphsNext.edge_labels(g)) == [("cel", "eli"), ("eli", "lia"), ("myc", "yce"), ("yce", "cel")]
end

# repeat with `Random.randstring()` with simulated observations

# repeat with `Mycelia.rand_ascii_greek_string()` with simulated observations

# repeat with `Mycelia.rand_latin1_string()` with simulated observations

# repeat with `Mycelia.rand_bmp_printable_string()` with simulated observations

# repeat with `Mycelia.rand_printable_unicode_string()` with simulated observations

Test.@testset "NGram Assembly" begin
    ## Example reads (possibly with simulated errors)
    alphabet = ['A', 'C', 'G', 'T']
    reads = ["ACGT", "CGTA", "GTAC", "TACG"]
    k = 3

    Test.@testset "Graph Construction" begin
        graph = build_ngram_graph(reads, k)
        Test.@test !isnothing(graph)
        ## Further tests on graph structure, nodes, edges, etc.
    end

    Test.@testset "Connected Components" begin
        components = find_connected_components(graph)
        ## Test expected number of components, contents, etc.
    end

    Test.@testset "Collapse Unbranching Paths" begin
        collapsed_graph = collapse_unbranching_paths(graph)
        ## Test properties of collapsed graph
    end

    Test.@testset "String Assembly" begin
        assembled = assemble_strings(graph)
        ## Test that assembled strings match expected sequences
    end
end
