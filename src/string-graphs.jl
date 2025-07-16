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
        weight_function=x->x,
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
"""
function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    underlying = MetaGraphsNext.underlying_graph(graph)
    return Graphs.connected_components(underlying)
end

""" Specialized method for directed graphs using weakly connected components. """
function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph})
    underlying = MetaGraphsNext.underlying_graph(graph)
    return Graphs.weakly_connected_components(underlying)
end

""" Specialized method for undirected graphs. """
function find_connected_components(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph})
    underlying = MetaGraphsNext.underlying_graph(graph)
    return Graphs.connected_components(underlying)
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
function assemble_strings(graph::MetaGraphsNext.MetaGraph{<:Graphs.AbstractSimpleGraph})
    assemblies = String[]
    for comp in find_connected_components(graph)
        seqs = _assemble_component(graph, comp)
        append!(assemblies, seqs)
    end
    return assemblies
end

function _assemble_component(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleDiGraph}, comp::Vector{Int})
    g = MetaGraphsNext.underlying_graph(graph)
    assemblies = String[]
    sources = [v for v in comp if Graphs.indegree(g, v) == 0]
    isempty(sources) && (sources = [comp[1]])
    for s in sources
        visited = Set{Int}()
        parts = String[]
        function dfs(v)
            v in visited && return
            push!(visited, v)
            data = MetaGraphsNext.get_vertex_data(graph, v)
            haskey(data, :sequence) && push!(parts, data[:sequence])
            for n in Graphs.outneighbors(g, v)
                dfs(n)
            end
        end
        dfs(s)
        !isempty(parts) && push!(assemblies, join(parts, ""))
    end
    return assemblies
end

function _assemble_component(graph::MetaGraphsNext.MetaGraph{<:Graphs.SimpleGraph}, comp::Vector{Int})
    g = MetaGraphsNext.underlying_graph(graph)
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

