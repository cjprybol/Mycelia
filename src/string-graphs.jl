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