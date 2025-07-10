import Pkg
Pkg.activate("..")

using Test
import Mycelia
import Graphs
import MetaGraphsNext

@testset "ngrams" begin
    s = "banana"
    n = 2
    expected = ["ba", "an", "na", "an", "na"]
    @test Mycelia.ngrams(s, n) == expected
end

@testset "string_to_ngram_graph" begin
    s = "banana"
    n = 2
    g = Mycelia.string_to_ngram_graph(; s, n)
    @test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    @test collect(MetaGraphsNext.labels(g)) == ["an", "ba", "na"]
    @test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
    comps = Mycelia.find_connected_components(g)
    @test length(comps) == 1 && sort(comps[1]) == [1, 2, 3]
end

@testset "collapse and assemble" begin
    g = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=Dict{Symbol,Any},
        edge_data_type=Int,
        weight_function=x->x,
        default_weight=0,
    )
    g["A"] = Dict(:sequence => "A")
    g["B"] = Dict(:sequence => "B")
    g["C"] = Dict(:sequence => "C")
    g["A", "B"] = 1
    g["B", "C"] = 1

    assembled = Mycelia.assemble_strings(g)
    @test assembled == ["ABC"]

    collapsed = Mycelia.collapse_unbranching_paths(g)
    @test sort(MetaGraphsNext.labels(collapsed)) == ["A", "C"]
    @test MetaGraphsNext.has_edge(collapsed, "A", "C")
    assembled_collapsed = Mycelia.assemble_strings(collapsed)
    @test assembled_collapsed == ["ABC"]
end

