import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import Graphs
import MetaGraphsNext

Test.@testset "ngrams" begin
    s = "banana"
    n = 2
    expected = ["ba", "an", "na", "an", "na"]
    Test.@test Mycelia.ngrams(s, n) == expected
end

Test.@testset "string_to_ngram_graph" begin
    s = "banana"
    n = 2
    g = Mycelia.string_to_ngram_graph(; s, n)
    Test.@test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    Test.@test collect(MetaGraphsNext.labels(g)) == ["an", "ba", "na"]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
    comps = Mycelia.find_connected_components(g)
    Test.@test length(comps) == 1 && sort(comps[1]) == [1, 2, 3]
end

Test.@testset "collapse and assemble" begin
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
    Test.@test assembled == ["ABC"]

    collapsed = Mycelia.collapse_unbranching_paths(g)
    Test.@test sort(MetaGraphsNext.labels(collapsed)) == ["A", "C"]
    Test.@test MetaGraphsNext.has_edge(collapsed, "A", "C")
    assembled_collapsed = Mycelia.assemble_strings(collapsed)
    Test.@test assembled_collapsed == ["ABC"]
end

