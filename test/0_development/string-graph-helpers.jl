# 
# ```bash
# julia --project=. -e 'include("test/0_development/unicode-graph-assembly.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run from the Mycelia base directory:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/0_development/string-graph-helpers.jl", "test/0_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
import Test
import Mycelia
import Graphs
import MetaGraphsNext

# Test string graph helpers
Test.@testset "ngrams" begin
    s = "banana"
    n = 2
    expected = ["ba", "an", "na", "an", "na"]
    Test.@test Mycelia.ngrams(s, n) == expected
end

# Test n-gram graph construction
Test.@testset "string_to_ngram_graph" begin
    s = "banana"
    n = 2
    g = Mycelia.string_to_ngram_graph(; s, n)
    Test.@test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    Test.@test collect(MetaGraphsNext.labels(g)) == ["an", "ba", "na"]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
    comps = Graphs.connected_components(g)
    Test.@test length(comps) == 1 && sort(comps[1]) == [1, 2, 3]
end

# Test n-gram graph assembly
Test.@testset "collapse and assemble" begin
    g = MetaGraphsNext.MetaGraph(
        Graphs.SimpleDiGraph(),
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

