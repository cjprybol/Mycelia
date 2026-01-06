# 
# ```bash
# julia --project=. --color=yes -e 'include("test/4_assembly/string-graph-helpers.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run from the Mycelia base directory:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/string-graph-helpers.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import Graphs
import MetaGraphsNext
import StatsBase

function assert_graph_roundtrip(graph, temp_graph_file)
    MetaGraphsNext.savemg(temp_graph_file, graph)
    loaded = MetaGraphsNext.loadmg(temp_graph_file)
    Test.@test Set(MetaGraphsNext.labels(loaded)) == Set(MetaGraphsNext.labels(graph))
    for label in MetaGraphsNext.labels(graph)
        Test.@test Mycelia.Rhizomorph.count_evidence(loaded[label]) ==
            Mycelia.Rhizomorph.count_evidence(graph[label])
    end
    Test.@test Set(MetaGraphsNext.edge_labels(loaded)) == Set(MetaGraphsNext.edge_labels(graph))
    for edge_label in MetaGraphsNext.edge_labels(graph)
        Test.@test Mycelia.Rhizomorph.count_evidence(loaded[edge_label...]) ==
            Mycelia.Rhizomorph.count_evidence(graph[edge_label...])
    end
    rm(temp_graph_file)
    return nothing
end

Test.@testset "banana ngrams of n=1" begin
    s = "banana"
    n = 1
    Test.@test Mycelia.ngrams(s, n) == ["b", "a", "n", "a", "n", "a"]
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    g = Mycelia.Rhizomorph.build_ngram_graph([s], n; dataset_id="test")
    Test.@test Set(MetaGraphsNext.labels(g)) == Set(keys(sorted_ngram_counts))
    for (ngram, count) in sorted_ngram_counts
        Test.@test Mycelia.Rhizomorph.count_evidence(g[ngram]) == count
    end
    Test.@test Set(MetaGraphsNext.edge_labels(g)) == Set([("a", "n"), ("b", "a"), ("n", "a")])
    comps = Graphs.connected_components(g)
    Test.@test length(comps) == 1 && sort.(comps) == [[1, 2, 3]]

    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    assert_graph_roundtrip(g, temp_graph_file)

    # dfs_tree_result = Graphs.dfs_tree(g, first(sources))
    # Graphs.dfs_tree_result
end

Test.@testset "banana ngrams of n=2" begin
    s = "banana"
    n = 2
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    Test.@test Mycelia.ngrams(s, n) == ["ba", "an", "na", "an", "na"]
    g = Mycelia.Rhizomorph.build_ngram_graph([s], n; dataset_id="test")
    Test.@test Set(MetaGraphsNext.labels(g)) == Set(keys(sorted_ngram_counts))
    for (ngram, count) in sorted_ngram_counts
        Test.@test Mycelia.Rhizomorph.count_evidence(g[ngram]) == count
    end
    Test.@test Set(MetaGraphsNext.edge_labels(g)) == Set([("an", "na"), ("ba", "an"), ("na", "an")])
    comps = Graphs.connected_components(g)
    Test.@test length(comps) == 1 && sort.(comps) == [[1, 2, 3]]
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    assert_graph_roundtrip(g, temp_graph_file)
end

Test.@testset "banana ngrams of n=3" begin
    s = "banana"
    n = 3
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    g = Mycelia.Rhizomorph.build_ngram_graph([s], n; dataset_id="test")
    Test.@test Set(MetaGraphsNext.labels(g)) == Set(keys(sorted_ngram_counts))
    for (ngram, count) in sorted_ngram_counts
        Test.@test Mycelia.Rhizomorph.count_evidence(g[ngram]) == count
    end
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    assert_graph_roundtrip(g, temp_graph_file)
    comps = Graphs.connected_components(g)
    Test.@test Mycelia.ngrams(s, n) == ["ban", "ana", "nan", "ana"]
    Test.@test Set(MetaGraphsNext.edge_labels(g)) == Set([("ana", "nan"), ("ban", "ana"), ("nan", "ana")])
    Test.@test length(comps) == 1 
    Test.@test sort.(comps) == [[1, 2, 3]]
end

Test.@testset "banana ngrams of n=4" begin
    s = "banana"
    n = 4
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    g = Mycelia.Rhizomorph.build_ngram_graph([s], n; dataset_id="test")
    Test.@test Set(MetaGraphsNext.labels(g)) == Set(keys(sorted_ngram_counts))
    for (ngram, count) in sorted_ngram_counts
        Test.@test Mycelia.Rhizomorph.count_evidence(g[ngram]) == count
    end
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    assert_graph_roundtrip(g, temp_graph_file)
    comps = Graphs.connected_components(g)
    Test.@test Mycelia.ngrams(s, n) == ["bana", "anan", "nana"]
    Test.@test Set(MetaGraphsNext.edge_labels(g)) == Set([("anan", "nana"), ("bana", "anan")])
    Test.@test length(comps) == 1 
    Test.@test sort.(comps) == [[1, 2, 3]]
end
