# 
# ```bash
# julia --project=. --color=yes -e 'include("test/0_development/string-graph-helpers.jl")'
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
## using Revise
import Test
import Mycelia
import Graphs
import MetaGraphsNext
import StatsBase

Test.@testset "banana ngrams of n=1" begin
    s = "banana"
    n = 1
    Test.@test Mycelia.ngrams(s, n) == ["b", "a", "n", "a", "n", "a"]
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    g = Mycelia.string_to_ngram_graph(; s, n)
    Test.@test collect(MetaGraphsNext.labels(g)) == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.label_for(g, i) for i in 1:Graphs.nv(g)] == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.code_for(g, l) for l in MetaGraphsNext.labels(g)] == collect(1:Graphs.nv(g))
    Test.@test [g[l] for l in MetaGraphsNext.labels(g)] == collect(values(sorted_ngram_counts))
    Test.@test MetaGraphsNext.default_weight(g) == 0
    Test.@test MetaGraphsNext.get_weight_function(g) == identity
    Test.@test [collect(MetaGraphsNext.all_neighbor_labels(g, l)) for l in MetaGraphsNext.labels(g)] == [["b", "n"], ["a"], ["a"]]
    Test.@test [collect(MetaGraphsNext.inneighbor_labels(g, l)) for l in MetaGraphsNext.labels(g)] == [["b", "n"], String[], ["a"]]
    Test.@test [collect(MetaGraphsNext.outneighbor_labels(g, l)) for l in MetaGraphsNext.labels(g)] == [["n"], ["a"], ["a"]]
    Test.@test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 2 0 0]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("a", "n"), ("b", "a"), ("n", "a")]
    comps = Graphs.connected_components(g)
    Test.@test length(comps) == 1 && sort.(comps) == [[1, 2, 3]]
    
    sources = findall(isempty, MetaGraphsNext.inneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    Test.@test length(sources) == 1
    destinations = findall(isempty, MetaGraphsNext.outneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    Test.@test length(destinations) == 0
    source = first(sources)
    Test.@test [MetaGraphsNext.label_for(g, vertex) for vertex in Graphs.randomwalk(g, source, Graphs.nv(g)*2)] == Mycelia.ngrams(s, n)
    
    Test.@test Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)) == 
    Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)*2) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)*2) == [2, 1, 3]
    
    Test.@test MetaGraphsNext.weighttype(g) == Int
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    MetaGraphsNext.savemg(temp_graph_file, g)
    g2 = MetaGraphsNext.loadmg(temp_graph_file)
    Test.@test g2 == g
    rm(temp_graph_file)

    # dfs_tree_result = Graphs.dfs_tree(g, first(sources))
    # Graphs.dfs_tree_result
end

Test.@testset "banana ngrams of n=2" begin
    s = "banana"
    n = 2
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    Test.@test Mycelia.ngrams(s, n) == ["ba", "an", "na", "an", "na"]
    g = Mycelia.string_to_ngram_graph(; s, n)
    Test.@test collect(MetaGraphsNext.labels(g)) == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.label_for(g, i) for i in 1:Graphs.nv(g)] == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.code_for(g, l) for l in MetaGraphsNext.labels(g)] == collect(1:Graphs.nv(g))
    Test.@test [g[l] for l in MetaGraphsNext.labels(g)] == collect(values(sorted_ngram_counts))
    Test.@test MetaGraphsNext.default_weight(g) == 0
    Test.@test MetaGraphsNext.get_weight_function(g) == identity
    Test.@test [collect(MetaGraphsNext.all_neighbor_labels(g, l)) for l in MetaGraphsNext.labels(g)] == [["ba", "na"], ["an"], ["an"]]
    Test.@test [collect(MetaGraphsNext.inneighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["ba", "na"], String[], ["an"]]
    Test.@test [collect(MetaGraphsNext.outneighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["na"], ["an"], ["an"]]
    Test.@test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
    comps = Graphs.connected_components(g)
    Test.@test length(comps) == 1 && sort.(comps) == [[1, 2, 3]]
    sources = findall(isempty, MetaGraphsNext.inneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    Test.@test length(sources) == 1
    destinations = findall(isempty, MetaGraphsNext.outneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    Test.@test length(destinations) == 0
    source = first(sources)
    Test.@test [MetaGraphsNext.label_for(g, vertex) for vertex in Graphs.randomwalk(g, source, Graphs.nv(g)*2-1)] == Mycelia.ngrams(s, n)
    Test.@test Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)) == 
    Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)*2) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)*2) == [2, 1, 3]
    Test.@test MetaGraphsNext.weighttype(g) == Int
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    MetaGraphsNext.savemg(temp_graph_file, g)
    g2 = MetaGraphsNext.loadmg(temp_graph_file)
    Test.@test g2 == g
    rm(temp_graph_file)
end

Test.@testset "banana ngrams of n=3" begin
    s = "banana"
    n = 3
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    g = Mycelia.string_to_ngram_graph(; s, n)
    Test.@test MetaGraphsNext.weighttype(g) == Int
    Test.@test collect(MetaGraphsNext.labels(g)) == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.label_for(g, i) for i in 1:Graphs.nv(g)] == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.code_for(g, l) for l in MetaGraphsNext.labels(g)] == collect(1:Graphs.nv(g))
    Test.@test [g[l] for l in MetaGraphsNext.labels(g)] == collect(values(sorted_ngram_counts))
    Test.@test MetaGraphsNext.default_weight(g) == 0
    Test.@test MetaGraphsNext.get_weight_function(g) == identity
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    MetaGraphsNext.savemg(temp_graph_file, g)
    g2 = MetaGraphsNext.loadmg(temp_graph_file)
    Test.@test g2 == g
    rm(temp_graph_file)
    comps = Graphs.connected_components(g)
    sources = findall(isempty, MetaGraphsNext.inneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    destinations = findall(isempty, MetaGraphsNext.outneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    Test.@test Mycelia.ngrams(s, n) == ["ban", "ana", "nan", "ana"]
    Test.@test [collect(MetaGraphsNext.all_neighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["ban", "nan"], ["ana"], ["ana"]]
    Test.@test [collect(MetaGraphsNext.inneighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["ban", "nan"], String[], ["ana"]]
    Test.@test [collect(MetaGraphsNext.outneighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["nan"], ["ana"], ["ana"]]
    Test.@test collect(Graphs.weights(g)) == [0 0 1; 1 0 0; 1 0 0]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("ana", "nan"), ("ban", "ana"), ("nan", "ana")]
    Test.@test length(comps) == 1 
    Test.@test sort.(comps) == [[1, 2, 3]]
    Test.@test length(sources) == 1
    Test.@test length(destinations) == 0
    source = first(sources)
    Test.@test [MetaGraphsNext.label_for(g, vertex) for vertex in Graphs.randomwalk(g, source, Graphs.nv(g)*2-2)] == Mycelia.ngrams(s, n)
    Test.@test Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)) == 
    Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)*2) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)*2) == [2, 1, 3]
end

Test.@testset "banana ngrams of n=4" begin
    s = "banana"
    n = 4
    sorted_ngram_counts = sort(StatsBase.countmap(Mycelia.ngrams(s, n)))
    g = Mycelia.string_to_ngram_graph(; s, n)
    Test.@test MetaGraphsNext.weighttype(g) == Int
    Test.@test collect(MetaGraphsNext.labels(g)) == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.label_for(g, i) for i in 1:Graphs.nv(g)] == collect(keys(sorted_ngram_counts))
    Test.@test [MetaGraphsNext.code_for(g, l) for l in MetaGraphsNext.labels(g)] == collect(1:Graphs.nv(g))
    Test.@test [g[l] for l in MetaGraphsNext.labels(g)] == collect(values(sorted_ngram_counts))
    Test.@test MetaGraphsNext.default_weight(g) == 0
    Test.@test MetaGraphsNext.get_weight_function(g) == identity
    temp_graph_file = "test.$(s).$(n).string-graph.jld2"
    MetaGraphsNext.savemg(temp_graph_file, g)
    g2 = MetaGraphsNext.loadmg(temp_graph_file)
    Test.@test g2 == g
    rm(temp_graph_file)
    comps = Graphs.connected_components(g)
    sources = findall(isempty, MetaGraphsNext.inneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    destinations = findall(isempty, MetaGraphsNext.outneighbor_labels(g, l) for l in MetaGraphsNext.labels(g))
    Test.@test Mycelia.ngrams(s, n) == ["bana", "anan", "nana"]
    Test.@test [collect(MetaGraphsNext.all_neighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["bana", "nana"], ["anan"], ["anan"]]
    Test.@test [collect(MetaGraphsNext.inneighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["bana"], String[], ["anan"]]
    Test.@test [collect(MetaGraphsNext.outneighbor_labels(g, l)) for l = MetaGraphsNext.labels(g)] == [["nana"], ["anan"], String[]]
    Test.@test collect(Graphs.weights(g)) == [0 0 1; 1 0 0; 0 0 0]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("anan", "nana"), ("bana", "anan")]
    Test.@test length(comps) == 1 
    Test.@test sort.(comps) == [[1, 2, 3]]
    Test.@test length(sources) == 1
    Test.@test length(destinations) == 1
    source = first(sources)
    destination = first(destinations)
    Test.@test [MetaGraphsNext.label_for(g, vertex) for vertex in Graphs.randomwalk(g, source, Graphs.nv(g))] == 
                [MetaGraphsNext.label_for(g, vertex) for vertex in Graphs.randomwalk(g, source, Graphs.nv(g)*100)] == Mycelia.ngrams(s, n)
    Test.@test Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)) == 
    Graphs.non_backtracking_randomwalk(g, source, Graphs.nv(g)*2) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)) ==
    Graphs.self_avoiding_walk(g, source, Graphs.nv(g)*2) == 
    Graphs.dag_longest_path(g) == Graphs.enumerate_paths(Graphs.dijkstra_shortest_paths(g, source), destination) == [2, 1, 3]
    Test.@test Graphs.a_star(g, source, destination) == [Graphs.Edge(2, 1), Graphs.Edge(1, 3)]
    Test.@test Graphs.yen_k_shortest_paths(g, source, destination, Graphs.weights(g)).paths == [[2, 1, 3]]
end
