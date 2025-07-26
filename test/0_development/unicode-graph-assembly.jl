# From the Mycelia base directory, run the tests from the root directory with:
# 
# ```bash
# julia --project=. -e 'include("test/0_development/unicode-graph-assembly.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run from the Mycelia base directory:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/0_development/unicode-graph-assembly.jl", "test/0_development", execute=false)'
# ````


## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
import Test
import Mycelia
import Random
import Graphs
import Graphs
import GraphPlot
import StatsBase
import MetaGraphsNext

# to implement


# create basic directed graph, find connected components, collapse unbranching paths, return minimal graph

# create basic undirectected graph, find connected components, collapse unbranching paths, return minimal graph

# add longer n's, find connected components, collapse unbranching paths, return minimal string graph
Test.@testset "banana string graph" begin
    s = "banana"
    n = 2
    g = Mycelia.string_to_ngram_graph(;s,n)
    display(Mycelia.plot_ngram_graph(g))
    Test.@test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    Test.@test collect(MetaGraphsNext.labels(g)) == ["an", "ba", "na"]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
end

# add longer n's, find connected components, collapse unbranching paths, return minimal string graph

Test.@testset "Mycelia string graph" begin
    s = "mycelia"
    n = 3
    g = Mycelia.string_to_ngram_graph(;s,n)
    display(Mycelia.plot_ngram_graph(g))
    Test.@test collect(Graphs.weights(g)) == [0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 1; 1 0 0 0 0]
    Test.@test collect(MetaGraphsNext.labels(g)) == ["cel", "eli", "lia", "myc", "yce"]
    Test.@test collect(MetaGraphsNext.edge_labels(g)) == [("cel", "eli"), ("eli", "lia"), ("myc", "yce"), ("yce", "cel")]
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
