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

@testset "Mycelia string graph" begin
    s = "mycelia"
    n = 3
    g = Mycelia.string_to_ngram_graph(;s,n)
    display(Mycelia.plot_ngram_graph(g))
    @test collect(Graphs.weights(g)) == [0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 1; 1 0 0 0 0]
    @test collect(MetaGraphsNext.labels(g)) == ["cel", "eli", "lia", "myc", "yce"]
    @test collect(MetaGraphsNext.edge_labels(g)) == [("cel", "eli"), ("eli", "lia"), ("myc", "yce"), ("yce", "cel")]
end

@testset "banana string graph" begin
    s = "banana"
    n = 2
    g = Mycelia.string_to_ngram_graph(;s,n)
    display(Mycelia.plot_ngram_graph(g))
    @test collect(Graphs.weights(g)) == [0 0 2; 1 0 0; 1 0 0]
    @test collect(MetaGraphsNext.labels(g)) == ["an", "ba", "na"]
    @test collect(MetaGraphsNext.edge_labels(g)) == [("an", "na"), ("ba", "an"), ("na", "an")]
end

# Random.randstring()
# Mycelia.rand_ascii_greek_string(100)
