import Test
import Mycelia
import Graphs

Test.@testset "TDA metrics (graph invariants)" begin
    Test.@testset "tda_betti_numbers" begin
        path_graph = Graphs.SimpleGraph(5)
        for (u, v) in ((1, 2), (2, 3), (3, 4), (4, 5))
            Graphs.add_edge!(path_graph, u, v)
        end

        betti0, betti1 = Mycelia.tda_betti_numbers(path_graph)
        Test.@test betti0 == 1
        Test.@test betti1 == 0

        cycle_graph = Graphs.SimpleGraph(5)
        for (u, v) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 1))
            Graphs.add_edge!(cycle_graph, u, v)
        end

        betti0, betti1 = Mycelia.tda_betti_numbers(cycle_graph)
        Test.@test betti0 == 1
        Test.@test betti1 == 1

        disconnected_cycles = Graphs.SimpleGraph(6)
        for (u, v) in ((1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 4))
            Graphs.add_edge!(disconnected_cycles, u, v)
        end

        betti0, betti1 = Mycelia.tda_betti_numbers(disconnected_cycles)
        Test.@test betti0 == 2
        Test.@test betti1 == 2

        empty_graph = Graphs.SimpleGraph(0)
        betti0, betti1 = Mycelia.tda_betti_numbers(empty_graph)
        Test.@test betti0 == 0
        Test.@test betti1 == 0
    end

    Test.@testset "tda_on_graph filtration sanity" begin
        cycle_graph = Graphs.SimpleGraph(5)
        for (u, v) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 1))
            Graphs.add_edge!(cycle_graph, u, v)
        end

        weights = Float64[1, 2, 3, 4, 5]
        cfg = Mycelia.TDAConfig(thresholds = [0, 3, 6])
        summary = Mycelia.tda_on_graph(cycle_graph, cfg; vertex_weights = weights)

        Test.@test summary.metrics.thresholds == Float64[0, 3, 6]
        Test.@test summary.metrics.betti0 == [1, 1, 0]
        Test.@test summary.metrics.betti1 == [1, 0, 0]
        Test.@test summary.metrics.persistence === nothing
    end

    Test.@testset "tda_graph_score" begin
        cfg = Mycelia.TDAConfig(thresholds = [0])
        g = Graphs.SimpleGraph(4)
        for (u, v) in ((1, 2), (2, 3), (3, 4), (4, 1))
            Graphs.add_edge!(g, u, v)
        end

        summary = Mycelia.tda_on_graph(g, cfg; vertex_weights = ones(Float64, Graphs.nv(g)))
        expected = Float64(maximum(summary.metrics.betti1) + maximum(summary.metrics.betti0))
        Test.@test Mycelia.tda_graph_score(summary.metrics) == expected
    end
end

