# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/tda_metrics_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/tda_metrics_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import Graphs
import DataFrames
import MetaGraphsNext

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
        expected = Float64(maximum(summary.metrics.betti1) +
                           maximum(summary.metrics.betti0))
        Test.@test Mycelia.tda_graph_score(summary.metrics) == expected
    end

    Test.@testset "Rhizomorph bubble metric extraction table" begin
        graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
            edge_data_type = Mycelia.Rhizomorph.StringEdgeData
        )

        for label in ["A", "B", "C", "D"]
            graph[label] = Mycelia.Rhizomorph.StringVertexData(label)
        end
        graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
        graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
        graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
        graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

        cfg = Mycelia.TDAConfig(thresholds = [1.0, 2.5, 4.0])
        vertex_weights = Dict("A" => 3.0, "B" => 2.0, "C" => 2.0, "D" => 3.0)

        table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = vertex_weights,
            graph_id = "diamond_bubble",
            filtration = :coverage_min,
            weight_name = :coverage,
            provenance = (assembly_id = "fixture", graph_family = :bubble)
        )

        Test.@test table isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(table) == 3
        Test.@test table.graph_id == fill("diamond_bubble", 3)
        Test.@test table.threshold_index == [1, 2, 3]
        Test.@test table.threshold == [1.0, 2.5, 4.0]
        Test.@test table.betti0 == [1, 2, 0]
        Test.@test table.betti1 == [1, 0, 0]
        Test.@test table.nv == fill(4, 3)
        Test.@test table.ne == fill(4, 3)
        Test.@test table.directed == fill(true, 3)
        Test.@test table.filtration == fill(:coverage_min, 3)
        Test.@test table.weight_name == fill(:coverage, 3)
        Test.@test table.weight_min == fill(2.0, 3)
        Test.@test table.weight_max == fill(3.0, 3)
        Test.@test table.weight_mean == fill(2.5, 3)
        Test.@test table.provenance[1].backend == :graph_betti
        Test.@test table.provenance[1].thresholds == [1.0, 2.5, 4.0]
        Test.@test table.provenance[1].filtration == :coverage_min
        Test.@test table.provenance[1].assembly_id == "fixture"
        Test.@test table.provenance[1].graph_family == :bubble
    end
end
