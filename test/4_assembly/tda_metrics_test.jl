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

function _tda_cycle_graph()
    graph = Graphs.SimpleGraph(4)
    for (u, v) in ((1, 2), (2, 3), (3, 4), (4, 1))
        Graphs.add_edge!(graph, u, v)
    end
    return graph
end

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

    Test.@testset "TDA table paths on plain graphs" begin
        graph = _tda_cycle_graph()
        cfg = Mycelia.TDAConfig(thresholds = [2.0, 0.0])
        vertex_weights = [1, 2, 3, 4]

        vector_table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = vertex_weights,
            graph_id = "plain_cycle"
        )

        Test.@test vector_table isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(vector_table) == 2
        Test.@test vector_table.threshold == [0.0, 2.0]
        Test.@test vector_table.betti0 == [1, 1]
        Test.@test vector_table.betti1 == [1, 0]
        Test.@test vector_table.weight_name == fill(:vertex_weight, 2)
        Test.@test vector_table.weight_min == fill(1.0, 2)
        Test.@test vector_table.weight_max == fill(4.0, 2)
        Test.@test vector_table.weight_mean == fill(2.5, 2)
        Test.@test vector_table.score == fill(2.0, 2)

        dict_table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4)
        )
        Test.@test dict_table.betti0 == vector_table.betti0
        Test.@test dict_table.betti1 == vector_table.betti1
        Test.@test dict_table.weight_mean == vector_table.weight_mean

        uniform_table = Mycelia.extract_tda_metrics(
            graph,
            Mycelia.TDAConfig(thresholds = Float64[]);
            graph_id = "uniform_cycle"
        )
        Test.@test uniform_table.threshold == [-Inf]
        Test.@test uniform_table.betti0 == [1]
        Test.@test uniform_table.betti1 == [1]
        Test.@test uniform_table.weight_name == [:uniform]
        Test.@test uniform_table.weight_min == [1.0]
        Test.@test uniform_table.weight_max == [1.0]
        Test.@test uniform_table.weight_mean == [1.0]

        summary = Mycelia.tda_on_graph(
            graph,
            cfg;
            vertex_weights = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4)
        )
        rows = Mycelia.tda_metric_rows(
            summary;
            graph_id = "plain_cycle",
            weight_name = :coverage,
            weight_stats = (weight_min = 1.0, weight_max = 4.0, weight_mean = 2.5)
        )
        table = Mycelia.tda_metric_table(
            summary;
            graph_id = "plain_cycle",
            weight_name = :coverage,
            weight_stats = (weight_min = 1.0, weight_max = 4.0, weight_mean = 2.5)
        )

        Test.@test rows isa Vector{NamedTuple}
        Test.@test length(rows) == 2
        Test.@test rows[1].graph_id == "plain_cycle"
        Test.@test rows[1].weight_name == :coverage
        Test.@test table isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(table) == 2
    end

    Test.@testset "plain graph metric extraction weight modes" begin
        graph = Graphs.SimpleGraph(4)
        for (u, v) in ((1, 2), (1, 3), (2, 4), (3, 4))
            Graphs.add_edge!(graph, u, v)
        end

        cfg = Mycelia.TDAConfig(thresholds = [4.0, 1.0, 2.5])
        dict_table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = Dict(1 => 3.0, 2 => 2.0, 3 => 2.0, 4 => 3.0),
            graph_id = "plain_diamond"
        )

        Test.@test dict_table.threshold == [1.0, 2.5, 4.0]
        Test.@test dict_table.betti0 == [1, 2, 0]
        Test.@test dict_table.betti1 == [1, 0, 0]
        Test.@test dict_table.weight_name == fill(:vertex_weight, 3)
        Test.@test dict_table.weight_min == fill(2.0, 3)
        Test.@test dict_table.weight_max == fill(3.0, 3)
        Test.@test dict_table.weight_mean == fill(2.5, 3)

        summary = Mycelia.tda_on_graph(graph, cfg; vertex_weights = [3.0, 2.0, 2.0, 3.0])
        rows = Mycelia.tda_metric_rows(
            summary;
            graph_id = "rows_plain_diamond",
            weight_name = :coverage,
            weight_stats = (weight_min = 2.0, weight_max = 3.0, weight_mean = 2.5)
        )
        table = Mycelia.tda_metric_table(
            summary;
            graph_id = "table_plain_diamond",
            weight_name = :coverage,
            weight_stats = (weight_min = 2.0, weight_max = 3.0, weight_mean = 2.5)
        )

        Test.@test length(rows) == 3
        Test.@test rows[1].graph_id == "rows_plain_diamond"
        Test.@test rows[1].betti1 == 1
        Test.@test table isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(table) == 3
        Test.@test table.graph_id == fill("table_plain_diamond", 3)
    end

    Test.@testset "uniform weights and empty graph statistics" begin
        path_graph = Graphs.SimpleGraph(3)
        Graphs.add_edge!(path_graph, 1, 2)
        Graphs.add_edge!(path_graph, 2, 3)

        metrics = Mycelia.tda_betti_curves(path_graph; thresholds = [0.0])
        Test.@test metrics.thresholds == [0.0]
        Test.@test metrics.betti0 == [1]
        Test.@test metrics.betti1 == [0]

        default_table = Mycelia.extract_tda_metrics(
            path_graph,
            Mycelia.TDAConfig(thresholds = Float64[])
        )
        Test.@test default_table.threshold == [-Inf]
        Test.@test default_table.weight_name == [:uniform]
        Test.@test default_table.weight_min == [1.0]
        Test.@test default_table.weight_max == [1.0]
        Test.@test default_table.weight_mean == [1.0]

        empty_table = Mycelia.extract_tda_metrics(
            Graphs.SimpleGraph(0),
            Mycelia.TDAConfig(thresholds = [0.0])
        )
        Test.@test empty_table.betti0 == [0]
        Test.@test empty_table.betti1 == [0]
        Test.@test isnan(empty_table.weight_min[1])
        Test.@test isnan(empty_table.weight_max[1])
        Test.@test isnan(empty_table.weight_mean[1])

        empty_stats = Mycelia._tda_weight_stats(Float64[])
        Test.@test isnan(empty_stats.weight_min)
        Test.@test isnan(empty_stats.weight_max)
        Test.@test isnan(empty_stats.weight_mean)

        weight_stats = Mycelia._tda_weight_stats([1, 2, 4])
        Test.@test weight_stats.weight_min == 1.0
        Test.@test weight_stats.weight_max == 4.0
        Test.@test weight_stats.weight_mean ≈ 7 / 3

        repeated_stats = Mycelia._tda_weight_stats([2, 2, 5])
        Test.@test repeated_stats.weight_min == 2.0
        Test.@test repeated_stats.weight_max == 5.0
        Test.@test repeated_stats.weight_mean == 3.0
    end

    Test.@testset "TDA helper validation and edge filtrations" begin
        graph = _tda_cycle_graph()
        cfg = Mycelia.TDAConfig(thresholds = [0.0, 1.0])

        Test.@test Mycelia._tda_vertex_weights(graph, nothing) == ones(Float64, 4)
        Test.@test Mycelia._tda_vertex_weights(graph, [1, 2, 3, 4]) == [1.0, 2.0, 3.0, 4.0]
        Test.@test Mycelia._tda_vertex_weights(
            graph,
            Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4)
        ) == [1.0, 2.0, 3.0, 4.0]
        Test.@test_throws ArgumentError Mycelia._tda_vertex_weights(graph, [1.0, 2.0])
        Test.@test_throws ArgumentError Mycelia._tda_vertex_weights(
            graph,
            Dict(1 => 1.0, 2 => 2.0, 4 => 4.0)
        )

        single_vertex_table = Mycelia.extract_tda_metrics(
            Graphs.SimpleGraph(1),
            Mycelia.TDAConfig(thresholds = [0.0])
        )
        Test.@test single_vertex_table.betti0 == [1]
        Test.@test single_vertex_table.betti1 == [0]

        zero_weight_table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = zeros(Float64, 4)
        )
        Test.@test zero_weight_table.threshold == [0.0, 1.0]
        Test.@test zero_weight_table.betti0 == [1, 0]
        Test.@test zero_weight_table.betti1 == [1, 0]
        Test.@test zero_weight_table.weight_min == fill(0.0, 2)
        Test.@test zero_weight_table.weight_max == fill(0.0, 2)
        Test.@test zero_weight_table.weight_mean == fill(0.0, 2)
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

        vector_table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = [3, 2, 2, 3],
            graph_id = "diamond_bubble"
        )
        Test.@test vector_table.betti0 == table.betti0
        Test.@test vector_table.betti1 == table.betti1
        Test.@test vector_table.weight_min == fill(2.0, 3)
        Test.@test vector_table.weight_max == fill(3.0, 3)
        Test.@test vector_table.weight_mean == fill(2.5, 3)

        conflict_table = Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = vertex_weights,
            graph_id = "diamond_bubble",
            filtration = :coverage_min,
            weight_name = :coverage,
            provenance = (
                backend = :caller_supplied,
                thresholds = [99.0],
                filtration = :caller_supplied,
                assembly_id = "fixture"
            )
        )

        Test.@test conflict_table.provenance[1].backend == :graph_betti
        Test.@test conflict_table.provenance[1].thresholds == [1.0, 2.5, 4.0]
        Test.@test conflict_table.provenance[1].filtration == :coverage_min
        Test.@test conflict_table.provenance[1].assembly_id == "fixture"

        Test.@test_throws ArgumentError Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = [1.0, 2.0]
        )
        Test.@test_throws ArgumentError Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = Dict("A" => 3.0, "B" => 2.0, "D" => 3.0)
        )
        Test.@test_throws TypeError Mycelia.extract_tda_metrics(
            graph,
            cfg;
            vertex_weights = "invalid"
        )
        Test.@test_throws ArgumentError Mycelia._tda_graph("invalid")
        Test.@test_throws ArgumentError Mycelia._tda_vertex_weights(graph, "invalid")

        Test.@test_throws MethodError Mycelia.extract_tda_metrics("invalid", cfg)

        plain_graph = Graphs.SimpleGraph(4)
        for (u, v) in ((1, 2), (1, 3), (2, 4), (3, 4))
            Graphs.add_edge!(plain_graph, u, v)
        end

        Test.@test_throws ArgumentError Mycelia.extract_tda_metrics(
            plain_graph,
            cfg;
            vertex_weights = [1.0, 2.0]
        )
        Test.@test_throws ArgumentError Mycelia.extract_tda_metrics(
            plain_graph,
            cfg;
            vertex_weights = Dict(1 => 3.0, 2 => 2.0, 4 => 3.0)
        )
    end
end
