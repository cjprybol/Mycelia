import Test
import Mycelia
import MetaGraphsNext
import Graphs

Test.@testset "Repeat Detection" begin
    Test.@testset "Repeat region identification" begin
        graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
            edge_data_type = Mycelia.Rhizomorph.StringEdgeData
        )

        repeat_center = "R"
        in_nodes = ["A", "B", "C"]
        out_nodes = ["D", "E", "F"]
        incoming_sources = ["X1", "X2", "X3"]
        outgoing_targets = ["Y1", "Y2", "Y3"]

        for label in
            vcat([repeat_center], in_nodes, out_nodes, incoming_sources, outgoing_targets)
            graph[label] = Mycelia.Rhizomorph.StringVertexData(label)
        end

        for label in in_nodes
            graph[label, repeat_center] = Mycelia.Rhizomorph.StringEdgeData(1)
        end
        for label in out_nodes
            graph[repeat_center, label] = Mycelia.Rhizomorph.StringEdgeData(1)
        end
        for (src, dst) in zip(incoming_sources, in_nodes)
            graph[src, dst] = Mycelia.Rhizomorph.StringEdgeData(1)
        end
        for (src, dst) in zip(out_nodes, outgoing_targets)
            graph[src, dst] = Mycelia.Rhizomorph.StringEdgeData(1)
        end

        for (idx, label) in enumerate(vcat([repeat_center], in_nodes, out_nodes))
            vertex_data = graph[label]
            for pos in 1:12
                Mycelia.Rhizomorph.add_evidence!(
                    vertex_data,
                    "dataset_01",
                    "obs_$(idx)",
                    Mycelia.Rhizomorph.EvidenceEntry(pos, Mycelia.Rhizomorph.Forward)
                )
            end
        end

        in_degrees, out_degrees = Mycelia.Rhizomorph.calculate_degrees(graph)
        Test.@test in_degrees[repeat_center] == 3
        Test.@test out_degrees[repeat_center] == 3

        candidates = Mycelia.Rhizomorph.find_repeat_candidates(in_degrees, out_degrees)
        Test.@test Set(candidates) == Set([repeat_center])

        local_vertices = Mycelia.Rhizomorph.get_local_subgraph(graph, repeat_center, 1)
        Test.@test repeat_center in local_vertices
        Test.@test all(label -> label in local_vertices, vcat(in_nodes, out_nodes))
        Test.@test all(label -> !(label in local_vertices), vcat(incoming_sources, outgoing_targets))

        repeats = Mycelia.Rhizomorph.resolve_repeats_next(graph; min_repeat_length = 1)
        Test.@test length(repeats) == 1

        region = repeats[1]
        Test.@test region.repeat_type == :interspersed
        Test.@test length(region.incoming_edges) == 3
        Test.@test length(region.outgoing_edges) == 3
        Test.@test region.copy_number_estimate >= 1.0
        Test.@test 0.0 <= region.confidence <= 1.0
    end

    Test.@testset "Repeat region merging" begin
        region1 = Mycelia.Rhizomorph.RepeatRegion(
            ["A", "B"],
            Tuple{String, String}[],
            Tuple{String, String}[],
            2.0,
            :tandem,
            0.6
        )
        region2 = Mycelia.Rhizomorph.RepeatRegion(
            ["B", "C"],
            Tuple{String, String}[],
            Tuple{String, String}[],
            3.0,
            :tandem,
            0.8
        )

        Test.@test Mycelia.Rhizomorph.regions_overlap(region1, region2)

        merged = Mycelia.Rhizomorph.merge_overlapping_repeats([region1, region2])
        Test.@test length(merged) == 1
        Test.@test Set(merged[1].repeat_vertices) == Set(["A", "B", "C"])
        Test.@test merged[1].repeat_type == :tandem
    end

    Test.@testset "Repeat region validation" begin
        Test.@test_throws AssertionError Mycelia.Rhizomorph.RepeatRegion(
            ["A"],
            Tuple{String, String}[],
            Tuple{String, String}[],
            1.0,
            :invalid,
            0.5
        )
        Test.@test_throws AssertionError Mycelia.Rhizomorph.RepeatRegion(
            ["A"],
            Tuple{String, String}[],
            Tuple{String, String}[],
            1.0,
            :tandem,
            2.0
        )
    end
end
