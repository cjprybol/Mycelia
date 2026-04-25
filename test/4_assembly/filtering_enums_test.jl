import Test
import Mycelia
import MetaGraphsNext
import Graphs
import FASTX
import Kmers

struct CoverageVertexData
    coverage::Vector{Int}
end

struct PlainVertexData
    label::String
end

function make_test_graph(vertex_pairs, edge_pairs; vertex_data_type)
    graph = MetaGraphsNext.MetaGraph(
        Graphs.Graph();
        label_type = String,
        vertex_data_type = vertex_data_type,
        edge_data_type = Int
    )

    for (label, data) in vertex_pairs
        graph[label] = data
    end

    for (src, dst) in edge_pairs
        graph[src, dst] = 1
    end

    return graph
end

Test.@testset "Rhizomorph filtering and enums" begin
    Test.@testset "Enums expose expected values" begin
        Test.@test instances(Mycelia.Rhizomorph.StrandOrientation) ==
                   (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
        Test.@test Bool(Mycelia.Rhizomorph.Forward)
        Test.@test !Bool(Mycelia.Rhizomorph.Reverse)
        Test.@test string(Mycelia.Rhizomorph.Forward) == "Forward"
        Test.@test string(Mycelia.Rhizomorph.Reverse) == "Reverse"

        Test.@test instances(Mycelia.Rhizomorph.GraphMode) ==
                   (Mycelia.Rhizomorph.SingleStrand, Mycelia.Rhizomorph.DoubleStrand)
        Test.@test Int(Mycelia.Rhizomorph.SingleStrand) == 0
        Test.@test Int(Mycelia.Rhizomorph.DoubleStrand) == 1
        Test.@test string(Mycelia.Rhizomorph.SingleStrand) == "SingleStrand"
        Test.@test string(Mycelia.Rhizomorph.DoubleStrand) == "DoubleStrand"
    end

    Test.@testset "filter_largest_components" begin
        graph = make_test_graph(
            [
                "A" => PlainVertexData("alpha"),
                "B" => PlainVertexData("beta"),
                "C" => PlainVertexData("gamma"),
                "D" => PlainVertexData("delta"),
                "E" => PlainVertexData("epsilon"),
                "F" => PlainVertexData("zeta")
            ],
            [("A", "B"), ("B", "C"), ("D", "E")];
            vertex_data_type = PlainVertexData
        )

        largest_only = Mycelia.Rhizomorph.filter_largest_components(graph; max_components = 1)
        Test.@test Set(MetaGraphsNext.labels(largest_only)) == Set(["A", "B", "C"])
        Test.@test length(collect(MetaGraphsNext.edge_labels(largest_only))) == 2

        threshold_kept = Mycelia.Rhizomorph.filter_largest_components(graph; min_nodes = 2)
        Test.@test Set(MetaGraphsNext.labels(threshold_kept)) == Set(["A", "B", "C", "D", "E"])
        Test.@test length(collect(MetaGraphsNext.edge_labels(threshold_kept))) == 3

        fallback_largest = Mycelia.Rhizomorph.filter_largest_components(
            graph;
            min_fraction = 0.9,
            min_nodes = 7
        )
        Test.@test Set(MetaGraphsNext.labels(fallback_largest)) == Set(["A", "B", "C"])

        empty_graph = make_test_graph(
            Pair{String, PlainVertexData}[],
            Tuple{String, String}[];
            vertex_data_type = PlainVertexData
        )
        empty_result = Mycelia.Rhizomorph.filter_largest_components(empty_graph)
        Test.@test isempty(collect(MetaGraphsNext.labels(empty_result)))
        Test.@test isempty(collect(MetaGraphsNext.edge_labels(empty_result)))
    end

    Test.@testset "filter_top_coverage" begin
        evidence_graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [
                FASTX.FASTA.Record("read1", "ATCA"),
                FASTX.FASTA.Record("read2", "ATCA"),
                FASTX.FASTA.Record("read3", "ATCG")
            ],
            3;
            dataset_id = "coverage_test"
        )
        top_evidence = Mycelia.Rhizomorph.filter_top_coverage(evidence_graph; top_fraction = 0.34)
        Test.@test Set(MetaGraphsNext.labels(top_evidence)) ==
                   Set([Kmers.DNAKmer{3}("ATC"), Kmers.DNAKmer{3}("TCA")])
        Test.@test length(collect(MetaGraphsNext.edge_labels(top_evidence))) == 1

        coverage_graph = make_test_graph(
            [
                "A" => CoverageVertexData([1, 2, 3]),
                "B" => CoverageVertexData([1, 2]),
                "C" => CoverageVertexData([1]),
                "D" => CoverageVertexData([1, 2, 3, 4])
            ],
            [("A", "D"), ("B", "C")];
            vertex_data_type = CoverageVertexData
        )
        top_coverage = Mycelia.Rhizomorph.filter_top_coverage(coverage_graph; top_fraction = 0.5)
        Test.@test Set(MetaGraphsNext.labels(top_coverage)) == Set(["A", "D"])
        top_coverage_edge = only(collect(MetaGraphsNext.edge_labels(top_coverage)))
        Test.@test Set(top_coverage_edge) == Set(["A", "D"])

        fallback_graph = make_test_graph(
            [
                "x" => PlainVertexData("ex"),
                "y" => PlainVertexData("why"),
                "z" => PlainVertexData("zee")
            ],
            [("x", "y"), ("y", "z")];
            vertex_data_type = PlainVertexData
        )
        fallback_result = Mycelia.Rhizomorph.filter_top_coverage(fallback_graph; top_fraction = 0.1)
        Test.@test Set(MetaGraphsNext.labels(fallback_result)) == Set(["x", "y", "z"])
        Test.@test length(collect(MetaGraphsNext.edge_labels(fallback_result))) == 2

        empty_graph = make_test_graph(
            Pair{String, CoverageVertexData}[],
            Tuple{String, String}[];
            vertex_data_type = CoverageVertexData
        )
        empty_result = Mycelia.Rhizomorph.filter_top_coverage(empty_graph)
        Test.@test isempty(collect(MetaGraphsNext.labels(empty_result)))
        Test.@test isempty(collect(MetaGraphsNext.edge_labels(empty_result)))
    end
end
