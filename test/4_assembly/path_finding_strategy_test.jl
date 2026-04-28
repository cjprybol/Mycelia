import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

function _strategy_source_sink_labels(graph)
    labels = collect(MetaGraphsNext.labels(graph))
    in_degree = Dict(label => 0 for label in labels)
    out_degree = Dict(label => 0 for label in labels)

    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        out_degree[src] += 1
        in_degree[dst] += 1
    end

    sources = [label for label in labels if in_degree[label] == 0 && out_degree[label] > 0]
    sinks = [label for label in labels if out_degree[label] == 0 && in_degree[label] > 0]
    return first(sources), first(sinks)
end

function _strategy_weighted_bubble_graph()
    records = [
        FASTX.FASTA.Record("path1", BioSequences.dna"ACGTCCTGCA"),
        FASTX.FASTA.Record("path2", BioSequences.dna"ACGTAATGCA")
    ]
    graph = Mycelia.Rhizomorph.build_kmer_graph(
        records, 4; dataset_id = "strategy_bubble", mode = :singlestrand
    )
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
    source, target = _strategy_source_sink_labels(weighted)
    return weighted, source, target
end

function _strategy_vertices(path)
    return [step.vertex_label for step in path.steps]
end

Test.@testset "Path Finding Strategies" begin
    weighted, source, target = _strategy_weighted_bubble_graph()

    Test.@testset "constructors validate bounds" begin
        Test.@test Mycelia.Rhizomorph.GreedyWalkStrategy(1) isa
                   Mycelia.Rhizomorph.PathFindingStrategy
        Test.@test Mycelia.Rhizomorph.ProbabilisticWalkStrategy(1) isa
                   Mycelia.Rhizomorph.PathFindingStrategy
        Test.@test Mycelia.Rhizomorph.ShortestProbabilityPathStrategy() isa
                   Mycelia.Rhizomorph.PathFindingStrategy
        Test.@test Mycelia.Rhizomorph.KShortestPathsStrategy(1) isa
                   Mycelia.Rhizomorph.PathFindingStrategy

        Test.@test_throws ArgumentError Mycelia.Rhizomorph.GreedyWalkStrategy(-1)
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.ProbabilisticWalkStrategy(-1)
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.KShortestPathsStrategy(-1)
    end

    Test.@testset "K-shortest strategy matches direct API" begin
        strategy = Mycelia.Rhizomorph.KShortestPathsStrategy(5)
        via_strategy = Mycelia.Rhizomorph.find_paths(
            weighted, strategy; source = source, target = target
        )
        direct = Mycelia.Rhizomorph.k_shortest_paths(weighted, source, target, 5)

        Test.@test length(via_strategy) == 2
        Test.@test length(via_strategy) == length(direct)
        for (strategy_path, direct_path) in zip(via_strategy, direct)
            Test.@test _strategy_vertices(strategy_path) == _strategy_vertices(direct_path)
            Test.@test strategy_path.total_probability == direct_path.total_probability
        end
    end

    Test.@testset "shortest-probability strategy uses common return shape" begin
        strategy = Mycelia.Rhizomorph.ShortestProbabilityPathStrategy()
        paths = Mycelia.Rhizomorph.find_paths(
            weighted, strategy; source = source, target = target
        )

        Test.@test length(paths) == 1
        Test.@test first(first(paths).steps).vertex_label == source
        Test.@test last(first(paths).steps).vertex_label == target
    end

    Test.@testset "walk strategies use start vertex" begin
        greedy_paths = Mycelia.Rhizomorph.find_paths(
            weighted,
            Mycelia.Rhizomorph.GreedyWalkStrategy(4);
            start_vertex = source
        )
        random_paths = Mycelia.Rhizomorph.find_paths(
            weighted,
            Mycelia.Rhizomorph.ProbabilisticWalkStrategy(4; seed = 42);
            source = source
        )

        Test.@test length(greedy_paths) == 1
        Test.@test length(random_paths) == 1
        Test.@test first(first(greedy_paths).steps).vertex_label == source
        Test.@test first(first(random_paths).steps).vertex_label == source
    end

    Test.@testset "missing strategy arguments throw" begin
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.find_paths(
            weighted,
            Mycelia.Rhizomorph.GreedyWalkStrategy(1)
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.find_paths(
            weighted,
            Mycelia.Rhizomorph.ShortestProbabilityPathStrategy();
            source = source
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.find_paths(
            weighted,
            Mycelia.Rhizomorph.KShortestPathsStrategy(1);
            target = target
        )
    end
end
