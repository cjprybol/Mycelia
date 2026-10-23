# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/genomic_graph_algorithms_test.jl")'
# ```

# Bidirectional genomic Dijkstra must agree with unidirectional genomic Dijkstra.
#
# Edge-weight definition shared by both searches: the cost of edge (u -> v) is
# `1 + (1 - coverage(kmers[v]) / max_coverage)`, i.e. it depends on the
# DESTINATION k-mer only. Path cost = sum over the path's edges.

import Test
import Mycelia
import Graphs
import Logging
import StableRNGs

Test.@testset "bidirectional_genomic_dijkstra agrees with genomic_dijkstra" begin
    build_case = function (n_vertices::Int, edges::Vector{Tuple{Int, Int}},
            coverages::Vector{Int})
        graph = Graphs.SimpleDiGraph(n_vertices)
        for (u, v) in edges
            Graphs.add_edge!(graph, u, v)
        end
        kmers = ["kmer$(i)" for i in 1:n_vertices]
        kmer_coverage = Dict(kmers[i] => coverages[i] for i in 1:n_vertices)
        return graph, kmers, kmer_coverage
    end

    destination_weighted_cost = function (path_indices::Vector{Int},
            kmers::Vector{String}, kmer_coverage::Dict{String, Int})
        max_coverage = maximum(values(kmer_coverage))
        total = 0.0
        for v in path_indices[2:end]
            total += 1.0 + (1.0 - kmer_coverage[kmers[v]] / max_coverage)
        end
        return total
    end

    check_agreement = function (graph, kmers, kmer_coverage, source::Int, target::Int)
        ## Both functions log every query; keep test output readable.
        expected,
        result = Logging.with_logger(Logging.NullLogger()) do
            (Mycelia.genomic_dijkstra(graph, source, target, kmers, kmer_coverage),
                Mycelia.bidirectional_genomic_dijkstra(
                    graph, source, target, kmers, kmer_coverage))
        end

        if isinf(expected.total_cost)
            Test.@test isinf(result.total_cost)
            Test.@test isempty(result.path)
            Test.@test result.path_length == 0
            return result
        end

        Test.@test result.total_cost ≈ expected.total_cost
        Test.@test !isempty(result.path)
        Test.@test result.path_length == length(result.path)

        kmer_index = Dict(kmer => i for (i, kmer) in enumerate(kmers))
        path_indices = [kmer_index[kmer] for kmer in result.path]
        Test.@test first(path_indices) == source
        Test.@test last(path_indices) == target
        Test.@test all(Graphs.has_edge(graph, path_indices[i], path_indices[i + 1])
        for i in 1:(length(path_indices) - 1))
        ## The returned path must actually cost what is reported.
        Test.@test destination_weighted_cost(path_indices, kmers, kmer_coverage) ≈
                   expected.total_cost
        return result
    end

    Test.@testset "(a) first meeting vertex is not on the shortest path" begin
        ## Direct edge 1 -> 3 (cost 1.99) beats 1 -> 2 -> 3 (cost 1.0 + 1.99).
        ## Vertex 2 is cheap, so both searches settle it first and the
        ## first-meeting rule would return the longer detour through it.
        graph, kmers, coverage = build_case(3, [(1, 2), (1, 3), (2, 3)], [1, 100, 1])
        result = check_agreement(graph, kmers, coverage, 1, 3)
        Test.@test result.path == ["kmer1", "kmer3"]
        Test.@test result.total_cost ≈ 1.99
    end

    Test.@testset "(b) one queue empties while a path still exists" begin
        ## Reverse side exhausts first: target 3 has a single ancestor chain
        ## 1 -> 2 -> 3, while the source fans out into a long cheap branch.
        graph, kmers,
        coverage = build_case(
            9,
            [(1, 2), (2, 3), (1, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)],
            [3, 50, 10, 100, 100, 100, 100, 100, 100]
        )
        check_agreement(graph, kmers, coverage, 1, 3)

        ## Forward side exhausts first: source 1 has a single descendant chain
        ## 1 -> 2 -> 3, while the target has a long cheap ancestor branch.
        graph, kmers,
        coverage = build_case(
            9,
            [(1, 2), (2, 3), (4, 3), (5, 4), (6, 5), (7, 6), (8, 7), (9, 8)],
            [3, 50, 10, 100, 100, 100, 100, 100, 100]
        )
        check_agreement(graph, kmers, coverage, 1, 3)
    end

    Test.@testset "(c) no path exists" begin
        graph, kmers, coverage = build_case(4, [(1, 2), (3, 4), (4, 3)], [5, 5, 5, 5])
        result = check_agreement(graph, kmers, coverage, 1, 4)
        Test.@test isinf(result.total_cost)
    end

    Test.@testset "(d) source == target" begin
        graph, kmers, coverage = build_case(3, [(1, 2), (2, 1), (2, 3)], [5, 50, 1])
        result = check_agreement(graph, kmers, coverage, 2, 2)
        Test.@test result.path == ["kmer2"]
        Test.@test result.total_cost == 0.0
    end

    Test.@testset "seeded random digraphs" begin
        rng = StableRNGs.StableRNG(20261023)
        for _ in 1:200
            n_vertices = rand(rng, 2:9)
            edges = Tuple{Int, Int}[]
            for u in 1:n_vertices, v in 1:n_vertices

                if u != v && rand(rng) < 0.3
                    push!(edges, (u, v))
                end
            end
            coverages = [rand(rng, (1, 3, 10, 50, 100)) for _ in 1:n_vertices]
            graph, kmers, coverage = build_case(n_vertices, edges, coverages)
            source = rand(rng, 1:n_vertices)
            target = rand(rng, 1:n_vertices)
            check_agreement(graph, kmers, coverage, source, target)
        end
    end
end
