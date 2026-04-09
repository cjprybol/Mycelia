import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "K-Shortest Paths" begin
    Test.@testset "Single path graph returns K=1" begin
        # Linear sequence ATCGATCG with k=4 has exactly one path
        sequence = BioSequences.dna"ATCGATCG"
        records = [FASTX.FASTA.Record("linear", sequence)]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "linear_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        labels = collect(MetaGraphsNext.labels(weighted))
        source = first(labels)
        target = last(labels)

        paths = Mycelia.Rhizomorph.k_shortest_paths(weighted, source, target, 5)

        # Only 1 path exists in a linear graph
        Test.@test length(paths) >= 1
        Test.@test length(paths) <= 1

        # The single path should be valid
        path = first(paths)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test length(path.steps) > 0
        Test.@test path.total_probability > 0.0
    end

    Test.@testset "Excluded vertices block paths" begin
        seq1 = BioSequences.dna"ATCGTTGA"
        seq2 = BioSequences.dna"ATCAATGA"
        records = [
            FASTX.FASTA.Record("path1", seq1),
            FASTX.FASTA.Record("path2", seq2)
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "bubble_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        labels = collect(MetaGraphsNext.labels(weighted))

        all_excluded = Set(labels[2:end])
        source = labels[1]
        target = labels[end]

        result = Mycelia.Rhizomorph._shortest_path_excluding(
            weighted, source, target, all_excluded, Set{Tuple{
                eltype(labels), eltype(labels)}}()
        )
        Test.@test result === nothing
    end

    Test.@testset "Bubble graph returns 2 paths" begin
        seq1 = BioSequences.dna"ATCGTTGA"
        seq2 = BioSequences.dna"ATCAATGA"
        records = [
            FASTX.FASTA.Record("path1", seq1),
            FASTX.FASTA.Record("path2", seq2)
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "bubble2_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        labels = collect(MetaGraphsNext.labels(weighted))
        in_edges = Dict(l => 0 for l in labels)
        out_edges = Dict(l => 0 for l in labels)
        for (src, dst) in MetaGraphsNext.edge_labels(weighted)
            out_edges[src] = get(out_edges, src, 0) + 1
            in_edges[dst] = get(in_edges, dst, 0) + 1
        end
        sources = [l for l in labels if get(in_edges, l, 0) == 0]
        sinks = [l for l in labels if get(out_edges, l, 0) == 0]

        if !isempty(sources) && !isempty(sinks)
            paths = Mycelia.Rhizomorph.k_shortest_paths(
                weighted, first(sources), first(sinks), 5
            )
            Test.@test length(paths) >= 1
            Test.@test length(paths) <= 2
            if length(paths) >= 2
                Test.@test paths[1].total_probability >= paths[2].total_probability
            end
            for path in paths
                Test.@test path isa Mycelia.Rhizomorph.GraphPath
                Test.@test path.total_probability > 0.0
            end
        end
    end
end
