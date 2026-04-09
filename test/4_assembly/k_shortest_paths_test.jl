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
end
