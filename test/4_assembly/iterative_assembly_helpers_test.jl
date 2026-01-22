import Test
import Mycelia
import MetaGraphsNext
import Graphs

Test.@testset "Iterative Assembly Helpers" begin
    Test.@testset "Sparsity and singleton detection" begin
        reads = [
            Mycelia.fastq_record(identifier="r1", sequence="AAAAAA", quality_scores=fill(30, 6)),
            Mycelia.fastq_record(identifier="r2", sequence="AT", quality_scores=fill(30, 2))
        ]
        sparsity = Mycelia.calculate_sparsity(reads, 2)
        Test.@test isapprox(sparsity, 0.875; atol=1e-6)
        Test.@test Mycelia.errors_are_singletons(reads, 2)
    end

    Test.@testset "K selection and memory estimates" begin
        reads = [Mycelia.fastq_record(identifier="r1", sequence="AAAAAA", quality_scores=fill(30, 6))]
        Test.@test Mycelia.find_initial_k(reads; k_range=[2, 3], sparsity_threshold=0.0) == 2
        Test.@test Mycelia.next_prime_k(5; max_k=10) == 7
        Test.@test Mycelia.next_prime_k(11; max_k=12) == 11

        estimate = Mycelia.estimate_memory_usage(10, 3)
        Test.@test estimate == 2436
    end

    Test.@testset "Label helpers" begin
        Test.@test Mycelia._label_k_length("ATG") == 3
        Test.@test Mycelia._label_k_length(Mycelia.Kmers.DNAKmer{3}("ATG")) == 3
        Test.@test Mycelia._label_k_length(Mycelia.BioSequences.LongDNA{4}("ATG")) == 3
    end

    Test.@testset "Memory limit checks" begin
        empty_graph = Graphs.SimpleGraph(0)
        Test.@test Mycelia.check_memory_limits(empty_graph, 1)

        records = [Mycelia.FASTX.FASTA.Record("seq1", "ATGC")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; mode=:singlestrand)
        labels = collect(MetaGraphsNext.labels(graph))
        expected = Mycelia.estimate_memory_usage(length(labels), 3)
        Test.@test Mycelia.check_memory_limits(graph, expected + 1)
        Test.@test !Mycelia.check_memory_limits(graph, expected - 1)
    end
end
