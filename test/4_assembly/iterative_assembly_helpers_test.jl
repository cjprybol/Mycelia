import Test
import Mycelia
import MetaGraphsNext
import Graphs
import Primes

Test.@testset "Iterative Assembly Helpers" begin
    Test.@testset "Sparsity and singleton detection" begin
        reads = [
            Mycelia.fastq_record(identifier = "r1", sequence = "AAAAAA", quality_scores = fill(30, 6)),
            Mycelia.fastq_record(identifier = "r2", sequence = "AT", quality_scores = fill(30, 2))
        ]
        sparsity = Mycelia.calculate_sparsity(reads, 2)
        Test.@test isapprox(sparsity, 0.875; atol = 1e-6)
        Test.@test Mycelia.errors_are_singletons(reads, 2)
    end

    Test.@testset "K selection and memory estimates" begin
        reads = [Mycelia.fastq_record(identifier = "r1", sequence = "AAAAAA", quality_scores = fill(30, 6))]
        Test.@test Mycelia.find_initial_k(reads; k_range = [2, 3], sparsity_threshold = 0.0) ==
                   2
        Test.@test Mycelia.next_prime_k(5; max_k = 10) == 7
        Test.@test Mycelia.next_prime_k(11; max_k = 12) == 11

        estimate = Mycelia.estimate_memory_usage(10, 3)
        Test.@test estimate == 2436
    end

    Test.@testset "build_k_ladder snaps rungs to primes (not odd)" begin
        # The :scalable corrector builds its k-ladder via n_k_rungs. Rungs must be
        # PRIME, not merely odd, to avoid period-p aliasing (9=3x3, 21=3x7 collapse
        # period-3/7 tandem repeats to self-overlapping k-mers). find_initial_k
        # already draws from Primes.primes, so a prime initial_k + this fix makes
        # the whole scalable walk prime end-to-end.

        # The canonical scalable config (initial_k=13, max_k=21, 3 rungs). The old
        # odd-snapping produced [13, 17, 21] with 21=3x7 COMPOSITE; the fix pins the
        # top rung to the largest prime <= max_k (19) and snaps the middle to a prime.
        ladder = Mycelia.build_k_ladder(13, 21; n_k_rungs = 3)
        Test.@test ladder == [13, 17, 19]
        Test.@test !(21 in ladder)                 # the composite top rung is gone
        Test.@test all(Primes.isprime, ladder)

        # Every rung prime across representative scalable configs with a prime initial_k.
        for (ik, mk, nr) in ((11, 31, 3), (3, 51, 5), (13, 41, 4), (17, 23, 3))
            l = Mycelia.build_k_ladder(ik, mk; n_k_rungs = nr)
            Test.@test all(Primes.isprime, l)       # all rungs prime
            Test.@test first(l) == ik               # first rung pinned to initial_k
            Test.@test last(l) == Primes.prevprime(mk)  # last rung = largest prime <= max_k
            Test.@test issorted(l)
            Test.@test allunique(l)
        end

        # Degenerate: initial_k >= max_k returns [initial_k] unchanged.
        Test.@test Mycelia.build_k_ladder(19, 19; n_k_rungs = 3) == [19]

        # Explicit k_ladder branch is user-provided and NOT prime-snapped (unchanged).
        Test.@test Mycelia.build_k_ladder(3, 30; k_ladder = [9, 15, 21]) == [9, 15, 21]

        # No-ladder mode still returns nothing (legacy prime-by-prime walk).
        Test.@test Mycelia.build_k_ladder(13, 21) === nothing
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
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; mode = :singlestrand)
        labels = collect(MetaGraphsNext.labels(graph))
        expected = Mycelia.estimate_memory_usage(length(labels), 3)
        Test.@test Mycelia.check_memory_limits(graph, expected + 1)
        Test.@test !Mycelia.check_memory_limits(graph, expected - 1)
    end
end
