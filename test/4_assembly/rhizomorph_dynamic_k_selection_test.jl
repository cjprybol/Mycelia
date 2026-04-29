import Test
import BioSequences
import FASTX
import Mycelia

Test.@testset "Rhizomorph Dynamic K Selection" begin
    Test.@testset "Short observations clamp the search space" begin
        observations = [
            FASTX.FASTA.Record("short1", "ATGATGA"),
            FASTX.FASTA.Record("short2", "ATGATG")
        ]

        plan = Mycelia.Rhizomorph.select_dynamic_kmer_plan(
            observations;
            min_k = 11,
            max_k = 31
        )

        Test.@test plan.initial_k == 5
        Test.@test plan.search_space == [5]
        Test.@test plan.candidate_ks == [5]
        Test.@test plan.min_sequence_length == 6
        Test.@test plan.max_candidate_k == 6
    end

    Test.@testset "Low-complexity observations fall back to the first feasible k" begin
        observations = [
            FASTX.FASTQ.Record("rep1", "AAAAAAAAAAAA", "IIIIIIIIIIII"),
            FASTX.FASTQ.Record("rep2", "AAAAAAAAAAAA", "IIIIIIIIIIII"),
            FASTX.FASTQ.Record("rep3", "AAAAAAAAAAAA", "IIIIIIIIIIII")
        ]

        plan = Mycelia.Rhizomorph.select_dynamic_kmer_plan(
            observations;
            min_k = 3,
            max_k = 7
        )

        Test.@test plan.initial_k == 3
        Test.@test plan.search_space == [3, 5, 7]
        Test.@test plan.candidate_ks == [3, 5]
        Test.@test !plan.singleton_separation_by_k[3]
    end

    Test.@testset "Noisy observations promote a larger starting k" begin
        observations = [
            "ATGATGATGATG",
            "ATGATGATGATG",
            "ATGATGATGATG",
            "ATGATGATGATG",
            "ATGATCATGATG"
        ]

        plan = Mycelia.Rhizomorph.select_dynamic_kmer_plan(
            observations;
            min_k = 3,
            max_k = 11,
            sparsity_threshold = 0.95,
            singleton_threshold = 1
        )

        Test.@test plan.initial_k == 5
        Test.@test plan.search_space == [3, 5, 7, 11]
        Test.@test plan.candidate_ks == [5, 7, 11]
        Test.@test plan.singleton_separation_by_k[5]
        Test.@test plan.sparsity_by_k[3] < 0.95
        Test.@test plan.sparsity_by_k[5] > plan.sparsity_by_k[3]
        Test.@test Set(keys(plan.sparsity_by_k)) == Set(plan.search_space)
        Test.@test plan.sequence_count == 5
        Test.@test plan.median_sequence_length == 12.0
    end

    Test.@testset "Prime progression stays deterministic" begin
        Test.@test Mycelia.Rhizomorph.dynamic_k_prime_pattern(5; max_k = 11) == [5, 7, 11]
        Test.@test Mycelia.Rhizomorph.dynamic_k_prime_pattern(11; max_k = 31) == [11, 13, 17, 23, 31]
    end

    Test.@testset "BioSequence observations are accepted directly" begin
        observations = [
            BioSequences.LongDNA{4}("ATGATGATGATG"),
            BioSequences.LongDNA{4}("ATGATGATGATG"),
            BioSequences.LongDNA{4}("ATGATCATGATG")
        ]

        plan = Mycelia.Rhizomorph.select_dynamic_kmer_plan(
            observations;
            min_k = 3,
            max_k = 11,
            sparsity_threshold = 0.95,
            singleton_threshold = 1
        )

        Test.@test plan.initial_k == 5
        Test.@test plan.candidate_ks == [5, 7, 11]
        Test.@test plan.sequence_count == 3
        Test.@test plan.singleton_separation_by_k[5]
    end

    Test.@testset "Unicode string observations are indexed safely" begin
        observations = [
            "αβγαβγαβγ",
            "αβγαβγαβγ",
            "αβγαβδαβγ"
        ]

        plan = Mycelia.Rhizomorph.select_dynamic_kmer_plan(
            observations;
            min_k = 3,
            max_k = 7,
            sparsity_threshold = 0.9,
            singleton_threshold = 1
        )

        Test.@test plan.initial_k == 3
        Test.@test plan.search_space == [3, 5, 7]
        Test.@test plan.sequence_count == 3
        Test.@test haskey(plan.sparsity_by_k, 3)
    end
end
