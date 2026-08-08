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
        Test.@test Mycelia.Rhizomorph.dynamic_k_prime_pattern(11; max_k = 31) ==
                   [11, 13, 17, 23, 31]
    end

    Test.@testset "dynamic_k_prime_pattern is defined exactly once, with keywords" begin
        # A second `dynamic_k_prime_pattern` used to live in
        # `src/development/intelligent-assembly.jl` with a POSITIONAL signature,
        # while the canonical one here takes KEYWORDS.
        #
        # This check is SOURCE-LEVEL on purpose, and the reason is the whole
        # point of the test. `src/development/` is not included by `Mycelia.jl`,
        # so the duplicate never entered any method table — a runtime oracle
        # against `Mycelia.Rhizomorph`'s methods is asking a different generic
        # function about a change that never touched it, and stays green with the
        # duplicate restored. Verified by mutation: re-adding the deleted method
        # to top-level `Mycelia` leaves a runtime check passing.
        #
        # Only a source-level oracle can observe a definition in an unloaded file.
        # NOTE the `m` flag. Without it Julia anchors `^` to the start of the
        # whole STRING, not each line, so reading a file wholesale and matching
        # `^\s*function` finds nothing — the oracle then reports 0 definitions
        # and goes red for a reason unrelated to what it is testing. That is how
        # the first version of this test "passed" its own positive control.
        definition_pattern = r"^\s*function\s+dynamic_k_prime_pattern\("m
        definition_files = String[]
        source_root = joinpath(pkgdir(Mycelia), "src")
        for (root, _, files) in walkdir(source_root), file in files

            endswith(file, ".jl") || continue
            path = joinpath(root, file)
            occursin(definition_pattern, read(path, String)) &&
                push!(definition_files, path)
        end

        canonical_suffix = joinpath("rhizomorph", "algorithms", "dynamic-k-selection.jl")
        Test.@test length(definition_files) == 1
        # `all` rather than `only`, so a second definition FAILS cleanly instead of
        # throwing ArgumentError from `only` — a red test that errors is harder to
        # read than one that reports the assertion it broke.
        Test.@test all(path -> endswith(path, canonical_suffix), definition_files)

        # Pin the METHOD COUNT too, so reintroducing a positional method inside
        # Rhizomorph itself (which a source check alone would also catch, but
        # only if it is spelled as a top-level `function`) fails here as well.
        Test.@test length(methods(Mycelia.Rhizomorph.dynamic_k_prime_pattern)) == 2

        # Assert the VALUE, not the type. The function carries a
        # `::Vector{Int}` return annotation, so Julia converts the result and an
        # `isa Vector{Int}` check can only fail if the call throws — it stayed
        # green under mutants that corrupted every default and the element type.
        Test.@test Mycelia.Rhizomorph.dynamic_k_prime_pattern() ==
                   [11, 13, 17, 23, 31, 41, 53, 67, 83, 101]
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
