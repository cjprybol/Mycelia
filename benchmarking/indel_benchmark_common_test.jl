# Unit test for benchmarking/indel_benchmark_common.jl.
#
# Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/indel_benchmark_common_test.jl
#
# The fixture assertions load Mycelia (they exercise the real read simulator),
# so this is not a millisecond-scale test like calibration_metrics_test.jl. The
# golden digests below are the whole point: they are the acceptance oracle that
# the shared `indel_bench_simulate_reads` is a VERBATIM transplant of the
# per-script fixture helpers it replaced, not merely an equivalent-looking one.

import BioSequences
import FASTX
import Random
import SHA
import Test

include(joinpath(@__DIR__, "indel_benchmark_common.jl"))

# Captured on the pre-consolidation base commit (origin/master 250f7f26) from
# `_indel_toy_make_fixture` in indel_frontier_fixed_toy_proof.jl, using the
# fixed-toy constants below. `Mycelia.observe` draws from the GLOBAL RNG, so
# these digests are sensitive to the exact interleaving of local and global RNG
# draws. If they move, the fixture transplant is wrong — revert it rather than
# re-pinning the constants.
const INDEL_BENCH_TEST_GENOME_LENGTH = 2_000
const INDEL_BENCH_TEST_SOURCE_READ_LENGTH = 1_200
const INDEL_BENCH_TEST_COVERAGE = 8
const INDEL_BENCH_TEST_ERROR_RATE = 0.05
const INDEL_BENCH_TEST_SEED = 42
const INDEL_BENCH_TEST_N_READS = 14
const INDEL_BENCH_TEST_MIN_READ_LENGTH = 1_190
const INDEL_BENCH_TEST_MAX_READ_LENGTH = 1_214
const INDEL_BENCH_TEST_REFERENCE_LENGTH = 2_000
const INDEL_BENCH_TEST_REFERENCE_SHA256 = "e504bccd87f51781ab718bad4fff6e0f222db36271d57c443616500153e75284"
const INDEL_BENCH_TEST_FIXTURE_SHA256 = "1044dc420ca6ad3a883df1bb689c757f45140c1eda9f98c7de26310dbddb2f5d"

"""
    _indel_bench_test_fixture_sha256(reads) -> String

Digest of the simulated fixture: per record in order, identifier, sequence, and
quality string, each newline-terminated. Matches the capture procedure used for
`INDEL_BENCH_TEST_FIXTURE_SHA256`.
"""
function _indel_bench_test_fixture_sha256(
        reads::Vector{FASTX.FASTQ.Record}
)::String
    buffer = IOBuffer()
    for record in reads
        print(buffer, FASTX.identifier(record))
        print(buffer, "\n")
        print(buffer, FASTX.sequence(String, record))
        print(buffer, "\n")
        print(buffer, String(collect(FASTX.quality(record))))
        print(buffer, "\n")
    end
    return Base.bytes2hex(SHA.sha256(Base.take!(buffer)))
end

Test.@testset "indel benchmark common" begin
    Test.@testset "golden fixture byte identity" begin
        reads,
        reference = indel_bench_simulate_reads(
            genome_length = INDEL_BENCH_TEST_GENOME_LENGTH,
            source_read_length = INDEL_BENCH_TEST_SOURCE_READ_LENGTH,
            coverage = INDEL_BENCH_TEST_COVERAGE,
            error_rate = INDEL_BENCH_TEST_ERROR_RATE,
            seed = INDEL_BENCH_TEST_SEED,
            tech = :nanopore
        )
        observed_lengths = length.(FASTX.sequence.(String, reads))
        Test.@test length(reads) == INDEL_BENCH_TEST_N_READS
        Test.@test minimum(observed_lengths) == INDEL_BENCH_TEST_MIN_READ_LENGTH
        Test.@test maximum(observed_lengths) == INDEL_BENCH_TEST_MAX_READ_LENGTH
        Test.@test length(reference) == INDEL_BENCH_TEST_REFERENCE_LENGTH
        Test.@test Base.bytes2hex(
            SHA.sha256(Vector{UInt8}(string(reference)))
        ) == INDEL_BENCH_TEST_REFERENCE_SHA256
        Test.@test _indel_bench_test_fixture_sha256(reads) ==
                   INDEL_BENCH_TEST_FIXTURE_SHA256
        # Identifier scheme is part of the digested bytes.
        Test.@test FASTX.identifier(first(reads)) == "nanopore_read_1"
    end

    Test.@testset "best reference alignment metric semantics" begin
        # A short contig against a much longer reference — the shape the
        # td-jt7r fixed-toy proof actually produces. These assertions are the
        # evidence for the `identity` -> `best_contig_reference_coverage`
        # rename: they pin BOTH that the coverage ratio is contiguity and that
        # the global match count it is built from is saturated. The reference
        # is drawn from a LOCAL RNG so it is aperiodic; no digest is pinned
        # here, only structure.
        reference_rng = Random.MersenneTwister(20_260_727)
        reference_string = String(rand(reference_rng, "ACGT", 2_000))
        contig = reference_string[401:600]
        reference = BioSequences.LongDNA{4}(reference_string)

        alignment = indel_bench_best_reference_alignment([contig], reference)

        # The global alignment pads the short contig out to the reference
        # length, so `aligned_bases` is the REFERENCE length.
        Test.@test alignment.orientation === :forward
        Test.@test alignment.matches == 200
        Test.@test alignment.aligned_bases == 2_000
        Test.@test alignment.contig_length == 200

        # Coverage therefore degenerates to contig_length / reference_length —
        # the exact shape of the published 0.1 (200 bp) and 0.0655 (131 bp)
        # figures on a 2 kb reference.
        Test.@test alignment.best_contig_reference_coverage ≈ 200 / 2_000
        Test.@test alignment.best_contig_reference_coverage ≈ 0.1

        # The contig is flawless, and the fit alignment says so. One assembly,
        # 0.1 by the coverage ratio and 1.0 by the accuracy ratio: this pair is
        # why the coverage ratio must not be called `identity`.
        Test.@test alignment.best_contig_fit_identity ≈ 1.0
        Test.@test alignment.best_contig_fit_matches == 200
        Test.@test alignment.best_contig_fit_edit_distance == 0

        # SATURATION. A contig carrying a substitution still scores a full 200
        # global matches, because the global alignment is free to match the
        # wrong base somewhere else in the reference. The global `matches` and
        # `edit_distance` columns therefore carry no accuracy information at
        # this contig/reference length ratio; only the fit alignment moves.
        substituted_base = contig[100] == 'A' ? 'C' : 'A'
        mutated = string(contig[1:99], substituted_base, contig[101:end])
        mutated_alignment = indel_bench_best_reference_alignment(
            [mutated], reference
        )
        Test.@test mutated_alignment.matches == 200
        Test.@test mutated_alignment.best_contig_reference_coverage ≈ 0.1
        Test.@test mutated_alignment.best_contig_fit_matches == 199
        Test.@test mutated_alignment.best_contig_fit_edit_distance == 1
        Test.@test mutated_alignment.best_contig_fit_identity ≈ 199 / 200

        # Indels inside the contig cost exactly their length in the fit
        # alignment, so it reports unit-cost (Levenshtein) identity.
        deleted = string(contig[1:99], contig[103:end])
        deleted_alignment = indel_bench_best_reference_alignment(
            [deleted], reference
        )
        Test.@test deleted_alignment.best_contig_fit_edit_distance == 3
        Test.@test deleted_alignment.best_contig_fit_matches == 197
        Test.@test deleted_alignment.best_contig_fit_identity ≈ 197 / 200

        # Saturation extends to STRAND: a reverse-strand contig ties its own
        # reverse complement at the coverage maximum, so the global
        # alignment records `:forward` for a contig that is on the reverse
        # strand. The fit identity must not inherit that arbitrary choice —
        # it is recomputed over both orientations, so it still reports 1.0.
        reverse_contig = string(
            BioSequences.reverse_complement(BioSequences.LongDNA{4}(contig))
        )
        reverse_alignment = indel_bench_best_reference_alignment(
            [reverse_contig], reference
        )
        Test.@test reverse_alignment.orientation === :forward
        Test.@test reverse_alignment.matches == 200
        Test.@test reverse_alignment.best_contig_reference_coverage ≈ 0.1
        Test.@test reverse_alignment.best_contig_fit_identity ≈ 1.0
        Test.@test reverse_alignment.best_contig_fit_edit_distance == 0

        # Selection ranks contiguity: a longer flawless contig wins on coverage
        # while both are 1.0 accurate. This is the long-standing selection
        # rule, unchanged by the rename.
        both = indel_bench_best_reference_alignment(
            [contig, reference_string[401:900]], reference
        )
        Test.@test both.contig_length == 500
        Test.@test both.best_contig_reference_coverage ≈ 0.25
        Test.@test both.best_contig_fit_identity ≈ 1.0

        # Empty assembly returns the all-zero `:none` sentinel, with no
        # division by zero in either ratio.
        empty_alignment = indel_bench_best_reference_alignment(
            String[], reference
        )
        Test.@test empty_alignment.orientation === :none
        Test.@test empty_alignment.best_contig_reference_coverage == 0.0
        Test.@test empty_alignment.best_contig_fit_identity == 0.0
        Test.@test empty_alignment.contig_length == 0

        # The fit helper is well defined on an empty contig too.
        Test.@test indel_bench_best_contig_fit_alignment(
            "", reference_string
        ).identity == 0.0
    end

    Test.@testset "n50" begin
        Test.@test indel_bench_n50(String[]) == 0
        Test.@test indel_bench_n50(["AAAA"]) == 4
        # lengths 5,4,3,2,1 => total 15, half 8 => 5+4=9 >= 8 => N50 = 4
        Test.@test indel_bench_n50([
            "AAAAA", "CCCC", "GGG", "TT", "A"
        ]) == 4
    end

    Test.@testset "rung accessors" begin
        symbol_rung = Dict(:requested => 3, :attempted => 0)
        string_rung = Dict("requested" => 3)
        tuple_rung = (requested = 3, attempted = 0)

        Test.@test indel_bench_rung_value(symbol_rung, :requested, -1) == 3
        Test.@test indel_bench_rung_value(string_rung, :requested, -1) == 3
        Test.@test indel_bench_rung_value(tuple_rung, :requested, -1) == 3
        Test.@test indel_bench_rung_value(symbol_rung, :missing_key, -1) == -1
        Test.@test indel_bench_rung_value(42, :requested, -1) == -1

        Test.@test indel_bench_rung_has_key(symbol_rung, :requested)
        Test.@test indel_bench_rung_has_key(string_rung, :requested)
        Test.@test indel_bench_rung_has_key(tuple_rung, :attempted)
        Test.@test !indel_bench_rung_has_key(symbol_rung, :missing_key)
        Test.@test !indel_bench_rung_has_key(42, :requested)

        Test.@test indel_bench_rung_counter(symbol_rung, :requested) == 3
        Test.@test indel_bench_rung_counter(symbol_rung, :attempted) == 0
        Test.@test isnothing(indel_bench_rung_counter(symbol_rung, :missing_key))
        # Non-Int and negative counters must not be laundered into totals.
        Test.@test isnothing(
            indel_bench_rung_counter(Dict(:requested => 3.0), :requested)
        )
        Test.@test isnothing(
            indel_bench_rung_counter(Dict(:requested => -1), :requested)
        )
    end

    Test.@testset "hashing and dependency provenance" begin
        scratch = Base.Filesystem.mktempdir()
        try
            present = joinpath(scratch, "present.txt")
            write(present, "indel-bench")
            Test.@test indel_bench_file_sha256(present) ==
                       Base.bytes2hex(SHA.sha256(Vector{UInt8}("indel-bench")))
            Test.@test indel_bench_optional_dependency_sha256(present) ==
                       indel_bench_file_sha256(present)
            Test.@test indel_bench_optional_dependency_sha256(
                joinpath(scratch, "absent.txt")
            ) == INDEL_BENCH_MISSING_DEPENDENCY_SENTINEL
        finally
            Base.rm(scratch; recursive = true, force = true)
        end

        dependency = indel_bench_dependency_provenance()
        Test.@test isfile(dependency.active_project_path)
        Test.@test length(dependency.project_toml_sha256) == 64
    end

    Test.@testset "code environment fingerprint" begin
        environment = indel_bench_code_environment(@__FILE__)
        Test.@test length(environment.code_environment_fingerprint) == 64
        Test.@test length(environment.code_sha) == 40
        Test.@test environment.benchmark_source_sha256 ==
                   indel_bench_file_sha256(@__FILE__)
        Test.@test environment.julia_version == string(VERSION)
        # Recomputing without changing anything must agree.
        Test.@test isnothing(
            indel_bench_assert_environment_unchanged(
            environment, @__FILE__; context = "self-check"
        ),
        )
        # A tampered fingerprint must be rejected, and the context must surface.
        tampered = merge(
            environment, (code_environment_fingerprint = repeat("0", 64),)
        )
        Test.@test_throws ErrorException indel_bench_assert_environment_unchanged(
            tampered, @__FILE__; context = "tamper-check"
        )
    end

    Test.@testset "publish preserves prior generation on failure" begin
        artifact_names = ("data.csv", "checks.csv", "manifest.csv")
        output_dir = Base.Filesystem.mktempdir()
        try
            prior = Dict(
                "data.csv" => "prior data\n",
                "checks.csv" => "prior checks\n",
                "manifest.csv" => "prior manifest\n"
            )
            for (name, contents) in prior
                write(joinpath(output_dir, name), contents)
            end
            prior_digests = Dict(
                name => indel_bench_file_sha256(joinpath(output_dir, name))
            for name in artifact_names
            )

            # A run that aborts mid-staging: stage two of three artifacts, throw
            # before publish, clean the staging dir in `finally` exactly as the
            # benchmark scripts do.
            staging_dir = Base.Filesystem.mktempdir(
                output_dir; prefix = ".unit-staging-"
            )
            aborted = false
            try
                write(joinpath(staging_dir, "data.csv"), "new data\n")
                write(joinpath(staging_dir, "checks.csv"), "new checks\n")
                # The real scripts throw here when an acceptance check fails or
                # a decode errors — before the manifest is staged and long
                # before `indel_bench_publish_artifacts` is reached.
                error("simulated mid-staging abort before the manifest")
            catch exception
                aborted = exception isa ErrorException
            finally
                Base.rm(staging_dir; recursive = true, force = true)
            end
            Test.@test aborted
            for name in artifact_names
                path = joinpath(output_dir, name)
                Test.@test isfile(path)
                Test.@test indel_bench_file_sha256(path) == prior_digests[name]
            end

            # An incomplete staging dir handed directly to publish must error
            # BEFORE touching the prior generation.
            partial_dir = Base.Filesystem.mktempdir(
                output_dir; prefix = ".unit-staging-"
            )
            try
                write(joinpath(partial_dir, "data.csv"), "new data\n")
                Test.@test_throws ErrorException indel_bench_publish_artifacts(
                    partial_dir, output_dir, artifact_names
                )
            finally
                Base.rm(partial_dir; recursive = true, force = true)
            end
            for name in artifact_names
                path = joinpath(output_dir, name)
                Test.@test isfile(path)
                Test.@test indel_bench_file_sha256(path) == prior_digests[name]
            end

            # A complete staging dir replaces every artifact.
            complete_dir = Base.Filesystem.mktempdir(
                output_dir; prefix = ".unit-staging-"
            )
            try
                for name in artifact_names
                    write(joinpath(complete_dir, name), "new $(name)\n")
                end
                indel_bench_publish_artifacts(
                    complete_dir, output_dir, artifact_names
                )
            finally
                Base.rm(complete_dir; recursive = true, force = true)
            end
            for name in artifact_names
                Test.@test read(joinpath(output_dir, name), String) ==
                           "new $(name)\n"
            end
        finally
            Base.rm(output_dir; recursive = true, force = true)
        end
    end

    Test.@testset "remove prior artifacts" begin
        output_dir = Base.Filesystem.mktempdir()
        try
            artifact_names = ("data.csv", "manifest.csv")
            for name in artifact_names
                write(joinpath(output_dir, name), "x")
            end
            indel_bench_remove_prior_artifacts(output_dir, artifact_names)
            Test.@test !any(
                isfile(joinpath(output_dir, name)) for name in artifact_names
            )
            # Idempotent on an already-clean directory, and creates it.
            fresh_dir = joinpath(output_dir, "nested", "fresh")
            indel_bench_remove_prior_artifacts(fresh_dir, artifact_names)
            Test.@test isdir(fresh_dir)
            # A directory squatting on an artifact name is refused, not deleted.
            Base.Filesystem.mkpath(joinpath(output_dir, "data.csv"))
            Test.@test_throws ErrorException indel_bench_remove_prior_artifacts(
                output_dir, artifact_names
            )
            Test.@test isdir(joinpath(output_dir, "data.csv"))
        finally
            Base.rm(output_dir; recursive = true, force = true)
        end
    end

    Test.@testset "signed rank" begin
        # All-positive n=5: V = 15, exact two-sided p = 2/2^5 = 0.0625 (matches
        # R's wilcox.test(1:5) exact result).
        all_positive = indel_bench_signed_rank([1.0, 2.0, 3.0, 4.0, 5.0])
        Test.@test all_positive.n_nonzero == 5
        Test.@test all_positive.statistic == 15.0
        Test.@test all_positive.method == :exact
        Test.@test all_positive.p_two_sided ≈ 0.0625

        # Sign symmetry: negating every delta flips V to the complementary rank
        # sum and leaves the two-sided p-value unchanged.
        all_negative = indel_bench_signed_rank([-1.0, -2.0, -3.0, -4.0, -5.0])
        Test.@test all_negative.statistic == 0.0
        Test.@test all_negative.p_two_sided ≈ all_positive.p_two_sided

        # Zeros are dropped, not counted.
        with_zeros = indel_bench_signed_rank(
            [0.0, 1.0, 0.0, 2.0, 3.0, 4.0, 5.0]
        )
        Test.@test with_zeros.n_nonzero == 5
        Test.@test with_zeros.statistic == 15.0
        Test.@test with_zeros.p_two_sided ≈ all_positive.p_two_sided

        # Mixed signs, n=3: V = 5, exact two-sided p = 0.5 (R: V=5, p=0.5).
        mixed = indel_bench_signed_rank([-1.0, 2.0, 3.0])
        Test.@test mixed.n_nonzero == 3
        Test.@test mixed.statistic == 5.0
        Test.@test mixed.p_two_sided ≈ 0.5

        # Tied magnitudes get midranks; every attainable rank sum is at least as
        # extreme as the observed one here, so p = 1.
        tied = indel_bench_signed_rank([1.0, 1.0, -1.0])
        Test.@test tied.n_nonzero == 3
        Test.@test tied.statistic == 4.0
        Test.@test tied.method == :exact
        Test.@test tied.p_two_sided ≈ 1.0

        # No nonzero differences => undefined, not a spurious p-value.
        empty_result = indel_bench_signed_rank(Float64[])
        Test.@test empty_result.n_nonzero == 0
        Test.@test empty_result.method == :undefined
        Test.@test isnan(empty_result.p_two_sided)
        all_zero = indel_bench_signed_rank([0.0, 0.0, 0.0])
        Test.@test all_zero.n_nonzero == 0
        Test.@test all_zero.method == :undefined
        Test.@test isnan(all_zero.p_two_sided)

        # Exact/approximate switch is at n_nonzero == 20.
        Test.@test indel_bench_signed_rank(collect(1.0:20.0)).method == :exact
        approximate = indel_bench_signed_rank(collect(1.0:21.0))
        Test.@test approximate.n_nonzero == 21
        Test.@test approximate.statistic == 231.0
        Test.@test approximate.method == :normal_approximation
        Test.@test 0.0 < approximate.p_two_sided < 1.0e-3

        # A balanced large sample is nowhere near significant.
        balanced = indel_bench_signed_rank(
            Float64[iseven(index) ? index : -index for index in 1:30]
        )
        Test.@test balanced.method == :normal_approximation
        Test.@test balanced.p_two_sided > 0.5

        # p-values stay in [0, 1] for both branches.
        Test.@test 0.0 <= indel_bench_signed_rank([1.0, -1.0]).p_two_sided <= 1.0
        Test.@test 0.0 <=
                   indel_bench_signed_rank(fill(1.0, 25)).p_two_sided <= 1.0

        # Non-finite deltas are rejected rather than silently dropped.
        Test.@test_throws ArgumentError indel_bench_signed_rank([NaN, 1.0])
        Test.@test_throws ArgumentError indel_bench_signed_rank([Inf, 1.0])
    end
end
