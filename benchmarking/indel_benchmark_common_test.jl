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
import DataFrames
import FASTX
import Random
import SHA
import Test

# indel_benchmark_common.jl is DELIBERATELY not included directly here: the
# runtime include below already pulls it in (indel_frontier_runtime.jl:41), and
# the shared file has no include guard, so a direct include would re-evaluate
# every `const` in it in this same scope with a freshly allocated value. Bindings
# like INDEL_BENCH_MISSING_VALUE, the `@__FILE__` path, and the run-status
# filename are not `===` across two evaluations, so Julia emits a
# const-redefinition warning for each — noise that can bury a genuine warning in
# the test output. Every `indel_bench_*` symbol asserted below is defined
# transitively; do not re-add the direct include.
#
# `_indel_frontier_label_offsets` is a pure, deterministic figure-layout helper
# that lives in indel_frontier_runtime.jl and had no coverage. That script's
# `main()` is guarded by the `PROGRAM_FILE` check, so including it defines
# constants and functions and runs nothing; the packages it imports are already
# loaded by Mycelia. Covering it here avoids a second test entry point that
# nothing would invoke.
include(joinpath(@__DIR__, "indel_frontier_runtime.jl"))
# The fixed-toy proof's gate/advisory split and exit decision had no coverage at
# all, and the evidence offered for downgrading its wall-clock check to an
# advisory was an 18/18 HEALTHY run — the one observation that cannot tell "the
# wall-clock check stopped being fatal" apart from "the correctness detectors
# were switched off". Including the script defines its constants and helpers and
# runs nothing, because its `main()` sits behind the same `PROGRAM_FILE` guard
# relied on for the runtime include above.
#
# It goes in a module because both scripts define `main`; a bare second include
# would silently overwrite the first one's method. Qualifying the calls also
# keeps it obvious which entry points are under test.
module IndelToyProof
include(joinpath(@__DIR__, "indel_frontier_fixed_toy_proof.jl"))
end

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

    Test.@testset "simulated read identifiers follow tech" begin
        # `identifier_prefix` used to default to a hardcoded "nanopore_read"
        # independently of `tech`, so an :illumina fixture silently carried
        # nanopore identifiers into the digested bytes.
        illumina_reads, _ = indel_bench_simulate_reads(
            genome_length = 200,
            source_read_length = 100,
            coverage = 2,
            error_rate = 0.01,
            seed = 7,
            tech = :illumina
        )
        Test.@test !isempty(illumina_reads)
        Test.@test all(
            startswith(FASTX.identifier(record), "illumina_read_")
            for record in illumina_reads
        )

        # The default is derived, not removed: an explicit prefix still wins.
        explicit_reads, _ = indel_bench_simulate_reads(
            genome_length = 200,
            source_read_length = 100,
            coverage = 2,
            error_rate = 0.01,
            seed = 7,
            tech = :illumina,
            identifier_prefix = "custom"
        )
        Test.@test all(
            startswith(FASTX.identifier(record), "custom_")
            for record in explicit_reads
        )

        # The :nanopore default is unchanged, so the golden digest above holds.
        nanopore_reads, _ = indel_bench_simulate_reads(
            genome_length = 200,
            source_read_length = 100,
            coverage = 2,
            error_rate = 0.01,
            seed = 7
        )
        Test.@test startswith(
            FASTX.identifier(first(nanopore_reads)), "nanopore_read_"
        )
    end

    Test.@testset "start-of-run precheck and run status sidecar" begin
        artifact_names = ("data.csv", "manifest.csv")
        output_dir = Base.Filesystem.mktempdir()
        try
            # Nothing recorded yet.
            Test.@test isnothing(indel_bench_read_run_status(output_dir))

            # begin_run marks the directory as having a run in flight.
            indel_bench_begin_run(
                output_dir, artifact_names; context = "unit-run"
            )
            running = indel_bench_read_run_status(output_dir)
            Test.@test !isnothing(running)
            Test.@test running["status"] == "running"
            Test.@test running["context"] == "unit-run"
            Test.@test running["generation_id"] == "pending"

            # Publishing advances the marker to complete and records the
            # generation, so an aborted run is distinguishable from a finished
            # one by inspection of the output directory alone.
            staging_dir = Base.Filesystem.mktempdir(
                output_dir; prefix = ".unit-staging-"
            )
            try
                for name in artifact_names
                    write(joinpath(staging_dir, name), "new $(name)\n")
                end
                indel_bench_publish_artifacts(
                    staging_dir,
                    output_dir,
                    artifact_names;
                    context = "unit-run",
                    generation_id = "abc123"
                )
            finally
                Base.rm(staging_dir; recursive = true, force = true)
            end
            complete = indel_bench_read_run_status(output_dir)
            Test.@test complete["status"] == "complete"
            Test.@test complete["generation_id"] == "abc123"
            Test.@test indel_bench_prior_generation_id(complete) == "abc123"

            # A second run in flight must not erase the identifier of the
            # generation currently on disk — that identifier is the whole point
            # of the publish-time overwrite report.
            indel_bench_begin_run(
                output_dir, artifact_names; context = "unit-rerun"
            )
            rerunning = indel_bench_read_run_status(output_dir)
            Test.@test rerunning["generation_id"] == "pending"
            Test.@test rerunning["previous_generation_id"] == "abc123"
            Test.@test indel_bench_prior_generation_id(rerunning) == "abc123"
            Test.@test indel_bench_prior_generation_id(nothing) == "none"
            Test.@test all(
                isfile(joinpath(output_dir, name)) for name in artifact_names
            )

            # The squatting-directory guard is now a START-of-run failure, and
            # it is non-destructive.
            Base.rm(joinpath(output_dir, "data.csv"))
            Base.Filesystem.mkpath(joinpath(output_dir, "data.csv"))
            Test.@test_throws ErrorException indel_bench_assert_publishable(
                output_dir, artifact_names
            )
            Test.@test_throws ErrorException indel_bench_begin_run(
                output_dir, artifact_names; context = "unit-run"
            )
            Test.@test isdir(joinpath(output_dir, "data.csv"))
            Test.@test isfile(joinpath(output_dir, "manifest.csv"))
        finally
            Base.rm(output_dir; recursive = true, force = true)
        end

        # A publish with BOTH descriptive keywords omitted is a supported call
        # shape (the failure-durability testset above uses it). It must still
        # advance the marker: the status field records whether the generation
        # published, and leaving it at `running` would make a directory holding a
        # complete generation indistinguishable from one whose run aborted --
        # defeating the sidecar's only purpose. The descriptive fields fall back
        # to "unknown" rather than gating the advance.
        keywordless_dir = Base.Filesystem.mktempdir()
        try
            indel_bench_begin_run(
                keywordless_dir, artifact_names; context = "unit-keywordless"
            )
            Test.@test indel_bench_read_run_status(
                keywordless_dir
            )["status"] == "running"
            staging_dir = Base.Filesystem.mktempdir(
                keywordless_dir; prefix = ".unit-staging-"
            )
            try
                for name in artifact_names
                    write(joinpath(staging_dir, name), "new $(name)\n")
                end
                indel_bench_publish_artifacts(
                    staging_dir, keywordless_dir, artifact_names
                )
            finally
                Base.rm(staging_dir; recursive = true, force = true)
            end
            published = indel_bench_read_run_status(keywordless_dir)
            Test.@test published["status"] == "complete"
            Test.@test published["context"] == "unknown"
            Test.@test published["generation_id"] == "unknown"
            Test.@test haskey(published, "completed_at")
        finally
            Base.rm(keywordless_dir; recursive = true, force = true)
        end

        # A never-created output directory is publishable by definition.
        Test.@test isnothing(
            indel_bench_assert_publishable(
            joinpath(Base.Filesystem.mktempdir(), "absent"), artifact_names
        ),
        )
    end

    Test.@testset "frontier figure label offsets" begin
        base = INDEL_FRONTIER_LABEL_BASE_OFFSET
        step = INDEL_FRONTIER_LABEL_LEVEL_OFFSET

        Test.@test _indel_frontier_label_offsets(
            Float64[], Float64[], Int[], Bool[]
        ) == Int[]

        # Mismatched input lengths are rejected rather than silently truncated.
        Test.@test_throws ArgumentError _indel_frontier_label_offsets(
            [1.0, 2.0], [1.0], [4, 4], [true, true]
        )

        # Two identical anchors with identical wide labels overlap: the first
        # placed stays at the base offset, the second is raised by whole levels.
        collided = _indel_frontier_label_offsets(
            [10.0, 10.0], [1.0, 1.0], [20, 20], [true, true]
        )
        Test.@test collided[1] == base
        Test.@test collided[2] > collided[1]
        Test.@test (collided[2] - base) % step == 0
        # Exactly ONE level, not two. A level bump must move a label by at least
        # its own collision height, or the first bump never clears the test and
        # every collision silently costs two levels of vertical space.
        Test.@test collided[2] == base + step

        # Well-separated anchors with short labels do not collide.
        separated = _indel_frontier_label_offsets(
            [10.0, 10_000.0], [1.0, 1.0], [4, 4], [true, true]
        )
        Test.@test separated == [base, base]

        # Alignment is load-bearing: a right-aligned label extends LEFT from its
        # anchor, so opposite alignments at one anchor occupy disjoint boxes.
        opposed = _indel_frontier_label_offsets(
            [10.0, 10.0], [1.0, 1.0], [20, 20], [false, true]
        )
        Test.@test opposed == [base, base]

        # Placement is top-down: of two colliding labels the HIGHER anchor is
        # placed first and keeps the base offset, and the lower one is pushed up
        # into empty space rather than the higher one being pushed onto it. The
        # third, far-above point sets the y span so the first two land close
        # enough in normalized space to collide.
        stacked = _indel_frontier_label_offsets(
            [10.0, 10.0, 10.0], [1.0, 1.01, 100.0], [20, 20, 20],
            [true, true, true]
        )
        Test.@test stacked[3] == base
        Test.@test stacked[2] == base
        Test.@test stacked[1] > base

        # Deterministic — the figure digest depends on it.
        anchors_x = Float64[10.0^(index % 4) for index in 1:12]
        anchors_y = Float64[1.0 + (index % 3) for index in 1:12]
        lengths = Int[10 + index for index in 1:12]
        alignments = Bool[isodd(index) for index in 1:12]
        first_pass = _indel_frontier_label_offsets(
            anchors_x, anchors_y, lengths, alignments
        )
        Test.@test first_pass == _indel_frontier_label_offsets(
            anchors_x, anchors_y, lengths, alignments
        )
        Test.@test all(offset >= base for offset in first_pass)

        # Levels are capped, so a degenerate pile cannot run away.
        piled = _indel_frontier_label_offsets(
            fill(10.0, 30), fill(1.0, 30), fill(24, 30), fill(true, 30)
        )
        Test.@test maximum(piled) <=
                   base + INDEL_FRONTIER_LABEL_MAX_LEVELS * step
    end
end

# --- Positive control for the fixed-toy gate/advisory split -----------------
#
# The synthetic arms below never assemble anything: `_indel_toy_checks` reads
# only counters, telemetry rows, and byte strings off the arm NamedTuples, so
# the whole control runs in milliseconds and can therefore live in the ordinary
# unit-test entry point rather than behind the expensive proof.

"""
    _indel_toy_test_reads() -> Vector{FASTX.FASTQ.Record}

Two reads long enough to clear `INDEL_TOY_MIN_REQUIRED_READ_LENGTH`. Only the
lengths are read by the checks, so the bases are arbitrary.
"""
function _indel_toy_test_reads()::Vector{FASTX.FASTQ.Record}
    sequence = repeat("A", 1_200)
    quality = repeat("I", 1_200)
    return FASTX.FASTQ.Record[
        FASTX.FASTQ.Record("toy_read_$(index)", sequence, quality)
        for index in 1:2
    ]
end

"""
    _indel_toy_test_rung(k; requested, attempted, completed, truncated, engaged)

One per-rung telemetry record in the shape the proof's validator expects: exact
nonnegative `Int` counters keyed by `Symbol`.
"""
function _indel_toy_test_rung(
        k::Int;
        requested::Int,
        attempted::Int,
        completed::Int,
        truncated::Int,
        engaged::Int
)::Dict{Symbol, Any}
    return Dict{Symbol, Any}(
        :k => k,
        :requested => requested,
        :attempted => attempted,
        :completed => completed,
        :truncated => truncated,
        :engaged => engaged
    )
end

"""
    _indel_toy_test_arm(; ...) -> NamedTuple

Minimal stand-in for `_indel_toy_run_arm`'s return value carrying exactly the
fields `_indel_toy_checks` reads.
"""
function _indel_toy_test_arm(;
        label::String,
        wall_seconds::Float64,
        telemetry::Vector{Dict{Symbol, Any}},
        requested::Int,
        attempted::Int,
        completed::Int,
        truncated::Int,
        engaged::Int,
        reference_coverage::Float64,
        assembly_bytes::Vector{UInt8}
)::NamedTuple
    return (
        label = label,
        wall_seconds = wall_seconds,
        best_contig_reference_coverage = reference_coverage,
        best_contig_fit_identity = 0.99,
        best_contig_length = 200,
        n_contigs = 1,
        k_progression = [11, 17],
        telemetry = telemetry,
        indel_requested = requested,
        indel_attempted = attempted,
        indel_completed = completed,
        indel_truncated = truncated,
        indel_engaged = engaged,
        trace_contract_errors = 0,
        window_divergences = 0,
        assembly_bytes = assembly_bytes
    )
end

"""
    _indel_toy_test_outcome(; wall_seconds, oracle_matches) -> NamedTuple

Run the real `_indel_toy_checks` + `_indel_toy_evaluate_checks` over a synthetic
fixture in which every correctness gate passes, then perturb exactly two inputs:
the nanopore wall clock (the check that was downgraded to an advisory) and the
pre-wiring oracle digest (a correctness gate that must NOT have been affected).

`oracle_matches = false` is injected through the `prewiring_sha256` seam rather
than by mutating the assembly bytes, because sha256 is not invertible: a
fabricated byte stream can never satisfy the committed golden digest, so a
healthy baseline is otherwise unconstructible.
"""
function _indel_toy_test_outcome(;
        wall_seconds::Float64,
        oracle_matches::Bool
)::NamedTuple
    reads = _indel_toy_test_reads()
    substitution_only_bytes = Vector{UInt8}(">toy_contig\nACGTACGTACGT\n")
    nanopore = _indel_toy_test_arm(
        label = "nanopore",
        wall_seconds = wall_seconds,
        telemetry = Dict{Symbol, Any}[
            _indel_toy_test_rung(
                11; requested = 0, attempted = 0, completed = 0,
                truncated = 0, engaged = 0
            ),
            _indel_toy_test_rung(
                17; requested = 4, attempted = 4, completed = 4,
                truncated = 0, engaged = 2
            )
        ],
        requested = 4,
        attempted = 4,
        completed = 4,
        truncated = 0,
        engaged = 2,
        reference_coverage = 0.1,
        assembly_bytes = Vector{UInt8}(">toy_contig\nACGTACGTACGTACGT\n")
    )
    zero_rung = Dict{Symbol, Any}[
        _indel_toy_test_rung(
            11; requested = 0, attempted = 0, completed = 0,
            truncated = 0, engaged = 0
        )
    ]
    illumina = _indel_toy_test_arm(
        label = "illumina",
        wall_seconds = 1.0,
        telemetry = zero_rung,
        requested = 0,
        attempted = 0,
        completed = 0,
        truncated = 0,
        engaged = 0,
        reference_coverage = 0.0655,
        assembly_bytes = substitution_only_bytes
    )
    oracle = _indel_toy_test_arm(
        label = "default_illumina_oracle",
        wall_seconds = 1.0,
        telemetry = zero_rung,
        requested = 0,
        attempted = 0,
        completed = 0,
        truncated = 0,
        engaged = 0,
        reference_coverage = 0.0655,
        assembly_bytes = substitution_only_bytes
    )
    prewiring_sha256 = oracle_matches ?
                       Base.bytes2hex(SHA.sha256(substitution_only_bytes)) :
                       repeat("0", 64)
    checks = IndelToyProof._indel_toy_checks(
        reads, nanopore, illumina, oracle; prewiring_sha256 = prewiring_sha256
    )
    return (
        checks = checks,
        outcome = IndelToyProof._indel_toy_evaluate_checks(checks)
    )
end

Test.@testset "fixed-toy severity split and exit decision" begin
    Test.@testset "healthy run: all gates pass, nothing breached" begin
        result = _indel_toy_test_outcome(
            wall_seconds = 79.7, oracle_matches = true
        )
        Test.@test result.outcome.passed
        Test.@test isempty(result.outcome.failed)
        Test.@test isempty(result.outcome.breached)
        Test.@test result.outcome.gate_total ==
                   IndelToyProof.INDEL_TOY_EXPECTED_GATE_COUNT
        Test.@test result.outcome.advisory_total == 1
    end

    # The change this branch makes, and the ONLY case whose verdict it moves.
    # 141.2 s is the wall clock the reviewer actually measured on a loaded host.
    Test.@testset "wall breached, correctness intact: PASSES" begin
        result = _indel_toy_test_outcome(
            wall_seconds = 141.2, oracle_matches = true
        )
        Test.@test result.outcome.passed
        Test.@test isempty(result.outcome.failed)
        Test.@test result.outcome.breached == ["nanopore_under_120_seconds"]
    end

    # The downgrade must not have blinded any correctness detector.
    Test.@testset "oracle hash mismatch: still FAILS" begin
        result = _indel_toy_test_outcome(
            wall_seconds = 79.7, oracle_matches = false
        )
        Test.@test !result.outcome.passed
        Test.@test result.outcome.failed ==
                   ["illumina_byte_identical_to_prewiring_oracle"]
        Test.@test isempty(result.outcome.breached)
    end

    # The discriminating case: an advisory breach must not be able to absorb a
    # gate failure that happens at the same time, which is the failure mode a
    # severity split introduces and neither of the single-fault cases can see.
    Test.@testset "advisory breach cannot mask a concurrent gate failure" begin
        result = _indel_toy_test_outcome(
            wall_seconds = 141.2, oracle_matches = false
        )
        Test.@test !result.outcome.passed
        Test.@test result.outcome.failed ==
                   ["illumina_byte_identical_to_prewiring_oracle"]
        Test.@test result.outcome.breached == ["nanopore_under_120_seconds"]
    end

    Test.@testset "advisory set is pinned to the wall-clock check alone" begin
        result = _indel_toy_test_outcome(
            wall_seconds = 79.7, oracle_matches = true
        )
        checks = result.checks
        advisories = String[
            row.check for row in DataFrames.eachrow(checks)
            if row.severity == "advisory"
        ]
        Test.@test advisories == ["nanopore_under_120_seconds"]
        Test.@test IndelToyProof.INDEL_TOY_ADVISORY_CHECKS ==
                   ("nanopore_under_120_seconds",)
        Test.@test Set(checks.severity) == Set(["gate", "advisory"])
        Test.@test DataFrames.nrow(checks) ==
                   IndelToyProof.INDEL_TOY_EXPECTED_GATE_COUNT + 1
    end
end
