import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_stage2_toy_benchmark.jl"))

function _test_stage2_benchmark_error(
        f::Function,
        exception_type::Type{<:Exception},
        expected_message::AbstractString,
)::Nothing
    caught = try
        f()
        nothing
    catch error
        error
    end
    message = caught === nothing ? "" : sprint(showerror, caught)
    Test.@test caught isa exception_type
    Test.@test occursin(expected_message, message)
    return nothing
end

Test.@testset "Stage-2 dnadiff full-discrepancy metric" begin
    report = """
    [Bases]
    TotalBases 8000 8030
    AlignedBases 8000(100.00%) 8000(99.63%)
    [Alignments]
    AvgIdentity 99.50 99.50
    """
    metrics = parse_dnadiff_discrepancy(report)
    Test.@test metrics.aligned_pct_ref == 100.0
    Test.@test metrics.avg_identity == 99.5
    Test.@test metrics.total_errors_per_100kbp ≈ 875.0

    mktemp() do path, stream
        write(stream, report)
        close(stream)
        Test.@test parse_dnadiff_discrepancy_file(path) == metrics
    end

    _test_stage2_benchmark_error(
        ErrorException,
        "dnadiff report lacks reference aligned percentage",
    ) do
        parse_dnadiff_discrepancy(
            "TotalBases 8000 8000\nAvgIdentity 100.0")
    end
end

Test.@testset "Stage-2 streaming artifact hash" begin
    mktemp() do path, stream
        write(stream, "abc")
        close(stream)
        Test.@test file_sha256(path) == bytes2hex(SHA.sha256("abc"))
    end
end

Test.@testset "Stage-2 fixture complexity constraints" begin
    random_sequence = make_primary_sequence(
        StableRNGs.StableRNG(1);
        fixture_length = 128,
        fixture_complexity = :random,
    )
    Test.@test length(random_sequence) == 128

    short_repeat_message = try
        make_primary_sequence(
            StableRNGs.StableRNG(1);
            fixture_length = STAGE2_REPEAT_RICH_MINIMUM_LENGTH - 1,
            fixture_complexity = :repeat_rich,
        )
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin(
        "repeat-rich Stage-2 fixture length must be at least 5399",
        short_repeat_message,
    )
    repeat_rich_sequence = make_primary_sequence(
        StableRNGs.StableRNG(1);
        fixture_length = STAGE2_REPEAT_RICH_MINIMUM_LENGTH,
        fixture_complexity = :repeat_rich,
    )
    Test.@test length(repeat_rich_sequence) ==
               STAGE2_REPEAT_RICH_MINIMUM_LENGTH
    Test.@test STAGE2_TOY_EVIDENCE_STATUS == "unattested_contract_only"
end

Test.@testset "Stage-2 nested attempt provenance" begin
    source = Dict{String, Any}(
        "stage2_attempt_id" => "attempt-inner",
        "stage2_attempt_relative_path" => "attempt-inner",
        "stage2_attempt_path_semantics" =>
            "relative-to-config-output-dir",
        "nested" => Dict("status" => "experimental"),
    )
    payload = serialized_stage2_deconvolution_provenance(source)
    Test.@test payload["stage2_attempt_id"] == "attempt-inner"
    Test.@test payload["stage2_attempt_relative_path"] == "attempt-inner"
    Test.@test payload["stage2_attempt_path_semantics"] ==
               "relative-to-config-output-dir"
    payload["nested"]["status"] = "mutated"
    Test.@test source["nested"]["status"] == "experimental"

    _test_stage2_benchmark_error(
        ErrorException,
        "lacks an attempt ID",
    ) do
        serialized_stage2_deconvolution_provenance(Dict{String, Any}())
    end
    _test_stage2_benchmark_error(
        ErrorException,
        "path is not bound to its attempt ID",
    ) do
        serialized_stage2_deconvolution_provenance(Dict{String, Any}(
            "stage2_attempt_id" => "attempt-inner",
            "stage2_attempt_relative_path" => "attempt-other",
            "stage2_attempt_path_semantics" =>
                "relative-to-config-output-dir",
        ))
    end
    _test_stage2_benchmark_error(
        ErrorException,
        "invalid path semantics",
    ) do
        serialized_stage2_deconvolution_provenance(Dict{String, Any}(
            "stage2_attempt_id" => "attempt-inner",
            "stage2_attempt_relative_path" => "attempt-inner",
            "stage2_attempt_path_semantics" => "absolute",
        ))
    end
end

Test.@testset "Stage-2 staging ownership and promotion" begin
    mktempdir() do directory
        final = joinpath(directory, "final")
        paths = stage2_output_paths(final)
        attempt = prepare_stage2_staging!(final; attempt_id = "attempt-one")
        staging = attempt.staging_root
        Test.@test staging == paths.staging
        Test.@test attempt.lock_path == paths.lock
        Test.@test isfile(paths.lock)
        Test.@test isfile(stage2_staging_sentinel_path(staging))
        sentinel = TOML.parsefile(stage2_staging_sentinel_path(staging))
        Test.@test sentinel["schema"] == STAGE2_STAGING_SENTINEL_SCHEMA
        Test.@test sentinel["attempt_id"] == "attempt-one"
        Test.@test sentinel["lock_path"] == paths.lock
        Test.@test sentinel["lock_state"] == "active"
        write(joinpath(staging, "stale.txt"), "stale")

        active_lock_message = try
            prepare_stage2_staging!(final; attempt_id = "attempt-two")
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("active attempt lock", active_lock_message)
        Test.@test read(joinpath(staging, "stale.txt"), String) == "stale"

        release_stage2_staging_attempt!(attempt)
        Test.@test !attempt.active
        Test.@test !ispath(paths.lock)
        replacement = prepare_stage2_staging!(
            final; attempt_id = "attempt-two")
        Test.@test replacement.staging_root == staging
        Test.@test !isfile(joinpath(staging, "stale.txt"))
        Test.@test isfile(stage2_staging_sentinel_path(staging))
        replacement_sentinel = TOML.parsefile(
            stage2_staging_sentinel_path(staging))
        Test.@test replacement_sentinel["attempt_id"] == "attempt-two"
        Test.@test replacement_sentinel["lock_state"] == "active"

        write(joinpath(staging, "artifact.txt"), "complete")
        promoted = promote_stage2_staging!(replacement)
        Test.@test promoted == abspath(final)
        Test.@test read(joinpath(final, "artifact.txt"), String) == "complete"
        Test.@test !ispath(staging)
        Test.@test !replacement.active
        Test.@test !ispath(paths.lock)
        final_sentinel = TOML.parsefile(
            stage2_staging_sentinel_path(final))
        Test.@test final_sentinel["lock_state"] == "promoted"
        final_exists_message = try
            prepare_stage2_staging!(final; attempt_id = "attempt-three")
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin(
            "final output root already exists", final_exists_message)
        Test.@test !ispath(paths.lock)
        Test.@test read(joinpath(final, "artifact.txt"), String) == "complete"
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = stage2_output_paths(final).staging
        mkpath(staging)
        unowned = joinpath(staging, "unowned.txt")
        write(unowned, "preserve")
        unowned_message = try
            prepare_stage2_staging!(final; attempt_id = "unowned-attempt")
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("unowned Stage-2 staging root", unowned_message)
        Test.@test read(unowned, String) == "preserve"
        Test.@test !ispath(stage2_output_paths(final).lock)
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        attempt = prepare_stage2_staging!(
            final; attempt_id = "tampered-path-attempt")
        staging = attempt.staging_root
        release_stage2_staging_attempt!(attempt)
        sentinel = stage2_staging_sentinel_path(staging)
        table = TOML.parsefile(sentinel)
        table["final_output_root"] = joinpath(directory, "different-final")
        open(sentinel, "w") do stream
            TOML.print(stream, table)
        end
        protected = joinpath(staging, "protected.txt")
        write(protected, "preserve")
        tampered_path_message = try
            prepare_stage2_staging!(final; attempt_id = "next-attempt")
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin(
            "bound to a different final path", tampered_path_message)
        Test.@test read(protected, String) == "preserve"
        Test.@test !ispath(stage2_output_paths(final).lock)
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        attempt = prepare_stage2_staging!(
            final; attempt_id = "promotion-conflict")
        staging = attempt.staging_root
        write(joinpath(staging, "artifact.txt"), "staged")
        mkpath(final)
        final_artifact = joinpath(final, "existing.txt")
        write(final_artifact, "preserve")
        promotion_message = try
            promote_stage2_staging!(attempt)
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin(
            "final output root already exists", promotion_message)
        Test.@test read(final_artifact, String) == "preserve"
        Test.@test isdir(staging)
        Test.@test !attempt.active
        Test.@test !ispath(stage2_output_paths(final).lock)
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        attempt = prepare_stage2_staging!(
            final; attempt_id = "bound-attempt")
        sentinel_path = stage2_staging_sentinel_path(attempt.staging_root)
        sentinel = TOML.parsefile(sentinel_path)
        sentinel["attempt_id"] = "different-attempt"
        open(sentinel_path, "w") do stream
            TOML.print(stream, sentinel)
        end
        wrong_attempt_message = try
            promote_stage2_staging!(attempt)
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin(
            "attempt ID differs from the active attempt",
            wrong_attempt_message,
        )
        Test.@test isdir(attempt.staging_root)
        Test.@test !attempt.active
        Test.@test !ispath(attempt.lock_path)
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = stage2_output_paths(final).staging
        target = joinpath(directory, "target")
        mkpath(target)
        symlink(target, staging)
        symlink_message = try
            prepare_stage2_staging!(final; attempt_id = "symlink-attempt")
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("must not be a symbolic link", symlink_message)
        Test.@test isdir(target)
        Test.@test !ispath(stage2_output_paths(final).lock)
    end
end
