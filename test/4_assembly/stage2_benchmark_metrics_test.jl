import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_stage2_toy_benchmark.jl"))

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

    Test.@test_throws ErrorException parse_dnadiff_discrepancy(
        "TotalBases 8000 8000\nAvgIdentity 100.0")
end

Test.@testset "Stage-2 streaming artifact hash" begin
    mktemp() do path, stream
        write(stream, "abc")
        close(stream)
        Test.@test file_sha256(path) == bytes2hex(SHA.sha256("abc"))
    end
end

Test.@testset "Stage-2 staging ownership and promotion" begin
    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = prepare_stage2_staging!(final)
        Test.@test staging == stage2_output_paths(final).staging
        Test.@test isfile(stage2_staging_sentinel_path(staging))
        write(joinpath(staging, "stale.txt"), "stale")

        replaced = prepare_stage2_staging!(final)
        Test.@test replaced == staging
        Test.@test !isfile(joinpath(staging, "stale.txt"))
        Test.@test isfile(stage2_staging_sentinel_path(staging))

        write(joinpath(staging, "artifact.txt"), "complete")
        promoted = promote_stage2_staging!(final)
        Test.@test promoted == abspath(final)
        Test.@test read(joinpath(final, "artifact.txt"), String) == "complete"
        Test.@test !ispath(staging)
        Test.@test_throws ErrorException prepare_stage2_staging!(final)
        Test.@test read(joinpath(final, "artifact.txt"), String) == "complete"
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = stage2_output_paths(final).staging
        mkpath(staging)
        unowned = joinpath(staging, "unowned.txt")
        write(unowned, "preserve")
        Test.@test_throws ErrorException prepare_stage2_staging!(final)
        Test.@test read(unowned, String) == "preserve"
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = prepare_stage2_staging!(final)
        sentinel = stage2_staging_sentinel_path(staging)
        table = TOML.parsefile(sentinel)
        table["final_output_root"] = joinpath(directory, "different-final")
        open(sentinel, "w") do stream
            TOML.print(stream, table)
        end
        protected = joinpath(staging, "protected.txt")
        write(protected, "preserve")
        Test.@test_throws ErrorException prepare_stage2_staging!(final)
        Test.@test read(protected, String) == "preserve"
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = prepare_stage2_staging!(final)
        write(joinpath(staging, "artifact.txt"), "staged")
        mkpath(final)
        final_artifact = joinpath(final, "existing.txt")
        write(final_artifact, "preserve")
        Test.@test_throws ErrorException promote_stage2_staging!(final)
        Test.@test read(final_artifact, String) == "preserve"
        Test.@test isdir(staging)
    end

    mktempdir() do directory
        final = joinpath(directory, "final")
        staging = stage2_output_paths(final).staging
        target = joinpath(directory, "target")
        mkpath(target)
        symlink(target, staging)
        Test.@test_throws ErrorException prepare_stage2_staging!(final)
        Test.@test isdir(target)
    end
end
