# Pure-Julia staleness and non-destructive staging tests for taxonkit taxdump setup.
#
# `setup_taxonkit_taxonomy` previously warned and returned without updating when
# an existing `~/.taxonkit` directory lacked `.last_updated`, had an unreadable
# `.last_updated`, or exceeded `max_age_days`. It also removed the existing
# directory before attempting to download `taxdump.tar.gz`, leaving an empty
# directory if the download failed.
#
# These tests run strictly in temporary directories with synthetic populate
# callbacks: they never touch `~/.taxonkit`, never require network access, and
# never invoke `wget`.
#
# Run:
#   julia --project=. test/1_data_acquisition/taxonkit_staleness_test.jl

module TaxonkitStalenessTest

import Dates
import Logging
import Mycelia
import Test

function populate_dummy_taxdump!(
        dir::AbstractString;
        marker::AbstractString = "v1",
        timestamp::Union{Nothing, AbstractString} = nothing
)
    mkpath(dir)
    write(joinpath(dir, "nodes.dmp"), "1\t|\t1\t|\tno rank\t|\t$(marker)\n")
    write(joinpath(dir, "names.dmp"), "1\t|\troot\t|\t\t|\tscientific name\t|\t$(marker)\n")
    if timestamp !== nothing
        write(joinpath(dir, ".last_updated"), timestamp)
    end
    return dir
end

function assert_throws_with_message(
        f::Function,
        ::Type{E},
        expected_fragment::AbstractString
) where {E <: Exception}
    try
        f()
        Test.@test false
    catch err
        Test.@test err isa E
        Test.@test occursin(expected_fragment, sprint(showerror, err))
    end
    return nothing
end

Test.@testset "Taxonkit taxonomy staleness and staging" begin
    fixed_now = Dates.DateTime(2026, 11, 16, 12, 0, 0)

    Test.@testset "_taxonkit_update_reason classification" begin
        mktempdir() do root
            # Absent directory -> :missing
            missing_dir = joinpath(root, "absent_taxonkit")
            Test.@test Mycelia._taxonkit_update_reason(
                missing_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :missing

            # Empty directory -> :empty
            empty_dir = joinpath(root, "empty_taxonkit")
            mkpath(empty_dir)
            Test.@test Mycelia._taxonkit_update_reason(
                empty_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :empty

            # (a) Populated directory with no .last_updated -> :no_timestamp
            undated_dir = populate_dummy_taxdump!(
                joinpath(root, "undated_taxonkit");
                marker = "legacy"
            )
            Test.@test Mycelia._taxonkit_update_reason(
                undated_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :no_timestamp

            # Garbage / blank .last_updated -> :unreadable_timestamp
            garbage_dir = populate_dummy_taxdump!(
                joinpath(root, "garbage_taxonkit");
                timestamp = "not-an-iso-timestamp"
            )
            Test.@test Mycelia._taxonkit_update_reason(
                garbage_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :unreadable_timestamp

            blank_ts_dir = populate_dummy_taxdump!(
                joinpath(root, "blank_timestamp_taxonkit");
                timestamp = "   \n"
            )
            Test.@test Mycelia._taxonkit_update_reason(
                blank_ts_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :unreadable_timestamp

            # (b) Populated directory with timestamp older than max_age_days -> :stale
            stale_dir = populate_dummy_taxdump!(
                joinpath(root, "stale_taxonkit");
                timestamp = string(fixed_now - Dates.Day(31))
            )
            Test.@test Mycelia._taxonkit_update_reason(
                stale_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :stale

            just_over_dir = populate_dummy_taxdump!(
                joinpath(root, "just_over_taxonkit");
                timestamp = string(fixed_now - Dates.Day(30) - Dates.Second(1))
            )
            Test.@test Mycelia._taxonkit_update_reason(
                just_over_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :stale

            # Exact max_age_days boundary and recent timestamps -> :fresh
            exact_boundary_dir = populate_dummy_taxdump!(
                joinpath(root, "exact_boundary_taxonkit");
                timestamp = string(fixed_now - Dates.Day(30))
            )
            Test.@test Mycelia._taxonkit_update_reason(
                exact_boundary_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :fresh

            fresh_dir = populate_dummy_taxdump!(
                joinpath(root, "fresh_taxonkit");
                timestamp = string(fixed_now - Dates.Day(5)) * "\n"
            )
            Test.@test Mycelia._taxonkit_update_reason(
                fresh_dir;
                force_update = false,
                max_age_days = 30,
                now = fixed_now
            ) === :fresh
            Test.@test Mycelia._taxonkit_update_reason(
                fresh_dir;
                force_update = false,
                max_age_days = 3,
                now = fixed_now
            ) === :stale

            # force_update=true overrides even a fresh installation -> :forced
            Test.@test Mycelia._taxonkit_update_reason(
                fresh_dir;
                force_update = true,
                max_age_days = 30,
                now = fixed_now
            ) === :forced

            # Negative max_age_days is rejected with a descriptive ArgumentError
            assert_throws_with_message(ArgumentError, "max_age_days") do
                Mycelia._taxonkit_update_reason(
                    fresh_dir;
                    force_update = false,
                    max_age_days = -1,
                    now = fixed_now
                )
            end
        end
    end

    Test.@testset "(a) populated dir without .last_updated triggers update instead of warning" begin
        mktempdir() do root
            taxonkit_dir = populate_dummy_taxdump!(
                joinpath(root, "undated");
                marker = "old-undated"
            )
            write(joinpath(taxonkit_dir, "obsolete_extra.dmp"), "should disappear on swap\n")

            update_calls = Ref(0)
            fake_populate = staging_dir -> begin
                update_calls[] += 1
                populate_dummy_taxdump!(staging_dir; marker = "refreshed-from-undated")
                return nothing
            end

            result = Test.@test_logs (
                :info, r"missing \.last_updated"
            ) (
                :info, r"updated successfully"
            ) begin
                Mycelia._setup_taxonkit_taxonomy(
                    fake_populate;
                    taxonkit_dir = taxonkit_dir,
                    max_age_days = 30,
                    now = fixed_now
                )
            end

            Test.@test result === nothing
            Test.@test update_calls[] == 1
            Test.@test occursin(
                "refreshed-from-undated",
                read(joinpath(taxonkit_dir, "nodes.dmp"), String)
            )
            Test.@test !isfile(joinpath(taxonkit_dir, "obsolete_extra.dmp"))
            Test.@test read(joinpath(taxonkit_dir, ".last_updated"), String) == string(fixed_now)
            Test.@test readdir(root) == ["undated"]

            # Immediately re-running at the same timestamp is a no-op (:fresh)
            logger = Test.TestLogger()
            Logging.with_logger(logger) do
                Mycelia._setup_taxonkit_taxonomy(
                    fake_populate;
                    taxonkit_dir = taxonkit_dir,
                    max_age_days = 30,
                    now = fixed_now
                )
            end
            Test.@test update_calls[] == 1
            Test.@test isempty(logger.logs)
        end
    end

    Test.@testset "(b) populated dir older than max_age_days triggers update instead of warning" begin
        mktempdir() do root
            stale_timestamp = string(fixed_now - Dates.Day(45))
            taxonkit_dir = populate_dummy_taxdump!(
                joinpath(root, "stale");
                marker = "old-stale",
                timestamp = stale_timestamp
            )

            update_calls = Ref(0)
            fake_populate = staging_dir -> begin
                update_calls[] += 1
                populate_dummy_taxdump!(staging_dir; marker = "refreshed-from-stale")
                return nothing
            end

            result = Test.@test_logs (
                :info, r"exceeds max_age_days=30"
            ) (
                :info, r"updated successfully"
            ) begin
                Mycelia._setup_taxonkit_taxonomy(
                    fake_populate;
                    taxonkit_dir = taxonkit_dir,
                    max_age_days = 30,
                    now = fixed_now
                )
            end

            Test.@test result === nothing
            Test.@test update_calls[] == 1
            Test.@test occursin(
                "refreshed-from-stale",
                read(joinpath(taxonkit_dir, "nodes.dmp"), String)
            )
            Test.@test read(joinpath(taxonkit_dir, ".last_updated"), String) == string(fixed_now)
            Test.@test readdir(root) == ["stale"]
        end
    end

    Test.@testset "unreadable timestamp, empty, missing, and forced paths update; fresh public call skips download" begin
        mktempdir() do root
            # Unreadable timestamp updates
            unreadable_dir = populate_dummy_taxdump!(
                joinpath(root, "unreadable");
                marker = "old",
                timestamp = "corrupt"
            )
            unreadable_calls = Ref(0)
            Test.@test_logs (
                :info, r"unreadable \.last_updated"
            ) (
                :info, r"updated successfully"
            ) begin
                Mycelia._setup_taxonkit_taxonomy(
                    staging_dir -> begin
                        unreadable_calls[] += 1
                        populate_dummy_taxdump!(staging_dir; marker = "fixed")
                    end;
                    taxonkit_dir = unreadable_dir,
                    max_age_days = 30,
                    now = fixed_now
                )
            end
            Test.@test unreadable_calls[] == 1
            Test.@test occursin("fixed", read(joinpath(unreadable_dir, "nodes.dmp"), String))

            # Empty directory updates
            empty_dir = joinpath(root, "empty_existing")
            mkpath(empty_dir)
            empty_calls = Ref(0)
            Test.@test_logs (
                :info, r"directory is empty"
            ) (
                :info, r"updated successfully"
            ) begin
                Mycelia._setup_taxonkit_taxonomy(
                    staging_dir -> begin
                        empty_calls[] += 1
                        populate_dummy_taxdump!(staging_dir; marker = "filled")
                    end;
                    taxonkit_dir = empty_dir,
                    max_age_days = 30,
                    now = fixed_now
                )
            end
            Test.@test empty_calls[] == 1
            Test.@test occursin("filled", read(joinpath(empty_dir, "nodes.dmp"), String))

            # Unexpected reason symbol is rejected cleanly
            assert_throws_with_message(ArgumentError, "Unexpected taxonkit update reason") do
                Mycelia._taxonkit_update_message(:invalid_reason, empty_dir; max_age_days = 30)
            end

            # Missing directory updates (including trailing slash handling)
            missing_dir = joinpath(root, "newly_created")
            missing_calls = Ref(0)
            Test.@test_logs (
                :info, r"directory does not exist"
            ) (
                :info, r"updated successfully"
            ) begin
                Mycelia._setup_taxonkit_taxonomy(
                    staging_dir -> begin
                        missing_calls[] += 1
                        populate_dummy_taxdump!(staging_dir; marker = "created")
                    end;
                    taxonkit_dir = missing_dir * "/",
                    max_age_days = 30,
                    now = fixed_now
                )
            end
            Test.@test missing_calls[] == 1
            Test.@test isfile(joinpath(missing_dir, "nodes.dmp"))
            Test.@test read(joinpath(missing_dir, ".last_updated"), String) == string(fixed_now)

            # force_update=true updates even when fresh
            forced_calls = Ref(0)
            Test.@test_logs (
                :info, r"force_update=true"
            ) (
                :info, r"updated successfully"
            ) begin
                Mycelia._setup_taxonkit_taxonomy(
                    staging_dir -> begin
                        forced_calls[] += 1
                        populate_dummy_taxdump!(staging_dir; marker = "forced-refresh")
                    end;
                    force_update = true,
                    taxonkit_dir = missing_dir,
                    max_age_days = 30,
                    now = fixed_now
                )
            end
            Test.@test forced_calls[] == 1
            Test.@test occursin("forced-refresh", read(joinpath(missing_dir, "nodes.dmp"), String))

            # Public setup_taxonkit_taxonomy on a fresh temp dir returns without invoking wget
            live_fresh_dir = populate_dummy_taxdump!(
                joinpath(root, "public_fresh");
                marker = "keep-me",
                timestamp = string(Dates.now())
            )
            Test.@test Mycelia.setup_taxonkit_taxonomy(
                taxonkit_dir = live_fresh_dir,
                max_age_days = 30
            ) === nothing
            Test.@test occursin("keep-me", read(joinpath(live_fresh_dir, "nodes.dmp"), String))
        end
    end

    Test.@testset "non-destructive staging preserves existing taxdump on failure" begin
        mktempdir() do root
            original_ts = string(fixed_now - Dates.Day(60))
            taxonkit_dir = populate_dummy_taxdump!(
                joinpath(root, "existing");
                marker = "precious-existing-data",
                timestamp = original_ts
            )

            # Failure during download/extraction after partial write to staging dir
            assert_throws_with_message(ErrorException, "simulated download failure") do
                Mycelia._stage_and_swap_taxonkit_dir(
                    taxonkit_dir;
                    now = fixed_now
                ) do staging_dir
                    write(joinpath(staging_dir, "partial.dmp"), "incomplete")
                    error("simulated download failure")
                end
            end

            Test.@test isdir(taxonkit_dir)
            Test.@test occursin(
                "precious-existing-data",
                read(joinpath(taxonkit_dir, "nodes.dmp"), String)
            )
            Test.@test !isfile(joinpath(taxonkit_dir, "partial.dmp"))
            Test.@test read(joinpath(taxonkit_dir, ".last_updated"), String) == original_ts
            Test.@test readdir(root) == ["existing"]

            # Empty extraction is also rejected without touching existing data
            assert_throws_with_message(ErrorException, "empty directory") do
                Mycelia._stage_and_swap_taxonkit_dir(
                    taxonkit_dir;
                    now = fixed_now
                ) do _staging_dir
                    return nothing
                end
            end

            Test.@test occursin(
                "precious-existing-data",
                read(joinpath(taxonkit_dir, "nodes.dmp"), String)
            )
            Test.@test read(joinpath(taxonkit_dir, ".last_updated"), String) == original_ts
            Test.@test readdir(root) == ["existing"]
        end
    end
end

end  # module
