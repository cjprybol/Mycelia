import Test
import Mycelia

Test.@testset "Checkpointing Tests" begin
    # Use a temporary directory for all checkpoint tests
    mktempdir() do checkpoint_dir
        Test.@testset "cached_stage basic round-trip (no input_files)" begin
            call_count = Ref(0)
            result = Mycelia.cached_stage("test_basic", checkpoint_dir) do
                call_count[] += 1
                Dict("key" => "value", "n" => 42)
            end
            Test.@test result == Dict("key" => "value", "n" => 42)
            Test.@test call_count[] == 1

            # Second call should load from cache, not recompute
            result2 = Mycelia.cached_stage("test_basic", checkpoint_dir) do
                call_count[] += 1
                Dict("key" => "different")
            end
            Test.@test result2 == Dict("key" => "value", "n" => 42)
            Test.@test call_count[] == 1  # not incremented
        end

        Test.@testset "cached_stage with input_files" begin
            mktempdir() do input_dir
                input_file = joinpath(input_dir, "data.txt")
                write(input_file, "original content")

                call_count = Ref(0)
                result = Mycelia.cached_stage("test_inputs", checkpoint_dir;
                    input_files = [input_file]) do
                    call_count[] += 1
                    "computed_from_original"
                end
                Test.@test result == "computed_from_original"
                Test.@test call_count[] == 1

                # Second call with same input should use cache
                result2 = Mycelia.cached_stage("test_inputs", checkpoint_dir;
                    input_files = [input_file]) do
                    call_count[] += 1
                    "should_not_see_this"
                end
                Test.@test result2 == "computed_from_original"
                Test.@test call_count[] == 1
            end
        end

        Test.@testset "cache invalidation when input file changes" begin
            mktempdir() do input_dir
                input_file = joinpath(input_dir, "data.txt")
                write(input_file, "v1")

                call_count = Ref(0)
                result1 = Mycelia.cached_stage("test_invalidation", checkpoint_dir;
                    input_files = [input_file]) do
                    call_count[] += 1
                    "result_v1"
                end
                Test.@test result1 == "result_v1"
                Test.@test call_count[] == 1

                # Modify the input file (changes size, and mtime via sleep)
                sleep(1.1)  # ensure mtime differs on 1-second resolution filesystems
                write(input_file, "version 2 with different length content")

                result2 = Mycelia.cached_stage("test_invalidation", checkpoint_dir;
                    input_files = [input_file]) do
                    call_count[] += 1
                    "result_v2"
                end
                Test.@test result2 == "result_v2"
                Test.@test call_count[] == 2  # recomputed due to file change
            end
        end

        Test.@testset "stale cache cleanup on recompute" begin
            mktempdir() do input_dir
                input_file = joinpath(input_dir, "data.txt")
                write(input_file, "v1")

                Mycelia.cached_stage("test_stale", checkpoint_dir;
                    input_files = [input_file]) do
                    "result_v1"
                end
                # Count hash-suffixed files for this stage
                stale_files_before = filter(
                    f -> startswith(f, "test_stale_") && endswith(f, ".jld2"),
                    readdir(checkpoint_dir))
                Test.@test length(stale_files_before) == 1

                # Change input and recompute
                sleep(1.1)
                write(input_file, "v2 longer content")

                Mycelia.cached_stage("test_stale", checkpoint_dir;
                    input_files = [input_file]) do
                    "result_v2"
                end
                # Old cache should be cleaned up, only new one remains
                stale_files_after = filter(
                    f -> startswith(f, "test_stale_") && endswith(f, ".jld2"),
                    readdir(checkpoint_dir))
                Test.@test length(stale_files_after) == 1
                Test.@test stale_files_before != stale_files_after
            end
        end

        Test.@testset "clear_stage" begin
            Mycelia.cached_stage("to_clear", checkpoint_dir) do
                "data"
            end
            Test.@test Mycelia.clear_stage("to_clear", checkpoint_dir) == true
            Test.@test Mycelia.clear_stage("nonexistent", checkpoint_dir) == false
        end

        Test.@testset "clear_stage with input_files" begin
            mktempdir() do input_dir
                input_file = joinpath(input_dir, "data.txt")
                write(input_file, "content for clear test")

                Mycelia.cached_stage("to_clear_hashed", checkpoint_dir;
                    input_files = [input_file]) do
                    "hashed_data"
                end

                # Base name clear should NOT find the hashed file
                Test.@test Mycelia.clear_stage("to_clear_hashed", checkpoint_dir) == false

                # Clearing with correct input_files should work
                Test.@test Mycelia.clear_stage("to_clear_hashed", checkpoint_dir;
                    input_files = [input_file]) == true

                # Should be gone now
                Test.@test Mycelia.clear_stage("to_clear_hashed", checkpoint_dir;
                    input_files = [input_file]) == false
            end
        end

        Test.@testset "list_stages" begin
            Mycelia.cached_stage("stage_a", checkpoint_dir) do
                1
            end
            Mycelia.cached_stage("stage_b", checkpoint_dir) do
                2
            end
            stages = Mycelia.list_stages(checkpoint_dir)
            Test.@test "stage_a" in stages
            Test.@test "stage_b" in stages

            # Non-existent directory returns empty
            Test.@test Mycelia.list_stages("/nonexistent/path") == String[]
        end

        Test.@testset "clear_all_stages" begin
            mktempdir() do tmpdir
                Mycelia.cached_stage("s1", tmpdir) do
                    1
                end
                Mycelia.cached_stage("s2", tmpdir) do
                    2
                end
                Test.@test Mycelia.clear_all_stages(tmpdir) == 2
                Test.@test Mycelia.list_stages(tmpdir) == String[]

                # Non-existent directory returns 0
                Test.@test Mycelia.clear_all_stages("/nonexistent/path") == 0
            end
        end

        Test.@testset "checkpoint_info" begin
            mktempdir() do tmpdir
                Mycelia.cached_stage("info_test", tmpdir) do
                    collect(1:1000)
                end
                info = Mycelia.checkpoint_info(tmpdir)
                Test.@test size(info, 1) == 1
                Test.@test info.name[1] == "info_test"
                Test.@test info.size_mb[1] >= 0.0

                # Non-existent directory returns empty DataFrame
                empty_info = Mycelia.checkpoint_info("/nonexistent/path")
                Test.@test size(empty_info, 1) == 0
            end
        end

        Test.@testset "_input_hash determinism" begin
            mktempdir() do input_dir
                f1 = joinpath(input_dir, "a.txt")
                f2 = joinpath(input_dir, "b.txt")
                write(f1, "hello")
                write(f2, "world")

                h1 = Mycelia._input_hash([f1, f2])
                h2 = Mycelia._input_hash([f2, f1])  # reversed order
                Test.@test h1 == h2  # sorted internally, order-independent
                Test.@test length(h1) == 16
            end
        end

        Test.@testset "_input_hash errors on missing file" begin
            Test.@test_throws ErrorException Mycelia._input_hash(["/nonexistent/file.txt"])
        end

        Test.@testset "alternative signature passes input_files through" begin
            mktempdir() do input_dir
                input_file = joinpath(input_dir, "alt.txt")
                write(input_file, "test data")

                call_count = Ref(0)
                compute_fn = () -> begin
                    call_count[] += 1
                    "alt_result"
                end
                result = Mycelia.cached_stage("test_alt", checkpoint_dir, compute_fn;
                    input_files = [input_file])
                Test.@test result == "alt_result"
                Test.@test call_count[] == 1
            end
        end
    end
end
