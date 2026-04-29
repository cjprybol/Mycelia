import Test
import Mycelia
import JLD2

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

    Test.@testset "cached_map" begin
    # Counters are `Threads.Atomic{Int}` because `cached_map` runs `fn` via
    # `Threads.@threads`; non-atomic read-modify-write (`Ref{Int}` / `x[] += 1`)
    # loses increments under multi-threaded CI.
    Test.@testset "first call computes all inputs" begin
        mktempdir() do dir
            call_count = Threads.Atomic{Int}(0)
            inputs = [1, 2, 3]
            results = Mycelia.cached_map("first_call", dir, inputs) do x
                Threads.atomic_add!(call_count, 1)
                2x
            end
            Test.@test results == [2, 4, 6]
            Test.@test call_count[] == 3
            Test.@test isfile(joinpath(dir, "first_call.jld2"))
            Test.@test !isfile(joinpath(dir, "first_call.jld2.lock"))
        end
    end

    Test.@testset "second call with same inputs: no recomputation" begin
        mktempdir() do dir
            call_count = Threads.Atomic{Int}(0)
            inputs = [10, 20, 30]
            fn = x -> begin
                Threads.atomic_add!(call_count, 1)
                x + 1
            end
            first = Mycelia.cached_map("memoized", dir, inputs, fn)
            Test.@test first == [11, 21, 31]
            Test.@test call_count[] == 3

            second = Mycelia.cached_map("memoized", dir, inputs, fn)
            Test.@test second == [11, 21, 31]
            Test.@test call_count[] == 3  # unchanged: all hits
        end
    end

    Test.@testset "incremental growth only computes new inputs" begin
        mktempdir() do dir
            call_count = Threads.Atomic{Int}(0)
            fn = x -> begin
                Threads.atomic_add!(call_count, 1)
                2x
            end
            first = Mycelia.cached_map("grow", dir, [1, 2, 3], fn)
            Test.@test first == [2, 4, 6]
            Test.@test call_count[] == 3

            second = Mycelia.cached_map("grow", dir, [1, 2, 3, 4, 5], fn)
            Test.@test second == [2, 4, 6, 8, 10]
            Test.@test call_count[] == 5  # only 4 and 5 were new
        end
    end

    Test.@testset "custom keyfn supports non-string inputs" begin
        mktempdir() do dir
            call_count = Threads.Atomic{Int}(0)
            # Tuple inputs: key by first element only, second element is metadata
            inputs = [("a", "meta-1"), ("b", "meta-2"), ("c", "meta-3")]
            fn = x -> begin
                Threads.atomic_add!(call_count, 1)
                uppercase(x[1]) * "|" * x[2]
            end
            first = Mycelia.cached_map("tuples", dir, inputs, fn; keyfn = x -> x[1])
            Test.@test first == ["A|meta-1", "B|meta-2", "C|meta-3"]
            Test.@test call_count[] == 3

            # Add a 4th tuple; existing keys hit cache, new key computes
            inputs2 = [("a", "meta-1"), ("b", "meta-2"), ("c", "meta-3"), ("d", "meta-4")]
            second = Mycelia.cached_map("tuples", dir, inputs2, fn; keyfn = x -> x[1])
            Test.@test second == ["A|meta-1", "B|meta-2", "C|meta-3", "D|meta-4"]
            Test.@test call_count[] == 4
        end
    end

    Test.@testset "force=true ignores cache and recomputes" begin
        mktempdir() do dir
            call_count = Threads.Atomic{Int}(0)
            fn = x -> begin
                Threads.atomic_add!(call_count, 1)
                x * 10
            end
            Mycelia.cached_map("forced", dir, [1, 2, 3], fn)
            Test.@test call_count[] == 3

            Mycelia.cached_map("forced", dir, [1, 2, 3], fn; force = true)
            Test.@test call_count[] == 6  # all three recomputed
        end
    end

    Test.@testset "corrupted cache file rebuilds gracefully" begin
        mktempdir() do dir
            cache_file = joinpath(dir, "corrupt.jld2")
            write(cache_file, "not a valid jld2 file, just garbage bytes")

            call_count = Threads.Atomic{Int}(0)
            # @test_logs asserts the warning is emitted and captures no crash
            results = Test.@test_logs (:warn,) match_mode = :any begin
                Mycelia.cached_map("corrupt", dir, [1, 2, 3]) do x
                    Threads.atomic_add!(call_count, 1)
                    x + 100
                end
            end
            Test.@test results == [101, 102, 103]
            Test.@test call_count[] == 3

            # On-disk cache must be a valid JLD2 file after the rebuild —
            # otherwise a later run would re-trigger the rebuild path and
            # recompute everything, defeating the whole point of the cache.
            reloaded = JLD2.load(cache_file, "data")
            Test.@test reloaded isa Dict
            Test.@test length(reloaded) == 3

            # Subsequent call should now hit the rebuilt cache
            results2 = Mycelia.cached_map("corrupt", dir, [1, 2, 3]) do x
                Threads.atomic_add!(call_count, 1)
                error("should not be called")
            end
            Test.@test results2 == [101, 102, 103]
            Test.@test call_count[] == 3
        end
    end

    Test.@testset "empty inputs returns empty result" begin
        mktempdir() do dir
            call_count = Threads.Atomic{Int}(0)
            results = Mycelia.cached_map("empty", dir, Int[]) do x
                Threads.atomic_add!(call_count, 1)
                x
            end
            Test.@test isempty(results)
            Test.@test call_count[] == 0
        end
    end

    Test.@testset "partial progress preserved when fn throws" begin
        # This is the core invariant of cached_map: if `fn` throws partway
        # through a long sweep, every successfully-computed entry must be on
        # disk so the next run picks up where we left off. Without the
        # exception handler inside cached_map, the final _atomic_save_dict
        # would be skipped on exception and all in-memory progress would be lost.
        mktempdir() do dir
            inputs = [1, 2, 3, 4, 5]
            # Throw on input 4. With `Threads.@threads`, work is partitioned
            # into per-thread chunks: when one worker throws mid-chunk, the
            # remaining items in that chunk never run. So the exact set of
            # cached entries depends on thread scheduling — we can only
            # guarantee that the failing key is absent and that *some* progress
            # was flushed.
            Test.@test_throws Exception Mycelia.cached_map(
                "partial", dir, inputs; save_every = 2) do x
                x == 4 ? error("simulated failure on input $x") : x * 100
            end

            cache_file = joinpath(dir, "partial.jld2")
            Test.@test isfile(cache_file)
            reloaded = JLD2.load(cache_file, "data")
            Test.@test reloaded isa Dict
            # The failing input must never be cached (its `fn` threw before
            # the assignment `cached[key] = result`).
            Test.@test !haskey(reloaded, "4")
            # The finally block must flush whatever *did* land in the dict —
            # at minimum one successful entry, at most four.
            Test.@test 1 <= length(reloaded) <= 4
            # Every cached entry must be correct (no partial/garbage values).
            for (k, v) in reloaded
                Test.@test v == parse(Int, k) * 100
            end

            # A subsequent successful run with a non-throwing fn fills in
            # whatever was missing and returns a complete result.
            retry_count = Threads.Atomic{Int}(0)
            missing_before_retry = 5 - length(reloaded)
            results = Mycelia.cached_map("partial", dir, inputs) do x
                Threads.atomic_add!(retry_count, 1)
                x * 100
            end
            Test.@test results == [100, 200, 300, 400, 500]
            # Only the missing entries are recomputed (at minimum input 4,
            # plus any inputs whose chunk-mate threw before they ran).
            Test.@test retry_count[] == missing_before_retry
        end
    end
end
end
