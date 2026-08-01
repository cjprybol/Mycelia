# OPT4 CONFIGURATION GUARDS (td-jbjd, review fixes C2 + C3)
# =======================================================================
#
# Two review-flagged gaps in the opt4 (frozen-read skip) kwarg surface:
#
#   C2 — checkpoint resume-config gap: `_iterative_checkpoint_configuration`
#     did not capture `skip_frozen_reads`/`freeze_streak_threshold`/
#     `freeze_across_rungs`, so resuming a checkpointed run under DIFFERENT
#     opt4 settings silently switched an exact run into an approximate one (or
#     vice versa) instead of tripping `_validate_checkpoint_configuration`'s
#     exact-equality mismatch guard the way every other behavior-affecting
#     kwarg (skip_solid, hard_window, soft_em, cheap_correct, ...) already
#     does. Fixed by adding all three to the captured config dict + its call
#     site in `mycelia_iterative_assemble`.
#
#   C3 — `freeze_streak_threshold=0` would make `freeze_streaks[i] >=
#     freeze_streak_threshold` true for EVERY read from pass 1 onward (streaks
#     start all-zero), silently disabling the corrector entirely whenever
#     `skip_frozen_reads=true`. Fixed with an `ArgumentError` guard requiring
#     `freeze_streak_threshold >= 1`.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_config_guards_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _OPT4CG_BASES = ['A', 'C', 'G', 'T']

function _opt4cg_reads(rng, ref; n_reads = 40, readlen = 40, n_err = 8)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _OPT4CG_BASES))
        end
        push!(records,
            FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

Test.@testset "opt4 configuration guards (td-jbjd, C2 + C3)" begin
    Test.@testset "C2: checkpoint resume rejects a skip_frozen_reads flip" begin
        rng = Random.MersenneTwister(31337)
        ref = join(rand(rng, _OPT4CG_BASES, 400))
        reads = _opt4cg_reads(rng, ref; n_reads = 40, readlen = 40, n_err = 8)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)
        output_directory = joinpath(tmp, "run")

        common = (; max_k = 13, max_iterations_per_k = 1,
            improvement_threshold = 0.0, stop_on_no_change = false,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            verbose = false, enable_checkpointing = true,
            checkpoint_interval = 1, output_dir = output_directory)

        # First invocation lays down a checkpoint with skip_frozen_reads=false.
        first_result = Mycelia.mycelia_iterative_assemble(
            fq; common..., skip_frozen_reads = false)
        Test.@test isfile(joinpath(
            output_directory, "checkpoints", "latest_checkpoint.json"))
        Test.@test first_result[:metadata][:corrector_errors][
            :frozen_reads_skipped] == 0

        # Re-invoking against the SAME output_dir with every other kwarg held
        # identical but skip_frozen_reads flipped to true must be rejected --
        # NOT silently proceed as if freezing had been the setting all along
        # (which would retroactively switch an already-exact checkpointed run
        # into an approximate one).
        observed_exception = try
            Mycelia.mycelia_iterative_assemble(
                fq; common..., skip_frozen_reads = true)
            nothing
        catch exception
            exception
        end
        Test.@test observed_exception isa ArgumentError
        Test.@test occursin(
            "checkpoint resume_configuration does not match this invocation",
            sprint(showerror, observed_exception))
        Test.@test occursin("skip_frozen_reads", sprint(showerror, observed_exception))

        # The mirror direction (checkpoint saved WITH freezing, resumed
        # WITHOUT) must be rejected the same way -- the guard is symmetric,
        # not a one-directional "only tighten" check.
        output_directory_2 = joinpath(tmp, "run2")
        common2 = merge(common, (; output_dir = output_directory_2))
        Mycelia.mycelia_iterative_assemble(
            fq; common2..., skip_frozen_reads = true, freeze_streak_threshold = 2)
        observed_exception_2 = try
            Mycelia.mycelia_iterative_assemble(
                fq; common2..., skip_frozen_reads = false)
            nothing
        catch exception
            exception
        end
        Test.@test observed_exception_2 isa ArgumentError
        Test.@test occursin(
            "checkpoint resume_configuration does not match this invocation",
            sprint(showerror, observed_exception_2))
    end

    Test.@testset "C3: freeze_streak_threshold must be >= 1" begin
        rng = Random.MersenneTwister(99)
        ref = join(rand(rng, _OPT4CG_BASES, 300))
        reads = _opt4cg_reads(rng, ref; n_reads = 20, readlen = 30, n_err = 5)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        for bad_threshold in (0, -1, -5)
            observed_exception = try
                Mycelia.mycelia_iterative_assemble(
                    fq; max_k = 13, max_iterations_per_k = 1,
                    graph_mode = :canonical, verbose = false,
                    enable_checkpointing = false,
                    output_dir = joinpath(tmp, "bad_$(bad_threshold)"),
                    skip_frozen_reads = true,
                    freeze_streak_threshold = bad_threshold)
                nothing
            catch exception
                exception
            end
            Test.@test observed_exception isa ArgumentError
            Test.@test occursin(
                "freeze_streak_threshold must be at least 1",
                sprint(showerror, observed_exception))
        end

        # The guard fires independent of skip_frozen_reads too -- a caller that
        # later flips skip_frozen_reads=true without revisiting an already-bad
        # threshold=0 must not be silently broken.
        observed_exception_disabled = try
            Mycelia.mycelia_iterative_assemble(
                fq; max_k = 13, max_iterations_per_k = 1,
                graph_mode = :canonical, verbose = false,
                enable_checkpointing = false,
                output_dir = joinpath(tmp, "disabled_zero"),
                skip_frozen_reads = false, freeze_streak_threshold = 0)
            nothing
        catch exception
            exception
        end
        Test.@test observed_exception_disabled isa ArgumentError

        # Positive control: threshold=1 (the minimum valid value) does not
        # raise.
        result = Mycelia.mycelia_iterative_assemble(
            fq; max_k = 13, max_iterations_per_k = 1,
            graph_mode = :canonical, verbose = false,
            enable_checkpointing = false,
            output_dir = joinpath(tmp, "min_valid"),
            skip_frozen_reads = true, freeze_streak_threshold = 1)
        Test.@test isfile(result[:metadata][:final_fastq_file])
    end
end
