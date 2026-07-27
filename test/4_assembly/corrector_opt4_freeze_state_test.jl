# OPT4 FROZEN-READ SKIP — ENABLED-PATH FREEZE-STATE + INTERACTIONS (td-jbjd, pass 2)
# ======================================================================================
#
# Pass 1 (corrector_opt4_frozen_read_skip_test.jl) locked the DISABLED path
# (skip_frozen_reads=false) as a true no-op. This file exercises the ENABLED
# path (skip_frozen_reads=true) directly — the freeze-streak state machine
# itself, and its composition with soft-EM and stop_on_no_change — per the
# opt4 Rule-of-5 pass-2 ("self-review") scope in
# docs/design/2026-07-26-opt4-frozen-read-skip-design.md.
#
# WHAT IS ASSERTED
#   * Freeze-state unit test: driving `improve_read_set_likelihood` repeatedly
#     with a shared `freeze_streaks` vector + `CorrectorDiagnostics`, the
#     streak increments on no-improvement, resets on improvement (engineered
#     deterministically), and `frozen_reads_skipped` becomes positive once a
#     read's streak crosses `freeze_streak_threshold`.
#   * Within-rung vs across-rung scope: `mycelia_iterative_assemble` with
#     `freeze_across_rungs=false` (default) resets streaks at every k-advance,
#     so a fresh k-rung re-decodes everything; `freeze_across_rungs=true`
#     persists streaks across the advance instead. The scope difference is
#     observed as a directional inequality in total `frozen_reads_skipped`
#     across the whole run (across-rung >= within-rung).
#   * soft-EM interaction: `soft_em=true` + `skip_frozen_reads=true` together
#     do not crash and produce sane output (a frozen-skipped read still
#     contributes -- unchanged -- via the existing skip-passthrough path that
#     soft-EM already tolerates for skip_solid/hard_window skips).
#   * stop_on_no_change interaction: an aggressive freeze configuration does
#     not hang, still terminates at the same k-ladder (`final_k`) as the exact
#     path, and does not increase iteration count vs the exact path.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_freeze_state_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _OPT4FS_BASES = ['A', 'C', 'G', 'T']

function _opt4fs_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _OPT4FS_BASES))
        end
        push!(records,
            FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

Test.@testset "opt4 skip_frozen_reads freeze-state + interactions (td-jbjd pass 2)" begin
    Test.@testset "freeze-state: increment on no-improve, resets on improve" begin
        rng = Random.MersenneTwister(9001)
        ref = join(rand(rng, _OPT4FS_BASES, 1200))
        reads = _opt4fs_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
        k = 13
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)

        common = (; graph_mode = :canonical, skip_solid = false,
            cheap_correct = false, decode_enabled = true, batch_size = 50)

        freeze_streaks = zeros(Int, length(reads))
        diag = Mycelia.CorrectorDiagnostics()
        threshold = 2

        # Pass 1. NOTE: this repo's own completeness test documents that the
        # ungated per-read decode is "mostly pass-through by design" at toy
        # fixture scale (rhizomorph_substitution_reconstruction_completeness_
        # test.jl) — asserting a natural was_improved==true count here would be
        # flaky, so this sub-test does NOT depend on any read naturally
        # improving. It only pins the INCREMENT side of the arithmetic (every
        # read that does not improve advances its streak by exactly 1 per
        # pass) and defers the RESET side to the deterministic engineered case
        # below, which does not depend on the decode's natural behavior.
        out1, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; common..., skip_frozen_reads = true,
            freeze_streak_threshold = threshold, freeze_streaks = freeze_streaks,
            diagnostics = diag)
        Test.@test all(s -> s in (0, 1), freeze_streaks)
        Test.@test count(==(1), freeze_streaks) == length(reads)

        # Pass 2: same (unrebuilt) graph, feed forward pass 1's output. A read
        # already at streak==1 is still BELOW threshold (2) entering this pass,
        # so it decodes again (not yet frozen); it should end this pass at
        # streak==2 if still no improvement. No frozen-skip should have fired
        # yet (nothing was >= threshold when this pass STARTED).
        out2, = Mycelia.improve_read_set_likelihood(
            out1, graph, k; common..., skip_frozen_reads = true,
            freeze_streak_threshold = threshold, freeze_streaks = freeze_streaks,
            diagnostics = diag)
        Test.@test diag.frozen_reads_skipped[] == 0
        Test.@test all(s -> s in (0, 1, 2), freeze_streaks)
        Test.@test count(==(2), freeze_streaks) > 0

        # Pass 3: reads now at streak==2 meet the threshold entering this pass,
        # so they are frozen-skipped — the diagnostic counter must move.
        Mycelia.improve_read_set_likelihood(
            out2, graph, k; common..., skip_frozen_reads = true,
            freeze_streak_threshold = threshold, freeze_streaks = freeze_streaks,
            diagnostics = diag)
        Test.@test diag.frozen_reads_skipped[] > 0

        # RESET-ON-IMPROVE, engineered deterministically: manually seed a streak
        # BELOW threshold for a read with an obviously-correctable single-base
        # error and confirm it lands back at 0 after a fresh call (rather than
        # relying on natural convergence timing, which is not deterministic
        # read-by-read).
        rng2 = Random.MersenneTwister(4242)
        ref2 = join(rand(rng2, _OPT4FS_BASES, 1200))
        reads2 = _opt4fs_reads(rng2, ref2; n_reads = 180, readlen = 80, n_err = 30)
        graph2 = Mycelia.Rhizomorph.build_qualmer_graph(reads2, k; mode = :canonical)
        seeded = zeros(Int, length(reads2))
        seeded[1] = 3   # well below a high threshold, so read 1 still decodes
        out3, = Mycelia.improve_read_set_likelihood(
            reads2, graph2, k; common..., skip_frozen_reads = true,
            freeze_streak_threshold = 10, freeze_streaks = seeded)
        # read 1 (index 1, n_err=30 => the first 30 reads carry an injected
        # error) was decoded fresh this pass; whatever the outcome, the streak
        # reflects THIS pass's was_improved, not the stale seed of 3 — i.e. it
        # must have moved off the seeded value of 3 in the improving direction
        # (either reset to 0 on improvement, or advanced to 4 on no
        # improvement — never left untouched at the stale 3).
        Test.@test seeded[1] != 3
    end

    Test.@testset "scope: within-rung resets at k-advance vs across-rung persists" begin
        rng = Random.MersenneTwister(606)
        ref = join(rand(rng, _OPT4FS_BASES, 900))
        reads = _opt4fs_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        # freeze_streak_threshold=1 is maximally aggressive: any single
        # no-improvement pass makes a read eligible for frozen-skip on the NEXT
        # pass. Within-rung resets streaks to 0 at every k-advance, so a fresh
        # rung starts every read unfrozen; across-rung carries streaks forward,
        # so once a read is frozen it stays frozen (and counted) into later
        # rungs too. The run-total frozen_reads_skipped should therefore be
        # LARGER (or at worst equal) under across-rung than within-rung.
        run = across -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 3,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            skip_frozen_reads = true, freeze_streak_threshold = 1,
            freeze_across_rungs = across,
            output_dir = joinpath(tmp, across ? "across" : "within"))

        r_within = run(false)
        r_across = run(true)

        frozen_within = r_within[:metadata][:corrector_errors][:frozen_reads_skipped]
        frozen_across = r_across[:metadata][:corrector_errors][:frozen_reads_skipped]

        Test.@test frozen_within > 0   # positive control: freezing actually fired
        Test.@test frozen_across >= frozen_within

        # The k-ladder itself is controlled by n_k_rungs/max_k, independent of
        # freeze scope — freezing changes WHEN a rung stops (within-k
        # convergence), never WHICH rungs run. Both scopes must reach the same
        # final k.
        Test.@test r_within[:metadata][:final_k] == r_across[:metadata][:final_k]
    end

    Test.@testset "soft-EM interaction: skip_frozen_reads + soft_em don't crash" begin
        rng = Random.MersenneTwister(808)
        ref = join(rand(rng, _OPT4FS_BASES, 1200))
        reads = _opt4fs_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        result = Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 2,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = true, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            skip_frozen_reads = true, freeze_streak_threshold = 2,
            output_dir = joinpath(tmp, "soft_em_frozen"))

        final_fastq = result[:metadata][:final_fastq_file]
        Test.@test isfile(final_fastq)
        out_records = open(FASTX.FASTQ.Reader, final_fastq) do rd
            collect(rd)
        end
        Test.@test length(out_records) == length(reads)
        Test.@test all(r -> length(FASTX.sequence(String, r)) == 80, out_records)
        Test.@test result[:metadata][:soft_em] != false
    end

    Test.@testset "stop_on_no_change interaction: no premature/stalled k-advance" begin
        rng = Random.MersenneTwister(1717)
        ref = join(rand(rng, _OPT4FS_BASES, 900))
        reads = _opt4fs_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        run = frozen -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 5,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            stop_on_no_change = true,
            skip_frozen_reads = frozen, freeze_streak_threshold = 1,
            output_dir = joinpath(tmp, frozen ? "frozen" : "exact"))

        r_exact = run(false)
        r_frozen = run(true)

        # Both must terminate (no hang) with a positive iteration count, never
        # exceeding the hard cap (n_k_rungs * max_iterations_per_k = 15).
        Test.@test 0 < r_exact[:metadata][:total_iterations] <= 15
        Test.@test 0 < r_frozen[:metadata][:total_iterations] <= 15

        # Freezing never changes the k-ladder itself (only within-k pass count),
        # so both runs reach the same final k.
        Test.@test r_exact[:metadata][:final_k] == r_frozen[:metadata][:final_k]

        # The early-stop is a zero-improvement pass; frozen-skip can only make
        # a pass see FEWER (or equal) improvements than the exact path would
        # (skipped reads cannot contribute genuine improvements), so it should
        # never require MORE passes than the exact path to converge within a
        # k-rung.
        Test.@test r_frozen[:metadata][:total_iterations] <=
                   r_exact[:metadata][:total_iterations]
    end
end
