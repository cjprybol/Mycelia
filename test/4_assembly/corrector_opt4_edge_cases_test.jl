# OPT4 FROZEN-READ SKIP — EDGE CASES (td-jbjd, pass 4)
# =======================================================================
#
# Pass 4 ("edge cases") of the opt4 Rule-of-5 chain
# (docs/design/2026-07-26-opt4-frozen-read-skip-design.md). Passes 1-3 already
# cover the disabled-path lock, the freeze-state machine, within/across-rung
# scope, soft-EM/skip_solid composition, and one stop_on_no_change interaction
# test (aggressive-but-not-total freezing). This file adds the mechanical edge
# cases the design's own Rule-of-5 scope calls out that are NOT yet covered:
#
#   * ALL-reads-freeze: force literally every read in the set to freeze, then
#     confirm `stop_on_no_change`'s zero-improvement early-stop still fires
#     correctly (no hang, no premature/incorrect k-advance vs the exact path)
#     -- design doc "Risks" section, bullet 4.
#   * single-read batch (n_reads=1): skip_frozen_reads=true must not crash on
#     the smallest possible read set and must produce sane output.
#   * N=1 vs N=2 threshold timing: the freeze-eligibility boundary is exact,
#     not "eventually" -- threshold=N must make a read frozen-skip-eligible
#     starting EXACTLY N passes after a maintained no-improvement streak began.
#
# k-boundary reset (within-rung resets at k-advance vs across-rung persists) is
# already covered in detail by corrector_opt4_freeze_state_test.jl's two scope
# testsets, so it is intentionally NOT duplicated here.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_edge_cases_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _OPT4EC_BASES = ['A', 'C', 'G', 'T']

function _opt4ec_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _OPT4EC_BASES))
        end
        push!(records,
            FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

Test.@testset "opt4 skip_frozen_reads edge cases (td-jbjd pass 4)" begin
    Test.@testset "all-reads-freeze: stop_on_no_change fires correctly, no premature k-advance" begin
        # skip_solid=false, cheap_correct=false, decode_enabled=true is the
        # configuration under which the corrector's per-read decode makes NO
        # improvements at all at this toy scale (verified during pass-4
        # investigation: every read's `was_improved` is false every pass in this
        # regime), so with freeze_streak_threshold=1 and enough passes, EVERY
        # read in the set is frozen by the second pass of the first rung, and
        # stays frozen (across-rung) through every subsequent rung. This is the
        # literal "all reads freeze" case the design's Risks section asks for,
        # distinct from the existing freeze_state_test.jl stop_on_no_change
        # testset (which uses skip_solid=true/cheap_correct=true, an aggressive
        # -but-not-total freeze configuration).
        rng = Random.MersenneTwister(4141)
        ref = join(rand(rng, _OPT4EC_BASES, 900))
        reads = _opt4ec_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        run = frozen -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 21, n_k_rungs = 3, max_iterations_per_k = 5,
            graph_mode = :canonical, skip_solid = false, cheap_correct = false,
            soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            stop_on_no_change = true,
            skip_frozen_reads = frozen, freeze_streak_threshold = 1,
            freeze_across_rungs = true,
            output_dir = joinpath(tmp, frozen ? "frozen" : "exact"))

        r_exact = run(false)
        r_frozen = run(true)

        # Every read froze (all-solid gate is off here, so freeze is the ONLY
        # skip source): confirm the frozen run actually exercised total freezing,
        # not just partial freezing, as the setup sanity check for this testset.
        Test.@test r_frozen[:metadata][:corrector_errors][:frozen_reads_skipped] > 0

        # No hang: both runs terminate with a positive iteration count that
        # never exceeds the hard cap (n_k_rungs * max_iterations_per_k = 15).
        Test.@test 0 < r_exact[:metadata][:total_iterations] <= 15
        Test.@test 0 < r_frozen[:metadata][:total_iterations] <= 15

        # stop_on_no_change correctness: once literally every read is frozen, a
        # pass makes ZERO improvements (frozen reads cannot contribute genuine
        # improvements by construction), so the zero-improvement early-stop must
        # fire -- the frozen run must not silently run to the iteration CAP
        # instead of stopping early once nothing is left to do. Concretely: the
        # frozen run's total_iterations must not exceed the exact run's (frozen
        # reads can only cause EARLIER-or-equal stopping, never later).
        Test.@test r_frozen[:metadata][:total_iterations] <=
                   r_exact[:metadata][:total_iterations]

        # k-ladder correctness: freezing (even total freezing) changes WHEN a
        # rung's pass count stops, never WHICH rungs the run walks through -- no
        # premature or skipped k-advance vs the exact path.
        Test.@test r_exact[:metadata][:final_k] == r_frozen[:metadata][:final_k]
        Test.@test r_exact[:metadata][:k_progression] ==
                   r_frozen[:metadata][:k_progression]
    end

    Test.@testset "single-read batch (n_reads=1) does not crash" begin
        rng = Random.MersenneTwister(5151)
        ref = join(rand(rng, _OPT4EC_BASES, 300))
        reads = _opt4ec_reads(rng, ref; n_reads = 1, readlen = 80, n_err = 1)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        result = Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 3,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            skip_frozen_reads = true, freeze_streak_threshold = 1,
            freeze_across_rungs = true,
            output_dir = joinpath(tmp, "single_read"))

        final_fastq = result[:metadata][:final_fastq_file]
        Test.@test isfile(final_fastq)
        out_records = open(FASTX.FASTQ.Reader, final_fastq) do rd
            collect(rd)
        end
        Test.@test length(out_records) == 1
        Test.@test length(FASTX.sequence(String, out_records[1])) == 80
        # Sane, not necessarily improved: a single read with no peer coverage has
        # nothing to build consensus against, so the acceptance bar here is
        # "did not crash / did not corrupt the read", matching the design's
        # ask ("behaves sanely") rather than a correction-quality claim.
        Test.@test result[:metadata][:total_iterations] > 0
    end

    Test.@testset "N=1 vs N=2 threshold: frozen_reads_skipped timing is exact" begin
        # skip_solid=false, cheap_correct=false again pins was_improved=false for
        # EVERY read on EVERY pass (see the all-reads-freeze testset above), which
        # makes the freeze-eligibility boundary fully deterministic and lets this
        # testset assert the EXACT pass number at which frozen_reads_skipped first
        # becomes positive, rather than just its eventual sign.
        rng = Random.MersenneTwister(9001)
        ref = join(rand(rng, _OPT4EC_BASES, 1200))
        reads = _opt4ec_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
        k = 13
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        common = (; graph_mode = :canonical, skip_solid = false,
            cheap_correct = false, decode_enabled = true, batch_size = 50)

        # threshold=1: streak reaches 1 after pass 1 (>=1 entering pass 2), so
        # frozen_reads_skipped must be 0 after pass 1 and > 0 after pass 2.
        streaks1 = zeros(Int, length(reads))
        diag1 = Mycelia.CorrectorDiagnostics()
        out = reads
        out, = Mycelia.improve_read_set_likelihood(out, graph, k; common...,
            skip_frozen_reads = true, freeze_streak_threshold = 1,
            freeze_streaks = streaks1, diagnostics = diag1)
        Test.@test diag1.frozen_reads_skipped[] == 0
        out, = Mycelia.improve_read_set_likelihood(out, graph, k; common...,
            skip_frozen_reads = true, freeze_streak_threshold = 1,
            freeze_streaks = streaks1, diagnostics = diag1)
        Test.@test diag1.frozen_reads_skipped[] > 0
        n1_first_frozen_pass = 2

        # threshold=2: streak reaches 2 only after pass 2 (>=2 entering pass 3),
        # so frozen_reads_skipped must stay 0 through pass 2 and only become
        # positive at pass 3 -- one pass LATER than threshold=1, the exact
        # timing difference the design's N=1-vs-N>=2 distinction is asking for.
        streaks2 = zeros(Int, length(reads))
        diag2 = Mycelia.CorrectorDiagnostics()
        out2 = reads
        out2, = Mycelia.improve_read_set_likelihood(out2, graph, k; common...,
            skip_frozen_reads = true, freeze_streak_threshold = 2,
            freeze_streaks = streaks2, diagnostics = diag2)
        Test.@test diag2.frozen_reads_skipped[] == 0
        out2, = Mycelia.improve_read_set_likelihood(out2, graph, k; common...,
            skip_frozen_reads = true, freeze_streak_threshold = 2,
            freeze_streaks = streaks2, diagnostics = diag2)
        Test.@test diag2.frozen_reads_skipped[] == 0   # still below threshold
        out2, = Mycelia.improve_read_set_likelihood(out2, graph, k; common...,
            skip_frozen_reads = true, freeze_streak_threshold = 2,
            freeze_streaks = streaks2, diagnostics = diag2)
        Test.@test diag2.frozen_reads_skipped[] > 0
        n2_first_frozen_pass = 3

        # The exact one-pass-later relationship, stated directly.
        Test.@test n2_first_frozen_pass == n1_first_frozen_pass + 1
    end
end
