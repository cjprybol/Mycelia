# OPT4 FROZEN-READ SKIP — DISABLED-PATH IDENTITY LOCK (td-jbjd, pass 1)
# =======================================================================
#
# opt4 is the corrector campaign's first APPROXIMATE optimization: skipping
# re-decode of "frozen" (converged, N-consecutive-pass-stable) reads trades a
# bounded amount of accuracy for skipped work, so it is opt-in and ships
# default OFF (`skip_frozen_reads=false`). Unlike opt5/opt1/opt2, the shipped
# deliverable for opt4 is an accuracy sign-off, not a byte-identity lock — but
# the DISABLED path (the default for every existing caller) MUST remain a true
# no-op: identical output whether the kwarg is passed explicitly as `false` or
# omitted entirely. That disabled-path guarantee is what this test locks.
#
# See docs/design/2026-07-26-opt4-frozen-read-skip-design.md.
#
# WHAT IS ASSERTED
#   * Unit (improve_read_set_likelihood): a run with `skip_frozen_reads=false`
#     passed explicitly is byte-identical to a run where the kwarg (and its
#     siblings freeze_streak_threshold/freeze_across_rungs) are omitted.
#   * Integration (mycelia_iterative_assemble), multi-pass/multi-k so the
#     freeze-streak state-threading machinery actually runs its bookkeeping
#     paths (just gated off by the flag): omitted-kwarg vs explicit-false
#     produce byte-identical final corrected FASTQ.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_frozen_read_skip_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _OPT4_BASES = ['A', 'C', 'G', 'T']

# Reads sampled from ref, with a single substitution injected into the first
# n_err of them, so Stage 0 + the per-read decode both have work across
# multiple batches/passes (exercising the new freeze-streak bookkeeping sites
# even though skip_frozen_reads stays off).
function _opt4_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _OPT4_BASES))
        end
        push!(records,
            FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

# Full-record fingerprint (identifier, sequence, quality).
function _opt4_recs(records)
    [(String(FASTX.identifier(r)), FASTX.sequence(String, r),
         String(FASTX.quality(r))) for r in records]
end

Test.@testset "opt4 skip_frozen_reads disabled-path identity lock (td-jbjd)" begin
    Test.@testset "unit: improve_read_set_likelihood, explicit-false == omitted" begin
        rng = Random.MersenneTwister(4242)
        ref = join(rand(rng, _OPT4_BASES, 1200))
        reads = _opt4_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
        k = 13
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # batch_size < n_reads => multiple batches, exercising both the serial
        # and (with enable_parallel) the parallel collect-results freeze-streak
        # bookkeeping sites — all of which must be no-ops here.
        common = (; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true,
            batch_size = 50)

        out_omitted, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; common...)
        out_explicit, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; common...,
            skip_frozen_reads = false, freeze_streak_threshold = 2,
            freeze_across_rungs = false)

        Test.@test _opt4_recs(out_omitted) == _opt4_recs(out_explicit)
    end

    Test.@testset "integration: mycelia_iterative_assemble, explicit-false == omitted" begin
        rng = Random.MersenneTwister(909)
        ref = join(rand(rng, _OPT4_BASES, 900))
        reads = _opt4_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        # Multi-k, multi-iteration so freeze-streak allocation/reset (k-advance)
        # and per-read update sites all execute at least once — still a no-op
        # since skip_frozen_reads is never true.
        run_common = out -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 3,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false, output_dir = out)

        r_omitted = run_common(joinpath(tmp, "omitted"))
        r_explicit = Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 3,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            output_dir = joinpath(tmp, "explicit"),
            skip_frozen_reads = false, freeze_streak_threshold = 2,
            freeze_across_rungs = false)

        final_recs = res -> _opt4_recs(
            open(FASTX.FASTQ.Reader, res[:metadata][:final_fastq_file]) do rd
            collect(rd)
        end)
        Test.@test final_recs(r_omitted) == final_recs(r_explicit)
    end
end
