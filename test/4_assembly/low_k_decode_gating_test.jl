# LOW-K DECODE GATING (td-9h5r)
# =============================
#
# The per-read graph-Viterbi decode is the dominant runtime term of the :scalable
# corrector (57-78% of every iteration in the #370 profile). At LOW k the graph is
# dense — nearly every k-mer sits on a bubble/repeat vertex — so the hard-window
# gate cannot discriminate and flags essentially every read for a full whole-read
# decode, correction Stage 0's cheap linear pass already largely did. This gates
# that wasted low-k decode volume WITHOUT touching the decode algorithm itself, via
# two composable levers, BOTH gated behind the :scalable tier (hard_window=true) so
# :exhaustive stays byte-identical:
#
#   (1) EXPLICIT floor  — `min_decode_k`: skip the decode for every rung k < floor
#       (deterministic; "start the graph-decode ladder at a higher k").
#   (2) ADAPTIVE gate   — `decode_gate_density`: skip the decode for a pass whose
#       post-Stage-0 hard-window decode fraction would be >= the threshold (the gate
#       is non-discriminating). A genuinely selective low-k decode stays ON.
#
# WHAT IS ASSERTED
#   * `_read_set` explicit floor: decode_enabled=false ⇒ the whole pass is gated
#     (skip_fraction == 1.0, decode_gated == true), but Stage 0 cheap corrections
#     still happen (error correction is not disabled, only the Viterbi is deferred).
#   * `_read_set` adaptive gate: a NON-discriminating hard set (covers ~all reads)
#     with a density threshold below 1.0 gates the pass; a SELECTIVE hard set (few
#     reads) does NOT gate even at the same threshold (the load-bearing decode is
#     preserved).
#   * Integration: mycelia_iterative_assemble with hard_window=true + a high
#     min_decode_k gates the low rungs and surfaces `decode_gated_rungs`.
#   * :exhaustive byte-identical: with hard_window=false the gating knobs are
#     IGNORED — an explicit min_decode_k / decode_gate_density produces identical
#     corrected output to the ungated run.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/low_k_decode_gating_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _LK_BASES = ['A', 'C', 'G', 'T']

function _lk_reads(rng, ref; n_reads = 120, readlen = 80)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = ref[s:(s + readlen - 1)]
        push!(records, FASTX.FASTQ.Record("r$i", seq, String(fill('I', readlen))))
    end
    return records
end

Test.@testset "low-k decode gating (td-9h5r)" begin
    R = Mycelia.Rhizomorph

    Test.@testset "explicit floor: decode_enabled=false gates the whole pass" begin
        rng = Random.MersenneTwister(101)
        ref = join(rand(rng, _LK_BASES, 1000))
        # Inject single-substitution errors so Stage 0 has simple errors to fix.
        reads = _lk_reads(rng, ref; n_reads = 150)
        for i in 1:20
            s = rand(rng, 1:921)
            seq = collect(ref[s:(s + 79)])
            seq[40] = rand(rng, filter(!=(seq[40]), _LK_BASES))
            push!(reads, FASTX.FASTQ.Record("e$i", String(seq), String(fill('I', 80))))
        end
        k = 13
        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # decode_enabled=false ⇒ every read skips the decode; Stage 0 still runs.
        _out, imp, skip, cheap, gated = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = false)
        Test.@test gated == true
        Test.@test skip == 1.0                 # nothing decoded
        Test.@test cheap >= 0                  # Stage 0 still allowed to run
        Test.@test imp == cheap                # all improvements are cheap-correction

        # decode_enabled=true (no adaptive threshold) ⇒ normal, decode runs.
        _o2, _i2, skip2, _c2, gated2 = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true)
        Test.@test gated2 == false
        Test.@test skip2 < 1.0                 # at least some reads reached the decode
    end

    Test.@testset "adaptive gate: fires on a non-discriminating hard set only" begin
        rng = Random.MersenneTwister(202)
        ref = join(rand(rng, _LK_BASES, 800))
        reads = _lk_reads(rng, ref; n_reads = 120)
        k = 13
        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        all_labels = collect(R.MetaGraphsNext.labels(graph))

        # NON-discriminating hard set: EVERY vertex is "hard" ⇒ every read with a
        # k-mer would decode ⇒ natural decode fraction ~1.0. A density threshold
        # below that must gate the pass.
        hard_all = Set(all_labels)
        _o, _i, skip_all, _c, gated_all = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = false,
            cheap_correct = false, hard_vertices = hard_all,
            decode_gate_density = 0.90)
        Test.@test gated_all == true
        Test.@test skip_all == 1.0

        # Same non-discriminating hard set but NO threshold ⇒ NOT gated (decodes).
        _o2, _i2, skip_none, _c2, gated_none = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = false,
            cheap_correct = false, hard_vertices = hard_all,
            decode_gate_density = nothing)
        Test.@test gated_none == false
        Test.@test skip_none < 1.0

        # SELECTIVE hard set: only a couple of vertices are hard ⇒ few reads decode
        # ⇒ natural decode fraction is LOW ⇒ the adaptive gate must NOT fire even at
        # the same 0.90 threshold (a load-bearing selective decode is preserved).
        hard_few = Set(all_labels[1:min(2, length(all_labels))])
        _o3, _i3, skip_few, _c3, gated_few = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = false,
            cheap_correct = false, hard_vertices = hard_few,
            decode_gate_density = 0.90)
        Test.@test gated_few == false
        Test.@test skip_few > 0.0              # most reads skipped (few are hard)
    end

    Test.@testset "adaptive gate requires an active hard-window gate" begin
        # With no hard set (hard_vertices=nothing) there is no discrimination signal,
        # so the adaptive gate is inert regardless of the threshold — the decode runs
        # as before (this keeps the gate from firing on a skip_solid-only pass).
        rng = Random.MersenneTwister(303)
        ref = join(rand(rng, _LK_BASES, 800))
        reads = _lk_reads(rng, ref; n_reads = 100)
        k = 13
        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        _o, _i, _skip, _c, gated = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = false,
            cheap_correct = false, hard_vertices = nothing,
            decode_gate_density = 0.10)
        Test.@test gated == false
    end

    Test.@testset "integration: explicit floor gates low rungs (:scalable)" begin
        rng = Random.MersenneTwister(404)
        ref = join(rand(rng, _LK_BASES, 1500))
        reads = _lk_reads(rng, ref; n_reads = 200, readlen = 100)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        # Force the graph-decode ladder to start at the top rung: a high min_decode_k
        # (>= the largest rung) gates every lower rung. hard_window=true ⇒ :scalable.
        res = Mycelia.mycelia_iterative_assemble(fq;
            max_k = 21, skip_solid = true, graph_mode = :doublestrand,
            n_k_rungs = 3, max_iterations_per_k = 2, hard_window = true,
            soft_em = false, cheap_correct = true, beam_width = nothing,
            min_decode_k = 1000, verbose = false, enable_checkpointing = false,
            output_dir = joinpath(tmp, "out"))
        meta = res[:metadata]
        Test.@test meta[:min_decode_k] == 1000
        # Every processed rung is below 1000 ⇒ all are gated.
        Test.@test !isempty(meta[:decode_gated_rungs])
        Test.@test Set(meta[:decode_gated_rungs]) == Set(res[:k_progression])
        # Every pass skipped the whole decode ⇒ skip fraction 1.0 throughout.
        Test.@test all(==(1.0), meta[:skip_fraction_per_pass])
    end

    Test.@testset ":exhaustive is byte-identical (gating knobs ignored)" begin
        # hard_window=false ⇒ the :exhaustive tier. The low-k gating knobs must be
        # NO-OPS: a run with an aggressive min_decode_k + density threshold must
        # produce byte-identical corrected reads to the ungated run.
        rng = Random.MersenneTwister(505)
        ref = join(rand(rng, _LK_BASES, 600))
        reads = _lk_reads(rng, ref; n_reads = 60, readlen = 80)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        run_corrector = (; min_decode_k, decode_gate_density, out) ->
            Mycelia.mycelia_iterative_assemble(fq;
                max_k = 17, skip_solid = false, graph_mode = :canonical,
                n_k_rungs = 3, max_iterations_per_k = 2, hard_window = false,
                soft_em = false, cheap_correct = false, beam_width = typemax(Int),
                min_decode_k = min_decode_k, decode_gate_density = decode_gate_density,
                verbose = false, enable_checkpointing = false, output_dir = out)

        r_plain = run_corrector(; min_decode_k = nothing, decode_gate_density = nothing,
            out = joinpath(tmp, "plain"))
        r_gated = run_corrector(; min_decode_k = 1000, decode_gate_density = 0.10,
            out = joinpath(tmp, "gated"))

        seqs(res) = [FASTX.sequence(String, r) for r in
                     open(FASTX.FASTQ.Reader, res[:metadata][:final_fastq_file]) do rd
                         collect(rd)
                     end]
        Test.@test seqs(r_plain) == seqs(r_gated)
        # And the exhaustive tier records NO gating (byte-identical passthrough).
        Test.@test r_gated[:metadata][:min_decode_k] === nothing
        Test.@test isempty(r_gated[:metadata][:decode_gated_rungs])
    end
end
