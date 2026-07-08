# Stage 3c (td-nn6l): per-hard-region WINDOWED decode. Instead of decoding a hard
# read end-to-end, `improve_read_likelihood_windowed` decodes ONLY the boundary-
# constrained hard sub-window(s) around the read's hard vertices (each <= 500 bp)
# and splices the corrected window(s) back. This test asserts:
#   (a) BOUNDING — the windowed decode touches only a tiny fraction of a long
#       read's bases (O(window) not O(read_length)), the #375 super-linear term;
#   (b) SPEED — per-hard-read decode time drops vs whole-read decode on long
#       (nanopore-regime ~3 kb) reads;
#   (c) CORRECTNESS — the windowed correction matches the whole-read decode on the
#       hard region (any divergence is surfaced and asserted to be zero);
#   (d) LENGTH PRESERVATION — the spliced read has the same length as the input
#       (the window ML path is length-preserving; divergent windows are dropped);
#   (e) PASSTHROUGH — a read that touches no hard vertex is returned unchanged.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/windowed_decode_test.jl")'

import Test
import Mycelia
import FASTX
import Random
import Logging

const _WD_BASES = ['A', 'C', 'G', 'T']

# Build a long-read fixture (nanopore regime) over a reference, with a fixed-locus
# substitution injected via several reads so the graph carries a hard region
# (bubble / weak-k-mer neighborhood) that the windowed decode targets.
function _wd_long_read_fixture(; seed = 42, reflen = 5000, readlen = 3000,
        n_clean = 60, n_err = 8, err_locus = 2500)
    rng = Random.MersenneTwister(seed)
    ref = join(rand(rng, _WD_BASES, reflen))
    reads = FASTX.FASTQ.Record[]
    for i in 1:n_clean
        s = rand(rng, 1:(reflen - readlen + 1))
        push!(reads, FASTX.FASTQ.Record("r$i", ref[s:(s + readlen - 1)],
            String(fill('I', readlen))))
    end
    err_reads = FASTX.FASTQ.Record[]
    for i in 1:n_err
        s = err_locus - div(readlen, 2)              # read spans the error locus
        off = div(readlen, 2) + 1
        seq = collect(ref[s:(s + readlen - 1)])
        seq[off] = first(filter(!=(seq[off]), _WD_BASES))   # deterministic substitution
        push!(err_reads, FASTX.FASTQ.Record("e$i", String(seq), String(fill('I', readlen))))
    end
    return (ref = ref, reads = vcat(reads, err_reads), err_reads = err_reads)
end

Test.@testset "windowed decode (td-nn6l Stage 3c)" begin
    R = Mycelia.Rhizomorph
    k = 13
    mode = :doublestrand   # the mode the :scalable tier actually runs (td-nt69)

    fx = _wd_long_read_fixture()
    graph = R.build_qualmer_graph(fx.reads, k; mode = mode)
    hard = Mycelia._hard_vertex_set(graph, k)
    Test.@test !isempty(hard)

    # Every injected long error read touches the hard region ⇒ is a decode target.
    for r in fx.err_reads
        Test.@test Mycelia.should_decode_read(r, k, hard) == true
    end

    Test.@testset "bounding: windows cover a tiny fraction of a long read" begin
        total_read_bases = 0
        total_window_bases = 0
        for r in fx.err_reads
            rlen = length(FASTX.sequence(r))
            wins = Mycelia._hard_window_ranges(r, k, hard; pad = k, max_window = 500)
            Test.@test !isempty(wins)                       # hard read ⇒ >= 1 window
            for w in wins
                Test.@test first(w) >= 1 && last(w) <= rlen  # in-bounds
                Test.@test length(w) <= 500                  # capped
            end
            total_read_bases += rlen
            total_window_bases += sum(length.(wins); init = 0)
        end
        # The windowed decode only touches the hard-window bases — far less than the
        # whole-read volume (the point of windowing: O(window) not O(read_length)).
        Test.@test total_window_bases < total_read_bases
        Test.@test total_window_bases < 0.25 * total_read_bases
        @info "windowed-decode bounding" total_read_bases total_window_bases fraction = round(
            total_window_bases / total_read_bases, digits = 4)
    end

    Test.@testset "correctness: windowed matches whole-read on the hard region" begin
        # For every hard read, the windowed decode must agree with the whole-read
        # decode ON EACH HARD WINDOW (the region it actually decodes), be
        # length-preserving, and report zero divergent (dropped) windows. Both use
        # exact ML (beam_width=typemax) to isolate the windowing effect.
        n_whole_improved = 0
        n_win_improved = 0
        Logging.with_logger(Logging.NullLogger()) do
            for r in fx.err_reads
                orig = FASTX.sequence(String, r)
                wins = Mycelia._hard_window_ranges(r, k, hard; pad = k, max_window = 500)

                rec_w, imp_w = Mycelia.improve_read_likelihood(
                    r, graph, k; graph_mode = mode, beam_width = typemax(Int))
                rec_win, imp_win, ndec, ndiv = Mycelia.improve_read_likelihood_windowed_detail(
                    r, graph, k, hard; graph_mode = mode, beam_width = typemax(Int))

                sw = FASTX.sequence(String, rec_w)
                swin = FASTX.sequence(String, rec_win)

                # Length preserved end-to-end; no window dropped for a length mismatch.
                Test.@test length(swin) == length(orig)
                Test.@test ndec >= 1
                Test.@test ndiv == 0

                # On each hard window the windowed correction equals the whole-read
                # correction (match-or-approximate contract; here it MATCHES).
                for w in wins
                    lo, hi = first(w), last(w)
                    Test.@test swin[lo:hi] == sw[lo:hi]
                end

                # Outside the hard windows the windowed decode leaves the read
                # untouched (only hard regions are rewritten).
                covered = falses(length(orig))
                for w in wins
                    covered[first(w):last(w)] .= true
                end
                for i in 1:length(orig)
                    covered[i] || Test.@test swin[i] == orig[i]
                end

                n_whole_improved += imp_w ? 1 : 0
                n_win_improved += imp_win ? 1 : 0
            end
        end
        # Windowing must not silently drop corrections the whole-read decode makes:
        # whenever whole-read improves reads, windowed improves at least as many
        # (it targets the same hard regions, exactly).
        Test.@test n_win_improved >= n_whole_improved
        @info "windowed-decode correctness" n_whole_improved n_win_improved
    end

    Test.@testset "speed: windowed decode is faster per hard read (long reads)" begin
        probe = fx.err_reads[1]
        Logging.with_logger(Logging.NullLogger()) do
            # Warm up both paths (exclude compilation from the timing).
            Mycelia.improve_read_likelihood(
                probe, graph, k; graph_mode = mode, beam_width = typemax(Int))
            Mycelia.improve_read_likelihood_windowed(
                probe, graph, k, hard; graph_mode = mode, beam_width = typemax(Int))

            t_whole = @elapsed for r in fx.err_reads
                Mycelia.improve_read_likelihood(
                    r, graph, k; graph_mode = mode, beam_width = typemax(Int))
            end
            t_win = @elapsed for r in fx.err_reads
                Mycelia.improve_read_likelihood_windowed(
                    r, graph, k, hard; graph_mode = mode, beam_width = typemax(Int))
            end
            n = length(fx.err_reads)
            @info "windowed-decode per-hard-read time" readlen = length(FASTX.sequence(probe)) whole_per_read_ms = round(
                t_whole / n * 1000, digits = 1) windowed_per_read_ms = round(
                t_win / n * 1000, digits = 1) speedup = round(t_whole / t_win, digits = 2)
            # Bounding a ~3 kb read's decode to its <=500 bp hard windows is a large,
            # robust margin (~4-5x here); assert strict improvement.
            Test.@test t_win < t_whole
        end
    end

    Test.@testset "passthrough: a non-hard read is returned unchanged" begin
        # A clean read pulled from a well-covered region far from the hard locus
        # touches no hard vertex ⇒ no windows ⇒ returned unchanged, not improved.
        clean_probe = FASTX.FASTQ.Record(
            "clean_probe", fx.ref[100:(100 + 199)], String(fill('I', 200)))
        Test.@test Mycelia.should_decode_read(clean_probe, k, hard) == false
        rec, improved = Mycelia.improve_read_likelihood_windowed(
            clean_probe, graph, k, hard; graph_mode = mode, beam_width = typemax(Int))
        Test.@test improved == false
        Test.@test FASTX.sequence(String, rec) == fx.ref[100:(100 + 199)]
    end

    Test.@testset "end-to-end :scalable assemble runs windowed decode" begin
        # The provenance flag is honest: :scalable now decodes hard reads
        # window-by-window (windowed_decode=true), and still produces a real
        # assembly with a non-trivial skip fraction.
        res = R.assemble_genome(fx.reads; k = k, corrector = :iterative, strategy = :scalable)
        Test.@test res isa R.AssemblyResult
        Test.@test !isempty(res.contigs)
        Test.@test res.assembly_stats["windowed_decode"] == true
        Test.@test res.assembly_stats["hard_window"] == true
        skip = res.assembly_stats["skip_fraction"]
        Test.@test 0.0 < skip <= 1.0
    end
end
