# PARALLEL SOFT-EM BYTE-IDENTITY (opt1)
# Corrected reads AND soft-EM weights must be identical between parallel and
# serial decode. CI runs single-threaded, so this test relaunches itself under
# `julia -t 4` to exercise genuine concurrency (otherwise @threads runs serially
# and proves nothing). RED signal before the fix: the guard demotes
# enable_parallel=true+soft-EM to serial and logs "ignoring enable_parallel".
import Test
import Mycelia
import FASTX
import Random
import Logging

if Threads.nthreads() == 1 && get(ENV, "MYCELIA_PSI_CHILD", "0") != "1"
    # Self-hoist to real threads.
    proj = Base.active_project()
    thisfile = @__FILE__
    cmd = `$(Base.julia_cmd()) --project=$(proj) --threads=4 -e "include(\"$(thisfile)\")"`
    Test.@testset "parallel soft-EM byte-identity (hoisted to -t4)" begin
        Test.@test success(pipeline(setenv(cmd, "MYCELIA_PSI_CHILD" => "1"; dir = pwd());
            stdout = stdout, stderr = stderr))
    end
else
    const _PSI_BASES = ['A', 'C', 'G', 'T']

    function _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
        reflen = length(ref)
        records = FASTX.FASTQ.Record[]
        for i in 1:n_reads
            s = rand(rng, 1:(reflen - readlen + 1))
            seq = collect(ref[s:(s + readlen - 1)])
            if i <= n_err
                p = rand(rng, 1:readlen)
                seq[p] = rand(rng, filter(!=(seq[p]), _PSI_BASES))
            end
            push!(records,
                FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
        end
        return records
    end

    _psi_seqs(records) = [FASTX.sequence(String, r) for r in records]
    _psi_weight_pairs(acc) = sort!(collect(acc.weights); by = p -> repr(p[1]))

    function _psi_run(reads, graph, k, hard; parallel::Bool)
        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        out, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true,
            batch_size = 50, enable_parallel = parallel, soft_weights = acc)
        return _psi_seqs(out), _psi_weight_pairs(acc)
    end

    Test.@testset "parallel soft-EM byte-identity (opt1)" begin
        Test.@test Threads.nthreads() > 1     # guard: real threads
        rng = Random.MersenneTwister(2024)
        ref = join(rand(rng, _PSI_BASES, 1500))
        reads = _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
        k = 13
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # (a) parallel path taken: no "ignoring enable_parallel" warning
        logger = Test.TestLogger()
        seqs_par, wts_par = Logging.with_logger(logger) do
            _psi_run(reads, graph, k, hard; parallel = true)
        end
        Test.@test !any(occursin("ignoring enable_parallel", r.message)
                        for r in logger.logs)

        # (b) byte-identity vs serial: reads AND soft-EM weights
        seqs_ser, wts_ser = _psi_run(reads, graph, k, hard; parallel = false)
        Test.@test seqs_par == seqs_ser
        Test.@test wts_par == wts_ser
    end

    # ---------------------------------------------------------------------------
    # WINDOWED-DECODE byte-identity — the per-window ordered-capture divergence.
    #
    # Production `:scalable` runs windowed_decode=true + soft_em=true. The windowed
    # decode merges a fresh staged accumulator ONCE PER WINDOW. A per-read
    # accumulator (the pre-fix parallel branch) pre-groups a read's windows as
    # (w1+w2) before folding onto the running total R — giving R+(w1+w2) where the
    # serial path computes (R+w1)+w2. Under float non-associativity these differ
    # whenever a SHARED edge carrying a FRACTIONAL responsibility is hit by >=2
    # windows of one read that prior reads also touched.
    #
    # Construction that forces all three conditions at once:
    #   * TWO tandem arrays of the SAME low-complexity unit, separated by a unique
    #     spacer. Each array raises its own hard window; both windows decode the
    #     SAME repeat sequence, so they accumulate the SAME canonical edges
    #     (within-read cross-window repeat). Every read spans both arrays, so those
    #     edges are also touched by prior reads (R != 0).
    #   * A period-3 low-complexity unit makes the repeat graph a tight, genuinely
    #     ambiguous cycle, so competing paths split responsibility FRACTIONALLY
    #     over those shared edges (integer 1.0 responsibilities would sum
    #     associatively and hide the bug).
    # This case DIVERGES under the old per-read fold (verified: the weight Dict
    # differs) and is byte-identical only under the per-window flat replay.
    # Reads must exceed the windowed long-read threshold (n_obs > 1024).
    function _psi_two_array_ref(rng; flank = 540, unit_len = 44, per_array = 5,
            spacer = 160)
        motif = "ACG"
        unit = (motif^(unit_len ÷ length(motif) + 1))[1:unit_len]
        left = join(rand(rng, _PSI_BASES, flank))
        mid = join(rand(rng, _PSI_BASES, spacer))
        right = join(rand(rng, _PSI_BASES, flank))
        arr = unit^per_array
        ref = left * arr * mid * arr * right
        array1_start = flank + 1
        array2_start = flank + length(arr) + spacer + 1
        return ref, array1_start, array2_start, length(arr)
    end

    function _psi_two_array_reads(rng, ref, array1_start, array2_start, arr_len;
            n_reads = 120, readlen = 1180)
        reflen = length(ref)
        max_start = reflen - readlen + 1
        records = FASTX.FASTQ.Record[]
        for i in 1:n_reads
            s = rand(rng, 1:max(1, max_start))
            seq = collect(ref[s:(s + readlen - 1)])
            # One error inside each array so both raise a hard window.
            for array_start in (array1_start, array2_start)
                off = rand(rng, 5:(arr_len - 5))
                read_pos = array_start + off - s + 1
                if 1 <= read_pos <= readlen
                    seq[read_pos] =
                        rand(rng, filter(!=(seq[read_pos]), _PSI_BASES))
                end
            end
            push!(records,
                FASTX.FASTQ.Record("lr$i", String(seq), String(fill('I', readlen))))
        end
        return records
    end

    function _psi_run_windowed(reads, graph, k, hard; parallel::Bool)
        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        out, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true,
            windowed_decode = true, batch_size = 32,
            enable_parallel = parallel, soft_weights = acc)
        return _psi_seqs(out), _psi_weight_pairs(acc)
    end

    Test.@testset "parallel soft-EM byte-identity — windowed decode (opt1)" begin
        Test.@test Threads.nthreads() > 1     # guard: real threads
        rng = Random.MersenneTwister(11)
        ref, array1_start, array2_start, arr_len = _psi_two_array_ref(rng)
        reads = _psi_two_array_reads(
            rng, ref, array1_start, array2_start, arr_len;
            n_reads = 120, readlen = 1180)
        k = 13

        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # Windowed path must actually engage: >=1 read is "long" AND splits into
        # >=2 hard windows (otherwise the divergence case is not exercised).
        multi_window_reads = count(reads) do r
            Mycelia._windowed_decode_read_is_long(r, k) || return false
            length(Mycelia._hard_window_ranges(
                r, k, hard; graph_mode = :canonical)) >= 2
        end
        Test.@test multi_window_reads >= 1

        seqs_par, wts_par = _psi_run_windowed(reads, graph, k, hard; parallel = true)
        seqs_ser, wts_ser = _psi_run_windowed(reads, graph, k, hard; parallel = false)
        Test.@test seqs_par == seqs_ser
        Test.@test wts_par == wts_ser
    end
end
