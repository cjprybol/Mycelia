# PARALLEL SOFT-EM BYTE-IDENTITY (opt1)
# Corrected reads AND soft-EM weights must be identical between parallel and
# serial decode. CI's default job is single-threaded, so this file relaunches
# itself under `julia -t 4` when `Threads.nthreads() == 1`. The parent re-asserts
# each child result line (`MYCELIA_PSI_RESULTS`) so a failure names the check
# instead of collapsing every child assertion into one `success(pipeline(...))`.
import Test
import Mycelia
import FASTX
import Random

# Parent-side contract: one name per child `_psi_assert`, in emit order.
# A renamed or reordered child check fails here instead of disappearing into
# a single process-success boolean.
const _PSI_EXPECTED_ASSERTIONS = [
    "repr_collision_same_print",
    "repr_collision_sort_keys_differ",
    "nthreads",
    "parallel_batches_taken",
    "serial_batches_zero",
    "weights_nonempty",
    "seqs_equal",
    "weights_equal",
    "single_batch_count",
    "single_batch_weights_nonempty",
    "single_batch_seqs_equal",
    "single_batch_weights_equal",
    "windowed_nthreads",
    "windowed_multi_reads",
    "windowed_divergence_witness",
    "windowed_fractional_shared_weight",
    "windowed_weights_nonempty",
    "windowed_seqs_equal",
    "windowed_weights_equal",
    "gc_weights_nonempty",
    "gc_seqs_equal",
    "gc_weights_equal",
]

if Threads.nthreads() == 1 && get(ENV, "MYCELIA_PSI_CHILD", "0") != "1"
    # Self-hoist to real threads.
    proj = Base.active_project()
    thisfile = @__FILE__
    results_file, results_io = mktemp()
    close(results_io)
    cmd = `$(Base.julia_cmd()) --project=$(proj) --threads=4 -e "include(\"$(thisfile)\")"`
    child_env = copy(ENV)
    child_env["MYCELIA_PSI_CHILD"] = "1"
    child_env["MYCELIA_PSI_RESULTS"] = results_file
    Test.@testset "parallel soft-EM byte-identity (hoisted to -t4)" begin
        proc = run(pipeline(
            setenv(cmd, child_env; dir = pwd());
            stdout = stdout, stderr = stderr); wait = false)
        wait(proc)
        lines = readlines(results_file)
        parsed = Tuple{String, String}[]
        for line in lines
            parts = split(line, '\t', limit = 2)
            if length(parts) == 2
                push!(parsed, (String(parts[2]), String(parts[1])))
            else
                push!(parsed, (line, "MALFORMED"))
            end
        end
        for i in eachindex(_PSI_EXPECTED_ASSERTIONS)
            expected = _PSI_EXPECTED_ASSERTIONS[i]
            if i <= length(parsed)
                name, status = parsed[i]
                Test.@test (name, status) == (expected, "PASS")
            else
                Test.@test (expected, "MISSING") == (expected, "PASS")
            end
        end
        if length(parsed) > length(_PSI_EXPECTED_ASSERTIONS)
            for i in (length(_PSI_EXPECTED_ASSERTIONS) + 1):length(parsed)
                name, status = parsed[i]
                Test.@test (name, status) == ("", "ABSENT")
            end
        end
        Test.@test success(proc)
    end
else
    const _PSI_BASES = ['A', 'C', 'G', 'T']
    _psi_failed = Ref(false)

    function _psi_assert(name::AbstractString, cond::Bool)
        path = get(ENV, "MYCELIA_PSI_RESULTS", "")
        if !isempty(path)
            open(path, "a") do io
                println(io, (cond ? "PASS" : "FAIL"), "\t", name)
                flush(io)
            end
        end
        cond || (_psi_failed[] = true)
        Test.@test cond
        return cond
    end

    # The hoisted child must keep going after a failed `@testset` so later
    # result lines still reach the parent. In-process (the 4-thread CI job)
    # rethrows, so a failure fails this file instead of `exit` aborting the
    # rest of the suite.
    function _psi_guarded(body::Function)
        try
            body()
            return nothing
        catch err
            if err isa Test.TestSetException &&
               !isempty(get(ENV, "MYCELIA_PSI_RESULTS", ""))
                return nothing
            end
            rethrow()
        end
    end

    # `repr` is not injective: two values can print the same and still be
    # unequal. A sort `by = repr(key)` then keeps Dict iteration order for the
    # tie, so equal accumulators can compare unequal. The sort key adds the
    # concrete type and `hash` (equal keys hash equal) and refuses a non-unique
    # key instead of falling through to iteration order.
    function _psi_edge_sort_key(edge)
        return (string(typeof(edge)), repr(edge), hash(edge))
    end

    struct _PsiReprCollideA
        x::Int
    end
    struct _PsiReprCollideB
        x::Int
    end
    Base.show(io::IO, a::_PsiReprCollideA) = print(io, "edge:", a.x)
    Base.show(io::IO, b::_PsiReprCollideB) = print(io, "edge:", b.x)

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

    function _psi_weight_pairs(acc)
        pairs = collect(acc.weights)
        sort!(pairs; by = p -> _psi_edge_sort_key(p[1]))
        sort_keys = [_psi_edge_sort_key(p[1]) for p in pairs]
        allunique(sort_keys) || error(
            "soft-EM weight sort key is not injective " *
            "($(length(sort_keys) - length(unique(sort_keys))) collision(s)); " *
            "refusing to compare weights in Dict iteration order")
        return pairs
    end

    function _psi_run(reads, graph, k, hard; parallel::Bool, diagnostics = nothing,
            batch_size::Int = 50)
        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        out, = Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true,
            batch_size = batch_size, enable_parallel = parallel, soft_weights = acc,
            diagnostics = diagnostics)
        return _psi_seqs(out), _psi_weight_pairs(acc)
    end

    # Same edge from >= 2 windows of one read, with a fractional per-window
    # responsibility, after an earlier read already put nonzero mass on it.
    # That is the float-non-associativity case a per-read fold groups differently
    # from the serial per-window fold. `multi_window_reads >= 1` does not imply
    # this: integer responsibilities associate, and a single window cannot.
    function _psi_windowed_divergence_witness(reads, graph, k, hard)
        prior = Set{Any}()
        for read in reads
            Mycelia._windowed_decode_read_is_long(read, k) || continue
            sink = Vector{Mycelia.Rhizomorph.SoftEdgeWeightAccumulator}()
            Mycelia.improve_read_likelihood_windowed(
                read, graph, k, hard;
                graph_mode = :canonical,
                soft_weights_sink = sink)
            hits = Dict{Any, Int}()
            fractional = Dict{Any, Bool}()
            for staged in sink
                for (edge, weight) in staged.weights
                    hits[edge] = get(hits, edge, 0) + 1
                    if isfinite(weight) && !isinteger(weight)
                        fractional[edge] = true
                    end
                end
            end
            for (edge, n_windows) in hits
                if n_windows >= 2 && get(fractional, edge, false) && (edge in prior)
                    return true
                end
            end
            for staged in sink
                for (edge, weight) in staged.weights
                    weight == 0.0 && continue
                    push!(prior, edge)
                end
            end
        end
        return false
    end

    _psi_guarded() do
        Test.@testset "soft-EM weight sort key survives a repr collision" begin
            a = _PsiReprCollideA(1)
            b = _PsiReprCollideB(1)
            _psi_assert("repr_collision_same_print", repr(a) == repr(b))
            _psi_assert("repr_collision_sort_keys_differ",
                _psi_edge_sort_key(a) != _psi_edge_sort_key(b))
        end
    end

    _psi_guarded() do
        Test.@testset "parallel soft-EM byte-identity (opt1)" begin
            _psi_assert("nthreads", Threads.nthreads() > 1)
            rng = Random.MersenneTwister(2024)
            ref = join(rand(rng, _PSI_BASES, 1500))
            reads = _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
            k = 13
            graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
            hard = Mycelia._hard_vertex_set(graph, k)

            # (a) parallel path ACTUALLY taken — assert via the opt1 actuation counter.
            # The old "ignoring enable_parallel" warning was DELETED with the guard, so
            # grepping the logs for it was vacuous. The counter increments once per
            # parallel decode batch; a serial run leaves it at 0, so the counter
            # DISCRIMINATES parallel from serial. A silent revert-to-serial regression
            # makes the >0 assertion fail even though byte-identity still holds.
            diag = Mycelia.CorrectorDiagnostics()
            seqs_par, wts_par = _psi_run(reads, graph, k, hard; parallel = true,
                diagnostics = diag)
            _psi_assert("parallel_batches_taken", diag.parallel_decode_batches[] > 0)
            serial_diag = Mycelia.CorrectorDiagnostics()
            _psi_run(reads, graph, k, hard; parallel = false, diagnostics = serial_diag)
            _psi_assert("serial_batches_zero", serial_diag.parallel_decode_batches[] == 0)

            # (b) byte-identity vs serial: reads AND soft-EM weights
            seqs_ser, wts_ser = _psi_run(reads, graph, k, hard; parallel = false)
            _psi_assert("weights_nonempty", !isempty(wts_par))
            _psi_assert("seqs_equal", seqs_par == seqs_ser)
            _psi_assert("weights_equal", wts_par == wts_ser)
        end
    end

    # One batch (`batch_size >= n_reads`). The multi-batch cases above never take
    # the `batch_end == total_reads` on the only iteration, so a fold bug that
    # showed up only when every read is in the single parallel batch would pass.
    _psi_guarded() do
        Test.@testset "parallel soft-EM byte-identity — single batch (opt1)" begin
            rng = Random.MersenneTwister(77)
            ref = join(rand(rng, _PSI_BASES, 800))
            n_reads = 40
            reads = _psi_reads(rng, ref; n_reads = n_reads, readlen = 80, n_err = 10)
            k = 13
            graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
            hard = Mycelia._hard_vertex_set(graph, k)
            batch_size = length(reads)
            diag = Mycelia.CorrectorDiagnostics()
            seqs_par, wts_par = _psi_run(
                reads, graph, k, hard; parallel = true, diagnostics = diag,
                batch_size = batch_size)
            _psi_assert("single_batch_count",
                batch_size >= length(reads) && diag.parallel_decode_batches[] == 1)
            seqs_ser, wts_ser = _psi_run(
                reads, graph, k, hard; parallel = false, batch_size = batch_size)
            _psi_assert("single_batch_weights_nonempty", !isempty(wts_par))
            _psi_assert("single_batch_seqs_equal", seqs_par == seqs_ser)
            _psi_assert("single_batch_weights_equal", wts_par == wts_ser)
        end
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
    # The witness assertion below fails if this fixture stops producing that
    # case, so byte-identity cannot pass vacuously on associative weights.
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

    _psi_guarded() do
        Test.@testset "parallel soft-EM byte-identity — windowed decode (opt1)" begin
            _psi_assert("windowed_nthreads", Threads.nthreads() > 1)
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
            _psi_assert("windowed_multi_reads", multi_window_reads >= 1)
            _psi_assert("windowed_divergence_witness",
                _psi_windowed_divergence_witness(reads, graph, k, hard))

            seqs_par, wts_par = _psi_run_windowed(reads, graph, k, hard; parallel = true)
            seqs_ser, wts_ser = _psi_run_windowed(reads, graph, k, hard; parallel = false)
            _psi_assert("windowed_fractional_shared_weight",
                any(p -> isfinite(p[2]) && p[2] > 1 && !isinteger(p[2]), wts_ser))
            _psi_assert("windowed_weights_nonempty", !isempty(wts_par))
            _psi_assert("windowed_seqs_equal", seqs_par == seqs_ser)
            _psi_assert("windowed_weights_equal", wts_par == wts_ser)
        end
    end

    _psi_guarded() do
        Test.@testset "parallel x GC on/off byte-identity (opt1 + opt5)" begin
            rng = Random.MersenneTwister(555)
            ref = join(rand(rng, _PSI_BASES, 1500))
            reads = _psi_reads(rng, ref; n_reads = 200, readlen = 80, n_err = 40)
            k = 13
            graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
            hard = Mycelia._hard_vertex_set(graph, k)
            runit = gc -> begin
                acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
                out, = Mycelia.improve_read_set_likelihood(
                    reads, graph, k; graph_mode = :canonical, skip_solid = true,
                    cheap_correct = true, hard_vertices = hard, decode_enabled = true,
                    batch_size = 50, enable_parallel = true, gc_between_batches = gc,
                    soft_weights = acc)
                (_psi_seqs(out), _psi_weight_pairs(acc))
            end
            off = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
                runit(false)
            end
            on = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
                runit(true)
            end
            _psi_assert("gc_weights_nonempty", !isempty(off[2]))
            _psi_assert("gc_seqs_equal", off[1] == on[1])
            _psi_assert("gc_weights_equal", off[2] == on[2])
        end
    end

    if !isempty(get(ENV, "MYCELIA_PSI_RESULTS", "")) && _psi_failed[]
        exit(1)
    end
end
