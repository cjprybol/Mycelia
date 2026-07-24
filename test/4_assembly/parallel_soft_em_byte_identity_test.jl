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
end
