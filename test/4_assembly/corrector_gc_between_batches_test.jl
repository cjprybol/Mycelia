# GC-BETWEEN-BATCHES KEYWORD (td-3xob / td-jbjd opt5)
# ==================================================
#
# The :scalable corrector previously forced a stop-the-world GC.gc() between every
# read batch, parking all @threads decode workers (the mechanism most consistent
# with the observed ~286% CPU on a flat/bounded memory plateau) while buying nothing
# once memory is already bounded. opt5 makes that GC opt-in via the
# `gc_between_batches` keyword (primary) with a MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES
# env-var fallback for memory-constrained hosts.
#
# GC is OUTPUT-NEUTRAL, so corrected reads MUST be byte-identical whether the
# between-batch GC fires or not. This is the opt5 byte-identity lock. The keyword
# must also gate the GC path INDEPENDENT of the env var.
#
# WHAT IS ASSERTED
#   * Unit (improve_read_set_likelihood), multi-batch (batch_size < n_reads):
#     gc_between_batches=false vs =true => byte-identical corrected reads.
#   * Env independence: with the env var force-unset AND with it force-"true",
#     the keyword value does not change corrected output (output-neutral either way).
#   * Integration (mycelia_iterative_assemble), forced small batch_size:
#     gc_between_batches=false vs =true => byte-identical final corrected fastq,
#     proving the keyword threads through public entry -> wrapper -> impl.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_gc_between_batches_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _GC_BASES = ['A', 'C', 'G', 'T']

# Reads sampled from ref, with a single substitution injected into the first
# n_err of them, so Stage 0 + the per-read decode both have work across batches.
function _gc_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _GC_BASES))
        end
        push!(records,
            FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

_seqs(records) = [FASTX.sequence(String, r) for r in records]

Test.@testset "corrector gc_between_batches keyword (td-3xob opt5)" begin
    R = Mycelia.Rhizomorph

    Test.@testset "unit: multi-batch decode is byte-identical GC on/off" begin
        rng = Random.MersenneTwister(4242)
        ref = join(rand(rng, _GC_BASES, 1200))
        reads = _gc_reads(rng, ref; n_reads = 180, readlen = 80, n_err = 30)
        k = 13
        graph = R.build_qualmer_graph(reads, k; mode = :canonical)
        hard = Mycelia._hard_vertex_set(graph, k)

        # batch_size < n_reads => multiple batches => the between-batch GC site
        # fires when enabled (batch_end < total_reads on batches 1..n-1).
        run_unit = gc -> Mycelia.improve_read_set_likelihood(
            reads, graph, k; graph_mode = :canonical, skip_solid = true,
            cheap_correct = true, hard_vertices = hard, decode_enabled = true,
            batch_size = 50, gc_between_batches = gc)

        out_off, = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            run_unit(false)
        end
        out_on, = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            run_unit(true)
        end
        # byte-identical: GC is output-neutral
        Test.@test _seqs(out_off) == _seqs(out_on)

        # Env-var independence: env force-"true" with keyword=false is still
        # output-identical to the keyword=false / env-unset baseline.
        out_env, = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => "true") do
            run_unit(false)
        end
        Test.@test _seqs(out_env) == _seqs(out_off)
    end

    Test.@testset "integration: mycelia_iterative_assemble threads the keyword" begin
        rng = Random.MersenneTwister(909)
        ref = join(rand(rng, _GC_BASES, 900))
        reads = _gc_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        run_full = (; gc, out) -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 2,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            gc_between_batches = gc, verbose = false,
            enable_checkpointing = false, output_dir = out)

        r_off = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            run_full(; gc = false, out = joinpath(tmp, "off"))
        end
        r_on = Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            run_full(; gc = true, out = joinpath(tmp, "on"))
        end

        final_seqs = res -> _seqs(
            open(FASTX.FASTQ.Reader, res[:metadata][:final_fastq_file]) do rd
                collect(rd)
            end)
        Test.@test final_seqs(r_off) == final_seqs(r_on)
    end

    # WIRING LOCK: byte-identity alone can't prove the keyword gates the GC
    # (the GC is output-invisible). Assert the pure gate predicate directly.
    Test.@testset "gate predicate truth table (opt5 wiring lock)" begin
        gc_enabled = Mycelia._gc_between_batches_enabled
        # keyword is primary: true regardless of env
        Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => nothing) do
            Test.@test gc_enabled(true) == true
            Test.@test gc_enabled(false) == false
        end
        Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => "false") do
            Test.@test gc_enabled(true) == true      # keyword wins over env
            Test.@test gc_enabled(false) == false
        end
        # env fallback enables when the keyword is false
        for tok in ("1", "true", "yes")
            Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => tok) do
                Test.@test gc_enabled(false) == true
            end
        end
        # env rejects everything else (documented case-sensitive allow-list)
        for tok in ("0", "no", "TRUE", "on", "", "garbage")
            Base.withenv("MYCELIA_CORRECTOR_GC_BETWEEN_BATCHES" => tok) do
                Test.@test gc_enabled(false) == false
            end
        end
    end
end
