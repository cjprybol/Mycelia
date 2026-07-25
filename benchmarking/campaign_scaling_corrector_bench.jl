#!/usr/bin/env julia
# corrector_bench.jl — standalone cross-checkout corrector wall-clock timer.
#
# Uses ONLY the stable public `Mycelia.mycelia_iterative_assemble` API present
# byte-identically (by signature) in both the pre-campaign baseline (5b7b495e)
# and the post-campaign (opt5+opt1+opt2) checkout. Does NOT pass
# gc_between_batches or MYCELIA_RCA_BATCH_SIZE (opt5-only knobs). `enable_parallel`
# IS passed explicitly (via BENCH_ENABLE_PARALLEL) because its *auto-on-when-
# nthreads>1* default lives one layer up, in assemble_genome's :scalable
# `_corrector_strategy_knobs` (introduced by opt1, PR #434) — NOT in
# mycelia_iterative_assemble's own keyword default, which is `false` in BOTH
# checkouts. Passing enable_parallel explicitly exercises opt1's thread-safe
# decode path directly through the stable API without touching any
# opt5/opt1-only kwarg name.
#
# ENV:
#   BENCH_PROJECT           path to the checkout being timed (informational only)
#   BENCH_FASTQ             input FASTQ (from gen_fixture.jl)
#   BENCH_OUTDIR            output_dir for mycelia_iterative_assemble
#   BENCH_MAX_K             (default 21)
#   BENCH_N_K_RUNGS         (default 3)
#   BENCH_MAX_ITER          (default 2)
#   BENCH_ENABLE_PARALLEL   "true"/"false" (default "false")
#   BENCH_WARMUP_GENOME_LEN tiny warm-up genome length (default 300)
#   BENCH_RESULT_JSON       path to write a small JSON result summary

import Mycelia
import Random
import FASTX
import SHA
import Dates

function tiny_fastq(path::String; genome_len::Int = 300, n_reads::Int = 40, seed::Int = 1)
    rng = Random.MersenneTwister(seed)
    bases = ('A', 'C', 'G', 'T')
    genome = join(rand(rng, collect(bases), genome_len))
    readlen = min(80, genome_len)
    open(path, "w") do io
        writer = FASTX.FASTQ.Writer(io)
        for i = 1:n_reads
            start = rand(rng, 1:(genome_len-readlen+1))
            frag = collect(genome[start:(start+readlen-1)])
            for j in eachindex(frag)
                if rand(rng) < 0.01
                    alts = filter(!=(frag[j]), collect(bases))
                    frag[j] = alts[rand(rng, 1:length(alts))]
                end
            end
            qual = String(fill(Char(20 + 33), length(frag)))
            write(writer, FASTX.FASTQ.Record("warm_$(i)", String(frag), qual))
        end
        close(writer)
    end
    return path
end

function sha256_of(path::String)
    return open(io -> bytes2hex(SHA.sha256(io)), path, "r")
end

function run_corrector(
    fastq::String,
    outdir::String;
    max_k::Int,
    n_k_rungs::Int,
    max_iterations_per_k::Int,
    enable_parallel::Bool,
)
    isdir(outdir) && rm(outdir; recursive = true, force = true)
    return Mycelia.mycelia_iterative_assemble(
        fastq;
        max_k = max_k,
        n_k_rungs = n_k_rungs,
        max_iterations_per_k = max_iterations_per_k,
        graph_mode = :canonical,
        skip_solid = true,
        cheap_correct = true,
        hard_window = true,
        soft_em = true,
        enable_parallel = enable_parallel,
        verbose = false,
        enable_checkpointing = false,
        output_dir = outdir,
    )
end

max_k = parse(Int, get(ENV, "BENCH_MAX_K", "21"))
n_k_rungs = parse(Int, get(ENV, "BENCH_N_K_RUNGS", "3"))
max_iter = parse(Int, get(ENV, "BENCH_MAX_ITER", "2"))
enable_parallel = lowercase(get(ENV, "BENCH_ENABLE_PARALLEL", "false")) == "true"
fastq = ENV["BENCH_FASTQ"]
outdir = ENV["BENCH_OUTDIR"]
result_json = get(ENV, "BENCH_RESULT_JSON", "")

println("=== corrector_bench.jl ===")
println("julia threads      : $(Threads.nthreads())")
println("project             : $(get(ENV, "BENCH_PROJECT", "?"))")
println("fastq               : $fastq")
println("max_k/n_k_rungs/iter: $max_k/$n_k_rungs/$max_iter")
println("enable_parallel     : $enable_parallel")

# --- Warm-up (compilation only; NOT timed) ----------------------------------
warm_dir = mktempdir()
warm_fastq = tiny_fastq(joinpath(warm_dir, "warm.fastq"))
println("warming up on tiny fixture ($(warm_fastq)) ...")
t_warm = time()
run_corrector(
    warm_fastq,
    joinpath(warm_dir, "warm_out");
    max_k = min(max_k, 11),
    n_k_rungs = 1,
    max_iterations_per_k = 1,
    enable_parallel = enable_parallel,
)
println(
    "warm-up wall time (compile-inclusive, excluded from measured runtime): $(round(time() - t_warm, digits=2))s",
)

# --- Timed real run ----------------------------------------------------------
GC.gc()
t0 = time()
result = run_corrector(
    fastq,
    outdir;
    max_k = max_k,
    n_k_rungs = n_k_rungs,
    max_iterations_per_k = max_iter,
    enable_parallel = enable_parallel,
)
runtime_s = time() - t0

final_fastq = result[:metadata][:final_fastq_file]
final_sha = sha256_of(final_fastq)
n_out_reads = count(_ -> true, FASTX.FASTQ.Reader(open(final_fastq)))

println(
    "RESULT runtime_s=$(round(runtime_s, digits=3)) final_fastq=$final_fastq sha256=$final_sha n_out_reads=$n_out_reads",
)

if !isempty(result_json)
    open(result_json, "w") do io
        write(
            io,
            """
  {
    "timestamp": "$(Dates.now())",
    "project": "$(get(ENV, "BENCH_PROJECT", "?"))",
    "nthreads": $(Threads.nthreads()),
    "enable_parallel": $(enable_parallel),
    "max_k": $max_k,
    "n_k_rungs": $n_k_rungs,
    "max_iterations_per_k": $max_iter,
    "fastq": "$fastq",
    "runtime_s": $(runtime_s),
    "final_fastq": "$final_fastq",
    "final_fastq_sha256": "$final_sha",
    "n_out_reads": $n_out_reads
  }
  """,
        )
    end
end
