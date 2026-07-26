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
#   BENCH_FASTQ             input FASTQ (from gen_fixture.jl)   (required)
#   BENCH_OUTDIR            output_dir for mycelia_iterative_assemble (required)
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

function tiny_fastq(path::String; genome_len::Int = 300, n_reads::Int = 40, seed::Int = 1)::String
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

# Resolve + validate the output dir OUTSIDE the timed region. NEVER rm -rf a
# populated pre-existing directory: a typo in the caller-provided BENCH_OUTDIR
# (e.g. /tmp, or the directory holding BENCH_FASTQ) must not destroy unrelated
# files. The target must be fresh (non-existent) or empty — safe by
# construction, not by blocklist. The corrector creates the dir if absent.
function prepare_outdir(outdir::String)::String
    isempty(outdir) && error("BENCH_OUTDIR is empty")
    ap = abspath(outdir)
    if ap == "/" || ap == dirname(ap) || ap == abspath(homedir()) || ap == abspath(pwd())
        error("BENCH_OUTDIR ($ap) is unsafe; use a dedicated benchmark subdirectory")
    end
    if ispath(ap) && !isempty(readdir(ap))
        error("BENCH_OUTDIR ($ap) already exists and is non-empty; point at a fresh or empty dedicated directory (refusing to delete pre-existing content)")
    end
    return ap
end

# Minimal JSON string escaping — paths can contain quotes/backslashes.
json_str(s) = "\"" * replace(
    string(s),
    "\\" => "\\\\",
    "\"" => "\\\"",
    "\n" => "\\n",
    "\r" => "\\r",
    "\t" => "\\t",
) * "\""

function run_corrector(
    fastq::String,
    outdir::String;
    max_k::Int,
    n_k_rungs::Int,
    max_iterations_per_k::Int,
    enable_parallel::Bool,
)
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
_ep_raw = lowercase(strip(get(ENV, "BENCH_ENABLE_PARALLEL", "false")))
enable_parallel =
    _ep_raw == "true" ? true :
    _ep_raw == "false" ? false :
    error("BENCH_ENABLE_PARALLEL must be \"true\" or \"false\", got: $(repr(_ep_raw))")
fastq = ENV["BENCH_FASTQ"]
outdir = ENV["BENCH_OUTDIR"]
result_json = get(ENV, "BENCH_RESULT_JSON", "")

println("=== corrector_bench.jl ===")
println("julia threads      : $(Threads.nthreads())")
println("project             : $(get(ENV, "BENCH_PROJECT", "?"))")
println("fastq               : $fastq")
println("max_k/n_k_rungs/iter: $max_k/$n_k_rungs/$max_iter")
println("enable_parallel     : $enable_parallel")
if enable_parallel && Threads.nthreads() == 1
    @warn "enable_parallel=true but Julia started single-threaded (-t1); the parallel decode path has no worker pool, so this run is effectively serial."
end

# --- Warm-up (compilation only; NOT timed) ----------------------------------
warm_dir = mktempdir()
try
    warm_fastq = tiny_fastq(
        joinpath(warm_dir, "warm.fastq");
        genome_len = parse(Int, get(ENV, "BENCH_WARMUP_GENOME_LEN", "300")),
    )
    println("warming up on tiny fixture ($(warm_fastq)) ...")
    t_warm = time()
    run_corrector(
        warm_fastq,
        prepare_outdir(joinpath(warm_dir, "warm_out"));
        max_k = min(max_k, 11),
        n_k_rungs = 1,
        max_iterations_per_k = 1,
        enable_parallel = enable_parallel,
    )
    println(
        "warm-up wall time (compile-inclusive, excluded from measured runtime): $(round(time() - t_warm, digits=2))s",
    )
finally
    rm(warm_dir; recursive = true, force = true)
end

# --- Input provenance + outdir prep (BEFORE timing) --------------------------
clean_outdir = prepare_outdir(outdir)
input_sha = sha256_of(fastq)
n_in_reads = open(FASTX.FASTQ.Reader, fastq) do r
    count(_ -> true, r)
end
println("input_fastq_sha256  : $input_sha  ($n_in_reads reads)")

# --- Timed real run ----------------------------------------------------------
GC.gc()
t0 = time()
result = run_corrector(
    fastq,
    clean_outdir;
    max_k = max_k,
    n_k_rungs = n_k_rungs,
    max_iterations_per_k = max_iter,
    enable_parallel = enable_parallel,
)
runtime_s = time() - t0

final_fastq = result[:metadata][:final_fastq_file]
final_sha = sha256_of(final_fastq)
n_out_reads = open(FASTX.FASTQ.Reader, final_fastq) do r
    count(_ -> true, r)
end

# A corrector that returns normally but emits an empty/degenerate FASTQ
# would still hash cleanly and read as a "successful" benchmark point.
n_out_reads > 0 ||
    error("corrector produced $n_out_reads output reads (degenerate result) — refusing to report a runtime/SHA for a failed run")
if n_out_reads < n_in_reads ÷ 2
    @warn "output read count ($n_out_reads) is <50% of input ($n_in_reads) — verify this is expected before trusting the timing."
end

println(
    "RESULT runtime_s=$(round(runtime_s, digits=3)) final_fastq=$final_fastq sha256=$final_sha n_out_reads=$n_out_reads n_in_reads=$n_in_reads input_sha256=$input_sha",
)

if !isempty(result_json)
    open(result_json, "w") do io
        write(
            io,
            """
  {
    "timestamp": $(json_str(Dates.now())),
    "project": $(json_str(get(ENV, "BENCH_PROJECT", "?"))),
    "nthreads": $(Threads.nthreads()),
    "enable_parallel": $(enable_parallel),
    "max_k": $max_k,
    "n_k_rungs": $n_k_rungs,
    "max_iterations_per_k": $max_iter,
    "fastq": $(json_str(fastq)),
    "input_fastq_sha256": $(json_str(input_sha)),
    "n_in_reads": $n_in_reads,
    "runtime_s": $(runtime_s),
    "final_fastq": $(json_str(final_fastq)),
    "final_fastq_sha256": $(json_str(final_sha)),
    "n_out_reads": $n_out_reads
  }
  """,
        )
    end
end
