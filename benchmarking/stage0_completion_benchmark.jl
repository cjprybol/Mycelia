# Stage 0 completion benchmark (td-1do7): does the iterative corrector COMPLETE
# once all-solid reads are skipped, and what skip fraction does it report?
#
# The iterative pipeline never completed at real scale because it decoded EVERY
# read (td-ve02). Stage 0 skips reads with no weak k-mer; this benchmark measures
# whether that lets the pipeline finish, and the skip fraction (the work saved).
#
# Toy-scale defaults so it completes locally; env vars scale it up:
#   MYCELIA_S0C_REFLEN / COVERAGE / READLEN / ERR / MAXK / SEED
#   MYCELIA_S0C_SKIP (true|false) — toggle the skip to compare

import Mycelia
import Random

getenv(k, d) = haskey(ENV, k) ? parse(typeof(d), ENV[k]) : d
const REFLEN   = getenv("MYCELIA_S0C_REFLEN", 3000)
const COVERAGE = getenv("MYCELIA_S0C_COVERAGE", 20)
const READLEN  = getenv("MYCELIA_S0C_READLEN", 100)
const ERR      = getenv("MYCELIA_S0C_ERR", 0.01)
const MAXK     = getenv("MYCELIA_S0C_MAXK", 19)
const SEED     = getenv("MYCELIA_S0C_SEED", 42)
const SKIP     = get(ENV, "MYCELIA_S0C_SKIP", "true") == "true"
const BASES = ['A', 'C', 'G', 'T']

function simulate_fastq(path)
    rng = Random.MersenneTwister(SEED)
    ref = join(rand(rng, BASES, REFLEN))
    n_reads = max(1, round(Int, COVERAGE * REFLEN / READLEN))
    open(path, "w") do io
        for i in 1:n_reads
            s = rand(rng, 1:(REFLEN - READLEN + 1))
            seq = collect(ref[s:(s + READLEN - 1)])
            quals = fill('I', READLEN)
            for j in 1:READLEN
                if rand(rng) < ERR
                    seq[j] = rand(rng, filter(!=(seq[j]), BASES))
                    quals[j] = Char(clamp(round(Int, 22 + 8 * randn(rng)), 2, 40) + 33)
                end
            end
            println(io, "@r$i"); println(io, join(seq)); println(io, "+"); println(io, join(quals))
        end
    end
    return n_reads
end

function main()
    outdir = mktempdir()
    fastq = joinpath(outdir, "reads.fastq")
    n_reads = simulate_fastq(fastq)
    println("=== Stage 0 completion benchmark ===")
    println("reflen=$REFLEN cov=$COVERAGE readlen=$READLEN err=$ERR max_k=$MAXK seed=$SEED skip_solid=$SKIP reads=$n_reads")
    t0 = time()
    Mycelia.mycelia_iterative_assemble(fastq;
        max_k = MAXK,
        skip_solid = SKIP,
        verbose = true,
        enable_parallel = true,
        enable_checkpointing = false,
        output_dir = joinpath(outdir, "asm"))
    elapsed = time() - t0
    println("=== COMPLETED in $(round(elapsed; digits = 1))s (skip_solid=$SKIP) ===")
    println("(skip fraction is reported per iteration in the verbose output above)")
end

main()
