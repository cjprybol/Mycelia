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
# Regression guard (td-gz51): a FLOOR on the Stage 0 skip fraction. The skip is
# what lets the iterative corrector complete at scale — if a regression silently
# disables it (classifier returns all-weak, inverted no-evidence fallback), the
# pipeline still finishes, just far slower, so a wall-clock benchmark would NOT
# catch it. This asserts that at every k >= 5 the reported skip fraction stays
# above the floor, turning a silent skip regression into a benchmark FAILURE.
const SKIP_FLOOR = getenv("MYCELIA_S0C_SKIP_FLOOR", 30.0)  # percent
# k below this is exempt: at tiny k almost every read has a weak k-mer, so a low
# skip fraction there is expected, not a regression.
const SKIP_FLOOR_MIN_K = getenv("MYCELIA_S0C_SKIP_FLOOR_MIN_K", 5)
# Self-test hook: inject a synthetic k=7 skip=0/100 record to prove the floor
# assertion actually fires (verification only; off by default, non-breaking).
const SELFTEST_BAD = get(ENV, "MYCELIA_S0C_SELFTEST_BAD", "false") in ("1", "true")
const BASES = ['A', 'C', 'G', 'T']

"""Run `f` with stdout captured line-by-line while still echoing live to the
real stdout, and return the captured lines (so the verbose Stage 0 output can be
parsed for the skip-fraction guard without silencing the benchmark)."""
function run_capturing_stdout(f)
    real_out = stdout
    rd, wr = redirect_stdout()
    lines = String[]
    reader = @async for ln in eachline(rd)
        println(real_out, ln)   # live passthrough
        push!(lines, ln)
    end
    try
        f()
    finally
        redirect_stdout(real_out)
        close(wr)
        wait(reader)
    end
    return lines
end

"""Parse captured verbose output into (k, skipped, total) records, associating
each "Stage 0 skipped" line with the most recently announced k."""
function parse_skip_records(lines::Vector{String})
    records = Tuple{Int, Int, Int}[]
    current_k = 0
    for ln in lines
        mk = match(r"for k=(\d+)", ln)
        mk === nothing || (current_k = parse(Int, mk.captures[1]))
        mb = match(r"Building qualmer graph with k=(\d+)", ln)
        mb === nothing || (current_k = parse(Int, mb.captures[1]))
        ms = match(r"Stage 0 skipped \(all-solid, no decode\): (\d+)/(\d+)", ln)
        ms === nothing && continue
        push!(records, (current_k, parse(Int, ms.captures[1]), parse(Int, ms.captures[2])))
    end
    return records
end

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
    captured = run_capturing_stdout() do
        Mycelia.mycelia_iterative_assemble(fastq;
            max_k = MAXK,
            skip_solid = SKIP,
            verbose = true,
            enable_parallel = true,
            enable_checkpointing = false,
            output_dir = joinpath(outdir, "asm"))
    end
    elapsed = time() - t0
    println("=== COMPLETED in $(round(elapsed; digits = 1))s (skip_solid=$SKIP) ===")

    # ---- Skip-fraction floor guard (td-gz51) --------------------------------
    # Only meaningful when the skip was actually requested; with skip off the
    # source prints no "Stage 0 skipped" lines and there is nothing to guard.
    records = parse_skip_records(captured)
    if SELFTEST_BAD
        push!(records, (7, 0, 100))
        println("[SELFTEST] injected synthetic k=7 skip=0/100 record to prove the floor assertion fires")
    end
    if SKIP
        guarded = filter(r -> r[1] >= SKIP_FLOOR_MIN_K && r[3] > 0, records)
        println("\n--- Stage 0 skip-fraction guard (floor = $(SKIP_FLOOR)% for k >= $(SKIP_FLOOR_MIN_K)) ---")
        violations = String[]
        if isempty(guarded)
            push!(violations,
                "no Stage 0 skip records observed for k >= $(SKIP_FLOOR_MIN_K) while skip_solid=true " *
                "(the skip may be silently disabled)")
        else
            for (k, skipped, total) in guarded
                frac = skipped / total * 100
                status = frac >= SKIP_FLOOR ? "OK " : "LOW"
                println("  [$status] k=$k skip=$skipped/$total ($(round(frac; digits = 1))%)")
                frac < SKIP_FLOOR && push!(violations,
                    "k=$k skip fraction $(round(frac; digits = 1))% < floor $(SKIP_FLOOR)%")
            end
        end
        if !isempty(violations)
            println("\n=== GUARD FAILED: Stage 0 skip fraction below floor ===")
            for v in violations
                println("  - $v")
            end
            error("Stage 0 skip-fraction floor guard failed: " * join(violations, "; "))
        end
        println("=== GUARD PASSED: all k >= $(SKIP_FLOOR_MIN_K) skip fractions >= $(SKIP_FLOOR)% ===")
    else
        println("(skip_solid=false — skip-fraction guard not applicable this run)")
    end
end

main()
