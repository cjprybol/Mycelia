# Focused micro-benchmark isolating the soft-EM competing-paths E-step (C5c,
# td-e70t) and its scaling with graph size. Builds a canonical qualmer graph at
# the dense k=9 rung from simulated reads at increasing genome sizes, then times
# `accumulate_competing_paths!` over all reads (the exact per-read E-step the
# :scalable corrector runs). Reports per-size seconds + a log-log alpha so the
# super-linear residual is visible in isolation (no full corrector run needed).
#
# Run: LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. benchmarking/softem_competing_paths_micro.jl

import Mycelia
import FASTX
import BioSequences
import Random

const READLEN  = 150
const COVERAGE = 20.0
const ERR_RATE = 0.01
const K        = parse(Int, get(ENV, "SEM_K", "9"))
const SEED     = 42
const SIZES    = [parse(Int, s) for s in split(get(ENV, "SEM_SIZES", "5000,10000,20000"), ",")]
const REPEATS  = parse(Int, get(ENV, "SEM_REPEATS", "3"))

function simulate_reads(refseq, readlen, coverage, error_rate, rng)
    glen = length(refseq)
    effective_readlen = min(readlen, glen)
    n_reads = max(1, ceil(Int, coverage * glen / effective_readlen))
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        start = glen == effective_readlen ? 1 : rand(rng, 1:(glen - effective_readlen + 1))
        frag = refseq[start:(start + effective_readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs_seq, quals = Mycelia.observe(frag; error_rate = error_rate, tech = :illumina)
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records
end

function build_fixture(glen)
    rng = Random.MersenneTwister(SEED)
    ref_rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = glen)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, ref_rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR_RATE, rng)
    return reads
end

# Time the E-step over all reads: mirrors the pipeline's per-read
# `accumulate_competing_paths!` call, using a fresh accumulator each repeat.
# `band`/`bound` typemax => unbounded (pre-fix behavior); finite => bounded fix.
function time_estep(reads, graph, k, band, bound)
    best = Inf
    for _ in 1:REPEATS
        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        t = @elapsed for r in reads
            Mycelia.accumulate_competing_paths!(
                acc, r, graph, k; graph_mode = :canonical,
                walk_band = band, successor_bound = bound)
        end
        best = min(best, t)
    end
    return best
end

function report_alpha(label, secs)
    xs = log.(Float64.(SIZES)); ys = log.(secs)
    xbar = sum(xs)/length(xs); ybar = sum(ys)/length(ys)
    alpha = sum((xs .- xbar) .* (ys .- ybar)) / sum((xs .- xbar).^2)
    alpha_lg = (ys[end] - ys[end-1]) / (xs[end] - xs[end-1])
    println("$(label): alpha(loglog)=$(round(alpha, digits=3))  alpha_large=$(round(alpha_lg, digits=3))")
end

function main()
    band = Mycelia._soft_em_walk_band(K)
    bound = Mycelia._SOFT_EM_ALT_SUCCESSOR_BOUND
    println("soft-EM competing-paths micro-benchmark (C5c isolation, k=$(K), band=$(band), succ_bound=$(bound))")
    println("size_bp\tnreads\tnvert\tbefore_s\tafter_s\tspeedup")
    before = Float64[]; after = Float64[]
    for glen in SIZES
        reads = build_fixture(glen)
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, K; mode = :canonical)
        nv = length(collect(Mycelia.MetaGraphsNext.labels(graph)))
        # warmup (compile both paths)
        Mycelia.accumulate_competing_paths!(
            Mycelia.Rhizomorph.SoftEdgeWeightAccumulator(), reads[1], graph, K;
            graph_mode = :canonical, walk_band = band, successor_bound = bound)
        tb = time_estep(reads, graph, K, typemax(Int), typemax(Int))
        ta = time_estep(reads, graph, K, band, bound)
        push!(before, tb); push!(after, ta)
        println("$(glen)\t$(length(reads))\t$(nv)\t$(round(tb, digits=4))\t$(round(ta, digits=4))\t$(round(tb/ta, digits=2))x")
    end
    length(SIZES) >= 2 || return
    report_alpha("BEFORE (unbounded)", before)
    report_alpha("AFTER  (bounded)  ", after)
end

main()
