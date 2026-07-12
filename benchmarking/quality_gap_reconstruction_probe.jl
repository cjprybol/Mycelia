# Reconstruction-fallback probe for the quality-gap diagnostic.
# =============================================================
#
# Counts how many `path_to_sequence: no (k-1) overlap for canonical k-mer ...`
# reconstruction-fallback warnings the corrector emits under :canonical vs
# :doublestrand graph_mode on the SAME 1 kb reads. A high canonical count is the
# mechanism by which :scalable (which forces :canonical) mangles corrected reads.
#
# Run:
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#     benchmarking/quality_gap_reconstruction_probe.jl

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end
import Mycelia
import FASTX
import BioSequences
import Random
import Logging

const GENOME_LEN = parse(Int, get(ENV, "MYCELIA_QGD_GENOME_LEN", "1000"))
const READLEN    = parse(Int, get(ENV, "MYCELIA_QGD_READLEN", "100"))
const COVERAGE   = parse(Float64, get(ENV, "MYCELIA_QGD_COVERAGE", "20"))
const ERR        = parse(Float64, get(ENV, "MYCELIA_QGD_ERR", "0.01"))
const K          = parse(Int, get(ENV, "MYCELIA_QGD_K", "21"))
const SEED       = parse(Int, get(ENV, "MYCELIA_QGD_SEED", "42"))

# A logger that counts (and swallows) warnings whose message mentions the
# canonical reconstruction fallback; everything else is discarded silently.
mutable struct FallbackCounter <: Logging.AbstractLogger
    count::Int
end
Logging.min_enabled_level(::FallbackCounter) = Logging.Warn
Logging.shouldlog(::FallbackCounter, args...) = true
function Logging.handle_message(l::FallbackCounter, level, message, args...; kwargs...)
    if occursin("no (k-1) overlap", string(message))
        l.count += 1
    end
    return nothing
end

function simulate_reads(refseq, readlen, coverage, error_rate, rng)
    glen = length(refseq)
    n_reads = max(1, ceil(Int, coverage * glen / readlen))
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        start = rand(rng, 1:(glen - readlen + 1))
        frag = refseq[start:(start + readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs_seq, quals = Mycelia.observe(frag; error_rate = error_rate, tech = :illumina)
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records
end

function count_fallbacks(reads, graph_mode::Symbol)
    input_dir = mktempdir(); output_dir = mktempdir()
    temp_fastq = joinpath(input_dir, "input.fastq")
    open(temp_fastq, "w") do io
        w = FASTX.FASTQ.Writer(io)
        for r in reads; write(w, r); end
        close(w)
    end
    counter = FallbackCounter(0)
    Logging.with_logger(counter) do
        Mycelia.mycelia_iterative_assemble(temp_fastq;
            max_k = max(K, 13), graph_mode = graph_mode, skip_solid = true,
            n_k_rungs = 3, max_iterations_per_k = 2, hard_window = true,
            soft_em = true, cheap_correct = true, beam_width = nothing,
            verbose = false, enable_checkpointing = false, output_dir = output_dir)
    end
    rm(input_dir; force = true, recursive = true)
    rm(output_dir; force = true, recursive = true)
    return counter.count
end

function main()
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = GENOME_LEN)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR, Random.MersenneTwister(SEED))
    println("Reads: $(length(reads))  (1 kb genome, err=$ERR, ~$(COVERAGE)x)")
    c_canon = count_fallbacks(reads, :canonical)
    c_dbl   = count_fallbacks(reads, :doublestrand)
    println("Reconstruction-fallback warnings under :canonical    = $c_canon")
    println("Reconstruction-fallback warnings under :doublestrand = $c_dbl")
end

main()
