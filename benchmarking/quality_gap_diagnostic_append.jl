# Supplement to quality_gap_diagnostic.jl — runs the two remaining FAST configs
# and APPENDS them to benchmarking/results/quality_gap_diagnostic.csv.
# ================================================================================
#
# The main grid's canonical FINE-LADDER config (scalable, n_k_rungs=nothing) is
# prohibitively slow (>5 min on a single 1 kb config: the prime-by-prime canonical
# walk from a low initial k reconstructs through the same broken canonical path as
# the coarse-ladder canonical config, so it predictably lands at ~197 contigs like
# configs 1 & 3 — its extreme runtime is itself a symptom of canonical's
# reconstruction pathology, not a new quality signal). This supplement therefore
# runs only the two remaining fast, informative configs:
#
#   * config 5 — scalable, STAGE-0 GATES OFF (cheap_correct + hard_window off),
#                still canonical + coarse ladder. Tests the stage-0/hard-window
#                factor. ~50 s.
#   * config 6 — DoubleStrand full-corrector: doublestrand, no skip, FINE ladder,
#                no stage-0 gates, 3 iters, auto-beam. Serves double duty: the
#                historical-winner reference AND the k-ladder-at-DoubleStrand test
#                (coarse ladder = the main grid's config 2 at 16 contigs, vs this
#                fine ladder) — if both are ~16, the ladder is not the culprit.
#
# Reuses the exact read simulation + correct-and-reassemble pipeline from the main
# script by `include`-ing it is not possible (that file calls main() at load), so
# the needed helpers are duplicated verbatim here. Read-only w.r.t. src/.

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end
import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import Random
import Logging
Logging.disable_logging(Logging.Warn)

const GENOME_LEN = parse(Int, get(ENV, "MYCELIA_QGD_GENOME_LEN", "1000"))
const READLEN    = parse(Int, get(ENV, "MYCELIA_QGD_READLEN", "100"))
const COVERAGE   = parse(Float64, get(ENV, "MYCELIA_QGD_COVERAGE", "20"))
const ERR        = parse(Float64, get(ENV, "MYCELIA_QGD_ERR", "0.01"))
const K          = parse(Int, get(ENV, "MYCELIA_QGD_K", "21"))
const SEED       = parse(Int, get(ENV, "MYCELIA_QGD_SEED", "42"))

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

function run_config(name, reads; graph_mode, skip_solid, n_k_rungs, hard_window,
        cheap_correct, soft_em, max_iterations_per_k, beam_width)
    t0 = time(); input_dir = mktempdir(); output_dir = mktempdir()
    try
        temp_fastq = joinpath(input_dir, "input.fastq")
        open(temp_fastq, "w") do io
            w = FASTX.FASTQ.Writer(io); for r in reads; write(w, r); end; close(w)
        end
        rd = Mycelia.mycelia_iterative_assemble(temp_fastq;
            max_k = max(K, 13), graph_mode = graph_mode, skip_solid = skip_solid,
            n_k_rungs = n_k_rungs, max_iterations_per_k = max_iterations_per_k,
            hard_window = hard_window, soft_em = soft_em, cheap_correct = cheap_correct,
            beam_width = beam_width, verbose = false, enable_checkpointing = false,
            output_dir = output_dir)
        cf = get(rd[:metadata], :final_fastq_file, nothing)
        (cf === nothing || !isfile(cf)) && error("no final_fastq_file")
        creads = open(FASTX.FASTQ.Reader, cf) do r; collect(r); end
        length(creads) == 0 && error("0 corrected reads")
        asm = Mycelia.Rhizomorph.assemble_genome(creads; k = K,
            graph_mode = Mycelia.Rhizomorph.DoubleStrand, corrector = :none, verbose = false)
        contigs = asm.contigs
        cp = joinpath(output_dir, "contigs.fasta")
        open(cp, "w") do io
            for (i, c) in enumerate(contigs); println(io, ">contig_$(i)"); println(io, c); end
        end
        n50 = 0; largest = 0
        if isfile(cp) && filesize(cp) > 0
            m = Mycelia.assembly_metrics(cp)
            m !== nothing && (n50 = m.n50; largest = m.largest_contig)
        end
        return (name = name, ok = true, corrected_reads = length(creads),
            n_contigs = length(contigs), n50 = n50, largest = largest,
            total_length = sum(length.(contigs); init = 0), runtime_s = round(time() - t0; digits = 2))
    catch e
        @warn "config failed" name exception = (e, catch_backtrace())
        return (name = name, ok = false, corrected_reads = 0, n_contigs = -1,
            n50 = -1, largest = -1, total_length = -1, runtime_s = round(time() - t0; digits = 2))
    finally
        rm(input_dir; force = true, recursive = true); rm(output_dir; force = true, recursive = true)
    end
end

function main()
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = GENOME_LEN)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR, Random.MersenneTwister(SEED))

    out_csv = joinpath(@__DIR__, "results", "quality_gap_diagnostic.csv")
    df = CSV.read(out_csv, DataFrames.DataFrame)

    function record!(r)
        push!(df, (r.name, r.ok, r.corrected_reads, r.n_contigs, r.n50, r.largest,
            r.total_length, r.runtime_s))
        CSV.write(out_csv, df)
        println(stderr, "[done] $(rpad(r.name, 46)) ok=$(r.ok) contigs=$(r.n_contigs) " *
                "n50=$(r.n50) largest=$(r.largest) $(r.runtime_s)s"); flush(stderr)
    end

    println(stderr, "[run ] scalable, stage0 gates OFF ..."); flush(stderr)
    record!(run_config("scalable, stage0 gates OFF", reads;
        graph_mode = :canonical, skip_solid = true, n_k_rungs = 3,
        hard_window = false, cheap_correct = false, soft_em = true,
        max_iterations_per_k = 2, beam_width = nothing))

    # NOTE (k-ladder factor): the fine-ladder DoubleStrand config below was also
    # found to be too slow to complete quickly (the prime-by-prime walk from a low
    # initial k is expensive even on DoubleStrand). It is left here, DISABLED, for
    # reproducibility. It is NOT needed to isolate the k-ladder factor: configs 1
    # (canonical, COARSE ladder = 197 contigs) and 2 (DoubleStrand, COARSE ladder =
    # 16 contigs) hold the ladder CONSTANT and differ only in graph_mode, so the
    # 197->16 collapse is attributable to graph_mode alone, independent of ladder.
    # Enable this to additionally confirm fine-vs-coarse at DoubleStrand.
    #
    # println(stderr, "[run ] DoubleStrand full-corrector (fine ladder) ..."); flush(stderr)
    # record!(run_config("DoubleStrand full-corrector (fine ladder)", reads;
    #     graph_mode = :doublestrand, skip_solid = false, n_k_rungs = nothing,
    #     hard_window = false, cheap_correct = false, soft_em = false,
    #     max_iterations_per_k = 3, beam_width = nothing))

    println("\n=== FINAL RESULTS ===")
    show(stdout, df; allrows = true, allcols = true, truncate = 0)
    println("\nWrote $out_csv")
end

main()
