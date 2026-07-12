# Quality-Gap Diagnostic: :scalable (canonical, ~145 contigs) vs the earlier
# DoubleStrand full-corrector (~21 contigs, N50 903) on a controlled 1 kb fixture
# =============================================================================
#
# GOAL (read-only w.r.t. src/): isolate WHICH factor of the :scalable corrector
# tier costs assembly quality relative to the earlier winning DoubleStrand route.
# This script changes NO corrector code; it drives the existing public API
# (`Mycelia.mycelia_iterative_assemble` + `Mycelia.Rhizomorph.assemble_genome`)
# across a small grid and tabulates final-assembly quality per config.
#
# Pipeline per config (mirrors src/rhizomorph/assembly.jl::
# _assemble_with_iterative_corrector):
#   1. simulate reads on a fixed 1 kb genome (err=0.01, ~20x, 100 bp) — the SAME
#      read-simulation primitive the correction-validation sweep uses
#      (Mycelia.observe on both strands).
#   2. run `mycelia_iterative_assemble` with the config's knobs -> corrected reads.
#   3. re-assemble the corrected reads via
#      `assemble_genome(corrected; k, corrector=:none)` on DoubleStrand (exactly
#      what the production route does) -> real contigs.
#   4. record contigs / N50 / largest / runtime.
#
# The ONLY thing varied across configs is the CORRECTION step; the re-assembly is
# held fixed (DoubleStrand, same k), so any quality delta is attributable to the
# correction knobs. A naive (no-correction) baseline anchors the top of the table.
#
# Factors (each varied one-at-a-time off the current :scalable baseline):
#   * graph_mode      :canonical (what :scalable forces) vs :doublestrand
#   * skip_solid      true (skip already-solid reads) vs false
#   * k-ladder        n_k_rungs=3 (coarse) vs nothing (prime-by-prime fine walk)
#   * stage0 gates    cheap_correct + hard_window on vs off
#
# Run:
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#     benchmarking/quality_gap_diagnostic.jl
#
# Env knobs (defaults chosen to match the task fixture):
#   MYCELIA_QGD_GENOME_LEN  genome length in bp        (default 1000)
#   MYCELIA_QGD_READLEN     read length in bp          (default 100)
#   MYCELIA_QGD_COVERAGE    target fold coverage       (default 20)
#   MYCELIA_QGD_ERR         per-base error rate        (default 0.01)
#   MYCELIA_QGD_K           re-assembly / corrector k  (default 21)
#   MYCELIA_QGD_SEED        RNG seed                    (default 42)

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end

import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import Dates
import Random
import Logging

# The :canonical corrector floods stderr with `path_to_sequence: no (k-1)
# overlap` reconstruction-fallback warnings — thousands per pass. That I/O flood
# dominates runtime and is itself a symptom (documented separately). Silence
# Warn/Info here so the grid completes; the reconstruction-fallback COUNT is
# captured in a dedicated one-shot probe (see :probe env below).
Logging.disable_logging(Logging.Warn)

# --- config ---------------------------------------------------------------

const GENOME_LEN = parse(Int, get(ENV, "MYCELIA_QGD_GENOME_LEN", "1000"))
const READLEN    = parse(Int, get(ENV, "MYCELIA_QGD_READLEN", "100"))
const COVERAGE   = parse(Float64, get(ENV, "MYCELIA_QGD_COVERAGE", "20"))
const ERR        = parse(Float64, get(ENV, "MYCELIA_QGD_ERR", "0.01"))
const K          = parse(Int, get(ENV, "MYCELIA_QGD_K", "21"))
const SEED       = parse(Int, get(ENV, "MYCELIA_QGD_SEED", "42"))

# --- read simulation (mirrors rhizomorph_correction_validation_sweep.jl) ---

"""
Simulate fixed-length reads by sampling random fragments across both strands of
`refseq` and passing each through `Mycelia.observe` at the given error rate
(Illumina model). Returns a Vector{FASTX.FASTQ.Record}.
"""
function simulate_reads(refseq::BioSequences.LongDNA{4}, readlen::Int,
        coverage::Real, error_rate::Float64, rng::Random.AbstractRNG)
    glen = length(refseq)
    effective_readlen = min(readlen, glen)
    n_reads = max(1, ceil(Int, coverage * glen / effective_readlen))
    records = Vector{FASTX.FASTQ.Record}()
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

# --- one correction+reassembly config -------------------------------------

"""
Run one config: correct `reads` with `mycelia_iterative_assemble` under the given
knobs, then re-assemble the corrected reads on DoubleStrand (corrector=:none) and
measure the resulting contigs. Returns a NamedTuple of metrics. On failure the
row is recorded with ok=false rather than aborting the whole grid.
"""
function run_config(name::String, reads;
        graph_mode::Symbol, skip_solid::Bool, n_k_rungs::Union{Int, Nothing},
        hard_window::Bool, cheap_correct::Bool, soft_em::Bool,
        max_iterations_per_k::Int, beam_width::Union{Int, Nothing})
    t0 = time()
    input_dir = mktempdir()
    output_dir = mktempdir()
    try
        temp_fastq = joinpath(input_dir, "input.fastq")
        open(temp_fastq, "w") do io
            w = FASTX.FASTQ.Writer(io)
            for r in reads
                write(w, r)
            end
            close(w)
        end

        result_dict = Mycelia.mycelia_iterative_assemble(
            temp_fastq;
            max_k = max(K, 13),
            graph_mode = graph_mode,
            skip_solid = skip_solid,
            n_k_rungs = n_k_rungs,
            max_iterations_per_k = max_iterations_per_k,
            hard_window = hard_window,
            soft_em = soft_em,
            cheap_correct = cheap_correct,
            beam_width = beam_width,
            verbose = false,
            enable_checkpointing = false,
            output_dir = output_dir)

        corrected_fastq = get(result_dict[:metadata], :final_fastq_file, nothing)
        if corrected_fastq === nothing || !isfile(corrected_fastq)
            error("corrector produced no :final_fastq_file")
        end
        corrected_reads = open(FASTX.FASTQ.Reader, corrected_fastq) do reader
            collect(reader)
        end
        n_corr = length(corrected_reads)
        n_corr == 0 && error("corrector produced 0 reads")

        # Re-assemble on DoubleStrand (auto-detected DNA mode), corrector=:none —
        # exactly the production re-assembly path. This is held fixed across every
        # config so quality deltas isolate the correction knobs.
        assembly = Mycelia.Rhizomorph.assemble_genome(
            corrected_reads; k = K,
            graph_mode = Mycelia.Rhizomorph.DoubleStrand,
            corrector = :none, verbose = false)

        contigs = assembly.contigs
        n_contigs = length(contigs)

        # Write contigs and compute N50 / largest via the shared metric.
        contigs_path = joinpath(output_dir, "contigs.fasta")
        open(contigs_path, "w") do io
            for (i, c) in enumerate(contigs)
                println(io, ">contig_$(i) length=$(length(c))")
                println(io, c)
            end
        end
        n50 = 0
        largest = 0
        total_len = sum(length.(contigs); init = 0)
        if isfile(contigs_path) && filesize(contigs_path) > 0
            m = Mycelia.assembly_metrics(contigs_path)
            if m !== nothing
                n50 = m.n50
                largest = m.largest_contig
            end
        end

        return (name = name, ok = true, corrected_reads = n_corr,
            n_contigs = n_contigs, n50 = n50, largest = largest,
            total_length = total_len, runtime_s = round(time() - t0; digits = 2))
    catch e
        @warn "config failed" name exception = (e, catch_backtrace())
        return (name = name, ok = false, corrected_reads = 0,
            n_contigs = -1, n50 = -1, largest = -1, total_length = -1,
            runtime_s = round(time() - t0; digits = 2))
    finally
        rm(input_dir; force = true, recursive = true)
        rm(output_dir; force = true, recursive = true)
    end
end

"""
Naive (no-correction) baseline: assemble the raw simulated reads directly on
DoubleStrand. Anchors the top of the table.
"""
function run_naive(name::String, reads)
    t0 = time()
    output_dir = mktempdir()
    try
        assembly = Mycelia.Rhizomorph.assemble_genome(
            reads; k = K, graph_mode = Mycelia.Rhizomorph.DoubleStrand,
            corrector = :none, verbose = false)
        contigs = assembly.contigs
        contigs_path = joinpath(output_dir, "contigs.fasta")
        open(contigs_path, "w") do io
            for (i, c) in enumerate(contigs)
                println(io, ">contig_$(i) length=$(length(c))")
                println(io, c)
            end
        end
        n50 = 0; largest = 0
        if isfile(contigs_path) && filesize(contigs_path) > 0
            m = Mycelia.assembly_metrics(contigs_path)
            if m !== nothing
                n50 = m.n50; largest = m.largest_contig
            end
        end
        return (name = name, ok = true, corrected_reads = length(reads),
            n_contigs = length(contigs), n50 = n50, largest = largest,
            total_length = sum(length.(contigs); init = 0),
            runtime_s = round(time() - t0; digits = 2))
    catch e
        @warn "naive failed" exception = (e, catch_backtrace())
        return (name = name, ok = false, corrected_reads = 0, n_contigs = -1,
            n50 = -1, largest = -1, total_length = -1,
            runtime_s = round(time() - t0; digits = 2))
    finally
        rm(output_dir; force = true, recursive = true)
    end
end

# --- main -----------------------------------------------------------------

function main()
    println("=== Quality-Gap Diagnostic ===")
    println("Start        : $(Dates.now())")
    println("Genome       : $(GENOME_LEN) bp (synthetic, seed=$SEED)")
    println("Reads        : $(READLEN) bp, ~$(COVERAGE)x, err=$(ERR), Illumina")
    println("Re-assembly  : DoubleStrand, k=$K, corrector=:none (held fixed)")
    println()

    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = GENOME_LEN)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    rng = Random.MersenneTwister(SEED)
    reads = simulate_reads(refseq, READLEN, COVERAGE, ERR, rng)
    println("Simulated $(length(reads)) reads.\n")

    # Baseline :scalable knobs (from _corrector_strategy_knobs(:scalable)).
    sc = (graph_mode = :canonical, skip_solid = true, n_k_rungs = 3,
        hard_window = true, cheap_correct = true, soft_em = true,
        max_iterations_per_k = 2, beam_width = nothing)

    df = DataFrames.DataFrame(
        config = String[], ok = Bool[], corrected_reads = Int[],
        contigs = Int[], n50 = Int[], largest = Int[],
        total_length = Int[], runtime_s = Float64[])

    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    out_csv = joinpath(results_dir, "quality_gap_diagnostic.csv")

    # Incremental runner: run one config, append its row to the DataFrame, and
    # re-write the CSV after EACH config so partial results survive a slow/hung
    # later config. Progress is printed to stderr (flushed) for live visibility.
    function record!(r)
        push!(df, (r.name, r.ok, r.corrected_reads, r.n_contigs, r.n50,
            r.largest, r.total_length, r.runtime_s))
        CSV.write(out_csv, df)
        println(stderr, "[done] $(rpad(r.name, 46)) ok=$(r.ok) contigs=$(r.n_contigs) " *
                "n50=$(r.n50) largest=$(r.largest) $(r.runtime_s)s")
        flush(stderr)
        return nothing
    end
    _progress(name) = (println(stderr, "[run ] $name ..."); flush(stderr))

    # 0. naive baseline (no correction).
    _progress("naive (no correction)")
    record!(run_naive("naive (no correction)", reads))

    # 1. current :scalable (canonical + skip + coarse ladder + stage0 gates).
    _progress("scalable (canonical, all gates)")
    record!(run_config("scalable (canonical, all gates)", reads; sc...))

    # 2. flip graph_mode -> doublestrand (isolate the HYPOTHESIS).
    _progress("scalable but DoubleStrand")
    record!(run_config("scalable but DoubleStrand", reads;
        sc..., graph_mode = :doublestrand))

    # 3. flip skip_solid off.
    _progress("scalable, skip_solid OFF")
    record!(run_config("scalable, skip_solid OFF", reads;
        sc..., skip_solid = false))

    # NOTE: the k-ladder factor and the stage-0-gate factor, plus the DoubleStrand
    # full-corrector reference, are run by the companion script
    # `quality_gap_diagnostic_append.jl` and APPENDED to the same CSV. The
    # canonical FINE-LADDER config (scalable, n_k_rungs=nothing, :canonical) was
    # found to be prohibitively slow (>5 min on this single 1 kb config: the
    # prime-by-prime canonical walk from a low initial k reconstructs through the
    # same broken canonical path — its runtime is itself a symptom of the canonical
    # pathology, and it predictably lands at ~197 contigs like configs 1 & 3). The
    # k-ladder factor is instead isolated at DoubleStrand (coarse = config 2 above
    # at 16 contigs, vs the append script's fine-ladder DoubleStrand config), which
    # is fast and equally decisive. See the append script's header for detail.

    println("\n=== RESULTS ===")
    show(stdout, df; allrows = true, allcols = true, truncate = 0)
    println()
    println("\nWrote $out_csv")
    return df
end

main()
