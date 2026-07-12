# e2e_phase_profile.jl
#
# End-to-end PHASE breakdown for the Rhizomorph correction-validation sweep's
# ITERATIVE arm — see bead td-9q84.
#
# WHY: after the O(V^2)->O(V) fixes, the 5 kb :scalable sweep COMPLETES
# (~22 min, was >3 h), producing a near-complete assembly. Before choosing the
# NEXT optimization target we need to know WHERE the ~22 min actually goes:
# is it still the CORRECTOR (the thing we have been optimizing), or is it the
# benchmark's re-assembly / scoring scaffolding around it?
#
# WHAT THE SWEEP TIMES (benchmarking/rhizomorph_correction_validation_sweep.jl):
#   run_arm() sets t0, calls assemble_genome(reads; corrector=:iterative), then
#   captures runtime = time() - t0. assembly_metrics() + genome_fraction run
#   AFTER that, so the reported iterative `runtime_s` is EXACTLY
#   assemble_genome(:iterative). Under the hood
#   (src/rhizomorph/assembly.jl::_assemble_with_iterative_corrector) that is:
#     (a) CORRECTION  — Mycelia.mycelia_iterative_assemble(..., strategy :scalable knobs)
#     (b) RE-ASSEMBLY — assemble_genome(corrected_reads; corrector=:none)
#   The default strategy is :scalable (n_k_rungs=3, max_iter=2, skip_solid,
#   hard_window, soft_em, cheap_correct, graph_mode=:doublestrand).
#
# This harness reconstructs the iterative arm from its parts so each phase is
# timed SEPARATELY, using the EXACT knobs _assemble_with_iterative_corrector
# feeds the corrector, so (a)+(b) reproduces the sweep's iterative `runtime_s`.
# It also times the shared scaffolding — read simulation, the naive arm, and the
# scoring (assembly_metrics + genome_fraction, and optional QUAST) — so the
# report can say what fraction of end-to-end wall time is corrector vs harness.
#
# METHOD (non-invasive — corrector/assembler behavior UNCHANGED): each phase is
# a plain @elapsed around the REAL library call the sweep makes. No library
# function is modified or monkey-patched. A warmup pass on a tiny fixture forces
# compilation so the reported seconds are ~pure runtime.
#
# RUN (from repo root or this worktree):
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. benchmarking/e2e_phase_profile.jl
#
# Optional env (defaults chosen so 1/2/3 kb finish in a few minutes total):
#   MYCELIA_E2E_GENOME_LENS  comma list of genome lengths, bp (default 1000,2000,3000)
#   MYCELIA_E2E_READLEN      read length, bp                  (default 150)
#   MYCELIA_E2E_ERR          per-base error rate              (default 0.05)
#   MYCELIA_E2E_COVERAGE     fold coverage                    (default 30, matches sweep)
#   MYCELIA_E2E_K            assembly k-mer size              (default 21, matches sweep)
#   MYCELIA_E2E_SEED         RNG seed                         (default 42)
#   MYCELIA_E2E_RUN_QUAST    truthy -> also time run_quast    (default false)

import Mycelia
import FASTX
import BioSequences
import Random
import Printf
import Dates

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
_parse_int_list(s, default) = isempty(strip(s)) ? default :
    Int[parse(Int, strip(x)) for x in split(s, ",") if !isempty(strip(x))]
_truthy(s) = lowercase(strip(s)) in ("1", "true", "yes", "on")

const GENOME_LENS = _parse_int_list(get(ENV, "MYCELIA_E2E_GENOME_LENS", ""), [1000, 2000, 3000])
const READLEN     = parse(Int, get(ENV, "MYCELIA_E2E_READLEN", "150"))
const ERR_RATE    = parse(Float64, get(ENV, "MYCELIA_E2E_ERR", "0.05"))
const COVERAGE    = parse(Float64, get(ENV, "MYCELIA_E2E_COVERAGE", "30"))
const K           = parse(Int, get(ENV, "MYCELIA_E2E_K", "21"))
const SEED        = parse(Int, get(ENV, "MYCELIA_E2E_SEED", "42"))
const RUN_QUAST   = _truthy(get(ENV, "MYCELIA_E2E_RUN_QUAST", "false"))
# readlen<=500 -> short/illumina, else long/nanopore (mirrors the sweep's regime_for_readlen)
const TECH        = READLEN <= 500 ? :illumina : :nanopore
const REGIME      = READLEN <= 500 ? "short-low-error" : "long-high-error"

# ---------------------------------------------------------------------------
# Read simulation — byte-for-byte the sweep's simulate_regime_reads recipe so the
# profiled workload matches the validated :scalable regime.
# ---------------------------------------------------------------------------
function simulate_regime_reads(refseq::BioSequences.LongDNA{4}, readlen::Int,
        coverage::Real, error_rate::Float64, tech::Symbol, rng::Random.AbstractRNG)
    glen = length(refseq)
    effective_readlen = min(readlen, glen)
    n_reads = max(1, ceil(Int, coverage * glen / effective_readlen))
    records = Vector{FASTX.FASTQ.Record}()
    for i in 1:n_reads
        start = glen == effective_readlen ? 1 : rand(rng, 1:(glen - effective_readlen + 1))
        frag = refseq[start:(start + effective_readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs_seq, quals = Mycelia.observe(frag; error_rate = error_rate, tech = tech)
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records
end

function write_fastq_records(reads)
    dir = mktempdir()
    path = joinpath(dir, "e2e_input.fastq")
    open(FASTX.FASTQ.Writer, path) do w
        for r in reads
            write(w, r)
        end
    end
    return path
end

function write_contigs_fasta(result, path::String)
    open(path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# Phase timing for ONE (genome_len) cell.
#
# Phases (each a plain @elapsed around the exact library call the sweep makes):
#   d1 sim_reads   — read simulation (shared: both arms use the same reads)
#   d2 naive_arm   — assemble_genome(reads; corrector=:none) (the "naive" arm +
#                    the reference cost of a naive assembly; the iterative arm's
#                    re-assembly is a naive assembly on corrected reads)
#   a  correction  — mycelia_iterative_assemble(...) with the :scalable knobs
#                    _assemble_with_iterative_corrector uses (the CORRECTOR)
#   b  reassembly  — assemble_genome(corrected_reads; corrector=:none)
#   c  scoring     — write contigs FASTA + assembly_metrics + genome_fraction
#                    (+ optional run_quast) — the harness scoring scaffolding
#
# iterative-arm wall (what the sweep reports as the iterative runtime_s) = a + b.
# ---------------------------------------------------------------------------
function profile_cell(glen::Int)
    rng = Random.MersenneTwister(SEED)
    ref_rec = Mycelia.random_fasta_record(moltype = :DNA, seed = SEED, L = glen)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, ref_rec)
    ref_path = joinpath(mktempdir(), "ref_$(glen).fasta")
    open(FASTX.FASTA.Writer, ref_path) do w
        write(w, ref_rec)
    end

    # k must not exceed the shortest read (the sweep clamps + warns identically).
    eff_readlen = min(READLEN, glen)
    k = K > eff_readlen ? eff_readlen : K

    # --- d1: read simulation ---
    local reads
    t_sim = @elapsed (reads = simulate_regime_reads(refseq, READLEN, COVERAGE, ERR_RATE, TECH, rng))
    n_reads = length(reads)

    # --- d2: naive arm ---
    local naive_asm
    t_naive = @elapsed (naive_asm = Mycelia.Rhizomorph.assemble_genome(
        reads; k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand,
        corrector = :none, verbose = false))

    # --- a: CORRECTION (mycelia_iterative_assemble with the :scalable knobs the
    #        corrector=:iterative route feeds it). max_k = max(k, 13). ---
    knobs = Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable)
    corrector_graph_mode = knobs.graph_mode === nothing ?
        Mycelia.Rhizomorph._graph_mode_symbol(Mycelia.Rhizomorph.DoubleStrand) : knobs.graph_mode
    max_k = max(k, 13)
    input_fastq = write_fastq_records(reads)
    corr_out = mktempdir()
    local result_dict
    t_correct = @elapsed (result_dict = Mycelia.mycelia_iterative_assemble(
        input_fastq;
        max_k = max_k,
        skip_solid = knobs.skip_solid,
        graph_mode = corrector_graph_mode,
        n_k_rungs = knobs.n_k_rungs,
        max_iterations_per_k = knobs.max_iterations_per_k,
        hard_window = knobs.hard_window,
        soft_em = knobs.soft_em,
        cheap_correct = knobs.cheap_correct,
        beam_width = knobs.beam_width,
        verbose = false,
        enable_checkpointing = false,
        output_dir = corr_out))

    corrected_fastq = get(result_dict[:metadata], :final_fastq_file, nothing)
    corrected_reads = open(FASTX.FASTQ.Reader, corrected_fastq) do reader
        collect(reader)
    end
    n_corrected = length(corrected_reads)

    # --- b: RE-ASSEMBLY of corrected reads (auto graph mode, corrector=:none) ---
    reassembly_k = k
    local iter_asm
    t_reassemble = @elapsed (iter_asm = Mycelia.Rhizomorph.assemble_genome(
        corrected_reads; k = reassembly_k, corrector = :none, verbose = false))

    # --- c: SCORING (write FASTA + assembly_metrics + genome_fraction [+ QUAST]) ---
    score_dir = mktempdir()
    iter_contigs = joinpath(score_dir, "iterative_contigs.fasta")
    naive_contigs = joinpath(score_dir, "naive_contigs.fasta")
    local n50, largest, gfrac
    t_score = @elapsed begin
        write_contigs_fasta(iter_asm, iter_contigs)
        write_contigs_fasta(naive_asm, naive_contigs)
        m = Mycelia.assembly_metrics(iter_contigs)
        n50 = m === nothing ? 0 : m.n50
        largest = m === nothing ? 0 : m.largest_contig
        total_len = sum(length.(iter_asm.contigs); init = 0)
        gfrac = glen == 0 ? 0.0 : round(total_len / glen * 100; digits = 1)
    end

    t_quast = 0.0
    if RUN_QUAST
        try
            t_quast = @elapsed Mycelia.run_quast(
                String[iter_contigs, naive_contigs];
                outdir = joinpath(score_dir, "quast"),
                reference = ref_path,
                min_contig = max(50, glen ÷ 10))
        catch e
            @warn "run_quast unavailable / failed — reporting t_quast=0" exception = e
        end
    end

    return (
        glen = glen, k = k, n_reads = n_reads, n_corrected = n_corrected,
        n_contigs_iter = length(iter_asm.contigs), n50 = n50, largest = largest, gfrac = gfrac,
        t_sim = t_sim, t_naive = t_naive, t_correct = t_correct,
        t_reassemble = t_reassemble, t_score = t_score, t_quast = t_quast,
    )
end

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
function print_cell_breakdown(io, r)
    iter_arm = r.t_correct + r.t_reassemble                 # == sweep iterative runtime_s
    e2e = r.t_sim + r.t_naive + iter_arm + r.t_score + r.t_quast

    println(io, "\n" * "-"^78)
    println(io, Printf.@sprintf(
        "CELL glen=%d bp  k=%d  reads=%d  corrected=%d  iter_contigs=%d  N50=%d  frac=%.1f%%",
        r.glen, r.k, r.n_reads, r.n_corrected, r.n_contigs_iter, r.n50, r.gfrac))
    println(io, "-"^78)

    # (1) ITERATIVE-ARM breakdown: % of the number the sweep reports as runtime_s.
    println(io, "Iterative arm (== sweep-reported iterative runtime_s = correction + re-assembly):")
    println(io, Printf.@sprintf("  %-28s %10s %10s", "phase", "seconds", "% arm"))
    for (name, secs) in (("(a) correction  [CORRECTOR]", r.t_correct),
                         ("(b) re-assembly  [scaffold]", r.t_reassemble))
        pct = iter_arm == 0 ? 0.0 : 100 * secs / iter_arm
        println(io, Printf.@sprintf("  %-28s %10.2f %9.1f%%", name, secs, pct))
    end
    println(io, Printf.@sprintf("  %-28s %10.2f %9.1f%%", "iterative arm TOTAL", iter_arm, 100.0))

    # (2) FULL end-to-end (per-cell) including the shared/scoring scaffolding.
    println(io, "\nFull end-to-end per cell (adds the shared read-sim + naive arm + scoring):")
    println(io, Printf.@sprintf("  %-28s %10s %10s", "phase", "seconds", "% e2e"))
    rows = (("(d1) read simulation  [shared]", r.t_sim),
            ("(d2) naive arm        [shared]", r.t_naive),
            ("(a)  correction  [CORRECTOR]  ", r.t_correct),
            ("(b)  re-assembly [scaffold]   ", r.t_reassemble),
            ("(c)  scoring     [scaffold]   ", r.t_score))
    for (name, secs) in rows
        pct = e2e == 0 ? 0.0 : 100 * secs / e2e
        println(io, Printf.@sprintf("  %-28s %10.2f %9.1f%%", name, secs, pct))
    end
    if RUN_QUAST
        pct = e2e == 0 ? 0.0 : 100 * r.t_quast / e2e
        println(io, Printf.@sprintf("  %-28s %10.2f %9.1f%%", "(c') QUAST      [scaffold]    ", r.t_quast, pct))
    end
    println(io, Printf.@sprintf("  %-28s %10.2f %9.1f%%", "end-to-end TOTAL", e2e, 100.0))
end

function print_scaling(io, results)
    println(io, "\n" * "="^78)
    println(io, "SCALING across genome sizes (seconds): correction vs re-assembly vs scoring")
    println(io, "="^78)
    println(io, Printf.@sprintf("  %-8s %10s %10s %10s %10s %10s %8s",
        "glen", "sim", "naive", "correct", "reassemb", "score", "corr%arm"))
    for r in results
        iter_arm = r.t_correct + r.t_reassemble
        corr_pct = iter_arm == 0 ? 0.0 : 100 * r.t_correct / iter_arm
        println(io, Printf.@sprintf("  %-8d %10.2f %10.2f %10.2f %10.2f %10.2f %7.1f%%",
            r.glen, r.t_sim, r.t_naive, r.t_correct, r.t_reassemble, r.t_score, corr_pct))
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    io = stdout
    println(io, "="^78)
    println(io, "e2e phase profile — Rhizomorph correction-validation sweep iterative arm (td-9q84)")
    println(io, "  started : $(Dates.now())")
    println(io, "  regime  : $(REGIME) ($(TECH)), readlen=$(READLEN), err=$(ERR_RATE), cov=$(COVERAGE)x, k=$(K)")
    println(io, "  sizes   : $(GENOME_LENS) bp")
    println(io, "  strategy: :scalable (default) — knobs = $(Mycelia.Rhizomorph._corrector_strategy_knobs(:scalable))")
    println(io, "  QUAST   : $(RUN_QUAST ? "enabled" : "disabled (set MYCELIA_E2E_RUN_QUAST=true)")")
    println(io, "="^78)

    # --- Warmup on a tiny genome so timed cells are ~pure runtime, not compile. ---
    println(io, "\n[warmup] compiling all phases on a 300 bp fixture ...")
    _w = profile_cell(300)
    println(io, "[warmup] done (warmup correction=$(Printf.@sprintf("%.2f", _w.t_correct))s, discarded)")

    results = Vector{Any}()
    for glen in GENOME_LENS
        println(io, "\n[cell] profiling glen=$(glen) bp ...")
        r = profile_cell(glen)
        push!(results, r)
        print_cell_breakdown(io, r)
    end

    print_scaling(io, results)

    # --- Verdict ---
    println(io, "\n" * "="^78)
    println(io, "DOMINANT-PHASE READOUT")
    println(io, "="^78)
    biggest = results[end]
    iter_arm = biggest.t_correct + biggest.t_reassemble
    e2e = biggest.t_sim + biggest.t_naive + iter_arm + biggest.t_score + biggest.t_quast
    corr_pct_arm = iter_arm == 0 ? 0.0 : 100 * biggest.t_correct / iter_arm
    corr_pct_e2e = e2e == 0 ? 0.0 : 100 * biggest.t_correct / e2e
    println(io, Printf.@sprintf(
        "  Largest cell (glen=%d): iterative arm = %.1fs; correction = %.1fs (%.1f%% of arm, %.1f%% of e2e).",
        biggest.glen, iter_arm, biggest.t_correct, corr_pct_arm, corr_pct_e2e))
    println(io, Printf.@sprintf(
        "  Re-assembly = %.1fs (%.1f%% of arm); scoring = %.1fs; read-sim+naive = %.1fs.",
        biggest.t_reassemble, 100 - corr_pct_arm, biggest.t_score, biggest.t_sim + biggest.t_naive))
    if corr_pct_arm >= 60
        println(io, "  => CORRECTION dominates the iterative arm. The corrector is still the")
        println(io, "     right optimization target (not the benchmark's re-assembly/scoring).")
    elseif biggest.t_reassemble > biggest.t_correct
        println(io, "  => RE-ASSEMBLY dominates. The corrector is already fast; the next win is")
        println(io, "     the naive re-assembly pass on corrected reads, not the corrector.")
    else
        println(io, "  => Mixed. Correction and scaffolding are comparable — see the table above.")
    end
    println(io, "  finished: $(Dates.now())")
    println(io, "="^78)
end

main()
