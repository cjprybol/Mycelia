# Rhizomorph REAL-genome validation benchmark (td-wlfq).
#
# Every Rhizomorph corrector result to date used simulated reads drawn from a
# RANDOM reference. Random references have no repeats, no biological k-mer
# structure, and trivially unique k-mers, so they cannot answer the core
# UNVALIDATED question:
#
#   1. Does the iterative corrector improve assembly quality over the naive
#      (uncorrected) baseline on a REAL genome?
#   2. Is the all-solid-read skip (skip_solid=true) LOSSLESS relative to
#      correcting every read (skip_solid=false)?
#
# This benchmark answers both against a real reference (Enterobacteria phage
# lambda, NCBI NC_001416, ~48 kb) with reference-based QUAST metrics
# (NGA50, genome fraction, misassemblies, duplication ratio).
#
# Public API under test (merged in PRs #348/#350):
#   Mycelia.Rhizomorph.assemble_genome(reads;
#       k, corrector=:none|:iterative, skip_solid=false, graph_mode=...)
#       -> AssemblyResult (with .contigs)
#
# Three arms, all fed the SAME simulated read set:
#   naive           corrector=:none      graph_mode=DoubleStrand
#   iterative       corrector=:iterative skip_solid=false graph_mode=Canonical
#   iterative+skip  corrector=:iterative skip_solid=true  graph_mode=Canonical
#
# graph_mode rationale (both are hard constraints of the current code, not a
# design choice of this benchmark):
#   * The naive graph path emits INVALID contigs under Canonical (undirected
#     canonical traversal is not orientation-aware -- see the @warn in
#     src/rhizomorph/assembly.jl::_assemble_kmer_graph), so the naive arm uses
#     the valid DoubleStrand reconstruction.
#   * skip_solid is only honored for graph_mode=:canonical (the corrector
#     disables the skip and warns on any other mode -- see
#     src/iterative-assembly.jl), so BOTH iterative arms use Canonical. The
#     corrector path does not use the invalid Canonical graph reconstruction;
#     it emits error-corrected reads read from its final FASTQ.
#
# Interpretation note: the v0 iterative corrector returns error-corrected READS
# as its `contigs` (one corrected read per contig), NOT graph-assembled
# unitigs. The naive arm returns graph-assembled contigs. QUAST contextualizes
# both against the reference. The lossless question (arm 2 vs arm 3) is the
# cleanest and most defensible comparison because the corrector is identical
# and ONLY skip_solid differs.
#
# Scaling: toy-runnable defaults finish locally in minutes on a real subregion
# of lambda; env vars scale to the full genome + higher coverage for a full run
# (Lovelace). QUAST is wrapped via Conda.jl and is skipped gracefully (with a
# clear message) when its conda env cannot be provisioned locally.
#
#   MYCELIA_RGV_ACCESSION   NCBI accession                 (default NC_001416)
#   MYCELIA_RGV_REGION_LEN  bp of real genome to use; 0=full   (default 4000)
#   MYCELIA_RGV_COVERAGE    fold read coverage                 (default 15)
#   MYCELIA_RGV_READLEN     simulated read length (bp)         (default 300)
#   MYCELIA_RGV_ERR         per-base substitution error rate   (default 0.02)
#   MYCELIA_RGV_K           k-mer size                         (default 21)
#   MYCELIA_RGV_SEED        RNG seed                           (default 42)
#   MYCELIA_RGV_MIN_CONTIG  QUAST --min-contig                 (default 100)
#   MYCELIA_RGV_QUAST       run QUAST: auto|on|off             (default auto)
#   MYCELIA_RGV_SKIPFRAC    measure skip fraction: true|false  (default true)
#   MYCELIA_RGV_SYNTHETIC   allow synthetic ref if offline     (default false)
#
# Usage:
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. benchmarking/rhizomorph_real_genome_validation.jl
#   MYCELIA_RGV_REGION_LEN=0 MYCELIA_RGV_COVERAGE=30 ... (full run, Lovelace)

import Mycelia
import Random
import Statistics
import Printf
import Dates
import FASTX

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
getenv(k, d::Int) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
getenv(k, d::Float64) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
getenv(k, d::String) = get(ENV, k, d)

const ACCESSION  = getenv("MYCELIA_RGV_ACCESSION", "NC_001416")
const REGION_LEN = getenv("MYCELIA_RGV_REGION_LEN", 4000)   # 0 => full genome
const COVERAGE   = getenv("MYCELIA_RGV_COVERAGE", 15)
const READLEN    = getenv("MYCELIA_RGV_READLEN", 300)
const ERR        = getenv("MYCELIA_RGV_ERR", 0.02)
const K          = getenv("MYCELIA_RGV_K", 21)
const SEED       = getenv("MYCELIA_RGV_SEED", 42)
const MIN_CONTIG = getenv("MYCELIA_RGV_MIN_CONTIG", 100)
const QUAST_MODE = getenv("MYCELIA_RGV_QUAST", "auto")      # auto|on|off
const DO_SKIPFRAC = getenv("MYCELIA_RGV_SKIPFRAC", "true") == "true"
const ALLOW_SYNTHETIC = getenv("MYCELIA_RGV_SYNTHETIC", "false") == "true"
const BASES = ('A', 'C', 'G', 'T')

# ---------------------------------------------------------------------------
# 1. Reference acquisition (real genome from NCBI; graceful offline handling)
# ---------------------------------------------------------------------------
"""
Acquire the reference genome sequence for `accession`. Uses Mycelia's
`download_genome_by_accession` helper (NCBI efetch via `get_sequence`).
Returns the uppercase DNA sequence String, or `nothing` if acquisition fails
(e.g. no network) and no synthetic fallback is permitted.
"""
function acquire_reference(accession::String, outdir::String)
    try
        fna = Mycelia.download_genome_by_accession(;
            accession = accession, outdir = outdir, compressed = false)
        if isfile(fna) && filesize(fna) > 0
            rec = first(collect(FASTX.FASTA.Reader(open(fna))))
            seq = uppercase(FASTX.FASTA.sequence(String, rec))
            if occursin(r"^[ACGTN]+$", seq) && length(seq) > 1000
                return seq
            end
        end
        @warn "Downloaded reference for $(accession) was empty or malformed."
    catch e
        @warn "Reference acquisition failed (offline?)" exception = e
    end
    if ALLOW_SYNTHETIC
        @warn "Falling back to a SYNTHETIC random reference (MYCELIA_RGV_SYNTHETIC=true). " *
              "This is a smoke test ONLY -- results do NOT constitute real-genome validation."
        rng = Random.MersenneTwister(SEED)
        return join(rand(rng, BASES, REGION_LEN == 0 ? 48000 : REGION_LEN))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# 2. Read simulation (pure Julia: window sampling + per-base substitution).
#    Conda-free and deterministic, so the toy path runs anywhere. Badread-based
#    simulators (Mycelia.simulate_nanopore_reads / simulate_pacbio_reads) are
#    the higher-fidelity option for a full Lovelace run but require a conda env.
# ---------------------------------------------------------------------------
"""
Simulate `coverage`-fold reads of length `readlen` from `ref` with per-base
substitution rate `err`. Returns a Vector of FASTX.FASTQ.Record (both strands
sampled so DoubleStrand/Canonical assembly is exercised).
"""
function simulate_reads(ref::String; coverage::Int, readlen::Int, err::Float64, seed::Int)
    rng = Random.MersenneTwister(seed)
    reflen = length(ref)
    rl = min(readlen, reflen)
    n_reads = max(1, round(Int, coverage * reflen / rl))
    records = Vector{FASTX.FASTQ.Record}(undef, n_reads)
    revcomp(c) = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'N'
    for i in 1:n_reads
        start = rand(rng, 1:(reflen - rl + 1))
        seq = collect(ref[start:(start + rl - 1)])
        # sample ~half the reads from the reverse strand
        if rand(rng, Bool)
            seq = [revcomp(c) for c in reverse(seq)]
        end
        quals = fill('I', rl)  # Phred+33 'I' == Q40
        for j in 1:rl
            if rand(rng) < err
                seq[j] = rand(rng, filter(!=(seq[j]), BASES))
                quals[j] = Char(clamp(round(Int, 22 + 8 * randn(rng)), 2, 40) + 33)
            end
        end
        records[i] = FASTX.FASTQ.Record("r$(i)", join(seq), join(quals))
    end
    return records
end

# ---------------------------------------------------------------------------
# Helpers: write sequences to FASTA, parse QUAST report.tsv
# ---------------------------------------------------------------------------
function write_fasta(path::String, contigs::AbstractVector{<:AbstractString}; prefix = "contig")
    open(path, "w") do io
        for (i, c) in enumerate(contigs)
            isempty(c) && continue
            println(io, ">$(prefix)_$(i)")
            println(io, c)
        end
    end
    return path
end

"""
Parse a QUAST `report.tsv` into Dict{metric_name => Dict{assembly_name => value}}.
QUAST's TSV has assemblies as columns; the first column holds metric names.
"""
function parse_quast_report(tsv::String)
    lines = readlines(tsv)
    isempty(lines) && return Dict{String, Dict{String, String}}()
    header = split(lines[1], '\t')
    assemblies = header[2:end]
    metrics = Dict{String, Dict{String, String}}()
    for line in lines[2:end]
        cols = split(line, '\t')
        length(cols) < 2 && continue
        metric = cols[1]
        d = Dict{String, String}()
        for (j, a) in enumerate(assemblies)
            d[a] = (j + 1) <= length(cols) ? cols[j + 1] : "-"
        end
        metrics[metric] = d
    end
    return metrics
end

# ---------------------------------------------------------------------------
# 3. Assembly arms via the public assemble_genome API
# ---------------------------------------------------------------------------
"""
Run one assembly arm; return (contigs::Vector{String}, runtime_seconds::Float64).
"""
function run_arm(reads; kwargs...)
    t0 = time()
    result = Mycelia.Rhizomorph.assemble_genome(reads; kwargs...)
    elapsed = time() - t0
    return (collect(String, result.contigs), elapsed)
end

# ---------------------------------------------------------------------------
# Skip-fraction measurement. assemble_genome hardcodes verbose=false, so the
# corrector's "Stage 0 skipped ... (X%)" line is never emitted through the
# public API. To surface the metric we make ONE supplementary direct call to
# mycelia_iterative_assemble(verbose=true) and parse the reported fraction.
# ---------------------------------------------------------------------------
function measure_skip_fraction(reads_fastq::String, k::Int, outdir::String)
    buf = IOBuffer()
    try
        redirect_stdout(buf) do
            Mycelia.mycelia_iterative_assemble(reads_fastq;
                max_k = max(k, 13),
                skip_solid = true,
                graph_mode = :canonical,
                verbose = true,
                enable_checkpointing = false,
                output_dir = joinpath(outdir, "skipfrac_probe"))
        end
    catch e
        @warn "Skip-fraction probe failed" exception = e
        return nothing
    end
    s = String(take!(buf))
    fracs = Float64[]
    for m in eachmatch(r"Stage 0 skipped[^()]*\(([\d.]+)%\)", s)
        push!(fracs, parse(Float64, m.captures[1]))
    end
    return isempty(fracs) ? nothing : maximum(fracs)
end

# ---------------------------------------------------------------------------
# 4/5. QUAST + reporting
# ---------------------------------------------------------------------------
function try_run_quast(assembly_fastas, ref_path, outdir)
    if QUAST_MODE == "off"
        println("\nQUAST disabled (MYCELIA_RGV_QUAST=off).")
        return nothing
    end
    try
        quast_dir = Mycelia.run_quast(assembly_fastas;
            outdir = joinpath(outdir, "quast"),
            reference = ref_path,
            min_contig = MIN_CONTIG)
        tsv = joinpath(quast_dir, "report.tsv")
        if isfile(tsv)
            return parse_quast_report(tsv)
        end
        println("QUAST produced no report.tsv.")
        return nothing
    catch e
        msg = "QUAST unavailable locally (conda env could not be provisioned, " *
              "or execution failed). The reference-based metrics (NGA50, genome " *
              "fraction, misassemblies, duplication ratio) require a full run on " *
              "Lovelace where the QUAST conda env is available."
        if QUAST_MODE == "on"
            @error msg exception = e
        else
            println("\n" * msg)
            println("(reason: $(e))")
        end
        return nothing
    end
end

function fmt(x)
    return x === nothing ? "n/a" : string(x)
end

function main()
    println("=" ^ 78)
    println("Rhizomorph REAL-genome validation benchmark (td-wlfq)")
    println("=" ^ 78)
    println("timestamp   : $(Dates.now())")
    println("julia threads: $(Threads.nthreads())")

    workdir = mktempdir()

    # -- 1. reference ------------------------------------------------------
    println("\n[1/5] Acquiring reference $(ACCESSION) ...")
    full_seq = acquire_reference(ACCESSION, workdir)
    if full_seq === nothing
        println("\nSKIP: could not acquire the reference genome and synthetic " *
                "fallback is disabled (set MYCELIA_RGV_SYNTHETIC=true for an " *
                "offline smoke test). Nothing to validate; exiting cleanly.")
        return
    end
    ref_seq = REGION_LEN == 0 ? full_seq : full_seq[1:min(REGION_LEN, length(full_seq))]
    ref_path = write_fasta(joinpath(workdir, "reference.fasta"), [ref_seq]; prefix = "ref")
    println("  reference length : $(length(full_seq)) bp (using $(length(ref_seq)) bp)")

    # -- 2. reads ----------------------------------------------------------
    println("\n[2/5] Simulating reads (coverage=$(COVERAGE)x, readlen=$(READLEN), err=$(ERR)) ...")
    reads = simulate_reads(ref_seq; coverage = COVERAGE, readlen = READLEN, err = ERR, seed = SEED)
    println("  simulated reads  : $(length(reads))")
    reads_fastq = joinpath(workdir, "reads.fastq")
    open(reads_fastq, "w") do io
        w = FASTX.FASTQ.Writer(io)
        for r in reads
            write(w, r)
        end
        close(w)
    end

    # -- 3. three assembly arms via assemble_genome ------------------------
    println("\n[3/5] Assembling three arms via Rhizomorph.assemble_genome ...")
    println("  arm 1/3 naive (corrector=:none, DoubleStrand) ...")
    naive_contigs, naive_t = run_arm(reads;
        k = K, corrector = :none, graph_mode = Mycelia.Rhizomorph.DoubleStrand)

    println("  arm 2/3 iterative (skip_solid=false, Canonical) ...")
    iter_contigs, iter_t = run_arm(reads;
        k = K, corrector = :iterative, skip_solid = false,
        graph_mode = Mycelia.Rhizomorph.Canonical)

    println("  arm 3/3 iterative+skip (skip_solid=true, Canonical) ...")
    skip_contigs, skip_t = run_arm(reads;
        k = K, corrector = :iterative, skip_solid = true,
        graph_mode = Mycelia.Rhizomorph.Canonical)

    naive_fa = write_fasta(joinpath(workdir, "naive.fasta"), naive_contigs)
    iter_fa  = write_fasta(joinpath(workdir, "iterative.fasta"), iter_contigs)
    skip_fa  = write_fasta(joinpath(workdir, "iterative_skip.fasta"), skip_contigs)

    # -- skip fraction (supplementary direct verbose call) -----------------
    skip_fraction = nothing
    if DO_SKIPFRAC
        println("\n  measuring skip fraction (direct verbose corrector probe) ...")
        skip_fraction = measure_skip_fraction(reads_fastq, K, workdir)
    end

    # -- lossless equivalence (arm 2 vs arm 3) -----------------------------
    iter_set = Set(iter_contigs)
    skip_set = Set(skip_contigs)
    lossless_exact = iter_set == skip_set
    only_iter = length(setdiff(iter_set, skip_set))
    only_skip = length(setdiff(skip_set, iter_set))

    # -- 4. QUAST ----------------------------------------------------------
    println("\n[4/5] Running QUAST vs the reference ...")
    quast = try_run_quast([naive_fa, iter_fa, skip_fa], ref_path, workdir)

    # match QUAST column names back to arms (basename minus extension)
    arm_files = Dict("naive" => "naive", "iterative" => "iterative", "iterative+skip" => "iterative_skip")
    function q(metric, colkey)
        quast === nothing && return nothing
        haskey(quast, metric) || return nothing
        col = arm_files[colkey]
        d = quast[metric]
        # QUAST may sanitize names; try exact then prefix match
        haskey(d, col) && return d[col]
        for (k2, v) in d
            startswith(k2, col) && return v
        end
        return nothing
    end

    # -- 5. report ---------------------------------------------------------
    println("\n[5/5] Results")
    println("=" ^ 78)
    arms = ["naive", "iterative", "iterative+skip"]
    contigcounts = Dict("naive" => length(naive_contigs),
        "iterative" => length(iter_contigs), "iterative+skip" => length(skip_contigs))
    runtimes = Dict("naive" => naive_t, "iterative" => iter_t, "iterative+skip" => skip_t)

    hdr = Printf.@sprintf("%-16s %10s %10s %12s %14s %14s %12s", "arm", "runtime_s",
        "#contigs", "NGA50", "genome_frac%", "misassembl.", "dup_ratio")
    println(hdr)
    println("-" ^ length(hdr))
    for a in arms
        row = Printf.@sprintf("%-16s %10.2f %10d %12s %14s %14s %12s",
            a, runtimes[a], contigcounts[a],
            fmt(q("NGA50", a)),
            fmt(q("Genome fraction (%)", a)),
            fmt(q("# misassemblies", a)),
            fmt(q("Duplication ratio (%)", a)))
        println(row)
    end

    println("\nSkip diagnostics")
    println("-" ^ 40)
    println("  reported skip fraction : $(skip_fraction === nothing ? "n/a (probe skipped/failed)" : string(skip_fraction) * "%")")
    println("  arm2 vs arm3 exact-equal contig sets : $(lossless_exact)")
    println("  contigs only in iterative (no-skip)  : $(only_iter)")
    println("  contigs only in iterative+skip       : $(only_skip)")

    # -- verdicts ----------------------------------------------------------
    println("\nVERDICTS")
    println("=" ^ 78)

    # Q1: iterative vs naive
    if quast !== nothing
        gf_naive = tryparse(Float64, something(q("Genome fraction (%)", "naive"), ""))
        gf_iter  = tryparse(Float64, something(q("Genome fraction (%)", "iterative"), ""))
        nga_naive = tryparse(Float64, something(q("NGA50", "naive"), ""))
        nga_iter  = tryparse(Float64, something(q("NGA50", "iterative"), ""))
        print("Q1 iterative vs naive: ")
        if gf_naive !== nothing && gf_iter !== nothing
            better_gf = gf_iter >= gf_naive
            gfmsg = better_gf ? "iterative genome fraction >= naive" : "iterative genome fraction < naive"
            ngmsg = (nga_naive !== nothing && nga_iter !== nothing) ?
                (nga_iter >= nga_naive ? "; NGA50 >= naive" : "; NGA50 < naive") : ""
            verdict = better_gf ? "IMPROVES (or matches)" : "DOES NOT IMPROVE"
            println("$(verdict) -- $(gfmsg)$(ngmsg).")
        else
            println("INCONCLUSIVE -- QUAST metrics not both parseable.")
        end
    else
        println("Q1 iterative vs naive: PENDING QUAST (needs Lovelace). " *
                "Contig counts: naive=$(contigcounts["naive"]), iterative=$(contigcounts["iterative"]).")
    end

    # Q2: skip lossless
    print("Q2 skip lossless (iterative+skip vs iterative): ")
    if lossless_exact
        println("LOSSLESS -- corrected contig sets are byte-identical.")
    elseif quast !== nothing
        gf_iter = tryparse(Float64, something(q("Genome fraction (%)", "iterative"), ""))
        gf_skip = tryparse(Float64, something(q("Genome fraction (%)", "iterative+skip"), ""))
        if gf_iter !== nothing && gf_skip !== nothing && isapprox(gf_iter, gf_skip; atol = 0.5)
            println("EFFECTIVELY LOSSLESS -- contig sets differ ($(only_iter)/$(only_skip)) " *
                    "but genome fraction matches within 0.5% ($(gf_iter) vs $(gf_skip)).")
        else
            println("NOT LOSSLESS -- contig sets differ ($(only_iter)/$(only_skip)) and " *
                    "genome fraction diverges ($(fmt(gf_iter)) vs $(fmt(gf_skip))).")
        end
    else
        println("contig sets differ ($(only_iter)/$(only_skip)); QUAST needed on Lovelace " *
                "to judge whether the difference is quality-neutral.")
    end

    # -- persist a CSV summary --------------------------------------------
    resultsdir = joinpath(@__DIR__, "results")
    mkpath(resultsdir)
    stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv = joinpath(resultsdir, "rhizomorph_real_genome_validation_$(stamp).csv")
    open(csv, "w") do io
        println(io, "arm,runtime_s,num_contigs,NGA50,genome_fraction_pct,misassemblies,duplication_ratio_pct")
        for a in arms
            println(io, join([a,
                Printf.@sprintf("%.3f", runtimes[a]),
                string(contigcounts[a]),
                fmt(q("NGA50", a)),
                fmt(q("Genome fraction (%)", a)),
                fmt(q("# misassemblies", a)),
                fmt(q("Duplication ratio (%)", a))], ","))
        end
        println(io, "# accession=$(ACCESSION) region_len=$(length(ref_seq)) coverage=$(COVERAGE) " *
                "readlen=$(READLEN) err=$(ERR) k=$(K) seed=$(SEED) skip_fraction_pct=$(fmt(skip_fraction)) " *
                "lossless_exact=$(lossless_exact) quast=$(quast !== nothing)")
    end
    println("\nwrote $(csv)")
    println("=" ^ 78)
end

main()
