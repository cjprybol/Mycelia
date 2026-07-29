# OPT4 FROZEN-READ SKIP — ASSEMBLY-LEVEL (dnadiff/ANI) ACCURACY CHECK
# =============================================================================
# td-jbjd, Rule-of-5 pass 5 follow-on. Design doc:
# docs/design/2026-07-26-opt4-frozen-read-skip-design.md ("Accuracy sign-off
# methodology") lists TWO checks: (1) per-base accuracy against injected-error
# ground truth (done by benchmarking/opt4_frozen_read_accuracy_signoff.jl,
# results in benchmarking/results/opt4_frozen_read_skip_signoff.md) and (2) an
# assembly-level dnadiff/ANI check (re-assemble the corrected reads and diff
# the resulting contigs against the reference). This script does (2), on the
# SAME fixture (Lambda phage NC_001416, seed 42, :scalable production knobs)
# pass 5 used.
#
# THREE ARMS, each producing an assembly diffed against the reference:
#   1. raw       -- the uncorrected simulated reads (sanity baseline: confirms
#                    the assembler works and correction helps at all)
#   2. exact      -- corrected reads from skip_frozen_reads=false
#   3. freeze_across -- corrected reads from skip_frozen_reads=true,
#                    threshold=2, across=true (the non-default, non-maximal
#                    scope pass 5 found gives ~15-17% decode-volume savings for
#                    a -0.02% relative recall cost at 21x per-base)
#
# CRITICAL mechanics note (why this actually tests something): each arm's
# CORRECTED reads are assembled WITHOUT re-running the iterative corrector --
# assemble_genome(...; corrector=:none, use_quality_scores=false) routes
# straight to the plain k-mer graph (_assemble_kmer_graph), never touching
# mycelia_iterative_assemble. If the corrected reads were re-corrected
# identically here, the exact-vs-freeze_across difference already present in
# the input reads would be washed out before it ever reached the assembler.
# Same k / graph_mode / dedup for every arm (only the input reads differ), so
# any contig-level delta between exact and freeze_across is attributable to
# the frozen-read-skip mechanism, not to a downstream assembly wobble.
#
# Reuses correct_reads_with_freeze / acquire_reference / simulate_substitution_
# reads from opt4_frozen_read_accuracy_signoff.jl (safe to include: guarded by
# abspath(PROGRAM_FILE) == @__FILE__, so including only DEFINES functions).
#
# Run directly (network required to download the reference once; mummer
# bioconda env must be present -- verified pre-run via `Mycelia.add_bioconda_
# env("mummer")` which is a no-op if already installed):
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/opt4_frozen_read_dnadiff_signoff.jl
#
# Config override env vars (same names/defaults as the per-base script):
#   MYCELIA_OPT4_ACCESSION   default NC_001416 (Lambda phage, 48502 bp)
#   MYCELIA_OPT4_COVERAGE    default 8.0 (set to 21.0 for VERDICT-tier parity
#                             with the per-base sign-off's Tier B)
#   MYCELIA_OPT4_ERR         default 0.02 (substitution rate)
#   MYCELIA_OPT4_READLEN     default 150
#   MYCELIA_OPT4_K           default 21
#   MYCELIA_OPT4_BATCH_SIZE  default 500
#   MYCELIA_OPT4_SEED        default 42

import Random
import FASTX
import BioSequences
import Dates
import CSV
import DataFrames

include(joinpath(@__DIR__, "opt4_frozen_read_accuracy_signoff.jl"))

"""
Assemble `records` (a Vector of FASTX.FASTQ.Record) into contigs WITHOUT
running the iterative corrector: `corrector=:none` (the AssemblyConfig
default) routes straight to the plain k-mer / qualmer graph, never calling
`mycelia_iterative_assemble`. `use_quality_scores=false` forces the
deterministic k-mer graph path (`_assemble_kmer_graph`) instead of the
qualmer path, so no placeholder quality injection or quality-weighted
traversal is in play -- the assembly is a pure function of the input
sequences. Returns the `AssemblyResult`.
"""
function assemble_no_correction(
        records::Vector{FASTX.FASTQ.Record}; k::Int, graph_mode::Symbol)
    gm = graph_mode == :doublestrand ? Mycelia.Rhizomorph.DoubleStrand :
         graph_mode == :canonical ? Mycelia.Rhizomorph.Canonical :
         Mycelia.Rhizomorph.SingleStrand
    return Mycelia.Rhizomorph.assemble_genome(
        records; corrector = :none, k = k, graph_mode = gm,
        use_quality_scores = false, verbose = false)
end

"""
Write `result.contigs` (Vector{String}) to a FASTA file at `path` and return
`path`.
"""
function write_contigs_fasta(result, path::String; label_prefix::String)
    recs = FASTX.FASTA.Record[]
    for (i, seq) in enumerate(result.contigs)
        push!(recs, FASTX.FASTA.Record("$(label_prefix)_contig_$(i)", seq))
    end
    open(FASTX.FASTA.Writer, path) do w
        for r in recs
            write(w, r)
        end
    end
    return path
end

function run_dnadiff_check()
    accession = get(ENV, "MYCELIA_OPT4_ACCESSION", "NC_001416")
    coverage = parse(Float64, get(ENV, "MYCELIA_OPT4_COVERAGE", "21.0"))
    err = parse(Float64, get(ENV, "MYCELIA_OPT4_ERR", "0.02"))
    readlen = parse(Int, get(ENV, "MYCELIA_OPT4_READLEN", "150"))
    k = parse(Int, get(ENV, "MYCELIA_OPT4_K", "21"))
    batch_size = parse(Int, get(ENV, "MYCELIA_OPT4_BATCH_SIZE", "500"))
    seed = parse(Int, get(ENV, "MYCELIA_OPT4_SEED", "42"))
    assigned_q = 20

    println("=== OPT4 frozen-read-skip ASSEMBLY-LEVEL (dnadiff/ANI) check ===")
    println("Start: $(Dates.now())")
    println(
        "Config: accession=$accession coverage=$(coverage)x err=$err readlen=$readlen k=$k batch_size=$batch_size seed=$seed")

    workdir = mktempdir(prefix = "opt4_dnadiff_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)

    # Verify mummer env up front so a missing dependency fails fast, before
    # spending the correction/assembly compute budget.
    try
        Mycelia.add_bioconda_env("mummer")
    catch e
        println("BLOCKER: mummer bioconda env unavailable: $e")
        rethrow(e)
    end

    refseq, ref_path,
    ref_label = acquire_reference(
        smoke = false, accession = accession, smoke_len = 0, seed = seed,
        workdir = workdir)
    glen = length(refseq)
    println("Reference: $ref_label ($glen bp)")

    rng = Random.MersenneTwister(seed)
    records, truth_by_id,
    observed_by_id,
    injected_total,
    sampled_bases = simulate_substitution_reads(
        refseq, readlen, coverage, err, rng; assigned_q = assigned_q)
    eff_cov = round(sampled_bases / glen; digits = 2)
    println(
        "Simulated $(length(records)) reads, effective coverage $(eff_cov)x, injected $injected_total substitutions")

    # --- Correction step: exact + freeze_across only (this check's two arms) ---
    correction_arms = [
        (label = "exact", skip_frozen = false, threshold = 2, across = false),
        (label = "freeze_across", skip_frozen = true, threshold = 2, across = true)
    ]

    corrected_records = Dict{String, Vector{FASTX.FASTQ.Record}}()
    corrected_records["raw"] = records  # uncorrected baseline

    for arm in correction_arms
        println("\n--- Correcting: $(arm.label) ---")
        t0 = time()
        corrected_by_id,
        result = correct_reads_with_freeze(
            records, k; skip_frozen = arm.skip_frozen, threshold = arm.threshold,
            across = arm.across, batch_size = batch_size,
            output_dir = joinpath(workdir, "correct_$(arm.label)"))
        runtime = time() - t0
        frozen_skipped = result[:metadata][:corrector_errors][:frozen_reads_skipped]
        println(
            "  correction runtime=$(round(runtime; digits=1))s frozen_reads_skipped=$frozen_skipped")
        # Rebuild FASTQ records in the ORIGINAL read order/id set from the
        # corrected sequence dict, reusing each original record's quality
        # string length-for-length (corrected sequences preserve read length
        # per the corrector's length contract). This keeps record count/order
        # identical to `records` so the assembler sees the same read set,
        # differing only in sequence content.
        recs = Vector{FASTX.FASTQ.Record}(undef, length(records))
        for (i, r) in enumerate(records)
            id = FASTX.identifier(r)
            corrected_seq = corrected_by_id[id]
            qual = FASTX.quality(r)
            recs[i] = FASTX.FASTQ.Record(id, corrected_seq, qual)
        end
        corrected_records[arm.label] = recs
    end

    # --- Assembly step: assemble each arm's reads WITHOUT re-correcting ---
    dnadiff_rows = DataFrames.DataFrame(
        arm = String[], n_contigs = Int[], total_contig_bases = Int[],
        largest_contig = Int[], avg_identity = Float64[],
        aligned_pct_ref = Float64[], aligned_pct_query = Float64[],
        total_snps = Float64[], total_indels = Float64[])

    assembly_results = Dict{String, Any}()
    for arm_label in ("raw", "exact", "freeze_across")
        println("\n--- Assembling (corrector=:none): $(arm_label) ---")
        t0 = time()
        result = assemble_no_correction(
            corrected_records[arm_label]; k = k, graph_mode = :doublestrand)
        runtime = time() - t0
        n_contigs = length(result.contigs)
        total_bases = sum(length.(result.contigs); init = 0)
        largest = isempty(result.contigs) ? 0 : maximum(length.(result.contigs))
        println(
            "  assembly runtime=$(round(runtime; digits=1))s n_contigs=$n_contigs total_bases=$total_bases largest=$largest")
        assembly_results[arm_label] = result

        query_fasta = joinpath(workdir, "$(arm_label)_contigs.fasta")
        write_contigs_fasta(result, query_fasta; label_prefix = arm_label)

        if isempty(result.contigs)
            println("  WARNING: zero contigs produced for arm=$arm_label; skipping dnadiff")
            push!(dnadiff_rows,
                (arm = arm_label, n_contigs = 0, total_contig_bases = 0,
                    largest_contig = 0, avg_identity = NaN, aligned_pct_ref = NaN,
                    aligned_pct_query = NaN, total_snps = 0, total_indels = 0))
            continue
        end

        dnadiff_outdir = joinpath(workdir, "dnadiff_$(arm_label)")
        dnadiff_paths = Mycelia.run_dnadiff(;
            reference = ref_path, query = query_fasta, outdir = dnadiff_outdir,
            prefix = "dnadiff_$(arm_label)")
        report = Mycelia.parse_dnadiff_report(dnadiff_paths.report)
        # Copy the raw report out (mktempdir workdir is auto-removed on exit) so
        # numbers are independently verifiable after the run.
        cp(dnadiff_paths.report,
            joinpath(results_dir, "opt4_dnadiff_$(arm_label).report"); force = true)
        # parse_dnadiff_report returns (; summary, distance, distance_aligned,
        # raw_sections). `summary` is a NamedTuple keyed by lowercased dnadiff
        # field names; extract by fuzzy match to be robust to exact naming.
        smry = Dict{Symbol, Any}(pairs(report.summary))
        for (_sec, _mets) in report.raw_sections
            for (_mk, _mv) in _mets
                get!(smry, _mk, _mv)  # summary wins on key collision
            end
        end
        println("  dnadiff summary: ",
            join(["$(k)=$(smry[k])" for k in sort(collect(keys(smry)))], " "))
        function _find(subs::Vector{String})
            for k in keys(smry)
                lk = lowercase(string(k))
                if any(occursin(s, lk) for s in subs)
                    return smry[k]
                end
            end
            return nothing
        end
        function _num(x)
            x === nothing && return NaN
            x isa Number && return Float64(x)
            s = string(x)
            m = match(r"\(([0-9.]+)%\)", s)  # "48000(98.9%)" -> 98.9
            m !== nothing && return something(tryparse(Float64, m.captures[1]), NaN)
            return something(tryparse(Float64, replace(s, "%" => "")), NaN)
        end
        avg_identity = _num(_find(["avgidentity", "avg_identity", "averageidentity"]))
        # parse_dnadiff_report's summary drops TotalSNPs/TotalIndels/AlignedBases
        # (it matches the bare tokens SNPs/Indels, but dnadiff emits TotalSNPs/
        # TotalIndels). Parse those fields DIRECTLY from the .report file, which
        # is committed alongside this CSV as reproducible evidence. Report line
        # format: `FieldName   <ref-col>   <query-col>` (ref col first).
        function _report_field(report_path::String, field::String)
            for line in eachline(report_path)
                parts = split(strip(line))
                if !isempty(parts) && parts[1] == field
                    return length(parts) >= 2 ? parts[2] : nothing
                end
            end
            return nothing
        end
        # AlignedBases ref col is like "48469(99.93%)" -> _num extracts 99.93.
        aligned_pct_ref = _num(_report_field(dnadiff_paths.report, "AlignedBases"))
        aligned_pct_query = _num(_find(["alignedbasesquery", "aligned_bases_query", "alignedquery"]))
        snps = _num(_report_field(dnadiff_paths.report, "TotalSNPs"))
        indels = _num(_report_field(dnadiff_paths.report, "TotalIndels"))

        println(
            "  dnadiff: avg_identity=$avg_identity aligned_pct_ref=$aligned_pct_ref " *
            "aligned_pct_query=$aligned_pct_query snps=$snps indels=$indels")

        push!(dnadiff_rows,
            (arm = arm_label, n_contigs = n_contigs, total_contig_bases = total_bases,
                largest_contig = largest, avg_identity = avg_identity,
                aligned_pct_ref = aligned_pct_ref, aligned_pct_query = aligned_pct_query,
                total_snps = snps, total_indels = indels))
    end

    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv_path = joinpath(
        results_dir, "opt4_frozen_read_skip_dnadiff_$(timestamp).csv")
    CSV.write(csv_path, dnadiff_rows)
    println("\nCSV written: $csv_path")
    println(dnadiff_rows)
    println("\n=== dnadiff check complete: $(Dates.now()) ===")
    return (
        rows = dnadiff_rows, csv_path = csv_path, glen = glen, ref_label = ref_label,
        n_reads = length(records), eff_cov = eff_cov,
        config = (; accession, coverage, err, readlen, k, batch_size, seed))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_dnadiff_check()
end
