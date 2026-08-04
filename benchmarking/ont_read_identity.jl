# Empirical ONT read-identity measurement (td-4e19d.29)
#
# WHY THIS EXISTS
# ---------------
# `Mycelia.simulate_nanopore_reads` calls Badread with NO --error_model,
# --qscore_model or --identity flags (src/simulation.jl:1091), so the error
# process it produces is whatever the INSTALLED Badread compiles in as its
# defaults. Its docstring asserts "Oxford Nanopore R10.4.1 / nanopore2023", but
# that is documentation ABOUT Badread's defaults, not a parameter the wrapper
# pins — a Badread upgrade that changed a default would silently change every
# ONT benchmark in this repo without touching a line of Mycelia.
#
# The per-base error rate matters because it fixes the arithmetic that governs
# whether k-mer assembly on these reads can work at all:
#
#     P(error-free k-mer) = (1 - e)^k
#     error-free k-mer coverage = raw coverage x (1 - e)^k
#
# so this script MEASURES e rather than inheriting it from a model name. Two
# independent estimates are reported:
#
#   * ALIGNMENT-MEASURED identity — reads mapped back to the reference with
#     minimap2, identity recomputed from CIGAR + NM. This is the operative
#     number: it is what a downstream consumer of these reads actually sees,
#     and it absorbs everything Badread layers on top of the identity
#     distribution (adapters, junk reads, random reads, chimeras, glitches).
#   * BADREAD-DECLARED identity — the `read_identity=` field Badread stamps
#     into each FASTQ header. This is what Badread TARGETED for that read.
#
# Reporting both is the cross-check. If they agree, the identity distribution is
# behaving as documented. If the alignment-measured value is materially worse,
# the gap is the contribution of the non-identity error sources, and reads that
# fail to align at all are counted separately rather than being dropped from a
# mean (which would bias the estimate optimistic).
#
# Two identity conventions are computed, because they differ substantially on
# nanopore data and quoting the wrong one misstates e:
#
#   * BLAST identity          — matches / (matches + mismatches + inserted +
#                               deleted bases). Every indel BASE counts.
#   * gap-compressed identity — each indel RUN counts once regardless of length.
#
# BLAST identity is the one that governs k-mer survival, because a k-mer is
# destroyed by every erroneous base it covers, not by every error EVENT. It is
# therefore the estimate used for the (1-e)^k table.
#
# Usage:
#   julia --project=. benchmarking/ont_read_identity.jl
#   julia --project=. benchmarking/ont_read_identity.jl --coverage 30 --seed 42
#   julia --project=. benchmarking/ont_read_identity.jl --reads path/to/ont.fq.gz
#   julia --project=. benchmarking/ont_read_identity.jl --output-dir /scratch/ident

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import CodecZlib
import CSV
import DataFrames
import Dates
import Statistics

const ACCESSION = "NC_001416"   # Lambda, matching the sweep and the Track-A pilot
const ORGANISM = "Lambda"
const K_LADDER = [11, 15, 21, 31]
const COVERAGE_LADDER = [10, 30, 50, 100]

arg_value(flag) =
    let i = findfirst(==(flag), ARGS)
        (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
    end

const COVERAGE = something(tryparse(Int, something(arg_value("--coverage"), "30")), 30)
const SEED = something(tryparse(Int, something(arg_value("--seed"), "42")), 42)
const READS_OVERRIDE = arg_value("--reads")
const OUTPUT_DIR = something(arg_value("--output-dir"),
    joinpath(@__DIR__, "results", "ont_read_identity"))

# === CIGAR / SAM parsing ===

"""
    parse_cigar(cigar) -> (; aligned_columns, inserted_bases, deleted_bases,
                            insertion_events, deletion_events, has_ops)

Tally a SAM CIGAR string.

`aligned_columns` sums M/=/X — reference-and-query-consuming columns, i.e.
matches PLUS mismatches. S/H (clipping) and N/P are deliberately excluded: soft-
clipped bases are not part of the alignment, so counting them would depress
identity for a read that merely extends past the reference end.

`has_ops` is false for "*" or an unparseable CIGAR, which lets the caller drop
the record rather than silently score it as a perfect zero-length alignment.
"""
function parse_cigar(cigar::AbstractString)
    aligned_columns = 0
    inserted_bases = 0
    deleted_bases = 0
    insertion_events = 0
    deletion_events = 0
    has_ops = false
    number_start = 1
    for (index, character) in enumerate(cigar)
        isdigit(character) && continue
        length_field = tryparse(Int, cigar[number_start:(index - 1)])
        number_start = index + 1
        length_field === nothing && return (; aligned_columns, inserted_bases,
            deleted_bases, insertion_events, deletion_events, has_ops = false)
        has_ops = true
        if character in ('M', '=', 'X')
            aligned_columns += length_field
        elseif character == 'I'
            inserted_bases += length_field
            insertion_events += 1
        elseif character == 'D'
            deleted_bases += length_field
            deletion_events += 1
        end
        # S/H/N/P intentionally ignored — see docstring.
    end
    return (; aligned_columns, inserted_bases, deleted_bases,
        insertion_events, deletion_events, has_ops)
end

"""
    identities_from_sam(sam_path) -> (; rows, n_primary, n_unmapped, n_skipped)

Per-read identity for every PRIMARY alignment in `sam_path`.

Secondary (flag 0x100) and supplementary (flag 0x800) records are excluded so
each read contributes exactly once; without that filter a chimeric read would be
counted twice and weight the distribution toward whichever fragment aligned
better. Unmapped records (flag 0x4) are COUNTED but carry no identity — a read
that could not be aligned is a real property of the read set, and dropping such
reads silently would bias the mean identity optimistic.

Identity is reconstructed from CIGAR + the NM tag (edit distance):

    mismatches = NM - inserted_bases - deleted_bases
    blast    = (aligned_columns - mismatches) /
               (aligned_columns + inserted_bases + deleted_bases)
    gap_comp = (aligned_columns - mismatches) /
               (aligned_columns - mismatches + mismatches +
                insertion_events + deletion_events)

A record with no NM tag is skipped rather than assumed error-free.
"""
function identities_from_sam(sam_path::AbstractString)
    rows = NamedTuple[]
    n_primary = 0
    n_unmapped = 0
    n_skipped = 0
    for line in eachline(sam_path)
        (isempty(line) || startswith(line, '@')) && continue
        fields = split(line, '\t')
        length(fields) >= 11 || (n_skipped += 1; continue)
        flag = tryparse(Int, fields[2])
        flag === nothing && (n_skipped += 1; continue)
        # 0x100 secondary, 0x800 supplementary
        (flag & 0x100 != 0 || flag & 0x800 != 0) && continue
        if flag & 0x4 != 0
            n_unmapped += 1
            push!(rows,
                (read_id = String(fields[1]), mapped = false,
                    blast_identity = NaN, gap_compressed_identity = NaN,
                    aligned_columns = 0, mismatches = 0, inserted_bases = 0,
                    deleted_bases = 0))
            continue
        end

        nm = nothing
        for field in fields[12:end]
            if startswith(field, "NM:i:")
                nm = tryparse(Int, field[6:end])
                break
            end
        end
        nm === nothing && (n_skipped += 1; continue)

        cigar = parse_cigar(fields[6])
        cigar.has_ops || (n_skipped += 1; continue)

        mismatches = nm - cigar.inserted_bases - cigar.deleted_bases
        # A negative mismatch count means NM and the CIGAR disagree (NM smaller
        # than the indel bases the CIGAR declares). That is a malformed record,
        # not a 100%-identity read, so it is skipped rather than clamped.
        mismatches < 0 && (n_skipped += 1; continue)

        matches = cigar.aligned_columns - mismatches
        blast_denominator = cigar.aligned_columns + cigar.inserted_bases +
                            cigar.deleted_bases
        gap_denominator = matches + mismatches + cigar.insertion_events +
                          cigar.deletion_events
        (blast_denominator <= 0 || gap_denominator <= 0) && (n_skipped += 1; continue)

        n_primary += 1
        push!(rows,
            (
                read_id = String(fields[1]), mapped = true,
                blast_identity = matches / blast_denominator,
                gap_compressed_identity = matches / gap_denominator,
                aligned_columns = cigar.aligned_columns, mismatches = mismatches,
                inserted_bases = cigar.inserted_bases,
                deleted_bases = cigar.deleted_bases))
    end
    return (; rows, n_primary, n_unmapped, n_skipped)
end

"""
    declared_identities(fastq_gz) -> Vector{Float64}

Badread stamps `read_identity=<pct>%` into each FASTQ header. Extract it as a
FRACTION so it is directly comparable with the alignment-measured values.

This is Badread reporting its own intent, not an independent measurement — it is
here purely as a cross-check on the alignment pipeline. Reads whose header lacks
the field (none, for Badread output) are skipped.
"""
function declared_identities(fastq_gz::AbstractString)
    values = Float64[]
    open(fastq_gz) do raw
        stream = CodecZlib.GzipDecompressorStream(raw)
        line_number = 0
        for line in eachline(stream)
            line_number += 1
            line_number % 4 == 1 || continue   # header lines only
            m = match(r"read_identity=([0-9.]+)%", line)
            m === nothing && continue
            parsed = tryparse(Float64, m.captures[1])
            parsed === nothing || push!(values, parsed / 100)
        end
        close(stream)
    end
    return values
end

# === Reporting helpers ===

function quantile_summary(values::AbstractVector{<:Real})
    isempty(values) && return nothing
    sorted = sort(collect(values))
    return (
        n = length(sorted),
        min = first(sorted),
        q05 = Statistics.quantile(sorted, 0.05),
        q25 = Statistics.quantile(sorted, 0.25),
        median = Statistics.median(sorted),
        mean = Statistics.mean(sorted),
        q75 = Statistics.quantile(sorted, 0.75),
        q95 = Statistics.quantile(sorted, 0.95),
        max = last(sorted)
    )
end

fmt(x) = string(round(x; digits = 4))

# === Main ===

if abspath(PROGRAM_FILE) == @__FILE__
    println("=== ONT read-identity measurement (td-4e19d.29) ===")
    println("Start: $(Dates.now())")
    mkpath(OUTPUT_DIR)

    # --- Badread provenance: version + the defaults it will actually apply --------
    # Recorded from the INSTALLED binary, never from the wrapper's docstring.
    badread_version = try
        strip(read(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread --version`,
            String))
    catch e
        @warn "could not read badread --version" exception = e
        "unknown"
    end
    badread_help = try
        read(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n badread badread simulate --help`,
            String)
    catch e
        @warn "could not read badread simulate --help" exception = e
        ""
    end
    # Strip ANSI styling before matching: the help text is emitted with escape codes.
    plain_help = replace(badread_help, r"\e\[[0-9;]*m" => "")
    function extract_default(label)
        let m = match(Regex("$(label)[^)]*default: ([^)]+)\\)"), plain_help)
            m === nothing ? "unknown" : strip(m.captures[1])
        end
    end
    default_identity = extract_default("Sequencing identity distribution")
    default_error_model = extract_default("or a model filename")
    default_length = extract_default("Fragment length distribution")

    println("\n--- Badread provenance (installed binary) ---")
    println("  version:                $(badread_version)")
    println("  default --identity:     $(default_identity)")
    println("  default --error_model:  $(default_error_model)")
    println("  default --length:       $(default_length)")

    # --- Reference + reads --------------------------------------------------------
    refs_dir = joinpath(OUTPUT_DIR, "refs")
    mkpath(refs_dir)
    ref_path = Mycelia.download_genome_by_accession(
        accession = ACCESSION, outdir = refs_dir, compressed = false)
    println("\nReference: $(ref_path) ($(Mycelia.total_fasta_size(ref_path)) bp)")

    reads_path = if READS_OVERRIDE !== nothing
        isfile(READS_OVERRIDE) || error("--reads path does not exist: $(READS_OVERRIDE)")
        READS_OVERRIDE
    else
        reads_dir = joinpath(OUTPUT_DIR, "reads")
        mkpath(reads_dir)
        # Same call path as the sweep and the pilot: no model flags, so this
        # inherits exactly the defaults reported above.
        Mycelia.simulate_nanopore_reads(
            fasta = ref_path, quantity = "$(COVERAGE)x",
            outfile = joinpath(reads_dir, "ont_$(COVERAGE)x_seed$(SEED).fq.gz"),
            seed = SEED, quiet = true)
    end
    println("Reads:     $(reads_path)")

    # --- Align ---------------------------------------------------------------------
    sam_path = joinpath(OUTPUT_DIR, "ont_$(COVERAGE)x_seed$(SEED).sam")
    if !isfile(sam_path) || filesize(sam_path) == 0
        Mycelia.add_bioconda_env("minimap2")
        println("\nAligning with minimap2 -ax map-ont ...")
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2
                      -ax map-ont --secondary=no $(ref_path) $(reads_path)`,
            stdout = sam_path, stderr = devnull))
    else
        println("\nReusing existing SAM: $(sam_path)")
    end

    alignment = identities_from_sam(sam_path)
    declared = declared_identities(reads_path)

    measured_blast = [r.blast_identity for r in alignment.rows if r.mapped]
    measured_gap = [r.gap_compressed_identity for r in alignment.rows if r.mapped]

    blast_summary = quantile_summary(measured_blast)
    gap_summary = quantile_summary(measured_gap)
    declared_summary = quantile_summary(declared)

    n_reads_total = length(declared)
    n_mapped = length(measured_blast)

    println("\n--- Observed read identity (Lambda / ONT / $(COVERAGE)x / seed $(SEED)) ---")
    println("  reads in FASTQ:            $(n_reads_total)")
    println("  primary alignments scored: $(n_mapped)")
    println("  unmapped reads:            $(alignment.n_unmapped)")
    println("  records skipped:           $(alignment.n_skipped)")

    for (label, summary) in (("BLAST identity (alignment-measured)", blast_summary),
        ("gap-compressed identity (alignment-measured)", gap_summary),
        ("Badread-declared identity (header)", declared_summary))
        summary === nothing && continue
        println("\n  $(label):")
        println("    n=$(summary.n)  mean=$(fmt(summary.mean))  median=$(fmt(summary.median))")
        println("    min=$(fmt(summary.min))  q05=$(fmt(summary.q05))  q25=$(fmt(summary.q25))  " *
                "q75=$(fmt(summary.q75))  q95=$(fmt(summary.q95))  max=$(fmt(summary.max))")
    end

    # --- The arithmetic this measurement exists to fix ------------------------------
    #
    # e is taken from the MEAN alignment-measured BLAST identity: k-mer survival is
    # destroyed per erroneous BASE, which is what BLAST identity counts.
    if blast_summary !== nothing
        error_rate = 1 - blast_summary.mean
        println("\n--- Implied k-mer survival at measured e = $(fmt(error_rate)) ---")
        println("  (P(error-free k-mer) = (1-e)^k; error-free coverage = raw coverage x P)")
        header = rpad("k", 5) * rpad("P(clean)", 12) *
                 join([rpad("$(c)x", 10) for c in COVERAGE_LADDER])
        println("  " * header)
        ladder_rows = NamedTuple[]
        for k in K_LADDER
            p = (1 - error_rate)^k
            println("  " * rpad(k, 5) * rpad(fmt(p), 12) *
                    join([rpad(fmt(c * p), 10) for c in COVERAGE_LADDER]))
            for c in COVERAGE_LADDER
                push!(ladder_rows,
                    (k = k, p_error_free_kmer = p, coverage = c,
                        error_free_kmer_coverage = c * p))
            end
        end
        CSV.write(joinpath(OUTPUT_DIR, "kmer_survival_ladder.tsv"),
            DataFrames.DataFrame(ladder_rows); delim = '\t')
    end

    # --- Persist ---------------------------------------------------------------------
    CSV.write(joinpath(OUTPUT_DIR, "per_read_identity.tsv"),
        DataFrames.DataFrame(alignment.rows); delim = '\t')

    summary_rows = NamedTuple[]
    for (source, summary) in (("blast_identity_aligned", blast_summary),
        ("gap_compressed_identity_aligned", gap_summary),
        ("badread_declared_identity", declared_summary))
        summary === nothing && continue
        push!(summary_rows,
            (
                organism = ORGANISM, accession = ACCESSION, technology = "ont",
                coverage = COVERAGE, seed = SEED, badread_version = badread_version,
                badread_default_identity = default_identity,
                badread_default_error_model = default_error_model,
                source = source, n = summary.n, min = summary.min, q05 = summary.q05,
                q25 = summary.q25, median = summary.median, mean = summary.mean,
                q75 = summary.q75, q95 = summary.q95, max = summary.max,
                n_reads_total = n_reads_total, n_unmapped = alignment.n_unmapped))
    end
    CSV.write(joinpath(OUTPUT_DIR, "read_identity_summary.tsv"),
        DataFrames.DataFrame(summary_rows); delim = '\t')

    println("\nWrote: $(OUTPUT_DIR)")
    println("End: $(Dates.now())")
end  # PROGRAM_FILE guard
