# Hybrid-OLC Engineering Validation — Rhizomorph on REAL genomes
# ==============================================================
# INTERIM ENGINEERING VALIDATION ONLY — this is not the manuscript H5 result.
# Compare the native naive and :scalable paths against the short-read hybrid-OLC
# route on the same reads, asking whether Stage-1 correction followed by MEGAHIT
# preserves reference accuracy while providing an established layout path.
#
# Data provenance (documented, honest):
#   * REFERENCE  — REAL genome downloaded from NCBI RefSeq by accession.
#   * READS      — realistic Illumina reads SIMULATED FROM THE REAL REFERENCE via
#                  ART (HS25 = HiSeq 2500 empirical error/quality model). SRA
#                  acquisition (prefetch/fasterq-dump) is unavailable in this
#                  environment (sra-tools not installed), so we use ART's
#                  platform-realistic error model applied to the REAL reference.
#                  This exercises the real genome's repeat/GC structure — the
#                  property a random SMOKE genome cannot test — while keeping read
#                  error characteristics matched to a real Illumina platform.
#
# Arms:
#   * naive      — assemble_genome(reads; corrector=:none, k=21)
#   * scalable   — assemble_genome(reads; corrector=:iterative,
#                                  strategy=:scalable, k=21)
#   * hybrid-olc — assemble_genome(reads; corrector=:iterative,
#                                  strategy=:scalable, layout=:olc,
#                                  olc_tool=:megahit, k=21)
#
# ART simulates paired-end reads, but the current public assembly API receives
# the concatenated R1-then-R2 records as one unpaired stream for every arm. A
# later paired-end validation is tracked separately.
#
# Reference-free structural metrics (#contigs, N50, largest contig) use the same
# k+1 minimum contig length for every arm. Reference-aligned genome fraction,
# identity, mismatches, and indels come from MUMmer dnadiff. Genome fraction is
# alignment based (not total_length/reference_length).
# Runtime is recorded for diagnostics, not as a fair performance comparison:
# the interim hybrid arm is pinned to one MEGAHIT thread (see its config below).
#
# Read-only w.r.t. src/ — this script only consumes the public Mycelia API.
#
# Usage:
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#       benchmarking/real_data_corrector_validation.jl
#
# Env knobs:
#   MYCELIA_RDV_TARGETS   comma list of names: phix174,lambda  (default: phix174,lambda)
#   MYCELIA_RDV_COVERAGE  fold coverage for ART                (default: 50)
#   MYCELIA_RDV_K         k-mer size                           (default: 21)
#   MYCELIA_RDV_SEED      ART rndSeed                          (default: 42)
#   MYCELIA_RDV_SMOKE     "true" -> phix174 only, coverage 30  (default: false)

import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import Dates
import JSON
import SHA

# --- Config -----------------------------------------------------------------

function _truthy(value::AbstractString)::Bool
    return lowercase(strip(value)) in ("1", "true", "yes", "on")
end

const SMOKE = _truthy(get(ENV, "MYCELIA_RDV_SMOKE", "false"))

# name => (versioned_accession, expected_size_bp)
const REGISTRY = Dict(
    "phix174" => ("NC_001422.1", 5386),   # phiX174 — classic Illumina control
    "lambda" => ("NC_001416.1", 48502),  # Enterobacteria phage lambda
)

const DEFAULT_TARGETS = ["phix174", "lambda"]
const TARGETS = let
    raw = get(ENV, "MYCELIA_RDV_TARGETS", SMOKE ? "phix174" : "phix174,lambda")
    [String(strip(x)) for x in split(raw, ",") if !isempty(strip(x))]
end

const COVERAGE = parse(Int, get(ENV, "MYCELIA_RDV_COVERAGE", SMOKE ? "30" : "50"))
const K = parse(Int, get(ENV, "MYCELIA_RDV_K", "21"))
const SEED = parse(Int, get(ENV, "MYCELIA_RDV_SEED", "42"))
const ARMS = (:naive, :scalable, :hybrid_olc)
const VALIDATION_SCOPE = "interim_engineering_validation"
const ART_PROFILE = "HS25"
const READ_LENGTH = 150
const FRAGMENT_MEAN = 300
const FRAGMENT_SD = 10
const SIMULATED_LAYOUT = "paired_end"
const ASSEMBLY_INPUT_LAYOUT = "single_end_unpaired"
const MIN_CONTIG_LENGTH = K + 1
const PROMOTION_ELIGIBLE = !SMOKE && TARGETS == DEFAULT_TARGETS &&
                           COVERAGE == 50 && K == 21 && SEED == 42
const MIN_HYBRID_IDENTITY = 99.9
const MIN_HYBRID_GENOME_FRACTION = 95.0
const SCALABLE_N50_RETENTION = 0.99
const SCALABLE_GENOME_FRACTION_TOLERANCE = 1.0
const IDENTITY_TOLERANCE = 0.02

# --- Helpers ----------------------------------------------------------------

function arm_name(arm::Symbol)::String
    arm in ARMS || throw(ArgumentError(
        "unknown validation arm :$(arm); expected one of $(ARMS)"))
    return arm == :hybrid_olc ? "hybrid-olc" : String(arm)
end

"""Load a FASTQ file into a `Vector{FASTX.FASTQ.Record}`."""
function load_fastq(path::String)::Vector{FASTX.FASTQ.Record}
    return open(FASTX.FASTQ.Reader, path) do reader
        return collect(reader)
    end
end

function sha256_file(path::String)::String
    return bytes2hex(SHA.sha256(read(path)))
end

function sha256_files(paths::Tuple{Vararg{String}})::String
    contents = UInt8[]
    for path in paths
        append!(contents, read(path))
    end
    return bytes2hex(SHA.sha256(contents))
end

function conda_package_record(environment::String, package::String)::String
    output = read(
        `$(Mycelia.CONDA_RUNNER) list -n $(environment) $(package) --json`,
        String,
    )
    records = filter(
        record -> record["name"] == package,
        JSON.parse(output),
    )
    length(records) == 1 || error(
        "expected one $(package) record in conda environment $(environment), " *
        "found $(length(records))",
    )
    record = only(records)
    return "$(record["name"])-$(record["version"])-" *
           "$(record["build_string"])@$(record["platform"])"
end

"""
Parse a MUMmer 3.23 `dnadiff` .report into the reference-based metrics we need.
The shared `Mycelia.parse_dnadiff_report` does not decode MUMmer's
`count(pct%)` single-token format for the `[Bases]` block, so we parse it here
(read-only w.r.t. src/). Returns (aligned_ref_bases, genome_fraction_pct,
avg_identity_pct, total_snps, total_indels); absent fields remain non-finite or
`nothing` so artifact validation fails closed. The single `AlignedBases` entry
is the M-to-M reference coverage; the first `AvgIdentity` is the 1-to-1 identity.

Report layout (columns are REF then QRY):
    AlignedBases   48487(99.97%)   48502(100.00%)
    AvgIdentity    99.99           99.99            # first occ = 1-to-1
    TotalSNPs      12              12
    TotalIndels    3               3
"""
function parse_dnadiff_local(report::String)::NamedTuple
    aln_ref = 0
    gf = NaN
    ident = NaN
    snps = nothing
    indels = nothing
    seen_aligned = false
    seen_ident = false

    function count_pct(
            tok::AbstractString)::Union{Nothing, Tuple{Int, Float64}}
        m = match(r"^(\d+)\(([\d.]+)%\)", tok)
        if m === nothing
            return nothing
        end
        return (parse(Int, m.captures[1]), parse(Float64, m.captures[2]))
    end

    for line in eachline(report)
        f = split(strip(line))
        isempty(f) && continue
        key = f[1]
        if key == "AlignedBases" && length(f) >= 2
            if seen_aligned
                aln_ref = 0
                gf = NaN
                continue
            end
            seen_aligned = true
            count_and_pct = count_pct(f[2])
            if count_and_pct !== nothing
                c, p = count_and_pct
                aln_ref = c
                gf = p
            end
        elseif key == "AvgIdentity" && !seen_ident
            seen_ident = true
            if length(f) >= 2
                v = tryparse(Float64, f[2])
                v !== nothing && (ident = v)
            end
        elseif key == "TotalSNPs" && length(f) >= 2
            v = tryparse(Int, f[2])
            v !== nothing && (snps = v)
        elseif key == "TotalIndels" && length(f) >= 2
            v = tryparse(Int, f[2])
            v !== nothing && (indels = v)
        end
    end
    return (aligned_ref_bases = aln_ref, genome_fraction = gf,
        avg_identity = ident, total_snps = snps, total_indels = indels)
end

"""Write assembly contigs to a FASTA file."""
function write_contigs(
        contigs::AbstractVector{<:AbstractString}, path::String)::String
    open(path, "w") do io
        for (i, contig) in enumerate(contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end
    return path
end

"""
Run one assembly arm and evaluate it against the reference with MUMmer `dnadiff`.
Returns a NamedTuple of metrics.

Reference-based metrics come from dnadiff (QUAST's old bioconda build does not
solve on osx-arm64; dnadiff/MUMmer installs cleanly and yields analogous
reference-aligned quantities with the definitions below):
  * genome fraction (%) = M-to-M AlignedBases % of the reference
  * avg identity (%)    = first 1-to-1 AvgIdentity
  * mismatches/100 kbp  = TotalSNPs / M-to-M aligned reference bases * 1e5
  * indels/100 kbp      = TotalIndels / M-to-M aligned reference bases * 1e5

The normalized error rates are harness diagnostics, not claimed to be identical
to QUAST's rates.
"""
function run_arm(
        name::String,
        reads::AbstractVector{<:FASTX.FASTQ.Record},
        arm::Symbol,
        reference::String,
        outdir::String)::NamedTuple
    tag = arm_name(arm)
    mkpath(outdir)
    contigs_path = joinpath(outdir, "$(name)_$(tag)_contigs.fasta")

    t0 = time()
    local result
    try
        if arm == :naive
            result = Mycelia.Rhizomorph.assemble_genome(
                reads;
                k = K,
                graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                corrector = :none,
                verbose = false,
            )
        elseif arm == :scalable
            result = Mycelia.Rhizomorph.assemble_genome(
                reads;
                k = K,
                graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                corrector = :iterative,
                strategy = :scalable,
                sequencing_tech = :illumina,
                verbose = false,
            )
        else
            result = Mycelia.Rhizomorph.assemble_genome(
                reads;
                k = K,
                graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                corrector = :iterative,
                strategy = :scalable,
                sequencing_tech = :illumina,
                layout = :olc,
                olc_tool = :megahit,
                # Pin one thread for a stable local engineering artifact. On
                # these corrected validation reads, the Bioconda osx-arm64
                # MEGAHIT 1.2.9 k-mer counter crashed at three threads; one thread
                # succeeds. Manuscript H5 performance tuning is explicitly out of
                # scope, so runtime is diagnostic rather than a fair comparison.
                olc_options = (;
                    k_list = string(K),
                    min_contig_len = MIN_CONTIG_LENGTH,
                    threads = 1,
                ),
                verbose = false,
            )
        end
    catch e
        e isa InterruptException && rethrow()
        @warn "Assembly arm failed" name tag exception = (e, catch_backtrace())
        return (name = name, arm = tag,
            validation_scope = VALIDATION_SCOPE,
            ok = false, runtime_s = round(time() - t0; digits = 2),
            raw_n_contigs = 0, n_contigs = 0, total_length = 0,
            largest_contig = 0, n50 = 0,
            genome_fraction = NaN, avg_identity = NaN,
            mismatches_per_100kbp = NaN, indels_per_100kbp = NaN)
    end
    runtime = time() - t0
    raw_n_contigs = length(result.contigs)
    contigs = filter(contig -> length(contig) >= MIN_CONTIG_LENGTH, result.contigs)
    write_contigs(contigs, contigs_path)
    n_contigs = length(contigs)
    total_length = sum(length.(contigs); init = 0)

    # Reference-free structural metrics
    n50 = 0
    largest = 0
    am = Mycelia.assembly_metrics(contigs_path)
    if am !== nothing
        n50 = am.n50
        largest = am.largest_contig
    end

    # Reference-based metrics via MUMmer dnadiff
    gf = NaN
    ident = NaN
    mm = NaN
    indels_100k = NaN
    if n_contigs > 0 && total_length > 0
        dd_out = joinpath(outdir, "$(name)_$(tag)_dnadiff")
        try
            paths = Mycelia.run_dnadiff(reference = reference, query = contigs_path,
                outdir = dd_out, prefix = "dnadiff", force = true)
            p = parse_dnadiff_local(paths.report)
            gf = p.genome_fraction
            ident = p.avg_identity
            if p.aligned_ref_bases > 0 &&
                    p.total_snps !== nothing &&
                    p.total_indels !== nothing
                mm = round(p.total_snps / p.aligned_ref_bases * 1e5; digits = 1)
                indels_100k = round(p.total_indels / p.aligned_ref_bases * 1e5; digits = 1)
            end
        catch e
            e isa InterruptException && rethrow()
            @warn "dnadiff failed" name tag exception = (e, catch_backtrace())
        end
    end

    metrics_complete = all(isfinite, (gf, ident, mm, indels_100k))
    return (name = name, arm = tag,
        validation_scope = VALIDATION_SCOPE,
        ok = metrics_complete, runtime_s = round(runtime; digits = 2),
        raw_n_contigs = raw_n_contigs, n_contigs = n_contigs,
        total_length = total_length,
        largest_contig = largest, n50 = n50,
        genome_fraction = gf, avg_identity = ident,
        mismatches_per_100kbp = mm, indels_per_100kbp = indels_100k)
end

"""Apply the interim hybrid-vs-naive and hybrid-vs-scalable engineering gate."""
function validate_comparative_gate(
        df::DataFrames.DataFrame,
        targets::AbstractVector{<:AbstractString})::Nothing
    for target in targets
        target_rows = df[df.name .== String(target), :]
        naive = only(DataFrames.eachrow(target_rows[target_rows.arm .== "naive", :]))
        scalable = only(DataFrames.eachrow(
            target_rows[target_rows.arm .== "scalable", :],
        ))
        hybrid = only(DataFrames.eachrow(
            target_rows[target_rows.arm .== "hybrid-olc", :],
        ))

        hybrid.n50 > naive.n50 || error(
            "hybrid-olc N50 must exceed naive for $(target)",
        )
        hybrid.n_contigs < naive.n_contigs || error(
            "hybrid-olc contig count must be below naive for $(target)",
        )
        hybrid.avg_identity >= MIN_HYBRID_IDENTITY || error(
            "hybrid-olc identity is below $(MIN_HYBRID_IDENTITY)% for $(target)",
        )
        hybrid.avg_identity >= naive.avg_identity - IDENTITY_TOLERANCE || error(
            "hybrid-olc identity materially regressed from naive for $(target)",
        )
        hybrid.genome_fraction >= MIN_HYBRID_GENOME_FRACTION || error(
            "hybrid-olc genome fraction is below $(MIN_HYBRID_GENOME_FRACTION)% " *
            "for $(target)",
        )
        hybrid.mismatches_per_100kbp <= naive.mismatches_per_100kbp || error(
            "hybrid-olc mismatch rate exceeds naive for $(target)",
        )
        hybrid.indels_per_100kbp <= naive.indels_per_100kbp || error(
            "hybrid-olc indel rate exceeds naive for $(target)",
        )

        hybrid.n50 >= SCALABLE_N50_RETENTION * scalable.n50 || error(
            "hybrid-olc retains less than $(100 * SCALABLE_N50_RETENTION)% of " *
            "scalable N50 for $(target)",
        )
        hybrid.avg_identity >= scalable.avg_identity - IDENTITY_TOLERANCE || error(
            "hybrid-olc identity materially regressed from scalable for $(target)",
        )
        hybrid.genome_fraction >= scalable.genome_fraction -
                                  SCALABLE_GENOME_FRACTION_TOLERANCE || error(
            "hybrid-olc genome fraction regressed by more than " *
            "$(SCALABLE_GENOME_FRACTION_TOLERANCE) point from scalable for $(target)",
        )
        hybrid.mismatches_per_100kbp <= scalable.mismatches_per_100kbp || error(
            "hybrid-olc mismatch rate exceeds scalable for $(target)",
        )
        hybrid.indels_per_100kbp <= scalable.indels_per_100kbp || error(
            "hybrid-olc indel rate exceeds scalable for $(target)",
        )
    end
    return nothing
end

"""Fail closed unless every requested target and arm has complete evidence."""
function validate_results(
        df::DataFrames.DataFrame,
        targets::AbstractVector{<:AbstractString},
        arms::Tuple{Vararg{Symbol}})::Nothing
    length(unique(targets)) == length(targets) ||
        error("validation targets must be unique: $(targets)")
    length(unique(arms)) == length(arms) ||
        error("validation arms must be unique: $(arms)")

    arm_names = arm_name.(arms)
    expected_pairs = Set((String(target), arm_name)
    for target in targets for arm_name in arm_names)
    actual_pairs = Set(zip(String.(df.name), String.(df.arm)))
    actual_pairs == expected_pairs || error(
        "incomplete target-arm matrix: expected $(sort!(collect(expected_pairs))), " *
        "got $(sort!(collect(actual_pairs)))")
    DataFrames.nrow(df) == length(expected_pairs) || error(
        "duplicate validation rows: expected $(length(expected_pairs)), " *
        "got $(DataFrames.nrow(df))")
    all(df.validation_scope .== VALIDATION_SCOPE) ||
        error("validation_scope must be $(VALIDATION_SCOPE) for every row")
    all(df.ok) || error("one or more validation arms failed or lack dnadiff metrics")

    finite_columns = [
        :runtime_s,
        :genome_fraction,
        :avg_identity,
        :mismatches_per_100kbp,
        :indels_per_100kbp,
    ]
    for column in finite_columns
        all(isfinite, df[!, column]) || error("non-finite metric in column $(column)")
    end

    all(df.runtime_s .>= 0.0) || error("runtime_s must be non-negative")
    all(df.raw_n_contigs .>= df.n_contigs) ||
        error("raw_n_contigs must be at least n_contigs")
    all(df.n_contigs .> 0) || error("every arm must produce at least one contig")
    all(df.total_length .> 0) || error("every arm must produce positive total length")
    all(df.largest_contig .> 0) || error("every arm must produce a positive largest contig")
    all(df.n50 .> 0) || error("every arm must produce a positive N50")
    all(df.n50 .<= df.largest_contig) || error("n50 must not exceed largest_contig")
    all(df.largest_contig .<= df.total_length) ||
        error("largest_contig must not exceed total_length")
    all((0.0 .< df.genome_fraction) .& (df.genome_fraction .<= 100.0)) ||
        error("genome_fraction must be in (0, 100]")
    all((0.0 .< df.avg_identity) .& (df.avg_identity .<= 100.0)) ||
        error("avg_identity must be in (0, 100]")
    all(df.mismatches_per_100kbp .>= 0.0) ||
        error("mismatches_per_100kbp must be non-negative")
    all(df.indels_per_100kbp .>= 0.0) ||
        error("indels_per_100kbp must be non-negative")

    all(!isempty, df.art_profile) || error("ART profile must not be empty")
    all(df.requested_coverage .> 0) || error("requested coverage must be positive")
    all(df.read_length .> 0) || error("read length must be positive")
    all(df.fragment_mean .> 0) || error("fragment mean must be positive")
    all(df.fragment_sd .>= 0) || error("fragment standard deviation must be non-negative")
    all(df.art_seed .>= 0) || error("ART seed must be non-negative")
    all(df.k .> 0) || error("k must be positive")
    all(df.min_contig_length .== df.k .+ 1) ||
        error("all arms must use the common k+1 minimum contig length")
    all(!isempty, df.simulated_layout) || error("simulated layout must not be empty")
    all(!isempty, df.assembly_input_layout) ||
        error("assembly input layout must not be empty")
    all(!isempty, df.input_read_order) || error("input read order must not be empty")
    length(unique(df.promotion_eligible)) == 1 ||
        error("promotion eligibility must be identical across the run")
    all(df.read_count .> 0) || error("read_count must be positive")
    all(df.read_bases .> 0) || error("read_bases must be positive")
    all(value -> occursin(r"^[0-9a-f]{64}$", value), df.reference_sha256) ||
        error("invalid reference SHA-256")
    all(value -> occursin(r"^[0-9a-f]{64}$", value), df.input_reads_sha256) ||
        error("invalid input-read SHA-256")

    provenance_columns = [
        :reference_accession,
        :reference_length,
        :reference_sha256,
        :art_profile,
        :requested_coverage,
        :read_length,
        :fragment_mean,
        :fragment_sd,
        :art_seed,
        :k,
        :min_contig_length,
        :simulated_layout,
        :assembly_input_layout,
        :input_read_order,
        :read_count,
        :read_bases,
        :input_reads_sha256,
    ]
    for target in targets
        target_rows = df[df.name .== String(target), :]
        for column in provenance_columns
            length(unique(target_rows[!, column])) == 1 || error(
                "$(column) must be identical across arms for $(target)",
            )
        end
        accession, expected_length = REGISTRY[String(target)]
        all(target_rows.reference_accession .== accession) ||
            error("unexpected reference accession for $(target)")
        all(target_rows.reference_length .== expected_length) ||
            error("unexpected reference length for $(target)")
    end

    for column in [:art_build, :megahit_build, :mummer_build, :gfatools_build]
        length(unique(df[!, column])) == 1 ||
            error("$(column) must identify one tool build for the run")
        all(!isempty, df[!, column]) || error("$(column) must not be empty")
    end
    all(!isempty, df.run_platform) || error("run_platform must not be empty")
    all(!isempty, df.julia_version) || error("julia_version must not be empty")

    validate_comparative_gate(df, targets)
    return nothing
end

# --- Main -------------------------------------------------------------------

function main()::Nothing
    println("=== Hybrid-OLC Engineering Validation ===")
    println("Claim scope  : INTERIM ENGINEERING VALIDATION ONLY (not manuscript H5)")
    println("Start        : $(Dates.now())")
    println("Targets      : $(join(TARGETS, ", "))")
    println("Coverage     : $(COVERAGE)x  |  k=$(K)  |  ART seed=$(SEED)")
    println("Read model   : ART $(ART_PROFILE), $(READ_LENGTH) bp paired-end")
    println("Assembly I/O : R1 then R2 as one unpaired stream for every arm")
    println("Contig floor : $(MIN_CONTIG_LENGTH) bp (k+1) for every arm")
    println("Promotable   : $(PROMOTION_ELIGIBLE)")

    workdir = mktempdir(prefix = "rdv_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    println("Workdir      : $(workdir)")

    rows = NamedTuple[]

    for name in TARGETS
        haskey(REGISTRY, name) || (@warn "Unknown target, skipping" name; continue)
        accession, expected_length = REGISTRY[name]
        println("\n" * "="^70)
        println("Target: $(name)  ($(accession), $(expected_length) bp)")
        println("="^70)

        target_dir = joinpath(workdir, name)
        mkpath(target_dir)

        # 1. Download the version-pinned real reference from NCBI RefSeq.
        println("[1/3] Downloading reference $(accession) ...")
        reference = Mycelia.download_genome_by_accession(
            accession = accession,
            outdir = target_dir,
            compressed = false,
        )
        if !isfile(reference) || filesize(reference) == 0
            @warn "Reference download failed; skipping target" name accession
            continue
        end
        ref_records = open(FASTX.FASTA.Reader, reference) do reader
            return collect(reader)
        end
        length(ref_records) == 1 || error(
            "expected one reference record for $(accession), found $(length(ref_records))",
        )
        ref_len = length(FASTX.sequence(BioSequences.LongDNA{4}, only(ref_records)))
        ref_len == expected_length || error(
            "reference length mismatch for $(accession): expected $(expected_length), " *
            "got $(ref_len)",
        )
        reference_sha256 = sha256_file(reference)
        println("      Reference: $(reference)  ($(ref_len) bp, $(reference_sha256))")

        # 2. Simulate paired Illumina reads, then preserve one R1-then-R2 record
        # stream for the current single-input assembly API used by every arm.
        println("[2/3] Simulating Illumina reads " *
                "(ART $(ART_PROFILE), $(COVERAGE)x, $(READ_LENGTH) bp PE) ...")
        outbase = joinpath(target_dir, "$(name)_illumina")
        art = Mycelia.simulate_illumina_reads(
            fasta = reference,
            coverage = COVERAGE,
            outbase = outbase,
            read_length = READ_LENGTH,
            mflen = FRAGMENT_MEAN,
            sdev = FRAGMENT_SD,
            seqSys = ART_PROFILE,
            paired = true,
            errfree = false,
            rndSeed = SEED,
            quiet = true,
        )
        r1 = replace(art.forward_reads, ".gz" => "")
        r2 = replace(art.reverse_reads, ".gz" => "")
        for (gz, fq) in ((art.forward_reads, r1), (art.reverse_reads, r2))
            (isfile(gz) && filesize(gz) > 0) || error("ART did not produce $(gz)")
            isfile(fq) || run(`gunzip -k -f $(gz)`)
            (isfile(fq) && filesize(fq) > 0) || error("gunzip did not produce $(fq)")
        end
        reads = vcat(load_fastq(r1), load_fastq(r2))
        read_bases = sum(length(FASTX.sequence(String, record)) for record in reads)
        input_reads_sha256 = sha256_files((r1, r2))
        println("      Reads: $(length(reads)) ($(read_bases) bp, " *
                "~$(round(read_bases / ref_len; digits = 1))x effective, " *
                "$(input_reads_sha256))")

        provenance = (
            reference_accession = accession,
            reference_length = ref_len,
            reference_sha256 = reference_sha256,
            art_profile = ART_PROFILE,
            requested_coverage = COVERAGE,
            read_length = READ_LENGTH,
            fragment_mean = FRAGMENT_MEAN,
            fragment_sd = FRAGMENT_SD,
            art_seed = SEED,
            k = K,
            min_contig_length = MIN_CONTIG_LENGTH,
            simulated_layout = SIMULATED_LAYOUT,
            assembly_input_layout = ASSEMBLY_INPUT_LAYOUT,
            input_read_order = "R1_then_R2",
            read_count = length(reads),
            read_bases = read_bases,
            input_reads_sha256 = input_reads_sha256,
            promotion_eligible = PROMOTION_ELIGIBLE,
        )

        # 3. Every arm receives this exact `reads` vector. The common k+1
        # post-filter is also passed explicitly to MEGAHIT in `run_arm`.
        println("[3/3] Assembling (naive vs scalable vs hybrid-olc, k=$(K)) + " *
                "dnadiff ...")
        target_rows = NamedTuple[]
        for arm in ARMS
            tag = arm_name(arm)
            row = merge(
                run_arm(
                    name,
                    reads,
                    arm,
                    reference,
                    joinpath(target_dir, tag),
                ),
                provenance,
            )
            push!(rows, row)
            push!(target_rows, row)
        end

        for row in target_rows
            println("   [$(row.arm)] contigs=$(row.n_contigs)/$(row.raw_n_contigs) " *
                    "total=$(row.total_length) N50=$(row.n50) " *
                    "largest=$(row.largest_contig) GF=$(row.genome_fraction)% " *
                    "ident=$(row.avg_identity)% " *
                    "mm/100kb=$(row.mismatches_per_100kbp) " *
                    "indels/100kb=$(row.indels_per_100kbp) t=$(row.runtime_s)s")
        end
    end

    isempty(rows) && error("no validation results produced")

    df = DataFrames.DataFrame(rows)
    n_rows = DataFrames.nrow(df)
    df.art_build = fill(conda_package_record("art", "art"), n_rows)
    df.megahit_build = fill(conda_package_record("megahit", "megahit"), n_rows)
    df.mummer_build = fill(conda_package_record("mummer", "mummer"), n_rows)
    df.gfatools_build = fill(conda_package_record("gfatools", "gfatools"), n_rows)
    df.run_platform = fill(Sys.MACHINE, n_rows)
    df.julia_version = fill(string(VERSION), n_rows)

    validate_results(df, TARGETS, ARMS)
    stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS_sss")
    csv_path = joinpath(
        results_dir,
        "real_data_corrector_validation_$(stamp)_p$(getpid()).csv",
    )
    CSV.write(csv_path, df)

    println("\n" * "="^70)
    println("SUMMARY — INTERIM ENGINEERING VALIDATION ONLY")
    println("Real references, same ART-simulated reads; not manuscript H5")
    println("="^70)
    show(df; allrows = true, allcols = true)
    println("\n\nCandidate CSV: $(csv_path)")
    println("Registration : UNREGISTERED until the manifest pointer is updated explicitly")
    println("Promotable   : $(PROMOTION_ELIGIBLE)")
    println("Done         : $(Dates.now())")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
