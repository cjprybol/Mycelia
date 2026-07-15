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
# identity, mismatches, and gap-affected indel bases come from MUMmer dnadiff.
# Genome fraction is alignment based (not total_length/reference_length).
# Assembly-call wall time is recorded for diagnostics, not as a fair performance
# comparison: the interim hybrid arm is pinned to one MEGAHIT thread (see its
# config below). It excludes contig filtering, FASTA writing, and dnadiff.
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

# Pure dnadiff-derived contiguity metrics (NGA50 + misassembly proxy).
include(joinpath(@__DIR__, "dnadiff_contiguity_metrics.jl"))

# --- Config -----------------------------------------------------------------

function _truthy(value::AbstractString)::Bool
    normalized = lowercase(strip(value))
    normalized in ("1", "true", "yes", "on") && return true
    normalized in ("0", "false", "no", "off") && return false
    throw(ArgumentError(
        "invalid boolean value $(repr(value)); expected true/false, yes/no, on/off, or 1/0",
    ))
end

const SMOKE = _truthy(get(ENV, "MYCELIA_RDV_SMOKE", "false"))

# name => (versioned_accession, expected_size_bp)
const REGISTRY = Dict(
    "phix174" => ("NC_001422.1", 5386),   # phiX174 — classic Illumina control
    "lambda" => ("NC_001416.1", 48502),  # Enterobacteria phage lambda
    # E. coli K-12 MG1655 — a REPEAT-BEARING reference (7 near-identical rRNA
    # operons ~5 kb each) where short-read assemblers fragment, so arm-to-arm
    # NGA50/misassembly actually spreads. Repeat-free phiX/lambda make "parity"
    # trivial; this fixture makes it a meaningful claim. Heavy (4.64 Mb @ 50x),
    # so it is NOT in DEFAULT_TARGETS — select it explicitly on HPC via
    # MYCELIA_RDV_TARGETS=ecoli_k12. The ref_len == expected_length assertion in
    # main() fails closed if this size is wrong.
    "ecoli_k12" => ("NC_000913.3", 4641652)
)

const DEFAULT_TARGETS = ["phix174", "lambda"]
const TARGETS = let
    raw = get(ENV, "MYCELIA_RDV_TARGETS", SMOKE ? "phix174" : "phix174,lambda")
    [String(strip(x)) for x in split(raw, ",") if !isempty(strip(x))]
end

const COVERAGE = parse(Int, get(ENV, "MYCELIA_RDV_COVERAGE", SMOKE ? "30" : "50"))
const K = parse(Int, get(ENV, "MYCELIA_RDV_K", "21"))
const SEED = parse(Int, get(ENV, "MYCELIA_RDV_SEED", "42"))
# The standalone SOTA assembler arm (SPAdes on RAW reads) is an independent
# contiguity anchor for the parity comparison. It is OPT-IN (MYCELIA_RDV_SOTA)
# so the fast phix/lambda default and SMOKE path do not pay a full SPAdes run or
# require the spades conda env; enable it for the repeat-bearing parity run.
const ENABLE_SOTA = _truthy(get(ENV, "MYCELIA_RDV_SOTA", "false"))
const _DEFAULT_ARMS = ENABLE_SOTA ?
                      (:naive, :scalable, :hybrid_olc, :sota_spades) :
                      (:naive, :scalable, :hybrid_olc)
# MYCELIA_RDV_ARMS lets a run select a subset of arms, e.g. to skip the :naive
# arm, whose no-prefilter qualmer-evidence store does not scale to bacterial
# genomes (see td-ck03) — the parity comparison itself only needs :hybrid_olc
# vs :sota_spades, with :scalable as the memory-bounded reference. Valid names:
# naive, scalable, hybrid_olc, sota_spades.
const _VALID_ARMS = (:naive, :scalable, :hybrid_olc, :sota_spades)
const ARMS = let
    raw = strip(get(ENV, "MYCELIA_RDV_ARMS", ""))
    if isempty(raw)
        _DEFAULT_ARMS
    else
        selected = Tuple(Symbol(strip(x)) for x in split(raw, ",") if !isempty(strip(x)))
        for a in selected
            a in _VALID_ARMS ||
                throw(ArgumentError("unknown arm :$(a) in MYCELIA_RDV_ARMS; " *
                                    "valid: $(_VALID_ARMS)"))
        end
        selected
    end
end
const VALIDATION_SCOPE = "interim_engineering_validation"
const ART_PROFILE = "HS25"
const READ_LENGTH = 150
const FRAGMENT_MEAN = 300
const FRAGMENT_SD = 10
const SIMULATED_LAYOUT = "paired_end"
const ASSEMBLY_INPUT_LAYOUT = "single_end_unpaired"
const INPUT_READS_SHA256_SEMANTICS = "sha256_raw_R1_then_raw_R2_bytes"
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
    arm == :hybrid_olc && return "hybrid-olc"
    arm == :sota_spades && return "sota-spades"
    return String(arm)
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
        String
    )
    records = filter(
        record -> record["name"] == package,
        JSON.parse(output)
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
avg_identity_pct, total_snps, total_indel_bases); absent fields remain non-finite or
`nothing` so artifact validation fails closed. MUMmer 3.23 emits one
`AlignedBases`, two `AvgIdentity` rows (1-to-1 then M-to-M), one `TotalSNPs`,
and one `TotalIndels`; any missing, malformed, or duplicate record is rejected.
`TotalIndels` is interpreted as the number of gap-affected bases, not events.

Report layout (columns are REF then QRY):
    AlignedBases   48487(99.97%)   48502(100.00%)
    AvgIdentity    99.99           99.99            # first occ = 1-to-1
    TotalSNPs      12              12
    TotalIndels    3               3
"""
function parse_dnadiff_local(report::String)::NamedTuple
    aligned_values = Union{Nothing, Tuple{Int, Float64}}[]
    alignment_markers = String[]
    identity_modes = String[]
    identity_values = Union{Nothing, Tuple{Float64, Float64}}[]
    snp_values = Union{Nothing, Tuple{Int, Int}}[]
    indel_values = Union{Nothing, Tuple{Int, Int}}[]
    section = ""
    alignment_mode = ""

    function count_pct(
            tok::AbstractString)::Union{Nothing, Tuple{Int, Float64}}
        m = match(r"^(\d+)\((\d+(?:\.\d+)?)%\)$", tok)
        if m === nothing
            return nothing
        end
        count = parse(Int, m.captures[1])
        percentage = parse(Float64, m.captures[2])
        0.0 <= percentage <= 100.0 || return nothing
        return (count, percentage)
    end

    function float_pair(
            fields::AbstractVector{<:AbstractString},
    )::Union{Nothing, Tuple{Float64, Float64}}
        length(fields) == 3 || return nothing
        ref_value = tryparse(Float64, fields[2])
        qry_value = tryparse(Float64, fields[3])
        (ref_value === nothing || qry_value === nothing) && return nothing
        all(isfinite, (ref_value, qry_value)) || return nothing
        all(value -> 0.0 <= value <= 100.0, (ref_value, qry_value)) ||
            return nothing
        ref_value == qry_value || return nothing
        return (ref_value, qry_value)
    end

    function int_pair(
            fields::AbstractVector{<:AbstractString},
    )::Union{Nothing, Tuple{Int, Int}}
        length(fields) == 3 || return nothing
        ref_value = tryparse(Int, fields[2])
        qry_value = tryparse(Int, fields[3])
        (ref_value === nothing || qry_value === nothing) && return nothing
        (ref_value >= 0 && qry_value >= 0 && ref_value == qry_value) ||
            return nothing
        return (ref_value, qry_value)
    end

    function alignment_marker_valid(fields::AbstractVector{<:AbstractString})::Bool
        length(fields) == 3 || return false
        ref_value = tryparse(Int, fields[2])
        qry_value = tryparse(Int, fields[3])
        return ref_value !== nothing && qry_value !== nothing &&
               ref_value >= 0 && qry_value >= 0 && ref_value == qry_value
    end

    for line in eachline(report)
        stripped = strip(line)
        if startswith(stripped, "[") && endswith(stripped, "]")
            section = stripped
            alignment_mode = ""
            continue
        end
        fields = split(stripped)
        isempty(fields) && continue
        key = fields[1]
        if section == "[Alignments]" && key in ("1-to-1", "M-to-M")
            alignment_mode = alignment_marker_valid(fields) ? String(key) : ""
            push!(alignment_markers, alignment_mode)
        elseif section == "[Bases]" && key == "AlignedBases"
            if length(fields) == 3
                ref_value = count_pct(fields[2])
                qry_value = count_pct(fields[3])
                push!(
                    aligned_values,
                    ref_value === nothing || qry_value === nothing ?
                    nothing : ref_value
                )
            else
                push!(aligned_values, nothing)
            end
        elseif section == "[Alignments]" && key == "AvgIdentity"
            push!(identity_modes, alignment_mode)
            push!(identity_values, float_pair(fields))
        elseif section == "[SNPs]" && key == "TotalSNPs"
            push!(snp_values, int_pair(fields))
        elseif section == "[SNPs]" && key == "TotalIndels"
            push!(indel_values, int_pair(fields))
        end
    end

    aligned = length(aligned_values) == 1 ? only(aligned_values) : nothing
    identities_valid = length(identity_values) == 2 &&
                       alignment_markers == ["1-to-1", "M-to-M"] &&
                       identity_modes == ["1-to-1", "M-to-M"] &&
                       all(value -> value !== nothing, identity_values)
    snp_pair = length(snp_values) == 1 ? only(snp_values) : nothing
    indel_pair = length(indel_values) == 1 ? only(indel_values) : nothing

    return (
        aligned_ref_bases = aligned === nothing ? 0 : first(aligned),
        genome_fraction = aligned === nothing ? NaN : last(aligned),
        avg_identity = identities_valid ? first(identity_values[1]) : NaN,
        total_snps = snp_pair === nothing ? nothing : first(snp_pair),
        total_indel_bases = indel_pair === nothing ? nothing : first(indel_pair)
    )
end

"""Write a CSV through a same-directory temporary file and atomic rename."""
function write_csv_atomically(
        path::String,
        df::DataFrames.AbstractDataFrame
)::String
    ispath(path) && throw(ArgumentError("refusing to overwrite existing CSV: $(path)"))
    directory = dirname(path)
    mkpath(directory)
    temporary_path, io = mktemp(directory; cleanup = false)
    close(io)
    try
        CSV.write(temporary_path, df)
        mv(temporary_path, path; force = false)
    catch
        ispath(temporary_path) && rm(temporary_path; force = true)
        rethrow()
    end
    return path
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
  * indel bases/100 kbp = TotalIndels / M-to-M aligned reference bases * 1e5

The normalized error rates are harness diagnostics, not claimed to be identical
to QUAST's rates. `assembly_runtime_s` measures only the `assemble_genome` call;
it excludes filtering, FASTA writing, and reference evaluation.
"""
function run_arm(
        name::String,
        reads::AbstractVector{<:FASTX.FASTQ.Record},
        arm::Symbol,
        reference::String,
        reference_length::Integer,
        outdir::String,
        report_staging_path::String,
        report_relative_path::String
)::NamedTuple
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
                verbose = false
            )
        elseif arm == :scalable
            result = Mycelia.Rhizomorph.assemble_genome(
                reads;
                k = K,
                graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                corrector = :iterative,
                strategy = :scalable,
                sequencing_tech = :illumina,
                verbose = false
            )
        elseif arm == :hybrid_olc
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
                # scope, so assembly runtime is diagnostic rather than a fair
                # comparison.
                olc_options = (;
                    k_list = string(K),
                    min_contig_len = MIN_CONTIG_LENGTH,
                    threads = 1
                ),
                verbose = false
            )
        elseif arm == :sota_spades
            # Independent SOTA contiguity anchor: SPAdes on the RAW (uncorrected)
            # reads, using the same single-stream input as every other arm. Run
            # as a standalone assembler baseline — NOT routed through
            # assemble_genome, whose :olc layout requires the Stage-1 corrector,
            # so a raw (uncorrected) external assembly cannot go through it. This
            # answers "are we at the field's contiguity level?", distinct from
            # the hybrid arm's accurization ablation.
            raw_fastq = joinpath(outdir, "$(name)_raw.fastq")
            Mycelia.write_fastq(records = reads, filename = raw_fastq)
            spades = Mycelia.run_spades(
                fastq1 = raw_fastq,
                outdir = joinpath(outdir, "spades")
            )
            contig_strings = String[]
            open(FASTX.FASTA.Reader, spades.contigs) do reader
                for record in reader
                    push!(contig_strings, FASTX.sequence(String, record))
                end
            end
            result = (; contigs = contig_strings)
        else
            throw(ArgumentError("unsupported validation arm :$(arm)"))
        end
    catch e
        e isa InterruptException && rethrow()
        @warn "Assembly arm failed" name tag exception = (e, catch_backtrace())
        return (name = name, arm = tag,
            validation_scope = VALIDATION_SCOPE,
            ok = false, assembly_runtime_s = round(time() - t0; digits = 2),
            raw_n_contigs = 0, n_contigs = 0, total_length = 0,
            largest_contig = 0, n50 = 0,
            nga50 = 0, misassembly_proxy = -1,
            genome_fraction = NaN, avg_identity = NaN,
            aligned_reference_bases = 0, total_snps = -1,
            total_indel_bases = -1,
            mismatches_per_100kbp = NaN,
            indel_bases_per_100kbp = NaN,
            dnadiff_report_path = report_relative_path,
            dnadiff_report_sha256 = "")
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
    aligned_reference_bases = 0
    total_snps = -1
    total_indel_bases = -1
    indel_bases_100k = NaN
    dnadiff_report_sha256 = ""
    # Contiguity signals beyond reference-free N50: NGA50 (from dnadiff .1coords
    # aligned blocks vs the reference length) and a misassembly proxy (from the
    # .report [Feature Estimates] block). Both dnadiff-derived — QUAST does not
    # solve on osx-arm64. See benchmarking/dnadiff_contiguity_metrics.jl.
    nga50_value = 0
    misassembly_value = -1
    if n_contigs > 0 && total_length > 0
        dd_out = joinpath(outdir, "$(name)_$(tag)_dnadiff")
        try
            paths = Mycelia.run_dnadiff(reference = reference, query = contigs_path,
                outdir = dd_out, prefix = "dnadiff", force = true)
            p = parse_dnadiff_local(paths.report)
            cp(paths.report, report_staging_path; force = false)
            dnadiff_report_sha256 = sha256_file(report_staging_path)
            report_text = read(paths.report, String)
            misassembly_value = misassembly_proxy(parse_dnadiff_feature_estimates(report_text))
            coords_path = joinpath(dd_out, "dnadiff.1coords")
            if isfile(coords_path)
                nga50_value = nga50(
                    aligned_ref_block_lengths(read(coords_path, String)),
                    reference_length
                )
            end
            gf = p.genome_fraction
            ident = p.avg_identity
            aligned_reference_bases = p.aligned_ref_bases
            total_snps = p.total_snps === nothing ? -1 : p.total_snps
            total_indel_bases = p.total_indel_bases === nothing ?
                                -1 : p.total_indel_bases
            if p.aligned_ref_bases > 0 &&
               p.total_snps !== nothing &&
               p.total_indel_bases !== nothing
                mm = round(p.total_snps / p.aligned_ref_bases * 1e5; digits = 1)
                indel_bases_100k = round(
                    p.total_indel_bases / p.aligned_ref_bases * 1e5;
                    digits = 1
                )
            end
        catch e
            e isa InterruptException && rethrow()
            @warn "dnadiff failed" name tag exception = (e, catch_backtrace())
        end
    end

    metrics_complete = all(isfinite, (gf, ident, mm, indel_bases_100k)) &&
                       aligned_reference_bases > 0 &&
                       total_snps >= 0 &&
                       total_indel_bases >= 0 &&
                       occursin(r"^[0-9a-f]{64}$", dnadiff_report_sha256)
    return (name = name, arm = tag,
        validation_scope = VALIDATION_SCOPE,
        ok = metrics_complete, assembly_runtime_s = round(runtime; digits = 2),
        raw_n_contigs = raw_n_contigs, n_contigs = n_contigs,
        total_length = total_length,
        largest_contig = largest, n50 = n50,
        nga50 = nga50_value, misassembly_proxy = misassembly_value,
        genome_fraction = gf, avg_identity = ident,
        aligned_reference_bases = aligned_reference_bases,
        total_snps = total_snps,
        total_indel_bases = total_indel_bases,
        mismatches_per_100kbp = mm,
        indel_bases_per_100kbp = indel_bases_100k,
        dnadiff_report_path = report_relative_path,
        dnadiff_report_sha256 = dnadiff_report_sha256)
end

"""Apply the interim hybrid-vs-naive and hybrid-vs-scalable engineering gate."""
function validate_comparative_gate(
        df::DataFrames.DataFrame,
        targets::AbstractVector{<:AbstractString})::Nothing
    # This interim gate asserts hybrid-olc vs naive AND vs scalable. If a run
    # selected a subset of arms that omits any of the three (e.g. MYCELIA_RDV_ARMS
    # drops :naive because it does not scale — td-ck03), the gate does not apply;
    # skip it rather than error on a missing arm. validate_results still enforces
    # per-arm evidence completeness for whatever arms did run.
    required = Set(["naive", "scalable", "hybrid-olc"])
    if !issubset(required, Set(df.arm))
        @info "validate_comparative_gate skipped: run does not include all of " *
              "naive/scalable/hybrid-olc" arms = sort(unique(df.arm))
        return nothing
    end
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
        hybrid.indel_bases_per_100kbp <= naive.indel_bases_per_100kbp || error(
            "hybrid-olc indel-base rate exceeds naive for $(target)",
        )

        hybrid.n50 >= SCALABLE_N50_RETENTION * scalable.n50 || error(
            "hybrid-olc retains less than $(100 * SCALABLE_N50_RETENTION)% of " *
            "scalable N50 for $(target)",
        )
        hybrid.avg_identity >= scalable.avg_identity - IDENTITY_TOLERANCE || error(
            "hybrid-olc identity materially regressed from scalable for $(target)",
        )
        hybrid.genome_fraction >=
        scalable.genome_fraction -
        SCALABLE_GENOME_FRACTION_TOLERANCE || error(
            "hybrid-olc genome fraction regressed by more than " *
            "$(SCALABLE_GENOME_FRACTION_TOLERANCE) point from scalable for $(target)",
        )
        hybrid.mismatches_per_100kbp <= scalable.mismatches_per_100kbp || error(
            "hybrid-olc mismatch rate exceeds scalable for $(target)",
        )
        hybrid.indel_bases_per_100kbp <= scalable.indel_bases_per_100kbp || error(
            "hybrid-olc indel-base rate exceeds scalable for $(target)",
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
        :assembly_runtime_s,
        :genome_fraction,
        :avg_identity,
        :mismatches_per_100kbp,
        :indel_bases_per_100kbp
    ]
    for column in finite_columns
        all(isfinite, df[!, column]) || error("non-finite metric in column $(column)")
    end

    all(df.assembly_runtime_s .>= 0.0) ||
        error("assembly_runtime_s must be non-negative")
    all(df.raw_n_contigs .>= df.n_contigs) ||
        error("raw_n_contigs must be at least n_contigs")
    all(df.n_contigs .> 0) || error("every arm must produce at least one contig")
    all(df.total_length .> 0) || error("every arm must produce positive total length")
    all(df.largest_contig .> 0) || error("every arm must produce a positive largest contig")
    all(df.n50 .> 0) || error("every arm must produce a positive N50")
    all(df.n50 .>= df.min_contig_length) ||
        error("n50 must be at least min_contig_length")
    all(df.largest_contig .>= df.min_contig_length) ||
        error("largest_contig must be at least min_contig_length")
    all(df.n50 .<= df.largest_contig) || error("n50 must not exceed largest_contig")
    all(df.largest_contig .<= df.total_length) ||
        error("largest_contig must not exceed total_length")
    all(df.total_length .>= df.n_contigs .* df.min_contig_length) ||
        error("total_length must cover every contig at the common minimum length")
    all(df.total_length .<= df.n_contigs .* df.largest_contig) ||
        error("total_length exceeds the maximum possible from n_contigs and largest_contig")
    all((0.0 .< df.genome_fraction) .& (df.genome_fraction .<= 100.0)) ||
        error("genome_fraction must be in (0, 100]")
    all((0.0 .< df.avg_identity) .& (df.avg_identity .<= 100.0)) ||
        error("avg_identity must be in (0, 100]")
    all(df.mismatches_per_100kbp .>= 0.0) ||
        error("mismatches_per_100kbp must be non-negative")
    all(df.indel_bases_per_100kbp .>= 0.0) ||
        error("indel_bases_per_100kbp must be non-negative")
    all((0 .< df.aligned_reference_bases) .&
        (df.aligned_reference_bases .<= df.reference_length)) ||
        error("aligned_reference_bases must be in (0, reference_length]")
    all(df.total_snps .>= 0) || error("total_snps must be non-negative")
    all(df.total_indel_bases .>= 0) ||
        error("total_indel_bases must be non-negative")
    all(df.nga50 .>= 0) || error("nga50 must be non-negative")
    all(df.nga50 .<= df.reference_length) ||
        error("nga50 must not exceed reference_length")
    all(df.misassembly_proxy .>= 0) ||
        error("misassembly_proxy must be non-negative")
    all(value -> occursin(r"^[0-9a-f]{64}$", value), df.dnadiff_report_sha256) ||
        error("invalid dnadiff report SHA-256")
    all(value -> !isempty(value) && !isabspath(value), df.dnadiff_report_path) ||
        error("dnadiff report paths must be non-empty repository-relative paths")
    all(
        value -> all(
            part -> !isempty(part) && part != "." && part != "..",
            split(value, "/"; keepempty = true)
        ),
        df.dnadiff_report_path
    ) || error("dnadiff report paths must use canonical path components")
    report_path_pattern = r"^benchmarking/results/(real_data_corrector_validation_\d{8}_\d{6}_\d{3}_p\d+)/dnadiff_reports/([^/]+\.report)$"
    report_path_matches = match.(Ref(report_path_pattern), df.dnadiff_report_path)
    all(value -> value !== nothing, report_path_matches) ||
        error("dnadiff report paths must use the canonical validation artifact layout")
    valid_report_matches = something.(report_path_matches)
    length(unique(df.dnadiff_report_path)) == DataFrames.nrow(df) ||
        error("dnadiff report paths must be unique per validation row")
    report_artifact_names = [value.captures[1] for value in valid_report_matches]
    length(unique(report_artifact_names)) == 1 ||
        error("dnadiff reports must identify one retained bundle for the run")
    expected_report_names = Set(
        "$(String(target))_$(arm_name)_dnadiff.report"
    for target in targets for arm_name in arm_names
    )
    report_names = [value.captures[2] for value in valid_report_matches]
    Set(report_names) == expected_report_names ||
        error("dnadiff report paths do not match the target-arm matrix")
    expected_report_names_by_row = ["$(row.name)_$(row.arm)_dnadiff.report"
                                    for row in DataFrames.eachrow(df)]
    report_names == expected_report_names_by_row ||
        error("dnadiff report filenames must match their target-arm rows")
    expected_mismatches = round.(
        df.total_snps ./ df.aligned_reference_bases .* 1e5;
        digits = 1
    )
    expected_indel_bases = round.(
        df.total_indel_bases ./ df.aligned_reference_bases .* 1e5;
        digits = 1
    )
    expected_genome_fraction = round.(
        df.aligned_reference_bases ./ df.reference_length .* 100.0;
        digits = 2
    )
    all(isapprox.(
        df.genome_fraction,
        expected_genome_fraction;
        atol = 1.0e-8,
        rtol = 0.0
    )) || error("genome_fraction does not match aligned reference bases")
    all(isapprox.(
        df.mismatches_per_100kbp,
        expected_mismatches;
        atol = 1.0e-8,
        rtol = 0.0
    )) || error("mismatches_per_100kbp does not match raw dnadiff counts")
    all(isapprox.(
        df.indel_bases_per_100kbp,
        expected_indel_bases;
        atol = 1.0e-8,
        rtol = 0.0
    )) || error("indel_bases_per_100kbp does not match raw dnadiff counts")

    all(df.art_profile .== ART_PROFILE) ||
        error("ART profile must match the configured profile")
    all(df.requested_coverage .== COVERAGE) ||
        error("requested coverage must match the configured coverage")
    all(df.read_length .== READ_LENGTH) ||
        error("read length must match the configured read length")
    all(df.fragment_mean .== FRAGMENT_MEAN) ||
        error("fragment mean must match the configured fragment mean")
    all(df.fragment_sd .== FRAGMENT_SD) ||
        error("fragment standard deviation must match the configured value")
    all(df.art_seed .== SEED) || error("ART seed must match the configured seed")
    all(df.k .== K) || error("k must match the configured k-mer size")
    all(df.min_contig_length .== df.k .+ 1) ||
        error("all arms must use the common k+1 minimum contig length")
    all(df.simulated_layout .== SIMULATED_LAYOUT) ||
        error("simulated layout must match the configured layout")
    all(df.assembly_input_layout .== ASSEMBLY_INPUT_LAYOUT) ||
        error("assembly input layout must match the configured layout")
    all(df.input_read_order .== "R1_then_R2") ||
        error("input read order must be R1_then_R2")
    all(df.input_reads_sha256_semantics .== INPUT_READS_SHA256_SEMANTICS) ||
        error("unexpected input-read SHA-256 semantics")
    all(df.promotion_eligible .== PROMOTION_ELIGIBLE) ||
        error("promotion eligibility does not match the configured run")
    all(df.read_count .> 0) || error("read_count must be positive")
    all(df.read_bases .> 0) || error("read_bases must be positive")
    all(df.read_bases .== df.read_count .* df.read_length) ||
        error("read_bases must equal read_count times read_length")
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
        :input_reads_sha256_semantics
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
    length(unique(df.run_platform)) == 1 ||
        error("run_platform must identify one platform for the run")
    all(!isempty, df.run_platform) || error("run_platform must not be empty")
    length(unique(df.julia_version)) == 1 ||
        error("julia_version must identify one Julia version for the run")
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
    stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS_sss")
    artifact_id = "$(stamp)_p$(getpid())"
    artifact_name = "real_data_corrector_validation_$(artifact_id)"
    artifact_relative = "benchmarking/results/$(artifact_name)"
    artifact_path = joinpath(results_dir, artifact_name)
    ispath(artifact_path) && error(
        "refusing to overwrite existing validation artifact: $(artifact_path)",
    )
    artifact_staging_dir = mktempdir(
        results_dir;
        prefix = ".rdv_artifact_staging_",
        cleanup = true
    )
    report_staging_dir = joinpath(artifact_staging_dir, "dnadiff_reports")
    mkpath(report_staging_dir)
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
            compressed = false
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
            quiet = true
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
            input_reads_sha256_semantics = INPUT_READS_SHA256_SEMANTICS,
            promotion_eligible = PROMOTION_ELIGIBLE
        )

        # 3. Every arm receives this exact `reads` vector. The common k+1
        # post-filter is also passed explicitly to MEGAHIT in `run_arm`.
        println("[3/3] Assembling (naive vs scalable vs hybrid-olc, k=$(K)) + " *
                "dnadiff ...")
        target_rows = NamedTuple[]
        for arm in ARMS
            tag = arm_name(arm)
            report_filename = "$(name)_$(tag)_dnadiff.report"
            row = merge(
                run_arm(
                    name,
                    reads,
                    arm,
                    reference,
                    ref_len,
                    joinpath(target_dir, tag),
                    joinpath(report_staging_dir, report_filename),
                    "$(artifact_relative)/dnadiff_reports/$(report_filename)"
                ),
                provenance
            )
            push!(rows, row)
            push!(target_rows, row)
        end

        for row in target_rows
            println("   [$(row.arm)] contigs=$(row.n_contigs)/$(row.raw_n_contigs) " *
                    "total=$(row.total_length) N50=$(row.n50) " *
                    "NGA50=$(row.nga50) misassembly=$(row.misassembly_proxy) " *
                    "largest=$(row.largest_contig) GF=$(row.genome_fraction)% " *
                    "ident=$(row.avg_identity)% " *
                    "mm/100kb=$(row.mismatches_per_100kbp) " *
                    "indel-bases/100kb=$(row.indel_bases_per_100kbp) " *
                    "assembly_t=$(row.assembly_runtime_s)s")
        end
    end

    isempty(rows) && error("no validation results produced")

    df = DataFrames.DataFrame(rows)
    n_rows = DataFrames.nrow(df)
    df.art_build = fill(conda_package_record("art", "art"), n_rows)
    df.megahit_build = fill(conda_package_record("megahit", "megahit"), n_rows)
    if :sota_spades in ARMS
        # Only query the spades env when the SOTA arm actually ran (otherwise the
        # env may not exist and conda_package_record would error).
        df.spades_build = fill(conda_package_record("spades", "spades"), n_rows)
    end
    df.mummer_build = fill(conda_package_record("mummer", "mummer"), n_rows)
    df.gfatools_build = fill(conda_package_record("gfatools", "gfatools"), n_rows)
    df.run_platform = fill(Sys.MACHINE, n_rows)
    df.julia_version = fill(string(VERSION), n_rows)

    validate_results(df, TARGETS, ARMS)
    staged_csv_path = joinpath(artifact_staging_dir, "comparison.csv")
    csv_path = joinpath(artifact_path, "comparison.csv")
    write_csv_atomically(staged_csv_path, df)
    mv(artifact_staging_dir, artifact_path; force = false)

    println("\n" * "="^70)
    println("SUMMARY — INTERIM ENGINEERING VALIDATION ONLY")
    println("Real references, same ART-simulated reads; not manuscript H5")
    println("="^70)
    show(df; allrows = true, allcols = true)
    println("\n\nCandidate CSV: $(csv_path)")
    println("dnadiff bundle: $(joinpath(artifact_path, "dnadiff_reports"))")
    println("Registration : UNREGISTERED until the manifest pointer is updated explicitly")
    println("Promotable   : $(PROMOTION_ELIGIBLE)")
    println("Done         : $(Dates.now())")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
