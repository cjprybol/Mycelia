# Track A baseline benchmark — greedy Viterbi on the viral tier (rhizomorph-paper td-bblmi)
#
# Runs the CURRENT default Rhizomorph assembler (greedy path-finding) across the
# viral tier to establish the baseline for all H1-H7 comparisons and to confirm
# the pre-registration power analysis (assumed NGA50 CV ~ 0.15) before the
# registration is locked.
#
# Matrix (288 cells): 4 organisms x 3 technologies x 4 coverages x 3 seeds x 2 decoder arms.
#   - decoder arm "qualmer": simulated FASTQ passed as-is -> assemble_genome auto-enables
#     quality scores -> qualmer graph (the assembler's real default on real reads).
#   - decoder arm "kmer":    quality stripped (FASTQ -> FASTA records) -> plain k-mer graph,
#     matching the existing FASTA-based benchmark and the future DP arm's graph type.
#
# Each cell: simulate reads -> assemble (k=31) -> QUAST vs reference -> parse NGA50 /
# misassemblies / genome fraction / duplication ratio. Per-cell JSON checkpoint enables
# crash-safe resume. A final step computes NGA50 CV per (organism x tech x coverage x arm)
# and writes a pass/fail power-analysis summary.
#
# Usage:
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl              # full 288-cell run
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --smoke      # 1 cell (Lambda/illumina/30x/seed42/qualmer)
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --organisms Lambda,T4 --arms kmer
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --coverages 30,100 --seeds 42 --technologies illumina,ont
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --output-dir /scratch/track_a
#
# Shard flags (--organisms/--technologies/--coverages/--seeds/--arms) take comma-separated
# values and compose, so an HPC array job can split the matrix and share one results tree.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import DataFrames
import CSV
import JSON
import Dates
import Statistics

# Optional provenance writer (git SHA + tool versions). Degrade gracefully if unavailable.
const HAVE_ARTIFACT_WRITER = try
    include(joinpath(@__DIR__, "benchmark_artifacts.jl"))
    true
catch e
    @warn "benchmark_artifacts.jl unavailable; skipping provenance bundle" exception = e
    false
end

# === Configuration ===

# (display_name, NCBI accession, expected_size_bp). Lambda/T4/SARS-CoV-2 accessions match
# benchmarking/real_genome_benchmark.jl; phi29 = Bacillus phage phi29 (NC_011048, ~19.3 kb).
const ORGANISMS = [
    ("Lambda", "NC_001416", 48_502),
    ("T4", "NC_000866", 168_903),
    ("phi29", "NC_011048", 19_282),
    ("SARS-CoV-2", "NC_045512", 29_903)
]
const TECHNOLOGIES = ["illumina", "pacbio", "ont"]
const COVERAGES = [10, 30, 50, 100]
const SEEDS = [42, 123, 456]
const DECODER_ARMS = ["qualmer", "kmer"]
const K = 31
const CV_THRESHOLD = 0.15  # assumed NGA50 coefficient of variation in the power analysis

# Canonical row schema (fixed order so in-memory and JSON-reloaded rows align in the DataFrame).
const ROW_KEYS = (
    :organism, :accession, :technology, :coverage, :seed, :decoder_arm, :k,
    :n_reads, :n_contigs, :NGA50, :misassemblies, :genome_fraction,
    :duplication_ratio, :largest_contig, :wall_seconds, :peak_rss_bytes, :status
)
const INT_KEYS = (
    :coverage, :seed, :k, :n_reads, :n_contigs, :largest_contig, :peak_rss_bytes)
const FLOAT_KEYS = (
    :NGA50, :misassemblies, :genome_fraction, :duplication_ratio, :wall_seconds)
const STR_KEYS = (:organism, :accession, :technology, :decoder_arm, :status)

# === Argument parsing ===

arg_value(flag) =
    let i = findfirst(==(flag), ARGS)
        (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
    end
arg_list(flag) =
    let v = arg_value(flag)
        v === nothing ? nothing : String.(split(v, ","))
    end

const SMOKE = "--smoke" in ARGS

organisms = ORGANISMS
technologies = TECHNOLOGIES
coverages = COVERAGES
seeds = SEEDS
arms = DECODER_ARMS

if SMOKE
    organisms = ORGANISMS[1:1]      # Lambda (accession verified in-repo)
    technologies = ["illumina"]
    coverages = [30]
    seeds = [42]
    arms = ["qualmer"]
else
    # NOTE: assign directly (no `let`). A `let` block introduces a new scope, so
    # `coverages = ...` inside it would create a local shadow and silently leave
    # the global matrix at its default — i.e. the shard/filter flags would be
    # ignored. `if/else` bodies do not introduce scope, so these reach the globals.
    _f = arg_list("--organisms")
    _f !== nothing && (organisms = filter(o -> o[1] in _f, ORGANISMS))
    _f = arg_list("--technologies")
    _f !== nothing && (technologies = _f)
    _f = arg_list("--coverages")
    _f !== nothing && (coverages = parse.(Int, _f))
    _f = arg_list("--seeds")
    _f !== nothing && (seeds = parse.(Int, _f))
    _f = arg_list("--arms")
    _f !== nothing && (arms = _f)
end

const OUTPUT_DIR = let v = arg_value("--output-dir")
    v === nothing ? joinpath(@__DIR__, "results", "track_a_baseline") : v
end

const N_CELLS = length(organisms) * length(technologies) * length(coverages) *
                length(seeds) * length(arms)

println("=== Track A baseline benchmark ===")
println("Start: $(Dates.now())")
println("Smoke mode: $SMOKE")
println("Organisms: $(join((o[1] for o in organisms), ", "))")
println("Technologies: $(join(technologies, ", "))")
println("Coverages: $(join(coverages, ", "))x")
println("Seeds: $(join(seeds, ", "))")
println("Decoder arms: $(join(arms, ", "))")
println("Cells to run: $N_CELLS")
println("Output dir: $OUTPUT_DIR")

# === Metrics + row helpers ===

function empty_metrics()
    (NGA50 = 0.0, misassemblies = 0.0, genome_fraction = 0.0,
        duplication_ratio = 0.0, largest_contig = 0.0)
end

function cell_row(org, acc, tech, cov, seed, arm; n_reads, n_contigs,
        wall_seconds, peak_rss_bytes, metrics, status)
    return (
        organism = String(org), accession = String(acc), technology = String(tech),
        coverage = Int(cov), seed = Int(seed), decoder_arm = String(arm), k = K,
        n_reads = Int(n_reads), n_contigs = Int(n_contigs),
        NGA50 = Float64(metrics.NGA50), misassemblies = Float64(metrics.misassemblies),
        genome_fraction = Float64(metrics.genome_fraction),
        duplication_ratio = Float64(metrics.duplication_ratio),
        largest_contig = Int(round(Float64(metrics.largest_contig))),
        wall_seconds = round(Float64(wall_seconds); digits = 3),
        peak_rss_bytes = Int(peak_rss_bytes), status = String(status)
    )
end

# Rebuild a canonical, type-coerced NamedTuple from a parsed JSON dict (resumed cells).
function canonical(d::AbstractDict)
    vals = map(ROW_KEYS) do key
        v = d[String(key)]
        if key in INT_KEYS
            v isa Integer ? Int(v) : Int(round(Float64(v)))
        elseif key in FLOAT_KEYS
            Float64(v)
        else
            String(v)
        end
    end
    return (; (ROW_KEYS .=> vals)...)
end

function save_cell_json(path, row)
    open(path, "w") do io
        JSON.print(io, Dict(string(k) => v for (k, v) in pairs(row)), 2)
    end
    return path
end

# === QUAST report parsing ===

# report.tsv: column 1 = metric label, column 2 = this assembly's value. NGA50 and
# misassembly rows appear only when QUAST is run with a reference.
function parse_quast_report(report_tsv)
    isfile(report_tsv) || return empty_metrics()
    df = CSV.read(report_tsv, DataFrames.DataFrame; delim = '\t', header = true)
    DataFrames.ncol(df) >= 2 || return empty_metrics()
    labelcol, valcol = names(df)[1], names(df)[2]
    getval(metric) =
        let i = findfirst(==(metric), df[!, labelcol])
            i === nothing ? missing : df[i, valcol]
        end
    tonum(x) = x === missing ? 0.0 :
               x isa Number ? Float64(x) :
               something(tryparse(Float64, strip(string(x))), 0.0)
    return (
        NGA50 = tonum(getval("NGA50")),
        misassemblies = tonum(getval("# misassemblies")),
        genome_fraction = tonum(getval("Genome fraction (%)")),
        duplication_ratio = tonum(getval("Duplication ratio")),
        largest_contig = tonum(getval("Largest contig"))
    )
end

# === Read simulation (per-technology adapter; APIs are non-uniform) ===

function simulate_reads(tech, ref_fasta, cov, seed, reads_dir)
    mkpath(reads_dir)
    if tech == "illumina"
        outbase = joinpath(reads_dir, "illumina_$(cov)x")
        res = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta, coverage = cov, outbase = outbase,
            rndSeed = seed, paired = true, quiet = true)
        recs = FASTX.FASTQ.Record[]
        for p in (res.forward_reads, res.reverse_reads)
            p === nothing && continue
            reader = Mycelia.open_fastx(p)
            append!(recs, collect(reader))
            close(reader)
        end
        return recs
    elseif tech == "pacbio"
        fq = Mycelia.simulate_pacbio_reads(
            fasta = ref_fasta, quantity = "$(cov)x",
            outfile = joinpath(reads_dir, "pacbio_$(cov)x.fq.gz"),
            seed = seed, quiet = true)
        reader = Mycelia.open_fastx(fq)
        recs = collect(reader)
        close(reader)
        return recs
    elseif tech == "ont"
        fq = Mycelia.simulate_nanopore_reads(
            fasta = ref_fasta, quantity = "$(cov)x",
            outfile = joinpath(reads_dir, "ont_$(cov)x.fq.gz"),
            seed = seed, quiet = true)
        reader = Mycelia.open_fastx(fq)
        recs = collect(reader)
        close(reader)
        return recs
    else
        error("unknown technology: $tech")
    end
end

# Quality-on arm keeps FASTQ records; quality-off arm strips quality to FASTA records.
function shape_for_arm(records, arm)
    if arm == "qualmer"
        return records
    elseif arm == "kmer"
        return [FASTX.FASTA.Record(FASTX.identifier(r), FASTX.sequence(r)) for r in records]
    else
        error("unknown decoder arm: $arm")
    end
end

# === Per-cell execution ===

function run_cell(org, acc, ref, tech, cov, seed, arm, cell_dir)
    records = simulate_reads(tech, ref, cov, seed, joinpath(cell_dir, "reads"))
    n_reads = length(records)
    asm_input = shape_for_arm(records, arm)

    GC.gc()
    rss0 = Sys.maxrss()
    timed = @timed Mycelia.Rhizomorph.assemble_genome(asm_input; k = K, verbose = false)
    result = timed.value
    wall_seconds = timed.time
    peak_rss_bytes = max(0, Sys.maxrss() - rss0)

    contigs_path = joinpath(cell_dir, "contigs.fasta")
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # Rhizomorph contigs are String
        end
    end
    n_contigs = length(result.contigs)

    metrics = if n_contigs > 0 && filesize(contigs_path) > 0
        try
            quast_dir = joinpath(cell_dir, "quast")
            Mycelia.run_quast([contigs_path]; outdir = quast_dir, reference = ref, min_contig = 500)
            parse_quast_report(joinpath(quast_dir, "report.tsv"))
        catch e
            @warn "QUAST failed" cell=basename(cell_dir) exception=(e, catch_backtrace())
            empty_metrics()
        end
    else
        empty_metrics()
    end

    status = n_contigs == 0 ? "empty_assembly" : "ok"
    return cell_row(org, acc, tech, cov, seed, arm;
        n_reads, n_contigs, wall_seconds, peak_rss_bytes, metrics, status)
end

# === Aggregation ===

function write_aggregate(root, rows)
    df = DataFrames.DataFrame(rows)
    CSV.write(joinpath(root, "track_a_results.tsv"), df; delim = '\t')
    return df
end

function write_power_analysis(root, df)
    cv_rows = NamedTuple[]
    for g in DataFrames.groupby(df, [:organism, :technology, :coverage, :decoder_arm])
        nga = Float64.(g.NGA50)
        m = Statistics.mean(nga)
        s = length(nga) > 1 ? Statistics.std(nga; corrected = true) : NaN
        cv = (m == 0 || isnan(s)) ? NaN : s / m
        push!(cv_rows,
            (
                organism = g.organism[1], technology = g.technology[1],
                coverage = g.coverage[1], decoder_arm = g.decoder_arm[1],
                n = length(nga), mean_nga50 = round(m; digits = 1),
                sd_nga50 = round(s; digits = 1), cv_nga50 = round(cv; digits = 4),
                passes = (isfinite(cv) && cv <= CV_THRESHOLD)))
    end
    cv_df = DataFrames.DataFrame(cv_rows)
    CSV.write(joinpath(root, "power_analysis_cv.tsv"), cv_df; delim = '\t')

    evaluable = filter(r -> isfinite(r.cv_nga50), cv_rows)
    n_eval = length(evaluable)
    n_pass = count(r -> r.passes, evaluable)
    max_cv = isempty(evaluable) ? NaN : maximum(r.cv_nga50 for r in evaluable)
    verdict = (n_eval > 0 && n_pass == n_eval) ? "supported" :
              n_eval == 0 ? "indeterminate (no evaluable cells)" : "NOT fully supported"

    open(joinpath(root, "power_analysis_summary.md"), "w") do io
        println(io, "# Track A power-analysis check — NGA50 CV vs assumed $(CV_THRESHOLD)\n")
        println(io, "Generated: $(Dates.now())\n")
        println(io, "- Evaluable cells (organism×tech×coverage×arm, ≥2 seeds, nonzero NGA50): $n_eval")
        println(io, "- Cells with CV ≤ $(CV_THRESHOLD): $n_pass / $n_eval")
        println(io, "- Max CV observed: $(round(max_cv; digits = 4))")
        println(io, "- **Verdict: assumed CV ≈ $(CV_THRESHOLD) is $verdict.**\n")
        fails = filter(r -> isfinite(r.cv_nga50) && !r.passes, cv_rows)
        if !isempty(fails)
            println(io, "## Cells exceeding CV $(CV_THRESHOLD)\n")
            for r in fails
                println(io,
                    "- $(r.organism) / $(r.technology) / $(r.coverage)x / $(r.decoder_arm): " *
                    "CV=$(r.cv_nga50) (mean NGA50 $(r.mean_nga50), n=$(r.n))")
            end
            println(io)
        end
        println(io,
            "_Caveat: n=3 seeds makes each CV estimate noisy; treat this as directional " *
            "evidence for the power analysis, not a definitive variance estimate._")
    end
    return cv_df
end

# === Main ===

mkpath(OUTPUT_DIR)
refs_dir = joinpath(OUTPUT_DIR, "refs")
cells_dir = joinpath(OUTPUT_DIR, "cells")
mkpath(refs_dir)
mkpath(cells_dir)

println("\n--- Phase 1: download reference genomes ---")
ref_paths = Dict{String, String}()
for (org, acc, expected) in organisms
    haskey(ref_paths, org) && continue
    print("  $org ($acc) ... ")
    rp = Mycelia.download_genome_by_accession(accession = acc, outdir = refs_dir, compressed = false)
    if !isfile(rp) || filesize(rp) == 0
        error("download failed for $org ($acc): $rp")
    end
    actual = try
        Mycelia.total_fasta_size(rp)
    catch
        -1
    end
    if actual > 0 && abs(actual - expected) > 0.2 * expected
        @warn "reference size mismatch" organism=org expected=expected actual=actual
    end
    ref_paths[org] = rp
    println("$(actual > 0 ? actual : "?") bp")
end

println("\n--- Phase 2: matrix ($N_CELLS cells) ---")
rows = NamedTuple[]
cell_index = 0
for (org, acc, _expected) in organisms,
    tech in technologies,
    cov in coverages,
    seed in seeds,
    arm in arms
    global cell_index += 1
    cell_id = "$(org)__$(tech)__$(cov)x__seed$(seed)__$(arm)"
    cell_dir = joinpath(cells_dir, cell_id)
    ckpt = joinpath(cell_dir, "cell_result.json")

    if isfile(ckpt)
        println("  [$cell_index/$N_CELLS] $cell_id — cached, skipping")
        push!(rows, canonical(JSON.parsefile(ckpt)))
        continue
    end

    mkpath(cell_dir)
    print("  [$cell_index/$N_CELLS] $cell_id ... ")
    row = try
        run_cell(org, acc, ref_paths[org], tech, cov, seed, arm, cell_dir)
    catch e
        @warn "cell failed" cell_id exception = (e, catch_backtrace())
        cell_row(org, acc, tech, cov, seed, arm;
            n_reads = 0, n_contigs = 0, wall_seconds = 0.0, peak_rss_bytes = 0,
            metrics = empty_metrics(), status = "error")
    end
    save_cell_json(ckpt, row)
    push!(rows, row)
    write_aggregate(OUTPUT_DIR, rows)  # rewrite aggregate each cell (cheap at this scale)
    println("$(row.status): $(row.n_contigs) contigs, NGA50=$(row.NGA50), " *
            "GF=$(row.genome_fraction)%, $(round(row.wall_seconds; digits = 1))s")
end

println("\n--- Phase 3: aggregate + power analysis ---")
results_df = write_aggregate(OUTPUT_DIR, rows)
cv_df = write_power_analysis(OUTPUT_DIR, results_df)

if HAVE_ARTIFACT_WRITER
    try
        write_benchmark_artifacts(
            ["track_a_results" => results_df, "track_a_power_analysis_cv" => cv_df];
            output_dir = joinpath(OUTPUT_DIR, "artifacts"),
            run_id = "track_a_baseline_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))",
            scale = SMOKE ? "smoke" : "full",
            dataset_ids = [o[2] for o in organisms],
            command_args = ARGS,
            metadata = Dict("benchmark" => "track_a_baseline", "k" => K,
                "cv_threshold" => CV_THRESHOLD)
        )
    catch e
        @warn "provenance artifact bundle failed (results TSVs still written)" exception = e
    end
end

n_ok = count(r -> r.status == "ok", rows)
println("\nDone: $(length(rows)) cells, $n_ok ok. Results in $OUTPUT_DIR")
println("End: $(Dates.now())")
