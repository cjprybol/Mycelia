"""
Generate coverage baseline from Julia .cov files.

Parses all .cov files under src/, merges multiple PIDs per source file,
computes per-module and per-tier coverage stats, and outputs:
- coverage/baseline.json
- coverage/baseline-summary.md
- coverage/exclusions.json
"""

import JSON
import Dates

const COVERAGE_DIR = dirname(@__FILE__)
const SRC_DIR = joinpath(dirname(COVERAGE_DIR), "src")

# Tier assignments from campaign spec
const TIER_MAP = Dict{String, Tuple{Int, String, Int}}(
    # filename => (tier_number, tier_name, target_pct)
    # Tier 1: Data Acquisition (90%)
    "reference-databases.jl" => (1, "Data Acquisition", 90),
    "simulation.jl" => (1, "Data Acquisition", 90),
    # Tier 2: Preprocessing/QC (90%)
    "alphabets.jl" => (2, "Preprocessing/QC", 90),
    "constants.jl" => (2, "Preprocessing/QC", 90),
    "fastx.jl" => (2, "Preprocessing/QC", 90),
    # Tier 3: K-mer Analysis (95%)
    "kmer-analysis.jl" => (3, "K-mer Analysis", 95),
    "kmer-saturation-analysis.jl" => (3, "K-mer Analysis", 95),
    "qualmer-analysis.jl" => (3, "K-mer Analysis", 95),
    "distance-metrics.jl" => (3, "K-mer Analysis", 95),
    # Tier 4: Assembly Core (95%)
    "graph-cleanup.jl" => (4, "Assembly Core", 95),
    "iterative-assembly.jl" => (4, "Assembly Core", 95),
    # Tier 5: Validation (80%)
    "quality-control-and-benchmarking.jl" => (5, "Validation", 80),
    # Tier 6: Annotation (80%)
    "annotation.jl" => (6, "Annotation", 80),
    "genome-features.jl" => (6, "Annotation", 80),
    # Tier 7: Comparative (85%)
    "pangenome-analysis.jl" => (7, "Comparative", 85),
    "sequence-comparison.jl" => (7, "Comparative", 85),
    # Tier 8: Tool Integration (75%)
    "alignments-and-mapping.jl" => (8, "Tool Integration", 75),
    "classification.jl" => (8, "Tool Integration", 75),
    "binning.jl" => (8, "Tool Integration", 75),
    # Tier 9: Analysis/Viz (85%)
    "clustering.jl" => (9, "Analysis/Viz", 85),
    "dimensionality-reduction.jl" => (9, "Analysis/Viz", 85),
    "plotting-and-visualization.jl" => (9, "Analysis/Viz", 85),
    # Tier 10: Infrastructure (90%)
    "utility-functions.jl" => (10, "Infrastructure", 90),
    "checkpointing.jl" => (10, "Infrastructure", 90),
    "slurm-sbatch.jl" => (10, "Infrastructure", 90)
)

# Rhizomorph submodules all go to Tier 4
const RHIZOMORPH_FILES = [
    "contigs.jl", "error-correction.jl", "generation.jl", "io.jl",
    "metrics.jl", "path-finding.jl", "repeats.jl", "simplification.jl",
    "information-theory.jl", "sequence-quality.jl", "assembly.jl",
    "edge-data.jl", "evidence-functions.jl", "evidence-structures.jl",
    "graph-construction.jl", "graph-query.jl", "graph-type-conversions.jl",
    "quality-functions.jl", "vertex-data.jl", "kmer-graphs.jl",
    "ngram-graphs.jl", "qualmer-graphs.jl", "rhizomorph.jl",
    "fasta-graphs.jl", "fastq-graphs.jl", "strand-conversions.jl",
    "string-graphs.jl"
]

# Excluded from coverage targets (but still measured)
const EXCLUDED_FILES = Set([
    "Mycelia.jl",              # module loader
    "precompile_workload.jl",  # precompilation statements
    "testing-utilities.jl",    # test helpers
    "bioconda.jl",             # external tool installation (network)
    "execution.jl",            # shell execution wrappers
    "slurm-templates.jl"      # HPC template generation
])

function parse_cov_file(path::String)::Dict{Int, Int}
    """Parse a .cov file, returning line_number => hit_count for executable lines."""
    hits = Dict{Int, Int}()
    for (i, line) in enumerate(eachline(path))
        stripped = lstrip(line)
        if startswith(stripped, '-') || isempty(stripped)
            continue
        end
        # First token is the hit count
        parts = split(stripped, limit = 2)
        count = tryparse(Int, parts[1])
        if count !== nothing
            hits[i] = count
        end
    end
    return hits
end

function find_cov_files(src_dir::String)::Dict{String, Vector{String}}
    """Find all .cov files and group by source file path."""
    grouped = Dict{String, Vector{String}}()
    for (root, _, files) in walkdir(src_dir)
        for f in files
            m = match(r"^(.+\.jl)\.\d+\.cov$", f)
            if m !== nothing
                source_name = m.captures[1]
                rel_path = joinpath(root, source_name)
                rel_path = relpath(rel_path, dirname(src_dir))
                cov_path = joinpath(root, f)
                if !haskey(grouped, rel_path)
                    grouped[rel_path] = String[]
                end
                push!(grouped[rel_path], cov_path)
            end
        end
    end
    return grouped
end

function merge_cov_files(paths::Vector{String})::Dict{Int, Int}
    """Merge multiple .cov files for the same source, taking max hit count per line."""
    merged = Dict{Int, Int}()
    for path in paths
        hits = parse_cov_file(path)
        for (line, count) in hits
            merged[line] = max(get(merged, line, 0), count)
        end
    end
    return merged
end

function get_tier(filename::String, rel_path::String)::Tuple{Int, String, Int}
    # Check rhizomorph submodules
    if contains(rel_path, "rhizomorph")
        return (4, "Assembly Core", 95)
    end
    return get(TIER_MAP, filename, (0, "Unassigned", 0))
end

function find_uncovered_sources(src_dir::String, covered::Set{String})::Vector{Tuple{
        String, Int}}
    """Find .jl source files with no .cov output and count their executable lines."""
    uncovered = Tuple{String, Int}[]
    for (root, _, files) in walkdir(src_dir)
        for f in files
            if endswith(f, ".jl") && !endswith(f, ".cov")
                rel_path = relpath(joinpath(root, f), dirname(src_dir))
                if rel_path ∉ covered
                    # Count executable lines (non-blank, non-comment)
                    n_exec = 0
                    for line in eachline(joinpath(root, f))
                        stripped = strip(line)
                        if !isempty(stripped) && !startswith(stripped, '#') &&
                           !startswith(stripped, "\"\"\"") && stripped != "\"\"\""
                            n_exec += 1
                        end
                    end
                    push!(uncovered, (rel_path, n_exec))
                end
            end
        end
    end
    return uncovered
end

function main()
    println("Scanning for .cov files in $(SRC_DIR)...")
    grouped = find_cov_files(SRC_DIR)
    println("Found $(length(grouped)) source files with coverage data")

    modules = Dict{String, Any}[]
    tier_stats = Dict{Int, Dict{String, Any}}()
    covered_paths = Set{String}(keys(grouped))

    for (rel_path, cov_paths) in sort(collect(grouped))
        merged = merge_cov_files(cov_paths)
        lines_total = length(merged)
        lines_hit = count(v -> v > 0, values(merged))
        coverage_pct = lines_total > 0 ?
                       round(100.0 * lines_hit / lines_total; digits = 1) : 0.0

        filename = basename(rel_path)
        is_excluded = filename in EXCLUDED_FILES
        tier_num, tier_name, target_pct = get_tier(filename, rel_path)

        mod = Dict{String, Any}(
            "file" => rel_path,
            "filename" => filename,
            "lines_total" => lines_total,
            "lines_hit" => lines_hit,
            "lines_missed" => lines_total - lines_hit,
            "coverage_pct" => coverage_pct,
            "tier" => tier_num,
            "tier_name" => tier_name,
            "target_pct" => target_pct,
            "excluded" => is_excluded,
            "cov_file_count" => length(cov_paths)
        )
        push!(modules, mod)

        # Aggregate tier stats (skip excluded and unassigned)
        if !is_excluded && tier_num > 0
            if !haskey(tier_stats, tier_num)
                tier_stats[tier_num] = Dict{String, Any}(
                    "tier" => tier_num,
                    "name" => tier_name,
                    "target_pct" => target_pct,
                    "lines_total" => 0,
                    "lines_hit" => 0,
                    "modules" => String[]
                )
            end
            ts = tier_stats[tier_num]
            ts["lines_total"] += lines_total
            ts["lines_hit"] += lines_hit
            push!(ts["modules"], filename)
        end
    end

    # Add 0% coverage modules (source files with no .cov output)
    uncovered = find_uncovered_sources(SRC_DIR, covered_paths)
    println("Found $(length(uncovered)) source files with zero coverage")
    for (rel_path, n_exec) in uncovered
        filename = basename(rel_path)
        is_excluded = filename in EXCLUDED_FILES
        tier_num, tier_name, target_pct = get_tier(filename, rel_path)

        mod = Dict{String, Any}(
            "file" => rel_path,
            "filename" => filename,
            "lines_total" => n_exec,
            "lines_hit" => 0,
            "lines_missed" => n_exec,
            "coverage_pct" => 0.0,
            "tier" => tier_num,
            "tier_name" => tier_name,
            "target_pct" => target_pct,
            "excluded" => is_excluded,
            "cov_file_count" => 0
        )
        push!(modules, mod)

        if !is_excluded && tier_num > 0 && n_exec > 0
            if !haskey(tier_stats, tier_num)
                tier_stats[tier_num] = Dict{String, Any}(
                    "tier" => tier_num,
                    "name" => tier_name,
                    "target_pct" => target_pct,
                    "lines_total" => 0,
                    "lines_hit" => 0,
                    "modules" => String[]
                )
            end
            ts = tier_stats[tier_num]
            ts["lines_total"] += n_exec
            push!(ts["modules"], filename)
        end
    end

    # Compute tier coverage percentages
    tiers = Dict{String, Any}[]
    for tier_num in sort(collect(keys(tier_stats)))
        ts = tier_stats[tier_num]
        ts["coverage_pct"] = ts["lines_total"] > 0 ?
                             round(100.0 * ts["lines_hit"] / ts["lines_total"]; digits = 1) :
                             0.0
        ts["lines_missed"] = ts["lines_total"] - ts["lines_hit"]
        ts["gap_pct"] = round(ts["target_pct"] - ts["coverage_pct"]; digits = 1)
        push!(tiers, ts)
    end

    # Aggregate totals
    total_lines = sum(m["lines_total"] for m in modules if !m["excluded"])
    total_hit = sum(m["lines_hit"] for m in modules if !m["excluded"])
    overall_pct = total_lines > 0 ? round(100.0 * total_hit / total_lines; digits = 1) : 0.0

    baseline = Dict{String, Any}(
        "generated_at" => string(Dates.now()),
        "overall" => Dict{String, Any}(
            "lines_total" => total_lines,
            "lines_hit" => total_hit,
            "lines_missed" => total_lines - total_hit,
            "coverage_pct" => overall_pct
        ),
        "tiers" => tiers,
        "modules" => sort(modules; by = m -> (m["tier"], m["filename"]))
    )

    # Write JSON
    json_path = joinpath(COVERAGE_DIR, "baseline.json")
    open(json_path, "w") do io
        JSON.print(io, baseline, 2)
    end
    println("Wrote $(json_path)")

    # Write markdown summary
    md_path = joinpath(COVERAGE_DIR, "baseline-summary.md")
    open(md_path, "w") do io
        println(io, "# Mycelia Coverage Baseline")
        println(io, "")
        println(io, "Generated: $(Dates.now())")
        println(io, "")
        println(io, "## Overall: $(overall_pct)% ($(total_hit)/$(total_lines) executable lines)")
        println(io, "")
        println(io, "## Tier Summary")
        println(io, "")
        println(io, "| Tier | Name | Coverage | Target | Gap | Lines | Modules |")
        println(io, "|------|------|----------|--------|-----|-------|---------|")
        for t in tiers
            println(io,
                "| $(t["tier"]) | $(t["name"]) | $(t["coverage_pct"])% | $(t["target_pct"])% | $(t["gap_pct"])pp | $(t["lines_hit"])/$(t["lines_total"]) | $(length(t["modules"])) |")
        end
        println(io, "")
        println(io, "## Module Detail")
        println(io, "")
        println(io, "| Module | Tier | Coverage | Target | Lines Hit | Lines Total | Excluded |")
        println(io, "|--------|------|----------|--------|-----------|-------------|----------|")
        for m in sort(modules; by = m -> (m["tier"], m["filename"]))
            ex = m["excluded"] ? "yes" : ""
            println(io,
                "| $(m["filename"]) | $(m["tier_name"]) | $(m["coverage_pct"])% | $(m["target_pct"])% | $(m["lines_hit"]) | $(m["lines_total"]) | $(ex) |")
        end
    end
    println("Wrote $(md_path)")

    # Write exclusions
    exclusions = Dict{String, Any}(
        "policy" => "Modules excluded from coverage targets but still measured",
        "exclusions" => [Dict("file" => f, "reason" => _exclusion_reason(f))
                         for f in sort(collect(EXCLUDED_FILES))]
    )
    exc_path = joinpath(COVERAGE_DIR, "exclusions.json")
    open(exc_path, "w") do io
        JSON.print(io, exclusions, 2)
    end
    println("Wrote $(exc_path)")

    # Print summary
    println("\n=== Coverage Baseline Summary ===")
    println("Overall: $(overall_pct)% ($(total_hit)/$(total_lines))")
    println("")
    for t in tiers
        status = t["coverage_pct"] >= t["target_pct"] ? "✓" : "△"
        println("  $(status) Tier $(t["tier"]) $(t["name"]): $(t["coverage_pct"])% / $(t["target_pct"])% target ($(length(t["modules"])) modules)")
    end
end

function _exclusion_reason(filename::String)::String
    reasons = Dict(
        "Mycelia.jl" => "Module loader — no testable logic",
        "precompile_workload.jl" => "Precompilation statements — not runtime code",
        "testing-utilities.jl" => "Test infrastructure helpers — tested transitively",
        "bioconda.jl" => "External tool installation — network-dependent",
        "execution.jl" => "Shell execution wrappers — system-dependent",
        "slurm-templates.jl" => "HPC template generation — environment-dependent"
    )
    return get(reasons, filename, "Excluded per campaign policy")
end

main()
