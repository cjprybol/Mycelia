"""
Generate a coverage baseline for the package's loaded source files.

This script:
0. Activates `coverage/Project.toml` so it can run standalone via `julia coverage/generate_baseline.jl`
1. Optionally runs the full Mycelia test suite with coverage enabled
2. Resolves the actual package source scope from `include(...)` statements
3. Uses Coverage.jl to merge `.cov` files per source file
4. Writes:
   - coverage/baseline.json
   - coverage/baseline-summary.md
   - coverage/exclusions.json

Set `MYCELIA_COVERAGE_COLLECT=false` to skip the test run and summarize the
existing coverage artifacts already present under `src/`.

Set `MYCELIA_COVERAGE_MODE=user` or `MYCELIA_COVERAGE_MODE=all` to opt out of
the default package-scoped coverage mode. The default passes
`--code-coverage=@src` into a direct `test/runtests.jl` subprocess under the
dedicated `coverage/` environment so the standalone collector instruments this
worktree's source tree instead of unrelated Julia/Base paths.

Set `MYCELIA_COVERAGE_NO_COMPILED_MODULES=false` to opt out of the default
no-precompile coverage mode. The default matches the HPC CI coverage
invocation so standalone runs collect `src/*.cov` artifacts reliably.

Set `MYCELIA_COVERAGE_NO_PKGIMAGES=false` to opt out of the default package-image
disablement for the standalone coverage process. The default keeps the
instrumented test run on source-evaluated code paths instead of cached package
images, which is necessary for reliable `src/*.cov` emission.

Set `MYCELIA_COVERAGE_WARM_PRECOMPILE=false` to skip the default warm-up
`Pkg.precompile()` pass before coverage collection. The warm-up keeps the
standalone no-precompile test invocation from spending most of its time in
first-run dependency compilation.

Set `MYCELIA_COVERAGE_SKIP_AQUA=false` to include Aqua.jl during standalone
coverage collection. The default skips Aqua because it validates package QA
policy rather than executable source behavior and dominates no-precompile
coverage runtime without contributing useful line hits.

Set `MYCELIA_COVERAGE_SKIP_JET=false` to include JET.jl during standalone
coverage collection. The default skips JET because it is a static-analysis gate
rather than a runtime line-coverage workload.

Set `MYCELIA_COVERAGE_SKIP_TEST_FILES=` to override the comma-separated list of
test files skipped only during standalone coverage collection. By default the
collector skips `hpc_assembly_wrappers_test.jl` and
`megahit_phix_workflow.jl`, the external validation integration suites
(`checkm_tools.jl`, `coverm_wrappers.jl`, `coverm_integration_extended.jl`,
`mosdepth_coverage_qc.jl`, `coverage_taxonomy_integration.jl`,
`quast_busco_wrappers_test.jl`), plus the `third_party_assemblers*.jl`
integration family, because those tests depend on provisioned conda
environments and large external tool stacks that are orthogonal to baseline
generation on a local workstation.

Coverage collection shells out to `test/runtests.jl` under the dedicated
coverage project so the standalone runner exercises Mycelia directly from this
worktree and emits `.cov` files beside the real `src/*.jl` files that the
baseline parser summarizes.
"""

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

import Coverage
import Dates
import JSON

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const COVERAGE_DIR = @__DIR__
const SRC_DIR = joinpath(ROOT_DIR, "src")
const TEST_DIR = joinpath(ROOT_DIR, "test")
const ROOT_RUNTESTS = joinpath(TEST_DIR, "runtests.jl")
const ROOT_ENTRYPOINT = joinpath(SRC_DIR, "Mycelia.jl")
const PACKAGE_COVERAGE_SCOPE = string('@', SRC_DIR)
const COVERAGE_RUN_ALL = lowercase(get(ENV, "MYCELIA_RUN_ALL", "true")) == "true"
const COVERAGE_RUN_EXTERNAL =
    lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "true")) == "true"
const COVERAGE_COLLECT =
    lowercase(get(ENV, "MYCELIA_COVERAGE_COLLECT", "true")) != "false"
const COVERAGE_NO_COMPILED_MODULES =
    lowercase(get(ENV, "MYCELIA_COVERAGE_NO_COMPILED_MODULES", "true")) == "true"
const COVERAGE_NO_PKGIMAGES =
    lowercase(get(ENV, "MYCELIA_COVERAGE_NO_PKGIMAGES", "true")) == "true"
const COVERAGE_WARM_PRECOMPILE =
    lowercase(get(ENV, "MYCELIA_COVERAGE_WARM_PRECOMPILE", "true")) == "true"
const COVERAGE_SKIP_AQUA =
    lowercase(get(ENV, "MYCELIA_COVERAGE_SKIP_AQUA", "true")) == "true"
const COVERAGE_SKIP_AQUA_AMBIGUITIES =
    lowercase(get(ENV, "MYCELIA_COVERAGE_SKIP_AQUA_AMBIGUITIES", "true")) == "true"
const COVERAGE_SKIP_JET =
    lowercase(get(ENV, "MYCELIA_COVERAGE_SKIP_JET", "true")) == "true"
const DEFAULT_COVERAGE_SKIP_TEST_FILES = [
    "hpc_assembly_wrappers_test.jl",
    "megahit_phix_workflow.jl",
    "checkm_tools.jl",
    "coverm_wrappers.jl",
    "coverm_integration_extended.jl",
    "mosdepth_coverage_qc.jl",
    "coverage_taxonomy_integration.jl",
    "quast_busco_wrappers_test.jl",
    "third_party_assemblers*.jl",
]
const COVERAGE_SKIP_TEST_FILES = let raw_skip_files = split(
        get(
            ENV,
            "MYCELIA_COVERAGE_SKIP_TEST_FILES",
            join(DEFAULT_COVERAGE_SKIP_TEST_FILES, ',')
        ),
        ','
    )
    sort!(
        unique!(
            filter!(
                !isempty,
                map(file -> strip(file), raw_skip_files)
            )
        )
    )
end
const COVERAGE_SKIP_HPC_ASSEMBLY_WRAPPERS =
    "hpc_assembly_wrappers_test.jl" in COVERAGE_SKIP_TEST_FILES
const COVERAGE_MODE = let coverage_mode = get(ENV, "MYCELIA_COVERAGE_MODE", "package")
    normalized_mode = lowercase(coverage_mode)
    if normalized_mode in ("package", "project", "default", "true")
        PACKAGE_COVERAGE_SCOPE
    elseif normalized_mode in ("none", "false")
        false
    elseif normalized_mode in ("user", "all") || startswith(coverage_mode, "@")
        coverage_mode
    else
        error(
            "Unsupported MYCELIA_COVERAGE_MODE=$(repr(coverage_mode)); use package, user, all, none, or an explicit @path"
        )
    end
end
const COVERAGE_MODE_LABEL = if COVERAGE_MODE === PACKAGE_COVERAGE_SCOPE
    "package ($(PACKAGE_COVERAGE_SCOPE))"
elseif COVERAGE_MODE isa Bool
    "none"
else
    COVERAGE_MODE
end

const TIER_RULES = Dict{String, Tuple{Int, String, Int}}(
    "reference-databases.jl" => (1, "Data Acquisition", 90),
    "simulation.jl" => (1, "Data Acquisition", 90),
    "ncbi-datasets-cli.jl" => (1, "Data Acquisition", 90),
    "busco-datasets.jl" => (1, "Data Acquisition", 90),
    "protein-databases.jl" => (1, "Data Acquisition", 90),

    "alphabets.jl" => (2, "Preprocessing/QC", 90),
    "constants.jl" => (2, "Preprocessing/QC", 90),
    "fastx.jl" => (2, "Preprocessing/QC", 90),
    "read-quality-control.jl" => (2, "Preprocessing/QC", 90),
    "xam.jl" => (2, "Preprocessing/QC", 90),

    "kmer-analysis.jl" => (3, "K-mer Analysis", 95),
    "kmer-saturation-analysis.jl" => (3, "K-mer Analysis", 95),
    "qualmer-analysis.jl" => (3, "K-mer Analysis", 95),
    "distance-metrics.jl" => (3, "K-mer Analysis", 95),
    "coverage-clustering.jl" => (3, "K-mer Analysis", 95),

    "assembly.jl" => (4, "Assembly Core", 95),
    "autocycler.jl" => (4, "Assembly Core", 95),
    "bcalm.jl" => (4, "Assembly Core", 95),
    "ggcat.jl" => (4, "Assembly Core", 95),
    "graph-cleanup.jl" => (4, "Assembly Core", 95),
    "iterative-assembly.jl" => (4, "Assembly Core", 95),
    "viterbi-next.jl" => (4, "Assembly Core", 95),

    "quality-control-and-benchmarking.jl" => (5, "Validation", 80),
    "variant-analysis.jl" => (5, "Validation", 80),
    "viterbi-polishing-and-error-correction.jl" => (5, "Validation", 80),

    "annotation.jl" => (6, "Annotation", 80),
    "genome-features.jl" => (6, "Annotation", 80),
    "codon-optimization.jl" => (6, "Annotation", 80),

    "pangenome-analysis.jl" => (7, "Comparative", 85),
    "sequence-comparison.jl" => (7, "Comparative", 85),
    "taxonomy-and-trees.jl" => (7, "Comparative", 85),
    "metagraph.jl" => (7, "Comparative", 85),
    "relational-matrices.jl" => (7, "Comparative", 85),

    "alignments-and-mapping.jl" => (8, "Tool Integration", 75),
    "binning.jl" => (8, "Tool Integration", 75),
    "bioconda.jl" => (8, "Tool Integration", 75),
    "classification.jl" => (8, "Tool Integration", 75),
    "foldseek.jl" => (8, "Tool Integration", 75),
    "metagenomic-classification.jl" => (8, "Tool Integration", 75),
    "pantools.jl" => (8, "Tool Integration", 75),
    "prokrustean.jl" => (8, "Tool Integration", 75),
    "rclone.jl" => (8, "Tool Integration", 75),
    "sentencepiece.jl" => (8, "Tool Integration", 75),

    "amino-acid-analysis.jl" => (9, "Analysis/Viz", 85),
    "clustering.jl" => (9, "Analysis/Viz", 85),
    "dimensionality-reduction.jl" => (9, "Analysis/Viz", 85),
    "plotting-and-visualization.jl" => (9, "Analysis/Viz", 85),
    "tda.jl" => (9, "Analysis/Viz", 85),

    "Mycelia.jl" => (10, "Infrastructure", 90),
    "checkpointing.jl" => (10, "Infrastructure", 90),
    "execution.jl" => (10, "Infrastructure", 90),
    "precompile_workload.jl" => (10, "Infrastructure", 90),
    "slurm-sbatch.jl" => (10, "Infrastructure", 90),
    "slurm-templates.jl" => (10, "Infrastructure", 90),
    "testing-utilities.jl" => (10, "Infrastructure", 90),
    "utility-functions.jl" => (10, "Infrastructure", 90)
)

const EXCLUSIONS = Dict{String, String}(
    "src/Mycelia.jl" => "Module loader and include graph root",
    "src/bioconda.jl" => "External environment management and installation side effects",
    "src/execution.jl" => "Shell execution wrappers depend on the host environment",
    "src/precompile_workload.jl" => "Precompile workload definitions rather than runtime logic",
    "src/slurm-templates.jl" => "Template generation for HPC environments",
    "src/testing-utilities.jl" => "Test helper infrastructure exercised transitively"
)

function parse_includes(path::String, seen::Set{String})::Vector{String}
    absolute_path = normpath(abspath(path))
    absolute_path in seen && return String[]
    push!(seen, absolute_path)

    included_paths = String[absolute_path]
    include_pattern = r"^include\(\"([^\"]+\.jl)\"\)"

    for line in eachline(absolute_path)
        stripped = strip(line)
        startswith(stripped, '#') && continue
        match_result = match(include_pattern, stripped)
        match_result === nothing && continue
        child_path = normpath(joinpath(dirname(absolute_path), only(match_result.captures)))
        append!(included_paths, parse_includes(child_path, seen))
    end

    return included_paths
end

function parse_direct_includes(path::String)::Vector{String}
    include_pattern = r"^include\(\"([^\"]+\.jl)\"\)"
    include_paths = String[]

    for line in eachline(path)
        stripped = strip(line)
        startswith(stripped, '#') && continue
        match_result = match(include_pattern, stripped)
        match_result === nothing && continue
        push!(
            include_paths,
            normpath(joinpath(dirname(path), only(match_result.captures)))
        )
    end

    return include_paths
end

function collect_source_scope()::NamedTuple{(:module_paths, :loaded_files), Tuple{Vector{String}, Vector{String}}}
    module_roots = parse_direct_includes(ROOT_ENTRYPOINT)
    module_paths = vcat([ROOT_ENTRYPOINT], module_roots)
    loaded_files = String[]

    for module_root in module_roots
        append!(loaded_files, parse_includes(module_root, Set{String}()))
    end

    return (
        module_paths = sort!(unique(module_paths)),
        loaded_files = sort!(unique(vcat([ROOT_ENTRYPOINT], loaded_files)))
    )
end

function relative_to_root(path::String)::String
    return relpath(path, ROOT_DIR)
end

function coverage_artifact_count(path::String)::Int
    folder = dirname(path)
    prefix = string(path, '.')
    return count(file -> startswith(file, prefix) && endswith(file, ".cov"), readdir(folder; join = true))
end

function total_coverage_artifact_count(paths::Vector{String})::Int
    return sum(coverage_artifact_count, paths)
end

function tier_for(rel_path::String)::Tuple{Int, String, Int}
    startswith(rel_path, "src/rhizomorph/") && return (4, "Assembly Core", 95)
    return get(TIER_RULES, basename(rel_path), (0, "Unassigned", 0))
end

function exclusion_reason(rel_path::String)::Union{Nothing, String}
    return get(EXCLUSIONS, rel_path, nothing)
end

function coverage_percent(lines_hit::Int, lines_total::Int)::Float64
    if lines_total == 0
        return 0.0
    end
    return round(100 * lines_hit / lines_total; digits = 1)
end

function coverage_summary(paths::Vector{String})::Tuple{Int, Int}
    lines_hit = 0
    lines_total = 0

    for path in paths
        file_coverage = Coverage.process_file(path)
        file_lines_hit, file_lines_total = Coverage.get_summary(file_coverage)
        lines_hit += file_lines_hit
        lines_total += file_lines_total
    end

    return (lines_hit, lines_total)
end

function module_record(path::String)::Dict{String, Any}
    rel_path = relative_to_root(path)
    source_files = if path == ROOT_ENTRYPOINT
        [path]
    else
        sort!(unique(parse_includes(path, Set{String}())))
    end
    lines_hit, lines_total = coverage_summary(source_files)
    tier_number, tier_name, target_pct = tier_for(rel_path)
    exclusion = exclusion_reason(rel_path)

    return Dict{String, Any}(
        "file" => rel_path,
        "filename" => basename(rel_path),
        "lines_total" => lines_total,
        "lines_hit" => lines_hit,
        "lines_missed" => lines_total - lines_hit,
        "coverage_pct" => coverage_percent(lines_hit, lines_total),
        "tier" => tier_number,
        "tier_name" => tier_name,
        "target_pct" => target_pct,
        "excluded" => exclusion !== nothing,
        "exclusion_reason" => exclusion,
        "cov_file_count" => total_coverage_artifact_count(source_files),
        "source_file_count" => length(source_files),
        "source_files" => map(relative_to_root, source_files)
    )
end

function aggregate_tiers(modules::Vector{Dict{String, Any}})::Vector{Dict{String, Any}}
    tier_stats = Dict{Int, Dict{String, Any}}()

    for module_record in modules
        module_record["excluded"] && continue
        tier_number = module_record["tier"]
        tier_number == 0 && continue

        if !haskey(tier_stats, tier_number)
            tier_stats[tier_number] = Dict{String, Any}(
                "tier" => tier_number,
                "name" => module_record["tier_name"],
                "target_pct" => module_record["target_pct"],
                "lines_total" => 0,
                "lines_hit" => 0,
                "modules" => String[]
            )
        end

        tier_record = tier_stats[tier_number]
        tier_record["lines_total"] += module_record["lines_total"]
        tier_record["lines_hit"] += module_record["lines_hit"]
        push!(tier_record["modules"], module_record["filename"])
    end

    tiers = Dict{String, Any}[]
    for tier_number in sort!(collect(keys(tier_stats)))
        tier_record = tier_stats[tier_number]
        lines_hit = tier_record["lines_hit"]
        lines_total = tier_record["lines_total"]
        tier_record["lines_missed"] = lines_total - lines_hit
        tier_record["coverage_pct"] = coverage_percent(lines_hit, lines_total)
        tier_record["gap_pct"] = round(tier_record["target_pct"] - tier_record["coverage_pct"]; digits = 1)
        push!(tiers, tier_record)
    end

    return tiers
end

function write_json(path::String, value)
    open(path, "w") do io
        JSON.print(io, value, 2)
    end
end

function write_summary_markdown(path::String, baseline::Dict{String, Any})
    modules = baseline["modules"]
    tiers = baseline["tiers"]
    overall = baseline["overall"]
    source_scope = baseline["source_scope"]

    open(path, "w") do io
        println(io, "# Mycelia Coverage Baseline")
        println(io, "")
        println(io, "Generated: $(baseline["generated_at"])")
        println(io, "")
        println(
            io,
            "Scope: $(source_scope["module_count"]) top-level modules spanning $(source_scope["loaded_file_count"]) loaded source files"
        )
        println(io, "")
        println(io, "## Overall: $(overall["coverage_pct"])% ($(overall["lines_hit"])/$(overall["lines_total"]) executable lines)")
        println(io, "")
        println(io, "## Tier Summary")
        println(io, "")
        println(io, "| Tier | Name | Coverage | Target | Gap | Lines | Modules |")
        println(io, "|------|------|----------|--------|-----|-------|---------|")

        for tier_record in tiers
            println(io,
                "| $(tier_record["tier"]) | $(tier_record["name"]) | $(tier_record["coverage_pct"])% | $(tier_record["target_pct"])% | $(tier_record["gap_pct"])pp | $(tier_record["lines_hit"])/$(tier_record["lines_total"]) | $(length(tier_record["modules"])) |")
        end

        println(io, "")
        println(io, "## Module Detail")
        println(io, "")
        println(io, "| Module | Tier | Coverage | Target | Lines Hit | Lines Total | Source Files | Excluded |")
        println(io, "|--------|------|----------|--------|-----------|-------------|--------------|----------|")

        for module_record in modules
            excluded = module_record["excluded"] ? "yes" : ""
            println(io,
                "| $(module_record["file"]) | $(module_record["tier_name"]) | $(module_record["coverage_pct"])% | $(module_record["target_pct"])% | $(module_record["lines_hit"]) | $(module_record["lines_total"]) | $(module_record["source_file_count"]) | $(excluded) |")
        end
    end
end

function baseline_payload(
        modules::Vector{Dict{String, Any}},
        module_paths::Vector{String},
        loaded_files::Vector{String}
    )::Dict{String, Any}
    tiers = aggregate_tiers(modules)

    in_scope_modules = filter(module_record -> !module_record["excluded"], modules)
    lines_total = sum(module_record["lines_total"] for module_record in in_scope_modules)
    lines_hit = sum(module_record["lines_hit"] for module_record in in_scope_modules)

    return Dict{String, Any}(
        "generated_at" => Dates.format(Dates.now(), Dates.DateFormat("yyyy-mm-ddTHH:MM:SS.sss")),
        "source_scope" => Dict{String, Any}(
            "entrypoint" => relative_to_root(ROOT_ENTRYPOINT),
            "module_count" => length(module_paths),
            "modules" => map(relative_to_root, module_paths),
            "loaded_file_count" => length(loaded_files),
            "loaded_files" => map(relative_to_root, loaded_files)
        ),
        "overall" => Dict{String, Any}(
            "lines_total" => lines_total,
            "lines_hit" => lines_hit,
            "lines_missed" => lines_total - lines_hit,
            "coverage_pct" => coverage_percent(lines_hit, lines_total)
        ),
        "tiers" => tiers,
        "modules" => modules
    )
end

function exclusions_payload(modules::Vector{Dict{String, Any}})::Dict{String, Any}
    excluded_modules = filter(module_record -> module_record["excluded"], modules)

    return Dict{String, Any}(
        "policy" => "Modules excluded from baseline targets but still measured within the package include graph",
        "exclusions" => [
            Dict{String, Any}(
                "file" => module_record["file"],
                "reason" => module_record["exclusion_reason"]
            ) for module_record in excluded_modules
        ]
    )
end

function clear_existing_coverage_artifacts()::Int
    removed = 0
    for coverage_root in (SRC_DIR, TEST_DIR)
        for (root, _, files) in walkdir(coverage_root)
            for file in files
                endswith(file, ".cov") || continue
                rm(joinpath(root, file); force = true)
                removed += 1
            end
        end
    end
    return removed
end

function coverage_artifact_paths()::Vector{String}
    artifact_paths = String[]
    for (root, _, files) in walkdir(SRC_DIR)
        for file in files
            endswith(file, ".cov") || continue
            push!(artifact_paths, joinpath(root, file))
        end
    end
    return sort!(artifact_paths)
end

function instantiate_coverage_project()
    command = `$(Base.julia_cmd()) --project=$(COVERAGE_DIR) -e "import Pkg; Pkg.instantiate()"`
    run(command)
end

function warm_coverage_project_precompile()
    command = `$(Base.julia_cmd()) --project=$(COVERAGE_DIR) -e "import Pkg; Pkg.precompile()"`
    run(command)
end

function coverage_collection_julia_args()::Vector{String}
    julia_args = String[]
    if COVERAGE_NO_COMPILED_MODULES
        push!(julia_args, "--compiled-modules=no")
    end
    if COVERAGE_NO_PKGIMAGES
        push!(julia_args, "--pkgimages=no")
    end
    if COVERAGE_MODE !== false
        push!(julia_args, "--code-coverage=$(COVERAGE_MODE)")
    end
    return julia_args
end

function coverage_test_command()::Cmd
    julia_command = [
        Base.julia_cmd().exec...,
        "--project=$(COVERAGE_DIR)",
        coverage_collection_julia_args()...,
        ROOT_RUNTESTS,
    ]
    return Cmd(julia_command)
end

function collect_coverage_artifacts()
    removed = clear_existing_coverage_artifacts()
    println("Removed $(removed) stale coverage artifact(s)")
    println("Coverage mode: run_all=$(COVERAGE_RUN_ALL), run_external=$(COVERAGE_RUN_EXTERNAL)")
    println("Coverage mode: code_coverage=$(COVERAGE_MODE_LABEL)")
    println("Coverage startup: pkg_test_subprocess_julia_args=$(join(coverage_collection_julia_args(), ' '))")
    println("Coverage startup: aqua=$(COVERAGE_SKIP_AQUA ? "skipped" : "enabled")")
    println("Coverage startup: aqua_ambiguities=$(COVERAGE_SKIP_AQUA_AMBIGUITIES ? "skipped" : "enabled")")
    println("Coverage startup: jet=$(COVERAGE_SKIP_JET ? "skipped" : "enabled")")
    println(
        "Coverage startup: skip_test_files=$(isempty(COVERAGE_SKIP_TEST_FILES) ? "none" : join(COVERAGE_SKIP_TEST_FILES, ", "))"
    )

    println("Ensuring coverage runner environment is instantiated")
    instantiate_coverage_project()

    if COVERAGE_WARM_PRECOMPILE
        println("Warming package precompile cache in the coverage runner project")
        warm_coverage_project_precompile()
    end

    environment = [
        "MYCELIA_RUN_ALL" => string(COVERAGE_RUN_ALL),
        "MYCELIA_RUN_EXTERNAL" => string(COVERAGE_RUN_EXTERNAL),
        "MYCELIA_SKIP_AQUA" => string(COVERAGE_SKIP_AQUA),
        "MYCELIA_SKIP_AQUA_AMBIGUITIES" => string(COVERAGE_SKIP_AQUA_AMBIGUITIES),
        "MYCELIA_SKIP_JET" => string(COVERAGE_SKIP_JET),
        "MYCELIA_SKIP_TEST_FILES" => join(COVERAGE_SKIP_TEST_FILES, ','),
        "MYCELIA_SKIP_HPC_ASSEMBLY_WRAPPERS" => string(COVERAGE_SKIP_HPC_ASSEMBLY_WRAPPERS),
        "JULIA_PKG_PRECOMPILE_AUTO" => (COVERAGE_NO_COMPILED_MODULES ? "0" : "1"),
    ]

    println("Running full test suite via test/runtests.jl under coverage/Project.toml")
    try
        run(setenv(coverage_test_command(), environment))
    catch error
        throw(ErrorException("Coverage collection failed: $(error)"))
    end

    generated = length(coverage_artifact_paths())
    println("Collected $(generated) coverage artifact(s)")
end

function main()
    source_scope = collect_source_scope()
    if COVERAGE_COLLECT
        collect_coverage_artifacts()
    else
        println("Skipping coverage collection; processing existing artifacts only")
    end

    artifact_count = total_coverage_artifact_count(source_scope.loaded_files)
    artifact_count == 0 && error(
        "No .cov files were found for the package include graph. Run the coverage suite first."
    )

    println(
        "Processing coverage for $(length(source_scope.module_paths)) top-level modules across $(length(source_scope.loaded_files)) loaded source files"
    )
    println("Found $(artifact_count) coverage artifact(s)")

    modules = sort!(
        map(module_record, source_scope.module_paths);
        by = module_record -> (module_record["tier"], module_record["file"])
    )

    baseline = baseline_payload(modules, source_scope.module_paths, source_scope.loaded_files)
    exclusions = exclusions_payload(modules)

    baseline_path = joinpath(COVERAGE_DIR, "baseline.json")
    summary_path = joinpath(COVERAGE_DIR, "baseline-summary.md")
    exclusions_path = joinpath(COVERAGE_DIR, "exclusions.json")

    write_json(baseline_path, baseline)
    write_summary_markdown(summary_path, baseline)
    write_json(exclusions_path, exclusions)

    overall = baseline["overall"]
    println("Wrote $(baseline_path)")
    println("Wrote $(summary_path)")
    println("Wrote $(exclusions_path)")
    println("Overall coverage: $(overall["coverage_pct"])% ($(overall["lines_hit"])/$(overall["lines_total"]))")
end

main()
