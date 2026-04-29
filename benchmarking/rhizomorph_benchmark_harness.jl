#!/usr/bin/env julia

import DataFrames
import TOML

const RHIZOMORPH_BENCHMARK_MANIFEST =
    joinpath(@__DIR__, "rhizomorph_benchmark_manifest.toml")

const RHIZOMORPH_BENCHMARK_SCALE_RANK = Dict(
    "ci" => 1,
    "full" => 2,
    "candidate" => 3
)

"""
    load_rhizomorph_benchmark_manifest(manifest_path=RHIZOMORPH_BENCHMARK_MANIFEST; validate=true)

Load the Rhizomorph benchmark dataset and H1-H7 slice manifest.
"""
function load_rhizomorph_benchmark_manifest(
        manifest_path::AbstractString = RHIZOMORPH_BENCHMARK_MANIFEST;
        validate::Bool = true)
    manifest = TOML.parsefile(manifest_path)
    if validate
        validate_rhizomorph_benchmark_manifest(manifest)
    end
    return manifest
end

"""
    validate_rhizomorph_benchmark_manifest(manifest)

Validate required fields, unique identifiers, scale labels, and slice dataset references.
"""
function validate_rhizomorph_benchmark_manifest(manifest::AbstractDict)
    _require_keys(manifest, ["schema_version", "schema", "datasets", "hypothesis_slices"], "manifest")
    schema = manifest["schema"]
    _require_keys(schema, ["dataset_required_fields", "hypothesis_required_fields"], "schema")

    dataset_required_fields = schema["dataset_required_fields"]
    hypothesis_required_fields = schema["hypothesis_required_fields"]
    allowed_suitability = Set(get(schema, "ci_suitability_values", collect(keys(RHIZOMORPH_BENCHMARK_SCALE_RANK))))
    canonical_suitability = Set(keys(RHIZOMORPH_BENCHMARK_SCALE_RANK))
    unsupported_suitability = setdiff(allowed_suitability, canonical_suitability)
    if !isempty(unsupported_suitability)
        error("Unsupported ci_suitability_values in schema: $(join(sort(collect(unsupported_suitability)), ", "))")
    end

    dataset_ids = Set{String}()
    for dataset in manifest["datasets"]
        _require_keys(dataset, dataset_required_fields, "dataset")
        dataset_id = string(dataset["id"])
        if dataset_id in dataset_ids
            error("Duplicate Rhizomorph benchmark dataset id: $(dataset_id)")
        end
        push!(dataset_ids, dataset_id)

        suitability = string(dataset["ci_suitability"])
        if !(suitability in allowed_suitability)
            error("Dataset $(dataset_id) has invalid ci_suitability=$(suitability)")
        end
        _require_keys(dataset["provenance"], ["kind", "source", "accessions"], "dataset $(dataset_id) provenance")
    end

    slice_ids = Set{String}()
    for slice in manifest["hypothesis_slices"]
        _require_keys(slice, hypothesis_required_fields, "hypothesis slice")
        slice_id = string(slice["id"])
        if slice_id in slice_ids
            error("Duplicate Rhizomorph benchmark hypothesis id: $(slice_id)")
        end
        push!(slice_ids, slice_id)

        for dataset_id in slice["dataset_ids"]
            if !(string(dataset_id) in dataset_ids)
                error("Hypothesis slice $(slice_id) references unknown dataset: $(dataset_id)")
            end
        end
    end

    return true
end

"""
    list_rhizomorph_benchmark_datasets(; manifest_path=RHIZOMORPH_BENCHMARK_MANIFEST)

Return a table of benchmark datasets, provenance, expected outputs, and CI/full suitability.
"""
function list_rhizomorph_benchmark_datasets(;
        manifest_path::AbstractString = RHIZOMORPH_BENCHMARK_MANIFEST)
    manifest = load_rhizomorph_benchmark_manifest(manifest_path)
    rows = NamedTuple[]

    for dataset in manifest["datasets"]
        provenance = dataset["provenance"]
        push!(rows, (
            id = string(dataset["id"]),
            name = string(dataset["name"]),
            category = string(dataset["category"]),
            ci_suitability = string(dataset["ci_suitability"]),
            materialization = string(dataset["materialization"]),
            provenance_kind = string(provenance["kind"]),
            provenance_source = string(provenance["source"]),
            accessions = join(string.(provenance["accessions"]), ","),
            expected_outputs = join(string.(dataset["expected_outputs"]), "; ")
        ))
    end

    return DataFrames.DataFrame(rows)
end

"""
    list_rhizomorph_benchmark_slices(; manifest_path=RHIZOMORPH_BENCHMARK_MANIFEST)

Return a table describing the H1-H7 benchmark slices and their harness entry points.
"""
function list_rhizomorph_benchmark_slices(;
        manifest_path::AbstractString = RHIZOMORPH_BENCHMARK_MANIFEST)
    manifest = load_rhizomorph_benchmark_manifest(manifest_path)
    rows = NamedTuple[]

    for slice in manifest["hypothesis_slices"]
        push!(rows, (
            id = string(slice["id"]),
            title = string(slice["title"]),
            status = string(get(slice, "status", "stub")),
            dataset_ids = join(string.(slice["dataset_ids"]), ","),
            entrypoint = string(slice["entrypoint"]),
            expected_outputs = join(string.(slice["expected_outputs"]), "; ")
        ))
    end

    return DataFrames.DataFrame(rows)
end

"""
    build_rhizomorph_benchmark_plan(; scale="ci", hypothesis_ids=nothing, dataset_ids=nothing, manifest_path=...)

Build a dry-run plan that maps H1-H7 benchmark slices to datasets and expected outputs.

`scale="ci"` includes only CI-suitable datasets, `scale="full"` includes CI and
full public-reference datasets, and `scale="candidate"` also includes planned
heterogeneous candidates.
"""
function build_rhizomorph_benchmark_plan(;
        scale::AbstractString = "ci",
        hypothesis_ids = nothing,
        dataset_ids = nothing,
        manifest_path::AbstractString = RHIZOMORPH_BENCHMARK_MANIFEST)
    _validate_scale(scale)
    manifest = load_rhizomorph_benchmark_manifest(manifest_path)
    dataset_by_id = Dict(string(dataset["id"]) => dataset for dataset in manifest["datasets"])
    selected_hypothesis_ids = _selected_id_set(hypothesis_ids)
    selected_dataset_ids = _selected_id_set(dataset_ids)
    rows = NamedTuple[]

    for slice in manifest["hypothesis_slices"]
        slice_id = string(slice["id"])
        if selected_hypothesis_ids !== nothing && !(slice_id in selected_hypothesis_ids)
            continue
        end

        for dataset_id_value in slice["dataset_ids"]
            dataset_id = string(dataset_id_value)
            dataset = dataset_by_id[dataset_id]
            if !_dataset_allowed_for_scale(dataset, scale)
                continue
            end
            if selected_dataset_ids !== nothing && !(dataset_id in selected_dataset_ids)
                continue
            end

            push!(rows, (
                scale = string(scale),
                hypothesis_id = slice_id,
                hypothesis_title = string(slice["title"]),
                dataset_id = dataset_id,
                dataset_category = string(dataset["category"]),
                ci_suitability = string(dataset["ci_suitability"]),
                materialization = string(dataset["materialization"]),
                entrypoint = string(slice["entrypoint"]),
                expected_inputs = join(string.(slice["expected_inputs"]), "; "),
                expected_outputs = join(string.(slice["expected_outputs"]), "; "),
                implemented = string(get(slice, "status", "stub")) != "stub"
            ))
        end
    end

    _validate_requested_ids("hypothesis", selected_hypothesis_ids, string.(getindex.(manifest["hypothesis_slices"], "id")))
    _validate_requested_ids("dataset", selected_dataset_ids, collect(keys(dataset_by_id)))
    return DataFrames.DataFrame(rows)
end

"""
    run_rhizomorph_benchmark_harness(; dry_run=true, kwargs...)

Return the dry-run benchmark plan. Executing benchmark slices is intentionally a
stub until follow-on issues implement each H1-H7 runner.
"""
function run_rhizomorph_benchmark_harness(; dry_run::Bool = true, kwargs...)
    plan = build_rhizomorph_benchmark_plan(; kwargs...)
    if dry_run
        return plan
    end
    error("Rhizomorph benchmark slice execution is not implemented yet; use dry_run=true to inspect the plan.")
end

function _require_keys(table::AbstractDict, keys, context::AbstractString)
    missing_keys = String[]
    for key in keys
        if !haskey(table, key)
            push!(missing_keys, string(key))
        end
    end
    if !isempty(missing_keys)
        error("Missing required $(context) key(s): $(join(missing_keys, ", "))")
    end
    return nothing
end

function _validate_scale(scale::AbstractString)
    if !haskey(RHIZOMORPH_BENCHMARK_SCALE_RANK, string(scale))
        error("Unknown Rhizomorph benchmark scale: $(scale). Use ci, full, or candidate.")
    end
    return nothing
end

function _selected_id_set(ids)
    if ids === nothing
        return nothing
    end
    return Set(string.(ids))
end

function _dataset_allowed_for_scale(dataset::AbstractDict, scale::AbstractString)
    dataset_scale = string(dataset["ci_suitability"])
    return RHIZOMORPH_BENCHMARK_SCALE_RANK[dataset_scale] <=
           RHIZOMORPH_BENCHMARK_SCALE_RANK[string(scale)]
end

function _validate_requested_ids(kind::AbstractString, requested_ids, known_ids)
    if requested_ids === nothing
        return nothing
    end

    known_id_set = Set(string.(known_ids))
    unknown_ids = sort(collect(setdiff(requested_ids, known_id_set)))
    if !isempty(unknown_ids)
        error("Unknown Rhizomorph benchmark $(kind) id(s): $(join(unknown_ids, ", "))")
    end
    return nothing
end

function _flag_values(args, flag::AbstractString)
    values = String[]
    index = 1
    while index <= length(args)
        if args[index] == flag
            if index == length(args)
                error("Missing value after $(flag)")
            end
            if startswith(args[index + 1], "-")
                error("Missing value after $(flag)")
            end
            push!(values, args[index + 1])
            index += 2
        else
            index += 1
        end
    end
    return values
end

function _flag_value(args, flag::AbstractString, default::AbstractString)
    values = _flag_values(args, flag)
    return isempty(values) ? default : values[end]
end

function print_rhizomorph_benchmark_usage()
    println("Rhizomorph benchmark harness")
    println()
    println("Usage:")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --list-datasets")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --list-slices")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --plan --scale ci")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --slice H1 --slice H7 --scale full")
    println()
    println("Scales: ci, full, candidate")
    println("Execution is currently stubbed; the script emits dry-run plans for follow-on runners.")
    return nothing
end

function main(args = ARGS)
    known_flags = Set([
        "--help",
        "-h",
        "--list-datasets",
        "--list-slices",
        "--plan",
        "--scale",
        "--slice",
        "--dataset",
        "--manifest",
        "--execute"
    ])
    for arg in args
        if startswith(arg, "-") && !(arg in known_flags)
            error("Unknown flag: $(arg)")
        end
    end

    if "--help" in args || "-h" in args
        print_rhizomorph_benchmark_usage()
        return nothing
    end

    manifest_path = _flag_value(args, "--manifest", RHIZOMORPH_BENCHMARK_MANIFEST)
    if "--list-datasets" in args
        show(list_rhizomorph_benchmark_datasets(manifest_path = manifest_path); allrows = true, allcols = true)
        println()
        return nothing
    end
    if "--list-slices" in args
        show(list_rhizomorph_benchmark_slices(manifest_path = manifest_path); allrows = true, allcols = true)
        println()
        return nothing
    end

    scale = _flag_value(args, "--scale", "ci")
    hypothesis_ids = _flag_values(args, "--slice")
    dataset_ids = _flag_values(args, "--dataset")
    dry_run = !("--execute" in args)
    plan = run_rhizomorph_benchmark_harness(
        dry_run = dry_run,
        scale = scale,
        hypothesis_ids = isempty(hypothesis_ids) ? nothing : hypothesis_ids,
        dataset_ids = isempty(dataset_ids) ? nothing : dataset_ids,
        manifest_path = manifest_path
    )
    show(plan; allrows = true, allcols = true)
    println()
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
