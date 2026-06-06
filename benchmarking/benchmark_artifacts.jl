import CSV
import DataFrames
import Dates
import JSON
import OrderedCollections
import TOML

const BENCHMARK_ARTIFACT_SCHEMA_VERSION = 1
const BENCHMARK_ARTIFACT_SUBDIRS = ("tables", "plots", "logs", "provenance")

"""
    benchmark_artifact_layout(output_dir; create=true)

Return the public benchmark artifact directory layout.

The layout is stable across CI-sized and full benchmark runs:
`tables/`, `plots/`, `logs/`, and `provenance/`.
"""
function benchmark_artifact_layout(output_dir::AbstractString; create::Bool = true)
    root = abspath(output_dir)
    paths = Dict{Symbol, String}(:root => root)
    if create
        mkpath(root)
    end

    for subdir in BENCHMARK_ARTIFACT_SUBDIRS
        path = joinpath(root, subdir)
        paths[Symbol(subdir)] = path
        if create
            mkpath(path)
        end
    end

    return (
        root = paths[:root],
        tables = paths[:tables],
        plots = paths[:plots],
        logs = paths[:logs],
        provenance = paths[:provenance]
    )
end

"""
    collect_benchmark_provenance(; kwargs...)

Collect reproducibility metadata for public benchmark artifacts.

The returned dictionary records schema version, run identity, command arguments,
dataset identifiers, Julia/Mycelia versions, git metadata when available, and
caller-supplied metadata.
"""
function collect_benchmark_provenance(;
        run_id::AbstractString,
        scale::AbstractString,
        dataset_ids = String[],
        command_args = String[],
        tool_versions = Dict{String, Any}(),
        metadata = Dict{String, Any}(),
        generated_at = nothing,
        repo_root::AbstractString = abspath(joinpath(@__DIR__, "..")))
    project_versions = _benchmark_project_versions(repo_root)
    merged_tool_versions = Dict{String, Any}("julia" => string(VERSION))
    for (key, value) in project_versions
        merged_tool_versions[string(key)] = value
    end
    for (key, value) in tool_versions
        merged_tool_versions[string(key)] = value
    end

    return OrderedCollections.OrderedDict{String, Any}(
        "schema_version" => BENCHMARK_ARTIFACT_SCHEMA_VERSION,
        "run_id" => string(run_id),
        "scale" => string(scale),
        "generated_at" => isnothing(generated_at) ? _utc_timestamp() : string(generated_at),
        "dataset_ids" => sort(string.(dataset_ids)),
        "command_args" => string.(command_args),
        "tool_versions" => merged_tool_versions,
        "git" => _benchmark_git_metadata(repo_root),
        "metadata" => Dict(string(key) => value for (key, value) in metadata)
    )
end

"""
    write_benchmark_table(table, table_name; layout, provenance, row_context_columns=Dict())

Write a CSV table under `layout.tables` and a sidecar provenance JSON file
under `layout.provenance`.

Stable `benchmark_*` provenance columns are prepended to the table before
writing. `row_context_columns` maps provenance column names to source columns
already present in the table, for example
`"benchmark_dataset_id" => "dataset_id"`.
"""
function write_benchmark_table(
        table,
        table_name::AbstractString;
        layout,
        provenance::AbstractDict,
        row_context_columns = Dict{String, String}())
    table_basename = _artifact_basename(table_name)
    output_table = DataFrames.DataFrame(table)
    _insert_table_provenance_columns!(output_table, provenance, row_context_columns)

    table_path = joinpath(layout.tables, "$(table_basename).csv")
    CSV.write(table_path, output_table)

    table_provenance = copy(Dict{String, Any}(provenance))
    table_provenance["table_name"] = table_basename
    table_provenance["table_path"] = relpath(table_path, layout.root)
    table_provenance["row_count"] = DataFrames.nrow(output_table)
    table_provenance["columns"] = names(output_table)

    provenance_path = joinpath(layout.provenance, "$(table_basename).provenance.json")
    _write_benchmark_json(provenance_path, table_provenance)

    return (
        table = table_path,
        provenance = provenance_path,
        rows = DataFrames.nrow(output_table),
        columns = names(output_table)
    )
end

"""
    write_benchmark_artifacts(tables; output_dir, run_id, scale, kwargs...)

Write public benchmark tables plus run-level provenance.

`tables` is an iterable of `name => table` pairs. The function returns paths to
the root layout, run provenance, and each table artifact.
"""
function write_benchmark_artifacts(
        tables;
        output_dir::AbstractString,
        run_id::AbstractString,
        scale::AbstractString,
        dataset_ids = String[],
        command_args = String[],
        tool_versions = Dict{String, Any}(),
        metadata = Dict{String, Any}(),
        generated_at = nothing,
        table_context_columns = Dict{String, Dict{String, String}}())
    layout = benchmark_artifact_layout(output_dir)
    provenance = collect_benchmark_provenance(
        run_id = run_id,
        scale = scale,
        dataset_ids = dataset_ids,
        command_args = command_args,
        tool_versions = tool_versions,
        metadata = metadata,
        generated_at = generated_at
    )

    run_provenance_path = joinpath(layout.provenance, "run.provenance.json")
    _write_benchmark_json(run_provenance_path, provenance)

    table_paths = OrderedCollections.OrderedDict{String, Any}()
    normalized_tables = sort(
        [(string(table_name), table) for (table_name, table) in tables];
        by = first
    )
    for (normalized_name, table) in normalized_tables
        row_context_columns = get(table_context_columns, normalized_name, Dict{String, String}())
        table_paths[normalized_name] = write_benchmark_table(
            table,
            normalized_name;
            layout = layout,
            provenance = provenance,
            row_context_columns = row_context_columns
        )
    end

    index = OrderedCollections.OrderedDict{String, Any}(
        "schema_version" => BENCHMARK_ARTIFACT_SCHEMA_VERSION,
        "run_id" => string(run_id),
        "scale" => string(scale),
        "directories" => OrderedCollections.OrderedDict{String, Any}(
            "tables" => relpath(layout.tables, layout.root),
            "plots" => relpath(layout.plots, layout.root),
            "logs" => relpath(layout.logs, layout.root),
            "provenance" => relpath(layout.provenance, layout.root)
        ),
        "tables" => OrderedCollections.OrderedDict{String, Any}(
            name => Dict(
                "table" => relpath(paths.table, layout.root),
                "provenance" => relpath(paths.provenance, layout.root),
                "rows" => paths.rows,
                "columns" => paths.columns
            )
            for (name, paths) in table_paths
        )
    )
    index_path = joinpath(layout.root, "artifact-index.json")
    _write_benchmark_json(index_path, index)

    return (
        root = layout.root,
        tables = table_paths,
        provenance = run_provenance_path,
        index = index_path,
        layout = layout
    )
end

function _insert_table_provenance_columns!(
        table::DataFrames.DataFrame,
        provenance::AbstractDict,
        row_context_columns)
    rows = DataFrames.nrow(table)
    git = get(provenance, "git", Dict{String, Any}())
    dataset_ids = get(provenance, "dataset_ids", String[])
    column_values = Pair{String, Any}[
        "benchmark_schema_version" => get(provenance, "schema_version", BENCHMARK_ARTIFACT_SCHEMA_VERSION),
        "benchmark_run_id" => get(provenance, "run_id", ""),
        "benchmark_scale" => get(provenance, "scale", ""),
        "benchmark_dataset_ids" => join(string.(dataset_ids), ","),
        "benchmark_git_commit" => _string_or_empty(get(git, "commit", nothing))
    ]

    for (provenance_column, source_column) in sort(collect(row_context_columns); by = first)
        if source_column in names(table)
            push!(column_values, string(provenance_column) => string.(table[!, source_column]))
        else
            push!(column_values, string(provenance_column) => "")
        end
    end

    for (column_name, value) in reverse(column_values)
        if column_name in names(table)
            continue
        end
        values = value isa AbstractVector ? value : fill(value, rows)
        DataFrames.insertcols!(table, 1, Symbol(column_name) => values)
    end

    return table
end

function _benchmark_project_versions(repo_root::AbstractString)
    project_toml = joinpath(repo_root, "Project.toml")
    if !isfile(project_toml)
        return Dict{String, Any}()
    end

    project = TOML.parsefile(project_toml)
    package_name = string(get(project, "name", "project"))
    package_version = string(get(project, "version", "unknown"))
    return Dict{String, Any}(package_name => package_version)
end

function _benchmark_git_metadata(repo_root::AbstractString)
    return OrderedCollections.OrderedDict{String, Any}(
        "commit" => _git_output(repo_root, ["rev-parse", "HEAD"]),
        "branch" => _git_output(repo_root, ["rev-parse", "--abbrev-ref", "HEAD"]),
        "remote_url" => _git_output(repo_root, ["config", "--get", "remote.origin.url"]),
        "is_dirty" => _git_dirty(repo_root)
    )
end

function _git_dirty(repo_root::AbstractString)
    status = _git_output(repo_root, ["status", "--porcelain"])
    return isnothing(status) ? nothing : !isempty(status)
end

function _git_output(repo_root::AbstractString, args::Vector{String})
    command = Cmd(vcat(["git", "-C", repo_root], args))
    try
        return chomp(read(command, String))
    catch
        return nothing
    end
end

function _write_benchmark_json(path::AbstractString, data)
    open(path, "w") do io
        JSON.print(io, _ordered_benchmark_json(data), 2)
        println(io)
    end
    return path
end

function _ordered_benchmark_json(data::AbstractDict)
    ordered = OrderedCollections.OrderedDict{String, Any}()
    entries = sort([(string(key), value) for (key, value) in data]; by = first)
    for (key, value) in entries
        ordered[key] = _ordered_benchmark_json(value)
    end
    return ordered
end

function _ordered_benchmark_json(data::AbstractVector)
    return [_ordered_benchmark_json(value) for value in data]
end

function _ordered_benchmark_json(data)
    return data
end

function _artifact_basename(name::AbstractString)
    sanitized = replace(lowercase(string(name)), r"[^a-z0-9_.-]+" => "_")
    sanitized = strip(sanitized, ['_', '.', '-'])
    if isempty(sanitized)
        error("Artifact name must contain at least one alphanumeric character.")
    end
    return sanitized
end

function _string_or_empty(value)
    return isnothing(value) ? "" : string(value)
end

function _utc_timestamp()
    return string(Dates.now(Dates.UTC)) * "Z"
end
