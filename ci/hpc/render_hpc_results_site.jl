#!/usr/bin/env julia

import Dates
import JSON

function load_json(path::String)
    return JSON.parsefile(path)
end

function ensure_dict(value)
    return value isa AbstractDict ? Dict{String, Any}(string(k) => v for (k, v) in value) : Dict{String, Any}()
end

function parse_duration(value)::Int
    if value isa Integer
        return Int(value)
    elseif value isa AbstractFloat
        return round(Int, value)
    elseif value isa AbstractString
        try
            return parse(Int, value)
        catch
            return 0
        end
    else
        return 0
    end
end

function endpoint_color(status::AbstractString)::String
    normalized = lowercase(strip(status))
    if normalized in ("pass", "passed", "success", "ok")
        return "brightgreen"
    elseif normalized in ("skipped", "not run", "unknown")
        return "lightgrey"
    else
        return "red"
    end
end

function status_message(section::Dict{String, Any})::String
    status = replace(lowercase(string(get(section, "status", "unknown"))), "_" => " ")
    duration = parse_duration(get(section, "duration_seconds", 0))
    return duration > 0 ? "$(status) $(duration)s" : status
end

function write_json(path::String, object)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, object, 4)
        write(io, '\n')
    end
end

function coverage_string(coverage)
    if coverage === nothing || !(coverage isa AbstractDict)
        return "n/a"
    end
    percent = get(coverage, "coverage_percent", nothing)
    if percent === nothing
        return "n/a"
    elseif percent isa Integer
        return string(percent) * "%"
    elseif percent isa AbstractFloat
        return string(round(percent; digits = 1)) * "%"
    else
        return string(percent)
    end
end

function build_badge(label::String, message::String, color::String)
    return Dict(
        "schemaVersion" => 1,
        "label" => label,
        "message" => message,
        "color" => color,
    )
end

function readme_text(result::Dict{String, Any}, repository_slug::String, branch_name::String)
    commit = string(get(result, "commit", "unknown"))
    source_branch = string(get(result, "branch", "unknown"))
    generated_at = string(get(result, "generated_at", "unknown"))
    hpc = ensure_dict(get(result, "hpc", Dict{String, Any}()))
    tests = ensure_dict(get(result, "tests", Dict{String, Any}()))
    benchmarks = ensure_dict(get(result, "benchmarks", Dict{String, Any}()))
    coverage = get(result, "coverage", nothing)

    raw_base = isempty(repository_slug) ? "" :
        "https://raw.githubusercontent.com/$(repository_slug)/$(branch_name)"
    tests_badge = isempty(raw_base) ? "" :
        "[![HPC Tests](https://img.shields.io/endpoint?url=$(raw_base)/latest-tests.json)]($(raw_base)/latest-hpc-results.json)"
    benchmarks_badge = isempty(raw_base) ? "" :
        "[![HPC Benchmarks](https://img.shields.io/endpoint?url=$(raw_base)/latest-benchmarks.json)]($(raw_base)/latest-hpc-results.json)"

    badge_block = isempty(tests_badge) ? "" : string(tests_badge, "\n", benchmarks_badge, "\n\n")
    cluster = string(get(hpc, "cluster", "unknown"))
    job_id = string(get(hpc, "job_id", ""))
    node = string(get(hpc, "node", ""))

    return """
# Mycelia HPC CI Results

$(badge_block)Latest published HPC CI summary for Mycelia.

| Field | Value |
| --- | --- |
| Commit | `$(commit)` |
| Source branch | `$(source_branch)` |
| Generated at (UTC) | `$(generated_at)` |
| Cluster | `$(cluster)` |
| Job ID | `$(isempty(job_id) ? "n/a" : job_id)` |
| Node | `$(isempty(node) ? "n/a" : node)` |
| Tests | `$(status_message(tests))` |
| Benchmarks | `$(status_message(benchmarks))` |
| Coverage | `$(coverage_string(coverage))` |

## Files

- [`latest-hpc-results.json`](latest-hpc-results.json)
- [`latest-tests.json`](latest-tests.json)
- [`latest-benchmarks.json`](latest-benchmarks.json)
- [`latest-meta.json`](latest-meta.json)
- [`$(commit)/hpc-results.json`]($(commit)/hpc-results.json)
- [`$(commit)/tests.json`]($(commit)/tests.json)
- [`$(commit)/benchmarks.json`]($(commit)/benchmarks.json)
- [`$(commit)/meta.json`]($(commit)/meta.json)
"""
end

function main()
    if length(ARGS) != 2
        error("Usage: render_hpc_results_site.jl <input-json> <output-dir>")
    end

    input_json = abspath(ARGS[1])
    output_dir = abspath(ARGS[2])
    result = ensure_dict(load_json(input_json))

    commit = string(get(result, "commit", "unknown"))
    branch_name = get(ENV, "HPC_RESULTS_BRANCH", "hpc-results")
    repository_slug = get(ENV, "HPC_RESULTS_REPOSITORY_SLUG", "")
    tests = ensure_dict(get(result, "tests", Dict{String, Any}()))
    benchmarks = ensure_dict(get(result, "benchmarks", Dict{String, Any}()))
    coverage = get(result, "coverage", nothing)
    hpc = ensure_dict(get(result, "hpc", Dict{String, Any}()))

    archive_dir = joinpath(output_dir, commit)
    mkpath(archive_dir)

    latest_result_path = joinpath(output_dir, "latest-hpc-results.json")
    archived_result_path = joinpath(archive_dir, "hpc-results.json")
    cp(input_json, latest_result_path; force = true)
    cp(input_json, archived_result_path; force = true)

    tests_endpoint = build_badge("HPC tests", status_message(tests), endpoint_color(string(get(tests, "status", "unknown"))))
    benchmarks_endpoint = build_badge(
        "HPC bench",
        status_message(benchmarks),
        endpoint_color(string(get(benchmarks, "status", "unknown"))),
    )
    meta = Dict(
        "schema_version" => get(result, "schema_version", 1),
        "commit" => commit,
        "branch" => string(get(result, "branch", "unknown")),
        "generated_at" => string(get(result, "generated_at", Dates.format(Dates.now(Dates.UTC), Dates.dateformat"yyyy-mm-ddTHH:MM:SSZ"))),
        "cluster" => string(get(hpc, "cluster", "unknown")),
        "job_id" => string(get(hpc, "job_id", "")),
        "node" => string(get(hpc, "node", "")),
        "tests_status" => string(get(tests, "status", "unknown")),
        "benchmarks_status" => string(get(benchmarks, "status", "unknown")),
        "coverage_percent" => coverage isa AbstractDict ? get(coverage, "coverage_percent", nothing) : nothing,
    )

    write_json(joinpath(output_dir, "latest-tests.json"), tests_endpoint)
    write_json(joinpath(output_dir, "latest-benchmarks.json"), benchmarks_endpoint)
    write_json(joinpath(output_dir, "latest-meta.json"), meta)
    write_json(joinpath(archive_dir, "tests.json"), tests_endpoint)
    write_json(joinpath(archive_dir, "benchmarks.json"), benchmarks_endpoint)
    write_json(joinpath(archive_dir, "meta.json"), meta)

    open(joinpath(output_dir, "README.md"), "w") do io
        write(io, readme_text(result, repository_slug, branch_name))
        write(io, '\n')
    end

    open(joinpath(output_dir, ".nojekyll"), "w") do io
        write(io, "")
    end
end

main()
