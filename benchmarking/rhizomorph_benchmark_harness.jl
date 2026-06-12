#!/usr/bin/env julia

import BioSequences
import DataFrames
import FASTX
import Mycelia
import TOML

include("benchmark_artifacts.jl")
include("standard_assembler_fixtures.jl")

const RHIZOMORPH_BENCHMARK_MANIFEST =
    joinpath(@__DIR__, "rhizomorph_benchmark_manifest.toml")

const RHIZOMORPH_BENCHMARK_SCALE_RANK = Dict(
    "ci" => 1,
    "full" => 2,
    "candidate" => 3
)

const RHIZOMORPH_EXECUTABLE_HYPOTHESES = Set(["H1", "H2", "H7"])
const RHIZOMORPH_DEFAULT_ASSEMBLERS = ["Rhizomorph", "MEGAHIT", "metaSPAdes"]

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

Return the dry-run benchmark plan, or execute implemented CI-safe slices and
write public-record artifacts when `dry_run=false`.
"""
function run_rhizomorph_benchmark_harness(;
        dry_run::Bool = true,
        output_dir::AbstractString = joinpath("results", "rhizomorph-benchmark"),
        run_id::AbstractString = "rhizomorph_benchmark_suite",
        command_args = String[],
        assemblers = RHIZOMORPH_DEFAULT_ASSEMBLERS,
        run_external::Bool = lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true",
        threads::Int = _benchmark_threads(),
        generated_at = nothing,
        kwargs...)
    plan = build_rhizomorph_benchmark_plan(; kwargs...)
    if dry_run
        return plan
    end
    scale = haskey(kwargs, :scale) ? string(kwargs[:scale]) : "ci"
    return execute_rhizomorph_benchmark_plan(
        plan;
        output_dir = output_dir,
        run_id = run_id,
        scale = scale,
        command_args = command_args,
        assemblers = string.(assemblers),
        run_external = run_external,
        threads = threads,
        generated_at = generated_at
    )
end

"""
    execute_rhizomorph_benchmark_plan(plan; output_dir, run_id, scale, ...)

Execute implemented Rhizomorph benchmark slices and write a reproducible
artifact bundle. H1 records graph-construction scalability/efficiency, H2
records Rhizomorph assembly accuracy, and H7 records assembler comparison rows.
Unsupported datasets or external tools are represented as skipped rows.
"""
function execute_rhizomorph_benchmark_plan(
        plan::DataFrames.DataFrame;
        output_dir::AbstractString,
        run_id::AbstractString,
        scale::AbstractString,
        command_args = String[],
        assemblers = RHIZOMORPH_DEFAULT_ASSEMBLERS,
        run_external::Bool = false,
        threads::Int = _benchmark_threads(),
        generated_at = nothing)
    mkpath(output_dir)
    work_dir = joinpath(output_dir, "work")
    mkpath(work_dir)

    graph_metrics = _run_graph_construction_slice(plan, work_dir, scale)
    assembly_metrics = _run_assembly_accuracy_slice(plan, work_dir)
    assembler_metrics = _run_assembler_comparison_slice(
        plan,
        work_dir;
        assemblers = string.(assemblers),
        run_external = run_external,
        threads = threads
    )
    summary = _summarize_rhizomorph_benchmark_tables(
        graph_metrics,
        assembly_metrics,
        assembler_metrics,
        plan
    )

    tables = [
        "rhizomorph_benchmark_plan" => plan,
        "graph_construction_metrics" => graph_metrics,
        "assembly_accuracy_metrics" => assembly_metrics,
        "assembler_comparison_metrics" => assembler_metrics,
        "benchmark_suite_summary" => summary
    ]
    table_context_columns = Dict{String, Dict{String, String}}(
        "rhizomorph_benchmark_plan" => Dict(
            "benchmark_dataset_id" => "dataset_id",
            "benchmark_hypothesis_id" => "hypothesis_id"
        ),
        "graph_construction_metrics" => Dict(
            "benchmark_dataset_id" => "dataset_id",
            "benchmark_hypothesis_id" => "hypothesis_id"
        ),
        "assembly_accuracy_metrics" => Dict(
            "benchmark_dataset_id" => "dataset_id",
            "benchmark_hypothesis_id" => "hypothesis_id"
        ),
        "assembler_comparison_metrics" => Dict(
            "benchmark_dataset_id" => "dataset_id",
            "benchmark_hypothesis_id" => "hypothesis_id"
        ),
        "benchmark_suite_summary" => Dict(
            "benchmark_hypothesis_id" => "hypothesis_id"
        )
    )

    return write_benchmark_artifacts(
        tables,
        output_dir = output_dir,
        run_id = run_id,
        scale = scale,
        dataset_ids = unique(string.(plan.dataset_id)),
        command_args = command_args,
        tool_versions = _assembler_tool_versions(assemblers),
        generated_at = generated_at,
        metadata = Dict{String, Any}(
            "artifact_kind" => "rhizomorph_benchmark_suite",
            "implemented_hypotheses" => sort(collect(RHIZOMORPH_EXECUTABLE_HYPOTHESES)),
            "run_external" => run_external,
            "threads" => threads,
            "work_dir" => abspath(work_dir)
        ),
        table_context_columns = table_context_columns
    )
end

"""
    write_rhizomorph_benchmark_plan_artifacts(; output_dir, scale="ci", kwargs...)

Write smoke-testable public-record artifacts for the Rhizomorph benchmark plan.

The emitted layout is stable for CI, full, and candidate scales:
`tables/`, `plots/`, `logs/`, and `provenance/`. The table schemas do not
change with scale; only row membership changes.
"""
function write_rhizomorph_benchmark_plan_artifacts(;
        output_dir::AbstractString,
        scale::AbstractString = "ci",
        hypothesis_ids = nothing,
        dataset_ids = nothing,
        manifest_path::AbstractString = RHIZOMORPH_BENCHMARK_MANIFEST,
        run_id::AbstractString = "rhizomorph_benchmark_plan",
        command_args = String[],
        tool_versions = Dict{String, Any}(),
        generated_at = nothing)
    plan = build_rhizomorph_benchmark_plan(
        scale = scale,
        hypothesis_ids = hypothesis_ids,
        dataset_ids = dataset_ids,
        manifest_path = manifest_path
    )
    datasets = list_rhizomorph_benchmark_datasets(manifest_path = manifest_path)
    slices = list_rhizomorph_benchmark_slices(manifest_path = manifest_path)
    selected_dataset_ids = unique(string.(plan.dataset_id))
    selected_hypothesis_ids = unique(string.(plan.hypothesis_id))

    tables = [
        "rhizomorph_benchmark_plan" => plan,
        "rhizomorph_benchmark_datasets" => datasets,
        "rhizomorph_benchmark_slices" => slices
    ]
    table_context_columns = Dict{String, Dict{String, String}}(
        "rhizomorph_benchmark_plan" => Dict(
            "benchmark_dataset_id" => "dataset_id",
            "benchmark_hypothesis_id" => "hypothesis_id"
        ),
        "rhizomorph_benchmark_datasets" => Dict(
            "benchmark_dataset_id" => "id"
        ),
        "rhizomorph_benchmark_slices" => Dict(
            "benchmark_hypothesis_id" => "id"
        )
    )

    return write_benchmark_artifacts(
        tables,
        output_dir = output_dir,
        run_id = run_id,
        scale = scale,
        dataset_ids = selected_dataset_ids,
        command_args = command_args,
        tool_versions = tool_versions,
        generated_at = generated_at,
        metadata = Dict{String, Any}(
            "manifest_path" => abspath(manifest_path),
            "selected_hypothesis_ids" => selected_hypothesis_ids,
            "artifact_kind" => "rhizomorph_benchmark_plan"
        ),
        table_context_columns = table_context_columns
    )
end

function _run_graph_construction_slice(plan::DataFrames.DataFrame, work_dir::AbstractString, scale::AbstractString)
    rows = NamedTuple[]
    slice_rows = DataFrames.filter(row -> row.hypothesis_id == "H1", plan)
    k_values = _benchmark_k_values(scale)
    modes = _benchmark_graph_modes(scale)
    memory_profiles = scale == "ci" ? [:ultralight] : [:ultralight, :full]

    for row in DataFrames.eachrow(slice_rows)
        dataset_id = string(row.dataset_id)
        if !_is_standard_assembler_fixture(dataset_id)
            push!(rows, _graph_construction_row(
                dataset_id = dataset_id,
                fixture_kind = string(row.dataset_category),
                graph_type = "kmer",
                mode = "",
                memory_profile = "",
                k = missing,
                input_records = missing,
                total_reference_bases = missing,
                input_bytes = missing,
                status = "skipped",
                runtime_seconds = missing,
                allocated_bytes = missing,
                gc_time_seconds = missing,
                vertices = missing,
                edges = missing,
                total_observations = missing,
                error = "No CI-safe materializer is registered for $(dataset_id)."
            ))
            continue
        end

        fixture_dir = joinpath(work_dir, "H1", dataset_id)
        fixture = materialize_standard_assembler_fixture(dataset_id; outdir = fixture_dir, emit_reads = false)
        records = collect(Mycelia.open_fastx(fixture.reference_fasta))
        input_bytes = filesize(fixture.reference_fasta)

        for k in k_values
            for mode in modes
                for memory_profile in memory_profiles
                    timed = nothing
                    graph = nothing
                    status = "ok"
                    error_message = missing
                    try
                        timed = @timed Mycelia.Rhizomorph.build_kmer_graph(
                            records,
                            k;
                            dataset_id = dataset_id,
                            mode = mode,
                            type_hint = :DNA,
                            memory_profile = memory_profile
                        )
                        graph = timed.value
                    catch error_value
                        status = "failed"
                        error_message = sprint(showerror, error_value)
                    end

                    stats = status == "ok" ? Mycelia.Rhizomorph.get_kmer_statistics(graph) : Dict()
                    push!(rows, _graph_construction_row(
                        dataset_id = dataset_id,
                        fixture_kind = string(row.dataset_category),
                        graph_type = "kmer",
                        mode = string(mode),
                        memory_profile = string(memory_profile),
                        k = k,
                        input_records = length(records),
                        total_reference_bases = fixture.total_reference_bases,
                        input_bytes = input_bytes,
                        status = status,
                        runtime_seconds = status == "ok" ? timed.time : missing,
                        allocated_bytes = status == "ok" ? timed.bytes : missing,
                        gc_time_seconds = status == "ok" ? _timed_property(timed, :gctime, missing) : missing,
                        vertices = status == "ok" ? Int(get(stats, :num_vertices, 0)) : missing,
                        edges = status == "ok" ? Int(get(stats, :num_edges, 0)) : missing,
                        total_observations = status == "ok" ? Int(get(stats, :total_observations, 0)) : missing,
                        error = error_message
                    ))
                end
            end
        end
    end

    return DataFrames.DataFrame(rows)
end

function _run_assembly_accuracy_slice(plan::DataFrames.DataFrame, work_dir::AbstractString)
    rows = NamedTuple[]
    slice_rows = DataFrames.filter(row -> row.hypothesis_id == "H2", plan)

    for row in DataFrames.eachrow(slice_rows)
        dataset_id = string(row.dataset_id)
        if !_is_standard_assembler_fixture(dataset_id)
            push!(rows, _assembly_accuracy_row(
                dataset_id = dataset_id,
                fixture_kind = string(row.dataset_category),
                assembler = "Rhizomorph",
                status = "skipped",
                runtime_seconds = missing,
                allocated_bytes = missing,
                n_contigs = missing,
                total_length = missing,
                n50 = missing,
                l50 = missing,
                longest_contig = missing,
                length_recovery = missing,
                kmer_precision = missing,
                kmer_recall = missing,
                contigs_path = missing,
                error = "No CI-safe read materializer is registered for $(dataset_id)."
            ))
            continue
        end

        fixture_dir = joinpath(work_dir, "H2", dataset_id, "fixture")
        fixture = materialize_standard_assembler_fixture(dataset_id; outdir = fixture_dir, emit_reads = true)
        run_result = _run_rhizomorph_fixture_assembly(
            fixture,
            joinpath(work_dir, "H2", dataset_id, "rhizomorph")
        )
        push!(rows, _assembly_accuracy_row_from_result(
            dataset_id,
            string(row.dataset_category),
            "Rhizomorph",
            fixture,
            run_result
        ))
    end

    return DataFrames.DataFrame(rows)
end

function _run_assembler_comparison_slice(
        plan::DataFrames.DataFrame,
        work_dir::AbstractString;
        assemblers,
        run_external::Bool,
        threads::Int)
    rows = NamedTuple[]
    slice_rows = DataFrames.filter(row -> row.hypothesis_id == "H7", plan)

    for row in DataFrames.eachrow(slice_rows)
        dataset_id = string(row.dataset_id)
        if !_is_standard_assembler_fixture(dataset_id)
            for assembler in assemblers
                push!(rows, _assembler_comparison_row(
                    dataset_id = dataset_id,
                    fixture_kind = string(row.dataset_category),
                    assembler = assembler,
                    status = "skipped",
                    runtime_seconds = missing,
                    allocated_bytes = missing,
                    n_contigs = missing,
                    total_length = missing,
                    n50 = missing,
                    l50 = missing,
                    longest_contig = missing,
                    length_recovery = missing,
                    kmer_precision = missing,
                    kmer_recall = missing,
                    contigs_path = missing,
                    tool_available = false,
                    run_external = run_external,
                    error = "No CI-safe read materializer is registered for $(dataset_id)."
                ))
            end
            continue
        end

        fixture_dir = joinpath(work_dir, "H7", dataset_id, "fixture")
        fixture = materialize_standard_assembler_fixture(dataset_id; outdir = fixture_dir, emit_reads = true)
        for assembler in assemblers
            assembler_dir = joinpath(work_dir, "H7", dataset_id, _artifact_basename(assembler))
            run_result = _run_benchmark_assembler(
                assembler,
                fixture,
                assembler_dir;
                run_external = run_external,
                threads = threads
            )
            push!(rows, _assembler_comparison_row_from_result(
                dataset_id,
                string(row.dataset_category),
                assembler,
                fixture,
                run_result,
                run_external
            ))
        end
    end

    return DataFrames.DataFrame(rows)
end

function _graph_construction_row(; kwargs...)
    values = Dict{Symbol, Any}(kwargs)
    return (
        hypothesis_id = "H1",
        dataset_id = values[:dataset_id],
        fixture_kind = values[:fixture_kind],
        graph_type = values[:graph_type],
        mode = values[:mode],
        memory_profile = values[:memory_profile],
        k = values[:k],
        input_records = values[:input_records],
        total_reference_bases = values[:total_reference_bases],
        input_bytes = values[:input_bytes],
        status = values[:status],
        runtime_seconds = values[:runtime_seconds],
        allocated_bytes = values[:allocated_bytes],
        gc_time_seconds = values[:gc_time_seconds],
        vertices = values[:vertices],
        edges = values[:edges],
        total_observations = values[:total_observations],
        error = values[:error]
    )
end

function _assembly_accuracy_row(; kwargs...)
    values = Dict{Symbol, Any}(kwargs)
    return (
        hypothesis_id = "H2",
        dataset_id = values[:dataset_id],
        fixture_kind = values[:fixture_kind],
        assembler = values[:assembler],
        status = values[:status],
        runtime_seconds = values[:runtime_seconds],
        allocated_bytes = values[:allocated_bytes],
        n_contigs = values[:n_contigs],
        total_length = values[:total_length],
        n50 = values[:n50],
        l50 = values[:l50],
        longest_contig = values[:longest_contig],
        length_recovery = values[:length_recovery],
        kmer_precision = values[:kmer_precision],
        kmer_recall = values[:kmer_recall],
        contigs_path = values[:contigs_path],
        error = values[:error]
    )
end

function _assembler_comparison_row(; kwargs...)
    values = Dict{Symbol, Any}(kwargs)
    return (
        hypothesis_id = "H7",
        dataset_id = values[:dataset_id],
        fixture_kind = values[:fixture_kind],
        assembler = values[:assembler],
        status = values[:status],
        runtime_seconds = values[:runtime_seconds],
        allocated_bytes = values[:allocated_bytes],
        n_contigs = values[:n_contigs],
        total_length = values[:total_length],
        n50 = values[:n50],
        l50 = values[:l50],
        longest_contig = values[:longest_contig],
        length_recovery = values[:length_recovery],
        kmer_precision = values[:kmer_precision],
        kmer_recall = values[:kmer_recall],
        contigs_path = values[:contigs_path],
        tool_available = values[:tool_available],
        run_external = values[:run_external],
        error = values[:error]
    )
end

function _assembly_accuracy_row_from_result(
        dataset_id::AbstractString,
        fixture_kind::AbstractString,
        assembler::AbstractString,
        fixture,
        run_result)
    if run_result.status != "ok"
        return _assembly_accuracy_row(
            dataset_id = dataset_id,
            fixture_kind = fixture_kind,
            assembler = assembler,
            status = run_result.status,
            runtime_seconds = getproperty(run_result, :runtime_seconds),
            allocated_bytes = getproperty(run_result, :allocated_bytes),
            n_contigs = missing,
            total_length = missing,
            n50 = missing,
            l50 = missing,
            longest_contig = missing,
            length_recovery = missing,
            kmer_precision = missing,
            kmer_recall = missing,
            contigs_path = getproperty(run_result, :contigs),
            error = getproperty(run_result, :error)
        )
    end

    metrics = _contig_metrics(run_result.contigs, fixture.total_reference_bases)
    accuracy = _kmer_accuracy_metrics(fixture.reference_fasta, run_result.contigs)
    return _assembly_accuracy_row(
        dataset_id = dataset_id,
        fixture_kind = fixture_kind,
        assembler = assembler,
        status = "ok",
        runtime_seconds = run_result.runtime_seconds,
        allocated_bytes = run_result.allocated_bytes,
        n_contigs = metrics.n_contigs,
        total_length = metrics.total_length,
        n50 = metrics.n50,
        l50 = metrics.l50,
        longest_contig = metrics.longest_contig,
        length_recovery = metrics.length_recovery,
        kmer_precision = accuracy.precision,
        kmer_recall = accuracy.recall,
        contigs_path = metrics.contigs_path,
        error = missing
    )
end

function _assembler_comparison_row_from_result(
        dataset_id::AbstractString,
        fixture_kind::AbstractString,
        assembler::AbstractString,
        fixture,
        run_result,
        run_external::Bool)
    base_row = _assembly_accuracy_row_from_result(
        dataset_id,
        fixture_kind,
        assembler,
        fixture,
        run_result
    )
    return _assembler_comparison_row(
        dataset_id = base_row.dataset_id,
        fixture_kind = base_row.fixture_kind,
        assembler = assembler,
        status = base_row.status,
        runtime_seconds = base_row.runtime_seconds,
        allocated_bytes = base_row.allocated_bytes,
        n_contigs = base_row.n_contigs,
        total_length = base_row.total_length,
        n50 = base_row.n50,
        l50 = base_row.l50,
        longest_contig = base_row.longest_contig,
        length_recovery = base_row.length_recovery,
        kmer_precision = base_row.kmer_precision,
        kmer_recall = base_row.kmer_recall,
        contigs_path = base_row.contigs_path,
        tool_available = _assembler_tool_available(assembler),
        run_external = run_external,
        error = base_row.error
    )
end

function _run_benchmark_assembler(
        assembler::AbstractString,
        fixture,
        outdir::AbstractString;
        run_external::Bool,
        threads::Int)
    if assembler == "Rhizomorph"
        return _run_rhizomorph_fixture_assembly(fixture, outdir)
    end

    if !run_external
        return (
            status = "skipped",
            runtime_seconds = missing,
            allocated_bytes = missing,
            contigs = missing,
            error = "Set MYCELIA_RUN_EXTERNAL=true or pass --run-external to run $(assembler)."
        )
    end

    # MEGAHIT citation: Li et al. 2015, DOI 10.1093/bioinformatics/btv033.
    # metaSPAdes citation: Nurk et al. 2017, DOI 10.1101/gr.213959.116.
    mkpath(outdir)
    result = nothing
    try
        timed = @timed begin
            if assembler == "MEGAHIT"
                result = Mycelia.run_megahit(
                    fastq1 = fixture.fastq1,
                    fastq2 = fixture.fastq2,
                    outdir = joinpath(outdir, "megahit"),
                    k_list = "21",
                    min_contig_len = 100,
                    threads = threads
                )
            elseif assembler == "metaSPAdes"
                result = Mycelia.run_metaspades(
                    fastq1 = fixture.fastq1,
                    fastq2 = fixture.fastq2,
                    outdir = joinpath(outdir, "metaspades"),
                    k_list = "21",
                    threads = threads
                )
            else
                error("Unsupported assembler: $(assembler)")
            end
        end
        return (
            status = "ok",
            runtime_seconds = timed.time,
            allocated_bytes = timed.bytes,
            contigs = result.contigs,
            error = missing
        )
    catch error_value
        return (
            status = "failed",
            runtime_seconds = missing,
            allocated_bytes = missing,
            contigs = missing,
            error = sprint(showerror, error_value)
        )
    end
end

function _run_rhizomorph_fixture_assembly(fixture, outdir::AbstractString)
    mkpath(outdir)
    contigs_fasta = joinpath(outdir, "rhizomorph.contigs.fasta")
    try
        records = FASTX.FASTQ.Record[]
        append!(records, collect(Mycelia.open_fastx(fixture.fastq1)))
        append!(records, collect(Mycelia.open_fastx(fixture.fastq2)))

        result = nothing
        timed = @timed begin
            result = Mycelia.Rhizomorph.assemble_genome(records; k = 21)
            _write_rhizomorph_contigs(result, contigs_fasta)
        end
        return (
            status = "ok",
            runtime_seconds = timed.time,
            allocated_bytes = timed.bytes,
            contigs = contigs_fasta,
            error = missing
        )
    catch error_value
        return (
            status = "failed",
            runtime_seconds = missing,
            allocated_bytes = missing,
            contigs = contigs_fasta,
            error = sprint(showerror, error_value)
        )
    end
end

function _write_rhizomorph_contigs(result, contigs_fasta::AbstractString)
    records = FASTX.FASTA.Record[]
    for (index, contig_sequence) in enumerate(result.contigs)
        contig_name = index <= length(result.contig_names) ? result.contig_names[index] : "contig_$(index)"
        push!(records, FASTX.FASTA.Record(contig_name, BioSequences.LongDNA{4}(contig_sequence)))
    end
    if isempty(records)
        open(contigs_fasta, "w") do io
            write(io, "")
        end
    else
        Mycelia.write_fasta(outfile = contigs_fasta, records = records, gzip = false)
    end
    return contigs_fasta
end

function _contig_metrics(contigs_fasta, expected_total_length::Integer)
    if contigs_fasta === missing || !isfile(contigs_fasta)
        return (
            contigs_path = contigs_fasta,
            n_contigs = missing,
            total_length = missing,
            n50 = missing,
            l50 = missing,
            longest_contig = missing,
            length_recovery = missing
        )
    end

    n_contigs, total_length, n50, l50 = Mycelia.assess_assembly_quality(contigs_fasta)
    longest_contig = 0
    for record in Mycelia.open_fastx(contigs_fasta)
        longest_contig = max(longest_contig, length(FASTX.sequence(record)))
    end

    return (
        contigs_path = contigs_fasta,
        n_contigs = n_contigs,
        total_length = total_length,
        n50 = n50,
        l50 = l50,
        longest_contig = longest_contig,
        length_recovery = expected_total_length == 0 ? missing : total_length / expected_total_length
    )
end

function _kmer_accuracy_metrics(reference_fasta::AbstractString, contigs_fasta, k::Int = 21)
    if contigs_fasta === missing || !isfile(contigs_fasta)
        return (precision = missing, recall = missing)
    end

    reference_kmers = _fasta_kmer_set(reference_fasta, k)
    contig_kmers = _fasta_kmer_set(contigs_fasta, k)
    if isempty(reference_kmers) || isempty(contig_kmers)
        return (precision = missing, recall = missing)
    end

    shared = length(intersect(reference_kmers, contig_kmers))
    return (
        precision = shared / length(contig_kmers),
        recall = shared / length(reference_kmers)
    )
end

function _fasta_kmer_set(fasta_path::AbstractString, k::Int)
    kmers = Set{String}()
    for record in Mycelia.open_fastx(fasta_path)
        sequence = string(FASTX.sequence(record))
        if length(sequence) < k
            continue
        end
        for index in 1:(lastindex(sequence) - k + 1)
            push!(kmers, sequence[index:(index + k - 1)])
        end
    end
    return kmers
end

function _summarize_rhizomorph_benchmark_tables(
        graph_metrics::DataFrames.DataFrame,
        assembly_metrics::DataFrames.DataFrame,
        assembler_metrics::DataFrames.DataFrame,
        plan::DataFrames.DataFrame)
    rows = NamedTuple[]
    for hypothesis_id in sort(unique(string.(plan.hypothesis_id)))
        table = if hypothesis_id == "H1"
            graph_metrics
        elseif hypothesis_id == "H2"
            assembly_metrics
        elseif hypothesis_id == "H7"
            assembler_metrics
        else
            DataFrames.DataFrame()
        end
        push!(rows, (
            hypothesis_id = hypothesis_id,
            implemented = hypothesis_id in RHIZOMORPH_EXECUTABLE_HYPOTHESES,
            planned_rows = DataFrames.nrow(DataFrames.filter(row -> row.hypothesis_id == hypothesis_id, plan)),
            result_rows = DataFrames.nrow(table),
            ok_rows = "status" in names(table) ? count(==("ok"), table.status) : 0,
            skipped_rows = "status" in names(table) ? count(==("skipped"), table.status) : 0,
            failed_rows = "status" in names(table) ? count(==("failed"), table.status) : 0
        ))
    end
    return DataFrames.DataFrame(rows)
end

function _is_standard_assembler_fixture(dataset_id::AbstractString)
    return haskey(STANDARD_ASSEMBLER_FIXTURE_SPECS, string(dataset_id))
end

function _benchmark_k_values(scale::AbstractString)
    if string(scale) == "ci"
        return [15, 21]
    end
    return [15, 21, 31]
end

function _benchmark_graph_modes(scale::AbstractString)
    if string(scale) == "ci"
        return [:singlestrand, :canonical]
    end
    return [:singlestrand, :canonical, :doublestrand]
end

function _benchmark_threads()
    parsed = tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_THREADS", "2"))
    return clamp(something(parsed, 2), 1, 64)
end

function _timed_property(timed, property::Symbol, default)
    return property in propertynames(timed) ? getproperty(timed, property) : default
end

function _assembler_tool_available(assembler::AbstractString)
    if assembler == "Rhizomorph"
        return true
    elseif assembler == "MEGAHIT"
        return Sys.which("megahit") !== nothing
    elseif assembler == "metaSPAdes"
        return Sys.which("metaspades.py") !== nothing || Sys.which("spades.py") !== nothing
    end
    return false
end

function _assembler_tool_versions(assemblers)
    versions = Dict{String, Any}()
    for assembler in string.(assemblers)
        versions[assembler] = Dict(
            "available_on_path" => _assembler_tool_available(assembler)
        )
    end
    return versions
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
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --plan --scale ci --write-artifacts --output-dir results/public-record")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --execute --slice H1 --scale ci --output-dir results/rhizomorph-ci")
    println()
    println("Scales: ci, full, candidate")
    println("Implemented executable slices: H1 graph construction, H2 Rhizomorph accuracy, H7 assembler comparison.")
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
        "--execute",
        "--write-artifacts",
        "--output-dir",
        "--run-id",
        "--assembler",
        "--run-external"
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
    assemblers = _flag_values(args, "--assembler")
    dry_run = !("--execute" in args)
    output_dir = _flag_value(args, "--output-dir", joinpath("results", "public-record"))
    run_id = _flag_value(args, "--run-id", dry_run ? "rhizomorph_benchmark_plan" : "rhizomorph_benchmark_suite")

    if dry_run
        plan = run_rhizomorph_benchmark_harness(
            dry_run = true,
            scale = scale,
            hypothesis_ids = isempty(hypothesis_ids) ? nothing : hypothesis_ids,
            dataset_ids = isempty(dataset_ids) ? nothing : dataset_ids,
            manifest_path = manifest_path
        )
        show(plan; allrows = true, allcols = true)
        println()
        if "--write-artifacts" in args
            artifacts = write_rhizomorph_benchmark_plan_artifacts(
                output_dir = output_dir,
                scale = scale,
                hypothesis_ids = isempty(hypothesis_ids) ? nothing : hypothesis_ids,
                dataset_ids = isempty(dataset_ids) ? nothing : dataset_ids,
                manifest_path = manifest_path,
                run_id = run_id,
                command_args = args
            )
            println("Artifacts written to: $(artifacts.root)")
            println("Artifact index: $(artifacts.index)")
            println("Run provenance: $(artifacts.provenance)")
        end
    else
        artifacts = run_rhizomorph_benchmark_harness(
            dry_run = false,
            output_dir = output_dir,
            run_id = run_id,
            command_args = args,
            assemblers = isempty(assemblers) ? RHIZOMORPH_DEFAULT_ASSEMBLERS : assemblers,
            run_external = "--run-external" in args ||
                           lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true",
            scale = scale,
            hypothesis_ids = isempty(hypothesis_ids) ? nothing : hypothesis_ids,
            dataset_ids = isempty(dataset_ids) ? nothing : dataset_ids,
            manifest_path = manifest_path
        )
        println("Benchmark artifacts written to: $(artifacts.root)")
        println("Artifact index: $(artifacts.index)")
        println("Run provenance: $(artifacts.provenance)")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
