#!/usr/bin/env julia

import DataFrames
import TOML

include("benchmark_artifacts.jl")

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

struct H1SyntheticEdge
    source::String
    target::String
    probability::Float64
end

struct H1SyntheticFixture
    id::String
    name::String
    source::String
    sink::String
    edges::Vector{H1SyntheticEdge}
    truth_path::Vector{String}
    expected_greedy_path::Vector{String}
    data_type::String
    coverage::String
    seed::Int
    ambiguity_margin::Union{Missing, Float64}
    observation_classes::Union{Nothing, Vector{String}}
    vertex_classes::Dict{String, String}
    expected_greedy_failure_code::String
end

"""
    h1_viterbi_dp_greedy_fixtures()

Return the pre-registered clean H1-G0 through H1-G4 synthetic fixtures from
`rhizomorph-paper/planning/PLAN-2026-06-02-h1-viterbi-dp-greedy-benchmark.md`.
This is still a local synthetic path-metric gate, not the real-data H1 assembly
quality decision rule.
"""
function h1_viterbi_dp_greedy_fixtures()::Vector{H1SyntheticFixture}
    return H1SyntheticFixture[
        _h1_fixture(
            "H1-G0",
            "linear control",
            H1SyntheticEdge[
                H1SyntheticEdge("S", "A1", 0.99),
                H1SyntheticEdge("A1", "A2", 0.99),
                H1SyntheticEdge("A2", "A3", 0.99),
                H1SyntheticEdge("A3", "T", 0.99)
            ],
            ["S", "A1", "A2", "A3", "T"],
            ["S", "A1", "A2", "A3", "T"]
        ),
        _h1_fixture(
            "H1-G1",
            "local-trap bubble",
            H1SyntheticEdge[
                H1SyntheticEdge("S", "A1", 0.99),
                H1SyntheticEdge("A1", "A2", 0.50),
                H1SyntheticEdge("A2", "T", 0.50),
                H1SyntheticEdge("S", "B1", 0.90),
                H1SyntheticEdge("B1", "B2", 0.90),
                H1SyntheticEdge("B2", "T", 0.90)
            ],
            ["S", "B1", "B2", "T"],
            ["S", "A1", "A2", "T"]
        ),
        _h1_fixture(
            "H1-G2",
            "delayed-exit bubble",
            H1SyntheticEdge[
                H1SyntheticEdge("S", "A1", 0.98),
                H1SyntheticEdge("A1", "A2", 0.98),
                H1SyntheticEdge("A2", "A3", 0.98),
                H1SyntheticEdge("A3", "T", 0.20),
                H1SyntheticEdge("S", "B1", 0.80),
                H1SyntheticEdge("B1", "B2", 0.80),
                H1SyntheticEdge("B2", "B3", 0.80),
                H1SyntheticEdge("B3", "T", 0.80)
            ],
            ["S", "B1", "B2", "B3", "T"],
            ["S", "A1", "A2", "A3", "T"]
        ),
        _h1_fixture(
            "H1-G3",
            "repeat-copy ambiguity",
            H1SyntheticEdge[
                H1SyntheticEdge("S", "L1", 0.93),
                H1SyntheticEdge("L1", "R1", 0.96),
                H1SyntheticEdge("R1", "L2", 0.60),
                H1SyntheticEdge("L2", "R2", 0.96),
                H1SyntheticEdge("R2", "T", 0.93),
                H1SyntheticEdge("R1", "T", 0.88),
                H1SyntheticEdge("S", "C1", 0.89),
                H1SyntheticEdge("C1", "C2", 0.89),
                H1SyntheticEdge("C2", "C3", 0.89),
                H1SyntheticEdge("C3", "T", 0.89)
            ],
            ["S", "L1", "R1", "L2", "R2", "T"],
            ["S", "L1", "R1", "T"];
            observation_classes = ["L", "R", "L", "R"],
            vertex_classes = Dict(
                "L1" => "L",
                "L2" => "L",
                "R1" => "R",
                "R2" => "R",
                "C1" => "C",
                "C2" => "C",
                "C3" => "C"
            ),
            expected_greedy_failure_code = "length_mismatch"
        ),
        _h1_fixture(
            "H1-G4",
            "near-tie bubble",
            H1SyntheticEdge[
                H1SyntheticEdge("S", "A1", 0.900),
                H1SyntheticEdge("A1", "A2", 0.900),
                H1SyntheticEdge("A2", "T", 0.900),
                H1SyntheticEdge("S", "B1", 0.899),
                H1SyntheticEdge("B1", "B2", 0.899),
                H1SyntheticEdge("B2", "T", 0.899)
            ],
            ["S", "A1", "A2", "T"],
            ["S", "A1", "A2", "T"];
            ambiguity_margin = abs(3 * (log(0.900) - log(0.899))) / 3
        )
    ]
end

function _h1_fixture(
        id::AbstractString,
        name::AbstractString,
        edges::Vector{H1SyntheticEdge},
        truth_path::Vector{String},
        expected_greedy_path::Vector{String};
        data_type::AbstractString = "clean",
        coverage::AbstractString = "fixture",
        seed::Int = 0,
        ambiguity_margin::Union{Missing, Float64} = missing,
        observation_classes::Union{Nothing, Vector{String}} = nothing,
        vertex_classes::Dict{String, String} = Dict{String, String}(),
        expected_greedy_failure_code::AbstractString = "none")::H1SyntheticFixture
    return H1SyntheticFixture(
        string(id),
        string(name),
        "S",
        "T",
        edges,
        truth_path,
        expected_greedy_path,
        string(data_type),
        string(coverage),
        seed,
        ambiguity_margin,
        observation_classes,
        vertex_classes,
        string(expected_greedy_failure_code)
    )
end

"""
    run_h1_viterbi_dp_greedy_smoke()

Run the clean H1-G0 through H1-G4 Viterbi-objective-vs-greedy smoke benchmark. The
harness-level exhaustive oracle enumerates all source-to-sink paths and selects
the maximum shared log-likelihood. It is a tiny-fixture oracle for the Viterbi DP
objective, not a call into the production Rhizomorph ViterbiDP strategy. The
harness-level `GreedyViterbi` baseline selects the locally best outgoing
transition at each step with the same edge probabilities and deterministic
tie-breaking.
"""
function run_h1_viterbi_dp_greedy_smoke()::DataFrames.DataFrame
    rows = NamedTuple[]
    for fixture in h1_viterbi_dp_greedy_fixtures()
        dp_result = _h1_select_viterbi_dp_path(fixture)
        greedy_result = _h1_select_greedy_viterbi_path(fixture)
        oracle_log_probability = _h1_path_log_probability(fixture, fixture.truth_path)
        log_likelihood_gap = dp_result.log_probability - greedy_result.log_probability
        if dp_result.failure_code != "none" || dp_result.path != fixture.truth_path
            error("H1 smoke DP path expectation failed for $(fixture.id).")
        end
        if greedy_result.failure_code != fixture.expected_greedy_failure_code ||
           greedy_result.path != fixture.expected_greedy_path
            error("H1 smoke greedy path expectation failed for $(fixture.id).")
        end

        for result in (dp_result, greedy_result)
            exact_match = result.path == fixture.truth_path
            expected_path = result.algorithm == "greedy" ?
                            fixture.expected_greedy_path : fixture.truth_path
            expected_match = result.path == expected_path
            push!(rows, (
                hypothesis_id = "H1",
                dataset_id = "rhizomorph_graph_unit_fixtures",
                fixture_id = fixture.id,
                fixture_name = fixture.name,
                organism = "synthetic",
                data_type = fixture.data_type,
                coverage = fixture.coverage,
                seed = fixture.seed,
                ambiguity_margin = fixture.ambiguity_margin,
                algorithm = result.algorithm,
                strategy_name = result.strategy_name,
                implementation_scope = result.implementation_scope,
                path_vertices = join(result.path, ","),
                truth_vertices = join(fixture.truth_path, ","),
                expected_strategy_vertices = join(expected_path, ","),
                exact_path_match = exact_match,
                expected_strategy_path_match = expected_match,
                path_accuracy = exact_match ? 1.0 : 0.0,
                sequence_identity = exact_match ? 1.0 : 0.0,
                normalized_edit_distance = exact_match ? 0.0 : 1.0,
                repeat_copy_number_error = _h1_repeat_copy_number_error(fixture, result.path),
                path_confidence = _h1_path_confidence(result.log_probability, oracle_log_probability),
                sequence_metric_note = "path_match_proxy_no_sequence",
                log_probability = result.log_probability,
                oracle_log_probability = oracle_log_probability,
                oracle_log_probability_gap = oracle_log_probability - result.log_probability,
                dp_log_probability = dp_result.log_probability,
                greedy_log_probability = greedy_result.log_probability,
                delta_log_probability = log_likelihood_gap,
                log_likelihood_gap_dp_minus_greedy = log_likelihood_gap,
                dp_failure_code = dp_result.failure_code,
                greedy_failure_code = greedy_result.failure_code,
                pair_failure_code = _h1_pair_failure_code(dp_result.failure_code, greedy_result.failure_code),
                tie_breaking = "deterministic_lexicographic_vertex_strand_edge",
                tie_breaking_exercised = false,
                emission_model = "neutral",
                runtime_s = result.runtime_s,
                peak_rss_mib = _current_peak_rss_mib(),
                failure_code = result.failure_code
            ))
        end
    end

    return DataFrames.DataFrame(rows)
end

"""
    write_h1_viterbi_dp_greedy_artifacts(; output_dir, kwargs...)

Write the H1-G0 through H1-G4 smoke path-metrics CSV and provenance artifacts.
"""
function write_h1_viterbi_dp_greedy_artifacts(;
        output_dir::AbstractString,
        run_id::AbstractString = "h1_viterbi_dp_greedy_smoke",
        command_args = String[],
        generated_at = nothing)
    path_metrics = run_h1_viterbi_dp_greedy_smoke()
    tables = ["h1_viterbi_dp_greedy_path_metrics" => path_metrics]
    table_context_columns = Dict{String, Dict{String, String}}(
        "h1_viterbi_dp_greedy_path_metrics" => Dict(
            "benchmark_dataset_id" => "dataset_id",
            "benchmark_hypothesis_id" => "hypothesis_id"
        )
    )

    return write_benchmark_artifacts(
        tables,
        output_dir = output_dir,
        run_id = run_id,
        scale = "local-smoke",
        dataset_ids = unique(string.(path_metrics.dataset_id)),
        command_args = command_args,
        generated_at = generated_at,
        metadata = Dict{String, Any}(
            "artifact_kind" => "h1_viterbi_dp_greedy_path_metrics",
            "plan_path" => "rhizomorph-paper/planning/PLAN-2026-06-02-h1-viterbi-dp-greedy-benchmark.md",
            "fixtures" => "H1-G0,H1-G1,H1-G2,H1-G3,H1-G4",
            "scope" => "clean synthetic local-expanded smoke; no real-data H1 decision rule",
            "non_finite_metric_semantics" => "Inf/-Inf are permitted only in paired likelihood diagnostic columns when dp_failure_code, greedy_failure_code, or pair_failure_code documents an invalid or length-incompatible path; they are diagnostic sentinel values, not finite effect sizes.",
            "peak_rss_mib_semantics" => "peak_rss_mib records process-level Sys.maxrss() smoke-run provenance sampled per result row; it is not per-algorithm incremental memory.",
            "empty_layout_directory_semantics" => "plots/ and logs/ are reserved public-record layout directories and may be empty for this table-only smoke artifact."
        ),
        table_context_columns = table_context_columns
    )
end

"""
    run_rhizomorph_benchmark_harness(; dry_run=true, kwargs...)

Return the dry-run benchmark plan by default. With `dry_run=false`, execute the
implemented H1-G0 through H1-G4 Viterbi-DP-vs-greedy synthetic smoke when `H1` is the
requested slice.
"""
function run_rhizomorph_benchmark_harness(; dry_run::Bool = true, kwargs...)
    if haskey(kwargs, :scale)
        _validate_scale(string(kwargs[:scale]))
    end
    if dry_run
        return build_rhizomorph_benchmark_plan(; kwargs...)
    end

    selected_hypothesis_ids = _selected_id_set(get(kwargs, :hypothesis_ids, nothing))
    selected_dataset_ids = _selected_id_set(get(kwargs, :dataset_ids, nothing))
    if selected_hypothesis_ids != Set(["H1"])
        error("Only the H1 Viterbi DP vs greedy smoke runner is implemented for --execute; pass --slice H1.")
    end
    if selected_dataset_ids !== nothing &&
       selected_dataset_ids != Set(["rhizomorph_graph_unit_fixtures"])
        error("H1 Viterbi DP vs greedy smoke only supports dataset rhizomorph_graph_unit_fixtures.")
    end

    execute_plan = build_rhizomorph_benchmark_plan(; kwargs...)
    if DataFrames.nrow(execute_plan) != 1 ||
       only(execute_plan.hypothesis_id) != "H1" ||
       only(execute_plan.dataset_id) != "rhizomorph_graph_unit_fixtures" ||
       !only(execute_plan.implemented)
        error("H1 execute mode requires the manifest-backed implemented H1 rhizomorph_graph_unit_fixtures row.")
    end

    return run_h1_viterbi_dp_greedy_smoke()
end

function _h1_pair_failure_code(dp_failure_code::AbstractString, greedy_failure_code::AbstractString)::String
    if dp_failure_code == "none" && greedy_failure_code == "none"
        return "none"
    end
    return "dp=$(dp_failure_code);greedy=$(greedy_failure_code)"
end

function _current_peak_rss_mib()::Float64
    return Sys.maxrss() / 1024^2
end

function _h1_path_confidence(log_probability::Float64, oracle_log_probability::Float64)::Float64
    if !isfinite(log_probability) || !isfinite(oracle_log_probability)
        return 0.0
    end
    return exp(log_probability - oracle_log_probability)
end

function _h1_repeat_copy_number_error(fixture::H1SyntheticFixture, path::Vector{String})::Int
    if fixture.id != "H1-G3"
        return 0
    end
    expected_copies = count(vertex -> startswith(vertex, "L"), fixture.truth_path)
    observed_copies = count(vertex -> startswith(vertex, "L"), path)
    return abs(expected_copies - observed_copies)
end

function _h1_select_viterbi_dp_path(fixture::H1SyntheticFixture)
    start_time = time()
    paths = _h1_enumerate_paths(fixture)
    if isempty(paths)
        return (
            algorithm = "viterbi_objective_oracle",
            strategy_name = "ExhaustiveViterbiObjectiveOracle",
            implementation_scope = "tiny_fixture_exhaustive_oracle_not_production_strategy",
            path = String[],
            log_probability = -Inf,
            runtime_s = time() - start_time,
            failure_code = "no_path"
        )
    end

    ranked = sort(
        [(path, _h1_path_log_probability(fixture, path), _h1_path_tie_key(fixture, path)) for path in paths];
        by = item -> (-item[2], item[3])
    )
    best_path, best_log_probability, _ = first(ranked)
    return (
        algorithm = "viterbi_objective_oracle",
        strategy_name = "ExhaustiveViterbiObjectiveOracle",
        implementation_scope = "tiny_fixture_exhaustive_oracle_not_production_strategy",
        path = best_path,
        log_probability = best_log_probability,
        runtime_s = time() - start_time,
        failure_code = "none"
    )
end

function _h1_select_greedy_viterbi_path(fixture::H1SyntheticFixture)
    start_time = time()
    outgoing = _h1_outgoing_edges(fixture)
    path = [fixture.source]
    current = fixture.source
    max_steps = length(fixture.edges) + 1

    for _ in 1:max_steps
        if current == fixture.sink
            return (
                algorithm = "greedy",
                strategy_name = "GreedyViterbi",
                implementation_scope = "harness_local_greedy_baseline",
                path = path,
                log_probability = _h1_path_log_probability(fixture, path),
                runtime_s = time() - start_time,
                failure_code = _h1_path_failure_code(fixture, path)
            )
        end

        candidates = get(outgoing, current, H1SyntheticEdge[])
        if isempty(candidates)
            return (
                algorithm = "greedy",
                strategy_name = "GreedyViterbi",
                implementation_scope = "harness_local_greedy_baseline",
                path = path,
                log_probability = -Inf,
                runtime_s = time() - start_time,
                failure_code = "no_path"
            )
        end

        best_edge = first(sort(candidates; by = edge -> (-log(edge.probability), _h1_edge_tie_key(edge))))
        push!(path, best_edge.target)
        current = best_edge.target
    end

    return (
        algorithm = "greedy",
        strategy_name = "GreedyViterbi",
        implementation_scope = "harness_local_greedy_baseline",
        path = path,
        log_probability = -Inf,
        runtime_s = time() - start_time,
        failure_code = "cycle_limit"
    )
end

function _h1_enumerate_paths(fixture::H1SyntheticFixture)::Vector{Vector{String}}
    outgoing = _h1_outgoing_edges(fixture)
    paths = Vector{String}[]
    stack = [String[fixture.source]]
    max_vertices = length(fixture.edges) + 2

    while !isempty(stack)
        path = pop!(stack)
        current = last(path)
        if current == fixture.sink
            push!(paths, path)
            continue
        end
        if length(path) > max_vertices
            continue
        end

        for edge in reverse(sort(get(outgoing, current, H1SyntheticEdge[]); by = _h1_edge_tie_key))
            if edge.target in path
                continue
            end
            push!(stack, vcat(path, [edge.target]))
        end
    end

    return paths
end

function _h1_path_log_probability(fixture::H1SyntheticFixture, path::Vector{String})::Float64
    edge_probability_by_key = Dict(_h1_edge_key(edge) => edge.probability for edge in fixture.edges)
    total = 0.0
    for index in 1:(length(path) - 1)
        probability = get(edge_probability_by_key, "$(path[index])->$(path[index + 1])", nothing)
        if probability === nothing
            return -Inf
        end
        total += log(probability)
    end
    emission_log_probability = _h1_path_emission_log_probability(fixture, path)
    if !isfinite(emission_log_probability)
        return -Inf
    end
    return total + emission_log_probability
end

function _h1_path_emission_log_probability(fixture::H1SyntheticFixture, path::Vector{String})::Float64
    if fixture.observation_classes === nothing
        return 0.0
    end
    nonterminal_path = path[2:(length(path) - 1)]
    if length(nonterminal_path) != length(fixture.observation_classes)
        return -Inf
    end
    total = 0.0
    for (vertex, observation_class) in zip(nonterminal_path, fixture.observation_classes)
        vertex_class = get(fixture.vertex_classes, vertex, nothing)
        if vertex_class === nothing
            return -Inf
        end
        total += vertex_class == observation_class ? 0.0 : -20.0
    end
    return total
end

function _h1_path_failure_code(fixture::H1SyntheticFixture, path::Vector{String})::String
    if fixture.observation_classes !== nothing
        nonterminal_path = path[2:(length(path) - 1)]
        if length(nonterminal_path) != length(fixture.observation_classes)
            return "length_mismatch"
        end
    end
    if !isfinite(_h1_path_log_probability(fixture, path))
        return "length_mismatch"
    end
    return "none"
end

function _h1_outgoing_edges(fixture::H1SyntheticFixture)::Dict{String, Vector{H1SyntheticEdge}}
    outgoing = Dict{String, Vector{H1SyntheticEdge}}()
    for edge in fixture.edges
        if !haskey(outgoing, edge.source)
            outgoing[edge.source] = H1SyntheticEdge[]
        end
        push!(outgoing[edge.source], edge)
    end
    return outgoing
end

function _h1_path_tie_key(fixture::H1SyntheticFixture, path::Vector{String})::String
    edge_probability_by_key = Dict(_h1_edge_key(edge) => edge.probability for edge in fixture.edges)
    edge_keys = String[]
    for index in 1:(length(path) - 1)
        edge_key = "$(path[index])->$(path[index + 1])"
        if !haskey(edge_probability_by_key, edge_key)
            return join(path, "|")
        end
        push!(edge_keys, "$(path[index + 1])|forward|$(edge_key)")
    end
    return join(edge_keys, "|")
end

function _h1_edge_tie_key(edge::H1SyntheticEdge)::String
    return "$(edge.target)|forward|$(_h1_edge_key(edge))"
end

function _h1_edge_key(edge::H1SyntheticEdge)::String
    return "$(edge.source)->$(edge.target)"
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
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --slice H1 --execute")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --slice H1 --execute --write-artifacts --output-dir benchmarking/results/h1_viterbi_dp_greedy_smoke")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --slice H1 --slice H7 --scale full")
    println("  julia --project=. benchmarking/rhizomorph_benchmark_harness.jl --plan --scale ci --write-artifacts --output-dir results/public-record")
    println()
    println("Scales: ci, full, candidate")
    println("Execution currently supports the H1-G0 through H1-G4 Viterbi-DP-vs-greedy synthetic smoke.")
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
        "--run-id"
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
    if "--write-artifacts" in args
        default_output_dir = dry_run ?
                             joinpath("results", "public-record") :
                             joinpath("benchmarking", "results", "h1_viterbi_dp_greedy_smoke")
        output_dir = _flag_value(args, "--output-dir", default_output_dir)
        run_id = _flag_value(
            args,
            "--run-id",
            dry_run ? "rhizomorph_benchmark_plan" : "h1_viterbi_dp_greedy_smoke"
        )
        artifacts = if dry_run
            write_rhizomorph_benchmark_plan_artifacts(
                output_dir = output_dir,
                scale = scale,
                hypothesis_ids = isempty(hypothesis_ids) ? nothing : hypothesis_ids,
                dataset_ids = isempty(dataset_ids) ? nothing : dataset_ids,
                manifest_path = manifest_path,
                run_id = run_id,
                command_args = args
            )
        else
            write_h1_viterbi_dp_greedy_artifacts(
                output_dir = output_dir,
                run_id = run_id,
                command_args = args
            )
        end
        println("Artifacts written to: $(artifacts.root)")
        println("Artifact index: $(artifacts.index)")
        println("Run provenance: $(artifacts.provenance)")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
