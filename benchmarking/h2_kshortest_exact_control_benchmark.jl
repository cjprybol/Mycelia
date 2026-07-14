# H2 exact-control benchmark: K-shortest-paths alternative-assembly discovery
# on two fully synthetic, hand-constructed toy graphs (a bubble graph and a
# repeat-choice graph), per the "exact-control tier" of
# rhizomorph-paper/planning/PLAN-2026-04-28-h2-kshortest-alternative-assemblies.md.
#
# Scope: ONLY the exact-control tier. The realistic diploid (SK1/S288C HiFi)
# and repeat-rich (phage T4) tiers named in the plan are explicitly out of
# scope for this script and are NOT run here.
#
# Usage:
#   julia --project=. benchmarking/h2_kshortest_exact_control_benchmark.jl
#   julia --project=. benchmarking/h2_kshortest_exact_control_benchmark.jl --output-dir /tmp/h2

import Mycelia
import MetaGraphsNext
import DataFrames
import Random

include(joinpath(@__DIR__, "benchmark_artifacts.jl"))

const H2_K_VALUES = [2, 5, 10, 50, 100]
const H2_PRIMARY_K = 10
const H2_BASELINE_SEEDS = [42, 123, 456]
const H2_WALKS_PER_SEED = 1000
const H2_WALK_MAX_STEPS = 12
const H2_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "h2_kshortest_exact_control_20260714")

const Rhizomorph = Mycelia.Rhizomorph

# ============================================================================
# Graph construction
# ============================================================================

# Direct MetaGraphsNext construction (not build_kmer_graph/FASTA-based), so
# named bubble/repeat arms can be hand-placed with exact edge weights. Follows
# the idiom in test/4_assembly/local_path_enumeration_test.jl:6-33.
function _new_named_graph()::MetaGraphsNext.MetaGraph
    return MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph();
        label_type = String,
        vertex_data_type = Any,
        edge_data_type = Rhizomorph.StrandWeightedEdgeData,
        weight_function = Rhizomorph.edge_data_weight,
        default_weight = 0.0
    )
end

function _add_edges!(
        graph::MetaGraphsNext.MetaGraph,
        edges::Vector{Tuple{String, String, Float64}},
        weighted::Bool
)::Nothing
    vertices = Set{String}()
    for (src, dst, _) in edges
        push!(vertices, src)
        push!(vertices, dst)
    end
    for vertex in vertices
        graph[vertex] = nothing
    end
    for (src, dst, weight) in edges
        effective_weight = weighted ? weight : 1.0
        graph[src, dst] = Rhizomorph.StrandWeightedEdgeData(
            effective_weight, Rhizomorph.Forward, Rhizomorph.Forward)
    end
    return nothing
end

# --- Bubble graph -----------------------------------------------------------
#
# Shared prefix S->P1->P2 (10.0/edge), a three-arm bubble at P2 (A1 highest at
# 9.0, A2 second at 6.0, decoy D1 lowest at 1.0), shared suffix Q1->Q2->T
# (10.0/edge). True answers: full A1 and A2 paths. Decoy: full D1 path.
function h2_build_bubble_graph(; weighted::Bool)::MetaGraphsNext.MetaGraph
    graph = _new_named_graph()
    edges = Tuple{String, String, Float64}[
        ("S", "P1", 10.0),
        ("P1", "P2", 10.0),
        ("P2", "A1_1", 9.0), ("A1_1", "A1_2", 9.0), ("A1_2", "Q1", 9.0),
        ("P2", "A2_1", 6.0), ("A2_1", "A2_2", 6.0), ("A2_2", "Q1", 6.0),
        ("P2", "D1_1", 1.0), ("D1_1", "D1_2", 1.0), ("D1_2", "Q1", 1.0),
        ("Q1", "Q2", 10.0),
        ("Q2", "T", 10.0)
    ]
    _add_edges!(graph, edges, weighted)
    return graph
end

# --- Repeat-choice graph -----------------------------------------------------
#
# Unique left flank S->F1 (10.0), two parallel VERTEX-DISJOINT acyclic
# repeat-resolution arms F1->F2 (R2 at 8.0/edge, R3 at 7.0/edge), plus a
# lower-weight error arm E1 (1.0/edge), unique right flank F2->T (10.0).
# k_shortest_paths (Yen's algorithm over Dijkstra shortest-path-excluding) is
# loopless on this DAG by construction: each arm is a disjoint acyclic chain,
# so no shared-vertex revisitation is required or possible.
function h2_build_repeat_choice_graph(; weighted::Bool)::MetaGraphsNext.MetaGraph
    graph = _new_named_graph()
    edges = Tuple{String, String, Float64}[
        ("S", "F1", 10.0),
        ("F1", "R2_1", 8.0), ("R2_1", "R2_2", 8.0), ("R2_2", "F2", 8.0),
        ("F1", "R3_1", 7.0), ("R3_1", "R3_2", 7.0), ("R3_2", "F2", 7.0),
        ("F1", "E1_1", 1.0), ("E1_1", "E1_2", 1.0), ("E1_2", "F2", 1.0),
        ("F2", "T", 10.0)
    ]
    _add_edges!(graph, edges, weighted)
    return graph
end

# Ground-truth arm definitions: label => the two arm-interior vertices whose
# joint presence in a path's vertex sequence uniquely identifies that arm.
# Because every arm is vertex-disjoint from every other arm by construction,
# substring/subset matching on these interior vertices is exact and sufficient
# (no SHA-256 / alignment-based truth matching needed for these toy graphs).
const H2_BUBBLE_TRUE_LABELS = ["A1", "A2"]
const H2_BUBBLE_ARM_VERTICES = Dict(
    "A1" => ["A1_1", "A1_2"],
    "A2" => ["A2_1", "A2_2"],
    "decoy" => ["D1_1", "D1_2"]
)
const H2_REPEAT_TRUE_LABELS = ["R2", "R3"]
const H2_REPEAT_ARM_VERTICES = Dict(
    "R2" => ["R2_1", "R2_2"],
    "R3" => ["R3_1", "R3_2"],
    "decoy" => ["E1_1", "E1_2"]
)

function _graph_spec(graph_id::AbstractString)
    if graph_id == "bubble"
        return (
            build = h2_build_bubble_graph,
            true_labels = H2_BUBBLE_TRUE_LABELS,
            arm_vertices = H2_BUBBLE_ARM_VERTICES,
            source = "S",
            target = "T"
        )
    elseif graph_id == "repeat_choice"
        return (
            build = h2_build_repeat_choice_graph,
            true_labels = H2_REPEAT_TRUE_LABELS,
            arm_vertices = H2_REPEAT_ARM_VERTICES,
            source = "S",
            target = "T"
        )
    end
    throw(ArgumentError("unknown graph_id: $graph_id"))
end

function _truth_label(path_vertices::Vector{String}, arm_vertices::Dict{String, Vector{String}})::String
    vertex_set = Set(path_vertices)
    for (label, arm_verts) in arm_vertices
        if all(v -> v in vertex_set, arm_verts)
            return label
        end
    end
    return "unlabeled"
end

# ============================================================================
# K-shortest-paths enumeration, re-derived scoring, dedup, ranking
# ============================================================================

struct H2CandidateRow
    rank::Int
    path_vertices::Vector{String}
    group_id::String
    score::Float64
    total_probability::Float64
    truth_label::String
    is_true_alternative::Bool
    is_decoy::Bool
    tie_broken::Bool
end

function _dedup_and_rank(
        raw_paths::Vector{<:Rhizomorph.GraphPath},
        true_labels::Vector{String},
        arm_vertices::Dict{String, Vector{String}}
)::Vector{H2CandidateRow}
    # Stage: exact string-group dedup keyed on the joined vertex-label
    # sequence. A full SHA-256 / near-duplicate-clustering pipeline (as the
    # manuscript plan's general deduplication section describes) is
    # unnecessary complexity for these two graphs: they are constructed so
    # that true and decoy arms are vertex-disjoint by design, so exact joined-
    # string equality is already a complete and correct dedup key here. This
    # is an explicit scope simplification relative to the full plan.
    seen = Dict{String, Int}()  # group_id => index into `entries`
    entries = NamedTuple[]
    for (orig_index, path) in enumerate(raw_paths)
        vertices = [step.vertex_label for step in path.steps]
        group_id = join(vertices, ",")
        if !haskey(seen, group_id)
            seen[group_id] = length(entries) + 1
            push!(entries, (
                orig_index = orig_index,
                path_vertices = vertices,
                group_id = group_id,
                total_probability = path.total_probability
            ))
        end
        # else: later (lower-ranked, since raw_paths is already best-first)
        # duplicate is dropped; the representative keeps the earliest
        # (highest-ranked) occurrence.
    end

    # Explicit re-derived ranking score, independent of k_shortest_paths's
    # internal ordering. score = -log(total_probability); lower is better.
    scored = [
        (
            score = -log(entry.total_probability),
            length = length(entry.path_vertices),
            group_id = entry.group_id,
            orig_index = entry.orig_index,
            path_vertices = entry.path_vertices,
            total_probability = entry.total_probability
        )
        for entry in entries
    ]

    # Tie-break: (score, path length ascending, joined-vertex string
    # lexicographic, original enumeration index) rather than trusting
    # k_shortest_paths's internal tie order.
    sort!(scored; by = e -> (e.score, e.length, e.group_id, e.orig_index))

    # Detect ties (within floating-point tolerance) for the tie_broken flag.
    tie_broken_flags = falses(length(scored))
    for i in eachindex(scored)
        for j in eachindex(scored)
            i == j && continue
            if isapprox(scored[i].score, scored[j].score; atol = 1e-12, rtol = 1e-9)
                tie_broken_flags[i] = true
            end
        end
    end

    rows = H2CandidateRow[]
    for (rank, entry) in enumerate(scored)
        truth_label = _truth_label(entry.path_vertices, arm_vertices)
        is_true_alternative = truth_label in true_labels
        is_decoy = truth_label == "decoy"
        push!(rows, H2CandidateRow(
            rank, entry.path_vertices, entry.group_id, entry.score,
            entry.total_probability, truth_label, is_true_alternative,
            is_decoy, tie_broken_flags[rank]
        ))
    end
    return rows
end

function _current_peak_rss_mib()::Float64
    return Sys.maxrss() / 1024^2
end

# ============================================================================
# Random-walk baseline
# ============================================================================
#
# IMPORTANT (empirically discovered while building this benchmark):
# `Rhizomorph.probabilistic_walk_next` internally computes correctly
# NORMALIZED transition probabilities (weight / total-outgoing-weight) to
# SAMPLE which edge to take at each step, but then stores the RAW, un-
# normalized edge weight (not the normalized probability) in each
# `WalkStep.probability`/`cumulative_probability` field. Consequently the
# returned `GraphPath.total_probability` is not a valid probability at all
# (e.g. ~7.29e6 for a 7-edge walk on this benchmark's bubble graph, versus the
# correct ~0.5625 that `k_shortest_paths` reports for the identical vertex
# sequence). This is a pre-existing library behavior, not something this
# benchmark script modifies; fixing `probabilistic_walk_next` is out of scope
# here. To keep the H2 decision rule's baseline comparison meaningful, this
# benchmark independently recomputes each walk's path probability from its
# realized vertex sequence using the SAME weight/total-outgoing-weight rule
# `k_shortest_paths` uses internally (see `_build_graph_path_from_vertices` /
# `_shortest_path_excluding` in path-finding.jl), rather than trusting
# `walk.total_probability`.
function _h2_recompute_path_probability(
        graph::MetaGraphsNext.MetaGraph,
        path_vertices::Vector{String}
)::Float64
    total_probability = 1.0
    for i in 1:(length(path_vertices) - 1)
        src = path_vertices[i]
        dst = path_vertices[i + 1]
        total_out = sum(
            Rhizomorph.edge_data_weight(graph[edge_src, edge_dst])
            for (edge_src, edge_dst) in MetaGraphsNext.edge_labels(graph)
            if edge_src == src;
            init = 0.0
        )
        edge_weight = Rhizomorph.edge_data_weight(graph[src, dst])
        total_probability *= edge_weight / total_out
    end
    return total_probability
end

struct H2BaselineRow
    graph_id::String
    seed::Int
    n_walks::Int
    n_walks_reaching_sink::Int
    best_walk_recovers_true_alt::Bool
    best_walk_score::Float64
end

function _walk_matches_true_alt(
        path_vertices::Vector{String},
        true_labels::Vector{String},
        arm_vertices::Dict{String, Vector{String}}
)::Bool
    return _truth_label(path_vertices, arm_vertices) in true_labels
end

function h2_random_walk_baseline(
        graph::MetaGraphsNext.MetaGraph,
        graph_id::AbstractString,
        source::AbstractString,
        target::AbstractString,
        true_labels::Vector{String},
        arm_vertices::Dict{String, Vector{String}}
)::Vector{H2BaselineRow}
    rows = H2BaselineRow[]
    for seed in H2_BASELINE_SEEDS
        n_reaching_sink = 0
        best_true_alt_score = Inf
        best_overall_score = Inf
        recovers_true_alt = false

        # probabilistic_walk_next calls Mycelia.Random.seed!(seed) internally,
        # mutating GLOBAL RNG state. Seed once per seed-block, then draw all
        # H2_WALKS_PER_SEED walks sequentially from the now-deterministic
        # global stream (do NOT reseed on every call, and never run these
        # concurrently).
        for walk_index in 1:H2_WALKS_PER_SEED
            walk = if walk_index == 1
                Rhizomorph.probabilistic_walk_next(
                    graph, source, H2_WALK_MAX_STEPS; seed = seed)
            else
                Rhizomorph.probabilistic_walk_next(
                    graph, source, H2_WALK_MAX_STEPS)
            end
            final_vertex = walk.steps[end].vertex_label
            if final_vertex == target
                n_reaching_sink += 1
                path_vertices = [step.vertex_label for step in walk.steps]
                # Do NOT use walk.total_probability here — see the module-
                # level note above _h2_recompute_path_probability: it is a
                # raw-weight product, not a normalized probability.
                recomputed_probability = _h2_recompute_path_probability(graph, path_vertices)
                score = -log(recomputed_probability)
                if score < best_overall_score
                    best_overall_score = score
                end
                if _walk_matches_true_alt(path_vertices, true_labels, arm_vertices)
                    recovers_true_alt = true
                    if score < best_true_alt_score
                        best_true_alt_score = score
                    end
                end
            end
        end

        # best_walk_score: the score of the best true-alt-matching walk when
        # one was recovered (feeds the H2 decision rule's "outrank baseline"
        # sub-condition directly); otherwise falls back to the best
        # sink-reaching walk's score as a diagnostic (no true-alt walk exists
        # to outrank, so the decision-rule sub-condition is trivially
        # satisfied regardless of this fallback value). Inf sentinel when no
        # walk reached the sink at all.
        best_score = recovers_true_alt ? best_true_alt_score : best_overall_score

        push!(rows, H2BaselineRow(
            string(graph_id), seed, H2_WALKS_PER_SEED, n_reaching_sink,
            recovers_true_alt, best_score
        ))
    end
    return rows
end

# ============================================================================
# Main per-graph, per-weighting, per-K sweep
# ============================================================================

function h2_run_kshortest_sweep()::NamedTuple
    graph_ids = ["bubble", "repeat_choice"]
    weightings = ["evidence_weighted", "uniform"]

    path_metric_rows = NamedTuple[]
    summary_rows = NamedTuple[]
    baseline_rows = NamedTuple[]

    for graph_id in graph_ids
        spec = _graph_spec(graph_id)

        # Baseline is evidence-weighted only, per the task spec.
        baseline_graph = spec.build(weighted = true)
        append!(
            baseline_rows,
            [
                (
                    graph_id = row.graph_id,
                    seed = row.seed,
                    n_walks = row.n_walks,
                    n_walks_reaching_sink = row.n_walks_reaching_sink,
                    best_walk_recovers_true_alt = row.best_walk_recovers_true_alt,
                    best_walk_score = row.best_walk_score
                )
                for row in h2_random_walk_baseline(
                    baseline_graph, graph_id, spec.source, spec.target,
                    spec.true_labels, spec.arm_vertices)
            ]
        )
        baseline_best_true_alt_score = let
            recovered_scores = [
                row.best_walk_score for row in baseline_rows
                if row.graph_id == graph_id && row.best_walk_recovers_true_alt
            ]
            isempty(recovered_scores) ? nothing : minimum(recovered_scores)
        end

        for weighting in weightings
            weighted = weighting == "evidence_weighted"
            graph = spec.build(weighted = weighted)

            for K in H2_K_VALUES
                elapsed = @elapsed begin
                    raw_paths = Rhizomorph.k_shortest_paths(
                        graph, spec.source, spec.target, K)
                end
                peak_rss = _current_peak_rss_mib()

                candidates = _dedup_and_rank(raw_paths, spec.true_labels, spec.arm_vertices)
                top_score = isempty(candidates) ? NaN : first(candidates).score
                top_prob = isempty(candidates) ? NaN : first(candidates).total_probability

                for candidate in candidates
                    push!(path_metric_rows, (
                        graph_id = graph_id,
                        weighting = weighting,
                        K = K,
                        rank = candidate.rank,
                        path_vertices = join(candidate.path_vertices, ","),
                        sha_or_string_group_id = candidate.group_id,
                        score = candidate.score,
                        total_probability = candidate.total_probability,
                        delta_log_prob_from_top = candidate.score - top_score,
                        probability_ratio_to_top = candidate.total_probability / top_prob,
                        truth_label = candidate.truth_label,
                        is_true_alternative = candidate.is_true_alternative,
                        is_decoy = candidate.is_decoy,
                        tie_broken = candidate.tie_broken,
                        runtime_s = elapsed,
                        peak_rss_mib = peak_rss
                    ))
                end

                # --- Summary metrics for this (graph_id, weighting, K) -----
                first_true_rank = Dict{String, Union{Missing, Int}}(
                    label => missing for label in ["A1", "A2", "R2", "R3"])
                for label in spec.true_labels
                    match = findfirst(c -> c.truth_label == label, candidates)
                    first_true_rank[label] = isnothing(match) ? missing : candidates[match].rank
                end

                n_recovered = count(label -> !ismissing(first_true_rank[label]), spec.true_labels)
                alt_recall_at_k = n_recovered / length(spec.true_labels)

                reciprocal_ranks = [
                    ismissing(first_true_rank[label]) ? 0.0 : 1.0 / first_true_rank[label]
                    for label in spec.true_labels
                ]
                mrr = sum(reciprocal_ranks) / length(reciprocal_ranks)

                is_primary_decision_row = (K == H2_PRIMARY_K) && weighted

                decision_rule_pass = false
                if is_primary_decision_row
                    all_recovered = n_recovered == length(spec.true_labels)
                    true_alt_ranks = [c.rank for c in candidates if c.is_true_alternative]
                    decoy_ranks = [c.rank for c in candidates if c.is_decoy]
                    ranked_above_decoys = isempty(true_alt_ranks) ? false :
                                          (isempty(decoy_ranks) || maximum(true_alt_ranks) < minimum(decoy_ranks))

                    best_true_alt_score = isempty(true_alt_ranks) ? NaN :
                        minimum(c.score for c in candidates if c.is_true_alternative)
                    outranks_baseline = if baseline_best_true_alt_score === nothing
                        true  # baseline never recovered a true alt: trivially outranked
                    else
                        !isempty(true_alt_ranks) && best_true_alt_score <= baseline_best_true_alt_score
                    end

                    decision_rule_pass = all_recovered && ranked_above_decoys && outranks_baseline
                end

                push!(summary_rows, (
                    graph_id = graph_id,
                    weighting = weighting,
                    K = K,
                    n_raw_paths = length(raw_paths),
                    n_dedup_candidates = length(candidates),
                    alt_recall_at_k = alt_recall_at_k,
                    first_true_rank_A1 = first_true_rank["A1"],
                    first_true_rank_A2 = first_true_rank["A2"],
                    first_true_rank_R2 = first_true_rank["R2"],
                    first_true_rank_R3 = first_true_rank["R3"],
                    mrr = mrr,
                    is_primary_decision_row = is_primary_decision_row,
                    decision_rule_pass = decision_rule_pass,
                    runtime_s = elapsed,
                    peak_rss_mib = peak_rss
                ))
            end
        end
    end

    return (
        path_metrics = DataFrames.DataFrame(path_metric_rows),
        summary = DataFrames.DataFrame(summary_rows),
        baseline = DataFrames.DataFrame(baseline_rows)
    )
end

# ============================================================================
# Artifact writing
# ============================================================================

function _write_summary_markdown(
        summary::DataFrames.DataFrame,
        baseline::DataFrames.DataFrame,
        output_dir::AbstractString
)::String
    primary = summary[summary.is_primary_decision_row .== true, :]
    lines = String[]
    push!(lines, "# H2 K-shortest-paths exact-control benchmark summary")
    push!(lines, "")
    push!(lines, "Scope: exact-control tier ONLY (bubble graph + repeat-choice graph),")
    push!(lines, "from PLAN-2026-04-28-h2-kshortest-alternative-assemblies.md. The")
    push!(lines, "realistic diploid (SK1/S288C HiFi) and repeat-rich (phage T4) tiers")
    push!(lines, "named in that plan were NOT run in this script.")
    push!(lines, "")
    push!(lines, "## Decision-rule verdicts (K=10, evidence-weighted)")
    push!(lines, "")
    for row in eachrow(primary)
        verdict = row.decision_rule_pass ? "PASS" : "FAIL"
        push!(lines, "- **$(row.graph_id)**: $(verdict) (alt_recall_at_k=$(row.alt_recall_at_k), mrr=$(round(row.mrr, digits = 4)))")
    end
    push!(lines, "")
    push!(lines, "## Random-walk baseline ($(H2_WALKS_PER_SEED) walks/seed, seeds $(H2_BASELINE_SEEDS))")
    push!(lines, "")
    for row in eachrow(baseline)
        push!(lines, "- $(row.graph_id) seed=$(row.seed): $(row.n_walks_reaching_sink)/$(row.n_walks) walks reached sink; recovers_true_alt=$(row.best_walk_recovers_true_alt); best_walk_score=$(round(row.best_walk_score, digits = 4))")
    end
    push!(lines, "")
    push!(lines, "## Scope caveats")
    push!(lines, "")
    push!(lines, "- Exact-control tier only; realistic diploid and repeat-rich tiers are unrun and out of scope for this script.")
    push!(lines, "- `Rhizomorph.probabilistic_walk_next` records the raw (un-normalized) edge weight in each `WalkStep`'s probability/cumulative_probability fields instead of the normalized transition probability it correctly uses internally for sampling; its `GraphPath.total_probability` is therefore not a valid probability (observed ~7.29e6 for a 7-edge bubble-graph walk instead of ~0.5625). This benchmark does not modify that library function; instead it independently recomputes each baseline walk's path probability from the realized vertex sequence using the same weight/total-outgoing-weight rule `k_shortest_paths` uses, so the baseline comparison in the decision rule is valid.")
    push!(lines, "- Deduplication uses exact joined-vertex-label string matching, not SHA-256 hashing plus 99.5%-identity near-duplicate clustering as in the full plan's general deduplication section. This is sufficient and exact here because both graphs are constructed so that true and decoy arms are vertex-disjoint by design.")
    push!(lines, "- Truth-matching uses exact vertex-subset matching against known arm-interior vertices, not sequence alignment (appropriate only for these hand-constructed toy graphs, per the plan's exact-control truth-matching plan).")
    push!(lines, "- `decision_rule_pass` is only meaningful (and only computed as non-trivially-false) on the K=10, evidence-weighted row per graph; see `is_primary_decision_row`.")

    path = joinpath(output_dir, "h2_kshortest_exact_control_summary.md")
    mkpath(output_dir)
    open(path, "w") do io
        println(io, join(lines, "\n"))
    end
    return path
end

function main(args::Vector{String} = ARGS)::Nothing
    output_dir = H2_DEFAULT_OUTPUT_DIR
    idx = findfirst(==("--output-dir"), args)
    if idx !== nothing && idx < length(args)
        output_dir = args[idx + 1]
    end

    results = h2_run_kshortest_sweep()

    artifacts = write_benchmark_artifacts(
        [
            "h2_kshortest_exact_control_path_metrics" => results.path_metrics,
            "h2_kshortest_exact_control_summary" => results.summary,
            "h2_kshortest_exact_control_baseline" => results.baseline
        ];
        output_dir = output_dir,
        run_id = "h2_kshortest_exact_control_20260714",
        scale = "local-smoke",
        dataset_ids = ["h2_bubble_graph_synthetic", "h2_repeat_choice_graph_synthetic"],
        command_args = ["julia", "--project=.",
            "benchmarking/h2_kshortest_exact_control_benchmark.jl"],
        metadata = Dict{String, Any}(
            "plan_path" => "rhizomorph-paper/planning/PLAN-2026-04-28-h2-kshortest-alternative-assemblies.md",
            "scope" => "exact-control tier ONLY (bubble graph + repeat-choice graph); realistic diploid and repeat-rich tiers are out of scope and unrun",
            "k_values" => H2_K_VALUES,
            "primary_k" => H2_PRIMARY_K,
            "baseline_seeds" => H2_BASELINE_SEEDS,
            "walks_per_seed" => H2_WALKS_PER_SEED,
            "walk_max_steps" => H2_WALK_MAX_STEPS,
            "dedup_simplification" => "exact joined-vertex-label string match only (no SHA-256 / near-duplicate clustering); valid here because true/decoy arms are vertex-disjoint by construction",
            "probabilistic_walk_next_finding" => "GraphPath.total_probability returned by Rhizomorph.probabilistic_walk_next is NOT a valid probability (it stores the raw un-normalized edge-weight product, not the normalized transition-probability product used internally for sampling); this benchmark recomputes baseline walk scores independently via _h2_recompute_path_probability rather than trusting that field",
            "peak_rss_mib_semantics" => "peak_rss_mib records process-level Sys.maxrss() sampled once per (graph_id, weighting, K) k_shortest_paths call; it is process-wide provenance, not per-algorithm incremental memory (same convention as the H1 Viterbi/DP/greedy smoke harness)",
            "score_definition" => "score = -log(total_probability); re-derived and re-sorted independently of k_shortest_paths's internal ordering, with tie-break (path length asc, joined-vertex string lexicographic, original enumeration index)",
            "decision_rule" => "at K=10, evidence-weighted: ALL true alternatives recovered in dedup top-K AND ranked strictly above every decoy candidate AND outrank the random-walk baseline's best recovered score for a true alternative (trivially satisfied if baseline never recovers one)"
        ),
        table_context_columns = Dict{String, Dict{String, String}}(
            "h2_kshortest_exact_control_path_metrics" => Dict("benchmark_graph_id" => "graph_id"),
            "h2_kshortest_exact_control_summary" => Dict("benchmark_graph_id" => "graph_id"),
            "h2_kshortest_exact_control_baseline" => Dict("benchmark_graph_id" => "graph_id")
        )
    )

    markdown_path = _write_summary_markdown(results.summary, results.baseline, output_dir)

    println("Wrote H2 K-shortest exact-control benchmark artifacts:")
    println("  root: $(artifacts.root)")
    println("  path_metrics: $(artifacts.tables["h2_kshortest_exact_control_path_metrics"].table)")
    println("  summary: $(artifacts.tables["h2_kshortest_exact_control_summary"].table)")
    println("  baseline: $(artifacts.tables["h2_kshortest_exact_control_baseline"].table)")
    println("  summary_markdown: $markdown_path")

    primary = results.summary[results.summary.is_primary_decision_row .== true, :]
    for row in eachrow(primary)
        println("  decision_rule_pass[$(row.graph_id)] = $(row.decision_rule_pass)")
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
