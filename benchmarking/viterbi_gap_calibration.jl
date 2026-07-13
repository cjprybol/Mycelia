# Tier-2 per-base Viterbi correction-confidence calibration (td-21eg).
#
# `collect_gap_truth` keeps a focused decoder-level contract test, while
# `run_gap_calibration` learns exclusively from telemetry emitted by the same
# iterative-corrector configuration used by the frontier. The logistic model is
# fit on proposed substitutions only and evaluated on a deterministic,
# read-grouped holdout so bases from one read never leak across the split.

import Mycelia
import Random
import FASTX

include(joinpath(@__DIR__, "calibration_metrics.jl"))
include(joinpath(@__DIR__, "benchmark_kmer_graph_fixture.jl"))

const CORRECTION_CONFIDENCE_FEATURES = (
    :raw_gap,
    :min_kmer_support,
    :all_stage0_solid,
    :competing_branch_support_ratio,
    :collapsed_frontier
)

function correction_calibration_serving_config(k::Int)::NamedTuple
    return (
        max_k = max(k, 13),
        skip_solid = false,
        graph_mode = :doublestrand,
        n_k_rungs = 3,
        max_iterations_per_k = 2,
        hard_window = false,
        soft_em = true,
        cheap_correct = true,
        beam_width = nothing
    )
end

function _contract_skip_fraction(contract_skips::Int, read_passes::Int,
        maximum_fraction::Float64)::Float64
    contract_skips >= 0 || throw(ArgumentError("contract_skips must be non-negative"))
    read_passes >= 0 || throw(ArgumentError("read_passes must be non-negative"))
    0.0 <= maximum_fraction <= 1.0 ||
        throw(ArgumentError("maximum contract-skip fraction must be in [0, 1]"))
    contract_skips <= read_passes ||
        throw(ArgumentError("contract skips cannot exceed read passes"))
    fraction = read_passes == 0 ? 0.0 : contract_skips / read_passes
    fraction <= maximum_fraction ||
        throw(ArgumentError("feature-contract skip fraction $(fraction) exceeds " *
                            "configured bound $(maximum_fraction)"))
    return fraction
end

"""
    collect_gap_truth(fixture, observed_sequence, truth_sequence;
        config, solid_kmers, quality_scores)

Decode one fixed-length observation and align `position_gaps[i]` with
`path.steps[i + 1]`. A k-mer transition contributes the newly emitted base at
read position `i + k`; its label is true iff that corrected base matches truth.
For each candidate base, emit the fixed feature vector
`(raw_gap, min_kmer_support, all_stage0_solid, competing_branch_support_ratio,
collapsed_frontier)` plus a `candidate_edit` flag. The support/solidity window is
`path.steps[i + 1:min(i + k, nsteps)]`, the decoded k-mers overlapping that
base. Collapsed-frontier `Inf` sentinels are retained as an explicit stratum.
"""
function collect_gap_truth(fixture::NamedTuple,
        observed_sequence::AbstractString,
        truth_sequence::AbstractString;
        config::Mycelia.ViterbiCorrectionConfig =
        Mycelia.ViterbiCorrectionConfig(record_position_gaps = true),
        solid_kmers::Union{Nothing, AbstractSet} = nothing,
        quality_scores::Union{Nothing, AbstractVector{<:Integer}} = nothing)::NamedTuple
    length(observed_sequence) == length(truth_sequence) ||
        throw(ArgumentError("observed and truth sequences must have equal length"))
    raw_observation = fixture.to_observation(observed_sequence)
    observation = if quality_scores === nothing
        raw_observation
    else
        length(quality_scores) == length(observed_sequence) ||
            throw(ArgumentError("quality scores must align with the observed sequence"))
        observation_k = length(String(first(raw_observation)))
        [Mycelia.QualityObservation(raw_observation[i],
             UInt8.(@view quality_scores[i:(i + observation_k - 1)]))
         for i in eachindex(raw_observation)]
    end
    result = Mycelia.correct_observations(fixture.graph, [observation]; config = config)
    # NB: `something(x, throw(e))` would raise EAGERLY (throw is an argument,
    # evaluated before `something` runs) — use a lazy nothing-check so we only
    # raise when the decode genuinely failed.
    path = result.paths[1].path
    path === nothing && throw(ArgumentError("decode produced no path"))
    gaps = result.paths[1].diagnostics[:position_gaps]
    length(gaps) == length(path.steps) - 1 ||
        throw(ArgumentError("position-gap/path alignment contract violated"))
    k = length(String(path.steps[1].vertex_label))
    corrected_sequence = String(Mycelia.Rhizomorph.path_to_sequence(path, fixture.graph))
    length(corrected_sequence) == length(observed_sequence) ||
        throw(ArgumentError("decoded sequence must preserve observation length"))
    solid = solid_kmers === nothing ? Mycelia._solid_kmer_set(fixture.graph) : solid_kmers
    scores = Float64[]
    labels = Bool[]
    positions = Int[]
    candidate_edits = Bool[]
    feature_rows = NamedTuple[]
    for i in eachindex(gaps)
        (isfinite(gaps[i]) || gaps[i] == Inf) ||
            throw(ArgumentError("position gaps must be finite or +Inf"))
        position = i + k
        position <= lastindex(truth_sequence) || continue
        shared_features = Mycelia.correction_confidence_features(
            fixture.graph, path, i, k, solid, Float64(gaps[i]))
        features = (
            raw_gap = shared_features.raw_gap,
            min_kmer_support = shared_features.min_kmer_support,
            all_stage0_solid = shared_features.all_stage0_solid,
            competing_branch_support_ratio =
            shared_features.competing_branch_support_ratio,
            collapsed_frontier = shared_features.raw_gap == Inf
        )
        push!(scores, features.raw_gap)
        push!(labels, corrected_sequence[position] == truth_sequence[position])
        push!(positions, position)
        push!(candidate_edits,
            corrected_sequence[position] != observed_sequence[position])
        push!(feature_rows, features)
    end
    features = Matrix{Float64}(undef, length(feature_rows),
        length(CORRECTION_CONFIDENCE_FEATURES))
    for (row_index, row) in enumerate(feature_rows)
        features[row_index, :] .= (
            row.collapsed_frontier ? 0.0 : row.raw_gap,
            row.min_kmer_support,
            row.all_stage0_solid ? 1.0 : 0.0,
            row.competing_branch_support_ratio,
            row.collapsed_frontier ? 1.0 : 0.0
        )
    end
    return (; scores, labels, positions, candidate_edits, features, feature_rows, result)
end

# The three calibration models compared head-to-head. Each is a (fit, predict)
# pair over the same (gap, label) sample; the runtime gate consumes only a raw
# gap threshold, so these differ only in how they turn the gap RANKING into
# calibrated probabilities (the ECE/Brier axis).
const CALIBRATION_MODELS = (
    (name = "isotonic", fit = fit_isotonic_map, predict = predict_isotonic),
    (name = "logistic", fit = fit_logistic_map, predict = predict_logistic),
    (name = "binned", fit = fit_binned_map, predict = predict_binned)
)

"""Fit+evaluate one (fit, predict) calibration model on finite gap/truth pairs."""
function calibrate_gap_probability(fit::Function, predict::Function,
        scores::AbstractVector{<:Real}, labels::AbstractVector{Bool};
        nbins::Int = 10)::NamedTuple
    model = fit(scores, labels)
    probabilities = predict(model, scores)
    return (
        model = model,
        probabilities = probabilities,
        ece = expected_calibration_error(probabilities, labels; nbins = nbins),
        brier = brier_score(probabilities, labels),
        reliability = reliability_bins(probabilities, labels; nbins = nbins)
    )
end

"""
    compare_gap_calibrations(scores, labels; nbins=10) -> NamedTuple

Fit all three calibration models on the same finite (gap, label) sample. Returns
`(auroc, rows)` where `auroc` is the shared raw-gap discrimination (identical
across models) and each row is `(model, ece, brier, cal)`. ECE/Brier isolate
CALIBRATION quality — the only axis on which the models differ.
"""
function compare_gap_calibrations(scores::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; nbins::Int = 10)::NamedTuple
    rows = [(model = m.name,
                cal = calibrate_gap_probability(
                    m.fit, m.predict, scores, labels; nbins = nbins))
            for m in CALIBRATION_MODELS]
    return (auroc = auroc(scores, labels),
        rows = [(model = r.model, ece = r.cal.ece, brier = r.cal.brier, cal = r.cal)
                for r in rows])
end

"""
    grouped_read_split(rows; holdout_fraction=0.2, seed=42) -> NamedTuple

Split candidate rows by `(error_rate, read_id)`, stratified by error rate. The
split is deterministic under `seed`; every error rate with at least two read
groups contributes at least one held-out and one training group.
"""
function grouped_read_split(rows::AbstractVector{<:NamedTuple};
        holdout_fraction::Float64 = 0.2, seed::Int = 42)::NamedTuple
    0.0 < holdout_fraction < 1.0 ||
        throw(ArgumentError("holdout_fraction must be in the open interval (0, 1)"))
    rng = Random.MersenneTwister(seed)
    heldout_groups = Set{Tuple{Float64, String}}()
    training_groups = Set{Tuple{Float64, String}}()
    for error_rate in sort!(unique(Float64(row.error_rate) for row in rows))
        read_ids = sort!(unique(String(row.read_id)
        for row in rows
        if row.error_rate == error_rate))
        Random.shuffle!(rng, read_ids)
        if length(read_ids) < 2
            foreach(read_id -> push!(training_groups, (error_rate, read_id)), read_ids)
            continue
        end
        n_heldout = clamp(round(Int, holdout_fraction * length(read_ids)),
            1, length(read_ids) - 1)
        foreach(read_id -> push!(heldout_groups, (error_rate, read_id)),
            @view read_ids[1:n_heldout])
        foreach(read_id -> push!(training_groups, (error_rate, read_id)),
            @view read_ids[(n_heldout + 1):end])
    end
    training_indices = [i
                        for (i, row) in enumerate(rows)
                        if (
        Float64(row.error_rate), String(row.read_id)) in training_groups]
    heldout_indices = [i
                       for (i, row) in enumerate(rows)
                       if (Float64(row.error_rate), String(row.read_id)) in heldout_groups]
    return (; training_indices, heldout_indices, training_groups, heldout_groups)
end

function _feature_matrix(rows::AbstractVector{<:NamedTuple})::Matrix{Float64}
    features = Matrix{Float64}(undef, length(rows),
        length(CORRECTION_CONFIDENCE_FEATURES))
    for (i, row) in enumerate(rows)
        features[i, :] .= (
            row.collapsed_frontier ? 0.0 : row.raw_gap,
            row.min_kmer_support,
            row.all_stage0_solid ? 1.0 : 0.0,
            row.competing_branch_support_ratio,
            row.collapsed_frontier ? 1.0 : 0.0
        )
    end
    return features
end

function _calibration_partition(rows::AbstractVector{<:NamedTuple},
        indices::AbstractVector{<:Integer})::NamedTuple
    candidate_indices = [i for i in indices if rows[i].candidate_edit]
    candidate_rows = rows[candidate_indices]
    return (
        n = length(candidate_rows),
        rows = candidate_rows,
        indices = candidate_indices,
        features = _feature_matrix(candidate_rows),
        scores = Float64[row.raw_gap for row in candidate_rows],
        collapsed_frontier = Bool[row.collapsed_frontier for row in candidate_rows],
        # The gap-only baseline still uses only gap information. Its second
        # column encodes the `Inf` sentinel state so beam-collapsed candidates
        # are modeled rather than dropped or treated as certainty.
        gap_features = hcat(
            Float64[row.collapsed_frontier ? 0.0 : row.raw_gap
                    for row in candidate_rows],
            Float64[row.collapsed_frontier ? 1.0 : 0.0
                    for row in candidate_rows]),
        labels = Bool[row.label for row in candidate_rows],
        groups = Tuple{Float64, String}[(Float64(row.error_rate), String(row.read_id))
                                        for row in candidate_rows]
    )
end

function _metric_row(model_name::String, scope::String,
        error_rate::Union{Nothing, Float64}, probabilities::Vector{Float64},
        labels::Vector{Bool}; nbins::Int)::NamedTuple
    return (
        scope = scope,
        error_rate = error_rate,
        model = model_name,
        n = length(labels),
        positive_frac = isempty(labels) ? NaN : count(labels) / length(labels),
        auroc = auroc(probabilities, labels),
        ece = expected_calibration_error(probabilities, labels; nbins = nbins),
        brier = brier_score(probabilities, labels),
        reliability = reliability_bins(probabilities, labels; nbins = nbins)
    )
end

function _heldout_metric_rows(multifeature_model::NamedTuple, gap_model::NamedTuple,
        heldout::NamedTuple; nbins::Int)::Vector{NamedTuple}
    rows = NamedTuple[]
    scopes = Tuple{String, Union{Nothing, Float64}, Vector{Int}}[
    (
        "pooled", nothing, collect(eachindex(heldout.labels))),
]
    for error_rate in sort!(unique(row.error_rate for row in heldout.rows))
        indices = findall(row -> row.error_rate == error_rate, heldout.rows)
        push!(scopes, ("error-rate", Float64(error_rate), indices))
    end
    for (scope, error_rate, indices) in scopes
        labels = heldout.labels[indices]
        (all(labels) || !any(labels)) &&
            throw(ArgumentError(
                "held-out $(scope) scope needs both correctness classes"))
        multifeature_probabilities = predict_logistic(
            multifeature_model, heldout.features[indices, :])
        gap_probabilities = predict_logistic(
            gap_model, heldout.gap_features[indices, :])
        push!(rows,
            _metric_row("multifeature-logistic", scope, error_rate,
                multifeature_probabilities, labels; nbins = nbins))
        push!(rows,
            _metric_row("gap-only-logistic", scope, error_rate,
                gap_probabilities, labels; nbins = nbins))
    end
    return rows
end

# Substitute each base with prob `rate` to a DIFFERENT base (length-preserving),
# so the observed read aligns 1:1 with its truth for per-position labeling.
function _inject_substitutions(seq::AbstractString, rate::Float64, rng)::String
    chars = collect(seq)
    alts = Dict('A' => "CGT", 'C' => "AGT", 'G' => "ACT", 'T' => "ACG")
    for i in eachindex(chars)
        # Only substitute ACGT bases, and always to a DIFFERENT base (the per-base
        # `alts` string excludes the original), so a "substitution" is never a no-op
        # that would mislabel an unchanged position as edited.
        haskey(alts, chars[i]) || continue
        if rand(rng) < rate
            opts = alts[chars[i]]
            chars[i] = opts[rand(rng, 1:lastindex(opts))]
        end
    end
    return String(chars)
end

"""
    build_coverage_fixture(reads, k; moltype=:DNA) -> NamedTuple

Build a `collect_gap_truth`-compatible fixture (`.graph`, `.to_observation`) from
a READ SET rather than a single truth sequence. The read-coverage k-mer graph
carries the error- and repeat-induced BRANCHES the decoder must navigate — the
only place the per-position Viterbi margin is FINITE. A truth-only fixture graph
is linear (unique k-mers, no competitor), so every gap is `Inf` and there is no
calibration signal; the coverage graph is what the production gate actually
decodes against.
"""
function build_coverage_fixture(reads::Vector{<:FASTX.FASTQ.Record}, k::Int;
        moltype::Symbol = :DNA, graph_mode::Symbol = :doublestrand)::NamedTuple
    graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k;
        dataset_id = "gapcal", mode = graph_mode)
    to_observation = seq -> _bench_convert_kmers(_bench_sequence_kmers(seq, k), k, moltype)
    return (graph = graph, to_observation = to_observation)
end

"""
    run_gap_calibration(; kwargs...) -> Union{Vector, NamedTuple}

At each error rate (through err=0.10), simulate a substitution-only Q20 read
cohort and run it through the exact iterative-corrector configuration used by
the frontier probability family: doublestrand graphs, `skip_solid=false`,
`cheap_correct=true`, `soft_em=true`, `hard_window=false`, three k rungs, two
iterations per rung, and the automatic beam. Candidate-substitution telemetry is
collected through `correction_feature_sink`, so training and serving share one
feature extractor and the same graph/decoder path. A pooled multi-feature
logistic and a gap-only logistic baseline are fit on identical candidate-bearing
read groups and evaluated on a deterministic 20% grouped holdout, pooled and by
error rate. Prints a head-to-head table and writes
`results/rhizomorph_gap_calibration.csv`. Deterministic under `seed`.

By default, returns the metric-row vector for backward compatibility. With
`return_artifact=true`, returns the fitted models, partitions, and all finite
rows needed by the correction-frontier benchmark.
"""
function run_gap_calibration(;
        k::Int = 11,
        genome_length::Int = 1200,
        readlen::Int = 120,
        coverage::Int = 30,
        error_rates::Vector{Float64} = [0.05, 0.08, 0.10],
        seed::Int = 42,
        nbins::Int = 10,
        holdout_fraction::Float64 = 0.2,
        replicates::Int = 5,
        max_contract_skip_fraction::Float64 = 0.25,
        results_dir::String = joinpath(@__DIR__, "results"),
        return_artifact::Bool = false)::Union{Vector{NamedTuple}, NamedTuple}
    isempty(error_rates) && throw(ArgumentError("error_rates must not be empty"))
    all(error_rate -> 0.0 <= error_rate <= 1.0, error_rates) ||
        throw(ArgumentError("error_rates must lie in [0, 1]"))
    length(unique(error_rates)) == length(error_rates) ||
        throw(ArgumentError("error_rates must be unique"))
    genome_length >= readlen ||
        throw(ArgumentError("genome_length must be at least readlen"))
    replicates >= 1 || throw(ArgumentError("replicates must be at least 1"))
    0.0 <= max_contract_skip_fraction <= 1.0 ||
        throw(ArgumentError("max_contract_skip_fraction must be in [0, 1]"))
    rng = Random.MersenneTwister(seed)
    n_reads = max(1, ceil(Int, coverage * genome_length / readlen))
    serving_config = correction_calibration_serving_config(k)
    mkpath(results_dir)
    dataset = NamedTuple[]
    contract_skips = Dict{Float64, Int}()
    contract_read_passes = Dict{Float64, Int}()
    println("=== Rhizomorph multi-feature correction-confidence calibration ===")
    println("genome=$(genome_length)  k=$(k)  readlen=$(readlen)  coverage=$(coverage)x  " *
            "n_reads=$(n_reads)  replicates=$(replicates)  seed=$(seed)")
    for er in error_rates
        contract_skips[er] = 0
        contract_read_passes[er] = 0
        rate_event_count = 0
        rate_groups = Set{String}()
        for replicate in 1:replicates
            truth_genome = String(rand(rng, ['A', 'C', 'G', 'T'], genome_length))
            reads = FASTX.FASTQ.Record[]
            truth_by_id = Dict{String, String}()
            for i in 1:n_reads
                start = rand(rng, 1:(genome_length - readlen + 1))
                clean = truth_genome[start:(start + readlen - 1)]
                observed = _inject_substitutions(clean, er, rng)
                read_id = "rep$(replicate)-err$(er)-r$(i)"
                push!(reads,
                    FASTX.FASTQ.Record(read_id, observed, String(fill('5', readlen))))
                truth_by_id[read_id] = clean
            end
            event_rows = NamedTuple[]
            function feature_sink(
                    event::Mycelia.CorrectionFeatureObservation)::Nothing
                haskey(truth_by_id, event.read_id) ||
                    throw(ArgumentError(
                        "feature event has unknown read id $(event.read_id)"))
                truth = truth_by_id[event.read_id]
                1 <= event.position <= lastindex(truth) ||
                    throw(ArgumentError(
                        "feature-event position is outside its truth read"))
                features = event.features
                push!(event_rows,
                    (
                        error_rate = er,
                        read_id = event.read_id,
                        k = event.k,
                        position = event.position,
                        raw_gap = features.raw_gap,
                        min_kmer_support = features.min_kmer_support,
                        all_stage0_solid = features.all_stage0_solid,
                        competing_branch_support_ratio =
                        features.competing_branch_support_ratio,
                        collapsed_frontier = features.raw_gap == Inf,
                        label = event.corrected_base == truth[event.position],
                        candidate_edit = event.corrected_base != event.observed_base
                    ))
                return nothing
            end
            input_dir = mktempdir(prefix = "gapcal_input_")
            output_dir = mktempdir(prefix = "gapcal_output_")
            input_fastq = joinpath(input_dir, "calibration.fastq")
            open(FASTX.FASTQ.Writer, input_fastq) do writer
                for read in reads
                    write(writer, read)
                end
            end
            Random.seed!(seed + 1_000_000 * replicate + round(Int, er * 1_000_000))
            correction_result = Mycelia.mycelia_iterative_assemble(
                input_fastq;
                serving_config...,
                correction_feature_sink = feature_sink,
                verbose = false,
                enable_checkpointing = false,
                output_dir = output_dir)
            contract_skips[er] += Int(get(
                correction_result[:metadata], :calibrated_feature_contract_skips, 0))
            contract_read_passes[er] +=
                n_reads * Int(correction_result[:metadata][:total_iterations])
            append!(dataset, event_rows)
            rate_event_count += length(event_rows)
            union!(rate_groups, (row.read_id for row in event_rows))
        end
        println("err=$(er): candidate_events=$(rate_event_count)  " *
                "candidate-bearing_reads=$(length(rate_groups))  " *
                "feature_contract_skips=$(contract_skips[er])")
        try
            _contract_skip_fraction(
                contract_skips[er], contract_read_passes[er],
                max_contract_skip_fraction)
        catch error
            error isa ArgumentError || rethrow()
            throw(ArgumentError("error rate $(er) $(error.msg)"))
        end
    end

    isempty(dataset) &&
        throw(ArgumentError("no candidate-substitution events available"))
    all(row.candidate_edit for row in dataset) ||
        throw(ArgumentError("correction feature sink emitted a non-candidate event"))
    split = grouped_read_split(dataset;
        holdout_fraction = holdout_fraction, seed = seed + 1)
    for error_rate in error_rates
        candidate_groups = Set((row.error_rate, row.read_id)
        for row in dataset
        if row.error_rate == error_rate)
        length(candidate_groups) >= 2 ||
            throw(ArgumentError("error rate $(error_rate) needs at least two " *
                                "candidate-bearing read groups"))
        any(group -> group[1] == error_rate, split.training_groups) ||
            throw(ArgumentError("error rate $(error_rate) has no training groups"))
        any(group -> group[1] == error_rate, split.heldout_groups) ||
            throw(ArgumentError("error rate $(error_rate) has no held-out groups"))
    end
    training = _calibration_partition(dataset, split.training_indices)
    heldout = _calibration_partition(dataset, split.heldout_indices)
    isempty(training.labels) &&
        throw(ArgumentError("no candidate edits in training groups"))
    isempty(heldout.labels) && throw(ArgumentError("no candidate edits in held-out groups"))
    any(training.collapsed_frontier) ||
        throw(ArgumentError(
            "training groups need collapsed-frontier candidate edits"))
    any(heldout.collapsed_frontier) ||
        throw(ArgumentError(
            "held-out groups need collapsed-frontier candidate edits"))
    any(!, training.collapsed_frontier) ||
        throw(ArgumentError("training groups need finite-gap candidate edits"))
    any(!, heldout.collapsed_frontier) ||
        throw(ArgumentError("held-out groups need finite-gap candidate edits"))
    (all(training.labels) || !any(training.labels)) &&
        throw(ArgumentError("candidate-edit training labels need both classes"))
    collapsed_training_indices = findall(training.collapsed_frontier)
    multifeature_model = fit_logistic_map(training.features, training.labels)
    gap_model = fit_logistic_map(training.gap_features, training.labels)
    out_rows = _heldout_metric_rows(
        multifeature_model, gap_model, heldout; nbins = nbins)
    println("train candidate edits=$(training.n)  held-out=$(heldout.n)  " *
            "train groups=$(length(unique(training.groups)))  " *
            "held-out groups=$(length(unique(heldout.groups)))")
    println("collapsed-frontier training edits=$(length(collapsed_training_indices))")
    for row in out_rows
        rate_label = row.error_rate === nothing ? "pooled" : "err=$(row.error_rate)"
        println("$(rpad(rate_label, 10)) $(rpad(row.model, 23)) " *
                "AUROC=$(round(row.auroc; digits = 3))  " *
                "ECE=$(round(row.ece; digits = 4))  " *
                "Brier=$(round(row.brier; digits = 4))  n=$(row.n)")
    end
    csv_path = joinpath(results_dir, "rhizomorph_gap_calibration.csv")
    open(csv_path, "w") do io
        println(io, "scope,error_rate,model,n,positive_frac,auroc,ece,brier")
        for r in out_rows
            error_rate = r.error_rate === nothing ? "" : string(r.error_rate)
            println(io,
                join(
                    (r.scope, error_rate, r.model, r.n, r.positive_frac,
                        r.auroc, r.ece, r.brier), ","))
        end
    end
    println("Wrote $(csv_path)")

    model_path = joinpath(results_dir, "rhizomorph_gap_calibration_model.csv")
    open(model_path, "w") do io
        println(io, "model,term,coefficient")
        println(io, "multifeature-logistic,intercept,$(multifeature_model.a)")
        for (feature_name, coefficient) in
            zip(CORRECTION_CONFIDENCE_FEATURES, multifeature_model.b)
            println(io,
                "multifeature-logistic,$(feature_name),$(coefficient)")
        end
        println(io, "gap-only-logistic,intercept,$(gap_model.a)")
        println(io, "gap-only-logistic,raw_gap,$(gap_model.b[1])")
        println(io,
            "gap-only-logistic,collapsed_frontier,$(gap_model.b[2])")
    end
    println("Wrote $(model_path)")

    manifest_path = joinpath(results_dir, "rhizomorph_gap_calibration_manifest.csv")
    open(manifest_path, "w") do io
        println(io, "key,value")
        println(io, "seed,$(seed)")
        println(io, "k,$(k)")
        println(io, "genome_length,$(genome_length)")
        println(io, "readlen,$(readlen)")
        println(io, "coverage,$(coverage)")
        println(io, "replicates,$(replicates)")
        println(io, "holdout_fraction,$(holdout_fraction)")
        println(io, "max_contract_skip_fraction,$(max_contract_skip_fraction)")
        println(io, "training_candidate_edits,$(training.n)")
        println(io, "heldout_candidate_edits,$(heldout.n)")
        println(io,
            "training_collapsed_frontier_edits,$(length(collapsed_training_indices))")
        println(io,
            "heldout_collapsed_frontier_edits,$(count(heldout.collapsed_frontier))")
        println(io, "error_rates,$(join(error_rates, ';'))")
        for error_rate in error_rates
            println(io,
                "contract_skips_$(error_rate),$(contract_skips[error_rate])")
            println(io,
                "contract_read_passes_$(error_rate)," *
                "$(contract_read_passes[error_rate])")
        end
    end
    println("Wrote $(manifest_path)")
    return return_artifact ?
           (
        rows = out_rows,
        multifeature_model = multifeature_model,
        gap_model = gap_model,
        feature_names = CORRECTION_CONFIDENCE_FEATURES,
        training = training,
        heldout = heldout,
        dataset = dataset,
        split = split,
        contract_skips = contract_skips,
        contract_read_passes = contract_read_passes,
        max_contract_skip_fraction = max_contract_skip_fraction,
        metrics_path = csv_path,
        model_path = model_path,
        manifest_path = manifest_path,
        serving_config = merge(serving_config, (assigned_q = 20,))
    ) : out_rows
end

# Run when invoked directly (julia benchmarking/viterbi_gap_calibration.jl).
if abspath(PROGRAM_FILE) == @__FILE__
    run_gap_calibration()
end
