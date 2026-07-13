import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

function _sink_record(id::AbstractString, sequence::AbstractString)::FASTX.FASTQ.Record
    return FASTX.FASTQ.Record(id, sequence, String(fill('I', length(sequence))))
end

function _sink_throws_message(f::Function, fragment::AbstractString)::Bool
    try
        f()
    catch error
        return occursin(fragment, sprint(showerror, error))
    end
    return false
end

mutable struct _FailingSolidKmers <: AbstractSet{Any}
    memberships::Int
    fail_after::Int
    failure::Exception
end

function Base.in(::Any, solid_kmers::_FailingSolidKmers)::Bool
    solid_kmers.memberships += 1
    solid_kmers.memberships > solid_kmers.fail_after &&
        throw(solid_kmers.failure)
    return true
end

mutable struct _CollectingFeatureSink
    observations::Vector{Mycelia.CorrectionFeatureObservation}
end

function (sink::_CollectingFeatureSink)(
        observation::Mycelia.CorrectionFeatureObservation)::Nothing
    push!(sink.observations, observation)
    return nothing
end

Test.@testset "correction feature serving sink" begin
    k = 7
    clean = "ATGCGTACGTACGTTAGCCGATACAGGTCA"
    reads = [_sink_record("clean$(i)", clean) for i in 1:10]
    errored = collect(clean)
    errored[15] = errored[15] == 'A' ? 'C' : 'A'
    push!(reads, _sink_record("err1", String(errored)))
    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k;
        dataset_id = "feature-sink", mode = :canonical, memory_profile = :full)

    observations = Mycelia.CorrectionFeatureObservation[]
    Mycelia.improve_read_set_likelihood(reads, graph, k;
        graph_mode = :canonical,
        correction_feature_sink = observation -> push!(observations, observation))

    Test.@test !isempty(observations)
    for observation in observations
        Test.@test observation.read_id in String.(FASTX.identifier.(reads))
        Test.@test observation.k == k
        Test.@test 1 <= observation.position <= length(clean)
        Test.@test observation.observed_base != observation.corrected_base
        Test.@test isfinite(observation.features.raw_gap) ||
                   observation.features.raw_gap == Inf
        Test.@test isfinite(observation.features.min_kmer_support)
        Test.@test 0.0 <= observation.features.competing_branch_support_ratio <= 1.0
    end

    requested_parallel = Mycelia.CorrectionFeatureObservation[]
    Mycelia.improve_read_set_likelihood(
        reads,
        graph,
        k;
        graph_mode = :canonical,
        enable_parallel = true,
        correction_feature_sink = observation ->
            push!(requested_parallel, observation),
    )
    Test.@test requested_parallel == observations
    # Test the parallel policy independently of the CI worker's actual thread count.
    Test.@test Mycelia._correction_pass_uses_parallel(true, 2, false, false)
    Test.@test !Mycelia._correction_pass_uses_parallel(true, 2, false, true)
    Test.@test !Mycelia._correction_pass_uses_parallel(true, 2, true, false)
    Test.@test !Mycelia._correction_pass_uses_parallel(true, 1, false, false)

    callable_diagnostics = Mycelia.CorrectorDiagnostics()
    callable_sink = _CollectingFeatureSink(Mycelia.CorrectionFeatureObservation[])
    Mycelia.improve_read_set_likelihood(
        reads, graph, k; graph_mode = :canonical, diagnostics = callable_diagnostics,
        correction_feature_sink = callable_sink)
    Test.@test callable_sink.observations == observations
    Test.@test callable_diagnostics.feature_events_emitted[] == length(observations)
    Test.@test callable_diagnostics.candidate_substitutions_evaluated[] >=
               callable_diagnostics.feature_events_emitted[]

    Test.@test _sink_throws_message(
        () -> Mycelia.improve_read_set_likelihood(reads, graph, k;
            graph_mode = :canonical,
            correction_feature_sink = _ -> error("sink sentinel")),
        "sink sentinel")

    multi_error_chars = collect(clean)
    for position in (10, 25)
        multi_error_chars[position] = multi_error_chars[position] == 'A' ? 'C' : 'A'
    end
    multi_error = _sink_record("multi-error", String(multi_error_chars))
    complete_observations = Mycelia.CorrectionFeatureObservation[]
    Mycelia.try_viterbi_path_improvement(
        multi_error,
        graph,
        k;
        graph_mode = :canonical,
        correction_feature_sink = observation ->
            push!(complete_observations, observation),
    )
    Test.@test length(complete_observations) >= 2

    typed_multi_error = FASTX.sequence(BioSequences.LongDNA{4}, multi_error)
    kmer_at_15 = collect(Kmers.UnambiguousDNAMers{k}(typed_multi_error))[15][1]
    hard_vertices = Set([BioSequences.canonical(kmer_at_15)])
    windows = Mycelia._hard_window_ranges(
        multi_error, k, hard_vertices; pad = k, max_window = 21)
    Test.@test only(windows) == 8:28
    window_observations = Mycelia.CorrectionFeatureObservation[]
    Mycelia.improve_read_likelihood_windowed_detail(
        multi_error,
        graph,
        k,
        hard_vertices;
        graph_mode = :canonical,
        pad = k,
        max_window = 21,
        correction_feature_sink = observation ->
            push!(window_observations, observation),
    )
    Test.@test !isempty(window_observations)

    # Re-run the exact decoded sub-window through the lower-level serving entry
    # point and splice its candidate into parent coordinates. Equality with the
    # window wrapper's events proves that it threaded both parent id and `lo - 1`.
    window = only(windows)
    original_sequence = FASTX.sequence(String, multi_error)
    sub_sequence = original_sequence[window]
    sub_read = _sink_record("sub-window", sub_sequence)
    direct_observations = Mycelia.CorrectionFeatureObservation[]
    direct_result = Mycelia.try_viterbi_path_improvement(
        sub_read,
        graph,
        k;
        graph_mode = :canonical,
        correction_read_id = "multi-error",
        correction_position_offset = first(window) - 1,
        correction_feature_sink = observation ->
            push!(direct_observations, observation),
    )
    Test.@test direct_result !== nothing
    candidate_chars = collect(original_sequence)
    candidate_chars[window] = collect(FASTX.sequence(String, first(direct_result)))
    candidate_sequence = String(candidate_chars)
    Test.@test window_observations == direct_observations
    for observation in window_observations
        Test.@test observation.read_id == "multi-error"
        Test.@test observation.position in window
        Test.@test observation.observed_base == original_sequence[observation.position]
        Test.@test observation.corrected_base == candidate_sequence[observation.position]
        Test.@test observation.observed_base != observation.corrected_base
    end

    # A later feature-extraction failure invalidates the whole read contract. No
    # observation buffered for an earlier edit may escape to the caller-owned sink.
    transactional_observations = Mycelia.CorrectionFeatureObservation[]
    transactional_diagnostics = Mycelia.CorrectorDiagnostics()
    Mycelia.try_viterbi_path_improvement(
        multi_error,
        graph,
        k;
        graph_mode = :canonical,
        diagnostics = transactional_diagnostics,
        correction_feature_sink = observation ->
            push!(transactional_observations, observation),
        solid_kmers = _FailingSolidKmers(
            0, k, ArgumentError("solid-kmer membership sentinel")),
    )
    Test.@test isempty(transactional_observations)
    Test.@test transactional_diagnostics.gate_skipped[] == 1

    interrupt_error = nothing
    try
        Mycelia.try_viterbi_path_improvement(
            multi_error,
            graph,
            k;
            graph_mode = :canonical,
            correction_feature_sink = _ -> nothing,
            solid_kmers = _FailingSolidKmers(0, 0, InterruptException()),
        )
    catch error
        interrupt_error = error
    end
    Test.@test interrupt_error isa InterruptException

    unexpected_error = nothing
    try
        Mycelia.try_viterbi_path_improvement(
            multi_error,
            graph,
            k;
            graph_mode = :canonical,
            correction_feature_sink = _ -> nothing,
            solid_kmers = _FailingSolidKmers(
                0, 0, ErrorException("unexpected feature sentinel")),
        )
    catch error
        unexpected_error = error
    end
    Test.@test unexpected_error isa ErrorException
    Test.@test occursin(
        "unexpected feature sentinel", sprint(showerror, unexpected_error))

    overflow_model = Mycelia.CorrectionConfidenceModel(
        0.0, floatmax(Float64), -floatmax(Float64), 0.0,
        floatmax(Float64), floatmax(Float64))
    serving_numeric_diagnostics = Mycelia.CorrectorDiagnostics()
    serving_numeric_error = nothing
    try
        Mycelia.try_viterbi_path_improvement(
            multi_error,
            graph,
            k;
            graph_mode = :canonical,
            diagnostics = serving_numeric_diagnostics,
            calibrated_probability_model = overflow_model,
        )
    catch error
        serving_numeric_error = error
    end
    Test.@test serving_numeric_error isa Mycelia._CorrectionConfidenceServingError
    Test.@test occursin(
        "logit must not be NaN", sprint(showerror, serving_numeric_error))
    Test.@test serving_numeric_diagnostics.structural_errors[] == 0
    Test.@test serving_numeric_diagnostics.gate_skipped[] == 0

    bounds_diagnostics = Mycelia.CorrectorDiagnostics()
    bounds_error = nothing
    try
        Mycelia.try_viterbi_path_improvement(
            multi_error,
            graph,
            k;
            graph_mode = :canonical,
            diagnostics = bounds_diagnostics,
            correction_feature_sink = _ -> nothing,
            solid_kmers = _FailingSolidKmers(
                0, 0, BoundsError([1], 2)),
        )
    catch error
        bounds_error = error
    end
    Test.@test bounds_error isa BoundsError
    Test.@test occursin("BoundsError", sprint(showerror, bounds_error))
    Test.@test bounds_diagnostics.structural_errors[] == 0
    Test.@test bounds_diagnostics.gate_skipped[] == 0

    numeric_model_error = nothing
    try
        Mycelia.correction_confidence_probability(
            Mycelia.CorrectionConfidenceModel(
                0.0, 1.0e308, -1.0e308, 0.0, 0.0, 0.0),
            Mycelia.CorrectionConfidenceFeatures(
                1.0e308, 1.0e308, false, 0.5))
    catch error
        numeric_model_error = error
    end
    Test.@test numeric_model_error isa ArgumentError
    Test.@test !(numeric_model_error isa Mycelia._CorrectionFeatureContractError)
    Test.@test occursin(
        "logit must not be NaN", sprint(showerror, numeric_model_error))

    ambiguous = _sink_record("ambiguous", repeat("N", length(clean)))
    ambiguous_diagnostics = Mycelia.CorrectorDiagnostics()
    ambiguous_result = Mycelia.find_optimal_sequence_path(
        ambiguous,
        graph,
        k;
        graph_mode = :canonical,
        diagnostics = ambiguous_diagnostics,
        correction_feature_sink = _ -> nothing,
    )
    Test.@test first(ambiguous_result) == ambiguous
    Test.@test ambiguous_diagnostics.unkmerizable_reads[] == 1
    Test.@test ambiguous_diagnostics.gate_skipped[] == 1

    extreme = Mycelia.CorrectionConfidenceModel(
        0.0, floatmax(Float64), -floatmax(Float64), 0.0, 0.0, 0.5)
    features = Mycelia.CorrectionConfidenceFeatures(
        floatmax(Float64), floatmax(Float64), false, 0.5)
    Test.@test _sink_throws_message(
        () -> Mycelia.correction_confidence_probability(extreme, features),
        "logit must not be NaN")

    collapsed_model = Mycelia.CorrectionConfidenceModel(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.2)
    collapsed_features = Mycelia.CorrectionConfidenceFeatures(
        Inf, 3.0, true, 0.8)
    Test.@test Mycelia.correction_confidence_probability(
        collapsed_model, collapsed_features) ≈ 1.0 / (1.0 + exp(-0.2))
    finite_zero_features = Mycelia.CorrectionConfidenceFeatures(
        0.0, 3.0, true, 0.8)
    Test.@test Mycelia.correction_confidence_probability(
        collapsed_model, finite_zero_features) == 0.5
end
