import FASTX
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

    Test.@test_throws InterruptException Mycelia.try_viterbi_path_improvement(
        multi_error,
        graph,
        k;
        graph_mode = :canonical,
        correction_feature_sink = _ -> nothing,
        solid_kmers = _FailingSolidKmers(0, 0, InterruptException()),
    )

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
end
