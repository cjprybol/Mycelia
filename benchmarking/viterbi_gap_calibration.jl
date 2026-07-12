# Tier-2 per-base Viterbi-gap calibration (td-4osf).

import Mycelia

include(joinpath(@__DIR__, "calibration_metrics.jl"))
include(joinpath(@__DIR__, "benchmark_kmer_graph_fixture.jl"))

"""
    collect_gap_truth(fixture, observed_sequence, truth_sequence; config)

Decode one fixed-length observation and align `position_gaps[i]` with
`path.steps[i + 1]`. A k-mer transition contributes the newly emitted base at
read position `i + k`; its label is true iff that corrected base matches truth.
Collapsed-frontier `Inf` sentinels are excluded before calibration.
"""
function collect_gap_truth(fixture::NamedTuple,
        observed_sequence::AbstractString,
        truth_sequence::AbstractString;
        config::Mycelia.ViterbiCorrectionConfig =
            Mycelia.ViterbiCorrectionConfig(record_position_gaps = true))::NamedTuple
    length(observed_sequence) == length(truth_sequence) ||
        throw(ArgumentError("observed and truth sequences must have equal length"))
    observation = fixture.to_observation(observed_sequence)
    result = Mycelia.correct_observations(fixture.graph, [observation]; config = config)
    path = something(result.paths[1].path, throw(ArgumentError("decode produced no path")))
    gaps = result.paths[1].diagnostics[:position_gaps]
    length(gaps) == length(path.steps) - 1 ||
        throw(ArgumentError("position-gap/path alignment contract violated"))
    k = length(String(path.steps[1].vertex_label))
    scores = Float64[]
    labels = Bool[]
    positions = Int[]
    for i in eachindex(gaps)
        isfinite(gaps[i]) || continue
        corrected_kmer = String(path.steps[i + 1].vertex_label)
        position = i + k
        position <= lastindex(truth_sequence) || continue
        push!(scores, gaps[i])
        push!(labels, corrected_kmer[end] == truth_sequence[position])
        push!(positions, position)
    end
    return (; scores, labels, positions, result)
end

"""Fit and evaluate an isotonic probability map on finite gap/truth pairs."""
function calibrate_gap_probability(scores::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; nbins::Int = 10)::NamedTuple
    model = fit_isotonic_map(scores, labels)
    probabilities = predict_isotonic(model, scores)
    return (
        model = model,
        probabilities = probabilities,
        auroc = auroc(scores, labels),
        ece = expected_calibration_error(probabilities, labels; nbins = nbins),
        brier = brier_score(probabilities, labels),
        reliability = reliability_bins(probabilities, labels; nbins = nbins),
    )
end
