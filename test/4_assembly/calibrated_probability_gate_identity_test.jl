# Multi-feature calibrated decode gate (td-21eg): typed-model validation,
# byte-identity when permissive, and revert-only behavior.

import FASTX
import Mycelia
import Random
import Test

_prob_rec(id::AbstractString,
    sequence::AbstractString)::FASTX.FASTQ.Record = FASTX.FASTQ.Record(
    id, sequence, String(fill('I', length(sequence))))
_prob_sequences(reads::AbstractVector{<:FASTX.FASTQ.Record})::Vector{String} = [FASTX.sequence(
                                                                                    String,
                                                                                    read)
                                                                                for read in
                                                                                    reads]
_prob_record_bytes(reads::AbstractVector{<:FASTX.FASTQ.Record})::Vector = [(
                                                                               identifier = String(FASTX.identifier(read)),
                                                                               sequence = FASTX.sequence(
                                                                                   String, read),
                                                                               quality = String(FASTX.quality(read))
                                                                           )
                                                                           for read in
                                                                               reads]
_prob_mismatches(a::AbstractString, b::AbstractString)::Int = count(
    pair -> pair[1] !=
            pair[2], zip(a, b))
function _prob_throws_message(f::Function, fragment::AbstractString)::Bool
    try
        f()
    catch error
        return error isa ArgumentError && occursin(fragment, error.msg)
    end
    return false
end

Test.@testset "calibrated probability gate" begin
    model = Mycelia.CorrectionConfidenceModel(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    Test.@test Mycelia.correction_confidence_probability(
        model, 1.0, 4.0, true, 0.75) == 0.5

    Test.@test _prob_throws_message(
        () -> Mycelia.CorrectionConfidenceModel(
            0.0, NaN, 0.0, 0.0, 0.0, 0.5),
        "coefficients must be finite")
    Test.@test _prob_throws_message(
        () -> Mycelia.CorrectionConfidenceModel(
            0.0, 0.0, 0.0, 0.0, 0.0, NaN),
        "coefficients must be finite")

    k = 7
    clean = "ATGCGTACGTACGTTAGCCGATACAGGTCA"
    reads = [_prob_rec("clean$(i)", clean) for i in 1:10]
    errored_chars = collect(clean)
    errored_chars[15] = errored_chars[15] == 'A' ? 'C' : 'A'
    push!(reads, _prob_rec("err1", String(errored_chars)))
    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k;
        dataset_id = "prob-gate", mode = :canonical, memory_profile = :full)

    function run_gate(model_arg::Union{Nothing, Mycelia.CorrectionConfidenceModel},
            threshold::Float64)::NamedTuple
        Random.seed!(11)
        diagnostics = Mycelia.CorrectorDiagnostics()
        corrected = FASTX.FASTQ.Record[]
        for read in reads
            result = Mycelia.try_viterbi_path_improvement(
                read, graph, k;
                graph_mode = :canonical,
                calibrated_probability_model = model_arg,
                calibrated_probability_threshold = threshold,
                diagnostics = diagnostics
            )
            push!(corrected, result === nothing ? read : first(result))
        end
        return (; corrected, diagnostics)
    end

    off = run_gate(nothing, 0.5)
    permissive = run_gate(model, 0.0)
    strict = run_gate(model, 1.0)

    # Engaging feature extraction and gap recording without rejecting any edit is
    # byte-identical to the default-off path, including identifiers and qualities.
    Test.@test _prob_record_bytes(off.corrected) ==
               _prob_record_bytes(permissive.corrected)
    Test.@test permissive.diagnostics.gate_skipped[] == 0
    Test.@test strict.diagnostics.gate_skipped[] == 0

    # A probability gate can only suppress proposed substitutions; it can never
    # fabricate a new edit relative to the observed input.
    observed = _prob_sequences(reads)
    for i in eachindex(reads)
        Test.@test _prob_mismatches(_prob_sequences(strict.corrected)[i], observed[i]) <=
                   _prob_mismatches(_prob_sequences(off.corrected)[i], observed[i])
    end
    off_edits = sum(_prob_mismatches(sequence, observed[i])
    for (i, sequence) in enumerate(_prob_sequences(off.corrected)))
    strict_edits = sum(_prob_mismatches(sequence, observed[i])
    for (i, sequence) in enumerate(_prob_sequences(strict.corrected)))
    Test.@test off_edits > 0
    Test.@test strict_edits < off_edits

    Test.@test _prob_throws_message(
        () -> Mycelia.improve_read_set_likelihood(
            reads, graph, k; calibrated_probability_model = model,
            calibrated_probability_threshold = 1.1),
        "must be in [0, 1]")
    Test.@test _prob_throws_message(
        () -> Mycelia.improve_read_set_likelihood(
            reads, graph, k; calibrated_gap_threshold = 1.0,
            calibrated_probability_model = model),
        "mutually exclusive")

    # Public record entry points validate configuration before data-dependent
    # early returns (short reads and un-k-merizable ambiguous reads).
    short_read = _prob_rec("short", "ATG")
    Test.@test _prob_throws_message(
        () -> Mycelia.improve_read_likelihood(
            short_read, graph, k;
            calibrated_gap_threshold = 1.0,
            calibrated_probability_model = model),
        "mutually exclusive")
    ambiguous_read = _prob_rec("ambiguous", repeat("N", length(clean)))
    Test.@test _prob_throws_message(
        () -> Mycelia.find_optimal_sequence_path(
            ambiguous_read, graph, k;
            calibrated_gap_threshold = 1.0,
            calibrated_probability_model = model),
        "mutually exclusive")
end
