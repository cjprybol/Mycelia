module IndelFrontierRuntimeContractTest

import Test

include(joinpath(
    @__DIR__,
    "..",
    "..",
    "benchmarking",
    "indel_frontier_runtime.jl",
))

function fake_runtime_result(
        diagnostics::Dict{Symbol, Any};
        graph_steps::Int = count(
            !=(:I), get(diagnostics, :move_trace, Symbol[])),
        corrected_labels::Union{Nothing, Vector{String}} = nothing,
)::Tuple{Mycelia.Rhizomorph.ViterbiDecodingResult, Vector{String}}
    labels = ["vertex_$(index)" for index in 1:graph_steps]
    steps = Mycelia.Rhizomorph.WalkStep{String}[
        Mycelia.Rhizomorph.WalkStep(
            label,
            Mycelia.Rhizomorph.Forward,
            1.0,
            1.0,
        )
        for label in labels
    ]
    path = Mycelia.Rhizomorph.GraphPath(steps)
    move_trace = get(diagnostics, :move_trace, Symbol[])
    get!(diagnostics, :path_length, graph_steps)
    get!(
        diagnostics,
        :move_counts,
        Dict{Symbol, Int}(
            move => count(==(move), move_trace)
            for move in (:M, :I, :D)
        ),
    )
    result = Mycelia.Rhizomorph.ViterbiDecodingResult(
        path, 0.0, diagnostics)
    corrected = corrected_labels === nothing ?
                copy(labels) : corrected_labels
    return result, corrected
end

function fake_runtime_row(
        replicate::Int,
        diagnostics::Dict{Symbol, Any},
        expected_columns::Int;
        graph_steps::Int = count(
            !=(:I), get(diagnostics, :move_trace, Symbol[])),
        corrected_labels::Union{Nothing, Vector{String}} = nothing,
)::NamedTuple
    path, corrected = fake_runtime_result(
        diagnostics;
        graph_steps = graph_steps,
        corrected_labels = corrected_labels,
    )
    return _indel_frontier_replicate_row(
        "measurement",
        replicate,
        1.0,
        path,
        corrected,
        expected_columns,
    )
end

Test.@testset "Runtime calibration requires a complete pair-HMM trace" begin
    missing_terminal = fake_runtime_row(
        1,
        Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :completed_columns => 3,
            :move_trace => Symbol[:M, :M, :M],
            :read_index_trace => Int[1, 2, 3],
        ),
        3,
    )
    Test.@test missing_terminal.pair_hmm_path_valid
    Test.@test !missing_terminal.terminal_contract_valid
    Test.@test missing_terminal.trace_complete
    Test.@test missing_terminal.path_trace_aligned
    Test.@test !missing_terminal.complete
    Test.@test !missing_terminal.full_decode

    malformed_trace = fake_runtime_row(
        2,
        Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => false,
            :completed_columns => 3,
            :decoded_read_index => 3,
            :move_trace => Symbol[:M, :M, :M],
            :read_index_trace => Int[1, 3, 3],
        ),
        3,
    )
    Test.@test malformed_trace.terminal_contract_valid
    Test.@test !malformed_trace.trace_complete
    Test.@test malformed_trace.path_trace_aligned
    Test.@test !malformed_trace.pair_hmm_valid

    complete = fake_runtime_row(
        3,
        Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => false,
            :completed_columns => 3,
            :decoded_read_index => 3,
            :move_trace => Symbol[:M, :D, :M, :I],
            :read_index_trace => Int[1, 1, 2, 3],
        ),
        3,
    )
    Test.@test complete.terminal_contract_valid
    Test.@test complete.trace_complete
    Test.@test complete.graph_step_trace_aligned
    Test.@test complete.path_length_aligned
    Test.@test complete.corrected_path_aligned
    Test.@test complete.move_counts_aligned
    Test.@test complete.path_trace_aligned
    Test.@test complete.pair_hmm_valid
    Test.@test complete.complete
    Test.@test complete.full_decode

    truncated = fake_runtime_row(
        4,
        Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => true,
            :completed_columns => 2,
            :decoded_read_index => 2,
            :move_trace => Symbol[:M, :M],
            :read_index_trace => Int[1, 2],
        ),
        3,
    )
    Test.@test truncated.truncated
    Test.@test !truncated.terminal_contract_valid
    Test.@test !truncated.trace_complete
    Test.@test truncated.path_trace_aligned
    Test.@test !truncated.full_decode

    path_mismatch = fake_runtime_row(
        5,
        Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => false,
            :completed_columns => 3,
            :decoded_read_index => 3,
            :move_trace => Symbol[:M, :D, :M, :I],
            :read_index_trace => Int[1, 1, 2, 3],
        ),
        3;
        graph_steps = 2,
    )
    Test.@test path_mismatch.terminal_contract_valid
    Test.@test path_mismatch.trace_complete
    Test.@test !path_mismatch.graph_step_trace_aligned
    Test.@test !path_mismatch.path_trace_aligned
    Test.@test !path_mismatch.full_decode

    corrected_mismatch = fake_runtime_row(
        6,
        Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => false,
            :completed_columns => 3,
            :decoded_read_index => 3,
            :move_trace => Symbol[:M, :D, :M, :I],
            :read_index_trace => Int[1, 1, 2, 3],
        ),
        3;
        corrected_labels = ["wrong_1", "wrong_2", "wrong_3"],
    )
    Test.@test corrected_mismatch.terminal_contract_valid
    Test.@test corrected_mismatch.trace_complete
    Test.@test corrected_mismatch.graph_step_trace_aligned
    Test.@test !corrected_mismatch.corrected_path_aligned
    Test.@test !corrected_mismatch.path_trace_aligned
    Test.@test !corrected_mismatch.full_decode
end

end
