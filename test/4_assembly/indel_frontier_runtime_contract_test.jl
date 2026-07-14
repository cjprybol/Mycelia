module IndelFrontierRuntimeContractTest

import Test

include(joinpath(
    @__DIR__,
    "..",
    "..",
    "benchmarking",
    "indel_frontier_runtime.jl",
))

function fake_runtime_path(diagnostics::Dict{Symbol, Any})::NamedTuple
    return (path = :available, diagnostics = diagnostics)
end

Test.@testset "Runtime calibration requires a complete pair-HMM trace" begin
    missing_terminal = _indel_frontier_replicate_row(
        "measurement",
        1,
        1.0,
        fake_runtime_path(Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :completed_columns => 3,
            :move_trace => Symbol[:M, :M, :M],
            :read_index_trace => Int[1, 2, 3],
        )),
        3,
    )
    Test.@test missing_terminal.pair_hmm_path_valid
    Test.@test !missing_terminal.terminal_contract_valid
    Test.@test missing_terminal.trace_complete
    Test.@test !missing_terminal.complete
    Test.@test !missing_terminal.full_decode

    malformed_trace = _indel_frontier_replicate_row(
        "measurement",
        2,
        1.0,
        fake_runtime_path(Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => false,
            :completed_columns => 3,
            :decoded_read_index => 3,
            :move_trace => Symbol[:M, :M, :M],
            :read_index_trace => Int[1, 3, 3],
        )),
        3,
    )
    Test.@test malformed_trace.terminal_contract_valid
    Test.@test !malformed_trace.trace_complete
    Test.@test !malformed_trace.pair_hmm_valid

    complete = _indel_frontier_replicate_row(
        "measurement",
        3,
        1.0,
        fake_runtime_path(Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => false,
            :completed_columns => 3,
            :decoded_read_index => 3,
            :move_trace => Symbol[:M, :D, :M, :I],
            :read_index_trace => Int[1, 1, 2, 3],
        )),
        3,
    )
    Test.@test complete.terminal_contract_valid
    Test.@test complete.trace_complete
    Test.@test complete.pair_hmm_valid
    Test.@test complete.complete
    Test.@test complete.full_decode

    truncated = _indel_frontier_replicate_row(
        "measurement",
        4,
        1.0,
        fake_runtime_path(Dict{Symbol, Any}(
            :algorithm => :viterbi_indel_pair_hmm,
            :truncated => true,
            :completed_columns => 2,
            :decoded_read_index => 2,
            :move_trace => Symbol[:M, :M],
            :read_index_trace => Int[1, 2],
        )),
        3,
    )
    Test.@test truncated.truncated
    Test.@test !truncated.terminal_contract_valid
    Test.@test !truncated.trace_complete
    Test.@test !truncated.full_decode
end

end
