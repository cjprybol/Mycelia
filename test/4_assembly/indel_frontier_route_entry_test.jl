# td-4mbg: the per-read decode route (`_frontier_decode_route`) through the REAL
# entry point, `improve_read_set_likelihood`, in BOTH decode loops.
#
# indel_substitution_reroute_test.jl pins the route decision by calling the
# helper directly. That cannot see the wiring: in round-1 review of the routing
# change, reverting either call site to the old pass-level test, or passing the
# pass-level indel params to the whole-read decode, left every existing test
# green. These two testsets close that gap, each in the serial and the parallel
# loop (`enable_parallel` selects the loop even on one thread;
# `parallel_decode_batches` proves the parallel one ran).
#
#   1. A short read whose windows are ALL rejected must decode exactly like the
#      substitution-only route (`windowed_decode = false`). The observable is the
#      soft-edge weight mass: a whole-read decode stages one unit per k-mer step
#      (length - k), a decode confined to the truncated hard window stages less.
#      Against the old routing the frontier run stages 497.0 against 997.0.
#   2. In a PARTIALLY admitted pass, a fully demoted short read must not reach
#      the pair-HMM. `attempted` must equal `admitted_windows`; passing the pass
#      params to the whole-read decode adds one attempt (3 vs 2).
#
# Testset 2 depends on a fixture whose per-window frontier work sits near the
# default 300_000 budget, so it guards its own preconditions (admitted, some
# windows rejected, one short read fully demoted) and fails loudly if a cost
# model change moves them, instead of passing vacuously.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/indel_frontier_route_entry_test.jl")'

import Test
import Mycelia
import FASTX
import MetaGraphsNext
import Logging
import Random

function _td4mbg_entry_fastq(identifier::String, sequence::String)::FASTX.FASTQ.Record
    return FASTX.FASTQ.Record(identifier, sequence, repeat("I", length(sequence)))
end

function _td4mbg_entry_params()::Mycelia.IndelDecodeParams
    profile = Mycelia.indel_error_profile(:nanopore)
    return Mycelia.IndelDecodeParams(
        profile.base_error_rate, profile.insertion_fraction,
        profile.deletion_fraction, profile.insertion_extend_probability,
        profile.deletion_extend_probability, 3, 3, 16)
end

"""
Order-`order` de Bruijn sequence over `alphabet` (Lyndon-word construction).
"""
function _td4mbg_entry_de_bruijn(order::Int, alphabet::Vector{Char})::String
    alphabet_size = length(alphabet)
    work = zeros(Int, alphabet_size * order + 1)
    indices = Int[]
    function visit(t::Int, period::Int)::Nothing
        if t > order
            order % period == 0 && append!(indices, @view work[2:(period + 1)])
            return nothing
        end
        work[t + 1] = work[t - period + 1]
        visit(t + 1, period)
        for index in (work[t - period + 1] + 1):(alphabet_size - 1)
            work[t + 1] = index
            visit(t + 1, t)
        end
        return nothing
    end
    visit(1, 1)
    return String([alphabet[index + 1] for index in indices])
end

function _td4mbg_entry_run(reads, graph, k, hard; windowed_decode::Bool,
        enable_parallel::Bool, soft_weights = nothing)
    diagnostics = Mycelia.CorrectorDiagnostics()
    telemetry = Dict{Symbol, Any}()
    result = Logging.with_logger(Logging.NullLogger()) do
        Mycelia.improve_read_set_likelihood(
            reads, graph, k;
            graph_mode = :singlestrand,
            beam_width = 16,
            hard_vertices = hard,
            windowed_decode = windowed_decode,
            indel_params = _td4mbg_entry_params(),
            indel_schedule = :frontier_budgeted,
            rung_telemetry = telemetry,
            soft_weights = soft_weights,
            diagnostics = diagnostics,
            enable_parallel = enable_parallel,
        )
    end
    return (sequences = [FASTX.sequence(String, record) for record in result[1]],
        telemetry = telemetry, diagnostics = diagnostics)
end

Test.@testset "td-4mbg route entry: fully rejected short read decodes whole" begin
    k = 3
    alphabet = collect("ACGT")
    # Every 4-mer: a fully branching k=3 graph in which no indel window is
    # affordable, so the frontier schedule rejects them all.
    graph_reads = FASTX.FASTQ.Record[
        _td4mbg_entry_fastq("branch_$(a)$(b)$(c)$(d)", string(a, b, c, d))
        for a in alphabet for b in alphabet for c in alphabet for d in alphabet]
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        graph_reads, k; mode = :singlestrand)
    hard = Set(MetaGraphsNext.labels(graph))
    cycle = _td4mbg_entry_de_bruijn(k, alphabet)
    sequence = collect(first(repeat(cycle, cld(1_000, length(cycle))), 1_000))
    for position in (100, 350, 620, 880)
        sequence[position] = sequence[position] == 'A' ? 'C' : 'A'
    end
    read = _td4mbg_entry_fastq("short_all_rejected", String(sequence))
    Test.@test !Mycelia._windowed_decode_read_is_long(read, k)

    for enable_parallel in (false, true)
        frontier_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        frontier = _td4mbg_entry_run([read], graph, k, hard;
            windowed_decode = true, enable_parallel = enable_parallel,
            soft_weights = frontier_weights)
        substitution_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        substitution = _td4mbg_entry_run([read], graph, k, hard;
            windowed_decode = false, enable_parallel = enable_parallel,
            soft_weights = substitution_weights)

        # Preconditions: the pass really is a forced frontier pass with every
        # window rejected, and the intended loop ran.
        Test.@test frontier.telemetry[:bounded_windowing_forced]
        Test.@test frontier.telemetry[:admitted_windows] == 0
        Test.@test frontier.telemetry[:rejected_windows] > 0
        Test.@test (frontier.diagnostics.parallel_decode_batches[] > 0) ==
                   enable_parallel

        # The route: identical to the substitution-only decode of the same read.
        Test.@test frontier_weights.weights == substitution_weights.weights
        Test.@test frontier.sequences == substitution.sequences
        # Coverage guard, so the equality above cannot pass as two empty decodes:
        # a whole-read decode stages one unit of mass per k-mer step.
        Test.@test sum(values(frontier_weights.weights); init = 0.0) ==
                   length(sequence) - k
        Test.@test frontier.diagnostics.indel_attempts[] == 0
    end
end

"""
Binary (A/C) de Bruijn sequence of order `order`, prefer-one greedy rule.
"""
function _td4mbg_entry_binary_de_bruijn(order::Int)::String
    seen = Set{String}()
    sequence = repeat("A", order)
    push!(seen, sequence)
    while true
        tail = sequence[(end - order + 2):end]
        if !(tail * "C" in seen)
            sequence *= "C"
            push!(seen, tail * "C")
        elseif !(tail * "A" in seen)
            sequence *= "A"
            push!(seen, tail * "A")
        else
            break
        end
    end
    return sequence
end

Test.@testset "td-4mbg route entry: demoted short read gets no pair-HMM in a partially admitted pass" begin
    k = 9
    branch_sequence = _td4mbg_entry_binary_de_bruijn(k + 1)
    rng = Random.MersenneTwister(7)
    reads = FASTX.FASTQ.Record[
        _td4mbg_entry_fastq("linear_short",
            join(rand(rng, ['G', 'T', 'A', 'C'], 1_000))),
        _td4mbg_entry_fastq("branch_short",
            first(repeat(branch_sequence, 2), 1_000)),
    ]
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        vcat(reads, [_td4mbg_entry_fastq("branch_full", branch_sequence)]), k;
        mode = :singlestrand)
    hard = Set(MetaGraphsNext.labels(graph))
    Test.@test !any(read -> Mycelia._windowed_decode_read_is_long(read, k), reads)

    # Fixture preconditions, measured on the scheduler itself: the pass is
    # admitted, some windows are rejected, and at least one short read has every
    # window demoted. If a cost-model change moves these, this fails here rather
    # than letting the assertion below pass vacuously.
    decision = Logging.with_logger(Logging.NullLogger()) do
        Mycelia._evaluate_indel_frontier_schedule(
            reads, graph, k, hard, fill(false, length(reads)),
            _td4mbg_entry_params(), :singlestrand)
    end
    Test.@test decision.admitted
    Test.@test decision.admitted_windows > 0
    Test.@test decision.rejected_windows > 0
    Test.@test any(
        sources -> !isempty(sources) &&
                   all(==(:substitution), values(sources)),
        values(decision.window_sources))

    for enable_parallel in (false, true)
        run = _td4mbg_entry_run(reads, graph, k, hard;
            windowed_decode = true, enable_parallel = enable_parallel)
        Test.@test (run.diagnostics.parallel_decode_batches[] > 0) ==
                   enable_parallel
        Test.@test run.telemetry[:admitted_windows] == decision.admitted_windows
        # Only admitted windows reach the pair-HMM. A demoted short read decoded
        # whole with the pass params would add an attempt.
        Test.@test run.telemetry[:attempted] == run.telemetry[:admitted_windows]
    end
end
