# td-4mbg: the per-read decode route (`_frontier_decode_route`) through the REAL
# entry point, `improve_read_set_likelihood`, in BOTH decode loops.
#
# indel_substitution_reroute_test.jl pins the route decision by calling the
# helper directly, which cannot see the wiring: reverting either call site to the
# pass-level test, or passing the pass-level indel params to the whole-read
# decode, leaves that file green. These two testsets cover the wiring, each in the
# serial and the parallel loop (`enable_parallel` selects the loop even on one
# thread; `parallel_decode_batches` proves the parallel one ran).
#
#   1. A short read whose windows are ALL rejected must decode exactly like the
#      substitution-only arm, checked against both the real `:illumina`
#      configuration (`indel_params = nothing`) and `windowed_decode = false`.
#      Soft-edge weight mass guards against two empty decodes comparing equal: a
#      whole-read decode stages one unit per k-mer step (length - k), a decode
#      confined to the truncated hard window stages less.
#   2. In a PARTIALLY admitted pass, each fully demoted short read must be decoded
#      whole (`demoted_whole_read_routes`) and must not reach the pair-HMM
#      (`attempted == admitted_windows`). The counter assertion fails if the
#      route stops demoting such reads; the attempt assertion catches the pass
#      params leaking into the whole-read decode. A call site that ignores the
#      route's `windowed` decision is caught by testset 1, not here.
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
        enable_parallel::Bool, soft_weights = nothing,
        indel_params::Union{Nothing, Mycelia.IndelDecodeParams} =
            _td4mbg_entry_params())
    diagnostics = Mycelia.CorrectorDiagnostics()
    telemetry = Dict{Symbol, Any}()
    result = Logging.with_logger(Logging.NullLogger()) do
        Mycelia.improve_read_set_likelihood(
            reads, graph, k;
            graph_mode = :singlestrand,
            beam_width = 16,
            hard_vertices = hard,
            windowed_decode = windowed_decode,
            indel_params = indel_params,
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
        # The real :illumina configuration: hard-window gate on, no indel params.
        illumina_weights = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        illumina = _td4mbg_entry_run([read], graph, k, hard;
            windowed_decode = true, enable_parallel = enable_parallel,
            soft_weights = illumina_weights, indel_params = nothing)

        # Preconditions: the pass really is a forced frontier pass with every
        # window rejected, and the intended loop ran.
        Test.@test frontier.telemetry[:bounded_windowing_forced]
        Test.@test frontier.telemetry[:admitted_windows] == 0
        Test.@test frontier.telemetry[:rejected_windows] > 0
        Test.@test (frontier.diagnostics.parallel_decode_batches[] > 0) ==
                   enable_parallel

        # The route: identical to the substitution-only decode of the same read,
        # in both the :illumina configuration and with windowing off.
        Test.@test frontier.diagnostics.demoted_whole_read_routes[] == 1
        Test.@test illumina.diagnostics.demoted_whole_read_routes[] == 0
        Test.@test frontier_weights.weights == illumina_weights.weights
        Test.@test frontier.sequences == illumina.sequences
        Test.@test frontier_weights.weights == substitution_weights.weights
        Test.@test frontier.sequences == substitution.sequences
        # Coverage guard, so the equalities above cannot pass as empty decodes:
        # a whole-read decode stages one unit of mass per k-mer step.
        Test.@test isapprox(sum(values(frontier_weights.weights); init = 0.0),
            length(sequence) - k; atol = 1e-9)
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
    fully_demoted = count(
        sources -> !isempty(sources) &&
                   all(==(:substitution), values(sources)),
        values(decision.window_sources))
    Test.@test fully_demoted > 0

    for enable_parallel in (false, true)
        run = _td4mbg_entry_run(reads, graph, k, hard;
            windowed_decode = true, enable_parallel = enable_parallel)
        Test.@test (run.diagnostics.parallel_decode_batches[] > 0) ==
                   enable_parallel
        Test.@test run.telemetry[:admitted_windows] == decision.admitted_windows
        # Each fully demoted short read is decoded whole. Counting them per read
        # also shows none was skipped, which would make the next check vacuous.
        Test.@test run.diagnostics.demoted_whole_read_routes[] == fully_demoted
        # Only admitted windows reach the pair-HMM. A demoted short read decoded
        # whole with the pass params would add an attempt.
        Test.@test run.telemetry[:attempted] == run.telemetry[:admitted_windows]
    end
end
