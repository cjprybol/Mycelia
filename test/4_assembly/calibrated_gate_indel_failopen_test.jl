# Calibrated gate under INDEL mode (td-4osf Stage 3, review convergent fix):
# the positional revert requires a stable base<->column map, which only a
# length-preserving SUBSTITUTION decode provides. A balanced indel decode (one
# insertion + one deletion) is net-length-neutral, so a length check alone would
# not exclude it — the gate is gated on `indel_params === nothing` instead.
#
# This test proves the exclusion behaviorally: with `indel_params` set, the
# `calibrated_gap_threshold` has NO effect on the corrected output (a large
# threshold that WOULD revert every finite-gap base in substitution mode is a
# no-op here). Contrast: the substitution-mode gate engagement is proven by
# rhizomorph_calibrated_gate_frontier_test.jl.
#
# Run:
#   julia --project=. test/4_assembly/calibrated_gate_indel_failopen_test.jl

import Test
import Mycelia
import FASTX
import Random

function _write_fastq(
        reads::AbstractVector{<:FASTX.FASTQ.Record}, dir::AbstractString)::String
    path = joinpath(dir, "in.fastq")
    open(FASTX.FASTQ.Writer, path) do w
        for r in reads
            write(w, r)
        end
    end
    return path
end

function _assemble_corrected(fastq::AbstractString, k::Int;
        threshold::Union{Nothing, Real},
        indel_params::Union{Nothing, Mycelia.IndelDecodeParams},
        outdir::AbstractString,
        probability_model::Union{Nothing, Mycelia.CorrectionConfidenceModel} = nothing,
        probability_threshold::Real = 0.5,
        correction_feature_sink::Any = nothing)::NamedTuple
    Random.seed!(42)
    result = Mycelia.mycelia_iterative_assemble(fastq;
        max_k = max(k, 13), skip_solid = false, graph_mode = :doublestrand,
        n_k_rungs = 1, max_iterations_per_k = 1, hard_window = false,
        soft_em = false, cheap_correct = true, beam_width = nothing,
        calibrated_gap_threshold = threshold, indel_params = indel_params,
        calibrated_probability_model = probability_model,
        calibrated_probability_threshold = probability_threshold,
        correction_feature_sink = correction_feature_sink,
        verbose = false, enable_checkpointing = false, output_dir = outdir)
    corrected_fastq = get(result[:metadata], :final_fastq_file, nothing)
    (corrected_fastq === nothing || !isfile(corrected_fastq)) &&
        error("corrector produced no :final_fastq_file")
    out = Dict{String, String}()
    open(FASTX.FASTQ.Reader, corrected_fastq) do reader
        for rec in reader
            out[FASTX.identifier(rec)] = FASTX.sequence(String, rec)
        end
    end
    return (sequences = out, metadata = result[:metadata])
end

Test.@testset "calibrated gate is excluded under indel mode (fail-open)" begin
    k = 11
    rng = Random.MersenneTwister(9)
    genome = String(rand(rng, ['A', 'C', 'G', 'T'], 500))
    alts = Dict('A' => "CGT", 'C' => "AGT", 'G' => "ACT", 'T' => "ACG")
    reads = FASTX.FASTQ.Record[]
    for i in 1:60
        st = rand(rng, 1:(500 - 120 + 1))
        clean = genome[st:(st + 119)]
        obs = String([rand(rng) < 0.06 ? alts[c][rand(rng, 1:3)] : c
                      for c in collect(clean)])
        push!(reads, FASTX.FASTQ.Record("r$i", obs, String(fill('I', 120))))
    end

    # An indel-capable decode profile (values mirror a modest nanopore-ish setting).
    indel = Mycelia.IndelDecodeParams(0.06, 0.4, 0.4, 0.1, 0.1, 3, 3, nothing)

    results = mktempdir() do directory
        fastq = _write_fastq(reads, directory)
        base = _assemble_corrected(fastq, k;
            threshold = nothing, indel_params = indel,
            outdir = joinpath(directory, "base"))
        # A huge threshold would revert every finite-gap correction if the gate ran.
        # Attach the sink to this same run so the strengthened provenance checks do
        # not add a fourth expensive indel assembly to the existing test.
        sink_calls = Ref(0)
        raw = _assemble_corrected(fastq, k;
            threshold = 1.0e6, indel_params = indel,
            outdir = joinpath(directory, "raw"),
            correction_feature_sink = _ -> (sink_calls[] += 1))
        probability_model = Mycelia.CorrectionConfidenceModel(
            0.0, 0.0, 0.0, 0.0, 0.0, 0.5)
        probability = _assemble_corrected(fastq, k;
            threshold = nothing, indel_params = indel,
            outdir = joinpath(directory, "probability"),
            probability_model = probability_model, probability_threshold = 1.0)
        (; base, raw, probability, sink_calls = sink_calls[])
    end

    # A real indel decode ran, but all substitution-only controls were disabled.
    Test.@test results.base.metadata[:corrector_errors][:indel_decodes] > 0
    Test.@test !isempty(results.base.sequences)
    Test.@test results.base.sequences == results.raw.sequences
    Test.@test results.base.sequences == results.probability.sequences

    for (result, kind) in ((results.raw, :raw_gap),
            (results.probability, :probability))
        metadata = result.metadata
        Test.@test metadata[:calibrated_gate_kind] == kind
        Test.@test metadata[:calibrated_gate_requested] === true
        Test.@test metadata[:calibrated_gate_effective] === false
        Test.@test metadata[:calibrated_gate_executed] === false
        Test.@test metadata[:calibrated_gate_disabled_reason] ==
                   "indel_params_requested"
        Test.@test metadata[:calibrated_gate_candidate_evaluations] == 0
        Test.@test metadata[:calibrated_gate_reverts] == 0
    end
    Test.@test results.sink_calls == 0
    Test.@test results.raw.metadata[:correction_feature_sink_requested] === true
    Test.@test results.raw.metadata[:correction_feature_sink_effective] === false
    Test.@test results.raw.metadata[:correction_feature_sink_executed] === false
    Test.@test results.raw.metadata[:correction_feature_sink_disabled_reason] ==
               "indel_params_requested"
    Test.@test results.raw.metadata[:calibrated_feature_events] == 0

    # Default-off assembly metadata remains free of opt-in gate/sink keys.
    Test.@test !haskey(results.base.metadata, :calibrated_gate_requested)
    Test.@test !haskey(results.base.metadata, :correction_feature_sink_requested)
end
