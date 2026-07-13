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
        reads::Vector{FASTX.FASTQ.Record}, dir::String)::String
    path = joinpath(dir, "in.fastq")
    open(FASTX.FASTQ.Writer, path) do w
        for r in reads
            write(w, r)
        end
    end
    return path
end

function _assemble_corrected(
        fastq::String,
        k::Int;
        threshold::Union{Float64, Nothing},
        indel_params::Union{Mycelia.IndelDecodeParams, Nothing},
        outdir::String,
)::NamedTuple
    Random.seed!(42)
    result = Mycelia.mycelia_iterative_assemble(fastq;
        max_k = max(k, 13), skip_solid = false, graph_mode = :doublestrand,
        n_k_rungs = 1, max_iterations_per_k = 1, hard_window = false,
        soft_em = false, cheap_correct = true, beam_width = nothing,
        calibrated_gap_threshold = threshold, indel_params = indel_params,
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
    return (
        sequences = out,
        bytes = read(corrected_fastq),
        metadata = result[:metadata],
    )
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

    corrected = mktempdir() do d
        base = _assemble_corrected(_write_fastq(reads, d), k;
            threshold = nothing, indel_params = indel, outdir = mktempdir())
        # A huge threshold would revert EVERY finite-gap correction if the gate ran.
        gated = _assemble_corrected(_write_fastq(reads, d), k;
            threshold = 1.0e6, indel_params = indel, outdir = mktempdir())
        substitution = _assemble_corrected(_write_fastq(reads, d), k;
            threshold = nothing, indel_params = nothing, outdir = mktempdir())
        (base = base, gated = gated, substitution = substitution)
    end

    # This fixture is above the measured indel budget, so raw nanopore intent must
    # stage OFF and take the byte-identical substitution fallback. PR #412's
    # positional calibrated gate remains excluded based on RAW indel intent.
    max_vertices = maximum(
        iteration[:memory_kmers]
        for iterations in values(corrected.base.metadata[:iteration_history])
        for iteration in iterations
    )
    Test.@test max_vertices > Mycelia._INDEL_DECODE_MAX_VERTICES
    staged = corrected.base.metadata[:staged_indel_rungs]
    dense_ks = Set(
        rung_k
        for (rung_k, counts) in corrected.base.metadata[:rung_vertex_counts]
        if maximum(counts) > Mycelia._INDEL_DECODE_MAX_VERTICES
    )
    Test.@test !isempty(dense_ks)
    Test.@test all(
        rung.n_vertices <= Mycelia._INDEL_DECODE_MAX_VERTICES &&
        !(rung.k in dense_ks)
        for rung in staged
    )
    Test.@test corrected.base.metadata[:corrector_errors][:indel_decodes] > 0
    Test.@test !isempty(corrected.base.sequences)
    Test.@test corrected.base.bytes == corrected.gated.bytes
    Test.@test corrected.base.bytes == corrected.substitution.bytes
end
