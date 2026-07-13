# Substitution-decode reconstruction COMPLETENESS (td-pw7p.1). Regression guard
# for the defect where the ungated full-read substitution corrector emitted
# TRUNCATED / length-divergent reads, which the length-preserving scoring contract
# then silently EXCLUDED — so only a small fraction of reads were scored (e.g.
# 5/200) and no frontier verdict was valid.
#
# The contract this test pins: in substitution mode EVERY corrected record is
# length-preserving OR the read fails open to its original sequence — either way
# it keeps its identifier and is scored at input length. So `reads_scored` must
# equal the read count and every corrected length must equal its input length.
#
# NON-VACUITY: err=0.10 at k=21 drives the ungated decode to dead-end at the first
# uncorrectable position, so WITHOUT the fail-open guard the length invariant below
# is violated for most reads (they come back truncated) and `reads_scored` collapses
# far below the read count. The positive control (`injected > 0`) rules out a
# trivially-passing null with nothing to correct.
#
# Run:
#   julia --project=. test/4_assembly/rhizomorph_substitution_reconstruction_completeness_test.jl

import Test
import Random
import BioSequences
import FASTX

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_correction_pr_curve.jl"))

# Substitution-only, ungated full-read decode — the exact operating point whose
# truncation the fix addresses. Single pass keeps the test CI-fast; the per-read
# Viterbi decode (and thus the truncation the fix guards) runs every pass.
const _COMPLETENESS_POINT = (
    name = "noskip+nogate", skip_solid = false, cheap_correct = true,
    hard_window = false, soft_em = false, n_k_rungs = 1, max_iterations_per_k = 1,
    beam_width = nothing, calibrated_gap_threshold = nothing)

function _run_completeness_case(; err::Float64, seed::Int)
    k = 21
    rng = Random.MersenneTwister(seed)
    genome = BioSequences.randdnaseq(rng, 500)
    records, truth_by_id, observed_by_id,
    injected, _ = simulate_substitution_reads(genome, 150, 10, err, rng)

    Test.@test injected > 0                     # positive control: errors exist to correct
    Test.@test !isempty(records)

    input_len_by_id = Dict(FASTX.identifier(r) => length(FASTX.sequence(String, r))
    for r in records)

    Random.seed!(seed)
    corrected = correct_reads_at_point(records, k, _COMPLETENESS_POINT)

    # (1) IDENTIFIER PRESERVATION: every input read is present in the output.
    for id in keys(input_len_by_id)
        Test.@test haskey(corrected, id)
    end

    # (2) LENGTH-PRESERVING CONTRACT (the non-vacuous core): every corrected
    # substitution record equals its input length — no truncated/length-divergent
    # reads survive to the output. Violated pre-fix for the truncated reads.
    n_length_matched = 0
    for (id, corrected_seq) in corrected
        if haskey(input_len_by_id, id)
            Test.@test length(corrected_seq) == input_len_by_id[id]
            n_length_matched += length(corrected_seq) == input_len_by_id[id] ? 1 : 0
        end
    end
    Test.@test n_length_matched == length(input_len_by_id)

    # (3) COMPLETENESS: with lengths preserved, per_base_metrics scores every read
    # (no length exclusions) — the property the frontier verdict depends on.
    m = per_base_metrics(truth_by_id, observed_by_id, corrected)
    Test.@test m.reads_scored == length(records)
end

Test.@testset "substitution reconstruction completeness (fail-open length contract)" begin
    # Both operating points named in the bead acceptance criteria.
    _run_completeness_case(err = 0.05, seed = 20260713)
    _run_completeness_case(err = 0.10, seed = 20260714)
end
