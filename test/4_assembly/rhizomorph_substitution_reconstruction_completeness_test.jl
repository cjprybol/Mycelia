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
# NON-VACUITY IS SELF-CERTIFIED (not seed-trusted): the test asserts the fail-open
# guard actually fired this run — `substitution_length_divergences > 0` — via the
# corrector's own counter. Without the guard those same reads come back truncated
# (length invariant violated) and `reads_scored` collapses far below the read
# count. So a seed/RNG/Julia-version drift that stopped truncating would fail the
# `> 0` assertion loudly rather than silently turning this into a green no-op on
# buggy code. `injected > 0` is the positive control (errors exist to correct).
#
# WHAT THIS TEST DOES NOT PIN (verified elsewhere / deferred):
#   * WHICH reads failed open, and that a failed-open read equals its ORIGINAL
#     sequence byte-for-byte, plus quality preservation — `correct_reads_at_point`
#     returns sequences only, not records, and exposes counts not per-read flags.
#     The fail-open path returning the original record verbatim (ID + qualities +
#     sequence) is verified by code review of `try_viterbi_path_improvement` ->
#     `find_optimal_sequence_path` (returns the original `read` on `nothing`).
#   * That completeness came from PRESERVED corrections vs TOTAL pass-through — at
#     this operating point the ungated decode is mostly pass-through by design
#     (that IS the truncation defect, tracked in a follow-up), so a `recall > 0`
#     floor would be flaky; we assert completeness + that the guard fired, not a
#     correction rate.
#   * Indel mode (deliberately exempt from the guard) — this harness is
#     substitution-only; indel-mode length-change survival is a separate test.
#   * One geometry (500 bp genome, 150 bp reads, 10x, k=21, single pass, two error
#     rates) — matches the bead's acceptance config; broader geometries not swept.
#
# Run:
#   julia --project=. test/4_assembly/rhizomorph_substitution_reconstruction_completeness_test.jl

import Test
import Random
import BioSequences
import FASTX
import Mycelia

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
    corrected,
    corrector_errors = correct_reads_at_point(
        records, k, _COMPLETENESS_POINT; return_corrector_errors = true)

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

    # (4) SELF-CERTIFIED NON-VACUITY: the fail-open guard actually fired this run.
    # If this is 0, the decode did NOT truncate at this seed/config, so (2)/(3)
    # would pass even on the pre-fix code — a green no-op. Asserting > 0 makes the
    # regression self-certifying against seed/RNG drift, and is the sole coverage
    # of the new counter.
    Test.@test get(corrector_errors, :substitution_length_divergences, 0) > 0

    # (5) TELEMETRY-DRIFT GUARD: the exported corrector-error dict must expose one
    # key per CorrectorDiagnostics field. If a future field is added to the struct
    # but not to the finalize export dict (a real, previously-latent omission —
    # gate_skipped), this count mismatch catches it.
    Test.@test length(corrector_errors) == fieldcount(Mycelia.CorrectorDiagnostics)
end

Test.@testset "substitution reconstruction completeness (fail-open length contract)" begin
    # Both operating points named in the bead acceptance criteria.
    _run_completeness_case(err = 0.05, seed = 20260713)
    _run_completeness_case(err = 0.10, seed = 20260714)
end
