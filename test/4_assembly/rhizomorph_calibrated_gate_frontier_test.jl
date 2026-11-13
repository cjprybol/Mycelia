# Calibrated gate frontier smoke (td-4osf Stage 3): at err=0.10, gating the
# ungated (noskip+nogate) corrector by a positive gap threshold moves along the
# precision-recall frontier in the intended direction — over-correction does not
# rise and recall is retained.
#
# This is the CI-fast counterpart to the full sweep + assess_frontier_dominance in
# benchmarking/rhizomorph_correction_pr_curve.jl (which runs at real scale and is
# not swept by Pkg.test()). To keep the monotone guarantee exact and the runtime
# small, it uses a SINGLE decode pass (n_k_rungs=1, max_iterations_per_k=1): the
# gate is then pure post-decode processing, so it can only REVERT corrections —
# never add an over-correction, never gain a true fix — relative to the ungated
# run on the identical reads + RNG seed.
#
# Run:
#   julia --project=. test/4_assembly/rhizomorph_calibrated_gate_frontier_test.jl

import Test
import Mycelia
import Random
import BioSequences

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_correction_pr_curve.jl"))

Test.@testset "calibrated gate frontier smoke (err=0.10, single pass)" begin
    k = 11
    rng = Random.MersenneTwister(7)
    genome = BioSequences.randdnaseq(rng, 600)
    records, truth_by_id, observed_by_id,
    injected, _ = simulate_substitution_reads(genome, 150, 10, 0.10, rng)
    Test.@test injected > 0                          # there ARE errors to correct

    base_point = (name = "noskip+nogate", skip_solid = false, cheap_correct = true,
        hard_window = false, soft_em = false, n_k_rungs = 1, max_iterations_per_k = 1,
        beam_width = nothing, calibrated_gap_threshold = nothing)
    # Threshold kept at 16 after the strand-label fix (mixed double-strand edges
    # are traversed as Forward, not as whichever evidence entry the set yields
    # first). Measured on this toy: the ungated corrector already has
    # over_rate == 0, and thresholds 0.25-32 leave recall, over-rate, and
    # correction rate unchanged. There is no over-correction left for the gate
    # to remove, so a strict reduction cannot be the check. Real-scale operating
    # points remain a benchmark question (assess_frontier_dominance).
    gated_point = (name = "calibrated-gap-16.0", skip_solid = false, cheap_correct = true,
        hard_window = false, soft_em = false, n_k_rungs = 1, max_iterations_per_k = 1,
        beam_width = nothing, calibrated_gap_threshold = 16.0)

    # Permissive limit: gate engaged (records gaps) but reverts nothing.
    permissive_point = merge(base_point,
        (name = "permissive", calibrated_gap_threshold = -Inf))

    run_point(pt) = begin
        Random.seed!(42)   # identical stochastic conditions per point (fair comparison)
        correct_reads_at_point(records, k, pt)
    end

    base_corr = run_point(base_point)          # gate OFF (threshold === nothing)
    perm_corr = run_point(permissive_point)    # gate ON, reverts nothing
    gated_corr = run_point(gated_point)        # gate ON, positive threshold

    # (1) PIPELINE BYTE-IDENTITY: recording gaps + a no-op gate never perturb the
    # decode, so gate-OFF and the permissive limit produce identical corrected reads.
    Test.@test base_corr == perm_corr

    base = per_base_metrics(truth_by_id, observed_by_id, base_corr)
    gated = per_base_metrics(truth_by_id, observed_by_id, gated_corr)

    Test.@test base.reads_scored > 0
    for m in (base, gated)
        Test.@test isfinite(m.recall) && 0.0 <= m.recall <= 1.0
        Test.@test m.over_rate >= 0.0
    end

    # (2) MONOTONE ENVELOPE (exact under a single post-decode pass): gating reverts
    # corrected bases toward observed, so it can only LOWER over-correction and the
    # true-fix count — never raise them.
    Test.@test gated.over_rate <= base.over_rate + 1e-9
    Test.@test gated.recall <= base.recall + 1e-9

    # (3) On this toy the ungated corrector has nothing for threshold 16 to
    # revert: over-rate is already 0 and the gate matches the ungated metrics.
    # A decoder change that reintroduces over-correction, or a gate that starts
    # reverting true fixes, breaks the equalities. Real-scale dominance is
    # assess_frontier_dominance, not this smoke.
    Test.@test base.over_rate == 0.0
    Test.@test gated.correction_rate == base.correction_rate
    Test.@test gated.over_rate == base.over_rate
    Test.@test gated.recall == base.recall
end
