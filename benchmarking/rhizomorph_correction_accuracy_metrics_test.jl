# Pure unit test for `per_base_metrics` — NO Mycelia dependency, runs in ms.
#
# Includes ONLY the pure metrics file, so the classification math is pinned
# without loading Mycelia (the simulator/corrector integration is tested in
# rhizomorph_correction_accuracy_benchmark_test.jl, which needs Mycelia).
#
# Run:
#   julia benchmarking/rhizomorph_correction_accuracy_metrics_test.jl

import Test

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_metrics.jl"))

Test.@testset "per_base_metrics — classification math" begin
    # --- Case A: a clean fix -------------------------------------------------
    m = per_base_metrics(Dict("r1" => "ACGT"), Dict("r1" => "AGGT"), Dict("r1" => "ACGT"))
    Test.@test m.injected == 1
    Test.@test m.total_edits == 1
    Test.@test m.tp == 1
    Test.@test m.mis_fixes == 0
    Test.@test m.over == 0
    Test.@test m.recall == 1.0
    Test.@test m.precision == 1.0
    Test.@test m.over_rate == 0.0
    Test.@test m.correction_rate == 1 / 4
    Test.@test m.reads_scored == 1
    Test.@test m.scored_fraction == 1.0
    Test.@test m.corrected_unjoined == 0
    # partition invariant
    Test.@test m.total_edits == m.tp + m.mis_fixes + m.over

    # --- Case B: an over-correction (correct base made wrong) ----------------
    m = per_base_metrics(Dict("r1" => "ACGT"), Dict("r1" => "ACGT"), Dict("r1" => "ACGA"))
    Test.@test m.injected == 0
    Test.@test m.total_edits == 1
    Test.@test m.tp == 0
    Test.@test m.mis_fixes == 0
    Test.@test m.over == 1
    Test.@test isnan(m.recall)          # 0 injected -> undefined recall
    Test.@test m.precision == 0.0
    Test.@test m.over_rate == 1 / 4
    Test.@test m.total_edits == m.tp + m.mis_fixes + m.over

    # --- Case C: a missed error (no edit at the error) -----------------------
    m = per_base_metrics(Dict("r1" => "ACGT"), Dict("r1" => "AGGT"), Dict("r1" => "AGGT"))
    Test.@test m.injected == 1
    Test.@test m.total_edits == 0
    Test.@test m.tp == 0
    Test.@test m.mis_fixes == 0
    Test.@test m.recall == 0.0
    Test.@test isnan(m.precision)       # 0 edits -> undefined precision
    Test.@test m.over_rate == 0.0

    # --- Case MIS: an error edited to a WRONG base (mis-fix) -----------------
    # pos-2 error C->G, corrector changes it to T (still wrong). Distinct from
    # both true_fix and over_correction; depresses recall AND precision.
    m = per_base_metrics(Dict("r1" => "ACGT"), Dict("r1" => "AGGT"), Dict("r1" => "ATGT"))
    Test.@test m.injected == 1
    Test.@test m.total_edits == 1
    Test.@test m.tp == 0
    Test.@test m.mis_fixes == 1
    Test.@test m.over == 0
    Test.@test m.recall == 0.0
    Test.@test m.precision == 0.0       # the one edit did not land on truth
    Test.@test m.over_rate == 0.0
    Test.@test m.total_edits == m.tp + m.mis_fixes + m.over

    # --- Case 3CLASS: all three edit classes in one read ---------------------
    # T=ACGTAC O=AGCTAC (p2,p3 errors) C=ATGTGC:
    #   p2 error edited to wrong base (mis), p3 error fixed (tp), p5 correct->G (over).
    m = per_base_metrics(Dict("r1" => "ACGTAC"), Dict("r1" => "AGCTAC"), Dict("r1" => "ATGTGC"))
    Test.@test m.injected == 2
    Test.@test m.tp == 1
    Test.@test m.mis_fixes == 1
    Test.@test m.over == 1
    Test.@test m.total_edits == 3
    Test.@test m.total_edits == m.tp + m.mis_fixes + m.over
    Test.@test m.recall == 1 / 2
    Test.@test m.precision == 1 / 3
    Test.@test m.correct_positions == 4
    Test.@test m.over_rate == 1 / 4

    # --- Case CLEAN: zero errors, zero edits (do-no-harm) --------------------
    m = per_base_metrics(Dict("r1" => "ACGT"), Dict("r1" => "ACGT"), Dict("r1" => "ACGT"))
    Test.@test m.injected == 0
    Test.@test m.total_edits == 0
    Test.@test m.over == 0
    Test.@test m.over_rate == 0.0
    Test.@test m.correction_rate == 0.0
    Test.@test isnan(m.recall)
    Test.@test isnan(m.precision)
    Test.@test m.reads_scored == 1

    # --- Case ALLERR: every position is an error (over_rate NaN branch) -------
    m = per_base_metrics(Dict("r1" => "AAAA"), Dict("r1" => "CCCC"), Dict("r1" => "AAAA"))
    Test.@test m.injected == 4
    Test.@test m.correct_positions == 0
    Test.@test m.tp == 4
    Test.@test m.recall == 1.0
    Test.@test m.precision == 1.0
    Test.@test isnan(m.over_rate)       # 0 correct positions -> undefined over-rate
    Test.@test m.correction_rate == 1.0

    # --- Case D: a dropped read (absent from corrected set) -------------------
    m = per_base_metrics(
        Dict("r1" => "ACGT", "r2" => "TTTT"),
        Dict("r1" => "AGGT", "r2" => "TTTT"),
        Dict("r1" => "ACGT"))
    Test.@test m.reads_dropped == 1
    Test.@test m.reads_scored == 1
    Test.@test m.tp == 1
    Test.@test m.scored_fraction == 1 / 2

    # --- Case E: length-mismatched corrected read (Viterbi indel path) --------
    m = per_base_metrics(Dict("r1" => "ACGT"), Dict("r1" => "AGGT"), Dict("r1" => "ACGTA"))
    Test.@test m.reads_excluded_len == 1
    Test.@test m.reads_scored == 0
    Test.@test m.total_bases == 0
    Test.@test isnan(m.correction_rate)  # 0 total_bases -> undefined
    Test.@test m.scored_fraction == 0.0

    # --- Case UNJOINED: corrected read whose id is not in truth ---------------
    # Signature of an ID-format mismatch; counted, not silently ignored.
    m = per_base_metrics(
        Dict("r1" => "ACGT"),
        Dict("r1" => "AGGT"),
        Dict("r1" => "ACGT", "r2" => "GGGG"))
    Test.@test m.corrected_unjoined == 1
    Test.@test m.reads_scored == 1

    # --- Case EMPTY: empty inputs are safe (all-NaN/zero) ---------------------
    m = per_base_metrics(Dict{String, String}(), Dict{String, String}(), Dict{
        String, String}())
    Test.@test m.n_reads == 0
    Test.@test m.reads_scored == 0
    Test.@test m.total_edits == 0
    Test.@test isnan(m.scored_fraction)
    Test.@test isnan(m.recall)
    Test.@test isnan(m.precision)
    Test.@test isnan(m.over_rate)
    Test.@test isnan(m.correction_rate)

    # --- Aggregation across reads --------------------------------------------
    # r1 clean fix (1 injected, 1 fixed); r2 over-correction (0 injected, 1 spurious).
    m = per_base_metrics(
        Dict("r1" => "ACGT", "r2" => "GGGG"),
        Dict("r1" => "AGGT", "r2" => "GGGG"),
        Dict("r1" => "ACGT", "r2" => "GGCG"))
    Test.@test m.injected == 1
    Test.@test m.total_edits == 2
    Test.@test m.tp == 1
    Test.@test m.over == 1
    Test.@test m.recall == 1.0
    Test.@test m.precision == 1 / 2
    Test.@test m.total_edits == m.tp + m.mis_fixes + m.over
end
