# Stage 2 gap calibration signal check (td-4osf): the per-position Viterbi gap
# must DISCRIMINATE a correct decoded base from a wrong one (AUROC > 0.5), and the
# gap<->path alignment contract must hold on real decodes.
#
# The finite-gap signal only exists on a BRANCHING graph — a read-coverage k-mer
# graph carries the error/repeat-induced competition the decoder navigates,
# whereas a truth-only (linear) graph has a single surviving state per depth and
# every gap is Inf. So this drives the coverage-graph harness in
# benchmarking/viterbi_gap_calibration.jl (not swept by Pkg.test()) at small scale.
#
# Run:
#   julia --project=. test/4_assembly/viterbi_gap_calibration_test.jl

import Test
import Mycelia
import FASTX

include(joinpath(@__DIR__, "..", "..", "benchmarking", "viterbi_gap_calibration.jl"))

Test.@testset "Tier-2 gap calibration signal (AUROC > 0.5)" begin
    # Small, deterministic coverage-graph run at err=0.08/0.10; aggregate finite
    # gaps + truth labels and fit all three calibration models.
    rows = mktempdir() do dir
        run_gap_calibration(; k = 11, genome_length = 500, readlen = 120, coverage = 25,
            error_rates = [0.08, 0.10], seed = 1, results_dir = dir)
    end

    Test.@test !isempty(rows)                       # class balance sufficient to fit
    # BOTH requested error rates must yield results — otherwise a silently-skipped
    # rate (insufficient class balance / regression) would pass unnoticed.
    Test.@test Set(r.error_rate for r in rows) == Set([0.08, 0.10])

    # AUROC is shared across models per error rate (property of the raw gap ranking).
    for er in unique(r.error_rate for r in rows)
        er_rows = filter(r -> r.error_rate == er, rows)
        aurocs = unique(r.auroc for r in er_rows)
        Test.@test length(aurocs) == 1
        au = only(aurocs)
        Test.@test !isnan(au)
        Test.@test au > 0.55                        # the gap is a usable signal (>> chance)
        # All three models present with finite calibration metrics.
        Test.@test Set(r.model for r in er_rows) == Set(["isotonic", "logistic", "binned"])
        for r in er_rows
            Test.@test isfinite(r.ece) && 0.0 <= r.ece <= 1.0
            Test.@test isfinite(r.brier) && 0.0 <= r.brier <= 1.0
        end
    end

    # --- Alignment contract + finite gaps on a real coverage-graph decode -------
    rng = Random.MersenneTwister(3)
    genome = String(rand(rng, ['A', 'C', 'G', 'T'], 500))
    k = 11
    reads = FASTX.FASTQ.Record[]
    truths = String[]
    for i in 1:40
        st = rand(rng, 1:(500 - 120 + 1))
        clean = genome[st:(st + 119)]
        obs = _inject_substitutions(clean, 0.08, rng)
        push!(reads, FASTX.FASTQ.Record("r$i", obs, String(fill('I', 120))))
        push!(truths, clean)
    end
    cov = build_coverage_fixture(reads, k)
    cfg = Mycelia.ViterbiCorrectionConfig(record_position_gaps = true, error_rate = 0.08)
    any_finite = false
    for (rec, clean) in zip(reads, truths)
        gt = try
            collect_gap_truth(cov, FASTX.sequence(String, rec), clean; config = cfg)
        catch
            continue
        end
        Test.@test length(gt.scores) == length(gt.labels)
        Test.@test all(isfinite, gt.scores)         # Inf sentinels already filtered
        Test.@test all(p -> p >= 1, gt.positions)
        any_finite |= !isempty(gt.scores)
    end
    Test.@test any_finite                            # the coverage graph yields finite gaps
end
