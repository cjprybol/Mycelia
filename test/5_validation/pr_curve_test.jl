import Test
import Mycelia

# Summarize a set of classifier operating points (a threshold GRID, i.e. an
# unordered cloud of (precision, recall) pairs) into a precision-recall curve and
# its area. Expected values below are hand-computed; see the doctstring of
# Mycelia.precision_recall_curve_summary for the all-points-interpolation AP.
Test.@testset "precision_recall_curve_summary" begin
    # Three non-dominated points: higher precision at lower recall.
    P = [0.9, 0.6, 0.3]
    R = [0.2, 0.5, 0.8]

    s = Mycelia.precision_recall_curve_summary(P, R)

    # Pareto frontier = all three, sorted by recall ascending.
    Test.@test s.frontier_recall == [0.2, 0.5, 0.8]
    Test.@test s.frontier_precision == [0.9, 0.6, 0.3]

    # All-points-interpolated average precision:
    #   (0.2-0)*0.9 + (0.5-0.2)*0.6 + (0.8-0.5)*0.3 = 0.18 + 0.18 + 0.09 = 0.45
    Test.@test isapprox(s.average_precision, 0.45; atol = 1e-9)

    # Trapezoid under the recall-sorted frontier:
    #   (0.5-0.2)*(0.9+0.6)/2 + (0.8-0.5)*(0.6+0.3)/2 = 0.225 + 0.135 = 0.36
    Test.@test isapprox(s.auc_trapezoid, 0.36; atol = 1e-9)

    # max F1 is at B=(0.6,0.5): F1 = 2*0.6*0.5/1.1 = 0.5454...
    Test.@test isapprox(s.max_f1, 2 * 0.6 * 0.5 / 1.1; atol = 1e-9)
    Test.@test s.best_index == 2
    Test.@test isapprox(s.best_precision, 0.6; atol = 1e-9)
    Test.@test isapprox(s.best_recall, 0.5; atol = 1e-9)
    Test.@test s.n_points == 3

    # A dominated point must not change AP (a worse point never raises the running
    # max precision), and must be excluded from the Pareto frontier.
    Pd = [0.9, 0.6, 0.3, 0.4]
    Rd = [0.2, 0.5, 0.8, 0.3]   # (0.4,0.3) is dominated by (0.6,0.5)
    sd = Mycelia.precision_recall_curve_summary(Pd, Rd)
    Test.@test isapprox(sd.average_precision, 0.45; atol = 1e-9)
    Test.@test sd.frontier_recall == [0.2, 0.5, 0.8]   # dominated point dropped
    Test.@test sd.n_points == 4

    # NaN points are dropped before any computation.
    Pn = [0.9, NaN, 0.6, 0.3]
    Rn = [0.2, 0.4, 0.5, 0.8]
    sn = Mycelia.precision_recall_curve_summary(Pn, Rn)
    Test.@test sn.n_points == 3
    Test.@test isapprox(sn.average_precision, 0.45; atol = 1e-9)

    # Single point: AP = area of the interpolated step from recall 0 to r at p;
    # trapezoid is undefined (<2 frontier points) -> NaN.
    s1 = Mycelia.precision_recall_curve_summary([0.8], [0.5])
    Test.@test isapprox(s1.average_precision, 0.5 * 0.8; atol = 1e-9)
    Test.@test isnan(s1.auc_trapezoid)
    Test.@test s1.best_index == 1

    # Empty (all dropped) -> NaN metrics, empty frontier, best_index 0.
    s0 = Mycelia.precision_recall_curve_summary(Float64[], Float64[])
    Test.@test s0.n_points == 0
    Test.@test isnan(s0.average_precision)
    Test.@test isempty(s0.frontier_recall)
    Test.@test s0.best_index == 0

    # Length mismatch is an error.
    Test.@test_throws Exception Mycelia.precision_recall_curve_summary([0.5, 0.6], [0.5])
end
