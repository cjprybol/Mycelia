# Unit test for the pre-registered paired-Wilcoxon analysis (beads td-59o7, td-9p91).
#
# This test carries the executability claim for td-59o7. The bead says explicitly:
# "confirm the pre-reg's paired test becomes executable — do not stop at 'column
# added'." So the central pair of assertions is:
#
#   * against a sweep CSV WITHOUT `seed`, the analysis RAISES and points at the
#     backfill tool (the pre-registered test is genuinely not runnable), and
#   * against the current schema WITH `seed`, it produces a real p-value, pairs
#     across seeds 42/123/456, and applies the pre-registration's decision rule.
#
# It also pins td-9p91's guard inside this analysis path, and cross-validates the
# statistics against scipy.stats.wilcoxon (reference values computed with
# scipy 1.16.3; see the comment at each fixture for the exact call).
#
# Dependency-free apart from CSV / DataFrames / Distributions / JSON (direct
# Mycelia deps). No network, no QUAST, no assembly.

import Test
import CSV
import DataFrames
import JSON

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rgv_paired_wilcoxon.jl"))

function _pwt_throws_with(f, needles::Vector{String})
    try
        f()
    catch e
        msg = sprint(showerror, e)
        for needle in needles
            Test.@test occursin(needle, msg)
        end
        return msg
    end
    Test.@test false
    return ""
end

# Build a sweep-shaped table. `with_seed=false` reproduces the pre-td-59o7 schema.
# `naive`/`iterative` give the per-cell metric values, one entry per (seed, cell).
function _pwt_frame(cells; with_seed = true, metric_source = "quast",
        min_contig = 500, metric = :quast_nga50, ok = true)
    rows = []
    for c in cells
        for (arm, value) in (("naive", c.naive), ("iterative", c.iterative))
            row = Dict{Symbol, Any}(
                :reference => get(c, :reference, "NC_001416"),
                :error_rate => get(c, :error_rate, 0.05),
                :readlen => get(c, :readlen, 150),
                :target_coverage => get(c, :target_coverage, 30.0),
                :k => 21,
                :arm => arm,
                :ok => get(c, :ok, ok),
                :metric_source => get(c, :metric_source, metric_source),
                :quast_min_contig => get(c, :min_contig, min_contig),
                metric => value
            )
            with_seed && (row[:seed] = get(c, :seed, 42))
            push!(rows, row)
        end
    end
    cols = collect(keys(rows[1]))
    return DataFrames.DataFrame([c => [r[c] for r in rows] for c in cols])
end

Test.@testset "RGV paired-Wilcoxon analysis" begin
    Test.@testset "Wilcoxon signed-rank vs scipy.stats.wilcoxon" begin
        # scipy: wilcoxon(A, method="exact") -> statistic 4.0, p 0.109375
        A = [1.5, -2.5, 3.5, 4.5, -0.5, 6.5, 7.5]
        ra = wilcoxon_signed_rank(A)
        Test.@test ra.n == 7
        Test.@test ra.n_zero_dropped == 0
        Test.@test ra.statistic == 4.0
        Test.@test occursin("exact", ra.method)
        Test.@test isapprox(ra.pvalue, 0.109375; atol = 1e-12)

        # scipy: wilcoxon(B, method="approx", correction=False)
        #        -> statistic 14.5, p 0.18371406577326088   (ties in |d| force approx)
        B = [1.0, -1.0, 2.0, 2.0, -3.0, 3.0, 4.0, -4.0, 5.0, 5.0]
        rb = wilcoxon_signed_rank(B)
        Test.@test rb.n == 10
        Test.@test rb.statistic == 14.5
        Test.@test occursin("normal approximation", rb.method)
        Test.@test occursin("tie correction", rb.method)
        Test.@test isapprox(rb.pvalue, 0.18371406577326088; atol = 1e-10)

        # scipy: wilcoxon(C, method="approx", zero_method="wilcox", correction=False)
        #        -> statistic 4.0, p 0.17177277759571152   (zeros dropped, ties present)
        C = [0.0, 2.0, -3.0, 4.0, 4.0, 0.0, -1.0, 5.0]
        rc = wilcoxon_signed_rank(C)
        Test.@test rc.n == 6
        Test.@test rc.n_zero_dropped == 2
        Test.@test rc.statistic == 4.0
        Test.@test isapprox(rc.pvalue, 0.17177277759571152; atol = 1e-10)

        # scipy: wilcoxon(D, method="exact") -> statistic 6.0, p 0.0068359375
        D = [3.1, -1.2, 4.7, 5.9, -0.3, 7.2, 8.8, 9.4, -2.6, 10.5, 11.1, 12.3]
        rd = wilcoxon_signed_rank(D)
        Test.@test rd.statistic == 6.0
        Test.@test isapprox(rd.pvalue, 0.0068359375; atol = 1e-12)

        # Degenerate inputs are reported as undefined, not as a null result.
        rz = wilcoxon_signed_rank([0.0, 0.0, 0.0])
        Test.@test rz.n == 0
        Test.@test rz.n_zero_dropped == 3
        Test.@test isnan(rz.pvalue)
        Test.@test occursin("undefined", rz.method)
        Test.@test isnan(wilcoxon_signed_rank(Float64[]).pvalue)

        # An all-positive set is the strongest one-directional signal available.
        rp = wilcoxon_signed_rank([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        Test.@test rp.w_minus == 0.0
        Test.@test isapprox(rp.pvalue, 2 / 2^6; atol = 1e-12)  # 2*(1/64)
    end

    Test.@testset "average_ranks handles ties as Wilcoxon requires" begin
        Test.@test average_ranks([10.0, 20.0, 30.0]) == [1.0, 2.0, 3.0]
        Test.@test average_ranks([5.0, 5.0, 9.0]) == [1.5, 1.5, 3.0]
        Test.@test average_ranks([7.0, 7.0, 7.0]) == [2.0, 2.0, 2.0]
        Test.@test average_ranks([3.0, 1.0, 2.0]) == [3.0, 1.0, 2.0]
    end

    Test.@testset "exact_signed_rank_pvalue endpoints" begin
        # n=1: W+ is 0 or 1 with probability 1/2 each, so the smallest attainable
        # two-sided p is 1.0 — a single pair can never be significant.
        Test.@test exact_signed_rank_pvalue(1, 1) == 1.0
        Test.@test exact_signed_rank_pvalue(1, 0) == 1.0
        # The distribution is symmetric about n(n+1)/4.
        for n in 2:8
            maxw = n * (n + 1) ÷ 2
            for w in 0:maxw
                Test.@test isapprox(exact_signed_rank_pvalue(n, w),
                    exact_signed_rank_pvalue(n, maxw - w); atol = 1e-12)
            end
        end
        # Extreme W+ for n=6: 2 * (1/64).
        Test.@test isapprox(exact_signed_rank_pvalue(6, 21), 2 / 64; atol = 1e-12)
    end

    Test.@testset "Benjamini-Hochberg matches the standard procedure" begin
        # Reference (statsmodels multipletests fdr_bh / the standard step-up):
        # [0.01, 0.04, 0.03, 0.2] -> [0.04, 0.05333333, 0.05333333, 0.2]
        adj = benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
        Test.@test isapprox(adj, [0.04, 0.4 / 7.5, 0.4 / 7.5, 0.2]; atol = 1e-8)
        Test.@test isapprox(adj[2], 0.05333333333; atol = 1e-8)
        # Monotone and never above 1.
        Test.@test all(adj .<= 1.0)
        Test.@test benjamini_hochberg([0.5])[1] == 0.5
        Test.@test benjamini_hochberg([0.9, 0.95]) == [0.95, 0.95]
        # NaN (undefined test) does not consume multiplicity budget.
        adj2 = benjamini_hochberg([0.01, NaN])
        Test.@test adj2[1] == 0.01      # m == 1, not 2
        Test.@test isnan(adj2[2])
        Test.@test all(isnan, benjamini_hochberg([NaN, NaN]))
    end

    Test.@testset "WITHOUT seed: the pre-registered test is NOT runnable" begin
        # This is the defect td-59o7 describes. Merged replicate rows are
        # indistinguishable, so pairs cannot be formed and the pre-registration's
        # paired-Wilcoxon rule over seeds 42/123/456 cannot be executed at all.
        df = _pwt_frame(
            [
                (naive = 1000.0, iterative = 1500.0),
                (naive = 1100.0, iterative = 1600.0)
            ];
            with_seed = false)
        Test.@test !("seed" in DataFrames.names(df))
        msg = _pwt_throws_with(() -> run_paired_analysis(df),
            ["no `seed` column",
                "NOT runnable",
                "42/123/456",
                "td-59o7",
                "rgv_seed_backfill.jl"])   # the message names the fix
        Test.@test occursin("paired", msg)
        # require_pairing_schema is the gate, and it also catches other key losses.
        _pwt_throws_with(
            () -> require_pairing_schema(DataFrames.select(
                _pwt_frame([(naive = 1.0, iterative = 2.0)]),
                DataFrames.Not(:reference))),
            ["missing pairing column", "reference"])
    end

    Test.@testset "WITH seed: the pre-registered test IS runnable" begin
        # Three pre-registered seeds x two error rates = 6 pairs, exactly the shape
        # the pre-registration's rule assumes. Seed is what separates the replicates.
        cells = []
        for (i, seed) in enumerate((42, 123, 456))
            for (j, err) in enumerate((0.01, 0.05))
                push!(cells,
                    (seed = seed, error_rate = err,
                        naive = 1000.0 + 10 * i + j,
                        iterative = 1400.0 + 13 * i + 2 * j))
            end
        end
        df = _pwt_frame(cells)
        Test.@test DataFrames.nrow(df) == 12          # 6 cells x 2 arms
        analysis = run_paired_analysis(df; metrics = (:quast_nga50,))

        # It ran, and it ran on the pre-registered seeds.
        Test.@test analysis.seeds == [42, 123, 456]
        r = analysis.results[1]
        Test.@test r.n_pairs == 6
        Test.@test r.n_dropped == 0
        Test.@test !isnan(r.test.pvalue)
        Test.@test !isnan(analysis.adjusted_p[1])
        # Every difference is positive here, so the exact two-sided p is 2/2^6.
        Test.@test isapprox(r.test.pvalue, 2 / 64; atol = 1e-12)
        Test.@test r.median_paired_difference > 0
        Test.@test r.median_relative_improvement > RGV_IMPROVEMENT_THRESHOLD
        Test.@test occursin("SUPPORTED", analysis.verdicts[1])
        Test.@test !occursin("PARTIALLY", analysis.verdicts[1])

        # The pairing key is what makes this work: drop seed from it and the six
        # cells collapse into duplicate rows per arm, which is refused.
        collapsed = DataFrames.select(df, DataFrames.Not(:seed))
        _pwt_throws_with(() -> run_paired_analysis(collapsed), ["no `seed` column"])
    end

    Test.@testset "decision rule follows the pre-registration" begin
        # Significant but small: PARTIALLY SUPPORTED (< 10% median improvement).
        small = _pwt_frame([(seed = s, error_rate = e,
                                naive = 1000.0, iterative = 1000.0 + d)
                            for (s, e, d) in [
            (42, 0.01, 10.0), (42, 0.05, 12.0), (123, 0.01, 11.0),
            (123, 0.05, 13.0), (456, 0.01, 9.0), (456, 0.05, 14.0)]])
        a_small = run_paired_analysis(small; metrics = (:quast_nga50,))
        Test.@test a_small.adjusted_p[1] < RGV_ALPHA
        Test.@test a_small.results[1].median_relative_improvement <
                   RGV_IMPROVEMENT_THRESHOLD
        Test.@test occursin("PARTIALLY SUPPORTED", a_small.verdicts[1])

        # Not significant: NOT SUPPORTED, reported as a null with effect size.
        mixed_signs = _pwt_frame([(seed = s, error_rate = e,
                                      naive = 1000.0, iterative = 1000.0 + d)
                                  for (s, e, d) in [
            (42, 0.01, 300.0), (42, 0.05, -280.0), (123, 0.01, 250.0),
            (123, 0.05, -320.0)]])
        a_null = run_paired_analysis(mixed_signs; metrics = (:quast_nga50,))
        Test.@test a_null.adjusted_p[1] >= RGV_ALPHA
        Test.@test occursin("NOT SUPPORTED", a_null.verdicts[1])
        Test.@test !isnan(a_null.results[1].median_paired_difference)  # effect size reported
    end

    Test.@testset "lower-is-better metrics are oriented correctly" begin
        # Misassemblies: FEWER is an improvement, so a negative paired difference
        # must yield a POSITIVE relative improvement.
        df = _pwt_frame(
            [(seed = s, error_rate = e, naive = 10.0, iterative = 4.0)
             for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
                (123, 0.05), (456, 0.01), (456, 0.05)]];
            metric = :quast_num_misassemblies)
        a = run_paired_analysis(df; metrics = (:quast_num_misassemblies,))
        r = a.results[1]
        Test.@test r.direction == :lower
        Test.@test r.median_paired_difference == -6.0            # fewer misassemblies
        Test.@test isapprox(r.median_relative_improvement, 0.6)  # a 60% improvement
        Test.@test occursin("SUPPORTED", a.verdicts[1])
    end

    Test.@testset "no directional claim where the pre-registration makes none" begin
        df = _pwt_frame(
            [(seed = s, error_rate = 0.01, naive = 1.0, iterative = 1.4)
             for s in (42, 123, 456)];
            metric = :quast_duplication_ratio)
        a = run_paired_analysis(df; metrics = (:quast_duplication_ratio,))
        Test.@test a.results[1].direction == :none
        Test.@test occursin("no directional claim", a.verdicts[1])
    end

    Test.@testset "td-9p91 guard fires inside THIS analysis path" begin
        # A realistic sweep CSV holds both QUAST-scored and degraded rows. Running
        # the analysis over them unfiltered must RAISE, not compare an
        # alignment-validated NGA50 against an internal size-ratio proxy.
        cells = [
            (seed = 42, error_rate = 0.01, naive = 1000.0, iterative = 1500.0,
                metric_source = "quast"),
            (seed = 123, error_rate = 0.01, naive = 900.0, iterative = 1400.0,
                metric_source = "internal:quast-failed"),
            (seed = 456, error_rate = 0.01, naive = 950.0, iterative = 1450.0,
                metric_source = "quast")
        ]
        df = _pwt_frame(cells)
        _pwt_throws_with(() -> run_paired_analysis(df; metrics = (:quast_nga50,)),
            ["MixedMetricDefinitionError",
                "RGV paired-Wilcoxon analysis",
                "internal:quast-failed", "quast"])

        # Selecting a definition explicitly is what makes the run legal.
        a = run_paired_analysis(df;
            metrics = (:quast_nga50,), metric_source = "quast")
        Test.@test a.results[1].n_pairs == 2
        Test.@test a.n_dropped_metric_source == 2
        Test.@test Dict(a.definition)[:metric_source] == ["quast"]

        # Mixed quast_min_contig is rejected on the same grounds: the T4 fix
        # changes the threshold, so rows either side of it are not comparable.
        thresholds = _pwt_frame([
            (seed = 42, error_rate = 0.01, naive = 1000.0, iterative = 1500.0,
                min_contig = 16_890),
            (seed = 123, error_rate = 0.01, naive = 900.0, iterative = 1400.0,
                min_contig = 5_000)
        ])
        _pwt_throws_with(
            () -> run_paired_analysis(thresholds; metrics = (:quast_nga50,)),
            ["quast_min_contig", "16890", "5000"])
    end

    Test.@testset "unpairable cells are dropped and reported, never imputed" begin
        base = [(seed = 42, error_rate = 0.01, naive = 1000.0, iterative = 1500.0),
            (seed = 123, error_rate = 0.01, naive = 900.0, iterative = 1400.0)]
        df = _pwt_frame(base)

        # A cell with only one arm cannot be paired.
        orphan = DataFrames.DataFrame(df[1:1, :])
        orphan[1, :seed] = 456
        with_orphan = vcat(df, orphan)
        a = run_paired_analysis(with_orphan; metrics = (:quast_nga50,))
        Test.@test a.results[1].n_pairs == 2
        Test.@test a.results[1].n_dropped == 1
        Test.@test occursin("no `iterative` row", a.results[1].dropped[1].reason)

        # A `missing` metric on one arm is dropped, not filled in.
        holey = DataFrames.DataFrame(df)
        holey[!, :quast_nga50] = Vector{Union{Missing, Float64}}(holey.quast_nga50)
        holey[2, :quast_nga50] = missing
        a2 = run_paired_analysis(holey; metrics = (:quast_nga50,))
        Test.@test a2.results[1].n_pairs == 1
        Test.@test occursin("`missing`", a2.results[1].dropped[1].reason)

        # Duplicate rows for one arm mean the pairing key is not unique — refused
        # for that cell rather than silently averaged.
        dup = vcat(df, DataFrames.DataFrame(df[1:1, :]))
        a3 = run_paired_analysis(dup; metrics = (:quast_nga50,))
        Test.@test a3.results[1].n_pairs == 1
        Test.@test occursin("duplicate rows", a3.results[1].dropped[1].reason)

        # ok=false rows are dropped by default and counted.
        notok = _pwt_frame(base; ok = false)
        _pwt_throws_with(() -> run_paired_analysis(notok; metrics = (:quast_nga50,)),
            ["no rows left after filtering"])
        a4 = run_paired_analysis(notok; metrics = (:quast_nga50,), keep_not_ok = true)
        Test.@test a4.results[1].n_pairs == 2
    end

    Test.@testset "zero-control pairs are excluded from % improvement and counted" begin
        # A zero baseline makes the RATIO undefined; dropping it silently would bias
        # the statistic the 10% threshold is applied to.
        df = _pwt_frame([
            (seed = 42, error_rate = 0.01, naive = 0.0, iterative = 1500.0),
            (seed = 123, error_rate = 0.01, naive = 1000.0, iterative = 1500.0),
            (seed = 456, error_rate = 0.01, naive = 1000.0, iterative = 1600.0)
        ])
        a = run_paired_analysis(df; metrics = (:quast_nga50,))
        r = a.results[1]
        Test.@test r.n_pairs == 3                                      # all pairs tested
        Test.@test r.n_zero_control_excluded_from_relative == 1
        Test.@test !isnan(r.median_relative_improvement)
        Test.@test r.median_paired_difference == 600.0                 # uses all 3
    end

    Test.@testset "multiple CSVs concatenate into one pairable table" begin
        # This is the parallelization payoff of having `seed`: seed-sharded runs
        # from different hosts merge into one table the paired test can use.
        mktempdir() do dir
            paths = String[]
            for (i, seed) in enumerate((42, 123, 456))
                df = _pwt_frame([(seed = seed, error_rate = 0.01,
                    naive = 1000.0 + i, iterative = 1500.0 + 2i)])
                p = joinpath(dir, "shard_seed$(seed).csv")
                CSV.write(p, df)
                push!(paths, p)
            end
            combined = load_sweep_csvs(paths)
            Test.@test DataFrames.nrow(combined) == 6
            a = run_paired_analysis(combined; metrics = (:quast_nga50,))
            Test.@test a.seeds == [42, 123, 456]
            Test.@test a.results[1].n_pairs == 3
            _pwt_throws_with(() -> load_sweep_csvs([joinpath(dir, "nope.csv")]),
                ["CSV not found"])
            _pwt_throws_with(() -> load_sweep_csvs(String[]), ["no CSV supplied"])
        end
    end

    Test.@testset "report and JSON record what was actually done" begin
        mktempdir() do dir
            df = _pwt_frame([(seed = s, error_rate = e, naive = 1000.0,
                                 iterative = 1600.0)
                             for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
                (123, 0.05), (456, 0.01), (456, 0.05)]])
            a = run_paired_analysis(df; metrics = (:quast_nga50,))
            report = write_paired_report(joinpath(dir, "report.md"), a;
                csv_paths = ["fixture.csv"])
            text = read(report, String)
            Test.@test occursin("paired Wilcoxon", text)
            Test.@test occursin("42, 123, 456", text)            # seeds used
            Test.@test occursin("`seed`", text)                   # pairing key shown
            Test.@test occursin("quast_nga50", text)
            Test.@test occursin("SUPPORTED", text)
            Test.@test occursin("metric_source_guard.jl", text)   # definition provenance

            payload = paired_analysis_json(a; csv_paths = ["fixture.csv"])
            json_path = joinpath(dir, "results.json")
            open(json_path, "w") do io
                JSON.print(io, payload, 2)
            end
            round = JSON.parsefile(json_path)
            Test.@test round["seeds"] == [42, 123, 456]
            Test.@test round["pair_keys"][end] == "seed"
            Test.@test round["metrics"][1]["metric"] == "quast_nga50"
            Test.@test round["metrics"][1]["n_pairs"] == 6
            Test.@test round["alpha"] == RGV_ALPHA
            Test.@test round["metric_definition"]["metric_source"] == ["quast"]
        end
    end

    Test.@testset "pairing key and metrics match the pre-registration" begin
        # seed must be part of the pairing key, or the pre-registered replicate
        # design is not represented in the analysis.
        Test.@test :seed in RGV_PAIR_KEYS
        # H1's primary metrics: NGA50 and misassembly count.
        Test.@test RGV_DEFAULT_METRICS == (:quast_nga50, :quast_num_misassemblies)
        Test.@test RGV_ALPHA == 0.05
        Test.@test RGV_IMPROVEMENT_THRESHOLD == 0.10
        Test.@test RGV_METRIC_DIRECTION[:quast_nga50] == :higher
        Test.@test RGV_METRIC_DIRECTION[:quast_num_misassemblies] == :lower
    end

    Test.@testset "an unknown metric column raises" begin
        df = _pwt_frame([(seed = 42, error_rate = 0.01, naive = 1.0, iterative = 2.0)])
        _pwt_throws_with(() -> run_paired_analysis(df; metrics = (:not_a_column,)),
            ["metric column `not_a_column` not present"])
    end

    Test.@testset "C16: a shard lacking `seed` must NOT create pseudo-replicates" begin
        # The worst failure available to this file, because it makes the result look
        # BETTER. `load_sweep_csvs` concatenates with `cols = :union`, so a shard
        # predating the `seed` column contributes rows with `seed = missing`. A
        # presence-only schema check passed, `build_pairs` stringified the key, and
        # `missing` became a distinct pseudo-seed that paired normally: the same
        # physical runs supplied twice reported n = 12 for 6 replicates and deflated
        # p by ~27x.
        cells = [(seed = 42, error_rate = e, naive = 1000.0 + i, iterative = 1300.0 + i)
                 for (i, e) in enumerate((0.01, 0.02, 0.03, 0.04, 0.05, 0.06))]
        withseed = _pwt_frame(cells)
        noseed = _pwt_frame(cells; with_seed = false)

        # Correct baseline: 6 real cells.
        a_ok = run_paired_analysis(withseed; metrics = (:quast_nga50,))
        Test.@test a_ok.results[1].n_pairs == 6

        merged = vcat(withseed, noseed; cols = :union)
        Test.@test "seed" in DataFrames.names(merged)      # column IS present...
        Test.@test count(ismissing, merged.seed) == 12     # ...but half is missing
        msg = _pwt_throws_with(
            () -> run_paired_analysis(merged; metrics = (:quast_nga50,)),
            ["missing", "pairing column", "`seed`",
                "pseudo-level", "inflating n", "deflating p"])
        Test.@test occursin("rgv_seed_backfill.jl", msg)

        # Every pairing key is checked, not just seed.
        holed = DataFrames.DataFrame(withseed)
        holed[!, :error_rate] = Vector{Union{Missing, Float64}}(holed.error_rate)
        holed[1, :error_rate] = missing
        _pwt_throws_with(() -> run_paired_analysis(holed; metrics = (:quast_nga50,)),
            ["`error_rate`"])

        # A `missing` arm belongs to neither side of a pair.
        badarm = DataFrames.DataFrame(withseed)
        badarm[!, :arm] = Vector{Union{Missing, String}}(badarm.arm)
        badarm[1, :arm] = missing
        _pwt_throws_with(() -> run_paired_analysis(badarm; metrics = (:quast_nga50,)),
            ["`missing` `arm`"])
    end

    Test.@testset "C9: a significant HARM is FALSIFIED, not partially supported" begin
        # The pre-registration's H1 row: "Falsified if p >= 0.05, THE DIRECTION IS
        # NEGATIVE, or ...". `decide` branched only on magnitude, so a treatment
        # worse on every pair reported as PARTIALLY SUPPORTED — the exact class of
        # misreport this PR exists to prevent.
        harm = _pwt_frame([(seed = s, error_rate = e, naive = 1000.0,
                               iterative = 1000.0 - d)
                           for (s, e, d) in [
            (42, 0.01, 400.0), (42, 0.05, 450.0), (123, 0.01, 500.0),
            (123, 0.05, 420.0), (456, 0.01, 470.0), (456, 0.05, 430.0),
            (42, 0.10, 460.0), (123, 0.10, 440.0), (456, 0.10, 480.0)]])
        a = run_paired_analysis(harm; metrics = (:quast_nga50,))
        Test.@test a.adjusted_p[1] < RGV_ALPHA                       # significant
        Test.@test a.results[1].median_relative_improvement < 0      # and WORSE
        Test.@test occursin("FALSIFIED", a.verdicts[1])
        Test.@test occursin("direction is", a.verdicts[1])
        Test.@test !occursin("SUPPORTED (", a.verdicts[1])
        # A positive-but-small effect is still PARTIALLY SUPPORTED, not falsified.
        small = _pwt_frame([(seed = s, error_rate = e, naive = 1000.0,
                                iterative = 1000.0 + d)
                            for (s, e, d) in [
            (42, 0.01, 10.0), (42, 0.05, 12.0), (123, 0.01, 11.0),
            (123, 0.05, 13.0), (456, 0.01, 9.0), (456, 0.05, 14.0)]])
        Test.@test occursin("PARTIALLY SUPPORTED",
            run_paired_analysis(small; metrics = (:quast_nga50,)).verdicts[1])
    end

    Test.@testset "C12: an undefined effect size is not silently '<= 10%'" begin
        # `NaN > 0.1` is false, so an undefined relative improvement fell through to
        # the "<= 10%" branch. Routine for misassembly counts, where a clean control
        # is zero on every pair.
        zc = _pwt_frame([(seed = s, error_rate = e, naive = 0.0, iterative = 3.0)
                         for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]])
        a = run_paired_analysis(zc; metrics = (:quast_nga50,))
        Test.@test a.adjusted_p[1] < RGV_ALPHA
        Test.@test isnan(a.results[1].median_relative_improvement)
        Test.@test a.results[1].n_zero_control_excluded_from_relative == 6
        Test.@test occursin("INDETERMINATE", a.verdicts[1])
        Test.@test !occursin("<= 10%", a.verdicts[1])
    end

    Test.@testset "C5: the mixed-src override is recorded in the deliverable" begin
        mixed = _pwt_frame([
            (seed = 42, error_rate = 0.01, naive = 1000.0, iterative = 1500.0,
                metric_source = "quast"),
            (seed = 123, error_rate = 0.01, naive = 900.0, iterative = 1400.0,
                metric_source = "internal:quast-failed"),
            (seed = 456, error_rate = 0.01, naive = 950.0, iterative = 1450.0,
                metric_source = "quast")])
        a = run_paired_analysis(mixed;
            metrics = (:quast_nga50,), allow_mixed_src = true)
        Test.@test a.allow_mixed_src
        Test.@test a.override_bound                       # it actually BOUND
        Test.@test "metric_source" in a.mixed_axes

        mktempdir() do dir
            text = read(write_paired_report(joinpath(dir, "r.md"), a;
                csv_paths = ["f.csv"]), String)
            # The artifact must say the guard was overridden...
            Test.@test occursin("GUARD OVERRIDDEN", text)
            Test.@test occursin("not validation-grade", lowercase(text))
            # ...and must NOT go on to assert the guard was enforced.
            Test.@test !occursin("enforced it — no override was used", text)
            payload = paired_analysis_json(a; csv_paths = ["f.csv"])
            Test.@test payload["metric_definition_override_bound"] == true
            Test.@test payload["allow_mixed_src"] == true
            Test.@test payload["exploratory"] == true
        end

        # An override that does NOT bind is reported differently — "available but
        # unused" is a different provenance claim from "load-bearing".
        clean = _pwt_frame([(seed = s, error_rate = 0.01, naive = 1000.0,
                                iterative = 1500.0) for s in (42, 123, 456)])
        b = run_paired_analysis(clean; metrics = (:quast_nga50,), allow_mixed_src = true)
        Test.@test b.allow_mixed_src
        Test.@test !b.override_bound
        mktempdir() do dir
            text = read(write_paired_report(joinpath(dir, "r.md"), b;
                csv_paths = ["f.csv"]), String)
            Test.@test occursin("did NOT bind", text)
            Test.@test !occursin("GUARD OVERRIDDEN", text)
        end
    end

    Test.@testset "C1/I3: the report does not claim to be a pre-registered test" begin
        df = _pwt_frame([(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                         for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]])
        a = run_paired_analysis(df; metrics = (:quast_nga50,))
        mktempdir() do dir
            text = read(write_paired_report(joinpath(dir, "r.md"), a;
                csv_paths = ["f.csv"]), String)
            # H1 is Viterbi-DP-vs-greedy; this compares correctors. Attributing it
            # to H1 would put an unregistered comparison into a manuscript as
            # confirmatory.
            Test.@test occursin("EXPLORATORY", text)
            Test.@test occursin("not a pre-registered test", text)
            Test.@test occursin("Viterbi DP vs greedy", text)
            Test.@test !occursin("(H1 test procedure + shared analysis rules)", text)
            # I3: the axes actually swept must sit next to the numbers.
            Test.@test occursin("Error rates swept", text)
            Test.@test occursin("0.01, 0.05", text)
            Test.@test occursin("~360 paired comparisons", text)
        end
    end
end
