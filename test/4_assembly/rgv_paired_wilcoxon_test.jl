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

# The WRITER of the provenance labels this file's reader consumes, loaded so the
# cross-file contract can be tested by round-trip rather than by grepping source.
# Namespaced into a module because both scripts define `main` at top level.
module _PWTBackfill
include(joinpath(@__DIR__, "..", "..", "benchmarking", "rgv_seed_backfill.jl"))
end

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
#
# `provenance=nothing` (the default) reproduces the shape of EVERY sweep CSV the
# documented backfill/analyse workflow produced before the provenance axis existed:
# definition columns present, no `*_provenance` siblings. That is a table whose
# provenance cannot be determined, NOT a table shown to be observed, so callers
# that want the measured case must ask for it explicitly.
function _pwt_frame(cells; with_seed = true, metric_source = "quast",
        min_contig = 500, metric = :quast_nga50, ok = true, provenance = nothing)
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
            if provenance !== nothing
                p = get(c, :provenance, provenance)
                row[:metric_source_provenance] = p
                row[:quast_min_contig_provenance] = p
            end
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

    Test.@testset "C6: multi-value flag parsing (glob and repeated forms)" begin
        # The fix has to be exercised, or the header's documented glob invocation
        # can silently regress to reading one shard again.
        saved = copy(ARGS)
        try
            empty!(ARGS)
            append!(ARGS, ["--csv", "a.csv", "--csv", "b.csv", "--metric-source", "quast"])
            Test.@test _pw_args("--csv") == ["a.csv", "b.csv"]          # repeated flag
            empty!(ARGS)
            append!(ARGS, ["--csv", "a.csv", "b.csv", "c.csv", "--metrics", "n50"])
            Test.@test _pw_args("--csv") == ["a.csv", "b.csv", "c.csv"] # shell glob
            Test.@test _pw_arg("--metrics") == "n50"                    # single-value intact
            empty!(ARGS)
            append!(ARGS, ["--metric-source", "quast"])
            Test.@test isempty(_pw_args("--csv"))                       # absent flag
        finally
            empty!(ARGS);
            append!(ARGS, saved)
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

        # This is the analysis that carries NaN into the payload, so it is the one
        # that must round-trip through JSON — the all-finite fixture elsewhere
        # cannot pin `_json_num`. Non-finite values must serialize as `null`
        # (which the figure script already reads as missing), not as a bare NaN
        # token that no strict JSON parser accepts.
        mktempdir() do dir
            payload = paired_analysis_json(a; csv_paths = ["f.csv"])
            Test.@test payload["metrics"][1]["median_relative_improvement"] === nothing
            jp = joinpath(dir, "results.json")
            open(jp, "w") do io
                JSON.print(io, payload, 2)
            end
            text = read(jp, String)
            Test.@test !occursin("NaN", text)
            round = JSON.parsefile(jp)                      # must parse back
            Test.@test round["metrics"][1]["median_relative_improvement"] === nothing
            # And the figure can consume it without a NaN blowing up the axes.
            Test.@test round["metrics"][1]["n_pairs"] == 6
        end
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
            text = read(
                write_paired_report(joinpath(dir, "r.md"), a;
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
            text = read(
                write_paired_report(joinpath(dir, "r.md"), b;
                    csv_paths = ["f.csv"]), String)
            Test.@test occursin("did NOT bind", text)
            Test.@test !occursin("GUARD OVERRIDDEN", text)
        end
    end

    Test.@testset "C: an absent `ok` column is not 'every row verified ok=true'" begin
        # `report.md` emitted the byte-identical line
        #   `- Rows: 12 read, 12 analyzed (0 dropped for ok=false, ...)`
        # both when the column was present and all-true, and when it was ABSENT
        # entirely. "0 dropped for ok=false" reads as "the check ran and every row
        # passed"; it was printed verbatim when the check could not run at all, and
        # `results.json` carried neither the count nor a presence flag, so no
        # consumer could recover the distinction from either artifact.
        # `metric_source_guard.jl` applies the opposite doctrine to the definition
        # columns: unverifiable is not the same as verified.
        cells = [(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                 for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]]
        with_ok = _pwt_frame(cells)                       # ok=true on every row
        no_ok = DataFrames.select(with_ok, DataFrames.Not(:ok))

        a_ok = run_paired_analysis(with_ok; metrics = (:quast_nga50,))
        a_no = run_paired_analysis(no_ok; metrics = (:quast_nga50,))
        # Identical numbers — the difference is what is KNOWN about them.
        Test.@test a_ok.n_rows_analyzed == a_no.n_rows_analyzed == 12
        Test.@test a_ok.ok_column_present && a_ok.ok_filter_applied
        Test.@test !a_no.ok_column_present && !a_no.ok_filter_applied

        mktempdir() do dir
            t_ok = read(
                write_paired_report(joinpath(dir, "ok.md"), a_ok;
                    csv_paths = ["f.csv"]), String)
            t_no = read(
                write_paired_report(joinpath(dir, "no.md"), a_no;
                    csv_paths = ["f.csv"]), String)
            Test.@test occursin("0 dropped for ok=false", t_ok)
            # Scoped to the `ok` clause: bare "ABSENT" also matches unrelated
            # disclosures elsewhere in the report and would silently stop testing
            # this one.
            Test.@test !occursin("`ok` column ABSENT", t_ok)
            # The unrunnable case must SAY it was unrunnable, and must not claim a
            # count of failures it never looked for.
            Test.@test occursin("`ok` column ABSENT", t_no)
            Test.@test occursin("could NOT run", t_no)
            Test.@test !occursin("dropped for ok=false", t_no)

            # And the distinction must survive into the machine-readable artifact.
            p_ok = paired_analysis_json(a_ok; csv_paths = ["f.csv"])
            p_no = paired_analysis_json(a_no; csv_paths = ["f.csv"])
            Test.@test p_ok["ok_column_present"] == true
            Test.@test p_no["ok_column_present"] == false
            Test.@test p_ok["ok_filter_applied"] != p_no["ok_filter_applied"]
            Test.@test p_ok["n_dropped_not_ok"] == 0
            Test.@test p_ok["n_dropped_metric_source"] == 0
        end

        # `--keep-not-ok` is a THIRD state: the column is present and the check was
        # deliberately disabled. That is not "0 failures" either.
        a_keep = run_paired_analysis(with_ok;
            metrics = (:quast_nga50,), keep_not_ok = true)
        Test.@test a_keep.ok_column_present
        Test.@test !a_keep.ok_filter_applied
        mktempdir() do dir
            t = read(
                write_paired_report(joinpath(dir, "k.md"), a_keep;
                    csv_paths = ["f.csv"]), String)
            Test.@test occursin("filter DISABLED", t)
            Test.@test !occursin("`ok` column ABSENT", t)
        end

        # Sub-finding: `coalesce(missing, false)` counted an UNKNOWN status as an
        # observed failure. "The harness said this run failed" and "we do not know
        # whether it succeeded" are different claims and are counted separately.
        holey = DataFrames.DataFrame(with_ok)
        holey[!, :ok] = Vector{Union{Missing, Bool}}(holey.ok)
        holey[1, :ok] = missing
        holey[3, :ok] = false
        a_h = run_paired_analysis(holey; metrics = (:quast_nga50,))
        Test.@test a_h.n_dropped_ok_missing == 1
        Test.@test a_h.n_dropped_not_ok == 1          # NOT 2
        Test.@test a_h.n_rows_analyzed == 10
        mktempdir() do dir
            t = read(
                write_paired_report(joinpath(dir, "h.md"), a_h;
                    csv_paths = ["f.csv"]), String)
            Test.@test occursin("1 dropped for ok=false", t)
            Test.@test occursin("1 dropped for ok=missing", t)
            p = paired_analysis_json(a_h; csv_paths = ["f.csv"])
            Test.@test p["n_dropped_not_ok"] == 1
            Test.@test p["n_dropped_ok_missing"] == 1
        end
    end

    Test.@testset "C: a non-boolean `ok` column is not a check that found failures" begin
        # `keep = coalesce.(okcol, false) .== true` is vacuously FALSE for a `String`
        # element, and the counter split then attributes every such row to
        # `n_dropped_not_ok` (only `ismissing` rows are subtracted out). So the
        # report asserted the harness had OBSERVED N failed runs whose cells
        # literally read "TRUE". Those runs succeeded: they are dropped, `n` shrinks,
        # the p-value moves, and the "(status unknown, not an observed failure)"
        # phrasing rules out the correct explanation.
        #
        # The trigger is ordinary. `CSV.read` infers a string type for an `ok` column
        # containing `true,false,NA` — one sentinel poisons the whole column — and
        # merging that shard with a Bool-typed shard under `cols = :union` produces
        # exactly the mixed column below.
        cells = [(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                 for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]]
        base = _pwt_frame(cells; provenance = "observed")
        n = DataFrames.nrow(base)

        # A merged two-shard table: the string shard's rows all read as successes.
        merged = DataFrames.DataFrame(base)
        merged[!, :ok] = Union{Bool, String}[i <= 6 ? "TRUE" : true for i in 1:n]
        a_mix = run_paired_analysis(merged; metrics = (:quast_nga50,))
        Test.@test a_mix.ok_column_present
        Test.@test !a_mix.ok_column_runnable
        Test.@test !a_mix.ok_filter_applied
        # The successful rows are NOT dropped, and nothing is reported as an
        # observed failure.
        Test.@test a_mix.n_rows_analyzed == n
        Test.@test a_mix.n_dropped_not_ok == 0
        Test.@test a_mix.n_dropped_ok_missing == 0

        # A wholly string-typed column is the same epistemic state. Under the old
        # filter every row dropped and the analysis died with "no rows left after
        # filtering", which at least failed loudly — the mixed case above is the one
        # that silently moved the p-value.
        stringy = DataFrames.DataFrame(base)
        stringy[!, :ok] = fill("TRUE", n)
        a_str = run_paired_analysis(stringy; metrics = (:quast_nga50,))
        Test.@test !a_str.ok_column_runnable
        Test.@test a_str.n_rows_analyzed == n

        mktempdir() do dir
            for (name, a) in (("mix", a_mix), ("str", a_str))
                text = read(
                    write_paired_report(joinpath(dir, "$name.md"), a;
                        csv_paths = ["f.csv"]), String)
                Test.@test occursin("NOT a boolean", text)
                Test.@test occursin("could NOT run", text)
                # It must never render as a count of observed failures.
                Test.@test !occursin("dropped for ok=false", text)
                Test.@test !occursin("`ok` column ABSENT", text)
                p = paired_analysis_json(a; csv_paths = ["f.csv"])
                Test.@test p["ok_column_present"] == true
                Test.@test p["ok_column_runnable"] == false
                Test.@test p["ok_filter_applied"] == false
                Test.@test occursin("String", p["ok_column_eltype"])
            end
        end

        # Positive control: an INTEGER 0/1 column really can express success, so it
        # must stay runnable and keep filtering. Widening the runnable predicate
        # must not quietly stop checking a schema that works today.
        ints = DataFrames.DataFrame(base)
        okints = Union{Missing, Int}[1 for _ in 1:n]
        okints[3] = 0
        okints[4] = missing
        ints[!, :ok] = okints
        a_int = run_paired_analysis(ints; metrics = (:quast_nga50,))
        Test.@test a_int.ok_column_runnable
        Test.@test a_int.ok_filter_applied
        Test.@test a_int.n_dropped_not_ok == 1
        Test.@test a_int.n_dropped_ok_missing == 1
        Test.@test a_int.n_rows_analyzed == n - 2
        # Bool remains runnable, and an ABSENT column is still its own state.
        a_bool = run_paired_analysis(base; metrics = (:quast_nga50,))
        Test.@test a_bool.ok_column_runnable
        Test.@test a_bool.ok_column_eltype == "Bool"
        a_absent = run_paired_analysis(
            DataFrames.select(base, DataFrames.Not(:ok)); metrics = (:quast_nga50,))
        Test.@test !a_absent.ok_column_present
        Test.@test !a_absent.ok_column_runnable
        Test.@test paired_analysis_json(a_absent;
            csv_paths = ["f.csv"])["ok_column_eltype"] === nothing
    end

    Test.@testset "B: an OPERATOR-ASSERTED definition cannot earn the guard assurance" begin
        # `rgv_seed_backfill.jl --metric-source quast` writes a definition on
        # operator say-so, with nothing to verify it against. The guard can then
        # only check that the asserted labels agree with EACH OTHER, and a
        # backfilled label agrees with itself by construction — so a legacy table
        # that genuinely mixed an alignment-validated NGA50 with a size-ratio proxy
        # passes, and the report went on to assert the guard had enforced
        # single-definition-ness. That is the same laundering the `--allow-mixed-src`
        # disclosure exists to prevent, reached by a different route.
        cells = [(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                 for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]]
        # Both tables must carry AFFIRMATIVE provenance on every definition axis:
        # the assurance is granted for observed provenance, not for the absence of a
        # contrary signal (see the "undeterminable" testset below).
        observed = _pwt_frame(cells; provenance = "observed")
        asserted = _pwt_frame(cells; provenance = "observed")
        asserted[!, :metric_source_provenance] = fill("operator-asserted (backfill)",
            DataFrames.nrow(asserted))

        a_obs = run_paired_analysis(observed; metrics = (:quast_nga50,))
        a_ast = run_paired_analysis(asserted; metrics = (:quast_nga50,))
        # The two tables produce IDENTICAL definition summaries — that is the point.
        Test.@test Dict(a_obs.definition) == Dict(a_ast.definition)
        Test.@test !a_obs.definition_operator_asserted
        Test.@test a_ast.definition_operator_asserted
        Test.@test a_ast.operator_asserted_axes == ["metric_source"]
        Test.@test a_ast.n_rows_operator_asserted == DataFrames.nrow(asserted)

        mktempdir() do dir
            t_obs = read(
                write_paired_report(joinpath(dir, "obs.md"), a_obs;
                    csv_paths = ["f.csv"]), String)
            t_ast = read(
                write_paired_report(joinpath(dir, "ast.md"), a_ast;
                    csv_paths = ["f.csv"]), String)
            # A measured definition still earns the unqualified assurance...
            Test.@test occursin("enforced it — no override was used", t_obs)
            Test.@test !occursin("OPERATOR-ASSERTED", t_obs)
            # ...an asserted one must not, and must say so before the numbers.
            Test.@test !occursin("enforced it — no override was used", t_ast)
            Test.@test occursin("OPERATOR-ASSERTED, NOT OBSERVED", t_ast)
            Test.@test occursin("operator assertion", t_ast)

            p_ast = paired_analysis_json(a_ast; csv_paths = ["f.csv"])
            Test.@test p_ast["metric_definition_operator_asserted"] == true
            Test.@test p_ast["operator_asserted_definition_axes"] == ["metric_source"]
            Test.@test p_ast["n_rows_operator_asserted_definition"] == 12
            p_obs = paired_analysis_json(a_obs; csv_paths = ["f.csv"])
            Test.@test p_obs["metric_definition_operator_asserted"] == false
            Test.@test isempty(p_obs["operator_asserted_definition_axes"])
        end

        # ONE asserted row is enough: the rows the backfill did not touch stay
        # labelled observed, but the aggregate still spans a claim nobody verified.
        partial = _pwt_frame(cells; provenance = "observed")
        prov = fill("observed", DataFrames.nrow(partial))
        prov[1] = "operator-asserted (backfill)"
        partial[!, :metric_source_provenance] = prov
        a_part = run_paired_analysis(partial; metrics = (:quast_nga50,))
        Test.@test a_part.definition_operator_asserted
        Test.@test a_part.n_rows_operator_asserted == 1

        # The claim is evaluated over the rows actually ANALYZED. An assertion on
        # rows the `--metric-source` filter removed does not taint what is left.
        filtered = _pwt_frame(
            [
                (seed = 42, error_rate = 0.01, naive = 1000.0, iterative = 1500.0,
                    metric_source = "quast"),
                (seed = 123, error_rate = 0.01, naive = 900.0, iterative = 1400.0,
                    metric_source = "internal:quast-failed"),
                (seed = 456, error_rate = 0.01, naive = 950.0, iterative = 1450.0,
                    metric_source = "quast")];
            provenance = "observed")
        filtered[!, :metric_source_provenance] = [src == "quast" ? "observed" :
                                                  "operator-asserted (backfill)"
                                                  for src in filtered.metric_source]
        a_filt = run_paired_analysis(filtered;
            metrics = (:quast_nga50,), metric_source = "quast")
        Test.@test a_filt.n_rows_operator_asserted == 0
        Test.@test !a_filt.definition_operator_asserted
    end

    Test.@testset "B: provenance that CANNOT be determined must not read as observed" begin
        # The two-state reader was a BLACKLIST: anything not containing
        # "operator-asserted" counted as observed. So it failed OPEN in exactly the
        # cases it exists for, and reported `metric_definition_operator_asserted:
        # false` — an affirmative claim the definition was MEASURED — over tables
        # that carry no evidence either way.
        cells = [(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                 for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]]

        # (1) NO sibling provenance column — the shape of every sweep CSV the
        # documented backfill -> analyse workflow produced before this PR.
        none = _pwt_frame(cells)
        a_none = run_paired_analysis(none; metrics = (:quast_nga50,))
        Test.@test a_none.definition_provenance_undeterminable
        Test.@test a_none.undeterminable_definition_axes ==
                   ["metric_source", "quast_min_contig"]
        # The existing boolean keeps its meaning — true iff a row was ASSERTED — so
        # the third state is added alongside it rather than overloaded onto it.
        Test.@test !a_none.definition_operator_asserted
        Test.@test a_none.n_rows_operator_asserted == 0

        # (2) A provenance column carrying a label the reader does not understand.
        # This is STRONGER evidence that something wrote provenance metadata than an
        # absent column is, yet a blacklist treated it more permissively.
        drifted = _pwt_frame(cells; provenance = "backfilled-by-operator")
        a_drift = run_paired_analysis(drifted; metrics = (:quast_nga50,))
        Test.@test a_drift.definition_provenance_undeterminable
        Test.@test !a_drift.definition_operator_asserted

        # (3) A provenance column that is present but `missing` on some rows.
        holey = _pwt_frame(cells; provenance = "observed")
        holey[!, :metric_source_provenance] = Vector{Union{Missing, String}}(
            holey.metric_source_provenance)
        holey[1, :metric_source_provenance] = missing
        a_holey = run_paired_analysis(holey; metrics = (:quast_nga50,))
        Test.@test a_holey.undeterminable_definition_axes == ["metric_source"]
        Test.@test !a_holey.definition_operator_asserted

        # (4) One axis affirmatively observed, the other with no sibling column: the
        # assurance is per-axis, so a partial disclosure does not earn it.
        partial = _pwt_frame(cells; provenance = "observed")
        DataFrames.select!(partial,
            DataFrames.Not(:quast_min_contig_provenance))
        a_partial = run_paired_analysis(partial; metrics = (:quast_nga50,))
        Test.@test a_partial.undeterminable_definition_axes == ["quast_min_contig"]

        # (5) The control: AFFIRMATIVELY observed on every present axis is the only
        # state that earns the unqualified assurance.
        observed = _pwt_frame(cells; provenance = "observed")
        a_obs = run_paired_analysis(observed; metrics = (:quast_nga50,))
        Test.@test !a_obs.definition_provenance_undeterminable
        Test.@test isempty(a_obs.undeterminable_definition_axes)

        mktempdir() do dir
            for (name, a) in (("none", a_none), ("drift", a_drift),
                ("holey", a_holey), ("partial", a_partial))
                text = read(
                    write_paired_report(joinpath(dir, "$name.md"), a;
                        csv_paths = ["f.csv"]), String)
                # The undeterminable state gets its own block, and the unqualified
                # assurance is withheld: "no contrary signal" is not evidence.
                Test.@test occursin("PROVENANCE UNDETERMINABLE", text)
                Test.@test occursin("cannot be determined from this table", text)
                Test.@test !occursin("enforced it — no override was used", text)
                Test.@test occursin(
                    "does NOT establish that the aggregate spans a single MEASURED definition",
                    text)
                p = paired_analysis_json(a; csv_paths = ["f.csv"])
                Test.@test p["metric_definition_provenance_undeterminable"] == true
                Test.@test p["undeterminable_definition_provenance_axes"] ==
                           a.undeterminable_definition_axes
                # NOT overloaded onto the asserted boolean.
                Test.@test p["metric_definition_operator_asserted"] == false
            end
            t_obs = read(
                write_paired_report(joinpath(dir, "obs.md"), a_obs;
                    csv_paths = ["f.csv"]), String)
            Test.@test occursin("enforced it — no override was used", t_obs)
            Test.@test !occursin("PROVENANCE UNDETERMINABLE", t_obs)
            p_obs = paired_analysis_json(a_obs; csv_paths = ["f.csv"])
            Test.@test p_obs["metric_definition_provenance_undeterminable"] == false
            Test.@test isempty(p_obs["undeterminable_definition_provenance_axes"])
        end
    end

    Test.@testset "B: the reader detects what the backfill WRITES (round-trip)" begin
        # A source grep certifies a coupling it cannot detect breaking. Leaving
        # `DEFINITION_PROVENANCE_ASSERTED` defined verbatim while changing only the
        # call site in `_mark_definition_provenance!` to emit another label keeps
        # both grep assertions passing — Julia does not warn on a defined-but-unused
        # const — while this reader stops detecting the writer's output entirely.
        # So the coupling is pinned BEHAVIOURALLY: the writer produces a frame, the
        # reader classifies it.
        #
        # Loaded into its own module because both scripts define `main`.
        df = DataFrames.DataFrame(
            arm = ["naive", "iterative"],
            metric_source = ["quast", missing])
        _PWTBackfill.insert_definition_column!(df, :metric_source, "quast")
        asserted, undeterminable, n_asserted = definition_provenance_axes(df)
        # The row the backfill FILLED is detected as an assertion...
        Test.@test asserted == ["metric_source"]
        Test.@test n_asserted == 1
        # ...and the row it left alone is detected as affirmatively OBSERVED, which
        # pins the other half of the contract: the reader must recognise the exact
        # label the writer uses for untouched rows, or every backfilled table also
        # reads as undeterminable.
        Test.@test isempty(undeterminable)

        # The writer's own constants must round-trip through the reader's two
        # recognised categories, so neither end can be renamed unilaterally.
        Test.@test occursin(RGV_ASSERTED_DEFINITION_MARKER,
            _PWTBackfill.DEFINITION_PROVENANCE_ASSERTED)
        Test.@test _PWTBackfill.DEFINITION_PROVENANCE_OBSERVED ==
                   RGV_OBSERVED_DEFINITION_LABEL
        # The column-naming convention the reader looks under, from the writer's own
        # function rather than a grep of its source.
        Test.@test _PWTBackfill.definition_provenance_column(:metric_source) ==
                   :metric_source_provenance
    end

    Test.@testset "B: the asserted marker matches what the backfill actually writes" begin
        # Retained as a cheap literal check ONLY. It cannot detect the call site
        # drifting away from the const it declares — that is what the round-trip
        # testset above is for.
        backfill = read(
            joinpath(@__DIR__, "..", "..", "benchmarking", "rgv_seed_backfill.jl"),
            String)
        Test.@test occursin(
            "DEFINITION_PROVENANCE_ASSERTED = \"$(RGV_ASSERTED_DEFINITION_MARKER)",
            backfill)
        # And the column-naming convention the reader looks under.
        Test.@test occursin("String(name) * \"_provenance\"", backfill)
    end

    Test.@testset "B: the operator-asserted qualifier does not depend on --allow-mixed-src" begin
        # This was the untested cell of the cross-product. A NON-BINDING
        # `--allow-mixed-src` used to take the `elseif` arm and suppress the
        # operator-asserted qualifier entirely, restoring an unqualified assurance
        # over a table whose definition was never observed. The two facts are
        # orthogonal, so the qualifier must appear in BOTH columns.
        cells = [(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                 for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]]
        asserted = _pwt_frame(cells; provenance = "observed")
        asserted[!, :metric_source_provenance] = fill("operator-asserted (backfill)",
            DataFrames.nrow(asserted))

        for allow in (false, true)
            a = run_paired_analysis(asserted;
                metrics = (:quast_nga50,), allow_mixed_src = allow)
            Test.@test a.definition_operator_asserted
            Test.@test !a.override_bound          # the flag is non-binding either way
            mktempdir() do dir
                text = read(
                    write_paired_report(joinpath(dir, "r.md"), a;
                        csv_paths = ["f.csv"]), String)
                Test.@test occursin("OPERATOR-ASSERTED, NOT OBSERVED", text)
                Test.@test occursin(
                    "does NOT establish that the aggregate spans a single MEASURED definition",
                    text)
                Test.@test !occursin("enforced it — no override was used", text)
                # The specific phrase the non-binding branch used to restore. It is
                # a claim about the VALUES and is not settled by whether a flag
                # bound, so the report must never make it.
                Test.@test !occursin("well-defined", text)
                # And the non-binding flag still reports itself accurately.
                Test.@test occursin("did NOT bind", text) == allow
            end
        end
    end

    Test.@testset "C1/I3: the report does not claim to be a pre-registered test" begin
        df = _pwt_frame([(seed = s, error_rate = e, naive = 1000.0, iterative = 1600.0)
                         for (s, e) in [(42, 0.01), (42, 0.05), (123, 0.01),
            (123, 0.05), (456, 0.01), (456, 0.05)]])
        a = run_paired_analysis(df; metrics = (:quast_nga50,))
        mktempdir() do dir
            text = read(
                write_paired_report(joinpath(dir, "r.md"), a;
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
            # Every outcome decide() can emit must be defined in the legend, or a
            # reader meets a verdict string the report never explains.
            for outcome in ("SUPPORTED", "PARTIALLY SUPPORTED", "FALSIFIED",
                "NOT SUPPORTED", "INDETERMINATE")
                Test.@test occursin("**$outcome**", text)
            end
        end
    end
end
