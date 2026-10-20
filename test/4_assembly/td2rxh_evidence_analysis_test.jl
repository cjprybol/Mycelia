module Td2rxhEvidenceAnalysisTest

# Contract tests for the td-2rxh evidence path: the `--cells td2rxh` grid and CSV
# schema of benchmarking/indel_reroute_evidence.jl, and the pre-registered
# decision rule (with amendment 1) in benchmarking/td2rxh_paired_analysis.jl. No
# assembly runs here; the harness runs on HPC.
#
# The decision-rule testsets pair a POSITIVE control (a clean win must be
# accepted) with NEGATIVE controls (a loss, a one-read-length win, a fit-identity
# regression, a raw-significant-but-not-after-Holm cell, a missing cell, an
# undefined cell, an unengaged treatment, duplicate and mixed-label shards), so a
# rule that always says yes, or always says no, fails. Several controls exist
# because a review-round mutant survived without them; each says which.

import Test
import CSV
import DataFrames
import StableRNGs
import Random

module Harness
include(joinpath(@__DIR__, "..", "..", "benchmarking", "indel_reroute_evidence.jl"))
end

module Analysis
include(joinpath(@__DIR__, "..", "..", "benchmarking", "td2rxh_paired_analysis.jl"))
end

const LEGACY_HEADER = [
    "label", "cell", "seed", "arm", "n_reads",
    "substitution_length_divergences", "window_divergences",
    "window_anchor_rejections", "trace_contract_errors",
    "indel_engaged", "indel_requested", "reassembly_k",
    "n_contigs", "largest_contig", "largest_ratio",
    "total_contig_bases", "wall_seconds"]
const NEW_COLUMNS = [
    "wall_seconds_per_read", "best_contig_reference_coverage", "best_contig_fit_identity"]
const VALID_LABEL = "td2rxh@0123456789ab+manifestabcdef012345"
const PREREGISTERED_CELLS = [
    "6000/1000/20x/e0.01", "6000/1000/20x/e0.02", "6000/1000/20x/e0.05",
    "6000/2000/20x/e0.01", "6000/2000/20x/e0.02", "6000/2000/20x/e0.05"]

"""
Synthetic two-arm table in the analysis input schema. `delta(c, s)` and
`fit_delta(c, s)` give nanopore-minus-illumina per (cell index, seed). Default
deltas are distinct and positive across seeds.
"""
function synthetic_table(; cells = PREREGISTERED_CELLS, seeds = 1:16,
        delta = (c, s) -> 0.001 * s, fit_delta = (c, s) -> 0.0,
        extra_seconds_per_read = 0.05, engaged = (c, s) -> 5, label = VALID_LABEL,
        rng = StableRNGs.StableRNG(7))
    rows = NamedTuple[]
    for (c, cell) in enumerate(cells), seed in seeds

        base_cov = round(0.3 + 0.1 * Random.rand(rng); digits = 6)
        read_length = Analysis.parse_cell_id(cell).read_length
        n_reads = ceil(Int, 20 * 6000 / read_length)
        for arm in ("illumina", "nanopore")
            is_ont = arm == "nanopore"
            push!(rows,
                (label = label, cell = cell, seed = seed, arm = arm,
                    n_reads = n_reads,
                    wall_seconds = 10.0 + (is_ont ? extra_seconds_per_read * n_reads : 0.0),
                    best_contig_reference_coverage = base_cov +
                                                     (is_ont ? delta(c, seed) : 0.0),
                    best_contig_fit_identity = 0.99 + (is_ont ? fit_delta(c, seed) : 0.0),
                    n_contigs = 100, largest_contig = 1000,
                    substitution_length_divergences = 0, window_anchor_rejections = 0,
                    indel_engaged = is_ont ? engaged(c, seed) : 0))
        end
    end
    return DataFrames.DataFrame(rows)
end

function analyse(df; n_seeds = 16)
    Analysis.validate_shards(df)
    paired, unpaired = Analysis.pair_arms(df; n_seeds = n_seeds)
    verdicts = Analysis.cell_verdicts(paired; n_seeds = n_seeds)
    return (paired = paired, unpaired = unpaired, verdicts = verdicts,
        h1 = Analysis.h1_verdict(verdicts), verdict = Analysis.acceptance(verdicts))
end

cell_index(cell) = findfirst(==(cell), PREREGISTERED_CELLS)
rate_of(c) = Analysis.parse_cell_id(PREREGISTERED_CELLS[c]).error_rate
readlen_of(c) = Analysis.parse_cell_id(PREREGISTERED_CELLS[c]).read_length

Test.@testset "td-2rxh evidence harness and pre-registered analysis" begin
    Test.@testset "td2rxh cell grid matches the pre-registration" begin
        cells = Harness.select_cells("td2rxh")
        Test.@test Harness.cell_id.(cells) == PREREGISTERED_CELLS
        Test.@test Analysis.TD2RXH_EXPECTED_CELLS == PREREGISTERED_CELLS
        Test.@test all(c -> c.kmax == 31, cells)
        Test.@test Analysis.TD2RXH_N_SEEDS == 16
        Test.@test Analysis.TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ == 0.30
        Test.@test Analysis.TD2RXH_ALPHA == 0.05
        Test.@test Harness.select_cells("all") === Harness.ALL_CELLS
        Test.@test Harness.select_cells("smoke") === Harness.SMOKE_CELLS
        Test.@test_throws ArgumentError Harness.select_cells("td2rhx")
    end

    Test.@testset "harness argument parsing is strict" begin
        Test.@test Harness.parse_args(["--cells", "td2rxh", "--seeds", "3"]) ==
                   Dict("cells" => "td2rxh", "seeds" => "3")
        Test.@test_throws ArgumentError Harness.parse_args(["--cell", "td2rxh"])
        Test.@test_throws ArgumentError Harness.parse_args(["--cells"])
        Test.@test_throws ArgumentError Harness.parse_args(["--cells", "--seeds"])
        Test.@test Harness.parse_bool_flag("warmup", "false") == false
        Test.@test_throws ArgumentError Harness.parse_bool_flag("warmup", "TRUE")
        Test.@test Harness.arm_order(1) == (:illumina, :nanopore)
        Test.@test Harness.arm_order(2) == (:nanopore, :illumina)
    end

    Test.@testset "CSV schema: prefix-compatible, new names and positions pinned" begin
        Test.@test Harness.CSV_HEADER == vcat(LEGACY_HEADER, NEW_COLUMNS)
        result = (substitution_length_divergences = 1, window_divergences = 2,
            window_anchor_rejections = 3, trace_contract_errors = 0,
            indel_engaged = missing, indel_requested = 5, reassembly_k = 31,
            n_contigs = 6, largest_contig = 700, largest_ratio = 0.1167,
            total_contig_bases = 5000, wall_seconds = 12.3456,
            wall_seconds_per_read = 0.1029, best_contig_reference_coverage = 0.1166667,
            best_contig_fit_identity = 0.9981234)
        fields = split(
            Harness.csv_row("label", first(Harness.TD2RXH_CELLS), 3,
                :nanopore, 120, result), ",")
        column(name) = fields[findfirst(==(name), Harness.CSV_HEADER)]
        Test.@test length(fields) == length(Harness.CSV_HEADER)
        # Distinct values per column: a swap of any two (review mutant C) fails.
        Test.@test column("wall_seconds") == "12.346"
        Test.@test column("wall_seconds_per_read") == "0.1029"
        Test.@test column("best_contig_reference_coverage") == "0.116667"
        Test.@test column("best_contig_fit_identity") == "0.998123"
        Test.@test column("indel_engaged") == ""
        Test.@test Harness.per_read_seconds(12.0, 120) ≈ 0.1
        Test.@test isnan(Harness.per_read_seconds(12.0, 0))
    end

    Test.@testset "exact sign-flip signed-rank test" begin
        # All 8 positive, distinct: 2 of 256 sign assignments are as extreme.
        Test.@test Analysis.signflip_signed_rank_test(collect(1.0:8.0)).pvalue ≈ 2 / 256
        # All tied: still exact (the asymptotic code path would leave the exact branch).
        tied = Analysis.signflip_signed_rank_test(fill(0.5, 8))
        Test.@test tied.pvalue ≈ 2 / 256
        Test.@test startswith(tied.method, "exact sign-flip")
        # Float noise below 1e-6 must not create distinct ranks or extra zeros.
        noisy = Analysis.signflip_signed_rank_test(
            [0.3333, 0.33329999999999993, 0.33330000000000004, 0.3334])
        Test.@test noisy.n == 4 && noisy.pvalue ≈ 2 / 16
        # The rounding is load-bearing: a float residue must become an exact zero.
        # Without it 0.3 - (0.1 + 0.2) = -5.6e-17 would count as a non-zero delta.
        residue = Analysis.signflip_signed_rank_test([0.3 - (0.1 + 0.2), 1.0, 2.0])
        Test.@test residue.n == 2 && residue.n_zero_dropped == 1
        zeros_dropped = Analysis.signflip_signed_rank_test([0.0, 0.0, 1.0, 2.0, 3.0])
        Test.@test zeros_dropped.n == 3 && zeros_dropped.n_zero_dropped == 2
        Test.@test zeros_dropped.pvalue ≈ 2 / 8
        undefined = Analysis.signflip_signed_rank_test(zeros(4))
        Test.@test undefined.n == 0 && isnan(undefined.pvalue)
        Test.@test_throws ArgumentError Analysis.signflip_signed_rank_test([1.0, NaN])
        # Symmetric: sign reversal leaves the two-sided p unchanged.
        Test.@test Analysis.signflip_signed_rank_test(-collect(1.0:8.0)).pvalue ≈ 2 / 256
    end

    Test.@testset "holm_adjust uses a fixed family; NaN counts as p = 1 (review mutant E)" begin
        Test.@test Analysis.holm_adjust([0.01, 0.04, 0.03]) ≈ [0.03, 0.06, 0.06]
        # NaN stays IN the family (m = 3), unlike the round-1 code (m = 2).
        Test.@test Analysis.holm_adjust([0.01, NaN, 0.03]) ≈ [0.03, 1.0, 0.06]
        Test.@test Analysis.holm_adjust([0.9, 0.8]) == [1.0, 1.0]
    end

    Test.@testset "POSITIVE control: a consistent win is accepted" begin
        r = analyse(synthetic_table())
        Test.@test DataFrames.nrow(r.paired) == 96 && DataFrames.nrow(r.unpaired) == 0
        Test.@test all(r.verdicts.status .== "COMPLETE")
        Test.@test all(r.verdicts.p_raw .≈ 2 / 2^16)
        Test.@test all(r.verdicts.p_holm .≈ 6 * 2 / 2^16)
        Test.@test all(r.verdicts.h2_holds)
        Test.@test all(r.verdicts.median_incremental_s_per_read .≈ 0.05)
        Test.@test r.h1 == "HOLDS"
        Test.@test r.verdict.verdict == "MET"
        Test.@test r.verdict.qualifying_error_rates == [0.01, 0.02, 0.05]
    end

    Test.@testset "NEGATIVE control: a consistent loss is rejected as NOT MET" begin
        r = analyse(synthetic_table(delta = (c, s) -> -0.001 * s))
        Test.@test !any(r.verdicts.h2_holds)
        Test.@test r.verdict.verdict == "NOT MET"
    end

    Test.@testset "at least ONE error rate suffices (review mutant B)" begin
        win_at_005 = (c, s) -> rate_of(c) == 0.05 ? 0.001 * s : -0.001 * s
        r = analyse(synthetic_table(delta = win_at_005))
        Test.@test r.verdict.qualifying_error_rates == [0.05]
        Test.@test r.verdict.verdict == "MET"
        # Fit identity is judged per qualifying rate, not across all cells.
        fit_loss_at_001 = (c, s) -> rate_of(c) == 0.01 ? -0.001 : 0.0
        r2 = analyse(synthetic_table(fit_delta = fit_loss_at_001))
        Test.@test r2.verdict.qualifying_error_rates == [0.02, 0.05]
    end

    Test.@testset "NEGATIVE control: win at only one read length is rejected" begin
        loses_at_2000 = (c, s) -> readlen_of(c) == 2000 ? -0.001 * s : 0.001 * s
        r = analyse(synthetic_table(delta = loses_at_2000))
        Test.@test count(r.verdicts.h2_holds) == 3
        Test.@test r.verdict.verdict == "NOT MET"
    end

    Test.@testset "NEGATIVE control: contiguity win with fit-identity loss" begin
        r = analyse(synthetic_table(fit_delta = (c, s) -> -0.001))
        Test.@test all(r.verdicts.h2_holds)
        Test.@test r.verdict.verdict == "NOT MET"
    end

    Test.@testset "Holm decides, not raw p (review mutant A)" begin
        # Negatives at the 7 smallest ranks (W- = 28 of 136): raw p is under alpha
        # but at least x2 under Holm pushes it over, so nothing may qualify.
        mixed = (c, s) -> s <= 7 ? -0.001 * s : 0.001 * s
        p = Analysis.signflip_signed_rank_test([mixed(1, s) for s in 1:16]).pvalue
        Test.@test p < 0.05 && 6 * p > 0.05 && 2 * p > 0.05
        r = analyse(synthetic_table(delta = mixed))
        Test.@test all(r.verdicts.p_raw .< 0.05)
        Test.@test !any(r.verdicts.h2_holds)
        Test.@test r.verdict.verdict == "NOT MET"
    end

    Test.@testset "a non-COMPLETE cell enters Holm at p = 1 even when its p is tiny" begin
        # Pick a negative-rank pattern whose exact p lies in (0.05/6, 0.05/5): it
        # passes Holm at multiplier 5 but not at 6. Then give cell 1 a tiny p but
        # no indel engagement. Entering the family at p = 1 (correct) keeps the
        # other five at multiplier 6, so nothing qualifies; letting its tiny p into
        # the family (review-round mutant F) would drop them to x5 and accept.
        function negatives_for(target)
            chosen, remaining = Int[], target
            for r in 16:-1:1
                r <= remaining && (push!(chosen, r); remaining -= r)
            end
            return chosen
        end
        pick = nothing
        for w in 1:60
            neg = negatives_for(w)
            p = Analysis.signflip_signed_rank_test([s in neg ? -s : s for s in 1:16]).pvalue
            if 5 * p < 0.05 < 6 * p
                pick = neg
                break
            end
        end
        Test.@test pick !== nothing
        pattern = (c, s) -> c == 1 ? 0.001 * s : (s in pick ? -0.001 * s : 0.001 * s)
        r = analyse(synthetic_table(delta = pattern, engaged = (c, s) -> c == 1 ? 0 : 5))
        Test.@test r.verdicts[1, :status] == "TREATMENT-NOT-APPLIED"
        Test.@test r.verdicts[1, :p_holm] == 1.0
        Test.@test all(r.verdicts[2:6, :p_holm] .> 0.05)
        Test.@test isempty(r.verdict.qualifying_error_rates)
    end

    Test.@testset "a missing cell stays in the Holm family and blocks NOT MET" begin
        present = filter(!=("6000/2000/20x/e0.05"), PREREGISTERED_CELLS)
        r = analyse(synthetic_table(cells = present))
        missing_row = r.verdicts[r.verdicts.cell .== "6000/2000/20x/e0.05", :]
        Test.@test DataFrames.nrow(r.verdicts) == 6
        Test.@test only(missing_row.status) == "MISSING"
        Test.@test only(missing_row.p_holm) == 1.0
        Test.@test count(==("no rows"), r.unpaired.reason) == 16
        # Family is still 6: the best cell is multiplied by 6, not 5.
        Test.@test minimum(r.verdicts.p_holm) ≈ 6 * 2 / 2^16
        Test.@test r.verdict.qualifying_error_rates == [0.01, 0.02]
        Test.@test r.h1 == "INCOMPLETE"
        loss = analyse(synthetic_table(cells = present, delta = (c, s) -> -0.001 * s))
        Test.@test loss.verdict.verdict == "INCOMPLETE"
    end

    Test.@testset "an incomplete cell (one seed's arm threw) cannot qualify" begin
        df = synthetic_table()
        drop = findfirst(
            r -> r.cell == "6000/1000/20x/e0.05" && r.seed == 3 &&
                 r.arm == "nanopore", eachrow(df))
        r = analyse(df[setdiff(1:DataFrames.nrow(df), [drop]), :])
        row = r.verdicts[r.verdicts.cell .== "6000/1000/20x/e0.05", :]
        Test.@test only(row.status) == "INCOMPLETE" && only(row.n_pairs) == 15
        Test.@test only(r.unpaired.reason) == "illumina only"
        Test.@test !only(row.h2_holds)
        Test.@test r.verdict.qualifying_error_rates == [0.01, 0.02]
    end

    Test.@testset "an all-zero cell is UNDEFINED, enters Holm at 1" begin
        r = analyse(synthetic_table(delta = (c, s) -> c == 1 ? 0.0 : 0.001 * s))
        row = r.verdicts[1, :]
        Test.@test row.status == "UNDEFINED" && row.p_holm == 1.0
        Test.@test r.verdict.qualifying_error_rates == [0.02, 0.05]
        # Amendment 2: an UNDEFINED cell is a measured null, so with no qualifying
        # rate and nothing missing the verdict is final.
        loss = analyse(synthetic_table(delta = (c, s) -> c == 1 ? 0.0 : -0.001 * s))
        Test.@test loss.verdict.verdict == "NOT MET"
        # ...and its timing still counts for H1' (it can refute).
        slow = analyse(synthetic_table(delta = (c, s) -> 0.0, extra_seconds_per_read = 2.5))
        Test.@test all(slow.verdicts.status .== "UNDEFINED")
        Test.@test all(slow.verdicts.h1_cell .== "exceeds")
        Test.@test slow.h1 == "REFUTED"
        Test.@test slow.verdict.verdict == "NOT MET"
        fast = analyse(synthetic_table(delta = (c, s) -> 0.0, extra_seconds_per_read = 0.05))
        Test.@test fast.h1 == "HOLDS"
    end

    Test.@testset "a TREATMENT-NOT-APPLIED cell keeps the verdict INCOMPLETE" begin
        loss = analyse(synthetic_table(delta = (c, s) -> -0.001 * s,
            engaged = (c, s) -> c == 1 ? 0 : 5))
        Test.@test loss.verdict.verdict == "INCOMPLETE"
        Test.@test loss.verdicts[1, :h1_cell] == "n/a"
        Test.@test loss.h1 == "INCOMPLETE"
    end

    Test.@testset "treatment never engaged -> TREATMENT-NOT-APPLIED" begin
        r = analyse(synthetic_table(engaged = (c, s) -> c == 1 ? 0 : 5))
        Test.@test r.verdicts[1, :status] == "TREATMENT-NOT-APPLIED"
        Test.@test !r.verdicts[1, :h2_holds]
        Test.@test Analysis._engaged("3") && !Analysis._engaged("") &&
                   !Analysis._engaged(missing) && !Analysis._engaged(0)
    end

    Test.@testset "empty assembly: NaN fit identity is excluded, not scored 0" begin
        df = synthetic_table()
        ill = (df.arm .== "illumina") .& (df.cell .== "6000/1000/20x/e0.01")
        df.best_contig_fit_identity = Float64.(df.best_contig_fit_identity)
        df[ill, :best_contig_fit_identity] .= NaN
        r = analyse(df)
        Test.@test r.verdicts[1, :n_fit_pairs] == 0
        Test.@test isnan(r.verdicts[1, :median_fit_identity_delta])
        Test.@test r.verdict.qualifying_error_rates == [0.02, 0.05]
        # One NaN pair must be EXCLUDED, not propagated: the finite median over the
        # other 15 pairs is 0, so the rate still qualifies.
        one = synthetic_table()
        one.best_contig_fit_identity = Float64.(one.best_contig_fit_identity)
        one[findfirst(ill .& (df.seed .== 1)), :best_contig_fit_identity] = NaN
        r1 = analyse(one)
        Test.@test r1.verdicts[1, :n_fit_pairs] == 15
        Test.@test r1.verdicts[1, :median_fit_identity_delta] == 0.0
        Test.@test r1.verdict.qualifying_error_rates == [0.01, 0.02, 0.05]
    end

    Test.@testset "H1' threshold and three states" begin
        Test.@test analyse(synthetic_table(extra_seconds_per_read = 0.29)).h1 == "HOLDS"
        Test.@test analyse(synthetic_table(extra_seconds_per_read = 0.31)).h1 == "REFUTED"
        # The threshold is inclusive: exactly 0.30 s/read holds.
        Test.@test analyse(synthetic_table(extra_seconds_per_read = 0.30)).h1 == "HOLDS"
        # Arms disagreeing on read count cannot be scored.
        df = synthetic_table()
        df[1, :n_reads] += 1
        Test.@test_throws ArgumentError analyse(df)
    end

    Test.@testset "shard integrity errors (review mutant D)" begin
        df = synthetic_table()
        Test.@test_throws ArgumentError Analysis.validate_shards(vcat(df, df[1:1, :]))
        mixed = copy(df)
        mixed[1:2, :label] .= "td2rxh@ba9876543210+manifestabcdef012345"
        Test.@test_throws ArgumentError Analysis.validate_shards(mixed)
        extra = copy(df)
        extra[1:2, :cell] .= "2000/1200/8x/e0.01"
        Test.@test_throws ArgumentError Analysis.validate_shards(extra)
        Test.@test_throws ArgumentError Analysis.validate_shards(
            DataFrames.select(df, DataFrames.Not(:indel_engaged)))
        stray_seed = copy(df)
        stray_seed[1:2, :seed] .= 17
        Test.@test_throws ArgumentError Analysis.validate_shards(stray_seed)
        bad_label = copy(df)
        bad_label.label .= "td2rxh@"
        Test.@test_throws ArgumentError Analysis.validate_shards(bad_label)
        no_label = copy(df)
        no_label.label = Vector{Union{Missing, String}}(missing, DataFrames.nrow(df))
        Test.@test_throws ArgumentError Analysis.validate_shards(no_label)
        Test.@test Analysis.validate_shards(df) === nothing
    end

    Test.@testset "analysis arguments are strict" begin
        Test.@test_throws ArgumentError Analysis.parse_analysis_args(
            ["--csv", "a.csv", "--output-dir"])
        Test.@test_throws ArgumentError Analysis.parse_analysis_args(
            ["--csv", "a.csv", "--output_dir", "x"])
        Test.@test_throws ArgumentError Analysis.parse_analysis_args(["--output-dir", "x"])
        Test.@test_throws ArgumentError Analysis.parse_analysis_args([
            "--csv", "--output-dir"])
    end

    Test.@testset "end to end: harness csv_row shards -> analysis main" begin
        dir = mktempdir()
        df = synthetic_table()
        # One empty illumina assembly: its NaN fit identity must survive the
        # csv_row -> CSV.read round trip and drop out of that cell's fit median (15 pairs), not score 0.
        empty_arm = findfirst(
            r -> r.cell == "6000/1000/20x/e0.01" && r.seed == 1 &&
                 r.arm == "illumina", eachrow(df))
        df.best_contig_fit_identity = Float64.(df.best_contig_fit_identity)
        df[empty_arm, :best_contig_fit_identity] = NaN
        paths = String[]
        for seed in 1:16
            path = joinpath(dir, "seed$(seed).csv")
            open(path, "w") do io
                println(io, join(Harness.CSV_HEADER, ","))
                for r in eachrow(df[df.seed .== seed, :])
                    cell = Harness.TD2RXH_CELLS[cell_index(r.cell)]
                    result = (substitution_length_divergences = 0, window_divergences = 0,
                        window_anchor_rejections = 0, trace_contract_errors = 0,
                        indel_engaged = r.indel_engaged, indel_requested = 1,
                        reassembly_k = 31, n_contigs = r.n_contigs,
                        largest_contig = r.largest_contig, largest_ratio = 0.1,
                        total_contig_bases = 1000, wall_seconds = r.wall_seconds,
                        wall_seconds_per_read = r.wall_seconds / r.n_reads,
                        best_contig_reference_coverage = r.best_contig_reference_coverage,
                        best_contig_fit_identity = r.best_contig_fit_identity)
                    println(io, Harness.csv_row(r.label, cell, r.seed,
                        Symbol(r.arm), r.n_reads, result))
                end
            end
            push!(paths, path)
        end
        out = joinpath(dir, "analysis")
        Test.@test Analysis.td2rxh_analysis_main(
            ["--csv", join(paths, ","), "--output-dir", out]) == 0
        report = read(joinpath(out, "report.md"), String)
        Test.@test occursin("td-2rxh ACCEPTANCE: MET", report)
        Test.@test occursin("H1': HOLDS", report)
        Test.@test occursin("Qualifying error rates: [0.01, 0.02, 0.05]", report)
        Test.@test isfile(joinpath(out, "td2rxh_paired_deltas.svg"))
        Test.@test isfile(joinpath(out, "td2rxh_paired_deltas.png"))
        verdict_table = CSV.read(joinpath(out, "td2rxh_cell_verdicts.csv"),
            DataFrames.DataFrame)
        Test.@test DataFrames.nrow(verdict_table) == 6
        Test.@test verdict_table[1, :n_fit_pairs] == 15
        # A duplicated shard is refused through the real entry point, and the
        # previous generation is gone, so no stale "MET" report survives the failure.
        Test.@test_throws ArgumentError Analysis.td2rxh_analysis_main(
            ["--csv", join(vcat(paths, paths[1:1]), ","), "--output-dir", out])
        Test.@test !isfile(joinpath(out, "report.md"))
    end
end

end
