module Td2rxhEvidenceAnalysisTest

# Contract tests for the td-2rxh evidence path: the `--cells td2rxh` grid and CSV
# schema of benchmarking/indel_reroute_evidence.jl, and the pre-registered
# decision rule in benchmarking/td2rxh_paired_analysis.jl. No assembly runs here;
# the harness runs on HPC. The decision-rule testsets include a POSITIVE control
# (a clean win must be accepted) and NEGATIVE controls (a loss, a fit-identity
# regression and a missing read length must each be rejected), so that a rule
# which always says yes, or always says no, fails.

import Test
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

"""
Synthetic two-arm table in the harness's CSV schema. `delta(cell_index, seed)`
returns the nanopore-minus-illumina coverage delta; deltas are distinct across
seeds so the exact Wilcoxon branch is exercised.
"""
function synthetic_table(; cells = Harness.TD2RXH_CELLS, seeds = 1:8,
        delta = (c, s) -> 0.01 * s, fit_delta = (c, s) -> 0.0,
        extra_seconds_per_read = 0.05, rng = StableRNGs.StableRNG(7))
    rows = NamedTuple[]
    for (c, cell) in enumerate(cells), seed in seeds

        base_cov = 0.3 + 0.1 * Random.rand(rng)
        n_reads = ceil(Int, cell.coverage * cell.genome_length / cell.read_length)
        for arm in ("illumina", "nanopore")
            is_ont = arm == "nanopore"
            push!(rows,
                (cell = Harness.cell_id(cell), seed = seed, arm = arm,
                    n_reads = n_reads, wall_seconds = 10.0 +
                                                      (is_ont ?
                                                       extra_seconds_per_read * n_reads :
                                                       0.0),
                    best_contig_reference_coverage = base_cov +
                                                     (is_ont ? delta(c, seed) : 0.0),
                    best_contig_fit_identity = 0.99 + (is_ont ? fit_delta(c, seed) : 0.0),
                    n_contigs = 100, largest_contig = 1000,
                    substitution_length_divergences = 0, window_anchor_rejections = 0,
                    indel_engaged = is_ont ? 5 : 0))
        end
    end
    return DataFrames.DataFrame(rows)
end

function analyse(df)
    paired, dropped = Analysis.pair_arms(df)
    verdicts = Analysis.cell_verdicts(paired)
    (paired = paired, dropped = dropped, verdicts = verdicts,
        verdict = Analysis.acceptance(verdicts))
end

Test.@testset "td-2rxh evidence harness and pre-registered analysis" begin
    Test.@testset "td2rxh cell grid matches the pre-registration" begin
        cells = Harness.select_cells("td2rxh")
        Test.@test length(cells) == 6
        Test.@test all(c -> c.genome_length == 6000 && c.coverage == 20 && c.kmax == 31, cells)
        Test.@test sort(unique(c.read_length for c in cells)) == [1000, 2000]
        Test.@test sort(unique(c.error_rate for c in cells)) == [0.01, 0.02, 0.05]
        Test.@test allunique(Harness.cell_id.(cells))
        Test.@test Harness.select_cells("all") === Harness.ALL_CELLS
        Test.@test Harness.select_cells("smoke") === Harness.SMOKE_CELLS
        Test.@test_throws ArgumentError Harness.select_cells("td2rhx")
    end

    Test.@testset "CSV schema is prefix-compatible and row width matches header" begin
        Test.@test Harness.CSV_HEADER[1:length(LEGACY_HEADER)] == LEGACY_HEADER
        result = (substitution_length_divergences = 1, window_divergences = 2,
            window_anchor_rejections = 3, trace_contract_errors = 0,
            indel_engaged = 4, indel_requested = 5, reassembly_k = 31,
            n_contigs = 6, largest_contig = 700, largest_ratio = 0.1167,
            total_contig_bases = 5000, wall_seconds = 12.3456,
            wall_seconds_per_read = 0.1029, best_contig_reference_coverage = 0.1167,
            best_contig_fit_identity = 0.99812)
        row = Harness.csv_row(
            "label", first(Harness.TD2RXH_CELLS), 3, :nanopore, 120, result)
        fields = split(row, ",")
        Test.@test length(fields) == length(Harness.CSV_HEADER)
        Test.@test fields[findfirst(==("wall_seconds"), Harness.CSV_HEADER)] == "12.346"
        Test.@test Harness.per_read_seconds(12.0, 120) ≈ 0.1
        Test.@test isnan(Harness.per_read_seconds(12.0, 0))
    end

    Test.@testset "holm_adjust" begin
        Test.@test Analysis.holm_adjust([0.01, 0.04, 0.03]) ≈ [0.03, 0.06, 0.06]
        adjusted = Analysis.holm_adjust([0.01, NaN, 0.03])
        Test.@test adjusted[1] ≈ 0.02 && adjusted[3] ≈ 0.03 && isnan(adjusted[2])
        Test.@test Analysis.holm_adjust([0.9, 0.8]) == [1.0, 1.0]
    end

    Test.@testset "parse_cell_id round-trips the harness id" begin
        for cell in Harness.TD2RXH_CELLS
            parsed = Analysis.parse_cell_id(Harness.cell_id(cell))
            Test.@test parsed.read_length == cell.read_length
            Test.@test parsed.error_rate == cell.error_rate
        end
        Test.@test_throws ArgumentError Analysis.parse_cell_id("6000/1000/20x")
    end

    Test.@testset "POSITIVE control: a consistent win is accepted" begin
        r = analyse(synthetic_table())
        Test.@test DataFrames.nrow(r.paired) == 48 && DataFrames.nrow(r.dropped) == 0
        Test.@test all(r.verdicts.wilcoxon_method .== "exact (null rank-sum distribution)")
        # n=8 all-positive exact two-sided p = 2/256; Holm x6 = 0.046875 < 0.05.
        Test.@test all(r.verdicts.p_raw .≈ 2 / 256)
        Test.@test all(r.verdicts.h2_holds)
        Test.@test all(r.verdicts.median_incremental_s_per_read .≈ 0.05)
        Test.@test r.verdict.h1_holds_all_cells
        Test.@test r.verdict.accepted
        Test.@test r.verdict.qualifying_error_rates == [0.01, 0.02, 0.05]
    end

    Test.@testset "NEGATIVE control: a consistent loss is rejected" begin
        r = analyse(synthetic_table(delta = (c, s) -> -0.01 * s))
        Test.@test !any(r.verdicts.h2_holds)
        Test.@test !r.verdict.accepted
    end

    Test.@testset "NEGATIVE control: win at only one read length is rejected" begin
        loses_at_2000 = (
            c, s) -> Harness.TD2RXH_CELLS[c].read_length == 2000 ?
                     -0.01 * s : 0.01 * s
        r = analyse(synthetic_table(delta = loses_at_2000))
        Test.@test count(r.verdicts.h2_holds) == 3
        Test.@test !r.verdict.accepted
    end

    Test.@testset "NEGATIVE control: contiguity win with fit-identity loss is rejected" begin
        r = analyse(synthetic_table(fit_delta = (c, s) -> -0.001))
        Test.@test all(r.verdicts.h2_holds)
        Test.@test !r.verdict.accepted
    end

    Test.@testset "H1' threshold and a missing arm" begin
        slow = analyse(synthetic_table(extra_seconds_per_read = 0.5))
        Test.@test !slow.verdict.h1_holds_all_cells
        df = synthetic_table()
        drop = findfirst(r -> r.seed == 3 && r.arm == "nanopore", eachrow(df))
        r = analyse(df[setdiff(1:DataFrames.nrow(df), [drop]), :])
        Test.@test DataFrames.nrow(r.dropped) == 1
        Test.@test DataFrames.nrow(r.paired) == 47
    end

    Test.@testset "NEGATIVE control: a missing read-length cell cannot qualify" begin
        df = synthetic_table(cells = filter(c -> c.read_length == 1000, Harness.TD2RXH_CELLS))
        r = analyse(df)
        Test.@test all(r.verdicts.h2_holds)
        Test.@test !r.verdict.accepted
    end
end

end
