module Td4mbgBaseFixAnalysisTest

# Contract tests for benchmarking/td4mbg_base_fix_analysis.jl, the pre-registered
# BASE vs FIX decision rule of bead td-4mbg (amendments 4 and 5). No assembly runs
# here; the shards and oracle CSVs are synthetic.
#
# Each gate gets a POSITIVE control (a clean FIX that moves the delta toward 0 is
# SUPPORTED) and NEGATIVE controls (a FIX that moves it away is REFUTED; an
# illumina change anywhere, or a nanopore change in a 2000 bp cell, VOIDS; a
# nanopore change in a 1000 bp cell does NOT void; a manifest or head mismatch
# VOIDS; a failed oracle gate or unequal oracle sha VOIDS while a failed advisory
# row does not; a missing shard is INCOMPLETE). The one-sided test is checked
# against hand enumeration.

import Test
import CSV
import DataFrames

module Analysis
include(joinpath(@__DIR__, "..", "..", "benchmarking", "td4mbg_base_fix_analysis.jl"))
end

const CELLS = [
    "6000/1000/20x/e0.01",
    "6000/1000/20x/e0.02",
    "6000/1000/20x/e0.05",
    "6000/2000/20x/e0.01",
    "6000/2000/20x/e0.02",
    "6000/2000/20x/e0.05",
]
const BASE_HEAD = "aaaaaaaaaaaa"
const FIX_HEAD = "bbbbbbbbbbbb"
const MANIFEST = "cccccccccccc"
const SHA_A = repeat("d", 64)
const SHA_B = repeat("e", 64)

"""
Synthetic run in the full shard schema. Illumina coverage is 0.9 everywhere;
nanopore coverage is `0.9 + delta(cell, seed)`. `wall` differs per run so the
identity checks must ignore the wall-time columns. `tweak(row)` may edit a row.
"""
function synthetic_run(;
    head,
    manifest = MANIFEST,
    delta = (c, s) -> -0.01 * s,
    wall = 1.0,
    tweak = r -> r,
)
    rows = NamedTuple[]
    for cell in CELLS, seed = 1:16, arm in ("illumina", "nanopore")
        cov = round(0.9 + (arm == "nanopore" ? delta(cell, seed) : 0.0); digits = 6)
        row = (
            label = "td2rxh@$(head)+manifest$(manifest)",
            cell = cell,
            seed = seed,
            arm = arm,
            n_reads = 120,
            substitution_length_divergences = 0,
            window_divergences = 0,
            window_anchor_rejections = 0,
            trace_contract_errors = 0,
            indel_engaged = arm == "nanopore" ? 3 : 0,
            indel_requested = 0,
            reassembly_k = 31,
            n_contigs = 100,
            largest_contig = 1000,
            largest_ratio = 0.5,
            total_contig_bases = 9000,
            wall_seconds = 10.0 * wall + seed,
            wall_seconds_per_read = 0.1 * wall,
            best_contig_reference_coverage = cov,
            best_contig_fit_identity = 1.0,
        )
        push!(rows, tweak(row))
    end
    return DataFrames.DataFrame(rows)
end

"""
FIX-run nanopore delta: `f(seed)` in the 1000 bp cells (where PR #487 may change
results), the BASE delta `-0.01 * seed` in the 2000 bp cells (where G0b requires
identity).
"""
fix_delta(f) = (c, s) -> startswith(c, "6000/1000/") ? f(s) : -0.01 * s

function write_run(dir, df; skip_seeds = Int[])
    mkpath(dir)
    for seed = 1:16
        seed in skip_seeds && continue
        CSV.write(joinpath(dir, "seed$(seed).csv"), df[df.seed .== seed, :])
    end
    return dir
end

function write_oracle(path; sha = SHA_A, fail_gate = false, fail_advisory = false)
    rows = [
        (
            check = "gate_$(i)",
            severity = "gate",
            passed = !(fail_gate && i == 1),
            detail = "ok",
        ) for i = 1:17
    ]
    push!(
        rows,
        (
            check = "nanopore_under_120_seconds",
            severity = "advisory",
            passed = !fail_advisory,
            detail = "wall",
        ),
    )
    push!(
        rows,
        (
            check = "illumina_byte_identical_to_prewiring_oracle",
            severity = "gate",
            passed = true,
            detail = "sha256=$(sha), origin_master=$(sha)",
        ),
    )
    CSV.write(path, DataFrames.DataFrame(rows))
    return path
end

function run_case(;
    base = synthetic_run(head = BASE_HEAD),
    fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s), wall = 2.0),
    base_sha = SHA_A,
    fix_sha = SHA_A,
    fix_fail_gate = false,
    fix_fail_advisory = false,
    skip_fix_seeds = Int[],
    fix_head = FIX_HEAD,
)
    root = mktempdir()
    return Analysis.td4mbg_analysis(;
        base_dir = write_run(joinpath(root, "base"), base),
        fix_dir = write_run(joinpath(root, "fix"), fix; skip_seeds = skip_fix_seeds),
        base_oracle = write_oracle(joinpath(root, "base_checks.csv"); sha = base_sha),
        fix_oracle = write_oracle(
            joinpath(root, "fix_checks.csv");
            sha = fix_sha,
            fail_gate = fix_fail_gate,
            fail_advisory = fix_fail_advisory,
        ),
        base_head = BASE_HEAD,
        fix_head = fix_head,
        manifest = MANIFEST,
        output_dir = joinpath(root, "out"),
    )
end

Test.@testset "td4mbg one-sided sign-flip test matches hand enumeration" begin
    t = Analysis.signflip_signed_rank_test_greater
    # +-1 +-2 +-3: sums 6,4,2,0,0,-2,-4,-6.
    Test.@test t([1.0, 2.0, 3.0]).pvalue == 1 / 8
    Test.@test t([1.0, -2.0, 3.0]).pvalue == 3 / 8      # T_obs = 2
    Test.@test t([-1.0, -2.0, -3.0]).pvalue == 1.0
    z = t([0.0, 1.0, 2.0, 3.0])
    Test.@test z.n == 3 && z.n_zero_dropped == 1 && z.pvalue == 1 / 8
    Test.@test t([1.0, 1.0]).pvalue == 1 / 4            # tied ranks 1.5, 1.5
    Test.@test t([1e-9, 2.0]).n == 1                    # rounds to zero, dropped
    Test.@test isnan(t([0.0, 0.0]).pvalue)
    Test.@test_throws ArgumentError t([NaN, 1.0])
    # Same statistic as the two-sided test it mirrors: two-sided |T| >= 2 is 6/8.
    Test.@test Analysis.signflip_signed_rank_test([1.0, -2.0, 3.0]).pvalue == 6 / 8
end

Test.@testset "td4mbg FIX closer to zero in every seed is SUPPORTED" begin
    r = run_case()
    Test.@test r.verdict == "SUPPORTED"
    Test.@test r.g2.n_pairs == 16 && r.g2.n_x_positive == 16
    Test.@test r.g2.pvalue_one_sided == 1 / 2^16
    Test.@test DataFrames.nrow(r.g0a) == 0 && DataFrames.nrow(r.g0b) == 0   # wall ignored
    Test.@test r.g0c_ok
    Test.@test r.g2.median_x ≈ 0.0425                   # x(s) = 0.005 s, median of 1..16
end

Test.@testset "td4mbg FIX farther from zero is REFUTED" begin
    r = run_case(fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.02 * s)))
    Test.@test r.verdict == "REFUTED"
    Test.@test r.g2.n_x_positive == 0
    same = run_case(fix = synthetic_run(head = FIX_HEAD))
    Test.@test same.verdict == "REFUTED"                # all x = 0: undefined, not significant
    Test.@test isnan(same.g2.pvalue_one_sided)
end

Test.@testset "td4mbg G0 identity scope" begin
    bump(cell, arm) =
        r ->
            (r.cell == cell && r.arm == arm && r.seed == 3) ? merge(r, (n_contigs = 101,)) :
            r
    fix_with(t) = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s), tweak = t)
    a = run_case(fix = fix_with(bump("6000/1000/20x/e0.05", "illumina")))
    Test.@test a.verdict == "VOID (G0 failed)"
    Test.@test only(a.g0a.column) == "n_contigs" && only(a.g0a.seed) == 3
    b = run_case(fix = fix_with(bump("6000/2000/20x/e0.02", "nanopore")))
    Test.@test b.verdict == "VOID (G0 failed)" && DataFrames.nrow(b.g0b) == 1
    # A 1000 bp nanopore row is what the fix is allowed to change.
    c = run_case(fix = fix_with(bump("6000/1000/20x/e0.02", "nanopore")))
    Test.@test c.verdict == "SUPPORTED"
end

Test.@testset "td4mbg run labels" begin
    m = run_case(
        fix = synthetic_run(
            head = FIX_HEAD,
            manifest = "dddddddddddd",
            delta = fix_delta(s -> -0.005 * s),
        ),
    )
    Test.@test m.verdict == "VOID (run labels)"
    Test.@test any(occursin("FIX manifest", p) for p in m.label_problems)
    h = run_case(fix_head = "ffffffffffff")
    Test.@test h.verdict == "VOID (run labels)"
end

Test.@testset "td4mbg oracle G0c" begin
    Test.@test run_case(fix_sha = SHA_B).verdict == "VOID (G0 failed)"
    Test.@test run_case(fix_fail_gate = true).verdict == "VOID (G0 failed)"
    adv = run_case(fix_fail_advisory = true)
    Test.@test adv.verdict == "SUPPORTED" && adv.oracle_fix.n_passed == 18
    Test.@test adv.oracle_fix.n_rows == 19
end

Test.@testset "td4mbg missing shard is INCOMPLETE" begin
    Test.@test startswith(run_case(skip_fix_seeds = [7]).verdict, "INCOMPLETE")
end

Test.@testset "td4mbg outputs are written and deterministic" begin
    root = mktempdir()
    args(out) = (;
        base_dir = write_run(joinpath(root, "base"), synthetic_run(head = BASE_HEAD)),
        fix_dir = write_run(
            joinpath(root, "fix"),
            synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s)),
        ),
        base_oracle = write_oracle(joinpath(root, "b.csv")),
        fix_oracle = write_oracle(joinpath(root, "f.csv")),
        base_head = BASE_HEAD,
        fix_head = FIX_HEAD,
        manifest = MANIFEST,
        output_dir = out,
        context_dir = joinpath(root, "base"),
    )
    r1 = Analysis.td4mbg_analysis(; args(joinpath(root, "out1"))...)
    Analysis.td4mbg_analysis(; args(joinpath(root, "out2"))...)
    for f in (
        "report.md",
        "td4mbg_g2_summary.csv",
        "td4mbg_g2_per_seed.csv",
        "td4mbg_g0_mismatches.csv",
        "td4mbg_oracle.csv",
        "td4mbg_descriptive_counters.csv",
        "td4mbg_context_drift.csv",
        "td4mbg_g2_paired.png",
        "td4mbg_g2_paired.svg",
    )
        Test.@test isfile(joinpath(root, "out1", f))
    end
    Test.@test read(joinpath(root, "out1", "report.md")) ==
               read(joinpath(root, "out2", "report.md"))
    Test.@test occursin("NOT CAPTURED", read(joinpath(root, "out1", "report.md"), String))
    Test.@test all(r1.drift.n_seeds_identical .== 16)   # context = BASE itself
    Test.@test_throws ArgumentError Analysis.td4mbg_analysis_main(["--base-dir", "x"])
end

end
