module Td4mbgBaseFixAnalysisTest

# Contract tests for benchmarking/td4mbg_base_fix_analysis.jl, the pre-registered
# BASE vs FIX decision rule for the demoted-read reroute (amendments 4 and 5). No
# assembly runs here; the shards and oracle CSVs are synthetic.
#
# Every gate is reached through the real entry point `td4mbg_analysis` and gets a
# POSITIVE control (a clean FIX that moves the delta toward 0 is SUPPORTED) and
# NEGATIVE controls, on BOTH runs where a guard reads both:
#   - G2: FIX moving away is REFUTED; an overshoot past 0 is REFUTED (pins
#     x = |d_B| - |d_F|, not d_F - d_B); only the e0.01 cell decides.
#   - G0a: a change to ANY output column of an illumina row, in ANY of the 6
#     cells, VOIDS; an on-disk-only difference ("1.00" vs "1.0") VOIDS.
#   - G0b: a nanopore change in each 2000 bp cell VOIDS; one in a 1000 bp cell
#     does not.
#   - Labels: a head or manifest mismatch on either run VOIDS.
#   - G0c: a failed gate on either oracle or unequal shas VOIDS; a failed
#     advisory row and a historical-sha mismatch shared by both commits do not.
#   - Completeness: a missing shard or row on either run is INCOMPLETE.
#   - Errors, not verdicts: malformed heads, BASE == FIX, same directory or
#     oracle file, truncated or mis-shaped oracle files, wrong oracle provenance,
#     an empty-assembly pair, bad read counts, and mismatched shard columns.
# The one-sided test is checked against hand enumeration.

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
    "6000/2000/20x/e0.05"
]
const CELLS_2000 = ["6000/2000/20x/e0.01", "6000/2000/20x/e0.02", "6000/2000/20x/e0.05"]
# Every column G0 must compare, written out here rather than read from the
# production constant, so narrowing the comparison cannot narrow the test.
const OUTPUT_COLUMNS = [
    "substitution_length_divergences",
    "window_divergences",
    "window_anchor_rejections",
    "trace_contract_errors",
    "indel_engaged",
    "indel_requested",
    "reassembly_k",
    "n_contigs",
    "largest_contig",
    "largest_ratio",
    "total_contig_bases",
    "best_contig_reference_coverage",
    "best_contig_fit_identity"
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
        tweak = r -> r
)
    rows = NamedTuple[]
    for cell in CELLS, seed in 1:16, arm in ("illumina", "nanopore")
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
            best_contig_fit_identity = 1.0
        )
        push!(rows, tweak(row))
    end
    return DataFrames.DataFrame(rows)
end

"""
FIX-run nanopore delta: `f(seed)` in the 1000 bp cells (where the PR may change
results) and the BASE delta `base(seed)` in the 2000 bp cells (where G0b requires
identity).
"""
fix_delta(f; base = s -> -0.01 * s) = (c, s) -> startswith(c, "6000/1000/") ? f(s) :
                                                base(s)

"""Edit one (cell, seed, arm) row's `column` with `f`."""
function bump(cell, seed, arm, column, f)
    r -> (r.cell == cell && r.seed == seed && r.arm == arm) ?
         merge(r, NamedTuple{(column,)}((f(r[column]),))) : r
end

nudge(v::Integer) = v + 1
nudge(v::AbstractFloat) = v + 0.001

function write_run(dir, df; skip_seeds = Int[])
    mkpath(dir)
    for seed in 1:16
        seed in skip_seeds && continue
        CSV.write(joinpath(dir, "seed$(seed).csv"), df[df.seed .== seed, :])
    end
    return dir
end

"""
Write `fixed_toy_acceptance_checks.csv` and `fixed_toy_run_manifest.csv` into
`dir`, as the oracle generator does: 17 ordinary gates, one advisory row, and the
prewiring-sha gate whose `passed` is whether the sha equals the historical one.
"""
function write_oracle(dir;
        head,
        manifest = MANIFEST,
        sha = SHA_A,
        historical = true,
        fail_gate = false,
        fail_advisory = false,
        n_ordinary_gates = 17,
        severity = "gate",
        passed_text = nothing,
        dirty = false,
        write_manifest = true,
        first_advisory = false,
        extra_advisory = 0,
        advisory_severity = "advisory")
    mkpath(dir)
    rows = [(check = "gate_$(i)",
                severity = i == 1 && first_advisory ? "advisory" : severity,
                passed = i == 1 && passed_text !== nothing ? passed_text :
                         string(!(fail_gate && i == 1)),
                detail = "ok") for i in 1:n_ordinary_gates]
    push!(rows,
        (check = "nanopore_under_120_seconds", severity = advisory_severity,
            passed = string(!fail_advisory), detail = "wall"))
    for i in 1:extra_advisory
        push!(rows, (check = "extra_$(i)", severity = "advisory", passed = "true",
            detail = "extra"))
    end
    push!(rows,
        (check = "illumina_byte_identical_to_prewiring_oracle", severity = "gate",
            passed = string(historical),
            detail = "sha256=$(sha), origin_master=$(SHA_A)"))
    path = joinpath(dir, "fixed_toy_acceptance_checks.csv")
    CSV.write(path, DataFrames.DataFrame(rows))
    if write_manifest
        CSV.write(joinpath(dir, "fixed_toy_run_manifest.csv"),
            DataFrames.DataFrame(
                code_sha = [head * repeat("0", 28)],
                manifest_toml_sha256 = [manifest * repeat("0", 52)],
                git_tracked_worktree_dirty = [string(dirty)]))
    end
    return path
end

function run_case(;
        base = synthetic_run(head = BASE_HEAD),
        fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s), wall = 2.0),
        base_oracle = (;),
        fix_oracle = (;),
        skip_base_seeds = Int[],
        skip_fix_seeds = Int[],
        base_head = BASE_HEAD,
        fix_head = FIX_HEAD,
        manifest = MANIFEST,
        postprocess = root -> nothing,
        same_dir = false,
        same_oracle = false,
        context = false)
    root = mktempdir()
    base_dir = write_run(joinpath(root, "base"), base; skip_seeds = skip_base_seeds)
    fix_dir = same_dir ? base_dir :
              write_run(joinpath(root, "fix"), fix; skip_seeds = skip_fix_seeds)
    # merge, not a repeated keyword: an override may replace `head`.
    base_checks = write_oracle(joinpath(root, "oracle_base");
        merge((head = base_head,), base_oracle)...)
    fix_checks = same_oracle ? base_checks :
                 write_oracle(joinpath(root, "oracle_fix");
        merge((head = fix_head,), fix_oracle)...)
    postprocess(root)
    return Analysis.td4mbg_analysis(;
        base_dir, fix_dir, base_oracle = base_checks, fix_oracle = fix_checks,
        base_head, fix_head, manifest, output_dir = joinpath(root, "out"),
        context_dir = context ? base_dir : nothing)
end

const VOID_G0 = "VOID (G0 failed)"

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

Test.@testset "td4mbg G2 direction and gated cell" begin
    r = run_case()
    Test.@test r.verdict == "SUPPORTED"
    Test.@test r.g2.cell == "6000/1000/20x/e0.01"
    Test.@test r.g2.n_pairs == 16 && r.g2.n_x_positive == 16
    Test.@test r.g2.pvalue_one_sided == 1 / 2^16
    Test.@test DataFrames.nrow(r.g0a) == 0 && DataFrames.nrow(r.g0b) == 0   # wall ignored
    Test.@test r.g0c_ok
    Test.@test r.g2.median_x ≈ 0.0425                   # x(s) = 0.005 s, median of 1..16

    away = run_case(fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.02 * s)))
    Test.@test away.verdict == "REFUTED" && away.g2.n_x_positive == 0
    same = run_case(fix = synthetic_run(head = FIX_HEAD))
    Test.@test same.verdict == "REFUTED (undefined test: every x = 0)"
    Test.@test isnan(same.g2.pvalue_one_sided)

    # Overshoot: BASE -0.01 s, FIX +0.02 s. x = |d_B| - |d_F| = -0.01 s, so
    # REFUTED; a signed x = d_F - d_B = +0.03 s would wrongly say SUPPORTED.
    over = run_case(fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> 0.02 * s)))
    Test.@test over.verdict == "REFUTED" && over.g2.n_x_positive == 0
    # Mirror: BASE +0.01 s, FIX -0.005 s. x = +0.005 s, SUPPORTED; signed x is
    # -0.015 s.
    up = s -> 0.01 * s
    mirror = run_case(base = synthetic_run(head = BASE_HEAD, delta = (c, s) -> up(s)),
        fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s; base = up)))
    Test.@test mirror.verdict == "SUPPORTED"

    # Only e0.01 decides: toward 0 there and away in e0.02/e0.05 is SUPPORTED,
    # and the reverse is REFUTED.
    split(toward,
        away) = (c,
        s) -> c == "6000/1000/20x/e0.01" ? toward(s) :
              startswith(c, "6000/1000/") ? away(s) : -0.01 * s
    a = run_case(fix = synthetic_run(
        head = FIX_HEAD, delta = split(s -> -0.005 * s, s -> -0.02 * s)))
    Test.@test a.verdict == "SUPPORTED"
    b = run_case(fix = synthetic_run(
        head = FIX_HEAD, delta = split(s -> -0.02 * s, s -> -0.005 * s)))
    Test.@test b.verdict == "REFUTED"
end

Test.@testset "td4mbg G0a compares every output column in every cell" begin
    fix_with(t) = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s),
        tweak = t)
    for column in OUTPUT_COLUMNS
        r = run_case(fix = fix_with(bump(
            "6000/2000/20x/e0.05", 5, "illumina", Symbol(column), nudge)))
        Test.@test r.verdict == VOID_G0
        Test.@test only(r.g0a.column) == column
    end
    # n_reads must change in both arms (the arms must agree on read count).
    nr = run_case(fix = fix_with(
        r -> (r.cell == "6000/1000/20x/e0.02" && r.seed == 9) ?
             merge(r, (n_reads = 121,)) : r))
    Test.@test nr.verdict == VOID_G0 && "n_reads" in nr.g0a.column
    for cell in CELLS
        r = run_case(fix = fix_with(bump(cell, 3, "illumina", :n_contigs, nudge)))
        Test.@test r.verdict == VOID_G0 && only(r.g0a.cell) == cell
    end
end

Test.@testset "td4mbg G0a compares on-disk text, not parsed values" begin
    # "1.00" parses equal to "1.0" but is a different field on disk.
    rewrite = root -> begin
        path = joinpath(root, "fix", "seed4.csv")
        lines = readlines(path)
        i = only(findall(l -> occursin("6000/2000/20x/e0.01,4,illumina", l), lines))
        endswith(lines[i], ",1.0") || error("fixture: unexpected row $(lines[i])")
        lines[i] = lines[i][1:(end - 4)] * ",1.00"
        write(path, join(lines, "\n") * "\n")
    end
    r = run_case(postprocess = rewrite)
    Test.@test r.verdict == VOID_G0
    Test.@test only(r.g0a.column) == "best_contig_fit_identity"
end

Test.@testset "td4mbg G0b scope" begin
    fix_with(t) = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s),
        tweak = t)
    for cell in CELLS_2000
        r = run_case(fix = fix_with(bump(cell, 2, "nanopore", :n_contigs, nudge)))
        Test.@test r.verdict == VOID_G0 && only(r.g0b.cell) == cell
    end
    # A 1000 bp nanopore row is what the fix is allowed to change.
    c = run_case(fix = fix_with(bump(
        "6000/1000/20x/e0.02", 3, "nanopore", :n_contigs, nudge)))
    Test.@test c.verdict == "SUPPORTED"
end

Test.@testset "td4mbg run labels" begin
    toward = fix_delta(s -> -0.005 * s)
    m = run_case(fix = synthetic_run(head = FIX_HEAD, manifest = "dddddddddddd", delta = toward))
    Test.@test m.verdict == "VOID (run labels)"
    Test.@test any(occursin("FIX manifest", p) for p in m.label_problems)
    bm = run_case(base = synthetic_run(head = BASE_HEAD, manifest = "dddddddddddd"))
    Test.@test bm.verdict == "VOID (run labels)"
    Test.@test any(occursin("BASE manifest", p) for p in bm.label_problems)
    h = run_case(fix_head = "ffffffffffff")
    Test.@test h.verdict == "VOID (run labels)"
    Test.@test any(occursin("FIX head", p) for p in h.label_problems)
    bh = run_case(base_head = "ffffffffffff")
    Test.@test bh.verdict == "VOID (run labels)"
    Test.@test any(occursin("BASE head", p) for p in bh.label_problems)
end

Test.@testset "td4mbg oracle G0c" begin
    Test.@test run_case(fix_oracle = (sha = SHA_B, historical = false)).verdict == VOID_G0
    Test.@test run_case(fix_oracle = (fail_gate = true,)).verdict == VOID_G0
    Test.@test run_case(base_oracle = (fail_gate = true,)).verdict == VOID_G0
    adv = run_case(fix_oracle = (fail_advisory = true,))
    Test.@test adv.verdict == "SUPPORTED" && adv.oracle_fix.n_passed == 18
    # Amendment 4: equality with the historical sha is reported, not gated. Both
    # commits drift to the same new sha (the generator then writes passed=false
    # on that gate row): G0c still passes.
    drift = run_case(base_oracle = (sha = SHA_B, historical = false),
        fix_oracle = (sha = SHA_B, historical = false))
    Test.@test drift.verdict == "SUPPORTED"
    Test.@test isempty(drift.oracle_fix.failed_gates)
end

Test.@testset "td4mbg oracle shape and provenance are errors, not verdicts" begin
    E = ArgumentError
    Test.@test_throws E run_case(fix_oracle = (n_ordinary_gates = 16,))   # truncated
    # The row-count and gate-count checks each need a case only they catch.
    Test.@test_throws E run_case(fix_oracle = (extra_advisory = 1,))      # 20 rows, 18 gates
    Test.@test_throws E run_case(fix_oracle = (first_advisory = true,))   # 19 rows, 17 gates
    # Row and gate counts both correct; only the severity vocabulary is wrong.
    Test.@test_throws E run_case(fix_oracle = (advisory_severity = "Advisory",))
    Test.@test_throws E run_case(fix_oracle = (severity = "Gate",))
    Test.@test_throws E run_case(fix_oracle = (passed_text = "",))
    Test.@test_throws E run_case(fix_oracle = (head = "eeeeeeeeeeee",))  # other commit
    Test.@test_throws E run_case(fix_oracle = (manifest = "dddddddddddd",))
    Test.@test_throws E run_case(fix_oracle = (dirty = true,))
    Test.@test_throws E run_case(fix_oracle = (write_manifest = false,))
    Test.@test_throws E run_case(same_oracle = true)
end

Test.@testset "td4mbg operator input is checked" begin
    E = ArgumentError
    Test.@test_throws E run_case(base_head = BASE_HEAD * repeat("0", 28))   # full sha
    Test.@test_throws E run_case(fix_head = uppercase(FIX_HEAD))
    Test.@test_throws E run_case(fix_head = BASE_HEAD,
        fix = synthetic_run(head = BASE_HEAD))                              # BASE == FIX
    Test.@test_throws E run_case(same_dir = true)
    extra = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s))
    extra.extra_column = fill(1, DataFrames.nrow(extra))
    Test.@test_throws E run_case(fix = extra)
    Test.@test_throws E Analysis.td4mbg_analysis_main(["--base-dir", "x"])
end

Test.@testset "td4mbg pipeline failures are errors (td2rxh amendment 3)" begin
    E = ArgumentError
    empty = r -> (r.cell == "6000/1000/20x/e0.01" && r.seed == 6) ?
                 merge(r, (n_contigs = 0,)) : r
    Test.@test_throws E run_case(base = synthetic_run(head = BASE_HEAD, tweak = empty),
        fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s),
            tweak = empty))
    Test.@test_throws E run_case(fix = synthetic_run(head = FIX_HEAD,
        delta = fix_delta(s -> -0.005 * s),
        tweak = bump("6000/1000/20x/e0.05", 8, "nanopore", :n_reads, nudge)))
    zero_reads = r -> (r.cell == "6000/2000/20x/e0.01" && r.seed == 1) ?
                      merge(r, (n_reads = 0,)) : r
    Test.@test_throws E run_case(
        base = synthetic_run(head = BASE_HEAD, tweak = zero_reads),
        fix = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s),
            tweak = zero_reads))
end

Test.@testset "td4mbg missing data is INCOMPLETE on either run" begin
    Test.@test startswith(run_case(skip_fix_seeds = [7]).verdict, "INCOMPLETE")
    Test.@test startswith(run_case(skip_base_seeds = [12]).verdict, "INCOMPLETE")
    partial = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s))
    drop = (partial.cell .== "6000/2000/20x/e0.02") .& (partial.seed .== 5) .&
           (partial.arm .== "nanopore")
    Test.@test startswith(run_case(fix = partial[.!drop, :]).verdict, "INCOMPLETE")
end

Test.@testset "td4mbg outputs are complete, readable and deterministic" begin
    root = mktempdir()
    base_dir = write_run(joinpath(root, "base"), synthetic_run(head = BASE_HEAD))
    # Every illumina row differs in total_contig_bases: 96 G0a fields, more than
    # any fixed print limit.
    many = synthetic_run(head = FIX_HEAD, delta = fix_delta(s -> -0.005 * s),
        tweak = r -> r.arm == "illumina" ? merge(r, (total_contig_bases = 9001,)) : r)
    fix_dir = write_run(joinpath(root, "fix"), many)
    args(out) = (; base_dir, fix_dir,
        base_oracle = write_oracle(joinpath(root, "ob"); head = BASE_HEAD),
        fix_oracle = write_oracle(joinpath(root, "of"); head = FIX_HEAD),
        base_head = BASE_HEAD, fix_head = FIX_HEAD, manifest = MANIFEST,
        output_dir = out, context_dir = base_dir)
    r1 = Analysis.td4mbg_analysis(; args(joinpath(root, "out1"))...)
    Analysis.td4mbg_analysis(; args(joinpath(root, "out2"))...)
    for f in ("report.md", "td4mbg_g2_summary.csv", "td4mbg_g2_per_seed.csv",
        "td4mbg_g0_mismatches.csv", "td4mbg_oracle.csv", "td4mbg_descriptive_counters.csv",
        "td4mbg_context_drift.csv", "td4mbg_g2_paired.png", "td4mbg_g2_paired.svg")
        Test.@test isfile(joinpath(root, "out1", f))
    end
    report = read(joinpath(root, "out1", "report.md"), String)
    Test.@test read(joinpath(root, "out1", "report.md")) ==
               read(joinpath(root, "out2", "report.md"))
    Test.@test r1.verdict == VOID_G0 && DataFrames.nrow(r1.g0a) == 96
    Test.@test count(l -> startswith(l, "| G0a |"), split(report, "\n")) == 96
    Test.@test occursin("descriptively only", report)
    Test.@test occursin("NOT CAPTURED", report)
    Test.@test occursin("| true |", report)               # Bools print as Bools
    Test.@test !occursin(r"\| G0a \| [^|]+ \| \d+\.0 \|", report)  # seeds as Ints
    Test.@test all(r1.drift.n_seeds_identical .== 16)     # context = BASE itself
end

end
