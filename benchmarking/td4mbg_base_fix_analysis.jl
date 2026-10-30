# td-4mbg pre-registered BASE vs FIX analysis of two `--cells td2rxh` array runs.
# ==============================================================================
#
# This file implements the decision rule recorded on bead td-4mbg BEFORE either
# array run existed: the original pre-registration as replaced by AMENDMENTS 1-4,
# plus AMENDMENT 5 (recorded before any shard was analyzed). The bead notes are
# the authoritative text; this is a condensed transcription. Summary:
#
#   RUNS: BASE = PR #487's merge base, FIX = PR #487's head. Both are
#       run_td2rxh_evidence_lrc.sbatch --cells td2rxh, seeds 1..16, on ONE shared
#       Manifest.toml. The two run labels must carry the expected heads and the
#       same manifest digest, otherwise the comparison is VOID before analysis.
#   G0 (controls; any failure VOIDS the comparison):
#       (a) the :illumina rows are identical between BASE and FIX in all 6 cells;
#       (b) the :nanopore rows are identical between BASE and FIX in the three
#           2000 bp cells (every read there is long, so PR #487 cannot reroute it);
#       (c) the fixed-toy oracle passes every gate row at BOTH commits, and its
#           illumina sha is equal between them. Equality with the historical
#           prewiring sha is reported, not gated.
#       "Identical" (amendment 5) = exact CSV field-string equality on every
#       column except `label` (embeds the commit) and the two wall-time columns.
#   G2 (outcome): at 6000/1000/20x/e0.01, per seed s, d(s) = nanopore - illumina
#       best_contig_reference_coverage in each run, and x(s) = |d_B(s)| - |d_F(s)|.
#       H1: x is shifted above 0 (FIX moves the paired delta toward 0). Exact
#       ONE-sided sign-flip test on the Wilcoxon signed-rank statistic, values
#       rounded to 1e-6, tie-averaged ranks, zeros dropped (as in
#       td2rxh_paired_analysis.jl), alpha 0.05. SUPPORTED iff p < alpha, else
#       REFUTED. The 1000 bp err-0.02 and err-0.05 cells get the same statistic,
#       reported without a gate.
#   DESCRIPTIVE ONLY: the fail-open counters (substitution_length_divergences,
#       window_divergences, window_anchor_rejections) and wall time per read, per
#       cell, arm and run, with paired FIX - BASE medians. demoted_whole_read_routes
#       is NOT CAPTURED: the harness does not export it (amendment 5).
#   CONTEXT ONLY (optional --context-dir): an earlier run on master, compared with
#       BASE cell by cell to measure dependency and master drift. It enters no gate.
#
# Order of evaluation: run labels -> completeness (every expected row in both
# runs) -> G0 -> G2. Missing rows give INCOMPLETE, never VOID or REFUTED.
#
# Usage:
#   julia --project=. benchmarking/td4mbg_base_fix_analysis.jl \
#       --base-dir <BASE shards> --fix-dir <FIX shards> \
#       --base-oracle <BASE fixed_toy_acceptance_checks.csv> \
#       --fix-oracle <FIX fixed_toy_acceptance_checks.csv> \
#       --base-head <12 hex> --fix-head <12 hex> --manifest <12 hex> \
#       --output-dir <dir> [--context-dir <earlier run shards>]

include(joinpath(@__DIR__, "td2rxh_paired_analysis.jl"))

const TD4MBG_G2_CELL = "6000/1000/20x/e0.01"
const TD4MBG_G2_CELLS =
    ["6000/1000/20x/e0.01", "6000/1000/20x/e0.02", "6000/1000/20x/e0.05"]
const TD4MBG_G0B_CELLS =
    ["6000/2000/20x/e0.01", "6000/2000/20x/e0.02", "6000/2000/20x/e0.05"]
const TD4MBG_ARMS = ("illumina", "nanopore")
const TD4MBG_NON_OUTPUT_COLUMNS = ("label", "wall_seconds", "wall_seconds_per_read")
const TD4MBG_DESCRIPTIVE_METRICS = (
    :substitution_length_divergences,
    :window_divergences,
    :window_anchor_rejections,
    :wall_seconds_per_read,
)
const TD4MBG_ORACLE_SHA_CHECK = "illumina_byte_identical_to_prewiring_oracle"
const TD4MBG_PREWIRING_SHA = "d36e3b6a10685346aa7b0238b48b4ab7fcefbed88f82cad7d959b0a831cdd311"
const TD4MBG_ORACLE_EXPECTED_ROWS = 19
const TD4MBG_RUN_LABEL_PATTERN = r"^td2rxh@([0-9a-f]{12})\+manifest([0-9a-f]{12})$"

"""
    signflip_signed_rank_test_greater(diffs) -> NamedTuple

Exact ONE-sided (H1: shifted above 0) sign-flip permutation test on the Wilcoxon
signed-rank statistic. Identical to `signflip_signed_rank_test` in rounding (to
`TD2RXH_DELTA_DIGITS`), zero handling (dropped and counted) and ranks
(tie-averaged), except that the p-value is the fraction of all `2^n` sign
assignments with `T >= T_obs`, where `T = sum(sign(d) * rank(|d|))`. With no
non-zero difference the test is undefined and `pvalue` is `NaN`.
"""
function signflip_signed_rank_test_greater(diffs::AbstractVector{<:Real})::NamedTuple
    rounded = round.(Float64.(diffs); digits = TD2RXH_DELTA_DIGITS)
    any(!isfinite, rounded) &&
        throw(ArgumentError("non-finite paired difference; drop or fix the pair first"))
    nonzero = filter(!iszero, rounded)
    n = length(nonzero)
    n_zero = length(rounded) - n
    n == 0 && return (
        n = 0,
        n_zero_dropped = n_zero,
        statistic = NaN,
        pvalue = NaN,
        method = "undefined (no non-zero differences)",
    )
    n > TD2RXH_MAX_EXACT_N && throw(
        ArgumentError("n = $(n) exceeds the exact enumeration limit $(TD2RXH_MAX_EXACT_N)"),
    )
    ranks = average_ranks(abs.(nonzero))
    observed = sum(sign.(nonzero) .* ranks)
    # Ranks are multiples of 1/2, so a tolerance far below 1/2 compares exactly.
    tolerance = 1e-9
    extreme = 0
    for mask = 0:((1<<n)-1)
        t = 0.0
        for i = 1:n
            t += ((mask >> (i - 1)) & 1) == 1 ? ranks[i] : -ranks[i]
        end
        t >= observed - tolerance && (extreme += 1)
    end
    return (
        n = n,
        n_zero_dropped = n_zero,
        statistic = observed,
        pvalue = extreme / (1 << n),
        method = "exact one-sided sign-flip ($(1 << n) assignments)",
    )
end

"""
    parse_run_label(label) -> NamedTuple

Split a `td2rxh@<12 hex head>+manifest<12 hex digest>` run label.
"""
function parse_run_label(label::AbstractString)::NamedTuple
    m = match(TD4MBG_RUN_LABEL_PATTERN, label)
    m === nothing && throw(
        ArgumentError("run label \"$(label)\" is not td2rxh@<12 hex>+manifest<12 hex>"),
    )
    return (head = String(m.captures[1]), manifest = String(m.captures[2]))
end

"""
    load_run(dir) -> NamedTuple

Read every `seed<N>.csv` shard in `dir` twice: typed (for statistics, checked by
`validate_shards`) and as strings (for the G0 field-identity checks, so a value
that parses equal but differs on disk is still a difference).
"""
function load_run(dir::AbstractString)::NamedTuple
    isdir(dir) || throw(ArgumentError("shard directory not found: $(dir)"))
    paths = sort(
        filter(p -> occursin(r"^seed\d+\.csv$", basename(p)), readdir(dir; join = true)),
    )
    isempty(paths) && throw(ArgumentError("no seed*.csv shards in $(dir)"))
    typed = DataFrames.vcat(
        [CSV.read(p, DataFrames.DataFrame) for p in paths]...;
        cols = :setequal,
    )
    strings = DataFrames.vcat(
        [CSV.read(p, DataFrames.DataFrame; types = String) for p in paths]...;
        cols = :setequal,
    )
    validate_shards(typed)
    label = only(unique(string.(typed.label)))
    paired, unpaired = pair_arms(typed)
    return (
        dir = dir,
        paths = paths,
        typed = typed,
        strings = strings,
        label = label,
        parsed_label = parse_run_label(label),
        paired = paired,
        unpaired = unpaired,
    )
end

"""
    label_problems(base, fix; base_head, fix_head, manifest) -> Vector{String}

Every way the two run labels fail amendment 4: wrong head on either run, or a
manifest digest that is not the expected shared one. Empty means the labels pass.
"""
function label_problems(
    base::NamedTuple,
    fix::NamedTuple;
    base_head::AbstractString,
    fix_head::AbstractString,
    manifest::AbstractString,
)::Vector{String}
    problems = String[]
    base.parsed_label.head == base_head ||
        push!(problems, "BASE head $(base.parsed_label.head) != expected $(base_head)")
    fix.parsed_label.head == fix_head ||
        push!(problems, "FIX head $(fix.parsed_label.head) != expected $(fix_head)")
    for (name, run) in (("BASE", base), ("FIX", fix))
        run.parsed_label.manifest == manifest || push!(
            problems,
            "$(name) manifest $(run.parsed_label.manifest) != expected $(manifest)",
        )
    end
    return problems
end

function _string_row(
    df::DataFrames.DataFrame,
    cell::AbstractString,
    seed::Integer,
    arm::AbstractString,
)
    idx = findall(
        isequal.(df.cell, cell) .& isequal.(df.seed, string(seed)) .& isequal.(df.arm, arm),
    )
    return length(idx) == 1 ? df[only(idx), :] : nothing
end

_field(x) = ismissing(x) ? "<empty>" : String(x)

"""
    field_mismatches(a, b; arm, cells) -> DataFrame

One row per (cell, seed, column) where the string-typed shards `a` and `b`
differ, over every column except `TD4MBG_NON_OUTPUT_COLUMNS`, for `arm` in
`cells`. A (cell, seed) absent from either side is one row with column `<row>`.
"""
function field_mismatches(
    a::DataFrames.DataFrame,
    b::DataFrames.DataFrame;
    arm::AbstractString,
    cells::Vector{String},
)::DataFrames.DataFrame
    sort(names(a)) == sort(names(b)) ||
        throw(ArgumentError("the two runs' shards have different columns"))
    columns = [c for c in names(a) if !(c in TD4MBG_NON_OUTPUT_COLUMNS)]
    out = DataFrames.DataFrame(
        cell = String[],
        seed = Int[],
        arm = String[],
        column = String[],
        a = String[],
        b = String[],
    )
    for cell in cells, seed = 1:TD2RXH_N_SEEDS

        ra = _string_row(a, cell, seed, arm)
        rb = _string_row(b, cell, seed, arm)
        if ra === nothing || rb === nothing
            push!(
                out,
                (
                    cell,
                    seed,
                    arm,
                    "<row>",
                    ra === nothing ? "absent" : "present",
                    rb === nothing ? "absent" : "present",
                ),
            )
            continue
        end
        for c in columns
            isequal(ra[c], rb[c]) ||
                push!(out, (cell, seed, arm, c, _field(ra[c]), _field(rb[c])))
        end
    end
    return out
end

"""
    oracle_summary(path) -> NamedTuple

Read a `fixed_toy_acceptance_checks.csv`. `gates_pass` is true iff every
`severity == "gate"` row has `passed == "true"`; the illumina sha is parsed from
the `$(TD4MBG_ORACLE_SHA_CHECK)` row's detail.
"""
function oracle_summary(path::AbstractString)::NamedTuple
    df = CSV.read(path, DataFrames.DataFrame; types = String)
    for c in ("check", "severity", "passed", "detail")
        c in names(df) || throw(ArgumentError("$(path) lacks column $(c)"))
    end
    passed = [isequal(x, "true") for x in df.passed]
    is_gate = [isequal(x, "gate") for x in df.severity]
    failed_gates =
        [string(df.check[i]) for i in eachindex(passed) if is_gate[i] && !passed[i]]
    row = findfirst(isequal(TD4MBG_ORACLE_SHA_CHECK), df.check)
    row === nothing && throw(ArgumentError("$(path) has no $(TD4MBG_ORACLE_SHA_CHECK) row"))
    m = match(r"sha256=([0-9a-f]{64})", coalesce(df.detail[row], ""))
    m === nothing &&
        throw(ArgumentError("$(path): no sha256 in the $(TD4MBG_ORACLE_SHA_CHECK) detail"))
    return (
        path = path,
        n_rows = DataFrames.nrow(df),
        n_passed = count(passed),
        n_gates = count(is_gate),
        failed_gates = failed_gates,
        gates_pass = isempty(failed_gates),
        illumina_sha = String(m.captures[1]),
    )
end

function _paired_row(paired::DataFrames.DataFrame, cell::AbstractString, seed::Integer)
    idx = findall((string.(paired.cell) .== cell) .& (Int.(paired.seed) .== seed))
    return length(idx) == 1 ? paired[only(idx), :] : nothing
end

function _coverage_delta(row)::Float64
    return round(
        Float64(row.best_contig_reference_coverage_nanopore) -
        Float64(row.best_contig_reference_coverage_illumina);
        digits = TD2RXH_DELTA_DIGITS,
    )
end

_median_or_nan(x) = isempty(x) ? NaN : Statistics.median(x)

"""
    g2_cell(base, fix, cell) -> (summary::NamedTuple, per_seed::DataFrame)

The G2 statistic for one cell: per seed d_B, d_F (nanopore - illumina coverage,
rounded) and x = |d_B| - |d_F|, then `signflip_signed_rank_test_greater(x)`.
`status` is COMPLETE only with all `TD2RXH_N_SEEDS` seeds paired in both runs.
"""
function g2_cell(base::NamedTuple, fix::NamedTuple, cell::AbstractString)::Tuple
    per_seed = DataFrames.DataFrame(
        cell = String[],
        seed = Int[],
        d_base = Float64[],
        d_fix = Float64[],
        x = Float64[],
        d_fix_minus_d_base = Float64[],
    )
    for seed = 1:TD2RXH_N_SEEDS
        rb = _paired_row(base.paired, cell, seed)
        rf = _paired_row(fix.paired, cell, seed)
        (rb === nothing || rf === nothing) && continue
        d_b = _coverage_delta(rb)
        d_f = _coverage_delta(rf)
        push!(
            per_seed,
            (
                cell,
                seed,
                d_b,
                d_f,
                round(abs(d_b) - abs(d_f); digits = TD2RXH_DELTA_DIGITS),
                round(d_f - d_b; digits = TD2RXH_DELTA_DIGITS),
            ),
        )
    end
    n_pairs = DataFrames.nrow(per_seed)
    test =
        n_pairs == 0 ?
        (n = 0, n_zero_dropped = 0, statistic = NaN, pvalue = NaN, method = "no pairs") :
        signflip_signed_rank_test_greater(per_seed.x)
    summary = (
        cell = String(cell),
        gated = cell == TD4MBG_G2_CELL,
        status = n_pairs == TD2RXH_N_SEEDS ? "COMPLETE" : "INCOMPLETE",
        n_pairs = n_pairs,
        n_nonzero = test.n,
        n_zero = test.n_zero_dropped,
        statistic = test.statistic,
        pvalue_one_sided = test.pvalue,
        median_x = _median_or_nan(per_seed.x),
        median_d_fix_minus_d_base = _median_or_nan(per_seed.d_fix_minus_d_base),
        n_x_positive = count(>(0), per_seed.x),
        median_d_base = _median_or_nan(per_seed.d_base),
        median_d_fix = _median_or_nan(per_seed.d_fix),
        method = test.method,
    )
    return summary, per_seed
end

function _metric_by_seed(
    df::DataFrames.DataFrame,
    cell,
    arm,
    metric::Symbol,
)::Dict{Int,Float64}
    out = Dict{Int,Float64}()
    metric_name = string(metric)
    metric_name in names(df) || return out
    for r in eachrow(df)
        (string(r.cell) == cell && string(r.arm) == arm) || continue
        v = r[metric_name]
        ismissing(v) || (out[Int(r.seed)] = Float64(v))
    end
    return out
end

"""
    descriptive_counters(base, fix) -> DataFrame

DESCRIPTIVE ONLY (amendment 2): per cell, arm and metric, the BASE and FIX
medians over seeds and the median paired FIX - BASE difference.
"""
function descriptive_counters(base::NamedTuple, fix::NamedTuple)::DataFrames.DataFrame
    out = DataFrames.DataFrame(
        cell = String[],
        arm = String[],
        metric = String[],
        median_base = Float64[],
        median_fix = Float64[],
        median_fix_minus_base = Float64[],
        n_pairs = Int[],
    )
    for cell in TD2RXH_EXPECTED_CELLS,
        arm in TD4MBG_ARMS,
        metric in TD4MBG_DESCRIPTIVE_METRICS

        b = _metric_by_seed(base.typed, cell, arm, metric)
        f = _metric_by_seed(fix.typed, cell, arm, metric)
        seeds = sort(collect(intersect(keys(b), keys(f))))
        push!(
            out,
            (
                cell,
                arm,
                string(metric),
                _median_or_nan(collect(values(b))),
                _median_or_nan(collect(values(f))),
                _median_or_nan([f[s] - b[s] for s in seeds]),
                length(seeds),
            ),
        )
    end
    return out
end

"""
    context_drift(context, base) -> DataFrame

CONTEXT ONLY: per cell and arm, how many seeds of an earlier run are field-
identical to BASE (same column rule as G0), and the median BASE - context
coverage difference.
"""
function context_drift(context::NamedTuple, base::NamedTuple)::DataFrames.DataFrame
    out = DataFrames.DataFrame(
        cell = String[],
        arm = String[],
        n_seeds_identical = Int[],
        n_seeds_compared = Int[],
        median_coverage_base_minus_context = Float64[],
    )
    for cell in TD2RXH_EXPECTED_CELLS, arm in TD4MBG_ARMS

        mismatches =
            field_mismatches(context.strings, base.strings; arm = arm, cells = [cell])
        b = _metric_by_seed(base.typed, cell, arm, :best_contig_reference_coverage)
        c = _metric_by_seed(context.typed, cell, arm, :best_contig_reference_coverage)
        seeds = sort(collect(intersect(keys(b), keys(c))))
        differing = unique(mismatches.seed)
        push!(
            out,
            (
                cell,
                arm,
                count(s -> !(s in differing), seeds),
                length(seeds),
                _median_or_nan([
                    round(b[s] - c[s]; digits = TD2RXH_DELTA_DIGITS) for s in seeds
                ]),
            ),
        )
    end
    return out
end

function td4mbg_decision(;
    labels_ok::Bool,
    complete::Bool,
    g0a_ok::Bool,
    g0b_ok::Bool,
    g0c_ok::Bool,
    g2::NamedTuple,
)::String
    labels_ok || return "VOID (run labels)"
    complete || return "INCOMPLETE (missing rows; analyze only a complete run)"
    (g0a_ok && g0b_ok && g0c_ok) || return "VOID (G0 failed)"
    g2.status == "COMPLETE" || return "INCOMPLETE (G2 cell)"
    return !isnan(g2.pvalue_one_sided) && g2.pvalue_one_sided < TD2RXH_ALPHA ? "SUPPORTED" :
           "REFUTED"
end

function _fmt(x)::String
    x isa AbstractFloat || return string(x)
    isnan(x) && return "NaN"
    return string(round(x; sigdigits = 4))
end

function write_td4mbg_figure(
    g2_results::Vector,
    output_dir::AbstractString,
    base_head::AbstractString,
    fix_head::AbstractString,
)::Nothing
    fig = CairoMakie.Figure(size = (1200, 520))
    for (i, (summary, per_seed)) in enumerate(g2_results)
        ax = CairoMakie.Axis(
            fig[1, i];
            xticks = ([1, 2], ["BASE $(base_head[1:8])", "FIX $(fix_head[1:8])"]),
            ylabel = i == 1 ? "nanopore − illumina best-contig reference coverage" : "",
            title = "$(summary.cell) ($(summary.gated ? "G2, gated" : "descriptive"))\n" *
                    "n=$(summary.n_pairs), x>0 in $(summary.n_x_positive), " *
                    "one-sided p=$(_fmt(summary.pvalue_one_sided))",
        )
        CairoMakie.hlines!(ax, [0.0]; color = :gray50, linestyle = :dash)
        for r in eachrow(per_seed)
            color = r.x > 0 ? :steelblue : (r.x < 0 ? :firebrick : :gray60)
            CairoMakie.lines!(ax, [1, 2], [r.d_base, r.d_fix]; color = color)
            CairoMakie.scatter!(
                ax,
                [1, 2],
                [r.d_base, r.d_fix];
                color = color,
                markersize = 6,
            )
        end
        CairoMakie.xlims!(ax, 0.7, 2.3)
    end
    CairoMakie.Label(
        fig[2, :],
        "One line per seed (6000 bp genome, 1000 bp reads, 20x). Blue: FIX moves the delta " *
        "toward 0 (x = |d_BASE| − |d_FIX| > 0); red: away; grey: unchanged. Exact one-sided " *
        "sign-flip signed-rank on x, alpha 0.05; only e0.01 is a pre-registered gate.";
        tellwidth = false,
        fontsize = 12,
    )
    CairoMakie.save(joinpath(output_dir, "td4mbg_g2_paired.png"), fig; px_per_unit = 2)
    CairoMakie.save(joinpath(output_dir, "td4mbg_g2_paired.svg"), fig)
    return nothing
end

function _table(io::IO, df::DataFrames.DataFrame)::Nothing
    cols = names(df)
    println(io, "| ", join(cols, " | "), " |")
    println(io, "|", repeat(" --- |", length(cols)))
    for r in eachrow(df)
        println(io, "| ", join([_fmt(r[c]) for c in cols], " | "), " |")
    end
    println(io)
    return nothing
end

"""
    td4mbg_analysis(; base_dir, fix_dir, base_oracle, fix_oracle, base_head, fix_head,
        manifest, output_dir, context_dir = nothing) -> NamedTuple

Apply the td-4mbg decision rule, write the report, CSVs and figure to
`output_dir`, and return the verdict with its components.
"""
function td4mbg_analysis(;
    base_dir::AbstractString,
    fix_dir::AbstractString,
    base_oracle::AbstractString,
    fix_oracle::AbstractString,
    base_head::AbstractString,
    fix_head::AbstractString,
    manifest::AbstractString,
    output_dir::AbstractString,
    context_dir::Union{Nothing,AbstractString} = nothing,
)::NamedTuple
    mkpath(output_dir)
    base = load_run(base_dir)
    fix = load_run(fix_dir)
    problems = label_problems(base, fix; base_head, fix_head, manifest)
    complete = DataFrames.nrow(base.unpaired) == 0 && DataFrames.nrow(fix.unpaired) == 0
    g0a = field_mismatches(
        base.strings,
        fix.strings;
        arm = "illumina",
        cells = TD2RXH_EXPECTED_CELLS,
    )
    g0b = field_mismatches(
        base.strings,
        fix.strings;
        arm = "nanopore",
        cells = TD4MBG_G0B_CELLS,
    )
    oracle_base = oracle_summary(base_oracle)
    oracle_fix = oracle_summary(fix_oracle)
    g0c_ok =
        oracle_base.gates_pass &&
        oracle_fix.gates_pass &&
        oracle_base.illumina_sha == oracle_fix.illumina_sha
    g2_results = [g2_cell(base, fix, cell) for cell in TD4MBG_G2_CELLS]
    g2 = first(only(filter(r -> first(r).gated, g2_results)))
    verdict = td4mbg_decision(;
        labels_ok = isempty(problems),
        complete,
        g0a_ok = DataFrames.nrow(g0a) == 0,
        g0b_ok = DataFrames.nrow(g0b) == 0,
        g0c_ok,
        g2,
    )
    counters = descriptive_counters(base, fix)
    drift = context_dir === nothing ? nothing : context_drift(load_run(context_dir), base)

    g0_rows = vcat(
        DataFrames.insertcols(g0a, 1, :gate => "G0a"),
        DataFrames.insertcols(g0b, 1, :gate => "G0b"),
    )
    g2_summary = DataFrames.DataFrame([first(r) for r in g2_results])
    g2_per_seed = vcat([last(r) for r in g2_results]...)
    oracle = DataFrames.DataFrame([
        (
            run = "BASE",
            rows_passed = "$(oracle_base.n_passed)/$(oracle_base.n_rows)",
            gates_pass = oracle_base.gates_pass,
            failed_gates = join(oracle_base.failed_gates, ";"),
            illumina_sha = oracle_base.illumina_sha,
            equals_historical = oracle_base.illumina_sha == TD4MBG_PREWIRING_SHA,
        ),
        (
            run = "FIX",
            rows_passed = "$(oracle_fix.n_passed)/$(oracle_fix.n_rows)",
            gates_pass = oracle_fix.gates_pass,
            failed_gates = join(oracle_fix.failed_gates, ";"),
            illumina_sha = oracle_fix.illumina_sha,
            equals_historical = oracle_fix.illumina_sha == TD4MBG_PREWIRING_SHA,
        ),
    ])
    CSV.write(joinpath(output_dir, "td4mbg_g0_mismatches.csv"), g0_rows)
    CSV.write(joinpath(output_dir, "td4mbg_oracle.csv"), oracle)
    CSV.write(joinpath(output_dir, "td4mbg_g2_summary.csv"), g2_summary)
    CSV.write(joinpath(output_dir, "td4mbg_g2_per_seed.csv"), g2_per_seed)
    CSV.write(joinpath(output_dir, "td4mbg_descriptive_counters.csv"), counters)
    drift === nothing || CSV.write(joinpath(output_dir, "td4mbg_context_drift.csv"), drift)
    write_td4mbg_figure(
        g2_results,
        output_dir,
        base.parsed_label.head,
        fix.parsed_label.head,
    )

    open(joinpath(output_dir, "report.md"), "w") do io
        println(io, "# td-4mbg BASE vs FIX: pre-registered analysis\n")
        println(io, "**Verdict: $(verdict)**\n")
        println(io, "## Inputs\n")
        println(
            io,
            "- BASE label `$(base.label)` ($(length(base.paths)) shards, $(DataFrames.nrow(base.typed)) rows)",
        )
        println(
            io,
            "- FIX label `$(fix.label)` ($(length(fix.paths)) shards, $(DataFrames.nrow(fix.typed)) rows)",
        )
        println(
            io,
            "- Expected heads BASE `$(base_head)`, FIX `$(fix_head)`, shared manifest `$(manifest)`",
        )
        println(io, "- Label problems: ", isempty(problems) ? "none" : join(problems, "; "))
        println(
            io,
            "- Unpaired (cell, seed): BASE $(DataFrames.nrow(base.unpaired)), FIX $(DataFrames.nrow(fix.unpaired))\n",
        )
        println(io, "## G0 controls\n")
        println(
            io,
            "- G0a, illumina identical in all 6 cells: ",
            DataFrames.nrow(g0a) == 0 ? "PASS" :
            "FAIL ($(DataFrames.nrow(g0a)) mismatching fields)",
        )
        println(
            io,
            "- G0b, nanopore identical in the 2000 bp cells: ",
            DataFrames.nrow(g0b) == 0 ? "PASS" :
            "FAIL ($(DataFrames.nrow(g0b)) mismatching fields)",
        )
        println(
            io,
            "- G0c, oracle gates pass at both commits and illumina sha equal: ",
            g0c_ok ? "PASS" : "FAIL",
        )
        println(
            io,
            "- Columns excluded from identity: ",
            join(TD4MBG_NON_OUTPUT_COLUMNS, ", "),
            "\n",
        )
        _table(io, oracle)
        DataFrames.nrow(g0_rows) == 0 || _table(io, first(g0_rows, 50))
        println(io, "## G2 outcome (x = |d_BASE| − |d_FIX|, H1: x > 0)\n")
        _table(io, DataFrames.select(g2_summary, DataFrames.Not(:method)))
        println(io, "Test: ", first(first(g2_results)).method, "; alpha $(TD2RXH_ALPHA).\n")
        println(io, "### Per seed, $(TD4MBG_G2_CELL)\n")
        _table(io, last(only(filter(r -> first(r).gated, g2_results))))
        println(io, "## Descriptive counters (not gates)\n")
        println(
            io,
            "demoted_whole_read_routes: NOT CAPTURED. The td2rxh harness does not export it (amendment 5).\n",
        )
        _table(io, counters)
        if drift !== nothing
            println(io, "## Context: earlier run vs BASE (drift; enters no gate)\n")
            _table(io, drift)
        end
    end
    return (
        verdict = verdict,
        label_problems = problems,
        complete = complete,
        g0a = g0a,
        g0b = g0b,
        g0c_ok = g0c_ok,
        oracle_base = oracle_base,
        oracle_fix = oracle_fix,
        g2 = g2,
        g2_results = g2_results,
        counters = counters,
        drift = drift,
    )
end

function td4mbg_analysis_main(args::Vector{String} = ARGS)::Int
    required = (
        "base-dir",
        "fix-dir",
        "base-oracle",
        "fix-oracle",
        "base-head",
        "fix-head",
        "manifest",
        "output-dir",
    )
    allowed = (required..., "context-dir")
    isodd(length(args)) && throw(ArgumentError("arguments must be --key value pairs"))
    opts = Dict{String,String}()
    for i = 1:2:length(args)
        key = args[i]
        startswith(key, "--") && key[3:end] in allowed ||
            throw(ArgumentError("unknown argument $(key)"))
        opts[key[3:end]] = args[i+1]
    end
    absent = [k for k in required if !haskey(opts, k)]
    isempty(absent) || throw(ArgumentError("missing --$(join(absent, ", --"))"))
    result = td4mbg_analysis(;
        base_dir = opts["base-dir"],
        fix_dir = opts["fix-dir"],
        base_oracle = opts["base-oracle"],
        fix_oracle = opts["fix-oracle"],
        base_head = opts["base-head"],
        fix_head = opts["fix-head"],
        manifest = opts["manifest"],
        output_dir = opts["output-dir"],
        context_dir = get(opts, "context-dir", nothing),
    )
    println(read(joinpath(opts["output-dir"], "report.md"), String))
    println("VERDICT: ", result.verdict)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(td4mbg_analysis_main())
end
