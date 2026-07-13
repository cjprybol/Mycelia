# Rhizomorph Correction PRECISION-RECALL Frontier
# ===============================================
#
# The per-base accuracy benchmark (rhizomorph_correction_accuracy_benchmark.jl)
# measures ONE operating point — the wired :scalable corrector — and finds
# precision=1.0 / over-correction=0 but conservative recall (0.99 -> 0.06 as error
# density rises). The recall shortfall is an EDIT-AGGRESSIVENESS property of the
# corrector's decision knobs, NOT the calibration map. This harness maps the
# precision-recall FRONTIER by re-running the SAME corrector at several
# aggressiveness operating points, so we can see whether recall can be bought back
# and at what precision / over-correction cost.
#
# Operating points (all exposed knobs — NO src changes):
#   scalable        skip_solid=true,  cheap_correct=true, hard_window=true,  iter=2   (the wired default)
#   noskip          skip_solid=FALSE, cheap_correct=true, hard_window=true,  iter=2   (isolates the skip_solid lever)
#   noskip+nogate   skip_solid=FALSE, cheap_correct=true, hard_window=FALSE, iter=2   (also drops the decode gate)
#
# A smoke run localized the recall lever to the hard_window DECODE GATE, not
# skip_solid: turning skip_solid off alone leaves recall unchanged; dropping the
# gate too lifts recall from ~0.24 to ~1.0 (at a precision cost that grows with
# error rate). Exact-beam / more-iteration points are dropped — the former is
# intractable on the dense low-k graph, the latter is a no-op once converged.
#
# Reuses simulate_substitution_reads + per_base_metrics from the merged accuracy
# harness. Reports (precision, recall, over_correction, mis_fixes, correction_rate,
# runtime) per (error_rate x operating_point).
#
# Usage:
#   MYCELIA_RPC_SMOKE=true LD_LIBRARY_PATH='' \
#     julia --project=. benchmarking/rhizomorph_correction_pr_curve.jl
#
# Env (RPC prefix): MYCELIA_RPC_SMOKE, MYCELIA_RPC_SMOKE_GENOME_LEN,
#   MYCELIA_RPC_ERR (default 0.05,0.10), MYCELIA_RPC_READLEN (default 150),
#   MYCELIA_RPC_COVERAGE (default 30; smoke default 15), MYCELIA_RPC_K (default 21),
#   MYCELIA_RPC_ASSIGNED_Q (default 20), MYCELIA_RPC_SEED (default 42),
#   MYCELIA_RPC_POINTS (comma list of point names; default all),
#   MYCELIA_RPC_GAP_THRESHOLDS (default 0,0.5,1,2,4),
#   MYCELIA_RPC_PROB_THRESHOLDS (dense near zero, where recall changes fastest)

import Dates
import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))
include(joinpath(@__DIR__, "viterbi_gap_calibration.jl"))

# === Operating points (named knob presets) ==================================

# Three informative points span the frontier (a smoke run localized the recall
# limiter to the hard_window decode gate, NOT skip_solid): the wired default, the
# skip_solid-off variant (isolates that lever), and both gates off. A prior
# exact-beam "exhaustive" point is dropped — it is intractable on the dense low-k
# graph (times out) and adds nothing once noskip+nogate reaches recall ~1.0. The
# auto-bounded Viterbi beam (256) applies to every point, so these stay tractable.
const PR_OPERATING_POINTS = [
    (name = "scalable", skip_solid = true, cheap_correct = true, hard_window = true,
        soft_em = true, n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing,
        calibrated_gap_threshold = nothing, calibrated_probability_model = nothing,
        calibrated_probability_threshold = 0.5),
    (name = "noskip", skip_solid = false, cheap_correct = true, hard_window = true,
        soft_em = true, n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing,
        calibrated_gap_threshold = nothing, calibrated_probability_model = nothing,
        calibrated_probability_threshold = 0.5),
    (name = "noskip+nogate", skip_solid = false,
        cheap_correct = true, hard_window = false,
        soft_em = true, n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing,
        calibrated_gap_threshold = nothing, calibrated_probability_model = nothing,
        calibrated_probability_threshold = 0.5)
]

function calibrated_gap_points()::Vector
    thresholds = _parse_float_list(
        get(ENV, "MYCELIA_RPC_GAP_THRESHOLDS", "0.0,0.5,1.0,2.0,4.0"), Float64[])
    return [(
                name = "calibrated-gap-$(threshold)", skip_solid = false,
                cheap_correct = true, hard_window = false, soft_em = true,
                n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing,
                calibrated_gap_threshold = threshold, calibrated_probability_model = nothing,
                calibrated_probability_threshold = 0.5
            ) for threshold in thresholds]
end

function calibrated_probability_points()::Vector
    thresholds = _parse_float_list(
        get(ENV, "MYCELIA_RPC_PROB_THRESHOLDS",
            "0.0,0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.5,0.75,0.9"),
        Float64[])
    all(threshold -> 0.0 <= threshold <= 1.0, thresholds) ||
        throw(ArgumentError("MYCELIA_RPC_PROB_THRESHOLDS values must be in [0, 1]"))
    return [(
                name = "calibrated-prob-$(threshold)", skip_solid = false,
                cheap_correct = true, hard_window = false, soft_em = true,
                n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing,
                calibrated_gap_threshold = nothing, calibrated_probability_model = nothing,
                calibrated_probability_threshold = threshold
            ) for threshold in thresholds]
end

function persist_runtime_probability_model(
        model::Mycelia.CorrectionConfidenceModel,
        results_dir::AbstractString)::String
    mkpath(results_dir)
    model_path = joinpath(results_dir, "runtime_probability_model.csv")
    coefficients = (
        intercept = model.intercept,
        raw_gap = model.gap_weight,
        min_kmer_support = model.min_kmer_support_weight,
        all_stage0_solid = model.stage0_solid_weight,
        competing_branch_support_ratio = model.branch_support_ratio_weight,
        collapsed_frontier = model.collapsed_frontier_weight
    )
    open(model_path, "w") do io
        println(io, "term,coefficient")
        for term in keys(coefficients)
            println(io, "$(term),$(getproperty(coefficients, term))")
        end
    end
    return model_path
end

"""
Fit the probability gate on an independent synthetic cohort and convert the
benchmark logistic coefficients to the runtime model's fixed feature order.
The calibration cohort uses a disjoint seed and never sees frontier truth.
"""
function fit_frontier_probability_model(errs::Vector{Float64}, readlen::Int,
        coverage::Float64, k::Int, seed::Int)::NamedTuple
    results_root = joinpath(@__DIR__, "results")
    mkpath(results_root)
    calibration_results_dir = mktempdir(
        results_root; prefix = "td21eg_calibration_", cleanup = false)
    artifact = run_gap_calibration(
        k = k,
        genome_length = max(1200, 10 * readlen),
        readlen = readlen,
        coverage = max(1, round(Int, coverage)),
        error_rates = errs,
        seed = seed + 10_000,
        results_dir = calibration_results_dir,
        return_artifact = true)
    expected_features = (
        :raw_gap,
        :min_kmer_support,
        :all_stage0_solid,
        :competing_branch_support_ratio,
        :collapsed_frontier
    )
    artifact.feature_names == expected_features ||
        throw(ArgumentError("calibration/runtime feature-order mismatch: " *
                            "expected $(expected_features), got $(artifact.feature_names)"))
    model = artifact.multifeature_model
    length(model.b) == length(expected_features) ||
        throw(ArgumentError("multi-feature logistic has wrong coefficient count"))
    runtime_model = Mycelia.CorrectionConfidenceModel(model.a, model.b...)
    println("Calibration cohort: multi-feature model fitted with grouped holdout; " *
            "seed=$(seed + 10_000), n_train=$(length(artifact.training.labels)), " *
            "n_heldout=$(length(artifact.heldout.labels)), " *
            "artifacts=$(calibration_results_dir)")
    println("Calibration coefficients: intercept=$(runtime_model.intercept), " *
            "gap=$(runtime_model.gap_weight), " *
            "min_support=$(runtime_model.min_kmer_support_weight), " *
            "stage0_solid=$(runtime_model.stage0_solid_weight), " *
            "competing_branch=$(runtime_model.branch_support_ratio_weight), " *
            "collapsed_frontier=$(runtime_model.collapsed_frontier_weight)")
    return (model = runtime_model, results_dir = calibration_results_dir)
end

"""
Correct `records` with an explicit knob preset `pt` (a PR_OPERATING_POINTS entry),
returning `corrected_by_id`. Mirrors correct_reads_scalable but with the point's
knobs instead of the hardwired :scalable set.
"""
function correct_reads_at_point(records::Vector{FASTX.FASTQ.Record}, k::Int,
        pt::NamedTuple)::Dict{String, String}
    input_dir = mktempdir()
    output_dir = mktempdir()
    temp_fastq = joinpath(input_dir, "corrector_input.fastq")
    open(FASTX.FASTQ.Writer, temp_fastq) do w
        for r in records
            write(w, r)
        end
    end
    probability_model = hasproperty(pt, :calibrated_probability_model) ?
                        pt.calibrated_probability_model : nothing
    probability_threshold = hasproperty(pt, :calibrated_probability_threshold) ?
                            pt.calibrated_probability_threshold : 0.5
    result = Mycelia.mycelia_iterative_assemble(
        temp_fastq;
        max_k = max(k, 13),
        skip_solid = pt.skip_solid,
        graph_mode = :doublestrand,
        n_k_rungs = pt.n_k_rungs,
        max_iterations_per_k = pt.max_iterations_per_k,
        hard_window = pt.hard_window,
        soft_em = pt.soft_em,
        cheap_correct = pt.cheap_correct,
        beam_width = pt.beam_width,
        calibrated_gap_threshold = pt.calibrated_gap_threshold,
        calibrated_probability_model = probability_model,
        calibrated_probability_threshold = probability_threshold,
        verbose = false,
        enable_checkpointing = false,
        output_dir = output_dir
    )
    corrected_fastq = get(result[:metadata], :final_fastq_file, nothing)
    (corrected_fastq === nothing || !isfile(corrected_fastq)) &&
        error("corrector produced no :final_fastq_file at point $(pt.name)")
    corrected_by_id = Dict{String, String}()
    open(FASTX.FASTQ.Reader, corrected_fastq) do reader
        for rec in reader
            corrected_by_id[FASTX.identifier(rec)] = String(FASTX.sequence(BioSequences.LongDNA{4}, rec))
        end
    end
    return corrected_by_id
end

function run_pr_curve(;
        probability_model::Union{Nothing, Mycelia.CorrectionConfidenceModel} = nothing)::NamedTuple
    smoke = _truthy(get(ENV, "MYCELIA_RPC_SMOKE", "false"))
    accession = get(ENV, "MYCELIA_RPC_ACCESSION", "NC_001416")
    smoke_len = parse(Int, get(ENV, "MYCELIA_RPC_SMOKE_GENOME_LEN", "2000"))
    errs = _parse_float_list(get(ENV, "MYCELIA_RPC_ERR", ""), [0.05, 0.10])
    readlens = _parse_int_list(get(ENV, "MYCELIA_RPC_READLEN", ""), [150])
    coverage = parse(Float64, get(ENV, "MYCELIA_RPC_COVERAGE", smoke ? "15" : "30"))
    k = parse(Int, get(ENV, "MYCELIA_RPC_K", "21"))
    assigned_q = parse(Int, get(ENV, "MYCELIA_RPC_ASSIGNED_Q", "20"))
    seed = parse(Int, get(ENV, "MYCELIA_RPC_SEED", "42"))
    want = strip(get(ENV, "MYCELIA_RPC_POINTS", ""))
    all_points = vcat(
        PR_OPERATING_POINTS, calibrated_gap_points(), calibrated_probability_points())
    points = isempty(want) ? all_points :
             filter(p -> p.name in split(want, ","), all_points)
    calibration_artifact = ""
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)

    println("=== Rhizomorph Correction Precision-Recall Frontier ===")
    println("Start: $(Dates.now())   error rates: $errs   coverage: $(coverage)x   k: $k")
    println("Operating points: $(join([p.name for p in points], ", "))")

    workdir = mktempdir(prefix = "rpc_")
    refseq, _rp,
    ref_label = acquire_reference(
        smoke = smoke, accession = accession, smoke_len = smoke_len, seed = seed, workdir = workdir)
    glen = length(refseq)
    println("Reference: $ref_label ($glen bp)")
    k = min(k, minimum(readlens), glen)
    if any(startswith(point.name, "calibrated-prob-") for point in points)
        effective_probability_model = if probability_model === nothing
            fit = fit_frontier_probability_model(
                errs, readlens[1], coverage, k, seed)
            calibration_artifact = fit.results_dir
            fit.model
        else
            supplied_dir = mktempdir(
                results_dir; prefix = "td21eg_caller_model_", cleanup = false)
            persist_runtime_probability_model(probability_model, supplied_dir)
            calibration_artifact = supplied_dir
            probability_model
        end
        points = [startswith(point.name, "calibrated-prob-") ?
                  merge(point,
                      (calibrated_probability_model = effective_probability_model,)) :
                  point for point in points]
    end
    rng = Random.MersenneTwister(seed)

    rows = DataFrames.DataFrame(
        error_rate = Float64[], point = String[], ok = Bool[], n_reads = Int[],
        reads_scored = Int[], injected = Int[], true_fixes = Int[], mis_fixes = Int[],
        over_corrections = Int[], recall = Float64[], precision = Float64[],
        over_correction_rate = Float64[], correction_rate = Float64[], runtime_s = Float64[],
        calibration_artifact = String[])

    for err in errs
        # SAME reads for every operating point at this error rate (fair comparison).
        records, truth_by_id,
        observed_by_id,
        injected_total,
        _sb = simulate_substitution_reads(
            refseq, readlens[1], coverage, err, rng; assigned_q = assigned_q)
        println("\n[err=$err] $(length(records)) reads, injected $injected_total substitutions")
        for pt in points
            # Re-seed the GLOBAL RNG before each point: mycelia_iterative_assemble's
            # statistical-resampling arm draws from it, so without this each point's
            # stochastic behavior would depend on how much prior points drew —
            # confounding the very knob-comparison this harness exists to make.
            Random.seed!(seed)
            t0 = time()
            local corrected
            ok = true
            try
                corrected = correct_reads_at_point(records, k, pt)
            catch e
                e isa InterruptException && rethrow()
                @warn "point failed" point=pt.name err exception=(e, catch_backtrace())
                ok = false
                corrected = Dict{String, String}()
            end
            rt = round(time() - t0; digits = 2)
            m = per_base_metrics(truth_by_id, observed_by_id, corrected)
            push!(rows,
                (error_rate = err, point = pt.name, ok = ok, n_reads = length(records),
                    reads_scored = m.reads_scored, injected = m.injected, true_fixes = m.tp,
                    mis_fixes = m.mis_fixes, over_corrections = m.over, recall = m.recall,
                    precision = m.precision, over_correction_rate = m.over_rate,
                    correction_rate = m.correction_rate, runtime_s = rt,
                    calibration_artifact = startswith(pt.name, "calibrated-prob-") ?
                                           calibration_artifact : ""))
            println("  $(rpad(pt.name, 15)) recall=$(rpad(round(m.recall; digits=4), 8)) " *
                    "precision=$(rpad(round(m.precision; digits=4), 8)) " *
                    "over_rate=$(rpad(round(m.over_rate; digits=6), 10)) " *
                    "mis=$(rpad(m.mis_fixes, 5)) corr_rate=$(rpad(round(m.correction_rate; digits=4), 8)) $(rt)s $(ok ? "" : "[FAILED]")")
        end
    end

    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv = joinpath(results_dir, "rhizomorph_correction_pr_curve_$(ts).csv")
    CSV.write(csv, rows)
    println("\nCSV: $csv")

    # Preserve the legacy raw-gap verdict for comparison.
    if 0.10 in errs && any(startswith(p.name, "calibrated-gap-") for p in points)
        dom = assess_frontier_dominance(rows; err = 0.10)
        if hasproperty(dom, :reason)
            println("\n[dominance] $(dom.reason)")
        elseif dom.dominates
            println("\n[dominance] err=0.10 CALIBRATED GATE DOMINATES noskip+nogate: " *
                    "$(dom.best.point) recall=$(round(dom.best.recall; digits = 4)) " *
                    "precision=$(round(dom.best.precision; digits = 4)) " *
                    "over_rate=$(round(dom.best.over_rate; digits = 6)) " *
                    "(baseline recall=$(round(dom.baseline.recall; digits = 4)) " *
                    "precision=$(round(dom.baseline.precision; digits = 4)) " *
                    "over_rate=$(round(dom.baseline.over_rate; digits = 6)))")
        else
            println("\n[dominance] err=0.10 NO calibrated point dominates noskip+nogate " *
                    "(baseline recall=$(round(dom.baseline.recall; digits = 4)) " *
                    "precision=$(round(dom.baseline.precision; digits = 4)) " *
                    "over_rate=$(round(dom.baseline.over_rate; digits = 6)))")
        end
    end
    # The multi-feature probability family is the td-21eg decision surface. Report
    # both requested error rates separately; only err=0.10 is the completion gate.
    if any(startswith(point.name, "calibrated-prob-") for point in points)
        for err in (0.05, 0.10)
            err in errs || continue
            dom = assess_frontier_dominance(
                rows; err = err, candidate_prefix = "calibrated-prob-")
            if hasproperty(dom, :reason)
                println("\n[probability dominance] $(dom.reason)")
            elseif dom.dominates
                println("\n[probability dominance] err=$(err) CALIBRATED PROBABILITY " *
                        "DOMINATES noskip+nogate: $(dom.best.point) " *
                        "recall=$(round(dom.best.recall; digits = 4)) " *
                        "precision=$(round(dom.best.precision; digits = 4)) " *
                        "over_rate=$(round(dom.best.over_rate; digits = 6)) " *
                        "(baseline recall=$(round(dom.baseline.recall; digits = 4)) " *
                        "precision=$(round(dom.baseline.precision; digits = 4)) " *
                        "over_rate=$(round(dom.baseline.over_rate; digits = 6)))")
            else
                println("\n[probability dominance] err=$(err) NO calibrated-prob point " *
                        "strictly dominates noskip+nogate (baseline recall=" *
                        "$(round(dom.baseline.recall; digits = 4)) precision=" *
                        "$(round(dom.baseline.precision; digits = 4)) over_rate=" *
                        "$(round(dom.baseline.over_rate; digits = 6)))")
            end
        end
    end
    println("\nRead the frontier per error rate: a point with HIGHER recall at ~equal precision")
    println("+ ~zero over_correction is a strict improvement (the corrector was over-conservative).")
    println("A point that trades precision/over-correction for recall exposes the real tradeoff.")
    println("=== done: $(Dates.now()) ===")
    return (; csv, rows)
end

"""
    assess_frontier_dominance(rows; err=0.10, recall_eps=0.02,
        candidate_prefix="calibrated-gap-") -> NamedTuple

Programmatic frontier-dominance check (previously eyeballed from the CSV). At
error rate `err`, decide whether any operating point selected by `candidate_prefix` DOMINATES
the `noskip+nogate` baseline. Dominance requires a STRICT up-and-right move, not a
tie: recall retained (>= baseline recall − `recall_eps`), over-correction
**strictly** below baseline, and precision not below baseline. The strict
over-correction criterion is load-bearing — each family's zero-threshold point
is by definition the ungated baseline (threshold 0 reverts nothing), so it TIES
every metric; a non-strict predicate would report it as "dominating" itself. Returns
`(dominates, baseline, best, candidates)`. A non-finite baseline or non-finite
candidate metrics (zero-denominator rows) are treated as non-dominating. A row
is eligible only when its run succeeded (`ok`) and scored every simulated read
(`reads_scored == n_reads`); partial or failed runs cannot establish dominance.
"""
function assess_frontier_dominance(rows::DataFrames.DataFrame;
        err::Float64 = 0.10, recall_eps::Float64 = 0.02,
        candidate_prefix::AbstractString = "calibrated-gap-")::NamedTuple
    at = rows[[isapprox(r, err; atol = 1e-9) for r in rows.error_rate], :]
    base_idx = findfirst(==("noskip+nogate"), at.point)
    base_idx === nothing && return (dominates = false,
        reason = "no noskip+nogate baseline at err=$(err)",
        baseline = nothing, best = nothing, candidates = NamedTuple[])
    base = at[base_idx, :]
    if !(base.ok && base.reads_scored == base.n_reads)
        return (dominates = false,
            reason = "failed or partial noskip+nogate baseline at err=$(err): " *
                     "ok=$(base.ok), reads_scored=$(base.reads_scored)/$(base.n_reads)",
            baseline = (recall = base.recall, precision = base.precision,
                over_rate = base.over_correction_rate),
            best = nothing, candidates = NamedTuple[])
    end
    if !(isfinite(base.recall) && isfinite(base.precision) &&
         isfinite(base.over_correction_rate))
        return (dominates = false,
            reason = "degenerate (non-finite) noskip+nogate baseline at err=$(err)",
            baseline = (recall = base.recall, precision = base.precision,
                over_rate = base.over_correction_rate),
            best = nothing, candidates = NamedTuple[])
    end
    winners = NamedTuple[]
    for row in DataFrames.eachrow(at)
        startswith(row.point, candidate_prefix) || continue
        (row.ok && row.reads_scored == row.n_reads) || continue
        (isfinite(row.recall) && isfinite(row.precision) &&
         isfinite(row.over_correction_rate)) || continue
        # STRICT improvement: retain recall, strictly reduce over-correction, keep
        # precision. The strict `<` on over-correction excludes the gap-0.0 tie.
        if row.recall >= base.recall - recall_eps &&
           row.over_correction_rate < base.over_correction_rate &&
           row.precision >= base.precision
            push!(winners,
                (point = row.point, recall = row.recall,
                    precision = row.precision, over_rate = row.over_correction_rate))
        end
    end
    best = isempty(winners) ? nothing :
           sort(winners; by = w -> (w.over_rate, -w.recall))[1]
    return (dominates = !isempty(winners),
        baseline = (recall = base.recall, precision = base.precision,
            over_rate = base.over_correction_rate),
        best = best, candidates = winners)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pr_curve()
end
