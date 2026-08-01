# OPT4 FROZEN-READ SKIP — REPRESENTATIVE-SCALE ACCURACY SIGN-OFF
# =============================================================================
# td-jbjd, Rule-of-5 pass 5 ("Verification"). Design doc:
# docs/design/2026-07-26-opt4-frozen-read-skip-design.md ("Accuracy sign-off
# methodology"). This is the load-bearing measurement that lets the maintainer
# decide whether `skip_frozen_reads` (opt-in, default off) is worth recommending
# to enable: it measures the accuracy/speed tradeoff on representative-scale
# data, complementing the toy-scale positive-control investigation in
# test/4_assembly/corrector_opt4_positive_control_test.jl (pass 4), whose
# header documents WHY toy scale could not exercise this question: Stage 0
# cheap correction (`_stage0_cheap_correct`, freeze-immune, runs unconditionally
# before the freeze gate) dominates toy-scale correction volume and the toy
# genome's dense graph trips the adaptive low-k decode gate
# (`decode_gate_density`, default 0.90) to zero out decode entirely. This
# script targets a scale/config (verified by `benchmarking/opt4_signoff_probe.jl`,
# a throwaway probe kept for reproducibility) where the per-read Viterbi decode
# demonstrably does real work (`decode_fraction_mean` > 0 and decode-attributable
# edits, i.e. `total_improvements - cheap_corrections_total`, > 0), so opt4's
# actual accuracy exposure -- confined to decode, per the pass-4 finding -- has
# something to measure against.
#
# DOES NOT reuse `rhizomorph_correction_accuracy_benchmark.jl`'s `run_accuracy_
# benchmark()` sweep driver (that harness has no `skip_frozen_reads` knob and is
# intentionally left unmutated -- design doc: "extend
# rhizomorph_correction_accuracy_benchmark.jl" was the METHODOLOGY, not a
# license to add opt4-specific plumbing to the shared harness). Instead this
# script `include`s it (safe: guarded by `abspath(PROGRAM_FILE) == @__FILE__`,
# so including only DEFINES functions) and reuses its
# `acquire_reference` / `simulate_substitution_reads` / `per_base_metrics`
# building blocks directly, adding only the `skip_frozen_reads` /
# `freeze_streak_threshold` / `freeze_across_rungs` passthrough that the shared
# harness's `correct_reads_scalable` does not expose.
#
# FOUR ARMS on the IDENTICAL read set (same simulated reads/errors, same seed):
#   1. exact             -- skip_frozen_reads=false (byte-identical-when-off path)
#   2. freeze_within      -- skip_frozen_reads=true, threshold=2, across=false (DEFAULT scope)
#   3. freeze_across       -- skip_frozen_reads=true, threshold=2, across=true (aggressive scope)
#   4. freeze_max_aggressive -- skip_frozen_reads=true, threshold=1, across=true
#      (the positive-control candidate: the most aggressive freeze setting this
#      design doc's Rule-of-5 pass 4 tried and could NOT make degrade at toy
#      scale -- here it is retried at a scale where decode does real work.)
#
# For each arm: per-base accuracy (recall/precision/over_correction_rate/
# correction_rate) via `per_base_metrics`, wall-clock runtime, and
# `frozen_reads_skipped` (from `result[:metadata][:corrector_errors]`).
#
# SCOPE NOTE: this script measures PER-BASE READ accuracy only (design doc
# check #1), not the assembly-level dnadiff/ANI check (design doc check #2) --
# deferred out of the ~40-minute compute budget for this pass; per-base
# accuracy against injected-error ground truth is the primary, most direct
# signal (it measures the corrector's own output, unmediated by a downstream
# assembler), and is explicitly what this pass was asked to reuse.
#
# Run directly (network required to download the reference once):
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/opt4_frozen_read_accuracy_signoff.jl
#
# Config override env vars (defaults chosen by benchmarking/opt4_signoff_probe.jl):
#   MYCELIA_OPT4_ACCESSION   default NC_001416 (Lambda phage, 48502 bp)
#   MYCELIA_OPT4_COVERAGE    default 8.0
#   MYCELIA_OPT4_ERR         default 0.02 (substitution rate)
#   MYCELIA_OPT4_READLEN     default 150
#   MYCELIA_OPT4_K           default 21
#   MYCELIA_OPT4_BATCH_SIZE  default 500
#   MYCELIA_OPT4_SEED        default 42

import Random
import FASTX
import Dates
import CSV
import DataFrames

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))

"""
Run the wired corrector (same knobs as `correct_reads_scalable`, plus the opt4
freeze passthrough) on `records`, at k, writing to a fresh `output_dir`. Returns
`(corrected_by_id, result)`.
"""
function correct_reads_with_freeze(records::Vector{FASTX.FASTQ.Record}, k::Int;
        skip_frozen::Bool, threshold::Int, across::Bool, batch_size::Int,
        output_dir::String)
    input_dir = mktempdir()
    temp_fastq = joinpath(input_dir, "in.fastq")
    open(FASTX.FASTQ.Writer, temp_fastq) do w
        for r in records
            write(w, r)
        end
    end
    max_k = max(k, 13)
    result = Mycelia.mycelia_iterative_assemble(
        temp_fastq;
        max_k = max_k, skip_solid = true, graph_mode = :doublestrand,
        n_k_rungs = 3, max_iterations_per_k = 2, batch_size = batch_size,
        hard_window = true, soft_em = true, cheap_correct = true,
        beam_width = nothing, verbose = false, enable_checkpointing = false,
        skip_frozen_reads = skip_frozen, freeze_streak_threshold = threshold,
        freeze_across_rungs = across, output_dir = output_dir)
    corrected_fastq = result[:metadata][:final_fastq_file]
    corrected_by_id = Dict{String, String}()
    open(FASTX.FASTQ.Reader, corrected_fastq) do rd
        for rec in rd
            corrected_by_id[FASTX.identifier(rec)] = FASTX.sequence(String, rec)
        end
    end
    return corrected_by_id, result
end

const ARMS = [
    (label = "exact", skip_frozen = false, threshold = 2, across = false),
    (label = "freeze_within", skip_frozen = true, threshold = 2, across = false),
    (label = "freeze_across", skip_frozen = true, threshold = 2, across = true),
    (
        label = "freeze_max_aggressive", skip_frozen = true, threshold = 1,
        across = true)
]

function run_signoff()
    accession = get(ENV, "MYCELIA_OPT4_ACCESSION", "NC_001416")
    coverage = parse(Float64, get(ENV, "MYCELIA_OPT4_COVERAGE", "8.0"))
    err = parse(Float64, get(ENV, "MYCELIA_OPT4_ERR", "0.02"))
    readlen = parse(Int, get(ENV, "MYCELIA_OPT4_READLEN", "150"))
    k = parse(Int, get(ENV, "MYCELIA_OPT4_K", "21"))
    batch_size = parse(Int, get(ENV, "MYCELIA_OPT4_BATCH_SIZE", "500"))
    seed = parse(Int, get(ENV, "MYCELIA_OPT4_SEED", "42"))
    assigned_q = 20

    println("=== OPT4 frozen-read-skip representative-scale accuracy sign-off ===")
    println("Start: $(Dates.now())")
    println(
        "Config: accession=$accession coverage=$(coverage)x err=$err readlen=$readlen k=$k batch_size=$batch_size seed=$seed")

    workdir = mktempdir(prefix = "opt4_signoff_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)

    refseq, ref_path,
    ref_label = acquire_reference(
        smoke = false, accession = accession, smoke_len = 0, seed = seed,
        workdir = workdir)
    glen = length(refseq)
    println("Reference: $ref_label ($glen bp)")

    rng = Random.MersenneTwister(seed)
    records, truth_by_id,
    observed_by_id,
    injected_total,
    sampled_bases = simulate_substitution_reads(
        refseq, readlen, coverage, err, rng; assigned_q = assigned_q)
    eff_cov = round(sampled_bases / glen; digits = 2)
    println(
        "Simulated $(length(records)) reads, effective coverage $(eff_cov)x, injected $injected_total substitutions")

    rows = DataFrames.DataFrame(
        arm = String[], skip_frozen_reads = Bool[], freeze_streak_threshold = Int[],
        freeze_across_rungs = Bool[], runtime_s = Float64[],
        frozen_reads_skipped = Int[], decode_fraction_mean = Float64[],
        decode_fraction_max = Float64[], total_improvements = Int[],
        cheap_corrections_total = Int[], decode_attributable_edits = Int[],
        reads_scored = Int[], injected_errors = Int[], total_edits = Int[],
        true_fixes = Int[], mis_fixes = Int[], over_corrections = Int[],
        recall = Float64[], precision = Float64[], over_correction_rate = Float64[],
        correction_rate = Float64[])

    arm_records = Dict{String, Any}()

    for arm in ARMS
        println("\n--- Arm: $(arm.label) ---")
        t0 = time()
        corrected_by_id,
        result = correct_reads_with_freeze(
            records, k; skip_frozen = arm.skip_frozen, threshold = arm.threshold,
            across = arm.across, batch_size = batch_size,
            output_dir = joinpath(workdir, arm.label))
        runtime = time() - t0
        md = result[:metadata]
        m = per_base_metrics(truth_by_id, observed_by_id, corrected_by_id)
        decode_edits = md[:total_improvements] - md[:cheap_corrections_total]
        frozen_skipped = md[:corrector_errors][:frozen_reads_skipped]

        push!(rows,
            (
                arm = arm.label, skip_frozen_reads = arm.skip_frozen,
                freeze_streak_threshold = arm.threshold,
                freeze_across_rungs = arm.across, runtime_s = round(runtime; digits = 3),
                frozen_reads_skipped = frozen_skipped,
                decode_fraction_mean = round(md[:decode_fraction_mean]; digits = 4),
                decode_fraction_max = round(md[:decode_fraction_max]; digits = 4),
                total_improvements = md[:total_improvements],
                cheap_corrections_total = md[:cheap_corrections_total],
                decode_attributable_edits = decode_edits, reads_scored = m.reads_scored,
                injected_errors = m.injected, total_edits = m.total_edits,
                true_fixes = m.tp, mis_fixes = m.mis_fixes, over_corrections = m.over,
                recall = round(m.recall; digits = 6), precision = round(m.precision; digits = 6),
                over_correction_rate = round(m.over_rate; digits = 6),
                correction_rate = round(m.correction_rate; digits = 6)))

        arm_records[arm.label] = (; corrected_by_id, result, m, runtime)
        println(
            "  runtime=$(round(runtime; digits=1))s frozen_reads_skipped=$frozen_skipped " *
            "decode_frac(mean/max)=$(round(md[:decode_fraction_mean]; digits=4))/$(round(md[:decode_fraction_max]; digits=4)) " *
            "decode_attributable_edits=$decode_edits")
        println(
            "  recall=$(round(m.recall; digits=4)) precision=$(round(m.precision; digits=4)) " *
            "over_rate=$(round(m.over_rate; digits=6)) correction_rate=$(round(m.correction_rate; digits=4)) " *
            "true_fixes=$(m.tp) mis_fixes=$(m.mis_fixes) over_corrections=$(m.over)")
    end

    # Byte-identity check between arms (belt-and-suspenders on top of the
    # per-base metrics, mirroring the pass-4 positive-control test's own check).
    exact_recs = sort(collect(arm_records["exact"].corrected_by_id))
    for label in ("freeze_within", "freeze_across", "freeze_max_aggressive")
        this_recs = sort(collect(arm_records[label].corrected_by_id))
        identical = exact_recs == this_recs
        n_diff = count(
            i -> exact_recs[i][2] != this_recs[i][2],
            1:min(length(exact_recs), length(this_recs)))
        println(
            "\n[identity] exact vs $label: byte-identical=$identical (reads differing: $n_diff / $(length(exact_recs)))")
    end

    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv_path = joinpath(results_dir, "opt4_frozen_read_skip_signoff_$(timestamp).csv")
    CSV.write(csv_path, rows)
    println("\nCSV written: $csv_path")
    println("\n=== Sign-off run complete: $(Dates.now()) ===")
    return (
        rows = rows, arm_records = arm_records, csv_path = csv_path, glen = glen,
        ref_label = ref_label, n_reads = length(records), eff_cov = eff_cov,
        injected_total = injected_total,
        config = (; accession, coverage, err, readlen, k, batch_size, seed))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_signoff()
end
