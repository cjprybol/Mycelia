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
#   scalable+iter   skip_solid=true,  cheap_correct=true, hard_window=true,  iter=5
#   noskip          skip_solid=FALSE, cheap_correct=true, hard_window=true,  iter=2   (decode reads scalable would skip)
#   noskip+nogate   skip_solid=FALSE, cheap_correct=true, hard_window=FALSE, iter=2   (also drop the decode gate)
#   exhaustive      skip_solid=FALSE, cheap_correct=false,hard_window=FALSE, iter=10, beam=exact  (max sensitivity)
#
# skip_solid is the primary recall lever: :scalable SKIPS "mostly-solid" reads, so
# errors on skipped reads are never corrected. Turning it off decodes them
# (recall up) at some precision / runtime cost — the core tradeoff under test.
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
#   MYCELIA_RPC_COVERAGE (default 15), MYCELIA_RPC_K (default 21),
#   MYCELIA_RPC_ASSIGNED_Q (default 20), MYCELIA_RPC_SEED (default 42),
#   MYCELIA_RPC_POINTS (comma list of point names; default all)

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))

# === Operating points (named knob presets) ==================================

# Three informative points span the frontier (a smoke run localized the recall
# limiter to the hard_window decode gate, NOT skip_solid): the wired default, the
# skip_solid-off variant (isolates that lever), and both gates off. A prior
# exact-beam "exhaustive" point is dropped — it is intractable on the dense low-k
# graph (times out) and adds nothing once noskip+nogate reaches recall ~1.0. The
# auto-bounded Viterbi beam (256) applies to every point, so these stay tractable.
const PR_OPERATING_POINTS = [
    (name = "scalable", skip_solid = true, cheap_correct = true, hard_window = true,
        soft_em = true, n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing),
    (name = "noskip", skip_solid = false, cheap_correct = true, hard_window = true,
        soft_em = true, n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing),
    (name = "noskip+nogate", skip_solid = false,
        cheap_correct = true, hard_window = false,
        soft_em = true, n_k_rungs = 3, max_iterations_per_k = 2, beam_width = nothing)
]

"""
Correct `records` with an explicit knob preset `pt` (a PR_OPERATING_POINTS entry),
returning `corrected_by_id`. Mirrors correct_reads_scalable but with the point's
knobs instead of the hardwired :scalable set.
"""
function correct_reads_at_point(records::Vector{FASTX.FASTQ.Record}, k::Int, pt)
    input_dir = mktempdir()
    output_dir = mktempdir()
    temp_fastq = joinpath(input_dir, "corrector_input.fastq")
    open(FASTX.FASTQ.Writer, temp_fastq) do w
        for r in records
            write(w, r)
        end
    end
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

function run_pr_curve()
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
    points = isempty(want) ? PR_OPERATING_POINTS :
             filter(p -> p.name in split(want, ","), PR_OPERATING_POINTS)

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
    rng = Random.MersenneTwister(seed)

    rows = DataFrames.DataFrame(
        error_rate = Float64[], point = String[], n_reads = Int[], reads_scored = Int[],
        injected = Int[], true_fixes = Int[], mis_fixes = Int[], over_corrections = Int[],
        recall = Float64[], precision = Float64[], over_correction_rate = Float64[],
        correction_rate = Float64[], runtime_s = Float64[])

    for err in errs
        # SAME reads for every operating point at this error rate (fair comparison).
        records, truth_by_id,
        observed_by_id,
        injected_total,
        _sb = simulate_substitution_reads(
            refseq, readlens[1], coverage, err, rng; assigned_q = assigned_q)
        println("\n[err=$err] $(length(records)) reads, injected $injected_total substitutions")
        for pt in points
            t0 = time()
            local corrected
            ok = true
            try
                corrected = correct_reads_at_point(records, k, pt)
            catch e
                @warn "point failed" point=pt.name err exception=(e, catch_backtrace())
                ok = false
                corrected = Dict{String, String}()
            end
            rt = round(time() - t0; digits = 2)
            m = per_base_metrics(truth_by_id, observed_by_id, corrected)
            push!(rows,
                (error_rate = err, point = pt.name, n_reads = length(records),
                    reads_scored = m.reads_scored, injected = m.injected, true_fixes = m.tp,
                    mis_fixes = m.mis_fixes, over_corrections = m.over, recall = m.recall,
                    precision = m.precision, over_correction_rate = m.over_rate,
                    correction_rate = m.correction_rate, runtime_s = rt))
            println("  $(rpad(pt.name, 15)) recall=$(rpad(round(m.recall; digits=4), 8)) " *
                    "precision=$(rpad(round(m.precision; digits=4), 8)) " *
                    "over_rate=$(rpad(round(m.over_rate; digits=6), 10)) " *
                    "mis=$(rpad(m.mis_fixes, 5)) corr_rate=$(rpad(round(m.correction_rate; digits=4), 8)) $(rt)s $(ok ? "" : "[FAILED]")")
        end
    end

    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv = joinpath(results_dir, "rhizomorph_correction_pr_curve_$(ts).csv")
    CSV.write(csv, rows)
    println("\nCSV: $csv")
    println("\nRead the frontier per error rate: a point with HIGHER recall at ~equal precision")
    println("+ ~zero over_correction is a strict improvement (the corrector was over-conservative).")
    println("A point that trades precision/over-correction for recall exposes the real tradeoff.")
    println("=== done: $(Dates.now()) ===")
    return (; csv, rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pr_curve()
end
