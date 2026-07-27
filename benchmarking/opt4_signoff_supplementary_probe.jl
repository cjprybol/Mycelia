# OPT4 PASS-5 SIGN-OFF — SUPPLEMENTARY PROBE (td-jbjd, throwaway script)
# =============================================================================
# The main sign-off (opt4_frozen_read_accuracy_signoff.jl) found
# freeze_within (skip_frozen_reads=true, threshold=2, across=false, the DEFAULT
# scope) skipped ZERO reads at the production `:scalable`-tier config
# (n_k_rungs=3, max_iterations_per_k=2). Mechanism: within-rung resets the
# per-read streak to 0 at every k-change, and a threshold=2 freeze needs 2
# CONSECUTIVE no-improvement passes at the SAME k before it can skip a 3rd
# pass -- but max_iterations_per_k=2 means there is no 3rd pass at any k, so
# the freeze condition can be satisfied only on the LAST iteration of a rung,
# by which point there is nothing left to skip. This probe raises
# max_iterations_per_k to check whether within-rung freezing CAN engage when
# given room, and if so what its accuracy/speed profile looks like -- purely
# to explain, not to change, the main sign-off's config choice (which
# deliberately matches the real `:scalable` production knobs).
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/opt4_signoff_supplementary_probe.jl

import Random
import FASTX

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))

function run_arm(records, truth_by_id, observed_by_id, k, workdir, label;
        skip_frozen, threshold, across, max_iterations_per_k, batch_size = 500)
    input_dir = mktempdir()
    temp_fastq = joinpath(input_dir, "in.fastq")
    open(FASTX.FASTQ.Writer, temp_fastq) do w
        for r in records
            write(w, r)
        end
    end
    max_k = max(k, 13)
    t0 = time()
    result = Mycelia.mycelia_iterative_assemble(
        temp_fastq;
        max_k = max_k, skip_solid = true, graph_mode = :doublestrand,
        n_k_rungs = 3, max_iterations_per_k = max_iterations_per_k,
        batch_size = batch_size, hard_window = true, soft_em = true,
        cheap_correct = true, beam_width = nothing, verbose = false,
        enable_checkpointing = false, skip_frozen_reads = skip_frozen,
        freeze_streak_threshold = threshold, freeze_across_rungs = across,
        output_dir = joinpath(workdir, label))
    runtime = time() - t0
    corrected_by_id = Dict{String, String}()
    open(FASTX.FASTQ.Reader, result[:metadata][:final_fastq_file]) do rd
        for rec in rd
            corrected_by_id[FASTX.identifier(rec)] = FASTX.sequence(String, rec)
        end
    end
    m = per_base_metrics(truth_by_id, observed_by_id, corrected_by_id)
    md = result[:metadata]
    println(
        "[$label] miPk=$max_iterations_per_k runtime=$(round(runtime; digits=1))s " *
        "frozen_reads_skipped=$(md[:corrector_errors][:frozen_reads_skipped]) " *
        "recall=$(round(m.recall; digits=4)) precision=$(round(m.precision; digits=4)) " *
        "over_rate=$(round(m.over_rate; digits=6)) true_fixes=$(m.tp)")
    return (; runtime, m, md, corrected_by_id)
end

function main()
    accession = "NC_001416"
    coverage = 8.0
    err = 0.02
    readlen = 150
    k = 21
    seed = 42
    workdir = mktempdir(prefix = "opt4_suppl_")
    refseq, ref_path,
    ref_label = acquire_reference(
        smoke = false, accession = accession, smoke_len = 0, seed = seed,
        workdir = workdir)
    rng = Random.MersenneTwister(seed)
    records, truth_by_id,
    observed_by_id,
    injected_total = simulate_substitution_reads(
        refseq, readlen, coverage, err, rng; assigned_q = 20)
    println("reads=$(length(records)) injected=$injected_total")

    for mipk in (4,)
        println("\n=== max_iterations_per_k=$mipk ===")
        run_arm(records, truth_by_id, observed_by_id, k, workdir, "exact_mipk$mipk";
            skip_frozen = false, threshold = 2, across = false,
            max_iterations_per_k = mipk)
        run_arm(
            records, truth_by_id, observed_by_id, k, workdir, "within_mipk$mipk";
            skip_frozen = true, threshold = 2, across = false,
            max_iterations_per_k = mipk)
    end
end

main()
