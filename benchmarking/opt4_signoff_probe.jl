# OPT4 PASS-5 SIGN-OFF — CONFIG PROBE (td-jbjd, pass 5, throwaway script)
# =============================================================================
# Establishes whether a given (accession, coverage, error_rate, k, batch_size)
# configuration makes the per-read Viterbi decode do REAL work (as opposed to
# the pass-4 finding: at toy scale Stage 0 dominates and the adaptive low-k
# decode gate disables decode entirely). Evidence checked per candidate:
#   - result[:metadata][:decode_fraction_mean] > 0  (decode was not gated to 0
#     on every pass -- see :2085 `effective_decode_gate_density` /
#     iterative-assembly.jl ~4831 `adaptive_gated`)
#   - total_improvements - cheap_corrections_total > 0 (decode, not just Stage
#     0, made edits)
# Run directly (network required, downloads the reference once):
#   LD_LIBRARY_PATH='' julia --project=. benchmarking/opt4_signoff_probe.jl

import Random
import FASTX

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))

function probe_config(; accession, coverage, err, readlen, k, batch_size, seed = 42)
    workdir = mktempdir(prefix = "opt4_probe_")
    refseq, ref_path,
    ref_label = acquire_reference(
        smoke = false, accession = accession, smoke_len = 0, seed = seed, workdir = workdir)
    glen = length(refseq)
    rng = Random.MersenneTwister(seed)
    records, truth_by_id,
    observed_by_id,
    injected_total,
    sampled_bases = simulate_substitution_reads(
        refseq, readlen, coverage, err, rng; assigned_q = 20)
    println("  genome=$ref_label ($glen bp) reads=$(length(records)) injected=$injected_total sampled_bases=$sampled_bases")

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
        n_k_rungs = 3, max_iterations_per_k = 2, batch_size = batch_size,
        hard_window = true, soft_em = true, cheap_correct = true,
        beam_width = nothing, verbose = false, enable_checkpointing = false,
        output_dir = joinpath(workdir, "run"))
    runtime = time() - t0
    md = result[:metadata]
    decode_edits = md[:total_improvements] - md[:cheap_corrections_total]
    println("  runtime=$(round(runtime; digits=1))s decode_fraction(mean/max)=" *
            "$(round(md[:decode_fraction_mean]; digits=4))/$(round(md[:decode_fraction_max]; digits=4)) " *
            "total_improvements=$(md[:total_improvements]) cheap_total=$(md[:cheap_corrections_total]) " *
            "decode_attributable_edits=$decode_edits final_k=$(md[:final_k])")
    return (; runtime, md, glen, n_reads = length(records), injected_total)
end

candidates = [
    (accession = "NC_001416", coverage = 8.0, err = 0.02,
        readlen = 150, k = 21, batch_size = 500),
    (accession = "NC_001416", coverage = 12.0, err = 0.03,
        readlen = 150, k = 31, batch_size = 500)
]

for (i, c) in enumerate(candidates)
    println("\n=== Candidate $i: $c ===")
    try
        probe_config(; c...)
    catch e
        println("  FAILED: $e")
    end
end
