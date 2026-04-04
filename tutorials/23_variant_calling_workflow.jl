# # Tutorial 23: Variant Calling Workflow
#
# This tutorial demonstrates a compact end-to-end small-variant workflow in
# Mycelia using:
# - `minimap2` for read alignment
# - `GATK HaplotypeCaller` for short-read variant calling
# - `Clair3` for an optional deep-learning caller walkthrough
# - `RTG vcfeval` for optional truth-set benchmarking
#
# The tutorial builds a tiny synthetic truth set so the control flow is easy to
# inspect and adapt for real projects.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand how to:
# - create a toy truth VCF from a reference sequence
# - simulate a perturbed haplotype and align reads back to the original reference
# - call variants with `Mycelia.run_gatk_haplotypecaller`
# - normalize and benchmark a callset with `normalize_vcf` and `evaluate_variant_calling_accuracy`
# - run the corresponding Clair3 workflow for Illumina-style data
# - update a reference FASTA with a VCF consensus
#
# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/23_variant_calling_workflow.jl", "tutorials/notebooks", execute=false)'
# ```
#
# To execute the external-tool sections, opt in with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. tutorials/23_variant_calling_workflow.jl
# ```
#
# To also run RTG vcfeval benchmarking:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true MYCELIA_RUN_VCFEVAL_TUTORIAL=true julia --project=. tutorials/23_variant_calling_workflow.jl
# ```
#
# To also run Clair3:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true MYCELIA_RUN_CLAIR3_TUTORIAL=true julia --project=. tutorials/23_variant_calling_workflow.jl
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import DataFrames
import Mycelia
import StableRNGs
import Test

const RUN_EXTERNAL = lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
const RUN_CLAIR3 = RUN_EXTERNAL &&
                   lowercase(get(ENV, "MYCELIA_RUN_CLAIR3_TUTORIAL", "false")) == "true"
const RUN_VCFEVAL = RUN_EXTERNAL &&
                    lowercase(get(ENV, "MYCELIA_RUN_VCFEVAL_TUTORIAL", "false")) == "true"

println("=== Variant Calling Workflow Tutorial ===")
println("External wrappers enabled: ", RUN_EXTERNAL)
println("Clair3 demo enabled: ", RUN_CLAIR3)
println("RTG vcfeval enabled: ", RUN_VCFEVAL)

# ## Part 1: Create a tiny synthetic truth set
#
# We start with a random reference, simulate a handful of variants, and write
# them as a VCF. This part is pure Julia and can run without external tools.

rng = StableRNGs.StableRNG(42)
workdir = mktempdir()
println("Workspace: ", workdir)

reference_record = Mycelia.random_fasta_record(
    moltype = :DNA,
    seed = rand(rng, 0:typemax(Int)),
    L = 20_000
)
reference_fasta = joinpath(workdir, "reference.fa")
Mycelia.write_fasta(outfile = reference_fasta, records = [reference_record])

truth_vcf = reference_fasta * ".vcf"
truth_table = Mycelia.simulate_variants(reference_record; n_variants = 12)
Mycelia.write_vcf_table(vcf_file = truth_vcf, vcf_table = truth_table, fasta_file = reference_fasta)

Test.@test isfile(reference_fasta)
Test.@test isfile(truth_vcf)
Test.@test DataFrames.nrow(truth_table) > 0

println("Reference FASTA: ", reference_fasta)
println("Truth VCF: ", truth_vcf)
println("Simulated truth variants: ", DataFrames.nrow(truth_table))

# ## Part 2: Materialize a mutated haplotype FASTA
#
# `update_fasta_with_vcf` applies the truth VCF to the original reference so we
# can simulate reads from the altered haplotype and map them back to the
# unmodified reference.

mutant_fasta = joinpath(workdir, "truth_haplotype.fa")

if RUN_EXTERNAL
    Mycelia.add_bioconda_env("samtools")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools faidx $(reference_fasta)`)

    Mycelia.update_fasta_with_vcf(
        in_fasta = reference_fasta,
        vcf_file = truth_vcf,
        out_fasta = mutant_fasta
    )
    Test.@test isfile(mutant_fasta)
    println("Mutated haplotype FASTA: ", mutant_fasta)
else
    println("Skipping consensus FASTA materialization.")
    println("Set MYCELIA_RUN_EXTERNAL=true to run bcftools/minimap2/GATK steps.")
    println("Add MYCELIA_RUN_VCFEVAL_TUTORIAL=true or MYCELIA_RUN_CLAIR3_TUTORIAL=true for optional sections.")
end

# ## Part 3: Simulate paired-end reads and align them with minimap2
#
# The same BAM can feed multiple callers. This is usually the first point where
# real projects branch into short-read versus long-read pipelines.

bam_file = ""

if RUN_EXTERNAL
    reads = Mycelia.simulate_illumina_reads(
        fasta = mutant_fasta,
        coverage = 40,
        read_length = 150,
        mflen = 350,
        sdev = 25,
        rndSeed = rand(rng, 0:typemax(Int)),
        quiet = true
    )

    map_result = Mycelia.minimap_map_paired_end(
        fasta = reference_fasta,
        forward = reads.forward_reads,
        reverse = reads.reverse_reads,
        threads = 2,
        outdir = joinpath(workdir, "alignment")
    )
    run(map_result.cmd)
    bam_file = map_result.outfile

    Test.@test isfile(bam_file)

    println("Forward reads: ", reads.forward_reads)
    println("Reverse reads: ", reads.reverse_reads)
    println("Aligned BAM: ", bam_file)
else
    println("Alignment section skipped because external execution is disabled.")
end

# ## Part 4: Call variants with GATK HaplotypeCaller
#
# Mycelia handles environment setup, BAM indexing, and GATK sequence dictionary
# generation. We add a FASTA index up front because downstream tools such as
# `bcftools norm` also rely on it.

gatk_vcf = joinpath(workdir, "gatk", "variants.gatk.vcf")
normalized_gatk_vcf = replace(gatk_vcf, ".vcf" => ".sorted.normalized.vcf.gz")
gatk_consensus_fasta = joinpath(workdir, "gatk", "consensus.fa")

if RUN_EXTERNAL
    Mycelia.run_gatk_haplotypecaller(
        bam_file,
        reference_fasta,
        gatk_vcf;
        ploidy = 1,
        memory_gb = 4,
        threads = 2
    )
    Test.@test isfile(gatk_vcf)

    normalized_gatk_vcf = Mycelia.normalize_vcf(
        reference_fasta = reference_fasta,
        vcf_file = gatk_vcf
    )
    Test.@test isfile(normalized_gatk_vcf)

    Mycelia.update_fasta_with_vcf(
        in_fasta = reference_fasta,
        vcf_file = gatk_vcf,
        out_fasta = gatk_consensus_fasta
    )
    Test.@test isfile(gatk_consensus_fasta)

    println("GATK VCF: ", gatk_vcf)
    println("Normalized GATK VCF: ", normalized_gatk_vcf)
    println("GATK consensus FASTA: ", gatk_consensus_fasta)
else
    println("GATK execution skipped.")
end

# ## Part 5: Benchmark against the truth set with RTG vcfeval
#
# `evaluate_variant_calling_accuracy` wraps the full vcfeval flow and returns
# parsed ROC tables plus summary metrics.

if RUN_VCFEVAL
    gatk_eval = Mycelia.evaluate_variant_calling_accuracy(
        truth_vcf,
        normalized_gatk_vcf,
        reference_fasta,
        joinpath(workdir, "gatk", "vcfeval");
        generate_plots = false,
        threads = 2,
        memory_gb = 4
    )

    println("GATK evaluation summary:")
    println(gatk_eval.summary_stats)
else
    println("vcfeval section skipped.")
end

# ## Part 6: Optional Clair3 workflow
#
# Clair3 is often the first long-read-oriented or hybrid benchmarking target
# users want beside GATK. The wrapper below uses the Illumina model because this
# tutorial simulates paired-end short reads.

clair3_output_dir = joinpath(workdir, "clair3")

if RUN_CLAIR3
    clair3_vcf = Mycelia.run_clair3(
        bam_file,
        reference_fasta,
        clair3_output_dir;
        platform = "ilmn",
        threads = 2,
        haploid_precise = true
    )
    Test.@test isfile(clair3_vcf)

    normalized_clair3_vcf = Mycelia.normalize_vcf(
        reference_fasta = reference_fasta,
        vcf_file = clair3_vcf
    )
    Test.@test isfile(normalized_clair3_vcf)

    println("Clair3 VCF: ", clair3_vcf)
    println("Normalized Clair3 VCF: ", normalized_clair3_vcf)

    if RUN_VCFEVAL
        clair3_eval = Mycelia.evaluate_variant_calling_accuracy(
            truth_vcf,
            normalized_clair3_vcf,
            reference_fasta,
            joinpath(clair3_output_dir, "vcfeval");
            generate_plots = false,
            threads = 2,
            memory_gb = 4
        )
        println("Clair3 evaluation summary:")
        println(clair3_eval.summary_stats)
    end
else
    println("Clair3 execution skipped.")
    println("Enable it with MYCELIA_RUN_CLAIR3_TUTORIAL=true after verifying the Clair3 model environment.")
end

# ## Part 7: Alternative callers
#
# The same alignment can also feed the other wrappers in `variant-analysis.jl`.
#
# ```julia
# freebayes_vcf = Mycelia.run_freebayes(
#     bam_file,
#     reference_fasta,
#     joinpath(workdir, "freebayes", "variants.freebayes.vcf");
#     ploidy = 1
# )
#
# bcftools_vcf = Mycelia.run_bcftools_call(
#     bam_file,
#     reference_fasta,
#     joinpath(workdir, "bcftools", "variants.bcftools.vcf");
#     threads = 2
# )
# ```
#
# For side-by-side benchmarking, use:
#
# ```julia
# comparison = Mycelia.run_variant_calling_comparison(
#     bam_file,
#     reference_fasta,
#     joinpath(workdir, "comparison");
#     callers = ["gatk", "clair3", "bcftools"],
#     platform = "ilmn",
#     baseline_vcf = truth_vcf
# )
# ```

# ## Citation Guidance
#
# If you publish results from this workflow, cite Mycelia and the wrapped tools
# you used, especially `minimap2`, `GATK`, `Clair3`, `BCFtools/HTSlib`, and
# `RTG vcfeval`. The project-wide citation checklist lives in
# `docs/src/references.md`.

println("Variant calling tutorial complete. Workspace: ", workdir)
