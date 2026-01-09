# # Tutorial 2: Quality Control and Preprocessing
#
# This tutorial demonstrates how to assess and improve the quality of genomic data
# before analysis. Quality control is essential for ensuring reliable downstream results.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - How to assess sequencing data quality using multiple read-level metrics
# - Common quality issues and their biological implications
# - Preprocessing techniques for improving read quality
# - How to use Mycelia's read-level QC tooling
# - Best practices for quality control in different sequencing data types

# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/02_quality_control.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import Statistics
import Random
import Plots

Random.seed!(42)

# ## Part 1: Understanding Quality Metrics
#
# Quality assessment involves multiple metrics that capture different aspects
# of read quality. Understanding these metrics helps identify problems and
# guide preprocessing decisions.

println("=== Quality Control Tutorial ===")

# ### Phred Quality Scores
#
# Phred scores represent the probability of base-calling errors:
# - Q10 = 10% error rate (1 in 10 bases wrong)
# - Q20 = 1% error rate (1 in 100 bases wrong)
# - Q30 = 0.1% error rate (1 in 1000 bases wrong)
# - Q40 = 0.01% error rate (1 in 10,000 bases wrong)

println("Quality Score Interpretation:")
println("Q10: 10% error rate (poor quality)")
println("Q20: 1% error rate (acceptable)")
println("Q30: 0.1% error rate (good quality)")
println("Q40: 0.01% error rate (excellent quality)")

# ## Part 2: Example Read Data
#
# Generate small, self-contained FASTQ datasets using Mycelia utilities.

println("\n=== Example Read Data ===")

output_dir = mkpath(joinpath(@__DIR__, "..", "results", "qc_tutorial"))

short_reference = Mycelia.random_fasta_record(moltype=:DNA, seed=1, L=150)
short_reads = Mycelia.create_test_reads(FASTX.sequence(short_reference), 200, 0.01)
short_fastq = joinpath(output_dir, "short_reads.fastq")
Mycelia.write_fastq(records=short_reads, filename=short_fastq)

forward_reads = Mycelia.create_test_reads(FASTX.sequence(short_reference), 200, 0.01)
reverse_reads = Mycelia.create_test_reads(FASTX.sequence(short_reference), 200, 0.01)
forward_fastq = joinpath(output_dir, "short_reads_R1.fastq")
reverse_fastq = joinpath(output_dir, "short_reads_R2.fastq")
Mycelia.write_fastq(records=forward_reads, filename=forward_fastq)
Mycelia.write_fastq(records=reverse_reads, filename=reverse_fastq)

long_reference = Mycelia.random_fasta_record(moltype=:DNA, seed=2, L=4000)
long_reads = Mycelia.create_test_reads(FASTX.sequence(long_reference), 30, 0.08)
long_fastq = joinpath(output_dir, "long_reads.fastq")
Mycelia.write_fastq(records=long_reads, filename=long_fastq)

println("Wrote example reads to $(output_dir)")
println("  Short reads: $(short_fastq)")
println("  Paired-end reads: $(forward_fastq), $(reverse_fastq)")
println("  Long reads: $(long_fastq)")

# ## Part 3: Read-Level Quality Assessment
#
# Use Mycelia's QC utilities to summarize read quality.

println("\n=== Read-Level Quality Assessment ===")

short_stats = Mycelia.analyze_fastq_quality(short_fastq)
println("Short-read summary:")
println("  Reads: $(short_stats.n_reads)")
println("  Mean quality: $(round(short_stats.mean_quality, digits=2))")
println("  Mean length: $(round(short_stats.mean_length, digits=1))")
println("  GC content: $(round(short_stats.gc_content, digits=1))%")
println("  Q30+ reads: $(round(short_stats.quality_distribution.q30_percent, digits=1))%")

per_read_quals = Mycelia.per_read_quality_scores(short_reads)
gc_contents = Mycelia.gc_content_per_read(short_reads)
read_lengths = Mycelia.read_length_distribution(short_reads)
dup_stats = Mycelia.duplication_stats(short_reads; min_fraction=0.05)

println("Read-level distribution metrics:")
println("  Mean per-read quality: $(round(Statistics.mean(per_read_quals), digits=2))")
println("  Mean GC: $(round(Statistics.mean(gc_contents) * 100, digits=1))%")
println("  Mean read length: $(round(Statistics.mean(read_lengths), digits=1)) bp")
println("  Duplicate fraction: $(round(dup_stats.duplicate_fraction * 100, digits=1))%")
println("  Overrepresented sequences: $(length(dup_stats.overrepresented))")

# ## Part 4: Automated Read-Level QC and Filtering
#
# Mycelia provides wrappers around standard read QC tools.

println("\n=== Automated Read-Level QC ===")

println("Mycelia will auto-install missing tools via conda when you run these wrappers.")

fastqc_dir = Mycelia.run_fastqc(fastq=short_fastq)
println("FastQC output: $(fastqc_dir)")

fastp_outputs = Mycelia.qc_filter_short_reads_fastp(
    forward_reads=forward_fastq,
    reverse_reads=reverse_fastq
)
println("fastp outputs:")
println("  Filtered R1: $(fastp_outputs.out_forward)")
println("  Filtered R2: $(fastp_outputs.out_reverse)")
println("  JSON report: $(fastp_outputs.json)")
println("  HTML report: $(fastp_outputs.html)")

trim_outputs = Mycelia.trim_galore_paired(
    forward_reads=forward_fastq,
    reverse_reads=reverse_fastq,
    outdir=output_dir
)
println("trim_galore outputs:")
println("  Trimmed R1: $(trim_outputs.trimmed_forward)")
println("  Trimmed R2: $(trim_outputs.trimmed_reverse)")

fastplong_outputs = Mycelia.qc_filter_long_reads_fastplong(
    in_fastq=long_fastq,
    min_length=1000
)
println("fastplong outputs:")
println("  Filtered reads: $(fastplong_outputs.out_fastq)")
println("  JSON report: $(fastplong_outputs.json_report)")
println("  HTML report: $(fastplong_outputs.html_report)")

# ## Part 5: Quality Visualization
#
# Create plots to visualize read-level metrics.

println("\n=== Quality Visualization ===")

quality_dist_plot = Plots.histogram(
    per_read_quals;
    bins=0:2:50,
    title="Per-read Quality Distribution",
    xlabel="Mean Phred Score",
    ylabel="Count",
    legend=false
)

gc_plot = Plots.histogram(
    gc_contents .* 100;
    bins=0:5:100,
    title="GC Content Distribution",
    xlabel="GC %",
    ylabel="Count",
    legend=false
)

length_plot = Plots.histogram(
    read_lengths;
    bins=10,
    title="Read Length Distribution",
    xlabel="Length (bp)",
    ylabel="Count",
    legend=false
)

per_base_plot = Mycelia.plot_per_base_quality(short_fastq)

if isinteractive()
    Plots.display(quality_dist_plot)
    Plots.display(gc_plot)
    Plots.display(length_plot)
    Plots.display(per_base_plot)
else
    Plots.savefig(quality_dist_plot, joinpath(output_dir, "quality_distribution.png"))
    Plots.savefig(gc_plot, joinpath(output_dir, "gc_distribution.png"))
    Plots.savefig(length_plot, joinpath(output_dir, "length_distribution.png"))
    Plots.savefig(per_base_plot, joinpath(output_dir, "per_base_quality.png"))
    println("Saved QC plots to $(output_dir)")
end

# ## Summary
println("\n=== Quality Control Summary ===")
println("✓ Interpreted Phred quality scores")
println("✓ Assessed read-level quality metrics")
println("✓ Applied automated read QC with fastp/fastplong/trim_galore")
println("✓ Visualized quality distributions")

println("\nNext: Tutorial 3 - K-mer Analysis and Feature Extraction")

nothing
