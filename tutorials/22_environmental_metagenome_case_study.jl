# # Tutorial 22: Environmental Metagenome Case Study
#
# This tutorial turns Mycelia's environmental metagenome capabilities into a
# concrete portfolio-ready story: a public freshwater metagenome cohort from
# Yukon River water at Pilot Station, Alaska, USA.
#
# ## Learning Objectives
# - Select a real public environmental metagenome cohort from bundled SRA metadata
# - Generate a cohort overview figure before downloading data
# - Run read-level QC summaries on per-sample FASTQ files
# - Compare samples with lightweight canonical k-mer Jaccard ordination
# - Save tables and figures for downstream presentation or portfolio use
#
# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/22_environmental_metagenome_case_study.jl", "tutorials/notebooks", execute=false)'
# ```
#
# ## Notes
# - Real SRA download is opt-in via `MYCELIA_DOWNLOAD_ENVIRONMENTAL_DATA=true`
# - When download is disabled, the tutorial generates small deterministic demo
#   FASTQ files so the QC and figure-generation path still runs locally

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import CSV
import DataFrames
import Mycelia
import Plots
import Random

Random.seed!(42)

const DOWNLOAD_ENVIRONMENTAL_DATA = lowercase(
    get(ENV, "MYCELIA_DOWNLOAD_ENVIRONMENTAL_DATA", "false")) == "true"
const ENVIRONMENTAL_RESULTS_DIR = mkpath(
    joinpath(@__DIR__, "..", "results", "environmental_metagenome_case_study"))

function synthesize_demo_fastqs(sample_names::AbstractVector{<:AbstractString}, outdir::AbstractString)
    mkpath(outdir)

    references = [
        repeat("ATGCGCGATCGGATCCGCGTATATCGCGAT", 20),
        repeat("ATATATGCGCGTATATATGCGCGGCGTATA", 20),
        repeat("GCGCGCGCATATATATATGCGCGCGATATA", 20)
    ]
    coverages = [80, 100, 120]
    error_rates = [0.01, 0.015, 0.02]

    fastqs = String[]
    for (i, sample_name) in enumerate(sample_names)
        reference = references[mod1(i, length(references))]
        coverage = coverages[mod1(i, length(coverages))]
        error_rate = error_rates[mod1(i, length(error_rates))]
        reads = Mycelia.create_test_reads(reference, coverage, error_rate)
        fastq_path = joinpath(outdir, "$(sample_name)_demo.fastq")
        Mycelia.write_fastq(records = reads, filename = fastq_path)
        push!(fastqs, fastq_path)
    end

    return fastqs
end

case_study = Mycelia.environmental_metagenome_case_study(max_runs = 3)
sample_table = case_study.samples

println("=== Environmental Metagenome Case Study ===")
println(case_study.title)
println("Habitat: $(case_study.habitat)")
println("Location: $(case_study.location)")
println("Selected runs:")
display(sample_table)

CSV.write(joinpath(ENVIRONMENTAL_RESULTS_DIR, "selected_runs.csv"), sample_table)

# ## Figure 1: Cohort Overview
#
# This figure is available immediately from the bundled SRA metadata.

metadata_plot = Plots.bar(
    sample_table.run,
    sample_table.bytes_gb;
    title = "Yukon River Environmental Metagenome Cohort",
    xlabel = "Run",
    ylabel = "Compressed FASTQ Size (GB)",
    legend = false,
    color = :steelblue,
    size = (900, 500),
    xrotation = 20
)
Plots.savefig(metadata_plot, joinpath(ENVIRONMENTAL_RESULTS_DIR, "cohort_overview.png"))

# ## Acquire Data
#
# When `MYCELIA_DOWNLOAD_ENVIRONMENTAL_DATA=true`, Mycelia downloads the real
# public FASTQs. Otherwise, we generate small deterministic demo FASTQs so the
# downstream QC and ordination steps remain runnable offline.

data_dir = mkpath(joinpath(ENVIRONMENTAL_RESULTS_DIR, "data"))
sample_names = collect(sample_table.run)
fastq_files = String[]

if DOWNLOAD_ENVIRONMENTAL_DATA
    println("\nDownloading public Yukon River SRA runs...")
    for sample_name in sample_names
        download_result = Mycelia.download_sra_data(sample_name; outdir = data_dir)
        push!(fastq_files, download_result.files[1])
    end
else
    println("\nMYCELIA_DOWNLOAD_ENVIRONMENTAL_DATA=false, generating demo FASTQs for local execution.")
    fastq_files = synthesize_demo_fastqs(sample_names, data_dir)
end

println("FASTQ inputs:")
for fastq_file in fastq_files
    println("  $(fastq_file)")
end

# ## Read-Level QC

quality_summary = Mycelia.environmental_fastq_quality_summary(
    fastq_files;
    sample_names = sample_names
)
display(quality_summary)
CSV.write(joinpath(ENVIRONMENTAL_RESULTS_DIR, "quality_summary.tsv"), quality_summary; delim = '\t')

quality_plot = Plots.scatter(
    quality_summary.gc_content,
    quality_summary.mean_quality;
    title = "Environmental Metagenome QC Summary",
    xlabel = "GC Content (%)",
    ylabel = "Mean Phred Quality",
    legend = false,
    markersize = 7,
    color = :darkgreen,
    size = (800, 500)
)

for row in DataFrames.eachrow(quality_summary)
    Plots.annotate!(quality_plot, row.gc_content, row.mean_quality, row.sample)
end

Plots.savefig(quality_plot, joinpath(ENVIRONMENTAL_RESULTS_DIR, "qc_gc_vs_quality.png"))

# ## K-mer Ordination
#
# We intentionally cap the number of reads used per sample so the comparison can
# run on a laptop while still capturing strong between-sample signal.

kmer_profiles = Mycelia.analyze_environmental_kmer_profiles(
    fastq_files;
    sample_names = sample_names,
    k = 21,
    max_reads = 2_000
)
display(kmer_profiles.spectrum_summary)
CSV.write(
    joinpath(ENVIRONMENTAL_RESULTS_DIR, "kmer_spectrum_summary.tsv"),
    kmer_profiles.spectrum_summary;
    delim = '\t'
)

pcoa_df = kmer_profiles.pcoa_df
if !DataFrames.hasproperty(pcoa_df, :PC2)
    pcoa_df[!, :PC2] = zeros(DataFrames.nrow(pcoa_df))
end

pcoa_plot = Plots.scatter(
    pcoa_df.PC1,
    pcoa_df.PC2;
    title = "Canonical 21-mer Jaccard PCoA",
    xlabel = "PC1",
    ylabel = "PC2",
    legend = false,
    markersize = 7,
    color = :darkorange,
    size = (800, 500)
)

for row in DataFrames.eachrow(pcoa_df)
    Plots.annotate!(pcoa_plot, row.PC1, row.PC2, row.sample)
end

Plots.savefig(pcoa_plot, joinpath(ENVIRONMENTAL_RESULTS_DIR, "kmer_jaccard_pcoa.png"))

spectrum_plot = Plots.bar(
    kmer_profiles.spectrum_summary.sample,
    kmer_profiles.spectrum_summary.unique_kmers;
    title = "Unique Canonical 21-mers per Sample",
    xlabel = "Run",
    ylabel = "Unique 21-mers",
    legend = false,
    color = :mediumpurple,
    size = (900, 500),
    xrotation = 20
)
Plots.savefig(spectrum_plot, joinpath(ENVIRONMENTAL_RESULTS_DIR, "unique_kmers_per_sample.png"))

# ## Summary

println("\nSaved outputs to $(ENVIRONMENTAL_RESULTS_DIR)")
println("Generated figures:")
println("  - cohort_overview.png")
println("  - qc_gc_vs_quality.png")
println("  - kmer_jaccard_pcoa.png")
println("  - unique_kmers_per_sample.png")

nothing
