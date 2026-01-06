# # Tutorial 5: Assembly Validation
#
# This tutorial covers comprehensive validation techniques for genome assemblies,
# including reference-based and reference-free approaches.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - Reference-based validation using alignment and comparison tools
# - Reference-free validation using k-mer and read-based approaches
# - BUSCO analysis for gene completeness assessment
# - Merqury for k-mer based quality assessment
# - Statistical validation and confidence intervals
# - Comparative validation across multiple assemblies

# ## Setup
import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import CSV
import DataFrames
import Plots
import Random
import Statistics

Random.seed!(42)

output_dir = mkpath(joinpath(@__DIR__, "..", "results", "assembly_validation"))

# ## Part 1: Validation Strategy Overview
#
# Assembly validation requires multiple complementary approaches to ensure
# comprehensive quality assessment.

println("=== Assembly Validation Tutorial ===")

println("Validation Approaches:")
println("1. Reference-based validation")
println("   - Alignment to known reference")
println("   - Synteny analysis")
println("   - Structural variant detection")
println()
println("2. Reference-free validation")
println("   - K-mer based approaches (Merqury)")
println("   - Read mapping validation")
println("   - Internal consistency checks")
println()
println("3. Functional validation")
println("   - Gene completeness (BUSCO)")
println("   - Annotation quality")
println("   - Comparative genomics")

# ## Part 2: Reference-Based Validation
#
# When reference genomes are available, direct comparison provides
# the most straightforward validation approach.

println("\n=== Reference-Based Validation ===")

# ### Preparing Test Data
#
# Create reference and assembly for validation demonstration

println("--- Preparing Test Data ---")

# Reference genome
reference_size = 100000  ## 100 kb
reference_genome = Mycelia.random_fasta_record(moltype=:DNA, seed=1, L=reference_size)
reference_file = joinpath(output_dir, "reference.fasta")
Mycelia.write_fasta(outfile=reference_file, records=[reference_genome])

# Simulated assembly with some differences
reference_seq = FASTX.FASTA.sequence(reference_genome)
contig_1 = FASTX.FASTA.Record("contig_1", reference_seq[1:40000])
contig_2_seq, _ = Mycelia.observe(reference_seq[40001:80000]; error_rate=0.005)
contig_2 = FASTX.FASTA.Record("contig_2", contig_2_seq)
contig_3 = FASTX.FASTA.Record("contig_3", reference_seq[80001:100000])
assembly_contigs = [contig_1, contig_2, contig_3]

assembly_file = joinpath(output_dir, "assembly.fasta")
Mycelia.write_fasta(outfile=assembly_file, records=assembly_contigs)

reads = Mycelia.create_test_reads(reference_seq, 200, 0.01)
reads_file = joinpath(output_dir, "reads.fastq")
Mycelia.write_fastq(records=reads, filename=reads_file)

println("Test data prepared:")
println("  Reference: $(reference_size) bp")
println("  Assembly: $(length(assembly_contigs)) contigs")
println("  Reads: $(reads_file)")

# ### Alignment-Based Validation
#
# Use alignment tools to compare assembly to reference

println("--- Alignment-Based Validation ---")

println("Mycelia will auto-install missing tools via conda when you run these wrappers.")

mummer_dir = Mycelia.run_mummer(
    reference_file,
    assembly_file;
    outdir=joinpath(output_dir, "mummer"),
    prefix="validation"
)
coords_file = joinpath(mummer_dir, "validation.coords")
coords_lines = filter(line -> !isempty(strip(line)) && occursin(r"^\d", strip(line)), readlines(coords_file))
coords_rows = [split(strip(line)) for line in coords_lines if length(split(strip(line))) >= 7]
aligned_lengths = [parse(Int, row[5]) for row in coords_rows]
percent_identities = [parse(Float64, row[7]) for row in coords_rows]
total_aligned = isempty(aligned_lengths) ? 0 : sum(aligned_lengths)
mean_identity = isempty(percent_identities) ? missing : Statistics.mean(percent_identities)
coverage_estimate = total_aligned == 0 ? 0.0 : total_aligned / reference_size * 100

quast_dir = Mycelia.run_quast(
    assembly_file;
    reference=reference_file,
    outdir=joinpath(output_dir, "quast")
)
quast_report = CSV.read(joinpath(quast_dir, "report.tsv"), DataFrames.DataFrame; delim='\t')
metric_col = names(quast_report)[1]
value_col = names(quast_report)[2]
quast_metrics = Dict(row[metric_col] => row[value_col] for row in DataFrames.eachrow(quast_report))
quast_misassemblies = get(quast_metrics, "Misassemblies", get(quast_metrics, "# misassemblies", missing))

alignment_stats = Dict(
    "total_aligned_bp" => total_aligned,
    "mean_identity_percent" => mean_identity,
    "estimated_coverage_percent" => coverage_estimate,
    "quast_genome_fraction_percent" => get(quast_metrics, "Genome fraction (%)", missing),
    "quast_avg_identity_percent" => get(quast_metrics, "Avg. % identity", missing),
    "quast_misassemblies" => quast_misassemblies
)

println("Alignment Statistics:")
for (metric, value) in alignment_stats
    println("  $metric: $value")
end

# ### Synteny Analysis
#
# Analyze conserved gene order and chromosomal structure

println("--- Synteny Analysis ---")

dotplot_file = Mycelia.run_mummer_plot(
    joinpath(mummer_dir, "validation.delta");
    outdir=mummer_dir,
    prefix="validation_dotplot",
    plot_type="png"
)
println("Synteny dot plot: $(dotplot_file)")

# ### Structural Variant Detection
#
# Identify large-scale differences between assembly and reference

println("--- Structural Variant Detection ---")

structural_variants = Dict(
    "misassemblies" => quast_misassemblies,
    "relocations" => get(quast_metrics, "Relocations", missing),
    "translocations" => get(quast_metrics, "Translocations", missing),
    "inversions" => get(quast_metrics, "Inversions", missing)
)

println("Structural variant summary (QUAST):")
for (metric, value) in structural_variants
    println("  $metric: $value")
end

# ## Part 3: Reference-Free Validation
#
# When no reference is available, use intrinsic data properties
# for validation.

println("\n=== Reference-Free Validation ===")

# ### K-mer Based Validation (Merqury)
#
# Use k-mer analysis to assess assembly quality without reference

println("--- K-mer Based Validation ---")

qv_results = Mycelia.assess_assembly_quality(
    assembly=assembly_file,
    observations=[reads_file],
    ks=[17, 21]
)
merqury_results = Dict(
    "mean_qv" => Statistics.mean(qv_results.qv),
    "best_qv" => maximum(qv_results.qv),
    "mean_js_divergence" => Statistics.mean(qv_results.js_divergence)
)

println("Merqury-style QV results:")
println(qv_results)
for (metric, value) in merqury_results
    println("  $metric: $(round(value, digits=3))")
end

# ### Read Mapping Validation
#
# Map original reads back to assembly for validation

println("--- Read Mapping Validation ---")

map_result = Mycelia.minimap_map(
    fasta=assembly_file,
    fastq=reads_file,
    mapping_type="map-hifi",
    threads=4,
    output_format="bam",
    sorted=true,
    quiet=false
)
run(map_result.cmd)
bam_file = map_result.outfile

flagstat_file = Mycelia.run_samtools_flagstat(bam_file)
flagstat_lines = readlines(flagstat_file)
total_reads = isempty(flagstat_lines) ? 0 : parse(Int, split(flagstat_lines[1])[1])
mapped_line_idx = findfirst(line -> occursin(" mapped", line), flagstat_lines)
mapped_reads = mapped_line_idx === nothing ? 0 : parse(Int, split(flagstat_lines[mapped_line_idx])[1])
mapped_percent = total_reads == 0 ? 0.0 : mapped_reads / total_reads * 100

mosdepth_outputs = Mycelia.run_mosdepth(bam_file; thresholds="1,10,30", no_per_base=true)
dist_df = Mycelia.parse_mosdepth_distribution(mosdepth_outputs.global_dist)
summary_df = Mycelia.parse_mosdepth_summary(mosdepth_outputs.summary)
coverage_qc = Mycelia.summarize_mosdepth_qc(dist_df; thresholds=[1, 10, 30])
coverage_total = DataFrames.subset(coverage_qc, :chromosome => x -> x .== "total")
fraction_1x = DataFrames.nrow(coverage_total) == 0 ? missing : coverage_total[1, :coverage_1X]
fraction_10x = DataFrames.nrow(coverage_total) == 0 ? missing : coverage_total[1, :coverage_10X]
zero_cov_fraction = ismissing(fraction_1x) ? missing : 1 - fraction_1x

chrom_col = names(summary_df)[1]
mean_col = names(summary_df)[4]
total_rows = summary_df[summary_df[!, chrom_col] .== "total", :]
mean_coverage = DataFrames.nrow(total_rows) == 0 ? missing : total_rows[1, mean_col]

qualimap_outputs = Mycelia.run_qualimap_bamqc(
    bam=bam_file,
    outdir=joinpath(output_dir, "qualimap")
)
contig_cov = Mycelia.parse_qualimap_contig_coverage(qualimap_outputs.report_txt)
coverage_ratios = contig_cov[!, "Standard Deviation"] ./ contig_cov[!, "Mean coverage"]
coverage_ratios = filter(isfinite, coverage_ratios)
coverage_cv = isempty(coverage_ratios) ? missing : Statistics.mean(coverage_ratios)

mapping_stats = Dict(
    "mapped_reads_percent" => mapped_percent,
    "mean_coverage" => mean_coverage,
    "coverage_uniformity_cv" => coverage_cv,
    "zero_coverage_fraction" => zero_cov_fraction
)

println("Read Mapping Statistics:")
for (metric, value) in mapping_stats
    println("  $metric: $value")
end

# ### Internal Consistency Validation
#
# Check assembly internal consistency

println("--- Internal Consistency ---")

contig_strings = [String(FASTX.FASTA.sequence(contig)) for contig in assembly_contigs]
contig_names = [String(FASTX.FASTA.identifier(contig)) for contig in assembly_contigs]
rhizomorph_result = Mycelia.Rhizomorph.AssemblyResult(contig_strings, contig_names)
structure_report = Mycelia.Rhizomorph.validate_assembly_structure(rhizomorph_result)
basic_report = Mycelia.Rhizomorph.validate_assembly(rhizomorph_result)

println("Rhizomorph structure report:")
println("  valid: $(structure_report["valid"])")
println("  issues: $(length(structure_report["issues"]))")
println("  warnings: $(length(structure_report["warnings"]))")
println("Rhizomorph summary metrics:")
for (metric, value) in basic_report
    println("  $metric: $value")
end

# ## Part 4: Functional Validation
#
# Validate assembly quality through functional analysis

println("\n=== Functional Validation ===")

# ### BUSCO Analysis
#
# Assess gene completeness using conserved orthologs

println("--- BUSCO Analysis ---")

busco_dir = Mycelia.run_busco(
    assembly_file;
    outdir=joinpath(output_dir, "busco"),
    auto_lineage=true
)
busco_summary_files = filter(file -> occursin("short_summary", file), readdir(busco_dir; join=true))
busco_summary_file = isempty(busco_summary_files) ? "" : first(busco_summary_files)
busco_text = isempty(busco_summary_file) ? "" : read(busco_summary_file, String)
busco_match = match(r"C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)", busco_text)
busco_results = busco_match === nothing ? Dict("summary_file" => busco_summary_file) : Dict(
    "complete" => parse(Float64, busco_match.captures[1]),
    "complete_single_copy" => parse(Float64, busco_match.captures[2]),
    "complete_duplicated" => parse(Float64, busco_match.captures[3]),
    "fragmented" => parse(Float64, busco_match.captures[4]),
    "missing" => parse(Float64, busco_match.captures[5]),
    "total_buscos" => parse(Int, busco_match.captures[6])
)

println("BUSCO Results:")
for (metric, value) in busco_results
    println("  $metric: $value")
end

# ### Gene Annotation Quality
#
# Validate through gene prediction and annotation

println("--- Gene Annotation Quality ---")

pyrodigal_assembly = Mycelia.run_pyrodigal(
    fasta_file=assembly_file,
    out_dir=joinpath(output_dir, "pyrodigal_assembly")
)
pyrodigal_reference = Mycelia.run_pyrodigal(
    fasta_file=reference_file,
    out_dir=joinpath(output_dir, "pyrodigal_reference")
)

assembly_gff = CSV.read(pyrodigal_assembly.gff, DataFrames.DataFrame; delim='\t', comment="#", header=false)
DataFrames.rename!(assembly_gff, [:seqid, :source, :feature, :start, :stop, :score, :strand, :phase, :attributes])
reference_gff = CSV.read(pyrodigal_reference.gff, DataFrames.DataFrame; delim='\t', comment="#", header=false)
DataFrames.rename!(reference_gff, [:seqid, :source, :feature, :start, :stop, :score, :strand, :phase, :attributes])

assembly_gene_lengths = abs.(assembly_gff.start .- assembly_gff.stop) .+ 1
reference_gene_lengths = abs.(reference_gff.start .- reference_gff.stop) .+ 1
frame_consistency = isempty(assembly_gene_lengths) ? missing : Statistics.mean(assembly_gene_lengths .% 3 .== 0)

annotation_stats = Dict(
    "assembly_gene_count" => DataFrames.nrow(assembly_gff),
    "reference_gene_count" => DataFrames.nrow(reference_gff),
    "assembly_mean_gene_length" => isempty(assembly_gene_lengths) ? missing : Statistics.mean(assembly_gene_lengths),
    "reference_mean_gene_length" => isempty(reference_gene_lengths) ? missing : Statistics.mean(reference_gene_lengths),
    "frame_consistency_fraction" => frame_consistency
)

println("Annotation quality summary:")
for (metric, value) in annotation_stats
    printable = value isa Number ? round(value, digits=3) : value
    println("  $metric: $(printable)")
end

# ## Part 5: Comparative Validation
#
# Compare multiple assemblies to identify best approach

println("\n=== Comparative Validation ===")

# ### Multi-Assembly Comparison
#
# Compare assemblies from different tools or parameters

println("--- Multi-Assembly Comparison ---")

alt_contig_1 = FASTX.FASTA.Record("alt1_contig_1", reference_seq[1:30000])
alt_contig_2 = FASTX.FASTA.Record("alt1_contig_2", reference_seq[30001:60000])
alt_contig_3 = FASTX.FASTA.Record("alt1_contig_3", reference_seq[60001:90000])
alt_contig_4 = FASTX.FASTA.Record("alt1_contig_4", reference_seq[90001:100000])
assembly_alt1 = joinpath(output_dir, "assembly_alt1.fasta")
Mycelia.write_fasta(outfile=assembly_alt1, records=[alt_contig_1, alt_contig_2, alt_contig_3, alt_contig_4])

alt2_seq, _ = Mycelia.observe(reference_seq; error_rate=0.01)
assembly_alt2 = joinpath(output_dir, "assembly_alt2.fasta")
Mycelia.write_fasta(outfile=assembly_alt2, records=[FASTX.FASTA.Record("alt2_contig_1", alt2_seq)])

assemblies = [assembly_file, assembly_alt1, assembly_alt2]
quast_multi_dir = Mycelia.run_quast(
    assemblies;
    reference=reference_file,
    outdir=joinpath(output_dir, "quast_multi")
)
quast_multi_report = CSV.read(joinpath(quast_multi_dir, "report.tsv"), DataFrames.DataFrame; delim='\t')
metric_col_multi = names(quast_multi_report)[1]
assembly_cols = names(quast_multi_report)[2:end]

comparison_table = DataFrames.DataFrame(
    assembly=String[],
    n50=Union{Float64, Missing}[],
    genome_fraction=Union{Float64, Missing}[],
    misassemblies=Union{Float64, Missing}[]
)
for col in assembly_cols
    n50_idx = findfirst(==("N50"), quast_multi_report[!, metric_col_multi])
    genome_idx = findfirst(==("Genome fraction (%)"), quast_multi_report[!, metric_col_multi])
    mis_idx = findfirst(==("Misassemblies"), quast_multi_report[!, metric_col_multi])
    n50_val = n50_idx === nothing ? missing : quast_multi_report[n50_idx, col]
    genome_val = genome_idx === nothing ? missing : quast_multi_report[genome_idx, col]
    mis_val = mis_idx === nothing ? missing : quast_multi_report[mis_idx, col]
    push!(comparison_table, (assembly=String(col), n50=n50_val, genome_fraction=genome_val, misassemblies=mis_val))
end

println("Assembly Comparison (QUAST):")
println(comparison_table)

# ### Statistical Validation
#
# Apply statistical tests to validation results

println("--- Statistical Validation ---")

contig_lengths = [length(FASTX.FASTA.sequence(contig)) for contig in assembly_contigs]
bootstrap_n50 = Float64[]
for _ in 1:200
    resampled = rand(contig_lengths, length(contig_lengths))
    sorted_lengths = sort(resampled, rev=true)
    cumulative = cumsum(sorted_lengths)
    target = sum(sorted_lengths) / 2
    idx = findfirst(x -> x >= target, cumulative)
    push!(bootstrap_n50, sorted_lengths[idx])
end
ci_low, ci_high = Statistics.quantile(bootstrap_n50, [0.025, 0.975])

qv_alt1 = Mycelia.assess_assembly_quality(assembly=assembly_alt1, observations=[reads_file], ks=[17, 21])
observed_diff = Statistics.mean(qv_results.qv) - Statistics.mean(qv_alt1.qv)
combined_qv = vcat(qv_results.qv, qv_alt1.qv)
n1 = length(qv_results.qv)
perm_diffs = Float64[]
for _ in 1:1000
    permuted = Random.shuffle(combined_qv)
    diff = Statistics.mean(permuted[1:n1]) - Statistics.mean(permuted[n1+1:end])
    push!(perm_diffs, diff)
end
p_value = Statistics.mean(abs.(perm_diffs) .>= abs(observed_diff))
bonferroni_p = min(p_value * 2, 1.0)
effect_size = observed_diff / Statistics.std(combined_qv)

println("Bootstrap N50 CI: ($(round(ci_low, digits=2)), $(round(ci_high, digits=2)))")
println("QV mean difference (assembly vs alt1): $(round(observed_diff, digits=3))")
println("Permutation p-value: $(round(p_value, digits=4))")
println("Bonferroni-adjusted p-value: $(round(bonferroni_p, digits=4))")
println("Effect size (Cohen's d): $(round(effect_size, digits=3))")

# ## Part 6: Validation Metrics Integration
#
# Combine multiple validation approaches for comprehensive assessment

println("\n=== Integrated Validation ===")

# ### Composite Quality Scores
#
# Combine multiple metrics into overall quality assessment

println("--- Composite Quality Scores ---")

n_contigs, total_length, n50, l50 = Mycelia.assess_assembly_quality(assembly_file)
busco_complete = get(busco_results, "complete", missing)
busco_norm = busco_complete isa Number ? busco_complete / 100 : 0.0
qv_norm = merqury_results["best_qv"] / 60
contiguity_norm = n50 / reference_size

weights = Dict("contiguity" => 0.4, "accuracy" => 0.35, "completeness" => 0.25)
overall_quality = 10 * (
    weights["contiguity"] * contiguity_norm +
    weights["accuracy"] * qv_norm +
    weights["completeness"] * busco_norm
)
composite_samples = [
    10 * (
        weights["contiguity"] * (n50_sample / reference_size) +
        weights["accuracy"] * qv_norm +
        weights["completeness"] * busco_norm
    )
    for n50_sample in bootstrap_n50
]
quality_ci = Statistics.quantile(composite_samples, [0.025, 0.975])

composite_score = Dict(
    "overall_quality" => overall_quality,
    "confidence_interval" => (quality_ci[1], quality_ci[2]),
    "primary_strengths" => ["Contiguity", "K-mer accuracy"],
    "primary_weaknesses" => ["BUSCO completeness sensitivity"]
)

println("Composite Quality Assessment:")
for (metric, value) in composite_score
    println("  $metric: $value")
end

# ### Validation Report Generation
#
# Create comprehensive validation reports

println("--- Validation Report ---")

report_rows = DataFrames.DataFrame(metric=String[], value=String[])
for (metric, value) in alignment_stats
    push!(report_rows, (metric="alignment_$(metric)", value=string(value)))
end
for (metric, value) in merqury_results
    push!(report_rows, (metric="merqury_$(metric)", value=string(value)))
end
for (metric, value) in mapping_stats
    push!(report_rows, (metric="mapping_$(metric)", value=string(value)))
end
for (metric, value) in busco_results
    push!(report_rows, (metric="busco_$(metric)", value=string(value)))
end
for (metric, value) in annotation_stats
    push!(report_rows, (metric="annotation_$(metric)", value=string(value)))
end
for (metric, value) in composite_score
    push!(report_rows, (metric="composite_$(metric)", value=string(value)))
end

report_tsv = joinpath(output_dir, "assembly_validation_report.tsv")
CSV.write(report_tsv, report_rows; delim='\t')
println("Saved validation report: $(report_tsv)")

# ## Part 7: Validation Visualization
#
# Create visualizations for validation results

println("\n=== Validation Visualization ===")

# ### Quality Metric Plots
#
# Visualize validation metrics

println("--- Quality Metric Plots ---")

contig_lengths_plot = Plots.histogram(
    contig_lengths;
    bins=10,
    xlabel="Contig length (bp)",
    ylabel="Count",
    title="Contig Length Distribution"
)
contig_lengths_path = joinpath(output_dir, "contig_lengths.png")
Plots.savefig(contig_lengths_plot, contig_lengths_path)

qv_alt2 = Mycelia.assess_assembly_quality(assembly=assembly_alt2, observations=[reads_file], ks=[17, 21])
qv_results[!, :assembler] .= "assembly"
qv_alt1[!, :assembler] .= "alt1"
qv_alt2[!, :assembler] .= "alt2"
qv_all = vcat(qv_results, qv_alt1, qv_alt2)
qv_heatmap = Mycelia.generate_qv_heatmap(qv_all; assembler_column="assembler")
qv_heatmap_path = joinpath(output_dir, "qv_heatmap.png")
Plots.savefig(qv_heatmap, qv_heatmap_path)

_, _, n50_alt1, _ = Mycelia.assess_assembly_quality(assembly_alt1)
_, _, n50_alt2, _ = Mycelia.assess_assembly_quality(assembly_alt2)
mean_qv_table = DataFrames.DataFrame(
    assembly=["assembly", "alt1", "alt2"],
    n50=[n50, n50_alt1, n50_alt2],
    mean_qv=[Statistics.mean(qv_results.qv), Statistics.mean(qv_alt1.qv), Statistics.mean(qv_alt2.qv)]
)
metric_scatter = Plots.scatter(
    mean_qv_table.n50,
    mean_qv_table.mean_qv;
    xlabel="N50 (bp)",
    ylabel="Mean QV",
    title="Contiguity vs Accuracy",
    series_annotations=mean_qv_table.assembly
)
metric_scatter_path = joinpath(output_dir, "n50_vs_qv.png")
Plots.savefig(metric_scatter, metric_scatter_path)

println("Saved plots:")
println("  Contig lengths: $(contig_lengths_path)")
println("  QV heatmap: $(qv_heatmap_path)")
println("  N50 vs QV: $(metric_scatter_path)")

# ### Genome Browser Integration
#
# Visualize validation results in genome browser context

println("--- Genome Browser Integration ---")

browser_tracks = Mycelia.run_mosdepth(
    bam_file;
    prefix=joinpath(output_dir, "assembly_coverage"),
    thresholds="1,10,30"
)
println("Genome browser tracks:")
println("  Per-base coverage: $(browser_tracks.per_base)")
println("  Thresholds BED: $(browser_tracks.thresholds_file)")
println("  Load these files into IGV or JBrowse for interactive inspection.")

# ## Part 8: Validation Best Practices
#
# Guidelines for effective assembly validation

println("\n=== Validation Best Practices ===")

println("Validation Strategy:")
println("- Use multiple complementary approaches")
println("- Always validate with original data")
println("- Compare with related genomes when available")
println("- Focus on metrics relevant to your research goals")
println()
println("Quality Thresholds:")
println("- BUSCO completeness: >90% for eukaryotes, >95% for prokaryotes")
println("- Assembly QV: >30 for high-quality assemblies")
println("- Read mapping: >95% of reads should map")
println("- Contig N50: Should be substantial fraction of chromosome size")
println()
println("Common Pitfalls:")
println("- Relying on single validation metric")
println("- Ignoring biological context")
println("- Not validating with original data")
println("- Accepting assemblies without proper validation")

# ## Part 9: Troubleshooting Assembly Issues
#
# Identify and address common assembly problems

println("\n=== Troubleshooting Assembly Issues ===")

# ### Common Problems and Solutions
#
# Systematic approach to assembly problem diagnosis

println("--- Common Assembly Problems ---")

problems_solutions = Dict(
    "Low contiguity" => [
        "Increase read length or coverage",
        "Optimize assembly parameters",
        "Use scaffolding approaches",
        "Check for contamination"
    ],
    "Poor gene completeness" => [
        "Check assembly coverage",
        "Examine repeat resolution",
        "Validate gene prediction parameters",
        "Consider alternative assemblers"
    ],
    "High error rate" => [
        "Increase polishing iterations",
        "Check read quality",
        "Validate assembly parameters",
        "Consider consensus approaches"
    ],
    "Missing sequences" => [
        "Check for contamination filtering",
        "Examine coverage bias",
        "Validate input data quality",
        "Consider hybrid approaches"
    ]
)

for (problem, solutions) in problems_solutions
    println("$problem:")
    for solution in solutions
        println("  - $solution")
    end
    println()
end

# ## Summary
println("=== Assembly Validation Summary ===")
println("✓ Understanding multiple validation approaches")
println("✓ Implementing reference-based validation techniques")
println("✓ Applying reference-free validation methods")
println("✓ Functional validation through gene completeness")
println("✓ Comparative validation across multiple assemblies")
println("✓ Statistical validation and confidence assessment")
println("✓ Integrated quality scoring and reporting")
println("✓ Troubleshooting common assembly issues")

println("Results saved in: $(output_dir)")
println("\nNext: Tutorial 6 - Gene Annotation")

nothing
