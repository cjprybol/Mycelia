# Visualization Gallery

This gallery showcases the plotting and visualization capabilities of Mycelia across different bioinformatics analysis workflows. Each example includes the code to generate the visualization and explanation of biological interpretation.

## Overview

Mycelia provides comprehensive visualization capabilities for:
- **Data Quality Assessment** - Evaluate sequencing data quality
- **Sequence Analysis** - K-mer analysis and composition
- **Assembly Visualization** - Assembly statistics and comparisons  
- **Annotation Plots** - Gene features and functional summaries
- **Comparative Analysis** - Phylogenetics and pangenome visualization
- **Performance Monitoring** - Benchmark and resource usage plots

## Data Quality & Preprocessing

### FASTQ Quality Assessment

#### Per-Base Quality Scores
```julia
# Plot quality scores across read positions
quality_data = Mycelia.analyze_fastq_quality("reads.fastq")
Mycelia.plot_base_quality_scores(quality_data)
```
*Shows quality degradation patterns typical in sequencing data*

#### Read Length Distribution
```julia
# Visualize read length characteristics
length_dist = Mycelia.calculate_read_lengths("reads.fastq")
Mycelia.plot_read_length_distribution(length_dist)
```
*Identifies sequencing platform characteristics and potential issues*

#### GC Content Analysis
```julia
# Analyze sequence composition bias
gc_data = Mycelia.calculate_gc_content("reads.fastq")
Mycelia.plot_gc_content_distribution(gc_data)
```
*Detects contamination and sequencing bias*

#### Coverage Uniformity
```julia
# Assess coverage evenness across genome
coverage_data = Mycelia.calculate_coverage("reads.bam", "reference.fasta")
Mycelia.plot_coverage_uniformity(coverage_data)
```
*Identifies coverage bias and potential assembly issues*

## Sequence Analysis & K-mers

### K-mer Frequency Spectra
```julia
# Generate k-mer frequency histograms
kmer_counts = Mycelia.count_kmers("reads.fastq", k=21)
spectrum = Mycelia.kmer_frequency_spectrum(kmer_counts)
Mycelia.plot_kmer_spectrum(spectrum, title="21-mer Frequency Spectrum")
```
*Enables genome size estimation and error detection*

### K-mer Composition Heatmaps
```julia
# Visualize k-mer composition patterns
kmer_matrix = Mycelia.build_kmer_composition_matrix(sequences, k=4)
Mycelia.plot_kmer_heatmap(kmer_matrix, title="Tetranucleotide Composition")
```
*Reveals sequence composition patterns and contamination*

### Genome Size Estimation
```julia
# Plot k-mer based genome size estimation
estimation_data = Mycelia.estimate_genome_size_from_kmers(kmer_counts)
Mycelia.plot_genome_size_estimation(estimation_data)
```
*Provides independent genome size validation*

## Assembly Visualization

### Assembly Statistics
```julia
# Comprehensive assembly quality metrics
assembly_stats = Mycelia.evaluate_assembly("contigs.fasta")
Mycelia.plot_assembly_statistics(assembly_stats)
```
*Multi-panel plot showing N50, contig sizes, and quality metrics*

### Contig Size Distributions
```julia
# Visualize assembly contiguity
contig_lengths = Mycelia.get_contig_lengths("contigs.fasta")
Mycelia.plot_contig_distribution(contig_lengths, log_scale=true)
```
*Shows assembly fragmentation and contiguity*

### Assembly Comparison Dot Plots
```julia
# Compare assemblies or validate against reference
comparison = Mycelia.align_assemblies("assembly1.fasta", "assembly2.fasta")
Mycelia.plot_assembly_dotplot(comparison)
```
*Identifies structural differences and misassemblies*

### Coverage vs Length Analysis
```julia
# Analyze relationship between contig size and coverage
coverage_length_data = Mycelia.calculate_coverage_per_contig("reads.bam", "contigs.fasta")
Mycelia.plot_coverage_vs_length(coverage_length_data)
```
*Detects potential contamination and assembly artifacts*

## Annotation & Gene Features

### Genome Browser Tracks
```julia
# Create genome browser-style visualization
annotations = Mycelia.load_gff3("annotations.gff3")
sequence = Mycelia.load_fasta("genome.fasta")
Mycelia.plot_genome_browser(sequence, annotations, region="chr1:1000-5000")
```
*Interactive genome browser with multiple annotation tracks*

### Gene Density Plots
```julia
# Visualize gene distribution across genome
gene_positions = Mycelia.extract_gene_positions("annotations.gff3")
Mycelia.plot_gene_density(gene_positions, window_size=10000)
```
*Identifies gene-rich and gene-poor regions*

### Functional Annotation Summary
```julia
# Summarize functional categories
functional_data = Mycelia.analyze_functional_annotations("annotations.gff3")
Mycelia.plot_functional_categories(functional_data, plot_type="pie")
```
*Overview of gene functional classifications*

### GO Term Enrichment
```julia
# Visualize gene ontology enrichment
enrichment_results = Mycelia.perform_go_enrichment(gene_list, background)
Mycelia.plot_go_enrichment(enrichment_results, top_n=20)
```
*Identifies overrepresented biological processes*

## Comparative Genomics

### Phylogenetic Trees
```julia
# Multiple tree layouts and styles
tree = Mycelia.build_phylogenetic_tree(core_genes)
Mycelia.plot_phylogenetic_tree(tree, layout="circular", show_support=true)
```
*Publication-ready phylogenetic trees with branch support*

### Pangenome Heatmaps
```julia
# Gene presence/absence visualization
pangenome_matrix = Mycelia.build_pangenome_matrix(genomes)
Mycelia.plot_pangenome_heatmap(pangenome_matrix, cluster_rows=true, cluster_cols=true)
```
*Shows gene distribution across genomes with clustering*

### Synteny Analysis
```julia
# Visualize conserved gene order
synteny_data = Mycelia.calculate_synteny(genome1, genome2)
Mycelia.plot_synteny_dotplot(synteny_data, min_block_size=1000)
```
*Identifies chromosomal rearrangements and conserved regions*

### Population Structure
```julia
# Principal component analysis of genomic variation
pca_data = Mycelia.perform_genomic_pca(snp_matrix)
Mycelia.plot_population_structure(pca_data, color_by="population")
```
*Reveals population stratification and admixture*

### Pangenome Accumulation Curves
```julia
# Model pangenome size vs number of genomes
accumulation_data = Mycelia.calculate_pangenome_accumulation(genomes)
Mycelia.plot_pangenome_curves(accumulation_data, fit_model=true)
```
*Determines if pangenome is open or closed*

## Performance & Benchmarking

### Performance Scaling
```julia
# Benchmark performance across different scales
benchmark_data = Mycelia.run_scaling_benchmark(data_sizes, n_replicates=5)
Mycelia.plot_performance_scaling(benchmark_data, metric="time")
```
*Shows computational scaling characteristics*

### Memory Usage Monitoring
```julia
# Track memory usage during analysis
memory_data = Mycelia.monitor_memory_usage(analysis_function, inputs)
Mycelia.plot_memory_usage(memory_data, show_peak=true)
```
*Identifies memory bottlenecks and optimization opportunities*

### Resource Utilization Dashboard
```julia
# Multi-metric performance dashboard
resource_data = Mycelia.collect_resource_metrics(benchmark_results)
Mycelia.plot_resource_dashboard(resource_data)
```
*Comprehensive view of CPU, memory, and I/O usage*

### Accuracy vs Performance Trade-offs
```julia
# Compare different algorithms or parameters
tradeoff_data = Mycelia.analyze_accuracy_performance_tradeoff(algorithms, datasets)
Mycelia.plot_accuracy_performance(tradeoff_data)
```

## Customization and Export

## Gallery Organization

The visualization gallery is organized to match typical bioinformatics workflows:

1. **Start with data quality** - Essential first step
2. **Progress through analysis** - Following logical workflow
3. **End with interpretation** - Publication-ready visualizations
4. **Include performance monitoring** - For optimization

Each visualization includes:
- **Code example** - Complete, runnable code
- **Biological interpretation** - What the plot tells you
- **Customization options** - How to modify for your needs
- **Best practices** - When and how to use effectively

## Interactive Examples

For hands-on exploration, see the [Tutorial Notebooks](tutorials.md) which include executable visualization examples with real data.