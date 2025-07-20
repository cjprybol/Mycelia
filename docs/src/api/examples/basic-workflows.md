# Basic Workflows

Complete examples of common bioinformatics analysis workflows using Mycelia functions. These examples show how to combine functions from different modules to accomplish typical research tasks.

## Overview

These workflows demonstrate:
- **Complete analysis pipelines** from start to finish
- **Function integration** across different modules
- **Parameter selection** for common use cases
- **Error handling** and quality control
- **Result interpretation** and next steps

## Workflow 1: Bacterial Genome Assembly

Complete bacterial genome analysis from raw reads to annotated assembly.

### Step 1: Data Acquisition
```julia
import Mycelia

# Download reference genome for comparison
reference = Mycelia.download_genome_by_accession("NC_000913.3")  # E. coli K-12

# Or use your own sequencing data
reads_file = "bacterial_reads.fastq"
```

### Step 2: Quality Control
```julia
# Assess initial data quality
println("=== Initial Quality Assessment ===")
initial_quality = Mycelia.analyze_fastq_quality(reads_file)
println("Total reads: $(initial_quality.n_reads)")
println("Mean quality: $(initial_quality.mean_quality)")
println("Mean length: $(initial_quality.mean_length)")

# Filter low-quality reads
println("\n=== Quality Filtering ===")
filtered_reads = Mycelia.filter_by_quality(
    reads_file,
    min_quality=25,      # Q25 threshold
    min_length=1000,     # Minimum 1kb reads
    max_n_percent=5      # Maximum 5% N's
)

Mycelia.write_fastq("filtered_reads.fastq", filtered_reads)
println("Filtered reads: $(length(filtered_reads))")
```

### Step 3: K-mer Analysis
```julia
println("\n=== K-mer Analysis ===")
# Count k-mers for genome size estimation
kmer_counts = Mycelia.count_kmers("filtered_reads.fastq", k=21)
spectrum = Mycelia.kmer_frequency_spectrum(kmer_counts)

# Estimate genome size
genome_size = Mycelia.estimate_genome_size_from_kmers(kmer_counts)
println("Estimated genome size: $(genome_size.size) bp")
println("Estimated coverage: $(genome_size.coverage)x")

# Plot k-mer spectrum
Mycelia.plot_kmer_spectrum(spectrum, title="21-mer Frequency Spectrum")
```

### Step 4: Genome Assembly
```julia
println("\n=== Genome Assembly ===")
# Assemble genome using hifiasm (for HiFi reads) or adjust for your data type
assembly_result = Mycelia.assemble_genome(
    "filtered_reads.fastq",
    assembler="hifiasm",     # Use "spades" for Illumina
    output_dir="assembly",
    threads=8,
    min_contig_length=500
)

println("Assembly completed: $(assembly_result.contigs)")
```

### Step 5: Assembly Validation
```julia
println("\n=== Assembly Validation ===")
# Calculate assembly statistics
assembly_stats = Mycelia.calculate_assembly_stats(assembly_result.contigs)
println("Assembly Statistics:")
println("  Contigs: $(assembly_stats.n_contigs)")
println("  Total length: $(assembly_stats.total_length) bp")
println("  N50: $(assembly_stats.n50)")
println("  Largest contig: $(assembly_stats.largest_contig)")

# Validate assembly quality
validation_result = Mycelia.validate_assembly(
    assembly_result.contigs,
    "filtered_reads.fastq",
    reference_genome=reference
)

println("Assembly Quality:")
println("  Completeness: $(validation_result.completeness)%")
println("  Accuracy: $(validation_result.accuracy)%")
```

### Step 6: Gene Annotation
```julia
println("\n=== Gene Annotation ===")
# Predict genes
predicted_genes = Mycelia.predict_genes(
    assembly_result.contigs,
    method="prodigal",
    genetic_code="standard"
)

println("Predicted genes: $(length(predicted_genes))")

# Functional annotation
functional_annotations = Mycelia.annotate_functions(
    predicted_genes,
    database="uniprot",
    evalue_threshold=1e-5
)

println("Functionally annotated genes: $(Mycelia.count_annotated(functional_annotations)))")

# Save annotations
Mycelia.write_gff3("bacterial_genome.gff3", predicted_genes, functional_annotations)
```

### Step 7: Results Summary
```julia
println("\n=== Analysis Summary ===")
println("Workflow completed successfully!")
println("Files generated:")
println("  - filtered_reads.fastq: Quality-controlled reads")
println("  - assembly/contigs.fasta: Assembled genome")
println("  - bacterial_genome.gff3: Gene annotations")
println("  - assembly_stats.json: Assembly statistics")
```

## Workflow 2: Comparative Genomics Analysis

Compare multiple bacterial genomes to build a pangenome and phylogenetic tree.

### Step 1: Data Preparation
```julia
# Download multiple related genomes
species_genomes = [
    "GCF_000005825.2",  # E. coli K-12 MG1655
    "GCF_000009605.1",  # Salmonella enterica
    "GCF_000027325.1",  # Yersinia pestis
    "GCF_000006945.2",  # Shigella flexneri
]

println("=== Downloading Genomes ===")
genome_files = []
for accession in species_genomes
    result = Mycelia.ncbi_genome_download_accession(accession)
    push!(genome_files, result.genome)
    println("Downloaded: $accession")
end
```

### Step 2: Gene Prediction
```julia
println("\n=== Gene Prediction ===")
all_genes = []
for (i, genome_file) in enumerate(genome_files)
    genes = Mycelia.predict_genes(genome_file, method="prodigal")
    push!(all_genes, genes)
    println("Genome $i: $(length(genes)) genes predicted")
end
```

### Step 3: Pangenome Construction
```julia
println("\n=== Pangenome Construction ===")
# Build pangenome from all genomes
pangenome = Mycelia.build_pangenome(
    all_genes,
    similarity_threshold=0.9,
    coverage_threshold=0.8
)

println("Pangenome Statistics:")
println("  Core genes: $(pangenome.core_genes)")
println("  Accessory genes: $(pangenome.accessory_genes)")
println("  Unique genes: $(pangenome.unique_genes)")
println("  Total gene families: $(pangenome.total_families)")

# Visualize pangenome
Mycelia.plot_pangenome_heatmap(pangenome, title="Gene Presence/Absence Matrix")
```

### Step 4: Phylogenetic Analysis
```julia
println("\n=== Phylogenetic Analysis ===")
# Extract core genes for phylogeny
core_gene_alignments = Mycelia.extract_core_gene_alignments(pangenome)

# Build phylogenetic tree
phylo_tree = Mycelia.build_phylogenetic_tree(
    core_gene_alignments,
    method="ml",
    model="GTR+G",
    bootstrap=100
)

println("Phylogenetic tree constructed with $(Mycelia.get_bootstrap_support(phylo_tree)) average support")

# Visualize tree
Mycelia.plot_phylogenetic_tree(phylo_tree, layout="rectangular", show_support=true)
```

### Step 5: Functional Analysis
```julia
println("\n=== Functional Analysis ===")
# Analyze functional categories
functional_analysis = Mycelia.analyze_pangenome_functions(pangenome)

println("Functional Distribution:")
for (category, count) in functional_analysis.categories
    println("  $category: $count genes")
end

# Identify core vs accessory functional differences
functional_comparison = Mycelia.compare_core_accessory_functions(pangenome)
```

## Workflow 3: Quality Control Pipeline

Comprehensive quality control workflow for different sequencing platforms.

### HiFi Sequencing Data
```julia
println("=== HiFi Quality Control ===")

# HiFi-specific quality assessment
hifi_quality = Mycelia.assess_hifi_quality("hifi_reads.fastq")
println("HiFi Quality Metrics:")
println("  Mean accuracy: $(hifi_quality.mean_accuracy)")
println("  Mean length: $(hifi_quality.mean_length)")
println("  Length N50: $(hifi_quality.length_n50)")

# HiFi-optimized filtering
hifi_filtered = Mycelia.filter_hifi_reads(
    "hifi_reads.fastq",
    min_accuracy=0.99,
    min_length=5000,
    max_length=30000
)

println("HiFi reads after filtering: $(length(hifi_filtered))")
```

### Illumina Sequencing Data
```julia
println("\n=== Illumina Quality Control ===")

# Paired-end Illumina data
illumina_quality = Mycelia.assess_illumina_quality("illumina_R1.fastq", "illumina_R2.fastq")

# Remove adapters and low-quality bases
cleaned_reads = Mycelia.preprocess_illumina_reads(
    "illumina_R1.fastq", "illumina_R2.fastq",
    adapter_removal=true,
    quality_trimming=true,
    min_quality=20,
    min_length=50
)

Mycelia.write_fastq("illumina_R1_clean.fastq", cleaned_reads.read1)
Mycelia.write_fastq("illumina_R2_clean.fastq", cleaned_reads.read2)
```

### Nanopore Sequencing Data
```julia
println("\n=== Nanopore Quality Control ===")

# Nanopore-specific assessment
nanopore_quality = Mycelia.assess_nanopore_quality("nanopore_reads.fastq")

# Filter based on quality and length
nanopore_filtered = Mycelia.filter_nanopore_reads(
    "nanopore_reads.fastq",
    min_quality=7,       # Lower threshold for Nanopore
    min_length=500,
    max_length=100000
)
```

## Workflow 4: Assembly Optimization

Optimize assembly parameters for best results.

### Parameter Testing
```julia
println("=== Assembly Parameter Optimization ===")

# Test different k-mer sizes
k_values = [19, 21, 25, 31, 35]
assembly_results = []

for k in k_values
    println("Testing k=$k...")
    result = Mycelia.assemble_genome(
        "reads.fastq",
        assembler="hifiasm",
        k=k,
        output_dir="assembly_k$k"
    )
    
    stats = Mycelia.calculate_assembly_stats(result.contigs)
    push!(assembly_results, (k=k, n50=stats.n50, contigs=stats.n_contigs))
end

# Find optimal k-mer size
optimal_k = Mycelia.find_optimal_assembly_parameters(assembly_results)
println("Optimal k-mer size: $(optimal_k.k)")
```

### Multi-Assembler Comparison
```julia
println("\n=== Multi-Assembler Comparison ===")

assemblers = ["hifiasm", "canu", "flye"]
assembler_results = []

for assembler in assemblers
    println("Running $assembler...")
    result = Mycelia.assemble_genome(
        "reads.fastq",
        assembler=assembler,
        output_dir="assembly_$assembler"
    )
    
    validation = Mycelia.validate_assembly(result.contigs, "reads.fastq")
    push!(assembler_results, (
        assembler=assembler,
        n50=validation.n50,
        completeness=validation.completeness
    ))
end

# Select best assembler
best_assembler = Mycelia.select_best_assembly(assembler_results)
println("Best assembler: $(best_assembler.assembler)")
```

## Workflow 5: Contamination Detection and Removal

Comprehensive contamination screening and removal pipeline.

### Multi-Source Contamination Screening
```julia
println("=== Contamination Screening ===")

# Screen for multiple contamination sources
contamination_results = Mycelia.screen_all_contamination(
    "reads.fastq",
    host_genome="human_genome.fasta",
    vector_db="vector_database.fasta",
    adapter_db="adapter_sequences.fasta"
)

println("Contamination Summary:")
println("  Host contamination: $(contamination_results.host_rate)%")
println("  Vector contamination: $(contamination_results.vector_rate)%")
println("  Adapter contamination: $(contamination_results.adapter_rate)%")

# Remove all contamination
clean_reads = Mycelia.remove_all_contamination(
    "reads.fastq",
    contamination_results
)

Mycelia.write_fastq("decontaminated_reads.fastq", clean_reads)
```

## Error Handling and Robustness

### Robust Pipeline with Error Handling
```julia
function robust_genome_analysis(reads_file::String)
    try
        # Quality control with validation
        if !isfile(reads_file)
            error("Input file not found: $reads_file")
        end
        
        quality_data = Mycelia.analyze_fastq_quality(reads_file)
        if quality_data.mean_quality < 15
            @warn "Low quality data detected (Q$(quality_data.mean_quality))"
        end
        
        # Assembly with retry logic
        assembly_result = nothing
        for attempt in 1:3
            try
                assembly_result = Mycelia.assemble_genome(reads_file)
                break
            catch e
                @warn "Assembly attempt $attempt failed: $e"
                if attempt == 3
                    rethrow(e)
                end
            end
        end
        
        # Validation with quality checks
        validation = Mycelia.validate_assembly(assembly_result.contigs, reads_file)
        if validation.completeness < 90
            @warn "Low assembly completeness: $(validation.completeness)%"
        end
        
        return assembly_result
        
    catch e
        @error "Analysis failed" exception=(e, catch_backtrace())
        return nothing
    end
end
```

## Performance Optimization

### Memory-Efficient Large File Processing
```julia
function process_large_dataset(large_file::String)
    # Check available memory
    available_memory = Mycelia.get_available_memory_gb()
    if available_memory < 8
        @warn "Limited memory available: $(available_memory)GB"
    end
    
    # Use streaming for large files
    file_size_gb = filesize(large_file) / 1024^3
    if file_size_gb > available_memory / 2
        println("Using streaming processing for $(file_size_gb)GB file")
        return Mycelia.stream_kmer_counting(large_file, k=21, chunk_size=50000)
    else
        return Mycelia.count_kmers(large_file, k=21)
    end
end
```

## Next Steps

After completing these basic workflows, you can:

1. **Explore advanced techniques** in [Advanced Usage](advanced-usage.md)
2. **Customize parameters** using the [Parameter Guide](../quick-reference/parameter-guide.md)
3. **Integrate external tools** following [Tool Integration](../workflows/visualization.md)
4. **Scale to HPC systems** with [Deployment Guide](../../../deployment-guide.md)

## See Also

- [Function Index](../quick-reference/function-index.md) - Complete function reference
- [Workflow-Specific Guides](../workflows/) - Detailed workflow documentation
- [Tutorials](../../tutorials.md) - Step-by-step learning materials
- [Troubleshooting Guide](../../troubleshooting.md) - Common issues and solutions