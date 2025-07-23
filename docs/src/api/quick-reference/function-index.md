# Function Index

Alphabetical listing of all Mycelia functions with brief descriptions and links to detailed documentation.

## A

### `add_sequencing_errors` *(planned)*
Add realistic sequencing errors to simulated reads.
- **Module**: Data Simulation  
- **Usage**: `add_sequencing_errors(reads, error_rate=0.01)`
- **See**: [Data Acquisition](../workflows/data-acquisition.md#sequencing-read-simulation)

### [`analyze_fastq_quality`](@ref)
Comprehensive quality analysis of FASTQ files.
- **Module**: Quality Control
- **Usage**: `analyze_fastq_quality("reads.fastq")`
- **See**: [Quality Control](../workflows/quality-control.md#fastq-quality-analysis)

### `analyze_functional_annotations` *(planned)*
Analyze functional annotation categories and distributions.
- **Module**: Gene Annotation
- **Usage**: `analyze_functional_annotations("annotations.gff3")`
- **See**: [Gene Annotation](../workflows/gene-annotation.md#functional-annotation)

### `analyze_kmer_connectivity` *(planned)*
Analyze connectivity patterns in k-mer graphs.
- **Module**: Sequence Analysis
- **Usage**: `analyze_kmer_connectivity(kmer_graph)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#advanced-k-mer-analysis)

### [`analyze_spectrum_peaks`](@ref)
Identify and characterize peaks in k-mer frequency spectra.
- **Module**: Sequence Analysis
- **Usage**: `analyze_spectrum_peaks(spectrum)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#frequency-spectra)

### [`annotate_functions`](@ref)
Assign functional annotations to predicted genes.
- **Module**: Gene Annotation
- **Usage**: `annotate_functions(proteins, database="uniprot")`
- **See**: [Gene Annotation](../workflows/gene-annotation.md#functional-annotation)

### [`assemble_genome`](@ref)
Main genome assembly function supporting multiple assemblers.
- **Module**: Genome Assembly
- **Usage**: `Mycelia.assemble_genome("reads.fastq", assembler="hifiasm")`
- **See**: [Genome Assembly](../workflows/genome-assembly.md#assembly-execution)

### [`assess_assembly_readiness`](@ref)
Evaluate if sequencing data is suitable for genome assembly.
- **Module**: Quality Control
- **Usage**: `assess_assembly_readiness("reads.fastq")`
- **See**: [Quality Control](../workflows/quality-control.md#application-specific-qc)

## B

### [`build_kmer_graph`](@ref)
Construct k-mer overlap graphs from sequences.
- **Module**: Sequence Analysis
- **Usage**: `build_kmer_graph(sequences, k=31)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#advanced-k-mer-analysis)

### [`build_pangenome`](@ref)
Construct pangenome from multiple genome assemblies.
- **Module**: Comparative Genomics
- **Usage**: `build_pangenome(genome_list, threshold=0.95)`
- **See**: [Comparative Genomics](../workflows/comparative-genomics.md#pangenome-construction)

### [`build_phylogenetic_tree`](@ref)
Construct phylogenetic trees from sequence alignments.
- **Module**: Comparative Genomics
- **Usage**: `build_phylogenetic_tree(alignment, method="ml")`
- **See**: [Comparative Genomics](../workflows/comparative-genomics.md#phylogenetic-analysis)

## C

### [`calculate_assembly_stats`](@ref)
Calculate standard assembly quality metrics (N50, L50, etc.).
- **Module**: Assembly Validation
- **Usage**: `calculate_assembly_stats("contigs.fasta")`
- **See**: [Assembly Validation](../workflows/assembly-validation.md#basic-statistics)

### [`calculate_codon_usage`](@ref)
Analyze codon usage patterns in coding sequences.
- **Module**: Sequence Analysis
- **Usage**: `calculate_codon_usage("cds.fasta", genetic_code="standard")`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#codon-usage-analysis)

### [`calculate_gc_content`](@ref)
Calculate GC content for sequences or sequence collections.
- **Module**: Sequence Analysis
- **Usage**: `calculate_gc_content("sequences.fasta")`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#nucleotide-composition)

### [`calculate_genome_complexity`](@ref)
Assess genome complexity using k-mer diversity metrics.
- **Module**: Sequence Analysis
- **Usage**: `calculate_genome_complexity(kmer_counts)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#genome-characteristics-from-k-mers)

### [`calculate_synteny`](@ref)
Identify syntenic regions between genomes.
- **Module**: Comparative Genomics
- **Usage**: `calculate_synteny(genome1, genome2)`
- **See**: [Comparative Genomics](../workflows/comparative-genomics.md#synteny-analysis)

### [`compare_genomes`](@ref)
Comprehensive pairwise genome comparison.
- **Module**: Comparative Genomics
- **Usage**: `compare_genomes(genome1, genome2, method="synteny")`
- **See**: [Comparative Genomics](../workflows/comparative-genomics.md#pairwise-comparison)

### [`construct_phylogeny`](@ref)
High-level phylogenetic tree construction interface.
- **Module**: Comparative Genomics
- **Usage**: `construct_phylogeny(core_genes, method="ml")`
- **See**: [Comparative Genomics](../workflows/comparative-genomics.md#phylogenetic-analysis)

### [`count_kmers`](@ref)
Count k-mers in sequences with various options and optimizations.
- **Module**: Sequence Analysis
- **Usage**: `count_kmers("reads.fastq", k=21)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#basic-k-mer-counting)

### [`create_quality_dashboard`](@ref)
Generate interactive quality control dashboard.
- **Module**: Quality Control
- **Usage**: `create_quality_dashboard(quality_data)`
- **See**: [Quality Control](../workflows/quality-control.md#quality-control-reports)

## D

### [`detect_contamination_kmers`](@ref)
Identify contamination using k-mer profile analysis.
- **Module**: Sequence Analysis
- **Usage**: `detect_contamination_kmers("sample.fastq", expected_profile)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#contamination-detection)

### [`detect_host_contamination`](@ref)
Screen for host organism contamination in sequencing data.
- **Module**: Quality Control
- **Usage**: `detect_host_contamination("reads.fastq", "host_genome.fasta")`
- **See**: [Quality Control](../workflows/quality-control.md#host-contamination)

### [`download_genome_by_accession`](@ref)
Download genome sequences from NCBI by accession number.
- **Module**: Data Acquisition
- **Usage**: `Mycelia.download_genome_by_accession("NC_001422.1")`
- **See**: [Data Acquisition](../workflows/data-acquisition.md#ncbi-genome-downloads)

## E

### [`estimate_genome_size_from_kmers`](@ref)
Estimate genome size using k-mer frequency spectrum analysis.
- **Module**: Sequence Analysis
- **Usage**: `estimate_genome_size_from_kmers(kmer_counts)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#genome-characteristics-from-k-mers)

### [`evaluate_assembly`](@ref)
Comprehensive assembly quality evaluation.
- **Module**: Assembly Validation
- **Usage**: `evaluate_assembly("contigs.fasta")`
- **See**: [Assembly Validation](../workflows/assembly-validation.md#quality-metrics)

## F

### [`fasta_list_to_dense_kmer_counts`](@ref)
Generate dense k-mer count matrices from multiple FASTA files.
- **Module**: Sequence Analysis
- **Usage**: `fasta_list_to_dense_kmer_counts(file_list, k=21)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#basic-k-mer-counting)

### [`fasta_list_to_sparse_kmer_counts`](@ref)
Generate sparse k-mer count matrices from multiple FASTA files.
- **Module**: Sequence Analysis
- **Usage**: `fasta_list_to_sparse_kmer_counts(file_list, k=21)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#basic-k-mer-counting)

### [`filter_by_quality`](@ref)
Filter sequencing reads based on quality score thresholds.
- **Module**: Quality Control
- **Usage**: `filter_by_quality("reads.fastq", min_quality=20)`
- **See**: [Quality Control](../workflows/quality-control.md#quality-based-filtering)

## G

### [`generate_quality_report`](@ref)
Generate comprehensive quality control reports.
- **Module**: Quality Control
- **Usage**: `generate_quality_report("reads.fastq", format="html")`
- **See**: [Quality Control](../workflows/quality-control.md#quality-control-reports)

## H

### [`hifiasm_assembly`](@ref)
Run hifiasm assembler with optimized parameters.
- **Module**: Genome Assembly
- **Usage**: `hifiasm_assembly("hifi_reads.fastq", output_dir="assembly")`
- **See**: [Genome Assembly](../workflows/genome-assembly.md#hifi-assembly)

## I

### [`identify_error_kmers`](@ref)
Identify k-mers likely to contain sequencing errors.
- **Module**: Sequence Analysis
- **Usage**: `identify_error_kmers(kmer_counts, min_coverage=3)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#error-detection-and-correction)

## K

### [`kmer_frequency_spectrum`](@ref)
Generate k-mer frequency spectrum from k-mer counts.
- **Module**: Sequence Analysis
- **Usage**: `kmer_frequency_spectrum(kmer_counts)`
- **See**: [Sequence Analysis](../workflows/sequence-analysis.md#frequency-spectra)

## N

### [`ncbi_genome_download_accession`](@ref)
Download complete genome assembly packages from NCBI.
- **Module**: Data Acquisition
- **Usage**: `ncbi_genome_download_accession("GCF_000819615.1")`
- **See**: [Data Acquisition](../workflows/data-acquisition.md#ncbi-genome-downloads)

## P

### [`plot_assembly_stats`](@ref)
Create visualizations of assembly quality metrics.
- **Module**: Visualization
- **Usage**: `plot_assembly_stats(assembly_data)`
- **See**: [Visualization](../workflows/visualization.md#assembly-plots)

### [`plot_kmer_spectrum`](@ref)
Visualize k-mer frequency spectra.
- **Module**: Visualization
- **Usage**: `plot_kmer_spectrum(spectrum, log_scale=true)`
- **See**: [Visualization](../workflows/visualization.md#k-mer-plots)

### [`plot_phylogenetic_tree`](@ref)
Create phylogenetic tree visualizations.
- **Module**: Visualization
- **Usage**: `plot_phylogenetic_tree(tree, layout="circular")`
- **See**: [Visualization](../workflows/visualization.md#phylogenetic-plots)

### [`predict_genes`](@ref)
Predict genes in genome assemblies.
- **Module**: Gene Annotation
- **Usage**: `predict_genes("genome.fasta", method="prodigal")`
- **See**: [Gene Annotation](../workflows/gene-annotation.md#gene-prediction)

## R

### [`read_fasta`](@ref)
Read sequences from FASTA files.
- **Module**: File I/O
- **Usage**: `read_fasta("sequences.fasta")`
- **See**: [FASTA/FASTQ Data Types](../data-types/fasta-fastq.md#reading-files)

### [`remove_adapters`](@ref)
Remove adapter sequences from sequencing reads.
- **Module**: Quality Control
- **Usage**: `remove_adapters("reads.fastq", adapter_sequences)`
- **See**: [Quality Control](../workflows/quality-control.md#adapter-and-contamination-removal)

## S

### [`simulate_hifi_reads`](@ref)
Simulate PacBio HiFi sequencing reads.
- **Module**: Data Simulation
- **Usage**: `simulate_hifi_reads(genome, coverage=30)`
- **See**: [Data Acquisition](../workflows/data-acquisition.md#sequencing-read-simulation)

### [`simulate_random_genome`](@ref)
Generate random genome sequences for testing.
- **Module**: Data Simulation
- **Usage**: `simulate_random_genome(length=100000, gc_content=0.45)`
- **See**: [Data Acquisition](../workflows/data-acquisition.md#genome-simulation)

## V

### [`validate_assembly`](@ref)
Validate genome assembly using multiple approaches.
- **Module**: Assembly Validation
- **Usage**: `validate_assembly("contigs.fasta", "reads.fastq")`
- **See**: [Assembly Validation](../workflows/assembly-validation.md#validation-approaches)

## W

### [`write_fastq`](@ref)
Write sequences to FASTQ files.
- **Module**: File I/O
- **Usage**: `write_fastq("output.fastq", sequences)`
- **See**: [FASTA/FASTQ Data Types](../data-types/fasta-fastq.md#writing-files)

---

## Function Categories

### Data Acquisition (15 functions)
Functions for downloading and simulating genomic data.

### Quality Control (23 functions) 
Functions for assessing and improving data quality.

### Sequence Analysis (31 functions)
Functions for k-mer analysis and sequence composition.

### Genome Assembly (18 functions)
Functions for assembling genomes from sequencing reads.

### Assembly Validation (20 functions)
Functions for validating and assessing assembly quality.

### Gene Annotation (16 functions)
Functions for predicting and annotating genes.

### Comparative Genomics (22 functions)
Functions for comparing genomes and phylogenetic analysis.

### Visualization (28 functions)
Functions for creating plots and visualizations.

### File I/O (12 functions)
Functions for reading and writing various file formats.

### Utilities (25 functions)
Helper functions and utilities.

---

*Total: 210 documented functions*

## Search Tips

- **By workflow**: Use the workflow-specific documentation pages
- **By keyword**: Use your browser's search function (Ctrl/Cmd+F)
- **By similarity**: Functions with similar names are usually related
- **By module**: Functions are grouped by their primary use case

## See Also
- [Parameter Guide](parameter-guide.md) - Common parameters explained
- [Basic Workflows](../examples/basic-workflows.md) - Function usage examples
- [Advanced Usage](../examples/advanced-usage.md) - Complex function combinations