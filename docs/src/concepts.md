# Bioinformatics Concepts and Tools

This guide explains the key bioinformatics concepts implemented in Mycelia and provides guidance on when to use each tool. It's designed for graduate-level researchers who need to understand both the biological context and computational approaches.

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Data Acquisition](#data-acquisition)
3. [Quality Control](#quality-control)
4. [Sequence Analysis](#sequence-analysis)
5. [Genome Assembly](#genome-assembly)
6. [Gene Annotation](#gene-annotation)
7. [Comparative Genomics](#comparative-genomics)
8. [Phylogenetics](#phylogenetics)
9. [Tool Selection Guide](#tool-selection-guide)

## Workflow Overview

Bioinformatics analyses typically follow a structured workflow. Understanding this flow helps you choose the right tools and interpret results correctly.

### Standard Genomic Analysis Pipeline

```
Raw Data → Quality Control → Feature Extraction → Assembly → Validation → Annotation → Comparative Analysis
```

### Decision Points

At each stage, you face key decisions:
- **Data Quality**: Is my data sufficient for analysis?
- **Method Selection**: Which algorithm fits my data type and research question?
- **Parameter Tuning**: How do I optimize for my specific dataset?
- **Result Interpretation**: What do the numbers mean biologically?

## Data Acquisition

### Concept: Data Sources and Types

**Biological Context**: Genomic data comes from various sources, each with different characteristics and limitations.

### Data Types

#### Sequencing Data
- **Illumina**: Short reads (150-300 bp), high accuracy, paired-end
- **PacBio/HiFi**: Long reads (10-25 kb), high accuracy, single-molecule
- **Oxford Nanopore**: Ultra-long reads (50+ kb), moderate accuracy, real-time

#### Reference Data
- **RefSeq**: Manually curated, high-quality reference sequences
- **GenBank**: Community-submitted sequences, variable quality
- **SRA**: Raw sequencing data from published studies

### When to Use Each

| Data Type | Use Case | Advantages | Limitations |
|-----------|----------|------------|-------------|
| Illumina | SNP calling, RNA-seq, metagenomics | High accuracy, cost-effective | Short reads, repetitive regions |
| PacBio HiFi | Genome assembly, structural variants | Long + accurate, span repeats | Higher cost, lower throughput |
| Nanopore | Real-time analysis, ultra-long reads | Longest reads, portable | Higher error rates |

### Mycelia Functions

```julia
# Download reference genomes
Mycelia.download_genome_by_accession(accession="NC_001422.1")

# Download complete assemblies
ncbi_genome_download_accession(accession="GCF_000819615.1")

# Simulate test data
simulate_random_genome(length=10000, gc_content=0.45)
simulate_hifi_reads(genome, coverage=20, error_rate=0.001)
```

## Quality Control

### Concept: Data Quality Assessment

**Biological Context**: Sequencing is not perfect. Raw data contains errors, biases, and artifacts that must be identified and corrected before analysis.

### Quality Metrics

#### Sequence Quality
- **Phred Scores**: Probability of base-calling error (Q20 = 1% error, Q30 = 0.1% error)
- **GC Content**: Should match expected organism profile
- **Length Distribution**: Indicates fragmentation or size selection

#### Coverage Metrics
- **Read Depth**: Number of reads covering each position
- **Uniformity**: Even coverage across genome
- **Duplication Rate**: PCR duplicates inflate coverage

### Quality Issues

#### Common Problems
- **Adapter Contamination**: Sequencing adapters in reads
- **Low Quality Ends**: Error-prone read ends
- **Overrepresented Sequences**: PCR bias or contamination
- **Length Bias**: Preferential sequencing of certain sizes

#### Solutions
- **Trimming**: Remove low-quality bases and adapters
- **Filtering**: Discard low-quality reads
- **Normalization**: Adjust for coverage bias
- **Deduplication**: Remove PCR duplicates

### When to Apply QC

| Analysis Type | QC Priority | Key Metrics |
|---------------|-------------|-------------|
| Genome Assembly | Critical | Coverage uniformity, read length |
| Variant Calling | High | Base quality, mapping quality |
| RNA-seq | High | rRNA contamination, 3' bias |
| Metagenomics | Medium | Host contamination, diversity |

### Mycelia Functions

```julia
# Analyze read quality
analyze_fastq_quality("reads.fastq")

# Calculate sequence statistics
calculate_sequence_stats(sequences)

# Quality visualization
plot_quality_distribution(quality_scores)
```

## Sequence Analysis

### Concept: K-mer Analysis

**Biological Context**: K-mers are subsequences of length k. They capture local sequence composition and are fundamental to many bioinformatics algorithms.

### K-mer Theory

#### Mathematical Foundation
- **K-mer Space**: 4^k possible DNA k-mers
- **Frequency Spectrum**: Distribution of k-mer frequencies
- **Coverage Estimation**: Genome size estimation from k-mer frequencies

#### Biological Interpretation
- **Repetitive Elements**: High-frequency k-mers indicate repeats
- **Sequencing Errors**: Low-frequency k-mers often represent errors
- **Genome Complexity**: Unique k-mers indicate complex regions

### K-mer Applications

#### Genome Assembly
- **Error Correction**: Remove k-mers with low frequency
- **Overlap Detection**: Find shared k-mers between reads
- **Graph Construction**: Build de Bruijn graphs

#### Quality Assessment
- **Coverage Estimation**: Estimate genome size and coverage
- **Contamination Detection**: Identify foreign DNA
- **Ploidy Estimation**: Detect polyploidy from k-mer frequencies

### Parameter Selection

#### Choosing K
- **Small K (k=15-21)**: Sensitive to errors, good for error correction
- **Medium K (k=25-31)**: Balanced sensitivity/specificity
- **Large K (k=35-51)**: Specific but may miss overlaps

#### Memory Considerations
- **Dense Matrices**: Store all possible k-mers (memory intensive)
- **Sparse Matrices**: Store only observed k-mers (memory efficient)
- **Counting Algorithms**: Bloom filters, Count-Min sketch

### When to Use K-mer Analysis

| Application | K-mer Size | Matrix Type | Use Case |
|-------------|------------|-------------|----------|
| Error Correction | 15-21 | Dense | Remove sequencing errors |
| Assembly | 25-31 | Sparse | Genome assembly |
| Repeat Detection | 31-51 | Sparse | Identify repetitive elements |
| Contamination | 21-31 | Sparse | Detect foreign DNA |

### Mycelia Functions

```julia
# Count k-mers
count_kmers("reads.fastq", k=21)

# Dense vs sparse counting
fasta_list_to_dense_kmer_counts(files, k=21)
fasta_list_to_sparse_kmer_counts(files, k=21)

# K-mer spectrum analysis
kmer_frequency_spectrum(kmer_counts)
```

## Genome Assembly

### Concept: Reconstructing Genomes

**Biological Context**: Sequencing breaks genomes into fragments. Assembly reconstructs the original genome by finding overlaps between fragments.

### Assembly Algorithms

#### Overlap-Layout-Consensus (OLC)
- **Overlap**: Find overlaps between all reads
- **Layout**: Order reads based on overlaps
- **Consensus**: Generate consensus sequence
- **Best For**: Long reads (PacBio, Nanopore)

#### De Bruijn Graph
- **K-mer Decomposition**: Break reads into k-mers
- **Graph Construction**: Connect overlapping k-mers
- **Path Finding**: Find Eulerian paths through graph
- **Best For**: Short reads (Illumina)

#### String Graph
- **Overlap Graph**: Nodes are reads, edges are overlaps
- **Transitivity Reduction**: Remove redundant edges
- **Contig Construction**: Find paths through graph
- **Best For**: Long reads with high accuracy

### Assembly Challenges

#### Repetitive Sequences
- **Problem**: Identical sequences confuse assembly
- **Solution**: Long reads that span repeats
- **Tools**: Repeat-aware assemblers, scaffolding

#### Heterozygosity
- **Problem**: Diploid genomes have two alleles
- **Solution**: Haplotype-aware assembly
- **Tools**: Trio binning, Hi-C scaffolding

#### Polyploidy
- **Problem**: Multiple copies of chromosomes
- **Solution**: Specialized polyploid assemblers
- **Tools**: Ploidy-aware algorithms

### Assembly Quality Metrics

#### Contiguity
- **N50**: Length where 50% of assembly is in contigs of this length or longer
- **L50**: Number of contigs containing 50% of assembly
- **Contig Count**: Total number of contigs (fewer is better)

#### Completeness
- **Genome Coverage**: Percentage of genome assembled
- **Gene Completeness**: Percentage of expected genes found
- **BUSCO Scores**: Conserved gene completeness

#### Accuracy
- **Base Accuracy**: Percentage of correct bases
- **Structural Accuracy**: Correct arrangement of sequences
- **Validation**: Comparison with reference genome

## Gene Annotation

### Concept: Identifying Functional Elements

**Biological Context**: Raw genome sequences are meaningless without annotation. Gene annotation identifies protein-coding genes, regulatory elements, and other functional sequences.

### Annotation Types

#### Structural Annotation
- **Gene Prediction**: Identify protein-coding genes
- **Exon-Intron Structure**: Define gene boundaries
- **Promoter Prediction**: Identify regulatory sequences
- **Repeat Annotation**: Classify repetitive elements

#### Functional Annotation
- **Protein Function**: Predict what proteins do
- **Pathway Mapping**: Assign genes to metabolic pathways
- **GO Terms**: Gene Ontology functional classification
- **Domain Annotation**: Identify protein domains

### Gene Prediction Methods

#### Ab Initio Prediction
- **Approach**: Use sequence signals (start codons, splice sites)
- **Advantages**: No external data required
- **Limitations**: Lower accuracy, misses non-canonical genes
- **Tools**: Prodigal, Pyrodigal, Augustus, GeneMark

#### Homology-Based Prediction
- **Approach**: Compare to known genes in databases
- **Advantages**: Higher accuracy for conserved genes
- **Limitations**: Misses novel genes, depends on database quality
- **Tools**: BLAST, DIAMOND, MMseqs2, MetaEuk

#### RNA-seq Guided Prediction
- **Approach**: Use transcriptome data to guide prediction
- **Advantages**: High accuracy, identifies expressed genes
- **Limitations**: Requires RNA-seq data, misses lowly expressed genes
- **Tools**: StringTie, Cufflinks, BRAKER

### Annotation Challenges

#### Eukaryotic Complexity
- **Alternative Splicing**: Multiple transcripts per gene
- **Pseudogenes**: Non-functional gene copies
- **Non-coding RNAs**: Functional RNAs that don't code for proteins

#### Prokaryotic Specifics
- **Operons**: Polycistronic transcripts
- **Overlapping Genes**: Genes sharing sequence
- **Horizontal Gene Transfer**: Genes from other species

### Quality Assessment

#### Annotation Completeness
- **Gene Density**: Number of genes per kb
- **Protein Completeness**: Percentage of complete proteins
- **Functional Coverage**: Percentage of genes with functional annotation

#### Annotation Accuracy
- **Validation**: Comparison with experimental data
- **Consistency**: Agreement between different methods
- **Benchmarking**: Comparison with reference annotations

### When to Use Different Approaches

| Genome Type | Method | Tools | Considerations |
|-------------|---------|-------|---------------|
| Bacterial | Ab initio | Prodigal / Pyrodigal | Simple gene structure |
| Viral | Specialized | Prodigal-gv | Viral coding patterns |
| Fungal | Hybrid | Augustus + BLAST | Moderate complexity |
| Eukaryotic metagenome | Homology-guided | MetaEuk | Fragmented contigs |
| Plant/Animal | RNA-seq guided | BRAKER | High complexity, alternative splicing |

## Comparative Genomics

### Concept: Comparing Genomes

**Biological Context**: Comparing genomes reveals evolutionary relationships, functional elements, and species-specific adaptations.

### Comparative Approaches

#### Pairwise Comparison
- **Synteny**: Conserved gene order between species
- **Orthology**: Corresponding genes in different species
- **Whole Genome Alignment**: Align entire genomes
- **Structural Variation**: Differences in genome structure

#### Pangenome Analysis
- **Core Genome**: Genes present in all individuals
- **Accessory Genome**: Genes present in some individuals
- **Unique Genes**: Genes present in single individuals
- **Gene Gain/Loss**: Evolution of gene content

### Pangenome Construction

#### Graph-Based Approaches
- **Gene Graphs**: Nodes are genes, edges are relationships
- **Sequence Graphs**: Nodes are sequences, edges are overlaps
- **Variation Graphs**: Represent genetic variation as graphs

#### Clustering Approaches
- **Sequence Similarity**: Group similar sequences
- **Synteny**: Group genes with conserved context
- **Functional Similarity**: Group genes with similar functions

### Evolutionary Analysis

#### Phylogenetic Trees
- **Species Trees**: Evolutionary relationships between species
- **Gene Trees**: Evolution of individual genes
- **Reconciliation**: Compare gene and species trees

#### Selection Analysis
- **Positive Selection**: Genes under selective pressure
- **Purifying Selection**: Genes with functional constraints
- **Neutral Evolution**: Genes evolving without selection

## Phylogenetics

### Concept: Evolutionary Relationships

**Biological Context**: Phylogenetics reconstructs evolutionary history from molecular data, revealing how species are related and when they diverged.

### Tree Construction Methods

#### Distance-Based Methods
- **UPGMA**: Assumes molecular clock
- **Neighbor-Joining**: Doesn't assume molecular clock
- **Minimum Evolution**: Finds tree with shortest total branch length

#### Character-Based Methods
- **Maximum Parsimony**: Minimizes evolutionary changes
- **Maximum Likelihood**: Finds most likely tree given data
- **Bayesian Inference**: Incorporates prior knowledge

### Molecular Evolution Models

#### Nucleotide Substitution Models
- **JC69**: All substitutions equally likely
- **K80**: Different rates for transitions/transversions
- **GTR**: General time-reversible model

#### Protein Evolution Models
- **Poisson**: All amino acid changes equally likely
- **JTT**: Jones-Taylor-Thornton model
- **LG**: Le-Gascuel model

### Tree Evaluation

#### Support Values
- **Bootstrap**: Resampling support for branches
- **Posterior Probability**: Bayesian support
- **SH-like aLRT**: Likelihood ratio test

#### Tree Comparison
- **Robinson-Foulds**: Topological distance
- **Likelihood Ratio**: Statistical comparison
- **Consensus Trees**: Combine multiple trees

### Applications

#### Taxonomy
- **Species Identification**: Place unknown species
- **Classification**: Revise taxonomic relationships
- **Diversity Assessment**: Measure evolutionary diversity

#### Epidemiology
- **Outbreak Tracking**: Trace disease transmission
- **Vaccine Design**: Understand pathogen evolution
- **Drug Resistance**: Track resistance evolution

### Choosing the Right Analysis

#### Data-Driven Decisions

**Consider Your Data Type**
- **Illumina**: Short read optimized tools
- **PacBio/Nanopore**: Long read specialized tools
- **Hybrid**: Tools supporting multiple data types

**Consider Your Organism**
- **Bacteria**: Simpler tools, less memory
- **Eukaryotes**: Complex tools, more resources
- **Viruses**: Specialized viral tools

**Consider Your Resources**
- **Local Machine**: Memory/CPU constraints
- **HPC Cluster**: Parallel processing available
- **Cloud**: Scalable but cost considerations

#### Question-Driven Decisions

**Assembly Quality**
- **Draft Assembly**: Fast tools, lower quality
- **Reference Quality**: Comprehensive tools, higher quality
- **Comparison**: Multiple assemblers, best result

**Annotation Depth**
- **Basic Annotation**: Fast gene prediction
- **Comprehensive**: Multiple evidence types
- **Comparative**: Cross-species evidence

**Analysis Scope**
- **Single Genome**: Individual analysis tools
- **Population**: Population genetics tools
- **Comparative**: Multi-genome tools

### Performance Considerations

#### Memory Requirements
- **K-mer Analysis**: Exponential with k size
- **Assembly**: Linear with genome size
- **Annotation**: Depends on database size

#### Time Complexity
- **Data Download**: Network dependent
- **Quality Control**: Linear with data size
- **Assembly**: Depends on algorithm and data
- **Annotation**: Depends on database searches

#### Accuracy Trade-offs
- **Speed vs Accuracy**: Fast tools may sacrifice accuracy
- **Sensitivity vs Specificity**: Detect more vs fewer false positives
- **Completeness vs Contiguity**: More complete vs fewer contigs

### Best Practices

#### Workflow Design
1. **Start Simple**: Use basic tools first
2. **Iterate**: Improve based on results
3. **Validate**: Check results at each step
4. **Document**: Record parameters and decisions

#### Quality Control
1. **Check Inputs**: Validate data quality
2. **Monitor Progress**: Track intermediate results
3. **Evaluate Outputs**: Assess final results
4. **Compare Methods**: Use multiple approaches

#### Reproducibility
1. **Version Control**: Track software versions
2. **Parameter Recording**: Document all settings
3. **Environment Management**: Use containers/environments
4. **Data Provenance**: Track data sources