"""
Mycelia Graph Type Tutorials

This tutorial demonstrates the usage of all 6 graph types in Mycelia:
1. N-gram Graphs - For unicode text analysis
2. K-mer Graphs - For DNA/RNA/protein sequence analysis
3. Qualmer Graphs - For quality-aware sequence analysis
4. String Graphs - For simplified N-gram analysis
5. FASTA Graphs - For variable-length BioSequence analysis
6. FASTQ Graphs - For quality-aware BioSequence analysis

Each section shows practical examples with real data.
"""

import Mycelia
import FASTX
import BioSequences

println("="^80)
println("MYCELIA GRAPH TYPE TUTORIALS")
println("="^80)

# =============================================================================
# Tutorial 1: N-gram Graphs - Unicode Text Analysis
# =============================================================================

println("\n1. N-GRAM GRAPHS - Unicode Text Analysis")
println("-"^50)

# Example: Analyzing text patterns
text = "The quick brown fox jumps over the lazy dog"
println("Input text: \"$text\"")

# Create N-gram graph
ngram_graph = Mycelia.string_to_ngram_graph(s=text, n=3)
println("Number of 3-grams: $(length(ngram_graph.vertex_labels))")

# Show some N-grams
ngrams = collect(values(ngram_graph.vertex_labels))
println("Example 3-grams: $(ngrams[1:min(5, length(ngrams))])")

# Use case: Text compression, pattern detection, linguistic analysis
println("Use cases: Text compression, pattern detection, linguistic analysis")

# =============================================================================
# Tutorial 2: K-mer Graphs - BioSequence Analysis
# =============================================================================

println("\n2. K-MER GRAPHS - BioSequence Analysis")
println("-"^50)

# Example: DNA sequence analysis
dna_sequence = "ATCGATCGATCGATCGTAGCTAGCTAGCT"
println("DNA sequence: $dna_sequence")

# Create FASTA record
dna_record = FASTX.FASTA.Record("sample_dna", dna_sequence)
dna_records = [dna_record]

# Create k-mer graph
dna_kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, dna_records)
dna_kmers = collect(values(dna_kmer_graph.vertex_labels))
println("DNA k-mers (k=5): $(length(dna_kmers)) unique k-mers")
println("Example k-mers: $(dna_kmers[1:min(3, length(dna_kmers))])")

# Get vertex data to see coverage information
first_kmer = first(dna_kmers)
vertex_data = dna_kmer_graph[first_kmer]
println("Coverage for $(first_kmer): $(length(vertex_data.coverage)) observations")

# Example: RNA sequence analysis
rna_sequence = "AUCGAUCGAUCGAUCGUAGCUAGCUAGCU"
println("\\nRNA sequence: $rna_sequence")

rna_record = FASTX.FASTA.Record("sample_rna", rna_sequence)
rna_records = [rna_record]

# Create RNA k-mer graph
rna_kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.RNAKmer{4}, rna_records)
rna_kmers = collect(values(rna_kmer_graph.vertex_labels))
println("RNA k-mers (k=4): $(length(rna_kmers)) unique k-mers")
println("Example k-mers: $(rna_kmers[1:min(3, length(rna_kmers))])")

# Example: Protein sequence analysis
protein_sequence = "ALAVALINEGLUTAMINELEUCINE"
println("\\nProtein sequence: $protein_sequence")

protein_record = FASTX.FASTA.Record("sample_protein", protein_sequence)
protein_records = [protein_record]

# Create protein k-mer graph
protein_kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.AAKmer{3}, protein_records)
protein_kmers = collect(values(protein_kmer_graph.vertex_labels))
println("Protein k-mers (k=3): $(length(protein_kmers)) unique k-mers")
println("Example k-mers: $(protein_kmers[1:min(3, length(protein_kmers))])")

# Use cases
println("Use cases: Genome assembly, repeat detection, sequence comparison")

# =============================================================================
# Tutorial 3: Qualmer Graphs - Quality-Aware Sequence Analysis
# =============================================================================

println("\n3. QUALMER GRAPHS - Quality-Aware Sequence Analysis")
println("-"^50)

# Example: High-quality sequencing data
high_quality_seq = "ATCGATCGATCGATCGTAGCTAGCT"
high_quality_scores = "HHHHHHHHHHHHHHHHHHHHHHHHH"  # PHRED 39
println("High-quality sequence: $high_quality_seq")

# Create FASTQ record
hq_record = FASTX.FASTQ.Record("high_quality", high_quality_seq, high_quality_scores)
hq_records = [hq_record]

# Create qualmer graph
hq_qualmer_graph = Mycelia.build_qualmer_graph(hq_records, k=5)
hq_kmers = collect(values(hq_qualmer_graph.vertex_labels))
println("High-quality k-mers: $(length(hq_kmers)) unique k-mers")

# Show quality information
first_hq_kmer = first(hq_kmers)
hq_vertex_data = hq_qualmer_graph[first_hq_kmer]
println("$(first_hq_kmer): joint prob = $(round(hq_vertex_data.joint_probability, digits=4)), mean quality = $(hq_vertex_data.mean_quality)")

# Example: Medium-quality sequencing data
medium_quality_seq = "ATCGATCGATCGATCGTAGCTAGCT"
medium_quality_scores = "?????????????????????????"  # PHRED 30
println("\\nMedium-quality sequence: $medium_quality_seq")

mq_record = FASTX.FASTQ.Record("medium_quality", medium_quality_seq, medium_quality_scores)
mq_records = [mq_record]

# Create qualmer graph
mq_qualmer_graph = Mycelia.build_qualmer_graph(mq_records, k=5)
mq_kmers = collect(values(mq_qualmer_graph.vertex_labels))
println("Medium-quality k-mers: $(length(mq_kmers)) unique k-mers")

first_mq_kmer = first(mq_kmers)
mq_vertex_data = mq_qualmer_graph[first_mq_kmer]
println("$(first_mq_kmer): joint prob = $(round(mq_vertex_data.joint_probability, digits=4)), mean quality = $(mq_vertex_data.mean_quality)")

# Quality comparison
println("\\nQuality comparison:")
println("High-quality joint probability: $(round(hq_vertex_data.joint_probability, digits=4))")
println("Medium-quality joint probability: $(round(mq_vertex_data.joint_probability, digits=4))")
println("Quality-aware assembly can make better decisions using this information!")

# Use cases
println("Use cases: Error correction, quality-aware assembly, variant calling")

# =============================================================================
# Tutorial 4: String Graphs - Simplified N-gram Analysis
# =============================================================================

println("\n4. STRING GRAPHS - Simplified N-gram Analysis")
println("-"^50)

# Example: Text processing with simplification
text_for_strings = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
println("Input text: $text_for_strings")

# Create N-gram graph first
string_ngram_graph = Mycelia.string_to_ngram_graph(s=text_for_strings, n=3)
original_ngrams = collect(values(string_ngram_graph.vertex_labels))
println("Original N-grams: $(length(original_ngrams)) 3-grams")

# String graphs are created by simplifying N-gram graphs
# (Path collapsing functionality may be work in progress)
try
    collapsed_paths = Mycelia.collapse_unbranching_paths(string_ngram_graph)
    println("Simplified paths: $(length(collapsed_paths)) paths")
catch e
    println("Path simplification: Work in progress - $(typeof(e))")
end

# Use cases
println("Use cases: Text compression, sequence simplification, path analysis")

# =============================================================================
# Tutorial 5: FASTA Graphs - Variable-Length BioSequence Analysis
# =============================================================================

println("\n5. FASTA GRAPHS - Variable-Length BioSequence Analysis")
println("-"^50)

# Example: Direct BioSequence graph construction
fasta_sequence = "ATCGATCGATCGATCGTAGCTAGCTAGCTGCGCGCGC"
println("FASTA sequence: $fasta_sequence")

fasta_record = FASTX.FASTA.Record("sample", fasta_sequence)
fasta_records = [fasta_record]

# Create BioSequence graph directly
bio_graph = Mycelia.build_biosequence_graph(fasta_records)
bio_sequences = collect(values(bio_graph.vertex_labels))
println("BioSequence graph: $(length(bio_sequences)) sequences")

if !isempty(bio_sequences)
    first_seq = first(bio_sequences)
    println("First sequence: $first_seq (length: $(length(first_seq)))")
    println("Sequence type: $(typeof(first_seq))")
end

# Example: K-mer to BioSequence conversion
println("\\nK-mer to BioSequence conversion:")
kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, fasta_records)
bio_converted = Mycelia.kmer_graph_to_biosequence_graph(kmer_graph)
converted_sequences = collect(values(bio_converted.vertex_labels))
println("Converted $(length(kmer_graph.vertex_labels)) k-mers to $(length(converted_sequences)) BioSequences")

# Use cases
println("Use cases: Genome assembly, contig construction, sequence alignment")

# =============================================================================
# Tutorial 6: FASTQ Graphs - Quality-Aware BioSequence Analysis
# =============================================================================

println("\n6. FASTQ GRAPHS - Quality-Aware BioSequence Analysis")
println("-"^50)

# Example: Quality-aware BioSequence graph
fastq_sequence = "ATCGATCGATCGATCGTAGCTAGCTAGCTGCGCGCGC"
fastq_quality = repeat("H", length(fastq_sequence))  # High quality
println("FASTQ sequence: $fastq_sequence")
println("Quality string: $fastq_quality")

fastq_record = FASTX.FASTQ.Record("quality_sample", fastq_sequence, fastq_quality)
fastq_records = [fastq_record]

# Create quality-aware BioSequence graph
quality_bio_graph = Mycelia.build_quality_biosequence_graph(fastq_records)
quality_sequences = collect(values(quality_bio_graph.vertex_labels))
println("Quality BioSequence graph: $(length(quality_sequences)) sequences")

if !isempty(quality_sequences)
    first_quality_seq = first(quality_sequences)
    quality_vertex_data = quality_bio_graph[first_quality_seq]
    println("First sequence: $first_quality_seq")
    println("Quality scores: $(quality_vertex_data.quality_scores[1:min(10, length(quality_vertex_data.quality_scores))]...)")
end

# Example: Qualmer to quality BioSequence conversion
println("\\nQualmer to quality BioSequence conversion:")
qualmer_graph = Mycelia.build_qualmer_graph(fastq_records, k=5)
quality_converted = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
quality_converted_sequences = collect(values(quality_converted.vertex_labels))
println("Converted $(length(qualmer_graph.vertex_labels)) qualmers to $(length(quality_converted_sequences)) quality BioSequences")

# Use cases
println("Use cases: Quality-aware assembly, error correction, variant calling")

# =============================================================================
# Tutorial 7: Graph Type Hierarchy and Conversions
# =============================================================================

println("\n7. GRAPH TYPE HIERARCHY AND CONVERSIONS")
println("-"^50)

# Demonstrate the hierarchy: Fixed-length → Variable-length
sequence_for_hierarchy = "ATCGATCGATCGATCGTAGCTAGCTAGCTGCGCGCGC"
records_for_hierarchy = [FASTX.FASTA.Record("hierarchy_demo", sequence_for_hierarchy)]
fastq_records_for_hierarchy = [FASTX.FASTQ.Record("hierarchy_demo", sequence_for_hierarchy, repeat("H", length(sequence_for_hierarchy)))]

println("Original sequence: $sequence_for_hierarchy")

# Fixed-length graphs
ngram_result = Mycelia.string_to_ngram_graph(s=sequence_for_hierarchy, n=3)
kmer_result = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, records_for_hierarchy)
qualmer_result = Mycelia.build_qualmer_graph(fastq_records_for_hierarchy, k=5)

println("\\nFixed-length graphs:")
println("  N-grams (n=3): $(length(ngram_result.vertex_labels)) vertices")
println("  K-mers (k=5): $(length(kmer_result.vertex_labels)) vertices")
println("  Qualmers (k=5): $(length(qualmer_result.vertex_labels)) vertices")

# Variable-length graphs
bio_result = Mycelia.build_biosequence_graph(records_for_hierarchy)
quality_bio_result = Mycelia.build_quality_biosequence_graph(fastq_records_for_hierarchy)

println("\\nVariable-length graphs:")
println("  BioSequences: $(length(bio_result.vertex_labels)) vertices")
println("  Quality BioSequences: $(length(quality_bio_result.vertex_labels)) vertices")

# Conversions
kmer_to_bio = Mycelia.kmer_graph_to_biosequence_graph(kmer_result)
qualmer_to_quality_bio = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_result)

println("\\nConversions:")
println("  K-mer → BioSequence: $(length(kmer_result.vertex_labels)) → $(length(kmer_to_bio.vertex_labels))")
println("  Qualmer → Quality BioSequence: $(length(qualmer_result.vertex_labels)) → $(length(qualmer_to_quality_bio.vertex_labels))")

# =============================================================================
# Tutorial 8: Practical Assembly Workflow
# =============================================================================

println("\n8. PRACTICAL ASSEMBLY WORKFLOW")
println("-"^50)

# Simulate a realistic assembly scenario
simulated_reads = [
    "ATCGATCGATCGATCGTAGC",
    "GATCGATCGATCGTAGCTAG",
    "TCGATCGATCGTAGCTAGCT",
    "GATCGATCGTAGCTAGCTGC",
    "TCGATCGTAGCTAGCTGCGC"
]

println("Simulated sequencing reads:")
for (i, read) in enumerate(simulated_reads)
    println("  Read $i: $read")
end

# Create records
fasta_reads = [FASTX.FASTA.Record("read_$i", read) for (i, read) in enumerate(simulated_reads)]
fastq_reads = [FASTX.FASTQ.Record("read_$i", read, repeat("H", length(read))) for (i, read) in enumerate(simulated_reads)]

# Method 1: K-mer graph assembly
println("\\nMethod 1: K-mer graph assembly")
assembly_kmer_graph = Mycelia.build_kmer_graph_next(Mycelia.Kmers.DNAKmer{5}, fasta_reads)
assembly_kmers = collect(values(assembly_kmer_graph.vertex_labels))
println("  K-mer graph: $(length(assembly_kmers)) unique k-mers")

# Method 2: Quality-aware Qualmer assembly
println("\\nMethod 2: Quality-aware Qualmer assembly")
assembly_qualmer_graph = Mycelia.build_qualmer_graph(fastq_reads, k=5)
assembly_qualmers = collect(values(assembly_qualmer_graph.vertex_labels))
println("  Qualmer graph: $(length(assembly_qualmers)) unique qualmers")

# Show quality information
if !isempty(assembly_qualmers)
    first_qualmer = first(assembly_qualmers)
    qualmer_vertex_data = assembly_qualmer_graph[first_qualmer]
    println("  Quality example: $(first_qualmer) with joint prob = $(round(qualmer_vertex_data.joint_probability, digits=4))")
end

# Method 3: BioSequence graph assembly
println("\\nMethod 3: BioSequence graph assembly")
assembly_bio_graph = Mycelia.build_biosequence_graph(fasta_reads)
assembly_sequences = collect(values(assembly_bio_graph.vertex_labels))
println("  BioSequence graph: $(length(assembly_sequences)) sequences")

println("\\nWorkflow summary:")
println("  1. Choose graph type based on data: FASTA → K-mer/BioSequence, FASTQ → Qualmer/Quality BioSequence")
println("  2. Consider quality: Use Qualmer graphs for quality-aware assembly")
println("  3. Use appropriate k-mer size: Smaller k for higher sensitivity, larger k for specificity")
println("  4. Convert between graph types as needed for different analysis steps")

# =============================================================================
# Summary
# =============================================================================

println("\n" * "="^80)
println("TUTORIAL SUMMARY")
println("="^80)
println("You've learned about all 6 graph types in Mycelia:")
println()
println("FIXED-LENGTH GRAPHS (Assembly Foundation):")
println("  1. N-gram Graphs - Unicode text analysis")
println("  2. K-mer Graphs - DNA/RNA/protein sequence analysis")
println("  3. Qualmer Graphs - Quality-aware sequence analysis")
println()
println("VARIABLE-LENGTH GRAPHS (Simplified Products):")
println("  4. String Graphs - Simplified N-gram analysis")
println("  5. FASTA Graphs - Variable-length BioSequence analysis")
println("  6. FASTQ Graphs - Quality-aware BioSequence analysis")
println()
println("KEY FEATURES:")
println("  ✓ Type-stable implementation with storage parameters")
println("  ✓ Quality-aware assembly with per-base quality scores")
println("  ✓ Multi-alphabet support (DNA, RNA, amino acids)")
println("  ✓ Strand-aware graph construction")
println("  ✓ Graph type conversions and hierarchy")
println("  ✓ MetaGraphsNext.jl foundation for performance")
println()
println("Choose the right graph type for your analysis:")
println("  - Use N-gram graphs for text analysis")
println("  - Use K-mer graphs for basic sequence analysis")
println("  - Use Qualmer graphs for quality-aware assembly")
println("  - Use BioSequence graphs for variable-length analysis")
println("  - Use Quality BioSequence graphs for quality-aware variable-length analysis")
println("="^80)