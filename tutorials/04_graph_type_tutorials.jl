"""
Rhizomorph Graph Type Tutorials

This tutorial demonstrates the usage of all 6 graph types in Mycelia's Rhizomorph submodule:
1. N-gram Graphs - For unicode text analysis
2. K-mer Graphs - For DNA/RNA/protein sequence analysis
3. Qualmer Graphs - For quality-aware sequence analysis
4. String Graphs - For variable-length string overlaps (OLC)
5. FASTA Graphs - For variable-length BioSequence overlaps (OLC)
6. FASTQ Graphs - For quality-aware BioSequence overlaps (OLC)

Each section shows practical examples with real data using Rhizomorph APIs.
"""

import Mycelia
import FASTX
import MetaGraphsNext
import Statistics

println("="^80)
println("RHIZOMORPH GRAPH TYPE TUTORIALS")
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
ngram_graph = Mycelia.Rhizomorph.build_ngram_graph([text], 3; dataset_id="text_demo")
ngram_labels = collect(MetaGraphsNext.labels(ngram_graph))
println("Number of 3-grams: $(length(ngram_labels))")

# Show some N-grams
ngrams = ngram_labels
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
dna_dataset_id = "dna_demo"
dna_kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, 5; dataset_id=dna_dataset_id)
dna_kmers = collect(MetaGraphsNext.labels(dna_kmer_graph))
println("DNA k-mers (k=5): $(length(dna_kmers)) unique k-mers")
println("Example k-mers: $(dna_kmers[1:min(3, length(dna_kmers))])")

# Get vertex data to see evidence counts
first_kmer = first(dna_kmers)
vertex_data = dna_kmer_graph[first_kmer]
coverage = Mycelia.Rhizomorph.count_evidence_entries(vertex_data)
println("Evidence entries for $(first_kmer): $(coverage) observations")

# Example: RNA sequence analysis
rna_sequence = "AUCGAUCGAUCGAUCGUAGCUAGCUAGCU"
println("\\nRNA sequence: $rna_sequence")

rna_record = FASTX.FASTA.Record("sample_rna", rna_sequence)
rna_records = [rna_record]

# Create RNA k-mer graph
rna_kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(rna_records, 4; dataset_id="rna_demo")
rna_kmers = collect(MetaGraphsNext.labels(rna_kmer_graph))
println("RNA k-mers (k=4): $(length(rna_kmers)) unique k-mers")
println("Example k-mers: $(rna_kmers[1:min(3, length(rna_kmers))])")

# Example: Protein sequence analysis
protein_sequence = "ALAVALINEGLUTAMINELEUCINE"
println("\\nProtein sequence: $protein_sequence")

protein_record = FASTX.FASTA.Record("sample_protein", protein_sequence)
protein_records = [protein_record]

# Create protein k-mer graph
protein_kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(protein_records, 3; dataset_id="protein_demo")
protein_kmers = collect(MetaGraphsNext.labels(protein_kmer_graph))
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
high_quality_scores = "HHHHHHHHHHHHHHHHHHHHHHHHH"  ## PHRED 39
println("High-quality sequence: $high_quality_seq")

# Create FASTQ record
hq_record = FASTX.FASTQ.Record("high_quality", high_quality_seq, high_quality_scores)
hq_records = [hq_record]

# Create qualmer graph
hq_dataset_id = "high_quality"
hq_qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(hq_records, 5; dataset_id=hq_dataset_id)
hq_kmers = collect(MetaGraphsNext.labels(hq_qualmer_graph))
println("High-quality k-mers: $(length(hq_kmers)) unique k-mers")

# Show quality information
first_hq_kmer = first(hq_kmers)
hq_vertex_data = hq_qualmer_graph[first_hq_kmer]
hq_joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(hq_vertex_data, hq_dataset_id)
hq_mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(hq_vertex_data, hq_dataset_id)
hq_joint_mean = Statistics.mean(hq_joint_quality)
hq_mean_mean = Statistics.mean(hq_mean_quality)
println("$(first_hq_kmer): mean joint Q = $(round(hq_joint_mean, digits=2)), mean Q = $(round(hq_mean_mean, digits=2))")

# Example: Medium-quality sequencing data
medium_quality_seq = "ATCGATCGATCGATCGTAGCTAGCT"
medium_quality_scores = "?????????????????????????"  ## PHRED 30
println("\\nMedium-quality sequence: $medium_quality_seq")

mq_record = FASTX.FASTQ.Record("medium_quality", medium_quality_seq, medium_quality_scores)
mq_records = [mq_record]

# Create qualmer graph
mq_dataset_id = "medium_quality"
mq_qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(mq_records, 5; dataset_id=mq_dataset_id)
mq_kmers = collect(MetaGraphsNext.labels(mq_qualmer_graph))
println("Medium-quality k-mers: $(length(mq_kmers)) unique k-mers")

first_mq_kmer = first(mq_kmers)
mq_vertex_data = mq_qualmer_graph[first_mq_kmer]
mq_joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(mq_vertex_data, mq_dataset_id)
mq_mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(mq_vertex_data, mq_dataset_id)
mq_joint_mean = Statistics.mean(mq_joint_quality)
mq_mean_mean = Statistics.mean(mq_mean_quality)
println("$(first_mq_kmer): mean joint Q = $(round(mq_joint_mean, digits=2)), mean Q = $(round(mq_mean_mean, digits=2))")

# Quality comparison
println("\\nQuality comparison:")
println("High-quality mean joint Q: $(round(hq_joint_mean, digits=2))")
println("Medium-quality mean joint Q: $(round(mq_joint_mean, digits=2))")
println("Quality-aware assembly can make better decisions using this information!")

# Use cases
println("Use cases: Error correction, quality-aware assembly, variant calling")

# =============================================================================
# Tutorial 4: String Graphs - Variable-Length Overlap Analysis
# =============================================================================

println("\n4. STRING GRAPHS - Variable-Length Overlap Analysis")
println("-"^50)

# Example: Overlap-aware text fragments
string_fragments = ["ABCDEF", "CDEFGH", "EFGHIJ", "GHIJKL"]
println("Input fragments: $(string_fragments)")

# Create variable-length string graph (OLC)
string_graph = Mycelia.Rhizomorph.build_string_graph(string_fragments; dataset_id="string_demo", min_overlap=3)
string_stats = Mycelia.Rhizomorph.get_string_graph_statistics(string_graph)
println("String graph: $(string_stats[:num_vertices]) strings, $(string_stats[:num_edges]) overlaps")

# Show example overlaps if present
if string_stats[:num_edges] > 0
    example_edge = first(MetaGraphsNext.edge_labels(string_graph))
    src, dst = example_edge
    overlap_len = string_graph[src, dst].overlap_length
    println("Example overlap: \"$src\" → \"$dst\" (overlap length: $overlap_len)")
end

# Use cases
println("Use cases: Text reconstruction, fragment assembly, overlap analysis")

# =============================================================================
# Tutorial 5: FASTA Graphs - Variable-Length BioSequence Analysis
# =============================================================================

println("\n5. FASTA GRAPHS - Variable-Length BioSequence Analysis")
println("-"^50)

# Example: Direct FASTA graph construction
fasta_sequence = "ATCGATCGATCGATCGTAGCTAGCTAGCTGCGCGCGC"
println("FASTA sequence: $fasta_sequence")

fasta_record = FASTX.FASTA.Record("sample", fasta_sequence)
fasta_records = [fasta_record]

# Create FASTA overlap graph directly (OLC)
bio_graph = Mycelia.Rhizomorph.build_fasta_graph(fasta_records; dataset_id="fasta_demo", min_overlap=5)
bio_sequences = collect(MetaGraphsNext.labels(bio_graph))
println("FASTA graph: $(length(bio_sequences)) sequences")

if !isempty(bio_sequences)
    first_seq = first(bio_sequences)
    println("First sequence: $first_seq (length: $(length(first_seq)))")
    println("Sequence type: $(typeof(first_seq))")
end

# Example: K-mer to FASTA graph conversion
println("\\nK-mer to FASTA graph conversion:")
kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(fasta_records, 5; dataset_id="kmer_demo")
bio_converted = Mycelia.Rhizomorph.convert_fixed_to_variable(kmer_graph)
converted_sequences = collect(MetaGraphsNext.labels(bio_converted))
kmer_labels = collect(MetaGraphsNext.labels(kmer_graph))
println("Converted $(length(kmer_labels)) k-mers to $(length(converted_sequences)) FASTA sequences")

# Use cases
println("Use cases: Genome assembly, contig construction, sequence alignment")

# =============================================================================
# Tutorial 6: FASTQ Graphs - Quality-Aware BioSequence Analysis
# =============================================================================

println("\n6. FASTQ GRAPHS - Quality-Aware BioSequence Analysis")
println("-"^50)

# Example: Quality-aware FASTQ graph
fastq_sequence = "ATCGATCGATCGATCGTAGCTAGCTAGCTGCGCGCGC"
fastq_quality = repeat("H", length(fastq_sequence))  ## High quality
println("FASTQ sequence: $fastq_sequence")
println("Quality string: $fastq_quality")

fastq_record = FASTX.FASTQ.Record("quality_sample", fastq_sequence, fastq_quality)
fastq_records = [fastq_record]

# Create quality-aware FASTQ graph
fastq_dataset_id = "fastq_demo"
quality_bio_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id=fastq_dataset_id, min_overlap=5)
quality_sequences = collect(MetaGraphsNext.labels(quality_bio_graph))
println("FASTQ graph: $(length(quality_sequences)) sequences")

if !isempty(quality_sequences)
    first_quality_seq = first(quality_sequences)
    quality_vertex_data = quality_bio_graph[first_quality_seq]
    println("First sequence: $first_quality_seq")
    joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(quality_vertex_data, fastq_dataset_id)
    println("Joint quality (first 10): $(joint_quality[1:min(10, length(joint_quality))])")
end

# Example: Qualmer to FASTQ graph conversion
println("\\nQualmer to FASTQ graph conversion:")
qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, 5; dataset_id=fastq_dataset_id)
quality_converted = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
quality_converted_sequences = collect(MetaGraphsNext.labels(quality_converted))
qualmer_labels = collect(MetaGraphsNext.labels(qualmer_graph))
println("Converted $(length(qualmer_labels)) qualmers to $(length(quality_converted_sequences)) FASTQ sequences")

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
ngram_result = Mycelia.Rhizomorph.build_ngram_graph([sequence_for_hierarchy], 3; dataset_id="hierarchy")
kmer_result = Mycelia.Rhizomorph.build_kmer_graph(records_for_hierarchy, 5; dataset_id="hierarchy")
qualmer_result = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records_for_hierarchy, 5; dataset_id="hierarchy")
ngram_labels = collect(MetaGraphsNext.labels(ngram_result))
kmer_labels = collect(MetaGraphsNext.labels(kmer_result))
qualmer_labels = collect(MetaGraphsNext.labels(qualmer_result))

println("\\nFixed-length graphs:")
println("  N-grams (n=3): $(length(ngram_labels)) vertices")
println("  K-mers (k=5): $(length(kmer_labels)) vertices")
println("  Qualmers (k=5): $(length(qualmer_labels)) vertices")

# Variable-length graphs
bio_result = Mycelia.Rhizomorph.build_fasta_graph(records_for_hierarchy; dataset_id="hierarchy", min_overlap=5)
quality_bio_result = Mycelia.Rhizomorph.build_fastq_graph(fastq_records_for_hierarchy; dataset_id="hierarchy", min_overlap=5)
bio_labels = collect(MetaGraphsNext.labels(bio_result))
quality_bio_labels = collect(MetaGraphsNext.labels(quality_bio_result))

println("\\nVariable-length graphs:")
println("  FASTA graph: $(length(bio_labels)) vertices")
println("  FASTQ graph: $(length(quality_bio_labels)) vertices")

# Conversions
kmer_to_bio = Mycelia.Rhizomorph.convert_fixed_to_variable(kmer_result)
qualmer_to_quality_bio = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_result)
kmer_to_bio_labels = collect(MetaGraphsNext.labels(kmer_to_bio))
qualmer_to_quality_labels = collect(MetaGraphsNext.labels(qualmer_to_quality_bio))

println("\\nConversions:")
println("  K-mer → FASTA graph: $(length(kmer_labels)) → $(length(kmer_to_bio_labels))")
println("  Qualmer → FASTQ graph: $(length(qualmer_labels)) → $(length(qualmer_to_quality_labels))")

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
assembly_kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(fasta_reads, 5; dataset_id="assembly_kmer")
assembly_kmers = collect(MetaGraphsNext.labels(assembly_kmer_graph))
println("  K-mer graph: $(length(assembly_kmers)) unique k-mers")

# Method 2: Quality-aware Qualmer assembly
println("\\nMethod 2: Quality-aware Qualmer assembly")
assembly_qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_reads, 5; dataset_id="assembly_qualmer")
assembly_qualmers = collect(MetaGraphsNext.labels(assembly_qualmer_graph))
println("  Qualmer graph: $(length(assembly_qualmers)) unique qualmers")

# Show quality information
if !isempty(assembly_qualmers)
    first_qualmer = first(assembly_qualmers)
    qualmer_vertex_data = assembly_qualmer_graph[first_qualmer]
    assembly_joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(qualmer_vertex_data, "assembly_qualmer")
    joint_mean = Statistics.mean(assembly_joint_quality)
    println("  Quality example: $(first_qualmer) with mean joint Q = $(round(joint_mean, digits=2))")
end

# Method 3: FASTA graph assembly
println("\\nMethod 3: FASTA graph assembly")
assembly_bio_graph = Mycelia.Rhizomorph.build_fasta_graph(fasta_reads; dataset_id="assembly_fasta", min_overlap=5)
assembly_sequences = collect(MetaGraphsNext.labels(assembly_bio_graph))
println("  FASTA graph: $(length(assembly_sequences)) sequences")

println("\\nWorkflow summary:")
println("  1. Choose graph type based on data: FASTA → K-mer/FASTA graph, FASTQ → Qualmer/FASTQ graph")
println("  2. Consider quality: Use Qualmer graphs for quality-aware assembly")
println("  3. Use appropriate k-mer size: Smaller k for higher sensitivity, larger k for specificity")
println("  4. Convert between graph types as needed for different analysis steps")

# =============================================================================
# Summary
# =============================================================================

println("\n" * "="^80)
println("TUTORIAL SUMMARY")
println("="^80)
println("You've learned about all 6 graph types in Mycelia Rhizomorph:")
println()
println("FIXED-LENGTH GRAPHS (Assembly Foundation):")
println("  1. N-gram Graphs - Unicode text analysis")
println("  2. K-mer Graphs - DNA/RNA/protein sequence analysis")
println("  3. Qualmer Graphs - Quality-aware sequence analysis")
println()
println("VARIABLE-LENGTH GRAPHS (OLC Products):")
println("  4. String Graphs - Variable-length string overlaps")
println("  5. FASTA Graphs - Variable-length BioSequence overlaps")
println("  6. FASTQ Graphs - Quality-aware BioSequence overlaps")
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
println("  - Use FASTA graphs for variable-length analysis")
println("  - Use FASTQ graphs for quality-aware variable-length analysis")
println("="^80)
