# # Tutorial 6: FASTQ Sequences to Qualmer Graphs and Back
#
# This tutorial demonstrates the round-trip workflow from FASTQ sequences through
# quality-aware k-mer (qualmer) graphs and back to reconstructed sequences.
# Qualmers incorporate both sequence information and quality scores for more
# accurate assembly decisions.

# ## Learning Objectives
# 
# By the end of this tutorial, you will understand:
# 1. How to create qualmer graphs from FASTQ data
# 2. How quality scores influence k-mer confidence
# 3. How to perform quality-aware assembly using package functions
# 4. How to reconstruct sequences while preserving quality information
# 5. The advantages of quality-aware vs coverage-only assembly

# ## Setup and Imports
# Following CLAUDE.md standards: only import top-level packages, use full namespacing

import Mycelia
import Test
import Statistics

# ## Part 1: Understanding Qualmers
#
# Qualmers are k-mers with associated quality information. When the same k-mer
# is observed multiple times with different quality scores, we calculate a
# joint probability that represents our confidence in that k-mer's correctness.

# ### Creating Sample FASTQ Data with Varying Quality

function create_sample_fastq_data()
    ## High-quality read
    hq_seq = "ATCGATCGATCGATCGATCG"
    hq_qual = "IIIIIIIIIIIIIIIIIIII"  ## Phred 40 (99.99% accuracy)
    
    ## Medium-quality read with overlap
    mq_seq = "GATCGATCGATCGATCGTAG"  
    mq_qual = "FFFFFFFFFFFFFFFFFF@@"  ## Phred 37 (99.98%) with lower end
    
    ## Low-quality read with errors
    lq_seq = "GATCGATCGATCGATCGTAC"  ## Error at end (C instead of G)
    lq_qual = "555555555555555555##"  ## Phred 20 (99%) with very low end
    
    ## Create FASTQ records
    records = [
        FASTX.FASTQ.Record("read1", hq_seq, hq_qual),
        FASTX.FASTQ.Record("read2", mq_seq, mq_qual),
        FASTX.FASTQ.Record("read3", lq_seq, lq_qual)
    ]
    
    return records
end

println("Creating sample FASTQ data with varying quality scores...")
fastq_records = create_sample_fastq_data()

# Display the records
for (i, record) in enumerate(fastq_records)
    println("\nRead $i:")
    println("  Sequence: ", String(FASTX.sequence(record)))
    println("  Quality:  ", String(FASTX.quality(record)))
    println("  Phred scores: ", FASTX.quality_scores(record))
end

# ## Part 2: Building Qualmer Graphs
#
# Qualmer graphs combine k-mer information with quality scores to make
# more informed assembly decisions.

# ### Build a qualmer graph with k=7
k = 7
println("\n" * "="^60)
println("Building qualmer graph with k=$k...")

qualmer_graph = Mycelia.build_qualmer_graph(fastq_records; k=k, graph_mode=Mycelia.DoubleStrand)

# Examine graph properties
println("\nQualmer graph statistics:")
println("  Number of vertices (unique k-mers): ", Graphs.nv(qualmer_graph))
println("  Number of edges: ", Graphs.ne(qualmer_graph))

# ### Inspect qualmer properties
println("\n" * "="^60)
println("Examining qualmer vertices and their properties:")

# Get first few vertices to examine
for v in Iterators.take(Graphs.vertices(qualmer_graph), 5)
    vertex_data = qualmer_graph[v]
    println("\nVertex $v:")
    println("  K-mer: ", vertex_data.kmer)
    println("  Coverage: ", vertex_data.coverage)
    println("  Mean quality: ", round(vertex_data.mean_quality, digits=2))
    println("  Joint probability: ", round(vertex_data.joint_probability, digits=4))
end

# ## Part 3: Quality-Aware vs Coverage-Only Assembly
#
# Let's compare how quality information affects assembly decisions.

# ### Find high-confidence paths using quality information
println("\n" * "="^60)
println("Finding high-confidence paths through the qualmer graph...")

# Get all vertices sorted by joint probability (confidence)
vertices_by_confidence = sort(collect(Graphs.vertices(qualmer_graph)), 
    by=v -> qualmer_graph[v].joint_probability, rev=true)

println("\nTop 5 most confident k-mers:")
for v in vertices_by_confidence[1:min(5, length(vertices_by_confidence))]
    vdata = qualmer_graph[v]
    println("  ", vdata.kmer, 
            " - Coverage: ", vdata.coverage,
            ", Joint prob: ", round(vdata.joint_probability, digits=4),
            ", Mean quality: ", round(vdata.mean_quality, digits=1))
end

# ### Compare with coverage-only approach
vertices_by_coverage = sort(collect(Graphs.vertices(qualmer_graph)), 
    by=v -> qualmer_graph[v].coverage, rev=true)

println("\nTop 5 highest coverage k-mers:")
for v in vertices_by_coverage[1:min(5, length(vertices_by_coverage))]
    vdata = qualmer_graph[v]
    println("  ", vdata.kmer, 
            " - Coverage: ", vdata.coverage,
            ", Joint prob: ", round(vdata.joint_probability, digits=4),
            ", Mean quality: ", round(vdata.mean_quality, digits=1))
end

# ## Part 4: Quality-Aware Path Finding
#
# Use the package function to find the most likely path through the graph.

# Find best starting vertex (highest confidence)
start_vertex = vertices_by_confidence[1]
quality_path = Mycelia.find_quality_weighted_path(qualmer_graph, start_vertex)

println("\n" * "="^60)
println("Quality-weighted path through graph:")
println("Path length: ", length(quality_path), " vertices")

# Show path k-mers and qualities
println("\nPath details:")
for (i, v) in enumerate(quality_path[1:min(10, length(quality_path))])
    vdata = qualmer_graph[v]
    println("  Step $i: ", vdata.kmer, 
            " (prob: ", round(vdata.joint_probability, digits=3), ")")
end

# ## Part 5: Converting to Quality-Aware BioSequence Graph
#
# Convert the qualmer graph to a variable-length quality-aware sequence graph.

println("\n" * "="^60)
println("Converting qualmer graph to quality-aware BioSequence graph...")

# Convert to FASTQ graph (variable-length with quality)
fastq_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph, k)

println("\nFASTQ graph statistics:")
println("  Number of vertices: ", Graphs.nv(fastq_graph))
println("  Number of edges: ", Graphs.ne(fastq_graph))

# Examine simplified vertices
println("\nExamining quality-aware sequence vertices:")
for v in Iterators.take(Graphs.vertices(fastq_graph), 3)
    vertex_data = fastq_graph[v]
    println("\nVertex $v:")
    println("  Sequence: ", vertex_data.sequence)
    println("  Length: ", length(vertex_data.sequence))
    println("  Quality scores: ", vertex_data.quality_scores[1:min(20, length(vertex_data.quality_scores))], "...")
    println("  Mean quality: ", round(Statistics.mean(vertex_data.quality_scores), digits=1))
end

# ## Part 6: Round-Trip Reconstruction
#
# Reconstruct FASTQ records from the quality-aware graph.

println("\n" * "="^60)
println("Reconstructing FASTQ records from the graph...")

# Extract paths and convert to FASTQ records
reconstructed_records = Mycelia.quality_biosequence_graph_to_fastq(fastq_graph, "reconstructed")

println("\nReconstructed ", length(reconstructed_records), " FASTQ records")

# Compare with original
println("\nComparison with original reads:")
for (i, (orig, recon)) in enumerate(zip(fastq_records[1:min(3, length(reconstructed_records))], 
                                        reconstructed_records[1:min(3, length(reconstructed_records))]))
    orig_seq = String(FASTX.sequence(orig))
    recon_seq = String(FASTX.sequence(recon))
    
    println("\nRead $i:")
    println("  Original:      ", orig_seq)
    println("  Reconstructed: ", recon_seq)
    println("  Match: ", orig_seq == recon_seq ? "✓" : "✗")
    
    ## Compare quality scores
    orig_qual = FASTX.quality_scores(orig)
    recon_qual = FASTX.quality_scores(recon)
    println("  Original quality range: ", minimum(orig_qual), "-", maximum(orig_qual))
    println("  Reconstructed quality range: ", minimum(recon_qual), "-", maximum(recon_qual))
end

# ## Part 7: Error Correction Using Quality Information
#
# Demonstrate how quality scores help identify and correct errors using package functions.

println("\n" * "="^60)
println("Demonstrating quality-aware error correction...")

# Create reads with a known error
error_reads = [
    FASTX.FASTQ.Record("good1", "ATCGATCGATCG", "IIIIIIIIIIII"),  ## High quality
    FASTX.FASTQ.Record("good2", "TCGATCGATCGA", "IIIIIIIIIIII"),  ## High quality
    FASTX.FASTQ.Record("error", "TCGATCTATCGA", "IIIIII##IIII"),  ## Error at low quality position
]

# Build qualmer graph
error_graph = Mycelia.build_qualmer_graph(error_reads; k=5, graph_mode=Mycelia.SingleStrand)

println("\nAnalyzing k-mers around error position:")
# The error creates k-mers: GATCT (wrong) vs GATCG (correct)
for v in Graphs.vertices(error_graph)
    vdata = error_graph[v]
    kmer_str = String(vdata.kmer)
    if occursin("GATC", kmer_str)
        println("  K-mer: ", kmer_str,
                " - Coverage: ", vdata.coverage,
                ", Joint prob: ", round(vdata.joint_probability, digits=4))
    end
end

# Use package function to identify potential errors
potential_errors = Mycelia.identify_potential_errors(error_graph)
println("\nPotential error k-mers identified: ", length(potential_errors))

for error_v in potential_errors
    vdata = error_graph[error_v]
    println("  Error k-mer: ", vdata.kmer,
            " - Coverage: ", vdata.coverage,
            ", Quality: ", round(vdata.mean_quality, digits=1),
            ", Confidence: ", round(vdata.joint_probability, digits=4))
end

# ## Part 8: Advanced Quality Metrics
#
# Calculate assembly quality metrics using the package function.

metrics = Mycelia.calculate_assembly_quality_metrics(qualmer_graph)

println("\n" * "="^60)
println("Assembly quality metrics:")
println("  Mean k-mer coverage: ", round(metrics.mean_coverage, digits=2))
println("  Mean quality score: ", round(metrics.mean_quality, digits=1))
println("  Mean k-mer confidence: ", round(metrics.mean_confidence, digits=4))
println("  Fraction of low-confidence k-mers: ", round(metrics.low_confidence_fraction, digits=3))
println("  Total unique k-mers: ", metrics.total_kmers)

# ## Part 9: Practical Example - Assembling a Small Genome Region
#
# Let's create a more realistic example with overlapping reads from a genome region.

function create_genome_region_reads()
    ## Simulate a 50bp genome region
    true_sequence = "ATCGATCGATCGTAGCTAGCTAGCTTGCATGCATGCATGCATGCATGCAT"
    
    ## Generate overlapping reads with varying quality
    reads = []
    
    ## High quality reads
    push!(reads, FASTX.FASTQ.Record("hq1", true_sequence[1:25], "I"^25))
    push!(reads, FASTX.FASTQ.Record("hq2", true_sequence[15:40], "I"^26))
    push!(reads, FASTX.FASTQ.Record("hq3", true_sequence[30:50], "I"^21))
    
    ## Medium quality reads with some errors
    read_mq1 = true_sequence[5:30]
    qual_mq1 = "FFFFFFFFFFFFFFFFFFFFFFFFFF"
    push!(reads, FASTX.FASTQ.Record("mq1", read_mq1, qual_mq1))
    
    ## Low quality read with error
    read_lq1 = true_sequence[20:45]
    read_lq1 = read_lq1[1:10] * "T" * read_lq1[12:end]  ## Error at position 11
    qual_lq1 = "AAAAAAAAAA#AAAAAAAAAAAAAA"  ## Low quality at error
    push!(reads, FASTX.FASTQ.Record("lq1", read_lq1, qual_lq1))
    
    return reads, true_sequence
end

genome_reads, true_seq = create_genome_region_reads()

println("\n" * "="^60)
println("Assembling genome region from overlapping reads...")
println("True sequence: ", true_seq)
println("Number of reads: ", length(genome_reads))

# Build qualmer graph
genome_graph = Mycelia.build_qualmer_graph(genome_reads; k=9, graph_mode=Mycelia.SingleStrand)

# Convert to sequence graph and extract contigs
seq_graph = Mycelia.qualmer_graph_to_quality_biosequence_graph(genome_graph, 9)

println("\nAssembly results:")
println("  Qualmer graph: ", Graphs.nv(genome_graph), " vertices, ", Graphs.ne(genome_graph), " edges")
println("  Sequence graph: ", Graphs.nv(seq_graph), " vertices, ", Graphs.ne(seq_graph), " edges")

# Use package function to find quality-weighted path
if Graphs.nv(genome_graph) > 0
    ## Start from highest confidence vertex
    confidence_sorted = sort(collect(Graphs.vertices(genome_graph)), 
        by=v -> genome_graph[v].joint_probability, rev=true)
    best_path = Mycelia.find_quality_weighted_path(genome_graph, confidence_sorted[1])
    
    println("\nBest quality-weighted path:")
    println("  Path length: ", length(best_path), " k-mers")
    
    ## Reconstruct sequence from path
    if length(best_path) > 1
        path_kmers = [String(genome_graph[v].kmer) for v in best_path]
        ## Simple reconstruction: first k-mer + last base of each subsequent k-mer
        reconstructed = path_kmers[1] * join([kmer[end] for kmer in path_kmers[2:end]])
        
        println("  Reconstructed length: ", length(reconstructed))
        println("  Reconstructed: ", reconstructed)
        
        ## Check accuracy
        if reconstructed == true_seq
            println("  ✓ Perfect reconstruction!")
        elseif occursin(reconstructed, true_seq)
            println("  ✓ Assembled sequence is a substring of true sequence")
        elseif occursin(true_seq, reconstructed)
            println("  ✓ True sequence is a substring of assembled sequence")
        else
            println("  ✗ Assembly differs from true sequence")
            println("  True:      ", true_seq)
            println("  Assembled: ", reconstructed)
        end
    end
end

# Find longest path (contig) from sequence graph
if Graphs.nv(seq_graph) > 0
    ## Get vertex with longest sequence
    longest_v = argmax(v -> length(seq_graph[v].sequence), Graphs.vertices(seq_graph))
    longest_seq = seq_graph[longest_v].sequence
    longest_qual = seq_graph[longest_v].quality_scores
    
    println("\nLongest contig from simplified graph:")
    println("  Length: ", length(longest_seq))
    println("  Sequence: ", longest_seq)
    println("  Mean quality: ", round(Statistics.mean(longest_qual), digits=1))
    
    ## Check accuracy
    if String(longest_seq) == true_seq
        println("  ✓ Perfect reconstruction!")
    else
        ## Find best alignment
        true_str = true_seq
        assembled_str = String(longest_seq)
        if occursin(assembled_str, true_str)
            println("  ✓ Assembled sequence is a substring of true sequence")
        elseif occursin(true_str, assembled_str)
            println("  ✓ True sequence is a substring of assembled sequence")
        else
            println("  ✗ Assembly differs from true sequence")
        end
    end
end

# Calculate final quality metrics for the genome assembly
final_metrics = Mycelia.calculate_assembly_quality_metrics(genome_graph)
println("\nFinal assembly quality metrics:")
println("  Mean coverage: ", round(final_metrics.mean_coverage, digits=2))
println("  Mean quality: ", round(final_metrics.mean_quality, digits=1))
println("  Mean confidence: ", round(final_metrics.mean_confidence, digits=4))
println("  Low confidence fraction: ", round(final_metrics.low_confidence_fraction, digits=3))

# ## Summary
#
# In this tutorial, we've demonstrated:
#
# 1. **Qualmer Construction**: Building quality-aware k-mer graphs from FASTQ data
# 2. **Joint Probability**: How multiple observations with different qualities are combined
# 3. **Quality vs Coverage**: The advantage of using quality scores over coverage alone
# 4. **Package Functions**: Using Mycelia's built-in functions for quality-weighted analysis:
#    - `find_quality_weighted_path()` for optimal path finding
#    - `calculate_assembly_quality_metrics()` for comprehensive quality assessment
#    - `identify_potential_errors()` for error detection
# 5. **Round-Trip Conversion**: Maintaining quality information through graph transformations
# 6. **Practical Assembly**: Using quality information for more accurate genome assembly
#
# Key advantages of qualmer graphs:
# - Better discrimination between true k-mers and errors
# - Quality-weighted path finding for more accurate assembly
# - Preservation of quality information for downstream analysis
# - Improved handling of repetitive regions with varying quality
# - Built-in error detection and quality assessment capabilities

println("\n" * "="^60)
println("Tutorial 6 completed!")
println("You've learned how to use quality-aware k-mer graphs for improved assembly accuracy.")
println("All analysis was performed using Mycelia's built-in qualmer analysis functions.")