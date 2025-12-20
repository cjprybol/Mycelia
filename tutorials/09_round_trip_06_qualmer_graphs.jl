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
import FASTX
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

qualmer_dataset_id = "qualmer_demo"
qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(
    fastq_records,
    k;
    dataset_id=qualmer_dataset_id,
    mode=:doublestrand
)
qualmer_labels = collect(Mycelia.MetaGraphsNext.labels(qualmer_graph))

# Examine graph properties
println("\nQualmer graph statistics:")
println("  Number of vertices (unique k-mers): ", length(qualmer_labels))
println("  Number of edges: ", Mycelia.MetaGraphsNext.ne(qualmer_graph))

# ### Inspect qualmer properties
println("\n" * "="^60)
println("Examining qualmer vertices and their properties:")

# Get first few vertices to examine
for label in Iterators.take(qualmer_labels, 5)
    vertex_data = qualmer_graph[label]
    joint_quality = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, qualmer_dataset_id)
    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, qualmer_dataset_id)
    joint_mean = isnothing(joint_quality) ? 0.0 : Statistics.mean(Float64.(joint_quality))
    mean_mean = isnothing(mean_quality) ? 0.0 : Statistics.mean(mean_quality)
    coverage = Mycelia.Rhizomorph.count_evidence_entries(vertex_data)
    println("\nVertex $label:")
    println("  K-mer: ", label)
    println("  Coverage: ", coverage)
    println("  Mean quality: ", round(mean_mean, digits=2))
    println("  Joint quality: ", round(joint_mean, digits=2))
end

# ## Part 3: Quality-Aware vs Coverage-Only Assembly
#
# Let's compare how quality information affects assembly decisions.

# ### Find high-confidence paths using quality information
println("\n" * "="^60)
println("Finding high-confidence paths through the qualmer graph...")

# Get all vertices sorted by joint probability (confidence)
vertices_by_confidence = sort(
    qualmer_labels,
    by=v -> Mycelia.Rhizomorph.mean_joint_quality(qualmer_graph, v, qualmer_dataset_id),
    rev=true
)

println("\nTop 5 most confident k-mers:")
for v in vertices_by_confidence[1:min(5, length(vertices_by_confidence))]
    vdata = qualmer_graph[v]
    coverage = Mycelia.Rhizomorph.count_evidence_entries(vdata)
    joint_mean = Mycelia.Rhizomorph.mean_joint_quality(qualmer_graph, v, qualmer_dataset_id)
    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vdata, qualmer_dataset_id)
    mean_mean = isnothing(mean_quality) ? 0.0 : Statistics.mean(mean_quality)
    println("  ", v,
            " - Coverage: ", coverage,
            ", Joint Q: ", round(joint_mean, digits=2),
            ", Mean Q: ", round(mean_mean, digits=1))
end

# ### Compare with coverage-only approach
vertices_by_coverage = sort(
    qualmer_labels,
    by=v -> Mycelia.Rhizomorph.count_evidence_entries(qualmer_graph[v]),
    rev=true
)

println("\nTop 5 highest coverage k-mers:")
for v in vertices_by_coverage[1:min(5, length(vertices_by_coverage))]
    vdata = qualmer_graph[v]
    coverage = Mycelia.Rhizomorph.count_evidence_entries(vdata)
    joint_mean = Mycelia.Rhizomorph.mean_joint_quality(qualmer_graph, v, qualmer_dataset_id)
    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vdata, qualmer_dataset_id)
    mean_mean = isnothing(mean_quality) ? 0.0 : Statistics.mean(mean_quality)
    println("  ", v,
            " - Coverage: ", coverage,
            ", Joint Q: ", round(joint_mean, digits=2),
            ", Mean Q: ", round(mean_mean, digits=1))
end

# ## Part 4: Quality-Aware Path Finding
#
# Use Rhizomorph path finding and quality scoring to select a likely path.

# Find candidate paths and select the highest-quality path
paths = Mycelia.Rhizomorph.find_eulerian_paths_next(qualmer_graph)
quality_path = []
best_path_score = 0.0
for path in paths
    score = Mycelia.Rhizomorph.mean_path_quality(qualmer_graph, path, qualmer_dataset_id)
    if score > best_path_score
        best_path_score = score
        quality_path = path
    end
end

println("\n" * "="^60)
println("Quality-weighted path through graph:")
println("Path length: ", length(quality_path), " vertices")

# Show path k-mers and qualities
println("\nPath details:")
for (i, v) in enumerate(quality_path[1:min(10, length(quality_path))])
    joint_mean = Mycelia.Rhizomorph.mean_joint_quality(qualmer_graph, v, qualmer_dataset_id)
    println("  Step $i: ", v,
            " (mean joint Q: ", round(joint_mean, digits=2), ")")
end

# ## Part 5: Converting to Quality-Aware BioSequence Graph
#
# Convert the qualmer graph to a variable-length quality-aware sequence graph.

println("\n" * "="^60)
println("Converting qualmer graph to quality-aware BioSequence graph...")

# Convert to FASTQ graph (variable-length with quality)
fastq_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(qualmer_graph)
fastq_labels = collect(Mycelia.MetaGraphsNext.labels(fastq_graph))

println("\nFASTQ graph statistics:")
println("  Number of vertices: ", length(fastq_labels))
println("  Number of edges: ", Mycelia.MetaGraphsNext.ne(fastq_graph))

# Examine simplified vertices
println("\nExamining quality-aware sequence vertices:")
for label in Iterators.take(fastq_labels, 3)
    vertex_data = fastq_graph[label]
    mean_quality = Statistics.mean(Mycelia.Rhizomorph.decode_quality_scores(vertex_data.quality_scores))
    println("\nVertex: ", label)
    println("  Sequence: ", string(vertex_data.sequence))
    println("  Length: ", length(vertex_data.sequence))
    println("  Quality scores: ", vertex_data.quality_scores[1:min(20, length(vertex_data.quality_scores))], "...")
    println("  Mean quality: ", round(mean_quality, digits=1))
end

# ## Part 6: Round-Trip Reconstruction
#
# Reconstruct FASTQ records from the quality-aware graph.

println("\n" * "="^60)
println("Reconstructing FASTQ records from the graph...")

# Extract paths and convert to FASTQ records
reconstructed_records = Mycelia.Rhizomorph.fastq_graph_to_records(fastq_graph, "reconstructed")

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
error_dataset_id = "error_demo"
error_graph = Mycelia.Rhizomorph.build_qualmer_graph(error_reads, 5; dataset_id=error_dataset_id, mode=:singlestrand)
error_labels = collect(Mycelia.MetaGraphsNext.labels(error_graph))

println("\nAnalyzing k-mers around error position:")
# The error creates k-mers: GATCT (wrong) vs GATCG (correct)
for label in error_labels
    vdata = error_graph[label]
    kmer_str = string(label)
    if occursin("GATC", kmer_str)
        coverage = Mycelia.Rhizomorph.count_evidence_entries(vdata)
        joint_mean = Mycelia.Rhizomorph.mean_joint_quality(error_graph, label, error_dataset_id)
        println("  K-mer: ", kmer_str,
                " - Coverage: ", coverage,
                ", Joint Q: ", round(joint_mean, digits=2))
    end
end

# Use quality filtering to identify potential errors
high_quality = Mycelia.Rhizomorph.find_high_quality_kmers(error_graph, 30; dataset_id=error_dataset_id)
potential_errors = [label for label in error_labels if !(label in high_quality)]
println("\nPotential error k-mers identified: ", length(potential_errors))

for error_label in potential_errors
    vdata = error_graph[error_label]
    coverage = Mycelia.Rhizomorph.count_evidence_entries(vdata)
    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vdata, error_dataset_id)
    mean_mean = isnothing(mean_quality) ? 0.0 : Statistics.mean(mean_quality)
    joint_mean = Mycelia.Rhizomorph.mean_joint_quality(error_graph, error_label, error_dataset_id)
    println("  Error k-mer: ", error_label,
            " - Coverage: ", coverage,
            ", Mean Q: ", round(mean_mean, digits=1),
            ", Joint Q: ", round(joint_mean, digits=2))
end

# ## Part 8: Advanced Quality Metrics
#
# Calculate assembly quality metrics using the package function.

metrics = Mycelia.Rhizomorph.get_qualmer_statistics(qualmer_graph; dataset_id=qualmer_dataset_id)

println("\n" * "="^60)
println("Assembly quality metrics:")
println("  Mean joint quality: ", round(metrics[:mean_joint_quality], digits=2))
println("  Min joint quality: ", metrics[:min_joint_quality])
println("  Max joint quality: ", metrics[:max_joint_quality])
println("  Total unique k-mers: ", metrics[:num_vertices])

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
genome_dataset_id = "genome_demo"
genome_graph = Mycelia.Rhizomorph.build_qualmer_graph(genome_reads, 9; dataset_id=genome_dataset_id, mode=:singlestrand)

# Convert to sequence graph and extract contigs
seq_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(genome_graph)
genome_labels = collect(Mycelia.MetaGraphsNext.labels(genome_graph))
seq_labels = collect(Mycelia.MetaGraphsNext.labels(seq_graph))

println("\nAssembly results:")
println("  Qualmer graph: ", length(genome_labels), " vertices, ", Mycelia.MetaGraphsNext.ne(genome_graph), " edges")
println("  Sequence graph: ", length(seq_labels), " vertices, ", Mycelia.MetaGraphsNext.ne(seq_graph), " edges")

# Use package function to find quality-weighted path
if !isempty(genome_labels)
    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(genome_graph)
    best_path = []
    best_score = 0.0
    for path in paths
        score = Mycelia.Rhizomorph.mean_path_quality(genome_graph, path, genome_dataset_id)
        if score > best_score
            best_score = score
            best_path = path
        end
    end
    
    println("\nBest quality-weighted path:")
    println("  Path length: ", length(best_path), " k-mers")
    
    ## Reconstruct sequence from path
    if length(best_path) > 1
        reconstructed = string(Mycelia.Rhizomorph.path_to_sequence(best_path, genome_graph))
        
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
if !isempty(seq_labels)
    contigs = Mycelia.Rhizomorph.find_contigs_next(seq_graph; min_contig_length=1)
    if isempty(contigs)
        longest_seq = seq_graph[seq_labels[argmax(length.(seq_labels))]].sequence
        longest_qual = seq_graph[seq_labels[argmax(length.(seq_labels))]].quality_scores
    else
        longest_contig = contigs[argmax(length.(getfield.(contigs, :sequence)))]
        longest_seq = longest_contig.sequence
        longest_qual = seq_graph[longest_contig.vertices[1]].quality_scores
    end
    
    println("\nLongest contig from simplified graph:")
    println("  Length: ", length(longest_seq))
    println("  Sequence: ", longest_seq)
    println("  Mean quality: ", round(Statistics.mean(Mycelia.Rhizomorph.decode_quality_scores(longest_qual)), digits=1))
    
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
final_metrics = Mycelia.Rhizomorph.get_qualmer_statistics(genome_graph; dataset_id=genome_dataset_id)
println("\nFinal assembly quality metrics:")
println("  Mean joint quality: ", round(final_metrics[:mean_joint_quality], digits=2))
println("  Min joint quality: ", final_metrics[:min_joint_quality])
println("  Max joint quality: ", final_metrics[:max_joint_quality])

# ## Summary
#
# In this tutorial, we've demonstrated:
#
# 1. **Qualmer Construction**: Building quality-aware k-mer graphs from FASTQ data
# 2. **Joint Probability**: How multiple observations with different qualities are combined
# 3. **Quality vs Coverage**: The advantage of using quality scores over coverage alone
# 4. **Package Functions**: Using Mycelia.Rhizomorph utilities for quality-weighted analysis:
#    - `find_eulerian_paths_next()` + `path_to_sequence()` for path-based reconstruction
#    - `get_qualmer_statistics()` for quality assessment
#    - `find_high_quality_kmers()` for error screening
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
