# # Tutorial 5: Direct FASTQ Sequence Graphs
#
# This tutorial demonstrates the direct workflow from FASTQ sequences to
# quality-aware sequence graphs without the intermediate qualmer step.
# This approach is useful when you want variable-length contigs with
# quality information from the start.

# ## Learning Objectives
# 
# By the end of this tutorial, you will understand:
# 1. How to create quality-aware sequence graphs directly from FASTQ data
# 2. The difference between qualmer-based and direct approaches
# 3. How to perform assembly with preserved per-base quality scores
# 4. When to use direct vs qualmer-mediated approaches
# 5. Quality-aware contig assembly and validation

# ## Setup and Imports
# Following CLAUDE.md standards: only import top-level packages, use full namespacing

import Mycelia
import Test
import Statistics

# ## Part 1: Direct Quality-Aware Sequence Graph Construction
#
# Unlike Tutorial 6 which goes FASTQ → Qualmer → FASTQ graphs,
# this tutorial demonstrates the direct FASTQ → FASTQ graphs approach.

# ### Creating Diverse FASTQ Data

function create_diverse_fastq_data()
    # Simulated reads from a 60bp region with overlaps
    true_sequence = "ATCGATCGATCGTAGCTAGCTAGCTTGCATGCATGCATGCATGCATGCATTAGCTAGC"
    
    reads = []
    
    # High-quality overlapping reads (perfect assembly case)
    push!(reads, FASTX.FASTQ.Record("read1", true_sequence[1:25], "I"^25))
    push!(reads, FASTX.FASTQ.Record("read2", true_sequence[20:45], "I"^26))
    push!(reads, FASTX.FASTQ.Record("read3", true_sequence[40:60], "I"^21))
    
    # Medium-quality reads with slight variations
    push!(reads, FASTX.FASTQ.Record("read4", true_sequence[10:35], "F"^26))
    push!(reads, FASTX.FASTQ.Record("read5", true_sequence[30:55], "F"^26))
    
    # Low-quality read with potential error
    error_seq = true_sequence[15:40]
    error_seq = error_seq[1:10] * "T" * error_seq[12:end]  # Introduce error
    push!(reads, FASTX.FASTQ.Record("read6", error_seq, "AAA###AAA" * "A"^17))
    
    # Short high-quality reads
    push!(reads, FASTX.FASTQ.Record("read7", true_sequence[5:20], "I"^16))
    push!(reads, FASTX.FASTQ.Record("read8", true_sequence[45:60], "I"^16))
    
    return reads, true_sequence
end

println("Creating diverse FASTQ dataset for direct assembly...")
fastq_reads, reference_seq = create_diverse_fastq_data()

println("Dataset overview:")
println("  Reference sequence: ", reference_seq)
println("  Number of reads: ", length(fastq_reads))
for (i, read) in enumerate(fastq_reads)
    seq = String(FASTX.sequence(read))
    qual_scores = FASTX.quality_scores(read)
    mean_qual = round(Statistics.mean(qual_scores), digits=1)
    println("  Read $i: length $(length(seq)), mean quality $mean_qual")
end

# ## Part 2: Direct Quality-Aware BioSequence Graph Construction
#
# Build a quality-aware sequence graph directly from FASTQ reads.

println("\n" * "="^60)
println("Building quality-aware BioSequence graph directly from FASTQ...")

# Use the direct build function
fastq_graph = Mycelia.build_quality_biosequence_graph(fastq_reads)

println("\nDirect FASTQ graph statistics:")
println("  Number of vertices: ", Graphs.nv(fastq_graph))
println("  Number of edges: ", Graphs.ne(fastq_graph))

# Examine the vertices (these should be variable-length sequences)
println("\nExamining quality-aware sequence vertices:")
for v in Iterators.take(Graphs.vertices(fastq_graph), min(5, Graphs.nv(fastq_graph)))
    vertex_data = fastq_graph[v]
    mean_qual = round(Statistics.mean(vertex_data.quality_scores), digits=1)
    println("\nVertex $v:")
    println("  Sequence: ", vertex_data.sequence)
    println("  Length: ", length(vertex_data.sequence))
    println("  Mean quality: ", mean_qual)
    println("  Quality range: ", minimum(vertex_data.quality_scores), "-", maximum(vertex_data.quality_scores))
end

# ## Part 3: Comparison with Qualmer-Mediated Approach
#
# Let's compare the direct approach with the qualmer-mediated approach.

println("\n" * "="^60)
println("Comparing direct vs qualmer-mediated approaches...")

# Build qualmer graph first, then convert to sequence graph
k = 15  # Use larger k for better comparison
qualmer_graph = Mycelia.build_qualmer_graph(fastq_reads; k=k, graph_mode=Mycelia.SingleStrand)
qualmer_to_fastq = Mycelia.qualmer_graph_to_quality_biosequence_graph(qualmer_graph, k)

println("\nApproach comparison:")
println("  Direct FASTQ graph:")
println("    Vertices: ", Graphs.nv(fastq_graph))
println("    Edges: ", Graphs.ne(fastq_graph))

println("  Qualmer-mediated graph (k=$k):")
println("    Qualmer vertices: ", Graphs.nv(qualmer_graph))
println("    Qualmer edges: ", Graphs.ne(qualmer_graph))
println("    Final FASTQ vertices: ", Graphs.nv(qualmer_to_fastq))
println("    Final FASTQ edges: ", Graphs.ne(qualmer_to_fastq))

# ## Part 4: Quality-Aware Assembly Analysis
#
# Analyze assembly quality using both approaches.

println("\n" * "="^60)
println("Analyzing assembly quality for both approaches...")

# For direct approach - analyze sequence lengths and qualities
direct_sequences = [fastq_graph[v].sequence for v in Graphs.vertices(fastq_graph)]
direct_qualities = [fastq_graph[v].quality_scores for v in Graphs.vertices(fastq_graph)]

if !isempty(direct_sequences)
    direct_lengths = [length(seq) for seq in direct_sequences]
    direct_mean_quals = [Statistics.mean(quals) for quals in direct_qualities]
    
    println("\nDirect approach analysis:")
    println("  Sequence count: ", length(direct_sequences))
    println("  Mean sequence length: ", round(Statistics.mean(direct_lengths), digits=1))
    println("  Longest sequence: ", maximum(direct_lengths), " bp")
    println("  Mean quality across all sequences: ", round(Statistics.mean(vcat(direct_mean_quals...)), digits=1))
end

# For qualmer approach - use the package quality metrics
if Graphs.nv(qualmer_graph) > 0
    qualmer_metrics = Mycelia.calculate_assembly_quality_metrics(qualmer_graph)
    
    println("\nQualmer approach analysis:")
    println("  K-mer count: ", qualmer_metrics.total_kmers)
    println("  Mean k-mer coverage: ", round(qualmer_metrics.mean_coverage, digits=1))
    println("  Mean k-mer quality: ", round(qualmer_metrics.mean_quality, digits=1))
    println("  Mean k-mer confidence: ", round(qualmer_metrics.mean_confidence, digits=4))
end

# ## Part 5: Contig Assembly and Reconstruction
#
# Extract contigs from the quality-aware sequence graph.

println("\n" * "="^60)
println("Extracting contigs and performing reconstruction...")

# Find the longest contigs from direct approach
if !isempty(direct_sequences)
    # Sort by length
    seq_length_pairs = [(i, length(seq)) for (i, seq) in enumerate(direct_sequences)]
    sort!(seq_length_pairs, by=x -> x[2], rev=true)
    
    println("\nTop contigs from direct approach:")
    for (rank, (vertex_idx, length)) in enumerate(seq_length_pairs[1:min(3, length(seq_length_pairs))])
        actual_vertex = collect(Graphs.vertices(fastq_graph))[vertex_idx]
        vertex_data = fastq_graph[actual_vertex]
        mean_qual = round(Statistics.mean(vertex_data.quality_scores), digits=1)
        
        println("  Contig $rank:")
        println("    Length: $length bp")
        println("    Mean quality: $mean_qual")
        println("    Sequence: ", vertex_data.sequence)
        
        # Check if this contig matches part of the reference
        contig_seq = String(vertex_data.sequence)
        if occursin(contig_seq, reference_seq)
            println("    ✓ Perfect match in reference")
        elseif occursin(reference_seq, contig_seq)
            println("    ✓ Contains entire reference")
        else
            # Check for partial matches
            best_match_len = 0
            for i in 1:(length(reference_seq) - length(contig_seq) + 1)
                if reference_seq[i:(i + length(contig_seq) - 1)] == contig_seq
                    best_match_len = length(contig_seq)
                    break
                end
            end
            if best_match_len > 0
                println("    ✓ Partial match ($best_match_len bp)")
            else
                println("    ✗ No direct match found")
            end
        end
    end
end

# ## Part 6: Round-Trip Validation
#
# Convert back to FASTQ records and validate quality preservation.

println("\n" * "="^60)
println("Performing round-trip validation...")

# Convert graph back to FASTQ records
reconstructed_fastq = Mycelia.quality_biosequence_graph_to_fastq(fastq_graph, "reconstructed")

println("\nRound-trip validation:")
println("  Original reads: ", length(fastq_reads))
println("  Reconstructed reads: ", length(reconstructed_fastq))

# Analyze quality preservation
original_qualities = []
reconstructed_qualities = []

for read in fastq_reads
    append!(original_qualities, FASTX.quality_scores(read))
end

for read in reconstructed_fastq
    append!(reconstructed_qualities, FASTX.quality_scores(read))
end

if !isempty(original_qualities) && !isempty(reconstructed_qualities)
    println("\nQuality score analysis:")
    println("  Original mean quality: ", round(Statistics.mean(original_qualities), digits=1))
    println("  Reconstructed mean quality: ", round(Statistics.mean(reconstructed_qualities), digits=1))
    println("  Original quality range: ", minimum(original_qualities), "-", maximum(original_qualities))
    println("  Reconstructed quality range: ", minimum(reconstructed_qualities), "-", maximum(reconstructed_qualities))
end

# Compare individual reads
println("\nIndividual read comparison:")
for i in 1:min(3, length(fastq_reads), length(reconstructed_fastq))
    orig = fastq_reads[i]
    recon = reconstructed_fastq[i]
    
    orig_seq = String(FASTX.sequence(orig))
    recon_seq = String(FASTX.sequence(recon))
    
    println("\nRead $i:")
    println("  Original:      ", orig_seq)
    println("  Reconstructed: ", recon_seq)
    println("  Length match: ", length(orig_seq) == length(recon_seq) ? "✓" : "✗")
    println("  Sequence match: ", orig_seq == recon_seq ? "✓" : "✗")
end

# ## Part 7: Error Detection and Quality Assessment
#
# Use quality information to identify potential assembly issues.

println("\n" * "="^60)
println("Quality-based error detection and assessment...")

# Analyze quality distribution across contigs
if !isempty(direct_qualities)
    all_quals = vcat(direct_qualities...)
    quality_stats = (
        mean = Statistics.mean(all_quals),
        median = Statistics.median(all_quals),
        std = Statistics.std(all_quals),
        min = minimum(all_quals),
        max = maximum(all_quals)
    )
    
    println("\nQuality distribution analysis:")
    println("  Mean: ", round(quality_stats.mean, digits=1))
    println("  Median: ", round(quality_stats.median, digits=1))
    println("  Std dev: ", round(quality_stats.std, digits=1))
    println("  Range: ", quality_stats.min, "-", quality_stats.max)
    
    # Identify low-quality regions
    low_quality_threshold = 20.0
    low_qual_count = count(q -> q < low_quality_threshold, all_quals)
    low_qual_fraction = low_qual_count / length(all_quals)
    
    println("\nLow-quality region analysis:")
    println("  Positions below Q$low_quality_threshold: ", low_qual_count)
    println("  Fraction of low-quality positions: ", round(low_qual_fraction, digits=3))
    
    if low_qual_fraction > 0.1
        println("  ⚠️  High fraction of low-quality positions detected")
    else
        println("  ✓ Good overall quality distribution")
    end
end

# ## Part 8: Performance and Use Case Analysis
#
# Discuss when to use direct vs qualmer-mediated approaches.

println("\n" * "="^60)
println("Performance and use case analysis...")

println("\nApproach comparison summary:")
println("\nDirect FASTQ → FASTQ graphs:")
println("  ✓ Faster construction (no intermediate k-mer step)")
println("  ✓ Variable-length sequences from start")
println("  ✓ Natural read-level quality preservation") 
println("  ✓ Good for high-quality, long reads")
println("  ✗ Less granular error detection")
println("  ✗ May struggle with complex repeat regions")

println("\nQualmer-mediated approach:")
println("  ✓ Fine-grained quality analysis at k-mer level")
println("  ✓ Better error detection and correction")
println("  ✓ Handles complex genomic features better")
println("  ✓ Quality-weighted assembly decisions")
println("  ✗ More computationally intensive")
println("  ✗ Requires k-mer size optimization")

println("\nRecommended use cases:")
println("  Direct approach:")
println("    - High-quality PacBio HiFi or ONT reads")
println("    - Simple genomes without complex repeats")
println("    - Rapid prototyping and initial assembly")
println("    - When computational resources are limited")
println("\n  Qualmer approach:")
println("    - Short reads (Illumina)")
println("    - Error-prone long reads")
println("    - Complex genomes with repeats")
println("    - When maximum accuracy is required")

# ## Part 9: Practical Assembly Example
#
# Demonstrate a complete assembly workflow using the direct approach.

println("\n" * "="^60)
println("Complete assembly workflow example...")

function create_realistic_reads()
    # Simulate a 100bp region with realistic read coverage
    true_seq = "ATCGATCGATCGTAGCTAGCTAGCTTGCATGCATGCATGCATGCATGCATTAGCTAGCATCGATCGTAGCTAGCTAGCTTGCATGCATGCAT"
    reads = []
    
    # Simulate 10x coverage with 25bp reads
    read_length = 25
    step_size = 10  # Overlap by 15bp
    
    for start in 1:step_size:(length(true_seq) - read_length + 1)
        end_pos = min(start + read_length - 1, length(true_seq))
        read_seq = true_seq[start:end_pos]
        
        # Vary quality based on position (simulate quality degradation)
        base_quality = 35
        qual_scores = [max(20, base_quality - abs(i - read_length÷2)) for i in 1:length(read_seq)]
        qual_string = String([Char(33 + q) for q in qual_scores])
        
        push!(reads, FASTX.FASTQ.Record("read_$(start)", read_seq, qual_string))
    end
    
    return reads, true_seq
end

realistic_reads, true_genome = create_realistic_reads()

println("\nRealistic assembly example:")
println("  True genome length: ", length(true_genome))
println("  Number of reads: ", length(realistic_reads))
println("  Expected coverage: ~10x")

# Build graph and assemble
realistic_graph = Mycelia.build_quality_biosequence_graph(realistic_reads)
println("  Graph vertices: ", Graphs.nv(realistic_graph))
println("  Graph edges: ", Graphs.ne(realistic_graph))

# Find best assembly
if Graphs.nv(realistic_graph) > 0
    sequences = [realistic_graph[v].sequence for v in Graphs.vertices(realistic_graph)]
    qualities = [realistic_graph[v].quality_scores for v in Graphs.vertices(realistic_graph)]
    
    # Find longest sequence
    longest_idx = argmax(length.(sequences))
    best_assembly = String(sequences[longest_idx])
    best_quality = Statistics.mean(qualities[longest_idx])
    
    println("\nBest assembly result:")
    println("  Assembled length: ", length(best_assembly))
    println("  True length: ", length(true_genome))
    println("  Mean quality: ", round(best_quality, digits=1))
    
    # Check assembly accuracy
    if best_assembly == true_genome
        println("  ✓ Perfect assembly!")
    elseif occursin(best_assembly, true_genome)
        println("  ✓ Assembly is subset of true genome")
        coverage = length(best_assembly) / length(true_genome)
        println("  Coverage: ", round(coverage * 100, digits=1), "%")
    elseif occursin(true_genome, best_assembly)
        println("  ✓ Assembly contains entire true genome")
    else
        println("  ⚠️  Assembly differs from true genome")
        println("  True:      ", true_genome[1:min(50, length(true_genome))], "...")
        println("  Assembled: ", best_assembly[1:min(50, length(best_assembly))], "...")
    end
end

# ## Summary
#
# In this tutorial, we've demonstrated:
#
# 1. **Direct Construction**: Building quality-aware sequence graphs directly from FASTQ
# 2. **Approach Comparison**: Direct vs qualmer-mediated assembly strategies
# 3. **Quality Preservation**: Maintaining per-base quality through assembly
# 4. **Contig Assembly**: Extracting variable-length contigs with quality scores
# 5. **Error Detection**: Using quality information for assembly validation
# 6. **Use Case Analysis**: When to choose each approach
# 7. **Practical Workflow**: Complete assembly pipeline demonstration
#
# Key insights:
# - Direct approach is faster and simpler for high-quality data
# - Qualmer approach provides better error handling for challenging data
# - Quality information is preserved throughout both workflows
# - Choice of approach depends on data quality and computational requirements
# - Both approaches support the complete graph hierarchy for downstream analysis

println("\n" * "="^60)
println("Tutorial 5 completed!")
println("You've learned to choose between direct and qualmer-mediated quality-aware assembly.")
println("Both approaches preserve quality information for downstream analysis.")