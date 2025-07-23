using Pkg
Pkg.activate(".")

import Mycelia
import Statistics

println("=== Basic Tutorial Test ===")

# Create synthetic data
synthetic_genome = "GGAAACCTGGAGCGAACTGGATCCCCGCCTCCTTTTGTGGGCCTCCGGCGCTGTGAGCTCTCTACGACCCGCCCAGCCAG"
println("Synthetic genome length: $(length(synthetic_genome))")

# Simulate reads
dna_reads = Mycelia._simulate_fastq_reads_from_sequence(
    synthetic_genome, 
    "test_viroid";
    coverage = 10,
    read_length = 30,
    error_rate = 0.01,
    sequence_type = "DNA"
)

println("Generated $(length(dna_reads)) DNA reads")

# Quality analysis
if !isempty(dna_reads)
    first_read = dna_reads[1]
    quality_scores = collect(Mycelia.FASTX.quality_scores(first_read))
    mean_qual = Statistics.mean(quality_scores)
    println("Mean quality of first read: $(round(mean_qual, digits=2))")
end

# Assembly
config = Mycelia.AssemblyConfig(k=8, use_quality_scores=true)
observations = [(read, i) for (i, read) in enumerate(dna_reads)]
result = Mycelia._assemble_qualmer_graph(observations, config)

println("Assembly results:")
println("  String contigs: $(length(result.contigs))")
println("  FASTQ contigs: $(length(result.fastq_contigs))")

# Check quality preservation
if haskey(result.assembly_stats, "quality_preserved")
    quality_preserved = result.assembly_stats["quality_preserved"] 
    println("  Quality preserved: $quality_preserved")
end

if haskey(result.assembly_stats, "mean_quality")
    mean_quality = result.assembly_stats["mean_quality"]
    println("  Mean assembly quality: $(round(mean_quality, digits=2))")
end

println("✓ Basic tutorial test completed successfully!")