"""
# Assembly in 5 Minutes with Mycelia

This tutorial will take you from zero to assembled genome in under 5 minutes.
No prior experience with genome assembly required!

## What You'll Learn
- How to install and load Mycelia
- How to perform your first genome assembly
- How to interpret the results
- How to save your assembly

## Prerequisites
- Julia installed (version 1.10 or higher)
- About 1GB of free disk space
- Internet connection (for package installation)

Let's get started!
"""

# %% [markdown]
# ## Step 1: Installation (One-time setup)
# 
# If you haven't installed Mycelia yet, run this cell:

# %%
# Uncomment the next two lines if you need to install Mycelia
# import Pkg
# Pkg.add(url="https://github.com/cjprybol/Mycelia.git")

# %% [markdown]
# ## Step 2: Load Mycelia
#
# Import the packages we'll need:

# %%
import Mycelia
import FASTX
import Pkg
import Statistics
import Kmers

println("✅ Mycelia loaded successfully!")

# %% [markdown]
# ## Step 3: Get Some Data to Assemble
#
# We'll use a small virus genome (phiX174) as our test case.
# It's only 5,386 base pairs, so it assembles quickly!

# %%
# Download the reference genome
println("📥 Downloading test genome (phiX174)...")
reference_file = Mycelia.download_genome_by_accession(accession="NC_001422.1")
println("✅ Downloaded to: $reference_file")

# Let's see what we downloaded
ref_record = first(Mycelia.open_fastx(reference_file))
println("\n📊 Reference genome info:")
println("  Name: ", FASTX.identifier(ref_record))
println("  Length: ", length(FASTX.sequence(ref_record)), " bp")

# %% [markdown]
# ## Step 4: Simulate Some Sequencing Reads
#
# For this tutorial, we'll simulate some reads. In real life, you'd use
# actual sequencing data from your experiment.

# %%
# Simulate PacBio-style long reads with realistic error rates
println("\n🧬 Simulating sequencing reads...")
reads_file = Mycelia.simulate_pacbio_reads(
    fasta=reference_file,
    quantity="20x"      # 20x coverage
)
println("✅ Created reads file: $reads_file")

# Check the reads
reads = collect(Mycelia.open_fastx(reads_file))
println("\n📊 Simulated reads info:")
println("  Number of reads: ", length(reads))
println("  Average length: ", round(Statistics.mean(length(FASTX.sequence(r)) for r in reads)), " bp")
println("  Total bases: ", sum(length(FASTX.sequence(r)) for r in reads), " bp")

# %% [markdown]
# ## Step 5: Run the Assembly! 🚀
#
# This is where the magic happens. Mycelia will automatically:
# - Detect that you have FASTQ data with quality scores
# - Choose optimal parameters
# - Build quality-aware graphs
# - Find the best assembly path

# %%
println("\n🔧 Starting probabilistic assembly...")
println("  Using intelligent assembly mode (automatic parameter selection)")

# Record start time
start_time = time()

# Run the assembly
assembly_result = Mycelia.mycelia_assemble(
    reads_file;
    output_dir="my_first_assembly",
    max_k=51,  # Maximum k-mer size to try
    memory_limit=4_000_000_000  # 4GB memory limit
)

# Calculate elapsed time
elapsed = round(time() - start_time, digits=1)
println("\n✅ Assembly completed in $elapsed seconds!")

# %% [markdown]
# ## Step 6: Check Your Results
#
# Let's see how well we did:

# %%
# Load the assembled contigs - check what files were created
println("📁 Assembly files created:")
for file in readdir("my_first_assembly")
    println("  - $file")
end

# Extract contigs from assembly result
# The assembly returns a dictionary with :final_assembly containing sequence strings
final_sequences = get(assembly_result, :final_assembly, String[])
println("✅ Assembly produced $(length(final_sequences)) unique sequences")

# Convert strings to FASTA records for further analysis
contigs = FASTX.FASTA.Record[]
for (i, seq) in enumerate(final_sequences)
    if !isempty(seq)  # Skip empty sequences
        record = FASTX.FASTA.Record("contig_$i", seq)
        push!(contigs, record)
    end
end
println("✅ Created $(length(contigs)) FASTA contigs from assembly")

println("\n📊 Assembly Results:")
println("  Number of contigs: ", length(contigs))

# Calculate basic statistics
contig_lengths = [length(FASTX.sequence(c)) for c in contigs]
total_length = sum(contig_lengths)
println("  Total length: ", total_length, " bp")
println("  Longest contig: ", maximum(contig_lengths), " bp")

# Compare to reference
ref_length = length(FASTX.sequence(ref_record))
completeness = round(100 * total_length / ref_length, digits=1)
println("\n🎯 Assembly completeness: $completeness%")

if length(contigs) == 1 && abs(total_length - ref_length) < 100
    println("🌟 Perfect assembly! The genome was assembled into a single contig!")
elseif completeness > 95
    println("🎉 Excellent assembly! Nearly complete genome recovered.")
elseif completeness > 80
    println("👍 Good assembly! Most of the genome was recovered.")
else
    println("🔧 Assembly needs optimization. Try adjusting parameters.")
end

# %% [markdown]
# ## Step 7: Save Your Assembly
#
# Save the assembly for further analysis:

# %%
# Save to a new file with a descriptive name
output_file = "phiX174_mycelia_assembly.fasta"
if !isempty(contigs)
    open(output_file, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for contig in contigs
            write(writer, contig)
        end
    end
    println("\n💾 Assembly saved to: $output_file")
else
    println("\n⚠️  No contigs found in assembly result")
end

# %% [markdown]
# ## Step 8: Visualize Assembly Quality (Optional)
#
# Let's create a simple visualization of our assembly:

# %%
# Assembly quality analysis
println("\n📈 Assembly Quality Analysis:")

# Basic coverage check
println("📊 Coverage Analysis:")
println("  Input reads: $(length(reads)) reads")
println("  Output contigs: $(length(contigs)) contigs")

if !isempty(contigs)
    total_assembled = sum(length(FASTX.sequence(c)) for c in contigs)
    println("  Total assembled length: $total_assembled bp")
    
    # Compare to reference length
    ref_length = length(FASTX.sequence(ref_record))
    recovery_rate = round(100 * total_assembled / ref_length, digits=1)
    println("  Genome recovery rate: $recovery_rate%")
    
    # K-mer analysis comparison (if assembly file was saved successfully)
    if isfile(output_file)
        println("\n🧬 K-mer Analysis:")
        
        # Count k-mers in reference
        ref_kmers = Mycelia.count_canonical_kmers(Kmers.DNAKmer{21}, reference_file)
        println("  Reference unique 21-mers: ", length(ref_kmers))
        
        # Count k-mers in assembly
        asm_kmers = Mycelia.count_canonical_kmers(Kmers.DNAKmer{21}, output_file)
        println("  Assembly unique 21-mers: ", length(asm_kmers))
        
        # Calculate k-mer recovery
        shared_kmers = length(intersect(keys(ref_kmers), keys(asm_kmers)))
        kmer_recovery = round(100 * shared_kmers / length(ref_kmers), digits=1)
        println("  K-mer recovery rate: $kmer_recovery%")
    else
        println("  ⚠️  K-mer analysis skipped - assembly file not available")
    end
else
    println("  ⚠️  No contigs produced - assembly may need debugging")
end

# %% [markdown]
# ## 🎉 Congratulations!
#
# You've just completed your first genome assembly with Mycelia! 
#
# ### What Just Happened?
# 
# Mycelia's intelligent assembly:
# 1. **Analyzed your reads** to determine optimal parameters
# 2. **Built quality-aware graphs** preserving base quality information
# 3. **Used probabilistic algorithms** to find the best assembly path
# 4. **Self-optimized** through multiple k-mer sizes
# 5. **Produced a high-quality assembly** automatically
#
# ### Next Steps
#
# Now that you've seen Mycelia in action, you can:
#
# 1. **Try with your own data** - Replace the simulated reads with your FASTQ file
# 2. **Explore different modes**:
#    - `method=:iterative` for maximum accuracy
#    - `method=:quality_aware` for speed with good quality data
# 3. **Learn more**:
#    - [Understanding the algorithms](../docs/src/theoretical-foundations.md)
#    - [Choosing assembly methods](../docs/src/assembly-method-selection.md)
#    - [Advanced tutorials](../tutorials/)
#
# ### Using Your Own Data
#
# To use your own sequencing data, simply replace the read simulation step:
#
# ```julia
# # Instead of simulating reads:
# reads_file = "path/to/your/reads.fastq"
# 
# # Then run assembly as before:
# assembly_result = Mycelia.mycelia_assemble(reads_file)
# ```
#
# ### Getting Help
#
# - Check the [FAQ](../docs/src/faq.md) for common questions and issues
# - Ask questions in [Discussions](https://github.com/cjprybol/Mycelia/discussions)
# - Report issues on [GitHub](https://github.com/cjprybol/Mycelia/issues)
#
# Happy assembling! 🧬

# %% [markdown]
# ## Appendix: Understanding the Output
#
# The assembly process creates several files in the output directory:
#
# - `final_assembly.fasta` - Your assembled contigs
# - `assembly_summary.txt` - Statistics about the assembly
# - `kmer_progression.txt` - How the assembly progressed through k-mer sizes
# - `graphs/` - Intermediate graph files (for debugging)
#
# The intelligent assembly typically tries several k-mer sizes (e.g., 21, 31, 41)
# and automatically stops when it finds the optimal assembly.