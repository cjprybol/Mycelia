# # Tutorial 1: Data Acquisition in Bioinformatics
#
# This tutorial demonstrates how to acquire genomic data from public databases
# and simulate synthetic datasets for bioinformatics analysis. Data acquisition
# is the critical first step in any bioinformatics workflow.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - How to download genome assemblies from NCBI databases
# - Different types of genomic data files and their purposes
# - How to simulate synthetic genomic data for testing
# - Best practices for data management and reproducibility
# - Quality control considerations for downloaded data

# ## Setup
#
# First, let's load the required packages and set up our environment.

# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/01_data_acquisition.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import Random
import Dates

# Set random seed for reproducibility
Random.seed!(42)

# ## Part 1: Understanding Genomic Data Types
#
# Before downloading data, it's important to understand the different types
# of genomic data files commonly used in bioinformatics:
#
# - **FASTA (.fna, .fasta)**: Raw sequence data (DNA/RNA/protein)
# - **FASTQ (.fastq)**: Sequencing reads with quality scores
# - **GFF3 (.gff3)**: Gene annotations and features
# - **Protein sequences (.faa)**: Translated protein sequences
# - **CDS sequences (.ffn)**: Coding sequences
# - **Sequence reports**: Metadata about sequences

println("=== Genomic Data Types ===")
println("FASTA: Raw sequences (genome assemblies)")
println("FASTQ: Sequencing reads with quality scores")
println("GFF3: Gene annotations and genomic features")
println("Protein: Translated amino acid sequences")
println("CDS: Coding DNA sequences")

# ## Part 2: Downloading Data from NCBI
#
# The National Center for Biotechnology Information (NCBI) is a primary source
# for genomic data. We'll demonstrate downloading a small, well-characterized
# viral genome: bacteriophage phiX174.
#
# ### Why phiX174?
# - Small genome (~5.4 kb) - fast to download and process
# - Well-characterized - extensively studied reference
# - Commonly used as sequencing control
# - Single-stranded DNA virus with overlapping genes

println("\n=== Downloading phiX174 Genome ===")

# ### Method 1: Download by Accession ID
#
# NCBI accession IDs uniquely identify sequence records.
# NC_001422.1 is the RefSeq accession for phiX174.

phix_accession = "NC_001422.1"
println("Downloading genome by accession: $phix_accession")

genome_file = Mycelia.download_genome_by_accession(accession=phix_accession)
println("Downloaded to: $genome_file")

# Verify the download
Test.@test isfile(genome_file)
Test.@test endswith(genome_file, ".fna.gz")
println("✓ Download successful")
println("File size: $(filesize(genome_file)) bytes")

# ### Method 2: Download Complete Assembly Package
#
# For more comprehensive analysis, download the complete assembly package
# which includes multiple file types (genome, proteins, annotations, etc.).

assembly_accession = "GCF_000819615.1"  ## Assembly accession for phiX174
println("\nDownloading complete assembly: $assembly_accession")

assembly_result = Mycelia.ncbi_genome_download_accession(accession=assembly_accession)
println("Assembly downloaded to directory: $(assembly_result.directory)")

# Examine what files were downloaded
println("\nDownloaded files:")
for (key, filepath) in pairs(assembly_result)
    if key != :directory && !isnothing(filepath)
        println("  $key: $(basename(filepath)) ($(filesize(filepath)) bytes)")
    end
end

# ### Examining Downloaded Data
#
# Let's examine the structure of the downloaded genome file.
# FASTA files start with header lines (beginning with '>') followed by sequence data.

println("\n=== Examining Downloaded Genome ===")

# Read the first few lines of the genome file
open(assembly_result.genome) do io
    header = readline(io)
    println("FASTA header: $header")
    
    ## Read first 100 characters of sequence
    sequence_start = read(io, String)[1:min(100, end)]
    println("Sequence start: $sequence_start...")
    
    ## Verify it's a valid FASTA format
    Test.@test startswith(header, '>')
    Test.@test all(c -> c in "ATCGN\n", sequence_start)
end

println("✓ Valid FASTA format confirmed")

# ## Part 3: Simulating Synthetic Data
#
# Synthetic data is crucial for:
# - Testing algorithms with known ground truth
# - Benchmarking performance
# - Developing new methods
# - Educational purposes

println("\n=== Simulating Synthetic Genomic Data ===")

# ### Simulating Random Genome Sequences
#
# Generate a random DNA sequence with specified characteristics.
# This is useful for testing algorithms and understanding their behavior.

println("Generating random DNA sequences...")

# Generate sequences with different lengths and properties
sequences = []

# Short sequence for quick testing
short_seq = Mycelia.random_fasta_record(moltype=:DNA, seed=1, L=100)
push!(sequences, ("Short (100 bp)", short_seq))

# Medium sequence for moderate testing
medium_seq = Mycelia.random_fasta_record(moltype=:DNA, seed=2, L=1000)
push!(sequences, ("Medium (1 kb)", medium_seq))

# Long sequence for performance testing
long_seq = Mycelia.random_fasta_record(moltype=:DNA, seed=3, L=10000)
push!(sequences, ("Long (10 kb)", long_seq))

# Display sequence characteristics
for (name, seq) in sequences
    seq_str = String(FASTX.sequence(seq))
    gc_content = count(c -> c in "GC", seq_str) / length(seq_str)
    println("$name: Length=$(length(seq_str)), GC=$(round(gc_content*100, digits=1))%")
end

# ### Writing Synthetic Data to Files
#
# Save synthetic sequences to files for use in downstream analyses.

println("\nWriting synthetic sequences to files...")

synthetic_files = String[]
for (i, (name, seq)) in enumerate(sequences)
    filename = "synthetic_sequence_$i.fasta"
    Mycelia.write_fasta(outfile=filename, records=[seq])
    push!(synthetic_files, filename)
    println("✓ Wrote $name to $filename")
end

# ## Part 4: Data Quality and Validation
#
# Always validate downloaded and simulated data before analysis.

println("\n=== Data Quality Validation ===")

# ### Validating Downloaded Data
#
# Check file integrity and format correctness.

println("Validating downloaded phiX174 genome...")

# Read and validate the genome sequence
genome_seq = open(assembly_result.genome) do io
    header = readline(io)
    replace(read(io, String), '\n' => "")
end

# Basic validation checks
Test.@test length(genome_seq) > 0
Test.@test all(c -> c in "ATCGN", genome_seq)
Test.@test count(c -> c == 'N', genome_seq) / length(genome_seq) < 0.01  ## < 1% N's

println("✓ Genome length: $(length(genome_seq)) bp")
println("✓ Valid DNA alphabet")
println("✓ Low N content: $(count(c -> c == 'N', genome_seq)) N's")

# Calculate basic statistics
at_content = count(c -> c in "AT", genome_seq) / length(genome_seq)
gc_content = count(c -> c in "GC", genome_seq) / length(genome_seq)

println("AT content: $(round(at_content*100, digits=1))%")
println("GC content: $(round(gc_content*100, digits=1))%")

# ### Validating Synthetic Data
#
# Verify synthetic sequences meet expected criteria.

println("\nValidating synthetic sequences...")

for (i, filename) in enumerate(synthetic_files)
    seq_record = first(Mycelia.open_fastx(filename))
    seq_str = String(FASTX.sequence(seq_record))
    
    ## Validate sequence properties
    Test.@test length(seq_str) > 0
    Test.@test all(c -> c in "ATCG", seq_str)  ## No N's in synthetic data
    
    println("✓ Synthetic sequence $i: $(length(seq_str)) bp, valid DNA")
end

# ## Part 5: Data Management Best Practices
#
# Proper data management ensures reproducibility and efficient analysis.

println("\n=== Data Management Best Practices ===")

# ### Organizing Data Files
#
# Create a structured directory layout for your project.

data_dir = "tutorial_data"
if !isdir(data_dir)
    mkdir(data_dir)
end

# Move downloaded data to organized structure
reference_dir = joinpath(data_dir, "reference")
synthetic_dir = joinpath(data_dir, "synthetic")
mkdir(reference_dir)
mkdir(synthetic_dir)

# Copy files to organized structure (keeping originals for now)
cp(assembly_result.genome, joinpath(reference_dir, "phiX174_genome.fasta"))
cp(assembly_result.gff3, joinpath(reference_dir, "phiX174_annotations.gff3"))

for (i, filename) in enumerate(synthetic_files)
    cp(filename, joinpath(synthetic_dir, filename))
end

println("✓ Data organized in structured directories:")
println("  - Reference data: $reference_dir")
println("  - Synthetic data: $synthetic_dir")

# ### Creating Data Manifest
#
# Document your data sources and processing steps.

manifest_file = joinpath(data_dir, "data_manifest.txt")
open(manifest_file, "w") do io
    println(io, "# Data Manifest - Tutorial 1: Data Acquisition")
    println(io, "# Generated: $(Dates.now())")
    println(io, "")
    println(io, "## Reference Data")
    println(io, "phiX174_genome.fasta - Downloaded from NCBI accession $phix_accession")
    println(io, "phiX174_annotations.gff3 - Downloaded from NCBI assembly $assembly_accession")
    println(io, "")
    println(io, "## Synthetic Data")
    for (i, filename) in enumerate(synthetic_files)
        println(io, "$filename - Random DNA sequence, seed=$(i), length varies")
    end
end

println("✓ Data manifest created: $manifest_file")

# ## Part 6: Advanced Data Acquisition Techniques
#
# For larger-scale analyses, consider these advanced approaches.

println("\n=== Advanced Techniques ===")

# ### Partial Downloads
#
# Download only specific file types to save time and space.

println("Demonstrating partial downloads...")

# Download only genome files (not proteins, annotations, etc.)
genome_only = Mycelia.ncbi_genome_download_accession(
    accession=assembly_accession, 
    include_string="genome"
)

println("✓ Genome-only download completed")
println("  Available files: $(keys(genome_only))")

# ### Batch Processing Considerations
#
# When downloading multiple genomes, consider:
# - Rate limiting to avoid overwhelming servers
# - Error handling for failed downloads
# - Progress tracking for large datasets
# - Disk space management

println("\nBatch processing considerations:")
println("- Implement rate limiting between downloads")
println("- Add retry logic for failed downloads")
println("- Monitor disk space usage")
println("- Use parallel processing judiciously")

# ## Summary and Next Steps
#
# In this tutorial, you learned:
# - How to download genomic data from NCBI databases
# - Different types of genomic data files and their purposes
# - How to simulate synthetic data for testing
# - Data validation and quality control practices
# - Best practices for data organization and management

println("\n=== Tutorial Summary ===")
println("✓ Downloaded phiX174 genome from NCBI")
println("✓ Generated synthetic sequences for testing")
println("✓ Validated data quality and format")
println("✓ Organized data in structured directories")
println("✓ Created data manifest for reproducibility")

# ### Cleanup
#
# Clean up temporary files (keeping organized data)

println("\nCleaning up temporary files...")

# Remove original downloaded files (we have copies in organized structure)
rm(genome_file, force=true)
rm(assembly_result.directory, recursive=true, force=true)

# Remove synthetic files (we have copies in organized structure)
for filename in synthetic_files
    rm(filename, force=true)
end

println("✓ Temporary files cleaned up")

# ### What's Next?
#
# Now that you have acquired genomic data, the next steps typically include:
# 1. **Quality Control**: Assess data quality and preprocess if needed
# 2. **Feature Extraction**: Analyze sequence composition and properties
# 3. **Comparative Analysis**: Compare sequences or genomes
# 4. **Functional Analysis**: Annotate genes and predict functions
#
# Continue with **Tutorial 2: Quality Control and Preprocessing** to learn
# how to assess and improve data quality before analysis.

println("\n=== Next Steps ===")
println("Continue with Tutorial 2: Quality Control and Preprocessing")
println("Data files available in: $data_dir")
println("Ready for downstream analysis!")

# ## Key Takeaways
#
# 1. **Data Acquisition Strategy**: Always start with small, well-characterized datasets
# 2. **Validation is Critical**: Never assume downloaded data is perfect
# 3. **Organization Matters**: Structure your data from the beginning
# 4. **Documentation**: Keep detailed records of data sources and processing
# 5. **Synthetic Data**: Valuable for testing and algorithm development
# 6. **File Formats**: Understand the purpose of different bioinformatics file types

nothing  ## Suppress output in notebook