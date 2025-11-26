# Testing Utilities Module
#
# This module provides utility functions specifically designed for testing workflows.
# These functions handle common testing scenarios like generating test genomes,
# downloading reference data with fallbacks, and creating test FASTQ records.
#
# Note: These functions are NOT exported. Use them via Mycelia.function_name()

"""
    Mycelia.create_test_genome_fasta(;length=5386, seed=42) -> String

Create a temporary FASTA file with a simulated genome for testing.

This function provides a reliable way to generate test genomic data without
requiring network access or external dependencies. Useful for unit tests that
need reproducible genome sequences.

# Arguments
- `length::Int=5386`: Length of the simulated genome in base pairs (default: phiX174 genome size)
- `seed::Int=42`: Random seed for reproducible sequence generation

# Returns
- `String`: Path to the created temporary FASTA file

# Details
Uses the canonical `random_fasta_record()` and `write_fasta()` functions to ensure
consistency with the rest of Mycelia's API.

# Example
```julia
# Create a 10kb test genome
test_genome = Mycelia.create_test_genome_fasta(length=10000, seed=123)

# Use it with read simulators
reads = Mycelia.simulate_illumina_reads(fasta=test_genome, coverage=30)

# Clean up when done
rm(test_genome)
```

# See Also
- [`Mycelia.random_fasta_record`](@ref): Generate random FASTA records
- [`Mycelia.get_test_genome_fasta`](@ref): Download real genome with simulated fallback
"""
function create_test_genome_fasta(;length::Int=5386, seed::Int=42)
    record = random_fasta_record(moltype=:DNA, seed=seed, L=length)
    temp_fasta = tempname() * ".fasta"
    write_fasta(outfile=temp_fasta, records=[record])
    return temp_fasta
end

"""
    Mycelia.get_test_genome_fasta(;use_ncbi=true, accession="GCF_000819615.1") -> NamedTuple

Get a FASTA file for testing, with automatic fallback to simulated genome.

This function attempts to download a real reference genome from NCBI (default: phiX174),
and automatically falls back to a simulated genome if the download fails. This provides
robust test infrastructure that works both with and without network access.

# Arguments
- `use_ncbi::Bool=true`: Whether to attempt NCBI download (set false to always use simulated)
- `accession::String="GCF_000819615.1"`: NCBI assembly accession to download (default: phiX174)

# Returns
A `NamedTuple` with the following fields:
- `fasta::String`: Path to the FASTA file
- `source::Symbol`: Either `:ncbi` (downloaded) or `:simulated` (generated locally)
- `cleanup::Function`: Function to call for cleanup (removes downloaded/generated files)
- `dataset_dir::Union{String,Nothing}`: Directory containing downloaded data (or `nothing` for simulated)

# Details
- Uses `ncbi_genome_download_accession()` with retry logic for robust downloads
- Automatically cleans up partial downloads before retrying
- Falls back to `create_test_genome_fasta()` if download fails
- Provides cleanup function for resource management

# Example
```julia
# Get a test genome (tries NCBI download, falls back to simulation)
genome_info = Mycelia.get_test_genome_fasta()
println("Using \$(genome_info.source) genome: \$(genome_info.fasta)")

# Use the genome
reads = Mycelia.simulate_illumina_reads(fasta=genome_info.fasta, coverage=50)

# Clean up when done
genome_info.cleanup()

# Force use of simulated genome (no network access)
sim_genome = Mycelia.get_test_genome_fasta(use_ncbi=false)
```

# Notes
- Default accession (GCF_000819615.1) is phiX174, a 5.4kb bacteriophage genome
- Simulated fallback uses the same length as the default phiX174 genome (5386 bp)
- Cleanup function is always provided, but checks if files exist before deletion

# See Also
- [`Mycelia.ncbi_genome_download_accession`](@ref): Download genomes from NCBI
- [`Mycelia.create_test_genome_fasta`](@ref): Create simulated test genomes
"""
function get_test_genome_fasta(;use_ncbi::Bool=true, accession::String="GCF_000819615.1")
    if use_ncbi
        try
            # Clean up any existing partial downloads
            if isdir(accession)
                rm(accession, recursive=true)
            end
            dataset = ncbi_genome_download_accession(
                accession=accession, 
                include_string="genome",
                max_attempts=3,
                initial_retry_delay=10.0
            )
            return (
                fasta = dataset.genome,
                source = :ncbi,
                cleanup = () -> isdir(accession) && rm(accession, recursive=true),
                dataset_dir = accession
            )
        catch e
            @warn "NCBI download failed after all retries, using simulated genome" exception=e
        end
    end
    
    # Fallback to simulated genome
    # Use same length as phiX174 for consistency
    temp_fasta = create_test_genome_fasta(length=5386, seed=42)
    return (
        fasta = temp_fasta,
        source = :simulated,
        cleanup = () -> isfile(temp_fasta) && rm(temp_fasta),
        dataset_dir = nothing
    )
end

"""
    Mycelia.create_test_reads(reference_sequence, coverage, error_rate) -> Vector{FASTX.FASTQ.Record}

Create FASTQ records from a DNA reference sequence with simulated sequencing errors.

This function generates multiple reads from a reference sequence by applying the
`observe()` function to introduce realistic sequencing errors and quality scores.

# Arguments
- `reference_sequence::Union{String, BioSequences.LongDNA{4}}`: Reference sequence to generate reads from
- `coverage::Int`: Number of reads to generate (simulates coverage depth)
- `error_rate::Float64`: Error rate for sequencing simulation (0.0 = no errors, 1.0 = 100% errors)

# Returns
- `Vector{FASTX.FASTQ.Record}`: Vector of FASTQ records with simulated errors and quality scores

# Details
- Uses `Mycelia.observe()` to introduce realistic sequencing errors
- Generates Phred+33 quality scores for each base
- Read IDs are numbered sequentially as "read_1", "read_2", etc.

# Example
```julia
# Generate 100 reads with 1% error rate from a reference sequence
ref_seq = BioSequences.dna"ACGTACGTACGTACGT"
reads = Mycelia.create_test_reads(ref_seq, 100, 0.01)

# Write to FASTQ file
Mycelia.write_fastq(records=reads, filename="test_reads.fastq")
```

# See Also
- [`Mycelia.observe`](@ref): Introduce sequencing errors and generate quality scores
- [`Mycelia.create_test_rna_reads`](@ref): RNA-specific version
- [`Mycelia.create_test_aa_reads`](@ref): Amino acid-specific version
"""
function create_test_reads(reference_sequence::Union{String, BioSequences.LongDNA{4}}, coverage::Int, error_rate::Float64)
    records = FASTX.FASTQ.Record[]

    for i in 1:coverage
        # Use the observe function to introduce errors
        bio_seq = reference_sequence isa String ? BioSequences.LongDNA{4}(reference_sequence) : reference_sequence
        observed_seq, quality_scores = observe(bio_seq, error_rate=error_rate)
        
        # Convert quality scores to string format (Phred+33)
        quality_string = String([Char(q + 33) for q in quality_scores])
        
        record = FASTX.FASTQ.Record("read_$i", string(observed_seq), quality_string)
        push!(records, record)
    end
    
    return records
end

"""
    Mycelia.create_test_rna_reads(reference_sequence, coverage, error_rate) -> Vector{FASTX.FASTQ.Record}

Create FASTQ records from an RNA reference sequence with simulated sequencing errors.

Similar to `create_test_reads()` but specifically for RNA sequences.

# Arguments
- `reference_sequence::Union{String, BioSequences.LongRNA{4}}`: RNA reference sequence
- `coverage::Int`: Number of reads to generate
- `error_rate::Float64`: Error rate for sequencing simulation

# Returns
- `Vector{FASTX.FASTQ.Record}`: Vector of FASTQ records with simulated errors and quality scores

# Example
```julia
# Generate RNA-seq reads
rna_seq = BioSequences.rna"ACGUACGUACGUACGU"
reads = Mycelia.create_test_rna_reads(rna_seq, 50, 0.01)
```

# See Also
- [`Mycelia.create_test_reads`](@ref): DNA-specific version
- [`Mycelia.observe`](@ref): Core sequencing error simulation
"""
function create_test_rna_reads(reference_sequence::Union{String, BioSequences.LongRNA{4}}, coverage::Int, error_rate::Float64)
    records = FASTX.FASTQ.Record[]
    
    for i in 1:coverage
        # Use the observe function to introduce errors
        bio_seq = reference_sequence isa String ? BioSequences.LongRNA{4}(reference_sequence) : reference_sequence
        observed_seq, quality_scores = observe(bio_seq, error_rate=error_rate)
        
        # Convert quality scores to string format (Phred+33)
        quality_string = String([Char(q + 33) for q in quality_scores])
        
        record = FASTX.FASTQ.Record("read_$i", string(observed_seq), quality_string)
        push!(records, record)
    end
    
    return records
end

"""
    Mycelia.create_test_aa_reads(reference_sequence, coverage, error_rate) -> Vector{FASTX.FASTQ.Record}

Create FASTQ records from an amino acid reference sequence with simulated sequencing errors.

Similar to `create_test_reads()` but specifically for protein/amino acid sequences.
Includes validation to ensure no termination characters (*) are present.

# Arguments
- `reference_sequence::Union{String, BioSequences.LongAA}`: Amino acid reference sequence
- `coverage::Int`: Number of reads to generate
- `error_rate::Float64`: Error rate for sequencing simulation

# Returns
- `Vector{FASTX.FASTQ.Record}`: Vector of FASTQ records with simulated errors and quality scores

# Throws
- `ErrorException`: If the observed sequence contains termination character '*'

# Example
```julia
# Generate protein sequencing reads
aa_seq = BioSequences.aa"ACDEFGHIKLMNPQRSTVWY"
reads = Mycelia.create_test_aa_reads(aa_seq, 100, 0.01)
```

# See Also
- [`Mycelia.create_test_reads`](@ref): DNA-specific version
- [`Mycelia.observe`](@ref): Core sequencing error simulation
"""
function create_test_aa_reads(reference_sequence::Union{String, BioSequences.LongAA}, coverage::Int, error_rate::Float64)
    records = FASTX.FASTQ.Record[]

    for i in 1:coverage
        # Use the observe function to introduce errors
        bio_seq = reference_sequence isa String ? BioSequences.LongAA(reference_sequence) : reference_sequence
        observed_seq, quality_scores = observe(bio_seq, error_rate=error_rate)

        # Validate that no termination characters (*) are present in amino acid sequences
        # This ensures compatibility with FASTQ format
        observed_seq_str = string(observed_seq)
        if occursin('*', observed_seq_str)
            error("Amino acid sequence contains termination character '*' which is invalid for FASTQ format: $observed_seq_str")
        end

        # Convert quality scores to string format (Phred+33)
        quality_string = String([Char(q + 33) for q in quality_scores])

        record = FASTX.FASTQ.Record("read_$i", observed_seq_str, quality_string)
        push!(records, record)
    end
    
    return records
end
