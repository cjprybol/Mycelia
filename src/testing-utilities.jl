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
    Mycelia.get_test_genome_fasta(;use_ncbi=true, accession="GCF_000819615.1", outdir=nothing, cleanup_at_exit=true) -> NamedTuple

Get a FASTA file for testing, with automatic fallback to simulated genome.

This function attempts to download a real reference genome from NCBI (default: phiX174),
and automatically falls back to a simulated genome if the download fails. This provides
robust test infrastructure that works both with and without network access.

# Arguments
- `use_ncbi::Bool=true`: Whether to attempt NCBI download (set false to always use simulated)
- `accession::String="GCF_000819615.1"`: NCBI assembly accession to download (default: phiX174)
- `outdir::Union{Nothing,AbstractString}=nothing`: Optional download directory (defaults to a temporary directory)
- `cleanup_at_exit::Bool=true`: Register cleanup to run at process exit

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
- When `cleanup_at_exit=true`, cleanup is also registered via `atexit` to run on test completion

# See Also
- [`Mycelia.ncbi_genome_download_accession`](@ref): Download genomes from NCBI
- [`Mycelia.create_test_genome_fasta`](@ref): Create simulated test genomes
"""
function get_test_genome_fasta(;
    use_ncbi::Bool=true,
    accession::String="GCF_000819615.1",
    outdir::Union{Nothing,AbstractString}=nothing,
    cleanup_at_exit::Bool=true
)
    if use_ncbi
        download_root = nothing
        created_tmp = false
        try
            created_tmp = outdir === nothing
            download_root = created_tmp ? mktempdir() : outdir
            outfolder = joinpath(download_root, accession)
            # Clean up any existing partial downloads
            if isdir(outfolder)
                rm(outfolder, recursive=true, force=true)
            end
            dataset = ncbi_genome_download_accession(
                accession=accession, 
                include_string="genome",
                outdir=download_root,
                max_attempts=3,
                initial_retry_delay=10.0
            )
            cleanup_target = created_tmp ? download_root : outfolder
            cleanup_fn = () -> isdir(cleanup_target) && rm(cleanup_target, recursive=true, force=true)
            if cleanup_at_exit
                atexit(cleanup_fn)
            end
            return (
                fasta = dataset.genome,
                source = :ncbi,
                cleanup = cleanup_fn,
                dataset_dir = outfolder
            )
        catch e
            if created_tmp && download_root !== nothing && isdir(download_root)
                rm(download_root, recursive=true, force=true)
            end
            @warn "NCBI download failed after all retries, using simulated genome" exception=e
        end
    end
    
    # Fallback to simulated genome
    # Use same length as phiX174 for consistency
    temp_fasta = create_test_genome_fasta(length=5386, seed=42)
    cleanup_fn = () -> isfile(temp_fasta) && rm(temp_fasta, force=true)
    if cleanup_at_exit
        atexit(cleanup_fn)
    end
    return (
        fasta = temp_fasta,
        source = :simulated,
        cleanup = cleanup_fn,
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

const _BINNING_TEST_INPUTS_CACHE = Ref{Union{Nothing,NamedTuple}}(nothing)

function _build_simulated_binning_contigs(;n_contigs::Int, contig_length::Int, seed::Int)
    @assert n_contigs > 0 "n_contigs must be positive"
    @assert contig_length > 0 "contig_length must be positive"
    rng = StableRNGs.StableRNG(seed)
    records = FASTX.FASTA.Record[]
    length_jitter = max(1, Int(round(contig_length * 0.1)))
    for i in 1:n_contigs
        base_len = contig_length + rand(rng, -length_jitter:length_jitter)
        seq_len = max(base_len, 5000)
        seq = BioSequences.randdnaseq(rng, seq_len)
        record_id = "genome_$(i)"
        push!(records, FASTX.FASTA.Record(record_id, seq))
    end
    return records
end

function _write_simple_gfa_from_fasta(fasta::String, gfa::String)
    mkpath(dirname(gfa))
    open(gfa, "w") do io
        println(io, "H\tVN:Z:1.0")
        for record in Mycelia.open_fastx(fasta)
            record_id = Mycelia.sanitize_fastx_identifier(String(FASTX.identifier(record)))
            sequence = FASTX.sequence(String, record)
            println(io, "S\t$(record_id)\t$(sequence)")
        end
    end
    return gfa
end

function _summarize_bam_depth_jgi(bam::String, depth_file::String)
    if !isfile(depth_file) || filesize(depth_file) == 0
        Mycelia.add_bioconda_env("metabat2")
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metabat2 jgi_summarize_bam_contig_depths --outputDepth $(depth_file) $(bam)`)
    end
    @assert isfile(depth_file) "JGI depth file not created: $(depth_file)"
    return depth_file
end

"""
    prepare_binning_test_inputs(;outdir=nothing, n_contigs=4, contig_length=50_000,
                                seed=7, depth_target=20.0, readset=:illumina_pe150,
                                threads=get_default_threads())

Generate simulated inputs for binning tool integration tests.

Returns a named tuple with contigs, depth table, coverage table, assembly graph,
mapping file, genomes list, and bin directories.
"""
function prepare_binning_test_inputs(;
    outdir::Union{Nothing,String}=nothing,
    n_contigs::Int=4,
    contig_length::Int=50_000,
    seed::Int=7,
    depth_target::Real=20.0,
    readset::Symbol=:illumina_pe150,
    threads::Int=get_default_threads()
)
    inputs_dir = outdir === nothing ? mktempdir() : outdir
    mkpath(inputs_dir)
    if outdir === nothing
        atexit(() -> isdir(inputs_dir) && rm(inputs_dir; recursive=true, force=true))
    end

    contigs_fasta = joinpath(inputs_dir, "contigs.fna")
    if !isfile(contigs_fasta)
        records = _build_simulated_binning_contigs(n_contigs=n_contigs, contig_length=contig_length, seed=seed)
        Mycelia.write_fasta(outfile=contigs_fasta, records=records, gzip=false)
    end

    contig_ids = String[]
    contig_lengths = Int[]
    for record in Mycelia.open_fastx(contigs_fasta)
        push!(contig_ids, String(FASTX.identifier(record)))
        push!(contig_lengths, length(FASTX.sequence(record)))
    end
    reference_table = DataFrames.DataFrame(sequence_id=contig_ids, length=contig_lengths)

    taxonomy_file = joinpath(inputs_dir, "taxonomy.tsv")
    if !isfile(taxonomy_file)
        open(taxonomy_file, "w") do io
            println(io, "contigs\tpredictions")
            for (idx, contig_id) in enumerate(contig_ids)
                println(io, "$(contig_id)\tk__Bacteria;p__Simulated;g__Sim$(idx)")
            end
        end
    end

    sim_outdir = joinpath(inputs_dir, "simulation")
    sim = Mycelia.simulate_metagenome_community(
        reference_fasta=contigs_fasta,
        reference_table=reference_table,
        n_organisms=length(contig_ids),
        depth_target=depth_target,
        abundance_profile=:log_normal,
        readset=readset,
        outdir=sim_outdir,
        selected_ids=contig_ids,
        rng=StableRNGs.StableRNG(seed),
        emit_truth_reads=false,
        cleanup=false
    )

    if readset != :illumina_pe150
        error("Binning test inputs currently require illumina_pe150 readset")
    end
    forward_reads = sim.reads.forward
    reverse_reads = sim.reads.reverse
    @assert forward_reads !== nothing "Forward reads not generated"
    @assert reverse_reads !== nothing "Reverse reads not generated"

    mem_gb = (Int(Sys.total_memory()) / 1e9 * 0.85)
    index_result = Mycelia.minimap_index(
        fasta=contigs_fasta,
        mapping_type="sr",
        mem_gb=mem_gb,
        threads=threads
    )
    if !isfile(index_result.outfile) || filesize(index_result.outfile) == 0
        run(index_result.cmd)
    end

    bam_path = joinpath(inputs_dir, "contigs.minimap2.sorted.bam")
    mapping = Mycelia.minimap_map_with_index(
        fasta=contigs_fasta,
        mapping_type="sr",
        fastq=(forward_reads, reverse_reads),
        index_file=index_result.outfile,
        mem_gb=mem_gb,
        outfile=bam_path,
        threads=threads,
        sorted=true
    )
    if !isfile(mapping.outfile) || filesize(mapping.outfile) == 0
        run(mapping.cmd)
    end

    depth_file = joinpath(inputs_dir, "jgi_depth.tsv")
    _summarize_bam_depth_jgi(mapping.outfile, depth_file)

    coverage_table = joinpath(inputs_dir, "coverm_contig.tsv")
    if !isfile(coverage_table) || filesize(coverage_table) == 0
        Mycelia.run_coverm_contig(
            bam_files=[mapping.outfile],
            output_tsv=coverage_table,
            threads=threads,
            quiet=true
        )
    end

    metacoag_abundance = joinpath(inputs_dir, "metacoag_abundance.tsv")
    if !isfile(metacoag_abundance) || filesize(metacoag_abundance) == 0
        open(coverage_table, "r") do io
            first_line = readline(io)
            first_fields = split(first_line, '\t')
            open(metacoag_abundance, "w") do out
                if isempty(first_fields) || lowercase(first_fields[1]) ∉ ("contig", "contigname")
                    write(out, first_line, '\n')
                end
                for line in eachline(io)
                    write(out, line, '\n')
                end
            end
        end
    end

    assembly_graph = joinpath(inputs_dir, "assembly_graph.gfa")
    if !isfile(assembly_graph)
        _write_simple_gfa_from_fasta(contigs_fasta, assembly_graph)
    end

    genomes = collect(values(sim.per_genome_fastas))
    bins_root = joinpath(inputs_dir, "bins")
    bins_a = joinpath(bins_root, "bins_a")
    bins_b = joinpath(bins_root, "bins_b")
    mkpath(bins_a)
    mkpath(bins_b)
    for (idx, genome) in enumerate(genomes)
        bin_name = "bin_$(idx).fna"
        cp(genome, joinpath(bins_a, bin_name); force=true)
        cp(genome, joinpath(bins_b, bin_name); force=true)
    end

    return (
        outdir=inputs_dir,
        contigs_fasta=contigs_fasta,
        depth_file=depth_file,
        coverage_table=coverage_table,
        marker_file=nothing,
        taxonomy_file=taxonomy_file,
        assembly_graph=assembly_graph,
        mapping_file=metacoag_abundance,
        genomes=genomes,
        bins_dirs=[bins_a, bins_b]
    )
end

"""
    get_binning_test_inputs(;force=false, kwargs...)

Cache and return binning integration test inputs.
"""
function get_binning_test_inputs(;force::Bool=false, kwargs...)
    if !force && _BINNING_TEST_INPUTS_CACHE[] !== nothing
        return _BINNING_TEST_INPUTS_CACHE[]
    end
    inputs = prepare_binning_test_inputs(;kwargs...)
    _BINNING_TEST_INPUTS_CACHE[] = inputs
    return inputs
end
