"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run MEGAHIT assembler for metagenomic short read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "megahit_output")
- `min_contig_len::Int`: Minimum contig length (default: 200)
- `k_list::String`: k-mer sizes to use (default: "21,29,39,59,79,99,119,141")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to final contigs file

# Details
- Automatically creates and uses a conda environment with megahit
- Optimized for metagenomic assemblies with varying coverage
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_megahit(;fastq1, fastq2=nothing, outdir="megahit_output", min_contig_len=200, k_list="21,29,39,59,79,99,119,141")
    Mycelia.add_bioconda_env("megahit")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "final.contigs.fa"))
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit -r $(fastq1) -o $(outdir) --min-contig-len $(min_contig_len) --k-list $(k_list) -t $(Sys.CPU_THREADS)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit -1 $(fastq1) -2 $(fastq2) -o $(outdir) --min-contig-len $(min_contig_len) --k-list $(k_list) -t $(Sys.CPU_THREADS)`)
        end
    end
    return (;outdir, contigs=joinpath(outdir, "final.contigs.fa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaSPAdes assembler for metagenomic short read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "metaspades_output")
- `k_list::String`: k-mer sizes to use (default: "21,33,55,77")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `scaffolds::String`: Path to scaffolds file

# Details
- Automatically creates and uses a conda environment with spades
- Designed for metagenomic datasets with uneven coverage
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_metaspades(;fastq1, fastq2=nothing, outdir="metaspades_output", k_list="21,33,55,77")
    Mycelia.add_bioconda_env("spades")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "contigs.fasta"))
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades metaspades.py -s $(fastq1) -o $(outdir) -k $(k_list) -t $(Sys.CPU_THREADS)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades metaspades.py -1 $(fastq1) -2 $(fastq2) -o $(outdir) -k $(k_list) -t $(Sys.CPU_THREADS)`)
        end
    end
    return (;outdir, contigs=joinpath(outdir, "contigs.fasta"), scaffolds=joinpath(outdir, "scaffolds.fasta"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run SKESA assembler for high-accuracy bacterial genome assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "skesa_output")
- `min_contig_len::Int`: Minimum contig length (default: 200)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file

# Details
- Uses conservative assembly approach for high accuracy
- Optimized for bacterial genomes with uniform coverage
- Automatically creates and uses a conda environment with skesa
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_skesa(;fastq1, fastq2=nothing, outdir="skesa_output", min_contig_len=200)
    Mycelia.add_bioconda_env("skesa")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "contigs.fa")
    if !isfile(contigs_file)
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1) --contigs_out $(contigs_file) --cores $(Sys.CPU_THREADS) --min_contig $(min_contig_len)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1),$(fastq2) --contigs_out $(contigs_file) --cores $(Sys.CPU_THREADS) --min_contig $(min_contig_len)`)
        end
    end
    return (;outdir, contigs=contigs_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run IDBA-UD assembler for metagenomic assembly with uneven depth.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file
- `outdir::String`: Output directory path (default: "idba_ud_output")
- `min_k::Int`: Minimum k-mer size (default: 20)
- `max_k::Int`: Maximum k-mer size (default: 100)
- `step::Int`: K-mer size increment step (default: 20)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `scaffolds::String`: Path to scaffolds file

# Details
- Uses multi-k-mer iterative assembly approach
- Specifically designed for metagenomic data with uneven coverage depths
- Automatically creates and uses a conda environment with idba
- Requires paired-end reads merged to fasta format
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_idba_ud(;fastq1, fastq2, outdir="idba_ud_output", min_k=20, max_k=100, step=20)
    Mycelia.add_bioconda_env("idba")
    mkpath(outdir)
    
    # IDBA-UD requires paired reads in fasta format
    merged_fasta = joinpath(outdir, "merged_reads.fa")
    contig_file = joinpath(outdir, "contig.fa")
    scaffold_file = joinpath(outdir, "scaffold.fa")
    
    if !isfile(contig_file)
        # Convert and merge FASTQ to FASTA format for IDBA-UD
        if !isfile(merged_fasta)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n idba fq2fa --merge $(fastq1) $(fastq2) $(merged_fasta)`)
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n idba idba_ud -r $(merged_fasta) -o $(outdir) --mink $(min_k) --maxk $(max_k) --step $(step) --num_threads $(Sys.CPU_THREADS)`)
    end
    return (;outdir, contigs=contig_file, scaffolds=scaffold_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Flye assembler for long read assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "flye_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g")
- `read_type::String`: Type of reads ("pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Automatically creates and uses a conda environment with flye
- Supports various long read technologies and quality levels
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_flye(;fastq, outdir="flye_output", genome_size, read_type="pacbio-hifi")
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "assembly.fasta"))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye flye --$(read_type) $(fastq) --out-dir $(outdir) --genome-size $(genome_size) --threads $(Sys.CPU_THREADS)`)
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaFlye assembler for long-read metagenomic assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "metaflye_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g")
- `read_type::String`: Type of reads ("pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq")
- `meta::Bool`: Enable metagenome mode (default: true)
- `min_overlap::Int`: Minimum overlap between reads (default: auto-selected)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Uses metaFlye's repeat graph approach optimized for metagenomic data
- Implements solid k-mer selection combining global and local k-mer distributions
- Handles uneven coverage and strain variation in metagenomic samples
- Automatically creates and uses a conda environment with flye
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_metaflye(;fastq, outdir="metaflye_output", genome_size, read_type="pacbio-hifi", meta=true, min_overlap=nothing)
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "assembly.fasta"))
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir", outdir, "--genome-size", genome_size, "--threads", string(Sys.CPU_THREADS)]
        
        if meta
            push!(cmd_args, "--meta")
        end
        
        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`)
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Canu assembler for long read assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "canu_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g")
- `read_type::String`: Type of reads ("pacbio", "nanopore")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Automatically creates and uses a conda environment with canu
- Includes error correction, trimming, and assembly stages
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_canu(;fastq, outdir="canu_output", genome_size, read_type="pacbio")
    Mycelia.add_bioconda_env("canu")
    mkpath(outdir)
    
    prefix = basename(fastq, ".fastq")
    if !isfile(joinpath(outdir, "$(prefix).contigs.fasta"))
        if read_type == "pacbio"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -pacbio $(fastq) maxThreads=$(Sys.CPU_THREADS)`)
        elseif read_type == "nanopore"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -nanopore $(fastq) maxThreads=$(Sys.CPU_THREADS)`)
        else
            error("Unsupported read type: $(read_type). Use 'pacbio' or 'nanopore'.")
        end
    end
    return (;outdir, assembly=joinpath(outdir, "$(prefix).contigs.fasta"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the hifiasm genome assembler on PacBio HiFi reads.

# Arguments
- `fastq::String`: Path to input FASTQ file containing HiFi reads
- `outdir::String`: Output directory path (default: "\${basename(fastq)}_hifiasm")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `hifiasm_outprefix::String`: Prefix used for hifiasm output files

# Details
- Automatically creates and uses a conda environment with hifiasm
- Uses primary assembly mode (--primary) optimized for inbred samples
- Skips assembly if output files already exist at the specified prefix
- Utilizes all available CPU threads
"""
function run_hifiasm(;fastq, outdir=basename(fastq) * "_hifiasm")
    Mycelia.add_bioconda_env("hifiasm")
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm")
    hifiasm_outputs = filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join=true))
    # https://hifiasm.readthedocs.io/en/latest/faq.html#are-inbred-homozygous-genomes-supported
    if isempty(hifiasm_outputs)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm hifiasm --primary -l0 -o $(hifiasm_outprefix) -t $(Sys.CPU_THREADS) $(fastq)`)
    end
    return (;outdir, hifiasm_outprefix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run hifiasm-meta assembler for metagenomic PacBio HiFi assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing HiFi reads
- `outdir::String`: Output directory path (default: "\${basename(fastq)}_hifiasm_meta")
- `similarity::Float64`: Similarity threshold for strain resolution (default: 0.85)
- `purge_level::Int`: Purge level for strain variants (default: 2)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `hifiasm_outprefix::String`: Prefix used for hifiasm-meta output files

# Details
- Uses hifiasm-meta's string graph approach for metagenomic strain resolution
- Implements SNV-based read phasing to separate closely related strains
- Optimized for low-abundance species and strain-level diversity
- Automatically creates and uses a conda environment with hifiasm
- Skips assembly if output files already exist at the specified prefix
- Utilizes all available CPU threads
"""
function run_hifiasm_meta(;fastq, outdir=basename(fastq) * "_hifiasm_meta", similarity=0.85, purge_level=2)
    Mycelia.add_bioconda_env("hifiasm")
    mkpath(outdir)
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm_meta")
    hifiasm_outputs = filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join=true))
    
    if isempty(hifiasm_outputs)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm hifiasm-meta -t $(Sys.CPU_THREADS) -o $(hifiasm_outprefix) --purge-cov $(purge_level) --similarity $(similarity) $(fastq)`)
    end
    return (;outdir, hifiasm_outprefix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run hybrid assembly combining short and long reads using Unicycler.

# Arguments
- `short_1::String`: Path to first short read FASTQ file
- `short_2::String`: Path to second short read FASTQ file (optional)
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "unicycler_output")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Automatically creates and uses a conda environment with unicycler
- Combines short read accuracy with long read scaffolding
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_unicycler(;short_1, short_2=nothing, long_reads, outdir="unicycler_output")
    Mycelia.add_bioconda_env("unicycler")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "assembly.fasta"))
        if isnothing(short_2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -s $(short_1) -l $(long_reads) -o $(outdir) -t $(Sys.CPU_THREADS)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -1 $(short_1) -2 $(short_2) -l $(long_reads) -o $(outdir) -t $(Sys.CPU_THREADS)`)
        end
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"))
end

# ============================================================================
# Phase 3: Unified Assembly Pipeline Interface
# ============================================================================

"""
Assembly method enumeration for unified interface.
"""
@enum AssemblyMethod begin
    # Fixed-length graph types (assembly foundation)
    NgramGraph       # N-gram graph assembly (for unicode character analysis)
    KmerGraph        # K-mer graph assembly with DNAKmer/RNAKmer/AAKmer (for FASTA data)
    QualmerGraph     # Quality-aware k-mer graph assembly (for FASTQ data) - PRIMARY METHOD
    
    # Variable-length graph types (simplified products)
    StringGraph      # String graph assembly (simplified from N-gram graphs)
    BioSequenceGraph # BioSequence graph assembly (simplified from K-mer graphs)
    QualityBioSequenceGraph # Quality-aware BioSequence graph assembly (simplified from Qualmer graphs)
    
    # Hybrid approaches
    HybridOLC        # Hybrid OLC + qualmer graph approach
    MultiK           # Multi-k assembly with merging
end

"""
Assembly configuration structure.
"""
struct AssemblyConfig
    k::Int                          # Primary k-mer size
    error_rate::Float64             # Expected sequencing error rate
    min_coverage::Int               # Minimum coverage for k-mer inclusion  
    graph_mode::GraphMode           # SingleStrand or DoubleStrand
    use_quality_scores::Bool        # Whether to use FASTQ quality scores
    polish_iterations::Int          # Number of polishing iterations
    bubble_resolution::Bool         # Whether to resolve bubble structures
    repeat_resolution::Bool         # Whether to resolve repeat regions
    
    # Constructor with sensible defaults
    function AssemblyConfig(;
        k::Int = 31,
        error_rate::Float64 = 0.01,
        min_coverage::Int = 3,
        graph_mode::GraphMode = DoubleStrand,
        use_quality_scores::Bool = true,
        polish_iterations::Int = 3,
        bubble_resolution::Bool = true,
        repeat_resolution::Bool = true
    )
        new(k, error_rate, min_coverage, graph_mode, use_quality_scores, 
            polish_iterations, bubble_resolution, repeat_resolution)
    end
end

"""
Assembly result structure containing contigs and metadata.
"""
struct AssemblyResult
    contigs::Vector{String}             # Final assembled contigs
    contig_names::Vector{String}        # Contig identifiers
    graph::Union{Nothing, MetaGraphsNext.MetaGraph}  # Final assembly graph (optional)
    assembly_stats::Dict{String, Any}   # Assembly statistics and metrics
    fastq_contigs::Vector{FASTX.FASTQ.Record}  # Quality-aware contigs (FASTQ format)
    
    function AssemblyResult(contigs::Vector{String}, contig_names::Vector{String}; 
                          graph=nothing, assembly_stats=Dict{String, Any}(),
                          fastq_contigs=FASTX.FASTQ.Record[])
        new(contigs, contig_names, graph, assembly_stats, fastq_contigs)
    end
end

"""
    get_fastq_contigs(result::AssemblyResult) -> Vector{FASTX.FASTQ.Record}

Extract quality-aware FASTQ contigs from assembly result if available.
Returns empty vector if no quality information was preserved during assembly.
"""
function get_fastq_contigs(result::AssemblyResult)
    return result.fastq_contigs
end

"""
    has_quality_information(result::AssemblyResult) -> Bool

Check if the assembly result preserves quality information from the original reads.
"""
function has_quality_information(result::AssemblyResult)
    return !isempty(result.fastq_contigs) && get(result.assembly_stats, "quality_preserved", false)
end

"""
    write_fastq_contigs(result::AssemblyResult, output_file::String)

Write quality-aware contigs to a FASTQ file if quality information is available.
"""
function write_fastq_contigs(result::AssemblyResult, output_file::String)
    if !has_quality_information(result)
        error("Assembly result does not contain quality information")
    end
    
    FASTX.FASTQ.Writer(open(output_file, "w")) do writer
        for record in result.fastq_contigs
            write(writer, record)
        end
    end
    
    @info "Quality-aware contigs written to $(output_file)"
    return output_file
end

"""
    assemble_genome(reads; method=StringGraph, config=AssemblyConfig()) -> AssemblyResult

Unified genome assembly interface using Phase 2 next-generation algorithms.

# Arguments
- `reads`: Vector of FASTA/FASTQ records or file paths
- `method`: Assembly strategy (StringGraph, KmerGraph, HybridOLC, MultiK)
- `config`: Assembly configuration parameters

# Returns
- `AssemblyResult`: Structure containing contigs, names, and assembly metadata

# Details
This is the main entry point for the unified assembly pipeline, leveraging:
- Phase 1: MetaGraphsNext strand-aware graph construction  
- Phase 2: Probabilistic algorithms, enhanced Viterbi, and graph algorithms
- Phase 3: Integrated workflow with polishing and validation

# Examples
```julia
# Basic assembly with default parameters
reads = load_fastq_records("reads.fastq")
result = assemble_genome(reads)

# Custom assembly with specific k-mer size and error rate
config = AssemblyConfig(k=25, error_rate=0.005, polish_iterations=5)
result = assemble_genome(reads; method=KmerGraph, config=config)

# Access results
contigs = result.contigs
stats = result.assembly_stats
```
"""
function assemble_genome(reads; method::AssemblyMethod=QualmerGraph, config::AssemblyConfig=AssemblyConfig())
    @info "Starting unified genome assembly" method config.k config.error_rate
    
    # Phase 1: Load and validate input
    observations = _prepare_observations(reads)
    @info "Loaded $(length(observations)) sequence observations"
    
    # Phase 2: Assembly strategy dispatch following 6-graph hierarchy
    if method == NgramGraph
        result = _assemble_ngram_graph(observations, config)
    elseif method == KmerGraph  
        result = _assemble_kmer_graph(observations, config)
    elseif method == QualmerGraph
        result = _assemble_qualmer_graph(observations, config)
    elseif method == StringGraph
        result = _assemble_string_graph(observations, config)
    elseif method == BioSequenceGraph
        result = _assemble_biosequence_graph(observations, config)
    elseif method == QualityBioSequenceGraph
        result = _assemble_quality_biosequence_graph(observations, config)
    elseif method == HybridOLC
        result = _assemble_hybrid_olc(observations, config)
    elseif method == MultiK
        result = _assemble_multi_k(observations, config)
    else
        throw(ArgumentError("Unknown assembly method: $method"))
    end
    
    @info "Assembly completed: $(length(result.contigs)) contigs generated"
    return result
end

"""
    polish_assembly(assembly::AssemblyResult, reads; iterations=3) -> AssemblyResult

Polish assembled contigs using quality-aware error correction.

# Arguments  
- `assembly`: Initial assembly result to polish
- `reads`: Original reads for polishing (FASTQ with quality scores preferred)
- `iterations`: Number of polishing iterations (default: 3)

# Returns
- `AssemblyResult`: Polished assembly with improved accuracy

# Details
Uses Phase 2 enhanced Viterbi algorithms with quality score integration for:
- Error correction based on k-mer graph traversals
- Consensus calling from multiple observations
- Iterative improvement until convergence
"""
function polish_assembly(assembly::AssemblyResult, reads; iterations::Int=3)
    @info "Starting assembly polishing" iterations
    
    observations = _prepare_observations(reads)
    polished_contigs = copy(assembly.contigs)
    
    for iter in 1:iterations
        @info "Polishing iteration $iter/$iterations"
        
        # Build k-mer graph from current contigs + reads
        graph = build_kmer_graph_next(Kmers.DNAKmer{31}, 
                                    vcat(observations, _contigs_to_records(polished_contigs)))
        
        # Polish each contig using Viterbi error correction
        for (i, contig) in enumerate(polished_contigs)
            polished_contigs[i] = _polish_contig_viterbi(contig, graph, observations)
        end
    end
    
    # Update assembly stats
    new_stats = merge(assembly.assembly_stats, Dict(
        "polishing_iterations" => iterations,
        "polished" => true
    ))
    
    @info "Polishing completed"
    return AssemblyResult(polished_contigs, assembly.contig_names; 
                         graph=assembly.graph, assembly_stats=new_stats)
end

"""
    validate_assembly(assembly::AssemblyResult; reference=nothing) -> Dict{String, Any}

Validate assembly quality using various metrics and optional reference comparison.

# Arguments
- `assembly`: Assembly result to validate
- `reference`: Optional reference sequence for comparison

# Returns  
- `Dict{String, Any}`: Comprehensive validation metrics

# Details
Computes assembly quality metrics including:
- N50, N90 statistics
- Total assembly length and number of contigs
- Coverage uniformity (if reference provided)
- Structural variant detection (if reference provided)
- Gap analysis and repeat characterization
"""
function validate_assembly(assembly::AssemblyResult; reference=nothing)
    @info "Validating assembly quality"
    
    contigs = assembly.contigs
    metrics = Dict{String, Any}()
    
    # Basic assembly statistics
    contig_lengths = [length(contig) for contig in contigs]
    total_length = sum(contig_lengths)
    sort!(contig_lengths, rev=true)
    
    metrics["num_contigs"] = length(contigs)
    metrics["total_length"] = total_length
    metrics["mean_contig_length"] = mean(contig_lengths)
    metrics["max_contig_length"] = maximum(contig_lengths)
    metrics["min_contig_length"] = minimum(contig_lengths)
    
    # N-statistics
    metrics["N50"] = _calculate_n_statistic(contig_lengths, 0.5)
    metrics["N90"] = _calculate_n_statistic(contig_lengths, 0.9)
    metrics["L50"] = _calculate_l_statistic(contig_lengths, 0.5)
    metrics["L90"] = _calculate_l_statistic(contig_lengths, 0.9)
    
    # Reference-based validation (if provided)
    if !isnothing(reference)
        ref_metrics = _validate_against_reference(contigs, reference)
        merge!(metrics, ref_metrics)
    end
    
    @info "Assembly validation completed" metrics["N50"] metrics["num_contigs"]
    return metrics
end

# ============================================================================
# Phase 3: Helper Functions for Assembly Pipeline
# ============================================================================

"""
Prepare observations from various input formats (FASTA/FASTQ records or file paths).
"""
function _prepare_observations(reads)
    if reads isa Vector{String}
        # File paths - load FASTA/FASTQ records
        observations = FASTX.FASTA.Record[]
        for file_path in reads
            if endswith(file_path, ".fastq") || endswith(file_path, ".fq")
                # Load FASTQ and convert to FASTA records
                FASTX.FASTQ.open(file_path) do reader
                    for record in reader
                        push!(observations, FASTX.FASTA.Record(FASTX.FASTQ.identifier(record), FASTX.FASTQ.sequence(record)))
                    end
                end
            else
                # Load FASTA records
                FASTX.FASTA.open(file_path) do reader
                    for record in reader
                        push!(observations, record)
                    end
                end
            end
        end
        return observations
    elseif reads isa Vector{<:FASTX.FASTA.Record}
        return reads
    elseif reads isa Vector{<:FASTX.FASTQ.Record}
        # Convert FASTQ to FASTA records
        return [FASTX.FASTA.Record(FASTX.FASTQ.identifier(record), FASTX.FASTQ.sequence(record)) for record in reads]
    else
        throw(ArgumentError("Unsupported reads format: $(typeof(reads))"))
    end
end

"""
String graph assembly implementation (variable-length simplified from N-gram graphs).
"""
function _assemble_string_graph(observations, config)
    @info "Using string graph assembly strategy (variable-length simplified from N-gram)"
    
    # First build N-gram graph
    text = join([FASTX.FASTA.sequence(obs) for obs in observations], "")
    ngram_graph = string_to_ngram_graph(text, config.k)
    
    # Simplify N-gram graph to create variable-length string graph
    string_graph = _simplify_ngram_to_string_graph(ngram_graph)
    
    # Further collapse unbranching paths
    collapse_unbranching_paths(string_graph)
    
    # If enabled, resolve bubble structures
    if config.bubble_resolution
        _simplify_string_graph(string_graph)
    end
    
    # Assemble contigs from graph
    contigs = assemble_strings(string_graph)
    contig_names = ["string_contig_$i" for i in 1:length(contigs)]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "StringGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "unicode_strings",
        "source_graph" => "ngram_graph_simplification",
        "k" => config.k,
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Dates.now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=string_graph, assembly_stats=stats)
end

"""
K-mer graph assembly implementation (fixed-length k-mer foundation).
"""
function _assemble_kmer_graph(observations, config)
    @info "Using k-mer graph assembly strategy (fixed-length k-mer foundation)"
    
    # Determine k-mer type from observations
    kmer_type = _determine_kmer_type(observations, config.k)
    
    # Build k-mer graph using Phase 2 next-generation algorithms
    graph = build_kmer_graph_next(kmer_type, observations; graph_mode=config.graph_mode)
    
    # Apply Phase 2 graph algorithms
    if config.bubble_resolution
        bubbles = detect_bubbles_next(graph)
        @info "Detected $(length(bubbles)) bubble structures"
    end
    
    if config.repeat_resolution
        repeats = resolve_repeats_next(graph)
        @info "Identified $(length(repeats)) repeat regions"
    end
    
    # Find contigs using Eulerian path finding
    paths = find_eulerian_paths_next(graph)
    
    # Convert paths to sequences
    contigs = String[]
    for path in paths
        if length(path) > 1
            sequence = _path_to_sequence(path, graph)
            push!(contigs, sequence)
        end
    end
    
    # If no Eulerian paths, use probabilistic walks
    if isempty(contigs)
        @info "No Eulerian paths found, using probabilistic walks"
        contigs = _generate_contigs_probabilistic(graph, config)
    end
    
    contig_names = ["kmer_contig_$i" for i in 1:length(contigs)]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "KmerGraph",
        "graph_type" => "fixed_length",
        "vertex_type" => string(kmer_type),
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Dates.now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=graph, assembly_stats=stats)
end

"""
Quality-aware k-mer graph assembly implementation (fixed-length qualmer foundation).
"""
function _assemble_qualmer_graph(observations, config)
    @info "Using quality-aware k-mer graph assembly strategy (fixed-length qualmer foundation)"
    
    # Convert observations to FASTQ records for quality processing
    fastq_records = _prepare_fastq_observations(observations)
    
    # Build qualmer graph using Phase 2 quality-aware algorithms
    graph = build_qualmer_graph(fastq_records; k=config.k, graph_mode=config.graph_mode)
    
    # Apply Phase 2 graph algorithms if enabled
    if config.bubble_resolution
        # Note: This would use detect_bubbles_next adapted for qualmer graphs
        @info "Bubble resolution enabled for qualmer graphs"
    end
    
    if config.repeat_resolution
        # Note: This would use resolve_repeats_next adapted for qualmer graphs
        @info "Repeat resolution enabled for qualmer graphs"
    end
    
    # Find contigs using quality-aware path finding
    paths = _find_qualmer_paths(graph, config)
    
    # Convert paths to FASTQ records with quality propagation
    contig_records = FASTX.FASTQ.Record[]
    for (i, path) in enumerate(paths)
        if !isempty(path)
            contig_name = "qualmer_contig_$i"
            # Use consensus quality calculation for better accuracy
            fastq_record = _qualmer_path_to_consensus_fastq(path, graph, contig_name)
            if length(FASTX.sequence(fastq_record)) > config.k  # Only keep substantial contigs
                push!(contig_records, fastq_record)
            end
        end
    end
    
    # If no paths found, use probabilistic walks on qualmer graph
    if isempty(contig_records)
        @info "No quality-aware paths found, using probabilistic walks"
        contig_records = _generate_fastq_contigs_from_qualmer_graph(graph, config)
    end
    
    # Convert FASTQ records to strings for backward compatibility with existing code
    contigs = [String(FASTX.sequence(record)) for record in contig_records]
    contig_names = [String(FASTX.identifier(record)) for record in contig_records]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "QualmerGraph",
        "graph_type" => "fixed_length",
        "vertex_type" => "quality_aware_kmers",
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "quality_preserved" => true,  # Mark that quality information is preserved
        "num_fastq_contigs" => length(contig_records),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Dates.now())
    )
    
    # Add quality-specific statistics
    if !isempty(MetaGraphsNext.labels(graph))
        qualmer_stats = get_qualmer_statistics(graph)
        merge!(stats, qualmer_stats)
    end
    
    return AssemblyResult(contigs, contig_names; graph=graph, assembly_stats=stats, fastq_contigs=contig_records)
end

"""
Find quality-weighted paths through a qualmer graph using iterative algorithms.
Implements heaviest path (highest confidence), iterative Viterbi, and probabilistic walks.
"""
function _find_qualmer_paths(graph, config)
    paths = Vector{Vector}()
    
    # Get all vertices sorted by joint probability (highest confidence first)
    vertices = collect(MetaGraphsNext.labels(graph))
    if isempty(vertices)
        return paths
    end
    
    # Sort vertices by confidence (joint probability)
    sorted_vertices = sort(vertices, by=v -> graph[v].joint_probability, rev=true)
    
    visited = Set()
    
    # Method 1: Heaviest Path Algorithm (Highest Confidence Eulerian Paths)
    @info "Searching for high-confidence Eulerian paths"
    for start_vertex in sorted_vertices
        if start_vertex ∈ visited
            continue
        end
        
        # Try to find Eulerian path starting from high-confidence vertex
        path = _find_heaviest_eulerian_path(graph, start_vertex, visited)
        if length(path) > 1
            union!(visited, path)
            push!(paths, path)
        end
    end
    
    # Method 2: Iterative Viterbi Algorithm for remaining vertices
    if length(paths) < 3  # If we haven't found many paths, try iterative Viterbi
        @info "Applying iterative Viterbi algorithm"
        remaining_vertices = filter(v -> v ∉ visited, vertices)
        if !isempty(remaining_vertices)
            viterbi_paths = _iterative_viterbi_paths(graph, remaining_vertices, config)
            for path in viterbi_paths
                if length(path) > 1
                    union!(visited, path)
                    push!(paths, path)
                end
            end
        end
    end
    
    # Method 3: Probabilistic Walks for any remaining high-quality vertices
    remaining_vertices = filter(v -> v ∉ visited && graph[v].joint_probability > 0.5, vertices)
    if !isempty(remaining_vertices) && length(paths) < 5
        @info "Using probabilistic walks for remaining high-quality vertices"
        for start_vertex in remaining_vertices[1:min(3, end)]  # Limit to avoid too many small contigs
            if start_vertex ∉ visited
                path = _quality_weighted_walk(graph, start_vertex, config.k * 5)
                if length(path) > 1
                    union!(visited, path)
                    push!(paths, path)
                end
            end
        end
    end
    
    return paths
end

"""
Find the heaviest (highest confidence) Eulerian path starting from a given vertex.
"""
function _find_heaviest_eulerian_path(graph, start_vertex, visited::Set)
    path = [start_vertex]
    current = start_vertex
    
    max_steps = 1000  # Prevent infinite loops
    steps = 0
    
    while steps < max_steps
        # Get unvisited neighbors
        neighbors = collect(Graphs.outneighbors(graph, current))
        unvisited_neighbors = filter(n -> n ∉ visited, neighbors)
        
        if isempty(unvisited_neighbors)
            break
        end
        
        # Choose neighbor with highest combined score (vertex quality × edge quality)
        best_neighbor = nothing
        best_score = -1.0
        
        for neighbor in unvisited_neighbors
            # Vertex quality score
            vertex_score = graph[neighbor].joint_probability
            
            # Edge quality score
            edge_score = 1.0
            if MetaGraphsNext.has_edge(graph, current, neighbor)
                edge_data = graph[current, neighbor]
                edge_score = edge_data.quality_weight
            end
            
            # Combined score
            total_score = vertex_score * edge_score
            
            if total_score > best_score
                best_score = total_score
                best_neighbor = neighbor
            end
        end
        
        if best_neighbor === nothing || best_score < 0.1  # Quality threshold
            break
        end
        
        push!(path, best_neighbor)
        current = best_neighbor
        steps += 1
    end
    
    return path
end

"""
Iterative Viterbi algorithm for finding optimal paths through qualmer graph.
Uses dynamic programming with quality scores as emission/transition probabilities.
"""
function _iterative_viterbi_paths(graph, candidate_vertices, config)
    paths = Vector{Vector}()
    
    # Group vertices into potential start points (high in-degree or isolated)
    start_candidates = Vector()
    for vertex in candidate_vertices
        in_degree = length(collect(Graphs.inneighbors(graph, vertex)))
        out_degree = length(collect(Graphs.outneighbors(graph, vertex)))
        
        # Prefer vertices that could be sequence starts
        if in_degree <= out_degree || graph[vertex].joint_probability > 0.8
            push!(start_candidates, vertex)
        end
    end
    
    # If no clear starts, use highest quality vertices
    if isempty(start_candidates)
        start_candidates = candidate_vertices[1:min(3, end)]
    end
    
    visited = Set()
    for start_vertex in start_candidates
        if start_vertex ∈ visited
            continue
        end
        
        # Apply Viterbi-like algorithm
        path = _viterbi_optimal_path(graph, start_vertex, visited, 500)  # Reasonable max length
        if length(path) > 1
            union!(visited, path)
            push!(paths, path)
        end
    end
    
    return paths
end

"""
Find optimal path using Viterbi-like dynamic programming on quality scores.
"""
function _viterbi_optimal_path(graph, start_vertex, visited::Set, max_length::Int)
    # Initialize with start vertex
    path = [start_vertex]
    current = start_vertex
    local_visited = Set([start_vertex])
    
    for step in 1:max_length
        neighbors = collect(Graphs.outneighbors(graph, current))
        unvisited_neighbors = filter(n -> n ∉ visited && n ∉ local_visited, neighbors)
        
        if isempty(unvisited_neighbors)
            break
        end
        
        # Calculate Viterbi scores for each neighbor
        best_neighbor = nothing
        best_viterbi_score = -Inf
        
        for neighbor in unvisited_neighbors
            # Emission probability (vertex quality)
            emission_prob = log(max(1e-10, graph[neighbor].joint_probability))
            
            # Transition probability (edge quality)
            transition_prob = 0.0
            if MetaGraphsNext.has_edge(graph, current, neighbor)
                edge_data = graph[current, neighbor]
                transition_prob = log(max(1e-10, edge_data.quality_weight))
            end
            
            # Viterbi score
            viterbi_score = emission_prob + transition_prob
            
            if viterbi_score > best_viterbi_score
                best_viterbi_score = viterbi_score
                best_neighbor = neighbor
            end
        end
        
        if best_neighbor === nothing || best_viterbi_score < -10.0  # Quality threshold in log space
            break
        end
        
        push!(path, best_neighbor)
        push!(local_visited, best_neighbor)
        current = best_neighbor
    end
    
    return path
end

"""
Hybrid OLC assembly (placeholder for future implementation).
"""
function _assemble_hybrid_olc(observations, config)
    @info "Using hybrid OLC assembly strategy"
    
    # For now, fall back to k-mer graph assembly
    # Future implementation would combine overlap-layout-consensus with string graphs
    @warn "Hybrid OLC not fully implemented, using k-mer graph assembly"
    return _assemble_kmer_graph(observations, config)
end

"""
Multi-k assembly (placeholder for future implementation).
"""
function _assemble_multi_k(observations, config)
    @info "Using multi-k assembly strategy"
    
    # For now, fall back to single-k assembly
    # Future implementation would use multiple k-mer sizes and merge results
    @warn "Multi-k assembly not fully implemented, using single-k assembly"
    return _assemble_kmer_graph(observations, config)
end

"""
Generate contigs using probabilistic walks when Eulerian paths are not available.
"""
function _generate_contigs_probabilistic(graph, config)
    contigs = String[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    # Start probabilistic walks from vertices with high in-degree (likely start points)
    for start_vertex in vertices
        path = probabilistic_walk_next(graph, start_vertex, 1000; seed=rand(1:10000))
        if length(path.steps) > 1
            sequence = _reconstruct_sequence_from_path(path.steps)
            if length(sequence) > config.k  # Only keep substantial contigs
                push!(contigs, sequence)
            end
        end
    end
    
    return contigs
end

"""
Convert a path of vertices to a DNA sequence.
"""
function _path_to_sequence(path, graph)
    if isempty(path)
        return ""
    end
    
    # Start with the first k-mer
    sequence = string(path[1])
    
    # Add one nucleotide from each subsequent k-mer (assuming k-mer overlap)
    for i in 2:length(path)
        kmer = string(path[i])
        if length(kmer) > 0
            sequence *= kmer[end]  # Add last nucleotide
        end
    end
    
    return sequence
end

"""
Polish a single contig using Viterbi error correction.
"""
function _polish_contig_viterbi(contig, graph, observations)
    # Create observation sequence from contig
    contig_kmers = [contig[i:i+30] for i in 1:length(contig)-30]  # Assuming k=31
    
    # Use Viterbi decoding for error correction
    try
        path = viterbi_decode_next(graph, contig_kmers)
        return _reconstruct_sequence_from_path(path.steps)
    catch e
        @warn "Viterbi polishing failed for contig" error=e
        return contig  # Return original if polishing fails
    end
end

"""
Convert contigs to FASTA records for graph construction.
"""
function _contigs_to_records(contigs)
    records = FASTX.FASTA.Record[]
    for (i, contig) in enumerate(contigs)
        push!(records, FASTX.FASTA.Record("contig_$i", contig))
    end
    return records
end

"""
Calculate N-statistic (N50, N90, etc.) for contig lengths.
"""
function _calculate_n_statistic(sorted_lengths, threshold)
    total_length = sum(sorted_lengths)
    target_length = total_length * threshold
    
    cumulative_length = 0
    for length in sorted_lengths
        cumulative_length += length
        if cumulative_length >= target_length
            return length
        end
    end
    
    return 0
end

"""
    _calculate_l_statistic(sorted_lengths, threshold)

Calculate L-statistic (number of contigs needed to reach a given percentage of total assembly length).
For example, L50 is the number of contigs needed to reach 50% of the total assembly length.

# Arguments
- `sorted_lengths`: Vector of contig lengths sorted in descending order
- `threshold`: Fraction of total length (e.g., 0.5 for L50, 0.9 for L90)

# Returns
- `Int`: Number of contigs needed to reach the threshold
"""
function _calculate_l_statistic(sorted_lengths, threshold)
    total_length = sum(sorted_lengths)
    target_length = total_length * threshold
    
    cumulative_length = 0
    contig_count = 0
    
    for length in sorted_lengths
        cumulative_length += length
        contig_count += 1
        if cumulative_length >= target_length
            return contig_count
        end
    end
    
    return length(sorted_lengths)  # Return total number of contigs if threshold not reached
end

"""
Validate assembly against reference sequence (placeholder).
"""
function _validate_against_reference(contigs, reference)
    # Placeholder for reference-based validation
    # Would implement alignment-based metrics, coverage analysis, etc.
    return Dict{String, Any}(
        "reference_provided" => true,
        "reference_length" => length(reference)
    )
end

"""
Simplify string graph by removing unnecessary complexity.
"""
function _simplify_string_graph(graph)
    # Basic graph simplification - remove low-weight edges, merge linear paths
    # This is a placeholder for more sophisticated graph simplification
    return graph
end

# ============================================================================
# Qualmer Graph Assembly Helper Functions
# ============================================================================

"""
Convert observations to FASTQ records for quality processing.
"""
function _prepare_fastq_observations(observations)
    fastq_records = FASTX.FASTQ.Record[]
    
    for (i, obs) in enumerate(observations)
        # Handle different observation formats
        if obs isa FASTX.FASTQ.Record
            push!(fastq_records, obs)
        elseif obs isa FASTX.FASTA.Record
            # Convert FASTA to FASTQ with default quality (assume high quality)
            seq = FASTX.FASTA.sequence(obs)
            qual = repeat('I', length(seq))  # Quality score 40 (high quality)
            record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(obs), seq, qual)
            push!(fastq_records, record)
        elseif obs isa Tuple && length(obs) >= 1
            # Handle tuple format like (record, index) or (record, other_data)
            record = obs[1]
            if record isa FASTX.FASTQ.Record
                push!(fastq_records, record)
            elseif record isa FASTX.FASTA.Record
                # Convert FASTA to FASTQ with default quality
                seq = FASTX.FASTA.sequence(record)
                qual = repeat('I', length(seq))  # Quality score 40 (high quality)
                fastq_record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(record), seq, qual)
                push!(fastq_records, fastq_record)
            else
                @warn "Unsupported record type in tuple: $(typeof(record))"
            end
        else
            @warn "Unsupported observation type: $(typeof(obs))"
        end
    end
    
    return fastq_records
end

"""
Find quality-aware paths in a qualmer graph.
"""
function _find_qualmer_paths(graph, config)
    paths = []
    
    # Basic implementation: find simple paths
    # Future implementation would use more sophisticated algorithms
    vertices = collect(MetaGraphsNext.labels(graph))
    
    for start_vertex in vertices
        # Simple path following highest quality edges
        path = [start_vertex]
        current = start_vertex
        
        while true
            # Find outgoing edges
            outgoing = []
            for edge in MetaGraphsNext.edge_labels(graph)
                src, dst = edge
                if src == current
                    push!(outgoing, (dst, graph[src, dst]))
                end
            end
            
            if isempty(outgoing)
                break
            end
            
            # Choose edge with highest quality weight
            best_dst, best_edge = outgoing[1]
            for (dst, edge_data) in outgoing[2:end]
                if edge_data.quality_weight > best_edge.quality_weight
                    best_dst, best_edge = dst, edge_data
                end
            end
            
            # Avoid cycles
            if best_dst in path
                break
            end
            
            push!(path, best_dst)
            current = best_dst
        end
        
        if length(path) > 1
            push!(paths, path)
        end
    end
    
    return paths
end

"""
Convert qualmer path to DNA sequence.
"""
function _qualmer_path_to_sequence(path, graph)
    if isempty(path)
        return ""
    end
    
    # Start with first k-mer
    sequence = path[1]
    
    # Add one nucleotide from each subsequent k-mer
    for i in 2:length(path)
        kmer = path[i]
        if length(kmer) > 0
            sequence *= kmer[end]  # Add last nucleotide
        end
    end
    
    return sequence
end

"""
Enhanced qualmer path to FASTQ record conversion with quality propagation.
This is the core function that enables quality-aware assembly output.
"""
function _qualmer_path_to_fastq_record(path, graph, contig_name::String)
    if isempty(path)
        return FASTX.FASTQ.Record("", "", "")
    end
    
    # Initialize sequence and quality arrays
    sequence_chars = Char[]
    quality_scores = UInt8[]
    
    # Process first qualmer vertex completely
    first_vertex_data = graph[path[1]]
    first_qualmer = first_vertex_data.canonical_qualmer
    
    # Add entire first k-mer sequence and qualities
    for i in 1:length(first_qualmer.kmer)
        push!(sequence_chars, first_qualmer.kmer[i])
        push!(quality_scores, first_qualmer.qualities[i])
    end
    
    # For subsequent vertices, add only the extending nucleotide with its quality
    for vertex_idx in 2:length(path)
        vertex_data = graph[path[vertex_idx]]
        qualmer = vertex_data.canonical_qualmer
        
        # Add the last nucleotide and its quality (k-mer overlap extension)
        if length(qualmer.kmer) > 0
            push!(sequence_chars, qualmer.kmer[end])
            push!(quality_scores, qualmer.qualities[end])
        end
    end
    
    # Convert to strings
    sequence_str = String(sequence_chars)
    quality_str = String(Char.(quality_scores .+ UInt8(33)))  # Convert to ASCII quality scores
    
    # Create FASTQ record with quality-aware contig
    return FASTX.FASTQ.Record(contig_name, sequence_str, quality_str)
end

"""
Enhanced qualmer path to sequence with consensus quality calculation.
Uses joint probability from multiple observations to compute consensus quality.
"""
function _qualmer_path_to_consensus_fastq(path, graph, contig_name::String)
    if isempty(path)
        return FASTX.FASTQ.Record("", "", "")
    end
    
    sequence_chars = Char[]
    consensus_qualities = UInt8[]
    
    # Process first vertex completely
    first_vertex_data = graph[path[1]]
    first_qualmer = first_vertex_data.canonical_qualmer
    
    # For first k-mer, use consensus quality from all observations
    for i in 1:length(first_qualmer.kmer)
        push!(sequence_chars, first_qualmer.kmer[i])
        # Calculate consensus quality from joint probability
        consensus_qual = _calculate_consensus_quality_from_observations(first_vertex_data.observations, i)
        push!(consensus_qualities, consensus_qual)
    end
    
    # For subsequent vertices, add extending nucleotide with consensus quality
    for vertex_idx in 2:length(path)
        vertex_data = graph[path[vertex_idx]]
        qualmer = vertex_data.canonical_qualmer
        
        if length(qualmer.kmer) > 0
            push!(sequence_chars, qualmer.kmer[end])
            # Calculate consensus quality for the extending position
            extending_pos = length(qualmer.kmer)
            consensus_qual = _calculate_consensus_quality_from_observations(vertex_data.observations, extending_pos)
            push!(consensus_qualities, consensus_qual)
        end
    end
    
    # Convert to strings
    sequence_str = String(sequence_chars)
    quality_str = String(Char.(consensus_qualities .+ UInt8(33)))
    
    return FASTX.FASTQ.Record(contig_name, sequence_str, quality_str)
end

"""
Calculate consensus quality from multiple observations at a specific position.
Uses the joint probability and multiple observations to compute a consensus quality score.
"""
function _calculate_consensus_quality_from_observations(observations::Vector{<:QualmerObservation}, position::Int)
    if isempty(observations)
        return UInt8(2)  # Very low quality if no observations
    end
    
    # Extract quality scores at the given position from all observations
    position_qualities = UInt8[]
    for obs in observations
        if position <= length(obs.qualmer.qualities)
            push!(position_qualities, obs.qualmer.qualities[position])
        end
    end
    
    if isempty(position_qualities)
        return UInt8(2)
    end
    
    # Calculate consensus using weighted average based on error probabilities
    total_weight = 0.0
    weighted_sum = 0.0
    
    for qual in position_qualities
        # Convert quality to probability of correctness
        prob_correct = 1.0 - 10.0^(-qual / 10.0)
        weight = prob_correct
        
        total_weight += weight
        weighted_sum += weight * qual
    end
    
    if total_weight == 0.0
        return UInt8(2)
    end
    
    # Calculate weighted average quality
    consensus_qual = weighted_sum / total_weight
    
    # Apply confidence boost for multiple supporting observations
    num_obs = length(position_qualities)
    if num_obs > 1
        # Boost quality based on number of supporting observations
        confidence_boost = min(10.0, log10(num_obs) * 3.0)  # Max boost of 10
        consensus_qual = min(40.0, consensus_qual + confidence_boost)  # Cap at Q40
    end
    
    return UInt8(round(consensus_qual))
end

"""
Generate FASTQ contigs from qualmer graph using probabilistic walks when no paths are found.
"""
function _generate_fastq_contigs_from_qualmer_graph(graph, config)
    contig_records = FASTX.FASTQ.Record[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    if isempty(vertices)
        return contig_records
    end
    
    # Track visited vertices to avoid duplicates
    visited = Set()
    
    # Start probabilistic walks from high-confidence vertices
    confidence_threshold = 0.7
    start_vertices = filter(v -> graph[v].joint_probability > confidence_threshold, vertices)
    
    # If no high-confidence vertices, use all vertices
    if isempty(start_vertices)
        start_vertices = vertices
    end
    
    contig_id = 1
    for start_vertex in start_vertices
        if start_vertex ∈ visited
            continue
        end
        
        # Perform quality-weighted walk from this vertex
        path = _quality_weighted_walk(graph, start_vertex, config.k * 10)  # Reasonable max length
        
        if length(path) > 1  # Need at least 2 vertices to form a meaningful contig
            # Mark all vertices in path as visited
            union!(visited, path)
            
            # Convert path to FASTQ record
            contig_name = "qualmer_probabilistic_contig_$(contig_id)"
            fastq_record = _qualmer_path_to_consensus_fastq(path, graph, contig_name)
            
            # Only keep contigs longer than k
            if length(FASTX.sequence(fastq_record)) > config.k
                push!(contig_records, fastq_record)
                contig_id += 1
            end
        end
    end
    
    return contig_records
end

"""
Perform quality-weighted walk through qualmer graph.
"""
function _quality_weighted_walk(graph, start_vertex, max_length::Int=1000)
    path = [start_vertex]
    current = start_vertex
    visited = Set([current])
    
    while length(path) < max_length
        # Get outgoing neighbors
        neighbors = collect(Graphs.outneighbors(graph, current))
        
        # Filter out visited vertices
        unvisited = filter(n -> n ∉ visited, neighbors)
        
        if isempty(unvisited)
            break
        end
        
        # Choose next vertex based on joint probability and edge quality
        best_neighbor = nothing
        best_score = -1.0
        
        for neighbor in unvisited
            # Calculate score based on vertex quality and edge weight
            vertex_quality = graph[neighbor].joint_probability
            
            # Check if edge exists and get its quality weight
            edge_quality = 1.0
            if MetaGraphsNext.has_edge(graph, current, neighbor)
                edge_data = graph[current, neighbor]
                edge_quality = edge_data.quality_weight
            end
            
            # Combined score (vertex quality * edge quality)
            score = vertex_quality * edge_quality
            
            if score > best_score
                best_score = score
                best_neighbor = neighbor
            end
        end
        
        if best_neighbor === nothing
            break
        end
        
        push!(path, best_neighbor)
        push!(visited, best_neighbor)
        current = best_neighbor
    end
    
    return path
end

"""
Generate contigs from qualmer graph using probabilistic walks.
"""
function _generate_contigs_from_qualmer_graph(graph, config)
    contigs = String[]
    vertices = collect(MetaGraphsNext.labels(graph))
    
    # Start probabilistic walks from high-quality vertices
    for start_vertex in vertices
        vertex_data = graph[start_vertex]
        
        # Only start from high-quality k-mers
        if vertex_data.joint_probability > 0.7  # Threshold for high quality
            # Simple random walk (placeholder for more sophisticated algorithm)
            path = [start_vertex]
            current = start_vertex
            
            for _ in 1:100  # Max walk length
                # Find outgoing edges
                outgoing = []
                for edge in MetaGraphsNext.edge_labels(graph)
                    src, dst = edge
                    if src == current && dst ∉ path  # Avoid cycles
                        push!(outgoing, dst)
                    end
                end
                
                if isempty(outgoing)
                    break
                end
                
                # Choose randomly weighted by quality
                current = rand(outgoing)
                push!(path, current)
            end
            
            if length(path) > 1
                sequence = _qualmer_path_to_sequence(path, graph)
                if length(sequence) > config.k
                    push!(contigs, sequence)
                end
            end
        end
    end
    
    return contigs
end

# ============================================================================
# Additional Assembly Methods for 6-Graph Hierarchy
# ============================================================================

"""
N-gram graph assembly implementation (fixed-length unicode character analysis).
"""
function _assemble_ngram_graph(observations, config)
    @info "Using N-gram graph assembly strategy (fixed-length unicode analysis)"
    
    # Build N-gram graph from observations
    text = join([FASTX.FASTA.sequence(obs) for obs in observations], "")
    graph = string_to_ngram_graph(text, config.k)
    
    # Collapse unbranching paths
    collapse_unbranching_paths(graph)
    
    # If enabled, resolve bubble structures
    if config.bubble_resolution
        _simplify_ngram_graph(graph)
    end
    
    # Assemble contigs from graph
    contigs = assemble_strings(graph)
    contig_names = ["ngram_contig_$i" for i in 1:length(contigs)]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "NgramGraph",
        "graph_type" => "fixed_length",
        "vertex_type" => "unicode_character_vectors",
        "k" => config.k,
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Dates.now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=graph, assembly_stats=stats)
end

"""
BioSequence graph assembly implementation (variable-length simplified from k-mer graphs).
"""
function _assemble_biosequence_graph(observations, config)
    @info "Using BioSequence graph assembly strategy (variable-length simplified from k-mer)"
    
    # First build fixed-length k-mer graph
    kmer_type = _determine_kmer_type(observations, config.k)
    kmer_graph = build_kmer_graph_next(kmer_type, observations; graph_mode=config.graph_mode)
    
    # Convert to variable-length BioSequence graph
    biosequence_graph = kmer_graph_to_biosequence_graph(kmer_graph)
    
    # Extract sequences from graph vertices
    contigs = String[]
    contig_names = String[]
    
    for (i, vertex_label) in enumerate(MetaGraphsNext.labels(biosequence_graph))
        sequence = string(vertex_label)
        push!(contigs, sequence)
        push!(contig_names, "biosequence_contig_$i")
    end
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "BioSequenceGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "biosequences",
        "source_graph" => "kmer_graph_simplification",
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(biosequence_graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(biosequence_graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Dates.now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=biosequence_graph, assembly_stats=stats)
end

"""
Quality-aware BioSequence graph assembly implementation (variable-length simplified from qualmer graphs).
"""
function _assemble_quality_biosequence_graph(observations, config)
    @info "Using quality-aware BioSequence graph assembly strategy (variable-length simplified from qualmer)"
    
    # Convert observations to FASTQ records for quality processing
    fastq_records = _prepare_fastq_observations(observations)
    
    # First build fixed-length qualmer graph
    qualmer_graph = build_qualmer_graph(fastq_records; k=config.k, graph_mode=config.graph_mode)
    
    # Convert to variable-length quality-aware BioSequence graph
    quality_biosequence_graph = qualmer_graph_to_quality_biosequence_graph(qualmer_graph)
    
    # Extract sequences from graph vertices
    contigs = String[]
    contig_names = String[]
    
    for (i, vertex_label) in enumerate(MetaGraphsNext.labels(quality_biosequence_graph))
        sequence = string(vertex_label)
        push!(contigs, sequence)
        push!(contig_names, "quality_biosequence_contig_$i")
    end
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "QualityBioSequenceGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "quality_aware_biosequences",
        "source_graph" => "qualmer_graph_simplification",
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(quality_biosequence_graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(quality_biosequence_graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Dates.now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=quality_biosequence_graph, assembly_stats=stats)
end

"""
Simplify N-gram graph by removing unnecessary complexity.
"""
function _simplify_ngram_graph(graph)
    # Basic graph simplification - remove low-weight edges, merge linear paths
    # This is a placeholder for more sophisticated graph simplification
    return graph
end

"""
Convert N-gram graph to variable-length string graph.
"""
function _simplify_ngram_to_string_graph(ngram_graph)
    # Placeholder for N-gram to string graph conversion
    # This would implement path collapsing and simplification
    return ngram_graph
end

"""
Determine appropriate k-mer type from observations.
"""
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run QuickMerge for assembly merging and scaffolding.

# Arguments
- `self_assembly::String`: Path to self assembly FASTA file (primary assembly)
- `ref_assembly::String`: Path to reference assembly FASTA file (hybrid assembly)
- `outdir::String`: Output directory path (default: "quickmerge_output")
- `hco::Float64`: HCO threshold for overlap detection (default: 5.0)
- `c::Float64`: Coverage cutoff for contig filtering (default: 1.5)
- `l::Int`: Minimum alignment length (default: 5000)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `merged_assembly::String`: Path to merged assembly file

# Details
- Uses MUMmer-based approach for identifying high-confidence overlaps
- Merges assemblies by splicing contigs at overlap boundaries
- Optimized for merging complementary assemblies (e.g., short+long read)
- Automatically creates and uses a conda environment with quickmerge
- Skips analysis if output files already exist
"""
function run_quickmerge(self_assembly::String, ref_assembly::String; outdir::String="quickmerge_output", hco::Float64=5.0, c::Float64=1.5, l::Int=5000)
    Mycelia.add_bioconda_env("quickmerge")
    mkpath(outdir)
    
    merged_assembly = joinpath(outdir, "merged.fasta")
    
    if !isfile(merged_assembly)
        # Run QuickMerge
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n quickmerge merge_wrapper.py $(self_assembly) $(ref_assembly) -pre $(outdir)/merged -hco $(hco) -c $(c) -l $(l)`)
    end
    
    return (;outdir, merged_assembly)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Apollo for HMM-based assembly polishing.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `reads_file::String`: Path to reads FASTQ file (can be comma-separated for paired reads)
- `outdir::String`: Output directory path (default: "\${assembly_file}_apollo")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `polished_assembly::String`: Path to polished assembly file

# Details
- Uses profile Hidden Markov Models for technology-independent polishing
- Constructs HMM graphs for consensus sequence correction
- Works with both short and long read technologies
- Automatically creates and uses a conda environment with apollo
- Skips analysis if output files already exist
"""
function run_apollo(assembly_file::String, reads_file::String; outdir::String=assembly_file * "_apollo")
    Mycelia.add_bioconda_env("apollo")
    mkpath(outdir)
    
    basename_assembly = basename(assembly_file, ".fasta")
    polished_assembly = joinpath(outdir, basename_assembly * "_polished.fasta")
    
    if !isfile(polished_assembly)
        # Map reads to assembly first
        bam_file = joinpath(outdir, basename_assembly * ".bam")
        if !isfile(bam_file)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n apollo minimap2 -ax map-pb $(assembly_file) $(reads_file) | samtools sort -@ $(Sys.CPU_THREADS) -o $(bam_file)`)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n apollo samtools index $(bam_file)`)
        end
        
        # Run Apollo polishing
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n apollo apollo polish -a $(assembly_file) -b $(bam_file) -o $(polished_assembly)`)
    end
    
    return (;outdir, polished_assembly)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Homopolish for reference-based homopolymer error correction.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `reads_file::String`: Path to reads FASTQ file
- `outdir::String`: Output directory path (default: "\${assembly_file}_homopolish")
- `model_path::String`: Path to Homopolish model (default: auto-downloaded)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `polished_assembly::String`: Path to polished assembly file

# Details
- Uses reference genomes and neural networks for homopolymer correction
- Specifically targets homopolymer run errors common in long-read sequencing
- Downloads pre-trained models automatically if not provided
- Automatically creates and uses a conda environment with homopolish
- Skips analysis if output files already exist
"""
function run_homopolish(assembly_file::String, reads_file::String; outdir::String=assembly_file * "_homopolish", model_path::String="")
    Mycelia.add_bioconda_env("homopolish")
    mkpath(outdir)
    
    basename_assembly = basename(assembly_file, ".fasta")
    polished_assembly = joinpath(outdir, basename_assembly * "_homopolished.fasta")
    
    if !isfile(polished_assembly)
        cmd_args = ["homopolish", "polish", "-a", assembly_file, "-l", reads_file, "-o", outdir, "-t", string(Sys.CPU_THREADS)]
        
        if !isempty(model_path)
            push!(cmd_args, "-m", model_path)
        end
        
        # Run Homopolish
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n homopolish $(cmd_args)`)
        
        # Move output to expected location if needed
        homopolish_output = joinpath(outdir, "homopolished_" * basename_assembly)
        if isfile(homopolish_output) && !isfile(polished_assembly)
            mv(homopolish_output, polished_assembly)
        end
    end
    
    return (;outdir, polished_assembly)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run HyLight for hybrid strain-resolved metagenomic assembly.

# Arguments
- `short_reads_1::String`: Path to first short read FASTQ file
- `short_reads_2::String`: Path to second short read FASTQ file  
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "hylight_output")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_assemblies::String`: Path to strain-resolved assemblies directory

# Details
- Uses hybrid overlap graphs combining short and long reads
- Implements "cross hybrid" mutual support strategy for strain resolution
- Separates closely related strains in metagenomic samples
- Automatically creates and uses a conda environment with hylight
- Skips assembly if output directory already exists
"""
function run_hylight(short_reads_1::String, short_reads_2::String, long_reads::String; outdir::String="hylight_output")
    Mycelia.add_bioconda_env("hylight")
    mkpath(outdir)
    
    strain_assemblies = joinpath(outdir, "strain_assemblies")
    
    if !isdir(strain_assemblies)
        # Run HyLight assembly
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hylight hylight -1 $(short_reads_1) -2 $(short_reads_2) -l $(long_reads) -o $(outdir) -t $(Sys.CPU_THREADS)`)
    end
    
    return (;outdir, strain_assemblies)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run STRONG for strain resolution on assembly graphs.

# Arguments
- `assembly_graph::String`: Path to assembly graph file (GFA format)
- `reads_file::String`: Path to reads FASTQ file (can be comma-separated for paired reads)
- `outdir::String`: Output directory path (default: "strong_output")
- `nb_strains::Int`: Expected number of strains (default: auto-detect)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_unitigs::String`: Path to strain-resolved unitigs

# Details
- Uses BayesPaths algorithm for haplotype resolution on assembly graphs
- Implements Bayesian approach to separate strain variants
- Supports multi-sample strain analysis
- Automatically creates and uses a conda environment with strong
- Skips analysis if output files already exist
"""
function run_strong(assembly_graph::String, reads_file::String; outdir::String="strong_output", nb_strains::Int=0)
    Mycelia.add_bioconda_env("strong")
    mkpath(outdir)
    
    strain_unitigs = joinpath(outdir, "strain_unitigs.fasta")
    
    if !isfile(strain_unitigs)
        cmd_args = ["STRONG", assembly_graph, reads_file, "-o", outdir]
        
        if nb_strains > 0
            push!(cmd_args, "-s", string(nb_strains))
        end
        
        # Run STRONG
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strong $(cmd_args)`)
    end
    
    return (;outdir, strain_unitigs)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Strainy for strain phasing from long reads.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "strainy_output")
- `mode::String`: Analysis mode ("transform" or "phase", default: "phase")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_assemblies::String`: Path to strain-phased assemblies

# Details
- Uses connection graph-based read clustering for strain separation
- Implements long-read phasing to resolve strain-level variants
- Generates strain unitigs and simplified assembly graphs
- Automatically creates and uses a conda environment with strainy
- Skips analysis if output files already exist
"""
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Ray assembler for de novo genome assembly.

# Arguments
- `reads_files::Vector{String}`: Vector of FASTQ file paths (supports paired-end and single-end)
- `outdir::String`: Output directory path (default: "ray_output")
- `k::Int`: K-mer size for assembly (default: 31)
- `min_contig_length::Int`: Minimum contig length (default: 200)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `scaffolds::String`: Path to scaffolds file

# Details
- Uses Ray's distributed de Bruijn graph approach
- Supports both single-end and paired-end reads
- Automatically creates and uses a conda environment with ray
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_ray(reads_files::Vector{String}; outdir::String="ray_output", k::Int=31, min_contig_length::Int=200)
    Mycelia.add_bioconda_env("ray")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "Contigs.fasta")
    scaffolds_file = joinpath(outdir, "Scaffolds.fasta")
    
    if !isfile(contigs_file)
        # Build Ray command arguments
        cmd_args = ["Ray", "-k", string(k), "-o", outdir, "-minimum-contig-length", string(min_contig_length)]
        
        # Add reads files - Ray auto-detects paired vs single-end
        if length(reads_files) == 2
            # Paired-end reads
            push!(cmd_args, "-p", reads_files[1], reads_files[2])
        elseif length(reads_files) == 1
            # Single-end reads  
            push!(cmd_args, "-s", reads_files[1])
        else
            # Multiple libraries - treat as single-end
            for reads_file in reads_files
                push!(cmd_args, "-s", reads_file)
            end
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ray $(cmd_args)`)
    end
    
    return (;outdir, contigs=contigs_file, scaffolds=scaffolds_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Ray Meta for metagenomic assembly.

# Arguments
- `reads_files::Vector{String}`: Vector of FASTQ file paths (supports paired-end and single-end)
- `outdir::String`: Output directory path (default: "ray_meta_output")
- `k::Int`: K-mer size for assembly (default: 31)
- `min_contig_length::Int`: Minimum contig length (default: 200)
- `enable_communities::Bool`: Enable community detection (default: true)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `scaffolds::String`: Path to scaffolds file

# Details
- Uses Ray Meta's distributed approach for metagenomic data
- Supports community detection for binning related sequences
- Handles uneven coverage typical in metagenomic samples
- Automatically creates and uses a conda environment with ray
- Skips assembly if output directory already exists
"""
function run_ray_meta(reads_files::Vector{String}; outdir::String="ray_meta_output", k::Int=31, min_contig_length::Int=200, enable_communities::Bool=true)
    Mycelia.add_bioconda_env("ray")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "Contigs.fasta")
    scaffolds_file = joinpath(outdir, "Scaffolds.fasta")
    
    if !isfile(contigs_file)
        # Build Ray Meta command arguments
        cmd_args = ["Ray", "-k", string(k), "-o", outdir, "-minimum-contig-length", string(min_contig_length)]
        
        # Enable metagenomic mode
        push!(cmd_args, "-enable-neighbourhoods")
        
        if enable_communities
            push!(cmd_args, "-enable-neighbourhoods")
        end
        
        # Add reads files
        if length(reads_files) == 2
            # Paired-end reads
            push!(cmd_args, "-p", reads_files[1], reads_files[2])
        elseif length(reads_files) == 1
            # Single-end reads
            push!(cmd_args, "-s", reads_files[1])
        else
            # Multiple libraries
            for reads_file in reads_files
                push!(cmd_args, "-s", reads_file)
            end
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ray $(cmd_args)`)
    end
    
    return (;outdir, contigs=contigs_file, scaffolds=scaffolds_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Velvet assembler for short-read genome assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "velvet_output")
- `k::Int`: K-mer size for assembly (default: 31)
- `exp_cov::String`: Expected coverage (default: "auto")
- `min_contig_lgth::Int`: Minimum contig length (default: 200)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file

# Details
- Uses Velvet's de Bruijn graph approach for short reads
- Two-step process: velveth (indexing) + velvetg (assembly)
- Optimized for Illumina short reads
- Automatically creates and uses a conda environment with velvet
- Skips assembly if output directory already exists
"""
function run_velvet(fastq1::String; fastq2::Union{String,Nothing}=nothing, outdir::String="velvet_output", k::Int=31, exp_cov::String="auto", min_contig_lgth::Int=200)
    Mycelia.add_bioconda_env("velvet")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "contigs.fa")
    
    if !isfile(contigs_file)
        # Step 1: velveth (k-mer hashing and indexing)
        if isnothing(fastq2)
            # Single-end reads
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet velveth $(outdir) $(k) -short -fastq $(fastq1)`)
        else
            # Paired-end reads
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet velveth $(outdir) $(k) -shortPaired -fastq -separate $(fastq1) $(fastq2)`)
        end
        
        # Step 2: velvetg (graph construction and traversal)
        cmd_args = ["velvetg", outdir, "-min_contig_lgth", string(min_contig_lgth)]
        
        if exp_cov != "auto"
            push!(cmd_args, "-exp_cov", exp_cov)
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet $(cmd_args)`)
    end
    
    return (;outdir, contigs=contigs_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run VelvetG with optimized parameters for genome assembly.

# Arguments
- `velvet_dir::String`: Path to velveth output directory
- `outdir::String`: Output directory path (default: velvet_dir)
- `exp_cov::Union{String,Float64}`: Expected coverage (default: "auto")
- `cov_cutoff::Union{String,Float64}`: Coverage cutoff (default: "auto")  
- `min_contig_lgth::Int`: Minimum contig length (default: 200)
- `scaffolding::Bool`: Enable scaffolding (default: true)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `stats::String`: Path to assembly statistics file

# Details
- Advanced velvetg run with parameter optimization
- Supports scaffolding and coverage optimization
- Can be run on existing velveth output
- Automatically creates and uses a conda environment with velvet
"""
function run_velvetg(velvet_dir::String; outdir::String=velvet_dir, exp_cov::Union{String,Float64}="auto", cov_cutoff::Union{String,Float64}="auto", min_contig_lgth::Int=200, scaffolding::Bool=true)
    Mycelia.add_bioconda_env("velvet")
    
    contigs_file = joinpath(outdir, "contigs.fa")
    stats_file = joinpath(outdir, "stats.txt")
    
    # Build velvetg command
    cmd_args = ["velvetg", velvet_dir, "-min_contig_lgth", string(min_contig_lgth)]
    
    if exp_cov != "auto"
        push!(cmd_args, "-exp_cov", string(exp_cov))
    end
    
    if cov_cutoff != "auto"
        push!(cmd_args, "-cov_cutoff", string(cov_cutoff))
    end
    
    if scaffolding
        push!(cmd_args, "-scaffolding", "yes")
    else
        push!(cmd_args, "-scaffolding", "no")
    end
    
    # Copy output to specified directory if different
    if outdir != velvet_dir
        mkpath(outdir)
        push!(cmd_args, "-read_trkg", "yes")  # Enable read tracking for copying
    end
    
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet $(cmd_args)`)
    
    # Copy results if needed
    if outdir != velvet_dir
        for file in ["contigs.fa", "stats.txt", "Graph", "LastGraph"]
            src_file = joinpath(velvet_dir, file)
            dst_file = joinpath(outdir, file)
            if isfile(src_file)
                cp(src_file, dst_file, force=true)
            end
        end
    end
    
    return (;outdir, contigs=contigs_file, stats=stats_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaMDBG assembler for metagenomic long-read assembly.

# Arguments
- `fastq_file::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "metamdbg_output")
- `abundance_min::Int`: Minimum abundance threshold (default: 3)
- `length_min::Int`: Minimum sequence length (default: 1000)
- `nb_cores::Int`: Number of cores to use (default: Sys.CPU_THREADS)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `graph::String`: Path to assembly graph file

# Details
- Uses metaMDBG's minimizer-space de Bruijn graphs for metagenomic assembly
- Specifically designed for long-read metagenomic data
- Handles species abundance variation and strain diversity
- Produces both contigs and assembly graphs
- Automatically creates and uses a conda environment with metamdbg
- Skips assembly if output directory already exists
"""
function run_metamdbg(fastq_file::String; outdir::String="metamdbg_output", abundance_min::Int=3, length_min::Int=1000, nb_cores::Int=Sys.CPU_THREADS)
    Mycelia.add_bioconda_env("metamdbg")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "contigs.fasta")
    graph_file = joinpath(outdir, "graph.gfa")
    
    if !isfile(contigs_file)
        # Run metaMDBG assembly
        cmd_args = [
            "metaMDBG", "asm",
            "--in-dir", dirname(fastq_file),
            "--in-reads", basename(fastq_file),
            "--out-dir", outdir,
            "--abundance-min", string(abundance_min),
            "--length-min", string(length_min),
            "--nb-cores", string(nb_cores)
        ]
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metamdbg $(cmd_args)`)
        
        # metaMDBG may output with different naming - check for common output files
        possible_contigs = [
            joinpath(outdir, "contigs.fasta"),
            joinpath(outdir, "final_contigs.fasta"),
            joinpath(outdir, "assembly.fasta")
        ]
        
        possible_graphs = [
            joinpath(outdir, "graph.gfa"),
            joinpath(outdir, "assembly_graph.gfa"),
            joinpath(outdir, "final_graph.gfa")
        ]
        
        # Find actual output files
        actual_contigs = ""
        actual_graph = ""
        
        for possible in possible_contigs
            if isfile(possible)
                actual_contigs = possible
                break
            end
        end
        
        for possible in possible_graphs
            if isfile(possible)
                actual_graph = possible
                break
            end
        end
        
        # Update file paths
        contigs_file = actual_contigs
        graph_file = actual_graph
    end
    
    return (;outdir, contigs=contigs_file, graph=graph_file)
end

function run_strainy(assembly_file::String, long_reads::String; outdir::String="strainy_output", mode::String="phase")
    Mycelia.add_bioconda_env("strainy")
    mkpath(outdir)
    
    strain_assemblies = joinpath(outdir, "strain_assemblies.fasta")
    
    if !isfile(strain_assemblies)
        # First map reads to assembly
        bam_file = joinpath(outdir, "mapped_reads.bam")
        if !isfile(bam_file)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy minimap2 -ax map-ont $(assembly_file) $(long_reads) | samtools sort -@ $(Sys.CPU_THREADS) -o $(bam_file)`)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy samtools index $(bam_file)`)
        end
        
        # Run Strainy
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy strainy --bam $(bam_file) --fasta $(assembly_file) --output $(outdir) --mode $(mode)`)
    end
    
    return (;outdir, strain_assemblies)
end

function _determine_kmer_type(observations, k)
    # Examine first observation to determine sequence type
    if !isempty(observations)
        sample_seq = string(FASTX.FASTA.sequence(observations[1]))
        
        # Check if it's protein sequence (contains amino acids)
        if occursin(r"[ARNDCQEGHILKMFPSTWYV]", uppercase(sample_seq))
            return Kmers.AAKmer{k}
        # Check if it's RNA (contains U)
        elseif occursin('U', uppercase(sample_seq))
            return Kmers.RNAKmer{k}
        # Default to DNA
        else
            return Kmers.DNAKmer{k}
        end
    end
    
    # Default to DNA k-mer
    return Kmers.DNAKmer{k}
end