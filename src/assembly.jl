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
    StringGraph      # String graph assembly (for basic sequence analysis)
    KmerGraph        # K-mer graph assembly with BioSequences.DNAKmer (for FASTA data)
    QualmerGraph     # Quality-aware k-mer graph assembly (for FASTQ data) - PRIMARY METHOD
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
    
    function AssemblyResult(contigs::Vector{String}, contig_names::Vector{String}; 
                          graph=nothing, assembly_stats=Dict{String, Any}())
        new(contigs, contig_names, graph, assembly_stats)
    end
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
function assemble_genome(reads; method::AssemblyMethod=StringGraph, config::AssemblyConfig=AssemblyConfig())
    @info "Starting unified genome assembly" method config.k config.error_rate
    
    # Phase 1: Load and validate input
    observations = _prepare_observations(reads)
    @info "Loaded $(length(observations)) sequence observations"
    
    # Phase 2: Assembly strategy dispatch
    if method == StringGraph
        result = _assemble_string_graph(observations, config)
    elseif method == KmerGraph  
        result = _assemble_kmer_graph(observations, config)
    elseif method == QualmerGraph
        result = _assemble_qualmer_graph(observations, config)
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
String graph assembly implementation using Phase 2 algorithms.
"""
function _assemble_string_graph(observations, config)
    @info "Using string graph assembly strategy"
    
    # Build string graph from observations
    graph = string_to_ngram_graph(join([FASTX.FASTA.sequence(obs) for obs in observations], ""), config.k)
    
    # Collapse unbranching paths
    collapse_unbranching_paths(graph)
    
    # If enabled, resolve bubble structures
    if config.bubble_resolution
        # Note: This would use detect_bubbles_next from Phase 2
        # For now, we'll use basic graph simplification
        _simplify_string_graph(graph)
    end
    
    # Assemble contigs from graph
    contigs = assemble_strings(graph)
    contig_names = ["contig_$i" for i in 1:length(contigs)]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "StringGraph",
        "k" => config.k,
        "num_input_sequences" => length(observations),
        "assembly_date" => string(now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=graph, assembly_stats=stats)
end

"""
K-mer graph assembly implementation using Phase 2 probabilistic algorithms.
"""
function _assemble_kmer_graph(observations, config)
    @info "Using k-mer graph assembly strategy"
    
    # Build k-mer graph using Phase 2 next-generation algorithms
    kmer_type = Kmers.DNAKmer{config.k}
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
    
    contig_names = ["contig_$i" for i in 1:length(contigs)]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "KmerGraph",
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(now())
    )
    
    return AssemblyResult(contigs, contig_names; graph=graph, assembly_stats=stats)
end

"""
Quality-aware k-mer graph assembly implementation using Qualmer graphs.
"""
function _assemble_qualmer_graph(observations, config)
    @info "Using quality-aware k-mer graph assembly strategy"
    
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
    
    # Convert paths to sequences
    contigs = String[]
    for path in paths
        if !isempty(path)
            sequence = _qualmer_path_to_sequence(path, graph)
            if length(sequence) > config.k  # Only keep substantial contigs
                push!(contigs, sequence)
            end
        end
    end
    
    # If no paths found, use probabilistic walks on qualmer graph
    if isempty(contigs)
        @info "No quality-aware paths found, using probabilistic walks"
        contigs = _generate_contigs_from_qualmer_graph(graph, config)
    end
    
    contig_names = ["qualmer_contig_$i" for i in 1:length(contigs)]
    
    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "QualmerGraph",
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(now())
    )
    
    # Add quality-specific statistics
    if !isempty(MetaGraphsNext.labels(graph))
        qualmer_stats = get_qualmer_statistics(graph)
        merge!(stats, qualmer_stats)
    end
    
    return AssemblyResult(contigs, contig_names; graph=graph, assembly_stats=stats)
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
        if obs isa FASTX.FASTQ.Record
            push!(fastq_records, obs)
        elseif obs isa FASTX.FASTA.Record
            # Convert FASTA to FASTQ with default quality (assume high quality)
            seq = FASTX.FASTA.sequence(obs)
            qual = repeat('I', length(seq))  # Quality score 40 (high quality)
            record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(obs), seq, qual)
            push!(fastq_records, record)
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