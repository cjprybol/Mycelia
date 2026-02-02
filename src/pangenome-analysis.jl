# Pangenome Analysis Module
# Functions for comparative genomics using existing k-mer infrastructure
# Includes integration with PGGB, Cactus, and vg toolkit

"""
    PangenomeAnalysisResult

Results of k-mer based pangenome analysis.
"""
struct PangenomeAnalysisResult
    genome_names::Vector{String}
    kmer_counts_by_genome::Dict{String, Dict}
    shared_kmers::Vector
    core_kmers::Vector  # Present in ALL genomes
    accessory_kmers::Vector  # Present in SOME genomes
    unique_kmers_by_genome::Dict{String, Vector}
    presence_absence_matrix::BitMatrix
    distance_matrix::Matrix{Float64}
    similarity_stats::NamedTuple
end

"""
    analyze_pangenome_kmers(genome_files::Vector{String}; kmer_type=Kmers.DNAKmer{21}, distance_metric=:jaccard)

Perform comprehensive k-mer based pangenome analysis using existing Mycelia infrastructure.

Leverages existing `count_canonical_kmers` and distance metric functions to analyze
genomic content across multiple genomes, identifying core, accessory, and unique regions.

# Arguments
- `genome_files`: Vector of FASTA file paths containing genome sequences
- `kmer_type`: K-mer type from Kmers.jl (default: `Kmers.DNAKmer{21}`)
- `distance_metric`: Distance metric (`:jaccard`, `:bray_curtis`, `:cosine`, `:js_divergence`)

# Returns
- `PangenomeAnalysisResult` with comprehensive pangenome statistics

# Example
```julia
genome_files = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
result = Mycelia.analyze_pangenome_kmers(genome_files, kmer_type=Kmers.DNAKmer{31})
println("Core k-mers: \$(length(result.core_kmers))")
println("Total pangenome size: \$(size(result.presence_absence_matrix, 1)) k-mers")
```
"""
function analyze_pangenome_kmers(genome_files::Vector{String};
        kmer_type = Kmers.DNAKmer{21}, distance_metric = :jaccard)
    if isempty(genome_files)
        error("No genome files provided")
    end

    # Validate files exist
    for file in genome_files
        if !isfile(file)
            error("Genome file does not exist: $(file)")
        end
    end

    genome_names = [basename(file) for file in genome_files]
    n_genomes = length(genome_files)

    # Count k-mers for each genome using existing infrastructure
    println("Counting k-mers for $(n_genomes) genomes...")
    kmer_counts_by_genome = Dict{String, Dict}()

    for (i, file) in enumerate(genome_files)
        println("  Processing $(genome_names[i])...")
        # Use existing count_canonical_kmers function
        kmer_counts = count_canonical_kmers(kmer_type, file)
        kmer_counts_by_genome[genome_names[i]] = kmer_counts
    end

    # Get all unique k-mers across all genomes
    all_kmers = Set()
    for kmer_counts in values(kmer_counts_by_genome)
        union!(all_kmers, keys(kmer_counts))
    end
    all_kmers = collect(all_kmers)
    n_kmers = length(all_kmers)

    println("Found $(n_kmers) unique k-mers across all genomes")

    # Build presence/absence matrix
    presence_absence_matrix = BitMatrix(undef, n_kmers, n_genomes)

    for (i, kmer) in enumerate(all_kmers)
        for (j, genome_name) in enumerate(genome_names)
            presence_absence_matrix[i, j] = haskey(kmer_counts_by_genome[genome_name], kmer)
        end
    end

    # Classify k-mers into categories
    core_kmers = []
    accessory_kmers = []
    shared_kmers = []
    unique_kmers_by_genome = Dict{String, Vector}()

    for genome_name in genome_names
        unique_kmers_by_genome[genome_name] = []
    end

    for (i, kmer) in enumerate(all_kmers)
        presence_count = sum(presence_absence_matrix[i, :])

        if presence_count == n_genomes
            # Present in all genomes - core
            push!(core_kmers, kmer)
            push!(shared_kmers, kmer)
        elseif presence_count == 1
            # Present in only one genome - unique
            genome_idx = findfirst(presence_absence_matrix[i, :])
            genome_name = genome_names[genome_idx]
            push!(unique_kmers_by_genome[genome_name], kmer)
        elseif presence_count > 1
            # Present in some genomes - accessory
            push!(accessory_kmers, kmer)
            push!(shared_kmers, kmer)
        end
    end

    # Calculate distance matrix using existing distance functions
    println("Calculating distance matrix...")
    distance_matrix = if distance_metric == :jaccard
        jaccard_distance(presence_absence_matrix)
    elseif distance_metric == :bray_curtis
        # Convert to count matrix for Bray-Curtis
        count_matrix = zeros(Int, n_kmers, n_genomes)
        for (i, kmer) in enumerate(all_kmers)
            for (j, genome_name) in enumerate(genome_names)
                count_matrix[i, j] = get(kmer_counts_by_genome[genome_name], kmer, 0)
            end
        end
        bray_curtis_distance(count_matrix)
    else
        error("Unsupported distance metric: $(distance_metric)")
    end

    # Calculate summary statistics
    core_size = length(core_kmers)
    accessory_size = length(accessory_kmers)
    pangenome_size = n_kmers
    unique_total = sum(length(kmers) for kmers in values(unique_kmers_by_genome))

    similarity_stats = (
        n_genomes = n_genomes,
        pangenome_size = pangenome_size,
        core_size = core_size,
        accessory_size = accessory_size,
        unique_total = unique_total,
        core_percentage = (core_size / pangenome_size) * 100,
        accessory_percentage = (accessory_size / pangenome_size) * 100,
        unique_percentage = (unique_total / pangenome_size) * 100,
        mean_pairwise_distance = Statistics.mean(distance_matrix[distance_matrix .> 0])
    )

    println("Pangenome Analysis Results:")
    println("  Total k-mers: $(pangenome_size)")
    println("  Core k-mers: $(core_size) ($(round(similarity_stats.core_percentage, digits=1))%)")
    println("  Accessory k-mers: $(accessory_size) ($(round(similarity_stats.accessory_percentage, digits=1))%)")
    println("  Unique k-mers: $(unique_total) ($(round(similarity_stats.unique_percentage, digits=1))%)")
    println("  Mean pairwise distance: $(round(similarity_stats.mean_pairwise_distance, digits=3))")

    return PangenomeAnalysisResult(
        genome_names,
        kmer_counts_by_genome,
        shared_kmers,
        core_kmers,
        accessory_kmers,
        unique_kmers_by_genome,
        presence_absence_matrix,
        distance_matrix,
        similarity_stats
    )
end

"""
    compare_genome_kmer_similarity(genome1_file::String, genome2_file::String; kmer_type=Kmers.DNAKmer{21}, metric=:js_divergence)

Compare two genomes using existing k-mer distance metrics.

Leverages existing distance metric functions to compare genomic similarity
between pairs of genomes using various distance measures.

# Arguments
- `genome1_file`: Path to first genome FASTA file
- `genome2_file`: Path to second genome FASTA file  
- `kmer_type`: K-mer type from Kmers.jl (default: `Kmers.DNAKmer{21}`)
- `metric`: Distance metric (`:js_divergence`, `:cosine`, `:jaccard`)

# Returns
- Named tuple with similarity/distance metrics and k-mer statistics

# Example
```julia
similarity = Mycelia.compare_genome_kmer_similarity(
    "genome1.fasta", "genome2.fasta", 
    kmer_type=Kmers.DNAKmer{31}, 
    metric=:js_divergence
)
println("JS divergence: \$(similarity.distance)")
println("Shared k-mers: \$(similarity.shared_kmers)")
```
"""
function compare_genome_kmer_similarity(genome1_file::String, genome2_file::String;
        kmer_type = Kmers.DNAKmer{21}, metric = :js_divergence)
    # Count k-mers using existing infrastructure
    kmer_counts1 = count_canonical_kmers(kmer_type, genome1_file)
    kmer_counts2 = count_canonical_kmers(kmer_type, genome2_file)

    # Calculate distance using existing functions
    distance = if metric == :js_divergence
        kmer_counts_to_js_divergence(kmer_counts1, kmer_counts2)
    elseif metric == :cosine
        kmer_counts_to_cosine_similarity(kmer_counts1, kmer_counts2)
    elseif metric == :jaccard
        # Calculate Jaccard manually for k-mer dicts
        shared_kmers = intersect(keys(kmer_counts1), keys(kmer_counts2))
        union_kmers = union(keys(kmer_counts1), keys(kmer_counts2))
        1.0 - length(shared_kmers) / length(union_kmers)
    else
        error("Unsupported metric: $(metric)")
    end

    # Calculate additional statistics
    shared_kmers = length(intersect(keys(kmer_counts1), keys(kmer_counts2)))
    total_kmers = length(union(keys(kmer_counts1), keys(kmer_counts2)))
    genome1_unique = length(setdiff(keys(kmer_counts1), keys(kmer_counts2)))
    genome2_unique = length(setdiff(keys(kmer_counts2), keys(kmer_counts1)))

    return (
        distance = distance,
        metric = metric,
        shared_kmers = shared_kmers,
        total_kmers = total_kmers,
        genome1_unique = genome1_unique,
        genome2_unique = genome2_unique,
        jaccard_similarity = shared_kmers / total_kmers
    )
end

"""
    build_genome_distance_matrix(genome_files::Vector{String}; kmer_type=Kmers.DNAKmer{21}, metric=:js_divergence)

Build a distance matrix between all genome pairs using existing distance metrics.

Creates a comprehensive pairwise distance matrix using established k-mer
distance functions, suitable for phylogenetic analysis and clustering.

# Arguments
- `genome_files`: Vector of genome FASTA file paths
- `kmer_type`: K-mer type from Kmers.jl (default: `Kmers.DNAKmer{21}`)
- `metric`: Distance metric (`:js_divergence`, `:cosine`, `:jaccard`)

# Returns
- Named tuple with distance matrix and genome names

# Example
```julia
genomes = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
result = Mycelia.build_genome_distance_matrix(genomes, kmer_type=Kmers.DNAKmer{31})
println("Distance matrix: \$(result.distance_matrix)")
```
"""
function build_genome_distance_matrix(genome_files::Vector{String};
        kmer_type = Kmers.DNAKmer{21}, metric = :js_divergence)
    n_genomes = length(genome_files)
    distance_matrix = zeros(Float64, n_genomes, n_genomes)
    genome_names = [basename(file) for file in genome_files]

    println("Building $(n_genomes)x$(n_genomes) distance matrix using $(metric)...")

    # Calculate pairwise distances
    for i in 1:n_genomes
        for j in (i + 1):n_genomes
            println("  Comparing $(genome_names[i]) vs $(genome_names[j])...")
            similarity = compare_genome_kmer_similarity(
                genome_files[i], genome_files[j];
                kmer_type = kmer_type, metric = metric
            )
            distance_matrix[i, j] = similarity.distance
            distance_matrix[j, i] = similarity.distance  # Symmetric
        end
    end

    return (
        distance_matrix = distance_matrix,
        genome_names = genome_names,
        metric = metric
    )
end

# PGGB Integration Functions

"""
    construct_pangenome_pggb(genome_files::Vector{String}, output_dir::String; 
                            threads::Int=2, segment_length::Int=5000, 
                            block_length::Int=3*segment_length, 
                            mash_kmer::Int=16, min_match_length::Int=19,
                            transclose_batch::Int=10000000, 
                            additional_args::Vector{String}=String[])

Construct a pangenome using PGGB (PanGenome Graph Builder).

Uses the PGGB tool to build pangenome graphs from multiple genome assemblies,
following the workflows established in the Mycelia-Dev benchmarking notebooks.

# Arguments
- `genome_files`: Vector of FASTA file paths to include in pangenome
- `output_dir`: Directory for PGGB output files
- `threads`: Number of threads for parallel processing (default: 2)
- `segment_length`: Segment length for mapping (default: 5000)
- `block_length`: Block length for mapping (default: 3*segment_length)
- `mash_kmer`: Kmer size for mash sketching (default: 16)
- `min_match_length`: Minimum match length (default: 19)
- `transclose_batch`: Batch size for transitive closure (default: 10000000)
- `additional_args`: Additional command line arguments

# Returns
- Path to the main GFA output file

# Example
```julia
genomes = ["reference.fasta", "assembly1.fasta", "assembly2.fasta"]
gfa_file = Mycelia.construct_pangenome_pggb(genomes, "pangenome_output")
```
"""
function construct_pangenome_pggb(genome_files::Vector{String}, output_dir::String;
        threads::Int = 2, segment_length::Int = 5000,
        block_length::Int = 3*segment_length,
        mash_kmer::Int = 16, min_match_length::Int = 19,
        transclose_batch::Int = 10000000,
        additional_args::Vector{String} = String[])

    # Validate input files
    for file in genome_files
        if !isfile(file)
            error("Genome file does not exist: $(file)")
        end
    end

    # Create joint FASTA file for PGGB input
    joint_fasta = joinpath(output_dir, "joint_genomes.fasta")
    mkpath(dirname(joint_fasta))

    # Concatenate genome files (don't merge to preserve identifiers)
    concatenate_files(files = genome_files, file = joint_fasta)

    # Index the joint FASTA if not already indexed
    if !isfile(joint_fasta * ".fai")
        samtools_index_fasta(fasta = joint_fasta)
    end

    # Ensure PGGB conda environment exists
    add_bioconda_env("pggb")

    # Build PGGB command
    pggb_args = [
        "-i", joint_fasta,
        "-o", output_dir,
        "-t", string(threads),
        "-s", string(segment_length),
        "-l", string(block_length),
        "-k", string(mash_kmer),
        "-G", string(min_match_length),
        "-B", string(transclose_batch),
        "-n", string(length(genome_files))
    ]

    # Add any additional arguments
    append!(pggb_args, additional_args)

    # Run PGGB
    println("Running PGGB pangenome construction with $(length(genome_files)) genomes...")
    cmd = `$(CONDA_RUNNER) run --live-stream -n pggb pggb $(pggb_args)`

    try
        run(cmd)
    catch e
        error("PGGB failed: $(e)")
    end

    # Find and return the main GFA output file
    gfa_files = filter(x -> endswith(x, ".gfa"), readdir(output_dir, join = true))
    if isempty(gfa_files)
        error("No GFA output file found in $(output_dir)")
    end

    main_gfa = first(sort(gfa_files, by = filesize, rev = true))  # Return largest GFA file
    println("PGGB completed. Main GFA file: $(main_gfa)")

    return main_gfa
end

"""
    call_variants_from_pggb_graph(gfa_file::String, reference_prefix::String; 
                                 threads::Int=2, ploidy::Int=1, output_file::String="")

Call variants from a PGGB-generated pangenome graph using vg deconstruct.

Uses the vg toolkit to extract variants from pangenome graphs, following
the methodology established in Mycelia-Dev benchmarking workflows.

# Arguments
- `gfa_file`: Path to GFA pangenome graph file from PGGB
- `reference_prefix`: Prefix for reference paths in the graph
- `threads`: Number of threads (default: 2)
- `ploidy`: Ploidy level (default: 1)
- `output_file`: Output VCF file path (default: gfa_file + ".vcf")

# Returns
- Path to output VCF file

# Example
```julia
gfa_file = "pangenome.gfa"
vcf_file = Mycelia.call_variants_from_pggb_graph(gfa_file, "reference")
```
"""
function call_variants_from_pggb_graph(gfa_file::String, reference_prefix::String;
        threads::Int = 2, ploidy::Int = 1, output_file::String = "")
    if !isfile(gfa_file)
        error("GFA file does not exist: $(gfa_file)")
    end

    # Set default output file name
    if isempty(output_file)
        output_file = gfa_file * ".vcf"
    end

    # Ensure vg conda environment exists
    add_bioconda_env("vg")

    # Build vg deconstruct command
    vg_args = [
        "deconstruct",
        "--path-prefix", reference_prefix,
        "--ploidy", string(ploidy),
        "--path-traversals",
        "--all-snarls",
        "--threads", string(threads),
        gfa_file
    ]

    println("Calling variants from pangenome graph: $(gfa_file)")
    cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg $(vg_args)`

    try
        run(pipeline(cmd, output_file))
    catch e
        error("vg deconstruct failed: $(e)")
    end

    println("Variant calling completed. VCF file: $(output_file)")

    return output_file
end

# Cactus Integration Functions

"""
    construct_pangenome_cactus(genome_files::Vector{String}, genome_names::Vector{String}, 
                              output_dir::String, reference_name::String;
                              max_cores::Int=8, max_memory_gb::Int=32,
                              output_formats::Vector{String}=["gbz", "gfa", "vcf", "odgi"])

Construct a pangenome using Cactus alignment-based approach.

Uses the Cactus pangenome pipeline via containerized execution to build
pangenome graphs from multiple genome assemblies using progressive alignment.

# Arguments
- `genome_files`: Vector of FASTA file paths
- `genome_names`: Vector of sample names corresponding to genome files
- `output_dir`: Directory for Cactus output
- `reference_name`: Name of the reference sample for pangenome construction
- `max_cores`: Maximum number of CPU cores (default: 8)
- `max_memory_gb`: Maximum memory in GB (default: 32)
- `output_formats`: Output formats to generate (default: ["gbz", "gfa", "vcf", "odgi"])

# Returns
- Dictionary with paths to output files for each format

# Example
```julia
genomes = ["ref.fasta", "asm1.fasta", "asm2.fasta"]
names = ["REFERENCE", "SAMPLE1", "SAMPLE2"]
outputs = Mycelia.construct_pangenome_cactus(genomes, names, "cactus_out", "REFERENCE")
```
"""
function construct_pangenome_cactus(
        genome_files::Vector{String}, genome_names::Vector{String},
        output_dir::String, reference_name::String;
        max_cores::Int = 8, max_memory_gb::Int = 32,
        output_formats::Vector{String} = ["gbz", "gfa", "vcf", "odgi"])

    # Validate inputs
    if length(genome_files) != length(genome_names)
        error("Number of genome files must match number of genome names")
    end

    for file in genome_files
        if !isfile(file)
            error("Genome file does not exist: $(file)")
        end
    end

    if !(reference_name in genome_names)
        error("Reference name '$(reference_name)' not found in genome names")
    end

    # Create output directory
    mkpath(output_dir)

    # Create Cactus configuration file
    config_table = DataFrames.DataFrame(
        samples = genome_names,
        file_paths = genome_files
    )

    config_file = joinpath(output_dir, "cactus_config.txt")
    uCSV.write(config_file, data = collect(DataFrames.eachcol(config_table)),
        header = missing, delim = '\t')

    # Set up job store and output names
    jobstore = joinpath(output_dir, "cactus-job-store")
    output_name = "pangenome"

    # Prepare output format flags
    format_flags = String[]
    for format in output_formats
        if format in ["gbz", "gfa", "vcf", "odgi"]
            push!(format_flags, "--$(format)")
        else
            @warn "Unknown output format: $(format)"
        end
    end

    # Build Cactus command using podman-hpc
    cactus_args = [
        "run", "-it",
        "-v", "$(output_dir):/app",
        "-w", "/app",
        "quay.io/comparative-genomics-toolkit/cactus:v2.8.1",
        "cactus-pangenome",
        "./$(basename(jobstore))",
        "./$(basename(config_file))",
        "--maxCores", string(max_cores),
        "--maxMemory", "$(max_memory_gb)Gb",
        "--outDir", ".",
        "--outName", output_name,
        "--reference", reference_name
    ]

    # Add format flags
    append!(cactus_args, format_flags)

    println("Running Cactus pangenome construction...")
    println("  Genomes: $(length(genome_files))")
    println("  Reference: $(reference_name)")
    println("  Output formats: $(join(output_formats, ", "))")

    cmd = `podman-hpc $(cactus_args)`

    # Set up logging
    log_file = joinpath(output_dir, "cactus.log")

    try
        run(pipeline(cmd, stdout = log_file, stderr = log_file))
    catch e
        @warn "Cactus may have failed. Check log file: $(log_file)"
        error("Cactus pangenome construction failed: $(e)")
    end

    # Find output files
    output_files = Dict{String, String}()
    for format in output_formats
        pattern = "$(output_name).$(format)"
        files = filter(x -> occursin(pattern, x), readdir(output_dir, join = true))
        if !isempty(files)
            output_files[format] = first(files)
        else
            @warn "Output file for format $(format) not found"
        end
    end

    println("Cactus pangenome construction completed.")
    println("Output files: $(length(output_files))")

    return output_files
end

# vg Toolkit Integration Functions

"""
    convert_gfa_to_vg_format(gfa_file::String; output_file::String="")

Convert a GFA pangenome graph to vg native format.

Converts GFA format pangenome graphs to vg's optimized binary format
for improved performance in downstream analysis.

# Arguments
- `gfa_file`: Input GFA file path
- `output_file`: Output vg file path (default: gfa_file with .vg extension)

# Returns
- Path to output vg file

# Example
```julia
vg_file = Mycelia.convert_gfa_to_vg_format("pangenome.gfa")
```
"""
function convert_gfa_to_vg_format(gfa_file::String; output_file::String = "")
    if !isfile(gfa_file)
        error("GFA file does not exist: $(gfa_file)")
    end

    if isempty(output_file)
        output_file = replace(gfa_file, r"\.gfa$" => ".vg")
    end

    # Ensure vg conda environment exists
    add_bioconda_env("vg")

    println("Converting GFA to vg format: $(gfa_file)")
    cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg convert -g $(gfa_file) -v`

    try
        run(pipeline(cmd, output_file))
    catch e
        error("GFA to vg conversion failed: $(e)")
    end

    println("Conversion completed. vg file: $(output_file)")

    return output_file
end

"""
    index_pangenome_graph(graph_file::String; index_types::Vector{String}=["xg", "gcsa"])

Create indexes for a pangenome graph to enable efficient querying and mapping.

Builds various index types for pangenome graphs to support different
analysis workflows including read mapping and path queries.

# Arguments
- `graph_file`: Input graph file (GFA or vg format)
- `index_types`: Types of indexes to build (default: ["xg", "gcsa"])

# Returns
- Dictionary mapping index types to their file paths

# Example
```julia
indexes = Mycelia.index_pangenome_graph("pangenome.vg", index_types=["xg", "gcsa", "snarls"])
```
"""
function index_pangenome_graph(graph_file::String; index_types::Vector{String} = [
        "xg", "gcsa"])
    if !isfile(graph_file)
        error("Graph file does not exist: $(graph_file)")
    end

    # Ensure vg conda environment exists
    add_bioconda_env("vg")

    index_files = Dict{String, String}()
    base_name = replace(graph_file, r"\.(gfa|vg)$" => "")

    for index_type in index_types
        index_file = "$(base_name).$(index_type)"

        if index_type == "xg"
            cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg index -x $(index_file) $(graph_file)`
        elseif index_type == "gcsa"
            # GCSA indexing requires pruning first
            pruned_file = "$(base_name).pruned.vg"
            prune_cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg prune $(graph_file)`
            run(pipeline(prune_cmd, pruned_file))

            cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg index -g $(index_file) $(pruned_file)`
        elseif index_type == "snarls"
            cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg snarls $(graph_file)`
        else
            @warn "Unknown index type: $(index_type)"
            continue
        end

        println("Building $(index_type) index for $(graph_file)")

        try
            if index_type == "snarls"
                run(pipeline(cmd, index_file))
            else
                run(cmd)
            end
            index_files[index_type] = index_file
            println("Index completed: $(index_file)")
        catch e
            @warn "Failed to build $(index_type) index: $(e)"
        end
    end

    return index_files
end
