# Pangenome Analysis Module
# Functions for comparative genomics using existing k-mer infrastructure

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
function analyze_pangenome_kmers(genome_files::Vector{String}; kmer_type=Kmers.DNAKmer{21}, distance_metric=:jaccard)
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
function compare_genome_kmer_similarity(genome1_file::String, genome2_file::String; kmer_type=Kmers.DNAKmer{21}, metric=:js_divergence)
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
function build_genome_distance_matrix(genome_files::Vector{String}; kmer_type=Kmers.DNAKmer{21}, metric=:js_divergence)
    n_genomes = length(genome_files)
    distance_matrix = zeros(Float64, n_genomes, n_genomes)
    genome_names = [basename(file) for file in genome_files]
    
    println("Building $(n_genomes)x$(n_genomes) distance matrix using $(metric)...")
    
    # Calculate pairwise distances
    for i in 1:n_genomes
        for j in (i+1):n_genomes
            println("  Comparing $(genome_names[i]) vs $(genome_names[j])...")
            similarity = compare_genome_kmer_similarity(
                genome_files[i], genome_files[j]; 
                kmer_type=kmer_type, metric=metric
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