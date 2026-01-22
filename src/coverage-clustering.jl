# Coverage-based clustering algorithms for automated error/signal separation
# Based on k-medoids analysis from Mycelia-Dev notebooks (2021-08-24-k-mediods-error-cluster-detection.ipynb)

"""
    KMedoidsCoverageResult

Results from k-medoids coverage clustering analysis.

# Fields
- `assignments::Vector{Int}`: Cluster assignments for each k-mer
- `medoids::Vector{Int}`: Indices of cluster medoids  
- `silhouettes::Vector{Float64}`: Silhouette scores for each k-mer
- `cost::Float64`: Total clustering cost
- `iterations::Int`: Number of iterations to convergence
- `coverage_thresholds::Vector{Float64}`: Coverage thresholds for each cluster
"""
struct KMedoidsCoverageResult
    assignments::Vector{Int}
    medoids::Vector{Int}
    silhouettes::Vector{Float64}
    cost::Float64
    iterations::Int
    coverage_thresholds::Vector{Float64}
end

"""
    CoverageClusteringStats

Statistical analysis of coverage-based clustering results.

# Fields
- `error_cluster_size::Int`: Size of detected error cluster
- `signal_cluster_size::Int`: Size of detected signal cluster  
- `separation_quality::Float64`: Quality of cluster separation (0-1)
- `optimal_threshold::Float64`: Recommended coverage threshold
- `false_positive_rate::Float64`: Estimated false positive rate
- `false_negative_rate::Float64`: Estimated false negative rate
"""
struct CoverageClusteringStats
    error_cluster_size::Int
    signal_cluster_size::Int
    separation_quality::Float64
    optimal_threshold::Float64
    false_positive_rate::Float64
    false_negative_rate::Float64
end

"""
    k_medoids_coverage_clustering(kmer_counts::Dict{T, Int}, k::Int=2; max_iterations::Int=100, distance_metric::Symbol=:euclidean) where T

Perform k-medoids clustering on k-mer coverage values to automatically separate errors from signal.

This algorithm clusters k-mers based on their coverage patterns to identify distinct populations
representing sequencing errors (low coverage, high variance) and true genomic signal (higher 
coverage, lower variance).

# Arguments
- `kmer_counts::Dict{T, Int}`: Dictionary mapping k-mers to their coverage counts
- `k::Int=2`: Number of clusters (typically 2 for error/signal separation)
- `max_iterations::Int=100`: Maximum iterations for convergence
- `distance_metric::Symbol=:euclidean`: Distance metric (:euclidean, :manhattan, :cosine)

# Returns
`KMedoidsCoverageResult`: Comprehensive clustering results with assignments and quality metrics

# Algorithm Details
Uses k-medoids clustering optimized for genomic coverage patterns:
1. **Coverage normalization**: Log-transforms coverage to handle exponential distributions
2. **Distance calculation**: Specialized metrics for coverage similarity
3. **Medoid selection**: Selects representative k-mers that minimize total dissimilarity
4. **Quality assessment**: Computes silhouette scores and separation metrics

# Examples
```julia
# Basic error/signal separation
kmer_counts = count_kmers(sequences, 21)
result = k_medoids_coverage_clustering(kmer_counts)

# Multi-cluster analysis with custom parameters  
result = k_medoids_coverage_clustering(kmer_counts, 3, max_iterations=200, distance_metric=:manhattan)

# Extract error cluster for filtering
error_cluster_id = argmin(result.coverage_thresholds)
error_kmers = [kmer for (kmer, assignment) in zip(keys(kmer_counts), result.assignments) 
               if assignment == error_cluster_id]
```

# References
Based on analysis in Mycelia-Dev notebook: `2021-08-24-k-mediods-error-cluster-detection.ipynb`
"""
function k_medoids_coverage_clustering(kmer_counts::Dict{T, Int}, k::Int=2; 
                                     max_iterations::Int=100, 
                                     distance_metric::Symbol=:euclidean) where T
    
    ## Convert to vectors for processing
    kmers = collect(keys(kmer_counts))
    coverages = collect(values(kmer_counts))
    n = length(kmers)
    
    if n < k
        throw(ArgumentError("Number of k-mers ($n) must be >= number of clusters ($k)"))
    end
    
    ## Log-transform coverages to handle exponential distribution
    log_coverages = log.(coverages .+ 1.0)  # Add 1 to avoid log(0)
    
    ## Create feature matrix (can be extended with additional coverage-based features)
    features = reshape(log_coverages, :, 1)
    
    ## Initialize medoids randomly
    medoid_indices = Random.randperm(n)[1:k]
    assignments = zeros(Int, n)
    previous_cost = Inf
    iterations = 0
    
    ## Main k-medoids loop
    for iter in 1:max_iterations
        iterations = iter
        
        ## Assignment step: assign each point to nearest medoid
        for i in 1:n
            min_dist = Inf
            best_cluster = 1
            
            for j in 1:k
                medoid_idx = medoid_indices[j]
                dist = _calculate_distance(features[i, :], features[medoid_idx, :], distance_metric)
                
                if dist < min_dist
                    min_dist = dist
                    best_cluster = j
                end
            end
            
            assignments[i] = best_cluster
        end
        
        ## Update step: find new medoids that minimize within-cluster cost
        new_medoids = copy(medoid_indices)
        
        for cluster_id in 1:k
            cluster_points = findall(x -> x == cluster_id, assignments)
            
            if isempty(cluster_points)
                continue  ## Keep current medoid if cluster is empty
            end
            
            best_medoid = cluster_points[1]
            min_total_cost = Inf
            
            ## Try each point in cluster as potential medoid
            for candidate in cluster_points
                total_cost = 0.0
                
                for point in cluster_points
                    total_cost += _calculate_distance(features[candidate, :], features[point, :], distance_metric)
                end
                
                if total_cost < min_total_cost
                    min_total_cost = total_cost
                    best_medoid = candidate
                end
            end
            
            new_medoids[cluster_id] = best_medoid
        end
        
        ## Calculate total cost
        current_cost = 0.0
        for i in 1:n
            medoid_idx = new_medoids[assignments[i]]
            current_cost += _calculate_distance(features[i, :], features[medoid_idx, :], distance_metric)
        end
        
        ## Check for convergence
        if abs(current_cost - previous_cost) < 1e-6
            medoid_indices = new_medoids
            break
        end
        
        medoid_indices = new_medoids
        previous_cost = current_cost
    end
    
    ## Calculate silhouette scores
    silhouettes = _calculate_silhouettes(features, assignments, medoid_indices, distance_metric)
    
    ## Calculate coverage thresholds for each cluster
    coverage_thresholds = zeros(Float64, k)
    for cluster_id in 1:k
        cluster_coverages = coverages[assignments .== cluster_id]
        coverage_thresholds[cluster_id] = isempty(cluster_coverages) ? 0.0 : Statistics.mean(cluster_coverages)
    end
    
    ## Final cost calculation
    final_cost = 0.0
    for i in 1:n
        medoid_idx = medoid_indices[assignments[i]]
        final_cost += _calculate_distance(features[i, :], features[medoid_idx, :], distance_metric)
    end
    
    return KMedoidsCoverageResult(
        assignments,
        medoid_indices,
        silhouettes,
        final_cost,
        iterations,
        coverage_thresholds
    )
end

"""
    analyze_coverage_clustering(result::KMedoidsCoverageResult, kmer_counts::Dict{T, Int}) where T

Analyze k-medoids clustering results to extract error/signal separation statistics.

# Arguments
- `result::KMedoidsCoverageResult`: Results from k-medoids clustering
- `kmer_counts::Dict{T, Int}`: Original k-mer coverage counts

# Returns
`CoverageClusteringStats`: Statistical analysis of clustering quality and separation
"""
function analyze_coverage_clustering(result::KMedoidsCoverageResult, kmer_counts::Dict{T, Int}) where T
    
    coverages = collect(values(kmer_counts))
    n_clusters = length(result.medoids)
    
    if n_clusters != 2
        @warn "Analysis optimized for 2-cluster error/signal separation, found $n_clusters clusters"
    end
    
    ## Identify error vs signal clusters based on coverage thresholds
    error_cluster_id = argmin(result.coverage_thresholds)
    signal_cluster_id = argmax(result.coverage_thresholds)
    
    error_cluster_size = count(x -> x == error_cluster_id, result.assignments)
    signal_cluster_size = count(x -> x == signal_cluster_id, result.assignments)
    
    ## Calculate separation quality using silhouette analysis
    separation_quality = Statistics.mean(result.silhouettes)
    
    ## Determine optimal threshold as midpoint between cluster means
    error_threshold = result.coverage_thresholds[error_cluster_id]
    signal_threshold = result.coverage_thresholds[signal_cluster_id]
    optimal_threshold = (error_threshold + signal_threshold) / 2.0
    
    ## Estimate false positive/negative rates based on threshold overlap
    error_coverages = [coverages[i] for i in 1:length(coverages) if result.assignments[i] == error_cluster_id]
    signal_coverages = [coverages[i] for i in 1:length(coverages) if result.assignments[i] == signal_cluster_id]
    
    ## False positive rate: fraction of error cluster above optimal threshold
    false_positive_rate = isempty(error_coverages) ? 0.0 : 
                         count(x -> x >= optimal_threshold, error_coverages) / length(error_coverages)
    
    ## False negative rate: fraction of signal cluster below optimal threshold  
    false_negative_rate = isempty(signal_coverages) ? 0.0 :
                         count(x -> x < optimal_threshold, signal_coverages) / length(signal_coverages)
    
    return CoverageClusteringStats(
        error_cluster_size,
        signal_cluster_size,
        separation_quality,
        optimal_threshold,
        false_positive_rate,
        false_negative_rate
    )
end

"""
    automatic_error_filtering(kmer_counts::Dict{T, Int}; 
                            separation_threshold::Float64=0.3,
                            max_error_rate::Float64=0.1) where T

Automatically filter error k-mers using k-medoids coverage clustering.

# Arguments
- `kmer_counts::Dict{T, Int}`: Dictionary of k-mer coverage counts
- `separation_threshold::Float64=0.3`: Minimum silhouette score for reliable separation
- `max_error_rate::Float64=0.1`: Maximum acceptable false positive rate

# Returns
- `filtered_kmers::Dict{T, Int}`: K-mers with errors removed
- `stats::CoverageClusteringStats`: Clustering analysis statistics
- `removed_kmers::Vector{T}`: K-mers identified as errors and removed
"""
function automatic_error_filtering(kmer_counts::Dict{T, Int}; 
                                 separation_threshold::Float64=0.3,
                                 max_error_rate::Float64=0.1) where T
    
    ## Perform k-medoids clustering
    result = k_medoids_coverage_clustering(kmer_counts, 2)
    stats = analyze_coverage_clustering(result, kmer_counts)
    
    ## Check if clustering provides reliable separation
    if stats.separation_quality < separation_threshold
        @warn "Poor cluster separation ($(stats.separation_quality) < $separation_threshold). " *
              "Consider manual threshold selection."
    end
    
    if stats.false_positive_rate > max_error_rate
        @warn "High false positive rate ($(stats.false_positive_rate) > $max_error_rate). " *
              "Results may remove valid k-mers."
    end
    
    ## Identify error cluster
    error_cluster_id = argmin(result.coverage_thresholds)
    
    ## Filter k-mers
    kmers = collect(keys(kmer_counts))
    removed_kmers = T[]
    filtered_kmers = Dict{T, Int}()
    
    for (i, kmer) in enumerate(kmers)
        if result.assignments[i] == error_cluster_id
            push!(removed_kmers, kmer)
        else
            filtered_kmers[kmer] = kmer_counts[kmer]
        end
    end
    
    @info "Removed $(length(removed_kmers)) error k-mers, retained $(length(filtered_kmers)) signal k-mers"
    
    return filtered_kmers, stats, removed_kmers
end

## Helper function for distance calculations
function _calculate_distance(x::Vector{Float64}, y::Vector{Float64}, metric::Symbol)
    if metric == :euclidean
        return sqrt(sum((x .- y).^2))
    elseif metric == :manhattan
        return sum(abs.(x .- y))
    elseif metric == :cosine
        dot_product = sum(x .* y)
        norm_x = sqrt(sum(x.^2))
        norm_y = sqrt(sum(y.^2))
        denom = norm_x * norm_y
        if denom == 0.0
            return (norm_x == 0.0 && norm_y == 0.0) ? 0.0 : 1.0
        end
        return 1.0 - (dot_product / denom)
    else
        throw(ArgumentError("Unsupported distance metric: $metric"))
    end
end

## Helper function for silhouette score calculation
function _calculate_silhouettes(features::Matrix{Float64}, assignments::Vector{Int}, 
                               medoid_indices::Vector{Int}, distance_metric::Symbol)
    n, _ = size(features)
    k = length(medoid_indices)
    silhouettes = zeros(Float64, n)
    
    for i in 1:n
        cluster_i = assignments[i]
        
        ## Calculate average distance to points in same cluster (a_i)
        same_cluster_distances = Float64[]
        for j in 1:n
            if i != j && assignments[j] == cluster_i
                dist = _calculate_distance(features[i, :], features[j, :], distance_metric)
                push!(same_cluster_distances, dist)
            end
        end
        
        a_i = isempty(same_cluster_distances) ? 0.0 : Statistics.mean(same_cluster_distances)
        
        ## Calculate minimum average distance to points in other clusters (b_i)
        b_i = Inf
        for other_cluster in 1:k
            if other_cluster == cluster_i
                continue
            end
            
            other_cluster_distances = Float64[]
            for j in 1:n
                if assignments[j] == other_cluster
                    dist = _calculate_distance(features[i, :], features[j, :], distance_metric)
                    push!(other_cluster_distances, dist)
                end
            end
            
            if !isempty(other_cluster_distances)
                avg_dist = Statistics.mean(other_cluster_distances)
                b_i = min(b_i, avg_dist)
            end
        end
        
        ## Calculate silhouette score
        if isinf(b_i) || (a_i == 0.0 && b_i == 0.0)
            silhouettes[i] = 0.0
        else
            silhouettes[i] = (b_i - a_i) / max(a_i, b_i)
        end
    end
    
    return silhouettes
end
