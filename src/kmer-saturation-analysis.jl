# K-mer saturation curve fitting for optimal k selection
# Based on analysis from Mycelia-Dev notebooks (2021-09-15-sequencing-saturation.ipynb)

"""
    SaturationCurveResult

Results from k-mer saturation curve analysis.

# Fields
- `optimal_k::Int`: Recommended optimal k-mer size
- `saturation_levels::Vector{Float64}`: Saturation level for each tested k
- `curve_parameters::Dict{Int, Vector{Float64}}`: Michaelis-Menten parameters [vmax, km] for each k
- `sampling_points::Vector{Int}`: Number of sequences sampled for curve fitting
- `unique_kmer_counts::Dict{Int, Vector{Int}}`: Unique k-mer counts at each sampling point for each k
- `convergence_quality::Dict{Int, Float64}`: R² values for curve fits
"""
struct SaturationCurveResult
    optimal_k::Int
    saturation_levels::Vector{Float64}
    curve_parameters::Dict{Int, Vector{Float64}}
    sampling_points::Vector{Int}
    unique_kmer_counts::Dict{Int, Vector{Int}}
    convergence_quality::Dict{Int, Float64}
end

"""
    ConnectivityAnalysisResult

Results from connectivity-based assembly threshold analysis.

# Fields
- `optimal_threshold::Int`: Recommended coverage threshold for assembly
- `connectivity_scores::Vector{Float64}`: Connectivity quality scores for each threshold
- `component_counts::Vector{Int}`: Number of connected components at each threshold
- `average_connectivity::Vector{Float64}`: Average node connectivity at each threshold
- `thresholds_tested::Vector{Int}`: Coverage thresholds that were evaluated
"""
struct ConnectivityAnalysisResult
    optimal_threshold::Int
    connectivity_scores::Vector{Float64}
    component_counts::Vector{Int}
    average_connectivity::Vector{Float64}
    thresholds_tested::Vector{Int}
end

"""
    analyze_kmer_saturation(sequences::Vector{T}, k_range::Vector{Int}; 
                           max_sampling_points::Int=20, 
                           min_sequences::Int=100) where T

Analyze k-mer saturation curves to determine optimal k-mer size for assembly.

This algorithm uses Michaelis-Menten curve fitting to model k-mer saturation patterns
and identifies the optimal k-mer size that provides sufficient diversity without 
excessive sparsity.

# Arguments
- `sequences::Vector{T}`: Vector of DNA sequences (BioSequences or strings)
- `k_range::Vector{Int}`: Range of k-mer sizes to test
- `max_sampling_points::Int=20`: Maximum number of sampling points for curve fitting
- `min_sequences::Int=100`: Minimum number of sequences required for analysis

# Returns
`SaturationCurveResult`: Comprehensive analysis with optimal k recommendation

# Algorithm Details
1. **Progressive Sampling**: Incrementally samples sequences to build saturation curves
2. **Michaelis-Menten Fitting**: Models k-mer saturation using enzyme kinetics equation: `v = (vmax * s) / (km + s)`
3. **Saturation Assessment**: Calculates saturation level as fraction of theoretical maximum k-mers
4. **Optimal Selection**: Chooses k with < 10% saturation for sufficient sparsity

# Mathematical Foundation
The Michaelis-Menten equation models k-mer discovery rate:
- `vmax`: Maximum number of unique k-mers (saturation point)
- `km`: Number of sequences at half-saturation
- Lower saturation indicates better discrimination power

# Examples
```julia
# Analyze k-mer saturation for assembly
sequences = load_fasta_sequences("genome.fasta")
k_range = [21, 31, 51, 71, 91]
result = analyze_kmer_saturation(sequences, k_range)

# Use optimal k for downstream analysis
optimal_k = result.optimal_k
kmer_counts = count_kmers(sequences, optimal_k)
```

# References
Based on Michaelis-Menten k-mer saturation analysis from: 
`2021-09-15-sequencing-saturation.ipynb`
"""
function analyze_kmer_saturation(sequences::Vector{T}, k_range::Vector{Int}; 
                                max_sampling_points::Int=20, 
                                min_sequences::Int=100) where T
    
    if length(sequences) < min_sequences
        throw(ArgumentError("Need at least $min_sequences sequences, got $(length(sequences))"))
    end
    
    if isempty(k_range)
        throw(ArgumentError("k_range cannot be empty"))
    end
    
    n_sequences = length(sequences)
    
    ## Create sampling points (logarithmic spacing for better curve resolution)
    sampling_points = unique([round(Int, 10^x) for x in range(log10(min_sequences), 
                                                             log10(n_sequences), 
                                                             length=max_sampling_points)])
    filter!(x -> x <= n_sequences, sampling_points)
    sort!(sampling_points)
    
    @info "Testing k-mer saturation with $(length(k_range)) k values and $(length(sampling_points)) sampling points"
    
    ## Initialize results storage
    unique_kmer_counts = Dict{Int, Vector{Int}}()
    curve_parameters = Dict{Int, Vector{Float64}}()
    convergence_quality = Dict{Int, Float64}()
    saturation_levels = Float64[]
    
    ## Analyze each k-mer size
    for k in k_range
        @info "Analyzing k=$k saturation..."
        
        kmer_counts_at_sampling = Int[]
        
        ## Progressive sampling to build saturation curve
        for n_sample in sampling_points
            sampled_sequences = sequences[1:n_sample]
            
            ## Count unique k-mers at this sampling point
            kmer_set = Set{String}()
            
            for seq in sampled_sequences
                seq_str = string(seq)
                if length(seq_str) >= k
                    for i in 1:(length(seq_str) - k + 1)
                        kmer = seq_str[i:(i + k - 1)]
                        ## Add canonical k-mer (lexicographically smaller of forward/reverse complement)
                        canonical_kmer = _get_canonical_kmer(kmer)
                        push!(kmer_set, canonical_kmer)
                    end
                end
            end
            
            push!(kmer_counts_at_sampling, length(kmer_set))
        end
        
        unique_kmer_counts[k] = kmer_counts_at_sampling
        
        ## Fit Michaelis-Menten curve: v = (vmax * s) / (km + s)
        try
            params, r_squared = _fit_michaelis_menten(sampling_points, kmer_counts_at_sampling)
            curve_parameters[k] = params
            convergence_quality[k] = r_squared
            
            ## Calculate saturation level
            vmax = params[1]
            theoretical_max = 4^k  ## Maximum possible k-mers for DNA alphabet
            saturation = vmax / theoretical_max
            push!(saturation_levels, saturation)
            
            @info "k=$k: vmax=$(round(Int, vmax)), saturation=$(round(100*saturation, digits=2))%, R²=$(round(r_squared, digits=3))"
            
        catch e
            @warn "Failed to fit curve for k=$k: $e"
            curve_parameters[k] = [NaN, NaN]
            convergence_quality[k] = 0.0
            push!(saturation_levels, 1.0)  ## Assume saturated if fitting fails
        end
    end
    
    ## Select optimal k (lowest saturation with good fit quality)
    optimal_k_index = 1
    best_score = Inf
    
    for (i, k) in enumerate(k_range)
        saturation = saturation_levels[i]
        r_squared = get(convergence_quality, k, 0.0)
        
        ## Scoring function: prefer low saturation and good fit quality
        score = saturation - 0.1 * r_squared  ## Weight R² positively
        
        if score < best_score && r_squared > 0.5  ## Require reasonable fit quality
            best_score = score
            optimal_k_index = i
        end
    end
    
    optimal_k = k_range[optimal_k_index]
    
    @info "Optimal k-mer size: $optimal_k (saturation: $(round(100*saturation_levels[optimal_k_index], digits=2))%)"
    
    return SaturationCurveResult(
        optimal_k,
        saturation_levels,
        curve_parameters,
        sampling_points,
        unique_kmer_counts,
        convergence_quality
    )
end

"""
    find_optimal_assembly_threshold(kmer_counts::Dict{T, Int}, k::Int;
                                   min_threshold::Int=1,
                                   max_threshold::Int=100) where T

Find optimal coverage threshold for assembly using connectivity analysis.

This algorithm tests different coverage thresholds and selects the one that optimizes
graph connectivity while maintaining reasonable assembly complexity.

# Arguments
- `kmer_counts::Dict{T, Int}`: K-mer coverage counts
- `k::Int`: K-mer size used
- `min_threshold::Int=1`: Minimum coverage threshold to test
- `max_threshold::Int=100`: Maximum coverage threshold to test

# Returns
`ConnectivityAnalysisResult`: Analysis results with optimal threshold recommendation

# Algorithm Details
1. **Threshold Testing**: Tests exponentially spaced coverage thresholds
2. **Graph Construction**: Builds k-mer overlap graphs at each threshold  
3. **Connectivity Analysis**: Measures connected components and average connectivity
4. **Optimization**: Selects threshold minimizing `components × |connectivity - 2|`

# Optimization Criterion
The optimal threshold minimizes graph fragmentation while targeting connectivity ≈ 2:
- **Connectivity = 2**: Indicates linear paths (ideal for assembly)
- **Higher connectivity**: Over-collapsed repetitive regions
- **Lower connectivity**: Under-collapsed, fragmented assembly

# Examples
```julia
# Find optimal threshold for assembly
kmer_counts = count_kmers(sequences, 31)
result = find_optimal_assembly_threshold(kmer_counts, 31)

# Filter k-mers using optimal threshold
filtered_kmers = filter(p -> p.second >= result.optimal_threshold, kmer_counts)
```
"""
function find_optimal_assembly_threshold(kmer_counts::Dict{T, Int}, k::Int;
                                       min_threshold::Int=1,
                                       max_threshold::Int=100) where T
    
    if isempty(kmer_counts)
        throw(ArgumentError("kmer_counts cannot be empty"))
    end
    
    ## Determine reasonable threshold range based on data
    max_coverage = maximum(values(kmer_counts))
    actual_max_threshold = min(max_threshold, max_coverage)
    
    ## Create exponentially spaced thresholds
    thresholds = unique([round(Int, 2^x) for x in range(0, log2(actual_max_threshold), length=15)])
    filter!(t -> t >= min_threshold && t <= actual_max_threshold, thresholds)
    sort!(thresholds)
    
    @info "Testing $(length(thresholds)) coverage thresholds from $min_threshold to $actual_max_threshold"
    
    connectivity_scores = Float64[]
    component_counts = Int[]
    average_connectivity = Float64[]
    
    ## Test each threshold
    for threshold in thresholds
        ## Filter k-mers by coverage threshold
        filtered_kmers = filter(p -> p.second >= threshold, kmer_counts)
        
        if length(filtered_kmers) < 10
            ## Too few k-mers remaining - assign poor score
            push!(connectivity_scores, Inf)
            push!(component_counts, 0)
            push!(average_connectivity, 0.0)
            continue
        end
        
        ## Build overlap graph from filtered k-mers
        kmer_list = collect(keys(filtered_kmers))
        overlap_graph = _build_kmer_overlap_graph(kmer_list, k)
        
        ## Calculate connectivity metrics
        n_components = _count_connected_components(overlap_graph)
        avg_connectivity = _calculate_average_connectivity(overlap_graph)
        
        ## Score function: minimize fragmentation while targeting connectivity ≈ 2
        connectivity_penalty = abs(avg_connectivity - 2.0)
        score = n_components * (1.0 + connectivity_penalty)
        
        push!(connectivity_scores, score)
        push!(component_counts, n_components)
        push!(average_connectivity, avg_connectivity)
        
        @debug "Threshold $threshold: $(length(filtered_kmers)) k-mers, $n_components components, " *
               "connectivity $(round(avg_connectivity, digits=2)), score $(round(score, digits=2))"
    end
    
    ## Find optimal threshold
    optimal_index = argmin(connectivity_scores)
    optimal_threshold = thresholds[optimal_index]
    
    @info "Optimal coverage threshold: $optimal_threshold " *
          "($(component_counts[optimal_index]) components, " *
          "connectivity $(round(average_connectivity[optimal_index], digits=2)))"
    
    return ConnectivityAnalysisResult(
        optimal_threshold,
        connectivity_scores,
        component_counts,
        average_connectivity,
        thresholds
    )
end

## Helper function for Michaelis-Menten curve fitting
function _fit_michaelis_menten(x::Vector{Int}, y::Vector{Int})
    ## Convert to Float64 for numerical stability
    x_float = Float64.(x)
    y_float = Float64.(y)
    
    ## Initial parameter estimates
    vmax_init = maximum(y_float) * 1.2  ## Slightly above maximum observed
    km_init = x_float[findmin(abs.(y_float .- vmax_init/2))[2]]  ## x at half-max
    
    ## Simple gradient descent fitting (could be replaced with LsqFit.jl for better performance)
    best_vmax, best_km = vmax_init, km_init
    best_sse = Inf
    
    ## Grid search around initial estimates
    for vmax in range(vmax_init * 0.8, vmax_init * 1.5, length=10)
        for km in range(km_init * 0.5, km_init * 2.0, length=10)
            ## Calculate predicted values
            y_pred = [(vmax * xi) / (km + xi) for xi in x_float]
            
            ## Sum of squared errors
            sse = sum((y_float .- y_pred).^2)
            
            if sse < best_sse
                best_sse = sse
                best_vmax, best_km = vmax, km
            end
        end
    end
    
    ## Calculate R-squared
    y_pred = [(best_vmax * xi) / (best_km + xi) for xi in x_float]
    ss_res = sum((y_float .- y_pred).^2)
    ss_tot = sum((y_float .- Statistics.mean(y_float)).^2)
    r_squared = 1.0 - (ss_res / (ss_tot + 1e-10))
    
    return [best_vmax, best_km], r_squared
end

## Helper function to build k-mer overlap graph
function _build_kmer_overlap_graph(kmers::Vector{T}, k::Int) where T
    kmer_strings = [string(kmer) for kmer in kmers]
    kmer_to_index = Dict(kmer => i for (i, kmer) in enumerate(kmer_strings))
    n = length(kmer_strings)
    
    ## Build adjacency list representation
    adjacency = [Int[] for _ in 1:n]
    
    for (i, kmer) in enumerate(kmer_strings)
        ## Find overlapping k-mers (suffix of length k-1 matches prefix of length k-1)
        suffix = kmer[2:end]
        
        for base in ['A', 'C', 'G', 'T']
            candidate = suffix * base
            if haskey(kmer_to_index, candidate) && candidate != kmer
                j = kmer_to_index[candidate]
                push!(adjacency[i], j)
            end
        end
    end
    
    return adjacency
end

## Helper function to count connected components using DFS
function _count_connected_components(adjacency::Vector{Vector{Int}})
    n = length(adjacency)
    visited = falses(n)
    components = 0
    
    for i in 1:n
        if !visited[i]
            components += 1
            ## DFS to mark all nodes in this component
            stack = [i]
            while !isempty(stack)
                node = pop!(stack)
                if !visited[node]
                    visited[node] = true
                    append!(stack, adjacency[node])
                end
            end
        end
    end
    
    return components
end

## Helper function to calculate average connectivity
function _calculate_average_connectivity(adjacency::Vector{Vector{Int}})
    if isempty(adjacency)
        return 0.0
    end
    
    total_degree = sum(length(neighbors) for neighbors in adjacency)
    return total_degree / length(adjacency)
end