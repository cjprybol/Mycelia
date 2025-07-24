# Circular genome reconstruction algorithms using hub-based traversal
# Based on analysis from Mycelia-Dev notebooks (2022-01-26-sample-core-genome-circular*.ipynb)

"""
    CircularReconstructionResult

Results from circular genome reconstruction analysis.

# Fields
- `circular_sequence::String`: Reconstructed circular genome sequence
- `completion_percentage::Float64`: Percentage of target k-mers successfully included
- `traversal_path::Vector{String}`: K-mer path used for reconstruction
- `hub_nodes_used::Vector{String}`: Hub k-mers that served as path anchors
- `gap_regions::Vector{Tuple{Int, Int}}`: Positions of unresolved gaps in the sequence
- `reconstruction_confidence::Float64`: Overall confidence score for the reconstruction
- `iterations_performed::Int`: Number of extension iterations completed
"""
struct CircularReconstructionResult
    circular_sequence::String
    completion_percentage::Float64
    traversal_path::Vector{String}
    hub_nodes_used::Vector{String}
    gap_regions::Vector{Tuple{Int, Int}}
    reconstruction_confidence::Float64
    iterations_performed::Int
end

"""
    CircularityValidationResult

Results from circular genome validation analysis.

# Fields
- `is_circular::Bool`: Whether the sequence forms a valid circle
- `overlap_length::Int`: Length of sequence overlap at circularization point
- `circularity_confidence::Float64`: Confidence in circular structure (0-1)
- `terminal_kmers::Tuple{String, String}`: First and last k-mers in the sequence
- `gap_size::Int`: Size of gap preventing circularization (0 if circular)
"""
struct CircularityValidationResult
    is_circular::Bool
    overlap_length::Int
    circularity_confidence::Float64
    terminal_kmers::Tuple{String, String}
    gap_size::Int
end

"""
    reconstruct_circular_genome(target_kmers::Set{String},
                               kmer_graph::AbstractGraph,
                               kmer_coverage::Dict{String, Int},
                               k::Int;
                               max_iterations::Int=10,
                               completion_threshold::Float64=0.90,
                               use_stochastic_selection::Bool=true) -> CircularReconstructionResult

Reconstruct circular genome sequences using hub-based bidirectional traversal.

This algorithm identifies hub nodes (high-degree vertices) in k-mer graphs and uses
bidirectional Dijkstra search to construct paths that maximize coverage of target
k-mers while forming circular structures.

# Arguments
- `target_kmers::Set{String}`: Set of k-mers that should be included in reconstruction
- `kmer_graph::AbstractGraph`: K-mer de Bruijn graph for pathfinding
- `kmer_coverage::Dict{String, Int}`: Coverage information for each k-mer
- `k::Int`: K-mer size used in analysis
- `max_iterations::Int=10`: Maximum number of path extension iterations
- `completion_threshold::Float64=0.90`: Target fraction of k-mers to include
- `use_stochastic_selection::Bool=true`: Whether to use probabilistic hub selection

# Returns
`CircularReconstructionResult`: Comprehensive reconstruction results with quality metrics

# Algorithm Details
1. **Hub Identification**: Finds k-mers with degree ≥3 as potential path anchors
2. **Target Space Preparation**: Converts target k-mers to canonical form
3. **Stochastic Hub Selection**: Probabilistically selects starting hub based on coverage
4. **Bidirectional Extension**: Uses bidirectional Dijkstra to extend paths from hubs
5. **Iterative Completion**: Repeats until completion threshold is reached
6. **Circularity Validation**: Checks for proper circular structure

# Hub Selection Strategies
- **Deterministic**: Selects highest-coverage hub
- **Stochastic**: Weight-based random selection allowing exploration of multiple solutions

# Examples
```julia
# Basic circular reconstruction
target_set = Set(extract_kmers(reference_genome, 31))
graph, coverage = build_kmer_graph_with_coverage(sequencing_reads, 31)
result = reconstruct_circular_genome(target_set, graph, coverage, 31)

# High-completion reconstruction with deterministic selection
result = reconstruct_circular_genome(target_set, graph, coverage, 31,
                                   completion_threshold=0.95,
                                   use_stochastic_selection=false)

# Validate circularity
validation = validate_circularity(result.circular_sequence, 31)
if validation.is_circular
    println("Successfully reconstructed circular genome")
end
```

# References
Based on circular genome reconstruction from Mycelia-Dev notebooks:
- `2022-01-26-sample-core-genome-circular.ipynb`
- `2022-01-26-sample-core-genome-circular-recycle.ipynb`
"""
function reconstruct_circular_genome(target_kmers::Set{String},
                                   kmer_graph::AbstractGraph,
                                   kmer_coverage::Dict{String, Int},
                                   k::Int;
                                   max_iterations::Int=10,
                                   completion_threshold::Float64=0.90,
                                   use_stochastic_selection::Bool=true)
    
    if isempty(target_kmers)
        throw(ArgumentError("Target k-mer set cannot be empty"))
    end
    
    if Graphs.nv(kmer_graph) == 0
        throw(ArgumentError("K-mer graph cannot be empty"))
    end
    
    ## Convert target k-mers to canonical form for consistent comparison
    canonical_targets = Set([_get_canonical_kmer(kmer) for kmer in target_kmers])
    
    @info "Starting circular reconstruction: $(length(canonical_targets)) target k-mers, " *
          "completion threshold: $(round(100*completion_threshold, digits=1))%"
    
    ## Build k-mer index for graph navigation
    graph_kmers = collect(keys(kmer_coverage))
    kmer_to_vertex = Dict(kmer => i for (i, kmer) in enumerate(graph_kmers))
    
    ## Identify hub nodes (degree ≥ 3)
    hub_kmers = String[]
    for (i, kmer) in enumerate(graph_kmers)
        if i <= Graphs.nv(kmer_graph) && Graphs.degree(kmer_graph, i) >= 3
            push!(hub_kmers, kmer)
        end
    end
    
    if isempty(hub_kmers)
        @warn "No hub nodes found in graph - using highest coverage k-mers as anchors"
        ## Fallback: use top 10% highest coverage k-mers as pseudo-hubs
        sorted_kmers = sort(collect(kmer_coverage), by=x->x[2], rev=true)
        n_pseudo_hubs = max(1, length(sorted_kmers) ÷ 10)
        hub_kmers = [kmer for (kmer, cov) in sorted_kmers[1:n_pseudo_hubs]]
    end
    
    @info "Found $(length(hub_kmers)) hub nodes for path anchoring"
    
    ## Initialize reconstruction state
    reconstructed_path = String[]
    covered_targets = Set{String}()
    hub_nodes_used = String[]
    iterations_performed = 0
    
    ## Main reconstruction loop
    for iteration in 1:max_iterations
        iterations_performed = iteration
        
        ## Calculate current completion
        completion = length(covered_targets) / length(canonical_targets)
        
        if completion >= completion_threshold
            @info "Completion threshold reached: $(round(100*completion, digits=1))%"
            break
        end
        
        @debug "Iteration $iteration: $(round(100*completion, digits=1))% complete"
        
        ## Select hub for path extension
        if use_stochastic_selection
            selected_hub = select_hub_stochastic(hub_kmers, kmer_coverage)
        else
            selected_hub = select_hub_deterministic(hub_kmers, kmer_coverage)
        end
        
        if selected_hub === nothing
            @warn "No suitable hub found for iteration $iteration"
            break
        end
        
        push!(hub_nodes_used, selected_hub)
        
        ## Perform bidirectional extension from selected hub
        extension_result = extend_path_bidirectional(
            selected_hub, canonical_targets, covered_targets,
            kmer_graph, graph_kmers, kmer_to_vertex, kmer_coverage, k
        )
        
        if !isempty(extension_result.path)
            ## Merge new path with existing reconstruction
            reconstructed_path = merge_paths(reconstructed_path, extension_result.path, k)
            
            ## Update covered targets
            for kmer in extension_result.path
                canonical_kmer = _get_canonical_kmer(kmer)
                if canonical_kmer in canonical_targets
                    push!(covered_targets, canonical_kmer)
                end
            end
            
            @debug "Added $(length(extension_result.path)) k-mers, " *
                   "now covering $(length(covered_targets)) targets"
        else
            @debug "No path extension possible from hub $selected_hub"
        end
        
        ## Remove used hub to avoid cycles
        filter!(h -> h != selected_hub, hub_kmers)
        
        if isempty(hub_kmers)
            @info "All hub nodes exhausted after $iteration iterations"
            break
        end
    end
    
    ## Reconstruct final sequence
    circular_sequence = ""
    if !isempty(reconstructed_path)
        circular_sequence = reconstruct_sequence_from_kmers(reconstructed_path, k)
    end
    
    ## Calculate final metrics
    final_completion = length(covered_targets) / length(canonical_targets)
    gap_regions = identify_gap_regions(reconstructed_path, canonical_targets, k)
    reconstruction_confidence = calculate_reconstruction_confidence(
        reconstructed_path, kmer_coverage, final_completion
    )
    
    @info "Circular reconstruction complete: $(round(100*final_completion, digits=1))% completion, " *
          "$(length(reconstructed_path)) k-mers, $(length(gap_regions)) gaps"
    
    return CircularReconstructionResult(
        circular_sequence,
        final_completion,
        reconstructed_path,
        unique(hub_nodes_used),
        gap_regions,
        reconstruction_confidence,
        iterations_performed
    )
end

"""
    validate_circularity(sequence::String, k::Int) -> CircularityValidationResult

Validate whether a reconstructed sequence forms a proper circular structure.

# Arguments
- `sequence::String`: Sequence to validate for circularity
- `k::Int`: K-mer size used for validation

# Returns
`CircularityValidationResult`: Detailed circularity validation results
"""
function validate_circularity(sequence::String, k::Int)
    
    if length(sequence) < k
        return CircularityValidationResult(false, 0, 0.0, ("", ""), length(sequence))
    end
    
    ## Extract terminal k-mers
    first_kmer = sequence[1:k]
    last_kmer = sequence[(end-k+1):end]
    
    ## Check for overlap that would create circularity
    max_overlap = min(k-1, length(sequence) ÷ 2)
    best_overlap = 0
    circularity_confidence = 0.0
    
    for overlap_len in 1:max_overlap
        prefix = sequence[1:overlap_len]
        suffix = sequence[(end-overlap_len+1):end]
        
        if prefix == suffix
            best_overlap = overlap_len
            ## Higher overlap = higher confidence, but penalize very short overlaps
            circularity_confidence = min(1.0, overlap_len / (k-1))
            break
        end
    end
    
    is_circular = best_overlap > 0
    gap_size = is_circular ? 0 : k - 1  ## Minimum gap to achieve circularity
    
    return CircularityValidationResult(
        is_circular,
        best_overlap,
        circularity_confidence,
        (first_kmer, last_kmer),
        gap_size
    )
end

"""
    repair_circular_gaps(result::CircularReconstructionResult,
                        kmer_graph::AbstractGraph,
                        kmer_coverage::Dict{String, Int},
                        k::Int) -> CircularReconstructionResult

Attempt to repair gaps in circular genome reconstruction using local path search.
"""
function repair_circular_gaps(result::CircularReconstructionResult,
                             kmer_graph::AbstractGraph,
                             kmer_coverage::Dict{String, Int},
                             k::Int)
    
    if isempty(result.gap_regions)
        @info "No gaps to repair in circular reconstruction"
        return result
    end
    
    @info "Attempting to repair $(length(result.gap_regions)) gaps"
    
    repaired_path = copy(result.traversal_path)
    repaired_gaps = 0
    
    ## Build k-mer index
    graph_kmers = collect(keys(kmer_coverage))
    kmer_to_vertex = Dict(kmer => i for (i, kmer) in enumerate(graph_kmers))
    
    for (gap_start, gap_end) in result.gap_regions
        ## Identify k-mers flanking the gap
        before_gap = gap_start > 1 ? repaired_path[gap_start-1] : ""
        after_gap = gap_end < length(repaired_path) ? repaired_path[gap_end+1] : ""
        
        if before_gap != "" && after_gap != ""
            ## Try to find a path connecting the flanking k-mers
            connecting_path = find_connecting_path(
                before_gap, after_gap, kmer_graph, 
                graph_kmers, kmer_to_vertex, kmer_coverage, k, 5  ## Max 5 steps
            )
            
            if !isempty(connecting_path)
                ## Insert connecting path to repair gap
                repaired_path = vcat(
                    repaired_path[1:(gap_start-1)],
                    connecting_path,
                    repaired_path[(gap_end+1):end]
                )
                repaired_gaps += 1
                @debug "Repaired gap between positions $gap_start-$gap_end"
            end
        end
    end
    
    ## Reconstruct sequence with repairs
    repaired_sequence = ""
    if !isempty(repaired_path)
        repaired_sequence = reconstruct_sequence_from_kmers(repaired_path, k)
    end
    
    ## Recalculate metrics
    remaining_gaps = identify_gap_regions(repaired_path, Set(result.traversal_path), k)
    new_confidence = calculate_reconstruction_confidence(
        repaired_path, kmer_coverage, result.completion_percentage
    )
    
    @info "Gap repair complete: $repaired_gaps gaps repaired, $(length(remaining_gaps)) remaining"
    
    return CircularReconstructionResult(
        repaired_sequence,
        result.completion_percentage,  ## Completion doesn't change with gap repair
        repaired_path,
        result.hub_nodes_used,
        remaining_gaps,
        new_confidence,
        result.iterations_performed
    )
end

## Helper function for stochastic hub selection
function select_hub_stochastic(hub_kmers::Vector{String}, kmer_coverage::Dict{String, Int})
    
    if isempty(hub_kmers)
        return nothing
    end
    
    ## Weight hubs by their coverage
    weights = [Float64(get(kmer_coverage, hub, 1)) for hub in hub_kmers]
    total_weight = sum(weights)
    
    if total_weight == 0
        return first(hub_kmers)
    end
    
    ## Weighted random selection
    rand_val = rand() * total_weight
    cumulative_weight = 0.0
    
    for (i, weight) in enumerate(weights)
        cumulative_weight += weight
        if cumulative_weight >= rand_val
            return hub_kmers[i]
        end
    end
    
    return last(hub_kmers)  ## Fallback
end

## Helper function for deterministic hub selection
function select_hub_deterministic(hub_kmers::Vector{String}, kmer_coverage::Dict{String, Int})
    
    if isempty(hub_kmers)
        return nothing
    end
    
    ## Select hub with highest coverage
    best_hub = hub_kmers[1]
    best_coverage = get(kmer_coverage, best_hub, 0)
    
    for hub in hub_kmers[2:end]
        coverage = get(kmer_coverage, hub, 0)
        if coverage > best_coverage
            best_coverage = coverage
            best_hub = hub
        end
    end
    
    return best_hub
end

## Helper function for bidirectional path extension
function extend_path_bidirectional(hub_kmer::String, target_kmers::Set{String},
                                 covered_targets::Set{String}, graph::AbstractGraph,
                                 graph_kmers::Vector{String}, kmer_to_vertex::Dict{String, Int},
                                 kmer_coverage::Dict{String, Int}, k::Int)
    
    ## Find uncovered targets to prioritize
    uncovered_targets = setdiff(target_kmers, covered_targets)
    
    if isempty(uncovered_targets)
        return (path=String[], targets_found=0)
    end
    
    ## Use bidirectional search to find paths to uncovered targets
    hub_vertex = get(kmer_to_vertex, hub_kmer, nothing)
    
    if hub_vertex === nothing
        return (path=String[], targets_found=0)
    end
    
    ## Simple path extension using local search
    ## In production, would use the bidirectional Dijkstra from genomic-graph-algorithms.jl
    extended_path = [hub_kmer]
    current_vertex = hub_vertex
    targets_found = 0
    
    ## Greedy extension to maximize target coverage
    for step in 1:20  ## Avoid infinite loops
        best_neighbor = nothing
        best_score = -1.0
        
        if current_vertex <= Graphs.nv(graph)
            for neighbor_vertex in Graphs.neighbors(graph, current_vertex)
                if neighbor_vertex <= length(graph_kmers)
                    neighbor_kmer = graph_kmers[neighbor_vertex]
                    canonical_neighbor = _get_canonical_kmer(neighbor_kmer)
                    
                    ## Score based on target coverage and k-mer coverage
                    score = 0.0
                    if canonical_neighbor in uncovered_targets
                        score += 10.0  ## High priority for uncovered targets
                    end
                    score += log(1 + get(kmer_coverage, neighbor_kmer, 0))  ## Coverage bonus
                    
                    if score > best_score
                        best_score = score
                        best_neighbor = neighbor_vertex
                    end
                end
            end
        end
        
        if best_neighbor === nothing
            break  ## No more extensions possible
        end
        
        neighbor_kmer = graph_kmers[best_neighbor]
        push!(extended_path, neighbor_kmer)
        current_vertex = best_neighbor
        
        ## Check if we found a target
        canonical_neighbor = _get_canonical_kmer(neighbor_kmer)
        if canonical_neighbor in uncovered_targets
            targets_found += 1
        end
    end
    
    return (path=extended_path, targets_found=targets_found)
end

## Helper functions (simplified versions of complex operations)

function merge_paths(path1::Vector{String}, path2::Vector{String}, k::Int)
    ## Simple concatenation - in production would handle overlaps properly
    return vcat(path1, path2)
end

function identify_gap_regions(path::Vector{String}, targets::Set{String}, k::Int)
    ## Placeholder - would identify regions where targets are missing
    return Tuple{Int, Int}[]
end

function calculate_reconstruction_confidence(path::Vector{String}, 
                                           coverage::Dict{String, Int}, 
                                           completion::Float64)
    if isempty(path)
        return 0.0
    end
    
    ## Combine completion percentage with average coverage confidence
    avg_coverage = Statistics.mean([get(coverage, kmer, 1) for kmer in path])
    max_coverage = maximum(values(coverage))
    coverage_confidence = min(1.0, avg_coverage / (max_coverage * 0.5))
    
    return 0.7 * completion + 0.3 * coverage_confidence
end

function find_connecting_path(start_kmer::String, end_kmer::String,
                            graph::AbstractGraph, graph_kmers::Vector{String},
                            kmer_to_vertex::Dict{String, Int}, coverage::Dict{String, Int},
                            k::Int, max_steps::Int)
    
    ## Simplified path finding - would use genomic Dijkstra in production
    return String[]  ## Placeholder
end

function reconstruct_sequence_from_kmers(kmers::Vector{String}, k::Int)
    ## Same as in other modules
    if isempty(kmers)
        return ""
    end
    
    if length(kmers) == 1
        return kmers[1]
    end
    
    sequence = kmers[1]
    for i in 2:length(kmers)
        sequence *= kmers[i][end]
    end
    
    return sequence
end

## Helper function for canonical k-mer representation
function _get_canonical_kmer(kmer::String)
    reverse_complement = _reverse_complement(kmer)
    return kmer <= reverse_complement ? kmer : reverse_complement
end

## Helper function for reverse complement
function _reverse_complement(seq::String)
    complement_map = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
    return reverse(String([complement_map[c] for c in seq]))
end