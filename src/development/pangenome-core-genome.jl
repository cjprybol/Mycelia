# # Hub-based core genome identification algorithms
# # Based on analysis from Mycelia-Dev notebooks (2022-01-23-sample-core-genome.ipynb and variants)

# """
#     CoreGenomeResult

# Results from hub-based core genome identification analysis.

# # Fields
# - `core_kmers::Set{T}`: K-mers identified as core genomic elements
# - `hub_nodes::Vector{Int}`: Indices of hub nodes (degree ≥3) in the graph
# - `core_coverage::Float64`: Fraction of input genomes containing core elements
# - `core_length_estimate::Int`: Estimated total length of core genome
# - `traversal_paths::Vector{Vector{T}}`: Representative paths through core regions
# - `hub_statistics::Dict{String, Float64}`: Statistical measures of hub connectivity
# """
# struct CoreGenomeResult{T}
#     core_kmers::Set{T}
#     hub_nodes::Vector{Int}
#     core_coverage::Float64
#     core_length_estimate::Int
#     traversal_paths::Vector{Vector{T}}
#     hub_statistics::Dict{String, Float64}
# end

# """
#     PangenomeStats

# Statistical summary of pangenome analysis results.

# # Fields
# - `total_kmers::Int`: Total unique k-mers in pangenome
# - `core_kmers::Int`: Number of core k-mers
# - `accessory_kmers::Int`: Number of accessory k-mers  
# - `unique_kmers::Int`: Number of strain-specific k-mers
# - `core_fraction::Float64`: Fraction of pangenome that is core
# - `alpha_diversity::Float64`: Pangenome α-diversity (within-genome diversity)
# - `beta_diversity::Float64`: Pangenome β-diversity (between-genome diversity)
# """
# struct PangenomeStats
#     total_kmers::Int
#     core_kmers::Int
#     accessory_kmers::Int
#     unique_kmers::Int
#     core_fraction::Float64
#     alpha_diversity::Float64
#     beta_diversity::Float64
# end

# """
#     identify_core_genome_hubs(kmer_graph::AbstractGraph, 
#                              kmer_coverage::Dict{T, Int},
#                              min_hub_degree::Int=3,
#                              core_coverage_threshold::Float64=0.95) where T

# Identify core genomic regions using hub-based graph topology analysis.

# This algorithm identifies core genome elements by detecting highly connected "hub" nodes 
# in k-mer de Bruijn graphs that represent conserved genomic regions present across 
# multiple genomes or high-coverage regions within a single genome.

# # Arguments
# - `kmer_graph::AbstractGraph`: K-mer de Bruijn graph
# - `kmer_coverage::Dict{T, Int}`: Coverage counts for each k-mer
# - `min_hub_degree::Int=3`: Minimum node degree to qualify as a hub
# - `core_coverage_threshold::Float64=0.95`: Minimum coverage fraction for core classification

# # Returns
# `CoreGenomeResult{T}`: Comprehensive analysis of core genomic elements

# # Algorithm Details
# 1. **Hub Detection**: Identifies nodes with degree ≥ `min_hub_degree`
# 2. **Coverage Analysis**: Filters hubs by coverage threshold
# 3. **Path Traversal**: Traces paths between hubs to identify core regions
# 4. **Core Classification**: Classifies k-mers based on connectivity and coverage patterns

# # Examples
# ```julia
# # Build k-mer graph from multiple genomes
# graph, kmers = build_kmer_graph(genome_sequences, 31)
# coverage = count_kmer_coverage(genome_sequences, kmers)

# # Identify core genome using hub analysis
# result = identify_core_genome_hubs(graph, coverage)

# # Extract core k-mers for downstream analysis
# core_sequences = extract_core_sequences(result.core_kmers, kmers)
# ```

# # References
# Based on hub-based core genome identification from Mycelia-Dev notebooks:
# - `2022-01-23-sample-core-genome.ipynb`
# - `2022-01-26-sample-core-genome-circular.ipynb`
# - `2022-01-30-sample-core-genome.ipynb`
# """
# function identify_core_genome_hubs(kmer_graph::AbstractGraph, 
#                                   kmer_coverage::Dict{T, Int},
#                                   min_hub_degree::Int=3,
#                                   core_coverage_threshold::Float64=0.95) where T
    
#     ## Extract k-mers and create index mapping
#     kmers = collect(keys(kmer_coverage))
#     kmer_to_index = Dict(kmer => i for (i, kmer) in enumerate(kmers))
#     n_vertices = Graphs.nv(kmer_graph)
    
#     if length(kmers) != n_vertices
#         throw(ArgumentError("Number of k-mers ($(length(kmers))) must match graph vertices ($n_vertices)"))
#     end
    
#     ## Identify hub nodes based on degree
#     hub_nodes = Int[]
#     for v in Graphs.vertices(kmer_graph)
#         if Graphs.degree(kmer_graph, v) >= min_hub_degree
#             push!(hub_nodes, v)
#         end
#     end
    
#     @info "Identified $(length(hub_nodes)) hub nodes with degree ≥ $min_hub_degree"
    
#     ## Calculate coverage statistics
#     coverages = collect(values(kmer_coverage))
#     coverage_percentile_95 = Statistics.quantile(coverages, core_coverage_threshold)
    
#     ## Filter hubs by coverage threshold
#     high_coverage_hubs = Int[]
#     for hub in hub_nodes
#         if kmer_coverage[kmers[hub]] >= coverage_percentile_95
#             push!(high_coverage_hubs, hub)
#         end
#     end
    
#     @info "$(length(high_coverage_hubs)) hubs meet coverage threshold (≥ $(round(coverage_percentile_95, digits=2)))"
    
#     ## Identify core k-mers using hub-based traversal
#     core_kmers = Set{T}()
#     traversal_paths = Vector{T}[]
    
#     ## Add high-coverage hubs as core k-mers
#     for hub in high_coverage_hubs
#         push!(core_kmers, kmers[hub])
#     end
    
#     ## Trace paths between hubs to identify core regions
#     for i in 1:length(high_coverage_hubs)
#         for j in (i+1):length(high_coverage_hubs)
#             source_hub = high_coverage_hubs[i]
#             target_hub = high_coverage_hubs[j]
            
#             ## Find paths between hubs using BFS (limited depth to avoid long paths)
#             paths = _find_hub_connecting_paths(kmer_graph, source_hub, target_hub, 
#                                              kmers, kmer_coverage, coverage_percentile_95)
            
#             ## Add high-coverage paths to core genome
#             for path in paths
#                 if length(path) <= 10  ## Limit path length to avoid including accessory regions
#                     for kmer in path
#                         if kmer_coverage[kmer] >= coverage_percentile_95
#                             push!(core_kmers, kmer)
#                         end
#                     end
#                     push!(traversal_paths, path)
#                 end
#             end
#         end
#     end
    
#     ## Calculate core coverage fraction
#     total_coverage = sum(values(kmer_coverage))
#     core_coverage_sum = sum(kmer_coverage[kmer] for kmer in core_kmers)
#     core_coverage_fraction = core_coverage_sum / total_coverage
    
#     ## Estimate core genome length (accounting for k-mer overlaps)
#     k = length(string(first(kmers)))  ## Assume all k-mers have same length
#     core_length_estimate = length(core_kmers) > 0 ? length(core_kmers) + k - 1 : 0
    
#     ## Calculate hub connectivity statistics
#     hub_degrees = [Graphs.degree(kmer_graph, hub) for hub in hub_nodes]
#     hub_statistics = Dict{String, Float64}(
#         "mean_hub_degree" => isempty(hub_degrees) ? 0.0 : Statistics.mean(hub_degrees),
#         "max_hub_degree" => isempty(hub_degrees) ? 0.0 : maximum(hub_degrees),
#         "hub_degree_std" => isempty(hub_degrees) ? 0.0 : Statistics.std(hub_degrees),
#         "hub_fraction" => length(hub_nodes) / n_vertices,
#         "high_coverage_hub_fraction" => length(high_coverage_hubs) / max(1, length(hub_nodes))
#     )
    
#     @info "Core genome identification complete: $(length(core_kmers)) core k-mers " *
#           "($(round(100*core_coverage_fraction, digits=1))% coverage)"
    
#     return CoreGenomeResult{T}(
#         core_kmers,
#         hub_nodes,
#         core_coverage_fraction,
#         core_length_estimate,
#         traversal_paths,
#         hub_statistics
#     )
# end

# """
#     analyze_pangenome_structure(genome_kmer_sets::Vector{Set{T}}) where T

# Analyze pangenome structure to classify k-mers as core, accessory, or unique.

# # Arguments
# - `genome_kmer_sets::Vector{Set{T}}`: K-mer sets for each genome in the pangenome

# # Returns
# `PangenomeStats`: Statistical summary of pangenome composition and diversity

# # Algorithm Details
# - **Core k-mers**: Present in ≥95% of genomes
# - **Accessory k-mers**: Present in 2-94% of genomes  
# - **Unique k-mers**: Present in only 1 genome
# - **α-diversity**: Average within-genome k-mer diversity
# - **β-diversity**: Between-genome k-mer dissimilarity
# """
# function analyze_pangenome_structure(genome_kmer_sets::Vector{Set{T}}) where T
    
#     if isempty(genome_kmer_sets)
#         throw(ArgumentError("At least one genome k-mer set required"))
#     end
    
#     n_genomes = length(genome_kmer_sets)
    
#     ## Count k-mer frequencies across genomes
#     kmer_counts = Dict{T, Int}()
#     for kmer_set in genome_kmer_sets
#         for kmer in kmer_set
#             kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
#         end
#     end
    
#     total_kmers = length(kmer_counts)
    
#     ## Classify k-mers by presence patterns
#     core_threshold = ceil(Int, 0.95 * n_genomes)  ## Present in ≥95% of genomes
    
#     core_kmers = 0
#     accessory_kmers = 0
#     unique_kmers = 0
    
#     for (kmer, count) in kmer_counts
#         if count >= core_threshold
#             core_kmers += 1
#         elseif count == 1
#             unique_kmers += 1
#         else
#             accessory_kmers += 1
#         end
#     end
    
#     core_fraction = core_kmers / total_kmers
    
#     ## Calculate α-diversity (within-genome diversity)
#     genome_sizes = [length(kmer_set) for kmer_set in genome_kmer_sets]
#     alpha_diversity = Statistics.mean(genome_sizes)
    
#     ## Calculate β-diversity (between-genome dissimilarity using Jaccard distance)
#     if n_genomes >= 2
#         jaccard_distances = Float64[]
        
#         for i in 1:(n_genomes-1)
#             for j in (i+1):n_genomes
#                 set1 = genome_kmer_sets[i]
#                 set2 = genome_kmer_sets[j]
                
#                 intersection_size = length(intersect(set1, set2))
#                 union_size = length(union(set1, set2))
                
#                 jaccard_distance = union_size > 0 ? 1.0 - (intersection_size / union_size) : 1.0
#                 push!(jaccard_distances, jaccard_distance)
#             end
#         end
        
#         beta_diversity = Statistics.mean(jaccard_distances)
#     else
#         beta_diversity = 0.0
#     end
    
#     return PangenomeStats(
#         total_kmers,
#         core_kmers,
#         accessory_kmers,
#         unique_kmers,
#         core_fraction,
#         alpha_diversity,
#         beta_diversity
#     )
# end

# """
#     extract_core_sequences(core_kmers::Set{T}, k::Int) where T

# Extract contiguous sequences from core k-mers using overlap-based assembly.

# # Arguments
# - `core_kmers::Set{T}`: Set of core k-mers identified by hub analysis
# - `k::Int`: K-mer size used in analysis

# # Returns
# `Vector{String}`: Assembled core genome sequences

# # Algorithm Details
# Uses greedy overlap extension to assemble core k-mers into longer sequences:
# 1. **Graph Construction**: Build overlap graph from core k-mers
# 2. **Path Finding**: Identify maximal non-branching paths
# 3. **Sequence Assembly**: Extend k-mers along paths to reconstruct sequences
# """
# function extract_core_sequences(core_kmers::Set{T}, k::Int) where T
    
#     if isempty(core_kmers)
#         return String[]
#     end
    
#     ## Convert k-mers to strings for processing
#     kmer_strings = [string(kmer) for kmer in core_kmers]
#     kmer_set = Set(kmer_strings)
    
#     ## Build overlap graph
#     overlap_graph = Dict{String, Vector{String}}()
    
#     for kmer in kmer_strings
#         overlap_graph[kmer] = String[]
        
#         ## Find overlapping k-mers (k-1 suffix matches k-1 prefix)
#         suffix = kmer[2:end]
        
#         for base in ['A', 'C', 'G', 'T']
#             candidate = suffix * base
#             if candidate in kmer_set && candidate != kmer
#                 push!(overlap_graph[kmer], candidate)
#             end
#         end
#     end
    
#     ## Find maximal non-branching paths
#     sequences = String[]
#     visited = Set{String}()
    
#     for start_kmer in kmer_strings
#         if start_kmer in visited
#             continue
#         end
        
#         ## Check if this is a potential start of a path (in-degree 0 or 1)
#         in_degree = sum(1 for (_, neighbors) in overlap_graph if start_kmer in neighbors)
        
#         if in_degree <= 1
#             path = _extend_core_path(start_kmer, overlap_graph, visited)
            
#             if length(path) >= 2  ## Only keep paths with multiple k-mers
#                 ## Assemble sequence from k-mer path
#                 sequence = path[1]  ## Start with first k-mer
#                 for i in 2:length(path)
#                     sequence *= path[i][end]  ## Add last nucleotide of each subsequent k-mer
#                 end
#                 push!(sequences, sequence)
#             end
#         end
#     end
    
#     ## Sort sequences by length (longest first)
#     sort!(sequences, by=length, rev=true)
    
#     @info "Assembled $(length(sequences)) core genome sequences from $(length(core_kmers)) core k-mers"
    
#     return sequences
# end

# ## Helper function to find paths between hubs with coverage filtering
# function _find_hub_connecting_paths(graph::AbstractGraph, source::Int, target::Int,
#                                    kmers::Vector{T}, kmer_coverage::Dict{T, Int},
#                                    min_coverage::Float64) where T
    
#     paths = Vector{T}[]
    
#     ## Use BFS to find paths with limited depth
#     max_depth = 5
#     queue = [(source, [kmers[source]], 0)]  ## (current_node, path, depth)
#     visited = Set{Int}()
    
#     while !isempty(queue)
#         current, path, depth = popfirst!(queue)
        
#         if current == target
#             push!(paths, path)
#             continue
#         end
        
#         if depth >= max_depth || current in visited
#             continue
#         end
        
#         push!(visited, current)
        
#         ## Explore neighbors with sufficient coverage
#         for neighbor in Graphs.neighbors(graph, current)
#             if kmer_coverage[kmers[neighbor]] >= min_coverage
#                 new_path = copy(path)
#                 push!(new_path, kmers[neighbor])
#                 push!(queue, (neighbor, new_path, depth + 1))
#             end
#         end
#     end
    
#     return paths
# end

# ## Helper function to extend core genome paths using greedy assembly
# function _extend_core_path(start_kmer::String, overlap_graph::Dict{String, Vector{String}}, 
#                           visited::Set{String})
    
#     path = [start_kmer]
#     push!(visited, start_kmer)
#     current = start_kmer
    
#     ## Extend forward as far as possible
#     while true
#         neighbors = overlap_graph[current]
#         unvisited_neighbors = [n for n in neighbors if !(n in visited)]
        
#         if length(unvisited_neighbors) == 1
#             ## Unique extension - continue path
#             next_kmer = unvisited_neighbors[1]
#             push!(path, next_kmer)
#             push!(visited, next_kmer)
#             current = next_kmer
#         else
#             ## Multiple or zero extensions - stop path
#             break
#         end
#     end
    
#     return path
# end