# # Protein clustering pipeline for functional pangenome analysis
# # Based on analysis from Mycelia-Dev notebook (2022-03-12-ecoli-tequatrovirus.ipynb)

# """
#     ProteinClusterResult

# Results from protein clustering analysis for pangenome studies.

# # Fields
# - `cluster_assignments::Dict{String, Int}`: Protein ID to cluster assignment mapping
# - `cluster_representatives::Dict{Int, String}`: Representative protein for each cluster
# - `cluster_sizes::Vector{Int}`: Number of proteins in each cluster
# - `functional_annotations::Dict{Int, String}`: Functional annotation for each cluster
# - `core_clusters::Set{Int}`: Cluster IDs identified as core functions
# - `accessory_clusters::Set{Int}`: Cluster IDs identified as accessory functions
# - `unique_clusters::Set{Int}`: Cluster IDs identified as unique functions
# """
# struct ProteinClusterResult
#     cluster_assignments::Dict{String, Int}
#     cluster_representatives::Dict{Int, String}
#     cluster_sizes::Vector{Int}
#     functional_annotations::Dict{Int, String}
#     core_clusters::Set{Int}
#     accessory_clusters::Set{Int}
#     unique_clusters::Set{Int}
# end

# """
#     FunctionalPangenomeStats

# Statistical summary of functional pangenome composition.

# # Fields
# - `total_proteins::Int`: Total number of proteins analyzed
# - `total_clusters::Int`: Total number of protein clusters identified
# - `core_functions::Int`: Number of core functional clusters
# - `accessory_functions::Int`: Number of accessory functional clusters
# - `unique_functions::Int`: Number of unique functional clusters
# - `core_function_fraction::Float64`: Fraction of clusters that are core
# - `functional_diversity::Float64`: Shannon diversity of functional categories
# - `pangenome_openness::Float64`: Estimated openness of the pangenome
# """
# struct FunctionalPangenomeStats
#     total_proteins::Int
#     total_clusters::Int
#     core_functions::Int
#     accessory_functions::Int
#     unique_functions::Int
#     core_function_fraction::Float64
#     functional_diversity::Float64
#     pangenome_openness::Float64
# end

# """
#     cluster_proteins_for_pangenome(protein_sequences::Dict{String, Vector{String}},
#                                   similarity_threshold::Float64=0.70;
#                                   clustering_method::Symbol=:mmseqs2,
#                                   annotation_source::Symbol=:pfam) -> ProteinClusterResult

# Cluster proteins across multiple genomes for functional pangenome analysis.

# This pipeline performs sequence-based clustering of proteins and classifies the resulting
# clusters as core, accessory, or unique based on their distribution across genomes.

# # Arguments
# - `protein_sequences::Dict{String, Vector{String}}`: Genome ID to protein sequences mapping
# - `similarity_threshold::Float64=0.70`: Sequence similarity threshold for clustering
# - `clustering_method::Symbol=:mmseqs2`: Clustering algorithm (:mmseqs2, :cd_hit, :diamond)
# - `annotation_source::Symbol=:pfam`: Functional annotation source (:pfam, :cog, :kegg)

# # Returns
# `ProteinClusterResult`: Comprehensive clustering results with functional classification

# # Algorithm Details
# 1. **Sequence Clustering**: Groups proteins by sequence similarity using specified threshold
# 2. **Functional Classification**: Assigns clusters as core (≥95% genomes), accessory (2-94%), or unique (1 genome)
# 3. **Representative Selection**: Selects longest sequence as cluster representative
# 4. **Annotation Integration**: Maps functional annotations to protein clusters

# # Clustering Methods
# - **MMseqs2**: Fast and sensitive sequence clustering (recommended)
# - **CD-HIT**: Classical clustering with adjustable word size
# - **DIAMOND**: Fast approximate clustering for large datasets

# # Core/Accessory Classification
# - **Core Functions**: Present in ≥95% of genomes (essential cellular functions)
# - **Accessory Functions**: Present in 2-94% of genomes (adaptive functions)
# - **Unique Functions**: Present in only 1 genome (strain-specific functions)

# # Examples
# ```julia
# # Basic protein clustering
# genome_proteins = Dict(
#     "genome1" => load_protein_sequences("genome1.faa"),
#     "genome2" => load_protein_sequences("genome2.faa")
# )
# result = cluster_proteins_for_pangenome(genome_proteins)

# # High-stringency clustering with KEGG annotations
# result = cluster_proteins_for_pangenome(genome_proteins, 0.90,
#                                        clustering_method=:mmseqs2,
#                                        annotation_source=:kegg)

# # Analyze core vs accessory functions
# core_fraction = length(result.core_clusters) / length(result.cluster_assignments)
# ```

# # References
# Based on protein clustering pipeline from:
# `2022-03-12-ecoli-tequatrovirus.ipynb`
# """
# function cluster_proteins_for_pangenome(protein_sequences::Dict{String, Vector{String}},
#                                        similarity_threshold::Float64=0.70;
#                                        clustering_method::Symbol=:mmseqs2,
#                                        annotation_source::Symbol=:pfam)

#     if isempty(protein_sequences)
#         throw(ArgumentError("At least one genome with protein sequences required"))
#     end

#     genome_ids = collect(keys(protein_sequences))
#     n_genomes = length(genome_ids)

#     @info "Clustering proteins from $n_genomes genomes using $clustering_method (threshold: $similarity_threshold)"

#     ## Combine all proteins with genome tracking
#     all_proteins = Dict{String, String}()  ## protein_id => sequence
#     protein_to_genome = Dict{String, String}()  ## protein_id => genome_id

#     protein_counter = 0
#     for (genome_id, sequences) in protein_sequences
#         for (i, seq) in enumerate(sequences)
#             protein_counter += 1
#             protein_id = "$(genome_id)_protein_$(i)"
#             all_proteins[protein_id] = seq
#             protein_to_genome[protein_id] = genome_id
#         end
#     end

#     @info "Total proteins to cluster: $(length(all_proteins))"

#     ## Perform sequence clustering
#     cluster_assignments = perform_sequence_clustering(all_proteins, similarity_threshold, clustering_method)

#     ## Analyze cluster distribution across genomes
#     cluster_genome_distribution = analyze_cluster_distribution(cluster_assignments, protein_to_genome, genome_ids)

#     ## Classify clusters as core/accessory/unique
#     core_threshold = ceil(Int, 0.95 * n_genomes)
#     core_clusters = Set{Int}()
#     accessory_clusters = Set{Int}()
#     unique_clusters = Set{Int}()

#     for (cluster_id, genomes_present) in cluster_genome_distribution
#         n_genomes_present = length(genomes_present)

#         if n_genomes_present >= core_threshold
#             push!(core_clusters, cluster_id)
#         elseif n_genomes_present == 1
#             push!(unique_clusters, cluster_id)
#         else
#             push!(accessory_clusters, cluster_id)
#         end
#     end

#     @info "Cluster classification: $(length(core_clusters)) core, " *
#           "$(length(accessory_clusters)) accessory, $(length(unique_clusters)) unique"

#     ## Select cluster representatives (longest sequence in each cluster)
#     cluster_representatives = select_cluster_representatives(cluster_assignments, all_proteins)

#     ## Calculate cluster sizes
#     cluster_sizes = calculate_cluster_sizes(cluster_assignments)

#     ## Annotate clusters functionally
#     functional_annotations = annotate_protein_clusters(cluster_representatives, annotation_source)

#     return ProteinClusterResult(
#         cluster_assignments,
#         cluster_representatives,
#         cluster_sizes,
#         functional_annotations,
#         core_clusters,
#         accessory_clusters,
#         unique_clusters
#     )
# end

# """
#     analyze_functional_pangenome(results::Vector{ProteinClusterResult}) -> FunctionalPangenomeStats

# Analyze functional pangenome composition and estimate openness from clustering results.

# # Arguments
# - `results::Vector{ProteinClusterResult}`: Results from incremental pangenome sampling

# # Returns
# `FunctionalPangenomeStats`: Comprehensive functional pangenome statistics
# """
# function analyze_functional_pangenome(results::Vector{ProteinClusterResult})

#     if isempty(results)
#         throw(ArgumentError("At least one clustering result required"))
#     end

#     final_result = results[end]  ## Use final (complete) clustering result

#     total_proteins = sum(values(length.(values(final_result.cluster_assignments))))
#     total_clusters = length(final_result.cluster_representatives)
#     core_functions = length(final_result.core_clusters)
#     accessory_functions = length(final_result.accessory_clusters)
#     unique_functions = length(final_result.unique_clusters)

#     core_function_fraction = core_functions / total_clusters

#     ## Calculate functional diversity using Shannon entropy
#     cluster_sizes = final_result.cluster_sizes
#     total_size = sum(cluster_sizes)
#     proportions = cluster_sizes ./ total_size
#     functional_diversity = -sum(p * log(p) for p in proportions if p > 0)

#     ## Estimate pangenome openness using power law fitting
#     if length(results) >= 3
#         pangenome_openness = estimate_pangenome_openness(results)
#     else
#         pangenome_openness = NaN
#     end

#     return FunctionalPangenomeStats(
#         total_proteins,
#         total_clusters,
#         core_functions,
#         accessory_functions,
#         unique_functions,
#         core_function_fraction,
#         functional_diversity,
#         pangenome_openness
#     )
# end

# """
#     extract_core_proteome(result::ProteinClusterResult) -> Dict{String, String}

# Extract core proteome sequences from clustering results.

# # Arguments
# - `result::ProteinClusterResult`: Protein clustering results

# # Returns
# `Dict{String, String}`: Core cluster ID to representative sequence mapping
# """
# function extract_core_proteome(result::ProteinClusterResult)

#     core_proteome = Dict{String, String}()

#     for cluster_id in result.core_clusters
#         if haskey(result.cluster_representatives, cluster_id)
#             representative_id = result.cluster_representatives[cluster_id]
#             cluster_name = "core_cluster_$(cluster_id)"

#             ## Add functional annotation if available
#             if haskey(result.functional_annotations, cluster_id)
#                 annotation = result.functional_annotations[cluster_id]
#                 cluster_name *= "_$(annotation)"
#             end

#             ## Note: Would need access to original sequences to include actual sequence
#             ## This is a placeholder structure
#             core_proteome[cluster_name] = representative_id
#         end
#     end

#     @info "Extracted $(length(core_proteome)) core protein families"

#     return core_proteome
# end

# ## Helper function for sequence clustering
# function perform_sequence_clustering(protein_sequences::Dict{String, String}, 
#                                    similarity_threshold::Float64,
#                                    method::Symbol)

#     protein_ids = collect(keys(protein_sequences))
#     n_proteins = length(protein_ids)

#     ## For demonstration, implement simple similarity-based clustering
#     ## In production, this would interface with MMseqs2, CD-HIT, or DIAMOND

#     cluster_assignments = Dict{String, Int}()
#     clusters = Vector{Vector{String}}()

#     @info "Performing $method clustering of $n_proteins proteins"

#     ## Simple greedy clustering algorithm
#     unassigned = Set(protein_ids)
#     cluster_id = 0

#     while !isempty(unassigned)
#         cluster_id += 1

#         ## Start new cluster with first unassigned protein
#         seed_protein = first(unassigned)
#         current_cluster = [seed_protein]
#         cluster_assignments[seed_protein] = cluster_id
#         delete!(unassigned, seed_protein)

#         ## Find similar proteins to add to this cluster
#         seed_sequence = protein_sequences[seed_protein]

#         to_remove = String[]
#         for protein_id in unassigned
#             similarity = calculate_sequence_similarity(seed_sequence, protein_sequences[protein_id])

#             if similarity >= similarity_threshold
#                 push!(current_cluster, protein_id)
#                 cluster_assignments[protein_id] = cluster_id
#                 push!(to_remove, protein_id)
#             end
#         end

#         ## Remove assigned proteins from unassigned set
#         for protein_id in to_remove
#             delete!(unassigned, protein_id)
#         end

#         push!(clusters, current_cluster)
#     end

#     @info "Created $cluster_id clusters from $n_proteins proteins"

#     return cluster_assignments
# end

# ## Helper function to analyze cluster distribution across genomes
# function analyze_cluster_distribution(cluster_assignments::Dict{String, Int},
#                                     protein_to_genome::Dict{String, String},
#                                     genome_ids::Vector{String})

#     cluster_genome_distribution = Dict{Int, Set{String}}()

#     for (protein_id, cluster_id) in cluster_assignments
#         genome_id = protein_to_genome[protein_id]

#         if !haskey(cluster_genome_distribution, cluster_id)
#             cluster_genome_distribution[cluster_id] = Set{String}()
#         end

#         push!(cluster_genome_distribution[cluster_id], genome_id)
#     end

#     return cluster_genome_distribution
# end

# ## Helper function to select cluster representatives
# function select_cluster_representatives(cluster_assignments::Dict{String, Int},
#                                       protein_sequences::Dict{String, String})

#     ## Group proteins by cluster
#     clusters = Dict{Int, Vector{String}}()
#     for (protein_id, cluster_id) in cluster_assignments
#         if !haskey(clusters, cluster_id)
#             clusters[cluster_id] = String[]
#         end
#         push!(clusters[cluster_id], protein_id)
#     end

#     ## Select longest sequence as representative
#     representatives = Dict{Int, String}()
#     for (cluster_id, protein_ids) in clusters
#         longest_protein = ""
#         longest_length = 0

#         for protein_id in protein_ids
#             seq_length = length(protein_sequences[protein_id])
#             if seq_length > longest_length
#                 longest_length = seq_length
#                 longest_protein = protein_id
#             end
#         end

#         representatives[cluster_id] = longest_protein
#     end

#     return representatives
# end

# ## Helper function to calculate cluster sizes
# function calculate_cluster_sizes(cluster_assignments::Dict{String, Int})

#     cluster_counts = Dict{Int, Int}()
#     for (protein_id, cluster_id) in cluster_assignments
#         cluster_counts[cluster_id] = get(cluster_counts, cluster_id, 0) + 1
#     end

#     max_cluster_id = maximum(keys(cluster_counts))
#     sizes = zeros(Int, max_cluster_id)

#     for (cluster_id, count) in cluster_counts
#         sizes[cluster_id] = count
#     end

#     return sizes
# end

# ## Helper function for functional annotation
# function annotate_protein_clusters(cluster_representatives::Dict{Int, String},
#                                   annotation_source::Symbol)

#     annotations = Dict{Int, String}()

#     ## Placeholder functional annotation
#     ## In production, this would interface with Pfam, COG, KEGG databases

#     for (cluster_id, representative_id) in cluster_representatives
#         ## Mock annotation based on cluster size and ID
#         if cluster_id <= 100
#             annotations[cluster_id] = "core_function_family_$(cluster_id)"
#         elseif cluster_id <= 500
#             annotations[cluster_id] = "accessory_function_family_$(cluster_id)"
#         else
#             annotations[cluster_id] = "hypothetical_protein_$(cluster_id)"
#         end
#     end

#     @info "Annotated $(length(annotations)) protein clusters using $annotation_source"

#     return annotations
# end

# ## Helper function for simple sequence similarity calculation
# function calculate_sequence_similarity(seq1::String, seq2::String)

#     ## Simple identity-based similarity for demonstration
#     ## In production, would use BLAST, DIAMOND, or other alignment tools

#     if length(seq1) == 0 || length(seq2) == 0
#         return 0.0
#     end

#     ## For demonstration, use Hamming distance for equal-length sequences
#     if length(seq1) == length(seq2)
#         matches = sum(seq1[i] == seq2[i] for i in 1:length(seq1))
#         return matches / length(seq1)
#     else
#         ## For different lengths, use a simple heuristic
#         shorter = min(length(seq1), length(seq2))
#         longer = max(length(seq1), length(seq2))

#         ## Penalize length differences
#         length_penalty = (longer - shorter) / longer

#         ## Calculate identity over shorter sequence
#         matches = sum(seq1[i] == seq2[i] for i in 1:shorter)
#         identity = matches / shorter

#         return identity * (1.0 - length_penalty)
#     end
# end

# ## Helper function to estimate pangenome openness
# function estimate_pangenome_openness(results::Vector{ProteinClusterResult})

#     ## Extract cluster counts as function of genome number
#     cluster_counts = [length(result.cluster_representatives) for result in results]
#     genome_numbers = collect(1:length(cluster_counts))

#     ## Fit power law: clusters = a * genomes^b
#     ## Openness parameter is b (b > 0 indicates open pangenome)

#     if length(cluster_counts) < 3
#         return NaN
#     end

#     ## Simple linear regression in log space: log(clusters) = log(a) + b * log(genomes)
#     log_genomes = log.(genome_numbers)
#     log_clusters = log.(cluster_counts)

#     ## Calculate regression coefficients
#     n = length(log_genomes)
#     sum_x = sum(log_genomes)
#     sum_y = sum(log_clusters)
#     sum_xy = sum(log_genomes .* log_clusters)
#     sum_x2 = sum(log_genomes .^ 2)

#     ## Slope coefficient (openness parameter)
#     b = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x^2)

#     return b
# end
