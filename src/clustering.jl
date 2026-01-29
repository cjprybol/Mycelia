"""
Cluster sequences based on pairwise identity threshold and return representative sequences.

This function treats the problem as a graph clustering task where:
- Each unique sequence is a node
- Edges exist between sequences with identity ≥ threshold
- Connected components represent clusters of similar sequences
- The first occurrence (by original order) in each cluster is retained as representative
"""
function deduplicate_sequences(df::DataFrames.DataFrame, threshold::Float64=99.5)
    # Get all unique sequences and create mapping
    all_sequences = unique([df.query; df.reference])
    seq_to_index = Dict(seq => i for (i, seq) in enumerate(all_sequences))
    n_sequences = length(all_sequences)
    
    # Create undirected graph
    graph = Graphs.SimpleGraph(n_sequences)
    
    # Add edges for sequence pairs meeting threshold
    for row in DataFrames.eachrow(df)
        if row."%_identity" >= threshold
            query_idx = seq_to_index[row.query]
            ref_idx = seq_to_index[row.reference]
            Graphs.add_edge!(graph, query_idx, ref_idx)
        end
    end
    
    # Find connected components (clusters)
    components = Graphs.connected_components(graph)
    
    # For each cluster, find the representative (earliest appearing sequence)
    representatives = String[]
    cluster_info = DataFrames.DataFrame(
        representative = String[],
        cluster_size = Int[],
        cluster_members = Vector{String}[]
    )
    
    for component in components
        # Get sequences in this cluster
        cluster_seqs = all_sequences[component]
        
        # Find the representative (first occurrence in original data)
        # Check both query and reference columns for first appearance
        first_positions = Int[]
        for seq in cluster_seqs
            query_positions = findall(x -> x == seq, df.query)
            ref_positions = findall(x -> x == seq, df.reference)
            all_positions = [query_positions; ref_positions]
            if !isempty(all_positions)
                push!(first_positions, minimum(all_positions))
            else
                push!(first_positions, typemax(Int))  # Should not happen
            end
        end
        
        # Representative is sequence with minimum first position
        min_pos_idx = argmin(first_positions)
        representative = cluster_seqs[min_pos_idx]
        
        push!(representatives, representative)
        push!(cluster_info, (
            representative = representative,
            cluster_size = length(cluster_seqs),
            cluster_members = cluster_seqs
        ))
    end
    
    return (;representatives, cluster_info)
end

"""
Filter original dataframe to only include rows involving representative sequences.
"""
function filter_to_representatives(df::DataFrames.DataFrame, representatives::Vector{String})
    rep_set = Set(representatives)
    mask = [row.query in rep_set && row.reference in rep_set for row in DataFrames.eachrow(df)]
    return df[mask, :]
end

"""
Generate summary statistics about the clustering results.
"""
function clustering_summary(original_df::DataFrames.DataFrame, cluster_info::DataFrames.DataFrame)
    n_original_sequences = length(unique([original_df.query; original_df.reference]))
    n_clusters = nrow(cluster_info)
    n_representatives = n_clusters
    
    cluster_sizes = cluster_info.cluster_size
    
    println("=== Sequence Deduplication Summary ===")
    println("Original sequences: $n_original_sequences")
    println("Representative sequences: $n_representatives")
    println("Sequences removed: $(n_original_sequences - n_representatives)")
    println("Reduction: $(round((1 - n_representatives/n_original_sequences) * 100, digits=2))%")
    println()
    println("Cluster size distribution:")
    size_counts = StatsBase.countmap(cluster_sizes)
    for size in sort(collect(keys(size_counts)))
        count = size_counts[size]
        println("  Size $size: $count clusters")
    end
    
    return (
        original_count = n_original_sequences,
        representative_count = n_representatives,
        reduction_percent = (1 - n_representatives/n_original_sequences) * 100
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Cluster protein or nucleotide sequences using MMseqs2 easy-cluster workflow.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequences to cluster
- `output::String`: Base path for output files (default: input path + ".mmseqs_easy_cluster")
- `tmp::String`: Path to temporary directory (default: auto-generated temp dir)

# Returns
- `String`: Path to the output cluster TSV file containing cluster assignments

# Details
Uses MMseqs2 with minimum sequence identity threshold of 50% (-min-seq-id 0.5) and 
minimum coverage threshold of 80% (-c 0.8). The output TSV file format contains 
tab-separated cluster representative and member sequences.
"""
# --cov-mode: coverage mode (0: coverage of query and target, 1: coverage of target, 2: coverage of query)
function mmseqs_easy_cluster(;fasta, output=fasta*".mmseqs_easy_cluster", tmp=mktempdir())
    outfile = "$(output)_cluster.tsv"
    if !isfile(outfile)
        Mycelia.add_bioconda_env("mmseqs2")
        # at least 50% equivalent
        # --min-seq-id 0.5 -c 0.8
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-cluster $(fasta) $(output) $(tmp)`)
    end
    rm(tmp, recursive=true)
    return "$(output)_cluster.tsv"
end

"""
Vertex data for hierarchical clustering graph nodes.
"""
struct HclustVertexData
    hclust_id::String
    height::Float64
    x::Float64
    y::Float64
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a hierarchical clustering tree into a directed metagraph representation.

# Arguments
- `hcl::Clustering.Hclust`: Hierarchical clustering result object

# Returns
- `MetaGraphsNext.MetaGraph`: Directed graph with metadata representing the clustering hierarchy

# Graph Properties
The resulting graph contains vertices with the following data:
- `hclust_id`: String identifier for each node
- `height`: Height/distance at each merge point (0.0 for leaves)
- `x`: Horizontal position for visualization (0-1 range)
- `y`: Vertical position based on normalized height
"""
function hclust_to_metagraph(hcl::Clustering.Hclust)
    # Create MetaGraphsNext graph
    mg = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph(),
        label_type=String,
        vertex_data_type=HclustVertexData,
        edge_data_type=Nothing,
        graph_data=hcl  # Store original clustering object as graph metadata
    )

    # Add leaf nodes
    for leaf_node in 1:Clustering.nnodes(hcl)
        hclust_id = string(-leaf_node)
        mg[hclust_id] = HclustVertexData(hclust_id, 0.0, 0.0, 0.0)  # x and y will be set later
    end

    # Set x positions for leaf nodes based on order
    for (i, ordered_leaf_node) in enumerate(hcl.order)
        hclust_id = string(-ordered_leaf_node)
        x = i / (length(hcl.order) + 1)
        vertex_data = mg[hclust_id]
        mg[hclust_id] = HclustVertexData(vertex_data.hclust_id, vertex_data.height, x, vertex_data.y)
    end

    # Add internal nodes
    for (i, (left, right)) in enumerate(eachrow(hcl.merges))
        hclust_id = string(i)
        mg[hclust_id] = HclustVertexData(hclust_id, 0.0, 0.0, 0.0)  # height, x, y will be set later
    end

    # Set up hierarchy with edges and positions
    for (i, ((left, right), height)) in enumerate(zip(eachrow(hcl.merges), hcl.heights))
        parent_id = string(i)
        left_child_id = string(left)
        right_child_id = string(right)

        # Update parent height
        parent_data = mg[parent_id]
        mg[parent_id] = HclustVertexData(parent_data.hclust_id, height, parent_data.x, parent_data.y)

        # Calculate parent x position as midpoint of children
        left_child_data = mg[left_child_id]
        right_child_data = mg[right_child_id]
        x = (left_child_data.x + right_child_data.x) / 2

        # Update parent x position
        parent_data = mg[parent_id]
        mg[parent_id] = HclustVertexData(parent_data.hclust_id, parent_data.height, x, parent_data.y)

        # Add edges from parent to children
        mg[parent_id, left_child_id] = nothing
        mg[parent_id, right_child_id] = nothing
    end

    # Calculate and set y positions based on normalized heights
    max_height = maximum(hcl.heights)
    for label in MetaGraphsNext.labels(mg)
        vertex_data = mg[label]
        relative_height = (max_height - vertex_data.height) / max_height
        mg[label] = HclustVertexData(vertex_data.hclust_id, vertex_data.height, vertex_data.x, relative_height)
    end

    return mg
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Performs hierarchical clustering on a distance matrix using Ward's linkage method.

# Arguments
- `distance_matrix::Matrix{<:Real}`: A symmetric distance/dissimilarity matrix

# Returns
- `HierarchicalCluster`: A hierarchical clustering object from Clustering.jl

# Details
Uses Ward's method (minimum variance) for clustering, which:
- Minimizes total within-cluster variance
- Produces compact, spherical clusters
- Works well for visualization in radial layouts
"""
function heirarchically_cluster_distance_matrix(distance_matrix)
    # Perform hierarchical clustering using Ward's method.
    # Ward's method minimizes the total within-cluster variance, which tends to
    # produce compact, spherical clusters. This can be beneficial for visualization,
    # especially in radial layouts, as it can prevent branches from being overly long
    # and overlapping.  Other linkage methods like 'average' and 'complete' are also
    # common choices, however 'single' linkage is prone to chaining effects.
    hcl = Clustering.hclust(distance_matrix, linkage=:ward, branchorder=:optimal)
    # this is equivalent to UPGMA
    # hcl = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)
    return hcl
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a CairoMakie plot of silhouette scores vs number of clusters.

Args:
    ks: Vector or range of k values tested.
    avg_silhouette_scores: Vector of average silhouette scores corresponding to ks.
    optimal_k: The optimal number of clusters identified.

Returns:
    A CairoMakie Figure object.
"""
function plot_silhouette_scores_makie(ks, avg_silhouette_scores, optimal_k)
    # Filter valid scores for y-axis limits
    valid_indices = findall(s -> !isnan(s) && s != Float64(-Inf), avg_silhouette_scores)
    if isempty(valid_indices)
        error("No valid silhouette scores to plot.")
    end
    
    valid_scores = avg_silhouette_scores[valid_indices]
    min_val = minimum(valid_scores)
    max_val = maximum(valid_scores)
    range_val = max_val - min_val
    
    # Set sensible y-axis limits
    ylims_lower = max(-1.05, min(0, min_val - range_val * 0.1))
    ylims_upper = min(1.05, max_val + range_val * 0.1)
    
    if isapprox(ylims_lower, ylims_upper)
        ylims_lower -= 0.1
        ylims_upper += 0.1
        ylims_lower = max(-1.05, ylims_lower)
        ylims_upper = min(1.05, ylims_upper)
    end
    
    # Create figure with higher resolution
    # fig = CairoMakie.Figure(resolution = (1200, 800), fontsize = 16)
    fig = CairoMakie.Figure(size = (900, 600), fontsize = 16)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Number of Clusters (k)",
        ylabel = "Average Silhouette Score",
        title = "Clustering Performance vs. Number of Clusters\n(higher is better)",
        titlesize = 20
    )
    
    # Plot scatter points
    CairoMakie.scatter!(
        ax,
        collect(ks),
        avg_silhouette_scores,
        color = :steelblue,
        markersize = 12,
        strokewidth = 1.5,
        strokecolor = :black,
        label = "Avg Score"
    )
    
    # Add line connecting points for better visualization
    CairoMakie.lines!(
        ax,
        collect(ks),
        avg_silhouette_scores,
        color = (:steelblue, 0.3),
        linewidth = 2
    )
    
    # Add vertical line at optimal k
    CairoMakie.vlines!(
        ax,
        [optimal_k],
        color = :red,
        linestyle = :dash,
        linewidth = 2.5,
        label = "Inferred optimum = $(optimal_k)"
    )
    
    # Add horizontal reference line at y=0
    CairoMakie.hlines!(
        ax,
        [0],
        color = :gray,
        linestyle = :dot,
        linewidth = 1.5,
        label = "Score = 0"
    )
    
    # Set axis limits
    CairoMakie.ylims!(ax, ylims_lower, ylims_upper)
    
    # Add legend
    CairoMakie.axislegend(ax, position = :rt)
    
    return fig
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a StatsPlots visualization of silhouette scores (legacy version).

Args:
    ks: Vector or range of k values tested.
    avg_silhouette_scores: Vector of average silhouette scores corresponding to ks.
    optimal_k: The optimal number of clusters identified.

Returns:
    A StatsPlots plot object.
"""
function plot_silhouette_scores_statsplots(ks, avg_silhouette_scores, optimal_k)
    # Filter valid scores for y-axis limits
    valid_indices = findall(s -> !isnan(s) && s != Float64(-Inf), avg_silhouette_scores)
    if isempty(valid_indices)
        error("No valid silhouette scores to plot.")
    end
    
    valid_scores = avg_silhouette_scores[valid_indices]
    min_val = minimum(valid_scores)
    max_val = maximum(valid_scores)
    range_val = max_val - min_val
    
    ylims_lower = max(-1.05, min(0, min_val - range_val * 0.1))
    ylims_upper = min(1.05, max_val + range_val * 0.1)
    
    if isapprox(ylims_lower, ylims_upper)
        ylims_lower -= 0.1
        ylims_upper += 0.1
        ylims_lower = max(-1.05, ylims_lower)
        ylims_upper = min(1.05, ylims_upper)
    end
    
    p = StatsPlots.scatter(
        ks,
        avg_silhouette_scores,
        title = "Clustering Performance vs. Number of Clusters\n(higher is better)",
        xlabel = "Number of Clusters (k)",
        ylabel = "Average Silhouette Score",
        label = "Avg Score",
        markersize = 5,
        markerstrokewidth = 0.5,
    )
    
    StatsPlots.ylims!(p, (ylims_lower, ylims_upper))
    StatsPlots.vline!(p, [optimal_k], color = :red, linestyle = :dash, label = "Inferred optimum = $(optimal_k)")
    
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identifies the optimal number of clusters using hierarchical clustering
and maximizing the average silhouette score, displaying progress.

Uses Clustering.clustering_quality for score calculation.

Args:
    distance_matrix: A square matrix of pairwise distances between items.
    ks: Optional pre-defined range or vector of k values to test. If provided, min_k and max_k are ignored.
    min_k: Minimum number of clusters to test (default: max(floor(log10(n)), 2)).
    max_k: Maximum number of clusters to test (default: min(n, ceil(sqrt(n)))).
    plot_backend: Either :cairomakie (default) or :statsplots for visualization.

Returns:
    A named tuple containing:
    - hcl: The hierarchical clustering result object.
    - optimal_number_of_clusters: The inferred optimal number of clusters (k).
    - assignments: Cluster assignments for the optimal k.
    - figure: The plot figure.
    - ks: The range/vector of k values tested.
    - silhouette_scores: Vector of average silhouette scores for each k.
"""
function identify_optimal_number_of_clusters(
        distance_matrix;
        ks = nothing,
        min_k = nothing,
        max_k = nothing,
        plot_backend = :cairomakie
    )
    # Validate that user doesn't pass both ks and min_k/max_k
    if !isnothing(ks) && (!isnothing(min_k) || !isnothing(max_k))
        error("Cannot specify both 'ks' and 'min_k'/'max_k'. Choose one approach.")
    end
    
    # Ensure the input is a square matrix
    if size(distance_matrix, 1) != size(distance_matrix, 2)
        error("Input distance_matrix must be square.")
    end
    
    n_items = size(distance_matrix, 1)
    if n_items < 2
        error("Need at least 2 items to perform clustering.")
    end

    # Compute hierarchical clustering
    println("Performing hierarchical clustering...")
    hcl = heirarchically_cluster_distance_matrix(distance_matrix)
    println("Hierarchical clustering complete.")

    N = Clustering.nnodes(hcl)
    if N != n_items
        @warn "Number of labels in Hclust object ($(N)) does not match distance matrix size ($(n_items)). Using Hclust labels count."
    end
    
    # Determine k range to test
    if isnothing(ks)
        # Use min_k and max_k with defaults
        default_min_k = max(Int(floor(log10(N))), 2)
        default_max_k = min(N, Int(ceil(sqrt(N))))
        
        actual_min_k = isnothing(min_k) ? default_min_k : min_k
        actual_max_k = isnothing(max_k) ? default_max_k : max_k
        
        if actual_max_k < 2
            println("Not enough items (N=$N) to test multiple cluster numbers (k >= 2). Cannot determine optimal k.")
            @warn "Returning trivial clustering with k=$N as optimal k could not be determined."
            return (hcl = hcl, optimal_number_of_clusters = N, assignments = Clustering.cutree(hcl, k=N), figure = nothing, ks = [N], silhouette_scores = [NaN])
        end
        
        ks = actual_min_k:actual_max_k
    else
        # Validate user-provided ks
        if isempty(ks)
            error("Provided 'ks' range/vector is empty.")
        end
        if any(k -> k < 2 || k > N, ks)
            error("All k values must be between 2 and $N (number of items).")
        end
    end
    
    k_vec = collect(ks)
    n_ks = length(k_vec)
    println("Evaluating average silhouette scores for k in range: $(minimum(k_vec)) to $(maximum(k_vec))...")

    # Preallocate array for silhouette scores
    avg_silhouette_scores = zeros(Float64, n_ks)

    # Progress meter setup
    pm = ProgressMeter.Progress(n_ks; dt=1, desc="Calculating Avg Silhouette Scores: ", barlen=50)

    # Parallel computation of silhouette scores
    Threads.@threads for i in 1:n_ks
        k = k_vec[i]
        avg_silhouette_scores[i] = Statistics.mean(Clustering.silhouettes(Clustering.cutree(hcl, k=k), distance_matrix))
        ProgressMeter.next!(pm)
    end

    # Validate results
    if all(s -> s == Float64(-Inf) || isnan(s), avg_silhouette_scores)
        error("Failed to calculate valid average silhouette scores for all values of k.")
    end

    # Find optimal k
    valid_indices = findall(s -> !isnan(s) && s != Float64(-Inf), avg_silhouette_scores)
    if isempty(valid_indices)
        error("No valid average silhouette scores were computed.")
    end

    max_score_valid, local_max_idx_valid = findmax(avg_silhouette_scores[valid_indices])
    max_idx_original = valid_indices[local_max_idx_valid]
    optimal_number_of_clusters = k_vec[max_idx_original]

    println("\nOptimal number of clusters inferred (max avg silhouette score): ", optimal_number_of_clusters)
    println("Maximum average silhouette score: ", max_score_valid)

    # Create visualization
    if plot_backend == :cairomakie
        fig = plot_silhouette_scores_makie(k_vec, avg_silhouette_scores, optimal_number_of_clusters)
    elseif plot_backend == :statsplots
        fig = plot_silhouette_scores_statsplots(k_vec, avg_silhouette_scores, optimal_number_of_clusters)
    else
        error("Unknown plot_backend: $(plot_backend). Choose :cairomakie or :statsplots.")
    end
    
    display(fig)

    return (
        hcl = hcl,
        optimal_number_of_clusters = optimal_number_of_clusters,
        assignments = Clustering.cutree(hcl, k=optimal_number_of_clusters),
        figure = fig,
        ks = k_vec,
        silhouette_scores = avg_silhouette_scores
    )
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function mmseqs_easy_linclust(;fasta, output=fasta*".mmseqs_easy_linclust", tmp=mktempdir())
#     Mycelia.add_bioconda_env("mmseqs2")   
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createdb $(fasta) $(fasta)_DB`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createindex --search-type 3 $(fasta)_DB $(tempdir())`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-linclust $(fasta)_DB $(output) $(tmp)`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta)_DB $(fasta)_DB $(output) $(output).tsv`)
#     return "$(output).tsv"
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Determines the optimal number of clusters for k-means clustering by maximizing the silhouette score.

# # Algorithm
# - Starts by evaluating the first 5 k values
# - Continues evaluation if optimal k is at the edge of evaluated range
# - Refines search by evaluating midpoints between k values around the current optimum
# - Iterates until convergence (optimal k remains stable)

# # Arguments
# - `distance_matrix::Matrix`: Square matrix of pairwise distances between points
# - `ks_to_try::Vector{Int}`: Vector of k values to evaluate. Defaults to [1, 2, ...] up to matrix size

# # Returns
# Named tuple containing:
# - `optimal_number_of_clusters::Int`: The k value giving highest silhouette score
# - `ks_assessed::Vector{Int}`: All k values that were evaluated
# - `within_cluster_sum_of_squares::Vector{Float64}`: WCSS for each k assessed
# - `silhouette_scores::Vector{Float64}`: Silhouette scores for each k assessed
# """
# function fit_optimal_number_of_clusters(distance_matrix, ks_to_try=[1, 2, Mycelia.ks(max=size(distance_matrix, 1))...])
#     # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
#     # N = size(distance_matrix, 1)
#     # ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
#     # Int(round(size(distance_matrix, 1)/2))
#     # insert!(ks_to_try, 2, Int(round(size(distance_matrix, 1)))))
#     # ks_to_try = vcat([2^i for i in 0:Int(floor(log2(size(distance_matrix, 1))))], size(distance_matrix, 1))
#     # @info "ks = $(ks_to_try)"
#     @show ks_to_try
    
#     # can calculate this for k >= 1
#     # within_cluster_sum_of_squares = Union{Float64, Missing}[]
#     within_cluster_sum_of_squares = Float64[]
#     # these are only valid for k >= 2 so set initial value to missing
#     # between_cluster_sum_of_squares = [missing, zeros(length(ks_to_try)-1)...]
#     # silhouette_scores = Union{Float64, Missing}[]
#     silhouette_scores = Float64[]

#     for k in ks_to_try[1:5]
#         @show k
#         @time this_clustering = Clustering.kmeans(distance_matrix, k)
#         push!(within_cluster_sum_of_squares, wcss(this_clustering))
#         if k == 1
#             this_silhouette_score = 0
#         else
#             this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
            
#         end
#         push!(silhouette_scores, this_silhouette_score)
#         @show this_silhouette_score
#     end      
#     optimal_silhouette, optimal_index = findmax(silhouette_scores)
#     while (optimal_index == length(silhouette_scores)) && (optimal_index != length(ks_to_try))
#         k = ks_to_try[optimal_index+1]
#         @show k
#         @time this_clustering = Clustering.kmeans(distance_matrix, k)
#         push!(within_cluster_sum_of_squares, wcss(this_clustering))
#         this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
#         push!(silhouette_scores, this_silhouette_score)
#         @show this_silhouette_score
#         optimal_silhouette, optimal_index = findmax(silhouette_scores)
#     end
#     previous_optimal_number_of_clusters = 0
#     optimal_number_of_clusters = ks_to_try[optimal_index]
#     done = false
#     while optimal_number_of_clusters != previous_optimal_number_of_clusters
#         @show optimal_number_of_clusters
#         if optimal_index == 1
#             window_of_focus = ks_to_try[optimal_index:optimal_index+1]
#             insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
#         elseif optimal_index == length(ks_to_try)
#             window_of_focus = ks_to_try[optimal_index-1:optimal_index]
#             insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
#         else
#             window_of_focus = ks_to_try[optimal_index-1:optimal_index+1]
#         end
#         # @show window_of_focus
#         midpoints = [
#             Int(round(Statistics.mean(window_of_focus[1:2]))),
#             Int(round(Statistics.mean(window_of_focus[2:3])))
#             ]
#         # @show sort(vcat(midpoints, window_of_focus))
#         @show midpoints
        
#         for k in midpoints
#             insertion_index = first(searchsorted(ks_to_try, k))
#             if ks_to_try[insertion_index] != k
#                 insert!(ks_to_try, insertion_index, k)
#                 @show k
#                 @time this_clustering = Clustering.kmeans(distance_matrix, k)
#                 insert!(within_cluster_sum_of_squares, insertion_index, wcss(this_clustering))
#                 this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
#                 insert!(silhouette_scores, insertion_index, this_silhouette_score)
#                 @show this_silhouette_score
#             end
#         end
        
#         previous_optimal_number_of_clusters = optimal_number_of_clusters
#         optimal_silhouette, optimal_index = findmax(silhouette_scores)
#         optimal_number_of_clusters = ks_to_try[optimal_index]
#     end
#     @assert length(within_cluster_sum_of_squares) == length(silhouette_scores)
#     ks_assessed = ks_to_try[1:length(within_cluster_sum_of_squares)]
#     return (;optimal_number_of_clusters, ks_assessed, within_cluster_sum_of_squares, silhouette_scores)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the document frequency of tokens across a collection of documents.

# Arguments
- `documents`: Collection of text documents where each document is a string

# Returns
- Dictionary mapping each unique token to the number of documents it appears in

# Description
Computes how many documents contain each unique token. Each document is tokenized 
by splitting on whitespace. Tokens are counted only once per document, regardless 
of how many times they appear within that document.
"""
function document_frequency(documents)
    document_tokens = Set(split(strip(first(documents))))
    countmap = StatsBase.countmap(document_tokens)
    for document in documents[2:end]
        document_tokens = Set(split(strip(document)))
        this_countmap = StatsBase.countmap(document_tokens)
        merge!(+, countmap, this_countmap)
    end
    return countmap
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Within-Cluster Sum of Squares (WCSS) for a clustering result.

# Arguments
- `clustering_result`: A clustering result object containing:
  - `counts`: Vector with number of points in each cluster
  - `assignments`: Vector of cluster assignments for each point
  - `costs`: Vector of distances/costs from each point to its cluster center

# Returns
- `Float64`: The total within-cluster sum of squared distances

# Description
WCSS measures the compactness of clusters by summing the squared distances 
between each data point and its assigned cluster center.
"""
function wcss(clustering_result)
    n_clusters = length(clustering_result.counts)
    total_squared_cost = 0.0
    for cluster_id in 1:n_clusters
        cluster_indices = clustering_result.assignments .== cluster_id
        total_squared_cost += sum(clustering_result.costs[cluster_indices] .^ 2)
    end
    return total_squared_cost
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Determine the optimal number of clusters using hierarchical clustering with iterative refinement.

# # Arguments
# - `distance_matrix::Matrix`: Square matrix of pairwise distances between observations
# - `ks_to_try::Vector{Int}`: Vector of cluster counts to evaluate (default: 1, 2, and sequence from `Mycelia.ks()`)

# # Returns
# Named tuple containing:
# - `optimal_number_of_clusters::Int`: Best number of clusters found
# - `ks_assessed::Vector{Int}`: All cluster counts that were evaluated
# - `silhouette_scores::Vector{Float64}`: Silhouette scores for each k assessed
# - `hclust_result`: Hierarchical clustering result object

# # Details
# Uses Ward's linkage method and silhouette scores to evaluate cluster quality. 
# Implements an adaptive search that focuses on promising regions between initially tested k values.
# For k=1, silhouette score is set to 0 as a special case.
# """
# function fit_optimal_number_of_clusters_hclust(distance_matrix, ks_to_try=[1, 2, Mycelia.ks(max=size(distance_matrix, 1))...])
#     # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
#     # N = size(distance_matrix, 1)
#     # ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
#     # Int(round(size(distance_matrix, 1)/2))
#     # insert!(ks_to_try, 2, Int(round(size(distance_matrix, 1)))))
#     # ks_to_try = vcat([2^i for i in 0:Int(floor(log2(size(distance_matrix, 1))))], size(distance_matrix, 1))
#     # @info "ks = $(ks_to_try)"
#     @show ks_to_try
    
#     # can calculate this for k >= 1
#     # within_cluster_sum_of_squares = Union{Float64, Missing}[]
#     # within_cluster_sum_of_squares = Float64[]
#     # these are only valid for k >= 2 so set initial value to missing
#     # between_cluster_sum_of_squares = [missing, zeros(length(ks_to_try)-1)...]
#     # silhouette_scores = Union{Float64, Missing}[]
#     silhouette_scores = Float64[]
#     # wcss_scores = Float64[]



#     @show "initial heirarchical clustering"
#     @time hclust_result = Clustering.hclust(distance_matrix, linkage=:ward, branchorder=:optimal)
#     # for k in ks_to_try[1:3]
#     for k in ks_to_try
#         @show k
#         this_clustering = Clustering.cutree(hclust_result, k=k)
#         # push!(within_cluster_sum_of_squares, wcss(this_clustering))
#         if k == 1
#             this_silhouette_score = 0
#             # this_wcss = Inf
#         else
#             this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
#             # this_wcss = wcss(this_clustering)
#         end
#         push!(silhouette_scores, this_silhouette_score)
#         # push!(within_cluster_sum_of_squares, this_wcss)
#         @show this_silhouette_score
#     end      
#     optimal_silhouette, optimal_index = findmax(silhouette_scores)
#     # while (optimal_index == length(silhouette_scores)) && (optimal_index != length(ks_to_try))
#     #     k = ks_to_try[optimal_index+1]
#     #     @show k
#     #     this_clustering = Clustering.cutree(hclust_result, k=k)
#     #     # push!(within_cluster_sum_of_squares, wcss(this_clustering))
#     #     this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
#     #     push!(silhouette_scores, this_silhouette_score)
#     #     @show this_silhouette_score
#     #     optimal_silhouette, optimal_index = findmax(silhouette_scores)
#     # end
#     previous_optimal_number_of_clusters = 0
#     optimal_number_of_clusters = ks_to_try[optimal_index]
#     done = false
#     while optimal_number_of_clusters != previous_optimal_number_of_clusters
#         @show optimal_number_of_clusters
#         if optimal_index == 1
#             window_of_focus = ks_to_try[optimal_index:optimal_index+1]
#             insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
#         elseif optimal_index == length(ks_to_try)
#             window_of_focus = ks_to_try[optimal_index-1:optimal_index]
#             insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
#         else
#             window_of_focus = ks_to_try[optimal_index-1:optimal_index+1]
#         end
#         # @show window_of_focus
#         midpoints = [
#             Int(round(Statistics.mean(window_of_focus[1:2]))),
#             Int(round(Statistics.mean(window_of_focus[2:3])))
#             ]
#         # @show sort(vcat(midpoints, window_of_focus))
#         @show midpoints
        
#         for k in midpoints
#             insertion_index = first(searchsorted(ks_to_try, k))
#             if ks_to_try[insertion_index] != k
#                 insert!(ks_to_try, insertion_index, k)
#                 @show k
#                 this_clustering = Clustering.cutree(hclust_result, k=k)
#                 this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
#                 insert!(silhouette_scores, insertion_index, this_silhouette_score)
#                 # this_wcss = wcss(this_clustering)
#                 # insert!(within_cluster_sum_of_squares, insertion_index, this_wcss)
#                 @show this_silhouette_score
#             end
#         end
        
#         previous_optimal_number_of_clusters = optimal_number_of_clusters
#         optimal_silhouette, optimal_index = findmax(silhouette_scores)
#         optimal_number_of_clusters = ks_to_try[optimal_index]
#     end
#     # @assert length(within_cluster_sum_of_squares) == length(silhouette_scores)
#     ks_assessed = ks_to_try[1:length(silhouette_scores)]
#     return (;optimal_number_of_clusters, ks_assessed, silhouette_scores, hclust_result)
# end

"""
    find_group_medoid(distance_matrix, indices, metric_type=:distance)

Find the medoid (most central point) within a group of indices.
For :distance, minimizes sum of distances. For :similarity, maximizes sum of similarities.
"""
function find_group_medoid(distance_matrix, indices; metric_type=:distance)
    if length(indices) == 1
        return indices[1]
    end
    
    best_score = metric_type == :distance ? Inf : -Inf
    best_idx = indices[1]
    
    for candidate in indices
        # Compute sum of pairwise metrics to all other members
        current_score = sum(distance_matrix[candidate, i] for i in indices if i != candidate)
        
        if metric_type == :distance
            if current_score < best_score
                best_score = current_score
                best_idx = candidate
            end
        else # :similarity
            if current_score > best_score
                best_score = current_score
                best_idx = candidate
            end
        end
    end
    
    return best_idx
end

"""
    find_second_most_diverse(distance_matrix, indices, medoid_idx, metric_type=:distance)

Find the point in the group that is most distant (or least similar) from the medoid.
"""
function find_second_most_diverse(distance_matrix, indices, medoid_idx; metric_type=:distance)
    if length(indices) <= 1
        return nothing
    end
    
    # We want to maximize distance or minimize similarity
    best_score = metric_type == :distance ? -Inf : Inf
    best_idx = nothing
    
    for candidate in indices
        if candidate == medoid_idx
            continue
        end
        
        val = distance_matrix[candidate, medoid_idx]
        
        if metric_type == :distance
            if val > best_score
                best_score = val
                best_idx = candidate
            end
        else # :similarity
            if val < best_score
                best_score = val
                best_idx = candidate
            end
        end
    end
    
    return best_idx
end

"""
    greedy_maxmin_diversity(distance_matrix, initial_selected, candidate_pool, 
                            n_to_select, weights, group_ids, max_per_group; 
                            metric_type=:distance)

Greedily select samples to maximize minimum pairwise distance (or minimize max similarity) to the selected set.
Uses `weights` (e.g. abundance) as a tie-breaker.
"""
function greedy_maxmin_diversity(
    distance_matrix, 
    initial_selected, 
    candidate_pool, 
    n_to_select,
    weights,
    group_ids,
    max_per_group;
    metric_type=:distance
)
    selected = copy(initial_selected)
    remaining = setdiff(candidate_pool, selected)
    
    # Track selection counts per group
    group_counts = Dict{Any, Int}()
    for idx in selected
        g = group_ids[idx]
        group_counts[g] = get(group_counts, g, 0) + 1
    end
    
    for _ in 1:n_to_select
        best_candidate = nothing
        
        # Tracking the "Min Distance" (we want to MAXIMIZE this)
        # Or "Max Similarity" (we want to MINIMIZE this)
        best_worst_pair_metric = metric_type == :distance ? -Inf : Inf
        
        best_weight = -Inf
        
        for candidate in remaining
            # Check constraints
            g = group_ids[candidate]
            # Handle max_per_group as Int or Dict
            limit = isa(max_per_group, Dict) ? get(max_per_group, g, typemax(Int)) : max_per_group
            
            if get(group_counts, g, 0) >= limit
                continue
            end
            
            # Calculate metric to nearest neighbor in selected set
            if isempty(selected)
                worst_pair_metric = metric_type == :distance ? Inf : -Inf
            elseif metric_type == :distance
                # Find the CLOSEST selected point (min distance)
                # We want the candidate where this value is LARGEST (MaxMin)
                worst_pair_metric = minimum(distance_matrix[candidate, s] for s in selected)
            else
                # Find the MOST SIMILAR selected point (max similarity)
                # We want the candidate where this value is SMALLEST (MinMax)
                worst_pair_metric = maximum(distance_matrix[candidate, s] for s in selected)
            end
            
            # Tie-breaking logic
            current_weight = isnothing(weights) ? 0.0 : weights[candidate]
            is_better = false
            
            if isempty(selected)
                if current_weight > best_weight
                    is_better = true
                end
            elseif metric_type == :distance
                if worst_pair_metric > best_worst_pair_metric * 1.0001 # slightly better distance
                    is_better = true
                elseif isapprox(worst_pair_metric, best_worst_pair_metric; rtol=1e-4) && current_weight > best_weight
                    is_better = true
                end
            else # similarity
                if worst_pair_metric < best_worst_pair_metric * 0.9999 # slightly lower max similarity (better)
                    is_better = true
                elseif isapprox(worst_pair_metric, best_worst_pair_metric; rtol=1e-4) && current_weight > best_weight
                    is_better = true
                end
            end
            
            if is_better
                best_worst_pair_metric = worst_pair_metric
                best_candidate = candidate
                best_weight = current_weight
            end
        end
        
        if isnothing(best_candidate)
            break
        end
        
        push!(selected, best_candidate)
        remaining = setdiff(remaining, [best_candidate])
        group_counts[group_ids[best_candidate]] = get(group_counts, group_ids[best_candidate], 0) + 1
    end
    
    return selected
end

"""
    select_diverse_representatives(distance_matrix, group_ids; 
                                   weights=nothing, 
                                   n_total=100, 
                                   max_per_group=2, 
                                   metric_type=:distance)

Select maximally diverse representatives using a multi-stage approach:
1. Select one medoid from each group (guarantees representation).
2. Select a second diverse sample from larger groups (if max_per_group >= 2).
3. Fill remaining slots using greedy MaxMin selection (or MinMax similarity).

# Arguments
- `distance_matrix`: Square matrix of pairwise metrics.
- `group_ids`: Vector of group identifiers (Strings, Ints, or Tuples) matching rows of the matrix.
- `weights`: Optional vector of importance weights (e.g., abundance) for tie-breaking.
- `n_total`: Target total number of samples.
- `max_per_group`: Integer limit or Dictionary mapping group_id => limit.
- `metric_type`: `:distance` (default, 0=identical) or `:similarity` (100% or 1.0=identical).
"""
function select_diverse_representatives(
    distance_matrix::AbstractMatrix,
    group_ids::AbstractVector;
    weights::Union{AbstractVector, Nothing}=nothing,
    n_total::Int=100,
    max_per_group::Union{Int, Dict}=2,
    metric_type::Symbol=:distance
)
    if size(distance_matrix, 1) != length(group_ids)
        error("Distance matrix dimension $(size(distance_matrix, 1)) does not match group_ids length $(length(group_ids))")
    end
    if !isnothing(weights) && length(weights) != length(group_ids)
        error("Weights vector length must match group_ids length")
    end

    # Group indices
    groups = Dict{Any, Vector{Int}}()
    for (i, g) in enumerate(group_ids)
        if !haskey(groups, g)
            groups[g] = Int[]
        end
        push!(groups[g], i)
    end
    
    selected_medoids = Int[]
    
    # Stage 1: Medoids & Second Diverse
    for (g, indices) in groups
        # 1. Medoid
        medoid_idx = find_group_medoid(distance_matrix, indices; metric_type=metric_type)
        push!(selected_medoids, medoid_idx)
        
        # 2. Second Diverse (if allowed)
        limit = isa(max_per_group, Dict) ? get(max_per_group, g, typemax(Int)) : max_per_group
        if limit >= 2 && length(indices) > 1
            second_idx = find_second_most_diverse(distance_matrix, indices, medoid_idx; metric_type=metric_type)
            if !isnothing(second_idx)
                push!(selected_medoids, second_idx)
            end
        end
    end
    
    # Stage 2: Fill remaining
    remaining_slots = n_total - length(selected_medoids)

    if remaining_slots > 0
        all_candidates = collect(1:size(distance_matrix, 1))
        final_selected = greedy_maxmin_diversity(
            distance_matrix,
            selected_medoids,
            all_candidates,
            remaining_slots,
            weights,
            group_ids,
            max_per_group;
            metric_type=metric_type
        )
        return final_selected
    else
        # If Stage 1 selected too many, we truncate based on simple logic or return all
        # Ideally, you might want to prioritize medoids of larger groups, but simply returning
        # the list truncated to n_total is a safe fallback to strictly enforce n_total.
        if length(selected_medoids) > n_total
            @warn "Stage 1 selected $(length(selected_medoids)) samples, which exceeds target $n_total. Returning truncated list."
            return selected_medoids[1:n_total]
        end
        return selected_medoids
    end
end


# =============================================================================
# Cluster Ranking and Relational Clustering Pipeline
# =============================================================================

"""
    ClusterRanking

Ranking of entities within a cluster, including medoid and backup selections.

# Fields
- `cluster_id::Int`: Cluster identifier
- `entity_indices::Vector{Int}`: All entity indices in this cluster
- `entity_ids::Vector{String}`: Entity identifiers corresponding to indices
- `medoid_index::Int`: Index of the cluster medoid (representative)
- `backup_index::Union{Int, Nothing}`: Second-most diverse entity (backup representative)
- `rankings::Vector{Int}`: Rank order of all entities (1 = medoid)
- `intra_cluster_distances::Vector{Float64}`: Distance from each entity to medoid
"""
struct ClusterRanking
    cluster_id::Int
    entity_indices::Vector{Int}
    entity_ids::Vector{String}
    medoid_index::Int
    backup_index::Union{Int, Nothing}
    rankings::Vector{Int}
    intra_cluster_distances::Vector{Float64}
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Rank members within each cluster by distance from the medoid.

Uses existing `find_group_medoid()` and `find_second_most_diverse()` functions
for medoid and backup selection.

# Arguments
- `distance_matrix::Matrix{Float64}`: Symmetric pairwise distance matrix
- `cluster_assignments::Vector{Int}`: Cluster ID for each entity (1-indexed)
- `entity_ids::Vector{String}`: Entity identifiers
- `metric_type::Symbol`: `:distance` (0=identical) or `:similarity` (1=identical). Default: `:distance`

# Returns
- `Vector{ClusterRanking}`: Rankings for each cluster, sorted by cluster_id

# Example
```julia
D = [0.0 0.1 0.5; 0.1 0.0 0.4; 0.5 0.4 0.0]
assignments = [1, 1, 2]
entity_ids = ["E1", "E2", "E3"]
rankings = Mycelia.rank_cluster_members(D, assignments, entity_ids)
```
"""
function rank_cluster_members(
    distance_matrix::Matrix{Float64},
    cluster_assignments::Vector{Int},
    entity_ids::Vector{String};
    metric_type::Symbol = :distance
)
    n_samples = size(distance_matrix, 1)

    if length(cluster_assignments) != n_samples
        throw(ArgumentError("Length of cluster_assignments ($(length(cluster_assignments))) must match matrix dimension ($n_samples)"))
    end

    if length(entity_ids) != n_samples
        throw(ArgumentError("Length of entity_ids ($(length(entity_ids))) must match matrix dimension ($n_samples)"))
    end

    unique_clusters = sort(unique(cluster_assignments))
    rankings = ClusterRanking[]

    for cluster_id in unique_clusters
        ## Get indices belonging to this cluster
        cluster_indices = findall(==(cluster_id), cluster_assignments)
        cluster_entity_ids = entity_ids[cluster_indices]

        if length(cluster_indices) == 1
            ## Single-member cluster: it's its own medoid
            push!(rankings, ClusterRanking(
                cluster_id,
                cluster_indices,
                cluster_entity_ids,
                cluster_indices[1],
                nothing,
                [1],
                [0.0]
            ))
            continue
        end

        ## Use existing find_group_medoid function
        medoid_idx = find_group_medoid(distance_matrix, cluster_indices; metric_type=metric_type)

        ## Use existing find_second_most_diverse function
        backup_idx = find_second_most_diverse(distance_matrix, cluster_indices, medoid_idx; metric_type=metric_type)

        ## Calculate distances from medoid to all cluster members
        distances_to_medoid = [distance_matrix[idx, medoid_idx] for idx in cluster_indices]

        ## Rank by distance (closest to medoid = rank 1)
        rank_order = sortperm(distances_to_medoid)
        rankings_vec = invperm(rank_order)

        push!(rankings, ClusterRanking(
            cluster_id,
            cluster_indices,
            cluster_entity_ids,
            medoid_idx,
            backup_idx,
            rankings_vec,
            distances_to_medoid
        ))
    end

    ## Sort by cluster_id for consistent output
    sort!(rankings, by=r -> r.cluster_id)

    return rankings
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert cluster rankings to a DataFrame for analysis and export.

# Arguments
- `rankings::Vector{ClusterRanking}`: Cluster ranking results
- `include_distances::Bool`: Include intra-cluster distances (default: true)

# Returns
- `DataFrames.DataFrame`: Table with columns:
  - `cluster_id`: Cluster identifier
  - `entity_id`: Entity identifier
  - `entity_index`: Original matrix index
  - `rank`: Rank within cluster (1 = medoid)
  - `is_medoid`: Whether this entity is the cluster medoid
  - `is_backup`: Whether this entity is the backup representative
  - `distance_to_medoid`: Distance from entity to medoid (if include_distances=true)
"""
function rankings_to_dataframe(
    rankings::Vector{ClusterRanking};
    include_distances::Bool = true
)
    rows = NamedTuple[]

    for cr in rankings
        for (i, entity_id) in enumerate(cr.entity_ids)
            entity_idx = cr.entity_indices[i]
            row = (
                cluster_id = cr.cluster_id,
                entity_id = entity_id,
                entity_index = entity_idx,
                rank = cr.rankings[i],
                is_medoid = entity_idx == cr.medoid_index,
                is_backup = !isnothing(cr.backup_index) && entity_idx == cr.backup_index
            )

            if include_distances
                row = merge(row, (distance_to_medoid = cr.intra_cluster_distances[i],))
            end

            push!(rows, row)
        end
    end

    return DataFrames.DataFrame(rows)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate summary statistics for cluster rankings.

# Arguments
- `rankings::Vector{ClusterRanking}`: Cluster ranking results
- `distance_matrix::Matrix{Float64}`: Original pairwise distance matrix

# Returns
- `DataFrames.DataFrame`: Summary with columns:
  - `cluster_id`: Cluster identifier
  - `n_members`: Number of entities in cluster
  - `medoid_id`: Entity ID of the medoid
  - `backup_id`: Entity ID of the backup (or missing)
  - `diameter`: Maximum intra-cluster distance
  - `avg_intra_distance`: Mean intra-cluster distance
  - `max_distance_to_medoid`: Maximum distance from any member to medoid
"""
function cluster_summary(
    rankings::Vector{ClusterRanking},
    distance_matrix::Matrix{Float64}
)
    rows = NamedTuple[]

    for cr in rankings
        indices = cr.entity_indices
        n_members = length(indices)

        ## Find medoid entity ID
        medoid_local_idx = findfirst(==(cr.medoid_index), indices)
        medoid_id = cr.entity_ids[medoid_local_idx]

        ## Find backup entity ID
        backup_id = if isnothing(cr.backup_index)
            missing
        else
            backup_local_idx = findfirst(==(cr.backup_index), indices)
            isnothing(backup_local_idx) ? missing : cr.entity_ids[backup_local_idx]
        end

        ## Compute intra-cluster distances
        intra_dists = Float64[]
        for i in 1:(n_members-1)
            for j in (i+1):n_members
                d = distance_matrix[indices[i], indices[j]]
                if !isnan(d)
                    push!(intra_dists, d)
                end
            end
        end

        diameter = isempty(intra_dists) ? 0.0 : maximum(intra_dists)
        avg_intra = isempty(intra_dists) ? 0.0 : Statistics.mean(intra_dists)

        push!(rows, (
            cluster_id = cr.cluster_id,
            n_members = n_members,
            medoid_id = medoid_id,
            backup_id = backup_id,
            diameter = diameter,
            avg_intra_distance = avg_intra,
            max_distance_to_medoid = maximum(cr.intra_cluster_distances)
        ))
    end

    return DataFrames.DataFrame(rows)
end


# =============================================================================
# Cluster Comparison Methods
# =============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a contingency matrix from two clusterings of the same samples.

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering

# Returns
- `NamedTuple` with fields:
  - `matrix::Matrix{Int}`: Contingency matrix (rows = labels1, cols = labels2)
  - `labels1::Vector`: Sorted unique labels from first clustering
  - `labels2::Vector`: Sorted unique labels from second clustering

# Example
```julia
labels_a = [1, 1, 2, 2, 3]
labels_b = [1, 1, 1, 2, 2]
result = Mycelia.contingency_matrix(labels_a, labels_b)
```
"""
function contingency_matrix(labels1::AbstractVector, labels2::AbstractVector)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length (got $(length(labels1)) and $(length(labels2)))"))
    end

    unique1 = sort(unique(labels1))
    unique2 = sort(unique(labels2))

    matrix = zeros(Int, length(unique1), length(unique2))
    idx1 = Dict(l => i for (i, l) in enumerate(unique1))
    idx2 = Dict(l => i for (i, l) in enumerate(unique2))

    for (l1, l2) in zip(labels1, labels2)
        matrix[idx1[l1], idx2[l2]] += 1
    end

    return (matrix = matrix, labels1 = unique1, labels2 = unique2)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Adjusted Rand Index (ARI) between two clusterings.

The ARI is a measure of similarity between two clusterings, adjusted for chance.
It ranges from -1 to 1, where:
- 1 indicates perfect agreement
- 0 indicates agreement no better than random chance
- Negative values indicate agreement worse than random

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering

# Returns
- `Float64`: The Adjusted Rand Index

# References
- Hubert, L. and Arabie, P. (1985). Comparing partitions. Journal of Classification.

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2]
labels_b = [1, 1, 2, 2, 2, 2]
ari = Mycelia.adjusted_rand_index(labels_a, labels_b)
```
"""
function adjusted_rand_index(labels1::AbstractVector, labels2::AbstractVector)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n < 2
        return 1.0  ## Trivial case
    end

    cont = contingency_matrix(labels1, labels2)
    nij = cont.matrix

    ## Sum of combinations C(n_ij, 2) for all cells
    sum_comb_nij = sum((binomial(Int(cell), 2) for cell in nij if cell >= 2), init=0)

    ## Row and column sums
    row_sums = vec(sum(nij, dims=2))
    col_sums = vec(sum(nij, dims=1))

    ## Sum of combinations for row sums
    sum_comb_ai = sum((binomial(Int(ai), 2) for ai in row_sums if ai >= 2), init=0)

    ## Sum of combinations for column sums
    sum_comb_bj = sum((binomial(Int(bj), 2) for bj in col_sums if bj >= 2), init=0)

    ## Total combinations
    total_comb = binomial(n, 2)

    if total_comb == 0
        return 1.0
    end

    ## Expected index under random permutation
    expected = (sum_comb_ai * sum_comb_bj) / total_comb

    ## Maximum index
    max_index = 0.5 * (sum_comb_ai + sum_comb_bj)

    ## Adjusted Rand Index
    if max_index == expected
        return 1.0  ## Perfect agreement when both sums are 0
    end

    ari = (sum_comb_nij - expected) / (max_index - expected)
    return ari
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Normalized Mutual Information (NMI) between two clusterings.

NMI is an information-theoretic measure of clustering similarity, normalized
by the geometric mean of the entropies. It ranges from 0 to 1, where:
- 1 indicates perfect agreement
- 0 indicates no mutual information

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering

# Returns
- `Float64`: The Normalized Mutual Information

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2]
labels_b = [1, 1, 2, 2, 2, 2]
nmi = Mycelia.normalized_mutual_information(labels_a, labels_b)
```
"""
function normalized_mutual_information(labels1::AbstractVector, labels2::AbstractVector)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n == 0
        return 1.0
    end

    ## Calculate cluster counts
    counts1 = StatsBase.countmap(labels1)
    counts2 = StatsBase.countmap(labels2)
    joint_counts = StatsBase.countmap(collect(zip(labels1, labels2)))

    ## Entropy of clustering 1
    H1 = 0.0
    for c in values(counts1)
        p = c / n
        if p > 0
            H1 -= p * log(p)
        end
    end

    ## Entropy of clustering 2
    H2 = 0.0
    for c in values(counts2)
        p = c / n
        if p > 0
            H2 -= p * log(p)
        end
    end

    ## Mutual information
    MI = 0.0
    for ((l1, l2), nij) in joint_counts
        if nij > 0
            p_ij = nij / n
            p_i = counts1[l1] / n
            p_j = counts2[l2] / n
            MI += p_ij * log(p_ij / (p_i * p_j))
        end
    end

    ## Normalize by geometric mean of entropies
    if H1 == 0.0 || H2 == 0.0
        return H1 == H2 ? 1.0 : 0.0
    end

    nmi = MI / sqrt(H1 * H2)
    return nmi
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Adjusted Mutual Information (AMI) between two clusterings.

AMI is the mutual information adjusted for chance, similar to how ARI adjusts
the Rand Index. It ranges from 0 to 1 (approximately), where:
- 1 indicates perfect agreement
- 0 indicates agreement no better than random chance

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering

# Returns
- `Float64`: The Adjusted Mutual Information

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2]
labels_b = [1, 1, 2, 2, 2, 2]
ami = Mycelia.adjusted_mutual_information(labels_a, labels_b)
```
"""
function adjusted_mutual_information(labels1::AbstractVector, labels2::AbstractVector)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n == 0
        return 1.0
    end

    cont = contingency_matrix(labels1, labels2)
    nij = cont.matrix

    ## Row and column sums
    a = vec(sum(nij, dims=2))  ## row sums
    b = vec(sum(nij, dims=1))  ## column sums

    R = length(a)  ## number of clusters in labels1
    C = length(b)  ## number of clusters in labels2

    ## Entropies
    H_U = -sum(ai/n * log(ai/n) for ai in a if ai > 0)
    H_V = -sum(bj/n * log(bj/n) for bj in b if bj > 0)

    ## Mutual information
    MI = 0.0
    for i in 1:R
        for j in 1:C
            if nij[i, j] > 0
                MI += (nij[i, j] / n) * log((n * nij[i, j]) / (a[i] * b[j]))
            end
        end
    end

    ## Expected mutual information (approximate)
    ## This is a simplified computation; exact EMI requires summing over hypergeometric
    EMI = 0.0
    for i in 1:R
        for j in 1:C
            for nij_val in max(1, a[i] + b[j] - n):min(a[i], b[j])
                ## Hypergeometric probability
                numerator = binomial(Int(a[i]), Int(nij_val)) * binomial(Int(n - a[i]), Int(b[j] - nij_val))
                denominator = binomial(Int(n), Int(b[j]))
                if denominator > 0
                    p_nij = numerator / denominator
                    if nij_val > 0 && p_nij > 0
                        EMI += (nij_val / n) * log((n * nij_val) / (a[i] * b[j])) * p_nij
                    end
                end
            end
        end
    end

    ## Adjusted Mutual Information
    max_MI = (H_U + H_V) / 2

    if max_MI == EMI
        return 1.0
    end

    ami = (MI - EMI) / (max_MI - EMI)
    return ami
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the V-measure between two clusterings.

V-measure is the harmonic mean of homogeneity and completeness:
- Homogeneity: each cluster contains only members of a single class
- Completeness: all members of a class are assigned to the same cluster

# Arguments
- `labels_true::AbstractVector`: Ground truth cluster assignments
- `labels_pred::AbstractVector`: Predicted cluster assignments

# Returns
- `NamedTuple` with fields:
  - `homogeneity::Float64`: Homogeneity score (0 to 1)
  - `completeness::Float64`: Completeness score (0 to 1)
  - `v_measure::Float64`: V-measure (harmonic mean, 0 to 1)

# Example
```julia
labels_true = [1, 1, 1, 2, 2, 2]
labels_pred = [1, 1, 2, 2, 2, 2]
result = Mycelia.v_measure(labels_true, labels_pred)
println("V-measure: ", result.v_measure)
```
"""
function v_measure(labels_true::AbstractVector, labels_pred::AbstractVector)
    n = length(labels_true)
    if length(labels_pred) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n == 0
        return (homogeneity = 1.0, completeness = 1.0, v_measure = 1.0)
    end

    counts_true = StatsBase.countmap(labels_true)
    counts_pred = StatsBase.countmap(labels_pred)
    joint_counts = StatsBase.countmap(collect(zip(labels_true, labels_pred)))

    ## Entropies
    H_true = -sum((c/n) * log(c/n) for c in values(counts_true) if c > 0)
    H_pred = -sum((c/n) * log(c/n) for c in values(counts_pred) if c > 0)

    ## Conditional entropy H(true | pred)
    H_true_given_pred = 0.0
    for (l_pred, n_pred) in counts_pred
        for l_true in keys(counts_true)
            n_joint = get(joint_counts, (l_true, l_pred), 0)
            if n_joint > 0
                H_true_given_pred -= (n_joint / n) * log(n_joint / n_pred)
            end
        end
    end

    ## Conditional entropy H(pred | true)
    H_pred_given_true = 0.0
    for (l_true, n_true) in counts_true
        for l_pred in keys(counts_pred)
            n_joint = get(joint_counts, (l_true, l_pred), 0)
            if n_joint > 0
                H_pred_given_true -= (n_joint / n) * log(n_joint / n_true)
            end
        end
    end

    ## Homogeneity: 1 - H(true | pred) / H(true)
    homogeneity = H_true == 0.0 ? 1.0 : 1.0 - H_true_given_pred / H_true

    ## Completeness: 1 - H(pred | true) / H(pred)
    completeness = H_pred == 0.0 ? 1.0 : 1.0 - H_pred_given_true / H_pred

    ## V-measure: harmonic mean
    if homogeneity + completeness == 0.0
        v = 0.0
    else
        v = 2.0 * homogeneity * completeness / (homogeneity + completeness)
    end

    return (homogeneity = homogeneity, completeness = completeness, v_measure = v)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Fowlkes-Mallows Index (FMI) between two clusterings.

FMI is the geometric mean of precision and recall computed from pair counting.
It ranges from 0 to 1, where 1 indicates perfect agreement.

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering

# Returns
- `Float64`: The Fowlkes-Mallows Index

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2]
labels_b = [1, 1, 2, 2, 2, 2]
fmi = Mycelia.fowlkes_mallows_index(labels_a, labels_b)
```
"""
function fowlkes_mallows_index(labels1::AbstractVector, labels2::AbstractVector)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n < 2
        return 1.0
    end

    ## Count pairs
    TP = 0  ## True positives: same cluster in both
    FP = 0  ## False positives: same in 1, different in 2
    FN = 0  ## False negatives: different in 1, same in 2

    for i in 1:(n-1)
        for j in (i+1):n
            same1 = labels1[i] == labels1[j]
            same2 = labels2[i] == labels2[j]

            if same1 && same2
                TP += 1
            elseif same1 && !same2
                FP += 1
            elseif !same1 && same2
                FN += 1
            end
        end
    end

    if TP == 0
        return 0.0
    end

    precision = TP / (TP + FP)
    recall = TP / (TP + FN)

    fmi = sqrt(precision * recall)
    return fmi
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute cluster purity between true and predicted clusterings.

Purity measures the extent to which clusters contain a single class.
For each predicted cluster, the majority true class is found, and
the proportion of correctly assigned samples is computed.

# Arguments
- `labels_true::AbstractVector`: Ground truth cluster assignments
- `labels_pred::AbstractVector`: Predicted cluster assignments

# Returns
- `NamedTuple` with fields:
  - `overall_purity::Float64`: Overall purity (0 to 1)
  - `cluster_purities::Dict`: Purity for each predicted cluster

# Example
```julia
labels_true = [1, 1, 1, 2, 2, 2]
labels_pred = [1, 1, 2, 2, 2, 2]
result = Mycelia.cluster_purity(labels_true, labels_pred)
println("Overall purity: ", result.overall_purity)
```
"""
function cluster_purity(labels_true::AbstractVector, labels_pred::AbstractVector)
    n = length(labels_true)
    if length(labels_pred) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n == 0
        return (overall_purity = 1.0, cluster_purities = Dict())
    end

    joint_counts = StatsBase.countmap(collect(zip(labels_pred, labels_true)))
    pred_clusters = unique(labels_pred)

    total_correct = 0
    cluster_purities = Dict{eltype(labels_pred), Float64}()

    for pred_c in pred_clusters
        ## Find counts of true labels for this predicted cluster
        true_counts = Dict{eltype(labels_true), Int}()
        for ((pc, tc), count) in joint_counts
            if pc == pred_c
                true_counts[tc] = count
            end
        end

        if !isempty(true_counts)
            max_count = maximum(values(true_counts))
            total_correct += max_count
            cluster_size = sum(values(true_counts))
            cluster_purities[pred_c] = max_count / cluster_size
        end
    end

    overall_purity = total_correct / n
    return (overall_purity = overall_purity, cluster_purities = cluster_purities)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Jaccard Index between two clusterings based on pair counting.

The Jaccard Index measures agreement by comparing pairs of samples that are
clustered together. It ranges from 0 to 1, where 1 indicates perfect agreement.

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering

# Returns
- `Float64`: The Jaccard Index

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2]
labels_b = [1, 1, 2, 2, 2, 2]
jaccard = Mycelia.jaccard_index(labels_a, labels_b)
```
"""
function jaccard_index(labels1::AbstractVector, labels2::AbstractVector)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    if n < 2
        return 1.0
    end

    ## Count pairs
    a = 0  ## Pairs together in both
    b = 0  ## Pairs together in 1 only
    c = 0  ## Pairs together in 2 only

    for i in 1:(n-1)
        for j in (i+1):n
            same1 = labels1[i] == labels1[j]
            same2 = labels2[i] == labels2[j]

            if same1 && same2
                a += 1
            elseif same1 && !same2
                b += 1
            elseif !same1 && same2
                c += 1
            end
        end
    end

    if a + b + c == 0
        return 1.0
    end

    return a / (a + b + c)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a comprehensive summary comparing two clusterings.

Computes multiple comparison metrics and returns them in a single report.

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering
- `name1::String`: Name for the first clustering (default: "Clustering1")
- `name2::String`: Name for the second clustering (default: "Clustering2")

# Returns
- `NamedTuple` with all comparison metrics:
  - `n_samples`: Number of samples
  - `n_clusters1`: Number of clusters in first clustering
  - `n_clusters2`: Number of clusters in second clustering
  - `adjusted_rand_index`: ARI score
  - `normalized_mutual_information`: NMI score
  - `adjusted_mutual_information`: AMI score
  - `v_measure`: V-measure results (homogeneity, completeness, v_measure)
  - `fowlkes_mallows_index`: FMI score
  - `jaccard_index`: Jaccard Index
  - `purity_1_to_2`: Purity treating labels1 as ground truth
  - `purity_2_to_1`: Purity treating labels2 as ground truth
  - `contingency`: Contingency matrix result

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2, 3, 3]
labels_b = [1, 1, 2, 2, 2, 3, 3, 3]
summary = Mycelia.clustering_comparison_summary(labels_a, labels_b)
println("ARI: ", summary.adjusted_rand_index)
println("NMI: ", summary.normalized_mutual_information)
```
"""
function clustering_comparison_summary(
    labels1::AbstractVector,
    labels2::AbstractVector;
    name1::String = "Clustering1",
    name2::String = "Clustering2"
)
    n = length(labels1)
    if length(labels2) != n
        throw(ArgumentError("Label vectors must have the same length"))
    end

    ## Compute all metrics
    ari = adjusted_rand_index(labels1, labels2)
    nmi = normalized_mutual_information(labels1, labels2)
    ami = adjusted_mutual_information(labels1, labels2)
    vm = v_measure(labels1, labels2)
    fmi = fowlkes_mallows_index(labels1, labels2)
    ji = jaccard_index(labels1, labels2)
    purity_1_to_2 = cluster_purity(labels1, labels2)
    purity_2_to_1 = cluster_purity(labels2, labels1)
    cont = contingency_matrix(labels1, labels2)

    return (
        name1 = name1,
        name2 = name2,
        n_samples = n,
        n_clusters1 = length(unique(labels1)),
        n_clusters2 = length(unique(labels2)),
        adjusted_rand_index = ari,
        normalized_mutual_information = nmi,
        adjusted_mutual_information = ami,
        v_measure = vm,
        fowlkes_mallows_index = fmi,
        jaccard_index = ji,
        purity_1_to_2 = purity_1_to_2,
        purity_2_to_1 = purity_2_to_1,
        contingency = cont
    )
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Print a formatted comparison report for two clusterings.

# Arguments
- `labels1::AbstractVector`: Cluster assignments from first clustering
- `labels2::AbstractVector`: Cluster assignments from second clustering
- `name1::String`: Name for the first clustering (default: "Clustering1")
- `name2::String`: Name for the second clustering (default: "Clustering2")
- `io::IO`: Output stream (default: stdout)

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2, 3, 3]
labels_b = [1, 1, 2, 2, 2, 3, 3, 3]
Mycelia.print_clustering_comparison(labels_a, labels_b; name1="K-means", name2="DBSCAN")
```
"""
function print_clustering_comparison(
    labels1::AbstractVector,
    labels2::AbstractVector;
    name1::String = "Clustering1",
    name2::String = "Clustering2",
    io::IO = stdout
)
    summary = clustering_comparison_summary(labels1, labels2; name1=name1, name2=name2)

    println(io, "")
    println(io, "=" ^ 60)
    println(io, "CLUSTERING COMPARISON REPORT")
    println(io, "=" ^ 60)
    println(io, "")
    println(io, "Datasets:")
    println(io, "  - $(name1): $(summary.n_clusters1) clusters")
    println(io, "  - $(name2): $(summary.n_clusters2) clusters")
    println(io, "  - Common samples: $(summary.n_samples)")
    println(io, "")
    println(io, "Agreement Metrics:")
    println(io, "  Adjusted Rand Index (ARI):       $(round(summary.adjusted_rand_index, digits=4))")
    println(io, "  Normalized Mutual Info (NMI):    $(round(summary.normalized_mutual_information, digits=4))")
    println(io, "  Adjusted Mutual Info (AMI):      $(round(summary.adjusted_mutual_information, digits=4))")
    println(io, "  V-measure:                       $(round(summary.v_measure.v_measure, digits=4))")
    println(io, "    - Homogeneity:                 $(round(summary.v_measure.homogeneity, digits=4))")
    println(io, "    - Completeness:                $(round(summary.v_measure.completeness, digits=4))")
    println(io, "  Fowlkes-Mallows Index (FMI):     $(round(summary.fowlkes_mallows_index, digits=4))")
    println(io, "  Jaccard Index:                   $(round(summary.jaccard_index, digits=4))")
    println(io, "")
    println(io, "Purity Scores:")
    println(io, "  $(name1) -> $(name2): $(round(summary.purity_1_to_2.overall_purity, digits=4))")
    println(io, "  $(name2) -> $(name1): $(round(summary.purity_2_to_1.overall_purity, digits=4))")
    println(io, "")
    println(io, "Interpretation Guidelines:")
    println(io, "  Score > 0.9:  Excellent agreement")
    println(io, "  Score 0.7-0.9: Good agreement")
    println(io, "  Score 0.5-0.7: Moderate agreement")
    println(io, "  Score 0.3-0.5: Weak agreement")
    println(io, "  Score < 0.3:  Poor agreement")
    println(io, "")

    return summary
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a clustering comparison summary to a DataFrame for easy export.

# Arguments
- `summary::NamedTuple`: Result from `clustering_comparison_summary`

# Returns
- `DataFrames.DataFrame`: Single-row DataFrame with all metrics

# Example
```julia
labels_a = [1, 1, 1, 2, 2, 2]
labels_b = [1, 1, 2, 2, 2, 2]
summary = Mycelia.clustering_comparison_summary(labels_a, labels_b)
df = Mycelia.comparison_summary_to_dataframe(summary)
```
"""
function comparison_summary_to_dataframe(summary::NamedTuple)
    return DataFrames.DataFrame(
        name1 = [summary.name1],
        name2 = [summary.name2],
        n_samples = [summary.n_samples],
        n_clusters1 = [summary.n_clusters1],
        n_clusters2 = [summary.n_clusters2],
        adjusted_rand_index = [summary.adjusted_rand_index],
        normalized_mutual_information = [summary.normalized_mutual_information],
        adjusted_mutual_information = [summary.adjusted_mutual_information],
        homogeneity = [summary.v_measure.homogeneity],
        completeness = [summary.v_measure.completeness],
        v_measure = [summary.v_measure.v_measure],
        fowlkes_mallows_index = [summary.fowlkes_mallows_index],
        jaccard_index = [summary.jaccard_index],
        purity_1_to_2 = [summary.purity_1_to_2.overall_purity],
        purity_2_to_1 = [summary.purity_2_to_1.overall_purity]
    )
end
