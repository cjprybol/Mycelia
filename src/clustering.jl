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
- `tmp::Union{String, Nothing}`: Path to temporary directory. If nothing, defaults to a subdirectory in output location.

# Returns
- `String`: Path to the output TSV file.

# Details
Uses MMseqs2 with minimum sequence identity threshold of 50% (-min-seq-id 0.5) and 
minimum coverage threshold of 80% (-c 0.8). The output TSV file format contains 
tab-separated cluster representative and member sequences.
"""
function mmseqs_easy_cluster(;fasta, output=fasta*".mmseqs_easy_cluster", tmp=nothing,
                             executor=:local, job_name="mmseqs_cluster", time="4-00:00:00", partition=nothing, account=nothing,
                             mem::Union{String, Nothing}=nothing, slurm_template::Union{Symbol, Nothing}=nothing, slurm_kwargs=Dict{Symbol,Any}()
                             )
    
    exec_instance = Mycelia.Execution.resolve_executor(executor)
    if exec_instance isa Mycelia.Execution.SlurmExecutor && !isnothing(slurm_template)
        exec_instance = Mycelia.Execution.SlurmExecutor(exec_instance.submit_script, slurm_template)
    end
    
    resolved_mem, resolved_threads = Mycelia.Execution.resolve_job_resources(
        mem, nothing, slurm_template, "64G", get_default_threads()
    )

    outfile = "$(output)_cluster.tsv"
    
    # Handle temp directory
    # If tmp is not provided, use a directory alongside output to ensure visibility on shared FS
    if isnothing(tmp)
        tmp_dir = "$(output)_tmp"
    else
        tmp_dir = tmp
    end
    
    Mycelia.Execution.with_executor(exec_instance) do
        run_analysis = true
        if exec_instance isa Mycelia.Execution.LocalExecutor
            if isfile(outfile)
               run_analysis = false 
            end
        end
        
        if run_analysis
            Mycelia.add_bioconda_env("mmseqs2")
            
            cmd_parts = String[]
            
            # Create tmp dir
            push!(cmd_parts, "mkdir -p $(tmp_dir)")
            
            # MMseqs command
            # --threads $(resolved_threads)
            # --min-seq-id 0.5 -c 0.8
            mmseqs_cmd = "$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-cluster $(fasta) $(output) $(tmp_dir) --threads $(resolved_threads) --min-seq-id 0.5 -c 0.8"
            push!(cmd_parts, mmseqs_cmd)
            
            # Cleanup tmp dir
            push!(cmd_parts, "rm -rf $(tmp_dir)")
            
            final_cmd = join(cmd_parts, "\n")
            
            Mycelia.Execution.submit(Mycelia.Execution.JobSpec(
                final_cmd;
                name=job_name,
                time=time,
                cpus=resolved_threads,
                mem=resolved_mem,
                partition=partition,
                account=account,
                extra=slurm_kwargs
            ))
        end
    end
    
    return outfile
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
    for leaf_node in hcl.labels
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
    fig = CairoMakie.Figure(resolution = (900, 600), fontsize = 16)
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

    N = length(hcl.labels)
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
    pm = ProgressMeter.Progress(n_ks, 1, "Calculating Avg Silhouette Scores: ", 50)

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
