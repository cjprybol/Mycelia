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
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a hierarchical clustering tree into a directed metagraph representation.

# Arguments
- `hcl::Clustering.Hclust`: Hierarchical clustering result object

# Returns
- `MetaGraphs.MetaDiGraph`: Directed graph with metadata representing the clustering hierarchy

# Graph Properties
The resulting graph contains the following vertex properties:
- `:hclust_id`: String identifier for each node
- `:height`: Height/distance at each merge point (0.0 for leaves)
- `:x`: Horizontal position for visualization (0-1 range)
- `:y`: Vertical position based on normalized height
- `:hcl`: Original clustering object (stored as graph property)
"""
function hclust_to_metagraph(hcl::Clustering.Hclust)
    total_nodes = length(hcl.order) + size(hcl.merges, 1)
    mg = MetaGraphs.MetaDiGraph(total_nodes)
    MetaGraphs.set_prop!(mg, :hcl, hcl)
    for leaf_node in hcl.labels
        MetaGraphs.set_prop!(mg, leaf_node, :hclust_id, string(-leaf_node))
        MetaGraphs.set_prop!(mg, leaf_node, :height, 0.0)
    end
    for (i, ordered_leaf_node) in enumerate(hcl.order)
        x = i / (length(hcl.order) + 1)
        MetaGraphs.set_prop!(mg, ordered_leaf_node, :x, x)
    end
    for (i, (left, right)) in enumerate(eachrow(hcl.merges))
        graph_i = length(hcl.order) + i
        hclust_id = string(i)
        MetaGraphs.set_prop!(mg, graph_i, :hclust_id, hclust_id)
    end
    MetaGraphs.set_indexing_prop!(mg, :hclust_id)

    for (i, ((left, right), height)) in enumerate(zip(eachrow(hcl.merges), hcl.heights))
        parent_vertex = mg[string(i), :hclust_id]
        MetaGraphs.set_prop!(mg, parent_vertex, :height, height)
        
        left_child_vertex = mg[string(left), :hclust_id]
        right_child_vertex = mg[string(right), :hclust_id]

        x = Statistics.middle(mg.vprops[left_child_vertex][:x], mg.vprops[right_child_vertex][:x])
        MetaGraphs.set_prop!(mg, parent_vertex, :x, x)

        e1 = Graphs.Edge(parent_vertex, left_child_vertex)
        e2 = Graphs.Edge(parent_vertex, right_child_vertex)
        Graphs.add_edge!(mg, e1)
        Graphs.add_edge!(mg, e2)
    end

    max_height = maximum(get.(values(mg.vprops), :height, missing))
    for v in Graphs.vertices(mg)
        relative_height = (max_height - mg.vprops[v][:height]) / max_height
        # y = 1 - relative_height
        MetaGraphs.set_prop!(mg, v, :y, relative_height)
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

Identifies the optimal number of clusters using hierarchical clustering
and maximizing the average silhouette score, displaying progress.

Uses Clustering.clustering_quality for score calculation.

Args:
    distance_matrix: A square matrix of pairwise distances between items.

Returns:
    A tuple containing:
    - hcl: The hierarchical clustering result object.
    - optimal_number_of_clusters: The inferred optimal number of clusters (k).
"""
function identify_optimal_number_of_clusters(
        distance_matrix;
        min_k = max(Int(floor(log(size(distance_matrix, 2)))), 2),
        max_k = min(size(distance_matrix, 2), Int(ceil(sqrt(size(distance_matrix, 2))))),
    )
    # Ensure the input is a square matrix
    if size(distance_matrix, 1) != size(distance_matrix, 2)
         error("Input distance_matrix must be square.")
    end
    n_items = size(distance_matrix, 1)
    if n_items < 2
        error("Need at least 2 items to perform clustering.")
    end

    # Compute hierarchical clustering using your custom function.
    println("Performing hierarchical clustering...")
    hcl = heirarchically_cluster_distance_matrix(distance_matrix)
    println("Hierarchical clustering complete.")

    # # Determine N based on hcl object or matrix size
    # local N::Int
    # if hasfield(typeof(hcl), :labels) && !isempty(hcl.labels)
    N = length(hcl.labels)
    if N != n_items
         @warn "Number of labels in Hclust object ($(N)) does not match distance matrix size ($(n_items)). Using Hclust labels count."
    end
    # else
    #      @warn "Hierarchical clustering result type $(typeof(hcl)) might not have a `.labels` field or it's empty. Using distance matrix size for N."
    #      N = n_items
    # end
    
    # Define the range of k values to test: 2 to min(N-1, ceil(sqrt(N)))
    # Cannot have more clusters than N items, and silhouette requires at least 1 item per cluster,
    # implicitly k <= N. Silhouettes are ill-defined for k=1. Need k < N for non-trivial clustering.
    if max_k < 2
        println("Not enough items (N=$N) to test multiple cluster numbers (k >= 2). Cannot determine optimal k.")
        # Optionally, plot a message or return a specific value indicating failure/trivial case
        p_trivial = StatsPlots.plot(title="Cannot optimize k (N=$N, max tested k=$(max_k))", xlabel="Number of Clusters", ylabel="Avg Silhouette Score")
        display(p_trivial)
        # Returning N clusters (each item its own) might be misleading as optimal
        # Consider returning nothing or throwing error, or returning k=N with a warning
        @warn "Returning trivial clustering with k=$N as optimal k could not be determined."
        return (hcl, N)
    end
    ks = min_k:max_k
    n_ks = length(ks)
    println("Evaluating average silhouette scores for k from $(min_k) to $(max_k)...")

    # Preallocate an array for average silhouette scores for each value of k.
    avg_silhouette_scores = zeros(Float64, n_ks)

    # --- ProgressMeter Setup ---
    pm = ProgressMeter.Progress(n_ks, 1, "Calculating Avg Silhouette Scores: ", 50)
    # --------------------------

    # Parallel loop: each iteration computes the average silhouette score for a given k.
    Threads.@threads for i in 1:n_ks
        k = ks[i]
        avg_silhouette_scores[i] = Statistics.mean(Clustering.silhouettes(Clustering.cutree(hcl, k=k), distance_matrix))
        ProgressMeter.next!(pm)
    end

    # Check if all scores resulted in errors
    if all(s -> s == Float64(-Inf) || isnan(s), avg_silhouette_scores)
        error("Failed to calculate valid average silhouette scores for all values of k.")
    end

    # Determine the optimal number of clusters by finding the MAXIMUM average score
    # Filter out -Inf/NaN before finding the maximum.
    valid_indices = findall(s -> !isnan(s) && s != Float64(-Inf), avg_silhouette_scores)
    if isempty(valid_indices)
        error("No valid average silhouette scores were computed.")
    end

    # Find maximum among valid scores
    max_score_valid, local_max_idx_valid = findmax(avg_silhouette_scores[valid_indices])
    # Map back to the original index in avg_silhouette_scores and ks
    max_idx_original = valid_indices[local_max_idx_valid]
    optimal_number_of_clusters = ks[max_idx_original]

    println("\nOptimal number of clusters inferred (max avg silhouette score): ", optimal_number_of_clusters)
    println("Maximum average silhouette score: ", max_score_valid)


    # Plotting the average silhouette scores.
    p = StatsPlots.scatter(
        ks,
        avg_silhouette_scores,
        title = "Clustering Performance vs. Number of Clusters\n(higher is better)",
        xlabel = "Number of Clusters (k)",
        ylabel = "Average Silhouette Score", # Updated label
        label = "Avg Score",
        markersize = 5,
        markerstrokewidth = 0.5,
        # legend=:outertopright
    )

    # Adjust y-limits based on valid scores, keeping [-1, 1] range in mind
    valid_scores = avg_silhouette_scores[valid_indices]
    min_val = minimum(valid_scores)
    max_val = maximum(valid_scores) # This is max_score_valid
    range = max_val - min_val
    # Set sensible limits slightly beyond observed range, but clamp to [-1.05, 1.05]
    ylims_lower = max(-1.05, min(0, min_val - range * 0.1))
    ylims_upper = min(1.05, max_val + range * 0.1)
    # Ensure lower < upper, handle case where range is 0
    if isapprox(ylims_lower, ylims_upper)
         ylims_lower -= 0.1
         ylims_upper += 0.1
         # Clamp again if needed
         ylims_lower = max(-1.05, ylims_lower)
         ylims_upper = min(1.05, ylims_upper)
    end
    StatsPlots.ylims!(p, (ylims_lower, ylims_upper))

    # Add reference line at y=0 (separation between reasonable/poor clusters)
    # StatsPlots.hline!(p, [0], color=:gray, linestyle=:dot, label="Score = 0")

    # Add vertical line for the optimum
    StatsPlots.vline!(p, [optimal_number_of_clusters], color = :red, linestyle = :dash, label = "Inferred optimum = $(optimal_number_of_clusters)")

    # Display the plot
    display(p)

    return (;hcl, optimal_number_of_clusters, assignments = Clustering.cutree(hcl, k=optimal_number_of_clusters), figure=p)
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the optimal number of clusters for k-means clustering by maximizing the silhouette score.

# Algorithm
- Starts by evaluating the first 5 k values
- Continues evaluation if optimal k is at the edge of evaluated range
- Refines search by evaluating midpoints between k values around the current optimum
- Iterates until convergence (optimal k remains stable)

# Arguments
- `distance_matrix::Matrix`: Square matrix of pairwise distances between points
- `ks_to_try::Vector{Int}`: Vector of k values to evaluate. Defaults to [1, 2, ...] up to matrix size

# Returns
Named tuple containing:
- `optimal_number_of_clusters::Int`: The k value giving highest silhouette score
- `ks_assessed::Vector{Int}`: All k values that were evaluated
- `within_cluster_sum_of_squares::Vector{Float64}`: WCSS for each k assessed
- `silhouette_scores::Vector{Float64}`: Silhouette scores for each k assessed
"""
function fit_optimal_number_of_clusters(distance_matrix, ks_to_try=[1, 2, Mycelia.ks(max=size(distance_matrix, 1))...])
    # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
    # N = size(distance_matrix, 1)
    # ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
    # Int(round(size(distance_matrix, 1)/2))
    # insert!(ks_to_try, 2, Int(round(size(distance_matrix, 1)))))
    # ks_to_try = vcat([2^i for i in 0:Int(floor(log2(size(distance_matrix, 1))))], size(distance_matrix, 1))
    # @info "ks = $(ks_to_try)"
    @show ks_to_try
    
    # can calculate this for k >= 1
    # within_cluster_sum_of_squares = Union{Float64, Missing}[]
    within_cluster_sum_of_squares = Float64[]
    # these are only valid for k >= 2 so set initial value to missing
    # between_cluster_sum_of_squares = [missing, zeros(length(ks_to_try)-1)...]
    # silhouette_scores = Union{Float64, Missing}[]
    silhouette_scores = Float64[]

    for k in ks_to_try[1:5]
        @show k
        @time this_clustering = Clustering.kmeans(distance_matrix, k)
        push!(within_cluster_sum_of_squares, wcss(this_clustering))
        if k == 1
            this_silhouette_score = 0
        else
            this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
            
        end
        push!(silhouette_scores, this_silhouette_score)
        @show this_silhouette_score
    end      
    optimal_silhouette, optimal_index = findmax(silhouette_scores)
    while (optimal_index == length(silhouette_scores)) && (optimal_index != length(ks_to_try))
        k = ks_to_try[optimal_index+1]
        @show k
        @time this_clustering = Clustering.kmeans(distance_matrix, k)
        push!(within_cluster_sum_of_squares, wcss(this_clustering))
        this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
        push!(silhouette_scores, this_silhouette_score)
        @show this_silhouette_score
        optimal_silhouette, optimal_index = findmax(silhouette_scores)
    end
    previous_optimal_number_of_clusters = 0
    optimal_number_of_clusters = ks_to_try[optimal_index]
    done = false
    while optimal_number_of_clusters != previous_optimal_number_of_clusters
        @show optimal_number_of_clusters
        if optimal_index == 1
            window_of_focus = ks_to_try[optimal_index:optimal_index+1]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        elseif optimal_index == length(ks_to_try)
            window_of_focus = ks_to_try[optimal_index-1:optimal_index]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        else
            window_of_focus = ks_to_try[optimal_index-1:optimal_index+1]
        end
        # @show window_of_focus
        midpoints = [
            Int(round(Statistics.mean(window_of_focus[1:2]))),
            Int(round(Statistics.mean(window_of_focus[2:3])))
            ]
        # @show sort(vcat(midpoints, window_of_focus))
        @show midpoints
        
        for k in midpoints
            insertion_index = first(searchsorted(ks_to_try, k))
            if ks_to_try[insertion_index] != k
                insert!(ks_to_try, insertion_index, k)
                @show k
                @time this_clustering = Clustering.kmeans(distance_matrix, k)
                insert!(within_cluster_sum_of_squares, insertion_index, wcss(this_clustering))
                this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
                insert!(silhouette_scores, insertion_index, this_silhouette_score)
                @show this_silhouette_score
            end
        end
        
        previous_optimal_number_of_clusters = optimal_number_of_clusters
        optimal_silhouette, optimal_index = findmax(silhouette_scores)
        optimal_number_of_clusters = ks_to_try[optimal_index]
    end
    @assert length(within_cluster_sum_of_squares) == length(silhouette_scores)
    ks_assessed = ks_to_try[1:length(within_cluster_sum_of_squares)]
    return (;optimal_number_of_clusters, ks_assessed, within_cluster_sum_of_squares, silhouette_scores)
end

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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determine the optimal number of clusters using hierarchical clustering with iterative refinement.

# Arguments
- `distance_matrix::Matrix`: Square matrix of pairwise distances between observations
- `ks_to_try::Vector{Int}`: Vector of cluster counts to evaluate (default: 1, 2, and sequence from `Mycelia.ks()`)

# Returns
Named tuple containing:
- `optimal_number_of_clusters::Int`: Best number of clusters found
- `ks_assessed::Vector{Int}`: All cluster counts that were evaluated
- `silhouette_scores::Vector{Float64}`: Silhouette scores for each k assessed
- `hclust_result`: Hierarchical clustering result object

# Details
Uses Ward's linkage method and silhouette scores to evaluate cluster quality. 
Implements an adaptive search that focuses on promising regions between initially tested k values.
For k=1, silhouette score is set to 0 as a special case.
"""
function fit_optimal_number_of_clusters_hclust(distance_matrix, ks_to_try=[1, 2, Mycelia.ks(max=size(distance_matrix, 1))...])
    # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
    # N = size(distance_matrix, 1)
    # ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
    # Int(round(size(distance_matrix, 1)/2))
    # insert!(ks_to_try, 2, Int(round(size(distance_matrix, 1)))))
    # ks_to_try = vcat([2^i for i in 0:Int(floor(log2(size(distance_matrix, 1))))], size(distance_matrix, 1))
    # @info "ks = $(ks_to_try)"
    @show ks_to_try
    
    # can calculate this for k >= 1
    # within_cluster_sum_of_squares = Union{Float64, Missing}[]
    # within_cluster_sum_of_squares = Float64[]
    # these are only valid for k >= 2 so set initial value to missing
    # between_cluster_sum_of_squares = [missing, zeros(length(ks_to_try)-1)...]
    # silhouette_scores = Union{Float64, Missing}[]
    silhouette_scores = Float64[]
    # wcss_scores = Float64[]



    @show "initial heirarchical clustering"
    @time hclust_result = Clustering.hclust(distance_matrix, linkage=:ward, branchorder=:optimal)
    # for k in ks_to_try[1:3]
    for k in ks_to_try
        @show k
        this_clustering = Clustering.cutree(hclust_result, k=k)
        # push!(within_cluster_sum_of_squares, wcss(this_clustering))
        if k == 1
            this_silhouette_score = 0
            # this_wcss = Inf
        else
            this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
            # this_wcss = wcss(this_clustering)
        end
        push!(silhouette_scores, this_silhouette_score)
        # push!(within_cluster_sum_of_squares, this_wcss)
        @show this_silhouette_score
    end      
    optimal_silhouette, optimal_index = findmax(silhouette_scores)
    # while (optimal_index == length(silhouette_scores)) && (optimal_index != length(ks_to_try))
    #     k = ks_to_try[optimal_index+1]
    #     @show k
    #     this_clustering = Clustering.cutree(hclust_result, k=k)
    #     # push!(within_cluster_sum_of_squares, wcss(this_clustering))
    #     this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
    #     push!(silhouette_scores, this_silhouette_score)
    #     @show this_silhouette_score
    #     optimal_silhouette, optimal_index = findmax(silhouette_scores)
    # end
    previous_optimal_number_of_clusters = 0
    optimal_number_of_clusters = ks_to_try[optimal_index]
    done = false
    while optimal_number_of_clusters != previous_optimal_number_of_clusters
        @show optimal_number_of_clusters
        if optimal_index == 1
            window_of_focus = ks_to_try[optimal_index:optimal_index+1]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        elseif optimal_index == length(ks_to_try)
            window_of_focus = ks_to_try[optimal_index-1:optimal_index]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        else
            window_of_focus = ks_to_try[optimal_index-1:optimal_index+1]
        end
        # @show window_of_focus
        midpoints = [
            Int(round(Statistics.mean(window_of_focus[1:2]))),
            Int(round(Statistics.mean(window_of_focus[2:3])))
            ]
        # @show sort(vcat(midpoints, window_of_focus))
        @show midpoints
        
        for k in midpoints
            insertion_index = first(searchsorted(ks_to_try, k))
            if ks_to_try[insertion_index] != k
                insert!(ks_to_try, insertion_index, k)
                @show k
                this_clustering = Clustering.cutree(hclust_result, k=k)
                this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
                insert!(silhouette_scores, insertion_index, this_silhouette_score)
                # this_wcss = wcss(this_clustering)
                # insert!(within_cluster_sum_of_squares, insertion_index, this_wcss)
                @show this_silhouette_score
            end
        end
        
        previous_optimal_number_of_clusters = optimal_number_of_clusters
        optimal_silhouette, optimal_index = findmax(silhouette_scores)
        optimal_number_of_clusters = ks_to_try[optimal_index]
    end
    # @assert length(within_cluster_sum_of_squares) == length(silhouette_scores)
    ks_assessed = ks_to_try[1:length(silhouette_scores)]
    return (;optimal_number_of_clusters, ks_assessed, silhouette_scores, hclust_result)
end