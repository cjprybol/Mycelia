"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function fit_optimal_number_of_clusters(distance_matrix)
    # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
    N = size(distance_matrix, 1)
    ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
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

    for k in ks_to_try[1:3]
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

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function plot_optimal_cluster_assessment_results(clustering_results)
    p1 = StatsPlots.plot(
        ks_assessed[1:length(within_cluster_sum_of_squares)],
        within_cluster_sum_of_squares,
        ylabel = "within cluster sum of squares\n(lower is better)",
        xlabel = "n clusters",
        legend=false
    )
    StatsPlots.vline!(p1, [optimal_number_of_clusters])
    p2 = StatsPlots.plot(
        ks_assessed[1:length(silhouette_scores)],
        silhouette_scores,
        ylabel = "silhouette scores\n(higher is better)",
        xlabel = "n clusters",
        title = "Optimal n clusters = $(optimal_number_of_clusters)",
        legend=false
    )
    StatsPlots.vline!(p2, [optimal_number_of_clusters])
    # TODO write me out
    display(p2)
    StatsPlots.savefig(p1, "$DIR/wcss.svg")
    display(p1)
    StatsPlots.savefig(p2, "$DIR/silhouette.svg")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
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

A short description of the function

```jldoctest
julia> 1 + 1
2
```
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

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function fit_optimal_number_of_clusters_hclust(distance_matrix)
    # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
    N = size(distance_matrix, 1)
    ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
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

    @show "initial heirarchical clustering"
    @time hclust_result = Clustering.hclust(distance_matrix)
    for k in ks_to_try[1:3]
        @show k
        this_clustering = Clustering.cutree(hclust_result, k=k)
        # push!(within_cluster_sum_of_squares, wcss(this_clustering))
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
        this_clustering = Clustering.cutree(hclust_result, k=k)
        # push!(within_cluster_sum_of_squares, wcss(this_clustering))
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
                this_clustering = Clustering.cutree(hclust_result, k=k)
                this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
                insert!(silhouette_scores, insertion_index, this_silhouette_score)
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