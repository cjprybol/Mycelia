function visualize_many_timeseries(time_series_data::Vector{Vector{Float64}};
                                  title::String="High-Density Time Series Visualization")
    # Create x-axis values (assuming equal length for all series)
    # If series have different lengths, you'll need to adjust this
    n_points = length(time_series_data[1])
    x_values = 1:n_points
    
    # Calculate statistics for the data
    all_values = vcat(time_series_data...)
    global_min, global_max = minimum(all_values), maximum(all_values)
    
    # Calculate mean and quantiles
    means = zeros(n_points)
    q10 = zeros(n_points)
    q25 = zeros(n_points)
    q75 = zeros(n_points)
    q90 = zeros(n_points)
    
    for i in 1:n_points
        point_values = [series[i] for series in time_series_data if length(series) >= i]
        means[i] = Statistics.mean(point_values)
        q10[i] = Statistics.quantile(point_values, 0.1)
        q25[i] = Statistics.quantile(point_values, 0.25)
        q75[i] = Statistics.quantile(point_values, 0.75)
        q90[i] = Statistics.quantile(point_values, 0.9)
    end
    
    # Setup the figure - use size instead of resolution
    fig = CairoMakie.Figure(size=(1200, 800), fontsize=18)
    
    # First visualization: Density heatmap approach
    ax1 = CairoMakie.Axis(fig[1, 1], 
               title="Density Heatmap", 
               xlabel="Time", 
               ylabel="Value",
               titlesize=22)
    
    # Create density heatmap - use CairoMakie.heatmap! instead of density_heatmap!
    # First prepare data for the heatmap - extract x and y coordinates separately
    x_coords = Float64[]
    y_coords = Float64[]
    
    for (series_idx, series) in enumerate(time_series_data)
        for (point_idx, value) in enumerate(series)
            push!(x_coords, point_idx)
            push!(y_coords, value)
        end
    end
    
    # Create 2D histogram manually
    bins_x = 100  # Number of bins in x direction
    bins_y = 100  # Number of bins in y direction
    
    x_range = (minimum(x_coords), maximum(x_coords))
    y_range = (minimum(y_coords), maximum(y_coords))
    
    x_edges = range(x_range[1], x_range[2], length=bins_x+1)
    y_edges = range(y_range[1], y_range[2], length=bins_y+1)
    
    histogram = zeros(bins_x, bins_y)
    
    for i in 1:length(x_coords)
        x, y = x_coords[i], y_coords[i]
        x_bin = max(1, min(bins_x, Int(floor((x - x_range[1]) / (x_range[2] - x_range[1]) * bins_x)) + 1))
        y_bin = max(1, min(bins_y, Int(floor((y - y_range[1]) / (y_range[2] - y_range[1]) * bins_y)) + 1))
        histogram[x_bin, y_bin] += 1
    end
    
    # Create heatmap
    hmap = CairoMakie.heatmap!(ax1, 
                    x_edges[1:end-1], 
                    y_edges[1:end-1], 
                    histogram, 
                    colormap=:viridis)
    
    # Add mean line
    CairoMakie.lines!(ax1, x_values, means, color=:white, linewidth=3)
    
    CairoMakie.Colorbar(fig[1, 2], hmap, label="Density")
    
    # Second visualization: Sample with transparency
    ax2 = CairoMakie.Axis(fig[2, 1], 
               title="Sample with Statistical Bands", 
               xlabel="Time", 
               ylabel="Value",
               titlesize=22)
    
    # Draw a random sample of the time series with high transparency
    sample_size = min(500, length(time_series_data))
    sample_indices = rand(1:length(time_series_data), sample_size)
    
    for idx in sample_indices
        CairoMakie.lines!(ax2, x_values[1:length(time_series_data[idx])], time_series_data[idx], 
               color=(ColorSchemes.viridis[rand()], 0.05), linewidth=0.5)
    end
    
    # Add statistical bands
    CairoMakie.band!(ax2, x_values, q10, q90, color=(:blue, 0.2))
    CairoMakie.band!(ax2, x_values, q25, q75, color=(:blue, 0.3))
    
    # Add mean line
    CairoMakie.lines!(ax2, x_values, means, color=:red, linewidth=3)
    
    # Add legend
    CairoMakie.Legend(fig[2, 2],
           [
            CairoMakie.LineElement(color=:red, linewidth=3),
            CairoMakie.PolyElement(color=(:blue, 0.2)),
            CairoMakie.PolyElement(color=(:blue, 0.3))
           ],
           ["Mean", "10-90% Range", "25-75% Range"]
    )
    
    # Link the axes
    CairoMakie.linkaxes!(ax1, ax2)
    
    # Add overall title
    CairoMakie.Label(fig[0, :], text=title, fontsize=26)
    
    return fig
end

# src/visualization.jl
"""
Plot embeddings with optional true and fitted cluster labels.

# Arguments
- `embeddings::Matrix{<:Real}`: 2D embedding matrix where each column is a data point
- `title::String`: Title of the plot
- `xlabel::String`: Label for the x-axis
- `ylabel::String`: Label for the y-axis
- `true_labels::Vector{<:Integer}`: Vector of true cluster labels (optional)
- `fit_labels::Vector{<:Integer}`: Vector of fitted cluster labels (optional)

# Returns
- `Plots.Plot`: Plot object that can be displayed or saved
"""
function plot_embeddings(embeddings; title="", xlabel="", ylabel="", true_labels=nothing, fit_labels=nothing)
    scatter(embeddings[1, :], embeddings[2, :],
           title=title,
           xlabel=xlabel,
           ylabel=ylabel,
           label="",
           legend=:topright)

    if true_labels !== nothing
        for i in unique(true_labels)
            idx = findall(x -> x == i, true_labels)
            scatter!(embeddings[1, idx], embeddings[2, idx],
                     label="True Cluster $i",
                     markershape=:star5,
                     legend=:topright)
        end
    end

    if fit_labels !== nothing
        for i in unique(fit_labels)
            idx = findall(x -> x == i, fit_labels)
            scatter!(embeddings[1, idx], embeddings[2, idx],
                     label="Fit Cluster $i",
                     legend=:bottomright)
        end
    end

    return plot!
end