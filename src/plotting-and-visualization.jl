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

"""
    plot_taxa_abundances(
        df::DataFrames.DataFrame, 
        taxa_level::String; 
        top_n::Int = 10,
        sample_id_col::String = "sample_id",
        filter_taxa::Union{Vector{Union{String, Missing}}, Nothing} = nothing,
        figure_width::Int = 1500,
        figure_height::Int = 1000,
        bar_width::Float64 = 0.7,
        x_rotation::Int = 45,
        sort_samples::Bool = true,
        color_seed::Union{Int, Nothing} = nothing,
        legend_fontsize::Float64 = 12.0,
        legend_itemsize::Float64 = 12.0,
        legend_padding::Float64 = 5.0,
        legend_rowgap::Float64 = 1.0,
        legend_labelwidth::Union{Nothing, Float64} = nothing,
        legend_titlesize::Float64 = 15.0,
        legend_nbanks::Int = 1
    )

Create a stacked bar chart showing taxa relative abundances for each sample.

# Arguments
- `df`: DataFrame with sample_id and taxonomic assignments at different levels
- `taxa_level`: Taxonomic level to analyze (e.g., "genus", "species")
- `top_n`: Number of top taxa to display individually, remainder grouped as "Other"
- `sample_id_col`: Column name containing sample identifiers
- `filter_taxa`: Taxa to exclude from visualization (default: nothing - no filtering)
- `figure_width`: Width of the figure in pixels
- `figure_height`: Height of the figure in pixels
- `bar_width`: Width of each bar (between 0 and 1)
- `x_rotation`: Rotation angle for x-axis labels in degrees
- `sort_samples`: Whether to sort samples alphabetically
- `color_seed`: Seed for reproducible color generation
- `legend_fontsize`: Font size for legend entries
- `legend_itemsize`: Size of the colored marker/icon in the legend
- `legend_padding`: Padding around legend elements
- `legend_rowgap`: Space between legend rows
- `legend_labelwidth`: Maximum width for legend labels (truncation)
- `legend_titlesize`: Font size for legend title
- `legend_nbanks`: Number of legend columns

# Returns
- `fig`: CairoMakie figure object
- `ax`: CairoMakie axis object
- `taxa_colors`: Dictionary mapping taxa to their assigned colors
"""
function plot_taxa_abundances(
    df::DataFrames.DataFrame, 
    taxa_level::String; 
    top_n::Int = 10,
    sample_id_col::String = "sample_id",
    filter_taxa::Union{Vector{Union{String, Missing}}, Nothing} = nothing,
    figure_width::Int = 1500,
    figure_height::Int = 1000,
    bar_width::Float64 = 0.7,
    x_rotation::Int = 45,
    sort_samples::Bool = true,
    color_seed::Union{Int, Nothing} = nothing,
    legend_fontsize::Float64 = 12.0,
    legend_itemsize::Float64 = 12.0,
    legend_padding::Float64 = 5.0,
    legend_rowgap::Float64 = 1.0,
    legend_labelwidth::Union{Nothing, Float64} = nothing,
    legend_titlesize::Float64 = 15.0,
    legend_nbanks::Int = 1
)
    # Validate inputs
    if !DataFrames.hasproperty(df, Symbol(sample_id_col))
        error("DataFrame must have a column named '$(sample_id_col)'")
    end
    if !DataFrames.hasproperty(df, Symbol(taxa_level))
        error("DataFrame must have a column named '$(taxa_level)'")
    end
    
    # Helper function to check if a value represents "missing" in either format
    is_missing_value = x -> x === missing || x == "missing"
    
    # Step 1: Group by sample_id and count taxa at the specified level
    samples = unique(df[:, sample_id_col])
    n_samples = length(samples)
    
    # Create ordered dictionary to store taxa counts for each sample
    sample_to_taxa_counts = OrderedCollections.OrderedDict{String, Dict{Any, Int}}()
    
    # Count occurrences of each taxon within each sample
    for sample in samples
        sample_data = df[df[:, sample_id_col] .== sample, :]
        sample_to_taxa_counts[string(sample)] = StatsBase.countmap(sample_data[:, taxa_level])
    end
    
    # Step 2: Calculate joint counts across all samples to identify top taxa
    joint_counts = Dict{Any, Int}()
    for (_, counts) in sample_to_taxa_counts
        for (taxon, count) in counts
            joint_counts[taxon] = get(joint_counts, taxon, 0) + count
        end
    end
    
    # Filter out specified taxa if filter_taxa is provided
    if !isnothing(filter_taxa)
        for taxon in filter_taxa
            delete!(joint_counts, taxon)
        end
    end
    
    # Collect all missing-like values (both actual missing and "missing" string)
    missing_keys = filter(is_missing_value, keys(joint_counts))
    missing_count = sum(joint_counts[k] for k in missing_keys)
    
    # Remove all missing-like values from joint_counts for sorting
    for k in missing_keys
        delete!(joint_counts, k)
    end
    
    # Sort taxa by joint counts (descending)
    sorted_taxa = sort(collect(keys(joint_counts)), by=x -> joint_counts[x], rev=true)
    
    # Identify top N taxa, the rest will be grouped as "Other"
    top_taxa = sorted_taxa[1:min(top_n, length(sorted_taxa))]
    
    # Step 3: Normalize counts to get relative abundances for each sample
    sample_to_taxa_relative_abundances = OrderedCollections.OrderedDict{String, Dict{String, Float64}}()
    
    for (sample, counts) in sample_to_taxa_counts
        # Initialize relative abundances
        abundances = Dict{String, Float64}()
        
        # Calculate total counts (excluding filtered taxa)
        if isnothing(filter_taxa)
            total_counts = sum(values(counts))
        else
            total_counts = sum([
                count for (taxon, count) in counts 
                if !(taxon in filter_taxa)
            ])
        end
        
        # Skip samples with no valid counts
        if total_counts == 0
            continue
        end
        
        # Calculate relative abundances for top taxa, "Other", and "Missing"
        other_abundance = 0.0
        missing_abundance = 0.0
        
        for (taxon, count) in counts
            if !isnothing(filter_taxa) && taxon in filter_taxa
                continue
            elseif is_missing_value(taxon)
                missing_abundance += count / total_counts
            elseif taxon in top_taxa
                abundances[string(taxon)] = count / total_counts
            else
                other_abundance += count / total_counts
            end
        end
        
        # Add "Other" category if needed
        if other_abundance > 0
            abundances["Other"] = other_abundance
        end
        
        # Add "Missing" category if needed
        if missing_abundance > 0
            abundances["Missing"] = missing_abundance
        end
        
        # Ensure all top taxa are represented in each sample (even if absent)
        for taxon in top_taxa
            if !haskey(abundances, string(taxon))
                abundances[string(taxon)] = 0.0
            end
        end
        
        sample_to_taxa_relative_abundances[sample] = abundances
    end
    
    # Step 4: Prepare data for visualization
    # Sort samples if requested
    if sort_samples
        samples = sort(collect(keys(sample_to_taxa_relative_abundances)))
    else
        samples = collect(keys(sample_to_taxa_relative_abundances))
    end
    
    # Prepare final set of taxa in desired order (most abundant at bottom)
    # This will hold our ordered taxa categories: top taxa, "Other", and "Missing" (in that order)
    final_taxa = copy(top_taxa)
    
    # Check if any sample has "Other" category
    has_other = any(haskey(abundances, "Other") for (_, abundances) in sample_to_taxa_relative_abundances)
    if has_other
        push!(final_taxa, "Other")
    end
    
    # Check if any sample has "Missing" category
    has_missing_category = any(haskey(abundances, "Missing") for (_, abundances) in sample_to_taxa_relative_abundances)
    if has_missing_category
        push!(final_taxa, "Missing")
    end
    
    # Create matrix for stacked bars
    # Dimensions: [taxa, samples]
    abundance_matrix = zeros(length(final_taxa), length(samples))
    
    # Fill the matrix with relative abundances
    for (col_idx, sample) in enumerate(samples)
        abundances = sample_to_taxa_relative_abundances[sample]
        for (row_idx, taxon) in enumerate(final_taxa)
            taxon_str = string(taxon)
            abundance_matrix[row_idx, col_idx] = get(abundances, taxon_str, 0.0)
        end
    end
    
    # Step 5: Create stacked bar chart visualization with CairoMakie
    
    # Set reproducible colors if seed is provided
    if !isnothing(color_seed)
        Random.seed!(color_seed)
    end
    
    # Generate distinctive colors, reserving specific colors for "Other" and "Missing"
    n_colors_needed = length(final_taxa)
    
    # Generate base colors for all taxa except special categories
    if n_colors_needed > 0
        colorscheme = Colors.distinguishable_colors(n_colors_needed, [Colors.RGB(1,1,1), Colors.RGB(0,0,0)], dropseed=true)
        colorscheme = reverse(colorscheme)  # Reverse to match original behavior
        
        # Use light gray for "Other" and dark gray for "Missing" if they exist
        if has_other
            other_idx = findfirst(x -> x == "Other", final_taxa)
            if !isnothing(other_idx)
                colorscheme[other_idx] = Colors.RGB(0.7, 0.7, 0.7)  # Light gray for "Other"
            end
        end
        
        if has_missing_category
            missing_idx = findfirst(x -> x == "Missing", final_taxa)
            if !isnothing(missing_idx)
                colorscheme[missing_idx] = Colors.RGB(0.4, 0.4, 0.4)  # Dark gray for "Missing"
            end
        end
    else
        colorscheme = Colors.RGB[]
    end
    
    # Create mapping of taxa to colors for legend
    taxa_colors = Dict(taxon => colorscheme[i] for (i, taxon) in enumerate(final_taxa))
    
    # Create the figure with appropriate size ratio for legend
    # Make the figure wider if using more banks
    adjusted_width = figure_width + (legend_nbanks > 1 ? 250 * (legend_nbanks - 1) : 0)
    fig = CairoMakie.Figure(size=(adjusted_width, figure_height))
    
    # Create the main axis for the plot
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="Sample",
        ylabel="Relative Abundance",
        title="$(titlecase(taxa_level)) Relative Abundance (top $(length(top_taxa)) classified)",
        xticks=(1:length(samples), samples),
        xticklabelrotation=x_rotation
    )
    
    # Plot stacked bars - CairoMakie style
    # In CairoMakie, we need to create stacked bars manually
    x_positions = 1:length(samples)
    
    # Start from the bottom of the stack
    previous_heights = zeros(length(samples))
    
    for (i, taxon) in enumerate(final_taxa)
        heights = abundance_matrix[i, :]
        
        # CairoMakie uses 'offset' parameter for stacking
        CairoMakie.barplot!(
            ax,
            x_positions,
            heights,
            offset = previous_heights,
            color = colorscheme[i],
            label = string(taxon),
            width = bar_width
        )
        
        # Update heights for next layer
        previous_heights .+= heights
    end
    
    # Add compact legend with smaller items and text
    fig[1, 2] = CairoMakie.Legend(
        fig,
        ax,
        "Taxa",
        framevisible = true,
        labelsize = legend_fontsize,        # Smaller text
        titlesize = legend_titlesize,       # Smaller title
        patchsize = (legend_itemsize, legend_itemsize),  # Smaller color patches
        padding = legend_padding,           # Less padding around elements
        rowgap = legend_rowgap,             # Less space between rows
        labelwidth = legend_labelwidth,     # Optional truncation of long labels
        nbanks = legend_nbanks              # Multiple columns if needed
    )
    
    # Return the figure, axis, and color mapping for further customization if needed
    return fig, ax, taxa_colors
end

"""
    generate_taxa_abundances_plot(
        joint_reads_to_taxon_lineage_table::DataFrames.DataFrame;
        taxa_level::String,
        top_n::Int = 30,
        kwargs...
    )

Convenience wrapper function to generate taxa abundance visualization
with default parameters and save to a file if requested.

# Arguments
- `joint_reads_to_taxon_lineage_table`: DataFrame with sample_id and taxonomic assignments
- `taxa_level`: Taxonomic level to analyze
- `top_n`: Number of top taxa to display individually
- `kwargs...`: Additional parameters to pass to plot_taxa_abundances

# Returns
- `fig`: CairoMakie figure object
- `ax`: CairoMakie axis object
- `taxa_colors`: Dictionary mapping taxa to their assigned colors
"""
function generate_taxa_abundances_plot(
    joint_reads_to_taxon_lineage_table::DataFrames.DataFrame;
    taxa_level::String,
    top_n::Int = 30,
    save_path::Union{String, Nothing} = nothing,
    adjust_legend_for_taxa_count::Bool = true,
    kwargs...
)
    # Automatically adjust legend parameters based on taxa count
    legend_params = Dict()
    
    if adjust_legend_for_taxa_count && top_n > 30
        # For large numbers of taxa, make legend more compact
        legend_params = Dict(
            :legend_fontsize => 10.0,
            :legend_itemsize => 10.0, 
            :legend_padding => 5.0,
            :legend_rowgap => 0.5,
            :legend_titlesize => 12.0,
            :legend_nbanks => top_n > 50 ? 2 : 1  # Use 2 columns for very large legends
        )
    end
    
    # Merge automatically determined legend parameters with any user-provided ones
    merged_kwargs = merge(legend_params, kwargs)
    
    fig, ax, taxa_colors = plot_taxa_abundances(
        joint_reads_to_taxon_lineage_table,
        taxa_level;
        top_n=top_n,
        merged_kwargs...
    )
    
    # Save figure if path is provided
    if !isnothing(save_path)
        CairoMakie.save(save_path, fig)
    end
    
    return fig, ax, taxa_colors
end

function visualize_xam_classifications(;xam, sample_id, accession_to_taxid_table)
    xam_classification_table_file = xam * ".classification_table.arrow"
    if !isfile(xam_classification_table_file)
        @time classification_table = Mycelia.classify_xam_with_blast_taxonomies(xam=xam, accession_to_taxid_table=accession_to_taxid_table)
        Arrow.write(xam_classification_table_file, classification_table)
    else
        classification_table = DataFrames.DataFrame(Arrow.Table(xam_classification_table_file))
    end

    unique_taxids = filter(x -> x > 0, sort(unique(classification_table[!, "final_assignment"])))
    taxid_lineage_table = Mycelia.taxids2taxonkit_summarized_lineage_table(unique_taxids)
    reads_to_taxon_lineage_table = DataFrames.leftjoin(classification_table, taxid_lineage_table, on="final_assignment" => "taxid")
    reads_to_taxon_lineage_table = Mycelia.assign_lowest_rank_to_reads_to_taxon_lineage_table(reads_to_taxon_lineage_table)
    p = Mycelia.sankey_visualize_reads_lowest_rank(reads_to_taxon_lineage_table[!, "lowest_rank"]; title="$sample_id")
    StatsPlots.savefig(p, xam * ".$(sample_id).sankey.png")
    StatsPlots.savefig(p, xam * ".$(sample_id).sankey.svg")
    StatsPlots.savefig(p, xam * ".$(sample_id).sankey.pdf")
    display(p)

    for taxa_level in ["domain", "family", "genus", "species"]
        taxa_counts = sort(collect(StatsBase.countmap(filter(!ismissing, reads_to_taxon_lineage_table[!, taxa_level]))), by=x->x[2], rev=true)
        taxa_count_pairs = Mycelia.streamline_counts(taxa_counts)

        fig = Mycelia.visualize_single_sample_taxa_count_pairs(taxa_count_pairs; title = "Distribution of Relative Abundance by Taxa:$(taxa_level) - sample_id:$(sample_id)")
        CairoMakie.save(xam * ".$(sample_id).$(taxa_level).counts.png", fig)
        CairoMakie.save(xam * ".$(sample_id).$(taxa_level).counts.svg", fig)
        CairoMakie.save(xam * ".$(sample_id).$(taxa_level).counts.pdf", fig)
        display(fig)

        fig = Mycelia.visualize_single_sample_taxa_count_pair_proportions(taxa_count_pairs, title = "Distribution of Relative Abundance by Taxa:$(taxa_level) - sample_id:$(sample_id)")
        CairoMakie.save(xam * ".$(sample_id).$(taxa_level).proportion.png", fig)
        CairoMakie.save(xam * ".$(sample_id).$(taxa_level).proportion.svg", fig)
        CairoMakie.save(xam * ".$(sample_id).$(taxa_level).proportion.pdf", fig)
        display(fig)
    end
end

function visualize_single_sample_taxa_count_pair_proportions(taxa_count_pairs; title = "Distribution of Relative Abundance by Taxa")
    df = DataFrames.DataFrame(
        taxa = first.(taxa_count_pairs),
        count = last.(taxa_count_pairs)
    )
    df[!, "proportion"] .= round.(df[!, "count"] ./ sum(df[!, "count"]), digits=3)

    # Create the plot
    fig = CairoMakie.Figure(size=(800, 500), fontsize=12)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Taxa",
        ylabel = "Proportion",
        title = title,
        xticks = (1:length(df.taxa), df.taxa),
        xticklabelrotation = π/2,
        xticklabelsize = max(15 - Int(floor(log(2, length(df.taxa)))), 6)
    )

    # Create bar plot
    bars = CairoMakie.barplot!(
        ax,
        1:length(df.taxa),
        df.proportion,
        # bar_labels = :y,
        # label_formatter = x -> string(x),
        strokewidth = 1,
        strokecolor = :black
    )

    # Add a bit of styling
    CairoMakie.ylims!(ax, 0, maximum(df.proportion) * 1.2)

    # Add count values on top of each bar as vertical labels
    for (i, proportion) in enumerate(df.proportion)
        CairoMakie.text!(
            ax,
            i,
            proportion + maximum(df.proportion) * 0.02,
            text = string(proportion),
            align = (:left, :center),
            fontsize = max(15 - Int(floor(log(2, length(df.taxa)))), 6),
            rotation = π/2  # Make the labels vertical (90 degrees)
        )
    end

    # Add grid lines for better readability
    ax.ygridvisible = true
    ax.ygridstyle = :dash

    # Add a subtle background color
    ax.backgroundcolor = (:grey90, 0.1)

    # Save the figure


    # Display the figure
    fig
end


function visualize_single_sample_taxa_count_pairs(taxa_count_pairs; title = "Distribution of Counts by Taxa")
    # Convert to DataFrame for easier plotting
    df = DataFrames.DataFrame(
        taxa = first.(taxa_count_pairs),
        count = last.(taxa_count_pairs)
    )

    # Create the plot
    fig = CairoMakie.Figure(size=(800, 500), fontsize=12)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel = "Taxa",
        ylabel = "Count",
        title = title,
        xticks = (1:length(df.taxa), df.taxa),
        xticklabelrotation = π/2,
        xticklabelsize = max(15 - Int(floor(log(2, length(df.taxa)))), 6)
    )

    # Create bar plot
    bars = CairoMakie.barplot!(
        ax,
        1:length(df.taxa),
        df.count,
        strokewidth = 1,
        strokecolor = :black
    )

    # Add a bit of styling
    CairoMakie.ylims!(ax, 0, maximum(df.count) * 1.2)

    # Add count values on top of each bar as vertical labels
    for (i, count) in enumerate(df.count)
        CairoMakie.text!(
            ax,
            i,
            count + maximum(df.count) * 0.02,
            text = string(count),
            align = (:left, :center),
            fontsize = max(15 - Int(floor(log(2, length(df.taxa)))), 6),
            rotation = π/2  # Make the labels vertical (90 degrees)
        )
    end

    # Add grid lines for better readability
    ax.ygridvisible = true
    ax.ygridstyle = :dash

    # Add a subtle background color
    ax.backgroundcolor = (:grey90, 0.1)

    # Display the figure
    fig
end

function sankey_visualize_reads_lowest_rank(lowest_ranks; title="")
    level_hits = [sum(lowest_ranks .>= i) for i in 1:9]
    connections = [(a, b, c) for (a, b, c) in zip(1:10, 2:11, level_hits)]
    push!(connections, (1, 11, sum(lowest_ranks .== 0)))
    labels = [
        "all reads",
        "domain",
        "realm",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "unclassified"
    ]

    total_count = length(lowest_ranks)

    counts = vcat(total_count, [c[3] for c in connections])

    percentages = round.(counts ./ total_count * 100, digits=1)

    labels = [label * "\n$(count)\n$(percentage)%" for (label, count, percentage) in zip(labels, counts, percentages)]

    p = SankeyPlots.sankey(
        title = title,
        [c[1] for c in connections],
        [c[2] for c in connections],
        [c[3] for c in connections],
        node_labels=labels,
        edge_color=:gradient,
        label_position=:bottom,
        label_size=7,
        compact=true,
        size=(1200,800),
        dpi=100
    )
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Plots a k-mer rarefaction curve from data stored in a TSV file.
The TSV file should contain two columns:
1. Number of FASTA files processed.
2. Cumulative unique k-mers observed at that point.

The plot is displayed and saved in PNG, PDF, and SVG formats.

# Arguments
- `rarefaction_data_path::AbstractString`: Path to the TSV file containing rarefaction data.
- `output_dir::AbstractString`: Directory where the output plots will be saved. Defaults to the directory of `rarefaction_data_path`.
- `output_basename::AbstractString`: Basename for the output plot files (without extension). Defaults to the basename of `rarefaction_data_path` without its original extension.
- `display_plot::Bool`: Whether to display the plot interactively. Defaults to `true`.

# Keyword Arguments
- `fig_size::Tuple{Int, Int}`: Size of the output figure, e.g., `(1000, 750)`.
- `title::AbstractString`: Title of the plot.
- `xlabel::AbstractString`: Label for the x-axis.
- `ylabel::AbstractString`: Label for the y-axis.
- `line_color`: Color of the plotted line.
- `line_style`: Style of the plotted line (e.g. `:dash`, `:dot`).
- `marker`: Marker style for points (e.g. `:circle`, `:xcross`).
- `markersize::Number`: Size of the markers.
- Any other keyword arguments will be passed to `Makie.Axis`.
"""
function plot_kmer_rarefaction(
    rarefaction_data_path::AbstractString;
    output_dir::AbstractString = dirname(rarefaction_data_path),
    output_basename::AbstractString = first(Base.Filesystem.splitext(basename(rarefaction_data_path))),
    display_plot::Bool = true,
    # Makie specific customizations
    fig_size::Tuple{Int, Int} = (1000, 750),
    title::AbstractString = "K-mer Rarefaction Curve",
    xlabel::AbstractString = "Number of FASTA Files Processed",
    ylabel::AbstractString = "Cumulative Unique K-mers",
    line_color = :blue,
    line_style = nothing,
    marker = nothing,
    markersize::Number = 10,
    axis_kwargs... # Capture other axis properties
)
    if !isfile(rarefaction_data_path)
        Base.@error "Rarefaction data file not found: $rarefaction_data_path"
        return nothing
    end

    data = try
        DelimitedFiles.readdlm(rarefaction_data_path, '\t', Int, header=false)
    catch e
        Base.@error "Failed to read rarefaction data from $rarefaction_data_path: $e"
        return nothing
    end

    if size(data, 2) != 2
        Base.@error "Rarefaction data file $rarefaction_data_path must have exactly two columns."
        return nothing
    end

    files_processed = data[:, 1]
    unique_kmers = data[:, 2]

    # Sort data by files_processed for a proper line plot,
    # as the input might be from unordered parallel processing.
    sort_indices = sortperm(files_processed)
    files_processed = files_processed[sort_indices]
    unique_kmers = unique_kmers[sort_indices]
    
    Base.mkpath(output_dir) # Ensure output directory exists

    fig = Makie.Figure(size = fig_size)
    ax = Makie.Axis(
        fig[1, 1],
        title = title,
        xlabel = xlabel,
        ylabel = ylabel;
        axis_kwargs... # Pass through other axis settings
    )

    Makie.lines!(ax, files_processed, unique_kmers, color = line_color, linestyle = line_style)
    if !isnothing(marker)
        Makie.scatter!(ax, files_processed, unique_kmers, color = line_color, marker = marker, markersize = markersize)
    end


    if display_plot
        Base.@info "Displaying rarefaction plot..."
        Makie.display(fig)
    end

    output_path_png = Base.Filesystem.joinpath(output_dir, output_basename * ".png")
    output_path_pdf = Base.Filesystem.joinpath(output_dir, output_basename * ".pdf")
    output_path_svg = Base.Filesystem.joinpath(output_dir, output_basename * ".svg")

    try
        Base.@info "Saving plot to $output_path_png"
        Makie.save(output_path_png, fig)
        Base.@info "Saving plot to $output_path_pdf"
        Makie.save(output_path_pdf, fig)
        Base.@info "Saving plot to $output_path_svg"
        Makie.save(output_path_svg, fig)
    catch e
        Base.@error "Failed to save plot: $e. Make sure a Makie backend (e.g., CairoMakie for saving, GLMakie for display) is active and correctly configured in your environment."
    end
    
    return fig # Return the figure object
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the optimal subsequence length based on error rate distribution.

# Arguments
- `error_rate`: Single error rate or array of error rates (between 0 and 1)
- `threshold`: Desired probability that a subsequence is error-free (default: 0.95)
- `sequence_length`: Maximum sequence length to consider for plotting
- `plot_result`: If true, returns a plot of probability vs. length

# Returns
- If `plot_result=false`: Integer representing optimal subsequence length
- If `plot_result=true`: Tuple of (optimal_length, plot)

# Examples
```julia
# Single error rate
optimal_subsequence_length(error_rate=0.01)

# Array of error rates
optimal_subsequence_length(error_rate=[0.01, 0.02, 0.01])

# With more stringent threshold
optimal_subsequence_length(error_rate=0.01, threshold=0.99)

# Generate plot
length, p = optimal_subsequence_length(error_rate=0.01, plot_result=true)
Plots.display(p)
```
"""
function optimal_subsequence_length(;error_rate::Union{Real, AbstractArray{<:Real}},
                                   threshold::Float64=0.95,
                                   sequence_length::Union{Nothing, Int}=nothing,
                                   plot_result::Bool=false)
    # Handle array input by calculating mean error rate
    avg_error_rate = isa(error_rate, AbstractArray) ? Statistics.mean(error_rate) : error_rate
    
    # Validate inputs
    if avg_error_rate <= 0
        optimal_length = typemax(Int)
    elseif avg_error_rate >= 1
        optimal_length = 1
    else
        # Calculate optimal length where P(error-free) >= threshold
        optimal_length = floor(Int, log(threshold) / log(1 - avg_error_rate))
        optimal_length = max(1, optimal_length)  # Ensure at least length 1
    end

    # Return early if no plot requested
    if !plot_result
        return optimal_length
    end
    
    # For plotting, determine sequence length to display
    max_length = isnothing(sequence_length) ? 2 * optimal_length : sequence_length
    
    # Calculate probabilities for different lengths
    lengths = 1:max_length
    probabilities = [(1 - avg_error_rate)^len for len in lengths]
    
    # Create DataFrame for plotting
    df = DataFrames.DataFrame(
        Length = collect(lengths),
        Probability = probabilities,
        Optimal = lengths .== optimal_length
    )

    quality_score = Mycelia.error_rate_to_q_value(error_rate)
    rounded_quality_score = Int(floor(quality_score))
    rounded_error_rate = round(avg_error_rate, digits=4)
    rounded_threshold_rate = round(threshold, digits=2)
    plot_title = 
    """
    Optimal Kmer Length Inference
    Error Rate: $(rounded_error_rate * 100)%≈Q$(rounded_quality_score)
    Threshold: $(rounded_threshold_rate * 100)% of kmers expected to be correct
    """
    
    # Create plot
    p = Plots.plot(
        df.Length, df.Probability,
        linewidth=2, 
        label="P(error-free)",
        xlabel="Subsequence Length",
        ylabel="Probability of Error-free Match",
        title=plot_title,
        grid=true,
        ylims=(0,1),
        alpha=0.8
    )
    
    # Add horizontal line for threshold
    Plots.hline!([threshold], linestyle=:dash, color=:red, label="Threshold")
    
    # Add vertical line for optimal length
    Plots.vline!([optimal_length], linestyle=:dash, color=:green, label="Optimal length: $optimal_length")
    
    # Highlight optimal point
    Plots.scatter!([optimal_length], [threshold], color=:orange, markersize=8, label="")
    
    return optimal_length, p
end