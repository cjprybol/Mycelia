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

"""
Plot embeddings with optional true and fitted cluster labels using Makie.jl, 
with legend outside, and color by fit labels, shape by true labels.

# Arguments
- `embeddings::Matrix{<:Real}`: 2D embedding matrix where each column is a data point
- `title::String`: Title of the plot
- `xlabel::String`: Label for the x-axis
- `ylabel::String`: Label for the y-axis
- `true_labels::Vector{<:Integer}`: Vector of true cluster labels (optional)
- `fit_labels::Vector{<:Integer}`: Vector of fitted cluster labels (optional)

# Returns
- `Makie.Figure`: Figure object that can be displayed or saved
"""
function plot_embeddings(embeddings; title="", xlabel="", ylabel="", true_labels=nothing, fit_labels=nothing)
    fig = CairoMakie.Figure(size=(900, 450))
    ax = CairoMakie.Axis(fig[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    npoints = size(embeddings, 2)

    fit_palette = fit_labels !== nothing ? Mycelia.n_maximally_distinguishable_colors(length(unique(fit_labels))) : [:gray]
    true_markers = true_labels !== nothing ? Mycelia.choose_top_n_markers(length(unique(true_labels))) : [:circle]
    fit_label_to_color = fit_labels !== nothing ? Dict(lbl => fit_palette[i] for (i, lbl) in enumerate(unique(fit_labels))) : Dict()
    true_label_to_marker = true_labels !== nothing ? Dict(lbl => true_markers[i] for (i, lbl) in enumerate(unique(true_labels))) : Dict()

    xs = embeddings[1, :]
    ys = embeddings[2, :]
    colors = [fit_labels !== nothing ? fit_label_to_color[fit_labels[i]] : :gray for i in 1:npoints]
    markers = [true_labels !== nothing ? true_label_to_marker[true_labels[i]] : :circle for i in 1:npoints]
    CairoMakie.scatter!(ax, xs, ys; color=colors, marker=markers, markersize=14, strokewidth=1, strokecolor=:black, label="")

    # Legend entries (as before) ...
    if fit_labels !== nothing
        for (i, lbl) in enumerate(unique(fit_labels))
            CairoMakie.scatter!(ax, [NaN], [NaN]; color=fit_palette[i], marker=:circle, markersize=14, label="Fit Cluster $lbl")
        end
    end
    if true_labels !== nothing
        for (i, lbl) in enumerate(unique(true_labels))
            CairoMakie.scatter!(ax, [NaN], [NaN]; color=:gray, marker=true_markers[i], markersize=14, label="True Cluster $lbl")
        end
    end

    # Place legend outside
    CairoMakie.axislegend(ax; position=(1.25, 0.5), nbanks=1)
    
    # Print ranges for debugging
    println("x range: ", minimum(xs), " to ", maximum(xs))
    println("y range: ", minimum(ys), " to ", maximum(ys))

    return fig
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add random noise to create a vector of jittered values.

Generates `n` values by adding random noise to the input value `x`. 
The noise is uniformly distributed between -1/3 and 1/3.

# Arguments
- `x`: Base value to add jitter to
- `n`: Number of jittered values to generate

# Returns
- Vector of length `n` containing jittered values around `x`
"""
function jitter(x, n)
    return [x + rand() / 3 * (ifelse(rand(Bool), 1, -1)) for i in 1:n]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a visualization of a genome assembly graph using Bandage.

# Arguments
- `gfa`: Path to input GFA (Graphical Fragment Assembly) file
- `img`: Optional output image path. Defaults to GFA filename with .png extension

# Returns
- Path to the generated image file
"""
function bandage_visualize(;gfa, img=gfa*".png")
    # run(`$(bandage) image --helpall`)
    bandage = Mycelia.download_bandage()
    if !isfile(img)
        run(`$(bandage) image $(gfa) $(img)`)
    end
    return img
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate `n` colors that are maximally distinguishable from each other.

# Arguments
- `n::Integer`: The number of distinct colors to generate

# Returns
A vector of `n` RGB colors that are optimized for maximum perceptual distinction,
using white (RGB(1,1,1)) and black (RGB(0,0,0)) as anchor colors.
"""
function n_maximally_distinguishable_colors(n)
    return Colors.distinguishable_colors(n, [Colors.RGB(1,1,1), Colors.RGB(0,0,0)], dropseed=true)
end

# https://www.giantfocal.com/toolkit/font-size-converter
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert pixel measurements to point measurements using the standard 4:3 ratio.

Points are the standard unit for typography (1 point = 1/72 inch), while pixels are 
used for screen measurements. This conversion uses the conventional 4:3 ratio where 
3 points equal 4 pixels.

# Arguments
- `pixels`: The number of pixels to convert

# Returns
- The equivalent measurement in points
"""
function pixels_to_points(pixels)
    return pixels / 4 * 3
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert typographic points to pixels using a 4:3 ratio (1 point = 4/3 pixels).

# Arguments
- `points`: Size in typographic points (pt)

# Returns
- Size in pixels (px)
"""
function points_to_pixels(points)
    return points / 3 * 4
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Draw a dendrogram visualization of hierarchical clustering results stored in a MetaDiGraph.

# Arguments
- `mg::MetaGraphs.MetaDiGraph`: Graph containing hierarchical clustering results. Must have
  `:hcl` in graph properties with clustering data and vertex properties containing `:x`, `:y` coordinates.

# Keywords
- `width::Integer=500`: Width of output image in pixels
- `height::Integer=500`: Height of output image in pixels 
- `fontsize::Integer=12`: Font size for node labels in points
- `margins::Float64`: Margin size in pixels, defaults to min(width,height)/20
- `mergenodesize::Float64=1`: Size of circular nodes at merge points
- `lineweight::Float64=1`: Thickness of dendrogram lines
- `filename::String`: Output filename, defaults to timestamp with .dendrogram.png extension

# Returns
Nothing, but saves dendrogram image to disk and displays preview.
"""
function draw_dendrogram_tree(
        mg::MetaGraphs.MetaDiGraph;
        width=500,
        height=500,
        fontsize=12,
        # margins=min(width, height)/25,
        margins=min(width, height)/20,
        # margins=20,
        mergenodesize=1,
        lineweight=1,
        filename=Dates.format(Dates.now(), "yyyymmddTHHMMSS") * ".dendrogram.png"
    )

    fontsizebuffer = points_to_pixels(fontsize) / 3
    available_width = width - 2 * margins
    available_height = height - 2 * margins

    leaf_nodes = mg.gprops[:hcl].labels

    # Create a new drawing
    Luxor.Drawing(width, height, filename)
    # Luxor.origin()  # Set the origin to the center of the drawing
    # Luxor.background("white")
    Luxor.sethue("black")
    Luxor.fontsize(fontsize)
    Luxor.setline(lineweight)

    for ordered_leaf_node in mg.gprops[:hcl].order
        x = mg.vprops[ordered_leaf_node][:x] * available_width + margins
        y = mg.vprops[ordered_leaf_node][:y] * available_height + margins
        Luxor.circle(x, y, mergenodesize, :fill)
        # Luxor.text(string(ordered_leaf_node), Luxor.Point(x, y), halign=:center, valign=:middle)
        Luxor.text(string(ordered_leaf_node), Luxor.Point(x, y + fontsizebuffer), halign=:center, valign=:top)
    end
    for (i, (left, right)) in enumerate(eachrow(mg.gprops[:hcl].merges))
        parent_vertex = mg[string(i), :hclust_id]
                                        
        x = mg.vprops[parent_vertex][:x] * available_width + margins
        y = mg.vprops[parent_vertex][:y] * available_height + margins
        # Luxor.text(string(ordered_leaf_node), Luxor.Point(x, y), halign=:center, valign=:middle)
        Luxor.circle(x, y, mergenodesize, :fill)

        left_child_vertex = mg[string(left), :hclust_id]
        left_child_x = mg.vprops[left_child_vertex][:x] * available_width + margins
        left_child_y = mg.vprops[left_child_vertex][:y] * available_height + margins
        # draw horizontal bar
        Luxor.line(Luxor.Point(left_child_x, y), Luxor.Point(x, y), action=:stroke)
        # draw vertical bar
        Luxor.line(Luxor.Point(left_child_x, left_child_y), Luxor.Point(left_child_x, y), action=:stroke)

        right_child_vertex = mg[string(right), :hclust_id]
        right_child_x = mg.vprops[right_child_vertex][:x] * available_width + margins
        right_child_y = mg.vprops[right_child_vertex][:y] * available_height + margins
        # draw horizontal bar
        Luxor.line(Luxor.Point(x, y), Luxor.Point(right_child_x, y), action=:stroke)
        # draw vertical bar
        Luxor.line(Luxor.Point(right_child_x, right_child_y), Luxor.Point(right_child_x, y), action=:stroke)
    end

    # Finish the drawing and save the file
    Luxor.finish()
    Luxor.preview()
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Draw a radial hierarchical clustering tree visualization and save it as an image file.

# Arguments
- `mg::MetaGraphs.MetaDiGraph`: A meta directed graph containing hierarchical clustering data
  with required graph properties `:hcl` containing clustering information.

# Keywords
- `width::Int=500`: Width of the output image in pixels
- `height::Int=500`: Height of the output image in pixels
- `fontsize::Int=12`: Font size for node labels
- `margins::Float64`: Margin size (automatically calculated as min(width,height)/20)
- `mergenodesize::Float64=1`: Size of the merge point nodes
- `lineweight::Float64=1`: Thickness of the connecting lines
- `filename::String`: Output filename (defaults to timestamp with ".radial.png" suffix)

# Details
The function creates a radial visualization of hierarchical clustering results where:
- Leaf nodes are arranged in a circle with labels
- Internal nodes represent merge points
- Connections show the hierarchical structure through arcs and lines

The visualization is saved as a PNG file and automatically previewed.

# Required Graph Properties
The input graph must have:
- `mg.gprops[:hcl].labels`: Vector of leaf node labels
- `mg.gprops[:hcl].order`: Vector of ordered leaf nodes
- `mg.gprops[:hcl].merges`: Matrix of merge operations
- `mg.vprops[v][:x]`: X coordinate for each vertex
- `mg.vprops[v][:y]`: Y coordinate for each vertex
"""
function draw_radial_tree(
        mg::MetaGraphs.MetaDiGraph;
        width=500,
        height=500,
        fontsize=12,
        # margins=min(width, height)/25,
        margins=min(width, height)/20,
        # margins=20,
        mergenodesize=1,
        lineweight=1,
        filename=Dates.format(Dates.now(), "yyyymmddTHHMMSS") * ".radial.png"
    )

    # check max label size against margin and update as necessary
    max_radius = (min(width, height) - (margins * 2)) / 2

    # fontsizebuffer = points_to_pixels(fontsize) / 3
    available_width = width - 2 * margins
    available_height = height - 2 * margins
    leaf_nodes = mg.gprops[:hcl].labels

    # Create a new drawing
    Luxor.Drawing(width, height, filename)
    Luxor.origin()  # Set the origin to the center of the drawing
    # Luxor.background("white")
    Luxor.sethue("black")
    Luxor.fontsize(fontsize)
    Luxor.setline(lineweight)

    for ordered_leaf_node in mg.gprops[:hcl].order
        original_x = mg.vprops[ordered_leaf_node][:x]
        original_y = mg.vprops[ordered_leaf_node][:y]
        polar_radian = original_x * 2 * MathConstants.pi
        adjusted_radius = original_y * max_radius
        x = cos(polar_radian) * adjusted_radius
        y = sin(polar_radian) * adjusted_radius
        Luxor.circle(x, y, mergenodesize, :fill)
        text_x = cos(polar_radian) * (max_radius + 10)
        text_y = sin(polar_radian) * (max_radius + 10)
        Luxor.text(string(ordered_leaf_node), Luxor.Point(text_x, text_y), angle=polar_radian, halign=:left, valign=:middle)
    end
    for (i, (left, right)) in enumerate(eachrow(mg.gprops[:hcl].merges))
        parent_vertex = mg[string(i), :hclust_id]

        original_x = mg.vprops[parent_vertex][:x]
        original_y = mg.vprops[parent_vertex][:y]

        polar_radian = original_x * 2 * MathConstants.pi
        adjusted_radius = original_y * max_radius
        x = cos(polar_radian) * adjusted_radius
        y = sin(polar_radian) * adjusted_radius
        Luxor.circle(x, y, mergenodesize, :fill)

        # left child
        left_child_vertex = mg[string(left), :hclust_id]
        left_child_original_x = mg.vprops[left_child_vertex][:x]
        left_child_original_y = mg.vprops[left_child_vertex][:y]
        left_child_polar_radian = left_child_original_x * 2 * MathConstants.pi
        left_child_adjusted_radius = left_child_original_y * max_radius
        left_child_x = cos(left_child_polar_radian) * left_child_adjusted_radius
        left_child_y = sin(left_child_polar_radian) * left_child_adjusted_radius
        # # draw horizontal bar
        Luxor.arc(Luxor.Point(0, 0), adjusted_radius, left_child_polar_radian, polar_radian, action=:stroke)
        # draw vertical bar
        join_x = cos(left_child_polar_radian) * adjusted_radius
        join_y = sin(left_child_polar_radian) * adjusted_radius
        Luxor.line(Luxor.Point(left_child_x, left_child_y), Luxor.Point(join_x, join_y), action=:stroke)

        # right child
        right_child_vertex = mg[string(right), :hclust_id]
        right_child_original_x = mg.vprops[right_child_vertex][:x]
        right_child_original_y = mg.vprops[right_child_vertex][:y]
        right_child_polar_radian = right_child_original_x * 2 * MathConstants.pi
        right_child_adjusted_radius = right_child_original_y * max_radius
        right_child_x = cos(right_child_polar_radian) * right_child_adjusted_radius
        right_child_y = sin(right_child_polar_radian) * right_child_adjusted_radius
        # # draw horizontal bar
        Luxor.arc(Luxor.Point(0, 0), adjusted_radius, polar_radian, right_child_polar_radian, action=:stroke)
        # draw vertical bar
        join_x = cos(right_child_polar_radian) * adjusted_radius
        join_y = sin(right_child_polar_radian) * adjusted_radius
        Luxor.line(Luxor.Point(right_child_x, right_child_y), Luxor.Point(join_x, join_y), action=:stroke)
    end

    # Finish the drawing and save the file
    Luxor.finish()
    Luxor.preview()
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a visualization of chromosome coverage data with statistical thresholds.

# Arguments
- `cdf::DataFrame`: Coverage data frame containing columns:
  - `index`: Chromosome position indices
  - `depth`: Coverage depth values
  - `chromosome`: Chromosome identifier
  - `mean_coverage`: Mean coverage value
  - `std_coverage`: Standard deviation of coverage
  - `3σ`: Boolean vector indicating +3 sigma regions
  - `-3σ`: Boolean vector indicating -3 sigma regions

# Returns
- A StatsPlots plot object showing:
  - Raw coverage data (black line)
  - Mean coverage and ±1,2,3σ thresholds (rainbow colors)
  - Highlighted regions exceeding ±3σ thresholds (red vertical lines)
"""
function chromosome_coverage_table_to_plot(cdf)
    p = StatsPlots.plot(
        xlims = extrema(cdf[!, "index"]),
        ylims=(1, maximum(cdf[!, "depth"]) * 1.1),
        title = cdf[1, "chromosome"],
        xlabel = "chromosome index",
        ylabel = "depth"
    )
    for r in find_true_ranges(cdf[!, "3σ"]; min_length=1000)
        range_mean = Statistics.mean(cdf[r[1]:r[2], :depth])
        StatsPlots.vline!(p,
            [r[1], r[2]],
            # [range_mean, range_mean],
            seriestype = :path,
            label="",
            c=:red,
            linewidth=3,
            alpha=0.1
        )
    end
    for r in find_true_ranges(cdf[!, "-3σ"]; min_length=1000)
        range_mean = Statistics.mean(cdf[r[1]:r[2], :depth])
        StatsPlots.vline!(p,
            [r[1], r[2]],
            # [range_mean, range_mean],
            seriestype = :path,
            label="",
            c=:red,
            linewidth=3,
            alpha=1/3
        )
    end

    color_vec = StatsPlots.cgrad(ColorSchemes.rainbow, 7, categorical = true)
    
    
    
    StatsPlots.plot!(p, equally_spaced_samples(cdf[!, "index"], 10_000), equally_spaced_samples(cdf[!, "depth"], 10_000), label="coverage", c=:black)
    # StatsPlots.plot!(p, equally_spaced_samples(cdf[!, "index"], 10_000), equally_spaced_samples(rolling_centered_avg(cdf[!, "depth"], window_size=1001), 10_000), label="101bp sliding window mean", c=:gray)
    mean_coverage = first(unique(cdf[!, "mean_coverage"]))
    stddev_coverage = first(unique(cdf[!, "std_coverage"]))
    StatsPlots.hline!(p, [mean_coverage + 3 * stddev_coverage], label="+3σ", c=color_vec[7])
    StatsPlots.hline!(p, [mean_coverage + 2 * stddev_coverage], label="+2σ", c=color_vec[6])
    StatsPlots.hline!(p, [mean_coverage + 1 * stddev_coverage], label="+σ", c=color_vec[5])
    StatsPlots.hline!(p, unique(cdf[!, "mean_coverage"]), label="mean_coverage", c=color_vec[4])
    StatsPlots.hline!(p, [mean_coverage + -1 * stddev_coverage], label="-σ", c=color_vec[3])
    StatsPlots.hline!(p, [mean_coverage + -2 * stddev_coverage], label="-2σ", c=color_vec[2])
    StatsPlots.hline!(p, [mean_coverage + -3 * stddev_coverage], label="-3σ", c=color_vec[1])
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a multi-panel visualization of genome coverage across chromosomes.

# Arguments
- `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# Returns
- `Plots.Figure`: A composite figure with coverage plots for each chromosome

# Details
Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
distribution across genomic positions for that chromosome.
"""
function visualize_genome_coverage(coverage_table)
    num_plots = length(unique(coverage_table[!, "chromosome"]))
    meta_figure = StatsPlots.plot(
        [chromosome_coverage_table_to_plot(cdf) for cdf in DataFrames.groupby(coverage_table, "chromosome")]...,
        layout = (num_plots, 1),
        size = (800, 600 * num_plots)) # Adjust size as needed
    return meta_figure
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Merge two colors by calculating their minimal color difference vector.

# Arguments
- `c1::Color`: First color input
- `c2::Color`: Second color input 

# Returns
- If colors are equal, returns the input color
- Otherwise returns the color difference vector (c1-c2 or c2-c1) with minimal RGB sum

# Details
Calculates two difference vectors:
- mix_a = c1 - c2 
- mix_b = c2 - c1
Returns the difference vector with the smallest sum of RGB components.
"""
function merge_colors(c1, c2)
    if c1 == c2
        return c1
    else
        mix_a = c1 - c2
        mix_b = c2 - c1
        mix_a_sum = mix_a.r + mix_a.g + mix_a.b
        mix_b_sum = mix_b.r + mix_b.g + mix_b.b
        min_value, min_index = findmin([mix_a_sum, mix_b_sum])
        mixed_color = [mix_a, mix_b][min_index]
        # return Colors.color_names["black"]
        return mixed_color
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Visualizes cluster assessment metrics and saves the resulting plots.

# Arguments
- `clustering_results`: A named tuple containing:
    * `ks_assessed`: Vector of k values tested
    * `within_cluster_sum_of_squares`: Vector of WCSS scores
    * `silhouette_scores`: Vector of silhouette scores
    * `optimal_number_of_clusters`: Integer indicating optimal k

# Details
Creates two plots:
1. Within-cluster sum of squares (WCSS) vs number of clusters
2. Silhouette scores vs number of clusters

Both plots include a vertical line indicating the optimal number of clusters.

# Outputs
Saves two SVG files in the project directory:
- `wcss.svg`: WCSS plot
- `silhouette.svg`: Silhouette scores plot
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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function my_plot(graph::KmerGraph)
#     graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))
#     filename = "/assets/images/$(graph_hash).svg"
#     p = plot_graph(graph)
#     Plots.savefig(p, dirname(pwd()) * filename)
#     display(p)
#     display("text/markdown", "![]($filename)")
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a visualization of a kmer graph where nodes represent kmers and their sizes reflect counts.

# Arguments
- `graph`: A MetaGraph where vertices have `:kmer` and `:count` properties

# Returns
- A Plots.jl plot object showing the graph visualization

# Details
- Node sizes are scaled based on kmer counts
- Plot dimensions scale logarithmically with number of vertices
- Each node is labeled with its kmer sequence
"""
function plot_graph(graph)
    
#     kmer_counts = MetaGraphs.get_prop(graph, :kmer_counts)
    kmers = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
    counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]
    scale = 150
    
    n = Graphs.nv(graph)
    p = GraphRecipes.graphplot(
        graph,
#         markersize = 1/log2(n),
        markersize = 1/2^2,
        size = (2 * scale * log(n), scale * log(n)),
        node_weights = counts,
        names = kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Plots a histogram of kmer counts against # of kmers with those counts

Returns the plot object for adding additional layers and saving

Creates a scatter plot visualizing the k-mer frequency spectrum - the relationship
between k-mer frequencies and how many k-mers occur at each frequency.

# Arguments
- `counts::AbstractVector{<:Integer}`: Vector of k-mer counts/frequencies
- `log_scale::Union{Function,Nothing} = log2`: Function to apply logarithmic scaling to both axes.
  Set to `nothing` to use linear scaling.
- `kwargs...`: Additional keyword arguments passed to `StatsPlots.plot()`

# Returns
- `Plots.Plot`: A scatter plot object that can be further modified or saved

# Details
The x-axis shows k-mer frequencies (how many times each k-mer appears),
while the y-axis shows how many distinct k-mers appear at each frequency.
Both axes are log-scaled by default using log2.
"""
function plot_kmer_frequency_spectra(counts; log_scale = log2, kwargs...)
    kmer_counts_hist = StatsBase.countmap(c for c in counts)
    xs = collect(keys(kmer_counts_hist))
    ys = collect(values(kmer_counts_hist))
    if isa(log_scale, Function)
        xs = log_scale.(xs)
        ys = log_scale.(ys)
    end
    
    p = StatsPlots.plot(
        xs,
        ys,
        xlims = (0, maximum(xs) + max(1, ceil(0.1 * maximum(xs)))),
        ylims = (0, maximum(ys) + max(1, ceil(0.1 * maximum(ys)))),
        seriestype = :scatter,
        legend = false,
        xlabel = isa(log_scale, Function) ? "$(log_scale)(observed frequency)" : "observed frequency",
        ylabel = isa(log_scale, Function) ? "$(log_scale)(# of kmers)" : "observed frequency",
        ;kwargs...
    )
    return p
end