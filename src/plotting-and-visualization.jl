"""
    plot_batch_qc_distributions(df::DataFrames.DataFrame; title="Batch QC Summary")

Visualize distributions of QC metrics across multiple samples.
Panel 1: Sorted Yield (Knee Plot).
Panel 2: Read Count Distribution.
Panel 3: Read Length Distributions.
Panel 4: Base Quality Distribution.
"""
function plot_batch_qc_distributions(df::DataFrames.DataFrame; title="Cohort QC Overview")
    # Filter for clean data only
    clean_df = filter(row -> row.stage == "after_filtering", df)
    
    n_samples = DataFrames.nrow(clean_df)
    if n_samples == 0
        @warn "No 'after_filtering' data found."
        return CairoMakie.Figure()
    end

    fig = CairoMakie.Figure(size = (1200, 900), fontsize=18)
    
    # UPDATED: Added halign=:center to center the title
    CairoMakie.Label(fig[0, :], "$title (N=$n_samples)", fontsize=24, font=:bold, halign=:center)

    # --- Panel 1: Yield Rank (Knee Plot) ---
    ax1 = CairoMakie.Axis(fig[1, 1], 
        title="Yield per Sample (Ranked)", 
        xlabel="Sample Rank", 
        ylabel="Yield (Gb)")
    
    sorted_yields = sort(clean_df.yield_gb)
    CairoMakie.barplot!(ax1, 1:length(sorted_yields), sorted_yields, 
        color=sorted_yields, colormap=:plasma, strokewidth=0.5)

    # --- Panel 2: Read Count Distribution ---
    ax2 = CairoMakie.Axis(fig[1, 2], 
        title="Read Count Distribution", 
        ylabel="Number of Reads",
        xticksvisible=false, xticklabelsvisible=false)
    
    CairoMakie.violin!(ax2, fill(1, n_samples), clean_df.total_reads, color=(:teal, 0.5), show_median=true)
    CairoMakie.boxplot!(ax2, fill(1, n_samples), clean_df.total_reads, color=:black, width=0.1)

    # --- Panel 3: Read Length Distributions ---
    ax3 = CairoMakie.Axis(fig[2, 1], 
        title="Read Length Distributions", 
        ylabel="Length (bp)", 
        xticks=([1, 2], ["Mean", "N50"]))
    
    CairoMakie.violin!(ax3, fill(1, n_samples), clean_df.mean_length, color=(:dodgerblue, 0.5), show_median=true)
    CairoMakie.boxplot!(ax3, fill(1, n_samples), clean_df.mean_length, color=:black, width=0.1)

    if "n50" in names(clean_df) && any(clean_df.n50 .> 0)
        CairoMakie.violin!(ax3, fill(2, n_samples), clean_df.n50, color=(:orange, 0.5), show_median=true)
        CairoMakie.boxplot!(ax3, fill(2, n_samples), clean_df.n50, color=:black, width=0.1)
    end

    # --- Panel 4: Quality Distribution ---
    ax4 = CairoMakie.Axis(fig[2, 2], 
        title="Base Quality Distribution", 
        ylabel="Percentage (%)", 
        limits=(nothing, nothing, 0, 100),
        xticks=([1, 2], ["Q20%", "Q30%"]))
    
    CairoMakie.violin!(ax4, fill(1, n_samples), clean_df.q20_percent, color=(:lightgreen, 0.5), show_median=true)
    CairoMakie.violin!(ax4, fill(2, n_samples), clean_df.q30_percent, color=(:forestgreen, 0.5), show_median=true)
    
    CairoMakie.boxplot!(ax4, fill(1, n_samples), clean_df.q20_percent, color=:black, width=0.1)
    CairoMakie.boxplot!(ax4, fill(2, n_samples), clean_df.q30_percent, color=:black, width=0.1)
    
    CairoMakie.hlines!(ax4, [90], color=:red, linestyle=:dash)

    return fig
end

"""
    visualize_fastplong_single(data; title=nothing)

Create a detailed QC report for a single sample.
Includes: Yield, Lengths, Quality (no GC), and Filtering Stats.
"""
function visualize_fastplong_single(data; title=nothing)
    if isnothing(data)
        return CairoMakie.Figure()
    end
    
    df = data.summary
    sample_name = isnothing(title) ? data.sample_id : title
    
    fig = CairoMakie.Figure(size = (1200, 800), fontsize=18)
    
    # UPDATED: Added halign=:center to center the title
    CairoMakie.Label(fig[0, :], "QC Report: $sample_name", fontsize=24, font=:bold, halign=:center)

    colors = [:grey80, :dodgerblue]
    stages = ["Raw", "Filtered"]

    # --- Panel 1: Reads & Yield ---
    ax1 = CairoMakie.Axis(fig[1, 1], title="Sequencing Yield", ylabel="Reads", xticks=(1:2, stages))
    CairoMakie.barplot!(ax1, 1:2, df.total_reads, color=colors, strokewidth=1)
    
    ax1_r = CairoMakie.Axis(fig[1, 1], yaxisposition=:right, ylabel="Yield (Gb)")
    CairoMakie.hidespines!(ax1_r); CairoMakie.hidexdecorations!(ax1_r)
    CairoMakie.scatter!(ax1_r, 1:2, df.yield_gb, color=:orange, markersize=15, label="Gb")
    CairoMakie.axislegend(ax1_r, position=:rt)

    # --- Panel 2: Read Lengths ---
    has_n50 = any(df.n50 .> 0)
    
    ax2 = CairoMakie.Axis(fig[1, 2], title="Read Lengths", ylabel="Base Pairs (bp)", xticks=(1:2, stages))
    
    if has_n50
        CairoMakie.barplot!(ax2, [0.85, 1.85], df.mean_length, width=0.3, color=colors, label="Mean")
        CairoMakie.barplot!(ax2, [1.15, 2.15], df.n50, width=0.3, color=colors, gap=0, hatch=:x, label="N50")
        CairoMakie.axislegend(ax2, position=:lt)
    else
        CairoMakie.barplot!(ax2, 1:2, df.mean_length, width=0.5, color=colors, label="Mean")
        CairoMakie.text!(ax2, 1.5, maximum(df.mean_length)*0.5, text="N50 not available", align=(:center, :center))
    end

    # --- Panel 3: Base Quality ---
    ax3 = CairoMakie.Axis(fig[2, 1], title="Base Quality", ylabel="Percentage (%)", limits=(nothing, nothing, 0, 100), xticks=(1:2, stages))
    CairoMakie.barplot!(ax3, [0.85, 1.85], df.q20_percent, width=0.3, color=colors, label="Q20%")
    CairoMakie.barplot!(ax3, [1.15, 2.15], df.q30_percent, width=0.3, color=colors, hatch=:/, label="Q30%")
    CairoMakie.hlines!(ax3, [90], color=:gray, linestyle=:dash)
    CairoMakie.axislegend(ax3, position=:rb)

    # --- Panel 4: Filtering Stats ---
    ax4 = CairoMakie.Axis(fig[2, 2], title="Filtering Reasons", ylabel="Reads Dropped", xticklabelrotation=π/4)
    
    if !isempty(data.filtering_stats)
        drop_reasons = filter(p -> p.first != "passed_filter_reads" && p.second > 0, data.filtering_stats)
        
        if !isempty(drop_reasons)
            reasons = collect(keys(drop_reasons))
            counts = collect(values(drop_reasons))
            labels = replace.(reasons, "_reads" => "", "_" => " ")
            
            CairoMakie.barplot!(ax4, 1:length(counts), counts, color=:firebrick)
            ax4.xticks = (1:length(counts), labels)
        else
            CairoMakie.text!(ax4, 0.5, 0.5, text="No reads dropped!", align=(:center, :center))
        end
    else
        CairoMakie.text!(ax4, 0.5, 0.5, text="No filtering data", align=(:center, :center))
    end

    return fig
end

"""
Generate a taxonomic rank consensus visualization from sampled data.

Returns a CairoMakie Figure object showing consensus agreement across taxonomic ranks.
Samples data before aggregation for improved performance on large datasets.
"""
function create_rank_consensus_plot(
    df::DataFrames.DataFrame;
    taxonomic_ranks::Vector{String} = [
        "domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"
    ],
    figure_size::Tuple{Int,Int} = (1200, 700),
    max_samples::Int = 10_000,
    perfect_threshold::Float64 = 0.999,
    title::String = "Mapping Alignment Consensus Agreement Across Taxonomic Ranks",
    sample_id::Union{String, Nothing} = nothing,
    title_fontsize::Int = 16,
    label_fontsize::Int = 16,
    tick_fontsize::Int = 16,
    legend_fontsize::Int = 16,
    include_n_observations::Bool = true
)

    function monotonicize_agreement(ys::Vector{<:Union{Missing, Float64}})
        n = length(ys)
        # Work from the end backwards
        for i in (n-1):-1:1
            if !ismissing(ys[i+1])
                if ismissing(ys[i])
                    ys[i] = ys[i+1]
                else
                    ys[i] = max(ys[i], ys[i+1])
                end
            end
        end
        return ys
    end

    function has_any_non_missing(ys::Vector{<:Union{Missing, Float64}})
        return any(!ismissing(y) for y in ys)
    end

    # Add diagnostics at the beginning
    # println("Input DataFrame: $(DataFrames.nrow(df)) total rows")
    
    unique_templates = unique(df.template)
    n_total_observations = length(unique_templates)
    println("Unique templates in input: $(n_total_observations)")

    unique_read_mapping_results = unique(df[!, ["template", "ismapped"]])
    total_reads_mapped = count(unique_read_mapping_results.ismapped)
    println("total_reads_mapped = $(total_reads_mapped)")
    println("total_reads = $(n_total_observations)")
    percent_reads_mapped = round(total_reads_mapped / n_total_observations * 100, digits=3)
    println("percentage_of_reads_mapped = $(percent_reads_mapped)")

    mapped_df = df[df.ismapped, :]
    unique_mapped_templates = unique(mapped_df.template)

    # println(count(df.ismapped) / DataFrames.nrow(df)
    
    # After sampling
    n_samples_to_take = min(length(unique_mapped_templates), max_samples)
    sampled_templates = Set(StatsBase.sample(unique_mapped_templates, n_samples_to_take, replace=false))
    println("Templates sampled: $(length(sampled_templates))")
    
    # sampled_df = DataFrames.filter(row -> row.template in sampled_templates, df)
    # println("Rows after template sampling: $(DataFrames.nrow(sampled_df))")
    
    # Filter dataframe to only include sampled templates
    sampled_df = DataFrames.filter(row -> row.template in sampled_templates, mapped_df)
    
    # NOW aggregate only the sampled data
    aggregated = aggregate_by_rank_nonmissing(sampled_df, taxonomic_ranks)
    
    # APPLY THE FIX: Get the maximum alignment proportion for each template-rank combination
    max_aggregated = DataFrames.combine(
        DataFrames.first,
        DataFrames.groupby(
            DataFrames.sort(
                aggregated,
                [:template, :rank, DataFrames.order(:total_relative_alignment_proportion, rev=true)]
            ),
            [:template, :rank]
        )
    )
    println("Rows after taking max per template-rank: $(DataFrames.nrow(max_aggregated))")
    
    # Group by template (this will be much smaller now)
    template_grouped_aggregated = DataFrames.groupby(max_aggregated, "template")
    
    numbered_ranks = ["$(i)_$(rank)" for (i, rank) in enumerate(taxonomic_ranks)]
    all_ys = []
    n_unmapped = 0
    
    # Process all groups (since we already sampled at the template level)
    for this_group in template_grouped_aggregated
        this_sorted_group = DataFrames.sort(this_group, :rank)
        # Now this dictionary creation is safe - each rank appears only once per template
        template_values_dict = Dict(row.rank => row.total_relative_alignment_proportion for row in DataFrames.eachrow(this_sorted_group))
        ys = [get(template_values_dict, this_rank, missing) for this_rank in numbered_ranks]
        # Make ys monotonically non-increasing
        ys = monotonicize_agreement(ys)
        
        # Only include reads that have at least one non-missing value
        if has_any_non_missing(ys)
            push!(all_ys, ys)
        else
            n_unmapped += 1
        end
    end

    # println("Final reads included: $(length(all_ys))")
    # println("Reads with no taxonomic data: $(n_unmapped)")

    # Compute aggregate consensus proportions
    n_ranks = length(numbered_ranks)
    perfect_consensus_props = Float64[]
    for rank_idx in 1:n_ranks
        # Get all non-missing consensus values for this rank
        rank_ys = [ys[rank_idx] for ys in all_ys if !ismissing(ys[rank_idx])]
        n_perfect = count(y -> y >= perfect_threshold, rank_ys)
        n_total = length(rank_ys)
        prop = n_total > 0 ? n_perfect / n_total : 0.0
        push!(perfect_consensus_props, prop)
    end

    # Build the title
    full_title = title
    if include_n_observations
        full_title = "$(full_title) (n=$(length(all_ys))/$(total_reads_mapped) mapped) [mapped = $(percent_reads_mapped)%]"
    end
    if sample_id !== nothing
        full_title = "$(full_title)\n$(sample_id)"
    end

    # Create the figure
    fig = CairoMakie.Figure(size=figure_size)

    # Title across the top
    fig[1, 1:2] = CairoMakie.Label(fig, full_title;
        fontsize=title_fontsize, font="bold", halign=:center)

    # Main plot
    ax = CairoMakie.Axis(
        fig[2, 1];
        xlabel="Taxonomic Rank (High to Low)",
        ylabel="Cumulative Relative Alignment Proportion",
        xticks=(1:length(taxonomic_ranks), taxonomic_ranks),
        yticks=0:0.1:1.0,
        limits=(0.2, length(taxonomic_ranks)+0.8, -0.01, 1.01),
        xlabelsize=label_fontsize,
        ylabelsize=label_fontsize,
        xticklabelsize=tick_fontsize,
        yticklabelsize=tick_fontsize,
        xticklabelrotation=π/4,
        xgridcolor=(:gray, 0.3),
        ygridcolor=(:gray, 0.3),
        xgridwidth=1.2,
        ygridwidth=1.2
    )

    # Barplot for perfect consensus
    CairoMakie.barplot!(
        ax,
        1:n_ranks,
        perfect_consensus_props;
        color=(:gray, 0.4),
        width=0.7,
        strokewidth=0,
        label="Perfect consensus"
    )

    # Individual read traces - plot first one with label for legend
    if !isempty(all_ys)
        CairoMakie.lines!(ax, all_ys[1]; color=:blue, linewidth=0.4, alpha=0.22, label="Individual reads")
        for ys in all_ys[2:end]
            CairoMakie.lines!(ax, ys; color=:blue, linewidth=0.4, alpha=0.22)
        end
    end

    # Only compute and plot summary statistics if we have data
    if !isempty(all_ys)
        median_ys = [Statistics.median([ys[i] for ys in all_ys if !ismissing(ys[i])]) for i in 1:n_ranks]
        CairoMakie.lines!(ax, median_ys; color=:orange, linewidth=3, label="Median trajectory")

        mean_ys = [Statistics.mean([ys[i] for ys in all_ys if !ismissing(ys[i])]) for i in 1:n_ranks]
        CairoMakie.lines!(ax, mean_ys; color=:green, linewidth=3, linestyle=:dash, label="Mean trajectory")
    end

    # Legend that automatically sizes itself
    fig[2, 2] = CairoMakie.Legend(
        fig, ax;
        tellwidth=true,  # Let it tell the layout how much width it needs
        tellheight=false,
        valign=:center,
        framevisible=true,
        labelsize=legend_fontsize
    )

    CairoMakie.resize!(fig.scene, figure_size...)
    
    return fig
end

"""
Generate and save taxonomic rank consensus visualizations from an Arrow file.

Loads data from the specified Arrow file, creates the visualization, and saves it as PNG and SVG files.
Output files are named based on the input file path with "_rank-level-confidences" suffix by default.
"""
function generate_rank_consensus_plots(
    inputfile_path::String;
    png_output_path::Union{String, Nothing} = nothing,
    svg_output_path::Union{String, Nothing} = nothing,
    force=false,
    kwargs...  # Pass through any arguments to create_rank_consensus_plot
)
    fig = nothing
    # Determine output paths
    base_path = splitext(inputfile_path)[1]  # Remove extension
    png_path = png_output_path !== nothing ? png_output_path : "$(base_path)_rank-level-confidences.png"
    svg_path = svg_output_path !== nothing ? svg_output_path : "$(base_path)_rank-level-confidences.svg"

    if (!isfile(png_path) && !isfile(svg_path)) || force
        # Load the data
        if occursin(r"\.arrow$", inputfile_path)
            df = DataFrames.DataFrame(Arrow.Table(inputfile_path))
        elseif occursin(r"\.tsv\.gz$", inputfile_path)
            df = Mycelia.read_tsvgz(inputfile_path)
        else
            error()
        end
        display(size(df))
        
        # Create the figure
        fig = create_rank_consensus_plot(df; sample_id=basename(inputfile_path), kwargs...)
        
        # Save the figures
        CairoMakie.save(png_path, fig)
        CairoMakie.save(svg_path, fig)
    
        println("Saved PNG: $(png_path)")
        println("Saved SVG: $(svg_path)")
    end
    
    return (;fig, png_path, svg_path)
end

function visualize_many_timeseries(
    time_series_data::Vector{Vector{Float64}};
    title::String = "High-Density Time Series Visualization",
    # many-trace panel controls
    many_trace_count::Int = 2000,
    trace_alpha::Float64 = 0.03,
    trace_color = :black,
    trace_linewidth::Float64 = 0.6,
    rng::Random.AbstractRNG = Random.default_rng()
)
    # Create x-axis values (assuming equal length for all series)
    n_points = length(time_series_data[1])
    x_values = 1:n_points

    # Calculate statistics for the data
    all_values = vcat(time_series_data...)
    global_min, global_max = minimum(all_values), maximum(all_values)

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

    # Figure with three rows now (density, bands, many-trace)
    fig = CairoMakie.Figure(size=(1200, 1200), fontsize=18)

    # 1) Density heatmap
    ax1 = CairoMakie.Axis(fig[1, 1],
        title = "Density Heatmap",
        xlabel = "Time",
        ylabel = "Value",
        titlesize = 22)

    x_coords = Float64[]
    y_coords = Float64[]
    for series in time_series_data
        for (point_idx, value) in enumerate(series)
            push!(x_coords, point_idx)
            push!(y_coords, value)
        end
    end

    bins_x = 100
    bins_y = 100
    x_range = (minimum(x_coords), maximum(x_coords))
    y_range = (minimum(y_coords), maximum(y_coords))
    x_edges = range(x_range[1], x_range[2], length = bins_x + 1)
    y_edges = range(y_range[1], y_range[2], length = bins_y + 1)

    histogram = zeros(bins_x, bins_y)
    for i in 1:length(x_coords)
        x = x_coords[i]
        y = y_coords[i]
        x_bin = max(1, min(bins_x, Int(floor((x - x_range[1]) / (x_range[2] - x_range[1]) * bins_x)) + 1))
        y_bin = max(1, min(bins_y, Int(floor((y - y_range[1]) / (y_range[2] - y_range[1]) * bins_y)) + 1))
        histogram[x_bin, y_bin] += 1
    end

    hmap = CairoMakie.heatmap!(ax1,
        x_edges[1:end-1],
        y_edges[1:end-1],
        histogram,
        colormap = :viridis)

    CairoMakie.lines!(ax1, x_values, means, color = :white, linewidth = 3)
    CairoMakie.Colorbar(fig[1, 2], hmap, label = "Density")

    # 2) Sample with statistical bands
    ax2 = CairoMakie.Axis(fig[2, 1],
        title = "Sample with Statistical Bands",
        xlabel = "Time",
        ylabel = "Value",
        titlesize = 22)

    sample_size = min(500, length(time_series_data))
    sample_indices = rand(rng, 1:length(time_series_data), sample_size)

    for idx in sample_indices
        series = time_series_data[idx]
        CairoMakie.lines!(ax2, x_values[1:length(series)], series,
            color = (ColorSchemes.viridis[rand(rng)], 0.05), linewidth = 0.5)
    end

    CairoMakie.band!(ax2, x_values, q10, q90, color = (:blue, 0.2))
    CairoMakie.band!(ax2, x_values, q25, q75, color = (:blue, 0.3))
    CairoMakie.lines!(ax2, x_values, means, color = :red, linewidth = 3)

    CairoMakie.Legend(fig[2, 2],
        [
            CairoMakie.LineElement(color = :red, linewidth = 3),
            CairoMakie.PolyElement(color = (:blue, 0.2)),
            CairoMakie.PolyElement(color = (:blue, 0.3))
        ],
        ["Mean", "10-90% Range", "25-75% Range"]
    )

    # 3) Many sampled line trace plot (thin, semi-transparent traces)
    ax3 = CairoMakie.Axis(fig[3, 1],
        title = "Many Sampled Traces",
        xlabel = "Time",
        ylabel = "Value",
        titlesize = 22)

    # sample without replacement, like your StatsBase.sample snippet
    sampled_series = StatsBase.sample(rng, time_series_data,
        min(many_trace_count, length(time_series_data));
        replace = false)

    for series in sampled_series
        if !isempty(series)
            CairoMakie.lines!(ax3,
                x_values[1:length(series)],
                series,
                color = (trace_color, trace_alpha),
                linewidth = trace_linewidth)
        end
    end

    # Optional: overlay mean for reference
    CairoMakie.lines!(ax3, x_values, means, color = (:black, 0.6), linewidth = 2)

    # Link axes across all three panels
    CairoMakie.linkaxes!(ax1, ax2, ax3)

    # Overall title
    CairoMakie.Label(fig[0, :], text = title, fontsize = 26)

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
    missing_count = 0
    if !isempty(missing_keys)
        missing_count = sum(joint_counts[k] for k in missing_keys)
    end
    
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
    
    # # Set reproducible colors if seed is provided
    # if !isnothing(color_seed)
    #     Random.seed!(color_seed)
    # end
    
    # Generate distinctive colors, reserving specific colors for "Other" and "Missing"
    n_colors_needed = length(final_taxa)
    
    # Generate base colors for all taxa except special categories
    if n_colors_needed > 0
        # colorscheme = Colors.distinguishable_colors(n_colors_needed, [Colors.RGB(1,1,1), Colors.RGB(0,0,0)], dropseed=true)
        colorscheme = Mycelia.n_maximally_distinguishable_colors(n_colors_needed)
        # colorscheme = reverse(colorscheme)  # Reverse to match original behavior
        
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
- `img`: Optional output image path. Defaults to `gfa * ".png"` or `format` extension
- `format`: Optional image format extension (e.g., `png`, `svg`) used when `img` is not provided
- `extra_args`: Additional CLI flags passed directly to Bandage (e.g., layout or label options)
- `force`: Re-render even if the output already exists

# Returns
- Path to the generated image file
"""
function bandage_visualize(; gfa, img=nothing, format::Union{Nothing, String}=nothing, extra_args::AbstractVector{<:AbstractString}=String[], force::Bool=false)
    # run(`$(bandage) image --helpall`)
    bandage = Mycelia.download_bandage()
    target_img = something(img, gfa * "." * (format === nothing ? "png" : format))
    if force || !isfile(target_img)
        cmd_parts = [bandage, "image", gfa, target_img]
        append!(cmd_parts, extra_args)
        run(Cmd(cmd_parts))
    end
    return target_img
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Reduce a FASTG to a simplified GFA using Bandage.

# Arguments
- `fastg`: Path to input FASTG file
- `gfa`: Optional output GFA path (defaults to `fastg * ".gfa"`)

# Returns
- Path to the reduced GFA file
"""
function bandage_reduce(;fastg, gfa=fastg*".gfa")
    bandage = Mycelia.download_bandage()
    if !isfile(gfa)
        run(`$(bandage) reduce $(fastg) $(gfa)`)
    end
    return gfa
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
    color_spectrum = Colors.distinguishable_colors(n+3, [Colors.RGB(1,1,1), Colors.RGB(0,0,0)], dropseed=false)
    filtered_color_spectrum = color_spectrum[[2, 5:length(color_spectrum)...]]
    return filtered_color_spectrum
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
- `mg::MetaGraphsNext.MetaGraph`: Graph containing hierarchical clustering results. Must have
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
        mg::MetaGraphsNext.MetaGraph;
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
- `mg::MetaGraphsNext.MetaGraph`: A meta directed graph containing hierarchical clustering data
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
        mg::MetaGraphsNext.MetaGraph;
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


# -----------------------------------------------------------------------------
# Radial Dendrogram Visualization
# -----------------------------------------------------------------------------

"""
Internal structure to represent a radial dendrogram node.
"""
mutable struct RadialNode
    left::Union{RadialNode, Nothing}
    right::Union{RadialNode, Nothing}
    height::Float64
    angle::Float64
    leaf_index::Int # -1 for internal nodes, >0 for original leaf index
end

"""
    hclust_to_radial_tree(hclust::Clustering.Hclust)

Convert a Clustering.Hclust object into a linked RadialNode tree structure.
"""
function hclust_to_radial_tree(hclust::Clustering.Hclust)
    n = length(hclust.order)
    merges = hclust.merges
    heights = hclust.heights
    
    nodes = Dict{Int, RadialNode}()
    
    # Initialize leaves
    for i in 1:n
        nodes[-i] = RadialNode(nothing, nothing, 0.0, 0.0, i)
    end
    
    # Build tree from merges
    root_node = nothing
    for i in 1:size(merges, 1)
        left_idx = merges[i, 1]
        right_idx = merges[i, 2]
        height = heights[i]
        
        left_node = nodes[left_idx]
        right_node = nodes[right_idx]
        
        parent_node = RadialNode(left_node, right_node, height, 0.0, -1)
        nodes[i] = parent_node
        root_node = parent_node
    end
    
    return root_node
end

"""
    assign_radial_angles!(node, leaf_order_map, angle_map, start_angle, end_angle)

Recursively assign angular positions to nodes based on the number of leaves in their subtrees.
"""
function assign_radial_angles!(
    node::RadialNode,
    leaf_order_map::Dict{Int, Int},
    angle_map::Dict{Int, Float64},
    start_angle::Float64,
    end_angle::Float64
)
    if isnothing(node.left) && isnothing(node.right)
        # It's a leaf
        leaf_pos = leaf_order_map[node.leaf_index]
        n_leaves = length(leaf_order_map)
        
        # Calculate angle based on position in the optimal ordering
        # Using (leaf_pos - 1) / (n_leaves - 1) spreads leaves exactly from start to end
        angle_fraction = n_leaves > 1 ? (leaf_pos - 1) / (n_leaves - 1) : 0.5
        node.angle = start_angle + angle_fraction * (end_angle - start_angle)
        angle_map[node.leaf_index] = node.angle
    else
        # It's an internal node
        left_leaves = count_leaves(node.left)
        right_leaves = count_leaves(node.right)
        total_leaves = left_leaves + right_leaves
        
        # Proportional split of the angular wedge
        mid_angle = start_angle + (left_leaves / total_leaves) * (end_angle - start_angle)
        
        assign_radial_angles!(node.left, leaf_order_map, angle_map, start_angle, mid_angle)
        assign_radial_angles!(node.right, leaf_order_map, angle_map, mid_angle, end_angle)
        
        # Parent angle is the midpoint of children angles
        node.angle = (node.left.angle + node.right.angle) / 2
    end
end

"""
Helper to count leaves in a RadialNode subtree.
"""
function count_leaves(node::RadialNode)
    if isnothing(node.left) && isnothing(node.right)
        return 1
    else
        return count_leaves(node.left) + count_leaves(node.right)
    end
end

"""
    plot_radial_dendrogram_branches!(ax, node, max_height, inner_radius, outer_radius; color, linewidth)

Recursively plot the lines (arcs and radials) of the dendrogram.
"""
function plot_radial_dendrogram_branches!(
    ax,
    node::RadialNode,
    max_height::Float64,
    inner_radius::Float64,
    outer_radius::Float64;
    line_color=(:black, 0.4),
    linewidth=1.0
)
    if isnothing(node.left) || isnothing(node.right)
        return
    end
    
    # Calculate radius r based on height (inverted: root is at inner_radius, leaves at outer_radius is standard,
    # but dendrograms usually have root at center (0 height) or outside. 
    # Standard dendrogram: leaves are at height 0, root at max_height.
    # Here we map height 0 (leaves) to outer_radius, and max_height (root) to inner_radius.
    
    # Normalize height: 0.0 (leaves) -> 1.0 (tips), max_height -> 0.0 (root)
    # Actually, hclust heights go from small (leaves merged) to large (root).
    # We want root at center (inner) and leaves at periphery (outer).
    
    function get_r(h)
        # Linear interpolation
        # h=0 (leaves) => outer_radius
        # h=max_height => inner_radius
        return outer_radius - (h / max_height) * (outer_radius - inner_radius)
    end

    parent_r = get_r(node.height)
    left_r = get_r(node.left.height)
    right_r = get_r(node.right.height)
    
    # Polar to Cartesian
    px, py = parent_r .* (cos(node.angle), sin(node.angle))
    lx, ly = left_r .* (cos(node.left.angle), sin(node.left.angle))
    rx, ry = right_r .* (cos(node.right.angle), sin(node.right.angle))
    
    # Draw arcs between children angles at parent radius
    # We need to draw the crossbar at the parent's radius spanning from left child angle to right child angle
    
    theta_left = node.left.angle
    theta_right = node.right.angle
    
    # Ensure we take the shortest path around circle? 
    # For dendrograms, we usually strictly follow the wedge.
    # Just ensure ordering for range
    if theta_left > theta_right
        theta_left, theta_right = theta_right, theta_left
    end
    
    n_arc = max(2, Int(floor(50 * abs(theta_right - theta_left))))
    arc_angles = range(theta_left, theta_right, length=n_arc)
    
    arc_x = [parent_r * cos(a) for a in arc_angles]
    arc_y = [parent_r * sin(a) for a in arc_angles]
    
    # Plot the "crossbar" arc
    CairoMakie.lines!(ax, arc_x, arc_y, color=line_color, linewidth=linewidth)
    
    # Plot the "drops" to children
    # Line from left child's radius to parent's radius at left child's angle
    # Wait, standard dendrograms usually:
    # 1. Draw stems from children up to parent height.
    # 2. Draw crossbar at parent height.
    
    # Stem for left child
    l_stem_start_x = left_r * cos(node.left.angle)
    l_stem_start_y = left_r * sin(node.left.angle)
    l_stem_end_x = parent_r * cos(node.left.angle)
    l_stem_end_y = parent_r * sin(node.left.angle)
    CairoMakie.lines!(ax, [l_stem_start_x, l_stem_end_x], [l_stem_start_y, l_stem_end_y], color=line_color, linewidth=linewidth)

    # Stem for right child
    r_stem_start_x = right_r * cos(node.right.angle)
    r_stem_start_y = right_r * sin(node.right.angle)
    r_stem_end_x = parent_r * cos(node.right.angle)
    r_stem_end_y = parent_r * sin(node.right.angle)
    CairoMakie.lines!(ax, [r_stem_start_x, r_stem_end_x], [r_stem_start_y, r_stem_end_y], color=line_color, linewidth=linewidth)
    
    # Recurse
    plot_radial_dendrogram_branches!(ax, node.left, max_height, inner_radius, outer_radius; line_color=line_color, linewidth=linewidth)
    plot_radial_dendrogram_branches!(ax, node.right, max_height, inner_radius, outer_radius; line_color=line_color, linewidth=linewidth)
end

"""
    plot_radial_dendrogram(hclust; rings, title, ...)

Generate a generalized radial dendrogram with optional concentric annotation rings.

# Arguments
- `hclust`: `Clustering.Hclust` object.
- `rings`: Vector of NamedTuples or Dicts defining annotation layers. 
   Each element should look like `(indices=[...], label="Name", color=:red, marker=:circle, size=10)`.
   Only `indices` is strictly required; defaults will be supplied for others.

# Keywords
- `title`: Plot title.
- `figure_size`: Tuple (width, height).
- `inner_radius`: Radius where the tree root sits (can be 0 for center).
- `outer_radius`: Radius where tree leaves end.
- `ring_start_radius`: Radius where the first annotation ring starts (default: `outer_radius + 0.05`).
- `ring_spacing`: Gap between annotation rings.
- `line_color`: Color of dendrogram branches.
- `line_width`: Width of dendrogram branches.

# Returns
- `CairoMakie.Figure`
"""
function plot_radial_dendrogram(
    hclust::Clustering.Hclust;
    rings::Vector = [], # Vector of NamedTuples/Dicts
    title::String = "",
    figure_size::Tuple{Int,Int} = (1000, 1000),
    inner_radius::Float64 = 0.0, # Tree root
    outer_radius::Float64 = 0.8, # Tree tips
    ring_start_radius::Union{Float64, Nothing} = nothing,
    ring_spacing::Float64 = 0.05,
    line_color = (:black, 0.5),
    line_width = 1.0,
    legend_title = "Legend"
)
    # 1. Build Radial Tree
    root = hclust_to_radial_tree(hclust)
    n_leaves = length(hclust.order)
    
    # 2. Assign Angles (Optimized)
    # Map leaf index to its position in the sort order (1..N)
    leaf_order_map = Dict{Int,Int}(idx => pos for (pos, idx) in enumerate(hclust.order))
    angle_map = Dict{Int, Float64}()
    
    start_angle = -π/2
    end_angle = start_angle + 2π
    
    assign_radial_angles!(root, leaf_order_map, angle_map, start_angle, end_angle)
    max_height = maximum(hclust.heights)

    # 3. Setup Figure
    fig = CairoMakie.Figure(size=figure_size)
    ax = CairoMakie.Axis(fig[1, 1], aspect=CairoMakie.DataAspect(), title=title)
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)

    # 4. Plot Tree
    plot_radial_dendrogram_branches!(
        ax, root, max_height, inner_radius, outer_radius; 
        line_color=line_color, linewidth=line_width
    )

    # 5. Plot Rings
    current_radius = isnothing(ring_start_radius) ? outer_radius + 0.02 : ring_start_radius
    
    legend_entries = []
    legend_labels = []
    
    # Default palettes
    default_colors = Mycelia.n_maximally_distinguishable_colors(max(length(rings), 3))
    
    for (i, ring_data) in enumerate(rings)
        # normalize inputs to ensure fields exist
        indices = get(ring_data, :indices, Int[])
        label = get(ring_data, :label, "Group $i")
        color = get(ring_data, :color, default_colors[mod1(i, length(default_colors))])
        marker = get(ring_data, :marker, :circle)
        msize = get(ring_data, :size, 10.0)
        
        xs = Float64[]
        ys = Float64[]
        
        indices_set = Set(indices)
        
        for leaf_idx in hclust.order
            if leaf_idx in indices_set
                ang = angle_map[leaf_idx]
                push!(xs, current_radius * cos(ang))
                push!(ys, current_radius * sin(ang))
            end
        end
        
        if !isempty(xs)
            CairoMakie.scatter!(ax, xs, ys; color=color, marker=marker, markersize=msize, strokewidth=0.5, strokecolor=:black)
            
            # Add to legend
            push!(legend_entries, CairoMakie.MarkerElement(color=color, marker=marker, markersize=15, strokecolor=:black, strokewidth=0.5))
            push!(legend_labels, label)
        end
        
        current_radius += ring_spacing
    end

    # 6. Add Legend
    if !isempty(legend_entries)
        CairoMakie.Legend(
            fig[1, 2],
            legend_entries,
            legend_labels,
            legend_title,
            framevisible=true
        )
    end
    
    # Adjust limits to fit everything
    limit_radius = current_radius + 0.05
    CairoMakie.xlims!(ax, -limit_radius, limit_radius)
    CairoMakie.ylims!(ax, -limit_radius, limit_radius)
    
    return fig
end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a visualization of chromosome coverage data with statistical thresholds.

# # Arguments
# - `cdf::DataFrame`: Coverage data frame containing columns:
#   - `index`: Chromosome position indices
#   - `depth`: Coverage depth values
#   - `chromosome`: Chromosome identifier
#   - `mean_coverage`: Mean coverage value
#   - `std_coverage`: Standard deviation of coverage
#   - `3σ`: Boolean vector indicating +3 sigma regions
#   - `-3σ`: Boolean vector indicating -3 sigma regions

# # Returns
# - A StatsPlots plot object showing:
#   - Raw coverage data (black line)
#   - Mean coverage and ±1,2,3σ thresholds (rainbow colors)
#   - Highlighted regions exceeding ±3σ thresholds (red vertical lines)
# """
# function chromosome_coverage_table_to_plot(cdf)

#     # mean_coverage = first(unique(cdf[!, "mean_coverage"]))
#     # stddev_coverage = first(unique(cdf[!, "std_coverage"]))

#     mean_coverage = Statistics.mean(cdf[!, "depth"])
#     median_coverage = Statistics.median(cdf[!, "depth"])
#     stddev_coverage = Statistics.std(cdf[!, "depth"])

#     p = StatsPlots.plot(
#         xlims = extrema(cdf[!, "index"]),
#         ylims=(0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#         title = cdf[1, "chromosome"],
#         xlabel = "chromosome index",
#         ylabel = "depth"
#     )

#     # cdf[!, "3σ"] = 
#     # cdf[!, "-3σ"] = 

#     # mean_coverage + 3 * stddev_coverage    

#     is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#     for r in find_true_ranges(is_3sigma_above; min_length=1000)
#         range_mean = Statistics.mean(cdf[r[1]:r[2], :depth])
#         StatsPlots.vline!(p,
#             [r[1], r[2]],
#             # [range_mean, range_mean],
#             seriestype = :path,
#             label="",
#             c=:red,
#             linewidth=3,
#             alpha=0.1
#         )
#     end
#     is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#     for r in find_true_ranges(is_3sigma_below; min_length=1000)
#         range_mean = Statistics.mean(cdf[r[1]:r[2], :depth])
#         StatsPlots.vline!(p,
#             [r[1], r[2]],
#             # [range_mean, range_mean],
#             seriestype = :path,
#             label="",
#             c=:red,
#             linewidth=3,
#             alpha=1/3
#         )
#     end

#     color_vec = StatsPlots.cgrad(ColorSchemes.rainbow, 8, categorical = true)
#     StatsPlots.plot!(p, equally_spaced_samples(cdf[!, "index"], 10_000), equally_spaced_samples(cdf[!, "depth"], 10_000), label="coverage", c=:black)
#     # StatsPlots.plot!(p, equally_spaced_samples(cdf[!, "index"], 10_000), equally_spaced_samples(rolling_centered_avg(cdf[!, "depth"], window_size=1001), 10_000), label="101bp sliding window mean", c=:gray)

#     StatsPlots.hline!(p, [mean_coverage + 3 * stddev_coverage], label="+3σ", c=color_vec[8])
#     StatsPlots.hline!(p, [mean_coverage + 2 * stddev_coverage], label="+2σ", c=color_vec[7])
#     StatsPlots.hline!(p, [mean_coverage + 1 * stddev_coverage], label="+σ", c=color_vec[6])
#     StatsPlots.hline!(p, [mean_coverage], label="mean_coverage", c=color_vec[5])
#     StatsPlots.hline!(p, [median_coverage], label="median", c=color_vec[4])
#     StatsPlots.hline!(p, [mean_coverage + -1 * stddev_coverage], label="-σ", c=color_vec[3])
#     StatsPlots.hline!(p, [mean_coverage + -2 * stddev_coverage], label="-2σ", c=color_vec[2])
#     StatsPlots.hline!(p, [mean_coverage + -3 * stddev_coverage], label="-3σ", c=color_vec[1])
#     return p
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a visualization of chromosome coverage data with statistical thresholds.

# # Arguments
# - `cdf::DataFrame`: Coverage data frame containing columns:
#   - `index`: Chromosome position indices
#   - `depth`: Coverage depth values
#   - `chromosome`: Chromosome identifier
#   - `mean_coverage`: Mean coverage value
#   - `std_coverage`: Standard deviation of coverage
#   - `3σ`: Boolean vector indicating +3 sigma regions
#   - `-3σ`: Boolean vector indicating -3 sigma regions

# # Returns
# - A CairoMakie figure showing:
#   - Raw coverage data (black line)
#   - Mean coverage and ±1,2,3σ thresholds (rainbow colors)
#   - Highlighted regions exceeding ±3σ thresholds (red vertical lines)
# """
# function chromosome_coverage_table_to_plot(cdf)
#     mean_coverage = Statistics.mean(cdf[!, "depth"])
#     median_coverage = Statistics.median(cdf[!, "depth"])
#     stddev_coverage = Statistics.std(cdf[!, "depth"])

#     # Create figure with proper spacing for labels and legend
#     fig = CairoMakie.Figure(
#         resolution = (1000, 600),
#         fontsize = 14,
#         figure_padding = (20, 80, 20, 20)  # left, right, bottom, top
#     )
    
#     # Create axis with proper margins
#     ax = CairoMakie.Axis(
#         fig[1, 1],
#         limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#         title = cdf[1, "chromosome"],
#         xlabel = "chromosome index",
#         ylabel = "depth",
#         titlesize = 16,
#         xlabelpadding = 10,
#         ylabelpadding = 10,
#         xgridvisible = true,
#         ygridvisible = true,
#         xgridcolor = (:gray, 0.3),
#         ygridcolor = (:gray, 0.3),
#         topspinevisible = false,
#         rightspinevisible = false
#     )

#     # Add highlighted regions for 3σ outliers
#     is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#     for r in find_true_ranges(is_3sigma_above; min_length=1000)
#         CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.1))
#     end
    
#     is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#     for r in find_true_ranges(is_3sigma_below; min_length=1000)
#         CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.33))
#     end

#     # Create rainbow color scheme
#     color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
    
#     # Plot main coverage line
#     CairoMakie.lines!(ax, 
#         equally_spaced_samples(cdf[!, "index"], 10_000), 
#         equally_spaced_samples(cdf[!, "depth"], 10_000), 
#         color = :black,
#         linewidth = 1.5,
#         label = "coverage"
#     )

#     # Add horizontal threshold lines
#     threshold_lines = [
#         (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
#         (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
#         (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
#         (mean_coverage, "mean_coverage", color_vec[5]),
#         (median_coverage, "median", color_vec[4]),
#         (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
#         (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
#         (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
#     ]
    
#     for (y_val, label, color) in threshold_lines
#         CairoMakie.hlines!(ax, y_val, 
#             color = color, 
#             linewidth = 2,
#             linestyle = :dash,
#             label = label
#         )
#     end

#     # Create external legend
#     CairoMakie.Legend(
#         fig[1, 2], 
#         ax, 
#         framevisible = true,
#         bgcolor = (:white, 0.9),
#         labelsize = 12,
#         rowgap = 5,
#         margin = (10, 10, 10, 10)
#     )

#     # Adjust layout to prevent label cutoff
#     CairoMakie.resize_to_layout!(fig)
    
#     return fig
# end

function visualize_genome_coverage(path_to_coverage_table::AbstractString; kwargs...)
    return visualize_genome_coverage(
        CSV.read(path_to_coverage_table, DataFrames.DataFrame, delim='\t', header=["chromosome", "index", "depth"]);
        kwargs...
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a multi-panel visualization of genome coverage across chromosomes.

# Arguments
- `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# Keyword Arguments
- `entity`: Optional string to specify the entity being analyzed. If provided, chromosome titles become "\$(entity): \$(chromosome)", otherwise just "\$(chromosome)"

# Returns
- `CairoMakie.Figure`: A composite figure with coverage plots for each chromosome

# Details
Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
distribution across genomic positions for that chromosome.
"""
function visualize_genome_coverage(coverage_table::DataFrames.AbstractDataFrame; entity::Union{String, Nothing} = nothing)
    chromosome_groups = DataFrames.groupby(coverage_table, "chromosome")
    num_plots = length(chromosome_groups)
    
    # Create main figure with appropriate size and spacing
    fig = CairoMakie.Figure(
        size = (1400, 400 * num_plots),
        fontsize = 12,
        figure_padding = (20, 40, 20, 20)  # left, right, bottom, top
    )
    
    for (i, cdf) in enumerate(chromosome_groups)
        mean_coverage = Statistics.mean(cdf[!, "depth"])
        median_coverage = Statistics.median(cdf[!, "depth"])
        stddev_coverage = Statistics.std(cdf[!, "depth"])

        # Determine subplot title
        chromosome = cdf[1, "chromosome"]
        subplot_title = isnothing(entity) ? chromosome : "$(entity): $(chromosome)"

        # Create axis for this chromosome
        ax = CairoMakie.Axis(
            fig[i, 1],
            limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
            title = subplot_title,
            xlabel = i == num_plots ? "chromosome index" : "",  # Only show xlabel on bottom plot
            ylabel = "depth",
            titlesize = 16,
            xlabelsize = 14,
            ylabelsize = 14,
            xlabelpadding = 10,
            ylabelpadding = 10,
            xgridvisible = false,
            ygridvisible = true,
            # xgridcolor = (:gray, 0.3),
            ygridcolor = (:gray, 0.3),
            topspinevisible = false,
            rightspinevisible = false
        )

        # Add highlighted regions for 3σ outliers
        is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
        for r in find_true_ranges(is_3sigma_above; min_length=1000)
            CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 01.0))
            println("$(r[1]) - $(r[2]) > +3σ")
        end
        
        is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
        for r in find_true_ranges(is_3sigma_below; min_length=1000)
            CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 1.0))
            println("$(r[1]) - $(r[2]) < -3σ")
        end

        # Create rainbow color scheme
        color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
        
        # Plot main coverage line
        CairoMakie.lines!(ax, 
            equally_spaced_samples(cdf[!, "index"], 10_000), 
            equally_spaced_samples(cdf[!, "depth"], 10_000), 
            color = :black,
            linewidth = 1.5,
            label = "coverage"
        )

        # Add horizontal threshold lines
        threshold_lines = [
            (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
            (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
            (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
            (mean_coverage, "mean_coverage", color_vec[5]),
            (median_coverage, "median", color_vec[4]),
            (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
            (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
            (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
        ]
        
        for (y_val, label, color) in threshold_lines
            CairoMakie.hlines!(ax, y_val, 
                color = color, 
                linewidth = 2,
                linestyle = :dash,
                label = label
            )
        end

        # Create individual legend for this subplot
        CairoMakie.Legend(
            fig[i, 2], 
            ax, 
            framevisible = true,
            backgroundcolor = (:white, 0.9),
            labelsize = 12,
            rowgap = 2,
            margin = (5, 5, 5, 5),
            valign = :center,
            halign = :left
        )
    end

    # Now that the layout is created, set column ratios
    CairoMakie.colsize!(fig.layout, 1, CairoMakie.Relative(0.75))  # 75% for plots
    CairoMakie.colsize!(fig.layout, 2, CairoMakie.Relative(0.25))  # 25% for legends

    # Adjust layout to prevent cutoff
    CairoMakie.resize_to_layout!(fig)
    
    return fig
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a multi-panel visualization of genome coverage across chromosomes.

# # Arguments
# - `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# # Keyword Arguments
# - `entity`: Optional string to specify the entity being analyzed. If provided, chromosome titles become "\$(entity): \$(chromosome)", otherwise just "\$(chromosome)"

# # Returns
# - `CairoMakie.Figure`: A composite figure with coverage plots for each chromosome

# # Details
# Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
# distribution across genomic positions for that chromosome.
# """
# function visualize_genome_coverage(coverage_table; entity::Union{String, Nothing} = nothing)
#     chromosome_groups = DataFrames.groupby(coverage_table, "chromosome")
#     num_plots = length(chromosome_groups)
    
#     # Create main figure with appropriate size and spacing
#     fig = CairoMakie.Figure(
#         size = (1400, 400 * num_plots),
#         fontsize = 12,
#         figure_padding = (20, 40, 20, 20)  # left, right, bottom, top
#     )
    
#     # Set column ratios - give much more space to plots than legends
#     CairoMakie.colsize!(fig.layout, 1, CairoMakie.Relative(0.75))  # 75% for plots
#     CairoMakie.colsize!(fig.layout, 2, CairoMakie.Relative(0.25))  # 25% for legends
    
#     for (i, cdf) in enumerate(chromosome_groups)
#         mean_coverage = Statistics.mean(cdf[!, "depth"])
#         median_coverage = Statistics.median(cdf[!, "depth"])
#         stddev_coverage = Statistics.std(cdf[!, "depth"])

#         # Determine subplot title
#         chromosome = cdf[1, "chromosome"]
#         subplot_title = isnothing(entity) ? chromosome : "$(entity): $(chromosome)"

#         # Create axis for this chromosome
#         ax = CairoMakie.Axis(
#             fig[i, 1],
#             limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#             title = subplot_title,
#             xlabel = i == num_plots ? "chromosome index" : "",  # Only show xlabel on bottom plot
#             ylabel = "depth",
#             titlesize = 14,
#             xlabelpadding = 10,
#             ylabelpadding = 10,
#             xgridvisible = true,
#             ygridvisible = true,
#             xgridcolor = (:gray, 0.3),
#             ygridcolor = (:gray, 0.3),
#             topspinevisible = false,
#             rightspinevisible = false
#         )

#         # Add highlighted regions for 3σ outliers
#         is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_above; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.1))
#         end
        
#         is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_below; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.33))
#         end

#         # Create rainbow color scheme
#         color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
        
#         # Plot main coverage line
#         CairoMakie.lines!(ax, 
#             equally_spaced_samples(cdf[!, "index"], 10_000), 
#             equally_spaced_samples(cdf[!, "depth"], 10_000), 
#             color = :black,
#             linewidth = 1.5,
#             label = "coverage"
#         )

#         # Add horizontal threshold lines
#         threshold_lines = [
#             (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
#             (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
#             (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
#             (mean_coverage, "mean_coverage", color_vec[5]),
#             (median_coverage, "median", color_vec[4]),
#             (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
#             (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
#             (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
#         ]
        
#         for (y_val, label, color) in threshold_lines
#             CairoMakie.hlines!(ax, y_val, 
#                 color = color, 
#                 linewidth = 2,
#                 linestyle = :dash,
#                 label = label
#             )
#         end

#         # Create individual legend for this subplot with smaller margins
#         CairoMakie.Legend(
#             fig[i, 2], 
#             ax, 
#             framevisible = true,
#             backgroundcolor = (:white, 0.9),
#             labelsize = 9,
#             rowgap = 2,
#             margin = (5, 5, 5, 5),  # Reduced margins
#             valign = :center,
#             halign = :left
#         )
#     end

#     # Adjust layout to prevent cutoff
#     CairoMakie.resize_to_layout!(fig)
    
#     return fig
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a multi-panel visualization of genome coverage across chromosomes.

# # Arguments
# - `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# # Keyword Arguments
# - `entity`: Optional string to specify the entity being analyzed. Default title is "Genome coverage", with entity it becomes "Genome coverage of \$(entity)"

# # Returns
# - `CairoMakie.Figure`: A composite figure with coverage plots for each chromosome

# # Details
# Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
# distribution across genomic positions for that chromosome.
# """
# function visualize_genome_coverage(coverage_table; entity::Union{String, Nothing} = nothing)
#     chromosome_groups = DataFrames.groupby(coverage_table, "chromosome")
#     num_plots = length(chromosome_groups)
    
#     # Determine figure title
#     fig_title = isnothing(entity) ? "Genome coverage" : "Genome coverage of $(entity)"
    
#     # Create main figure with appropriate size and spacing
#     fig = CairoMakie.Figure(
#         size = (1400, 400 * num_plots + 60),  # Extra space for main title
#         fontsize = 12,
#         figure_padding = (20, 120, 20, 40)  # left, right, bottom, top (increased top for title)
#     )
    
#     # Add main figure title
#     CairoMakie.Label(fig[0, :], fig_title, fontsize = 18, font = :bold)
    
#     for (i, cdf) in enumerate(chromosome_groups)
#         mean_coverage = Statistics.mean(cdf[!, "depth"])
#         median_coverage = Statistics.median(cdf[!, "depth"])
#         stddev_coverage = Statistics.std(cdf[!, "depth"])

#         # Create axis for this chromosome
#         ax = CairoMakie.Axis(
#             fig[i, 1],
#             limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#             title = cdf[1, "chromosome"],
#             xlabel = i == num_plots ? "chromosome index" : "",  # Only show xlabel on bottom plot
#             ylabel = "depth",
#             titlesize = 14,
#             xlabelpadding = 10,
#             ylabelpadding = 10,
#             xgridvisible = true,
#             ygridvisible = true,
#             xgridcolor = (:gray, 0.3),
#             ygridcolor = (:gray, 0.3),
#             topspinevisible = false,
#             rightspinevisible = false
#         )

#         # Add highlighted regions for 3σ outliers
#         is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_above; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.1))
#         end
        
#         is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_below; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.33))
#         end

#         # Create rainbow color scheme
#         color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
        
#         # Plot main coverage line
#         CairoMakie.lines!(ax, 
#             equally_spaced_samples(cdf[!, "index"], 10_000), 
#             equally_spaced_samples(cdf[!, "depth"], 10_000), 
#             color = :black,
#             linewidth = 1.5,
#             label = "coverage"
#         )

#         # Add horizontal threshold lines
#         threshold_lines = [
#             (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
#             (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
#             (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
#             (mean_coverage, "mean_coverage", color_vec[5]),
#             (median_coverage, "median", color_vec[4]),
#             (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
#             (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
#             (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
#         ]
        
#         for (y_val, label, color) in threshold_lines
#             CairoMakie.hlines!(ax, y_val, 
#                 color = color, 
#                 linewidth = 2,
#                 linestyle = :dash,
#                 label = label
#             )
#         end

#         # Create individual legend for this subplot
#         CairoMakie.Legend(
#             fig[i, 2], 
#             ax, 
#             framevisible = true,
#             backgroundcolor = (:white, 0.9),
#             labelsize = 10,
#             rowgap = 3,
#             margin = (10, 10, 10, 10),
#             valign = :center
#         )
#     end

#     # Adjust layout to prevent cutoff
#     CairoMakie.resize_to_layout!(fig)
    
#     return fig
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a multi-panel visualization of genome coverage across chromosomes.

# # Arguments
# - `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# # Keyword Arguments
# - `entity`: Optional string to specify the entity being analyzed. Default title is "Genome coverage", with entity it becomes "Genome coverage of \$(entity)"

# # Returns
# - `CairoMakie.Figure`: A composite figure with coverage plots for each chromosome

# # Details
# Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
# distribution across genomic positions for that chromosome.
# """
# function visualize_genome_coverage(coverage_table; entity::Union{String, Nothing} = nothing)
#     chromosome_groups = DataFrames.groupby(coverage_table, "chromosome")
#     num_plots = length(chromosome_groups)
    
#     # Determine figure title
#     fig_title = isnothing(entity) ? "Genome coverage" : "Genome coverage of $(entity)"
    
#     # Create main figure with appropriate size and spacing
#     fig = CairoMakie.Figure(
#         size = (1400, 400 * num_plots + 60),  # Extra space for main title
#         fontsize = 12,
#         figure_padding = (20, 120, 20, 40)  # left, right, bottom, top (increased top for title)
#     )
    
#     # Add main figure title
#     CairoMakie.Label(fig[0, :], fig_title, fontsize = 18, font = :bold)
    
#     for (i, cdf) in enumerate(chromosome_groups)
#         mean_coverage = Statistics.mean(cdf[!, "depth"])
#         median_coverage = Statistics.median(cdf[!, "depth"])
#         stddev_coverage = Statistics.std(cdf[!, "depth"])

#         # Create axis for this chromosome
#         ax = CairoMakie.Axis(
#             fig[i, 1],
#             limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#             title = cdf[1, "chromosome"],
#             xlabel = i == num_plots ? "chromosome index" : "",  # Only show xlabel on bottom plot
#             ylabel = "depth",
#             titlesize = 14,
#             xlabelpadding = 10,
#             ylabelpadding = 10,
#             xgridvisible = true,
#             ygridvisible = true,
#             xgridcolor = (:gray, 0.3),
#             ygridcolor = (:gray, 0.3),
#             topspinevisible = false,
#             rightspinevisible = false
#         )

#         # Add highlighted regions for 3σ outliers
#         is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_above; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.1))
#         end
        
#         is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_below; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.33))
#         end

#         # Create rainbow color scheme
#         color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
        
#         # Plot main coverage line
#         CairoMakie.lines!(ax, 
#             equally_spaced_samples(cdf[!, "index"], 10_000), 
#             equally_spaced_samples(cdf[!, "depth"], 10_000), 
#             color = :black,
#             linewidth = 1.5,
#             label = "coverage"
#         )

#         # Add horizontal threshold lines
#         threshold_lines = [
#             (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
#             (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
#             (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
#             (mean_coverage, "mean_coverage", color_vec[5]),
#             (median_coverage, "median", color_vec[4]),
#             (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
#             (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
#             (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
#         ]
        
#         for (y_val, label, color) in threshold_lines
#             CairoMakie.hlines!(ax, y_val, 
#                 color = color, 
#                 linewidth = 2,
#                 linestyle = :dash,
#                 label = label
#             )
#         end

#         # Create individual legend for this subplot
#         CairoMakie.Legend(
#             fig[i, 2], 
#             ax, 
#             framevisible = true,
#             backgroundcolor = (:white, 0.9),
#             labelsize = 10,
#             rowgap = 3,
#             margin = (10, 10, 10, 10),
#             valign = :center
#         )
#     end

#     # Adjust layout to prevent cutoff
#     CairoMakie.resize_to_layout!(fig)
    
#     return fig
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a multi-panel visualization of genome coverage across chromosomes.

# # Arguments
# - `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# # Returns
# - `CairoMakie.Figure`: A composite figure with coverage plots for each chromosome

# # Details
# Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
# distribution across genomic positions for that chromosome.
# """
# function visualize_genome_coverage(coverage_table)
#     chromosome_groups = DataFrames.groupby(coverage_table, "chromosome")
#     num_plots = length(chromosome_groups)
    
#     # Create main figure with appropriate size and spacing
#     fig = CairoMakie.Figure(
#         size = (1400, 400 * num_plots),
#         fontsize = 12,
#         figure_padding = (20, 120, 20, 20)  # left, right, bottom, top
#     )
    
#     for (i, cdf) in enumerate(chromosome_groups)
#         mean_coverage = Statistics.mean(cdf[!, "depth"])
#         median_coverage = Statistics.median(cdf[!, "depth"])
#         stddev_coverage = Statistics.std(cdf[!, "depth"])

#         # Create axis for this chromosome
#         ax = CairoMakie.Axis(
#             fig[i, 1],
#             limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#             title = cdf[1, "chromosome"],
#             xlabel = i == num_plots ? "chromosome index" : "",  # Only show xlabel on bottom plot
#             ylabel = "depth",
#             titlesize = 14,
#             xlabelpadding = 10,
#             ylabelpadding = 10,
#             xgridvisible = true,
#             ygridvisible = true,
#             xgridcolor = (:gray, 0.3),
#             ygridcolor = (:gray, 0.3),
#             topspinevisible = false,
#             rightspinevisible = false
#         )

#         # Add highlighted regions for 3σ outliers
#         is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_above; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.1))
#         end
        
#         is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_below; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.33))
#         end

#         # Create rainbow color scheme
#         color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
        
#         # Plot main coverage line
#         CairoMakie.lines!(ax, 
#             equally_spaced_samples(cdf[!, "index"], 10_000), 
#             equally_spaced_samples(cdf[!, "depth"], 10_000), 
#             color = :black,
#             linewidth = 1.5,
#             label = "coverage"
#         )

#         # Add horizontal threshold lines
#         threshold_lines = [
#             (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
#             (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
#             (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
#             (mean_coverage, "mean_coverage", color_vec[5]),
#             (median_coverage, "median", color_vec[4]),
#             (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
#             (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
#             (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
#         ]
        
#         for (y_val, label, color) in threshold_lines
#             CairoMakie.hlines!(ax, y_val, 
#                 color = color, 
#                 linewidth = 2,
#                 linestyle = :dash,
#                 label = label
#             )
#         end

#         # Create individual legend for this subplot
#         CairoMakie.Legend(
#             fig[i, 2], 
#             ax, 
#             framevisible = true,
#             backgroundcolor = (:white, 0.9),
#             labelsize = 10,
#             rowgap = 3,
#             margin = (10, 10, 10, 10),
#             valign = :center
#         )
#     end

#     # Adjust layout to prevent cutoff
#     CairoMakie.resize_to_layout!(fig)
    
#     return fig
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a multi-panel visualization of genome coverage across chromosomes.

# # Arguments
# - `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# # Returns
# - `CairoMakie.Figure`: A composite figure with coverage plots for each chromosome

# # Details
# Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
# distribution across genomic positions for that chromosome.
# """
# function visualize_genome_coverage(coverage_table)
#     chromosome_groups = DataFrames.groupby(coverage_table, "chromosome")
#     num_plots = length(chromosome_groups)
    
#     # Create main figure with appropriate size and spacing
#     fig = CairoMakie.Figure(
#         resolution = (1200, 400 * num_plots),
#         fontsize = 12,
#         figure_padding = (20, 100, 20, 20)  # left, right, bottom, top
#     )
    
#     # Store all axes for shared legend
#     all_axes = []
    
#     for (i, cdf) in enumerate(chromosome_groups)
#         mean_coverage = Statistics.mean(cdf[!, "depth"])
#         median_coverage = Statistics.median(cdf[!, "depth"])
#         stddev_coverage = Statistics.std(cdf[!, "depth"])

#         # Create axis for this chromosome
#         ax = CairoMakie.Axis(
#             fig[i, 1],
#             limits = (extrema(cdf[!, "index"])..., 0, max(maximum(cdf[!, "depth"]), mean_coverage + 3 * stddev_coverage) * 1.1),
#             title = cdf[1, "chromosome"],
#             xlabel = i == num_plots ? "chromosome index" : "",  # Only show xlabel on bottom plot
#             ylabel = "depth",
#             titlesize = 14,
#             xlabelpadding = 10,
#             ylabelpadding = 10,
#             xgridvisible = true,
#             ygridvisible = true,
#             xgridcolor = (:gray, 0.3),
#             ygridcolor = (:gray, 0.3),
#             topspinevisible = false,
#             rightspinevisible = false
#         )
        
#         push!(all_axes, ax)

#         # Add highlighted regions for 3σ outliers
#         is_3sigma_above = cdf[!, "depth"] .> mean_coverage + 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_above; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.1))
#         end
        
#         is_3sigma_below = cdf[!, "depth"] .< mean_coverage - 3 * stddev_coverage
#         for r in find_true_ranges(is_3sigma_below; min_length=1000)
#             CairoMakie.vspan!(ax, r[1], r[2], color = (:red, 0.33))
#         end

#         # Create rainbow color scheme
#         color_vec = CairoMakie.cgrad(:rainbow, 8, categorical = true)
        
#         # Plot main coverage line
#         CairoMakie.lines!(ax, 
#             equally_spaced_samples(cdf[!, "index"], 10_000), 
#             equally_spaced_samples(cdf[!, "depth"], 10_000), 
#             color = :black,
#             linewidth = 1.5,
#             label = "coverage"
#         )

#         # Add horizontal threshold lines
#         threshold_lines = [
#             (mean_coverage + 3 * stddev_coverage, "+3σ", color_vec[8]),
#             (mean_coverage + 2 * stddev_coverage, "+2σ", color_vec[7]),
#             (mean_coverage + 1 * stddev_coverage, "+σ", color_vec[6]),
#             (mean_coverage, "mean_coverage", color_vec[5]),
#             (median_coverage, "median", color_vec[4]),
#             (mean_coverage - 1 * stddev_coverage, "-σ", color_vec[3]),
#             (mean_coverage - 2 * stddev_coverage, "-2σ", color_vec[2]),
#             (mean_coverage - 3 * stddev_coverage, "-3σ", color_vec[1])
#         ]
        
#         for (y_val, label, color) in threshold_lines
#             CairoMakie.hlines!(ax, y_val, 
#                 color = color, 
#                 linewidth = 2,
#                 linestyle = :dash,
#                 label = label
#             )
#         end
#     end

#     # Create shared legend using the first axis (all have same elements)
#     CairoMakie.Legend(
#         fig[1:num_plots, 2], 
#         all_axes[1], 
#         framevisible = true,
#         bgcolor = (:white, 0.9),
#         labelsize = 11,
#         rowgap = 3,
#         margin = (10, 10, 10, 10),
#         valign = :top
#     )

#     # Adjust layout to prevent cutoff
#     CairoMakie.resize_to_layout!(fig)
    
#     return fig
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a multi-panel visualization of genome coverage across chromosomes.

# # Arguments
# - `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# # Returns
# - `Plots.Figure`: A composite figure with coverage plots for each chromosome

# # Details
# Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
# distribution across genomic positions for that chromosome.
# """
# function visualize_genome_coverage(coverage_table)
#     num_plots = length(unique(coverage_table[!, "chromosome"]))
#     meta_figure = StatsPlots.plot(
#         [chromosome_coverage_table_to_plot(cdf) for cdf in DataFrames.groupby(coverage_table, "chromosome")]...,
#         layout = (num_plots, 1),
#         size = (800, 600 * num_plots)) # Adjust size as needed
#     return meta_figure
# end

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
    
#     kmer_counts = MetaGraphsNext.get_prop(graph, :kmer_counts)
    kmers = [MetaGraphsNext.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
    counts = [MetaGraphsNext.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]
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

"""
    plot_per_base_quality(fastq_file::String; max_position::Union{Int,Nothing}=nothing, sample_size::Union{Int,Nothing}=nothing)

Create per-base quality boxplots for FASTQ data, similar to FastQC output.

# Arguments
- `fastq_file::String`: Path to FASTQ file to analyze
- `max_position::Union{Int,Nothing}=nothing`: Maximum read position to plot (default: auto-detect from data)
- `sample_size::Union{Int,Nothing}=nothing`: Number of reads to sample for analysis (default: use all reads)

# Returns
- `Plots.Plot`: Boxplot showing quality distribution at each base position

# Examples
```julia
# Basic per-base quality plot
p = Mycelia.plot_per_base_quality("reads.fastq")

# Plot first 100 positions only, sampling 10000 reads
p = Mycelia.plot_per_base_quality("reads.fastq", max_position=100, sample_size=10000)
```

# Notes
- Quality scores are displayed in Phred scale
- Green zone: Q>=30 (high quality)
- Yellow zone: Q20-29 (medium quality)  
- Red zone: Q<20 (low quality)
- For large files, consider using sample_size to improve performance
"""
function plot_per_base_quality(fastq_file::String; 
                              max_position::Union{Int,Nothing}=nothing, 
                              sample_size::Union{Int,Nothing}=nothing)
    
    if !isfile(fastq_file)
        error("FASTQ file does not exist: $(fastq_file)")
    end

    # Collect quality scores by position
    position_qualities = Dict{Int, Vector{Int}}()
    read_count = 0
    max_read_length = 0

    reader = FASTX.FASTQ.Reader(open(fastq_file, "r"))
    
    try
        for record in reader
            read_count += 1
            
            # Sample reads if requested
            if sample_size !== nothing && read_count > sample_size
                break
            end
            
            quality_scores = FASTX.quality_scores(record)
            read_length = length(quality_scores)
            max_read_length = max(max_read_length, read_length)
            
            # Apply position limit if specified
            end_pos = max_position !== nothing ? min(max_position, read_length) : read_length
            
            for (pos, qual) in enumerate(quality_scores[1:end_pos])
                if !haskey(position_qualities, pos)
                    position_qualities[pos] = Int[]
                end
                push!(position_qualities[pos], qual)
            end
        end
    finally
        close(reader)
    end
    
    if read_count == 0
        error("No reads found in FASTQ file: $(fastq_file)")
    end
    
    # Determine plotting range
    max_pos_to_plot = max_position !== nothing ? max_position : max_read_length
    positions = 1:max_pos_to_plot
    
    # Prepare data for boxplot
    plot_data = []
    plot_positions = []
    
    for pos in positions
        if haskey(position_qualities, pos) && !isempty(position_qualities[pos])
            quals = position_qualities[pos]
            append!(plot_data, quals)
            append!(plot_positions, fill(pos, length(quals)))
        end
    end
    
    if isempty(plot_data)
        error("No quality data found for specified position range")
    end
    
    # Create the boxplot
    p = StatsPlots.boxplot(
        plot_positions,
        plot_data,
        xlabel="Position in Read",
        ylabel="Quality Score (Phred)",
        title="Per-Base Sequence Quality",
        legend=false,
        fillalpha=0.7,
        linewidth=1.5,
        outliers=false  # Don't show outliers to reduce clutter
    )
    
    # Add quality zone background colors
    max_qual = maximum(plot_data)
    min_qual = minimum(plot_data)
    y_range = max_qual - min_qual
    
    # High quality zone (Q>=30) - green background
    if max_qual >= 30
        StatsPlots.plot!(p, [0, max_pos_to_plot+1], [30, 30], 
                        fillrange=[max_qual+1, max_qual+1], 
                        fillalpha=0.1, fillcolor=:green, 
                        linealpha=0, label=nothing)
    end
    
    # Medium quality zone (Q20-29) - yellow background  
    if max_qual >= 20
        y_top = min(29, max_qual)
        StatsPlots.plot!(p, [0, max_pos_to_plot+1], [20, 20], 
                        fillrange=[y_top, y_top], 
                        fillalpha=0.1, fillcolor=:yellow, 
                        linealpha=0, label=nothing)
    end
    
    # Low quality zone (Q<20) - red background
    if min_qual < 20
        y_top = min(19, max_qual)
        StatsPlots.plot!(p, [0, max_pos_to_plot+1], [min_qual-1, min_qual-1], 
                        fillrange=[y_top, y_top], 
                        fillalpha=0.1, fillcolor=:red, 
                        linealpha=0, label=nothing)
    end
    
    # Add quality threshold lines
    StatsPlots.hline!(p, [20], linestyle=:dash, color=:orange, linewidth=1, alpha=0.7, label=nothing)
    StatsPlots.hline!(p, [30], linestyle=:dash, color=:green, linewidth=1, alpha=0.7, label=nothing)
    
    # Set reasonable axis limits
    StatsPlots.plot!(p, 
                    xlims=(0.5, max_pos_to_plot + 0.5),
                    ylims=(max(0, min_qual - 2), max_qual + 2))
    
    # Add summary information
    avg_qual_per_pos = [Statistics.mean(position_qualities[pos]) for pos in positions if haskey(position_qualities, pos)]
    overall_avg_qual = Statistics.mean(avg_qual_per_pos)
    
    StatsPlots.annotate!(p, max_pos_to_plot * 0.7, max_qual * 0.95, 
                        Plots.text("Avg Quality: $(round(overall_avg_qual, digits=1))", 10, :black))
    StatsPlots.annotate!(p, max_pos_to_plot * 0.7, max_qual * 0.90, 
                        Plots.text("Reads: $(read_count)", 10, :black))
    
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Visualize biological diversity saturation (e.g., k-mer accumulation) against reference benchmarks.

Designed to track the accumulation of unique biological features (k-mers, genes, taxa) as more samples are sequenced or analyzed. 
This visualization helps assess "saturation"—whether sampling more data yields diminishing returns in novel diversity.

# Arguments
- `diversity_counts`: The primary data. Can be:
    - A single `Vector{<:Number}` (one series).
    - A `Vector{Vector{<:Number}}` (multiple series, e.g., different experimental cohorts).

# Keyword Arguments
## Data & Grouping
- `x_values`: X-axis values (e.g., "Genomes Processed"). Defaults to `1:length`.
- `labels`: Label(s) for the diversity count series (for the legend).
- `grouping_values`: A vector aligned with the x-axis indicating distinct groups (e.g., "Discovery Wave", "Year"). 
   If provided, the line color changes at group transitions, and vertical separators are added.
   *Note: If `grouping_values` is provided, gradient coloring applies to these segments.*

## Benchmarks
- `reference_values`: Vector of scalar numbers representing horizontal benchmark lines (e.g., "Total RefSeq Diversity").
- `reference_labels`: Labels corresponding to `reference_values`.
- `reference_colors`: Colors for `reference_values`. Defaults to dark gray.

## Aesthetics
- `grouping_colors`: Vector of colors for the segments defined by `grouping_values`. Defaults to a blue-gradient.
- `title`: Plot title.
- `xlabel`: X-axis label.
- `ylabel`: Y-axis label.
- `figure_size`: Tuple (width, height).
"""
function plot_diversity_saturation(
    diversity_counts::Union{AbstractVector{<:Real}, AbstractVector{<:AbstractVector{<:Real}}};
    # Data context
    x_values::Union{Nothing, AbstractVector} = nothing,
    labels::Union{Nothing, String, Vector{String}} = nothing,
    # Segmentation / Waves
    grouping_values::Union{Nothing, AbstractVector} = nothing,
    grouping_colors::Union{Nothing, AbstractVector} = nothing,
    # Reference Lines
    reference_values::AbstractVector{<:Real} = Float64[],
    reference_labels::AbstractVector{String} = String[],
    reference_colors::Union{Nothing, AbstractVector} = nothing,

    # Styling
    title::String = "Diversity Saturation Analysis",
    xlabel::String = "Samples Processed",
    ylabel::String = "Cumulative Unique Features",
    figure_size::Tuple{Int, Int} = (1200, 700),
    linewidth::Real = 3.0,
    base_color = :blue # Default color if no grouping is used
  )
  # Normalize input to vector of vectors
  series_list = eltype(diversity_counts) <: Real ? [diversity_counts] : diversity_counts
  n_series = length(series_list)
  # Normalize labels
  series_labels = if isnothing(labels)
      n_series == 1 ? ["Diversity"] : ["Series $i" for i in 1:n_series]
  elseif labels isa String
      [labels]
  else
      labels
  end

  # Determine X values (assume all series share the same X if not specified, or match length)
  # If x_values not provided, use 1:N of the first series
  xs = isnothing(x_values) ? collect(1:length(series_list[1])) : x_values

  # Setup Colors
  # Default wave colors: Progressive gradient from Drab/Gray (Wave 1) to Vivid Blue (Wave N)
  # If user didn't provide colors, we generate them based on the number of unique groups
  unique_groups = isnothing(grouping_values) ? [] : sort(unique(grouping_values))
  n_groups = length(unique_groups)

  segment_palette = if !isnothing(grouping_colors)
      grouping_colors
  elseif n_groups > 0
      # Generate a custom blue-ish gradient similar to the user's example
      # Interpolate between Gray-Blue and Vivid Blue
      c1 = Colors.RGB(0.65, 0.67, 0.70)
      c2 = Colors.RGB(0.00, 0.30, 0.90)
      Colors.range(c1, stop=c2, length=n_groups)
  else
      [base_color]
  end

  # Reference colors
  ref_palette = isnothing(reference_colors) ? 
      [Colors.RGB(0.4, 0.4, 0.4) for _ in 1:length(reference_values)] : 
      reference_colors

  # --- Plotting ---
  fig = CairoMakie.Figure(size = figure_size, fontsize = 14, backgroundcolor = :white)

  ax = CairoMakie.Axis(
      fig[1, 1],
      xlabel = xlabel,
      ylabel = ylabel,
      title = title,
      titlesize = 18,
      xlabelsize = 16,
      ylabelsize = 16,
      xticklabelsize = 12,
      yticklabelsize = 12,
      xgridvisible = false, # Turn OFF vertical grid lines (standard for saturation plots)
      ygridvisible = true,
      ygridcolor = Colors.RGB(0.9, 0.9, 0.9)
  )

  legend_elements = []
  legend_entries = String[]

  # 1. Plot Data Series
  for (si, y_data) in enumerate(series_list)
      current_x = xs[1:length(y_data)] # Handle if x is longer than y
      
      if !isnothing(grouping_values) && n_groups > 0
          # --- Segmented Plotting ---
          # We assume grouping_values corresponds to the X axis
          
          # Find transitions for vertical lines
          transitions = Int[]
          for i in 2:length(grouping_values)
              if grouping_values[i] != grouping_values[i-1]
                  push!(transitions, i)
              end
          end
          
          # Add vertical lines at transitions
          if !isempty(transitions)
              CairoMakie.vlines!(ax, current_x[transitions],
                  color = Colors.RGB(0.5, 0.5, 0.5),
                  linewidth = 1.5,
                  linestyle = :dash
              )
          end

          # Plot segments
          for (gi, group_val) in enumerate(unique_groups)
              mask = grouping_values .== group_val
              indices = findall(mask)
              
              if isempty(indices); continue; end

              # Extend segment to overlap with next for continuity (visual gap closing)
              if maximum(indices) < length(current_x)
                  # Only extend if the NEXT point exists in data
                  push!(indices, maximum(indices) + 1)
              end
              
              # Safety for color indexing
              c_idx = mod1(gi, length(segment_palette))
              color = segment_palette[c_idx]
              
              CairoMakie.lines!(ax, current_x[indices], y_data[indices],
                  color = color,
                  linewidth = linewidth
              )
              
              # Add to legend only once (for the first series, if multiple)
              if si == 1
                  push!(legend_elements, CairoMakie.LineElement(color = color, linewidth = linewidth))
                  push!(legend_entries, "Group: $group_val")
              end
          end
      else
          # --- Continuous Plotting ---
          # If no grouping, just plot the line
          color = n_series > 1 ? Mycelia.n_maximally_distinguishable_colors(n_series)[si] : segment_palette[1]
          
          CairoMakie.lines!(ax, current_x, y_data,
              color = color,
              linewidth = linewidth
          )
          
          push!(legend_elements, CairoMakie.LineElement(color = color, linewidth = linewidth))
          push!(legend_entries, series_labels[si])
      end
  end

  # 2. Plot References
  for (i, val) in enumerate(reference_values)
      if i > length(reference_labels); break; end
      
      c = ref_palette[mod1(i, length(ref_palette))]
      
      CairoMakie.hlines!(ax, [val],
          color = c,
          linewidth = linewidth,
          linestyle = :dashdot
      )
      
      push!(legend_elements, CairoMakie.LineElement(color = c, linewidth = linewidth, linestyle = :dashdot))
      push!(legend_entries, reference_labels[i])
  end

  # 3. Create Legend
  if !isempty(legend_elements)
      CairoMakie.Legend(
          fig[1, 2],
          legend_elements,
          legend_entries,
          framevisible = true,
          framecolor = Colors.RGB(0.8, 0.8, 0.8),
          backgroundcolor = Colors.RGB(0.98, 0.98, 0.98),
          labelsize = 12,
          rowgap = 5
      )
  end

  CairoMakie.resize_to_layout!(fig)
  return fig
end

import CairoMakie
import Statistics
import LinearAlgebra
import OrderedCollections

# --- Data Structures ---

"""
    PointSeries(indices, label, color; marker=:circle, size=10.0, z_order=1)

Defines a group of points that share a common visual style (e.g., "Environmental Isolates").
Used for coloring points and creating the source legend.
"""
struct PointSeries
    indices::Vector{Int}
    label::String
    color::Any             # Symbol (e.g., :red) or String hex or RGB object
    marker::Symbol         # Default marker for this group
    size::Float64
    z_order::Int           # Higher numbers plot on top of lower numbers
end

# Constructor with keywords
PointSeries(idx, lbl, col; marker=:circle, size=10.0, z_order=1) = 
    PointSeries(idx, lbl, col, marker, size, z_order)

"""
    EllipseConfig(indices, label, color; anchor=:top, align=(:center, :bottom), offset=(0.0, 5.0))

Defines a target for drawing a confidence ellipse and its label.
"""
struct EllipseConfig
    indices::Vector{Int}
    label::String          
    color::Any             
    label_anchor::Symbol   # :top, :bottom, :left, :right
    label_align::Tuple{Symbol, Symbol} 
    label_offset::Tuple{Float64, Float64}
end

# Constructor with keywords
function EllipseConfig(indices, label, color; 
                       anchor=:top, 
                       align=(:center, :bottom), 
                       offset=(0.0, 5.0))
    return EllipseConfig(indices, label, color, anchor, align, offset)
end


# --- Helper Functions ---

"""
    get_ellipse_points(x_coords, y_coords; scale=3.5)

Calculates the boundary points of a confidence ellipse for a set of 2D coordinates.
Returns a vector of `CairoMakie.Point2f` and a dictionary of anchor points for labeling.
"""
function get_ellipse_points(x_coords, y_coords; scale=3.5)
    # Need at least 3 points to define a covariance
    if length(x_coords) < 3
        return CairoMakie.Point2f[], Dict{Symbol, CairoMakie.Point2f}()
    end
    
    # Calculate Covariance and Mean
    cov_matrix = Statistics.cov(hcat(x_coords, y_coords))
    mean_vec = [Statistics.mean(x_coords), Statistics.mean(y_coords)]
    
    # Eigen decomposition for axis lengths and rotation
    eigen_decomp = LinearAlgebra.eigen(cov_matrix)
    vals, vecs = eigen_decomp.values, eigen_decomp.vectors
    
    # Calculate angle and axis lengths
    angle = atan(vecs[2, 2], vecs[1, 2])
    major_axis = scale * sqrt(vals[2])
    minor_axis = scale * sqrt(vals[1])
    
    # Generate parametric points
    t = range(0, 2π, length=100)
    ellipse_x = major_axis .* cos.(t)
    ellipse_y = minor_axis .* sin.(t)
    
    # Rotate and translate points
    rotated_x = ellipse_x .* cos(angle) .- ellipse_y .* sin(angle) .+ mean_vec[1]
    rotated_y = ellipse_x .* sin(angle) .+ ellipse_y .* cos(angle) .+ mean_vec[2]
    
    points = [CairoMakie.Point2f(x, y) for (x, y) in zip(rotated_x, rotated_y)]
    
    # Define anchors for text placement
    anchors = Dict(
        :top    => points[argmax(rotated_y)],
        :bottom => points[argmin(rotated_y)],
        :left   => points[argmin(rotated_x)],
        :right  => points[argmax(rotated_x)]
    )
    
    return points, anchors
end


# --- Main Plotting Function ---

"""
    plot_generalized_pcoa(x_all, y_all, series_list; ...)

A generalized visualizer for PCoA/PCA data within the Mycelia ecosystem.
Separates data logic from drawing logic using `PointSeries` and `EllipseConfig`.

# Arguments
- `x_all`, `y_all`: Vectors containing coordinates for all points in the dataset.
- `series_list`: Vector of `PointSeries`. Rules for how to color/group specific indices.

# Keywords
- `species_map`: Dict mapping indices to species names (used to assign shapes).
- `marker_map`: Dict mapping species names to Makie symbols (e.g. :rect, :circle).
- `ellipses`: Vector of `EllipseConfig`.
- `title`, `xlabel`, `ylabel`.
- `layout_target`: Optional Makie layout location (e.g. `fig[1,1]`) for subplots.
"""
function plot_generalized_pcoa(
    x_all::AbstractVector, 
    y_all::AbstractVector, 
    series_list::Vector{PointSeries};
    species_map::Dict = Dict(), # Index -> String
    marker_map::AbstractDict = Dict(), # String -> Symbol
    ellipses::Vector{EllipseConfig} = EllipseConfig[],
    title::String = "Principal Coordinate Analysis",
    xlabel::String = "PC1",
    ylabel::String = "PC2",
    output_file::String = "",
    figure_size = (1000, 750),
    layout_target = nothing
)
    # 1. Initialize Figure or Layout
    if isnothing(layout_target)
        fig = CairoMakie.Figure(size=figure_size)
        ax = CairoMakie.Axis(
            fig[1, 1], 
            title=title, 
            xlabel=xlabel, 
            ylabel=ylabel,
            titlesize=16, xlabelsize=14, ylabelsize=14
        )
        legend_layout = fig[1, 2] = CairoMakie.GridLayout()
    else
        # Subplot mode
        ax = CairoMakie.Axis(
            layout_target, 
            title=title, 
            xlabel=xlabel, 
            ylabel=ylabel
        )
        legend_layout = nothing # Legends handled externally or manually for subplots
    end

    # 2. Pre-allocate visual property vectors
    n = length(x_all)
    # Default styling for points not covered by any series
    plot_colors = Vector{Any}(fill(:gray90, n))
    plot_markers = Vector{Symbol}(fill(:circle, n))
    plot_sizes = Vector{Float64}(fill(6.0, n))
    plot_z_orders = Vector{Int}(fill(0, n))
    
    # 3. Apply Series Configurations
    # Last series in the list takes precedence if indices overlap
    for series in series_list
        for idx in series.indices
            if 1 <= idx <= n
                plot_colors[idx] = series.color
                plot_sizes[idx] = series.size
                plot_z_orders[idx] = series.z_order
                
                # Assign marker:
                # 1. Default to series marker
                plot_markers[idx] = series.marker
                
                # 2. If species logic exists, override it
                if !isempty(species_map) && haskey(species_map, idx)
                    sp = species_map[idx]
                    if haskey(marker_map, sp)
                        plot_markers[idx] = marker_map[sp]
                    end
                end
            end
        end
    end

    # 4. Sort indices by z-order to ensure correct layering
    # Stable sort preserves order of series definition for equal z-orders
    draw_order = sortperm(plot_z_orders)
    
    x_sorted = x_all[draw_order]
    y_sorted = y_all[draw_order]
    c_sorted = plot_colors[draw_order]
    m_sorted = plot_markers[draw_order]
    s_sorted = plot_sizes[draw_order]

    # 5. Draw the Scatter Plot
    CairoMakie.scatter!(
        ax, 
        x_sorted, 
        y_sorted, 
        color=c_sorted, 
        marker=m_sorted, 
        markersize=s_sorted, 
        strokewidth=0.5,
        strokecolor=:black
    )

    # 6. Draw Ellipses and Labels
    for ell in ellipses
        # Filter indices that exist in the current dataset
        valid_indices = filter(i -> 1 <= i <= n, ell.indices)
        if isempty(valid_indices) continue end
        
        ex = x_all[valid_indices]
        ey = y_all[valid_indices]
        
        pts, anchors = get_ellipse_points(ex, ey)
        
        if !isempty(pts)
            # Draw ellipse outline
            CairoMakie.poly!(
                ax, 
                pts, 
                color=:transparent, 
                strokecolor=ell.color, 
                strokewidth=1.5
            )
            
            # Draw label if anchor exists
            if haskey(anchors, ell.label_anchor)
                CairoMakie.text!(
                    ax, 
                    ell.label,
                    position = anchors[ell.label_anchor],
                    fontsize = 14,
                    color = ell.color,
                    font = :bold,
                    align = ell.label_align,
                    offset = ell.label_offset
                )
            end
        end
    end

    # 7. Set Axis Limits with Padding
    x_span = maximum(x_sorted) - minimum(x_sorted)
    y_span = maximum(y_sorted) - minimum(y_sorted)
    
    CairoMakie.xlims!(ax, minimum(x_sorted) - x_span*0.1, maximum(x_sorted) + x_span*0.1)
    CairoMakie.ylims!(ax, minimum(y_sorted) - y_span*0.1, maximum(y_sorted) + y_span*0.1)

    # 8. Generate Legends (Only if not in subplot mode)
    if !isnothing(legend_layout)
        # -- Source/Group Legend --
        group_elements = [
            CairoMakie.MarkerElement(color=s.color, marker=:circle, markersize=10, strokewidth=0.5, strokecolor=:black) 
            for s in series_list
        ]
        group_labels = [s.label for s in series_list]
        
        CairoMakie.Legend(
            legend_layout[1, 1], 
            group_elements, 
            group_labels, 
            "Source", 
            halign=:left, valign=:top, framevisible=false
        )
        
        # -- Species/Shape Legend (if used) --
        if !isempty(marker_map)
            # Create a sorted list of species for consistent legend order
            sorted_species = sort(collect(keys(marker_map)))
            
            shape_elements = [
                CairoMakie.MarkerElement(color=:black, marker=marker_map[s], markersize=10) 
                for s in sorted_species
            ]
            # Helper to abbreviate Genus
            abbrev_labels = [replace(s, r"^(\w)\w+\s" => s"\1. ") for s in sorted_species]
            
            CairoMakie.Legend(
                legend_layout[2, 1], 
                shape_elements, 
                abbrev_labels, 
                "Species", 
                halign=:left, valign=:top, framevisible=false
            )
        end
        
        if !isempty(marker_map)
            CairoMakie.rowgap!(legend_layout, 1, 20)
        end
    end

    # 9. Save and Return
    if !isempty(output_file) && isnothing(layout_target)
        CairoMakie.save(output_file, fig)
        println("Saved PCoA plot to: $output_file")
    end

    return isnothing(layout_target) ? fig : ax
end

import CairoMakie
import Clustering
import Statistics
import DataFrames
import OrderedCollections

# --- Tree Visualization Helpers ---

mutable struct ClusterNode
    left::Union{ClusterNode, Nothing}
    right::Union{ClusterNode, Nothing}
    height::Float64
    x::Float64
    leaf_index::Int
end

function hclust_to_tree(hclust::Clustering.Hclust)
    n = length(hclust.order)
    merges = hclust.merges
    heights = hclust.heights
    
    nodes = Dict{Int, ClusterNode}()
    
    for i in 1:n
        nodes[-i] = ClusterNode(nothing, nothing, 0.0, 0.0, i)
    end
    
    for i in 1:size(merges, 1)
        left_idx = merges[i, 1]
        right_idx = merges[i, 2]
        height = heights[i]
        
        left_node = nodes[left_idx]
        right_node = nodes[right_idx]
        
        parent_node = ClusterNode(left_node, right_node, height, 0.0, -1)
        nodes[i] = parent_node
    end
    
    return nodes[size(merges, 1)]
end

function assign_x_coordinates!(node::ClusterNode, order::Vector{Int}, position_map::Dict{Int, Float64})
    if isnothing(node.left) && isnothing(node.right)
        # It's a leaf
        leaf_pos = findfirst(==(node.leaf_index), order)
        node.x = Float64(leaf_pos)
        position_map[node.leaf_index] = node.x
    else
        assign_x_coordinates!(node.left, order, position_map)
        assign_x_coordinates!(node.right, order, position_map)
        node.x = (node.left.x + node.right.x) / 2
    end
end

function plot_dendrogram!(ax, node::ClusterNode, scale_func::Function, orientation::Symbol)
    if isnothing(node.left) || isnothing(node.right)
        return
    end
    
    parent_h = scale_func(node.height)
    left_h = scale_func(node.left.height)
    right_h = scale_func(node.right.height)
    
    if orientation == :vertical
        # Standard dendrogram (top/bottom)
        # Vertical lines
        CairoMakie.lines!(ax, [node.left.x, node.left.x], [left_h, parent_h], color=:black, linewidth=1)
        CairoMakie.lines!(ax, [node.right.x, node.right.x], [right_h, parent_h], color=:black, linewidth=1)
        # Horizontal bar
        CairoMakie.lines!(ax, [node.left.x, node.right.x], [parent_h, parent_h], color=:black, linewidth=1)
    else
        # Horizontal dendrogram (left/right)
        # Horizontal lines (height is x-axis now)
        CairoMakie.lines!(ax, [left_h, parent_h], [node.left.x, node.left.x], color=:black, linewidth=1)
        CairoMakie.lines!(ax, [right_h, parent_h], [node.right.x, node.right.x], color=:black, linewidth=1)
        # Vertical bar
        CairoMakie.lines!(ax, [parent_h, parent_h], [node.left.x, node.right.x], color=:black, linewidth=1)
    end
    
    plot_dendrogram!(ax, node.left, scale_func, orientation)
    plot_dendrogram!(ax, node.right, scale_func, orientation)
end

function plot_dendrogram_from_hclust!(ax, hclust::Clustering.Hclust; 
                                      scale_func=sqrt, 
                                      orientation=:vertical, 
                                      cut_height=nothing)
    root = hclust_to_tree(hclust)
    position_map = Dict{Int, Float64}()
    assign_x_coordinates!(root, hclust.order, position_map)
    
    # Ensure min height isn't negative for log scales or sqrt
    min_h = minimum(hclust.heights)
    offset = max(0.001, min_h * 0.1)
    safe_scale(h) = scale_func(h + offset)
    
    plot_dendrogram!(ax, root, safe_scale, orientation)
    
    n = length(hclust.order)
    
    # Set axis limits
    if orientation == :vertical
        CairoMakie.xlims!(ax, 0.5, n + 0.5)
        if !isnothing(cut_height)
            CairoMakie.hlines!(ax, [safe_scale(cut_height)], color=:red, linestyle=:dash, linewidth=2)
        end
    else
        CairoMakie.ylims!(ax, 0.5, n + 0.5)
        CairoMakie.hidedecorations!(ax) # Usually hide axis for side dendrogram
        CairoMakie.xlims!(ax, nothing, nothing) # Let autoscaling handle height width
        if !isnothing(cut_height)
            CairoMakie.vlines!(ax, [safe_scale(cut_height)], color=:red, linestyle=:dash, linewidth=2)
        end
    end
    
    return ax
end

# --- Main Generalized Visualization Function ---

"""
    visualize_hierarchical_heatmap(
        data_matrix::AbstractMatrix;
        
        # Axis Information
        x_labels::Union{Nothing, Vector{String}} = nothing,
        y_labels::Union{Nothing, Vector{String}} = nothing,
        
        # Clustering Information
        x_hclust::Union{Nothing, Clustering.Hclust} = nothing,
        x_cluster_assignments::Union{Nothing, Vector{Int}} = nothing,
        y_hclust::Union{Nothing, Clustering.Hclust} = nothing,
        y_cluster_assignments::Union{Nothing, Vector{Int}} = nothing,
        
        # Visualization Configuration
        color_normalization::Symbol = :row, # :row, :col, :global
        dendrogram_scale_func::Function = sqrt,
        color_func::Function = (val, max_val) -> CairoMakie.RGBAf(0.0, 0.0, 1.0, 0.2 + 0.8 * (val / max_val)),
        
        # Layout & Text
        title::String = "Cluster Heatmap",
        show_values::Bool = true,
        figure_size::Tuple{Int, Int} = (1600, 1000)
    )

Generalized function to visualize a heatmap with optional hierarchical clustering on X and Y axes.
Supports "collapsed" clusters where the heatmap cells represent groups of leaves from the tree.

# Arguments
- `data_matrix`: Matrix of values to plot. If clustering assignments are provided, 
  this should be aggregated (e.g. counts) matching the number of clusters.
  Dimensions should be (N_Y_Items, N_X_Items).
  
- `x_cluster_assignments`: If provided, maps the leaves of `x_hclust` to the columns of `data_matrix`.
  Must be length equal to number of leaves in `x_hclust`.
  
- `color_normalization`: Defines how the maximum value for color scaling is determined.
  - `:row`: Max value per row.
  - `:col`: Max value per column.
  - `:global`: Max value in entire matrix.
"""
function visualize_hierarchical_heatmap(
    data_matrix::AbstractMatrix;
    
    # X-axis (Columns) config
    x_labels::Union{Nothing, Vector{String}} = nothing,
    x_hclust::Union{Nothing, Clustering.Hclust} = nothing,
    x_cluster_assignments::Union{Nothing, Vector{Int}} = nothing,
    
    # Y-axis (Rows) config
    y_labels::Union{Nothing, Vector{String}} = nothing,
    y_hclust::Union{Nothing, Clustering.Hclust} = nothing,
    y_cluster_assignments::Union{Nothing, Vector{Int}} = nothing,
    
    # Options
    color_normalization::Symbol = :row,
    dendrogram_scale_func::Function = sqrt,
    color_func::Function = (val, max_val) -> CairoMakie.RGBAf(0.0, 0.0, 1.0, 0.2 + 0.8 * (val / max(1e-10, max_val))),
    cut_height_x::Union{Nothing, Real} = nothing,
    cut_height_y::Union{Nothing, Real} = nothing,
    
    title::String = "",
    show_values::Bool = true,
    figure_size::Tuple{Int, Int} = (1600, 1000)
)
    # --- 1. Compute Layout Geometries ---
    
    # Helper to calculate cell positions and sizes based on trees/clusters
    function calculate_axis_geometry(n_items, hclust, assignments)
        if isnothing(hclust)
            # No tree: Uniform spacing
            centers = collect(1.0:n_items)
            widths = ones(Float64, n_items)
            ordering = 1:n_items
            total_span = Float64(n_items)
        else
            # Tree exists
            if isnothing(assignments)
                # No collapsing: 1-to-1 mapping
                centers = collect(1.0:n_items)
                widths = ones(Float64, n_items)
                ordering = hclust.order # Heatmap follows tree order
                total_span = Float64(n_items)
            else
                # Collapsing: Clusters have variable widths
                # Identify unique clusters in tree order
                tree_ordered_assignments = assignments[hclust.order]
                
                # We assume the data_matrix is ordered such that column j corresponds to 
                # the j-th unique cluster encountered when traversing the tree leaves.
                # If your data_matrix is ordered differently, you must reorder it before calling this.
                unique_clusters_ordered = unique(tree_ordered_assignments)
                
                centers = Float64[]
                widths = Float64[]
                
                # Calculate span for each cluster
                for cluster_id in unique_clusters_ordered
                    # Find indices in the *ordered* list
                    indices = findall(==(cluster_id), tree_ordered_assignments)
                    min_pos = minimum(indices)
                    max_pos = maximum(indices)
                    
                    push!(centers, (min_pos + max_pos) / 2.0)
                    push!(widths, max_pos - min_pos + 1.0)
                end
                
                ordering = 1:length(unique_clusters_ordered) # 1..N_Clusters
                total_span = length(hclust.order)
            end
        end
        return centers, widths, ordering, total_span
    end

    n_rows, n_cols = size(data_matrix)
    
    # Calculate X geometries
    x_centers, x_widths, x_ordering, x_span = calculate_axis_geometry(
        n_cols, x_hclust, x_cluster_assignments
    )
    
    # Calculate Y geometries
    y_centers, y_widths, y_ordering, y_span = calculate_axis_geometry(
        n_rows, y_hclust, y_cluster_assignments
    )
    
    # --- 2. Setup Figure and Layout ---
    
    fig = CairoMakie.Figure(size = figure_size)
    
    # Layout definition
    # Col 1: Y-Dendrogram (optional)
    # Col 2: Heatmap
    # Col 3: Legend
    # Row 1: X-Dendrogram (optional)
    # Row 2: Heatmap
    
    hm_row, hm_col = 2, 2
    
    # --- 3. Draw Heatmap ---
    
    ax_hm = CairoMakie.Axis(
        fig[hm_row, hm_col],
        title = isnothing(x_hclust) && isnothing(y_hclust) ? title : "",
        titlesize = 24,
        xticklabelrotation = π/4,
        xgridvisible = false,
        ygridvisible = false
    )
    
    # Apply labels if provided
    if !isnothing(x_labels)
        ax_hm.xticks = (x_centers, x_labels)
    end
    
    if !isnothing(y_labels)
        ax_hm.yticks = (y_centers, y_labels)
    end
    
    # Set limits based on the total span of the trees (or count)
    CairoMakie.xlims!(ax_hm, 0.5, x_span + 0.5)
    CairoMakie.ylims!(ax_hm, 0.5, y_span + 0.5)
    
    # Draw Cells
    global_max = maximum(data_matrix)
    
    for r_idx in 1:n_rows
        # Get visual Y coordinates
        y_c = y_centers[r_idx]
        y_w = y_widths[r_idx]
        y_min = y_c - y_w/2
        y_max = y_c + y_w/2
        
        # Calculate normalization factor for this row
        row_max = maximum(data_matrix[r_idx, :])
        
        for c_idx in 1:n_cols
            # Get visual X coordinates
            x_c = x_centers[c_idx]
            x_w = x_widths[c_idx]
            x_min = x_c - x_w/2
            x_max = x_c + x_w/2
            
            val = data_matrix[r_idx, c_idx]
            
            # Determine max for normalization
            norm_max = if color_normalization == :row
                row_max
            elseif color_normalization == :col
                maximum(data_matrix[:, c_idx])
            else
                global_max
            end
            
            if val > 0
                color = color_func(val, norm_max)
                
                CairoMakie.poly!(ax_hm,
                    CairoMakie.Point2f[
                        (x_min, y_min), (x_max, y_min),
                        (x_max, y_max), (x_min, y_max)
                    ],
                    color = color,
                    strokecolor = :black,
                    strokewidth = 0.5
                )
                
                if show_values
                    # Contrast text color check
                    text_color = val > (norm_max * 0.5) ? :white : :black
                    CairoMakie.text!(ax_hm, x_c, y_c, 
                        text = string(val), 
                        align = (:center, :center), 
                        color = text_color,
                        fontsize = 12
                    )
                end
            end
        end
    end
    
    # --- 4. Draw X-Axis Dendrogram (Top) ---
    if !isnothing(x_hclust)
        ax_dendro_x = CairoMakie.Axis(
            fig[hm_row - 1, hm_col],
            ylabel = "Height",
            title = title,
            titlesize = 24,
            xgridvisible = false,
            ygridvisible = false,
            xticklabelsvisible = false,
            xticksvisible = false,
            bottomspinevisible = false
        )
        
        plot_dendrogram_from_hclust!(
            ax_dendro_x, x_hclust, 
            scale_func=dendrogram_scale_func, 
            orientation=:vertical,
            cut_height=cut_height_x
        )
        
        # Link X axis
        CairoMakie.linkxaxes!(ax_hm, ax_dendro_x)
        CairoMakie.rowsize!(fig.layout, hm_row - 1, CairoMakie.Relative(0.2))
    end
    
    # --- 5. Draw Y-Axis Dendrogram (Left) ---
    if !isnothing(y_hclust)
        ax_dendro_y = CairoMakie.Axis(
            fig[hm_row, hm_col - 1],
            xlabel = "Height",
            xgridvisible = false,
            ygridvisible = false,
            yticklabelsvisible = false,
            yticksvisible = false,
            rightspinevisible = false
        )
        
        plot_dendrogram_from_hclust!(
            ax_dendro_y, y_hclust, 
            scale_func=dendrogram_scale_func, 
            orientation=:horizontal,
            cut_height=cut_height_y
        )
        
        # Invert X axis for the side dendrogram so root is left/top depending on preference?
        # Usually roots are far away from the map.
        CairoMakie.xreverse!(ax_dendro_y) 
        
        # Link Y axis
        CairoMakie.linkyaxes!(ax_hm, ax_dendro_y)
        CairoMakie.colsize!(fig.layout, hm_col - 1, CairoMakie.Relative(0.15))
    end
    
    # --- 6. Legend ---
    # Construct a representative legend based on the normalization
    legend_labels = ["20%", "40%", "60%", "80%", "100%"]
    legend_vals = [0.2, 0.4, 0.6, 0.8, 1.0]
    
    # Dummy max for legend generation
    dummy_max = 100.0
    legend_elements = [
        CairoMakie.PolyElement(
            color = color_func(v * dummy_max, dummy_max), 
            strokecolor = :black
        ) for v in legend_vals
    ]
    
    legend_title = if color_normalization == :row
        "Relative Proportion\n(within row)"
    elseif color_normalization == :col
        "Relative Proportion\n(within col)"
    else
        "Relative Proportion\n(global)"
    end
    
    CairoMakie.Legend(
        fig[hm_row, hm_col + 1],
        legend_elements,
        legend_labels,
        legend_title,
        framevisible = true
    )

    return fig
end

# =============================================================================
# Unified Microbiome Visualization System
# =============================================================================
# These structs and functions provide a consistent, publication-quality
# visualization system for microbiome abundance data across large sample sets.
# Handles 300-600+ samples with adaptive sizing and automatic view selection.
# =============================================================================

"""
    AxisOrdering

Specification for ordering samples or taxa in microbiome plots.

# Fields
- `method::Symbol`: Ordering method - `:hclust`, `:preordered`, `:sort`, `:alphabetical`
- `preordered_labels::Union{Vector{String}, Nothing}`: Labels in desired order (for `:preordered`)
- `sort_by::Union{Symbol, Function, Nothing}`: Column or function for sorting (for `:sort`)
- `distance_metric::Symbol`: Distance metric for clustering - `:braycurtis`, `:euclidean`, `:cosine`
- `linkage::Symbol`: Linkage method for clustering - `:single`, `:complete`, `:average`, `:ward`

# Examples
```julia
# Hierarchical clustering with Bray-Curtis (default)
ordering = AxisOrdering()

# Pre-specified order
ordering = AxisOrdering(method=:preordered, preordered_labels=["Sample1", "Sample2", "Sample3"])

# Sort by mean abundance
ordering = AxisOrdering(method=:sort, sort_by=:mean_abundance)

# Alphabetical
ordering = AxisOrdering(method=:alphabetical)
```
"""
Base.@kwdef struct AxisOrdering
    method::Symbol = :hclust
    preordered_labels::Union{Vector{String}, Nothing} = nothing
    sort_by::Union{Symbol, Function, Nothing} = nothing
    distance_metric::Symbol = :braycurtis
    linkage::Symbol = :average
end

"""
    MicrobiomePlotConfig

Configuration for microbiome abundance visualizations.
Provides consistent defaults across all visualization functions.

# Figure Dimensions
- `max_width::Int`: Maximum figure width in pixels (default: 1600)
- `min_height::Int`: Minimum figure height in pixels (default: 800)
- `pixels_per_sample::Int`: Height per sample in portrait mode (default: 12)
- `orientation::Symbol`: `:auto`, `:landscape`, or `:portrait` (default: :auto)

# X-axis Labels
- `min_label_fontsize::Float64`: Minimum label font size (default: 6.0)
- `max_label_fontsize::Float64`: Maximum label font size (default: 12.0)
- `label_rotation::Float64`: Label rotation in degrees (default: 90.0)

# Taxa and Legend
- `top_n_taxa::Int`: Number of top taxa to show individually (default: 25)
- `legend_fontsize::Float64`: Legend text font size (default: 8.0)
- `legend_position::Symbol`: Legend position (default: :right)

# Margins (pixels)
- `left_margin::Int`: Left margin (default: 80)
- `bottom_margin::Int`: Bottom margin (default: 120)
- `right_margin::Int`: Right margin for legend (default: 150)
- `top_margin::Int`: Top margin (default: 40)

# Ordering
- `sample_ordering::AxisOrdering`: How to order samples (default: hierarchical clustering)
- `taxa_ordering::AxisOrdering`: How to order taxa (default: sort by mean abundance)

# Dendrograms
- `show_sample_dendrogram::Bool`: Show dendrogram for samples (default: false)
- `show_taxa_dendrogram::Bool`: Show dendrogram for taxa (default: false)
- `dendrogram_width::Int`: Dendrogram width in pixels (default: 100)
- `color_branches_by_cluster::Bool`: Color branches by cluster assignment (default: false)
- `n_clusters::Union{Int, Nothing}`: Number of clusters for branch coloring (default: nothing)

# Colors
- `color_palette::Symbol`: Color palette - `:maximally_distinguishable`, `:viridis`, etc.
- `custom_colors::Union{Vector, Nothing}`: Custom color vector (default: nothing)

# Large Dataset Handling
- `large_dataset_view::Symbol`: View type - `:auto`, `:tall`, `:paginated`, `:heatmap`, `:all`
- `samples_per_page::Int`: Samples per page in paginated view (default: 150)
- `heatmap_threshold::Int`: Switch to heatmap above this sample count (default: 300)

# Output
- `output_formats::Vector{Symbol}`: Output formats (default: [:png, :svg])
- `dpi::Int`: Output resolution (default: 300)

# Examples
```julia
# Default configuration
config = MicrobiomePlotConfig()

# Portrait orientation for many samples
config = MicrobiomePlotConfig(orientation=:portrait, pixels_per_sample=10)

# Show dendrograms with 5 clusters
config = MicrobiomePlotConfig(
    show_sample_dendrogram=true,
    color_branches_by_cluster=true,
    n_clusters=5
)
```
"""
Base.@kwdef struct MicrobiomePlotConfig
    # Figure dimensions
    max_width::Int = 1600
    min_height::Int = 800
    pixels_per_sample::Int = 12
    orientation::Symbol = :auto

    # X-axis label sizing
    min_label_fontsize::Float64 = 6.0
    max_label_fontsize::Float64 = 12.0
    label_rotation::Float64 = 90.0

    # Taxa/legend
    top_n_taxa::Int = 25
    legend_fontsize::Float64 = 8.0
    legend_position::Symbol = :right

    # Margins (pixels)
    left_margin::Int = 80
    bottom_margin::Int = 120
    right_margin::Int = 150
    top_margin::Int = 40

    # Ordering
    sample_ordering::AxisOrdering = AxisOrdering()
    taxa_ordering::AxisOrdering = AxisOrdering(method=:sort, sort_by=:mean_abundance)

    # Dendrogram options
    show_sample_dendrogram::Bool = false
    show_taxa_dendrogram::Bool = false
    dendrogram_width::Int = 100
    color_branches_by_cluster::Bool = false
    n_clusters::Union{Int, Nothing} = nothing

    # Color scheme
    color_palette::Symbol = :maximally_distinguishable
    custom_colors::Union{Vector, Nothing} = nothing

    # View type for large datasets
    large_dataset_view::Symbol = :auto
    samples_per_page::Int = 150
    heatmap_threshold::Int = 300

    # Output
    output_formats::Vector{Symbol} = [:png, :svg]
    dpi::Int = 300
end

"""
    adaptive_label_fontsize(n_samples::Int; config::MicrobiomePlotConfig=MicrobiomePlotConfig())

Calculate adaptive font size for x-tick labels based on sample count.
Prevents label overlap while maintaining readability.

Uses linear interpolation between `config.max_label_fontsize` (for ≤50 samples)
and `config.min_label_fontsize` (for ≥500 samples).

# Arguments
- `n_samples::Int`: Number of samples to display
- `config::MicrobiomePlotConfig`: Configuration with font size bounds

# Returns
- `Float64`: Recommended font size in points

# Examples
```julia
adaptive_label_fontsize(30)   # Returns 12.0 (max)
adaptive_label_fontsize(275)  # Returns 9.0 (midpoint)
adaptive_label_fontsize(600)  # Returns 6.0 (min)
```
"""
function adaptive_label_fontsize(n_samples::Int; config::MicrobiomePlotConfig=MicrobiomePlotConfig())
    if n_samples <= 50
        return config.max_label_fontsize
    elseif n_samples >= 500
        return config.min_label_fontsize
    else
        # Linear interpolation
        ratio = (n_samples - 50) / (500 - 50)
        return config.max_label_fontsize - ratio * (config.max_label_fontsize - config.min_label_fontsize)
    end
end

"""
    calculate_figure_size(n_samples::Int, n_taxa::Int; config::MicrobiomePlotConfig=MicrobiomePlotConfig())

Calculate figure dimensions for abundance plots.
Automatically selects orientation based on sample count when `config.orientation == :auto`.

# Arguments
- `n_samples::Int`: Number of samples to display
- `n_taxa::Int`: Number of taxa (affects legend height)
- `config::MicrobiomePlotConfig`: Configuration with dimension parameters

# Returns
- `Tuple{Int, Int}`: (width, height) in pixels

# Sizing Strategy
| Samples | Orientation | Width | Height |
|---------|-------------|-------|--------|
| 1-50 | Landscape | 800-1200 | 600-800 |
| 50-150 | Square | 1200-1400 | 800-1000 |
| 150-300 | Portrait | 1400-1600 | 1200-2000 |
| 300-600 | Portrait | 1600 | 2000-3600 |
| 600+ | Portrait | 1600 | 4000 (max) |

# Examples
```julia
calculate_figure_size(30, 15)   # (900, 800) landscape
calculate_figure_size(200, 25)  # (1600, 2400) portrait
calculate_figure_size(500, 30)  # (1600, 4000) portrait
```
"""
function calculate_figure_size(n_samples::Int, n_taxa::Int; config::MicrobiomePlotConfig=MicrobiomePlotConfig())
    orientation = config.orientation
    if orientation == :auto
        orientation = n_samples > 100 ? :portrait : :landscape
    end

    if orientation == :portrait
        # Portrait: fixed width, variable height based on samples
        width = config.max_width
        sample_height = n_samples * config.pixels_per_sample
        taxa_legend_height = n_taxa * 18  # Approximate height per taxa in legend
        height = max(config.min_height, sample_height, taxa_legend_height)
        height = min(height, 4000)  # Cap at reasonable maximum
    else
        # Landscape: fixed height, variable width based on samples
        height = config.min_height
        width = max(800, n_samples * 8 + config.right_margin + config.left_margin)
        width = min(width, 3000)  # Cap at reasonable maximum
    end

    return (width, height)
end

"""
    compute_axis_ordering(data_matrix::Matrix, labels::Vector{String}; ordering::AxisOrdering=AxisOrdering())

Order samples or taxa using the specified method.

# Arguments
- `data_matrix::Matrix`: Abundance matrix (rows=items to order, cols=features)
- `labels::Vector{String}`: Labels for each row
- `ordering::AxisOrdering`: Ordering specification

# Returns
- `Tuple{Vector{Int}, Union{Clustering.Hclust, Nothing}}`: (ordered_indices, hclust_result_or_nothing)

# Ordering Methods
- `:hclust`: Hierarchical clustering with specified distance metric and linkage
- `:preordered`: Use provided label order from `ordering.preordered_labels`
- `:sort`: Sort by computed values (e.g., mean abundance)
- `:alphabetical`: Alphabetical sorting of labels

# Examples
```julia
# Hierarchical clustering (default)
indices, hclust = compute_axis_ordering(matrix, labels)

# Alphabetical
indices, _ = compute_axis_ordering(matrix, labels, ordering=AxisOrdering(method=:alphabetical))

# Pre-specified order
ordered = AxisOrdering(method=:preordered, preordered_labels=["B", "A", "C"])
indices, _ = compute_axis_ordering(matrix, labels, ordering=ordered)
```
"""
function compute_axis_ordering(data_matrix::Matrix, labels::Vector{String}; ordering::AxisOrdering=AxisOrdering())
    if ordering.method == :preordered && !isnothing(ordering.preordered_labels)
        # Use provided order
        label_to_idx = Dict(l => i for (i, l) in enumerate(labels))
        ordered_indices = [label_to_idx[l] for l in ordering.preordered_labels if haskey(label_to_idx, l)]
        return (ordered_indices, nothing)

    elseif ordering.method == :alphabetical
        perm = sortperm(labels)
        return (perm, nothing)

    elseif ordering.method == :sort && !isnothing(ordering.sort_by)
        # Sort by computed values (e.g., mean abundance)
        if ordering.sort_by isa Symbol
            # Default: mean across rows
            values = vec(Statistics.mean(data_matrix, dims=2))
        else
            # Custom function
            values = [ordering.sort_by(row) for row in eachrow(data_matrix)]
        end
        perm = sortperm(values, rev=true)
        return (perm, nothing)

    else  # :hclust (default)
        dist_func = if ordering.distance_metric == :braycurtis
            Distances.BrayCurtis()
        elseif ordering.distance_metric == :cosine
            Distances.CosineDist()
        else
            Distances.Euclidean()
        end

        # Compute pairwise distance matrix
        dist_matrix = Distances.pairwise(dist_func, data_matrix, dims=1)

        # Handle NaN/Inf values that can occur with zero vectors
        dist_matrix = replace(dist_matrix, NaN => 0.0, Inf => 1.0)

        # Convert to symmetric matrix and perform hierarchical clustering
        linkage_method = if ordering.linkage == :single
            :single
        elseif ordering.linkage == :complete
            :complete
        elseif ordering.linkage == :ward
            :ward
        else
            :average
        end

        hcl = Clustering.hclust(dist_matrix, linkage=linkage_method)
        return (hcl.order, hcl)
    end
end

"""
    determine_view_types(n_samples::Int; config::MicrobiomePlotConfig=MicrobiomePlotConfig())

Determine which visualization types to generate based on sample count.

# Arguments
- `n_samples::Int`: Number of samples
- `config::MicrobiomePlotConfig`: Configuration with view preferences

# Returns
- `Vector{Symbol}`: List of view types to generate (`:barplot`, `:heatmap`, `:paginated`)

# Auto-selection Logic
| Samples | Views Generated |
|---------|-----------------|
| ≤150 | [:barplot] |
| 150-300 | [:barplot, :heatmap] |
| >300 | [:heatmap, :paginated] |

# Examples
```julia
determine_view_types(50)   # [:barplot]
determine_view_types(200)  # [:barplot, :heatmap]
determine_view_types(500)  # [:heatmap, :paginated]
```
"""
function determine_view_types(n_samples::Int; config::MicrobiomePlotConfig=MicrobiomePlotConfig())
    if config.large_dataset_view == :auto
        if n_samples <= 150
            return [:barplot]
        elseif n_samples <= config.heatmap_threshold
            return [:barplot, :heatmap]
        else
            return [:heatmap, :paginated]
        end
    elseif config.large_dataset_view == :all
        return [:barplot, :heatmap, :paginated]
    else
        return [config.large_dataset_view]
    end
end

"""
    calculate_tick_step(n_samples::Int; max_labels::Int=100)

Calculate step size for showing every Nth label to prevent overlap.

# Arguments
- `n_samples::Int`: Total number of samples
- `max_labels::Int`: Maximum number of labels to show (default: 100)

# Returns
- `Int`: Step size (1 means show all labels)

# Examples
```julia
calculate_tick_step(50)   # 1 (show all)
calculate_tick_step(200)  # 2 (show every 2nd)
calculate_tick_step(500)  # 5 (show every 5th)
```
"""
function calculate_tick_step(n_samples::Int; max_labels::Int=100)
    if n_samples <= max_labels
        return 1
    else
        return ceil(Int, n_samples / max_labels)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract short sample identifier from full sample name using a regex pattern.

Useful for making x-axis labels more readable when sample names contain
barcode sequences, UUIDs, or other long suffixes.

# Arguments
- `sample_id::AbstractString`: Full sample identifier
- `pattern::Regex`: Pattern with a capture group for the short ID portion.
  The first capture group is returned as the short ID.
- `max_length::Int`: Maximum length before truncation with ellipsis (default: 16)

# Returns
- `String`: Short sample identifier (first capture group, or truncated original)

# Examples
```julia
# Extract prefix before barcode sequence
extract_short_sample_id("Sample-01-ACGTACGT-TGCATGCA"; pattern=r"(Sample-\\d+).*")
# => "Sample-01"

# Extract numeric portion
extract_short_sample_id("experiment_12345_replicate_A"; pattern=r"experiment_(\\d+)")
# => "12345"

# Fallback truncation when pattern doesn't match
extract_short_sample_id("VeryLongSampleNameWithoutPattern"; pattern=r"(NoMatch)")
# => "VeryLongSampleN…"
```
"""
function extract_short_sample_id(sample_id::AbstractString;
                                  pattern::Regex=r"(.+)",
                                  max_length::Int=16)
    m = match(pattern, sample_id)
    if m !== nothing && length(m.captures) >= 1
        return m.captures[1]
    end
    # Fallback: truncate with ellipsis
    return length(sample_id) > max_length ? sample_id[1:max_length] * "…" : sample_id
end

# Broadcast-friendly vector version
extract_short_sample_id(ids::AbstractVector; kwargs...) = [extract_short_sample_id(id; kwargs...) for id in ids]

"""
    _prepare_abundance_data(
        abundance_df::DataFrames.DataFrame;
        sample_col::Symbol,
        taxon_col::Symbol,
        abundance_col::Symbol,
        top_n::Int,
        config::MicrobiomePlotConfig
    )

Internal function to prepare abundance data for visualization.
Converts long-format DataFrame to matrix form and computes ordering.

# Returns
Named tuple with:
- `samples`: Vector of ordered sample names
- `taxa`: Vector of ordered taxa names (including "Other", "Unclassified")
- `matrix`: Abundance matrix (taxa × samples)
- `colors`: Vector of colors for each taxon
- `sample_hclust`: Hierarchical clustering result for samples (or nothing)
- `taxa_hclust`: Hierarchical clustering result for taxa (or nothing)
"""
function _prepare_abundance_data(
    abundance_df::DataFrames.DataFrame;
    sample_col::Symbol,
    taxon_col::Symbol,
    abundance_col::Symbol,
    top_n::Int,
    config::MicrobiomePlotConfig
)
    # Get unique samples and taxa
    all_samples = unique(abundance_df[!, sample_col])
    all_taxa = unique(abundance_df[!, taxon_col])

    # Calculate total abundance per taxon to find top N
    taxa_totals = DataFrames.combine(
        DataFrames.groupby(abundance_df, taxon_col),
        abundance_col => sum => :total
    )
    DataFrames.sort!(taxa_totals, :total, rev=true)

    # Identify top taxa, treating missing/unclassified specially
    is_special = x -> ismissing(x) || lowercase(string(x)) in ["unclassified", "missing", "unknown", "other"]
    regular_taxa = filter(row -> !is_special(row[taxon_col]), taxa_totals)
    top_taxa = regular_taxa[1:min(top_n, DataFrames.nrow(regular_taxa)), taxon_col]

    # Build sample × taxon matrix for ordering
    n_samples = length(all_samples)
    sample_to_idx = Dict(s => i for (i, s) in enumerate(all_samples))

    # Create wide-format matrix for all taxa (for clustering)
    n_all_taxa = length(all_taxa)
    taxa_to_idx = Dict(t => i for (i, t) in enumerate(all_taxa))
    full_matrix = zeros(n_samples, n_all_taxa)

    for row in eachrow(abundance_df)
        si = sample_to_idx[row[sample_col]]
        ti = taxa_to_idx[row[taxon_col]]
        full_matrix[si, ti] = row[abundance_col]
    end

    # Order samples
    sample_order, sample_hclust = compute_axis_ordering(
        full_matrix,
        string.(all_samples),
        ordering=config.sample_ordering
    )
    ordered_samples = all_samples[sample_order]

    # Create final matrix with top taxa + Other + Unclassified
    final_taxa = Vector{String}()
    for t in top_taxa
        push!(final_taxa, string(t))
    end

    # Calculate "Other" (sum of non-top regular taxa)
    other_taxa = setdiff(regular_taxa[!, taxon_col], top_taxa)
    has_other = !isempty(other_taxa)
    if has_other
        push!(final_taxa, "Other")
    end

    # Check for unclassified/missing
    special_taxa = filter(row -> is_special(row[taxon_col]), taxa_totals)
    has_unclassified = DataFrames.nrow(special_taxa) > 0
    if has_unclassified
        push!(final_taxa, "Unclassified")
    end

    # Build final abundance matrix (taxa × samples)
    n_final_taxa = length(final_taxa)
    n_ordered_samples = length(ordered_samples)
    abundance_matrix = zeros(n_final_taxa, n_ordered_samples)

    for (col_idx, sample) in enumerate(ordered_samples)
        sample_data = filter(row -> row[sample_col] == sample, abundance_df)

        for row in eachrow(sample_data)
            taxon = row[taxon_col]
            abundance = row[abundance_col]

            if is_special(taxon)
                # Add to Unclassified
                if has_unclassified
                    tidx = findfirst(==("Unclassified"), final_taxa)
                    abundance_matrix[tidx, col_idx] += abundance
                end
            elseif string(taxon) in final_taxa[1:length(top_taxa)]
                # Add to specific taxon
                tidx = findfirst(==(string(taxon)), final_taxa)
                abundance_matrix[tidx, col_idx] += abundance
            else
                # Add to Other
                if has_other
                    tidx = findfirst(==("Other"), final_taxa)
                    abundance_matrix[tidx, col_idx] += abundance
                end
            end
        end
    end

    # Generate colors
    n_colors = length(final_taxa)
    if !isnothing(config.custom_colors) && length(config.custom_colors) >= n_colors
        colors = config.custom_colors[1:n_colors]
    elseif config.color_palette == :maximally_distinguishable
        colors = Mycelia.n_maximally_distinguishable_colors(n_colors)
    else
        # Use ColorSchemes
        colors = [ColorSchemes.get(ColorSchemes.colorschemes[config.color_palette], i / n_colors) for i in 1:n_colors]
    end

    # Override colors for special categories
    if has_other
        other_idx = findfirst(==("Other"), final_taxa)
        colors[other_idx] = Colors.RGB(0.7, 0.7, 0.7)  # Light gray
    end
    if has_unclassified
        unclass_idx = findfirst(==("Unclassified"), final_taxa)
        colors[unclass_idx] = Colors.RGB(0.4, 0.4, 0.4)  # Dark gray
    end

    return (
        samples = string.(ordered_samples),
        taxa = final_taxa,
        matrix = abundance_matrix,
        colors = colors,
        sample_hclust = sample_hclust,
        taxa_hclust = nothing  # Could add taxa clustering if needed
    )
end

"""
    _create_barplot(
        data::NamedTuple;
        config::MicrobiomePlotConfig,
        title::String,
        rank::String
    )

Create a stacked barplot of microbiome abundances with adaptive sizing.

# Arguments
- `data::NamedTuple`: Prepared data from `_prepare_abundance_data`
- `config::MicrobiomePlotConfig`: Visualization configuration
- `title::String`: Plot title
- `rank::String`: Taxonomic rank for axis labels

# Returns
- `CairoMakie.Figure`: The generated figure
"""
function _create_barplot(
    data::NamedTuple;
    config::MicrobiomePlotConfig,
    title::String = "",
    rank::String = "taxon"
)
    samples = data.samples
    taxa = data.taxa
    matrix = data.matrix
    colors = data.colors

    n_samples = length(samples)
    n_taxa = length(taxa)

    # Calculate figure size
    width, height = calculate_figure_size(n_samples, n_taxa, config=config)

    # Determine orientation
    orientation = config.orientation
    if orientation == :auto
        orientation = n_samples > 100 ? :portrait : :landscape
    end

    # Calculate adaptive font size and tick step
    label_fontsize = adaptive_label_fontsize(n_samples, config=config)
    tick_step = calculate_tick_step(n_samples)

    # Create figure
    fig = CairoMakie.Figure(size=(width, height), fontsize=12)

    # Determine plot title
    plot_title = if isempty(title)
        "$(titlecase(rank)) Relative Abundance (top $(length(taxa) - count(t -> t in ["Other", "Unclassified"], taxa)))"
    else
        title
    end

    if orientation == :portrait
        # Portrait: samples on Y axis, abundance on X axis (horizontal bars)
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Relative Abundance",
            ylabel = "Sample",
            title = plot_title,
            yticks = (1:tick_step:n_samples, samples[1:tick_step:n_samples]),
            yticklabelsize = label_fontsize,
            yreversed = true  # First sample at top
        )

        # Plot horizontal stacked bars
        x_positions = 1:n_samples
        previous_widths = zeros(n_samples)

        for (i, taxon) in enumerate(taxa)
            widths = matrix[i, :]

            CairoMakie.barplot!(
                ax,
                x_positions,
                widths,
                offset = previous_widths,
                color = colors[i],
                direction = :x,  # Horizontal bars
                label = taxon
            )

            previous_widths .+= widths
        end

        CairoMakie.xlims!(ax, 0, 1.05)

    else
        # Landscape: samples on X axis (traditional vertical bars)
        ax = CairoMakie.Axis(
            fig[1, 1],
            xlabel = "Sample",
            ylabel = "Relative Abundance",
            title = plot_title,
            xticks = (1:tick_step:n_samples, samples[1:tick_step:n_samples]),
            xticklabelrotation = deg2rad(config.label_rotation),
            xticklabelsize = label_fontsize
        )

        # Plot vertical stacked bars
        x_positions = 1:n_samples
        previous_heights = zeros(n_samples)

        for (i, taxon) in enumerate(taxa)
            heights = matrix[i, :]

            CairoMakie.barplot!(
                ax,
                x_positions,
                heights,
                offset = previous_heights,
                color = colors[i],
                label = taxon
            )

            previous_heights .+= heights
        end

        CairoMakie.ylims!(ax, 0, 1.05)
    end

    # Add legend
    legend_nbanks = n_taxa > 15 ? 2 : 1
    fig[1, 2] = CairoMakie.Legend(
        fig,
        ax,
        titlecase(rank),
        framevisible = true,
        labelsize = config.legend_fontsize,
        titlesize = config.legend_fontsize + 2,
        patchsize = (12, 12),
        nbanks = legend_nbanks
    )

    return fig
end

"""
    _create_heatmap_with_dendrograms(
        data::NamedTuple;
        config::MicrobiomePlotConfig,
        title::String,
        rank::String
    )

Create a heatmap visualization with optional dendrograms.

# Arguments
- `data::NamedTuple`: Prepared data from `_prepare_abundance_data`
- `config::MicrobiomePlotConfig`: Visualization configuration
- `title::String`: Plot title
- `rank::String`: Taxonomic rank for axis labels

# Returns
- `CairoMakie.Figure`: The generated figure
"""
function _create_heatmap_with_dendrograms(
    data::NamedTuple;
    config::MicrobiomePlotConfig,
    title::String = "",
    rank::String = "taxon"
)
    samples = data.samples
    taxa = data.taxa
    matrix = data.matrix
    sample_hclust = data.sample_hclust

    n_samples = length(samples)
    n_taxa = length(taxa)

    # Calculate figure size (samples on X-axis, taxa on Y-axis)
    width = max(800, n_samples * 4 + 200)
    width = min(width, config.max_width, 3000)
    height = max(config.min_height, n_taxa * 15 + 200)
    height = min(height, 1400)

    # Calculate tick step for sample labels
    tick_step = calculate_tick_step(n_samples, max_labels=50)
    label_fontsize = adaptive_label_fontsize(n_samples, config=config)

    # Determine grid layout based on dendrogram settings
    show_sample_dendro = config.show_sample_dendrogram && !isnothing(sample_hclust)

    # Create figure
    fig = CairoMakie.Figure(size=(width, height), fontsize=10)

    # Determine column positions (dendrogram is now above, not to the left)
    hm_col = 1
    legend_col = 2

    # Plot title
    plot_title = if isempty(title)
        "$(titlecase(rank)) Abundance Heatmap (n=$n_samples samples)"
    else
        title
    end

    # Create heatmap axis (samples on X, taxa on Y)
    ax_hm = CairoMakie.Axis(
        fig[1, hm_col],
        xlabel = "Sample",
        ylabel = titlecase(rank),
        title = plot_title,
        xticks = (1:tick_step:n_samples, samples[1:tick_step:n_samples]),
        xticklabelrotation = deg2rad(45),
        xticklabelsize = label_fontsize,
        yticks = (1:n_taxa, taxa),
        yticklabelsize = 8,
        yreversed = true
    )

    # Create heatmap (samples on X, taxa on Y)
    # CairoMakie.heatmap!(ax, x, y, values) expects values[i,j] at (x[i], y[j])
    # So we need shape (n_samples × n_taxa), but matrix is (n_taxa × n_samples)
    hm = CairoMakie.heatmap!(
        ax_hm,
        1:n_samples,
        1:n_taxa,
        matrix',  # Transpose: (n_samples × n_taxa) to match axis order
        colormap = :viridis
    )

    # Add colorbar
    CairoMakie.Colorbar(
        fig[1, legend_col],
        hm,
        label = "Relative Abundance"
    )

    # Add sample dendrogram if requested (above heatmap, samples on X-axis)
    if show_sample_dendro
        ax_dendro = CairoMakie.Axis(
            fig[0, hm_col],  # Row 0 = above heatmap
            ylabel = "Height",
            xgridvisible = false,
            ygridvisible = false,
            xticklabelsvisible = false,
            xticksvisible = false,
            bottomspinevisible = false
        )

        # Plot dendrogram using existing function if available
        if isdefined(Mycelia, :plot_dendrogram_from_hclust!)
            Mycelia.plot_dendrogram_from_hclust!(
                ax_dendro,
                sample_hclust,
                orientation = :vertical
            )
        end

        CairoMakie.yreverse!(ax_dendro)
        CairoMakie.linkxaxes!(ax_hm, ax_dendro)
        CairoMakie.rowsize!(fig.layout, 0, CairoMakie.Relative(0.15))
    end

    return fig
end

"""
    _create_paginated_barplots(
        data::NamedTuple;
        config::MicrobiomePlotConfig,
        title::String,
        rank::String
    )

Create paginated barplots for large datasets (600+ samples).

# Arguments
- `data::NamedTuple`: Prepared data from `_prepare_abundance_data`
- `config::MicrobiomePlotConfig`: Visualization configuration
- `title::String`: Plot title
- `rank::String`: Taxonomic rank for axis labels

# Returns
- `Vector{CairoMakie.Figure}`: Vector of figures, one per page
"""
function _create_paginated_barplots(
    data::NamedTuple;
    config::MicrobiomePlotConfig,
    title::String = "",
    rank::String = "taxon"
)
    samples = data.samples
    taxa = data.taxa
    matrix = data.matrix
    colors = data.colors

    n_samples = length(samples)
    n_taxa = length(taxa)
    samples_per_page = config.samples_per_page

    # Calculate number of pages
    n_pages = ceil(Int, n_samples / samples_per_page)

    figures = CairoMakie.Figure[]

    for page in 1:n_pages
        # Calculate sample range for this page
        start_idx = (page - 1) * samples_per_page + 1
        end_idx = min(page * samples_per_page, n_samples)
        page_samples = samples[start_idx:end_idx]
        page_matrix = matrix[:, start_idx:end_idx]

        n_page_samples = length(page_samples)

        # Create page-specific data
        page_data = (
            samples = page_samples,
            taxa = taxa,
            matrix = page_matrix,
            colors = colors,
            sample_hclust = nothing,
            taxa_hclust = nothing
        )

        # Create page title
        page_title = if isempty(title)
            "$(titlecase(rank)) Abundance (Page $page/$n_pages, samples $start_idx-$end_idx)"
        else
            "$title (Page $page/$n_pages)"
        end

        # Use portrait orientation for paginated views
        page_config = MicrobiomePlotConfig(
            max_width = config.max_width,
            min_height = config.min_height,
            pixels_per_sample = config.pixels_per_sample,
            orientation = :portrait,
            min_label_fontsize = config.min_label_fontsize,
            max_label_fontsize = config.max_label_fontsize,
            label_rotation = config.label_rotation,
            top_n_taxa = config.top_n_taxa,
            legend_fontsize = config.legend_fontsize,
            legend_position = config.legend_position,
            color_palette = config.color_palette,
            custom_colors = config.custom_colors
        )

        fig = _create_barplot(
            page_data,
            config = page_config,
            title = page_title,
            rank = rank
        )

        push!(figures, fig)
    end

    return figures
end

"""
    save_plot(
        fig,
        base_path::String;
        formats::Vector{Symbol} = [:png, :svg],
        dpi::Int = 300
    )

Save a plot in multiple formats. Works with both CairoMakie figures and
StatsPlots/Plots.jl plots.

# Arguments
- `fig`: Figure object (CairoMakie.Figure or Plots.Plot)
- `base_path::String`: Base path without extension (e.g., "output/genus_barplot")
- `formats::Vector{Symbol}`: Output formats (default: [:png, :svg])
- `dpi::Int`: Resolution in dots per inch (default: 300)

# Returns
- `Dict{Symbol, String}`: Mapping of format to saved file path

# Examples
```julia
# CairoMakie figure
paths = save_plot(makie_fig, "results/genus_abundance")

# StatsPlots/Plots.jl plot
paths = save_plot(statsplots_p, "results/barplot")

# Returns: Dict(:png => "results/genus_abundance.png", :svg => "results/genus_abundance.svg")
```
"""
function save_plot(
    fig,
    base_path::String;
    formats::Vector{Symbol} = [:png, :svg],
    dpi::Int = 300
)
    paths = Dict{Symbol, String}()

    # Create directory if it doesn't exist
    dir = dirname(base_path)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end

    # Detect figure type - check for CairoMakie.Figure directly or via .figure property
    # (FigureAxisPlot has a .figure field that contains a CairoMakie.Figure)
    is_makie = fig isa CairoMakie.Figure ||
               (hasproperty(fig, :figure) && fig.figure isa CairoMakie.Figure)
    is_plots = isdefined(Main, :Plots) && fig isa Main.Plots.Plot ||
               isdefined(Mycelia, :StatsPlots) && (
                   typeof(fig) <: StatsPlots.Plots.Plot ||
                   string(typeof(fig)) |> t -> occursin("Plot{", t)
               )

    for fmt in formats
        path = "$(base_path).$(fmt)"

        try
            if is_makie
                # CairoMakie figure
                if fmt == :png
                    CairoMakie.save(path, fig, px_per_unit=dpi/72)
                elseif fmt in [:svg, :pdf]
                    CairoMakie.save(path, fig)
                else
                    @warn "Unknown format for CairoMakie: $fmt, skipping"
                    continue
                end
            else
                # Assume StatsPlots/Plots.jl
                StatsPlots.savefig(fig, path)
            end

            paths[fmt] = path
            @info "Saved: $path"
        catch e
            @warn "Failed to save $path: $e"
        end
    end

    return paths
end

# Alias for backward compatibility
const save_microbiome_plot = save_plot

"""
    plot_microbiome_abundance(
        abundance_df::DataFrames.DataFrame;
        sample_col::Symbol = :sample,
        taxon_col::Symbol = :taxon,
        abundance_col::Symbol = :relative_abundance,
        rank::String = "genus",
        config::MicrobiomePlotConfig = MicrobiomePlotConfig(),
        title::String = "",
        output_dir::Union{String, Nothing} = nothing
    )

Create publication-quality microbiome abundance visualizations.
Automatically selects the best view type(s) based on sample count.

# Arguments
- `abundance_df::DataFrames.DataFrame`: Long-format abundance data
- `sample_col::Symbol`: Column containing sample identifiers
- `taxon_col::Symbol`: Column containing taxon names
- `abundance_col::Symbol`: Column containing abundance values (0-1 scale)
- `rank::String`: Taxonomic rank name for labels
- `config::MicrobiomePlotConfig`: Visualization configuration
- `title::String`: Optional plot title
- `output_dir::Union{String, Nothing}`: Directory to save plots (optional)

# Returns
- `Dict{Symbol, Any}`: Dictionary with visualization results:
  - `:barplot` => CairoMakie.Figure (if generated)
  - `:heatmap` => CairoMakie.Figure (if generated)
  - `:paginated` => Vector{CairoMakie.Figure} (if generated)
  - `:data` => Prepared data NamedTuple

# View Selection (when `config.large_dataset_view == :auto`)
| Sample Count | Views Generated |
|--------------|-----------------|
| ≤150 | barplot only |
| 150-300 | barplot + heatmap |
| >300 | heatmap + paginated |

# Examples
```julia
# Basic usage
results = plot_microbiome_abundance(df, sample_col=:sample_id, taxon_col=:genus)

# With custom configuration for many samples
config = MicrobiomePlotConfig(
    orientation = :portrait,
    show_sample_dendrogram = true,
    top_n_taxa = 20
)
results = plot_microbiome_abundance(df, config=config, output_dir="figures/")

# Access specific views
barplot_fig = results[:barplot]
heatmap_fig = results[:heatmap]
```
"""
function plot_microbiome_abundance(
    abundance_df::DataFrames.DataFrame;
    sample_col::Symbol = :sample,
    taxon_col::Symbol = :taxon,
    abundance_col::Symbol = :relative_abundance,
    rank::String = "genus",
    config::MicrobiomePlotConfig = MicrobiomePlotConfig(),
    title::String = "",
    output_dir::Union{String, Nothing} = nothing
)
    # Validate input columns exist
    for col in [sample_col, taxon_col, abundance_col]
        if !DataFrames.hasproperty(abundance_df, col)
            error("DataFrame must have column: $col")
        end
    end

    n_samples = length(unique(abundance_df[!, sample_col]))
    @info "Processing $n_samples samples for microbiome visualization"

    # Prepare data
    data = _prepare_abundance_data(
        abundance_df,
        sample_col = sample_col,
        taxon_col = taxon_col,
        abundance_col = abundance_col,
        top_n = config.top_n_taxa,
        config = config
    )

    # Determine which views to generate
    view_types = determine_view_types(n_samples, config=config)
    @info "Generating views: $(join(string.(view_types), ", "))"

    results = Dict{Symbol, Any}()
    results[:data] = data

    # Generate each view type
    for view_type in view_types
        if view_type == :barplot
            fig = _create_barplot(
                data,
                config = config,
                title = title,
                rank = rank
            )
            results[:barplot] = fig

            if !isnothing(output_dir)
                save_microbiome_plot(
                    fig,
                    joinpath(output_dir, "$(rank)_barplot"),
                    formats = config.output_formats,
                    dpi = config.dpi
                )
            end

        elseif view_type == :heatmap
            fig = _create_heatmap_with_dendrograms(
                data,
                config = config,
                title = title,
                rank = rank
            )
            results[:heatmap] = fig

            if !isnothing(output_dir)
                save_microbiome_plot(
                    fig,
                    joinpath(output_dir, "$(rank)_heatmap"),
                    formats = config.output_formats,
                    dpi = config.dpi
                )
            end

        elseif view_type == :paginated
            figures = _create_paginated_barplots(
                data,
                config = config,
                title = title,
                rank = rank
            )
            results[:paginated] = figures

            if !isnothing(output_dir)
                for (i, fig) in enumerate(figures)
                    save_microbiome_plot(
                        fig,
                        joinpath(output_dir, "$(rank)_barplot_page$(lpad(i, 2, '0'))"),
                        formats = config.output_formats,
                        dpi = config.dpi
                    )
                end
            end
        end
    end

    return results
end

# ============================================================================
# Coverage Integration Functions for Contig-Level Abundance
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Merge mosdepth coverage data with taxonomy assignments to create a unified
contig-level abundance table.

# Arguments
- `coverage_df::DataFrames.DataFrame`: DataFrame from `parse_mosdepth_summary()` with columns:
  - `chrom`: Contig/chromosome ID
  - `mean`: Mean coverage depth
  - Optional: `length`, `bases`, `min`, `max`
- `taxonomy_df::DataFrames.DataFrame`: DataFrame with taxonomy assignments, must have:
  - A contig ID column (specified by `contig_col`)
  - Taxonomy rank columns (species, genus, family, etc.)
- `contig_col::Symbol`: Column name in taxonomy_df containing contig IDs (default: :contig_id)
- `coverage_col::Symbol`: Column name in coverage_df containing coverage values (default: :mean)
- `coverage_contig_col::Symbol`: Column name in coverage_df containing contig IDs (default: :chrom)
- `min_coverage::Float64`: Minimum coverage threshold to include contigs (default: 0.0)
- `min_length::Int`: Minimum contig length to include (default: 0, requires `length` column)

# Returns
- `DataFrames.DataFrame`: Merged table with coverage and taxonomy information

# Example
```julia
coverage_df = Mycelia.parse_mosdepth_summary("sample.mosdepth.summary.txt")
taxonomy_df = Mycelia.parse_blast_report("sample.blast.txt") |> Mycelia.ensure_lineage_columns

merged = Mycelia.merge_coverage_with_taxonomy(
    coverage_df,
    taxonomy_df,
    contig_col = :query_id,  # BLAST query ID column
    min_coverage = 1.0
)
```
"""
function merge_coverage_with_taxonomy(
    coverage_df::DataFrames.DataFrame,
    taxonomy_df::DataFrames.DataFrame;
    contig_col::Symbol = :contig_id,
    coverage_col::Symbol = :mean,
    coverage_contig_col::Symbol = :chrom,
    min_coverage::Float64 = 0.0,
    min_length::Int = 0
)
    # Validate required columns in coverage_df
    if !DataFrames.hasproperty(coverage_df, coverage_contig_col)
        error("coverage_df must have column: $coverage_contig_col")
    end
    if !DataFrames.hasproperty(coverage_df, coverage_col)
        error("coverage_df must have column: $coverage_col")
    end

    # Validate required columns in taxonomy_df
    if !DataFrames.hasproperty(taxonomy_df, contig_col)
        error("taxonomy_df must have column: $contig_col")
    end

    # Apply coverage filter
    filtered_coverage = coverage_df[coverage_df[!, coverage_col] .>= min_coverage, :]

    # Apply length filter if column exists and min_length > 0
    if min_length > 0 && DataFrames.hasproperty(filtered_coverage, :length)
        filtered_coverage = filtered_coverage[filtered_coverage[!, :length] .>= min_length, :]
    end

    @info "Coverage filtering: $(DataFrames.nrow(coverage_df)) -> $(DataFrames.nrow(filtered_coverage)) contigs (min_coverage=$min_coverage, min_length=$min_length)"

    # Merge coverage with taxonomy
    merged = DataFrames.leftjoin(
        filtered_coverage,
        taxonomy_df,
        on = coverage_contig_col => contig_col,
        matchmissing = :notequal
    )

    # Count successful joins
    n_with_taxonomy = sum(!ismissing, merged[!, :domain])  # Assuming domain is a standard column
    if DataFrames.hasproperty(merged, :domain)
        n_with_taxonomy = sum(!ismissing, merged[!, :domain])
    else
        # Try to count based on any taxonomy column present
        taxonomy_cols = [:species, :genus, :family, :order, :class, :phylum, :kingdom, :domain]
        for col in taxonomy_cols
            if DataFrames.hasproperty(merged, col)
                n_with_taxonomy = sum(!ismissing, merged[!, col])
                break
            end
        end
    end

    @info "Taxonomy join: $(DataFrames.nrow(merged)) contigs, $n_with_taxonomy with taxonomy assignments"

    return merged
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute coverage-weighted relative abundances at a specified taxonomic rank.

Takes merged coverage + taxonomy data and aggregates by taxonomic rank,
weighting abundances by coverage depth. The output is formatted for use
with `plot_microbiome_abundance()`.

# Arguments
- `merged_df::DataFrames.DataFrame`: Output from `merge_coverage_with_taxonomy()`
- `sample_id::String`: Sample identifier to include in output
- `rank::Symbol`: Taxonomic rank to aggregate by (default: :genus)
  Options: :species, :genus, :family, :order, :class, :phylum, :kingdom, :domain
- `coverage_col::Symbol`: Column containing coverage values (default: :mean)
- `include_unclassified::Bool`: Whether to include contigs without taxonomy (default: true)
- `unclassified_label::String`: Label for unclassified contigs (default: "Unclassified")

# Returns
- `DataFrames.DataFrame`: Long-format abundance table with columns:
  - `sample`: Sample identifier
  - `taxon`: Taxon name at specified rank
  - `relative_abundance`: Coverage-weighted relative abundance (0-1)
  - `total_coverage`: Sum of coverage for this taxon
  - `n_contigs`: Number of contigs assigned to this taxon

# Example
```julia
merged = Mycelia.merge_coverage_with_taxonomy(coverage_df, taxonomy_df)
abundance = Mycelia.compute_coverage_weighted_abundance(
    merged,
    sample_id = "Sample_001",
    rank = :genus
)

# Use with visualization
results = Mycelia.plot_microbiome_abundance(
    abundance,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance
)
```
"""
function compute_coverage_weighted_abundance(
    merged_df::DataFrames.DataFrame,
    sample_id::String;
    rank::Symbol = :genus,
    coverage_col::Symbol = :mean,
    include_unclassified::Bool = true,
    unclassified_label::String = "Unclassified"
)
    # Validate rank column exists
    if !DataFrames.hasproperty(merged_df, rank)
        error("merged_df must have rank column: $rank")
    end
    if !DataFrames.hasproperty(merged_df, coverage_col)
        error("merged_df must have coverage column: $coverage_col")
    end

    # Create working copy with taxon column
    work_df = DataFrames.select(merged_df, coverage_col, rank)
    DataFrames.rename!(work_df, rank => :taxon)

    # Handle missing values
    if include_unclassified
        work_df.taxon = coalesce.(work_df.taxon, unclassified_label)
    else
        work_df = work_df[.!ismissing.(work_df.taxon), :]
    end

    # Convert to string type for consistency
    work_df.taxon = string.(work_df.taxon)

    # Aggregate by taxon
    aggregated = DataFrames.combine(
        DataFrames.groupby(work_df, :taxon),
        coverage_col => sum => :total_coverage,
        DataFrames.nrow => :n_contigs
    )

    # Compute relative abundance
    total_coverage = sum(aggregated.total_coverage)
    if total_coverage > 0
        aggregated.relative_abundance = aggregated.total_coverage ./ total_coverage
    else
        aggregated.relative_abundance = zeros(DataFrames.nrow(aggregated))
    end

    # Add sample column
    aggregated.sample = fill(sample_id, DataFrames.nrow(aggregated))

    # Reorder columns
    result = DataFrames.select(aggregated, :sample, :taxon, :relative_abundance, :total_coverage, :n_contigs)

    # Sort by abundance descending
    DataFrames.sort!(result, :relative_abundance, rev=true)

    @info "Aggregated to $(DataFrames.nrow(result)) taxa at $rank level for $sample_id"

    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Process multiple samples to create a combined coverage-weighted abundance table.

This is a convenience function that processes multiple mosdepth + taxonomy file
pairs and combines them into a single abundance table suitable for multi-sample
visualization.

# Arguments
- `sample_data::Vector{NamedTuple}`: Vector of named tuples, each containing:
  - `sample_id::String`: Sample identifier
  - `coverage_file::String`: Path to mosdepth summary file
  - `taxonomy_file::String`: Path to taxonomy file (BLAST, mmseqs, or Arrow)
  - Optional: `contig_col::Symbol`: Column name for contig IDs in taxonomy file
- `rank::Symbol`: Taxonomic rank to aggregate by (default: :genus)
- `min_coverage::Float64`: Minimum coverage threshold (default: 0.0)
- `min_length::Int`: Minimum contig length (default: 0)
- `include_unclassified::Bool`: Include unclassified contigs (default: true)
- `taxonomy_format::Symbol`: Format of taxonomy files (default: :auto)
  Options: :auto, :blast, :mmseqs, :arrow

# Returns
- `DataFrames.DataFrame`: Combined abundance table for all samples

# Example
```julia
samples = [
    (sample_id = "Sample_001",
     coverage_file = "sample1.mosdepth.summary.txt",
     taxonomy_file = "sample1.blast.txt"),
    (sample_id = "Sample_002",
     coverage_file = "sample2.mosdepth.summary.txt",
     taxonomy_file = "sample2.blast.txt"),
]

combined = Mycelia.process_samples_coverage_abundance(
    samples,
    rank = :family,
    min_coverage = 1.0
)

results = Mycelia.plot_microbiome_abundance(combined)
```
"""
function process_samples_coverage_abundance(
    sample_data::Vector{<:NamedTuple};
    rank::Symbol = :genus,
    min_coverage::Float64 = 0.0,
    min_length::Int = 0,
    include_unclassified::Bool = true,
    taxonomy_format::Symbol = :auto
)
    all_abundances = DataFrames.DataFrame()

    for (i, sample) in enumerate(sample_data)
        @info "Processing sample $(i)/$(length(sample_data)): $(sample.sample_id)"

        # Load coverage data
        coverage_df = parse_mosdepth_summary(sample.coverage_file)

        # Load taxonomy data based on format
        taxonomy_df = _load_taxonomy_file(
            sample.taxonomy_file,
            format = taxonomy_format
        )

        # Get contig column name (default or from sample tuple)
        contig_col = haskey(sample, :contig_col) ? sample.contig_col : _detect_contig_column(taxonomy_df)

        # Merge coverage with taxonomy
        merged = merge_coverage_with_taxonomy(
            coverage_df,
            taxonomy_df,
            contig_col = contig_col,
            min_coverage = min_coverage,
            min_length = min_length
        )

        # Compute abundance for this sample
        abundance = compute_coverage_weighted_abundance(
            merged,
            sample.sample_id,
            rank = rank,
            include_unclassified = include_unclassified
        )

        # Append to combined table
        all_abundances = DataFrames.vcat(all_abundances, abundance, cols = :union)
    end

    @info "Combined abundance table: $(DataFrames.nrow(all_abundances)) rows, $(length(unique(all_abundances.sample))) samples"

    return all_abundances
end

"""
Internal helper to load taxonomy files in various formats.
"""
function _load_taxonomy_file(
    filepath::String;
    format::Symbol = :auto
)
    # Auto-detect format
    if format == :auto
        if endswith(filepath, ".arrow")
            format = :arrow
        elseif occursin("_lca.tsv", filepath) || occursin("mmseqs", lowercase(filepath))
            format = :mmseqs
        else
            format = :blast
        end
    end

    if format == :arrow
        return DataFrames.DataFrame(Arrow.Table(filepath))
    elseif format == :mmseqs
        df = parse_mmseqs_easy_taxonomy_lca_tsv(filepath)
        # Ensure lineage columns exist
        if DataFrames.hasproperty(df, :taxon_id) && !DataFrames.hasproperty(df, :domain)
            taxids = collect(skipmissing(df.taxon_id))
            if !isempty(taxids)
                lineage = taxids2taxonkit_summarized_lineage_table(taxids)
                df = DataFrames.leftjoin(df, lineage, on = :taxon_id => :taxid)
            end
        end
        return df
    elseif format == :blast
        df = parse_blast_report(filepath)
        return ensure_lineage_columns(df)
    else
        error("Unknown taxonomy format: $format")
    end
end

"""
Internal helper to detect the contig ID column in taxonomy DataFrames.
"""
function _detect_contig_column(df::DataFrames.DataFrame)
    # Common contig column names in order of preference
    candidates = [:contig_id, :query_id, Symbol("query id"), :qseqid, :contig, :sequence_id]

    for col in candidates
        if DataFrames.hasproperty(df, col)
            return col
        end
    end

    # If none found, return first column as fallback
    @warn "Could not detect contig column, using first column: $(names(df)[1])"
    return Symbol(names(df)[1])
end

# ============================================================================
# BLAST-specific Coverage Integration
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute coverage-weighted abundance from BLAST results and mosdepth coverage.

This is a convenience function for the common workflow of combining BLAST taxonomy
assignments with mosdepth contig coverage to produce relative abundances.

# Arguments
- `blast_file::String`: Path to BLAST output file (outfmt 6 or 7)
- `mosdepth_summary::String`: Path to mosdepth summary file (.mosdepth.summary.txt)
- `sample_id::String`: Sample identifier for the output DataFrame
- `rank::Symbol`: Taxonomic rank to aggregate by (default: :genus)
- `min_coverage::Float64`: Minimum mean coverage to include contig (default: 1.0)
- `min_length::Int`: Minimum contig length to include (default: 0)
- `evalue_max::Float64`: Maximum e-value for BLAST hits (default: 1e-10)
- `include_unclassified::Bool`: Include contigs without BLAST hits (default: true)

# Returns
- `DataFrames.DataFrame`: Abundance table with columns:
  - `sample`: Sample identifier
  - `taxon`: Taxon name at specified rank
  - `relative_abundance`: Coverage-weighted relative abundance (0-1)
  - `total_coverage`: Sum of coverage for this taxon
  - `n_contigs`: Number of contigs assigned to this taxon

# Example
```julia
abundance = Mycelia.blast_coverage_abundance(
    "contigs.blast.txt",
    "contigs.mosdepth.summary.txt",
    "Sample_001",
    rank = :family,
    min_coverage = 3.0
)

# Combine multiple samples
all_abundances = DataFrames.vcat([
    Mycelia.blast_coverage_abundance(blast_files[i], mosdepth_files[i], sample_ids[i])
    for i in 1:length(sample_ids)
]...)

# Visualize
results = Mycelia.plot_microbiome_abundance(all_abundances, rank = "family")
```
"""
function blast_coverage_abundance(
    blast_file::String,
    mosdepth_summary::String,
    sample_id::String;
    rank::Symbol = :genus,
    min_coverage::Float64 = 1.0,
    min_length::Int = 0,
    evalue_max::Float64 = 1e-10,
    include_unclassified::Bool = true
)
    @info "Processing BLAST coverage abundance for $sample_id"

    # Load and parse BLAST results
    blast_df = parse_blast_report(blast_file)
    @info "Loaded $(DataFrames.nrow(blast_df)) BLAST hits"

    # Filter by e-value if column exists
    if DataFrames.hasproperty(blast_df, :evalue) || DataFrames.hasproperty(blast_df, Symbol("evalue"))
        evalue_col = DataFrames.hasproperty(blast_df, :evalue) ? :evalue : Symbol("evalue")
        blast_df = blast_df[blast_df[!, evalue_col] .<= evalue_max, :]
        @info "After e-value filter (≤$evalue_max): $(DataFrames.nrow(blast_df)) hits"
    end

    # Get best hit per query (by bit score or e-value)
    query_col = _detect_contig_column(blast_df)
    if DataFrames.hasproperty(blast_df, Symbol("bit score"))
        # Group by query and take best hit by bit score
        blast_df = DataFrames.combine(
            DataFrames.groupby(blast_df, query_col),
            sdf -> sdf[argmax(sdf[!, Symbol("bit score")]), :]
        )
        @info "Best hits per query: $(DataFrames.nrow(blast_df)) contigs"
    end

    # Add lineage columns
    taxonomy_df = ensure_lineage_columns(blast_df)

    # Load mosdepth coverage
    coverage_df = parse_mosdepth_summary(mosdepth_summary)

    # Merge coverage with taxonomy
    merged = merge_coverage_with_taxonomy(
        coverage_df,
        taxonomy_df,
        contig_col = query_col,
        min_coverage = min_coverage,
        min_length = min_length
    )

    # Compute coverage-weighted abundance
    abundance = compute_coverage_weighted_abundance(
        merged,
        sample_id,
        rank = rank,
        include_unclassified = include_unclassified
    )

    return abundance
end

# ============================================================================
# MMseqs2-specific Coverage Integration
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute coverage-weighted abundance from mmseqs2 easy-taxonomy results and mosdepth coverage.

This is a convenience function for combining mmseqs2 taxonomic assignments
(from UniRef databases) with mosdepth contig coverage.

# Arguments
- `mmseqs_lca_file::String`: Path to mmseqs2 LCA TSV file (from easy-taxonomy workflow)
- `mosdepth_summary::String`: Path to mosdepth summary file (.mosdepth.summary.txt)
- `sample_id::String`: Sample identifier for the output DataFrame
- `rank::Symbol`: Taxonomic rank to aggregate by (default: :genus)
- `min_coverage::Float64`: Minimum mean coverage to include contig (default: 1.0)
- `min_length::Int`: Minimum contig length to include (default: 0)
- `min_support::Float64`: Minimum -log(E-value) support for mmseqs assignments (default: 0.0)
- `include_unclassified::Bool`: Include contigs without taxonomy (default: true)

# Returns
- `DataFrames.DataFrame`: Abundance table with columns:
  - `sample`: Sample identifier
  - `taxon`: Taxon name at specified rank
  - `relative_abundance`: Coverage-weighted relative abundance (0-1)
  - `total_coverage`: Sum of coverage for this taxon
  - `n_contigs`: Number of contigs assigned to this taxon

# Example
```julia
abundance = Mycelia.mmseqs_coverage_abundance(
    "contigs_lca.tsv",
    "contigs.mosdepth.summary.txt",
    "Sample_001",
    rank = :family
)

# For UniRef50 vs UniRef90 comparison
uniref50_abundance = Mycelia.mmseqs_coverage_abundance(
    "contigs.uniref50_lca.tsv",
    "contigs.mosdepth.summary.txt",
    "Sample_001_UniRef50"
)

uniref90_abundance = Mycelia.mmseqs_coverage_abundance(
    "contigs.uniref90_lca.tsv",
    "contigs.mosdepth.summary.txt",
    "Sample_001_UniRef90"
)
```
"""
function mmseqs_coverage_abundance(
    mmseqs_lca_file::String,
    mosdepth_summary::String,
    sample_id::String;
    rank::Symbol = :genus,
    min_coverage::Float64 = 1.0,
    min_length::Int = 0,
    min_support::Float64 = 0.0,
    include_unclassified::Bool = true
)
    @info "Processing mmseqs2 coverage abundance for $sample_id"

    # Load and parse mmseqs LCA results
    mmseqs_df = parse_mmseqs_easy_taxonomy_lca_tsv(mmseqs_lca_file)
    @info "Loaded $(DataFrames.nrow(mmseqs_df)) mmseqs2 taxonomy assignments"

    # Filter by support if specified
    support_col = Symbol("support -log(E-value)")
    if min_support > 0 && DataFrames.hasproperty(mmseqs_df, support_col)
        mmseqs_df = mmseqs_df[mmseqs_df[!, support_col] .>= min_support, :]
        @info "After support filter (≥$min_support): $(DataFrames.nrow(mmseqs_df)) assignments"
    end

    # Add full lineage columns from taxon_id
    if DataFrames.hasproperty(mmseqs_df, :taxon_id) && !DataFrames.hasproperty(mmseqs_df, :domain)
        taxids = collect(skipmissing(mmseqs_df.taxon_id))
        if !isempty(taxids)
            @info "Looking up lineage for $(length(taxids)) unique taxids"
            lineage = taxids2taxonkit_summarized_lineage_table(taxids)
            mmseqs_df = DataFrames.leftjoin(mmseqs_df, lineage, on = :taxon_id => :taxid)
        end
    end

    # Load mosdepth coverage
    coverage_df = parse_mosdepth_summary(mosdepth_summary)

    # Merge coverage with taxonomy
    merged = merge_coverage_with_taxonomy(
        coverage_df,
        mmseqs_df,
        contig_col = :contig_id,
        min_coverage = min_coverage,
        min_length = min_length
    )

    # Compute coverage-weighted abundance
    abundance = compute_coverage_weighted_abundance(
        merged,
        sample_id,
        rank = rank,
        include_unclassified = include_unclassified
    )

    return abundance
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create violin plots for alpha diversity metrics comparison across methods or groups.

Generates a 2x2 grid showing Shannon index, Simpson's index, species richness,
and Pielou's evenness as violin plots with overlaid scatter points.

# Arguments
- `alpha_df::DataFrames.DataFrame`: DataFrame with alpha diversity metrics.
  Required columns: `sample`, `shannon`, `simpsons`, `richness`, `evenness`.
  Optional column: `method` or `group` for grouping (uses first found, or "All" if neither exists).
- `figsize::Tuple{Int,Int}`: Figure dimensions in pixels (default: (1000, 800))
- `color`: Color for violins and points (default: :teal)
- `group_col::Symbol`: Column to use for grouping (default: auto-detect :method or :group)

# Returns
- `CairoMakie.Figure`: Figure with 2x2 grid of violin plots

# Examples
```julia
# Single method/group
alpha_df = calculate_alpha_diversity(abundance_matrix, samples)
fig = plot_alpha_diversity_violins(alpha_df)

# Multiple methods for comparison
combined_df = vcat(
    insertcols!(alpha_df1, :method => "Method A"),
    insertcols!(alpha_df2, :method => "Method B")
)
fig = plot_alpha_diversity_violins(combined_df)
```
"""
function plot_alpha_diversity_violins(alpha_df::DataFrames.DataFrame;
                                       figsize::Tuple{Int,Int}=(1000, 800),
                                       color=:teal,
                                       group_col::Union{Symbol,Nothing}=nothing)
    fig = CairoMakie.Figure(size=figsize)

    metrics = [(:shannon, "Shannon Index"),
               (:simpsons, "Simpson's Index"),
               (:richness, "Species Richness"),
               (:evenness, "Pielou's Evenness")]

    # Auto-detect grouping column
    if group_col === nothing
        if DataFrames.hasproperty(alpha_df, :method)
            group_col = :method
        elseif DataFrames.hasproperty(alpha_df, :group)
            group_col = :group
        else
            # No grouping - add a dummy column
            alpha_df = DataFrames.copy(alpha_df)
            alpha_df[!, :_group] = fill("All", DataFrames.nrow(alpha_df))
            group_col = :_group
        end
    end

    groups = sort(unique(alpha_df[!, group_col]))
    n_groups = length(groups)

    for (i, (metric, title)) in enumerate(metrics)
        row = div(i - 1, 2) + 1
        col = mod(i - 1, 2) + 1
        ax = CairoMakie.Axis(fig[row, col],
            title = title,
            xlabel = n_groups > 1 ? string(titlecase(string(group_col))) : "",
            ylabel = title
        )

        for (j, grp) in enumerate(groups)
            group_data = DataFrames.filter(r -> r[group_col] == grp, alpha_df)
            values = group_data[!, metric]

            # Filter out missing/NaN values
            valid_values = collect(skipmissing(filter(!isnan, values)))

            if length(valid_values) > 1
                CairoMakie.violin!(ax, fill(j, length(valid_values)), valid_values,
                    color = (color, 0.5), strokewidth = 1)
            end
            # Always show individual points
            if !isempty(valid_values)
                CairoMakie.scatter!(ax, fill(j, length(valid_values)), valid_values,
                    color = color, markersize = 4, alpha = 0.6)
            end
        end

        ax.xticks = (1:n_groups, string.(groups))
    end

    return fig
end
