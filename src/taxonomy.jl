import DataFrames
import CairoMakie
import StatsBase
import OrderedCollections
import Colors
import Random

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






# function blastdb_accessions_to_taxid(;blastdb, outfile = blastdb * "." * Mycelia.normalized_current_date() * ".accession_to_taxid.arrow")
function blastdb_accessions_to_taxid(;blastdb, outfile = blastdb * ".accession_to_taxid.arrow")
    if !isfile(outfile)
        # Processing BLAST DB: 100%|██████████████████████████████| Time: 1:01:09
        #   completed:  112880307
        #   total:      112880307
        #   percent:    100.0
        # 4684.583511 seconds (6.22 G allocations: 306.554 GiB, 77.82% gc time, 0.29% compilation time: 2% of which was recompilation)
        @time accession_to_taxid_table = Mycelia.blastdb2table(blastdb = blastdb, ALL_FIELDS=false, accession=true, taxid=true)
        # 194.124952 seconds (3.76 M allocations: 6.231 GiB, 66.96% gc time, 10.46% compilation time: <1% of which was recompilation)
        accession_to_taxid_table[!, "taxid"] .= parse.(Int, accession_to_taxid_table[!, "taxid"])
        @time Arrow.write(outfile, accession_to_taxid_table)
    else
        accession_to_taxid_table = DataFrames.DataFrame(Arrow.Table(outfile))
        if !(eltype(accession_to_taxid_table[!, "taxid"]) <: Int)
            accession_to_taxid_table[!, "taxid"] .= parse.(Int, accession_to_taxid_table[!, "taxid"])
        end
    end
    return accession_to_taxid_table
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

function streamline_counts(counts; threshold=0.01, min_len=30)
    if length(counts) < min_len
        return counts
    end
    total_counts = sum(last, counts)
    
    # Determine if threshold is relative (float between 0-1) or absolute (integer)
    is_relative = isa(threshold, AbstractFloat) && 0.0 <= threshold <= 1.0
    
    new_counts = Vector{eltype(counts)}()
    other_counts = 0
    
    for (item, count) in counts
        if is_relative
            # For relative threshold, check if the item's proportion is >= threshold
            if count / total_counts >= threshold
                push!(new_counts, item => count)
            else
                other_counts += count
            end
        else
            # For absolute threshold, check if the count is >= threshold
            if count >= threshold
                push!(new_counts, item => count)
            else
                other_counts += count
            end
        end
    end
    
    if other_counts > 0
        push!(new_counts, "Other" => other_counts)
    end
    
    return new_counts
end

function assign_lowest_rank_to_reads_to_taxon_lineage_table(reads_to_taxon_lineage_table)
    reads_to_taxon_lineage_table[!, "lowest_rank"] .= 0
    for (i, row) in enumerate(DataFrames.eachrow(reads_to_taxon_lineage_table))
        if !ismissing(row["species"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 9
        elseif !ismissing(row["genus"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 8
        elseif !ismissing(row["family"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 7
        elseif !ismissing(row["order"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 6
        elseif !ismissing(row["class"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 5
        elseif !ismissing(row["phylum"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 4
        elseif !ismissing(row["kingdom"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 3
        elseif !ismissing(row["realm"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 2
        elseif !ismissing(row["domain"])
            reads_to_taxon_lineage_table[i, "lowest_rank"] = 1
        end
    end
    return reads_to_taxon_lineage_table
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

function classify_xam_with_blast_taxonomies(;xam, accession_to_taxid_table)
    xam_table_columns_of_interest = [
        "template",
        "ismapped",
        "isprimary",
        "flag",
        "reference",
        "position",
        "alignment_score",
    ]
    
    blast_tax_table_columns_of_interest = [
        "accession",
        "taxid",
    ]
    
    # 151.408048 seconds (317.18 k allocations: 460.065 MiB, 0.10% compilation time: 20% of which was recompilation)
    @time xam_table = Mycelia.xam_to_dataframe(xam)

    taxid_aware_xam_table = DataFrames.leftjoin(
        DataFrames.select(xam_table, xam_table_columns_of_interest),
        DataFrames.select(accession_to_taxid_table, blast_tax_table_columns_of_interest),
        on="reference" => "accession",
        matchmissing = :notequal
    )
    
    template_taxid_score_table = DataFrames.combine(DataFrames.groupby(taxid_aware_xam_table, [:template, :taxid]), :alignment_score => sum => :total_alignment_score)
    # replace missing (unclassified) with 0 (NCBI taxonomies start at 1, so 0 is a common, but technically non-standard NCBI taxon identifier
    template_taxid_score_table = DataFrames.coalesce.(template_taxid_score_table, 0)
    
    # For each template, identify top taxid and calculate score difference
    results_df = DataFrames.combine(DataFrames.groupby(template_taxid_score_table, :template)) do group
        # Sort scores in descending order
        sorted = DataFrames.sort(group, :total_alignment_score, rev=true)

        # Get top taxid and score
        top_taxid = sorted[1, :taxid]
        top_score = sorted[1, :total_alignment_score]

        # Calculate ratio with next best (if it exists)
        score_ratio = Inf
        if DataFrames.nrow(sorted) > 1
            second_score = sorted[2, :total_alignment_score]
            score_ratio = top_score/second_score
        end

        # Store all additional taxids and their scores (excluding the top one)
        additional_taxids = OrderedCollections.OrderedDict{Int, Float64}()
        if DataFrames.nrow(sorted) > 1
            for i in 2:DataFrames.nrow(sorted)
                additional_taxids[sorted[i, :taxid]] = sorted[i, :total_alignment_score]
            end
        end

        # Return a new row with the results
        return DataFrames.DataFrame(
            top_taxid = top_taxid,
            top_score = top_score,
            ratio_to_next_best_score = score_ratio,
            additional_taxids = [additional_taxids]  # Wrap in array to make it a single element
        )
    end
    classification_table = Mycelia.apply_conservative_taxonomy(results_df)
    return classification_table
end

function apply_conservative_taxonomy(results_df; ratio_threshold=2.0)
    # Initialize output columns with default values
    n_rows = DataFrames.nrow(results_df)
    final_assignment = copy(results_df.top_taxid)
    confidence_level = fill("high", n_rows)
    
    # Collect all sets of taxids that need LCA calculation
    lca_needed_indices = Int[]
    competing_taxids_list = Vector{Vector{Int}}()
    
    # First pass: identify which rows need LCA and prepare the taxid sets
    for i in 1:n_rows
        top_taxid = results_df.top_taxid[i]
        top_score = results_df.top_score[i]
        ratio = results_df.ratio_to_next_best_score[i]
        additional_dict = results_df.additional_taxids[i]
        
        # Skip if no competitors or ratio is high enough
        if isempty(additional_dict) || ratio >= ratio_threshold
            continue
        end
        
        # Collect taxids that are within the threshold
        competing_taxids = [top_taxid]
        for (taxid, score) in additional_dict
            if top_score / score < ratio_threshold
                push!(competing_taxids, taxid)
            end
        end
        
        # If we have multiple competing taxids, add to the batch
        if length(competing_taxids) > 1
            push!(lca_needed_indices, i)
            push!(competing_taxids_list, competing_taxids)
            confidence_level[i] = "lca"  # Mark as needing LCA
        end
    end
    
    # If we have any rows that need LCA calculation
    if !isempty(lca_needed_indices)
        # Batch calculate all LCAs
        lca_results = Mycelia.batch_taxids2lca(competing_taxids_list)
        
        # Apply LCA results to the appropriate rows
        for (idx, lca_idx) in enumerate(lca_needed_indices)
            final_assignment[lca_idx] = lca_results[idx]
        end
    end
    
    # Create a new dataframe with the original data plus our new columns
    updated_results = DataFrames.hcat(
        results_df,
        DataFrames.DataFrame(
            final_assignment = final_assignment,
            confidence_level = confidence_level
        )
    )
    
    return updated_results
end


# this is faster than NCBI version
# run(pipeline(
#         `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon 1 --children --as-json-lines`,
#         `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli dataformat tsv taxonomy --template tax-summary`
#     )
# )
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Retrieves and formats the complete NCBI taxonomy hierarchy into a structured DataFrame.

# Details
- Automatically sets up taxonkit environment and downloads taxonomy database if needed
- Starts from root taxid (1) and includes all descendant taxa
- Reformats lineage information into separate columns for each taxonomic rank

# Returns
DataFrame with columns:
- `taxid`: Taxonomy identifier
- `lineage`: Full taxonomic lineage string
- `taxid_lineage`: Lineage with taxonomy IDs
- Individual rank columns:
  - superkingdom, kingdom, phylum, class, order, family, genus, species
  - corresponding taxid columns (e.g., superkingdom_taxid)

# Dependencies
Requires taxonkit (installed automatically via Bioconda)
"""
function list_full_taxonomy()
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    p = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --add-prefix --fill-miss-rank --show-lineage-taxids --format '{k};{K};{p};{c};{o};{f};{g};{s}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "lineage", "taxid_lineage"]
    table = DataFrames.DataFrame(data, header)
    ranks = [
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]
    for rank in ranks
        table[!, rank] .= ""
        table[!, "$(rank)_taxid"] .= ""
    end

    for (i, row) in enumerate(DataFrames.eachrow(table))
        for (rank, x) in zip(ranks, split(row["lineage"], ';'))
            table[i, rank] = x
        end
        for (rank, x) in zip(ranks, split(row["taxid_lineage"], ';'))
            table[i, "$(rank)_taxid"] = x
        end
    end
    return table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return an ordered list of taxonomic ranks from highest (top) to lowest (species).

# Arguments
- `synonyms::Bool=false`: If true, includes alternative names for certain ranks (e.g. "domain" for "superkingdom")

# Returns
- `Vector{String}`: An array of taxonomic rank names in hierarchical order
"""
function list_ranks(;synonyms=false)
    if !synonyms
        return [
            "top",
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        ]
    else
        return [
            "top",
            "superkingdom/domain",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        ]
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a DataFrame containing the top-level taxonomic nodes.

The DataFrame has two fixed rows representing the most basic taxonomic classifications:
- taxid=0: "unclassified"
- taxid=1: "root"

Returns
-------
DataFrame
    Columns:
    - taxid::Int : Taxonomic identifier
    - name::String : Node name
"""
function list_toplevel()
    return DataFrames.DataFrame(taxid=[0, 1], name=["unclassified", "root"])
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

List all taxonomic entries at the specified rank level.

# Arguments
- `rank::String`: Taxonomic rank to query. Must be one of:
  - "top" (top level)
  - "superkingdom"/"domain"  
  - "kingdom"
  - "phylum" 
  - "class"
  - "order"
  - "family"
  - "genus"
  - "species"

# Returns
DataFrame with columns:
- `taxid`: NCBI taxonomy ID
- `name`: Scientific name at the specified rank
"""
function list_rank(rank)
    if rank == "top"
        return list_toplevel()
    else
        Mycelia.add_bioconda_env("taxonkit")
        if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
            setup_taxonkit_taxonomy()
        end
        ranks_to_shorthand = Dict(
            "superkingdom" => "k",
            "kingdom" => "K",
            "phylum" => "p",
            "class" => "c",
            "order" => "o",
            "family" => "f",
            "genus" => "g",
            "species" => "s"
        )
        shorthand = ranks_to_shorthand[rank]
        p = pipeline(
            `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit list --ids 1`,
            `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit filter --equal-to "$(rank)"`,
            `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format "{$(shorthand)}"`
        )
        data, header = uCSV.read(open(p), delim='\t')
        header = ["taxid", "name"]
        return DataFrames.DataFrame(data, header)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns an array of all taxonomic superkingdoms (e.g., Bacteria, Archaea, Eukaryota).

# Returns
- `Vector{String}`: Array containing names of all superkingdoms in the taxonomy database
"""
function list_superkingdoms()
    return list_rank("superkingdom")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists all taxonomic kingdoms in the database.

Returns a vector of kingdom names as strings. Kingdoms represent the highest
major taxonomic rank in biological classification.

# Returns
- `Vector{String}`: Array of kingdom names
"""
function list_kingdoms()
    return list_rank("kingdom")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted list of all unique phyla in the database.
"""
function list_phylums()
    return list_rank("phylum")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns an array of all taxonomic classes in the database.

Classes represent a major taxonomic rank between phylum and order in biological classification.

# Returns
- `Vector{String}`: Array of class names sorted alphabetically
"""
function list_classes()
    return list_rank("class")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists all orders in the taxonomic database.

Returns a vector of strings containing valid order names according to current mycological taxonomy.
Uses the underlying `list_rank()` function with rank="order".

# Returns
- `Vector{String}`: Alphabetically sorted list of order names
"""
function list_orders()
    return list_rank("order")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted vector of all family names present in the database.
"""
function list_families()
    return list_rank("family")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted vector of all genera names present in the database.
"""
function list_genera()
    return list_rank("genus")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted vector of all species names present in the database.
"""
function list_species()
    return list_rank("species")
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function list_subtaxa(taxid)
#     return parse.(Int, filter(!isempty, strip.(readlines(`conda run --no-capture-output -n taxonkit taxonkit list --ids 10239`))))
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns an array of Integer taxon IDs representing all sub-taxa under the specified taxonomic ID.

# Arguments
- `taxid`: NCBI taxonomy identifier for the parent taxon

# Returns
Vector{Int} containing all descendant taxon IDs

# Details
- Requires taxonkit to be installed via Bioconda
- Automatically sets up taxonkit database if not present
- Uses local taxonomy database in ~/.taxonkit/
"""
function list_subtaxa(taxid)
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    return parse.(Int, filter(!isempty, strip.(readlines(`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit list --ids $(taxid)`))))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert scientific name(s) to NCBI taxonomy ID(s) using taxonkit.

# Arguments
- `name::AbstractString`: Scientific name(s) to query. Can be a single name or multiple names separated by newlines.

# Returns
- `DataFrame` with columns:
  - `name`: Input scientific name
  - `taxid`: NCBI taxonomy ID
  - `rank`: Taxonomic rank (e.g., "species", "genus")

# Dependencies
Requires taxonkit package (installed automatically via Bioconda)
"""
function name2taxid(name)
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    p = pipeline(`echo $(name)`, `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit name2taxid --show-rank`)
    data, header = uCSV.read(open(p), delim='\t')
    header = ["name", "taxid", "rank"]
    return DataFrames.DataFrame(data, header)
end

# other ncbi-datasets reports that I didn't find as useful initially
# run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon 10114 --report names`))
# x = JSON.parse(open(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon "rattus norvegicus"`)))
# run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download taxonomy taxon 33554 --children`))

# more useful
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a vector of NCBI taxonomy IDs into a detailed taxonomy table using NCBI Datasets CLI.

# Arguments
- `taxids::AbstractVector{Int}`: Vector of NCBI taxonomy IDs to query

# Returns
- `DataFrame`: Table containing taxonomy information with columns including:
  - tax_id
  - species
  - genus
  - family
  - order
  - class
  - phylum
  - kingdom

# Dependencies
Requires ncbi-datasets-cli Conda package (automatically installed if missing)
"""
function taxids2ncbi_taxonomy_table(taxids::AbstractVector{Int})
    Mycelia.add_bioconda_env("ncbi-datasets-cli")
    # joint_table = DataFrames.DataFrame()
    # ProgressMeter.@showprogress for taxid in taxids
    #     cmd1 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon $(taxid) --as-json-lines`
    #     cmd2 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli dataformat tsv taxonomy --template tax-summary`
    #     io = open(pipeline(cmd1, cmd2))
    #     try
    #         append!(joint_table, CSV.read(io, DataFrames.DataFrame, delim='\t', header=1), promote=true)
    #     catch e
    #         error("unable to process taxid: $(taxid)\n$(e)")
    #     end
    # end
    temp_file = tempname() * ".taxonids.txt"
    unique_taxids = sort(unique(taxids))
    open(temp_file, "w") do io
        for taxid in unique_taxids
            println(io, taxid)
        end
    end
    cmd1 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon --inputfile $(temp_file) --as-json-lines`
    cmd2 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli dataformat tsv taxonomy --template tax-summary`
    io = open(pipeline(cmd1, cmd2))
    joint_table = CSV.read(io, DataFrames.DataFrame, delim='\t', header=1)
    rm(temp_file)    
    return joint_table
end

# more complete
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert NCBI taxonomic IDs to their complete taxonomic lineage information using taxonkit.

# Arguments
- `taxids::AbstractVector{Int}`: Vector of NCBI taxonomy IDs

# Returns
A DataFrame with columns:
- `taxid`: Original query taxonomy ID
- `lineage`: Full taxonomic lineage as semicolon-separated string
- `lineage-taxids`: Corresponding taxonomy IDs for each rank in lineage
- `lineage-ranks`: Taxonomic ranks for each level in lineage
"""
# function taxids2taxonkit_lineage_table(taxids::AbstractVector{Int})
function taxids2taxonkit_full_lineage_table(taxids::AbstractVector{Int})
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    f = tempname()
    open(f, "w") do io
        for taxid in taxids
            println(io, taxid)
        end
    end
    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n taxonkit taxonkit lineage --show-lineage-taxids --show-lineage-ranks $(f)`
    data, header = uCSV.read(open(pipeline(cmd)), delim='\t', header=false, typedetectrows=100)
    rm(f)
    header = ["taxid", "lineage", "lineage-taxids", "lineage-ranks"]
    return DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert taxonomic IDs to a structured lineage rank mapping.

Takes a vector of taxonomic IDs and returns a nested dictionary mapping each input taxid 
to its complete taxonomic lineage information. For each taxid, creates a dictionary where:
- Keys are taxonomic ranks (e.g., "species", "genus", "family")
- Values are NamedTuples containing:
  - `lineage::String`: The taxonomic name at that rank
  - `taxid::Union{Int, Missing}`: The corresponding taxonomic ID (if available)

Excludes "no rank" entries from the final output.

Returns:
    Dict{Int, Dict{String, NamedTuple{(:lineage, :taxid), Tuple{String, Union{Int, Missing}}}}}
"""
function taxids2taxonkit_taxid2lineage_ranks(taxids::AbstractVector{Int})
    table = taxids2taxonkit_full_lineage_table(taxids)
    # table = taxids2taxonkit_lineage_table(taxids)
    taxid_to_lineage_ranks = Dict{Int, Dict{String, @NamedTuple{lineage::String, taxid::Union{Int, Missing}}}}()
    for row in DataFrames.eachrow(table)
        lineage_ranks = String.(split(row["lineage-ranks"], ';'))
        lineage_taxids = [something(tryparse(Int, x), missing) for x in split(row["lineage-taxids"], ';')]
        # lineage_taxids = something.(tryparse.(Int, split(row["lineage-taxids"], ';')), missing)
        lineage = String.(split(row["lineage"], ';'))
        row_dict = Dict(rank => (;lineage, taxid) for (lineage, rank, taxid) in zip(lineage, lineage_ranks, lineage_taxids))
        delete!(row_dict, "no rank")
        taxid_to_lineage_ranks[row["taxid"]] = row_dict
    end
    return taxid_to_lineage_ranks
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a vector of taxonomy IDs to a summarized lineage table using taxonkit.

# Arguments
- `taxids::AbstractVector{Int}`: Vector of NCBI taxonomy IDs

# Returns
DataFrame with the following columns:
- `taxid`: Original input taxonomy ID
- `species_taxid`, `species`: Species level taxonomy ID and name
- `genus_taxid`, `genus`: Genus level taxonomy ID and name  
- `family_taxid`, `family`: Family level taxonomy ID and name
- `superkingdom_taxid`, `superkingdom`: Superkingdom level taxonomy ID and name

Missing values are used when a taxonomic rank is not available.
"""
function taxids2taxonkit_summarized_lineage_table(taxids::AbstractVector{Int})
    taxid_to_lineage_ranks = taxids2taxonkit_taxid2lineage_ranks(taxids)
    taxids_to_lineage_table = DataFrames.DataFrame()
    for (taxid, lineage_ranks) in taxid_to_lineage_ranks
        # 
        row = (
            taxid = taxid,
            species_taxid = haskey(lineage_ranks, "species") ? lineage_ranks["species"].taxid : missing,
            species = haskey(lineage_ranks, "species") ? lineage_ranks["species"].lineage : missing,
            
            genus_taxid = haskey(lineage_ranks, "genus") ? lineage_ranks["genus"].taxid : missing,
            genus = haskey(lineage_ranks, "genus") ? lineage_ranks["genus"].lineage : missing,
            
            family_taxid = haskey(lineage_ranks, "family") ? lineage_ranks["family"].taxid : missing,
            family = haskey(lineage_ranks, "family") ? lineage_ranks["family"].lineage : missing,
            
            order_taxid = haskey(lineage_ranks, "order") ? lineage_ranks["order"].taxid : missing,
            order = haskey(lineage_ranks, "order") ? lineage_ranks["order"].lineage : missing,

            class_taxid = haskey(lineage_ranks, "class") ? lineage_ranks["class"].taxid : missing,
            class = haskey(lineage_ranks, "class") ? lineage_ranks["class"].lineage : missing,
            
            phylum_taxid = haskey(lineage_ranks, "phylum") ? lineage_ranks["phylum"].taxid : missing,
            phylum = haskey(lineage_ranks, "phylum") ? lineage_ranks["phylum"].lineage : missing,
            
            kingdom_taxid = haskey(lineage_ranks, "kingdom") ? lineage_ranks["kingdom"].taxid : missing,
            kingdom = haskey(lineage_ranks, "kingdom") ? lineage_ranks["kingdom"].lineage : missing,
            
            realm_taxid = haskey(lineage_ranks, "realm") ? lineage_ranks["realm"].taxid : missing,
            realm = haskey(lineage_ranks, "realm") ? lineage_ranks["realm"].lineage : missing,
            
            domain_taxid = haskey(lineage_ranks, "domain") ? lineage_ranks["domain"].taxid : missing,
            domain = haskey(lineage_ranks, "domain") ? lineage_ranks["domain"].lineage : missing,

        )
        push!(taxids_to_lineage_table, row, promote=true)
    end
    return taxids_to_lineage_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Lowest Common Ancestor (LCA) taxonomic ID for a set of input taxonomic IDs.

# Arguments
- `ids::Vector{Int}`: Vector of NCBI taxonomic IDs

# Returns
- `Int`: The taxonomic ID of the lowest common ancestor

# Details
Uses taxonkit to compute the LCA. Automatically sets up the required taxonomy database 
if not already present in `~/.taxonkit/`.

# Dependencies
- Requires taxonkit (installed via Bioconda)
- Requires taxonomy database (downloaded automatically if missing)
"""
function taxids2lca(ids::Vector{Int})
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    # Convert the list of integers to a space-separated string
    input_str = join(ids, " ")

    # Pass the input string to the `taxonkit lca` command and capture the output
    output = read(pipeline(`echo $(input_str)`, `$(Mycelia.CONDA_RUNNER) run --live-stream -n taxonkit taxonkit lca`), String)

    # Split the output string and select the last item
    lca_id = split(chomp(output), "\t")[end]

    # Convert the LCA identifier to an integer and return it
    return parse(Int, lca_id)
end

function batch_taxids2lca(ids_list::Vector{Vector{Int}})
    # Ensure taxonkit environment and taxonomy are set up
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end

    # Write the queries to a temporary file: one query per line
    tmpfile = tempname()
    open(tmpfile, "w") do io
        for ids in ids_list
            # Join taxids with a space (TaxonKit lca accepts space-separated IDs)
            println(io, join(ids, " "))
        end
    end

    # Run taxonkit lca in batch mode, reading input from the temporary file.
    # The command here uses `cat tmpfile | taxonkit lca`
    cmd = pipeline(`cat $(tmpfile)`, `$(Mycelia.CONDA_RUNNER) run --live-stream -n taxonkit taxonkit lca`)
    output_str = read(cmd, String)

    # Clean up temporary file
    rm(tmpfile)

    # Parse the output: assume each line corresponds to one input query.
    # TaxonKit's output is expected to be tab-separated, with the last field being the LCA.
    lca_lines = split(chomp(output_str), "\n")
    result = Vector{Int}(undef, length(lca_lines))
    for (i, line) in enumerate(lca_lines)
        fields = split(line, "\t")
        result[i] = parse(Int, fields[end])
    end

    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a vector of species/taxon names to their corresponding NCBI taxonomy IDs.

# Arguments
- `names::AbstractVector{<:AbstractString}`: Vector of scientific names or common names

# Returns
- `Vector{Int}`: Vector of NCBI taxonomy IDs corresponding to the input names

Progress is displayed using ProgressMeter.
"""
function names2taxids(names::AbstractVector{<:AbstractString})
    results = []
    ProgressMeter.@showprogress for name in names
        push!(results, Mycelia.name2taxid(name))
    end
    return reduce(vcat, results)
end