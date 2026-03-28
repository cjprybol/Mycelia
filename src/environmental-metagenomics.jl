"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return a curated catalog of public environmental metagenome datasets that map
well onto Mycelia workflows.

The initial catalog is intentionally small and focuses on a single well-known
marine case study so that tutorials can stay reproducible while still pointing
to a real public dataset.
"""
function environmental_metagenome_dataset_catalog()
    return DataFrames.DataFrame(
        dataset_id = ["tara_oceans_surface_prokaryotes"],
        title = ["Tara Oceans surface prokaryote metagenomes"],
        public_accession = ["PRJEB1787"],
        environment = ["marine"],
        sample_scope = ["Surface, prokaryote-enriched 0.2-3 um fraction"],
        starter_runs = ["ERR598952,ERR598957"],
        source_url = ["https://www.ebi.ac.uk/ena/browser/view/PRJEB1787"],
        metadata_url = ["https://ocean-microbiome.embl.de/companion.html"],
        notes = ["Tutorial uses a bundled lightweight abundance table inspired by this study so the workflow stays runnable offline."]
    )
end

function _tara_oceans_surface_sample_metadata()
    return DataFrames.DataFrame(
        sample = [
            "TARA_NATL_001",
            "TARA_SATL_001",
            "TARA_NPAC_001",
            "TARA_SPAC_001",
            "TARA_IND_001",
            "TARA_SOCE_001"
        ],
        region = [
            "North Atlantic",
            "South Atlantic",
            "North Pacific",
            "South Pacific",
            "Indian Ocean",
            "Southern Ocean"
        ],
        biome = [
            "temperate gyre",
            "subtropical gyre",
            "temperate gyre",
            "subtropical gyre",
            "tropical gyre",
            "polar"
        ],
        depth_m = [5, 5, 5, 5, 5, 5],
        temperature_c = [18.4, 23.1, 16.8, 24.7, 27.2, 3.8],
        size_fraction = fill("0.2-3 um", 6),
        study_accession = fill("PRJEB1787", 6)
    )
end

function _tara_oceans_surface_abundance_table()
    taxa = [
        "Pelagibacter",
        "Prochlorococcus",
        "Synechococcus",
        "Roseobacter",
        "Flavobacteriaceae",
        "Alteromonas",
        "Oceanospirillales",
        "Thaumarchaeota"
    ]

    sample_names = [
        "TARA_NATL_001",
        "TARA_SATL_001",
        "TARA_NPAC_001",
        "TARA_SPAC_001",
        "TARA_IND_001",
        "TARA_SOCE_001"
    ]

    counts = [
        4200 4600 3900 4100 4300 3100
        1800 2500 1500 2700 2400 400
        1100 800 1400 700 900 1200
        900 700 850 650 780 950
        700 400 1200 350 500 1800
        500 350 650 300 420 900
        300 250 500 220 260 650
        120 80 140 70 90 300
    ]

    abundance_table = DataFrames.DataFrame(sample = String[], taxon = String[], count = Int[])
    for (taxon_index, taxon) in enumerate(taxa)
        for (sample_index, sample_name) in enumerate(sample_names)
            push!(
                abundance_table,
                (sample = sample_name, taxon = taxon, count = counts[taxon_index, sample_index])
            )
        end
    end

    return abundance_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load a small bundled environmental metagenome example that mirrors the shape of
an analysis-ready taxonomic abundance table.

The example is tied to a real public dataset entry from
`environmental_metagenome_dataset_catalog()` but ships as lightweight counts and
metadata so it can be executed in tests, docs, and tutorials without network
access.
"""
function environmental_metagenome_example_data(;
        dataset_id::String = "tara_oceans_surface_prokaryotes")
    catalog = environmental_metagenome_dataset_catalog()
    dataset_info = catalog[catalog.dataset_id .== dataset_id, :]

    if DataFrames.nrow(dataset_info) != 1
        error("Unknown environmental metagenome dataset id: $(dataset_id)")
    end

    if dataset_id == "tara_oceans_surface_prokaryotes"
        return (
            dataset_info = dataset_info,
            sample_metadata = _tara_oceans_surface_sample_metadata(),
            abundance_table = _tara_oceans_surface_abundance_table()
        )
    end

    error("No bundled example is available for dataset id: $(dataset_id)")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a long-form abundance table into a taxa x samples matrix suitable for
diversity calculations and ordination.
"""
function prepare_environmental_abundance_matrix(
        abundance_df::DataFrames.DataFrame;
        sample_col::Symbol = :sample,
        taxon_col::Symbol = :taxon,
        abundance_col::Symbol = :count
)
    sample_names = collect(unique(abundance_df[!, sample_col]))
    taxon_names = collect(unique(abundance_df[!, taxon_col]))

    sample_to_index = Dict(sample_name => i for (i, sample_name) in enumerate(sample_names))
    taxon_to_index = Dict(taxon_name => i for (i, taxon_name) in enumerate(taxon_names))

    abundance_matrix = zeros(Float64, length(taxon_names), length(sample_names))
    for row in DataFrames.eachrow(abundance_df)
        abundance_matrix[
            taxon_to_index[row[taxon_col]],
            sample_to_index[row[sample_col]]
        ] = Float64(row[abundance_col])
    end

    return (
        matrix = abundance_matrix,
        sample_names = sample_names,
        taxon_names = taxon_names
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert count abundances to per-sample relative abundances.
"""
function environmental_relative_abundance_table(
        abundance_df::DataFrames.DataFrame;
        sample_col::Symbol = :sample,
        abundance_col::Symbol = :count
)
    totals = DataFrames.combine(
        DataFrames.groupby(abundance_df, sample_col),
        abundance_col => sum => :sample_total
    )
    relative = DataFrames.leftjoin(abundance_df, totals, on = sample_col)
    relative[!, :relative_abundance] = relative[!, abundance_col] ./ relative.sample_total
    return DataFrames.select(relative, DataFrames.Not(:sample_total))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Summarize the most abundant taxa across all samples in an environmental
metagenome table.
"""
function summarize_environmental_taxa(
        abundance_df::DataFrames.DataFrame;
        taxon_col::Symbol = :taxon,
        abundance_col::Symbol = :count,
        top_n::Int = 8
)
    summarized = DataFrames.combine(
        DataFrames.groupby(abundance_df, taxon_col),
        abundance_col => sum => :total_count
    )
    DataFrames.sort!(summarized, :total_count, rev = true)
    return first(summarized, min(top_n, DataFrames.nrow(summarized)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run alpha and beta diversity analysis for an environmental metagenome abundance
table and attach sample metadata for plotting.
"""
function analyze_environmental_metagenome(
        abundance_df::DataFrames.DataFrame,
        sample_metadata::DataFrames.DataFrame;
        sample_col::Symbol = :sample,
        taxon_col::Symbol = :taxon,
        abundance_col::Symbol = :count,
        metadata_sample_col::Symbol = :sample,
        metric::Symbol = :bray_curtis
)
    prepared = prepare_environmental_abundance_matrix(
        abundance_df;
        sample_col = sample_col,
        taxon_col = taxon_col,
        abundance_col = abundance_col
    )

    alpha_diversity = calculate_alpha_diversity(prepared.matrix, prepared.sample_names)
    alpha_diversity = DataFrames.leftjoin(
        alpha_diversity,
        sample_metadata,
        on = :sample => metadata_sample_col
    )

    beta_diversity = beta_diversity_pcoa(
        prepared.matrix,
        prepared.sample_names;
        metric = metric,
        metadata = sample_metadata
    )

    return (
        abundance_matrix = prepared.matrix,
        sample_names = prepared.sample_names,
        taxon_names = prepared.taxon_names,
        relative_abundance = environmental_relative_abundance_table(
            abundance_df;
            sample_col = sample_col,
            abundance_col = abundance_col
        ),
        alpha_diversity = alpha_diversity,
        beta_diversity = beta_diversity,
        top_taxa = summarize_environmental_taxa(
            abundance_df;
            taxon_col = taxon_col,
            abundance_col = abundance_col
        )
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a scatter plot that summarizes alpha diversity across environmental
samples.
"""
function plot_environmental_alpha_diversity(
        alpha_diversity::DataFrames.DataFrame;
        group_col::Symbol = :region,
        title::String = "Environmental Metagenome Alpha Diversity"
)
    return StatsPlots.scatter(
        alpha_diversity.richness,
        alpha_diversity.shannon;
        group = alpha_diversity[!, group_col],
        xlabel = "Observed richness",
        ylabel = "Shannon diversity",
        title = title,
        legend = :outerright,
        markerstrokewidth = 0,
        markersize = 7
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a PC1 vs PC2 scatter plot from an environmental beta-diversity result.
"""
function plot_environmental_beta_diversity(
        beta_diversity_result;
        group_col::Symbol = :region,
        title::String = "Environmental Metagenome Beta Diversity"
)
    pcoa_df = beta_diversity_result.pcoa_df
    variance_explained = beta_diversity_result.variance_explained
    pc1_percent = round(100 * variance_explained[1], digits = 1)
    pc2_percent = round(100 * variance_explained[2], digits = 1)

    return StatsPlots.scatter(
        pcoa_df.PC1,
        pcoa_df.PC2;
        group = pcoa_df[!, group_col],
        xlabel = "PC1 ($(pc1_percent)%)",
        ylabel = "PC2 ($(pc2_percent)%)",
        title = title,
        legend = :outerright,
        markerstrokewidth = 0,
        markersize = 7
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate publication-ready figures for the bundled environmental metagenome
analysis and optionally save them to disk.
"""
function generate_environmental_metagenome_figures(
        analysis_result;
        title_prefix::String = "Environmental Metagenome",
        output_dir::Union{Nothing, String} = nothing,
        top_n_taxa::Int = 8
)
    abundance_views = plot_microbiome_abundance(
        analysis_result.relative_abundance;
        sample_col = :sample,
        taxon_col = :taxon,
        abundance_col = :relative_abundance,
        rank = "genus",
        title = "$(title_prefix): relative abundance",
        config = MicrobiomePlotConfig(
            top_n_taxa = top_n_taxa,
            sample_ordering = AxisOrdering(method = :alphabetical),
            output_formats = [:png, :svg]
        )
    )

    alpha_plot = plot_environmental_alpha_diversity(
        analysis_result.alpha_diversity;
        title = "$(title_prefix): alpha diversity"
    )

    beta_plot = plot_environmental_beta_diversity(
        analysis_result.beta_diversity;
        title = "$(title_prefix): Bray-Curtis PCoA"
    )

    saved_paths = Dict{Symbol, Dict{Symbol, String}}()
    if output_dir !== nothing
        saved_paths[:abundance_barplot] = save_plot(
            abundance_views[:barplot],
            joinpath(output_dir, "environmental_abundance")
        )
        saved_paths[:alpha_diversity] = save_plot(
            alpha_plot,
            joinpath(output_dir, "environmental_alpha_diversity")
        )
        saved_paths[:beta_diversity_pcoa] = save_plot(
            beta_plot,
            joinpath(output_dir, "environmental_beta_diversity_pcoa")
        )
    end

    return (
        abundance_views = abundance_views,
        alpha_diversity_plot = alpha_plot,
        beta_diversity_plot = beta_plot,
        saved_paths = saved_paths
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the full bundled environmental metagenome case study from dataset selection
through figure generation.
"""
function run_environmental_metagenome_case_study(;
        dataset_id::String = "tara_oceans_surface_prokaryotes",
        metric::Symbol = :bray_curtis,
        output_dir::Union{Nothing, String} = nothing,
        top_n_taxa::Int = 8)
    example_data = environmental_metagenome_example_data(; dataset_id = dataset_id)

    analysis_result = analyze_environmental_metagenome(
        example_data.abundance_table,
        example_data.sample_metadata;
        metric = metric
    )

    dataset_title = example_data.dataset_info.title[1]
    figures = generate_environmental_metagenome_figures(
        analysis_result;
        title_prefix = dataset_title,
        output_dir = output_dir,
        top_n_taxa = top_n_taxa
    )

    return (
        dataset_info = example_data.dataset_info,
        sample_metadata = example_data.sample_metadata,
        abundance_table = example_data.abundance_table,
        analysis = analysis_result,
        figures = figures
    )
end
