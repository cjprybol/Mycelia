function aggregate_by_rank_nonmissing(df::DataFrames.DataFrame, ranks::Vector{String}=[
        "domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"
    ])
    results = DataFrames.DataFrame()
    for (i, rank) in enumerate(ranks)
        rank_taxid = string(rank, "_taxid")
        if rank_taxid in names(df)
            # Drop rows where rank_taxid is missing
            subdf = DataFrames.filter(row -> !ismissing(row[rank_taxid]), df)
            if DataFrames.nrow(subdf) > 0
                grouped = DataFrames.groupby(subdf, ["template", rank_taxid])
                agg = DataFrames.combine(grouped,
                    :relative_alignment_proportion => sum => :total_relative_alignment_proportion,
                    :percent_identity => maximum => :max_percent_identity,
                    :mappingquality => maximum => :max_mappingquality,
                )
                DataFrames.rename!(agg, Dict(rank_taxid => :rank_taxid))
                agg[!, :rank] .= "$(i)_$(rank)"
                results = DataFrames.vcat(results, agg)
            end
        end
    end
    return results
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

"""
Merge XAM alignment data with taxonomic information and calculate alignment metrics.

This function processes XAM alignment data by:
1. Loading an accession-to-taxid mapping table
2. Left-joining alignment data with taxonomic IDs
3. Retrieving full taxonomic lineage information
4. Calculating percent identity scores
5. Calculating relative alignment score proportions per read
6. Writing results to cached files

# Arguments
- `xam`: Path to XAM file or XAM data structure
- `accession2taxid_file`: Path to accession-to-taxid mapping file (.tsv.gz)
- `output_prefix`: Prefix for output files (.tsv.gz). Defaults to "xam"
- `verbose::Bool=true`: Whether to print progress information
- `force_recalculate::Bool=false`: Whether to force recalculation even if cached files exist

# Returns
Paths to the output file: tsv_out

Used to return an Arrow format but Arrow integration is poor with union and custom datatypes in Julia
"""
function merge_xam_with_taxonomies(;
    xam, 
    accession2taxid_file::String,
    output_prefix::String = xam,
    verbose::Bool = true,
    force_recalculate::Bool = false
)
    # Setup output file paths
    tsv_out = output_prefix * ".taxonomy-aware.tsv.gz"
    # arrow_out = output_prefix * ".taxonomy-aware.arrow"
    
    # # Check if we can use cached results
    # use_cache = (!force_recalculate && 
    #              isfile(tsv_out) && 
    #              isfile(arrow_out) && 
    #              filesize(tsv_out) > 0 && 
    #              filesize(arrow_out) > 0)
    # Check if we can use cached results
    use_cache = (!force_recalculate && 
                 isfile(tsv_out) && 
                 filesize(tsv_out) > 0)
    
    if use_cache
        verbose && println("Using cached taxonomy-aware data from: $tsv_out")
        # return (tsv_out=tsv_out, arrow_out=arrow_out)
        return tsv_out
    end
    
    # Remove existing output files if they exist to avoid conflicts
    if isfile(tsv_out)
        verbose && println("Removing existing file: $tsv_out")
        rm(tsv_out)
    end
    # if isfile(arrow_out)
    #     verbose && println("Removing existing file: $arrow_out")
    #     rm(arrow_out)
    # end
    
    # Load accession to taxid mapping table
    verbose && println("Loading accession-to-taxid mapping table...")
    accession2taxid_table = if endswith(accession2taxid_file, ".arrow")
        DataFrames.DataFrame(Arrow.Table(accession2taxid_file))
    elseif endswith(accession2taxid_file, ".tsv.gz")
        Mycelia.read_tsvgz(accession2taxid_file)
    else
        error("Unsupported file format for accession2taxid_file. Expected .tsv.gz or .arrow")
    end
    
    # Convert XAM to DataFrame
    verbose && println("Converting XAM data to DataFrame...")
    xam_table = if verbose
        @time Mycelia.xam_to_dataframe(xam)
    else
        Mycelia.xam_to_dataframe(xam)
    end
    
    # Left join with taxonomic IDs
    verbose && println("Joining alignment data with taxonomic IDs...")
    taxid_aware_xam_table = if verbose
        @time DataFrames.leftjoin(
            xam_table,
            accession2taxid_table,
            on = "reference" => "accession",
            matchmissing = :notequal
        )
    else
        DataFrames.leftjoin(
            xam_table,
            accession2taxid_table,
            on = "reference" => "accession",
            matchmissing = :notequal
        )
    end
    
    # Get unique taxids and retrieve lineage information (filter out missing values)
    verbose && println("Retrieving taxonomic lineage information...")
    unique_taxids = convert(Vector{Int}, filter(!ismissing, unique(taxid_aware_xam_table[!, "taxid"])))
    
    if !isempty(unique_taxids)
        xam_taxid2lineage_table = if verbose
            @time Mycelia.taxids2taxonkit_summarized_lineage_table(unique_taxids)
        else
            Mycelia.taxids2taxonkit_summarized_lineage_table(unique_taxids)
        end
        
        # Add full taxonomic lineage
        verbose && println("Adding full taxonomic lineage information...")
        fully_taxonomy_aware_xam_table = if verbose
            @time DataFrames.leftjoin(
                taxid_aware_xam_table,
                xam_taxid2lineage_table,
                on="taxid",
                matchmissing = :notequal
            )
        else
            DataFrames.leftjoin(
                taxid_aware_xam_table,
                xam_taxid2lineage_table,
                on="taxid",
                matchmissing = :notequal
            )
        end
    else
        verbose && println("No valid taxids found, skipping lineage lookup...")
        fully_taxonomy_aware_xam_table = taxid_aware_xam_table
    end
    
    # Calculate percent identity
    verbose && println("Calculating percent identity scores...")
    if "alignlength" in names(fully_taxonomy_aware_xam_table) && "mismatches" in names(fully_taxonomy_aware_xam_table)
        DataFrames.transform!(fully_taxonomy_aware_xam_table, 
            [:alignlength, :mismatches] => 
            ((len, mis) -> begin
                # Handle missing values by returning missing for percent identity
                result = similar(len, Union{Float64, Missing})
                for i in eachindex(len, mis)
                    if ismissing(len[i]) || ismissing(mis[i])
                        result[i] = missing
                    else
                        result[i] = (len[i] - mis[i]) / len[i] * 100
                    end
                end
                result
            end) => 
            :percent_identity)
    else
        verbose && println("Warning: alignlength and/or mismatches columns not found. Skipping percent identity calculation.")
    end
    
    # Calculate relative alignment score proportions per read/template
    verbose && println("Calculating relative alignment score proportions...")
    if "template" in names(fully_taxonomy_aware_xam_table) && "alignment_score" in names(fully_taxonomy_aware_xam_table)
        # Group by template and calculate total scores per template
        template_totals = DataFrames.combine(
            DataFrames.groupby(fully_taxonomy_aware_xam_table, :template), 
            :alignment_score => sum => :total_template_score
        )
        
        # Merge back to get proportions
        fully_taxonomy_aware_xam_table = DataFrames.leftjoin(
            fully_taxonomy_aware_xam_table,
            template_totals,
            on = :template
        )
        
        # Calculate relative proportion
        DataFrames.transform!(fully_taxonomy_aware_xam_table,
            [:alignment_score, :total_template_score] =>
            ((score, total) -> begin
                # Handle missing values properly
                result = similar(score, Union{Float64, Missing})
                for i in eachindex(score, total)
                    if ismissing(score[i]) || ismissing(total[i]) || total[i] == 0
                        result[i] = missing
                    else
                        result[i] = score[i] / total[i]
                    end
                end
                result
            end) =>
            :relative_alignment_proportion
        )
        
        # Clean up the temporary total column
        DataFrames.select!(fully_taxonomy_aware_xam_table, DataFrames.Not(:total_template_score))
    else
        verbose && println("Warning: template and/or alignment_score columns not found. Skipping relative proportion calculation.")
    end

    # sanitize nothings before writing out
    verbose && println("Sanitizing `nothing` values and replacing with `missing`")
    fully_taxonomy_aware_xam_table = Mycelia.dataframe_replace_nothing_with_missing(fully_taxonomy_aware_xam_table)
    
    # Write output files
    # verbose && println("Writing results to $tsv_out and $arrow_out...")
    verbose && println("Writing results to $tsv_out...")
    
    if verbose
        @time Mycelia.write_tsvgz(df=fully_taxonomy_aware_xam_table, filename=tsv_out)
        GC.gc()  # Clean up after TSV write
        # @time Arrow.write(arrow_out, fully_taxonomy_aware_xam_table)
    else
        Mycelia.write_tsvgz(df=fully_taxonomy_aware_xam_table, filename=tsv_out)
        GC.gc()
        # Arrow.write(arrow_out, fully_taxonomy_aware_xam_table)
    end
    
    verbose && println("Taxonomy merging complete!")
    
    # return (tsv_out=tsv_out, arrow_out=arrow_out)
    return tsv_out
end

"""
Classify reads based on taxonomic alignments using individual alignment scoring.

This function takes taxonomy-aware alignment data and performs classification by:
1. Loading the taxonomy-aware alignment data
2. Analyzing alignment score distributions per read
3. Identifying dominant taxonomic assignments
4. Applying conservative taxonomy classification

# Arguments
- `taxonomy_aware_file`: Path to taxonomy-aware alignment file (.tsv.gz or .arrow)
- `min_relative_proportion::Float64=60.0`: Minimum relative proportion threshold for accepting a taxonomic assignment
- `verbose::Bool=true`: Whether to print progress information

# Returns
A DataFrame with taxonomic classification results including individual alignment metrics
"""
function classify_reads_by_taxonomy(;
    taxonomy_aware_file::String,
    min_relative_proportion::Float64 = 60.0,
    verbose::Bool = true
)
    # Load taxonomy-aware data
    verbose && println("Loading taxonomy-aware alignment data...")
    taxonomy_data = if endswith(taxonomy_aware_file, ".arrow")
        if verbose
            @time DataFrames.DataFrame(Arrow.Table(taxonomy_aware_file))
        else
            DataFrames.DataFrame(Arrow.Table(taxonomy_aware_file))
        end
    elseif endswith(taxonomy_aware_file, ".tsv.gz")
        if verbose
            @time Mycelia.read_tsvgz(taxonomy_aware_file)
        else
            Mycelia.read_tsvgz(taxonomy_aware_file)
        end
    else
        error("Unsupported file format. Expected .tsv.gz or .arrow")
    end
    
    # Check for required columns
    required_columns = ["template", "taxid", "alignment_score", "relative_alignment_proportion_percent"]
    missing_columns = setdiff(required_columns, names(taxonomy_data))
    if !isempty(missing_columns)
        error("Missing required columns: $(join(missing_columns, ", "))")
    end
    
    # Filter out rows with missing values in critical columns
    verbose && println("Filtering out rows with missing values in critical columns...")
    filtered_data = DataFrames.filter(
        row -> !ismissing(row.taxid) && !ismissing(row.alignment_score) && !ismissing(row.relative_alignment_proportion_percent),
        taxonomy_data
    )
    
    verbose && println("Retained $(DataFrames.nrow(filtered_data)) rows out of $(DataFrames.nrow(taxonomy_data)) total after filtering missing values")
    
    # Filter alignments that meet the relative proportion threshold
    verbose && println("Filtering alignments by relative proportion threshold (≥$min_relative_proportion%)...")
    high_confidence_alignments = DataFrames.filter(
        row -> row.relative_alignment_proportion_percent >= min_relative_proportion,
        filtered_data
    )
    
    verbose && println("Found $(DataFrames.nrow(high_confidence_alignments)) alignments above threshold out of $(DataFrames.nrow(filtered_data)) total")
    
    # For each template, analyze the taxonomic assignments
    verbose && println("Analyzing taxonomic assignments per read...")
    classification_results = DataFrames.combine(DataFrames.groupby(filtered_data, :template)) do group
        # Sort by relative proportion descending
        sorted = DataFrames.sort(group, :relative_alignment_proportion_percent, rev=true)
        
        # Get top assignment
        top_taxid = sorted[1, :taxid]
        top_proportion = sorted[1, :relative_alignment_proportion_percent]
        top_score = sorted[1, :alignment_score]
        
        # Check if top assignment meets threshold
        passes_threshold = top_proportion >= min_relative_proportion
        
        # Calculate proportion ratio with next best (if it exists)
        proportion_ratio = Inf
        if DataFrames.nrow(sorted) > 1
            second_proportion = sorted[2, :relative_alignment_proportion_percent]
            proportion_ratio = top_proportion / second_proportion
        end
        
        # Count total alignments for this template
        total_alignments = DataFrames.nrow(sorted)
        
        # Store all taxonomic assignments and their proportions
        all_assignments = OrderedCollections.OrderedDict{Int, Float64}()
        for i in 1:DataFrames.nrow(sorted)
            all_assignments[sorted[i, :taxid]] = sorted[i, :relative_alignment_proportion_percent]
        end
        
        return DataFrames.DataFrame(
            top_taxid = top_taxid,
            top_proportion_percent = top_proportion,
            top_alignment_score = top_score,
            passes_threshold = passes_threshold,
            proportion_ratio_to_next_best = proportion_ratio,
            total_alignments = total_alignments,
            all_taxonomic_assignments = [all_assignments]
        )
    end
    
    # Apply conservative taxonomy classification
    verbose && println("Applying conservative taxonomy classification...")
    
    # Prepare data for conservative taxonomy (adapting to expected format)
    conservative_input = DataFrames.DataFrame(
        template = classification_results.template,
        top_taxid = classification_results.top_taxid,
        top_score = classification_results.top_alignment_score,
        ratio_to_next_best_score = classification_results.proportion_ratio_to_next_best,
        additional_taxids = classification_results.all_taxonomic_assignments
    )
    
    final_classification = Mycelia.apply_conservative_taxonomy(conservative_input)
    
    # Merge back with our detailed results
    detailed_classification = DataFrames.leftjoin(
        classification_results,
        final_classification,
        on = :template
    )
    
    verbose && println("Read classification complete!")
    verbose && println("Summary:")
    verbose && println("  - Total reads: $(DataFrames.nrow(detailed_classification))")
    verbose && println("  - Reads passing threshold: $(sum(detailed_classification.passes_threshold))")
    verbose && println("  - Reads with multiple alignments: $(sum(detailed_classification.total_alignments .> 1))")
    
    return detailed_classification
end

# function assign_lowest_rank_to_reads_to_taxon_lineage_table(reads_to_taxon_lineage_table)
#     reads_to_taxon_lineage_table[!, "lowest_rank"] .= 0
#     for (i, row) in enumerate(DataFrames.eachrow(reads_to_taxon_lineage_table))
#         if !ismissing(row["species"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 9
#         elseif !ismissing(row["genus"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 8
#         elseif !ismissing(row["family"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 7
#         elseif !ismissing(row["order"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 6
#         elseif !ismissing(row["class"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 5
#         elseif !ismissing(row["phylum"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 4
#         elseif !ismissing(row["kingdom"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 3
#         elseif !ismissing(row["realm"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 2
#         elseif !ismissing(row["domain"])
#             reads_to_taxon_lineage_table[i, "lowest_rank"] = 1
#         end
#     end
#     return reads_to_taxon_lineage_table
# end

# function classify_xam_with_blast_taxonomies(;xam, accession_to_taxid_table)
#     xam_table_columns_of_interest = [
#         "template",
#         "ismapped",
#         "isprimary",
#         "flag",
#         "reference",
#         "position",
#         "alignment_score",
#     ]
    
#     blast_tax_table_columns_of_interest = [
#         "accession",
#         "taxid",
#     ]
    
#     # 151.408048 seconds (317.18 k allocations: 460.065 MiB, 0.10% compilation time: 20% of which was recompilation)
#     @time xam_table = Mycelia.xam_to_dataframe(xam)

#     taxid_aware_xam_table = DataFrames.leftjoin(
#         DataFrames.select(xam_table, xam_table_columns_of_interest),
#         DataFrames.select(accession_to_taxid_table, blast_tax_table_columns_of_interest),
#         on="reference" => "accession",
#         matchmissing = :notequal
#     )
    
#     template_taxid_score_table = DataFrames.combine(DataFrames.groupby(taxid_aware_xam_table, [:template, :taxid]), :alignment_score => sum => :total_alignment_score)
#     # replace missing (unclassified) with 0 (NCBI taxonomies start at 1, so 0 is a common, but technically non-standard NCBI taxon identifier
#     template_taxid_score_table = DataFrames.coalesce.(template_taxid_score_table, 0)
    
#     # For each template, identify top taxid and calculate score difference
#     results_df = DataFrames.combine(DataFrames.groupby(template_taxid_score_table, :template)) do group
#         # Sort scores in descending order
#         sorted = DataFrames.sort(group, :total_alignment_score, rev=true)

#         # Get top taxid and score
#         top_taxid = sorted[1, :taxid]
#         top_score = sorted[1, :total_alignment_score]

#         # Calculate ratio with next best (if it exists)
#         score_ratio = Inf
#         if DataFrames.nrow(sorted) > 1
#             second_score = sorted[2, :total_alignment_score]
#             score_ratio = top_score/second_score
#         end

#         # Store all additional taxids and their scores (excluding the top one)
#         additional_taxids = OrderedCollections.OrderedDict{Int, Float64}()
#         if DataFrames.nrow(sorted) > 1
#             for i in 2:DataFrames.nrow(sorted)
#                 additional_taxids[sorted[i, :taxid]] = sorted[i, :total_alignment_score]
#             end
#         end

#         # Return a new row with the results
#         return DataFrames.DataFrame(
#             top_taxid = top_taxid,
#             top_score = top_score,
#             ratio_to_next_best_score = score_ratio,
#             additional_taxids = [additional_taxids]  # Wrap in array to make it a single element
#         )
#     end
#     classification_table = Mycelia.apply_conservative_taxonomy(results_df)
#     return classification_table
# end

# function apply_conservative_taxonomy(results_df; ratio_threshold=2.0)
#     # Initialize output columns with default values
#     n_rows = DataFrames.nrow(results_df)
#     final_assignment = copy(results_df.top_taxid)
#     confidence_level = fill("high", n_rows)
    
#     # Collect all sets of taxids that need LCA calculation
#     lca_needed_indices = Int[]
#     competing_taxids_list = Vector{Vector{Int}}()
    
#     # First pass: identify which rows need LCA and prepare the taxid sets
#     for i in 1:n_rows
#         top_taxid = results_df.top_taxid[i]
#         top_score = results_df.top_score[i]
#         ratio = results_df.ratio_to_next_best_score[i]
#         additional_dict = results_df.additional_taxids[i]
        
#         # Skip if no competitors or ratio is high enough
#         if isempty(additional_dict) || ratio >= ratio_threshold
#             continue
#         end
        
#         # Collect taxids that are within the threshold
#         competing_taxids = [top_taxid]
#         for (taxid, score) in additional_dict
#             if top_score / score < ratio_threshold
#                 push!(competing_taxids, taxid)
#             end
#         end
        
#         # If we have multiple competing taxids, add to the batch
#         if length(competing_taxids) > 1
#             push!(lca_needed_indices, i)
#             push!(competing_taxids_list, competing_taxids)
#             confidence_level[i] = "lca"  # Mark as needing LCA
#         end
#     end
    
#     # If we have any rows that need LCA calculation
#     if !isempty(lca_needed_indices)
#         # Batch calculate all LCAs
#         lca_results = Mycelia.batch_taxids2lca(competing_taxids_list)
        
#         # Apply LCA results to the appropriate rows
#         for (idx, lca_idx) in enumerate(lca_needed_indices)
#             final_assignment[lca_idx] = lca_results[idx]
#         end
#     end
    
#     # Create a new dataframe with the original data plus our new columns
#     updated_results = DataFrames.hcat(
#         results_df,
#         DataFrames.DataFrame(
#             final_assignment = final_assignment,
#             confidence_level = confidence_level
#         )
#     )
    
#     return updated_results
# end


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
    setup_taxonkit_taxonomy()
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
        setup_taxonkit_taxonomy()
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
    setup_taxonkit_taxonomy()
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
    setup_taxonkit_taxonomy()
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
    setup_taxonkit_taxonomy()
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
    setup_taxonkit_taxonomy()
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
    setup_taxonkit_taxonomy()

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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_mmseqs_easy_taxonomy(;out_dir, query_fasta, target_database, outfile, force=false)
#     add_bioconda_env("mmseqs2")
#     out_dir = mkpath(joinpath(out_dir, "mmseqs_easy_taxonomy"))
#     outfile = joinpath(out_dir, outfile * ".mmseqs_easy_taxonomy." * basename(target_database) * ".txt")
#     # note I tried adjusting all of the following, and none of them improved overall runtime
#     # in any meaningful way
#     # -s FLOAT                         Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
#     # https://github.com/soedinglab/MMseqs2/issues/577#issuecomment-1191584081
#     # apparently orf-filter 1 speeds up by 50%!
#     # --orf-filter INT                 Prefilter query ORFs with non-selective search
#     #                               Only used during nucleotide-vs-protein classification
#     #                               NOTE: Consider disabling when classifying short reads [0]
#     # --lca-mode INT                   LCA Mode 1: single search LCA , 2/3: approximate 2bLCA, 4: top hit [3]
#     # --lca-search BOOL                Efficient search for LCA candidates [0]
#     # ^ this looks like it actually runs @ 1 with s=1.0 & --orf-filter=1
    
#     # 112 days to process 600 samples at this rate....
#     # 278 minutes or 4.5 hours for a single sample classification!!
#     # 16688.050696 seconds (1.43 M allocations: 80.966 MiB, 0.02% gc time, 0.00% compilation time)
#     # this is for default parameters
#     # lowering sensitivity and taking LCA
#     # 16590.725343 seconds (1.06 M allocations: 53.487 MiB, 0.00% compilation time)
#     # took just as long!
#     # difference was only 10 minutes
#     # 15903.218456 seconds (969.92 k allocations: 48.624 MiB, 0.01% gc time)
#     # use default parameters
    
#     if force || (!force && !isfile(outfile))
#         cmd = 
#         `$(CONDA_RUNNER) run --no-capture-output -n mmseqs2 mmseqs
#          easy-taxonomy
#          $(query_fasta)
#          $(target_database)
#          $(outfile)
#          $(joinpath(out_dir, "tmp"))
#         `
#         @time run(pipeline(cmd))
#     else
#         @info "target outfile $(outfile) already exists, remove it or set force=true to re-generate"
#     end
#     return outfile
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse an MMseqs2 easy-taxonomy tophit report into a structured DataFrame.

# Arguments
- `tophit_report::String`: Path to the MMseqs2 easy-taxonomy tophit report file (tab-delimited)

# Returns
- `DataFrame`: A DataFrame with columns:
  - `target_id`: Target sequence identifier
  - `number of sequences aligning to target`: Count of aligned sequences
  - `unique coverage of target`: Ratio of uniqueAlignedResidues to targetLength
  - `Target coverage`: Ratio of alignedResidues to targetLength
  - `Average sequence identity`: Mean sequence identity
  - `taxon_id`: Taxonomic identifier
  - `taxon_rank`: Taxonomic rank
  - `taxon_name`: Species name and lineage
"""
function parse_mmseqs_easy_taxonomy_tophit_report(tophit_report)
    data, header = uCSV.read(tophit_report, delim='\t')
    # tophit_report
    # (1) Target identifier 
    # (2) Number of sequences aligning to target
    # (3) Unique coverage of target uniqueAlignedResidues / targetLength
    # (4) Target coverage alignedResidues / targetLength
    # (5) Average sequence identity
    # (6) Taxonomical information identifier, species, lineage
    header = [
        "target_id",
        "number of sequences aligning to target",
        "unique coverage of target (uniqueAlignedResidues / targetLength)",
        "Target coverage (alignedResidues / targetLength)",
        "Average sequence identity",
        "taxon_id",
        "taxon_rank",
        "taxon_name"
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse the taxonomic Last Common Ancestor (LCA) TSV output from MMseqs2's easy-taxonomy workflow.

# Arguments
- `lca_tsv`: Path to the TSV file containing MMseqs2 taxonomy results

# Returns
DataFrame with columns:
- `contig_id`: Sequence identifier
- `taxon_id`: NCBI taxonomy identifier 
- `taxon_rank`: Taxonomic rank (e.g. species, genus)
- `taxon_name`: Scientific name
- `fragments_retained`: Number of fragments kept
- `fragments_taxonomically_assigned`: Number of fragments with taxonomy
- `fragments_in_agreement_with_assignment`: Fragments matching contig taxonomy
- `support -log(E-value)`: Statistical support score
"""
function parse_mmseqs_easy_taxonomy_lca_tsv(lca_tsv)
    data, header = uCSV.read(lca_tsv, delim='\t')
    # contig
    # (1) a single taxonomy numeric identifier
    # (2) a taxonomic rank column
    # (3) taxonomic name column
    # fragments retained
    # fragments taxonomically assigned
    # fragments in agreement with the contig label (i.e. same taxid or have it as an ancestor)
    # the support received -log(E-value)
    header = [
        "contig_id",
        "taxon_id",
        "taxon_rank",
        "taxon_name",
        "fragments_retained",
        "fragments_taxonomically_assigned",
        "fragments_in_agreement_with_assignment",
        "support -log(E-value)"
    ]
    return DataFrames.DataFrame(data, header)
end