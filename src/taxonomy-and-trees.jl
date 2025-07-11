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