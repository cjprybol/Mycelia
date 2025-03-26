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
    joint_table = DataFrames.DataFrame()
    ProgressMeter.@showprogress for taxid in taxids
        cmd1 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon $(taxid) --as-json-lines`
        cmd2 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli dataformat tsv taxonomy --template tax-summary`
        io = open(pipeline(cmd1, cmd2))
        append!(joint_table, DataFrames.DataFrame(uCSV.read(io, delim='\t', header=1)), promote=true)
    end
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
    data, header = uCSV.read(open(pipeline(cmd)), delim='\t', header=false)
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
            superkingdom_taxid = haskey(lineage_ranks, "superkingdom") ? lineage_ranks["superkingdom"].taxid : missing,
            superkingdom = haskey(lineage_ranks, "superkingdom") ? lineage_ranks["superkingdom"].lineage : missing,
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