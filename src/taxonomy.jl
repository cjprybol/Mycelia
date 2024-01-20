function list_full_taxonomy()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --add-prefix --fill-miss-rank --show-lineage-taxids --format '{k}\;{K}\;{p}\;{c}\;{o}\;{f}\;{g}\;{s}'`
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

function list_toplevel()
    return DataFrames.DataFrame(taxid=[0, 1], name=["unclassified", "root"])
end
    

"""
- top
- superkingdom/domain
- kingdom
- phylum
- class
- order
- family
- genus
- species
"""
function list_rank(rank)
    if rank == "top"
        return list_toplevel()
    else
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
            `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
            `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "$(rank)"`,
            `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format "{$(shorthand)}"`
        )
        data, header = uCSV.read(open(p), delim='\t')
        header = ["taxid", "name"]
        return DataFrames.DataFrame(data, header)
    end
end


function list_superkingdoms()
    return list_rank("superkingdom")
end

function list_kingdoms()
    return list_rank("kingdom")
end

function list_phylums()
    return list_rank("phylum")
end

function list_classes()
    return list_rank("class")
end

function list_orders()
    return list_rank("order")
end

function list_families()
    return list_rank("family")
end

function list_genera()
    return list_rank("genus")
end

function list_species()
    return list_rank("species")
end

# function list_subtaxa(taxid)
#     return parse.(Int, filter(!isempty, strip.(readlines(`conda run --no-capture-output -n taxonkit taxonkit list --ids 10239`))))
# end

function list_subtaxa(taxid)
    return parse.(Int, filter(!isempty, strip.(readlines(`conda run --no-capture-output -n taxonkit taxonkit list --ids $(taxid)`))))
end

function name2taxid(name)
    p = pipeline(`echo $(name)`, `conda run --no-capture-output -n taxonkit taxonkit name2taxid --show-rank`)
    data, header = uCSV.read(open(p), delim='\t')
    header = ["name", "taxid", "rank"]
    return DataFrames.DataFrame(data, header)
end

function taxids2lineage_name_and_rank(taxids::AbstractVector{Int})
    f = tempname()
    open(f, "w") do io
        for taxid in taxids
            println(io, taxid)
        end
    end
    data, header = uCSV.read(open(pipeline(`conda run --live-stream -n taxonkit taxonkit lineage --show-name --show-rank $(f)`)), delim='\t', header=false)
    rm(f)
    header = ["taxid", "lineage", "tax_name", "tax_rank"]
    return DataFrames.DataFrame(data, header)
end

# function lca(taxids)
# end


"""
Downloads and unpacks the desired .tar.gz prebuilt kraken index

Go to https://benlangmead.github.io/aws-indexes/k2 and identify the appropriate ".tar.gz" url download
"""
function download_kraken_index(;url, directory="$(homedir())/workspace/kraken")
    @assert occursin(r"\.tar\.gz$", url)
    filename = last(split(url, '/'))
    output_path = joinpath(directory, filename)
    # @show url
    # @show filename
    # @show output_path
    mkpath(directory)
    extracted_directory = replace(basename(output_path), ".tar.gz" => "")
    if !isdir(extracted_directory)
        mkpath(extracted_directory)
    end
    if isempty(readdir(extracted_directory))
        # download the file only if needed
        if !isfile(output_path)
            download(url, output_path)
            # run(`wget -P $(directory) $(url)`)
        end
        run(`tar -xvzf $(output_path) -C $(extracted_directory)`)
    end
    return extracted_directory
end



