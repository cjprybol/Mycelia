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

function list_rank(rank)
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