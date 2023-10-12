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
    for rank inranks
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
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{$(shorthand)}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end


function list_superkingdoms()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "superkingdom"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{k}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_kingdoms()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "kingdom"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{K}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_phylums()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "phylum"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{p}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_classes()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "class"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{c}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_orders()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "order"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{o}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_families()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "family"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{f}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_genera()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "genus"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{g}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end

function list_species()
    p = pipeline(
        `conda run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `conda run --no-capture-output -n taxonkit taxonkit filter --equal-to "species"`,
        `conda run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format '{s}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "name"]
    DataFrames.DataFrame(data, header)
end