# ```
# https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/

# GFF3 has 9 required fields, though not all are utilized (either blank or a default value of ‘.’).

#     Sequence ID
#     Source
#         Describes the algorithm or the procedure that generated this feature. Typically Genescane or Genebank, respectively.
#     Feature Type
#         Describes what the feature is (mRNA, domain, exon, etc.).
#         These terms are constrained to the [Sequence Ontology terms](http://www.sequenceontology.org/).
#     Feature Start
#     Feature End
#     Score
#         Typically E-values for sequence similarity and P-values for predictions.
#     Strand
#     Phase
#         Indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
#     Atributes
#         A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent . You can see the full list [here](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

# ```

# """
#     create_chromosome_genedata_table(chromosome)

# Take a chromosome from GenomicAnnotations.jl in GFF (and possibly genbank)
# and return the formed dataframe.
# """
# function create_chromosome_genedata_table(chromosome)
#     # genedata is already provided as a dataframe with all of the Attributes as columns
#     table = copy(chromosome.genedata)
    
#     # the rest of these need to be created
#     # I'm inserting them directly into their correct locations
#     DataFrames.insertcols!(table, 1, "sequence-id" => fill(chromosome.name, DataFrames.nrow(table)))
    
#     DataFrames.insertcols!(table, 3, "feature" => GenomicAnnotations.feature.(chromosome.genes))

#     loci = getproperty.(GenomicAnnotations.locus.(chromosome.genes), :position)
#     DataFrames.insertcols!(table, 4, "start" => first.(loci))
#     DataFrames.insertcols!(table, 5, "stop" => last.(loci))
    
#     DataFrames.insertcols!(table, 7, "strand" => .!GenomicAnnotations.iscomplement.(chromosome.genes))        
    
#     return table
# end

# function create_chromosome_genedata_table(chromosome)
#     table = chromosome.genedata
#     table[!, "chromosome"] .= chromosome.name
#     table[!, "feature"] = GenomicAnnotations.feature.(chromosome.genes)
#     table[!, "strand"] = .!GenomicAnnotations.iscomplement.(chromosome.genes)
#     table[!, "locus"] = getproperty.(GenomicAnnotations.locus.(chromosome.genes), :position)
#     table[!, "start"] = first.(chromosome.genedata[!, "locus"])
#     table[!, "stop"] = last.(chromosome.genedata[!, "locus"])
#     return table
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Update GFF annotations with protein descriptions from MMseqs2 search results.

# Arguments
- `gff_file::String`: Path to input GFF3 format file
- `mmseqs_file::String`: Path to MMseqs2 easy-search output file

# Returns
- `DataFrame`: Modified GFF table with updated attribute columns containing protein descriptions

# Details
Takes sequence matches from MMseqs2 and adds their descriptions as 'label' and 'product' 
attributes in the GFF file. Only considers top hits from MMseqs2 results. Preserves existing 
GFF attributes while prepending new annotations.
"""
function update_gff_with_mmseqs(gff_file, mmseqs_file)
    mmseqs_results = read_mmseqs_easy_search(mmseqs_file)
    top_hits = DataFrames.combine(DataFrames.groupby(mmseqs_results, "query"), group -> group[Base.argmax(group.bits), :])
    # id_to_product = Dict{String, String}()
    # for row in DataFrames.eachrow(top_hits)
    #     id = Dict(a => b for (a, b) in split.(split(last(split(row["qheader"], " # ")), ';'), '='))["ID"]
    #     product = replace(row["theader"], " " => "__")
    #     id_to_product[id] = product
    # end
    id_to_product = Dict{String, String}(row["query"] => replace(row["theader"], " " => "__") for row in DataFrames.eachrow(top_hits))

    gff_table = Mycelia.read_gff(gff_file)
    for (i, row) in enumerate(DataFrames.eachrow(gff_table))
        id = Dict(a => b for (a,b) in split.(filter(!isempty, split(row["attributes"], ';')), '='))["ID"]
        product = get(id_to_product, id, "")
        gff_table[i, "attributes"] = "label=\"$(product)\";product=\"$(product)\";" * row["attributes"]
    end
    return gff_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Opens a GFF (General Feature Format) file for reading.

# Arguments
- `path::String`: Path to GFF file. Can be:
    - Local file path
    - HTTP/FTP URL (FTP URLs are automatically converted to HTTP)
    - Gzipped file (automatically decompressed)

# Returns
- `IO`: An IO stream ready to read the GFF content
"""
function open_gff(path::String)
    if isfile(path)
        io = open(path)
    elseif occursin(r"^ftp", path) || occursin(r"^http", path)
        path = replace(path, r"^ftp:" => "http:")
        io = IOBuffer(HTTP.get(path).body)
    else
        error("unable to locate file $path")
    end
    path_base = basename(path)
    if occursin(r"\.gz$", path_base)
        io = CodecZlib.GzipDecompressorStream(io)
    end
    return io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Reads a GFF (General Feature Format) file and parses it into a DataFrame.

# Arguments
- `gff::AbstractString`: Path to the GFF file

# Returns
- `DataFrame`: A DataFrame containing the parsed GFF data with standard columns:
  seqid, source, type, start, end, score, strand, phase, and attributes
"""
function read_gff(gff::AbstractString)
    return read_gff(open_gff(gff))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Read a GFF (General Feature Format) file into a DataFrame.

# Arguments
- `gff_io`: An IO stream containing GFF formatted data

# Returns
- `DataFrame`: A DataFrame with standard GFF columns:
  - seqid: sequence identifier
  - source: feature source
  - type: feature type
  - start: start position (1-based)
  - end: end position
  - score: numeric score
  - strand: strand (+, -, or .)
  - phase: phase (0, 1, 2 or .)
  - attributes: semicolon-separated key-value pairs
"""
function read_gff(gff_io)
    data, header = uCSV.read(gff_io, delim='\t', header=0, comment='#')
    header = [
        "#seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
    types=[
        String,
        String,
        String,
        Int,
        Int,
        Int,
        String,
        Int,
        String
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Takes a GFF (General Feature Format) DataFrame and expands the attributes column into separate columns.

# Arguments
- `gff_df::DataFrame`: A DataFrame containing GFF data with an 'attributes' column formatted as key-value pairs
  separated by semicolons (e.g., "ID=gene1;Name=BRCA1;Type=gene")

# Returns
- `DataFrame`: The input DataFrame with additional columns for each attribute key found in the 'attributes' column
"""
function split_gff_attributes_into_columns(gff_df)
    # 1. Extract unique keys from the attribute column
    all_keys = Set{String}()
    for row in DataFrames.eachrow(gff_df)
        attributes = split(row["attributes"], ';')
        for attribute in attributes
            if !isempty(attribute)
                key = split(attribute, '=')[1]
                push!(all_keys, key)
            end
        end
    end

  # 2. Create new columns for each key
    for key in all_keys
        gff_df[!, key] = Vector{Union{Missing, String}}(missing, size(gff_df, 1))
    end

    # 3. Populate the new columns with values
    for (i, row) in enumerate(DataFrames.eachrow(gff_df))
        attributes = split(row["attributes"], ';')
        for attribute in attributes
            if !isempty(attribute)
            key_value = split(attribute, '=')
                if length(key_value) == 2
                    key, value = key_value
                    gff_df[i, key] = value
                end
            end
        end
    end
    return gff_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write GFF (General Feature Format) data to a tab-delimited file.

# Arguments
- `gff`: DataFrame/Table containing GFF formatted data
- `outfile`: String path where the output file should be written

# Returns
- `String`: Path to the written output file
"""
function write_gff(;gff, outfile)
    uCSV.write(outfile, gff, delim='\t')
    return outfile
end

# TODO: switch to using GenomicAnnotations if GFF3 package isn't updated
# PR -> https://github.com/BioJulia/GFF3.jl/pull/12
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

Retrieve GFF3 formatted genomic feature data from NCBI or direct FTP source.

# Arguments
- `db::String`: NCBI database to query ("nuccore" for DNA or "protein" for protein sequences)
- `accession::String`: NCBI accession number
- `ftp::String`: Direct FTP URL to GFF3 file (typically gzipped)

# Returns
- `IO`: IOBuffer containing uncompressed GFF3 data
"""
function get_gff(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=gff3&id=$(accession)"
        return IOBuffer(HTTP.get(url).body)
    elseif !isempty(ftp)
        return CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body))
    else
        @error "invalid call"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert FASTA sequence and GFF annotation files to GenBank format using EMBOSS seqret.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequence data
- `gff::String`: Path to input GFF file containing genomic features
- `genbank::String`: Path for output GenBank file

# Details
Requires EMBOSS toolkit (installed via Bioconda). The function will:
1. Create necessary output directories
2. Run seqret to combine sequence and features
3. Generate a GenBank format file at the specified location
"""
function fasta_and_gff_to_genbank(;fasta, gff, genbank=gff * ".genbank")
    add_bioconda_env("emboss")
    # https://www.insdc.org/submitting-standards/feature-table/
    genbank_directory = dirname(genbank)
    genbank_basename = basename(genbank)
    genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
    # genbank_extension = "genbank"
    # https://bioinformatics.stackexchange.com/a/11140
    # https://www.biostars.org/p/72220/#72272
    # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
#     -osname
    # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
    run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
    return genbank
end

# function gff_to_genbank(gff, genbank)
#     # https://www.insdc.org/submitting-standards/feature-table/
#     genbank_directory = dirname(genbank)
#     genbank_basename = basename(genbank)
#     genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
#     # genbank_extension = "genbank"
#     # https://bioinformatics.stackexchange.com/a/11140
#     # https://www.biostars.org/p/72220/#72272
#     # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
# #     -osname
#     # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
#     run(`seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
#     # return genbank
# end

# function genbank_to_fasta_and_gff(genbank_file)

# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Opens and parses a GenBank format file containing genomic sequence and annotation data.

# Arguments
- `genbank_file::AbstractString`: Path to the GenBank (.gb or .gbk) file

# Returns
- `Vector{GenomicAnnotations.Chromosome}`: Vector containing parsed chromosome data
"""
function open_genbank(genbank_file)
    return GenomicAnnotations.readgbk(genbank_file)
end