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

Update GFF annotations using a pre-loaded DataFrame of MMseqs2 results.

# Keyword Arguments
- `gff_file::String`: Path to input GFF3 format file.
- `mmseqs_results::DataFrames.DataFrame`: DataFrame of MMseqs2 easy-search results.
- `id_map::Union{Dict{String, String}, Nothing}=nothing`: Optional dictionary mapping original GFF IDs to normalized IDs.
- `no_hit_label::String="No Hit"`: Label to apply to annotations with no matching result in the mmseqs_results

# Returns
- `DataFrame`: Modified GFF table with updated attribute columns.

# Details
This function takes sequence matches from a `DataFrame` of MMseqs2 results and adds their
descriptions as 'label' and 'product' attributes in the GFF file.

The 'product' attribute is always set to the full 'theader' value. The 'label' is a
normalized version derived by parsing the cluster name from the 'theader' field.

If an `id_map` is provided, it links GFF protein IDs to normalized IDs in the `mmseqs_results`.
"""
function update_gff_with_mmseqs(
    ;
    gff_file::String,
    mmseqs_results::DataFrames.DataFrame,
    id_map::Union{Dict{String, String}, Nothing}=nothing,
    no_hit_label::String="No Hit",
    extract_cluster_name::Bool=false
)

    function normalize_header_label(theader::String)
        return replace(theader, ' ' => "__")
    end

    # This inner function extracts the cluster name and normalizes it for the 'label'
    function get_normalized_label(theader::String)
        if !extract_cluster_name
            return normalize_header_label(theader)
        end

        # Extract the part of the header before the sequence count (e.g., "n=5")
        match_result = Base.match(r"^(\S+\s+.*?)\s+n=\d+", theader)
        label_base = match_result !== nothing ? Base.String(match_result[1]) : theader

        # Split the base, take all words from the second onward, and join with an underscore
        parts = Base.split(label_base, ' ')
        # This check handles cases where there might be fewer than two words
        return Base.length(parts) > 1 ? normalize_header_label(Base.join(parts[2:end], "_")) : ""
    end

    # Create a dictionary mapping each query to its full 'theader' (for product)
    # and the generated normalized label
    query_to_info = Dict{String, Tuple{String, String}}(
        row.query => (normalize_header_label(row.theader), get_normalized_label(row.theader)) for row in DataFrames.eachrow(mmseqs_results)
    )

    gff_table = Mycelia.read_gff(gff_file)
    for i in 1:DataFrames.nrow(gff_table)
        row = gff_table[i, :]
        attributes_dict = Dict(a => b for (a,b) in Base.split.(Base.filter(!Base.isempty, Base.split(row.attributes, ';')), '='))
        original_id = Base.get(attributes_dict, "ID", "")

        # Use the id_map if provided, otherwise use the original ID for lookup
        lookup_key = if id_map !== nothing
            Base.get(id_map, original_id, nothing)
        else
            original_id
        end

        info = if lookup_key !== nothing
            Base.get(query_to_info, lookup_key, nothing)
        else
            nothing
        end
        
        if info !== nothing
            product, normalized_label = info
        else
            product = normalize_header_label(no_hit_label)
            normalized_label = product
        end
        # Prepend the new label and product attributes to the existing attributes
        gff_table[i, "attributes"] = "label=\"$(normalized_label)\";product=\"$(product)\";" * row.attributes
    end
    return gff_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Update GFF annotations with protein descriptions from an MMseqs2 search results file.

# Arguments
- `gff_file::String`: Path to input GFF3 format file
- `mmseqs_file::String`: Path to MMseqs2 easy-search output file

# Keyword Arguments
- `id_map::Union{Dict{String, String}, Nothing}=nothing`: Optional dictionary mapping original GFF IDs to normalized IDs.
- `extract_cluster_name::Bool=false`: If true, parse the `theader` column to extract only the cluster name.

# Returns
- `DataFrame`: Modified GFF table with updated attribute columns containing protein descriptions
"""
function update_gff_with_mmseqs(
    gff_file::String,
    mmseqs_file::String;
    id_map::Union{Dict{String, String}, Nothing}=nothing,
    extract_cluster_name::Bool=false
)
    mmseqs_results = Mycelia.read_mmseqs_easy_search(mmseqs_file)
    # Pass keywords through to the next call
    return update_gff_with_mmseqs(
        gff_file=gff_file,
        mmseqs_results=mmseqs_results,
        id_map=id_map,
        extract_cluster_name=extract_cluster_name
    )
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
    append_suffix_to_gff_id_final(gff_df::DataFrames.DataFrame; suffix::Union{String, Nothing}=nothing)

Correctly parses a GFF attributes column by splitting it into key-value pairs
before appending a suffix to the primary ID tag. This is robust against ID-like
substrings appearing inside quoted values.

# Arguments
- `gff_df::DataFrames.DataFrame`: DataFrame with GFF data, must have an "attributes" column.

# Keyword Arguments
- `suffix::Union{String, Nothing}=nothing`: Suffix to append. Defaults to row number if nothing.

# Returns
- `DataFrames.DataFrame`: A new DataFrame with correctly modified IDs.
"""
function append_suffix_to_gff_id(gff_df::DataFrames.DataFrame; suffix::Union{String, Nothing}=nothing)
    modified_df = deepcopy(gff_df)

    if !("attributes" in DataFrames.names(modified_df))
        error("Input DataFrame must contain a column named 'attributes'.")
    end

    for i in 1:DataFrames.nrow(modified_df)
        attributes_str = modified_df[i, :attributes]
        
        # Step 1: Split the string into individual "key=value" parts. This is the key fix.
        parts = split(attributes_str, ';')
        
        id_found_and_modified = false
        # Step 2: Iterate through the parts to find the correct ID tag.
        for (j, part) in enumerate(parts)
            # Use strip() to handle any potential whitespace around the key=value pair
            # Use startswith() to reliably find the primary ID tag for the feature.
            if startswith(strip(part), "ID=")
                
                # Determine the suffix to use
                current_suffix = suffix === nothing ? "_$(i)" : suffix
                
                # Step 3: Modify only this part of the array
                parts[j] = part * current_suffix
                
                id_found_and_modified = true
                
                # Once we find and modify the ID, we can stop searching this row's attributes.
                break 
            end
        end
        
        # Step 4: If we made a change, join the parts back together.
        if id_found_and_modified
            modified_df[i, :attributes] = join(parts, ';')
        else
            # Optional warning if a row that should have an ID doesn't
            if "type" in DataFrames.names(modified_df)
                row_type = modified_df[i, :type]
                if row_type in ["gene", "Gene", "CDS", "mRNA"]
                    @warn "Row $i (type: $row_type) did not contain a parsable 'ID=' tag."
                end
            end
        end
    end

    return modified_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert FASTA sequence and GFF annotation files to GenBank format using EMBOSS seqret.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequence data. Can be gzipped.
- `gff::String`: Path to input GFF file containing genomic features.
- `genbank::String`: Path for output GenBank file.

# Details
Requires EMBOSS toolkit (installed via Bioconda). The function will:
1. Create necessary output directories.
2. If the FASTA file is gzipped, create a temporary uncompressed version.
3. Run `seqret` to combine sequence and features.
4. Generate a GenBank format file at the specified location.
5. Clean up any temporary files.
"""
function fasta_and_gff_to_genbank(;fasta, gff, genbank=gff * ".genbank")
    add_bioconda_env("emboss")
    # https://www.insdc.org/submitting-standards/feature-table/
    genbank_directory = dirname(genbank)
    genbank_basename = basename(genbank)
    genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
    
    # Handle gzipped fasta
    if endswith(fasta, ".gz")
        temp_fasta = Mycelia.write_fasta(records = collect(Mycelia.open_fastx(fasta)))
        try
            run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret -sequence $(temp_fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
        finally
            rm(temp_fasta)
        end
    else
        # https://bioinformatics.stackexchange.com/a/11140
        # https://www.biostars.org/p/72220/#72272
        # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
        # -osname
        # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
        run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
    end
    
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

"""
    read_bed(file_path::String)

Reads a BED file from `file_path` into a `DataFrame`.

This function reads a standard BED file, which is a tab-delimited format without a header.
It assigns the standard BED column names and converts the appropriate columns to their
correct integer types for numerical analysis.

# Arguments
- `file_path::String`: The path to the BED file.

# Returns
- `DataFrames.DataFrame`: A DataFrame containing the BED file data with correctly typed columns.
"""
function read_bed(file_path::String)
    # Define column names for all possible 12 columns in a BED file.
    bed_column_names = [
        :chrom, :chromStart, :chromEnd, :name, :score, :strand,
        :thickStart, :thickEnd, :itemRgb, :blockCount, :blockSizes, :blockStarts
    ]

    # Read the raw data. We let CSV.jl detect types initially.
    # This is more efficient than reading everything as a string first.
    df = CSV.read(file_path, DataFrames.DataFrame, header=false, delim='\t', comment="#")

    # Rename the columns of the DataFrame based on how many were read.
    DataFrames.rename!(df, bed_column_names[1:size(df, 2)])

    # List of columns that should be converted to integers.
    # We will check for their existence before attempting conversion.
    int_columns = [:chromStart, :chromEnd, :score, :thickStart, :thickEnd, :blockCount]

    # Use a single `transform!` call for all type conversions.
    # This is more idiomatic and robust. It modifies the DataFrame in-place.
    # `ByRow` ensures the parsing function is applied to each element.
    # `AsTable` and `=>` are used to select and transform multiple columns.
    transform_pairs = []
    for col_name in int_columns
        if DataFrames.hasproperty(df, col_name)
            # Ensure the column is not already an integer type
            if eltype(df[!, col_name]) != Int
                push!(transform_pairs, col_name => DataFrames.ByRow(x -> parse(Int, string(x))) => col_name)
            end
        end
    end

    if !isempty(transform_pairs)
        DataFrames.transform!(df, transform_pairs...)
    end

    return df
end

"""
    bed_to_gff(bed_df::DataFrames.DataFrame; source::String="bed_to_gff", feature_type::String="feature")

Converts a BED `DataFrame` to a GFF3-compatible `DataFrame`.

This function creates a new DataFrame with columns corresponding to the GFF3 specification.
It converts coordinates from BED's 0-based system to GFF's 1-based system.
The 'attributes' column is constructed using the 'name' column from the BED data if available.

# Arguments
- `bed_df::DataFrames.DataFrame`: A DataFrame with BED data, typically from `read_bed`.
- `source::String`: The source value to use in the second column of the GFF DataFrame.
- `feature_type::String`: The feature type to use in the third column of the GFF DataFrame.

# Returns
- `DataFrames.DataFrame`: A new DataFrame representing the data in GFF format.
"""
function bed_to_gff(bed_df::DataFrames.DataFrame; source::String="bed_to_gff", feature_type::String="feature")
    # Initialize the GFF DataFrame with required columns, using Symbols for special names
    gff_df = DataFrames.DataFrame(
        Symbol("#seqid") => bed_df.chrom,
        :source => source,
        :type => feature_type,
        :start => bed_df.chromStart .+ 1, # Convert from 0-based (BED) to 1-based (GFF)
        :end => bed_df.chromEnd,          # Use :end as the column name for the stop coordinate
        :score => ".",
        :strand => ".",
        :phase => ".",
        :attributes => ""
    )

    # Populate optional columns if they exist in the source BED DataFrame
    if DataFrames.hasproperty(bed_df, :score)
        # Ensure score is a string or "."
        gff_df.score = [s === missing ? "." : string(s) for s in bed_df.score]
    end

    if DataFrames.hasproperty(bed_df, :strand)
        # Ensure strand is a string or "."
        gff_df.strand = [s === missing ? "." : s for s in bed_df.strand]
    end

    # Construct the attributes column (9th column)
    # GFF attributes are a semicolon-separated list of tag=value pairs.
    if DataFrames.hasproperty(bed_df, :name)
        # Use the 'name' column for the feature ID if it exists
        gff_df.attributes = "ID=" .* bed_df.name
    else
        # Otherwise, create a unique ID from the feature's coordinates
        gff_df.attributes = "ID=" .* gff_df[!, Symbol("#seqid")] .* ":" .* string.(gff_df.start) .* "-" .* string.(gff_df[!, :end])
    end
    
    return gff_df
end
