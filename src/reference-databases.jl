const NCBI_SUPERKINGDOMS = Dict(
    "Bacteria" => 2,
    "Archaea" => 2157,
    "Eukaryota" => 2759,
    "Viruses" => 10239,
    "Other sequences (-)" => 28384,
    "Unclassified (N/A)" => 12908
)

const UN_LANGUAGES = ["ar", "zh", "en", "fr", "ru", "es"]
const UN_PAIRS = [
    "ar-en", "ar-es", "ar-fr", "ar-ru", "ar-zh",
    "en-es", "en-fr", "en-ru", "en-zh",
    "es-fr", "es-ru", "es-zh",
    "fr-ru", "fr-zh",
    "ru-zh"
]
const UN_CORPUS_BASE_URL = "https://www.un.org/dgacm/sites/www.un.org.dgacm/files/files/UNCORPUS"
const UN_CORPUS_USER_AGENT = "Mycelia/UNCorpus"

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BLAST database (or a subset of it) to FASTA format, compressed with gzip.

If `outfile` already exists, is non-empty, and `force == false`, this function is a
no-op and simply returns the existing output path.

# Arguments
- `blastdb::AbstractString`: Name of the BLAST database to convert (e.g. "nr", "nt").
- `entries::AbstractVector{<:AbstractString}`: Optional list of specific sequence IDs
  to extract. Mutually exclusive with `taxids`.
- `taxids::AbstractVector{<:Integer}`: Optional list of taxonomic IDs to extract.
  Mutually exclusive with `entries`.
- `outfile::AbstractString`: Path for the output FASTA (gzip) file. If empty, a
  filename of the form `<blastdb>.<YYYY-MM-DD>.<ext>.gz` is inferred from the
  BLAST DB metadata, where `<YYYY-MM-DD>` is the date when the database was last
  updated and `<ext>` is `.fna` for nucleotide DBs and `.faa` for protein DBs.
- `force::Bool`: If `true`, regenerate the FASTA even if `outfile` already exists
  and is non-empty. If `false`, the function is a no-op in that case.
- `max_cores::Integer`: Maximum number of threads to use for `pigz` compression.

# Returns
- `outfile::String`: Path to the generated (or existing) FASTA file.
"""
function blastdb_to_fasta(; blastdb::AbstractString,
        entries::AbstractVector{<:AbstractString} = String[],
        taxids::AbstractVector{<:Integer} = Int[],
        outfile::AbstractString = "",
        force::Bool = true,
        max_cores::Integer = 16)

    # Validate mutually exclusive options
    if !isempty(entries) && !isempty(taxids)
        error("Can only specify entries OR taxids, not both")
    end

    # Infer output filename if not provided
    if isempty(outfile)
        blastdb_metadata = get_blastdb_metadata(blastdb = blastdb)

        dbtype = blastdb_metadata["dbtype"]
        if dbtype == "Nucleotide"
            extension = ".fna.gz"
        elseif dbtype == "Protein"
            extension = ".faa.gz"
        else
            error("Unsupported BLAST dbtype: $dbtype")
        end

        last_updated_raw = blastdb_metadata["last-updated"]
        last_updated_dt = Dates.DateTime(last_updated_raw)  # adjust if different format is needed
        last_update_date = Dates.format(last_updated_dt, "yyyy-mm-dd")

        outfile = "$(blastdb).$(last_update_date)" * extension
    end

    # No-op if output already exists and we are not forcing regeneration
    if isfile(outfile) && (filesize(outfile) > 0) && !force
        return String(outfile)
    end

    temp_file = ""
    blast_cmd = nothing

    # Build blastdbcmd command, optionally writing temp files for entries/taxids
    if isempty(entries) && isempty(taxids)
        blast_cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -entry all -outfmt %f`
    elseif !isempty(entries)
        temp_file = tempname() * ".entries.txt"
        open(temp_file, "w") do io
            for entry in entries
                println(io, entry)
            end
        end
        blast_cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -entry_batch $(temp_file) -outfmt %f`
    else
        # !isempty(taxids)
        temp_file = tempname() * ".taxids.txt"
        open(temp_file, "w") do io
            for taxid in taxids
                println(io, taxid)
            end
        end
        blast_cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -taxidlist $(temp_file) -outfmt %f`
    end

    @assert blast_cmd !== nothing "Internal error: blast_cmd was not defined"

    # Ensure cores is at least 1 and not more than Threads.nthreads()
    cores = max(1, min(Int(max_cores), Threads.nthreads()))

    # Create output directory if needed
    outdir = dirname(outfile)
    if !isempty(outdir) && outdir != "."
        mkpath(outdir)
    end

    # Build compression pipeline
    p = pipeline(blast_cmd, `$(CONDA_RUNNER) run --live-stream -n pigz pigz -c -p $(cores)`)

    # Only set up envs if we're actually going to run the commands
    Mycelia.add_bioconda_env("blast")
    Mycelia.add_bioconda_env("pigz")

    try
        run(pipeline(p, outfile))
    finally
        # Robust cleanup of temporary file, even on error
        if !isempty(temp_file) && isfile(temp_file)
            rm(temp_file; force = true)
        end
    end

    return String(outfile)
end

# not sure how to use this - doesn't seem super helpful?
# function get_blastdb_dups(;blastdb)
#     Mycelia.add_bioconda_env("blast")
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -get_dups -entry all`)
# end

function get_blastdb_info(; blastdb)
    Mycelia.add_bioconda_env("blast")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -info`)
end

function get_blastdb_metadata(; blastdb)
    Mycelia.add_bioconda_env("blast")
    # need to read this into a json object and return it
    JSON.parse(read(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -metadata`, String))
end

function get_blastdb_tax_info(; blastdb, entries = String[], taxids = Int[])
    Mycelia.add_bioconda_env("blast")
    if !isempty(entries) && !isempty(taxids)
        error("Can only specify entries OR taxids, not both")
    end
    # need to read this into a table and return it
    # %T means taxid
    # %L means common taxonomic name
    # %S means scientific name
    # %K means taxonomic super kingdom
    # %B means BLAST name
    # %n means num of seqs
    temp_file = ""
    col_names = ["taxid", "scientific name", "common taxonomic name",
        "taxonomic super kingdom", "BLAST name", "num of seqs"]
    if isempty(entries) && isempty(taxids)
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -tax_info -outfmt '%T	%S	%L	%K	%B	%n'`
    elseif !isempty(entries)
        temp_file = tempname() * ".entries.txt"
        open(temp_file, "w") do file
            for entry in entries
                println(file, entry)
            end
        end
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -tax_info -outfmt '%T	%S	%L	%K	%B	%n'`
    elseif !isempty(taxids)
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -tax_info -outfmt '%T	%S	%L	%K	%B	%n'`
    end
    df = CSV.read(
        open(cmd),
        DataFrames.DataFrame;
        header = col_names,
        delim = '\t',
        comment = "#"
    )
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Convert a BLAST database to Arrow format with sequence and taxonomy information.
# Uses a simple serial approach with direct Arrow streaming for minimal memory usage.

# # Arguments
# - `blastdb::String`: Path to the BLAST database
# - `outfile::String=""`: Output file path. If empty, generates name based on input database
# - `force::Bool=false`: Whether to overwrite existing output file

# # Returns
# - `String`: Path to the generated output file (.arrow)

# # Output Format
# Arrow file containing columns (in this order):
# - sequence SHA256
# - sequence hash
# - sequence id
# - accession
# - gi
# - sequence title
# - BLAST name
# - taxid
# - taxonomic super kingdom
# - scientific name
# - scientific names for leaf-node taxids
# - common taxonomic name
# - common taxonomic names for leaf-node taxids
# - leaf-node taxids
# - membership integer
# - ordinal id
# - PIG
# - sequence length
# - sequence
# """
# function blastdb2table(;blastdb, outfile="", force=false)
#     # Set up environment and validate database
#     Mycelia.add_bioconda_env("blast")
#     blast_db_info = Mycelia.local_blast_database_info()
#     filtered = blast_db_info[blast_db_info[!, "BLAST database path"] .== blastdb, :]
#     @assert DataFrames.nrow(filtered) == 1
#     blast_db_info = filtered[1, :]
#     @show blast_db_info

#     # Determine database type and file extension
#     if blast_db_info["BLAST database molecule type"] == "Protein"
#         extension = ".faa"
#     elseif blast_db_info["BLAST database molecule type"] == "Nucleotide"
#         extension = ".fna"
#     else
#         @show blast_db_info["BLAST database molecule type"]
#         error("unexpected blast database molecule type")
#     end

#     # Set output file name if not provided
#     if outfile == ""
#         outfile = blastdb * extension * ".arrow"
#     end

#     # Skip processing if output exists (unless forced)
#     if isfile(outfile) && filesize(outfile) > 0 && !force
#         @show Mycelia.filesize_human_readable(outfile)
#         return outfile
#     end

#     # Define the mapping from outfmt symbols to column names.
#     symbol_header_map = OrderedCollections.OrderedDict(
#         "%s" => "sequence",            # field 1
#         "%a" => "accession",           # field 2
#         "%g" => "gi",                  # field 3
#         "%o" => "ordinal id",          # field 4
#         "%i" => "sequence id",         # field 5
#         "%t" => "sequence title",      # field 6
#         "%l" => "sequence length",     # field 7
#         "%h" => "sequence hash",       # field 8
#         "%T" => "taxid",               # field 9
#         "%X" => "leaf-node taxids",    # field 10
#         "%e" => "membership integer",  # field 11
#         "%L" => "common taxonomic name",   # field 12
#         "%C" => "common taxonomic names for leaf-node taxids",  # field 13
#         "%S" => "scientific name",     # field 14
#         "%N" => "scientific names for leaf-node taxids",  # field 15
#         "%B" => "BLAST name",          # field 16
#         "%K" => "taxonomic super kingdom", # field 17
#         "%P" => "PIG"                  # field 18
#     )
#     outfmt_string = join(collect(keys(symbol_header_map)), '\t')

#     # Define the desired output column order.
#     header_order = [
#         "sequence SHA256",
#         "sequence hash",
#         "sequence id",
#         "accession",
#         "gi",
#         "sequence title",
#         "BLAST name",
#         "taxid",
#         "taxonomic super kingdom",
#         "scientific name",
#         "scientific names for leaf-node taxids",
#         "common taxonomic name",
#         "common taxonomic names for leaf-node taxids",
#         "leaf-node taxids",
#         "membership integer",
#         "ordinal id",
#         "PIG",
#         "sequence length"
#         # "sequence"
#     ]

#     total_sequences = blast_db_info["number of sequences"]
#     progress = ProgressMeter.Progress(total_sequences, desc="Converting BLAST DB to Arrow: ", dt=1.0)

#     cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -entry all -outfmt $(outfmt_string)`
#     open(Arrow.Writer, outfile) do writer
#         for (i, line) in enumerate(eachline(cmd))
#             fields = split(strip(line), '\t')
#             # Pad fields with empty strings if necessary
#             if length(fields) < length(symbol_header_map)
#                 fields = vcat(fields, fill("", length(symbol_header_map) - length(fields)))
#             end
#             # Process sequence: clean and compute SHA256 from field 1 (the raw sequence)
#             seq = uppercase(String(filter(x -> isvalid(Char, x), fields[1])))
#             seq_sha256 = Mycelia.seq2sha256(seq)

#             # Build a mapping from column names to the extracted field values.
#             mapped = Dict{String, String}()
#             idx = 1
#             for (_, colname) in symbol_header_map
#                 mapped[colname] = fields[idx]
#                 idx += 1
#             end

#             # Construct a NamedTuple which is Arrow-compatible
#             row = DataFrames.DataFrame(
#                 sequence_SHA256 = seq_sha256,
#                 sequence_hash = mapped["sequence hash"],
#                 sequence_id = mapped["sequence id"],
#                 accession = mapped["accession"],
#                 gi = mapped["gi"],
#                 sequence_title = mapped["sequence title"],
#                 BLAST_name = mapped["BLAST name"],
#                 taxid = mapped["taxid"],
#                 taxonomic_super_kingdom = mapped["taxonomic super kingdom"],
#                 scientific_name = mapped["scientific name"],
#                 scientific_names_for_leaf_node_taxids = mapped["scientific names for leaf-node taxids"],
#                 common_taxonomic_name = mapped["common taxonomic name"],
#                 common_taxonomic_names_for_leaf_node_taxids = mapped["common taxonomic names for leaf-node taxids"],
#                 leaf_node_taxids = mapped["leaf-node taxids"],
#                 membership_integer = mapped["membership integer"],
#                 ordinal_id = mapped["ordinal id"],
#                 PIG = mapped["PIG"],
#                 sequence_length = mapped["sequence length"],
#                 # sequence = mapped["sequence"]
#             )
#             Arrow.write(writer, row)

#             # Update progress meter
#             ProgressMeter.next!(progress; showvalues = [
#                 (:completed, i),
#                 (:total, total_sequences),
#                 (:percent, round(i/total_sequences*100, digits=1))
#             ])
#         end
#     end

#     ProgressMeter.finish!(progress)
#     println("Done! Output saved to $(outfile)")
#     return outfile
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Convert a BLAST database to an in-memory table with sequence and taxonomy information.

# # Arguments
# - `blastdb::String`: Path to the BLAST database
# - `outfile::String=""`: Optional output file path. If provided, results will be saved to this file
# - `force::Bool=false`: Whether to overwrite existing output file
# - `ALL_FIELDS::Bool=true`: If true, include all fields regardless of other flag settings
# - Field selection flags (default to false unless ALL_FIELDS is true):
#   - `sequence_sha256::Bool`: Include SHA256 hash of the sequence
#   - `sequence_hash::Bool`: Include sequence hash
#   - `sequence_id::Bool`: Include sequence ID
#   - `accession::Bool`: Include accession number
#   - `gi::Bool`: Include GI number
#   - `sequence_title::Bool`: Include sequence title
#   - `blast_name::Bool`: Include BLAST name
#   - `taxid::Bool`: Include taxid
#   - `taxonomic_super_kingdom::Bool`: Include taxonomic super kingdom
#   - `scientific_name::Bool`: Include scientific name
#   - `scientific_names_leaf_nodes::Bool`: Include scientific names for leaf-node taxids
#   - `common_taxonomic_name::Bool`: Include common taxonomic name
#   - `common_names_leaf_nodes::Bool`: Include common taxonomic names for leaf-node taxids
#   - `leaf_node_taxids::Bool`: Include leaf-node taxids
#   - `membership_integer::Bool`: Include membership integer
#   - `ordinal_id::Bool`: Include ordinal ID
#   - `pig::Bool`: Include PIG
#   - `sequence_length::Bool`: Include sequence length
#   - `sequence::Bool`: Include the full sequence

# # Returns
# - `DataFrame`: DataFrame containing the requested columns from the BLAST database
# """
# function blastdb2table(;
#     blastdb, 
#     outfile="", 
#     force=false,
#     # Master field selection flag
#     ALL_FIELDS=true,
#     # Individual field selection flags
#     sequence_sha256=false,
#     sequence_hash=false,
#     sequence_id=false,
#     accession=false,
#     gi=false,
#     sequence_title=false,
#     blast_name=false,
#     taxid=false,
#     taxonomic_super_kingdom=false,
#     scientific_name=false,
#     scientific_names_leaf_nodes=false,
#     common_taxonomic_name=false,
#     common_names_leaf_nodes=false,
#     leaf_node_taxids=false,
#     membership_integer=false,
#     ordinal_id=false,
#     pig=false,
#     sequence_length=false,
#     sequence=false
# )
#     # If ALL_FIELDS is true, override all field flags to true
#     if ALL_FIELDS
#         sequence_sha256 = true
#         sequence_hash = true
#         sequence_id = true
#         accession = true
#         gi = true
#         sequence_title = true
#         blast_name = true
#         taxid = true
#         taxonomic_super_kingdom = true
#         scientific_name = true
#         scientific_names_leaf_nodes = true
#         common_taxonomic_name = true
#         common_names_leaf_nodes = true
#         leaf_node_taxids = true
#         membership_integer = true
#         ordinal_id = true
#         pig = true
#         sequence_length = true
#         sequence = true
#     end

#     Mycelia.add_bioconda_env("blast")
#     blast_db_info = Mycelia.local_blast_database_info()
#     filtered = blast_db_info[blast_db_info[!, "BLAST database path"] .== blastdb, :]
#     @assert DataFrames.nrow(filtered) == 1
#     blast_db_info = filtered[1, :]
#     @show blast_db_info

#     # Determine database type
#     if blast_db_info["BLAST database molecule type"] == "Protein"
#         extension = ".faa"
#     elseif blast_db_info["BLAST database molecule type"] == "Nucleotide"
#         extension = ".fna"
#     else
#         @show blast_db_info["BLAST database molecule type"]
#         error("unexpected blast database molecule type")
#     end

#     # Check if output file exists (if specified)
#     if outfile != "" && isfile(outfile) && filesize(outfile) > 0 && !force
#         @show Mycelia.filesize_human_readable(outfile)
#         # Load and return the existing file
#         return Arrow.Table(outfile) |> DataFrames.DataFrame
#     end

#     # Check if sequence needs to be extracted for SHA256 calculation
#     needs_sequence = sequence || sequence_sha256

#     # Define the mapping from outfmt symbols to column names and their inclusion status
#     field_config = [
#         # format, column_name, include_flag, dependency (if any)
#         ("%s", "sequence", needs_sequence, nothing),  # Need sequence for SHA256
#         ("%a", "accession", accession, nothing),
#         ("%g", "gi", gi, nothing),
#         ("%o", "ordinal_id", ordinal_id, nothing),
#         ("%i", "sequence_id", sequence_id, nothing),
#         ("%t", "sequence_title", sequence_title, nothing),
#         ("%l", "sequence_length", sequence_length, nothing),
#         ("%h", "sequence_hash", sequence_hash, nothing),
#         ("%T", "taxid", taxid, nothing),
#         ("%X", "leaf_node_taxids", leaf_node_taxids, nothing),
#         ("%e", "membership_integer", membership_integer, nothing),
#         ("%L", "common_taxonomic_name", common_taxonomic_name, nothing),
#         ("%C", "common_names_leaf_nodes", common_names_leaf_nodes, nothing),
#         ("%S", "scientific_name", scientific_name, nothing),
#         ("%N", "scientific_names_leaf_nodes", scientific_names_leaf_nodes, nothing),
#         ("%B", "blast_name", blast_name, nothing),
#         ("%K", "taxonomic_super_kingdom", taxonomic_super_kingdom, nothing),
#         ("%P", "pig", pig, nothing)
#     ]

#     # Generate format string for fields we need to fetch
#     formats_to_fetch = [fmt for (fmt, _, include, _) in field_config if include]
#     outfmt_string = join(formats_to_fetch, '\t')

#     # Create field name mapping for processing
#     format_to_colname = Dict(fmt => colname for (fmt, colname, _, _) in field_config)

#     # Create list of field names to include in final output
#     output_columns = []
#     if sequence_sha256
#         push!(output_columns, "sequence_sha256")
#     end

#     for (_, colname, include, _) in field_config
#         if include && colname != "sequence"  # Handle sequence separately
#             push!(output_columns, colname)
#         end
#     end

#     if sequence
#         push!(output_columns, "sequence")
#     end

#     # Create DataFrame to hold results
#     result_df = DataFrames.DataFrame()

#     total_sequences = blast_db_info["number of sequences"]
#     progress = ProgressMeter.Progress(total_sequences, desc="Processing BLAST DB: ", dt=1.0)

#     # Run blastdbcmd to get data
#     cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -entry all -outfmt $(outfmt_string)`

#     # Process results and build DataFrame
#     rows = []

#     for (i, line) in enumerate(eachline(cmd))
#         fields = split(strip(line), '\t')
#         # Pad fields with empty strings if necessary
#         if length(fields) < length(formats_to_fetch)
#             fields = vcat(fields, fill("", length(formats_to_fetch) - length(fields)))
#         end

#         # Map fields to column names
#         field_map = Dict{String, String}()
#         for (j, fmt) in enumerate(formats_to_fetch)
#             field_map[format_to_colname[fmt]] = fields[j]
#         end

#         # Process sequence if needed (for SHA256 or to include in output)
#         seq_sha256 = ""
#         if needs_sequence
#             seq = uppercase(String(filter(x -> isvalid(Char, x), field_map["sequence"])))
#             if sequence_sha256
#                 seq_sha256 = Mycelia.seq2sha256(seq)
#             end
#         end

#         # Create a Dict for this row
#         row = Dict{String, String}()

#         # Add sequence_sha256 if requested
#         if sequence_sha256
#             row["sequence_sha256"] = seq_sha256
#         end

#         # Add other fields if requested
#         for (_, colname, include, _) in field_config
#             if include && colname != "sequence"  # Handle sequence separately
#                 row[colname] = field_map[colname]
#             end
#         end

#         # Add sequence if requested
#         if sequence
#             row["sequence"] = field_map["sequence"]
#         end

#         push!(rows, row)

#         # Update progress meter
#         ProgressMeter.next!(progress; showvalues = [
#             (:completed, i),
#             (:total, total_sequences),
#             (:percent, round(i/total_sequences*100, digits=1))
#         ])
#     end

#     # Convert rows to DataFrame
#     result_df = DataFrames.DataFrame(rows)

#     # Reorder columns if necessary to match expected order
#     final_columns = filter(col -> col in names(result_df), output_columns)
#     result_df = result_df[:, final_columns]

#     ProgressMeter.finish!(progress)
#     return result_df
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BLAST database to an in-memory table with sequence and taxonomy information.

# Arguments
- `blastdb::String`: Path to the BLAST database
- `outfile::String=""`: Optional output file path. If provided, results will be saved to this file
- `force::Bool=false`: Whether to overwrite existing output file
- `ALL_FIELDS::Bool=true`: If true, include all fields regardless of other flag settings
- Field selection flags (default to false unless ALL_FIELDS is true):
  - `sequence_sha256::Bool`: Include SHA256 hash of the sequence
  - `sequence_hash::Bool`: Include sequence hash
  - `sequence_id::Bool`: Include sequence ID
  - `accession::Bool`: Include accession number
  - `gi::Bool`: Include GI number
  - `sequence_title::Bool`: Include sequence title
  - `blast_name::Bool`: Include BLAST name
  - `taxid::Bool`: Include taxid
  - `taxonomic_super_kingdom::Bool`: Include taxonomic super kingdom
  - `scientific_name::Bool`: Include scientific name
  - `scientific_names_leaf_nodes::Bool`: Include scientific names for leaf-node taxids
  - `common_taxonomic_name::Bool`: Include common taxonomic name
  - `common_names_leaf_nodes::Bool`: Include common taxonomic names for leaf-node taxids
  - `leaf_node_taxids::Bool`: Include leaf-node taxids
  - `membership_integer::Bool`: Include membership integer
  - `ordinal_id::Bool`: Include ordinal ID
  - `pig::Bool`: Include PIG
  - `sequence_length::Bool`: Include sequence length
  - `sequence::Bool`: Include the full sequence

# Returns
- `DataFrame`: DataFrame containing the requested columns from the BLAST database
"""
function blastdb2table(;
        blastdb::String,
        blastdbs_dir::Union{Nothing, String} = nothing,
        # Master field selection flag
        ALL_FIELDS::Bool = true,
        # Individual field selection flags
        sequence_sha256::Bool = false,
        sequence_hash::Bool = false,
        sequence_id::Bool = false,
        accession::Bool = false,
        gi::Bool = false,
        sequence_title::Bool = false,
        blast_name::Bool = false,
        taxid::Bool = false,
        taxonomic_super_kingdom::Bool = false,
        scientific_name::Bool = false,
        scientific_names_leaf_nodes::Bool = false,
        common_taxonomic_name::Bool = false,
        common_names_leaf_nodes::Bool = false,
        leaf_node_taxids::Bool = false,
        membership_integer::Bool = false,
        ordinal_id::Bool = false,
        pig::Bool = false,
        sequence_length::Bool = false,
        sequence::Bool = false
)

    # If ALL_FIELDS is true, override all field flags to true
    if ALL_FIELDS
        sequence_sha256 = true
        sequence_hash = true
        sequence_id = true
        accession = true
        gi = true
        sequence_title = true
        blast_name = true
        taxid = true
        taxonomic_super_kingdom = true
        scientific_name = true
        scientific_names_leaf_nodes = true
        common_taxonomic_name = true
        common_names_leaf_nodes = true
        leaf_node_taxids = true
        membership_integer = true
        ordinal_id = true
        pig = true
        sequence_length = true
        sequence = true
    end

    Mycelia.add_bioconda_env("blast")
    # Check if sequence needs to be extracted for SHA256 calculation
    needs_sequence = sequence || sequence_sha256

    # Define the mapping from outfmt symbols to column names and their inclusion status
    field_config = [
        ("%s", "sequence", needs_sequence),  # Need sequence for SHA256
        ("%a", "accession", accession),
        ("%g", "gi", gi),
        ("%o", "ordinal_id", ordinal_id),
        ("%i", "sequence_id", sequence_id),
        ("%t", "sequence_title", sequence_title),
        ("%l", "sequence_length", sequence_length),
        ("%h", "sequence_hash", sequence_hash),
        ("%T", "taxid", taxid),
        ("%X", "leaf_node_taxids", leaf_node_taxids),
        ("%e", "membership_integer", membership_integer),
        ("%L", "common_taxonomic_name", common_taxonomic_name),
        ("%C", "common_names_leaf_nodes", common_names_leaf_nodes),
        ("%S", "scientific_name", scientific_name),
        ("%N", "scientific_names_leaf_nodes", scientific_names_leaf_nodes),
        ("%B", "blast_name", blast_name),
        ("%K", "taxonomic_super_kingdom", taxonomic_super_kingdom),
        ("%P", "pig", pig)
    ]

    # Generate format string for fields we need to fetch
    active_fields = [(fmt, colname) for (fmt, colname, include) in field_config if include]
    formats_to_fetch = [fmt for (fmt, _) in active_fields]
    columns_to_fetch = [colname for (_, colname) in active_fields]
    outfmt_string = join(formats_to_fetch, '\t')

    # Define output columns
    output_columns = String[]
    if sequence_sha256
        push!(output_columns, "sequence_sha256")
    end

    # Add the other active columns except "sequence" which is handled separately
    for (_, colname, include) in field_config
        if include && colname != "sequence"
            push!(output_columns, colname)
        end
    end

    if sequence
        push!(output_columns, "sequence")
    end

    # If no columns were requested, return an empty table early.
    if isempty(output_columns)
        @info "No BLAST fields requested; returning empty table"
        return DataFrames.DataFrame()
    end

    # Decide which directory to scan for BLAST databases. Prefer the provided
    # directory, otherwise use the directory component of `blastdb` when given,
    # and finally fall back to the default BLASTDB path.
    db_dir = if blastdbs_dir !== nothing
        blastdbs_dir
    elseif occursin('/', blastdb) || occursin('\\', blastdb)
        dirname(blastdb)
    else
        Mycelia.DEFAULT_BLASTDB_PATH
    end
    db_dir = abspath(db_dir)
    target_db = isabspath(blastdb) ? blastdb : joinpath(db_dir, basename(blastdb))

    blast_db_info = Mycelia.local_blast_database_info(; blastdbs_dir = db_dir)
    filtered = blast_db_info[blast_db_info[!, "BLAST database path"] .== target_db, :]
    if DataFrames.nrow(filtered) != 1
        @info "BLAST database not found in local listing; returning empty table" blastdb target_db DataFrames.nrow(filtered)
        return DataFrames.DataFrame(Dict(Symbol(col) => String[] for col in output_columns))
    end
    blast_db_info = filtered[1, :]

    # Determine database type
    if blast_db_info["BLAST database molecule type"] == "Protein"
        extension = ".faa"
    elseif blast_db_info["BLAST database molecule type"] == "Nucleotide"
        extension = ".fna"
    else
        @show blast_db_info["BLAST database molecule type"]
        error("unexpected blast database molecule type")
    end

    # Use dynamic storage for all columns
    column_data = OrderedCollections.OrderedDict{String, Vector{String}}()
    for col in output_columns
        column_data[col] = String[]
    end

    # For progress bar, try to get total sequences, else fallback to nothing
    total_sequences = get(blast_db_info, "number of sequences", nothing)
    if total_sequences !== nothing
        progress = ProgressMeter.Progress(total_sequences, desc = "Processing BLAST DB: ", dt = 1.0)
    else
        progress = ProgressMeter.Progress(desc = "Processing BLAST DB: ", dt = 1.0)
    end

    # Run blastdbcmd to get data
    cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -entry all -outfmt $(outfmt_string)`

    # Sequence index in fields (if needed)
    sequence_idx = findfirst(==("%s"), formats_to_fetch)

    # Process each line from blastdbcmd
    nseqs = 0
    for line in eachline(cmd)
        nseqs += 1
        fields = split(strip(line), '\t')
        # Pad fields with empty strings if necessary
        if length(fields) < length(formats_to_fetch)
            append!(fields, fill("", length(formats_to_fetch) - length(fields)))
        end

        # Process sequence if needed (for SHA256 or to include in output)
        seq = ""
        if needs_sequence && sequence_idx !== nothing
            seq_raw = fields[sequence_idx]
            seq = uppercase(String(filter(x -> isvalid(Char, x), seq_raw)))
            if sequence_sha256
                seq_sha256 = Mycelia.seq2sha256(seq)
                push!(column_data["sequence_sha256"], seq_sha256)
            end
        elseif sequence_sha256
            # If sequence is not available, push empty string for sha256
            push!(column_data["sequence_sha256"], "")
        end

        # Add other fields directly to column vectors
        for (j, colname) in enumerate(columns_to_fetch)
            # sequence is handled below
            if colname != "sequence"
                push!(column_data[colname], fields[j])
            end
        end

        # Add sequence field at the end if requested
        if sequence
            push!(column_data["sequence"], seq)
        end

        # Update progress meter
        if total_sequences !== nothing
            ProgressMeter.next!(progress;
                showvalues = [
                    (:completed, nseqs),
                    (:total, total_sequences),
                    (:percent, round(nseqs/total_sequences*100, digits = 1))
                ])
        else
            ProgressMeter.next!(progress)
        end
    end

    # Consistency check
    nrows = length(first(values(column_data)))
    for (col, vec) in column_data
        @assert length(vec) == nrows "Column $col has length $(length(vec)), expected $nrows"
    end

    # Construct DataFrame directly from column vectors
    result_df = DataFrames.DataFrame(column_data)

    ProgressMeter.finish!(progress)
    return result_df
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Convert a BLAST database to Arrow format with taxonomy information.

# # Arguments
# - `blastdb::String`: Path to the BLAST database
# - `outfile::String=""`: Output file path. If empty, generates name based on input database
# - `force::Bool=false`: Whether to overwrite existing output file

# # Returns
# - `String`: Path to the generated output file (.arrow)

# # Output Format
# Arrow file containing columns (in this order):
# - sequence hash
# - sequence id
# - accession
# - gi
# - sequence title
# - BLAST name
# - taxid
# - taxonomic super kingdom
# - scientific name
# - scientific names for leaf-node taxids
# - common taxonomic name
# - common taxonomic names for leaf-node taxids
# - leaf-node taxids
# """
# function blastdb2tax_table(; blastdb)
#     # Set up environment and validate database
#     Mycelia.add_bioconda_env("blast")
#     blast_db_info = Mycelia.local_blast_database_info()
#     filtered = blast_db_info[blast_db_info[!, "BLAST database path"] .== blastdb, :]
#     @assert DataFrames.nrow(filtered) == 1
#     blast_db_info = filtered[1, :]
#     @show blast_db_info

#     # Define the mapping from outfmt symbols to column names.
#     symbol_header_map = OrderedCollections.OrderedDict(
#         "%a" => "accession",           # field 2
#         "%g" => "gi",                  # field 3
#         "%i" => "sequence id",         # field 5
#         "%t" => "sequence title",      # field 6
#         "%h" => "sequence hash",       # field 8
#         "%T" => "taxid",               # field 9
#         "%X" => "leaf-node taxids",    # field 10
#         "%L" => "common taxonomic name",   # field 12
#         "%C" => "common taxonomic names for leaf-node taxids",  # field 13
#         "%S" => "scientific name",     # field 14
#         "%N" => "scientific names for leaf-node taxids",  # field 15
#         "%B" => "BLAST name",          # field 16
#         "%K" => "taxonomic super kingdom", # field 17
#     )
#     outfmt_string = join(collect(keys(symbol_header_map)), '\t')
#     cmd = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(blastdb) -entry all -outfmt $(outfmt_string)`
#     return CSV.read(open(cmd), DataFrames.DataFrame, delim='\t', header=collect(values(symbol_header_map)))
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads Sequence Read Archive (SRA) data using the prefetch tool from sra-tools.

# Arguments
- `SRR`: SRA accession number (e.g., "SRR12345678")
- `outdir`: Directory where the downloaded data will be saved. Defaults to current directory.

# Notes
- Requires sra-tools which will be installed in a Conda environment
- Downloads are saved in .sra format
- Internet connection required
"""
function prefetch(; SRR, outdir = pwd())
    Mycelia.add_bioconda_env("sra-tools")
    final_dir = joinpath(outdir, SRR)
    sra_archive = joinpath(final_dir, "$(SRR).sra")
    if !isfile(sra_archive)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n sra-tools prefetch $(SRR) -O $(outdir)`)
    else
        @info "SRA archive already present: $(sra_archive)"
    end
    return (directory = final_dir, archive = sra_archive)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Download and compress sequencing reads from the SRA database using fasterq-dump.

# Arguments
- `outdir::String=""`: Output directory for the FASTQ files. Defaults to current directory.
- `srr_identifier::String=""`: SRA run accession number (e.g., "SRR12345678")

# Returns
Named tuple containing paths to the generated files:
- `forward_reads`: Path to forward reads file (*_1.fastq.gz) or `missing`
- `reverse_reads`: Path to reverse reads file (*_2.fastq.gz) or `missing`
- `unpaired_reads`: Path to unpaired reads file (*.fastq.gz) or `missing`

# Outputs
Creates compressed FASTQ files in the output directory:
- `{srr_identifier}_1.fastq.gz`: Forward reads (for paired-end data)
- `{srr_identifier}_2.fastq.gz`: Reverse reads (for paired-end data)
- `{srr_identifier}.fastq.gz`: Unpaired reads (for single-end data)

# Dependencies
Requires:
- `fasterq-dump` from the SRA Toolkit (installed via Conda)
- `gzip` for compression

# Notes
- Skips download if output files already exist
- Uses up to 4 threads or system maximum, whichever is lower
- Allocates 1GB memory for processing
- Skips technical reads
- Handles both paired-end and single-end data automatically
"""
function fasterq_dump(; outdir = pwd(), srr_identifier = "")
    Mycelia.add_bioconda_env("sra-tools")
    prefetch_results = Mycelia.prefetch(SRR = srr_identifier, outdir = outdir)

    final_outdir = prefetch_results.directory

    forward_reads = joinpath(final_outdir, "$(srr_identifier)_1.fastq")
    reverse_reads = joinpath(final_outdir, "$(srr_identifier)_2.fastq")
    unpaired_reads = joinpath(final_outdir, "$(srr_identifier).fastq")

    forward_reads_gz = forward_reads * ".gz"
    reverse_reads_gz = reverse_reads * ".gz"
    unpaired_reads_gz = unpaired_reads * ".gz"

    forward_and_reverse_present = isfile(forward_reads_gz) && isfile(reverse_reads_gz)
    unpaired_present = isfile(unpaired_reads_gz)

    if !(forward_and_reverse_present || unpaired_present)
        # --progress doesn't work well for jupyter output
        fasterq_dump_cmd = `
            $(Mycelia.CONDA_RUNNER) run --live-stream -n sra-tools fasterq-dump
                --outdir $(final_outdir)
                --mem 1G
                --split-3
                --threads $(min(get_default_threads(), 4))
                --skip-technical
                $(final_outdir)`
        @time run(fasterq_dump_cmd)
        isfile(forward_reads) && run(`gzip $(forward_reads)`)
        isfile(reverse_reads) && run(`gzip $(reverse_reads)`)
        isfile(unpaired_reads) && run(`gzip $(unpaired_reads)`)
    else
        @info "$(forward_reads_gz) & $(reverse_reads_gz) already present"
    end
    return (
        forward_reads = isfile(forward_reads_gz) ? forward_reads_gz : missing,
        reverse_reads = isfile(reverse_reads_gz) ? reverse_reads_gz : missing,
        unpaired_reads = isfile(unpaired_reads_gz) ? unpaired_reads_gz : missing
    )
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Downloads and quality filters paired-end reads from the Sequence Read Archive (SRA).

# # Arguments
# - `outdir::String`: Output directory path for downloaded and processed files
# - `srr_identifier::String`: SRA run accession number (e.g., "SRR12345678")

# # Details
# 1. Downloads paired-end FASTQ files using fasterq-dump
# 2. Performs quality trimming using trim_galore
# 3. Removes intermediate compressed FASTQ files after processing

# # Returns
# Nothing, but creates the following files in `outdir`:
# - `trim_galore/[srr_identifier]_1_val_1.fq.gz`: Trimmed forward reads
# - `trim_galore/[srr_identifier]_2_val_2.fq.gz`: Trimmed reverse reads
# """
"""
Downloads sequencing reads from NCBI's Sequence Read Archive (SRA).

Downloads reads using fasterq-dump. The function automatically detects whether
the data is single-end or paired-end and returns appropriate file paths.
Users should apply quality control based on their knowledge of the data type.

# Arguments
- `srr_identifier`: SRA run identifier (e.g., "SRR1234567")
- `outdir`: Output directory for downloaded files (default: current directory)

# Returns
Named tuple with:
- `srr_id`: The SRA identifier
- `outdir`: Output directory path
- `files`: Vector of downloaded file paths (1 file for single-end, 2 for paired-end)
- `is_paired`: Boolean indicating if data is paired-end

# Example
```julia
# Download SRA data
result = Mycelia.download_sra_data("SRR1234567", outdir="./data")

# Apply appropriate QC based on data type
if result.is_paired
    # Paired-end data - use paired-end QC
    Mycelia.trim_galore_paired(forward_reads=result.files[1], reverse_reads=result.files[2])
else
    # Single-end data - use single-end QC
    Mycelia.qc_filter_short_reads_fastp(input=result.files[1])
end
```
"""
function download_sra_data(srr_identifier::String; outdir::String = pwd())
    if isempty(srr_identifier)
        error("SRA identifier cannot be empty")
    end

    @info "Downloading SRA data: $(srr_identifier)"
    fasterq_dump(outdir = outdir, srr_identifier = srr_identifier)

    # Check what files were created to determine data type
    forward_reads = joinpath(outdir, "$(srr_identifier)_1.fastq.gz")
    reverse_reads = joinpath(outdir, "$(srr_identifier)_2.fastq.gz")
    single_reads = joinpath(outdir, "$(srr_identifier).fastq.gz")

    if isfile(forward_reads) && isfile(reverse_reads)
        # Paired-end data
        return (
            srr_id = srr_identifier,
            outdir = outdir,
            files = [forward_reads, reverse_reads],
            is_paired = true
        )
    elseif isfile(single_reads)
        # Single-end data
        return (
            srr_id = srr_identifier,
            outdir = outdir,
            files = [single_reads],
            is_paired = false
        )
    else
        error("No FASTQ files found after download. Check SRA identifier: $(srr_identifier)")
    end
end

"""
Prefetches multiple SRA runs in parallel.

Downloads SRA run files (.sra) to local storage without converting to FASTQ.
Useful for batch downloading before processing with fasterq-dump.

# Arguments
- `srr_identifiers`: Vector of SRA run identifiers
- `outdir`: Output directory for prefetched files (default: current directory)
- `max_parallel`: Maximum number of parallel downloads (default: 4)

# Returns
Vector of named tuples with prefetch results for each SRA run

# Example
```julia
runs = ["SRR1234567", "SRR1234568", "SRR1234569"]
results = Mycelia.prefetch_sra_runs(runs, outdir="./sra_data")
```
"""
function prefetch_sra_runs(srr_identifiers::Vector{String}; outdir::String = pwd(), max_parallel::Int = 4)
    if isempty(srr_identifiers)
        error("No SRA identifiers provided")
    end

    results = []

    @info "Prefetching $(length(srr_identifiers)) SRA runs with max $(max_parallel) parallel downloads"

    # Process in chunks to respect max_parallel limit
    for chunk in Iterators.partition(srr_identifiers, max_parallel)
        chunk_tasks = []

        for srr_id in chunk
            task = Threads.@spawn begin
                try
                    @info "Prefetching $(srr_id)"
                    result = prefetch(SRR = srr_id, outdir = outdir)
                    @info "Completed prefetch for $(srr_id)"
                    (srr_id = srr_id, success = true, result = result, error = nothing)
                catch e
                    @error "Failed to prefetch $(srr_id): $(e)"
                    (srr_id = srr_id, success = false, result = nothing, error = string(e))
                end
            end
            push!(chunk_tasks, task)
        end

        # Wait for chunk to complete
        chunk_results = [fetch(task) for task in chunk_tasks]
        append!(results, chunk_results)
    end

    successful = sum([r.success for r in results])
    @info "Prefetch completed: $(successful)/$(length(srr_identifiers)) successful"

    return results
end

"""
Parallel FASTQ dump for multiple SRA files.

Converts multiple SRA files to FASTQ format in parallel. More efficient than
sequential processing for large batches.

# Arguments
- `srr_identifiers`: Vector of SRA run identifiers
- `outdir`: Output directory for FASTQ files (default: current directory)
- `max_parallel`: Maximum number of parallel conversions (default: 2)

# Returns
Vector of named tuples with conversion results

# Example
```julia
runs = ["SRR1234567", "SRR1234568"]
results = Mycelia.fasterq_dump_parallel(runs, outdir="./fastq_data")
```
"""
function fasterq_dump_parallel(srr_identifiers::Vector{String}; outdir::String = pwd(), max_parallel::Int = 2)
    if isempty(srr_identifiers)
        error("No SRA identifiers provided")
    end

    results = []

    @info "Converting $(length(srr_identifiers)) SRA runs to FASTQ with max $(max_parallel) parallel processes"

    # Process in chunks - fewer parallel for fasterq-dump as it's more resource intensive
    for chunk in Iterators.partition(srr_identifiers, max_parallel)
        chunk_tasks = []

        for srr_id in chunk
            task = Threads.@spawn begin
                try
                    @info "Converting $(srr_id) to FASTQ"
                    fasterq_dump(outdir = outdir, srr_identifier = srr_id)
                    @info "Completed FASTQ conversion for $(srr_id)"
                    (srr_id = srr_id, success = true, error = nothing)
                catch e
                    @error "Failed to convert $(srr_id): $(e)"
                    (srr_id = srr_id, success = false, error = string(e))
                end
            end
            push!(chunk_tasks, task)
        end

        # Wait for chunk to complete
        chunk_results = [fetch(task) for task in chunk_tasks]
        append!(results, chunk_results)
    end

    successful = sum([r.success for r in results])
    @info "FASTQ conversion completed: $(successful)/$(length(srr_identifiers)) successful"

    return results
end

function _un_request(url::AbstractString, method::AbstractString;
        headers::Dict{String, String} = Dict{String, String}(),
        max_redirects::Int = 5
)
    response = HTTP.request(method, url; headers = headers, status_exception = false)
    if response.status in (301, 302, 303, 307, 308)
        max_redirects > 0 || error("Too many redirects for $(url)")
        location = HTTP.header(response, "Location")
        location === nothing && error("Redirect without Location for $(url)")
        return _un_request(location, method; headers = headers, max_redirects = max_redirects -
                                                                                1)
    end
    return response
end

function _un_part_exists(url::AbstractString)
    base_headers = Dict("User-Agent" => UN_CORPUS_USER_AGENT)
    response = _un_request(url, "HEAD"; headers = base_headers)
    if response.status == 404
        return false
    elseif 200 <= response.status < 300
        return true
    elseif response.status == 405 || response.status == 403
        response = _un_request(
            url,
            "GET";
            headers = merge(base_headers, Dict("Range" => "bytes=0-0"))
        )
        if response.status == 404
            return false
        elseif response.status == 200 || response.status == 206
            return true
        elseif response.status == 403
            return true
        end
    end

    error("Unexpected HTTP status $(response.status) for $(url)")
end

function _un_download(url::AbstractString, outfile::AbstractString)
    with_retry() do
        Downloads.download(url, outfile; headers = Dict("User-Agent" => UN_CORPUS_USER_AGENT))
    end
    return outfile
end

function _merge_un_archive_parts(parts::Vector{String}, merged_archive::AbstractString)
    open(merged_archive, "w") do out_io
        buffer = Vector{UInt8}(undef, 1024 * 1024)
        for part in parts
            open(part, "r") do in_io
                while !eof(in_io)
                    n = readbytes!(in_io, buffer)
                    if n == 0
                        break
                    end
                    write(out_io, view(buffer, 1:n))
                end
            end
        end
    end
    return merged_archive
end

function _download_and_extract_split_archive(
        base_name::AbstractString,
        url_base::AbstractString,
        outdir::AbstractString;
        force::Bool = false,
        verbose_extract::Bool = false
)
    mkpath(outdir)
    merged_archive = joinpath(outdir, base_name)

    if isfile(merged_archive) && !force
        @info "Using existing merged archive: $(merged_archive)"
        Mycelia.tar_extract(tarchive = merged_archive, directory = outdir, verbose = verbose_extract)
        return outdir
    end

    parts = String[]
    index = 0

    while true
        suffix = lpad(string(index), 2, '0')
        filename = "$(base_name).$(suffix)"
        url = "$(url_base)/$(filename)"

        if !_un_part_exists(url)
            break
        end

        outfile = joinpath(outdir, filename)
        if force || !isfile(outfile)
            @info "Downloading $(url)"
            _un_download(url, outfile)
        else
            @info "Using existing part: $(outfile)"
        end
        push!(parts, outfile)
        index += 1
    end

    if isempty(parts)
        direct_url = "$(url_base)/$(base_name)"
        if _un_part_exists(direct_url)
            outfile = joinpath(outdir, base_name)
            if force || !isfile(outfile)
                @info "Downloading $(direct_url)"
                _un_download(direct_url, outfile)
            else
                @info "Using existing archive: $(outfile)"
            end

            extracted = false
            try
                Mycelia.tar_extract(tarchive = outfile, directory = outdir, verbose = verbose_extract)
                extracted = true
            finally
                if extracted
                    isfile(outfile) && rm(outfile; force = true)
                else
                    @warn "Extraction failed; keeping downloaded archive for inspection."
                end
            end
            return outdir
        end

        error("No archive parts found for $(base_name) at $(url_base)")
    end

    extracted = false
    try
        @info "Merging $(length(parts)) archive parts into $(merged_archive)"
        _merge_un_archive_parts(parts, merged_archive)
        Mycelia.tar_extract(tarchive = merged_archive, directory = outdir, verbose = verbose_extract)
        extracted = true
    finally
        if extracted
            for part in parts
                isfile(part) && rm(part; force = true)
            end
            isfile(merged_archive) && rm(merged_archive; force = true)
        else
            @warn "Extraction failed; keeping downloaded archive parts for inspection."
        end
    end

    return outdir
end

function _resolve_un_archives(subsets::Vector{String}, formats::Vector{String})
    if isempty(subsets)
        error("subsets cannot be empty")
    end
    if isempty(formats)
        error("formats cannot be empty")
    end

    normalized_formats = unique([lowercase(strip(f)) for f in formats if !isempty(strip(f))])
    allowed_formats = Set(["txt", "xml"])
    if isempty(normalized_formats)
        error("formats cannot be empty")
    end
    for fmt in normalized_formats
        if !(fmt in allowed_formats)
            error("Unsupported format: $(fmt). Use \"txt\" or \"xml\".")
        end
    end

    allowed_langs = Set(UN_LANGUAGES)
    allowed_pairs = Set(UN_PAIRS)
    normalized_subsets = unique([lowercase(strip(s)) for s in subsets if !isempty(strip(s))])
    if isempty(normalized_subsets)
        error("subsets cannot be empty")
    end

    subset_tokens = String[]
    explicit_tokens = String[]
    pairs = String[]
    all_requested = false

    for token in normalized_subsets
        if token == "6-way"
            token = "6way"
        end

        if token in ("all", "bitext", "tei", "6way")
            if token == "all"
                all_requested = true
            else
                push!(subset_tokens, token)
                push!(explicit_tokens, token)
            end
            continue
        end

        if occursin(r"^[a-z]{2}-[a-z]{2}$", token)
            langs = split(token, "-")
            if !all(in(allowed_langs), langs)
                error("Unsupported language pair: $(token). Allowed languages: $(join(UN_LANGUAGES, ", ")).")
            end
            pair = join(sort(langs), "-")
            if !(pair in allowed_pairs)
                error("Unsupported language pair: $(token).")
            end
            push!(pairs, pair)
            continue
        end

        error("Unsupported subset: $(token). Use \"all\", \"bitext\", \"tei\", \"6way\", or a pair like \"en-fr\".")
    end

    pairs = unique(pairs)
    subset_tokens = unique(subset_tokens)
    explicit_tokens = unique(explicit_tokens)

    if all_requested
        subset_tokens = unique(vcat(subset_tokens, ["bitext", "tei", "6way"]))
    end

    want_txt = all_requested ? true : ("txt" in normalized_formats)
    want_xml = all_requested ? true : ("xml" in normalized_formats)

    if !all_requested
        if (("bitext" in explicit_tokens) || ("6way" in explicit_tokens) ||
            !isempty(pairs)) && !want_txt
            error("Bitext/6-way downloads require formats to include \"txt\".")
        end
        if ("tei" in explicit_tokens) && !want_xml
            error("TEI downloads require formats to include \"xml\".")
        end
    end

    archives = String[]
    if ("bitext" in subset_tokens) && want_txt
        append!(archives, ["UNv1.0.$(pair).tar.gz" for pair in UN_PAIRS])
    end
    if !isempty(pairs)
        append!(archives, ["UNv1.0.$(pair).tar.gz" for pair in pairs])
    end
    if ("6way" in subset_tokens) && want_txt
        push!(archives, "UNv1.0.6way.tar.gz")
    end
    if ("tei" in subset_tokens) && want_xml
        append!(archives, ["UNv1.0-TEI.$(lang).tar.gz" for lang in UN_LANGUAGES])
    end

    archives = unique(archives)
    if isempty(archives)
        error("No archives selected. Check subsets and formats.")
    end

    return archives
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Download and extract the United Nations Parallel Corpus v1.0.

# Arguments
- `outdir::String`: Output directory (default: current directory). Files are extracted into
  `joinpath(outdir, "un_corpus")`.
- `subsets::Vector{String}`: List of subsets to download. Options:
  - `"all"`: All subsets and formats (ignores `formats`)
  - `"bitext"`: All pairwise plain-text files
  - `"tei"`: All XML/TEI files
  - `"6way"`: The fully aligned 6-way subcorpus
  - Specific pairs: e.g., `"en-fr"`, `"ru-zh"` (order-insensitive)
- `formats::Vector{String}`: Formats to download. `"txt"` covers bitext and 6-way; `"xml"` covers TEI.
- `base_url::String`: Base URL for downloads.
- `force::Bool`: Re-download archive parts even if they already exist.

# Details
Handles the split archive format (`.tar.gz.00`, `.tar.gz.01`, ...) by downloading all parts,
concatenating them, and extracting the merged archive.
If `"all"` is included in `subsets`, all formats are downloaded regardless of `formats`.

# Returns
- `String`: Path to the corpus directory.

# Example
```julia
# Download English-French bitexts
Mycelia.download_un_parallel_corpus(outdir="./data", subsets=["en-fr"])

# Download TEI/XML for all languages
Mycelia.download_un_parallel_corpus(outdir="./data", subsets=["tei"], formats=["xml"])

# Download all subsets and formats
Mycelia.download_un_parallel_corpus(outdir="./data", subsets=["all"], formats=["txt", "xml"])
```
"""
function download_un_parallel_corpus(;
        outdir::String = pwd(),
        subsets::Vector{String} = ["all"],
        formats::Vector{String} = ["txt"],
        base_url::String = UN_CORPUS_BASE_URL,
        force::Bool = false,
        verbose_extract::Bool = false
)
    archives = _resolve_un_archives(subsets, formats)

    corpus_dir = joinpath(outdir, "un_corpus")
    mkpath(corpus_dir)
    url_base = rstrip(base_url, '/')

    for archive in archives
        @info "Downloading UN Parallel Corpus archive: $(archive)"
        _download_and_extract_split_archive(
            archive,
            url_base,
            corpus_dir;
            force = force,
            verbose_extract = verbose_extract
        )
    end

    return corpus_dir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Download and extract the United Nations Parallel Corpus v1.0 test sets.

# Arguments
- `outdir::String`: Output directory (default: current directory). Files are extracted into
  `joinpath(outdir, "un_corpus_testsets")`.
- `base_url::String`: Base URL for downloads.
- `force::Bool`: Re-download archive even if it already exists.

# Returns
- `String`: Path to the test set directory.
"""
function download_un_parallel_corpus_testset(;
        outdir::String = pwd(),
        base_url::String = UN_CORPUS_BASE_URL,
        force::Bool = false,
        verbose_extract::Bool = false
)
    testset_dir = joinpath(outdir, "un_corpus_testsets")
    mkpath(testset_dir)
    url_base = rstrip(base_url, '/')
    archive = "UNv1.0.testsets.tar.gz"
    _download_and_extract_split_archive(
        archive,
        url_base,
        testset_dir;
        force = force,
        verbose_extract = verbose_extract
    )
    return testset_dir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Collect UN Parallel Corpus plain-text files by language.

# Arguments
- `corpus_root::String`: Root directory containing extracted UN corpus files.
- `subsets::Vector{String}`: Subset filters (e.g. `["all"]`, `["en-fr"]`, `["6way"]`).
- `strict_naming::Bool`: If `true` (default), require "UNv1.0" and pair naming ("-" or "6way")
  in filenames. Set to `false` for testset/devset files that use simpler naming conventions.

# Returns
- `Dict{String, Vector{String}}`: Map of language code to file paths.
"""
function collect_un_language_files(corpus_root::String; subsets::Vector{String}, strict_naming::Bool = true)
    language_files = Dict{String, Vector{String}}()
    use_all = "all" in subsets
    for (root, _dirs, files) in walkdir(corpus_root)
        for filename in files
            filepath = joinpath(root, filename)
            if endswith(filename, ".xml") || endswith(filename, ".gz") ||
               endswith(filename, ".tar")
                continue
            end
            if strict_naming
                if !occursin("UNv1.0", filename)
                    continue
                end
                if !occursin("-", filename) && !occursin("6way", filename)
                    continue
                end
            end
            parts = split(filename, '.')
            if length(parts) < 2
                continue
            end
            lang = parts[end]
            if length(lang) > 3
                continue
            end
            if !all(isletter, lang)
                continue
            end
            if !use_all
                matches_subset = any(subset -> occursin(subset, filename), subsets)
                if !matches_subset
                    continue
                end
            end
            push!(get!(Vector{String}, language_files, lang), filepath)
        end
    end
    return language_files
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads a genomic sequence from NCBI's nucleotide database by its accession number.

# Arguments
- `accession::String`: NCBI nucleotide accession number (e.g. "NC_045512")
- `outdir::String`: Output directory path. Defaults to current directory
- `compressed::Bool`: If true, compresses output file with gzip. Defaults to true

# Returns
- `String`: Path to the downloaded file (.fna or .fna.gz)
"""
function download_genome_by_accession(; accession, outdir = pwd(), compressed = true)
    temp_fasta = joinpath(outdir, accession * ".fna")
    if compressed
        outfile = temp_fasta * ".gz"
    else
        outfile = temp_fasta
    end
    if !isfile(outfile)
        try
            Mycelia.write_fasta(outfile = outfile,
                records = Mycelia.get_sequence(db = "nuccore", accession = accession))
        catch e
            println("An error occurred: ", e)
        end
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads a genome file from NCBI FTP server to the specified directory.

# Arguments
- `ftp::String`: NCBI FTP path for the genome (e.g. "ftp://ftp.ncbi.nlm.nih.gov/.../")
- `outdir::String`: Output directory path. Defaults to current working directory.

# Returns
- `String`: Path to the downloaded file

# Notes
- If the target file already exists, returns the existing file path without re-downloading
- Downloads the genomic.fna.gz version of the genome
"""
function download_genome_by_ftp(; ftp, outdir = pwd())
    url = Mycelia.ncbi_ftp_path_to_url(ftp_path = ftp, extension = "genomic.fna.gz")
    outfile = joinpath(outdir, basename(url))
    if !isfile(outfile)
        return Downloads.download(url, outfile)
    else
        return outfile
    end
end

function download_genomes_by_ftp(;
        ftp_paths,
        outdir = pwd()
)
    results = DataFrames.DataFrame(ftp_path = String[], fna_path = String[])
    n = length(ftp_paths)
    p = ProgressMeter.Progress(n; desc = "Downloading genomes: ", dt = 0.5)
    prog_lock = Threads.ReentrantLock()
    df_lock = Threads.ReentrantLock()

    Threads.@threads for i in 1:n
        ftp = ftp_paths[i]
        fna_file = Mycelia.download_genome_by_ftp(ftp = ftp, outdir = outdir)
        Threads.lock(df_lock) do
            push!(results, (ftp_path = ftp, fna_path = fna_file))
        end
        Threads.lock(prog_lock) do
            ProgressMeter.next!(p)
        end
    end

    return results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists available BLAST databases from the specified source.
"""
function list_blastdbs(; source::String = "")
    Mycelia.add_bioconda_env("blast")
    @assert source in Set(["", "ncbi", "aws", "gcp"])
    blast_table_data = readlines(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl --showall tsv`)
    # filter out "Connected to AWS" etc...
    if occursin(r"^Connected to"i, blast_table_data[1])
        io = IOBuffer(join(blast_table_data[2:end], '\n'))
    else
        io = IOBuffer(join(blast_table_data, '\n'))
    end
    data, header = uCSV.read(io, delim = '\t', typedetectrows = 100, comment = "#")
    header = ["NAME", "DESCRIPTION", "SIZE (GB)", "LAST_UPDATED"]
    blast_database_table = DataFrames.DataFrame(data, header)
    blast_database_table[!, "LAST_UPDATED"] = map(
        x -> Dates.Date(first(split(x, "T")), Dates.DateFormat("yyyy-mm-dd")),
        blast_database_table[!, "LAST_UPDATED"])
    return blast_database_table
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Display available BLAST databases from specified source.

# # Arguments
# - `source::String="ncbi"`: Database source (default: "ncbi")

# # Returns
# - `DataFrame`: Table of available BLAST databases with columns:
#   - NAME: Database name
#   - SIZE (GB): Database size
#   - LAST_UPDATED: Update date
#   - and other metadata fields
# """
# function showall_blastdbs(;source="ncbi")
#     try
#         run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
#     catch
#         run(`apt-get install ncbi-blast+ perl-doc -y`)
#     end
#     blast_table_header = filter(!isempty, split(readlines(`update_blastdb --source $(source) --showall pretty`)[2], "  "))
#     data, header = uCSV.read(IOBuffer(join(readlines(`update_blastdb --source $(source) --showall tsv`)[2:end], "\n")), delim="\t")
#     df = sort(DataFrames.DataFrame(data, blast_table_header), "SIZE (GB)", rev=true)
#     df[!, "LAST_UPDATED"] = map(dt -> Dates.Date(Dates.DateTime(first(split(dt, '.')), "yyyy-mm-ddTHH:MM:SS")), df[!, "LAST_UPDATED"])
#     return df
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Query information about local BLAST databases and return a formatted summary.

# Arguments
- `blastdbs_dir::String`: Directory containing BLAST databases (default: "~/workspace/blastdb")

# Returns
- `DataFrame` with columns:
  - BLAST database path
  - BLAST database molecule type
  - BLAST database title
  - date of last update
  - number of bases/residues
  - number of sequences
  - number of bytes
  - BLAST database format version
  - human readable size

# Dependencies
Requires NCBI BLAST+ tools. Will attempt to install via apt-get if not present.

# Side Effects
- May install system packages (ncbi-blast+, perl-doc) using sudo/apt-get
- Filters out numbered database fragments from results
"""
function local_blast_database_info(; blastdbs_dir = Mycelia.DEFAULT_BLASTDB_PATH)
    Mycelia.add_bioconda_env("blast")
    # try
    #     run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    # catch
    #     run(`apt-get install ncbi-blast+ perl-doc -y`)
    # end
    # %f means the BLAST database absolute file name path
    # %p means the BLAST database molecule type
    # %t means the BLAST database title
    # %d means the date of last update of the BLAST database
    # %l means the number of bases/residues in the BLAST database
    # %n means the number of sequences in the BLAST database
    # %U means the number of bytes used by the BLAST database
    # %v means the BLAST database format version
    symbol_header_map = OrderedCollections.OrderedDict(
        "%f" => "BLAST database path",
        "%p" => "BLAST database molecule type",
        "%t" => "BLAST database title",
        "%d" => "date of last update",
        "%l" => "number of bases/residues",
        "%n" => "number of sequences",
        "%U" => "number of bytes",
        "%v" => "BLAST database format version"
    )
    outfmt_string = join(collect(keys(symbol_header_map)), "\t")
    data,
    header = uCSV.read(
        open(`$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -list $(blastdbs_dir) -list_outfmt $(outfmt_string)`),
        delim = '\t')
    header = collect(values(symbol_header_map))
    if isempty(data)
        # Return empty DataFrame with expected columns to avoid DimensionMismatch when no DBs are present
        empty_cols = Dict(Symbol(h) => String[] for h in header)
        df = DataFrames.DataFrame(empty_cols)
    else
        df = DataFrames.DataFrame(data, header)
    end
    # remove numbered database fragments from summary results
    df = df[map(x -> !occursin(r"\.\d+$", x), df[!, "BLAST database path"]), :]
    df[!, "human readable size"] = Base.format_bytes.(df[!, "number of bytes"])
    return df
end

function _blastdb_url_reachable(url::AbstractString; connect_timeout::Integer = 5, max_time::Integer = 10)
    cmd = `$(CONDA_RUNNER) run -n blast curl -sfI --connect-timeout $(connect_timeout) --max-time $(max_time) $(url)`
    return success(cmd)
end

function _blastdb_auto_source()
    ncbi_metadata_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/blastdb-metadata-1-1.json"
    if _blastdb_url_reachable(ncbi_metadata_url)
        return ""
    end

    gcp_latest_dir_url = "https://storage.googleapis.com/blast-db/latest-dir"
    if _blastdb_url_reachable(gcp_latest_dir_url)
        @info "NCBI BLAST metadata unreachable; falling back to GCP source"
        return "gcp"
    end

    aws_latest_dir_url = "http://s3.amazonaws.com/ncbi-blast-databases/latest-dir"
    if _blastdb_url_reachable(aws_latest_dir_url)
        @info "NCBI BLAST metadata unreachable; falling back to AWS source"
        return "aws"
    end

    return ""
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Smart downloading of blast dbs depending on interactive, non interactive context

For a list of all available databases, run: `Mycelia.list_blastdbs()`

Downloads and sets up BLAST databases from various sources.

# Arguments
- `db`: Name of the BLAST database to download
- `dbdir`: Directory to store the downloaded database (default: "~/workspace/blastdb")
- `source`: Download source - one of ["", "aws", "gcp", "ncbi"]. Empty string auto-detects fastest source
- `wait`: Whether to wait for download completion (default: true)

# Returns
- String path to the downloaded database directory
"""
function download_blast_db(; db, dbdir = Mycelia.DEFAULT_BLASTDB_PATH, source = "", wait = true)
    Mycelia.add_bioconda_env("blast")
    # try
    #     run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    # catch
    #     run(`apt-get install ncbi-blast+ perl-doc -y`)
    # end
    @assert source in ["", "aws", "gcp", "ncbi"]
    mkpath(dbdir)
    current_directory = pwd()
    try
        cd(dbdir)
        if isempty(source)
            source = _blastdb_auto_source()
            if isempty(source)
                @info "source not provided, letting blast auto-detect fastest download option"
                cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress`
            else
                @info "downloading from source $(source)"
                cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source)`
            end
        else
            @info "downloading from source $(source)"
            if source == "ncbi"
                # --timeout 360 --passive no 
                cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source) --timeout 360 --passive no`
                # cmd = `update_blastdb --timeout 360 --passive no --decompress --source $(source) $(db)`
            else
                cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source)`
                # cmd = `update_blastdb --decompress --source $(source) $(db)`
            end
        end
        run(cmd, wait = wait)
    finally
        cd(current_directory)
    end
    return "$(dbdir)/$(db)"
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Load and parse the assembly summary metadata from NCBI's FTP server for either GenBank or RefSeq databases.

# # Arguments
# - `db::String`: Database source, must be either "genbank" or "refseq"

# # Returns
# - `DataFrame`: Parsed metadata table with properly typed columns including:
#   - Integer columns: taxid, species_taxid, genome metrics, and gene counts
#   - Float columns: gc_percent
#   - Date columns: seq_rel_date, annotation_date
#   - String columns: all other fields

# # Details
# Downloads the assembly summary file from NCBI's FTP server and processes it by:
# 1. Parsing the tab-delimited file with commented headers
# 2. Converting numeric strings to proper Integer/Float types
# 3. Parsing date strings to Date objects
# 4. Handling missing values throughout
# """
# function load_ncbi_metadata(db)
#     if !(db in ["genbank", "refseq"])
#         error()
#     end
#     ncbi_summary_url = "https://ftp.ncbi.nih.gov/genomes/$(db)/assembly_summary_$(db).txt"
#     # ncbi_summary_file = basename(ncbi_summary_url)
#     # if !isfile(ncbi_summary_file)
#     #     download(ncbi_summary_url, ncbi_summary_file)
#     # end
#     buffer = IOBuffer(HTTP.get(ncbi_summary_url).body)
#     # types=[]
#     # ncbi_summary_table = DataFrames.DataFrame(uCSV.read(ncbi_summary_file, comment = "## ", header=1, delim='\t', encodings=Dict("na" => missing), allowmissing=true, typedetectrows=100)...)
#     ncbi_summary_table = DataFrames.DataFrame(uCSV.read(buffer, comment = "## ", header=1, delim='\t', types=String)...)
#     ints = [
#         "taxid",
#         "species_taxid",
#         "genome_size",
#         "genome_size_ungapped",
#         "replicon_count",
#         "scaffold_count",
#         "contig_count",
#         "total_gene_count",
#         "protein_coding_gene_count",
#         "non_coding_gene_count"
#     ]
#     floats = ["gc_percent"]
#     dates = ["seq_rel_date", "annotation_date"]
#     for int in ints
#         ncbi_summary_table[!, int] = something.(tryparse.(Int, ncbi_summary_table[!, int]), missing)
#     end
#     for float in floats
#         ncbi_summary_table[!, float] = something.(tryparse.(Float64, ncbi_summary_table[!, float]), missing)
#     end
#     for date in dates
#         # ncbi_summary_table[!, date] = Dates.Date.(ncbi_summary_table[!, date], Dates.dateformat"yyyy/mm/dd")
#         parsed_dates = map(date_string -> tryparse(Dates.Date, date_string, Dates.dateformat"yyyy/mm/dd"), ncbi_summary_table[!, date])
#         ncbi_summary_table[!, date] = something.(parsed_dates, missing)
#     end
#     return ncbi_summary_table
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Load and parse the assembly summary metadata from NCBI's FTP server for either GenBank or RefSeq databases using CSV.jl for improved performance.

# # Arguments
# - `db::String`: Database source, must be either "genbank" or "refseq".

# # Returns
# - `DataFrames.DataFrame`: Parsed metadata table with columns typed according to NCBI specifications. Handles missing values represented by "na" or empty fields.

# # Details
# Downloads the assembly summary file directly from NCBI's FTP server via HTTP and parses it efficiently using `CSV.File`.
# 1.  Fetches data directly from the URL stream.
# 2.  Skips the initial comment line (`## ...`).
# 3.  Parses the tab-delimited file, recognizing the header line starting with `# assembly_accession...`.
# 4.  Uses `CSV.jl`'s type detection and specific type mapping for performance.
# 5.  Automatically handles missing values ("na", "", "NA").
# 6.  Parses specified date columns using the format "yyyy/mm/dd".
# """
# function load_ncbi_metadata(db::String)
#     # Validate database input
#     if !(db in ["genbank", "refseq"])
#         throw(ArgumentError("Invalid database specified: '$db'. Must be 'genbank' or 'refseq'."))
#     end

#     ncbi_summary_url = "https://ftp.ncbi.nih.gov/genomes/$(db)/assembly_summary_$(db).txt"
#     @info "Fetching NCBI assembly summary from $ncbi_summary_url"

#     try
#         # Fetch data using HTTP.get, response body is an IO stream
#         response = HTTP.get(ncbi_summary_url, status_exception = true) # Throw error for non-2xx status

#         # Define column types for direct parsing by CSV.jl
#         # Let CSV.jl infer non-specified columns (mostly String)
#         # Specify Int, Float, and Date types for known columns
#         types_dict = Dict(
#             :taxid => Int64,
#             :species_taxid => Int64,
#             :genome_size => Int64,
#             :genome_size_ungapped => Int64,
#             :replicon_count => Int64,
#             :scaffold_count => Int64,
#             :contig_count => Int64,
#             :total_gene_count => Int64,
#             :protein_coding_gene_count => Int64,
#             :non_coding_gene_count => Int64,
#             :gc_percent => Float64,
#             :seq_rel_date => Dates.Date,
#             :annotation_date => Dates.Date
#         )

#         # Parse the stream directly using CSV.File
#         # skipto=2: Skip the first line starting with "##"
#         # header=1: The first line *after* skipping is the header (starts with '#')
#         # delim='\t': Tab-separated values
#         # missingstrings: Define strings that represent missing data
#         # types: Apply specific types for efficiency; others inferred
#         # dateformat: Specify the format for date parsing
#         # pool=true: Can improve performance for string columns with repeated values
#         csv_file = CSV.File(
#             response.body;
#             skipto=2,
#             header=1,
#             delim='\t',
#             missingstrings=["na", "NA", ""],
#             types=types_dict,
#             dateformat="yyyy/mm/dd",
#             pool=true,
#             # normalizenames=true could be useful if header names have tricky characters,
#             # but NCBI headers seem clean after the initial '#'.
#             # CSV.jl handles the leading '#' in the header line automatically.
#         )

#         # Materialize the CSV.File into a DataFrame
#         ncbi_summary_table = DataFrames.DataFrame(csv_file)
#         @info "Successfully loaded and parsed NCBI assembly summary into a DataFrame."

#         return ncbi_summary_table

#     catch e
#         @error "Failed to download or parse NCBI metadata from $ncbi_summary_url" exception=(e, catch_backtrace())
#         # Depending on desired behavior, you might rethrow, return an empty DataFrame, or handle specific errors (e.g., HTTP.StatusError)
#         rethrow(e)
#     end
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Load and parse the assembly summary metadata from NCBI's FTP server for either GenBank or RefSeq databases using CSV.jl for improved performance.

# Handles files with initial comment lines (`##...`) and header lines starting with '#'.

# # Arguments
# - `db::String`: Database source, must be either "genbank" or "refseq".

# # Returns
# - `DataFrames.DataFrame`: Parsed metadata table with columns typed according to NCBI specifications. Handles missing values represented by "na" or empty fields.

# # Details
# Downloads the assembly summary file directly from NCBI's FTP server via HTTP and parses it efficiently using `CSV.File`.
# 1.  Fetches the entire data into an in-memory buffer.
# 2.  Manually reads the first line (comment `##`) and second line (header `# ...`).
# 3.  Cleans the extracted header names (removes leading '#', splits by tab).
# 4.  Parses the rest of the buffer using `CSV.File`, providing the cleaned header.
# 5.  Uses `CSV.jl`'s type detection and specific type mapping for performance.
# 6.  Automatically handles missing values ("na", "", "NA") using `missingstring`.
# 7.  Parses specified date columns using the format "yyyy/mm/dd".
# """
# function load_ncbi_metadata(db::String)
#     # Validate database input
#     if !(db in ["genbank", "refseq"])
#         throw(ArgumentError("Invalid database specified: '$db'. Must be 'genbank' or 'refseq'."))
#     end

#     ncbi_summary_url = "https://ftp.ncbi.nih.gov/genomes/$(db)/assembly_summary_$(db).txt"
#     @info "Fetching NCBI assembly summary from $ncbi_summary_url"

#     local response_body::Vector{UInt8}
#     try
#         # Fetch data using HTTP.get
#         response = HTTP.get(ncbi_summary_url, status_exception=true) # Throw error for non-2xx status
#         response_body = response.body
#     catch e
#         @error "Failed to download NCBI metadata from $ncbi_summary_url" exception=(e, catch_backtrace())
#         rethrow(e)
#     end

#     try
#         # Create an IOBuffer from the downloaded body
#         buffer = IOBuffer(response_body)

#         # Read and discard the first line (## comment)
#         readline(buffer)

#         # Read the second line (header starting with #)
#         header_line_raw = readline(buffer)

#         # Clean the header line: remove leading '#', strip whitespace, split by tab
#         header_string = lstrip(header_line_raw, ['#', ' ']) # Remove leading '#' and potential space
#         # Using Base Julia string split:
#         # header_names = String.(Base.split(Base.strip(header_string), '\t'))
#         # Using StringManipulation for potentially more robust splitting/cleaning:
#         header_names = split(strip(header_string), '\t')

#         # Convert header names to Symbols for CSV.jl and DataFrames
#         header_symbols = Symbol.(header_names)

#         # Define column types for direct parsing by CSV.jl
#         types_dict = Dict(
#             :taxid => Int,
#             :species_taxid => Int,
#             :genome_size => Int,
#             :genome_size_ungapped => Int,
#             :replicon_count => Int,
#             :scaffold_count => Int,
#             :contig_count => Int,
#             :total_gene_count => Int,
#             :protein_coding_gene_count => Int,
#             :non_coding_gene_count => Int,
#             :gc_percent => Float64,
#             :seq_rel_date => String,
#             :annotation_date => String
#         )

#         # Parse the rest of the buffer using CSV.File
#         # The buffer is now positioned at the start of the data (line 3)
#         # Provide the cleaned header symbols explicitly
#         # Use 'missingstring' (singular) instead of 'missingstrings' (plural)
#         csv_file = CSV.File(
#             buffer; # Pass the buffer, already positioned after the header
#             header=header_symbols, # Provide the cleaned header names
#             delim='\t',
#             missingstring=["na", "NA", ""], # Corrected keyword
#             types=types_dict,
#             pool=true,
#             # No need for skipto, datarow, or numerical header argument now
#         )

#         # Materialize the CSV.File into a DataFrame
#         ncbi_summary_table = DataFrames.DataFrame(csv_file) # copycols=true is default but explicit
#         @info "Successfully loaded and parsed NCBI assembly summary into a DataFrame."

#         return ncbi_summary_table

#     catch e
#         @error "Failed to parse NCBI metadata buffer from $ncbi_summary_url" exception=(e, catch_backtrace())
#         # Rethrow the parsing error
#         rethrow(e)
#     end
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Load and parse NCBI assembly summary metadata (GenBank/RefSeq), using a daily cache.

# Checks for `homedir()/workspace/.ncbi/YYYY-MM-DD.assembly_summary_{db}.txt`.
# Uses the cache if valid (exists, readable, not empty). Otherwise, downloads
# from NCBI, caches the result, and then parses.

# Handles NCBI's header format and uses CSV.jl for parsing. Requires necessary
# modules like Logging and Base.Filesystem to be in scope (e.g., via `using`).

# # Arguments
# - `db::String`: Database source ("genbank" or "refseq").

# # Returns
# - `DataFrames.DataFrame`: Parsed metadata table.

# # Errors
# - Throws `ArgumentError` for invalid `db`.
# - Throws error if cache directory cannot be created.
# - Throws error if data cannot be obtained from cache or download.
# - Rethrows errors from HTTP download or CSV parsing.
# """
# function load_ncbi_metadata(db::String)
#     # Validate database input
#     if !(db in ["genbank", "refseq"])
#         throw(ArgumentError("Invalid database specified: '$db'. Must be 'genbank' or 'refseq'."))
#     end

#     # --- Cache Path Setup ---
#     todays_date_str = Dates.format(Dates.today(), "yyyy-mm-dd")
#     original_filename = "assembly_summary_$(db).txt"
#     cache_dir = joinpath(homedir(), "workspace", ".ncbi")
#     cached_filename = "$(todays_date_str).$(original_filename)"
#     cached_filepath = joinpath(cache_dir, cached_filename)

#     # Ensure cache directory exists
#     try
#         mkpath(cache_dir)
#     catch e
#         @error "Failed to create cache directory at $cache_dir" exception=(e, catch_backtrace())
#         rethrow(e)
#     end

#     # --- Data Loading Logic ---
#     local response_body::Union{Vector{UInt8}, Nothing} = nothing
#     source_description = ""

#     # 1. Attempt to load from cache
#     if isfile(cached_filepath)
#         @info "Found cached file for today. Attempting to load: $cached_filepath"
#         try
#             content = read(cached_filepath)
#             if isempty(content)
#                 @warn "Cached file is empty: $cached_filepath. Will attempt download."
#             else
#                 response_body = content
#                 source_description = "cache file: $cached_filepath"
#                 @info "Successfully loaded non-empty data from cache."
#             end
#         catch e
#             @error "Failed to read cached file: $cached_filepath. Attempting download." exception=(e, catch_backtrace())
#         end
#     else
#         @info "No cached file found for today at: $cached_filepath. Attempting download."
#     end

#     # 2. If cache loading failed, download
#     if isnothing(response_body)
#         ncbi_summary_url = "https://ftp.ncbi.nih.gov/genomes/$(db)/assembly_summary_$(db).txt"
#         @info "Fetching NCBI assembly summary from $ncbi_summary_url"

#         try
#             response = HTTP.get(ncbi_summary_url; status_exception=true)
#             response_body = response.body
#             source_description = "NCBI URL: $ncbi_summary_url"

#             try
#                 open(cached_filepath, "w") do io
#                     write(io, response_body)
#                 end
#                 @info "Successfully cached downloaded data to: $cached_filepath"
#             catch e
#                 @warn "Failed to write cache file to $cached_filepath. Proceeding with in-memory data." exception=(e, catch_backtrace())
#             end
#         catch e
#             @error "Failed to download NCBI metadata from $ncbi_summary_url" exception=(e, catch_backtrace())
#         end
#     end

#     # 3. Final Check: Ensure data was loaded
#     if isnothing(response_body)
#         @error "Failed to obtain data from both cache and download for database '$db'."
#         error("Could not load NCBI metadata for '$db'. Check connection and cache permissions.")
#     end

#     # --- Parsing Logic ---
#     @info "Proceeding to parse data obtained from $(source_description)."
#     try
#         buffer = IOBuffer(response_body)
#         readline(buffer) # Skip ## comment line
#         header_line_raw = readline(buffer) # Read # header line

#         header_string = lstrip(header_line_raw, ['#', ' '])
#         header_names = split(strip(header_string), '\t')
#         header_symbols = Symbol.(header_names)

#         # --- CORRECTED and EXPANDED types_dict ---
#         types_dict = Dict(
#             # Strings (already present or added based on error)
#             :assembly_accession => String,
#             :bioproject => String,
#             :biosample => String,
#             :wgs_master => String,
#             :refseq_category => String,
#             :organism_name => String,
#             :infraspecific_name => String,
#             :isolate => String,
#             :version_status => String,
#             :assembly_level => String,
#             :release_type => String,
#             :genome_rep => String,
#             :seq_rel_date => String, # Keep as String unless specific parsing needed
#             :asm_name => String,
#             :asm_submitter => String,       # CORRECTED from :submitter
#             :gbrs_paired_asm => String,
#             :paired_asm_comp => String,
#             :ftp_path => String,
#             :excluded_from_refseq => String,
#             :relation_to_type_material => String,
#             :asm_not_live_date => String, # Keep as String unless specific parsing needed
#             :assembly_type => String,      # Added based on error output
#             :group => String,              # Added based on error output
#             :annotation_provider => String,# Added based on error output
#             :annotation_name => String,    # Added based on error output
#             :annotation_date => String,    # Added based on error output
#             :pubmed_id => String,          # Added based on error output (safer as String)

#             # Integers (already present or added based on error)
#             :taxid => Int,
#             :species_taxid => Int,
#             :genome_size => Int,           # Added based on error output
#             :genome_size_ungapped => Int,  # Added based on error output
#             :replicon_count => Int,        # Added based on error output
#             :scaffold_count => Int,        # Added based on error output
#             :contig_count => Int,          # Added based on error output
#             :total_gene_count => Int,      # Added based on error output
#             :protein_coding_gene_count => Int, # Added based on error output
#             :non_coding_gene_count => Int, # Added based on error output

#             # Floats (added based on error output)
#             :gc_percent => Float64         # Added based on error output
#         )
#         # --- End of types_dict ---

#         # Validate that headers derived match the keys expected (optional but good practice)
#         # This helps catch discrepancies early if NCBI changes format
#         if Set(header_symbols) != Set(keys(types_dict))
#              missing_in_dict = setdiff(Set(header_symbols), Set(keys(types_dict)))
#              extra_in_dict = setdiff(Set(keys(types_dict)), Set(header_symbols))
#              if !isempty(missing_in_dict)
#                  @warn "Headers found in data but not in types_dict: $missing_in_dict"
#              end
#              # Don't warn about extra_in_dict usually, as the error already caught the critical case.
#              # CSV.jl handles extra columns by inferring types.
#         end

#         # Parse using CSV.File
#         csv_file = CSV.File(
#             buffer;
#             header=header_symbols,
#             delim='\t',
#             missingstring=["na", "NA", ""],
#             types=types_dict,
#             pool=true,
#             # validate=false # Avoid using this; fixing types_dict is better
#         )

#         # Materialize into a DataFrame
#         ncbi_summary_table = DataFrames.DataFrame(csv_file; copycols=true)
#         @info "Successfully parsed NCBI assembly summary into a DataFrame."

#         return ncbi_summary_table

#     catch e
#         @error "Failed to parse NCBI metadata buffer from $(source_description)" exception=(e, catch_backtrace())
#         rethrow(e) # Rethrow parsing error
#     end
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse NCBI assembly summary metadata from a local file.

This function handles the NCBI file format (comment line, header line, then data)
and parses it into a DataFrame with appropriate column types.

# Arguments
- `filepath::String`: Path to the NCBI assembly summary file to parse.

# Returns
- `DataFrames.DataFrame`: Parsed metadata table.

# Errors
- Throws error if file doesn't exist or cannot be read.
- Rethrows errors from CSV parsing.
"""
function parse_ncbi_metadata_file(filepath::String)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    if filesize(filepath) == 0
        error("File is empty: $filepath")
    end

    @info "Parsing NCBI metadata from: $filepath"

    try
        # Manually read the first two lines to get the correct header
        header_symbols = Symbol[]
        open(filepath, "r") do io
            readline(io) # Skip ## comment line
            header_line_raw = readline(io) # Read # header line
            header_string = lstrip(header_line_raw, ['#', ' '])
            header_names = split(strip(header_string), '\t')
            header_symbols = Symbol.(header_names)
        end # File is automatically closed here

        # Define column types
        types_dict = Dict(
            :assembly_accession => String, :bioproject => String, :biosample => String,
            :wgs_master => String, :refseq_category => String, :organism_name => String,
            :infraspecific_name => String, :isolate => String, :version_status => String,
            :assembly_level => String, :release_type => String, :genome_rep => String,
            :seq_rel_date => String, :asm_name => String, :asm_submitter => String,
            :gbrs_paired_asm => String, :paired_asm_comp => String, :ftp_path => String,
            :excluded_from_refseq => String, :relation_to_type_material => String,
            :asm_not_live_date => String, :assembly_type => String, :group => String,
            :annotation_provider => String, :annotation_name => String, :annotation_date => String,
            :pubmed_id => String,
            :taxid => Int, :species_taxid => Int, :genome_size => Int,
            :genome_size_ungapped => Int, :replicon_count => Int, :scaffold_count => Int,
            :contig_count => Int, :total_gene_count => Int, :protein_coding_gene_count => Int,
            :non_coding_gene_count => Int,
            :gc_percent => Float64
        )

        # Validate that headers derived match the keys expected
        if Set(header_symbols) != Set(keys(types_dict))
            missing_in_dict = setdiff(Set(header_symbols), Set(keys(types_dict)))
            if !isempty(missing_in_dict)
                @warn "Headers found in data file but not explicitly typed in types_dict (will be inferred by CSV.jl): $missing_in_dict"
            end
        end

        # Parse using CSV.File directly from the filepath
        csv_file = CSV.File(
            filepath;
            skipto = 3,             # Skip the comment and header lines we already processed
            header = header_symbols, # Provide the headers we extracted
            delim = '\t',
            missingstring = ["na", "NA", ""],
            types = types_dict,
            pool = true,
            strict = true
        )

        # Materialize into a DataFrame
        ncbi_summary_table = DataFrames.DataFrame(csv_file; copycols = true)
        @info "Successfully parsed NCBI assembly summary into a DataFrame from $filepath"

        return ncbi_summary_table

    catch e
        @error "Failed to parse NCBI metadata file: $filepath" exception=(
            e, catch_backtrace())
        rethrow(e)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load and parse NCBI assembly summary metadata (GenBank/RefSeq), using a daily cache.

Checks for `homedir()/workspace/.ncbi/YYYY-MM-DD.assembly_summary_{db}.txt`.
Uses the cache if valid (exists, readable, not empty). Otherwise, downloads
from NCBI using `Downloads.download()`, caches the result (replacing any
previous version for the *same day*), and then parses the cached file.

# Arguments
- `db::String`: Database source ("genbank" or "refseq").
- `date::Union{Dates.Date,Nothing}=nothing`: Specific date to load from cache. 
  If `nothing`, uses today's date. If specified, only uses cache (no download).

# Returns
- `DataFrames.DataFrame`: Parsed metadata table.

# Errors
- Throws `ArgumentError` for invalid `db`.
- Throws error if cache directory cannot be created.
- Throws error if data cannot be obtained from cache or download.
- Throws error if specific date is requested but not found in cache.
- Rethrows errors from `Downloads.download` or CSV parsing.
"""
function load_ncbi_metadata(db::String; date::Union{Dates.Date, Nothing} = nothing)
    # Validate database input
    if !(db in ["genbank", "refseq"])
        throw(ArgumentError("Invalid database specified: '$db'. Must be 'genbank' or 'refseq'."))
    end

    # --- Cache Path Setup ---
    target_date = isnothing(date) ? Dates.today() : date
    date_str = Dates.format(target_date, "yyyy-mm-dd")
    original_filename = "assembly_summary_$(db).txt"
    cache_dir = joinpath(homedir(), "workspace", ".ncbi")
    cached_filename = "$(date_str).$(original_filename)"
    cached_filepath = joinpath(cache_dir, cached_filename)
    ncbi_summary_url = "https://ftp.ncbi.nih.gov/genomes/$(db)/assembly_summary_$(db).txt"

    # Ensure cache directory exists
    try
        mkpath(cache_dir)
    catch e
        @error "Failed to create cache directory at $cache_dir" exception=(
            e, catch_backtrace())
        rethrow(e)
    end

    # --- Data Acquisition Logic ---
    local source_description::String = ""
    data_needs_download = false

    # If a specific date was requested, only use cache (no download)
    if !isnothing(date)
        if isfile(cached_filepath) && filesize(cached_filepath) > 0
            @info "Found cached file for requested date ($date_str): $cached_filepath"
            source_description = "cache file: $cached_filepath"
        else
            error("No valid cache file found for requested date $date_str at: $cached_filepath. Cannot download when specific date is requested.")
        end
    else
        # Using today's date - can attempt download if needed
        # 1. Attempt to use existing cache
        if isfile(cached_filepath)
            if filesize(cached_filepath) > 0
                @info "Found valid cached file for today: $cached_filepath"
                source_description = "cache file: $cached_filepath"
            else
                @warn "Cached file exists but is empty: $cached_filepath. Will attempt download to replace it."
                data_needs_download = true
            end
        else
            @info "No cached file found for today at: $cached_filepath. Attempting download."
            data_needs_download = true
        end

        # 2. If cache is missing or empty, download
        if data_needs_download
            @info "Downloading NCBI assembly summary from $ncbi_summary_url to $cached_filepath"
            try
                Downloads.download(ncbi_summary_url, cached_filepath)
                # Verify download success
                if isfile(cached_filepath) && filesize(cached_filepath) > 0
                    @info "Successfully downloaded and cached data to: $cached_filepath"
                    source_description = "downloaded file: $cached_filepath (from $ncbi_summary_url)"
                else
                    @error "Download completed but resulted in an empty or missing file at $cached_filepath."
                end
            catch e
                @error "Failed to download NCBI metadata from $ncbi_summary_url to $cached_filepath" exception=(
                    e, catch_backtrace())
            end
        end

        # 3. Final Check: Ensure a valid data file exists at the cached path
        if !isfile(cached_filepath) || filesize(cached_filepath) == 0
            @error "Failed to obtain valid data for '$db'. Could not use cache and download failed or resulted in empty file."
            error("Could not load NCBI metadata for '$db' from cache or download. Check path '$cached_filepath', network connection, and permissions.")
        end
    end

    # --- Parse the file ---
    return parse_ncbi_metadata_file(cached_filepath)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Loads NCBI RefSeq metadata into a DataFrame. RefSeq is NCBI's curated collection 
of genomic, transcript and protein sequences.

# Arguments
- `date::Union{Dates.Date,Nothing}=nothing`: Specific date to load from cache.
  If `nothing`, uses today's date. If specified, only uses cache (no download).

# Returns
- `DataFrames.DataFrame`: Contains metadata columns including accession numbers, taxonomic information,
and sequence details from RefSeq.
"""
function load_refseq_metadata(; date::Union{Dates.Date, Nothing} = nothing)
    return load_ncbi_metadata("refseq"; date = date)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load metadata for GenBank sequences into a DataFrame.

This is a convenience wrapper around `load_ncbi_metadata("genbank")` that
specifically loads metadata from the GenBank database.

# Arguments
- `date::Union{Dates.Date,Nothing}=nothing`: Specific date to load from cache.
  If `nothing`, uses today's date. If specified, only uses cache (no download).

# Returns
- `DataFrames.DataFrame`: Contains metadata fields like accession numbers, taxonomy,
and sequence information from GenBank.
"""
function load_genbank_metadata(; date::Union{Dates.Date, Nothing} = nothing)
    return load_ncbi_metadata("genbank"; date = date)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a complete NCBI FTP URL by combining a base FTP path with a file extension.

# Arguments
- `ftp_path::String`: Base FTP directory path for the resource
- `extension::String`: File extension to append to the resource name

# Returns
- `String`: Complete FTP URL path to the requested resource

Extensions include:
- genomic.fna.gz
- genomic.gff.gz
- protein.faa.gz
- assembly_report.txt
- assembly_stats.txt
- cds_from_genomic.fna.gz
- feature_count.txt.gz
- feature_table.txt.gz
- genomic.gbff.gz
- genomic.gtf.gz
- protein.gpff.gz
- translated_cds.faa.gz
"""
function ncbi_ftp_path_to_url(; ftp_path, extension)
    # genomic.fna.gz
    # genomic.gff.gz
    # protein.faa.gz
    # assembly_report.txt
    # assembly_stats.txt
    # cds_from_genomic.fna.gz
    # feature_count.txt.gz
    # feature_table.txt.gz
    # genomic.gbff.gz
    # genomic.gtf.gz
    # protein.gpff.gz
    # translated_cds.faa.gz    
    f_name = basename(ftp_path) * "_" * extension
    new_path = joinpath(ftp_path, f_name)
    return new_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

Retrieve FASTA format sequences from NCBI databases or direct FTP URLs.

# Arguments
- `db::String`: NCBI database type ("nuccore" for DNA or "protein" for protein sequences)
- `accession::String`: NCBI sequence accession number
- `ftp::String`: Direct FTP URL to a FASTA file (alternative to db/accession pair)

# Returns
- `FASTX.FASTA.Reader`: Reader object containing the requested sequence(s)
"""
function get_sequence(; db::String = "", accession::String = "", ftp::String = "")
    if !isempty(db) && !isempty(accession)
        sleep(0.34)
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=fasta&id=$(accession)"
        resp = HTTP.get(url)
        if resp.status != 200
            throw(error("HTTP request failed (status $(resp.status)) for URL: $url"))
        end
        try
            return collect(FASTX.FASTA.Reader(IOBuffer(resp.body)))
        catch e
            response_string = read(IOBuffer(resp.body), String)
            throw(error("Failed to parse FASTA from NCBI: $(e)\n$(response_string)"))
        end
    elseif !isempty(ftp)
        resp = HTTP.get(ftp)
        if resp.status != 200
            throw(error("HTTP request failed (status $(resp.status)) for FTP URL: $ftp"))
        end
        try
            return collect(FASTX.FASTA.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(resp.body))))
        catch e
            throw(error("Failed to decompress or parse FASTA from FTP: $(e)"))
        end
    else
        throw(ArgumentError("Must provide either both `db` and `accession`, or `ftp`"))
    end
end

# function ncbi_datasets_download_by_taxon_id
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and constructs a MetaDiGraph representation of the NCBI taxonomy database.

# Arguments
- `path_to_taxdump`: Directory path where taxonomy files will be downloaded and extracted

# Returns
- `MetaGraphsNext.MetaGraph`: A directed graph where:
  - Vertex labels are `tax_id` values
  - Vertex metadata stores name, rank, division, and genetic code properties
  - Edges represent parent-child relationships in the taxonomy

# Dependencies
Requires internet connection for initial download. Uses DataFrames, MetaGraphsNext, and ProgressMeter.
"""
function load_ncbi_taxonomy(;
        path_to_taxdump = joinpath(Mycelia.DEFAULT_BLASTDB_PATH, "taxdump")
    # path_to_prebuilt_graph="$(path_to_taxdump)/ncbi_taxonomy.jld2"
)
    taxdump_url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    taxdump_local_tarball = joinpath(dirname(path_to_taxdump), basename(taxdump_url))
    taxdump_out = replace(taxdump_local_tarball, ".tar.gz" => "")
    # if isfile(path_to_prebuilt_graph) && filesize(path_to_prebuilt_graph) > 0
    #     println("Using prebuilt graph"
    #     ncbi_taxonomy = JLD2.load(path_to_prebuilt_graph, "ncbi_taxonomy")
    #     return (;ncbi_taxonomy, path_to_prebuilt_graph)
    # end
    if !isdir(taxdump_out)
        mkpath(taxdump_out)
        if !isfile(taxdump_local_tarball)
            download(taxdump_url, taxdump_local_tarball)
        end
        run(`tar -xf $(taxdump_local_tarball) -C $(taxdump_out)`)
    end

    names_dmp = DataFrames.DataFrame(
        tax_id = Int[],
        name_txt = String[],
        unique_name = String[],
        name_class = String[]
    )
    ProgressMeter.@showprogress for line in
                                    split(read(open("$(taxdump_out)/names.dmp"), String), "\t|\n")
        if isempty(line)
            continue
        else
            (tax_id_string, name_txt, unique_name, name_class) = split(line, "\t|\t")
            tax_id = parse(Int, tax_id_string)
            row = (; tax_id, name_txt, unique_name, name_class)
            push!(names_dmp, row)
        end
    end
    ncbi_taxonomy = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph();
        label_type = Int,
        vertex_data_type = Dict{Symbol, Any},
        edge_data_type = Nothing,
        graph_data = Dict{Symbol, Any}()
    )
    ProgressMeter.@showprogress for group in DataFrames.groupby(names_dmp, "tax_id")
        tax_id = group[1, "tax_id"]
        vertex_data = Dict{Symbol, Any}(:tax_id => tax_id)
        for row in DataFrames.eachrow(group)
            unique_name = isempty(row["unique_name"]) ? row["name_txt"] : row["unique_name"]
            # remove quotes since neo4j doesn't like them
            unique_name = replace(unique_name, '"' => "")
            # replace spaces and dashes with underscores
            name_class = Symbol(replace(replace(row["name_class"], r"\s+" => "-"), "-" => "_"))
            if haskey(vertex_data, name_class)
                current_value = vertex_data[name_class]
                if (current_value isa Vector) && !(unique_name in current_value)
                    vertex_data[name_class] = [current_value..., unique_name]
                elseif !(current_value isa Vector) && (current_value != unique_name)
                    vertex_data[name_class] = [current_value, unique_name]
                else
                    continue
                end
            else
                vertex_data[name_class] = unique_name
            end
        end
        ncbi_taxonomy[tax_id] = vertex_data
    end
    divisions = Dict()
    for line in split(read(open("$(taxdump_out)/division.dmp"), String), "\t|\n")
        if !isempty(line)
            (id_string, shorthand, full_name, notes) = split(line, "\t|\t")
            id = parse(Int, id_string)
            divisions[id] = Dict(:division_cde => shorthand, :division_name => full_name)
        end
    end
    divisions

    ProgressMeter.@showprogress for line in
                                    split(read(open("$(taxdump_out)/nodes.dmp"), String), "\t|\n")
        if isempty(line)
            continue
        else
            fields = split(line, "\t|\t")
            @assert length(fields) >= 9
            tax_id_string = fields[1]
            parent_tax_id_string = fields[2]
            rank = fields[3]
            embl_code = fields[4]
            division_id_string = fields[5]
            genetic_code_id_string = fields[7]
            mitochondrial_genetic_code_id_string = fields[9]

            division_id = parse(Int, division_id_string)
            genetic_code_id = parse(Int, genetic_code_id_string)
            mitochondrial_genetic_code_id = parse(Int, mitochondrial_genetic_code_id_string)

            tax_id = parse(Int, tax_id_string)
            parent_tax_id = parse(Int, parent_tax_id_string)
            if !haskey(ncbi_taxonomy, tax_id)
                ncbi_taxonomy[tax_id] = Dict{Symbol, Any}(:tax_id => tax_id)
            end
            if !haskey(ncbi_taxonomy, parent_tax_id)
                ncbi_taxonomy[parent_tax_id] = Dict{Symbol, Any}(:tax_id => parent_tax_id)
            end

            Graphs.add_edge!(ncbi_taxonomy, tax_id, parent_tax_id)
            vertex_data = ncbi_taxonomy[tax_id]
            vertex_data[:rank] = rank
            # these should probably be broken out as independent nodes!
            vertex_data[:division_id] = division_id
            vertex_data[:division_cde] = divisions[division_id][:division_cde]
            vertex_data[:division_name] = divisions[division_id][:division_name]
            vertex_data[:genetic_code_id] = genetic_code_id
            vertex_data[:mitochondrial_genetic_code_id] = mitochondrial_genetic_code_id
        end
    end
    # JLD2 graph killed a colab instance after 200Gb of size!
    # JLD2.save("$(homedir())/workspace/blastdb/taxdump/ncbi_taxonomy.jld2", "ncbi_taxonomy", ncbi_taxonomy)
    # return (;ncbi_taxonomy, path_to_prebuilt_graph)
    return ncbi_taxonomy
end

function _get_ncbi_tax_id_map(ncbi_taxonomy)
    graph_data = ncbi_taxonomy[]
    if (graph_data isa AbstractDict) && haskey(graph_data, :node_2_taxid_map)
        return graph_data[:node_2_taxid_map]
    end
    node_2_taxid_map = sort!(collect(MetaGraphsNext.labels(ncbi_taxonomy)))
    if graph_data isa AbstractDict
        graph_data[:node_2_taxid_map] = node_2_taxid_map
    end
    return node_2_taxid_map
end

function _ncbi_tax_id_to_vertex(ncbi_taxonomy, tax_id::Int)
    if haskey(ncbi_taxonomy, tax_id)
        return tax_id
    end
    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return the NCBI translation table ID for a taxonomy node.

# Arguments
- `ncbi_taxonomy`: Output from `load_ncbi_taxonomy`
- `tax_id`: NCBI taxonomy identifier

# Keywords
- `type::Symbol=:genomic`: `:genomic` or `:mitochondrial`
- `inherit::Bool=true`: Walk up the parent lineage if the code is missing/zero

# Returns
`Int` translation table ID, or `missing` if not found.
"""
function get_ncbi_genetic_code(ncbi_taxonomy, tax_id::Int; type::Symbol = :genomic, inherit::Bool = true)
    prop = if type === :genomic
        :genetic_code_id
    elseif type === :mitochondrial
        :mitochondrial_genetic_code_id
    else
        throw(ArgumentError("type must be :genomic or :mitochondrial"))
    end

    vertex = _ncbi_tax_id_to_vertex(ncbi_taxonomy, tax_id)
    if vertex === nothing
        return missing
    end

    vertex_data = ncbi_taxonomy[vertex]
    code = get(vertex_data, prop, missing)
    if !inherit
        return code
    end

    while code === missing || code == 0
        parents = collect(MetaGraphsNext.outneighbor_labels(ncbi_taxonomy, vertex))
        if isempty(parents)
            return code
        end
        vertex = first(parents)
        vertex_data = ncbi_taxonomy[vertex]
        code = get(vertex_data, prop, missing)
    end

    return code
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Download an accession using NCBI datasets command line tool

the .zip download output to outpath will be unzipped

returns the outfolder

ncbi's default include string is 
include_string = "gff3,rna,cds,protein,genome,seq-report"

Downloads and extracts a genome from NCBI using the datasets command line tool.

# Arguments
- `accession`: NCBI accession number for the genome
- `outdir`: Directory where files will be downloaded (defaults to current directory)
- `outpath`: Full path for the temporary zip file (defaults to `outdir/accession.zip`)
- `include_string`: Data types to download (defaults to all "gff3,rna,cds,protein,genome,seq-report").
- `max_attempts`: Maximum number of download retry attempts (default 3)
- `initial_retry_delay`: Initial delay in seconds between retries (default 10.0)
  
# Returns
- Path to the extracted genome data directory

# Notes
- Requires the ncbi-datasets-cli conda package (automatically installed if missing)
- Downloaded zip file is automatically removed after extraction
- If output folder already exists, download is skipped
- Data is extracted to `outdir/accession/ncbi_dataset/data/accession`
- Uses exponential backoff retry logic to handle transient NCBI API failures

# TODO: Apply similar retry logic to other NCBI download functions:
# - `download_genome_by_accession` in this file
# - `load_ncbi_metadata` in this file  
# - taxonomy functions in `taxonomy-and-trees.jl` that use ncbi-datasets-cli
"""
function ncbi_genome_download_accession(;
        accession,
        outdir = pwd(),
        outpath = joinpath(outdir, accession * ".zip"),
        include_string = "gff3,rna,cds,protein,genome,seq-report",
        max_attempts::Int = 3,
        initial_retry_delay::Float64 = 10.0
)
    outfolder = joinpath(outdir, accession)
    if !isdir(outfolder)
        add_bioconda_env("ncbi-datasets-cli")
        if isfile(outpath)
            @info "$(outpath) already exists, skipping download..."
        else
            mkpath(outdir)
            # Use retry logic to handle transient NCBI API failures
            with_retry(
                max_attempts = max_attempts,
                initial_delay = initial_retry_delay,
                on_retry = (attempt,
                    ex,
                    delay) -> @warn "NCBI download attempt $attempt/$max_attempts for $accession failed, retrying in $(delay)s..." exception=ex
            ) do
                run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome accession $(accession) --include $(include_string) --filename $(outpath) --no-progressbar`)
            end
        end
        run(`unzip -q -d $(outfolder) $(outpath)`)
    end
    final_outfolder = joinpath(outfolder, "ncbi_dataset", "data", accession)
    isfile(outpath) && rm(outpath)
    # Mapping include_string items to file paths
    include_items = split(include_string, ",")
    genome_value = nothing
    cds_value = nothing
    gff3_value = nothing
    protein_value = nothing
    seqreport_value = nothing
    for included_item in include_items
        if included_item == "genome"
            candidate_genomes = filter(
                x -> occursin(accession, basename(x)) &&
                     occursin(Mycelia.FASTA_REGEX, basename(x)),
                readdir(final_outfolder, join = true))
            @assert length(candidate_genomes) == 1
            genome_value = first(candidate_genomes)
            @assert isfile(genome_value)
            @assert filesize(genome_value) > 0
        elseif included_item == "cds"
            expected_cds_value = joinpath(final_outfolder, "cds_from_genomic.fna")
            if isfile(expected_cds_value) && filesize(expected_cds_value) > 0
                cds_value = expected_cds_value
            end
        elseif included_item == "gff3"
            expected_gff3_value = joinpath(final_outfolder, "genomic.gff")
            if isfile(expected_gff3_value) && filesize(expected_gff3_value) > 0
                gff3_value = expected_gff3_value
            end
        elseif included_item == "protein"
            expected_protein_value = joinpath(final_outfolder, "protein.faa")
            if isfile(expected_protein_value) && filesize(expected_protein_value) > 0
                protein_value = expected_protein_value
            end
        elseif included_item == "seq-report"
            expected_seqreport_value = joinpath(final_outfolder, "sequence_report.jsonl")
            if isfile(expected_seqreport_value) && filesize(expected_seqreport_value) > 0
                seqreport_value = expected_seqreport_value
            end
        end
    end
    return (
        directory = final_outfolder,
        genome = genome_value,
        cds = cds_value,
        gff3 = gff3_value,
        protein = protein_value,
        seqreport = seqreport_value
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

Retrieve GenBank records from NCBI or directly from an FTP site.

# Arguments
- `db::String`: NCBI database to query ("nuccore" for nucleotide or "protein" for protein sequences)
- `accession::String`: NCBI accession number for the sequence
- `ftp::String`: Direct FTP URL to a GenBank file (gzipped)

# Returns
- `GenomicAnnotations.GenBank.Reader`: A reader object containing the GenBank record

# Details
When using NCBI queries (`db` and `accession`), the function implements rate limiting 
(0.5s sleep) to comply with NCBI's API restrictions of max 2 requests per second.
"""
function get_genbank(; db = ""::String, accession = ""::String, ftp = ""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        # url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=genbank&id=$(accession)&rettype=text"
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=genbank&id=$(accession)&retmode=text"
        # readgbk can't read from an io buffer, so need to download to a temp file
        # outfile = tempname()
        # open(outfile, "w") do io
        #     write(io, HTTP.get(url).body)
        # end
        # genbank_data = GenomicAnnotations.readgbk(outfile)
        # rm(outfile)
        # return genbank_data
        return GenomicAnnotations.GenBank.Reader(IOBuffer(HTTP.get(url).body))
    elseif !isempty(ftp)
        return GenomicAnnotations.GenBank.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
    else
        @error "invalid call"
    end
end

# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and extracts the NCBI taxonomy database required for taxonkit operations.

Downloads `taxdump.tar.gz` from NCBI FTP server and extracts it to `~/.taxonkit/`.
This is a prerequisite for using taxonkit-based taxonomy functions.

# Arguments
- `force_update::Bool=false`: Force download even if taxdump already exists
- `max_age_days::Int=30`: Maximum age in days before warning about stale data

# Requirements
- Working internet connection
- Sufficient disk space (~100MB)
- `taxonkit` must be installed separately

# Returns
- Nothing

# Throws
- `SystemError` if download fails or if unable to create directory
- `ErrorException` if tar extraction fails
"""
function setup_taxonkit_taxonomy(; force_update::Bool = false, max_age_days::Int = 30)
    taxonkit_dir = joinpath(homedir(), ".taxonkit")
    timestamp_file = joinpath(taxonkit_dir, ".last_updated")

    # Check if we need to download/update
    needs_update = force_update

    if !needs_update
        # Check if directory exists and has files
        if !isdir(taxonkit_dir) || isempty(readdir(taxonkit_dir))
            needs_update = true
        else
            # Check age of installation
            if isfile(timestamp_file)
                try
                    last_updated = Dates.DateTime(read(timestamp_file, String))
                    age_days = Dates.value(Dates.now() - last_updated) /
                               (1000 * 60 * 60 * 24)

                    if age_days > max_age_days
                        @warn "Taxonkit taxonomy data is $(round(age_days, digits=1)) days old (> $max_age_days days). Consider updating with force_update=true."
                    end
                catch
                    @warn "Could not read timestamp file. Taxonkit data age unknown."
                end
            else
                @warn "No timestamp found for taxonkit data. Age unknown. Consider updating with force_update=true."
            end
        end
    end

    # Download and extract if needed
    if needs_update
        @info "Downloading NCBI taxonomy database..."

        # Clean up existing directory if it exists
        if isdir(taxonkit_dir)
            @info "Removing existing taxonkit directory..."
            rm(taxonkit_dir, recursive = true)
        end

        # Create fresh directory
        mkpath(taxonkit_dir)

        # Download in a temporary location first
        temp_file = tempname() * ".tar.gz"
        try
            run(`wget -q -O $temp_file ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`)

            # Extract to taxonkit directory
            Mycelia.tar_extract(tarchive = temp_file, directory = taxonkit_dir)

            # Write timestamp
            write(timestamp_file, string(Dates.now()))

            @info "Taxonkit taxonomy database updated successfully"

        finally
            # Clean up temporary file
            isfile(temp_file) && rm(temp_file)
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load and parse BV-BRC (Bacterial and Viral Bioinformatics Resource Center) genome metadata using a daily cache.

Checks for cached files at `homedir()/workspace/.bvbrc/YYYY-MM-DD.genome_summary.tsv` and
`homedir()/workspace/.bvbrc/YYYY-MM-DD.genome_metadata.tsv`. Uses the cache if valid 
(exists, readable, not empty). Otherwise, downloads from BV-BRC FTP server, caches the 
results, and then parses the cached files. The two tables (genome_summary and genome_metadata) 
are joined on the genome_id column.

# Keyword Arguments
- `date::Union{Dates.Date,Nothing}=nothing`: Specific date to load from cache. 
  If `nothing`, uses today's date. If specified, only uses cache (no download).
- `summary_url::String="ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary"`: URL for genome summary file.
- `metadata_url::String="ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata"`: URL for genome metadata file.

# Returns
- `DataFrames.DataFrame`: Inner join of genome summary and metadata tables on genome_id.

# Errors
- Throws error if cache directory cannot be created.
- Throws error if specific date is requested but not found in cache.
- Throws error if data cannot be obtained from cache or download.
- Rethrows errors from `Downloads.download` or CSV parsing.
"""
function load_bvbrc_genome_metadata(;
        date::Union{Dates.Date, Nothing} = nothing,
        summary_url::String = "ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary",
        metadata_url::String = "ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata"
)
    # --- Cache Path Setup ---
    cache_dir = joinpath(homedir(), "workspace", ".bvbrc")

    try
        mkpath(cache_dir)
    catch e
        @error "Failed to create cache directory at $cache_dir" exception=(
            e, catch_backtrace())
        rethrow(e)
    end

    target_date = isnothing(date) ? Dates.today() : date
    date_str = Dates.format(target_date, "yyyy-mm-dd")
    summary_cache = joinpath(cache_dir, "$(date_str).genome_summary.tsv")
    metadata_cache = joinpath(cache_dir, "$(date_str).genome_metadata.tsv")

    function is_valid_file(f)
        isfile(f) && filesize(f) > 0
    end

    # --- Data Acquisition Logic ---
    summary_file = summary_cache
    metadata_file = metadata_cache

    # If a specific date was requested, only use cache (no download)
    if !isnothing(date)
        if is_valid_file(summary_cache) && is_valid_file(metadata_cache)
            @info "Found cached BVBRC metadata for requested date ($date_str): $summary_cache, $metadata_cache"
        else
            error("No valid cache files found for requested date $date_str at: $summary_cache, $metadata_cache. Cannot download when specific date is requested.")
        end
    else
        # Using today's date - can attempt download if needed
        if is_valid_file(summary_cache) && is_valid_file(metadata_cache)
            @info "Using cached BVBRC metadata for today: $summary_cache, $metadata_cache"
        else
            # Try to download and cache
            try
                @info "Downloading BVBRC genome summary from $summary_url"
                Downloads.download(summary_url, summary_cache)
                @info "Downloading BVBRC genome metadata from $metadata_url"
                Downloads.download(metadata_url, metadata_cache)

                if is_valid_file(summary_cache) && is_valid_file(metadata_cache)
                    @info "Successfully downloaded and cached BVBRC metadata"
                else
                    @error "Downloaded files are missing or empty"
                    error("Download completed but resulted in empty or missing files")
                end
            catch e
                @warn "Failed to download BVBRC metadata" exception=(e, catch_backtrace())

                # Fallback: find most recent previous cache
                cached_summary_files = sort(
                    filter(f -> occursin(r"\d{4}-\d{2}-\d{2}\.genome_summary\.tsv", f),
                        readdir(cache_dir; join = true)),
                    rev = true
                )
                cached_metadata_files = sort(
                    filter(f -> occursin(r"\d{4}-\d{2}-\d{2}\.genome_metadata\.tsv", f),
                        readdir(cache_dir; join = true)),
                    rev = true
                )

                if !isempty(cached_summary_files) && !isempty(cached_metadata_files)
                    summary_file = cached_summary_files[1]
                    metadata_file = cached_metadata_files[1]
                    @info "Falling back to most recent cached BVBRC metadata: $summary_file, $metadata_file"
                else
                    @error "No valid BVBRC metadata cache available and download failed"
                    error("Could not load BVBRC metadata from cache or download. Check network connection and permissions.")
                end
            end
        end
    end

    # --- Parse the files ---
    try
        @info "Reading genome summary file: $summary_file"
        genome_summary = CSV.read(
            summary_file,
            DataFrames.DataFrame,
            delim = '\t',
            header = 1,
            types = Dict("genome_id" => String)
        )

        @info "Reading genome metadata file: $metadata_file"
        genome_metadata = CSV.read(
            metadata_file,
            DataFrames.DataFrame,
            delim = '\t',
            header = 1,
            types = Dict("genome_id" => String)
        )

        # Join the DataFrames
        @info "Joining genome summary and metadata on genome_id"
        bvbrc_genome_summary = DataFrames.innerjoin(
            genome_summary,
            genome_metadata,
            on = "genome_id",
            makeunique = true
        )

        @info "Successfully loaded and joined BVBRC genome metadata"
        return bvbrc_genome_summary

    catch e
        @error "Failed to parse BVBRC metadata files" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Export sequences from a BLAST database to a gzipped FASTA file.

# Arguments
- `path_to_db`: Path to the BLAST database
- `fasta`: Output path for the gzipped FASTA file (default: `path_to_db * ".fna.gz"`)

# Details
Uses conda's BLAST environment to extract sequences using `blastdbcmd`.
The output is automatically compressed using `pigz`.
If the output file already exists, the function will skip extraction.

"""
function export_blast_db(; path_to_db, fasta = path_to_db * ".fna.gz")
    Mycelia.add_bioconda_env("blast")
    if !isfile(fasta)
        # -long_seqids adds GI identifiers - these are cross-referenceable through other means so I'm dropping
        @time run(pipeline(
            pipeline(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd  -entry all -outfmt '%f' -db $(path_to_db)`,
                `pigz`),
            fasta))
    else
        @info "$(fasta) already present"
    end
end

# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/taxonomy/
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Retrieve taxonomic information for a given NCBI taxonomy ID.

# Arguments
- `taxa_id`: NCBI taxonomy identifier (integer)

# Returns
- `DataFrame`: Taxonomy summary containing fields like tax_id, rank, species, etc.
"""
function ncbi_taxon_summary(taxa_id)
    Mycelia.add_bioconda_env("ncbi-datasets")
    p = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets datasets summary taxonomy taxon $(taxa_id) --as-json-lines`,
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets dataformat tsv taxonomy --template tax-summary`
    )
    return DataFrames.DataFrame(uCSV.read(open(p), delim = '\t', header = 1))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Exports a taxonomy mapping table from a BLAST database in seqid2taxid format.

# Arguments
- `path_to_db::String`: Path to the BLAST database
- `outfile::String`: Output file path (defaults to input path + ".seqid2taxid.txt.gz")

# Returns
- `String`: Path to the created output file

# Details
Creates a compressed tab-delimited file mapping sequence IDs to taxonomy IDs.
Uses blastdbcmd without GI identifiers for better cross-referencing compatibility.
If the output file already exists, returns the path without regenerating.

# Dependencies
Requires BLAST+ tools installed via Bioconda.
"""
function export_blast_db_taxonomy_table(; path_to_db, outfile = path_to_db *
                                                                ".seqid2taxid.txt.gz")
    Mycelia.add_bioconda_env("blast")
    if !isfile(outfile)
        # -long_seqids adds GI identifiers - these are cross-referenceable through other means so I'm dropping
        @time run(pipeline(
            pipeline(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd  -entry all -outfmt "%a %T" -db $(path_to_db)`,
                `gzip`),
            outfile))
    else
        @info "$(outfile) already present"
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Loads a BLAST database taxonomy mapping table from a gzipped file into a DataFrame.

# Arguments
- `compressed_blast_db_taxonomy_table_file::String`: Path to a gzipped file containing BLAST taxonomy mappings

# Returns
- `DataFrame`: A DataFrame with columns `:sequence_id` and `:taxid` containing the sequence-to-taxonomy mappings

# Format
Input file should be a space-delimited text file (gzipped) with two columns:
1. sequence identifier
2. taxonomy identifier (taxid)
"""
function load_blast_db_taxonomy_table(compressed_blast_db_taxonomy_table_file)
    return CSV.read(
        CodecZlib.GzipDecompressorStream(open(compressed_blast_db_taxonomy_table_file)),
        delim = ' ', header = ["sequence_id", "taxid"], DataFrames.DataFrame)
    # data, header = uCSV.read(CodecZlib.GzipDecompressorStream(open(compressed_blast_db_taxonomy_table_file)), delim=' ')
    # header = ["sequence_id", "taxid"]
    # DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a reference table from a BLAST database with taxonomy lineage and sequence lengths.

# Arguments
- `blastdb::String`: BLAST database name or path prefix.
- `blastdbs_dir::String=Mycelia.DEFAULT_BLASTDB_PATH`: Directory where BLAST databases are stored.
- `download_if_missing::Bool=false`: Download the BLAST database if not found locally.
- `reference_fasta::String=""`: Optional path for exported FASTA (defaults to `<blastdb>.fna.gz` in `blastdbs_dir`).
- `taxonomy_map_out::String=""`: Optional path for seqid→taxid mapping output.
- `table_out::String=""`: Optional CSV output path for the reference table.
- `arrow_out::String=""`: Optional Arrow output path for the reference table.
- `force::Bool=false`: Regenerate outputs even if present.

# Returns
- `DataFrames.DataFrame` with sequence_id, taxid, lineage columns, and length.
"""
function prepare_blast_reference_table(;
        blastdb::String,
        blastdbs_dir::String = Mycelia.DEFAULT_BLASTDB_PATH,
        download_if_missing::Bool = false,
        reference_fasta::String = "",
        taxonomy_map_out::String = "",
        table_out::String = "",
        arrow_out::String = "",
        force::Bool = false
)
    db_path = if download_if_missing
        Mycelia.download_blast_db(db = blastdb, dbdir = blastdbs_dir)
    elseif occursin('/', blastdb) || occursin('\\', blastdb)
        blastdb
    else
        joinpath(blastdbs_dir, blastdb)
    end

    if isempty(reference_fasta)
        reference_fasta = joinpath(blastdbs_dir, "$(basename(db_path)).fna.gz")
    end
    Mycelia.blastdb_to_fasta(blastdb = db_path, outfile = reference_fasta, force = force)

    if isempty(taxonomy_map_out)
        taxonomy_map_out = db_path * ".seqid2taxid.txt.gz"
    end
    Mycelia.export_blast_db_taxonomy_table(path_to_db = db_path, outfile = taxonomy_map_out)
    accession_table = Mycelia.load_blast_db_taxonomy_table(taxonomy_map_out)
    lineage = Mycelia.taxids2taxonkit_summarized_lineage_table(unique(accession_table.taxid))
    reference_table = DataFrames.leftjoin(accession_table, lineage, on = :taxid)

    lengths = Dict{String, Int}()
    for record in Mycelia.open_fastx(reference_fasta)
        lengths[String(FASTX.identifier(record))] = length(FASTX.sequence(record))
    end
    reference_table[!, :length] = [get(lengths, id, missing)
                                   for id in reference_table.sequence_id]

    if !isempty(table_out)
        CSV.write(table_out, reference_table)
    end
    if !isempty(arrow_out)
        Arrow.write(arrow_out, reference_table)
    end

    return reference_table
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ncbi_datasets_genome(; kwargs...)

# Download and rehydrate a data package from [NCBI datasets genome tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)

# Specify the download using either

# - [taxon](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/)

# or

# - [accession(https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_accession/)

# # Arguments
# - `annotated::Bool=false`: Limit to annotated genomes.
# - `api_key::String=""`: Specify an NCBI API key.
# - `assembly_level::Array{String}=[]`: Limit to genomes at specific assembly levels (e.g., "chromosome", "complete", "contig", "scaffold"). Default is empty (no specific level).
# - `assembly_source::String="all"`: Limit to 'RefSeq' (GCF_) or 'GenBank' (GCA_) genomes. Default is "all".
# - `assembly_version::String=""`: Limit to 'latest' assembly accession version or include 'all' (latest + previous versions).
# - `chromosomes::Array{String}=[]`: Limit to a specified, comma-delimited list of chromosomes, or 'all' for all chromosomes.
# - `debug::Bool=false`: Emit debugging info.
# - `dehydrated::Bool=false`: Download a dehydrated zip archive including the data report and locations of data files (use the rehydrate command to retrieve data files).
# - `exclude_atypical::Bool=false`: Exclude atypical assemblies.
# - `filename::String=""`: Specify a custom file name for the downloaded data package. Default is "taxon.zip" or "accession.zip" if left blank.
# - `include::Array{String}=["genome"]`: Specify the data files to include (e.g., "genome", "rna", "protein"). Default includes genomic sequence files only.
# - `mag::String="all"`: Limit to metagenome assembled genomes (only) or remove them from the results (exclude). Default is "all".
# - `no_progressbar::Bool=false`: Hide the progress bar.
# - `preview::Bool=false`: Show information about the requested data package without downloading.
# - `reference::Bool=false`: Limit to reference genomes.
# - `released_after::String=""`: Limit to genomes released on or after a specified date (MM/DD/YYYY).
# - `released_before::String=""`: Limit to genomes released on or before a specified date (MM/DD/YYYY).
# - `search::Array{String}=[]`: Limit results to genomes with specified text in the searchable fields (e.g., species, assembly name).

# # Returns
# - The result of the API call.
# """

# function ncbi_datasets_genome(;
#         taxon=missing,
#         accession=missing,
#         annotated=false,
#         api_key="",
#         assembly_level=[],
#         assembly_source="all",
#         assembly_version="",
#         chromosomes=[],
#         debug=false,
#         dehydrated=false,
#         exclude_atypical=false,
#         filename="",
#         outdir="",
#         include=["genome"],
#         mag="all",
#         no_progressbar=false,
#         preview=false,
#         reference=false,
#         released_after="",
#         released_before="",
#         search=[])

#     # Base command
#     command = "$(CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome "

#     if !ismissing(taxon) && !ismissing(accession)
#         @error "can only provide taxon or accession, not both" taxon accession
#     elseif ismissing(taxon) && ismissing(accession)
#         @error "must provide either taxon or accession"
#     elseif !ismissing(taxon) && ismissing(accession)
#         command *= "taxon $(taxon) "
#         if isempty(filename)
#             filename = string(taxon) * ".zip"
#         end
#     elseif !ismissing(accession) && ismissing(taxon)
#         command *= "accession $(accession) "
#         if isempty(filename)
#             filename = string(accession) * ".zip"
#         end
#     end

#     @assert occursin(r"\.zip$", filename)

#     if !isempty(outdir)
#         filename = joinpath(outdir, filename)
#     end

#     annotated && (command *= "--annotated ")
#     !isempty(api_key) && (command *= "--api-key $api_key ")
#     !isempty(assembly_level) && (command *= "--assembly-level $(join(assembly_level, ',')) ")
#     command *= "--assembly-source $assembly_source "
#     !isempty(assembly_version) && (command *= "--assembly-version $assembly_version ")
#     !isempty(chromosomes) && (command *= "--chromosomes $(join(chromosomes, ',')) ")
#     debug && (command *= "--debug ")
#     dehydrated && (command *= "--dehydrated ")
#     exclude_atypical && (command *= "--exclude-atypical ")
#     command *= "--filename $filename "
#     !isempty(include) && (command *= "--include $(join(include, ',')) ")
#     command *= "--mag $mag "
#     no_progressbar && (command *= "--no-progressbar ")
#     preview && (command *= "--preview ")
#     reference && (command *= "--reference ")
#     !isempty(released_after) && (command *= "--released-after $released_after ")
#     !isempty(released_before) && (command *= "--released-before $released_before ")
#     for s in search
#         command *= "--search $s "
#     end

#     # Execute the command
#     println("Executing command: $command")
#     run(`$command`)
#     if dehydrated
#         @info "add code to rehydrate here"
#     end
#     return true
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and sets up MMseqs2 reference databases for sequence searching and analysis.

# Arguments
- `db::String`: Name of database to download (see table below)
- `dbdir::String`: Directory to store the downloaded database (default: "~/workspace/mmseqs")
- `force::Bool`: If true, force re-download even if database exists (default: false)
- `wait::Bool`: If true, wait for download to complete (default: true)

# Returns 
- Path to the downloaded database as a String

# Available Databases

| Database           | Type       | Taxonomy | Description                               |
|-------------------|------------|----------|-------------------------------------------|
| UniRef100         | Aminoacid  | Yes      | UniProt Reference Clusters - 100% identity|
| UniRef90          | Aminoacid  | Yes      | UniProt Reference Clusters - 90% identity |
| UniRef50          | Aminoacid  | Yes      | UniProt Reference Clusters - 50% identity |
| UniProtKB         | Aminoacid  | Yes      | Universal Protein Knowledge Base          |
| NR               | Aminoacid  | Yes      | NCBI Non-redundant proteins              |
| NT               | Nucleotide | No       | NCBI Nucleotide collection               |
| GTDB             | Aminoacid  | Yes      | Genome Taxonomy Database                  |
| PDB              | Aminoacid  | No       | Protein Data Bank structures             |
| Pfam-A.full      | Profile    | No       | Protein family alignments                |
| SILVA            | Nucleotide | Yes      | Ribosomal RNA database                   |

```
  Name                  Type            Taxonomy        Url                                                           
- UniRef100             Aminoacid            yes        https://www.uniprot.org/help/uniref
- UniRef90              Aminoacid            yes        https://www.uniprot.org/help/uniref
- UniRef50              Aminoacid            yes        https://www.uniprot.org/help/uniref
- UniProtKB             Aminoacid            yes        https://www.uniprot.org/help/uniprotkb
- UniProtKB/TrEMBL      Aminoacid            yes        https://www.uniprot.org/help/uniprotkb
- UniProtKB/Swiss-Prot  Aminoacid            yes        https://uniprot.org
- NR                    Aminoacid            yes        https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- NT                    Nucleotide             -        https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- GTDB                  Aminoacid            yes        https://gtdb.ecogenomic.org
- PDB                   Aminoacid              -        https://www.rcsb.org
- PDB70                 Profile                -        https://github.com/soedinglab/hh-suite
- Pfam-A.full           Profile                -        https://pfam.xfam.org
- Pfam-A.seed           Profile                -        https://pfam.xfam.org
- Pfam-B                Profile                -        https://xfam.wordpress.com/2020/06/30/a-new-pfam-b-is-released
- CDD                   Profile                -        https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
- eggNOG                Profile                -        http://eggnog5.embl.de
- VOGDB                 Profile                -        https://vogdb.org
- dbCAN2                Profile                -        http://bcb.unl.edu/dbCAN2
- SILVA                 Nucleotide           yes        https://www.arb-silva.de
- Resfinder             Nucleotide             -        https://cge.cbs.dtu.dk/services/ResFinder
- Kalamari              Nucleotide           yes        https://github.com/lskatz/Kalamari
```
"""
function download_mmseqs_db(; db, dbdir = "$(homedir())/workspace/mmseqs", force = false, wait = true)
    Mycelia.add_bioconda_env("mmseqs2")
    mkpath(dbdir)
    # sanitized_db = replace(db, "/" => "_")
    db_path = joinpath(dbdir, db)
    mkpath(dirname(db_path))
    if !isfile(db_path) || force
        cmd = `$(CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs databases --compressed 1 $(db) --remove-tmp-files 1 $(dbdir)/$(db) $(dbdir)/tmp`
        @time run(cmd, wait = wait)
    else
        @info "db $db @ $(db_path) already exists, set force=true to overwrite"
    end
    return db_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert NCBI sequence IDs to accession numbers.
"""
function ids_to_accessions(ids::Vector{String}, database::String)
    if isempty(ids)
        return String[]
    end

    # Batch convert IDs to accessions
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    id_string = join(ids, ",")

    params = [
        "db" => database,
        "id" => id_string,
        "rettype" => "acc",
        "retmode" => "text"
    ]

    sleep(0.34)  # Rate limiting

    try
        response = HTTP.get(fetch_url, query = params)
        content = String(response.body)

        # Split by newlines and filter out empty lines
        accessions = filter(!isempty, split(strip(content), '\n'))

        return accessions
    catch e
        @error "Failed to convert IDs to accessions: $e"
        return String[]
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Download a sequence from NCBI by accession number with flexible format options.

# Arguments
- `accession::String`: NCBI accession number
- `outdir::String`: Output directory
- `database::String`: NCBI database ("nuccore", "protein", etc.)
- `format::String`: Output format ("fasta", "fasta_cds_na", "fasta_cds_aa", etc.)
- `suffix::String`: File suffix to append to accession for filename
- `compressed::Bool`: Whether to gzip compress the output (default: true)

# Returns
- `String`: Path to the downloaded file
"""
function download_sequence_by_accession(; accession::String,
        outdir::String = pwd(),
        database::String = "nuccore",
        format::String = "fasta",
        suffix::String = ".fasta",
        compressed::Bool = true)

    # Construct output filename
    base_filename = accession * suffix
    if compressed && !endswith(suffix, ".gz")
        outfile = joinpath(outdir, base_filename * ".gz")
        temp_file = joinpath(outdir, base_filename)
    else
        outfile = joinpath(outdir, base_filename)
        temp_file = outfile
    end

    # Skip if file already exists
    if isfile(outfile)
        @info "File already exists: $outfile"
        return outfile
    end

    # Download using NCBI E-utilities
    try
        @info "Downloading $accession from $database (format: $format)"

        # Use efetch to get the sequence
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = [
            "db" => database,
            "id" => accession,
            "rettype" => format,
            "retmode" => "text"
        ]

        sleep(0.34)  # Rate limiting

        response = HTTP.get(fetch_url, query = params)
        content = String(response.body)

        # Write to file
        open(temp_file, "w") do io
            write(io, content)
        end

        # Compress if requested
        if compressed && temp_file != outfile
            run(`gzip $(temp_file)`)
        end

        # Verify file was created
        if isfile(outfile) && filesize(outfile) > 0
            @info "Successfully downloaded: $outfile"
            return outfile
        else
            @error "Download failed or resulted in empty file: $outfile"
            return ""
        end

    catch e
        @error "Failed to download $accession: $e"
        return ""
    end
end
