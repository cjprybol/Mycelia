# """
#     build_exact_match_table(uniref_df::DataFrames.DataFrame,
#                             observed_df::DataFrames.DataFrame;
#                             taxonomy_level::Symbol = :cluster,
#                             compute_bits::Bool = false)

# Constructs a DataFrame mimicking MMseqs2 easy-search output for exact (hash-based) matches.

# Inputs must both contain a `seq_hash` column. `observed_df` must contain at least
# `query` and `qheader`. If `compute_bits=true`, a `sequence` column in `observed_df`
# is preferred; otherwise representative sequence from UniRef is used.

# taxonomy_level:
#   :cluster        -> taxid = common_taxon_id, taxname = common_taxon
#   :representative -> taxid = rep_ncbi_taxonomy, taxname = common_taxon (fallback) or rep_protein_name

# The output columns are:
# query qheader target theader pident fident nident alnlen mismatch gapopen qstart qend qlen tstart tend tlen evalue bits taxid taxname

# Column `theader` is the full FASTA header for the target (target ID plus description),
# matching MMseqs2 style, e.g.:
# UniRef100_P99999 Cytochrome c n=5 Tax=Hominidae TaxID=9604 RepID=CYC_HUMAN
# """
# function build_exact_match_table(uniref_df::DataFrames.DataFrame,
#                                  observed_df::DataFrames.DataFrame;
#                                  taxonomy_level::Symbol = :cluster,
#                                  compute_bits::Bool = false)

#     @assert :seq_hash in propertynames(uniref_df) "uniref_df lacks seq_hash"
#     @assert :seq_hash in propertynames(observed_df) "observed_df lacks seq_hash"
#     @assert :query in propertynames(observed_df) "observed_df lacks query"
#     @assert :qheader in propertynames(observed_df) "observed_df lacks qheader"

#     joined = DataFrames.innerjoin(
#         observed_df,
#         uniref_df,
#         on = :seq_hash,
#         makeunique = true
#     )

#     has_obs_seq = :sequence in propertynames(observed_df)
#     rep_len_col = similar(joined.seq_hash, Int)
#     for (i, _) in enumerate(rep_len_col)
#         len_entry = joined.rep_sequence_length_attr[i]
#         rep_len = if len_entry !== missing
#             Int(len_entry)
#         elseif joined.rep_sequence[i] !== missing
#             length(joined.rep_sequence[i])
#         elseif has_obs_seq && joined.sequence[i] !== missing
#             length(joined.sequence[i])
#         else
#             0
#         end
#         rep_len_col[i] = rep_len
#     end

#     # Build full FASTA header (theader) including target ID at start
#     theader_col = Vector{String}(undef, DataFrames.nrow(joined))
#     for i in eachindex(theader_col)
#         cluster_name = joined.name[i] === missing ? "" : String(joined.name[i])
#         cluster_name = replace(cluster_name, r"^Cluster:\s*" => "")
#         parts = String[]
#         if !isempty(cluster_name)
#             push!(parts, cluster_name)
#         end
#         mc = joined.member_count[i]
#         if mc !== missing
#             push!(parts, "n=$(mc)")
#         end
#         taxname_cluster = joined.common_taxon[i]
#         if taxname_cluster !== missing
#             push!(parts, "Tax=$(taxname_cluster)")
#         end
#         taxid_cluster = joined.common_taxon_id[i]
#         if taxid_cluster !== missing
#             push!(parts, "TaxID=$(taxid_cluster)")
#         end
#         repid = joined.rep_db_id[i]
#         if repid !== missing
#             push!(parts, "RepID=$(repid)")
#         end
#         desc = join(parts, " ")
#         theader_col[i] = isempty(desc) ? String(joined.entry_id[i]) :
#             (String(joined.entry_id[i]) * " " * desc)
#     end

#     taxid_col = Vector{Union{Missing,String}}(undef, DataFrames.nrow(joined))
#     taxname_col = Vector{Union{Missing,String}}(undef, DataFrames.nrow(joined))
#     for i in eachindex(taxid_col)
#         if taxonomy_level == :cluster
#             taxid_col[i] = joined.common_taxon_id[i]
#             taxname_col[i] = joined.common_taxon[i]
#         elseif taxonomy_level == :representative
#             taxid_col[i] = joined.rep_ncbi_taxonomy[i]
#             taxname_col[i] = joined.common_taxon[i] !== missing ? joined.common_taxon[i] : joined.rep_protein_name[i]
#         else
#             error("taxonomy_level must be :cluster or :representative")
#         end
#     end

#     n = DataFrames.nrow(joined)
#     pident = fill(100.0, n)
#     fident = fill(1.0, n)
#     nident = rep_len_col
#     alnlen = rep_len_col
#     mismatch = fill(0, n)
#     gapopen = fill(0, n)
#     qstart = fill(1, n)
#     qend = rep_len_col
#     qlen = rep_len_col
#     tstart = fill(1, n)
#     tend = rep_len_col
#     tlen = rep_len_col
#     evalue = fill(0.0, n)

#     bits = Vector{Union{Missing,Float64}}(undef, n)
#     if compute_bits
#         for i in 1:n
#             seq = if has_obs_seq && joined.sequence[i] !== missing
#                 String(joined.sequence[i])
#             elseif joined.rep_sequence[i] !== missing
#                 String(joined.rep_sequence[i])
#             else
#                 ""
#             end
#             bits[i] = isempty(seq) ? missing : _compute_bits(seq)
#         end
#     else
#         fill!(bits, missing)
#     end

#     result = DataFrames.DataFrame(
#         query = joined.query,
#         qheader = joined.qheader,
#         target = joined.entry_id,
#         theader = theader_col,
#         pident = pident,
#         fident = fident,
#         nident = nident,
#         alnlen = alnlen,
#         mismatch = mismatch,
#         gapopen = gapopen,
#         qstart = qstart,
#         qend = qend,
#         qlen = qlen,
#         tstart = tstart,
#         tend = tend,
#         tlen = tlen,
#         evalue = evalue,
#         bits = bits,
#         taxid = taxid_col,
#         taxname = taxname_col
#     )

#     return result
# end

# function _compute_bits(seq::AbstractString)
#     raw = 0
#     for c in seq
#         raw += get(BLOSUM62_DIAG, c, 0)
#     end
#     return (BLAST_LAMBDA * raw - log(BLAST_K)) / log(2)
# end

# # --------------- Internal per-entry processing ----------------
# # This function remains the same as the previous working version.
# function _process_entry!(
#     entry_node::EzXML.Node,
#     entry_id_col,
#     updated_col,
#     name_col,
#     rep_db_type_col,
#     rep_db_id_col,
#     rep_uniprotkb_accession_col,
#     rep_uniparc_id_col,
#     rep_uniref90_id_col,
#     rep_uniref50_id_col,
#     rep_protein_name_col,
#     rep_ncbi_taxonomy_col,
#     rep_source_organism_col,
#     rep_is_seed_col,
#     rep_length_reported_col,
#     rep_sequence_col::Union{Nothing,Vector{Union{Missing,String}}},
#     rep_sequence_length_attr_col,
#     rep_sequence_checksum_col,
#     member_count_col,
#     common_taxon_col,
#     common_taxon_id_col,
#     entry_property_types_col,
#     entry_property_values_col,
#     rep_property_types_col,
#     rep_property_values_col,
#     member_db_ids_col,
#     member_uniprotkb_accessions_col,
#     member_uniparc_ids_col,
#     store_sequence::Bool
# )
#     # Attributes
#     eid_attr = EzXML.haskey(entry_node, "id") ? entry_node["id"] : ""
#     updated_attr = EzXML.haskey(entry_node, "updated") ? entry_node["updated"] : missing

#     # <name>
#     entry_name_node = EzXML.findfirst("name", entry_node)
#     entry_name = entry_name_node === nothing ? missing : strip(EzXML.nodecontent(entry_name_node))

#     # Entry-level properties
#     e_prop_types = String[]
#     e_prop_values = String[]
#     local_member_count = missing
#     local_common_taxon = missing
#     local_common_taxon_id = missing

#     for child in EzXML.eachelement(entry_node)
#         if EzXML.nodename(child) == "property"
#             t = EzXML.haskey(child, "type") ? child["type"] : nothing
#             v = EzXML.haskey(child, "value") ? child["value"] : nothing
#             if t !== nothing && v !== nothing
#                 push!(e_prop_types, t)
#                 push!(e_prop_values, v)
#                 if t == "member count"
#                     local_member_count = something(tryparse(Int, v), missing)
#                 elseif t == "common taxon"
#                     local_common_taxon = v
#                 elseif t == "common taxon ID"
#                     local_common_taxon_id = v
#                 end
#             end
#         end
#     end

#     # Representative member
#     rep_db_type, rep_db_id = missing, missing
#     rep_props_types, rep_props_values = String[], String[]
#     rep_uniprotkb_acc, rep_uniparc_id_val, rep_uniref90_id_val, rep_uniref50_id_val = missing, missing, missing, missing
#     rep_protein_name_val, rep_ncbi_taxonomy_val, rep_source_organism_val = missing, missing, missing
#     rep_is_seed_val, rep_length_prop = missing, missing
#     rep_sequence_text, rep_seq_len_attr, rep_seq_checksum = missing, missing, missing

#     rep_member_node = EzXML.findfirst("representativeMember", entry_node)
#     if rep_member_node !== nothing
#         db_ref_node = EzXML.findfirst("dbReference", rep_member_node)
#         if db_ref_node !== nothing
#             rep_db_type = EzXML.haskey(db_ref_node, "type") ? db_ref_node["type"] : missing
#             rep_db_id = EzXML.haskey(db_ref_node, "id") ? db_ref_node["id"] : missing

#             for p in EzXML.eachelement(db_ref_node)
#                 if EzXML.nodename(p) == "property"
#                     pt = EzXML.haskey(p, "type") ? p["type"] : nothing
#                     pv = EzXML.haskey(p, "value") ? p["value"] : nothing
#                     if pt !== nothing && pv !== nothing
#                         push!(rep_props_types, pt)
#                         push!(rep_props_values, pv)
#                         if pt == "UniProtKB accession"; rep_uniprotkb_acc = pv
#                         elseif pt == "UniParc ID"; rep_uniparc_id_val = pv
#                         elseif pt == "UniRef90 ID"; rep_uniref90_id_val = pv
#                         elseif pt == "UniRef50 ID"; rep_uniref50_id_val = pv
#                         elseif pt == "protein name"; rep_protein_name_val = pv
#                         elseif pt == "NCBI taxonomy"; rep_ncbi_taxonomy_val = pv
#                         elseif pt == "source organism"; rep_source_organism_val = pv
#                         elseif pt == "isSeed"; rep_is_seed_val = lowercase(pv) == "true"
#                         elseif pt == "length"; rep_length_prop = something(tryparse(Int, pv), missing)
#                         end
#                     end
#                 end
#             end
#         end

#         seq_node = EzXML.findfirst("sequence", rep_member_node)
#         if seq_node !== nothing
#             rep_sequence_text = replace(EzXML.nodecontent(seq_node), r"\s+" => "")
#             len_attr = EzXML.haskey(seq_node, "length") ? seq_node["length"] : nothing
#             rep_seq_len_attr = len_attr === nothing ? missing : something(tryparse(Int, len_attr), missing)
#             cks = EzXML.haskey(seq_node, "checksum") ? seq_node["checksum"] : nothing
#             rep_seq_checksum = cks === nothing ? missing : cks
#         end
#     end

#     # Other Members
#     member_ids_local, member_uniprot_accs, member_uniparc_ids_local = String[], String[], String[]
#     for child in EzXML.eachelement(entry_node)
#         if EzXML.nodename(child) == "member"
#             dbref = EzXML.findfirst("dbReference", child)
#             if dbref !== nothing
#                 mid = EzXML.haskey(dbref, "id") ? dbref["id"] : ""
#                 push!(member_ids_local, mid)
#                 local_up, local_upc = nothing, nothing
#                 for p in EzXML.eachelement(dbref)
#                     if EzXML.nodename(p) == "property"
#                         pt = EzXML.haskey(p, "type") ? p["type"] : nothing
#                         pv = EzXML.haskey(p, "value") ? p["value"] : nothing
#                         if pt == "UniProtKB accession"; local_up = pv
#                         elseif pt == "UniParc ID"; local_upc = pv
#                         end
#                     end
#                 end
#                 if local_up !== nothing; push!(member_uniprot_accs, local_up); end
#                 if local_upc !== nothing; push!(member_uniparc_ids_local, local_upc); end
#             end
#         end
#     end

#     # Normalize sequence length
#     if rep_sequence_text !== missing
#         actual_len = length(rep_sequence_text)
#         if rep_seq_len_attr !== actual_len
#             rep_seq_len_attr = actual_len
#         end
#     end

#     # Push data to columns
#     push!(entry_id_col, eid_attr)
#     push!(updated_col, updated_attr)
#     push!(name_col, entry_name)
#     push!(rep_db_type_col, rep_db_type)
#     push!(rep_db_id_col, rep_db_id)
#     push!(rep_uniprotkb_accession_col, rep_uniprotkb_acc)
#     push!(rep_uniparc_id_col, rep_uniparc_id_val)
#     push!(rep_uniref90_id_col, rep_uniref90_id_val)
#     push!(rep_uniref50_id_col, rep_uniref50_id_val)
#     push!(rep_protein_name_col, rep_protein_name_val)
#     push!(rep_ncbi_taxonomy_col, rep_ncbi_taxonomy_val)
#     push!(rep_source_organism_col, rep_source_organism_val)
#     push!(rep_is_seed_col, rep_is_seed_val)
#     push!(rep_length_reported_col, rep_length_prop)
#     if store_sequence; push!(rep_sequence_col, rep_sequence_text); end
#     push!(rep_sequence_length_attr_col, rep_seq_len_attr)
#     push!(rep_sequence_checksum_col, rep_seq_checksum)
#     push!(member_count_col, local_member_count)
#     push!(common_taxon_col, local_common_taxon)
#     push!(common_taxon_id_col, local_common_taxon_id)
#     push!(entry_property_types_col, e_prop_types)
#     push!(entry_property_values_col, e_prop_values)
#     push!(rep_property_types_col, rep_props_types)
#     push!(rep_property_values_col, rep_props_values)
#     push!(member_db_ids_col, member_ids_local)
#     push!(member_uniprotkb_accessions_col, member_uniprot_accs)
#     push!(member_uniparc_ids_col, member_uniparc_ids_local)

#     return nothing
# end

# Takes forever to run and uses too much memory, not recommended
# """
#     parse_uniref100_to_dataframe(xml_path::String;
#                                  limit::Union{Nothing,Int}=nothing,
#                                  store_sequence::Bool=true,
#                                  include_ids::Union{Nothing, Set{String}}=nothing)

# Parses a UniRef100 XML file into a DataFrame, with options for filtering and memory management.
# """
# function parse_uniref100_to_dataframe(xml_path::String;
#                                       limit::Union{Nothing,Int}=nothing,
#                                       store_sequence::Bool=true,
#                                       include_ids::Union{Nothing, Set{String}}=nothing)

#     # ---------------- Column storage ----------------
#     entry_id = String[]
#     updated = Vector{Union{Missing,String}}()
#     name = Vector{Union{Missing,String}}()
#     rep_db_type = Vector{Union{Missing,String}}()
#     rep_db_id = Vector{Union{Missing,String}}()
#     rep_uniprotkb_accession = Vector{Union{Missing,String}}()
#     rep_uniparc_id = Vector{Union{Missing,String}}()
#     rep_uniref90_id = Vector{Union{Missing,String}}()
#     rep_uniref50_id = Vector{Union{Missing,String}}()
#     rep_protein_name = Vector{Union{Missing,String}}()
#     rep_ncbi_taxonomy = Vector{Union{Missing,String}}()
#     rep_source_organism = Vector{Union{Missing,String}}()
#     rep_is_seed = Vector{Union{Missing,Bool}}()
#     rep_length_reported = Vector{Union{Missing,Int}}()
#     rep_sequence = store_sequence ? Vector{Union{Missing,String}}() : nothing
#     rep_sequence_length_attr = Vector{Union{Missing,Int}}()
#     rep_sequence_checksum = Vector{Union{Missing,String}}()
#     member_count = Vector{Union{Missing,Int}}()
#     common_taxon = Vector{Union{Missing,String}}()
#     common_taxon_id = Vector{Union{Missing,String}}()
#     entry_property_types = Vector{Vector{String}}()
#     entry_property_values = Vector{Vector{String}}()
#     rep_property_types = Vector{Vector{String}}()
#     rep_property_values = Vector{Vector{String}}()
#     member_db_ids = Vector{Vector{String}}()
#     member_uniprotkb_accessions = Vector{Vector{String}}()
#     member_uniparc_ids = Vector{Vector{String}}()

#     # ---------------- Open & detect gzip ----------------
#     raw_io = open(xml_path, "r")
#     stream = raw_io
#     try
#         if endswith(lowercase(xml_path), ".gz")
#             stream = CodecZlib.GzipDecompressorStream(raw_io)
#         end

#         # ---------------- Line-based entry extraction ----------------
#         buf = IOBuffer()
#         in_entry = false
#         records_read = 0
#         records_kept = 0

#         while !eof(stream)
#             line = readline(stream)

#             if !in_entry && occursin("<entry", line)
#                 in_entry = true
#                 truncate(buf, 0)
#             end

#             if in_entry
#                 write(buf, line)
#                 if occursin("</entry>", line)
#                     in_entry = false
#                     records_read += 1
                    
#                     entry_xml_string = String(take!(buf))
#                     doc = EzXML.parsexml(entry_xml_string)
#                     root = EzXML.root(doc)

#                     # --- MODIFICATION: Filter by ID before full processing ---
#                     if include_ids !== nothing
#                         entry_id_val = EzXML.haskey(root, "id") ? root["id"] : ""
#                         if entry_id_val ∉ include_ids
#                             continue # Skip this entry
#                         end
#                     end
                    
#                     _process_entry!(
#                         root,
#                         entry_id, updated, name,
#                         rep_db_type, rep_db_id,
#                         rep_uniprotkb_accession, rep_uniparc_id, rep_uniref90_id, rep_uniref50_id,
#                         rep_protein_name, rep_ncbi_taxonomy, rep_source_organism, rep_is_seed,
#                         rep_length_reported,
#                         rep_sequence, rep_sequence_length_attr, rep_sequence_checksum,
#                         member_count, common_taxon, common_taxon_id,
#                         entry_property_types, entry_property_values,
#                         rep_property_types, rep_property_values,
#                         member_db_ids, member_uniprotkb_accessions, member_uniparc_ids,
#                         store_sequence
#                     )
#                     records_kept += 1

#                     # --- MODIFICATION: Update progress meter ---
#                     if records_read % 1_000_000 == 0
#                         println("Records read from file: ", records_read)
#                     end

#                     if limit !== nothing && records_kept >= limit
#                         break
#                     end
#                 end
#             end
#         end

#         println("Finished parsing. Total records read: ", records_read, ". Total records kept: ", records_kept)
#     finally
#         close(raw_io)
#     end

#     # ---------------- Build DataFrame ----------------
#     df_cols = Dict{Symbol, Any}(
#         :entry_id => entry_id,
#         :updated => updated,
#         :name => name,
#         :rep_db_type => rep_db_type,
#         :rep_db_id => rep_db_id,
#         :rep_uniprotkb_accession => rep_uniprotkb_accession,
#         :rep_uniparc_id => rep_uniparc_id,
#         :rep_uniref90_id => rep_uniref90_id,
#         :rep_uniref50_id => rep_uniref50_id,
#         :rep_protein_name => rep_protein_name,
#         :rep_ncbi_taxonomy => rep_ncbi_taxonomy,
#         :rep_source_organism => rep_source_organism,
#         :rep_is_seed => rep_is_seed,
#         :rep_length_reported => rep_length_reported,
#         :rep_sequence_length_attr => rep_sequence_length_attr,
#         :rep_sequence_checksum => rep_sequence_checksum,
#         :member_count => member_count,
#         :common_taxon => common_taxon,
#         :common_taxon_id => common_taxon_id,
#         :entry_property_types => entry_property_types,
#         :entry_property_values => entry_property_values,
#         :rep_property_types => rep_property_types,
#         :rep_property_values => rep_property_values,
#         :member_db_ids => member_db_ids,
#         :member_uniprotkb_accessions => member_uniprotkb_accessions,
#         :member_uniparc_ids => member_uniparc_ids,
#     )

#     if store_sequence
#         df_cols[:rep_sequence] = rep_sequence
#     end

#     return DataFrames.DataFrame(df_cols)
# end

# # not working
# function run_phage_annotation(;
#     fasta::AbstractString,
#     pharroka_dbdir = joinpath(homedir(), "workspace", "pharroka"),
#     phold_dbdir = joinpath(homedir(), "workspace", "phold"),
#     phynteny_dbdir = joinpath(homedir(), "workspace", "phynteny"),
#     prefix = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
#     pharroka_outdir = replace(fasta, Mycelia.FASTA_REGEX => "_pharroka"),
#     phold_outdir = joinpath(pharroka_outdir, "phold"),
#     phynteny_outdir = joinpath(phold_outdir, "phynteny"),
#     threads = Sys.CPU_THREADS
# )
#     Mycelia.add_bioconda_env("pharokka")
#     Mycelia.add_bioconda_env("phold")
#     Mycelia.add_bioconda_env("phynteny_transformer")
#     if !isdir(pharroka_dbdir) || isempty(readdir(pharroka_dbdir))
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pharroka install_databases.py -o $(pharroka_dbdir)`)
#     end
#     if !isdir(pharroka_outdir) || isempty(readdir(pharroka_outdir))
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pharroka  pharokka.py -i $(fasta) -o $(pharroka_outdir) -d $(pharroka_dbdir) -t $(threads)  -m  -g prodigal-gv --force`)
#     end
#     if !isdir(phold_dbdir) || isempty(readdir(phold_dbdir))
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phold phold install --database $(phold_dbdir)`)
#     end
#     pharroko_gbk = joinpath(pharroka_outdir, prefix * ".gbk")
#     @assert isfile(pharroko_gbk)
#     if !isdir(phold_outdir) || isempty(readdir(phold_outdir))
#         # skip GPU integration at the moment
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phold phold run --input $(pharroko_gbk) --prefix $(prefix) --output $(phold_outdir) --database $(phold_dbdir) --threads $(threads) --cpu --force`)
#     end
#     if !isdir(phynteny_dbdir) || isempty(readddir(phynteny_dbdir))
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phynteny_transformer install_models -o $(phynteny_dbdir)`)
#     end
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phynteny_transformer  phynteny_transformer --help --advanced`)
# end

"""
    run_pgap(;
        fasta,
        organism,
        pgap_dir,
        outdir,
        as_string,
        force,
        threads,
        prefix,
        min_seq_length,
        filter_short,
        ignore_all_errors,
    )

Run the PGAP (Prokaryotic Genome Annotation Pipeline) tool on a FASTA file.

This function:
- Ensures PGAP is present/updated.
- Materializes a temporary FASTA for execution by reading the input with `Mycelia.open_fastx`,
  filtering out sequences shorter than `min_seq_length` (default 200 nt), and writing the
  result with `Mycelia.write_fasta`.
- Treats a run as successful if either a success marker (`.pgap_success`) is present OR
  all expected output files exist for the given `prefix`.
- If an existing `outdir` is complete, results are reused unless `force=true`.
  If incomplete, it is removed automatically so reruns do not require `force`.
- Cleans up partial output directories on failure so retries do not require `force`.
- Optionally appends the PGAP flag `--ignore-all-errors` (when `ignore_all_errors=true`) to
  continue despite QC errors and still obtain a draft annotation.

Arguments:
- `fasta::AbstractString`: Input FASTA (compressed or uncompressed).
- `organism::AbstractString`: Organism name for annotation.
- `pgap_dir`: Directory containing PGAP installation. Defaults to `~/workspace/pgap`.
- `outdir`: Output directory for results. Defaults to input filename with "_pgap" suffix.
- `as_string::Bool=false`: If true, returns the command string instead of executing it (no temp files created).
- `force::Bool=false`: If true, overwrites an existing completed output directory.
- `threads::Int`: CPU threads to use. Defaults to system CPU count.
- `prefix`: Prefix for output files. Defaults to input filename without extension.
- `min_seq_length::Int=200`: Minimum sequence length (nt) retained in the filtered FASTA.
- `filter_short::Bool=true`: If true, sequences shorter than `min_seq_length` are removed prior to PGAP.
- `ignore_all_errors::Bool=false`: If true, adds `--ignore-all-errors` to the PGAP run to ignore non-essential QC errors.

Returns:
- `String`: Path to `outdir` when `as_string=false`, or the bash command string when `as_string=true`.
"""
function run_pgap(;
    fasta::AbstractString,
    organism::AbstractString,
    pgap_dir = joinpath(homedir(), "workspace", "pgap"),
    outdir = replace(fasta, Mycelia.FASTA_REGEX => "_pgap"),
    as_string::Bool = false,
    force::Bool = false,
    threads::Int = Sys.CPU_THREADS,
    prefix = replace(replace(basename(fasta), Mycelia.FASTA_REGEX => ""), r"[^A-Za-z0-9_-]" => "_"),
    min_seq_length::Int = 200,
    filter_short::Bool = true,
    ignore_all_errors::Bool = false,
)
    # Expected outputs from a successful run (accept either prefix.fna or prefixprefix.fna)
    function _expected_complete(outdir_::AbstractString, prefix_::AbstractString)
        singles = String[
            joinpath(outdir_, "$(prefix_)_cds_from_genomic.fna"),
            joinpath(outdir_, "$(prefix_)_translated_cds.faa"),
            joinpath(outdir_, "$(prefix_)_with_genomic_fasta.gff"),
            joinpath(outdir_, "$(prefix_).faa"),
            joinpath(outdir_, "$(prefix_).gbk"),
            joinpath(outdir_, "$(prefix_).gff"),
            joinpath(outdir_, "$(prefix_).sqn"),
            joinpath(outdir_, "$(prefix_)ani-tax-report.txt"),
            joinpath(outdir_, "$(prefix_)ani-tax-report.xml"),
            joinpath(outdir_, "$(prefix_)checkm.txt"),
            joinpath(outdir_, "$(prefix_)cwltool.log"),
            joinpath(outdir_, "$(prefix_)fastaval.xml"),
            joinpath(outdir_, "$(prefix_)uuid.txt"),
        ]
        any_groups = Vector{Vector{String}}([
            [ joinpath(outdir_, "$(prefix_).fna"),
              joinpath(outdir_, "$(prefix_)$(prefix_).fna") ],
        ])
        if !all(isfile, singles)
            return false
        end
        return all(group -> any(isfile, group), any_groups)
    end

    success_marker = joinpath(outdir, ".pgap_success")
    _run_is_successful(outdir_::AbstractString, prefix_::AbstractString) =
        isfile(success_marker) || _expected_complete(outdir_, prefix_)

    # Ensure PGAP availability
    pgap_py_script = joinpath(pgap_dir, "pgap.py")
    if !isfile(pgap_py_script)
        mkpath(dirname(pgap_py_script))
        run(`wget -O $(pgap_py_script) https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py`)
        @assert isfile(pgap_py_script)
    end
    if isempty(filter(x -> occursin(r"input\-", x), readdir(pgap_dir)))
        @time run(setenv(`python3 $(pgap_py_script) --update`, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
    end

    # Handle pre-existing outdir
    if isdir(outdir)
        if force
            rm(outdir, recursive=true)
        else
            if _run_is_successful(outdir, prefix)
                # @info "Detected completed PGAP run in $outdir; reusing outputs"
                return outdir
            else
                @warn "$outdir exists but appears incomplete; removing it to allow a clean retry"
                rm(outdir, recursive=true)
            end
        end
    end
    @assert !isdir(outdir)

    # If only returning the command string, do not create temp files
    if as_string
        base_cmd = "PGAP_INPUT_DIR=$(pgap_dir) python3 $(pgap_py_script) --output $(outdir) --genome $(fasta) --report-usage-true --taxcheck --auto-correct-tax --organism '$(organism)' --cpu $(threads) --prefix $(prefix)"
        if ignore_all_errors
            base_cmd *= " --ignore-all-errors"
        end
        return base_cmd
    end

    # Build filtered temporary FASTA
    filtered_fasta::Union{Nothing,String} = nothing
    try
        reader = Mycelia.open_fastx(fasta)
        records = FASTX.FASTA.Record[]
        for record in reader
            seq = FASTX.sequence(record)
            if (!filter_short) || (length(seq) >= min_seq_length)
                # Make an owned copy (iterator may reuse buffers)
                push!(records, copy(record))
            end
        end
        if isempty(records)
            throw(ArgumentError("All sequences were filtered out (min_seq_length=$(min_seq_length) nt); nothing to annotate."))
        end
        filtered_fasta = Mycelia.write_fasta(; outfile = tempname() * ".fna", records = records, gzip = false)

        # Construct command with optional ignore-all-errors flag
        cmd = `python3 $(pgap_py_script) --output $(outdir) --genome $(filtered_fasta) --report-usage-true --taxcheck --auto-correct-tax --organism $(organism) --cpu $(threads) --prefix $(prefix)`
        if ignore_all_errors
            cmd = `$cmd --ignore-all-errors`
        end

        try
            @time run(setenv(cmd, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
            open(success_marker, "w") do io
                write(io, "ok\n")
            end
            return outdir
        catch e
            if isdir(outdir)
                try
                    rm(outdir, recursive=true)
                catch
                end
            end
            rethrow(e)
        end
    finally
        if !isnothing(filtered_fasta) && isfile(filtered_fasta)
            rm(filtered_fasta)
        end
    end
end

"""
    run_pgap_taxonomy_check(;
        fasta,
        organism,
        pgap_dir,
        outdir,
        as_string,
        prefix,
        no_internet,
        min_seq_length,
        filter_ambiguous,
    )

Run PGAP's taxonomy check (--taxcheck-only) to validate or obtain corrections for organism taxonomy.

This function runs only the taxonomy validation step of PGAP, which is much faster than a full
annotation. It verifies the provided organism name against NCBI taxonomy and can suggest corrections.

The taxonomy check accepts Genus only or Genus species names.

By default, this uses the same output directory as the full PGAP annotation (_pgap suffix), allowing
seamless reuse of existing taxonomy check results from previous annotation runs. If the taxonomy
report already exists, it will be reused. If the output directory exists but lacks the taxonomy
report, PGAP will be run with a temporary output directory and the taxonomy files will be moved
to the final location.

Sequences are filtered before submission to PGAP to remove short contigs and those containing
ambiguous nucleotides, which can cause PGAP errors.

Arguments:
- `fasta::AbstractString`: Input FASTA file (used to calculate ANI for taxonomy validation).
- `organism::AbstractString`: Organism name (Genus or Genus species) to validate.
- `pgap_dir`: Directory containing PGAP installation. Defaults to `~/workspace/pgap`.
- `outdir`: Output directory for results. Defaults to input filename with "_pgap" suffix (same as full annotation).
- `as_string::Bool=false`: If true, returns the command string instead of executing it.
- `prefix`: Prefix for output files. Defaults to input filename without extension with invalid characters replaced.
- `no_internet::Bool=false`: If true, skip all remote HTTP checks to avoid rate limiting.
- `min_seq_length::Int=200`: Minimum sequence length in base pairs. Shorter sequences are filtered out.
- `filter_ambiguous::Bool=true`: If true, filter out sequences containing ambiguous nucleotides (N, etc.).

Returns:
- `NamedTuple`: Contains `outdir`, `txt_report`, and `xml_report` paths.
"""
function run_pgap_taxonomy_check(;
    fasta::AbstractString,
    organism::AbstractString,
    pgap_dir = joinpath(homedir(), "workspace", "pgap"),
    outdir = replace(fasta, Mycelia.FASTA_REGEX => "_pgap"),
    as_string::Bool = false,
    prefix = replace(replace(basename(fasta), Mycelia.FASTA_REGEX => ""), r"[^A-Za-z0-9_-]" => "_"),
    no_internet::Bool = false,
    min_seq_length::Int = 200,
    filter_ambiguous::Bool = false,
)
    # Expected output from taxonomy check
    tax_report = joinpath(outdir, "$(prefix)ani-tax-report.txt")
    tax_report_xml = joinpath(outdir, "$(prefix)ani-tax-report.xml")
    
    # Check if taxonomy report already exists
    if isfile(tax_report)
        return (outdir = outdir, txt_report = tax_report, xml_report = tax_report_xml)
    end
    
    # Ensure PGAP availability
    pgap_py_script = joinpath(pgap_dir, "pgap.py")
    if !isfile(pgap_py_script)
        mkpath(dirname(pgap_py_script))
        Base.run(`wget -O $(pgap_py_script) https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py`)
        @assert isfile(pgap_py_script)
    end
    if !no_internet && isempty(filter(x -> occursin(r"input\-", x), readdir(pgap_dir)))
        @time Base.run(setenv(`python3 $(pgap_py_script) --update`, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
    end

    # If only returning the command string, return immediately
    if as_string
        if !no_internet
            return "PGAP_INPUT_DIR=$(pgap_dir) python3 $(pgap_py_script) --ignore-all-errors --report-usage-true --taxcheck-only --output $(outdir) --genome $(fasta) --organism '$(organism)' --prefix $(prefix)"
        else
            return "PGAP_INPUT_DIR=$(pgap_dir) python3 $(pgap_py_script) --ignore-all-errors --no-internet --report-usage-false --taxcheck-only --output $(outdir) --genome $(fasta) --organism '$(organism)' --prefix $(prefix)"
        end
    end

    # Build temporary uncompressed FASTA for PGAP with filtering
    temp_fasta::Union{Nothing,String} = nothing
    temp_outdir::Union{Nothing,String} = nothing
    
    try
        reader = Mycelia.open_fastx(fasta)
        records = FASTX.FASTA.Record[]
        filtered_count = 0
        
        for record in reader
            seq = FASTX.sequence(record)
            seq_length = length(seq)
            
            # Filter by length
            if seq_length <= min_seq_length
                filtered_count += 1
                continue
            end
            
            # Filter by ambiguous nucleotides using BioSequences
            if filter_ambiguous
                if BioSequences.hasambiguity(BioSequences.LongDNA{4}(seq))
                    filtered_count += 1
                    continue
                end
            end
            
            # Make an owned copy (iterator may reuse buffers)
            push!(records, copy(record))
        end
        
        if isempty(records)
            throw(ArgumentError("All sequences were filtered out (filtered: $(filtered_count), min_length: $(min_seq_length) bp, filter_ambiguous: $(filter_ambiguous)); nothing to analyze."))
        end
        
        if filtered_count > 0
            # @info "Filtered $(filtered_count) sequences; $(length(records)) sequences retained for taxonomy check"
        end
        
        temp_fasta = Mycelia.write_fasta(; outfile = tempname() * ".fna", records = records, gzip = false)

        # PGAP won't write to an existing directory, so we need to handle this
        # If outdir exists, use a temporary directory and move the results
        use_temp_dir = isdir(outdir)
        actual_outdir = if use_temp_dir
            # Create a temporary directory path without actually creating the directory
            temp_base = tempname()
            # Remove it if it exists
            if isdir(temp_base)
                rm(temp_base, recursive=true)
            end
            temp_outdir = temp_base
            temp_base
        else
            outdir
        end
        
        # Construct and run the taxonomy check command
        if !no_internet
            cmd = `python3 $(pgap_py_script) --ignore-all-errors --report-usage-true --taxcheck-only --output $(actual_outdir) --genome $(temp_fasta) --organism $(organism) --prefix $(prefix)`
        else
            cmd = `python3 $(pgap_py_script) --ignore-all-errors --no-internet --report-usage-false --taxcheck-only --output $(actual_outdir) --genome $(temp_fasta) --organism $(organism) --prefix $(prefix)`
        end
        
        temp_tax_report = joinpath(actual_outdir, "$(prefix)ani-tax-report.txt")
        temp_tax_report_xml = joinpath(actual_outdir, "$(prefix)ani-tax-report.xml")
        
        try
            @time Base.run(setenv(cmd, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
            
            # If we used a temporary directory, move the taxonomy files to the final location
            if use_temp_dir
                if !isdir(outdir)
                    mkpath(outdir)
                end
                if isfile(temp_tax_report)
                    mv(temp_tax_report, tax_report, force=true)
                end
                if isfile(temp_tax_report_xml)
                    mv(temp_tax_report_xml, tax_report_xml, force=true)
                end
                # Clean up temporary directory
                if !isnothing(temp_outdir) && isdir(temp_outdir)
                    rm(temp_outdir, recursive=true, force=true)
                end
            end
            
            return (outdir = outdir, txt_report = tax_report, xml_report = tax_report_xml)
        catch e
            # Clean up temporary directory if it exists
            if !isnothing(temp_outdir) && isdir(temp_outdir)
                try
                    rm(temp_outdir, recursive=true, force=true)
                catch
                end
            end
            # Only clean up the specific taxonomy files if we wrote directly to outdir
            if !use_temp_dir
                if isfile(tax_report)
                    try
                        rm(tax_report)
                    catch
                    end
                end
                if isfile(tax_report_xml)
                    try
                        rm(tax_report_xml)
                    catch
                    end
                end
            end
            rethrow(e)
        end
    finally
        if !isnothing(temp_fasta) && isfile(temp_fasta)
            rm(temp_fasta)
        end
    end
end

"""
    PGAPTaxonomyMatch

Represents a single ANI match from a PGAP taxonomy report.

Fields:
- `ani::Float64`: ANI value between the query assembly and this type
- `query_coverage::Float64`: Percentage of query assembly covered
- `subject_coverage::Float64`: Percentage of subject (type) assembly covered
- `new_seq::Int`: Count of bases best assigned to this type assembly
- `cntm_seq::Int`: Portion of new_seq allocated for contamination evaluation
- `assembly_id::Int`: Release ID of the type assembly
- `flags::String`: Type flags (C=contaminant, E=effectively published, T=trusted species)
- `organism::String`: Organism name of this type assembly
- `assembly_accession::String`: Assembly accession
- `assembly_name::String`: Assembly name
"""
struct PGAPTaxonomyMatch
    ani::Float64
    query_coverage::Float64
    subject_coverage::Float64
    new_seq::Int
    cntm_seq::Int
    assembly_id::Int
    flags::String
    organism::String
    assembly_accession::String
    assembly_name::String
end

"""
    PGAPTaxonomyReport

Parsed results from a PGAP taxonomy check report.

Fields:
- `assembly_name::String`: Name of the query assembly
- `submitted_organism::String`: Organism name as submitted
- `submitted_taxid::Union{Int,Nothing}`: Taxonomic ID of submitted organism
- `submitted_rank::Union{String,Nothing}`: Taxonomic rank of submitted organism
- `submitted_lineage::Union{String,Nothing}`: Full lineage of submitted organism
- `predicted_organism::Union{String,Nothing}`: Predicted organism name (if different)
- `predicted_taxid::Union{Int,Nothing}`: Taxonomic ID of predicted organism
- `predicted_rank::Union{String,Nothing}`: Taxonomic rank of predicted organism
- `predicted_lineage::Union{String,Nothing}`: Full lineage of predicted organism
- `best_match_organism::Union{String,Nothing}`: Best match organism (for CONFIRMED genus-level submissions)
- `best_match_taxid::Union{Int,Nothing}`: Best match taxonomic ID
- `best_match_rank::Union{String,Nothing}`: Best match taxonomic rank
- `best_match_lineage::Union{String,Nothing}`: Best match lineage
- `has_type::Bool`: Whether submitted organism has a type assembly
- `status::String`: One of CONFIRMED, MISASSIGNED, INCONCLUSIVE, or CONTAMINATED
- `confidence::String`: Confidence level (HIGH, MEDIUM, LOW)
- `matches::Vector{PGAPTaxonomyMatch}`: All ANI matches sorted by relevance
"""
struct PGAPTaxonomyReport
    assembly_name::String
    submitted_organism::String
    submitted_taxid::Union{Int,Nothing}
    submitted_rank::Union{String,Nothing}
    submitted_lineage::Union{String,Nothing}
    predicted_organism::Union{String,Nothing}
    predicted_taxid::Union{Int,Nothing}
    predicted_rank::Union{String,Nothing}
    predicted_lineage::Union{String,Nothing}
    best_match_organism::Union{String,Nothing}
    best_match_taxid::Union{Int,Nothing}
    best_match_rank::Union{String,Nothing}
    best_match_lineage::Union{String,Nothing}
    has_type::Bool
    status::String
    confidence::String
    matches::Vector{PGAPTaxonomyMatch}
end

"""
    parse_pgap_taxonomy_report(report_file::AbstractString)

Parse a PGAP taxonomy check report file into a structured format.

Arguments:
- `report_file::AbstractString`: Path to the ani-tax-report.txt file

Returns:
- `PGAPTaxonomyReport`: Structured representation of the taxonomy check results
"""
function parse_pgap_taxonomy_report(report_file::AbstractString)
    lines = readlines(report_file)
    
    # Initialize variables
    assembly_name = ""
    submitted_organism = ""
    submitted_taxid = nothing
    submitted_rank = nothing
    submitted_lineage = nothing
    predicted_organism = nothing
    predicted_taxid = nothing
    predicted_rank = nothing
    predicted_lineage = nothing
    best_match_organism = nothing
    best_match_taxid = nothing
    best_match_rank = nothing
    best_match_lineage = nothing
    has_type = false
    status = ""
    confidence = ""
    matches = PGAPTaxonomyMatch[]
    
    in_table = false
    
    for line in lines
        line = strip(line)
        
        # Skip empty lines
        if isempty(line)
            continue
        end
        
        # Parse header information
        if startswith(line, "ANI report for assembly:")
            assembly_name = strip(split(line, ":", limit=2)[2])
        elseif startswith(line, "Submitted organism:")
            # Parse: Organism name (taxid = X, rank = Y, lineage = Z)
            parts = split(line, ":", limit=2)[2]
            if contains(parts, "(taxid")
                org_match = match(r"^(.*?)\s*\(taxid\s*=\s*(\d+),\s*rank\s*=\s*([^,]+),\s*lineage\s*=\s*(.+)\)\s*$", strip(parts))
                if !isnothing(org_match)
                    submitted_organism = strip(org_match.captures[1])
                    submitted_taxid = parse(Int, org_match.captures[2])
                    submitted_rank = strip(org_match.captures[3])
                    submitted_lineage = strip(org_match.captures[4])
                end
            else
                submitted_organism = strip(parts)
            end
        elseif startswith(line, "Predicted organism:")
            parts = split(line, ":", limit=2)[2]
            if contains(parts, "(taxid")
                org_match = match(r"^(.*?)\s*\(taxid\s*=\s*(\d+),\s*rank\s*=\s*([^,]+),\s*lineage\s*=\s*(.+)\)\s*$", strip(parts))
                if !isnothing(org_match)
                    predicted_organism = strip(org_match.captures[1])
                    predicted_taxid = parse(Int, org_match.captures[2])
                    predicted_rank = strip(org_match.captures[3])
                    predicted_lineage = strip(org_match.captures[4])
                end
            end
        elseif startswith(line, "Best match:")
            parts = split(line, ":", limit=2)[2]
            if contains(parts, "(taxid")
                org_match = match(r"^(.*?)\s*\(taxid\s*=\s*(\d+),\s*rank\s*=\s*([^,]+),\s*lineage\s*=\s*(.+)\)\s*$", strip(parts))
                if !isnothing(org_match)
                    best_match_organism = strip(org_match.captures[1])
                    best_match_taxid = parse(Int, org_match.captures[2])
                    best_match_rank = strip(org_match.captures[3])
                    best_match_lineage = strip(org_match.captures[4])
                end
            end
        elseif startswith(line, "Submitted organism has type:")
            has_type = contains(lowercase(line), "yes")
        elseif startswith(line, "Status:")
            status = strip(split(line, ":", limit=2)[2])
        elseif startswith(line, "Confidence:")
            confidence = strip(split(line, ":", limit=2)[2])
        elseif startswith(line, "ANI") && contains(line, "(Coverages)")
            # Start of table header
            in_table = true
            continue
        elseif startswith(line, "---")
            # Table separator line
            continue
        elseif startswith(line, "Table legend:")
            # End of table
            in_table = false
            continue
        elseif in_table
            # Parse match lines
            # Format: ANI (query_cov subject_cov) NewSeq CntmSeq Assembly Flg Organism (accession, name)
            match_regex = r"^\s*(\d+\.\d+)\s+\(\s*(\d+\.\d+)\s+(\d+\.\d+)\)\s+(\d+)\s+(\d+)\s+(\d+)\s+([A-Z]*)\s+(.+?)\s+\(([^,]+),\s*(.+)\)\s*$"
            m = match(match_regex, line)
            if !isnothing(m)
                push!(matches, PGAPTaxonomyMatch(
                    parse(Float64, m.captures[1]),
                    parse(Float64, m.captures[2]),
                    parse(Float64, m.captures[3]),
                    parse(Int, m.captures[4]),
                    parse(Int, m.captures[5]),
                    parse(Int, m.captures[6]),
                    m.captures[7],
                    strip(m.captures[8]),
                    strip(m.captures[9]),
                    strip(m.captures[10])
                ))
            end
        end
    end
    
    return PGAPTaxonomyReport(
        assembly_name,
        submitted_organism,
        submitted_taxid,
        submitted_rank,
        submitted_lineage,
        predicted_organism,
        predicted_taxid,
        predicted_rank,
        predicted_lineage,
        best_match_organism,
        best_match_taxid,
        best_match_rank,
        best_match_lineage,
        has_type,
        status,
        confidence,
        matches
    )
end

function run_bakta(;
        fasta::AbstractString,
        outdir = replace(fasta, Mycelia.FASTA_REGEX => "_bakta"),
        baktadb_dir=joinpath(homedir(), "workspace", "bakta"),
        mode="full",
        threads=Sys.CPU_THREADS,
        as_string = false,
        force = false
    )
    @assert mode in ["full", "light"]
    @assert isfile(fasta)
    Mycelia.add_bioconda_env("bakta")
    bakta_db_path = joinpath(baktadb_dir, mode)
    current_db_path_contents = filter(x -> !occursin(r"^\.", x), readdir(bakta_db_path))
    if !isdir(bakta_db_path) || isempty(current_db_path_contents)
        @info "downloading bakta db"
        mkpath(bakta_db_path)
        @time run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bakta bakta_db download --output $(bakta_db_path) --type $(mode)`)
        @info "bakta db downloaded!"
    end
    bakta_db_path = joinpath(bakta_db_path, "db")
    @assert isdir(bakta_db_path)
    # --prefix ecoli123
    # --locus-tag eco634
    # --prodigal-tf eco.tf
    # --replicons replicon.tsv
    if !isdir(outdir) || isempty(filter(x -> !occursin(r"^\.", x), readdir(outdir))) || force
        if as_string
            bash_command = "$(Mycelia.CONDA_RUNNER) run --live-stream -n bakta bakta --verbose --db $(bakta_db_path) --output $(outdir) --threads $(threads) $(fasta) --force"
            return bash_command
        else
            @info "running bakta"
            @time run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bakta bakta --verbose --db $(bakta_db_path) --output $(outdir) --threads $(threads) $(fasta) --force`)
            return outdir
        end
    else
        @warn "$(outdir) exists and is non-empty - manually remove the directory before running or use `force=true` to overwrite existing outputs"
        return outdir
    end
end

"""
Run VirSorter2 viral sequence identification tool.

VirSorter2 identifies viral sequences in genomic and metagenomic data using
machine learning models and database comparisons.

# Arguments
- `input_fasta`: Path to input FASTA file
- `output_directory`: Output directory path
- `database_path`: Path to VirSorter2 database directory
- `include_groups`: Comma-separated viral groups to include
    - full set = `dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae`
    - tools original default set = `dsDNAphage,ssDNA`
    - Lavidaviridae = A family of small double‑stranded DNA "virophages" that parasitize the replication machinery of certain NCLDVs
    - NCLDV = Nucleocytoplasmic Large DNA Viruses = An informal clade of large double‑stranded DNA viruses that replicate (at least in part) in the cytoplasm of eukaryotic cells.
- `min_score`: Minimum score threshold for viral sequences
- `min_length`: Minimum sequence length threshold
- `threads`: Number of CPU threads to use
- `provirus_off`: Disable provirus detection
- `max_orf_per_seq`: Maximum ORFs per sequence
- `prep_for_dramv`: Prepare output for DRAMv annotation
- `label`: Label for output files
- `forceall`: Force rerun all steps
- `force`: Force rerun even if output files already exist

# Returns
NamedTuple containing paths to all generated output files and directories
"""
function run_virsorter2(;
    input_fasta,
    output_directory = input_fasta * "_virsorter2",
    database_path = mkpath(joinpath(homedir(), "workspace", "virsorter2_db")),
    include_groups = "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae", # full set is dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae
    min_score = 0.5,
    min_length = 1500,
    threads = Sys.CPU_THREADS,
    provirus_off = false,
    max_orf_per_seq = nothing,
    prep_for_dramv = false,
    label = nothing,
    forceall = false,
    force = false
)
    
    # Define output file paths based on VirSorter2 structure
    final_viral_combined = joinpath(output_directory, "final-viral-combined.fa")
    final_viral_score = joinpath(output_directory, "final-viral-score.tsv")
    final_viral_boundary = joinpath(output_directory, "final-viral-boundary.tsv")
    
    # Label-specific files if label was provided
    if label !== nothing
        labeled_viral_combined = joinpath(output_directory, "$(label)-final-viral-combined.fa")
        labeled_viral_score = joinpath(output_directory, "$(label)-final-viral-score.tsv")
        labeled_viral_boundary = joinpath(output_directory, "$(label)-final-viral-boundary.tsv")
    else
        labeled_viral_combined = nothing
        labeled_viral_score = nothing
        labeled_viral_boundary = nothing
    end
    
    # DRAMv output files if prep_for_dramv was enabled
    dramv_dir = joinpath(output_directory, "for-dramv")
    dramv_viral_combined = prep_for_dramv ? joinpath(dramv_dir, "final-viral-combined-for-dramv.fa") : nothing
    dramv_affi_contigs = prep_for_dramv ? joinpath(dramv_dir, "affi-contigs.tab") : nothing
    
    # Check if output files already exist (unless force is true)
    if !force
        # Build list of expected output files to check
        expected_files = [final_viral_combined, final_viral_score, final_viral_boundary]
        
        # Add label-specific files if label was provided
        if label !== nothing
            push!(expected_files, labeled_viral_combined, labeled_viral_score, labeled_viral_boundary)
        end
        
        # Add DRAMv files if prep_for_dramv is enabled
        if prep_for_dramv
            push!(expected_files, dramv_viral_combined, dramv_affi_contigs)
        end
        
        # Filter out nothing values and check if all files exist
        files_to_check = filter(x -> x !== nothing, expected_files)
        all_files_exist = all(Base.Filesystem.isfile, files_to_check)
        
        if all_files_exist
            @warn "All VirSorter2 output files already exist in $(output_directory). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            # Intermediate directories and files that might be useful
            iter_dir = joinpath(output_directory, "iter-0")
            checkpoints_dir = joinpath(output_directory, "checkpoints")
            logs_dir = joinpath(output_directory, "logs")
            
            # Configuration and temporary files
            config_yaml = joinpath(output_directory, "config.yaml")
            
            # Return comprehensive NamedTuple with all output paths (same as below)
            return (
                # Main output directory
                output_directory = output_directory,
                database_path = database_path,
                
                # Primary output files (most commonly used)
                final_viral_combined = final_viral_combined,
                final_viral_score = final_viral_score,
                final_viral_boundary = final_viral_boundary,
                
                # Label-specific outputs (if label was provided)
                labeled_viral_combined = labeled_viral_combined,
                labeled_viral_score = labeled_viral_score,
                labeled_viral_boundary = labeled_viral_boundary,
                
                # DRAMv compatibility outputs (if enabled)
                dramv_directory = prep_for_dramv ? dramv_dir : nothing,
                dramv_viral_combined = dramv_viral_combined,
                dramv_affi_contigs = dramv_affi_contigs,
                
                # Intermediate directories
                iter_directory = iter_dir,
                checkpoints_directory = checkpoints_dir,
                logs_directory = logs_dir,
                
                # Configuration
                config_file = config_yaml,
                
                # Environment information
                conda_environment = "vs2"
            )
        end
    end
    
    # Install VirSorter2 environment with specific version convention
    env_name = "vs2"
    
    # Check if environment exists, if not create it
    env_exists = try
        run(pipeline(`$(Mycelia.CONDA_RUNNER) env list`, `grep -q "^$(env_name) "`))
        true
    catch
        false
    end
    
    if !env_exists
        println("Creating VirSorter2 conda environment...")
        run(`$(Mycelia.CONDA_RUNNER) create -n $(env_name) -c conda-forge -c bioconda virsorter=2 -y`)
        println("Setting up VirSorter2 database (this may take a while)...")
        
        # Remove incomplete database directory if it exists
        if Base.Filesystem.isdir(database_path)
            Base.Filesystem.rm(database_path, recursive=true, force=true)
        end
        
        # Setup database
        setup_cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n $(env_name) virsorter setup -d $(database_path) -j $(threads)`
        
        # Run database setup
        process = run(setup_cmd)
        
        if !success(process)
            error("VirSorter2 database setup failed. You may need to download manually from https://osf.io/v46sc/download")
        end
    end
    
    # Setup database if it doesn't exist or is incomplete
    if !Base.Filesystem.isdir(database_path) || isempty(Base.Filesystem.readdir(database_path))
        @warn "VirSorter2 database not detected"
        @warn "to install, delete any existing directory at the output location, and then run"
        @warn "`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n $(env_name) virsorter setup -d $(database_path) -j $(threads)`"
        # # Remove incomplete database directory if it exists
        # if Base.Filesystem.isdir(database_path)
        #     Base.Filesystem.rm(database_path, recursive=true, force=true)
        # end
        
        # # Setup database
        # setup_cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n $(env_name) virsorter setup -d $(database_path) -j $(threads)`
        
        # # Run database setup
        # process = run(setup_cmd)
        
        # if !success(process)
        #     error("VirSorter2 database setup failed. You may need to download manually from https://osf.io/v46sc/download")
        # end
    end
    
    # Build VirSorter2 command arguments
    cmd_args = ["$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", env_name, "virsorter", "run"]
    
    # Add required arguments
    push!(cmd_args, "-w", output_directory)
    push!(cmd_args, "-i", input_fasta)
    push!(cmd_args, "-j", string(threads))
    
    # Add optional arguments
    push!(cmd_args, "--include-groups", include_groups)
    push!(cmd_args, "--min-score", string(min_score))
    push!(cmd_args, "--min-length", string(min_length))
    
    if provirus_off
        push!(cmd_args, "--provirus-off")
    end
    
    if max_orf_per_seq !== nothing
        push!(cmd_args, "--max-orf-per-seq", string(max_orf_per_seq))
    end
    
    if prep_for_dramv
        push!(cmd_args, "--prep-for-dramv")
    end
    
    if label !== nothing
        push!(cmd_args, "--label", label)
    end
    
    if forceall
        push!(cmd_args, "--forceall")
    end
    
    # Add positional argument (default is 'all')
    push!(cmd_args, "all")
    
    # Run VirSorter2
    println("Running VirSorter2...")
    run(Cmd(cmd_args))
    
    # Intermediate directories and files that might be useful
    iter_dir = joinpath(output_directory, "iter-0")
    checkpoints_dir = joinpath(output_directory, "checkpoints")
    logs_dir = joinpath(output_directory, "logs")
    
    # Configuration and temporary files
    config_yaml = joinpath(output_directory, "config.yaml")
    
    # Return comprehensive NamedTuple with all output paths
    return (
        # Main output directory
        output_directory = output_directory,
        database_path = database_path,
        
        # Primary output files (most commonly used)
        final_viral_combined = final_viral_combined,
        final_viral_score = final_viral_score,
        final_viral_boundary = final_viral_boundary,
        
        # Label-specific outputs (if label was provided)
        labeled_viral_combined = labeled_viral_combined,
        labeled_viral_score = labeled_viral_score,
        labeled_viral_boundary = labeled_viral_boundary,
        
        # DRAMv compatibility outputs (if enabled)
        dramv_directory = prep_for_dramv ? dramv_dir : nothing,
        dramv_viral_combined = dramv_viral_combined,
        dramv_affi_contigs = dramv_affi_contigs,
        
        # Intermediate directories
        iter_directory = iter_dir,
        checkpoints_directory = checkpoints_dir,
        logs_directory = logs_dir,
        
        # Configuration
        config_file = config_yaml,
        
        # Environment information
        conda_environment = env_name
    )
end

"""
Run geNomad mobile genetic element identification tool.

geNomad identifies viruses and plasmids in genomic and metagenomic data
using machine learning and database comparisons.

# Arguments
- `input_fasta`: Path to input FASTA file
- `output_directory`: Output directory path
- `genomad_dbpath`: Path to geNomad database directory
- `threads`: Number of CPU threads to use
- `cleanup`: Remove intermediate files after completion
- `splits`: Number of splits for memory management
- `force`: Force rerun even if output files already exist

# Returns
NamedTuple containing paths to all generated output files and directories
"""
function run_genomad(;
    input_fasta,
    output_directory=input_fasta * "_genomad",
    genomad_dbpath = mkpath(joinpath(homedir(), "workspace", "genomad")),
    threads = Sys.CPU_THREADS,
    cleanup = true,
    splits = nothing,
    force = false
)
    # Get the base name for output files (remove path and extensions)
    input_basename = splitext(basename(input_fasta))[1]
    if endswith(input_basename, ".fna") || endswith(input_basename, ".fa") || endswith(input_basename, ".fasta")
        input_basename = splitext(input_basename)[1]
    end
    
    # Define expected output paths
    summary_dir = joinpath(output_directory, "$(input_basename)_summary")
    
    # Core output files
    virus_summary = joinpath(summary_dir, "$(input_basename)_virus_summary.tsv")
    plasmid_summary = joinpath(summary_dir, "$(input_basename)_plasmid_summary.tsv")
    virus_fasta = joinpath(summary_dir, "$(input_basename)_virus.fna")
    plasmid_fasta = joinpath(summary_dir, "$(input_basename)_plasmid.fna")
    virus_proteins = joinpath(summary_dir, "$(input_basename)_virus_proteins.faa")
    plasmid_proteins = joinpath(summary_dir, "$(input_basename)_plasmid_proteins.faa")
    virus_genes = joinpath(summary_dir, "$(input_basename)_virus_genes.tsv")
    plasmid_genes = joinpath(summary_dir, "$(input_basename)_plasmid_genes.tsv")
    summary_json = joinpath(summary_dir, "$(input_basename)_summary.json")
    
    # Module directories
    marker_classification_dir = joinpath(output_directory, "$(input_basename)_marker_classification")
    nn_classification_dir = joinpath(output_directory, "$(input_basename)_nn_classification")
    find_proviruses_dir = joinpath(output_directory, "$(input_basename)_find_proviruses")
    annotate_dir = joinpath(output_directory, "$(input_basename)_annotate")
    aggregated_classification_dir = joinpath(output_directory, "$(input_basename)_aggregated_classification")
    
    # Log files
    marker_classification_log = joinpath(output_directory, "$(input_basename)_marker_classification.log")
    nn_classification_log = joinpath(output_directory, "$(input_basename)_nn_classification.log")
    find_proviruses_log = joinpath(output_directory, "$(input_basename)_find_proviruses.log")
    annotate_log = joinpath(output_directory, "$(input_basename)_annotate.log")
    aggregated_classification_log = joinpath(output_directory, "$(input_basename)_aggregated_classification.log")
    summary_log = joinpath(output_directory, "$(input_basename)_summary.log")
    
    # Check if output files already exist (unless force is true)
    if !force
        # Build list of expected core output files to check
        expected_files = [
            virus_summary,
            plasmid_summary,
            virus_fasta,
            plasmid_fasta,
            virus_proteins,
            plasmid_proteins,
            virus_genes,
            plasmid_genes,
            summary_json
        ]
        
        # Check if all core files exist
        all_files_exist = all(Base.Filesystem.isfile, expected_files)
        
        if all_files_exist
            @warn "All geNomad output files already exist in $(output_directory). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            # Return NamedTuple with all paths (same as below)
            return (
                # Main output directory
                output_directory = output_directory,
                summary_directory = summary_dir,
                
                # Summary files (most commonly used)
                virus_summary = virus_summary,
                plasmid_summary = plasmid_summary,
                virus_fasta = virus_fasta,
                plasmid_fasta = plasmid_fasta,
                virus_proteins = virus_proteins,
                plasmid_proteins = plasmid_proteins,
                virus_genes = virus_genes,
                plasmid_genes = plasmid_genes,
                summary_json = summary_json,
                
                # Module directories
                marker_classification_dir = marker_classification_dir,
                nn_classification_dir = nn_classification_dir,
                find_proviruses_dir = find_proviruses_dir,
                annotate_dir = annotate_dir,
                aggregated_classification_dir = aggregated_classification_dir,
                
                # Log files
                marker_classification_log = marker_classification_log,
                nn_classification_log = nn_classification_log,
                find_proviruses_log = find_proviruses_log,
                annotate_log = annotate_log,
                aggregated_classification_log = aggregated_classification_log,
                summary_log = summary_log
            )
        end
    end
    
    # Add genomad environment
    Mycelia.add_bioconda_env("genomad")
    
    # Download database if it doesn't exist or is empty
    final_genomad_dbpath = joinpath(genomad_dbpath, "genomad_db")
    if !Base.Filesystem.isdir(final_genomad_dbpath) || isempty(Base.Filesystem.readdir(final_genomad_dbpath))
        run(`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n genomad genomad download-database $(genomad_dbpath)`)
    end
    
    # Build command arguments
    cmd_args = ["$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "genomad", "genomad", "end-to-end"]
    
    if cleanup
        push!(cmd_args, "--cleanup")
    end
    
    push!(cmd_args, "--threads", string(threads))
    
    if splits !== nothing
        push!(cmd_args, "--splits", string(splits))
    end
    
    push!(cmd_args, input_fasta, output_directory, final_genomad_dbpath)
    
    # Run genomad
    run(Cmd(cmd_args))
    
    # Return NamedTuple with all paths
    return (
        # Main output directory
        output_directory = output_directory,
        summary_directory = summary_dir,
        
        # Summary files (most commonly used)
        virus_summary = virus_summary,
        plasmid_summary = plasmid_summary,
        virus_fasta = virus_fasta,
        plasmid_fasta = plasmid_fasta,
        virus_proteins = virus_proteins,
        plasmid_proteins = plasmid_proteins,
        virus_genes = virus_genes,
        plasmid_genes = plasmid_genes,
        summary_json = summary_json,
        
        # Module directories
        marker_classification_dir = marker_classification_dir,
        nn_classification_dir = nn_classification_dir,
        find_proviruses_dir = find_proviruses_dir,
        annotate_dir = annotate_dir,
        aggregated_classification_dir = aggregated_classification_dir,
        
        # Log files
        marker_classification_log = marker_classification_log,
        nn_classification_log = nn_classification_log,
        find_proviruses_log = find_proviruses_log,
        annotate_log = annotate_log,
        aggregated_classification_log = aggregated_classification_log,
        summary_log = summary_log
    )
end

"""
    run_phispy(input_file::String; output_dir::String="", 
           phage_genes::Int=2, color::Bool=false, prefix::String="",
           phmms::String="", threads::Int=1, metrics::Vector{String}=String[],
           expand_slope::Bool=false, window_size::Int=30, 
           min_contig_size::Int=5000, skip_search::Bool=false,
           output_choice::Int=3, training_set::String="", 
           prokka_args::NamedTuple=NamedTuple(), force::Bool=false)

Run PhiSpy to identify prophages in bacterial genomes.

PhiSpy identifies prophage regions in bacterial (and archaeal) genomes using
multiple approaches including gene composition, AT/GC skew, and optional HMM searches.

# Arguments
- `input_file::String`: Path to input file (FASTA or GenBank format)
- `output_dir::String`: Output directory (default: input_file * "_phispy")  
- `phage_genes::Int`: Minimum phage genes required per prophage region (default: 2, set to 0 for mobile elements)
- `color::Bool`: Add color annotations for CDS based on function (default: false)
- `prefix::String`: Prefix for output filenames (default: basename of input)
- `phmms::String`: Path to HMM database for additional phage gene detection
- `threads::Int`: Number of threads for HMM searches (default: 1)
- `metrics::Vector{String}`: Metrics to use for prediction (default: all standard metrics)
- `expand_slope::Bool`: Expand Shannon slope calculations (default: false)
- `window_size::Int`: Window size for calculations (default: 30)
- `min_contig_size::Int`: Minimum contig size to analyze (default: 5000)
- `skip_search::Bool`: Skip HMM search if already done (default: false)
- `output_choice::Int`: Bitmask for output files (default: 3 for coordinates + GenBank)
- `training_set::String`: Path to custom training set
- `prokka_args::NamedTuple`: Additional arguments to pass to Prokka if FASTA input is provided
- `force::Bool`: Force rerun even if output files already exist (default: false)

# Output Choice Codes (add values for multiple outputs)
- 1: prophage_coordinates.tsv
- 2: GenBank format output  
- 4: prophage and bacterial sequences
- 8: prophage_information.tsv
- 16: prophage.tsv
- 32: GFF3 format (prophages only)
- 64: prophage.tbl
- 128: test data used in random forest
- 256: GFF3 format (full genome)
- 512: all output files

# Returns
A NamedTuple with paths to generated output files (contents depend on output_choice):
- `prophage_coordinates`: prophage_coordinates.tsv file path
- `genbank_output`: Updated GenBank file with prophage annotations
- `prophage_sequences`: Prophage and bacterial sequence files
- `prophage_information`: prophage_information.tsv file path
- `prophage_simple`: prophage.tsv file path
- `gff3_prophages`: GFF3 file with prophage regions only
- `prophage_table`: prophage.tbl file path
- `test_data`: Random forest test data file path
- `gff3_genome`: GFF3 file with full genome annotations
- `output_dir`: Path to output directory
- `input_genbank`: Path to GenBank file used (original or generated by Prokka)
"""
function run_phispy(input_file::String; 
                output_dir::String="",
                phage_genes::Int=2,
                color::Bool=false,
                prefix::String="",
                phmms::String="",
                threads::Int=1,
                metrics::Vector{String}=String[],
                expand_slope::Bool=false,
                window_size::Int=30,
                min_contig_size::Int=5000,
                skip_search::Bool=false,
                output_choice::Int=512,
                training_set::String="",
                prokka_args::NamedTuple=NamedTuple(),
                force::Bool=false)
    
    # Validate input file exists
    if !Base.Filesystem.isfile(input_file)
        throw(ArgumentError("Input file does not exist: $input_file"))
    end
    
    # Set default output directory
    if Base.isempty(output_dir)
        output_dir = input_file * "_phispy"
    end
    
    # Set default prefix from input filename if not provided
    if Base.isempty(prefix)
        prefix = Base.Filesystem.splitext(Base.Filesystem.basename(input_file))[1]
    end
    
    # Build expected output file paths based on output_choice
    output_files = Dict{Symbol, String}()
    
    # Map output choice bits to file paths
    output_mapping = [
        (1, :prophage_coordinates, "prophage_coordinates.tsv"),
        (2, :genbank_output, prefix * ".gbk"),
        (4, :prophage_sequences, prefix * "_phage.fasta"), # Multiple files, returning one representative
        (8, :prophage_information, "prophage_information.tsv"),
        (16, :prophage_simple, "prophage.tsv"),
        (32, :gff3_prophages, prefix * "_prophage.gff3"),
        (64, :prophage_table, "prophage.tbl"),
        (128, :test_data, "test_data.tsv"),
        (256, :gff3_genome, prefix * ".gff3")
    ]
    
    # Check which outputs were requested and add paths
    expected_files = String[]
    for (bit, key, filename) in output_mapping
        if (output_choice & bit) != 0
            file_path = Base.Filesystem.joinpath(output_dir, filename)
            output_files[key] = file_path
            Base.push!(expected_files, file_path)
        end
    end
    
    # Check if output files already exist (unless force is true)
    if !force && Base.Filesystem.isdir(output_dir)
        # Check if all expected files exist
        # all_files_exist = Base.all(Base.Filesystem.isfile, expected_files)
        
        # if all_files_exist && !Base.isempty(expected_files)
        if !isempty(readdir(output_dir))
            @warn "PhiSpy output files already exist in $(output_dir). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            # Add standard output files
            output_files[:output_dir] = Base.Filesystem.abspath(output_dir)
            
            # For existing runs, we need to determine the input genbank path
            # Check if there's a Prokka temp directory
            # prokka_temp_dir = output_dir * "_prokka_temp"
            # if Base.Filesystem.isdir(prokka_temp_dir)
            #     prokka_prefix = prefix * "_prokka"
            #     potential_genbank = Base.Filesystem.joinpath(prokka_temp_dir, prokka_prefix * ".gbk")
            #     if Base.Filesystem.isfile(potential_genbank)
            #         output_files[:input_genbank] = Base.Filesystem.abspath(potential_genbank)
            #     else
            #         output_files[:input_genbank] = Base.Filesystem.abspath(input_file)
            #     end
            # else
            #     output_files[:input_genbank] = Base.Filesystem.abspath(input_file)
            # end
            
            return NamedTuple(output_files)
        end
    end
    
    # Determine input file type and prepare GenBank file
    input_genbank = ""
    file_extension = Base.lowercase(Base.Filesystem.splitext(input_file)[2])
    
    if file_extension in [".fasta", ".fa", ".fna", ".fas"]
        # Input is FASTA - need to run Prokka first
        @info "FASTA input detected. Running Prokka for annotation..."
        
        # Ensure bioconda environment is set up
        Mycelia.add_bioconda_env("phispy")
        
        prokka_output_dir = output_dir * "_prokka_temp"
        prokka_prefix = prefix * "_prokka"
        
        # Set up Prokka arguments
        prokka_kwargs = Dict{Symbol, Any}(
            :output_dir => prokka_output_dir,
            :prefix => prokka_prefix
        )
        
        # Add any additional Prokka arguments provided
        for (key, value) in Base.pairs(prokka_args)
            prokka_kwargs[key] = value
        end
        
        # Run Prokka
        prokka_results = Mycelia.run_prokka(input_file; prokka_kwargs...)
        input_genbank = prokka_results.gbk
        
        @info "Prokka completed. Using generated GenBank file: $input_genbank"
        
    elseif file_extension in [".gb", ".gbk", ".genbank", ".gbf"] || 
           (file_extension == ".gz" && Base.any(x -> Base.occursin(x, Base.lowercase(input_file)), [".gb", ".gbk", ".genbank", ".gbf"]))
        # Input is already GenBank format
        input_genbank = input_file
        @info "GenBank input detected: $input_genbank"
    else
        throw(ArgumentError("Unsupported file format. Please provide FASTA (.fasta, .fa, .fna, .fas) or GenBank (.gb, .gbk, .genbank, .gbf) files."))
    end
    
    # Ensure bioconda environment is set up
    Mycelia.add_bioconda_env("phispy")
    
    # Create output directory if it doesn't exist
    if !Base.Filesystem.isdir(output_dir)
        Base.Filesystem.mkdir(output_dir)
    end
    
    # Build PhiSpy command
    cmd_args = String[
        "$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "phispy", "PhiSpy.py",
        input_genbank,
        "-o", output_dir
    ]
    
    # Add optional arguments
    if phage_genes != 2
        Base.push!(cmd_args, "--phage_genes", Base.string(phage_genes))
    end
    
    if color
        Base.push!(cmd_args, "--color")
    end
    
    # if !Base.isempty(prefix)
    #     Base.push!(cmd_args, "--prefix", prefix)
    # end
    
    if !Base.isempty(phmms)
        Base.push!(cmd_args, "--phmms", phmms)
    end
    
    if threads != 1
        Base.push!(cmd_args, "--threads", Base.string(threads))
    end
    
    if !Base.isempty(metrics)
        Base.push!(cmd_args, "--metrics")
        Base.append!(cmd_args, metrics)
    end
    
    if expand_slope
        Base.push!(cmd_args, "--expand_slope")
    end
    
    if window_size != 30
        Base.push!(cmd_args, "--window_size", Base.string(window_size))
    end
    
    if min_contig_size != 5000
        Base.push!(cmd_args, "--min_contig_size", Base.string(min_contig_size))
    end
    
    if skip_search
        Base.push!(cmd_args, "--skip_search")
    end
    
    if output_choice != 3
        Base.push!(cmd_args, "--output_choice", Base.string(output_choice))
    end
    
    if !Base.isempty(training_set)
        Base.push!(cmd_args, "-t", training_set)
    end
    
    # Run PhiSpy
    try
        @info "Running PhiSpy with command: $(Base.join(cmd_args, " "))"
        Base.run(Base.Cmd(cmd_args))
    catch e
        throw(ErrorException("PhiSpy failed to run: $e"))
    end
    
    # Verify output directory exists
    if !Base.Filesystem.isdir(output_dir)
        throw(ErrorException("PhiSpy output directory was not created: $output_dir"))
    end
    
    # Add standard output files
    output_files[:output_dir] = Base.Filesystem.abspath(output_dir)
    # output_files[:input_genbank] = Base.Filesystem.abspath(input_genbank)
    
    # # Check if key output files exist and warn if missing
    # key_files = [
    #     Base.get(output_files, :prophage_coordinates, ""),
    #     Base.get(output_files, :genbank_output, ""),
    #     Base.get(output_files, :prophage_information, "")
    # ]
    
    # existing_files = Base.filter(f -> !Base.isempty(f) && Base.Filesystem.isfile(f), key_files)
    # missing_files = Base.filter(f -> !Base.isempty(f) && !Base.Filesystem.isfile(f), key_files)
    
    # if !Base.isempty(missing_files)
    #     @warn "Some expected output files were not created: $(Base.join(missing_files, ", "))"
    # end
    
    # if Base.isempty(existing_files)
    #     @warn "No standard PhiSpy output files were found. Check PhiSpy logs for errors."
    # else
    #     @info "PhiSpy completed successfully. Generated $(Base.length(existing_files)) output files."
    # end

    @assert isdir(output_files[:output_dir]) && !isempty(readdir(output_files[:output_dir]))
    
    return NamedTuple(output_files)
end

"""
    run_prokka(input_fasta::String; output_dir::String="", prefix::String="", 
           cpus::Int=0, kingdom::String="Bacteria", genus::String="", 
           species::String="", strain::String="", force_overwrite::Bool=false,
           addgenes::Bool=false, compliant::Bool=false, fast::Bool=false,
           evalue::Float64=1e-06, mincontiglen::Int=1, force::Bool=false)

Run Prokka for rapid prokaryotic genome annotation.

Prokka annotates bacterial, archaeal and viral genomes quickly and produces 
standards-compliant output files including GFF3, GenBank, and FASTA formats.

# Arguments
- `input_fasta::String`: Path to input FASTA file containing contigs
- `output_dir::String`: Output directory (default: input_fasta * "_prokka")
- `prefix::String`: Output file prefix (default: basename of input file)
- `cpus::Int`: Number of CPUs to use, 0 for all available (default: 0)
- `kingdom::String`: Annotation mode - "Bacteria", "Archaea", "Viruses", or "Mitochondria" (default: "Bacteria")
- `genus::String`: Genus name for annotation
- `species::String`: Species name for annotation  
- `strain::String`: Strain name for annotation
- `force_overwrite::Bool`: Force overwrite existing output directory (default: false)
- `addgenes::Bool`: Add 'gene' features for each 'CDS' feature (default: false)
- `compliant::Bool`: Force GenBank/ENA/DDJB compliance (default: false)
- `fast::Bool`: Fast mode - skip CDS product searching (default: false)
- `evalue::Float64`: Similarity e-value cut-off (default: 1e-06)
- `mincontiglen::Int`: Minimum contig size (default: 1, NCBI needs 200)
- `force::Bool`: Force rerun even if output files already exist (default: false)

# Returns
A NamedTuple with paths to all generated output files:
- `gff`: Master annotation in GFF3 format
- `gbk`: Standard GenBank file  
- `fna`: Nucleotide FASTA of input contigs
- `faa`: Protein FASTA of translated CDS sequences
- `ffn`: Nucleotide FASTA of all transcripts
- `sqn`: ASN1 Sequin file for GenBank submission
- `fsa`: Nucleotide FASTA for tbl2asn
- `tbl`: Feature table file
- `err`: NCBI discrepancy report
- `log`: Complete run log
- `txt`: Annotation statistics
- `tsv`: Tab-separated feature table
- `output_dir`: Path to output directory
"""
function run_prokka(input_fasta::String; 
                output_dir::String="",
                prefix::String="",
                cpus::Int=0,
                kingdom::String="Bacteria",
                genus::String="",
                species::String="", 
                strain::String="",
                force_overwrite::Bool=false,
                addgenes::Bool=false,
                compliant::Bool=false,
                fast::Bool=false,
                evalue::Float64=1e-06,
                mincontiglen::Int=1,
                force::Bool=false)
    
    # Validate input file exists
    if !Base.Filesystem.isfile(input_fasta)
        Base.throw(ArgumentError("Input FASTA file does not exist: $input_fasta"))
    end
    
    # Set default output directory
    if Base.isempty(output_dir)
        output_dir = input_fasta * "_prokka"
    end
    
    # Set default prefix from input filename if not provided
    if Base.isempty(prefix)
        prefix = Base.Filesystem.splitext(Base.Filesystem.basename(input_fasta))[1]
    end
    
    # Build output file paths
    output_files = (
        gff = Base.Filesystem.joinpath(output_dir, prefix * ".gff"),
        gbk = Base.Filesystem.joinpath(output_dir, prefix * ".gbk"), 
        fna = Base.Filesystem.joinpath(output_dir, prefix * ".fna"),
        faa = Base.Filesystem.joinpath(output_dir, prefix * ".faa"),
        ffn = Base.Filesystem.joinpath(output_dir, prefix * ".ffn"),
        sqn = Base.Filesystem.joinpath(output_dir, prefix * ".sqn"),
        fsa = Base.Filesystem.joinpath(output_dir, prefix * ".fsa"),
        tbl = Base.Filesystem.joinpath(output_dir, prefix * ".tbl"),
        err = Base.Filesystem.joinpath(output_dir, prefix * ".err"),
        log = Base.Filesystem.joinpath(output_dir, prefix * ".log"), 
        txt = Base.Filesystem.joinpath(output_dir, prefix * ".txt"),
        tsv = Base.Filesystem.joinpath(output_dir, prefix * ".tsv"),
        output_dir = Base.Filesystem.abspath(output_dir)
    )
    
    # Check if output files already exist (unless force is true)
    if !force
        # Build list of expected core output files to check
        expected_files = [
            output_files.gff,
            output_files.gbk,
            output_files.fna,
            output_files.faa,
            output_files.ffn,
            output_files.txt,
            output_files.tsv
        ]
        
        # Check if all core files exist
        all_files_exist = Base.all(Base.Filesystem.isfile, expected_files)
        
        if all_files_exist
            @warn "All Prokka output files already exist in $(output_dir). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            return output_files
        end
    end
    
    # Ensure bioconda environment is set up
    Mycelia.add_bioconda_env("prokka")
    
    # Build prokka command
    cmd_args = String[
        "$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "prokka", "prokka",
        "--outdir", output_dir,
        "--prefix", prefix,
        "--kingdom", kingdom
    ]
    
    # Add optional arguments
    if cpus > 0
        Base.push!(cmd_args, "--cpus", Base.string(cpus))
    elseif cpus == 0
        Base.push!(cmd_args, "--cpus", "0")  # Use all available CPUs
    end
    
    if !Base.isempty(genus)
        Base.push!(cmd_args, "--genus", genus)
    end
    
    if !Base.isempty(species) 
        Base.push!(cmd_args, "--species", species)
    end
    
    if !Base.isempty(strain)
        Base.push!(cmd_args, "--strain", strain)
    end
    
    if force_overwrite
        Base.push!(cmd_args, "--force")
    end
    
    if addgenes
        Base.push!(cmd_args, "--addgenes")  
    end
    
    if compliant
        Base.push!(cmd_args, "--compliant")
    end
    
    if fast
        Base.push!(cmd_args, "--fast")
    end
    
    if evalue != 1e-06
        Base.push!(cmd_args, "--evalue", Base.string(evalue))
    end
    
    if mincontiglen != 1
        Base.push!(cmd_args, "--mincontiglen", Base.string(mincontiglen))
    end
    
    # Add input file
    Base.push!(cmd_args, Base.Filesystem.abspath(input_fasta))
    
    # Run prokka
    try
        Base.run(Base.Cmd(cmd_args))
    catch e
        Base.throw(ErrorException("Prokka failed to run: $e"))
    end
    
    # Verify output directory was created
    if !Base.Filesystem.isdir(output_dir)
        Base.throw(ErrorException("Prokka output directory was not created: $output_dir"))
    end
    
    # Verify key output files exist
    key_files = [output_files.gff, output_files.gbk, output_files.fna]
    missing_files = Base.filter(f -> !Base.Filesystem.isfile(f), key_files)
    
    if !Base.isempty(missing_files)
        @warn "Some expected output files were not created: $(Base.join(missing_files, ", "))"
    end
    
    return output_files
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run AMRFinderPlus on FASTA input to identify antimicrobial resistance genes.

# Arguments
- `fasta::String`: Path to input FASTA file (must match Mycelia.FASTA_REGEX pattern)
- `output_dir::String`: Output directory path (default: input filename + "_amrfinderplus")
- `force::Bool`: Force rerun even if output files already exist (default: false)

# Returns
Path to the output directory containing AMRFinderPlus results

# Details
- For nucleotide FASTA files, automatically runs Mycelia.run_pyrodigal to generate protein sequences
- For protein FASTA files, runs AMRFinderPlus directly  
- Validates input file extension against Mycelia.FASTA_REGEX
- Creates output directory if it doesn't exist
- Skips processing if results already exist in output directory unless force=true
- Uses --plus flag for enhanced detection capabilities

# Files Generated
- `<basename>.amrfinderplus.tsv`: AMRFinderPlus results table
- For nucleotide inputs: intermediate pyrodigal outputs in subdirectory
"""
function run_amrfinderplus(;
        fasta::String,
        output_dir::String = fasta * "_amrfinderplus",
        force::Bool = false
    )
    
    # Validate input file extension
    if !Base.occursin(Mycelia.FASTA_REGEX, fasta)
        Base.error("Input file does not match FASTA format: $(fasta)")
    end
    
    if !Base.Filesystem.isfile(fasta)
        Base.error("Input FASTA file not found: $(fasta)")
    end
    
    # Get base filename for outputs
    base_name = Base.replace(Base.Filesystem.basename(fasta), Mycelia.FASTA_REGEX => "")
    amrfinder_output = Base.Filesystem.joinpath(output_dir, "$(base_name).amrfinderplus.tsv")
    
    # Check if output files already exist (unless force is true)
    if !force
        if Base.Filesystem.isfile(amrfinder_output)
            @warn "AMRFinderPlus output file already exists: $(amrfinder_output). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
                return (;output_dir, amrfinder_output)
        end
    end
    
    # Create output directory
    if !Base.Filesystem.isdir(output_dir)
        Base.Filesystem.mkpath(output_dir)
    end
    
    # Determine sequence type by checking file extension first, then by detection
    sequence_type = :unknown
    
    # Check common file extensions first for efficiency
    fasta_lower = Base.lowercase(fasta)
    if Base.endswith(fasta_lower, ".faa") || Base.endswith(fasta_lower, ".faa.gz")
        sequence_type = :protein
    elseif Base.endswith(fasta_lower, ".fna") || Base.endswith(fasta_lower, ".fna.gz")
        sequence_type = :nucleotide
    else
        # For generic extensions (.fa, .fasta, etc.), detect from sequence content
        @info "Generic FASTA extension detected, analyzing sequence content to determine type"
        for (i, record) in Base.enumerate(Mycelia.open_fastx(fasta))
            if i > 3 Base.break end  # Sample first 3 sequences
            seq_ext = Mycelia.detect_sequence_extension(record)
            if seq_ext == ".faa"
                sequence_type = :protein
                Base.break
            elseif seq_ext == ".fna"
                sequence_type = :nucleotide
                Base.break
            end
        end
        
        if sequence_type == :unknown
            Base.error("Could not determine sequence type for: $(fasta)")
        end
    end
    
    # Determine protein FASTA file to use based on sequence type
    protein_fasta = ""
    
    if sequence_type == :protein
        # Input is already protein - use directly
        protein_fasta = fasta
        @info "Using protein FASTA directly: $(fasta)"
    elseif sequence_type == :nucleotide
        # Input is nucleotide - need to run pyrodigal first
        @info "Nucleotide FASTA detected, running pyrodigal to generate protein sequences"
        pyrodigal_dir = Base.Filesystem.joinpath(output_dir, "pyrodigal")
        pyrodigal_results = Mycelia.run_pyrodigal(fasta_file=fasta, out_dir=pyrodigal_dir)
        protein_fasta = pyrodigal_results.faa
    else
        Base.error("Unsupported sequence type: $(sequence_type)")
    end
    
    Base.@assert Base.Filesystem.isfile(protein_fasta) "Protein FASTA file not found: $(protein_fasta)"
    
    # Run AMRFinderPlus
    @info "Running AMRFinderPlus on protein sequences: $(protein_fasta)"
    
    # Ensure AMRFinderPlus is available
    Mycelia.add_bioconda_env("ncbi-amrfinderplus")
    Base.run(`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n ncbi-amrfinderplus amrfinder -u`)
    
    cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n ncbi-amrfinderplus amrfinder
           -p $(protein_fasta)
           --plus
           --output $(amrfinder_output)`
    
    Base.run(cmd)
    
    Base.@assert Base.Filesystem.isfile(amrfinder_output) "AMRFinderPlus output not generated: $(amrfinder_output)"
    @info "AMRFinderPlus completed successfully. Results: $(amrfinder_output)"
    
    return (;output_dir, amrfinder_output)
end

# VIBRANT
function run_vibrant(;input_fasta, output_dir=input_fasta * "_vibrant", threads=Sys.CPU_THREADS)
    Mycelia._install_vibrant()
    # mkpath(output_dir)
    # run(`bash -lc "VIBRANT_run.py -i $input_fasta -folder $output_dir"`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n vibrant VIBRANT_run.py -i $input_fasta -folder $output_dir -t $threads`)
    # phages_circular.fna
    # integrated_prophage_coordinates.tsv
    # figure_PCA.pdf
    # figure_PCA.tsv
    # summary_normalized.tsv
    # for more information
    # https://github.com/AnantharamanLab/VIBRANT?tab=readme-ov-file#output-explanations--
    return output_dir
end


"""
    run_phageboost(input_fasta::AbstractString, output_dir::AbstractString; force_reinstall::Bool=false)

Run PhageBoost on the provided FASTA file, automatically handling conda environment setup.

This function will:
1. Check if the phageboost_env conda environment exists
2. Create and set up the environment if it doesn't exist
3. Validate that PhageBoost is properly installed
4. Run PhageBoost on the input FASTA file
5. Return the output directory path and list of generated files

# Arguments
- `input_fasta::AbstractString`: Path to the input FASTA file
- `output_dir::AbstractString`: Directory where PhageBoost outputs will be saved
- `force_reinstall::Bool=false`: If true, recreate the environment even if it exists

# Returns
- `NamedTuple` with fields:
  - `output_dir::String`: Path to the output directory
  - `files::Vector{String}`: List of files generated in the output directory
"""
function run_phageboost(;input_fasta::AbstractString, output_dir::AbstractString=input_fasta * "_phageboost", force_reinstall::Bool=false)
    # Check if input file exists
    if !Base.isfile(input_fasta)
        throw(ArgumentError("Input FASTA file does not exist: $input_fasta"))
    end
    
    # Setup environment if needed
    _setup_phageboost_environment(force_reinstall)
    
    # # Validate PhageBoost installation
    # _validate_phageboost_installation()
    
    # Ensure output directory exists
    Base.mkpath(output_dir)
    
    # Run PhageBoost
    println("Running PhageBoost on $input_fasta...")
    try
        # Base.run(`bash -lc "conda activate phageboost_env && PhageBoost -f $input_fasta -o $output_dir"`)
        cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n phageboost_env PhageBoost -f $input_fasta -o $output_dir`
        display(cmd)
        run(cmd)
        println("PhageBoost completed successfully")
    catch e
        throw(ErrorException("PhageBoost execution failed: $e"))
    end
    
    # Get list of output files
    output_files = _get_output_files(output_dir)
    
    return (output_dir=output_dir, files=output_files)
end


"""
    parallel_pyrodigal(normalized_fastas::Vector{String})

Runs Mycelia.run_pyrodigal on a list of FASTA files in parallel using Threads.

Args:
    normalized_fastas: A vector of strings, where each string is a path to a FASTA file.

Returns:
    A tuple containing two elements:
    1. successes (Vector{Tuple{String, Any}}): A vector of tuples, where each tuple contains the
       filename and the result returned by a successful Mycelia.run_pyrodigal call.
    2. failures (Vector{Tuple{String, String}}): A vector of tuples, where each tuple contains the
       filename and the error message string for a failed Mycelia.run_pyrodigal call.
"""
function parallel_pyrodigal(normalized_fastas::Vector{String})
    num_files = Base.length(normalized_fastas)
    Base.println("Processing $(num_files) FASTA files using $(Threads.nthreads()) threads...")

    # Create a Progress object for manual updates
    p = ProgressMeter.Progress(num_files, 1, "Running Pyrodigal: ", 50)

    # Use Channels to collect results and failures thread-safely
    # Channel{Tuple{Filename, ResultType}} - adjust ResultType if known
    successes = Base.Channel{Tuple{String, Any}}(num_files)
    failures = Base.Channel{Tuple{String, String}}(num_files)

    # Use Threads.@threads for parallel execution
    Threads.@threads for fasta_file in normalized_fastas
        result = nothing # Initialize result variable in the loop's scope
        try
            # --- Execute the function ---
            # Base.println("Thread $(Threads.threadid()) processing: $(fasta_file)") # Optional: for debugging
            result = Mycelia.run_pyrodigal(fasta_file = fasta_file) # Capture the result

            # --- Store success ---
            Base.put!(successes, (fasta_file, result))

        catch e
            # --- Store failure ---
            err_msg = Base.sprint(Base.showerror, e) # Get the error message as a string
            Base.println(Base.stderr, "ERROR processing $(fasta_file) on thread $(Threads.threadid()): $(err_msg)")
            Base.put!(failures, (fasta_file, err_msg))
        finally
            # --- Always update progress ---
            ProgressMeter.next!(p)
        end
    end

    # Close channels now that all threads are done writing
    Base.close(successes)
    Base.close(failures)

    # Collect results and failures from the channels
    successful_results = Base.collect(successes)
    failed_files = Base.collect(failures)

    # --- Report Summary ---
    Base.println("\n--- Pyrodigal Processing Summary ---")
    num_success = Base.length(successful_results)
    num_failed = Base.length(failed_files)
    Base.println("Successfully processed: $(num_success)")
    Base.println("Failed: $(num_failed)")

    if !Base.isempty(failed_files)
        Base.println("\nFailures:")
        for (file, err) in failed_files
            Base.println("- File: $(file)\n  Error: $(err)")
        end
    end
    Base.println("------------------------------------")

    return (;successful_results, failed_files) # Return both successes and failures
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Perform comprehensive annotation of a FASTA file including gene prediction, protein homology search,
# and terminator prediction.

# # Arguments
# - `fasta::String`: Path to input FASTA file
# - `identifier::String`: Unique identifier for output directory (default: FASTA filename without extension)
# - `basedir::String`: Base directory for output (default: current working directory)
# - `mmseqsdb::String`: Path to MMseqs2 UniRef50 database (default: "~/workspace/mmseqs/UniRef50")
# - `threads::Int`: Number of CPU threads to use (default: all available)

# # Processing Steps
# 1. Creates output directory and copies input FASTA
# 2. Runs Prodigal for gene prediction (nucleotide, amino acid, and GFF output)
# 3. Performs MMseqs2 homology search against UniRef50
# 4. Predicts terminators using TransTerm
# 5. Combines annotations into a unified GFF file
# 6. Generates GenBank format output

# # Returns
# - `String`: Path to the output directory containing all generated files

# # Files Generated
# - `.prodigal.fna`: Predicted genes (nucleotide)
# - `.prodigal.faa`: Predicted proteins
# - `.prodigal.gff`: Prodigal GFF annotations
# - `.gff`: Combined annotations
# - `.gff.genbank`: Final GenBank format
# """
# function annotate_fasta(;
#         fasta,
#         identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
#         basedir = pwd(),        
#         mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
#         threads=Sys.CPU_THREADS
#     )
#     # @show basedir
#     outdir = joinpath(basedir, identifier)
#     @assert outdir != fasta
#     # if !isdir(outdir)
#     #     @show isdir(outdir)
#     mkpath(outdir)
#     f = joinpath(outdir, basename(fasta))
#     # make this an rclone copy for portability
#     !isfile(f) && cp(fasta, f)
#     nucleic_acid_fasta = f * ".prodigal.fna"
#     amino_acid_fasta = f * ".prodigal.faa"
#     gff_file = f * ".prodigal.gff"
#     if !isfile(nucleic_acid_fasta) || !isfile(amino_acid_fasta) || !isfile(gff_file)
#         Mycelia.run_prodigal(fasta_file=f)
#     end
#     mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
#     mmseqs_gff_file = Mycelia.write_gff(gff = Mycelia.update_gff_with_mmseqs(gff_file, mmseqs_outfile), outfile = mmseqs_outfile * ".gff")
#     transterm_gff_file = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
#     joint_gff = Mycelia.write_gff(
#         gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file)), ["#seqid", "start", "end"]),
#         outfile=f * ".gff")
#     annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff, genbank = joint_gff * ".genbank")

#     transterm_gff_file_raw_fasta = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
#     joint_gff_raw_fasta = Mycelia.write_gff(
#         gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file_raw_fasta)), ["#seqid", "start", "end"]),
#         outfile=f * ".transterm_raw.gff")
#     annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff_raw_fasta, genbank = joint_gff_raw_fasta * ".genbank")
#     # else
#     #     @info "$(outdir) already present, skipping..."
#     # end
#     return outdir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform basic annotation of a FASTA file including gene prediction, protein homology search,
and terminator prediction applicable for phage and bacteria.

# Arguments
- `fasta::String`: Path to input FASTA file
- `identifier::String`: Unique identifier for output directory (default: FASTA filename without extension)
- `basedir::String`: Base directory for output (default: current working directory)
- `mmseqsdb::String`: Path to MMseqs2 UniRef50 database (default: `joinpath(homedir(), "workspace/mmseqs/UniRef50")`)
- `threads::Int`: Number of CPU threads to use (default: all available). Note: This argument is not explicitly used by Pyrodigal or MMseqs2 in this version of the function, they might use their own defaults or require modifications to `run_pyrodigal` or `run_mmseqs_easy_search` to respect it.

# Processing Steps
1. Creates output directory and copies input FASTA.
2. Runs Pyrodigal for gene prediction (nucleotide, amino acid, and GFF output).
3. Performs MMseqs2 homology search against UniRef50.
4. Predicts terminators using TransTerm.
5. Combines annotations into a unified GFF file.
6. Generates GenBank format output.

# Returns
- `String`: Path to the output directory containing all generated files.

# Files Generated (within the output directory specified by `identifier`)
- `(basename(fasta)).pyrodigal.fna`: Predicted genes (nucleotide) from Pyrodigal.
- `(basename(fasta)).pyrodigal.faa`: Predicted proteins from Pyrodigal.
- `(basename(fasta)).pyrodigal.gff`: Pyrodigal GFF annotations.
- `(basename(fasta)).gff`: Combined GFF annotations (MMseqs2 and TransTerm).
- `(basename(fasta)).gff.genbank`: Final GenBank format from the first combined GFF.
- `(basename(fasta)).transterm_raw.gff`: Combined GFF (MMseqs2 and a second TransTerm run).
- `(basename(fasta)).transterm_raw.gff.genbank`: Final GenBank format from the second combined GFF.
"""
function annotate_fasta(;
        fasta::String,
        # identifier::String = replace(basename(fasta), Mycelia.FASTA_REGEX => ""), # Assuming Mycelia.FASTA_REGEX is defined
        # basedir::String = pwd(),      
        mmseqsdb::String = joinpath(homedir(), "workspace/mmseqs/UniRef50"),
        threads::Int = Sys.CPU_THREADS,
        outdir::AbstractString = replace(fasta, Mycelia.FASTA_REGEX => "") * "_annotation"
    )
    
    # outdir = joinpath(basedir, identifier)
    # @assert outdir != fasta "Output directory cannot be the same as the input FASTA file path."

    if isdir(outdir) && !isempty(readdir(outdir))
        @warn "$outdir already exists, remove to regenerate"
        return outdir
    end
    mkpath(outdir)
    
    # Path to the FASTA file copied into the output directory
    f_in_outdir = joinpath(outdir, basename(fasta))
    # Copy input FASTA to output directory if it's not already there or needs update
    if !isfile(f_in_outdir) || mtime(fasta) > mtime(f_in_outdir)
        cp(fasta, f_in_outdir, force=true)
    end

    # --- Gene Prediction using Pyrodigal ---
    # Call run_pyrodigal, assuming it's part of the Mycelia module.
    # run_pyrodigal handles its own output file naming and existence checks.
    # We direct its output to be within our main `outdir`.
    pyrodigal_outputs = Mycelia.run_pyrodigal(fasta_file=f_in_outdir, out_dir=outdir)
    
    nucleic_acid_fasta = pyrodigal_outputs.fna
    amino_acid_fasta = pyrodigal_outputs.faa
    gff_file_pyrodigal = pyrodigal_outputs.gff # Renamed to avoid confusion with later gff_file variables
    # --- End of Pyrodigal section ---

    mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
    # Update GFF with MMseqs results, using Pyrodigal's GFF as base
    mmseqs_gff_file = Mycelia.write_gff(gff = Mycelia.update_gff_with_mmseqs(gff_file_pyrodigal, mmseqs_outfile), outfile = mmseqs_outfile * ".gff")
    
    # Predict terminators using TransTerm on the original sequence copy
    if occursin(r"\.gz$", f_in_outdir)
        @warn "transterm doesn't seem to work with gzip compressed fasta files"
    end
    transterm_results = Mycelia.run_transterm(fasta=f_in_outdir) # Assuming run_transterm returns path or object usable by transterm_output_to_gff
    transterm_gff_file = Mycelia.transterm_output_to_gff(transterm_results)
    
    # Combine MMseqs and TransTerm GFFs
    combined_gff_path = joinpath(outdir, basename(f_in_outdir) * ".gff")
    joint_gff_df = DataFrames.sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file)), ["#seqid", "start", "end"])
    Mycelia.write_gff(gff=joint_gff_df, outfile=combined_gff_path)
    
    genbank_path = combined_gff_path * ".genbank"
    annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f_in_outdir, gff=combined_gff_path, genbank=genbank_path)
    
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Annotate amino acid sequences in a FASTA file using MMseqs2 search against UniRef50 database.

# Arguments
- `fasta`: Path to input FASTA file containing amino acid sequences
- `identifier`: Name for the output directory (defaults to FASTA filename without extension)
- `basedir`: Base directory for output (defaults to current directory)
- `mmseqsdb`: Path to MMseqs2 formatted UniRef50 database (defaults to ~/workspace/mmseqs/UniRef50)
- `threads`: Number of CPU threads to use (defaults to system thread count)

# Returns
- Path to the output directory containing MMseqs2 search results

The function creates a new directory named by `identifier` under `basedir`, copies the input FASTA file,
and runs MMseqs2 easy-search against the specified database. If the output directory already exists,
the function skips processing and returns the directory path.
"""
function annotate_aa_fasta(;
        fasta,
        identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
        basedir = pwd(),
        mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
        threads=Sys.CPU_THREADS
    )
    # @show basedir
    outdir = joinpath(basedir, identifier)
    @assert outdir != fasta
    if !isdir(outdir)
        @show isdir(outdir)
        mkpath(outdir)
        f = joinpath(outdir, basename(fasta))
        # make this an rclone copy for portability
        cp(fasta, f, force=true)
        amino_acid_fasta = f

        mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
    else
        @info "$(outdir) already present, skipping..."
    end
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensure the `padloc` environment and database are installed.

Downloads the environment if missing and updates the padloc database.
"""
function setup_padloc()
    padloc_is_already_installed = check_bioconda_env_is_installed("padloc")
    if !padloc_is_already_installed
        Mycelia.add_bioconda_env("padlocbio::padloc")
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --db-update`)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the 'padloc' tool from the 'padlocbio' conda environment on a given FASTA file.

https://doi.org/10.1093/nar/gkab883
https://github.com/padlocbio/padloc

This function:
- Ensures that the 'padloc' environment is available.
- Attempts to update the 'padloc' database (via setup_padloc()).
- Decompresses compressed FASTA inputs (.gz, .bz2, .xz, .zip) into a temporary directory.
- Runs 'padloc' to produce a CSV plus associated prodigal and HMM output files.
- If the padloc command fails, it logs the error, cleans up temporary files, optionally removes the newly created output directory, and returns a tuple with all output paths set to missing plus the captured error object.
- On success, returns a tuple of discovered/expected files (csv may be missing if not produced) and error = nothing.

Parameters:
- fasta_file::AbstractString: Path to input FASTA (possibly compressed).
- outdir::AbstractString: Output directory (default: derived from FASTA filename with '_padloc' suffix).
- threads: Number of CPU threads for padloc.
- cleanup_on_failure::Bool: If true, and the output directory did not exist prior to this call, it will be removed on failure.
- return_error::Bool: If true, errors are captured and returned; if false, the error is rethrown after cleanup.
- verbose::Bool: If false (default), suppress stdout/stderr from the padloc command unless it fails.
"""
function run_padloc(; fasta_file,
                     outdir = replace(fasta_file, Mycelia.FASTA_REGEX => "") * "_padloc",
                     threads = Sys.CPU_THREADS,
                     cleanup_on_failure::Bool = true,
                     return_error::Bool = true,
                     verbose::Bool = false)

    padloc_outfile = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "") * "_padloc.csv")
    pre_existing_outdir = isdir(outdir)

    if isfile(padloc_outfile)
        if verbose
            @info "$(padloc_outfile) already present"
        end
        padloc_faa = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.faa"))
        padloc_gff = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.gff"))
        padloc_domtblout = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => ".domtblout"))
        return (csv = padloc_outfile,
                faa = padloc_faa,
                gff = padloc_gff,
                domtblout = padloc_domtblout,
                error = nothing)
    end

    setup_padloc()

    if !isdir(outdir)
        mkpath(outdir)
    end

    compressed_ext_regex = r"\.(gz|bz2|xz|zip)$"i
    is_compressed = occursin(compressed_ext_regex, fasta_file)

    temp_dir = ""
    temp_fasta = ""
    input_fasta = fasta_file

    if is_compressed
        temp_dir = mktempdir()
        temp_fname = replace(basename(fasta_file), compressed_ext_regex => "")
        temp_fasta = joinpath(temp_dir, temp_fname)
        if verbose
            @info "Detected compressed FASTA; creating temporary uncompressed copy at $(temp_fasta)"
        end

        open(temp_fasta, "w") do io
            lf = lowercase(fasta_file)
            if endswith(lf, ".gz")
                run(pipeline(`gzip -dc $fasta_file`, stdout=io))
            elseif endswith(lf, ".bz2")
                run(pipeline(`bzip2 -dc $fasta_file`, stdout=io))
            elseif endswith(lf, ".xz")
                run(pipeline(`xz -dc $fasta_file`, stdout=io))
            elseif endswith(lf, ".zip")
                run(pipeline(`unzip -p $fasta_file`, stdout=io))
            else
                run(pipeline(`gzip -dc $fasta_file`, stdout=io))
            end
        end
        input_fasta = temp_fasta
    end

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --fna $(input_fasta) --outdir $(outdir) --cpu $(threads)`
    run_error = nothing

    stdout_str = ""
    stderr_str = ""

    try
        if verbose
            run(cmd)
        else
            out_pipe = Pipe()
            err_pipe = Pipe()

            process = run(pipeline(cmd, stdout=out_pipe, stderr=err_pipe), wait=false)
            close(out_pipe.in)
            close(err_pipe.in)

            stdout_task = @async read(out_pipe, String)
            stderr_task = @async read(err_pipe, String)
            
            wait(process) # Throws ProcessFailedException on non-zero exit

            stdout_str = fetch(stdout_task)
            stderr_str = fetch(stderr_task)
        end
    catch err
        run_error = err
        
        # If we were capturing output, fetch it for the log.
        # This check is needed because tasks might not be defined if `verbose=true`.
        if @isdefined(stdout_task)
            stdout_str = fetch(stdout_task)
        end
        if @isdefined(stderr_task)
            stderr_str = fetch(stderr_task)
        end

        @error "padloc command failed" fasta_file=fasta_file outdir=outdir cmd=string(cmd) error=err stdout=stdout_str stderr=stderr_str

        if cleanup_on_failure && !pre_existing_outdir && isdir(outdir)
            try
                rm(outdir; force=true, recursive=true)
                if verbose
                    @info "Removed output directory after failure: $(outdir)"
                end
            catch cleanup_err
                @warn "Failed to remove output directory after failure" outdir=outdir cleanup_error=cleanup_err
            end
        end
        if !return_error
            rethrow()
        end
    finally
        if is_compressed && temp_dir != ""
            try
                rm(temp_dir; force=true, recursive=true)
                if verbose
                    @info "Removed temporary directory $(temp_dir)"
                end
            catch err
                @warn "Failed to remove temporary files at $(temp_dir)" error=err
            end
        end
    end

    if run_error !== nothing
        return (csv = missing,
                faa = missing,
                gff = missing,
                domtblout = missing,
                error = run_error)
    end

    padloc_faa = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.faa"))
    padloc_gff = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.gff"))
    padloc_domtblout = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => ".domtblout"))

    if !isfile(padloc_outfile)
        padloc_outfile = missing
    end

    return (csv = padloc_outfile,
            faa = padloc_faa,
            gff = padloc_gff,
            domtblout = padloc_domtblout,
            error = nothing)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Multi-Locus Sequence Typing (MLST) analysis on a genome assembly.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing the genome assembly

# Returns
- Path to the output file containing MLST results (`<input>.mlst.out`)

# Details
Uses the `mlst` tool from PubMLST to identify sequence types by comparing allelic 
profiles of housekeeping genes against curated MLST schemes.

# Dependencies
- Requires Bioconda and the `mlst` package
- Automatically sets up conda environment if not present
"""
function run_mlst(fasta_file)
    Mycelia.add_bioconda_env("mlst")    
    mlst_outfile = "$(fasta_file).mlst.out"
    @show mlst_outfile
    if !isfile(mlst_outfile)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mlst mlst $(fasta_file)`, mlst_outfile))
    else
        @info "$(mlst_outfile) already present"
    end
    return mlst_outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run ECTyper for serotyping E. coli genome assemblies.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing assembled genome(s)

# Returns
- `String`: Path to output directory containing ECTyper results
"""
function run_ectyper(fasta_file)
    Mycelia.add_bioconda_env("ectyper")
    outdir = fasta_file * "_ectyper"
    if !isdir(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ectyper ectyper -i $(fasta_file) -o $(outdir)`)
    else
        @info "$(outdir) already present"
    end
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run TransTermHP to predict rho-independent transcription terminators in DNA sequences.

# Arguments
- `fasta`: Path to input FASTA file containing DNA sequences. Can be gzip compressed.
- `gff`: Optional path to GFF annotation file. If provided, improves prediction accuracy

# Returns
- `String`: Path to output file containing TransTermHP predictions

# Details
- If the input `fasta` is gzip compressed, a temporary uncompressed copy will be used only when necessary.
- Uses Conda environment 'transtermhp' for execution.
- Automatically generates coordinate file from FASTA or GFF input.
- Removes temporary files after completion.
- Requires Mycelia's Conda setup.
"""
function run_transterm(;fasta, gff="")
    # Determine the base path for output files, removing .gz if present.
    base_path_for_outputs = isempty(gff) ? fasta : gff
    if endswith(lowercase(base_path_for_outputs), ".gz")
        base_path_for_outputs = base_path_for_outputs[1:end-3]
    end

    # Determine the final output file path from the base path.
    # The logic assumes coordinate files from GFFs will also follow a predictable pattern
    # that can be replicated to determine the final output name.
    transterm_calls_file = replace(base_path_for_outputs, r"\.(gff|gff3|fasta|fna|faa|fa)$" => ".transterm.txt")

    if isfile(transterm_calls_file)
        return transterm_calls_file
    end

    fasta_to_use = fasta
    temp_fasta_path = ""
    coordinates = ""

    try
        # Only decompress fasta if it's gzipped AND we are not using a GFF file.
        # generate_transterm_coordinates_from_fasta requires an uncompressed file.
        # The `transterm` binary itself needs the fasta, but it may handle compression,
        # so we only decompress when our Julia functions require it.
        if isempty(gff) && endswith(lowercase(fasta), ".gz")
            temp_fasta_path = write_fasta(records=collect(open_fastx(fasta)), gzip=false)
            fasta_to_use = temp_fasta_path
        end

        if isempty(gff)
            coordinates = generate_transterm_coordinates_from_fasta(fasta_to_use)
        else
            coordinates = generate_transterm_coordinates_from_gff(gff)
        end
        
        Mycelia.add_bioconda_env("transtermhp")

        conda_base = dirname(dirname(Mycelia.CONDA_RUNNER))
        dat_file = joinpath(conda_base, "envs", "transtermhp", "data", "expterm.dat")

        @assert isfile(dat_file) "expterm.dat not found in Conda environment"

        # Note: We pass the original `fasta` path to the command, assuming `transterm`
        # can handle it (even if compressed). We only needed `fasta_to_use` for the
        # Julia-side coordinate generation.
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n transtermhp transterm -p $(dat_file) $(fasta_to_use) $(coordinates)`
        stderr_buffer = IOBuffer()
        
        # The stdout of the process is redirected to `transterm_calls_file`.
        # The stderr is captured in `stderr_buffer`.
        process = run(pipeline(cmd, stdout=transterm_calls_file, stderr=stderr_buffer), wait=false)
        wait(process)

        if !success(process)
            stderr_output = String(take!(stderr_buffer))
            error("TransTermHP failed with exit code $(process.exitcode):\n$stderr_output")
        end
    finally
        # Cleanup temporary files
        if !isempty(coordinates) && isfile(coordinates)
            rm(coordinates)
        end
        if !isempty(temp_fasta_path) && isfile(temp_fasta_path)
            rm(temp_fasta_path)
        end
    end
    
    return transterm_calls_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate minimal coordinate files required for TransTermHP analysis from FASTA sequences.

Creates artificial gene annotations at sequence boundaries to enable TransTermHP to run
without real gene annotations. For each sequence in the FASTA file, generates two
single-base-pair "genes" at positions 1-2 and (L-1)-L, where L is sequence length.

# Arguments
- `fasta`: Path to input FASTA file containing sequences to analyze

# Returns
- Path to generated coordinate file (original path with ".coords" extension)

# Format
Generated coordinate file follows TransTermHP format:
gene_id start stop chromosome

where chromosome matches FASTA sequence identifiers.

See also: [`run_transterm`](@ref)
"""
function generate_transterm_coordinates_from_fasta(fasta)
    # 10. USING TRANSTERM WITHOUT GENOME ANNOTATIONS

    # TransTermHP uses known gene information for only 3 things: (1) tagging the
    # putative terminators as either "inside genes" or "intergenic," (2) choosing the
    # background GC-content percentage to compute the scores, because genes often
    # have different GC content than the intergenic regions, and (3) producing
    # slightly more readable output. Items (1) and (3) are not really necessary, and
    # (2) has no effect if your genes have about the same GC-content as your
    # intergenic regions.

    # Unfortunately, TransTermHP doesn't yet have a simple option to run without an
    # annotation file (either .ptt or .coords), and requires at least 2 genes to be
    # present. The solution is to create fake, small genes that flank each
    # chromosome. To do this, make a fake.coords file that contains only these two
    # lines:

    # 	fakegene1	1 2	chome_id
    # 	fakegene2	L-1 L	chrom_id

    # where L is the length of the input sequence and L-1 is 1 less than the length
    # of the input sequence. "chrom_id" should be the word directly following the ">"
    # in the .fasta file containing your sequence. (If, for example, your .fasta file
    # began with ">seq1", then chrom_id = seq1).

    # This creates a "fake" annotation with two 1-base-long genes flanking the
    # sequence in a tail-to-tail arrangement: --> <--. TransTermHP can then be run
    # with:

    # 	transterm -p expterm.dat sequence.fasta fake.coords

    # If the G/C content of your intergenic regions is about the same as your genes,
    # then this won't have too much of an effect on the scores terminators receive.
    # On the other hand, this use of TransTermHP hasn't been tested much at all, so
    # it's hard to vouch for its accuracy.

    coords_table = DataFrames.DataFrame(
        gene_id = String[],
        start = Int[],
        stop = Int[],
        chromosome = String[]
    )
    for record in Mycelia.open_fastx(fasta)
        row = (
            gene_id = FASTX.identifier(record) * "_start",
            start = 1,
            stop = 2,
            chromosome = FASTX.identifier(record)
        )
        push!(coords_table, row)

        row = (
            gene_id = FASTX.identifier(record) * "_stop",
            start = length(FASTX.sequence(record))-1,
            stop = length(FASTX.sequence(record)),
            chromosome = FASTX.identifier(record)
        )
        push!(coords_table, row)
    end
    transterm_coordinates_file = fasta * ".coords"
    uCSV.write(transterm_coordinates_file, data = collect(DataFrames.eachcol(coords_table)), delim="  ")
    return transterm_coordinates_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GFF file to a coordinates file compatible with TransTermHP format.

# Arguments
- `gff_file::String`: Path to input GFF file

# Processing
- Converts 1-based to 0-based coordinates
- Extracts gene IDs from the attributes field
- Retains columns: gene_id, start, end, seqid

# Returns
- Path to the generated coordinates file (original filename with '.coords' suffix)

# Output Format
Space-delimited file with columns: gene_id, start, end, seqid
"""
function generate_transterm_coordinates_from_gff(gff_file)
    raw_gff = Mycelia.read_gff(gff_file)
    # switch start to be zero index by subtracting one
    raw_gff[!, "start"] = raw_gff[!, "start"] .- 1
    raw_gff[!, "gene_id"] = last.(split.(first.(split.(raw_gff[!, "attributes"], ";")), '='))
    raw_gff = raw_gff[!, ["gene_id", "start", "end", "#seqid"]]
    transterm_coordinates_file = gff_file * ".coords"
    uCSV.write(transterm_coordinates_file, data = collect(DataFrames.eachcol(raw_gff)), delim="  ")
    return transterm_coordinates_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run tRNAscan-SE to identify and annotate transfer RNA genes in the provided sequence file.
Handles gzip-compressed input files by decompressing them to a temporary location.

# Arguments
- `fna_file::String`: Path to the input FASTA nucleotide file. Can be gzip-compressed.
- `outdir::String`: Path to the output directory (default: `replace(fna_file, Mycelia.FASTA_REGEX => "") * "_trnascan"`).
- `model::String`: The tRNA model to use. Defaults to `"G"`.
    - `"G"`: General model (all three domains).
    - `"E"`: Eukaryotic.
    - `"B"`: Bacterial.
    - `"A"`: Archaeal.
    - `"O"`: Other organellar.
    - `"M mammal"` or `"M vert"`: Mitochondrial models.
- `force::Bool`: If true, forces a rerun even if all output files already exist. Defaults to `false`.
- `verbose::Bool`: If true, prints info messages and command output. Defaults to `false`.

# Returns
- `NamedTuple`: A named tuple containing the paths to all generated files and the output directory: `(out, bed, fasta, structure, stats, log, outdir)`.

# Output Files
Creates the following files in `outdir`:
- `*.trnascan.out`: Main output with tRNA predictions.
- `*.trnascan.bed`: BED format coordinates of tRNAs.
- `*.trnascan.fasta`: FASTA sequences of predicted tRNAs.
- `*.trnascan.struct`: Secondary structure predictions.
- `*.trnascan.stats`: Summary statistics.
- `*.trnascan.log`: Program execution log.

# Notes
- Automatically sets up the tRNAscan-SE Conda environment via `Mycelia.add_bioconda_env`.
- Checks for the existence of expected output files to determine if the process should be skipped.
- Cleans up partially generated files if the `tRNAscan-SE` command fails.
"""
function run_trnascan(; fna_file::String, outdir::String=replace(fna_file, Mycelia.FASTA_REGEX => "") * "_trnascan", model::String="G", force::Bool=false, verbose::Bool=false)
    Mycelia.add_bioconda_env("trnascan-se")

    if !ispath(outdir)
        mkpath(outdir)
    end

    ID = replace(basename(fna_file), Mycelia.FASTA_REGEX => "")
    
    output_files = (
        out       = joinpath(outdir, "$(ID).trnascan.out"),
        bed       = joinpath(outdir, "$(ID).trnascan.bed"),
        fasta     = joinpath(outdir, "$(ID).trnascan.fasta"),
        structure = joinpath(outdir, "$(ID).trnascan.struct"),
        stats     = joinpath(outdir, "$(ID).trnascan.stats"),
        log       = joinpath(outdir, "$(ID).trnascan.log"),
    )

    all_files_exist = all(ispath, values(output_files))

    if all_files_exist && !force
        if verbose
            @info "All output files already exist in $(outdir). Skipping tRNAscan-SE run. Use `force=true` to rerun."
        end
        return (; output_files..., outdir=outdir)
    end
    
    input_to_process = fna_file
    temp_fna_path = ""

    try
        # Handle compressed input file
        if endswith(fna_file, ".gz")
            if verbose
                @info "Input file is gzipped. Decompressing to a temporary file."
            end
            temp_fna_path, temp_io = mktemp()
            close(temp_io) # Immediately close the stream from mktemp

            open(CodecZlib.GzipDecompressorStream, fna_file, "r") do f_in
                open(temp_fna_path, "w") do f_out
                    write(f_out, read(f_in))
                end
            end
            input_to_process = temp_fna_path
        end

        model_parts = split(model)
        model_flag = "-" * model_parts[1]
        
        cmd_parts = [
            "$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "trnascan-se", "tRNAscan-SE",
            model_flag
        ]
        if length(model_parts) > 1
            push!(cmd_parts, model_parts[2])
        end

        append!(cmd_parts, [
            "--output", output_files.out,
            "--bed", output_files.bed,
            "--fasta", output_files.fasta,
            "--struct", output_files.structure,
            "--stats", output_files.stats,
            "--log", output_files.log,
            input_to_process
        ])

        trnascan_cmd = Cmd(cmd_parts)
        if verbose
            @info "Running tRNAscan-SE command:" trnascan_cmd
        end

        try
            if verbose
                run(trnascan_cmd)
            else
                stderr_buffer = IOBuffer()
                # Redirect stdout to devnull and capture stderr
                process = run(pipeline(trnascan_cmd, stdout=devnull, stderr=stderr_buffer), wait=false)
                wait(process)
                
                if !success(process)
                    stderr_output = String(take!(stderr_buffer))
                    # Throw a more informative error including stderr
                    error("tRNAscan-SE failed with exit code $(process.exitcode). Stderr:\n$stderr_output")
                end
            end
        catch e
            @error "tRNAscan-SE failed. Cleaning up partial results."
            for file_path in values(output_files)
                if ispath(file_path)
                    rm(file_path)
                end
            end
            rethrow(e)
        end

    finally
        # Cleanup the temporary file if it was created
        if !isempty(temp_fna_path) && ispath(temp_fna_path)
            if verbose
                @info "Cleaning up temporary file: $(temp_fna_path)"
            end
            rm(temp_fna_path)
        end
    end

    return (; output_files..., outdir=outdir)
end

# function run_counterselection_spacer_detection(strain, out_dir, normalized_fasta_file)
#     counter_selection_dir = "$(out_dir)/counter-selection"
#     if !isdir(counter_selection_dir)
#         mkdir(counter_selection_dir)
#     end

#     if isempty(readdir(counter_selection_dir))
#         regex = BioSequences.biore"TTT[CG][ACGT]{25}"dna
#         k = 29
#         KMER_TYPE = BioSequences.BigDNAMer{k}

#         spacer_table = DataFrames.DataFrame(
#             ID = [],
#             contig = [],
#             PAM_and_spacer = [],
#             spacer = [],
#             strand = [],
#             start = [],
#             stop = [],
# #             free_energy = [],
# #             visualization_url = []
#         )

#         ProgressMeter.@showprogress for record in collect(FASTX.FASTA.Reader(open(normalized_fasta_file)))
#             for (i, kmer, reverse_complement_kmer) in BioSequences.each(KMER_TYPE, FASTX.FASTA.sequence(record))
#                 strand = missing
#                 if occursin(regex, kmer)
#                     strand = "+"
#                 elseif occursin(regex, reverse_complement_kmer)
#                     strand = "-"
#                     kmer = reverse_complement_kmer
#                 end
#                 if !ismissing(strand)
#                     spacer = BioSequences.DNAMer(kmer[i] for i in 5:length(kmer)) 
# #                     RNAfold_output = read(pipeline(`echo "$(string(spacer))"`, `RNAfold --noLP`), String)

# #                     rna_sequence, structure, free_energy = match(r"([ACGU]{25})\n([.()]{25})\s\(\s*(.*?)\)", RNAfold_output).captures
# #                     url = "http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence=$(rna_sequence)&structure=$(structure)"

#                     kmer_range = i:i+k-1

#                     t = DataFrames.DataFrame(
#                         ID = ID,
#                         contig = FASTX.FASTA.identifier(record),
#                         PAM_and_spacer = kmer,
#                         spacer = spacer,
#                         strand = strand,
#                         start = i,
#                         stop = i+k-1,
# #                         free_energy = free_energy,
# #                         visualization_url = url
#                     )
#                     spacer_table = vcat(spacer_table, t)
#                 end
#             end
#         end


#         uCSV.write(
#             "$(counter_selection_dir)/$(ID)-cpf1-spacers.tsv",
#             delim='\t',
#             data = collect(DataFrames.eachcol(spacer_table)),
#             header = DataFrames.names(spacer_table)
#         )
#         if isfile("$(out_dir)/rna.ps")
#             rm("$(out_dir)/rna.ps")
#         end
#     end
#     return counter_selection_dir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Prodigal gene prediction software on input FASTA file to identify protein-coding genes
in metagenomes or single genomes.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing genomic sequences
- `out_dir::String=dirname(fasta_file)`: Directory for output files. Defaults to input file's directory

# Returns
Named tuple containing paths to all output files:
- `fasta_file`: Input FASTA file path
- `out_dir`: Output directory path  
- `gff`: Path to GFF format gene predictions
- `gene_scores`: Path to all potential genes and their scores
- `fna`: Path to nucleotide sequences of predicted genes
- `faa`: Path to protein translations of predicted genes
- `std_out`: Path to captured stdout
- `std_err`: Path to captured stderr
"""
function run_prodigal(;fasta_file, out_dir=dirname(fasta_file))
    
    # if isempty(out_dir)
    #     prodigal_dir = mkpath("$(fasta_file)_prodigal")
    # else
    #     prodigal_dir = mkpath(out_dir)
    # end

    # $ prodigal
    # -------------------------------------
    # PRODIGAL v2.6.3 [February, 2016]         
    # Univ of Tenn / Oak Ridge National Lab
    # Doug Hyatt, Loren Hauser, et al.     
    # -------------------------------------

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    gff = "$(out_dir)/$(basename(fasta_file)).prodigal.gff"
    faa = "$(out_dir)/$(basename(fasta_file)).prodigal.faa"
    fna = "$(out_dir)/$(basename(fasta_file)).prodigal.fna"
    gene_scores = "$(out_dir)/$(basename(fasta_file)).prodigal.all_potential_gene_scores.txt"
    std_out = "$(out_dir)/$(basename(fasta_file)).prodigal.out"
    std_err = "$(out_dir)/$(basename(fasta_file)).prodigal.err"
    
    # I usually delete the rest, so don't reprocess if outputs of interest are present
    if (!isfile(gff) && !isfile(faa))
        add_bioconda_env("prodigal")
        cmd = 
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n prodigal prodigal
        -f gff
        -m
        -p meta
        -o $(gff)
        -i $(fasta_file)
        -a $(faa)
        -d $(fna)
        -s $(gene_scores)
        `
        p = pipeline(cmd, stdout=std_out, stderr=std_err)
        run(p)
    end
    return (;fasta_file, out_dir, gff, gene_scores, fna, faa, std_out, std_err)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Pyrodigal gene prediction on a FASTA file using the meta procedure optimized for metagenomic sequences.

Pyrodigal is a reimplementation of the Prodigal gene finder, which identifies protein-coding sequences in bacterial and archaeal genomes.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing genomic sequences
- `out_dir::String`: Output directory path (default: input filename + "_pyrodigal")

# Returns
Named tuple containing:
- `fasta_file`: Input FASTA file path
- `out_dir`: Output directory path
- `gff`: Path to GFF output file with gene predictions
- `faa`: Path to FASTA file with predicted protein sequences 
- `fna`: Path to FASTA file with nucleotide sequences

# Notes
- Uses metagenomic mode (`-p meta`) optimized for mixed communities
- Masks runs of N nucleotides (`-m` flag)
- Minimum gene length set to 33bp
- Maximum overlap between genes set to 31bp
- Requires Pyrodigal to be available in a Conda environment
- Skips processing if output files already exist
"""
function run_pyrodigal(;fasta_file, out_dir=fasta_file * "_pyrodigal")
    # https://pyrodigal.readthedocs.io/en/stable/guide/cli.html#command-line-interface

    # -a trans_file         Write protein translations to the selected file.
    # -c                    Closed ends. Do not allow genes to run off edges.
    # -d nuc_file           Write nucleotide sequences of genes to the selected file.
    # -f output_type        Select output format.
    # -g tr_table           Specify a translation table to use.
    # -i input_file         Specify FASTA input file.
    # -m                    Treat runs of N as masked sequence and don't build genes across them.
    # -n                    Bypass Shine-Dalgarno trainer and force a full motif scan.
    # -o output_file        Specify output file.
    # -p mode               Select procedure.
    # -s start_file         Write all potential genes (with scores) to the selected file.
    # -t training_file      Write a training file (if none exists); otherwise, read and use the specified training file.
    # -j jobs, --jobs jobs           The number of threads to use if input contains multiple sequences.
    # --min-gene MIN_GENE            The minimum gene length.
    # --min-edge-gene MIN_EDGE_GENE  The minimum edge gene length.
    # --max-overlap MAX_OVERLAP      The maximum number of nucleotides that can overlap between two genes on the same strand.
    #                             This must be lower or equal to the minimum gene length.
    # --no-stop-codon                Disables translation of stop codons into star characters (*) for complete genes.
    # --pool {thread,process}        The sort of pool to use to process genomes in parallel. Processes may be faster than
    #                             threads on some machines, refer to documentation. (default: thread)

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    gff = "$(out_dir)/$(basename(fasta_file)).pyrodigal.gff"
    faa = "$(out_dir)/$(basename(fasta_file)).pyrodigal.faa"
    fna = "$(out_dir)/$(basename(fasta_file)).pyrodigal.fna"
    gene_scores = "$(out_dir)/$(basename(fasta_file)).pyrodigal.gene_scores.txt"
    std_out = "$(out_dir)/$(basename(fasta_file)).pyrodigal.out"
    std_err = "$(out_dir)/$(basename(fasta_file)).pyrodigal.err"
    mkpath(out_dir)
    
    # I usually delete the rest, so don't reprocess if outputs of interest are present
    # `max_overlap` must be lower than `min_gene`
    if (!isfile(gff) || !isfile(faa))
        Mycelia.add_bioconda_env("pyrodigal")
        cmd = 
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n pyrodigal pyrodigal
        -f gff
        -m
        -p meta
        -o $(gff)
        -i $(fasta_file)
        -a $(faa)
        -d $(fna)
        -s $(gene_scores)
        --min-gene 33
        --max-overlap 31
        `
        p = pipeline(cmd, stdout=std_out, stderr=std_err)
        run(p)
        if fasta_file != "$(out_dir)/$(basename(fasta_file))"
            cp(fasta_file, "$(out_dir)/$(basename(fasta_file))", force=true)
        end
        if isfile(std_out) && (filesize(std_out) == 0)
            rm(std_out)
        end
        if isfile(std_err) && (filesize(std_err) == 0)
            rm(std_err)
        end
    end
    # return (;fasta_file, out_dir, gff, gene_scores, fna, faa, std_out, std_err)
    return (;fasta_file, out_dir, gff, faa, fna)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse TransTerm terminator prediction output into a structured DataFrame.

Takes a TransTerm output file path and returns a DataFrame containing parsed terminator predictions.
Each row represents one predicted terminator with the following columns:

- `chromosome`: Identifier of the sequence being analyzed
- `term_id`: Unique terminator identifier (e.g. "TERM 19")
- `start`: Start position of the terminator
- `stop`: End position of the terminator
- `strand`: Strand orientation ("+" or "-")
- `location`: Context type, where:
    * G/g = in gene interior (≥50bp from ends)
    * F/f = between two +strand genes
    * R/r = between two -strand genes
    * T = between ends of +strand and -strand genes
    * H = between starts of +strand and -strand genes
    * N = none of the above
    Lowercase indicates opposite strand from region
- `confidence`: Overall confidence score (0-100)
- `hairpin_score`: Hairpin structure score
- `tail_score`: Tail sequence score  
- `notes`: Additional annotations (e.g. "bidir")

# Arguments
- `transterm_output::AbstractString`: Path to TransTerm output file

# Returns
- `DataFrame`: Parsed terminator predictions with columns as described above

See TransTerm HP documentation for details on scoring and location codes.
"""
function parse_transterm_output(transterm_output)
    
   #     3. FORMAT OF THE TRANSTERM OUTPUT

#     The organism's genes are listed sorted by their end coordinate and terminators
#     are output between them. A terminator entry looks like this:

#         TERM 19  15310 - 15327  -      F     99      -12.7 -4.0 |bidir
#         (name)   (start - end)  (sense)(loc) (conf) (hp) (tail) (notes)

#     where 'conf' is the overall confidence score, 'hp' is the hairpin score, and
#     'tail' is the tail score. 'Conf' (which ranges from 0 to 100) is what you
#     probably want to use to assess the quality of a terminator. Higher is better.
#     The confidence, hp score, and tail scores are described in the paper cited
#     above.  'Loc' gives type of region the terminator is in:

#         'G' = in the interior of a gene (at least 50bp from an end),
#         'F' = between two +strand genes,
#         'R' = between two -strand genes,
#         'T' = between the ends of a +strand gene and a -strand gene,
#         'H' = between the starts of a +strand gene and a -strand gene,
#         'N' = none of the above (for the start and end of the DNA)

#     Because of how overlapping genes are handled, these designations are not
#     exclusive. 'G', 'F', or 'R' can also be given in lowercase, indicating that
#     the terminator is on the opposite strand as the region.  Unless the
#     --all-context option is given, only candidate terminators that appear to be in
#     an appropriate genome context (e.g. T, F, R) are output. 

#     Following the TERM line is the sequence of the hairpin and the 5' and 3'
#     tails, always written 5' to 3'.
    
    transterm_table = DataFrames.DataFrame()
    chromosome = ""
    for line in Iterators.filter(x -> occursin(r"^\s*(SEQUENCE|TERM)", x), eachline(transterm_output))
        line = strip(line)
        if occursin(r"^SEQUENCE", line)
            chromosome = split(line)[2]
        else
            transterm_regex = r"(TERM \d+)\s+(\d+) - (\d+)\s+(\S)\s+(\w+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\|(.*)"
            term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes = match(transterm_regex, line).captures
            notes = strip(notes)
            row = (;chromosome, term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes)
            push!(transterm_table, row)
        end
    end
    return transterm_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert TransTerm terminator predictions output to GFF3 format.

Parses TransTerm output and generates a standardized GFF3 file with the following transformations:
- Sets source field to "transterm"
- Sets feature type to "terminator"  
- Converts terminator IDs to GFF attributes
- Renames fields to match GFF3 spec

# Arguments
- `transterm_output::String`: Path to the TransTerm output file

# Returns
- `String`: Path to the generated GFF3 file (original filename with .gff extension)
"""
function transterm_output_to_gff(transterm_output)
    transterm_table = parse_transterm_output(transterm_output)
    
    final_gff_table = if DataFrames.nrow(transterm_table) == 0
        # If the input table is empty, create an empty DataFrame with the correct GFF columns and types.
        DataFrames.DataFrame(
            "#seqid" => String[],
            "source" => String[],
            "type" => String[],
            "start" => Int[],
            "end" => Int[],
            "score" => Float64[],
            "strand" => String[],
            "phase" => String[],
            "attributes" => String[]
        )
    else
        # If there is data, process it as before.
        transterm_table[!, "source"] .= "transterm"
        transterm_table[!, "type"] .= "terminator"
        transterm_table[!, "phase"] .= "."
        transterm_table[!, "attributes"] = map(x -> "label=" * replace(x, " " => "_"), transterm_table[!, "term_id"])
        DataFrames.rename!(transterm_table,
            "chromosome" => "#seqid",
            "stop" => "end",
            "confidence" => "score",
        )
        # Select and reorder columns for the final GFF structure.
        transterm_table[!, ["#seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    end

    transterm_gff_path = transterm_output * ".gff"
    # Write the resulting DataFrame (either with data or empty with correct columns) to a GFF file.
    # uCSV will correctly write just the header if the DataFrame is empty.
    uCSV.write(transterm_gff_path, final_gff_table, delim='\t')
    
    return transterm_gff_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a VirSorter score TSV file and return a DataFrame.

# Arguments
- `virsorter_score_tsv::String`: The file path to the VirSorter score TSV file.

# Returns
- `DataFrame`: A DataFrame containing the parsed data from the TSV file. If the file is empty, returns a DataFrame with the appropriate headers but no data.
"""
function parse_virsorter_score_tsv(virsorter_score_tsv)
    data, header = uCSV.read(virsorter_score_tsv, delim='\t', header=1)
    if length(data) == 0
        data = [[] for i in 1:length(header)]
    end
    return DataFrames.DataFrame(data, header)
end

"""
    add_trna_attributes!(gff_df::DataFrames.DataFrame)

Modifies a GFF `DataFrame` in-place to add informative `label` and `product`
attributes for tRNA features.

This function parses the tRNA type from the `ID` attribute, constructs
`label` and `product` fields, and appends them to the `attributes` column.
It is designed to work with attribute formats generated by tRNAscan-SE,
for example: `ID=...tRNA13-fMetCAT`.

# Arguments
- `gff_df::DataFrames.DataFrame`: A DataFrame containing GFF data, which will be modified.

# Returns
- `DataFrames.DataFrame`: The modified DataFrame with updated attributes.
"""
function add_trna_attributes(gff_df::DataFrames.DataFrame)
    # This helper function will be applied to each row of the DataFrame.
    function format_attributes(id_string::AbstractString)
        # 1. Extract the tRNA type string from the end of the ID.
        # Example: from "ID=...tRNA13-fMetCAT", get "fMetCAT".
        local trna_type
        try
            # Split the ID by the hyphen to isolate the type.
            trna_type = last(split(id_string, '-'))
        catch e
            # If split fails (e.g., no hyphen), return the original string.
            return id_string
        end

        # 2. Parse the amino acid and anticodon from the type string.
        local amino_acid, anticodon
        if startswith(trna_type, "fMet")
            amino_acid = "fMet"
            anticodon = trna_type[5:end] # The rest of the string after "fMet"
        elseif startswith(trna_type, "Undet")
            amino_acid = "Undet"
            anticodon = trna_type[6:end]
        else
            # Standard case: first 3 letters are the amino acid.
            amino_acid = trna_type[1:3]
            anticodon = trna_type[4:end]
        end

        # 3. Create the new label and product strings.
        # This format is descriptive in genome browsers.
        label = "tRNA-$(amino_acid)"
        product = "tRNA-$(amino_acid) ($(anticodon))"

        # 4. Construct the final attributes string and return it.
        return "label=\"$(label)\";product=\"$(product)\";$(id_string);"
    end

    # Apply the helper function to the 'attributes' column of the DataFrame.
    # `ByRow` ensures the function is applied to each element individually.
    return DataFrames.transform(gff_df, :attributes => DataFrames.ByRow(format_attributes) => :attributes)
end