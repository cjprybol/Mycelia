"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function prefetch(;SRR, outdir=pwd())
    Mycelia.add_bioconda_env("sra-tools")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n sra-tools prefetch $(SRR) -O $(outdir)`)
end

# function fastq_dump(SRR, outdir=SRR)
#     Mycelia.add_bioconda_env("sra-tools")
#     $(CONDA_RUNNER) run --live-stream -n sra-tools fastq-dump /home/jovyan/workspace/pacbio-metagenomic-datasets/ATCC-MSA-1003/SRR9328980 --split-3 --skip-technical --gzip --outdir /home/jovyan/workspace/pacbio-metagenomic-datasets/ATCC-MSA-1003/SRR9328980




# taxon=694009
# $(CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome taxon 694009 --dehydrated --filename /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009.zip
# unzip /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009.zip -d /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009
# $(CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets rehydrate --directory /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009



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

Download mmseqs databases

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
function download_mmseqs_db(;db, dbdir="$(homedir())/workspace/mmseqs", force=false, wait=true)
    Mycelia.add_bioconda_env("mmseqs2")
    mkpath(dbdir)
    # sanitized_db = replace(db, "/" => "_")
    db_path = joinpath(dbdir, db)
    mkpath(dirname(db_path))
    if !isfile(db_path) || force
        cmd = `$(CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs databases --compressed 1 $(db) --remove-tmp-files 1 $(dbdir)/$(db) $(dbdir)/tmp`
        @time run(cmd, wait=wait)
    else
        @info "db $db @ $(db_path) already exists, set force=true to overwrite"
    end
    return db_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function list_blastdbs(;source="ncbi")
    # Mycelia.add_bioconda_env("blast")
    # readlines(`$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl --showall`)
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    readlines(`update_blastdb --source $(source) --showall`)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function showall_blastdbs(;source="ncbi")
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    blast_table_header = filter(!isempty, split(readlines(`update_blastdb --source $(source) --showall pretty`)[2], "  "))
    data, header = uCSV.read(IOBuffer(join(readlines(`update_blastdb --source $(source) --showall tsv`)[2:end], "\n")), delim="\t")
    df = sort(DataFrames.DataFrame(data, blast_table_header), "SIZE (GB)", rev=true)
    df[!, "LAST_UPDATED"] = map(dt -> Dates.Date(Dates.DateTime(first(split(dt, '.')), "yyyy-mm-ddTHH:MM:SS")), df[!, "LAST_UPDATED"])
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function local_blast_database_info(;blastdbs_dir="$(homedir())/workspace/blastdb")
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
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
    data, header = uCSV.read(open(`blastdbcmd -list $(blastdbs_dir) -list_outfmt $(outfmt_string)`), delim='\t')
    header = collect(values(symbol_header_map))
    df = DataFrames.DataFrame(data, header)
    # remove numbered database fragments from summary results
    df = df[map(x -> !occursin(r"\.\d+$", x), df[!, "BLAST database path"]), :]
    df[!, "human readable size"] = Base.format_bytes.(df[!, "number of bytes"])
    return df
end

# function remove_non_ascii(input::String)
#     return String(filter(x -> x <= '\x7f', input))
# end

# function sanitize_string(input::String)
#     return String(filter(x -> isvalid(Char, x), input))
# end

function blastdb2table(;blastdb, outfile="", force=false)
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    blast_db_info = Mycelia.local_blast_database_info()
    # @info "local blast databases found"
    # display(blast_db_info)
    # @show blastdb
    filtered_blast_db_table = blast_db_info[blast_db_info[!, "BLAST database path"] .== blastdb, :]
    @assert DataFrames.nrow(filtered_blast_db_table) == 1
    blast_db_info = filtered_blast_db_table[1, :]
    n_sequences = blast_db_info["number of sequences"]
    if blast_db_info["BLAST database molecule type"] == "Protein"
        extension = ".faa"
    elseif blast_db_info["BLAST database molecule type"] == "Nucleotide"
        extension = ".fna"
    else
        @show blast_db_info["BLAST database molecule type"]
        error("unexpected blast database molecule type")
    end
    if outfile == ""
        outfile = blastdb * extension * ".tsv.gz"
    end
    lets_go = false
    if !isfile(outfile)
        lets_go = true
    elseif filesize(outfile) == 0
        lets_go = true
    elseif force
        lets_go = true
    else
        @show isfile(outfile)
        @show Mycelia.filesize_human_readable(outfile)
    end
    !lets_go && return outfile
    
    symbol_header_map = OrderedCollections.OrderedDict(
        "%s" => "sequence",
        "%a" => "accession",
        "%g" => "gi",
        "%o" => "ordinal id",
        "%i" => "sequence id",
        "%t" => "sequence title",
        "%l" => "sequence length",
        "%h" => "sequence hash",
        "%T" => "taxid",
        "%X" => "leaf-node taxids",
        "%e" => "membership integer",
        "%L" => "common taxonomic name",
        "%C" => "common taxonomic names for leaf-node taxids",
        "%S" => "scientific name",
        "%N" => "scientific names for leaf-node taxids",
        "%B" => "BLAST name",
        "%K" => "taxonomic super kingdom",
        "%P" => "PIG"
    )
    outfmt_string = join(collect(keys(symbol_header_map)), "\t")
    # @show outfmt_string
    outfile_io = CodecZlib.GzipCompressorStream(open(outfile, "w"))
    
    header = collect(values(symbol_header_map))
    header = ["sequence SHA256", header...]
    println(outfile_io, join(header, '\t'))
    p = ProgressMeter.Progress(n_sequences, desc="Processing $(n_sequences) records from Blast DB $(blastdb): ")
    # io = open(pipeline(`blastdbcmd -entry 'all' -db $(blastdb) -outfmt $(outfmt_string)`, `head`))
    io = `blastdbcmd -entry 'all' -db $(blastdb) -outfmt $(outfmt_string)`
    for line in eachline(io)
        split_line = split(line, '\t')
        # remove invalid characters
        seq = uppercase(String(filter(x -> isvalid(Char, x), split_line[1])))
        seq_sha256 = Mycelia.seq2sha256(seq)
        updated_line = join([seq_sha256, split_line...], "\t")
        println(outfile_io, updated_line)
        ProgressMeter.next!(p)
    end
    close(outfile_io)
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Smart downloading of blast dbs depending on interactive, non interactive context

For a list of all available databases, run: `Mycelia.list_blastdbs()`
"""
function download_blast_db(;db, dbdir="$(homedir())/workspace/blastdb", source="", wait=true)
    # Mycelia.add_bioconda_env("blast")
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    @assert source in ["", "aws", "gcp", "ncbi"]
    mkpath(dbdir)
    current_directory = pwd()
    cd(dbdir)
    if isempty(source)
        @info "source not provided, letting blast auto-detect fastest download option"
        # cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress`
        cmd = `update_blastdb --decompress $(db)`
    else
        @info "downloading from source $(source)"
        if source == "ncbi"
            # --timeout 360 --passive no 
            # cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source) --timeout 360 --passive no`
            cmd = `update_blastdb --timeout 360 --passive no --decompress --source $(source) $(db)`
        else
            # cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source)`
            cmd = `update_blastdb --decompress --source $(source) $(db)`
        end
    end
    run(cmd, wait=wait)
    cd(current_directory)
    return "$(dbdir)/$(db)"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function blastdb_to_fasta(;db, dbdir="$(homedir())/workspace/blastdb", compressed=true, outfile="$(dbdir)/$(db).$(string(Dates.today())).fasta.gz")
    # todo add more
    if db == "nr"
        outfile = replace(outfile, r"\.fasta\.gz" => ".faa.gz" )
    elseif db == "nt"
        outfile = replace(outfile, r"\.fasta\.gz" => ".fna.gz" )
    end
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    # p = pipeline(`$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(dbdir)/$(db) -entry all -outfmt %f`)
    p = pipeline(`blastdbcmd -db $(dbdir)/$(db) -entry all -outfmt %f`)
    if compressed
        p = pipeline(p, `gzip`)
    end
    run(pipeline(p, outfile))
    return outfile
end

# neo_import_dir = "/Users/cameronprybol/Library/Application Support/Neo4j Desktop/Application/relate-data/dbmss/dbms-8ab8baac-5dea-4137-bb24-e0b426447940/import"


# uploading over API is slow for remote and local connections
# Progress:   0%|▏                                        |  ETA: 8:13:39
# upload_nodes_over_api(graph, ADDRESS=local_neo4j_bolt_address, PASSWORD=local_neo4j_password)
# Progress:   0%|         

# # push to Neo4J Aura
# # run(`sudo touch /etc/neo4j/neo4j.conf`)
# run(`sudo neo4j stop`)
# # remote database needs to be running
# # needs to be big enough
# # leave off port from address
# # run(`neo4j-admin push-to-cloud --overwrite --verbose --bolt-uri=$(ADDRESS) --username=$(USERNAME) --password=$(PASSWORD)`)
# # run(`sudo neo4j-admin push-to-cloud --overwrite --verbose --dump-to "$(DIR)/test.db.dump" --bolt-uri=$(a) --username=$(USERNAME) --password=$(PASSWORD)`)
# run(`sudo neo4j-admin push-to-cloud --overwrite --verbose --bolt-uri=$(a) --username=$(USERNAME) --password=$(PASSWORD)`)


# cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.tax_id IS UNIQUE"
# @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.`scientific name` IS UNIQUE"
# @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.identifier IS UNIQUE"
# @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# parameters = ["$(n): row.$(n)" for n in filter(x -> x != "TYPE", names(node_table))]
# parameters = "{" * join(parameters, ", ") * "}"

# window_size = 10000
# V = DataFrames.nrow(node_table)
# windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
# ProgressMeter.@showprogress for (i, w) in enumerate(windows)
#     df_sub = node_table[w, :]
#     f = "node$i.tsv"
#     local_f_path = "$(temp_upload_dir)/$(f)"
#     uCSV.write(local_f_path, df_sub, delim='\t')
#     run(`chmod 777 $(local_f_path)`)
#     f_url = "file:///$(local_f_path)"
#     cmd =
#     """
#     LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
#     CREATE (node:$(NODE_TYPE) $(parameters))
#     """
#     cmd = rstrip(replace(cmd, '\n' => ' '))
#     cypher_cmd = Mycelia.cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd)
#     run(cypher_cmd) 
# end


# ProgressMeter.@showprogress for (i, w) in enumerate(windows)
#     df_sub = edge_table[w, :]
#     f = "edge$i.tsv"
#     local_f_path = "$(temp_upload_dir)/$(f)"
#     uCSV.write(local_f_path, df_sub, delim='\t')
#     run(`chmod 777 $(local_f_path)`)
#     f_url = "file:///$(local_f_path)"
#     cmd = 
#     """
#     LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
#     MATCH (src:$(src_type) {identifier: row.src})
#     MATCH (dst:$(dst_type) {identifier: row.dst})
#     MERGE (src)-[p:$(edge_type)]->(dst)
#     """
#     cmd = rstrip(replace(cmd, '\n' => ' '))
#     cypher_cmd = Mycelia.cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd)
#     run(cypher_cmd) 
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function upload_nodes_to_neo4j(;graph, address, username="neo4j", password, format="auto", database="neo4j", neo4j_import_directory)
    
    node_types = unique(MetaGraphs.props(graph, v)[:TYPE] for v in Graphs.vertices(graph))
    # node_type_strings = Mycelia.type_to_string.(node_types)
    
    for node_type in node_types
        @info "uploading node_type => $(Mycelia.type_to_string(node_type))..."
        node_type_table = node_type_to_dataframe(node_type=node_type, graph=graph)
        try
            upload_node_table(table=node_type_table, address=address, password=password, neo4j_import_dir=neo4j_import_directory)
        catch e
            showerror(stdout, e)
        end
    end
    
    @info "done!"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function node_type_to_dataframe(;node_type, graph)
    node_type_indices = filter(v -> MetaGraphs.props(graph, v)[:TYPE] == node_type, Graphs.vertices(graph))
    node_type_parameters = unique(reduce(vcat, map(v -> collect(keys(MetaGraphs.props(graph, v))), node_type_indices)))
    node_type_table = DataFrames.DataFrame(Dict(p => [] for p in node_type_parameters))
    for node_index in node_type_indices      
        push!(node_type_table, MetaGraphs.props(graph, node_index))
    end
    # normalize
    node_type_table[!, "TYPE"] = Mycelia.type_to_string.(node_type_table[!, "TYPE"])
    for column in names(node_type_table)
        T = Union{unique(typeof.(node_type_table[!, column]))...}
        if T <: AbstractDict
            node_type_table[!, column] = map(d -> JSON.json(string(JSON.json(d))), node_type_table[!, column])
        else
            node_type_table[!, column] = JSON.json.(string.(node_type_table[!, column]))
        end
    end  
    return node_type_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function upload_node_table(;table, window_size=1000, address, password, username="neo4j", database="neo4j", neo4j_import_dir)
    nrows = DataFrames.nrow(table)
    windows = (i:min(i+window_size-1,nrows) for i in 1:window_size:nrows)
    
    node_types = unique(table[!, "TYPE"])
    @assert length(node_types) == 1
    NODE_TYPE = Mycelia.type_to_string(first(node_types))
    parameters = ["$(n): row.$(n)" for n in filter(x -> !(x in ["TYPE"]), names(table))]
    parameters = "{" * join(parameters, ", ") * "}"

    ProgressMeter.@showprogress for (i, window) in enumerate(windows)
        df_sub = table[window, :]
        f = "node$i.tsv"
        local_f_path = "$(neo4j_import_dir)/$(f)"
        uCSV.write(local_f_path, df_sub, delim='\t')
        run(`chmod 777 $(local_f_path)`)
        f_url = "file:///$(f)"
        cmd =
        """
        LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
        CREATE (:`$(NODE_TYPE)` $(parameters))
        """
        cmd = rstrip(replace(cmd, '\n' => ' '))
        cypher_cmd = Mycelia.cypher(cmd, address = address, username = username, password = password, database = database)
        run(cypher_cmd) 
    end
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function upload_node_over_api(graph, v; ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j")
    node_type = MetaGraphs.props(graph, v)[:TYPE]
    node_identifier = MetaGraphs.props(graph, v)[:identifier]
    node_parameters = filter(x -> 
            !(x[1] in (:TYPE, :identifier)) && 
            !(ismissing(x[2]) || isempty(x[2])), 
        MetaGraphs.props(graph, v))
    params_string = join(["$(string(key)): \"$(string(value))\"" for (key, value) in node_parameters], ", ")
    node_type_string = Mycelia.type_to_string(node_type)
    node_identifier_string = string(node_identifier)
    cmd = 
    """
    MERGE (`$(node_identifier_string)`:`$(node_type_string)` {$(params_string)})
    """
    cmd = strip(cmd)
    cypher_cmd = Mycelia.cypher(cmd, address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE)
    run(cypher_cmd)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function upload_nodes_over_api(graph; ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j")
    ProgressMeter.@showprogress for v in Graphs.vertices(graph)
        upload_node_over_api(graph, v, ADDRESS=ADDRESS, USERNAME=USERNAME, PASSWORD=PASSWORD, DATABASE=DATABASE)
    end    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function create_node_constraints(graph; address, username="neo4j", password, database="neo4j")
    node_types = unique(MetaGraphs.props(graph, v)[:TYPE] for v in Graphs.vertices(graph))
    node_type_strings = map(t -> Mycelia.type_to_string(t), node_types)
    for t in node_type_strings
        cmd = "CREATE CONSTRAINT ON (t:`$(t)`) ASSERT t.identifier IS UNIQUE"
        try
            cypher = Mycelia.cypher(cmd, address = address, username = username, password = password, database = database)
            @show cypher
            @time run(cypher)
        catch
            continue
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
# function type_to_string(T::KMER_TYPE) where {KMER_TYPE <: Kmers.Kmer}
#     @show "here"
#     return 
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function type_to_string(T)
    if T <: Kmers.Kmer
        return "Kmers.DNAKmer{$(T.parameters[2])}"
    else
        return string(T)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function type_to_string(T::AbstractString)
    return string(T)
end


# function batch_upload_edge_type_over_url_from_graph(src_type, dst_type, edge_type, graph, ADDRESS, USERNAME, PASSWORD, DATABASE; window_size=100)    
#     src_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == src_type, Graphs.vertices(graph))
#     dst_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == dst_type, Graphs.vertices(graph))
#     edges_to_upload = []
#     for src_node in src_nodes
#         outneighbors = Graphs.outneighbors(graph, src_node)
#         outneighbors = filter(outneighbor -> outneighbor in dst_nodes, outneighbors)
#         for outneighbor in outneighbors
#             this_edge = Graphs.Edge(src_node, outneighbor)
#             @assert MetaGraphs.get_prop(graph, this_edge, :TYPE) == edge_type
#             push!(edges_to_upload, this_edge)
#         end
#     end
    
#     N = length(edges_to_upload)
#     windows = [i:min(i+window_size-1,N) for i in 1:window_size:N]
    
#     ProgressMeter.@showprogress for window in windows
#         cmds = []
#         for (i, e) in enumerate(edges_to_upload[window])
#             edge_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, e)))
#             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, e, param))'" for param in edge_params]
#             joint_params = join(params, ", ")
#             node_cmds = 
#             """
#             MERGE (src$(i):$(MetaGraphs.props(graph, e.src)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.src)[:identifier])'})
#             MERGE (dst$(i):$(MetaGraphs.props(graph, e.dst)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.dst)[:identifier])'})
#             """
#             if !isempty(joint_params)
#                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE]) {$(joint_params)}]->(dst$(i))"
#             else
#                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE])]->(dst$(i))"
#             end
#             cmd = node_cmds * relationship_cmd
#             cmd = replace(cmd, '\n' => ' ')
#             push!(cmds, cmd)
#         end
#         cmd = join(cmds, ' ')
#         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
#         run(cypher_cmd)
#     end    
# end

# function batch_upload_node_type_over_url_from_graph(node_type, graph, ADDRESS, USERNAME, PASSWORD, DATABASE, window_size=100)
#     node_type_params = Set{Symbol}()
#     vertices_of_type = [v for v in Graphs.vertices(graph) if (graph.vprops[v][:TYPE] == node_type)]
    
#     V = length(vertices_of_type)
#     windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
    
#     ProgressMeter.@showprogress for window in windows
#         cmds = []
#         for (i, v) in enumerate(vertices_of_type[window])
#             node_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, v)))
#             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, v, param))'" for param in node_params]
#             # params = ["$(string(param)):'$(escape_string(MetaGraphs.get_prop(graph, v, param)))'" for param in node_params]
#             joint_params = join(params, ", ")
#             cmd = "MERGE (node$(i):$(node_type) {$(joint_params)})"
#             push!(cmds, cmd)
#         end
#         cmd = join(cmds, ' ')
#         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
#         run(cypher_cmd)
#     end    
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function upload_node_type_over_url_from_graph(;node_type, graph, ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j", window_size=100)
    node_type_params = Set{Symbol}()
    vertices_of_type = [v for v in Graphs.vertices(graph) if (graph.vprops[v][:TYPE] == node_type)]
    
    V = length(vertices_of_type)
    windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
    
    ProgressMeter.@showprogress for window in windows
        cmds = []
        for (i, v) in enumerate(vertices_of_type[window])
            node_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, v)))
            params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, v, param))'" for param in node_params]
            # params = ["$(string(param)):'$(escape_string(MetaGraphs.get_prop(graph, v, param)))'" for param in node_params]
            joint_params = join(params, ", ")
            cmd = "MERGE (node$(i):$(node_type) {$(joint_params)})"
            push!(cmds, cmd)
        end
        cmd = join(cmds, ' ')
        cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
        run(cypher_cmd)
    end    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function upload_edge_type_over_url_from_graph(;src_type, dst_type, edge_type, graph, ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j", window_size=100)    
    src_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == src_type, Graphs.vertices(graph))
    dst_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == dst_type, Graphs.vertices(graph))
    edges_to_upload = []
    for src_node in src_nodes
        outneighbors = Graphs.outneighbors(graph, src_node)
        outneighbors = filter(outneighbor -> outneighbor in dst_nodes, outneighbors)
        for outneighbor in outneighbors
            this_edge = Graphs.Edge(src_node, outneighbor)
            @assert MetaGraphs.get_prop(graph, this_edge, :TYPE) == edge_type
            push!(edges_to_upload, this_edge)
        end
    end
    
    N = length(edges_to_upload)
    windows = (i:min(i+window_size-1,N) for i in 1:window_size:N)
    
    ProgressMeter.@showprogress for window in windows
        cmds = []
        for (i, e) in enumerate(edges_to_upload[window])
            edge_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, e)))
            params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, e, param))'" for param in edge_params]
            joint_params = join(params, ", ")
            node_cmds = 
            """
            MERGE (src$(i):$(MetaGraphs.props(graph, e.src)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.src)[:identifier])'})
            MERGE (dst$(i):$(MetaGraphs.props(graph, e.dst)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.dst)[:identifier])'})
            """
            if !isempty(joint_params)
                relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE]) {$(joint_params)}]->(dst$(i))"
            else
                relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE])]->(dst$(i))"
            end
            cmd = node_cmds * relationship_cmd
            cmd = replace(cmd, '\n' => ' ')
            push!(cmds, cmd)
        end
        cmd = join(cmds, ' ')
        cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
        run(cypher_cmd)
    end    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function cypher(cmd;
    address="neo4j://localhost:7687",
    username="neo4j",
    password="password",
    format="auto",
    database="neo4j"
    )
    # cmd = `cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
    cmd = `/home/jovyan/.local/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
#    cmd = `/home/jupyter-cjprybol/software/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
    return cmd
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function list_databases(;address, username="neo4j", password)
    cmd = "show databases"
    database = "system"
    cmd = cypher(cmd, address=address, username=username, password=password, database=database)
    return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, quotes='"', encodings=Dict("FALSE" => false, "TRUE" => true))...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function create_database(;database, address, username="neo4j", password)
    current_databases = list_databases(;address, username, password)
    if database in current_databases[!, "name"]
        return
    else
        f = run
        cmd = "create database $(database)"
        # switch database to system, so that we can create the user-specific database in the system
        database = "system"
        run(cypher(;address, username, password, database, cmd, f))
    end
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function cypher(;address, username, password, database, cmd)
#     return `cypher-shell --address $address --username $username --password $password --database $(database) --format auto $(cmd)`
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function load_ncbi_metadata(db)
    if !(db in ["genbank", "refseq"])
        error()
    end
    ncbi_summary_url = "https://ftp.ncbi.nih.gov/genomes/$(db)/assembly_summary_$(db).txt"
    # ncbi_summary_file = basename(ncbi_summary_url)
    # if !isfile(ncbi_summary_file)
    #     download(ncbi_summary_url, ncbi_summary_file)
    # end
    buffer = IOBuffer(HTTP.get(ncbi_summary_url).body)
    # types=[]
    # ncbi_summary_table = DataFrames.DataFrame(uCSV.read(ncbi_summary_file, comment = "## ", header=1, delim='\t', encodings=Dict("na" => missing), allowmissing=true, typedetectrows=100)...)
    ncbi_summary_table = DataFrames.DataFrame(uCSV.read(buffer, comment = "## ", header=1, delim='\t', types=String)...)
    ints = [
        "taxid",
        "species_taxid",
        "genome_size",
        "genome_size_ungapped",
        "replicon_count",
        "scaffold_count",
        "contig_count",
        "total_gene_count",
        "protein_coding_gene_count",
        "non_coding_gene_count"
    ]
    floats = ["gc_percent"]
    dates = ["seq_rel_date", "annotation_date"]
    for int in ints
        ncbi_summary_table[!, int] = something.(tryparse.(Int, ncbi_summary_table[!, int]), missing)
    end
    for float in floats
        ncbi_summary_table[!, float] = something.(tryparse.(Float64, ncbi_summary_table[!, float]), missing)
    end
    for date in dates
        # ncbi_summary_table[!, date] = Dates.Date.(ncbi_summary_table[!, date], Dates.dateformat"yyyy/mm/dd")
        parsed_dates = map(date_string -> tryparse(Dates.Date, date_string, Dates.dateformat"yyyy/mm/dd"), ncbi_summary_table[!, date])
        ncbi_summary_table[!, date] = something.(parsed_dates, missing)
    end
    return ncbi_summary_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function taxonomic_id_to_children(tax_id; DATABASE_ID, USERNAME="neo4j", PASSWORD)
    DATABASE = "neo4j"
    ADDRESS="neo4j+s://$(database_id).databases.neo4j.io:7687"
    
    # NOTE! *, or 0 distance (e.g. [*0..2]) step range will include source node!!!!
    cmd = "MATCH (n)<-[*]-(n2) WHERE n.tax_id IS NOT NULL AND n.tax_id = \"$(tax_id)\" RETURN DISTINCT n2.tax_id AS tax_id"
    println(cmd)
    
    cypher = cypher(cmd, address=ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE)
    tax_ids = readlines(open(cypher))[2:end]
    tax_ids = strip.(tax_ids, '"')
    tax_ids = parse.(Int, tax_ids)
    return unique(tax_ids)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function load_refseq_metadata()
    return load_ncbi_metadata("refseq")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Description

```jldoctest
julia> 1 + 1
2
```
"""
function load_genbank_metadata()
    return load_ncbi_metadata("genbank")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

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

```jldoctest
julia> 1 + 1
2
```
"""
function ncbi_ftp_path_to_url(;ftp_path, extension)
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

Description

```jldoctest
julia> 1 + 1
2
```
"""
function countmap_columns(table)
    for n in names(refseq_metadata)
        display(n)
        display(StatsBase.countmap(refseq_metadata[!, n]))
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

```jldoctest
julia> 1 + 1
2
```
"""
function get_sequence(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/3 second sleep to set max of 3 requests per second when looping
        sleep(0.34)
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=fasta&id=$(accession)"
        body = HTTP.get(url).body
        try
            return FASTX.FASTA.Reader(IOBuffer(body))
        catch e
            @error e body
        end
    elseif !isempty(ftp)
        return FASTX.FASTA.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
    else
        @error "invalid call"
    end
end

# function ncbi_datasets_download_by_taxon_id

function load_ncbi_taxonomy(;
        path_to_taxdump="$(homedir())/workspace/blastdb/taxdump"
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
    ProgressMeter.@showprogress for line in split(read(open("$(taxdump_out)/names.dmp"), String), "\t|\n")
        if isempty(line)
            continue
        else
            (tax_id_string, name_txt, unique_name, name_class) = split(line, "\t|\t")
            tax_id = parse(Int, tax_id_string)
            row = (;tax_id, name_txt, unique_name, name_class)
            push!(names_dmp, row)
        end
    end
    unique_tax_ids = sort(unique(names_dmp[!, "tax_id"]))

    ncbi_taxonomy = MetaGraphs.MetaDiGraph(length(unique_tax_ids))
    ProgressMeter.@showprogress for (index, group) in enumerate(collect(DataFrames.groupby(names_dmp, "tax_id")))
        MetaGraphs.set_prop!(ncbi_taxonomy, index, :tax_id, group[1, "tax_id"])
        for row in DataFrames.eachrow(group)
            unique_name = isempty(row["unique_name"]) ? row["name_txt"] : row["unique_name"]
            # remove quotes since neo4j doesn't like them
            unique_name = replace(unique_name, '"' => "")
            # replace spaces and dashes with underscores
            name_class = Symbol(replace(replace(row["name_class"], r"\s+" => "-"), "-" => "_"))
    #         name_class = Symbol(row["name_class"])
            if haskey(MetaGraphs.props(ncbi_taxonomy, index), name_class)
                current_value = MetaGraphs.get_prop(ncbi_taxonomy, index, name_class)
                if (current_value isa Array) && !(unique_name in current_value)
                    new_value = [current_value..., unique_name]
                    MetaGraphs.set_prop!(ncbi_taxonomy, index, name_class, new_value)
                elseif !(current_value isa Array) && (current_value != unique_name)
                    new_value = [current_value, unique_name]
                    MetaGraphs.set_prop!(ncbi_taxonomy, index, name_class, new_value)
                else
                    continue
                end
            else
                MetaGraphs.set_prop!(ncbi_taxonomy, index, name_class, unique_name)
            end
        end
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

    node_2_taxid_map = map(index -> ncbi_taxonomy.vprops[index][:tax_id], Graphs.vertices(ncbi_taxonomy))
    ProgressMeter.@showprogress for line in split(read(open("$(taxdump_out)/nodes.dmp"), String), "\t|\n")
        if isempty(line)
            continue
        else
            (tax_id_string, parent_tax_id_string, rank, embl_code, division_id_string) = split(line, "\t|\t")

            division_id = parse(Int, division_id_string)

            tax_id = parse(Int, tax_id_string)
            lightgraphs_tax_ids = searchsorted(node_2_taxid_map, tax_id)
            @assert length(lightgraphs_tax_ids) == 1
            lightgraphs_tax_id = first(lightgraphs_tax_ids)

            parent_tax_id = parse(Int, parent_tax_id_string)
            lightgraphs_parent_tax_ids = searchsorted(node_2_taxid_map, parent_tax_id)
            @assert length(lightgraphs_parent_tax_ids) == 1
            lightgraphs_parent_tax_id = first(lightgraphs_parent_tax_ids)

            Graphs.add_edge!(ncbi_taxonomy, lightgraphs_tax_id, lightgraphs_parent_tax_id)
            MetaGraphs.set_prop!(ncbi_taxonomy, lightgraphs_tax_id, :rank, rank)
            # these should probably be broken out as independent nodes!
            MetaGraphs.set_prop!(ncbi_taxonomy, lightgraphs_tax_id, :division_id, division_id)
            MetaGraphs.set_prop!(ncbi_taxonomy, lightgraphs_tax_id, :division_cde, divisions[division_id][:division_cde])
            MetaGraphs.set_prop!(ncbi_taxonomy, lightgraphs_tax_id, :division_name, divisions[division_id][:division_name])
        end
    end
    # JLD2 graph killed a colab instance after 200Gb of size!
    # JLD2.save("$(homedir())/workspace/blastdb/taxdump/ncbi_taxonomy.jld2", "ncbi_taxonomy", ncbi_taxonomy)
    # return (;ncbi_taxonomy, path_to_prebuilt_graph)
    return ncbi_taxonomy
end