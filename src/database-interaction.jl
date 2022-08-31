# neo_import_dir = "/Users/cameronprybol/Library/Application Support/Neo4j Desktop/Application/relate-data/dbmss/dbms-8ab8baac-5dea-4137-bb24-e0b426447940/import"

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
    database="system"
    )
    cmd = `cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
#    cmd = `/home/jupyter-cjprybol/software/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
    return cmd
end

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
    buffer = IOBuffer(HTTP.get(ncbi_summary_url).body)
    ncbi_summary_table = DataFrames.DataFrame(uCSV.read(buffer, comment = "#  ", header=1, delim='\t')...)
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
function taxonomic_id_to_children(tax_id; DATABASE_ID, USERNAME, PASSWORD)
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

Description

```jldoctest
julia> 1 + 1
2
```
"""
function ncbi_ftp_path_to_url(ftp_path, extension)
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
        return FASTX.FASTA.Reader(IOBuffer(HTTP.get(url).body))
    elseif !isempty(ftp)
        return FASTX.FASTA.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
    else
        @error "invalid call"
    end
end

# function ncbi_datasets_download_by_taxon_id