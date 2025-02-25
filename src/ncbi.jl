# function fastq_dump(SRR, outdir=SRR)
#     
#     $(CONDA_RUNNER) run --live-stream -n sra-tools fastq-dump /home/jovyan/workspace/pacbio-metagenomic-datasets/ATCC-MSA-1003/SRR9328980 --split-3 --skip-technical --gzip --outdir /home/jovyan/workspace/pacbio-metagenomic-datasets/ATCC-MSA-1003/SRR9328980


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
function prefetch(;SRR, outdir=pwd())
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
function fasterq_dump(;outdir="", srr_identifier="")
    Mycelia.add_bioconda_env("sra-tools")
    prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
    
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
                --threads $(min(Sys.CPU_THREADS, 4))
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
        reverse_reads = isfile(forward_reads_gz) ? reverse_reads_gz : missing,
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
# function download_and_filter_sra_reads(;outdir="", srr_identifier="")
#     forward_reads = joinpath(outdir, "$(srr_identifier)_1.fastq")
#     reverse_reads = joinpath(outdir, "$(srr_identifier)_2.fastq")
#     forward_reads_gz = forward_reads * ".gz"
#     reverse_reads_gz = reverse_reads * ".gz"
#     trimmed_forward_reads = joinpath(outdir, "trim_galore", "$(srr_identifier)_1_val_1.fq.gz")
#     trimmed_reverse_reads = joinpath(outdir, "trim_galore", "$(srr_identifier)_2_val_2.fq.gz")

#     if !(isfile(trimmed_forward_reads) && isfile(trimmed_reverse_reads))
#         @info "processing $(srr_identifier)"
#         fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
#         trim_galore(outdir=outdir, identifier=srr_identifier)
#     # else
#         # @info "$(srr_identifier) already processed..."
#     end
#     isfile(forward_reads_gz) && rm(forward_reads_gz)
#     isfile(reverse_reads_gz) && rm(reverse_reads_gz)
# end

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
function download_genome_by_accession(;accession, outdir=pwd(), compressed = true)
    temp_fasta = joinpath(outdir, accession * ".fna")
    if compressed
        outfile = temp_fasta * ".gz"
    else
        outfile = temp_fasta
    end
    if !isfile(outfile)
        try
            # pull the entire record so that if the download fails we don't leave an empty file
            fasta_records = collect(Mycelia.get_sequence(db = "nuccore", accession = accession))
            open(temp_fasta, "w") do io
                fastx_io = FASTX.FASTA.Writer(io)
                for fasta_record in fasta_records
                    write(fastx_io, fasta_record)
                end
                close(fastx_io)
                if compressed
                    run(`gzip $(temp_fasta)`)
                end
                @assert isfile(outfile)
            end
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
function download_genome_by_ftp(;ftp, outdir=pwd())
    url = Mycelia.ncbi_ftp_path_to_url(ftp_path=ftp, extension="genomic.fna.gz")
    outfile = joinpath(outdir, basename(url))
    if !isfile(outfile)
        return Downloads.download(url, outfile)
    else
        return outfile
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists available BLAST databases from the specified source.
"""
function list_blastdbs(;source::String="")
    Mycelia.add_bioconda_env("blast")
    @assert source in Set(["", "ncbi", "aws", "gcp"])
    data, header = uCSV.read(open(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl --showall tsv`), delim='\t')
    header = ["NAME", "DESCRIPTION", "SIZE (GB)", "LAST_UPDATED"]
    blast_database_table = DataFrames.DataFrame(data, header)
    blast_database_table[!, "LAST_UPDATED"] = map(x -> Dates.Date(first(split(x, "T")), Dates.DateFormat("yyyy-mm-dd")), blast_database_table[!, "LAST_UPDATED"])
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
function local_blast_database_info(;blastdbs_dir="$(homedir())/workspace/blastdb")
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
    data, header = uCSV.read(open(`blastdbcmd -list $(blastdbs_dir) -list_outfmt $(outfmt_string)`), delim='\t')
    header = collect(values(symbol_header_map))
    df = DataFrames.DataFrame(data, header)
    # remove numbered database fragments from summary results
    df = df[map(x -> !occursin(r"\.\d+$", x), df[!, "BLAST database path"]), :]
    df[!, "human readable size"] = Base.format_bytes.(df[!, "number of bytes"])
    return df
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
function download_blast_db(;db, dbdir="$(homedir())/workspace/blastdb", source="", wait=true)
    Mycelia.add_bioconda_env("blast")
    # try
    #     run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    # catch
    #     run(`apt-get install ncbi-blast+ perl-doc -y`)
    # end
    @assert source in ["", "aws", "gcp", "ncbi"]
    mkpath(dbdir)
    current_directory = pwd()
    cd(dbdir)
    if isempty(source)
        @info "source not provided, letting blast auto-detect fastest download option"
        cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress`
        # cmd = `update_blastdb --decompress $(db)`
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
    run(cmd, wait=wait)
    cd(current_directory)
    return "$(dbdir)/$(db)"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BLAST database to FASTA format.

# Arguments
- `db::String`: Name of the BLAST database to convert (e.g. "nr", "nt")
- `dbdir::String`: Directory containing the BLAST database files
- `compressed::Bool`: Whether to gzip compress the output file
- `outfile::String`: Path for the output FASTA file

# Returns
- Path to the generated FASTA file as String

# Notes
- For "nr" database, output extension will be .faa.gz (protein)
- For "nt" database, output extension will be .fna.gz (nucleotide)
- Requires ncbi-blast+ and perl-doc packages to be installed
"""
function blastdb_to_fasta(;db, dbdir="$(homedir())/workspace/blastdb", compressed=true, outfile="$(dbdir)/$(db).$(string(Dates.today())).fasta.gz")
    # todo add more
    if db == "nr"
        outfile = replace(outfile, r"\.fasta\.gz" => ".faa.gz" )
    elseif db == "nt"
        outfile = replace(outfile, r"\.fasta\.gz" => ".fna.gz" )
    end
    # try
    #     run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    # catch
    #     run(`apt-get install ncbi-blast+ perl-doc -y`)
    # end
    Mycelia.add_bioconda_env("blast")
    # p = pipeline(`$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(dbdir)/$(db) -entry all -outfmt %f`)
    p = pipeline(`$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(dbdir)/$(db) -entry all -outfmt %f`)
    if compressed
        p = pipeline(p, `gzip`)
    end
    run(pipeline(p, outfile))
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BLAST database to a tabular format with sequence and taxonomy information.

# Arguments
- `blastdb::String`: Path to the BLAST database
- `outfile::String=""`: Output file path. If empty, generates name based on input database
- `force::Bool=false`: Whether to overwrite existing output file

# Returns
- `String`: Path to the generated output file (.tsv.gz)

# Output Format
Tab-separated file containing columns:
- sequence SHA256
- sequence
- accession
- gi
- ordinal id
- sequence id
- sequence title
- sequence length
- sequence hash
- taxid
- leaf-node taxids
- membership integer
- common taxonomic name
- common taxonomic names for leaf-node taxids
- scientific name
- scientific names for leaf-node taxids
- BLAST name
- taxonomic super kingdom
- PIG
"""
function blastdb2table(;blastdb, outfile="", force=false)
    # try
    #     run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    # catch
    #     run(`apt-get install ncbi-blast+ perl-doc -y`)
    # end
    Mycelia.add_bioconda_env("blast")
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
    io = `$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -entry 'all' -db $(blastdb) -outfmt $(outfmt_string)`
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

Load and parse the assembly summary metadata from NCBI's FTP server for either GenBank or RefSeq databases.

# Arguments
- `db::String`: Database source, must be either "genbank" or "refseq"

# Returns
- `DataFrame`: Parsed metadata table with properly typed columns including:
  - Integer columns: taxid, species_taxid, genome metrics, and gene counts
  - Float columns: gc_percent
  - Date columns: seq_rel_date, annotation_date
  - String columns: all other fields

# Details
Downloads the assembly summary file from NCBI's FTP server and processes it by:
1. Parsing the tab-delimited file with commented headers
2. Converting numeric strings to proper Integer/Float types
3. Parsing date strings to Date objects
4. Handling missing values throughout
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

Loads NCBI RefSeq metadata into a DataFrame. RefSeq is NCBI's curated collection 
of genomic, transcript and protein sequences.

# Returns
- `DataFrame`: Contains metadata columns including accession numbers, taxonomic information,
and sequence details from RefSeq.
"""
function load_refseq_metadata()
    return load_ncbi_metadata("refseq")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load metadata for GenBank sequences into a DataFrame.

This is a convenience wrapper around `load_ncbi_metadata("genbank")` that
specifically loads metadata from the GenBank database.

# Returns
- `DataFrame`: Contains metadata fields like accession numbers, taxonomy,
and sequence information from GenBank.
"""
function load_genbank_metadata()
    return load_ncbi_metadata("genbank")
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
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and constructs a MetaDiGraph representation of the NCBI taxonomy database.

# Arguments
- `path_to_taxdump`: Directory path where taxonomy files will be downloaded and extracted

# Returns
- `MetaDiGraph`: A directed graph where:
  - Vertices represent taxa with properties:
    - `:tax_id`: NCBI taxonomy identifier
    - `:scientific_name`, `:common_name`, etc.: Name properties
    - `:rank`: Taxonomic rank
    - `:division_id`, `:division_cde`, `:division_name`: Division information
  - Edges represent parent-child relationships in the taxonomy

# Dependencies
Requires internet connection for initial download. Uses DataFrames, MetaGraphs, and ProgressMeter.
"""
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
  
# Returns
- Path to the extracted genome data directory

# Notes
- Requires the ncbi-datasets-cli conda package (automatically installed if missing)
- Downloaded zip file is automatically removed after extraction
- If output folder already exists, download is skipped
- Data is extracted to `outdir/accession/ncbi_dataset/data/accession`
"""
function ncbi_genome_download_accession(;
        accession,
        outdir = pwd(),
        outpath = joinpath(outdir, accession * ".zip"),
        include_string = "gff3,rna,cds,protein,genome,seq-report"
    )
    outfolder = joinpath(outdir, accession)
    if !isdir(outfolder)
        add_bioconda_env("ncbi-datasets-cli")
        if isfile(outpath)
            @info "$(outpath) already exists, skipping download..."
        else
            mkpath(outdir)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome accession $(accession) --include $(include_string) --filename $(outpath) --no-progressbar`)
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
            candidate_genomes = filter(x -> occursin(accession, basename(x)) && occursin(Mycelia.FASTA_REGEX, basename(x)), readdir(final_outfolder, join=true))
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

function get_ncbi_dataset_filename(item)
    filenames = Dict(
        "genome" => "GCF_000819615.1_ViralProj14015_genomic.fna",
        "cds" => "cds_from_genomic.fna",
        "gff3" => "genomic.gff",
        "protein" => "protein.faa",
        "seq-report" => "sequence_report.jsonl"
    )
    return get(filenames, item, "")
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
function get_genbank(;db=""::String, accession=""::String, ftp=""::String)
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
function setup_taxonkit_taxonomy()
    run(`wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`)
    Mycelia.tar_extract(tarchive="taxdump.tar.gz", directory=mkpath("$(homedir())/.taxonkit"))
end