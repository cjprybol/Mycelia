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