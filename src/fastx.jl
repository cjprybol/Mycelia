function run_fastqc(fastq; outdir = replace(fastq, Mycelia.FASTQ_REGEX => "_fastqc"))
    Mycelia.add_bioconda_env("fastqc")
    if !isdir(outdir)
        mkpath(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n fastqc fastqc --outdir $(outdir) $(fastq)`)
    else
        @warn "$outdir already exists"
    end
end

function run_fastqc(forward, reverse; outdir = Mycelia.find_matching_prefix(forward, reverse) * "_fastqc")
    Mycelia.add_bioconda_env("fastqc")
    if !isdir(outdir)
        mkpath(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n fastqc fastqc --outdir $(outdir) $(forward) $(reverse)`)
    else
        @warn "$outdir already exists"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyze sequence duplication rates in a FASTQ file.

This function processes a FASTQ file to quantify both exact sequence duplications and 
canonical duplications (considering sequences and their reverse complements as equivalent).
The function makes two passes through the file: first to count total records, then to
analyze unique sequences.

# Arguments
- `fastq::String`: Path to the input FASTQ file to analyze
- `results_table::String`: Optional. Path where the results will be saved as a tab-separated file.
  Defaults to the same path as the input file but with extension changed to ".duplication_rates.tsv"

# Returns
- `String`: Path to the results table file

# Output
Generates a tab-separated file containing the following metrics:
- `total_records`: Total number of sequence records in the file
- `total_unique_observations`: Count of unique sequence strings
- `total_unique_canonical_observations`: Count of unique canonical sequences 
  (after normalizing for reverse complements)
- `percent_unique_observations`: Percentage of sequences that are unique
- `percent_unique_canonical_observations`: Percentage of sequences that are unique after canonicalization
- `percent_duplication_rate`: Percentage of sequences that are duplicates (100 - percent_unique_observations)
- `percent_canonical_duplication_rate`: Percentage of sequences that are duplicates after canonicalization

# Notes
- If the specified results file already exists and is not empty, the function will
  return early without recomputing.
- Progress is displayed during processing with a progress bar showing speed.

# Example
```julia
# Analyze a FASTQ file and save results to default location
result_path = assess_duplication_rates("data/sample.fastq")

# Specify custom output path
result_path = assess_duplication_rates("data/sample.fastq", results_table="results/duplication_analysis.tsv")
```
"""
function assess_duplication_rates(fastq; results_table=replace(fastq, Mycelia.FASTQ_REGEX => ".duplication_rates.tsv"))
    # @show results_table
    if isfile(results_table) && (filesize(results_table) > 0)
        println("$results_table already exists.")
        display(DataFrames.DataFrame(uCSV.read(results_table, delim='\t', header=1)))
        return results_table
    end
    # First pass: count total records
    println("Counting total records in FASTQ file...")
    total_records = 0
    for _ in Mycelia.open_fastx(fastq)
        total_records += 1
    end
    println("Found $total_records total records.")
    
    # Initialize sets for unique sequences
    unique_observations = Set{String}()
    unique_canonical_observations = Set{String}()
    
    # Create progress meter
    println("Analyzing duplication rates...")
    prog = ProgressMeter.Progress(total_records, desc="Processing: ", showspeed=true)
    
    # Second pass: process records with progress meter
    for record in Mycelia.open_fastx(fastq)
        seq = FASTX.sequence(record)
        # converted_seq = Mycelia.convert_sequence(seq)
        converted_seq = BioSequences.LongDNA{4}(seq)
        push!(unique_observations, string(converted_seq))
        push!(unique_canonical_observations, string(BioSequences.canonical(converted_seq)))
        
        # Update progress meter
        ProgressMeter.next!(prog)
    end
    
    total_unique_observations = length(unique_observations)
    total_unique_canonical_observations = length(unique_canonical_observations)
    percent_unique_observations = (total_unique_observations / total_records) * 100
    percent_unique_canonical_observations = (total_unique_canonical_observations / total_records) * 100
    percent_duplication_rate = 100 - percent_unique_observations
    percent_canonical_duplication_rate = 100 - percent_unique_canonical_observations

    results_df = DataFrames.DataFrame(;total_records,
                total_unique_observations,
                total_unique_canonical_observations,
                percent_unique_observations,
                percent_unique_canonical_observations,
                percent_duplication_rate,
                percent_canonical_duplication_rate
            )
    display(results_df)
    uCSV.write(results_table, results_df, delim='\t')
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compare two FASTA sequences and calculate alignment statistics.

# Arguments
- `reference_fasta::String`: Path to the reference FASTA file
- `query_fasta::String`: Path to the query FASTA file

# Returns
DataFrame with the following columns:
- `alignment_percent_identity`: Percentage of matching bases in alignment
- `total_equivalent_bases`: Number of equivalent bases between sequences
- `total_alignment_length`: Length of the alignment
- `query_length`: Length of query sequence
- `total_variants`: Total number of variants (SNPs + indels)
- `total_snps`: Number of single nucleotide polymorphisms
- `total_indels`: Number of insertions and deletions
- `alignment_coverage_query`: Percentage of query sequence covered
- `alignment_coverage_reference`: Percentage of reference sequence covered
- `size_equivalence_to_reference`: Size ratio of query to reference (%)

# Notes
- Uses minimap2 with progressively relaxed settings (asm5→asm10→asm20)
- Returns empty string values for alignment statistics if no alignment is found
- Requires minimap2 to be installed and accessible in PATH
"""
# uses minimap
function pairwise_minimap_fasta_comparison(;reference_fasta, query_fasta)
    header = [
        "Query",
        "Query length",
        "Query start",
        "Query end",
        "Query strand",
        "Target",
        "Target length",
        "Target start",
        "Target end",
        "Matches",
        "Alignment length",
        "Mapping quality",
        "Cigar",
        "CS tag"]
    
#     asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    results5 = read(`minimap2 -x asm5 --cs -cL $reference_fasta $query_fasta`)
    if !isempty(results5)
        results = results5
    else
        @warn "no hit with asm5, trying asm10"
        results10 = read(`minimap2 -x asm10 --cs -cL $reference_fasta $query_fasta`)
        if !isempty(results10)
            results = results10
        else
            @warn "no hits with asm5 or asm10, trying asm20"
            results20 = read(`minimap2 -x asm20 --cs -cL $reference_fasta $query_fasta`)
            if !isempty(results20)
                results = results20
            end
        end
    end
    if !isempty(results)
        data =  DelimitedFiles.readdlm(IOBuffer(results), '\t')
        data_columns_of_interest = [collect(1:length(header)-2)..., collect(size(data, 2)-1:size(data, 2))...]
        minimap_results = DataFrames.DataFrame(data[:, data_columns_of_interest], header)

        equivalent_matches = reduce(vcat, map(x -> collect(eachmatch(r":([0-9]+)", replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_equivalent_bases = sum(map(match -> parse(Int, first(match.captures)), equivalent_matches))

        insertion_matches = reduce(vcat, map(x -> collect(eachmatch(r"\+([a-z]+)"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_inserted_bases = sum(map(match -> length(first(match.captures)), insertion_matches))
        deletion_matches = reduce(vcat, map(x -> collect(eachmatch(r"\-([a-z]+)"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_deleted_bases = sum(map(match -> length(first(match.captures)), deletion_matches))
        substitution_matches = reduce(vcat, map(x -> collect(eachmatch(r"\*([a-z]{2})"i, replace(x, "cs:Z:" => ""))), minimap_results[!, "CS tag"]))
        total_substituted_bases = length(substitution_matches)
        total_variants = length(insertion_matches) + length(deletion_matches) + length(substitution_matches)
        total_variable_bases = total_inserted_bases + total_deleted_bases + total_substituted_bases

        total_alignment_length = sum(minimap_results[!, "Alignment length"])
        total_matches = sum(minimap_results[!, "Matches"])
        
        alignment_percent_identity = round(total_matches / total_alignment_length * 100, digits=2)
        size_equivalence_to_reference = round(minimap_results[1, "Query length"]/minimap_results[1, "Target length"] * 100, digits=2)
        alignment_coverage_query = round(total_alignment_length / minimap_results[1, "Query length"] * 100, digits=2)
        alignment_coverage_reference = round(total_alignment_length / minimap_results[1, "Target length"] * 100, digits=2)

        results = DataFrames.DataFrame(
            alignment_percent_identity = alignment_percent_identity,
            total_equivalent_bases = total_equivalent_bases,
            total_alignment_length = total_alignment_length,
            query_length = minimap_results[1, "Query length"],
            total_variants = total_variants,
            total_snps = total_substituted_bases,
            total_indels = length(insertion_matches) + length(deletion_matches),
            alignment_coverage_query = alignment_coverage_query,
            alignment_coverage_reference = alignment_coverage_reference,
            size_equivalence_to_reference = size_equivalence_to_reference,
        )
    else
        query_length = length(FASTX.sequence(first(FASTX.FASTA.Reader(open(query_fasta)))))
        target_length = length(FASTX.sequence(first(FASTX.FASTA.Reader(open(reference_fasta)))))
        size_equivalence_to_reference = round(query_length/target_length * 100, digits=2)

        # unable to find any matches
        results = DataFrames.DataFrame(
            alignment_percent_identity = "",
            total_equivalent_bases = "",
            total_alignment_length = "",
            query_length = query_length,
            total_variants = "",
            total_snps = "",
            total_indels = "",
            alignment_coverage_query = 0,
            alignment_coverage_reference = 0,
            size_equivalence_to_reference = size_equivalence_to_reference
        )
    end
    return results
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Trim paired-end FASTQ reads using Trim Galore, a wrapper around Cutadapt and FastQC.

# # Arguments
# - `outdir::String`: Output directory containing input FASTQ files
# - `identifier::String`: Prefix for input/output filenames

# # Input files
# Expects paired FASTQ files in `outdir` named:
# - `{identifier}_1.fastq.gz` (forward reads)
# - `{identifier}_2.fastq.gz` (reverse reads)

# # Output files
# Creates trimmed reads in `outdir/trim_galore/`:
# - `{identifier}_1_val_1.fq.gz` (trimmed forward reads)
# - `{identifier}_2_val_2.fq.gz` (trimmed reverse reads)

# # Dependencies
# Requires trim_galore conda environment:
# """
# function trim_galore(;outdir="", identifier="")
    
#     trim_galore_dir = joinpath(outdir, "trim_galore")
    
#     forward_reads = joinpath(outdir, "$(identifier)_1.fastq.gz")
#     reverse_reads = joinpath(outdir, "$(identifier)_2.fastq.gz")
    
#     trimmed_forward_reads = joinpath(trim_galore_dir, "$(identifier)_1_val_1.fq.gz")
#     trimmed_reverse_reads = joinpath(trim_galore_dir, "$(identifier)_2_val_2.fq.gz")
    
#     # mamba create -n trim_galore -c bioconda trim_galore
#     if !isfile(trimmed_forward_reads) && !isfile(trimmed_reverse_reads)
#         cmd = `conda run -n trim_galore trim_galore --suppress_warn --cores $(min(Sys.CPU_THREADS, 4)) --output_dir $(trim_galore_dir) --paired $(forward_reads) $(reverse_reads)`
#         run(cmd)
#     else
#         @info "$(trimmed_forward_reads) & $(trimmed_reverse_reads) already present"
#     end
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Trim paired-end FASTQ reads using Trim Galore, a wrapper around Cutadapt and FastQC.

# Arguments
- `forward_reads::String`: Path to forward reads FASTQ file
- `reverse_reads::String`: Path to reverse reads FASTQ file
- `outdir::String`: Output directory for trimmed files

# Returns
- `Tuple{String, String}`: Paths to trimmed forward and reverse read files

# Dependencies
Requires trim_galore conda environment:
- `mamba create -n trim_galore -c bioconda trim_galore`
"""
function trim_galore_paired(;forward_reads::String, reverse_reads::String, outdir::String=pwd())
    # Create output directory if it doesn't exist
    trim_galore_dir = mkpath(joinpath(outdir, "trim_galore"))
    
    # Get base filename without path and extension for output naming
    forward_base = basename(forward_reads)
    reverse_base = basename(reverse_reads)
    
    # Construct output filenames according to trim_galore naming convention
    trimmed_forward = joinpath(trim_galore_dir, replace(forward_base, Mycelia.FASTQ_REGEX => "_val_1.fq.gz"))
    trimmed_reverse = joinpath(trim_galore_dir, replace(reverse_base, Mycelia.FASTQ_REGEX => "_val_2.fq.gz"))
    
    if !isfile(trimmed_forward) && !isfile(trimmed_reverse)
        cmd = `$(Mycelia.CONDA_RUNNER) run -n trim_galore trim_galore --suppress_warn --cores $(min(Sys.CPU_THREADS, 4)) --output_dir $(trim_galore_dir) --paired $(forward_reads) $(reverse_reads)`
        run(cmd)
    else
        @info "$(trimmed_forward) & $(trimmed_reverse) already present"
    end
    
    return (;trimmed_forward, trimmed_reverse, outdir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform quality control (QC) filtering and trimming on short-read FASTQ files using fastp.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output FASTQ file.
- `adapter_seq::String`: Adapter sequence to trim.
- `quality_threshold::Int`: Minimum phred score for trimming (default 20).
- `min_length::Int`: Minimum read length to retain (default 50).

# Returns
- `String`: Path to the filtered and trimmed FASTQ file.

# Details
This function uses fastp to remove adapter contamination, trim low‐quality bases from the 3′ end,
and discard reads shorter than `min_length`. It’s a simple wrapper that executes the external fastp command.
"""
function qc_filter_short_reads_fastp(in_fastq::String, out_fastq::String; adapter_seq::String, quality_threshold::Int=20, min_length::Int=50)
    if !isfile(out_fastq)
        Mycelia.add_bioconda_env("fastp")
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastp fastp --in1 $(in_fastq) --out1 $(out_fastq) --adapter_sequence $(adapter_seq) --qualified_quality_phred $(quality_threshold) --length_required $(min_length)`
        run(cmd)
    end
    return out_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform QC filtering on long-read FASTQ files using fastplong.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output FASTQ file.
- `quality_threshold::Int`: Minimum average quality to retain a read (default 10).
- `min_length::Int`: Minimum read length (default 1000).
- `max_length::Int=0`: Maximum read length (default 0, no maximum).

# Returns
- `String`: Path to the filtered FASTQ file.

# Details
This function uses fastplong to filter long reads based on quality and length criteria.
It is optimized for Oxford Nanopore, PacBio, or similar long-read datasets.
"""
function qc_filter_long_reads_fastplong(;
                            in_fastq::String,
                            report_title::String=in_fastq * " fastplong report",
                            out_fastq::String=Mycelia.replace(in_fastq, Mycelia.FASTQ_REGEX => ".fastplong.fq.gz"),
                            html_report::String=out_fastq * ".html",
                            json_report::String=out_fastq * ".json",
                            min_length::Int=1000,
                            max_length::Int=0)
    # Build command with required parameters
    if !isfile(out_fastq)
        Mycelia.add_bioconda_env("fastplong")
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastplong fastplong
                --in $(in_fastq)
                --out $(out_fastq)
                --report_title $(report_title)
                --html $(html_report)
                --json $(json_report)
                --length_required $(min_length)`
        # Add max length if specified
        if max_length > 0
            push!(cmd, "--length_limit")
            push!(cmd, string(max_length))
        end

        run(`$cmd`)
    end
    return out_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Filter and process long reads from a FASTQ file using Filtlong.

This function filters long sequencing reads based on quality and length criteria, 
then compresses the output using pigz.

# Arguments
- `in_fastq::String`: Path to the input FASTQ file.
- `out_fastq::String`: Path to the output filtered and compressed FASTQ file. 
   Defaults to the input filename with ".filtlong.fq.gz" appended.
- `min_mean_q::Int`: Minimum mean quality score for reads to be kept. Default is 20.
- `keep_percent::Int`: Percentage of reads to keep after filtering. Default is 95.

# Returns
- `out_fastq`

# Details
This function uses Filtlong to filter long reads and pigz for compression. It requires
the Bioconda environment for Filtlong to be set up, which is handled internally.
"""
function qc_filter_long_reads_filtlong(;
        in_fastq,
        out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".filtlong.fq.gz"),
        min_mean_q = 20,
        keep_percent = 95
    )
    if !isfile(out_fastq)
        Mycelia.add_bioconda_env("filtlong")
        Mycelia.add_bioconda_env("pigz")
        p1 = pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n filtlong filtlong --min_mean_q $(min_mean_q) --keep_percent $(keep_percent) $(in_fastq)`,
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz`
        )
        p2 = pipeline(p1, out_fastq)
        run(p2)
    end
    return out_fastq
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Perform QC filtering on long-read FASTQ files using chopper.

# # Arguments
# - `in_fastq::String`: Path to the input FASTQ file.
# - `out_fastq::String`: Path to the output FASTQ file.
# - `quality_threshold::Int`: Minimum average quality to retain a read (default 10).
# - `min_length::Int`: Minimum read length (default 1000).

# # Returns
# - `String`: Path to the filtered FASTQ file.

# # Details
# This function uses chopper to discard long reads that do not meet the minimum quality or length thresholds.
# It is intended for Oxford Nanopore or similar long-read datasets.

# # Dependencies
# Requires chopper to be installed via conda
# """
# function qc_filter_long_reads_chopper(in_fastq::String, out_fastq::String; quality_threshold::Int=20, min_length::Int=1000)
#     # chopper reads from STDIN and writes to STDOUT, so we pipe the file
#     Mycelia.add_bioconda_env("chopper")
#     cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n chopper chopper --input (in_fastq) --quality $(quality_threshold) --minlength $(min_length)`
#     p = pipeline(`cat $(in_fastq)`, cmd)
#     , out_fastq)
#     open(out_fastq, "w") do outf
#         run(pipeline(`cat $(in_fastq)`, cmd; stdout=outf))
#     end
#     return out_fastq
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate basic statistics for FASTQ/FASTA sequence files using seqkit.

# Arguments
- `fastq::String`: Path to input FASTQ/FASTA file

# Details
Automatically installs and uses seqkit from Bioconda to compute sequence statistics
including number of sequences, total bases, GC content, average length, etc.

# Dependencies
- Requires Conda and Bioconda channel
- Installs seqkit package if not present

# Returns
Returns a DataFrame of the table

https://bioinf.shenwei.me/seqkit/usage/#stats
"""
function fastx_stats(fastx)
    Mycelia.add_bioconda_env("seqkit")
    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit stats --N 90 --all --tabular $(fastx)`
    return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, delim='\t'))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a FASTX (FASTA/FASTQ) file into a normalized tab-separated table format with standardized sequence identifiers.

# Arguments
- `fastx_file::String`: Path to input FASTX file
- `outfile::String`: Path to output compressed TSV file (defaults to input filename + ".tsv.gz")
- `force::Bool=false`: If true, overwrites existing output file

# Returns
- `String`: Path to the created output file

# Output Format
Creates a gzipped TSV file with the following columns:
- fasta_identifier: Original FASTA filename
- sequence_sha256: SHA256 hash of the sequence
- sequence_identifier: Original sequence ID from FASTA
- sequence_description: Full sequence description from FASTA
- sequence: The actual sequence
"""
function fastx2normalized_table(fastx)

    # Mycelia.add_bioconda_env("seqkit")
    # cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit
    #     fx2tab
    # return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, delim='\t'))

    @assert isfile(fastx) && filesize(fastx) > 0
    normalized_table = DataFrames.DataFrame(
        record_identifier = String[],
        record_description = String[],
        record_sha256 = String[],
        record_quality = Union{Vector{Float64}, Missing}[],
        record_alphabet = String[],
        record_type = Symbol[],
        mean_record_quality = Union{Float64, Missing}[],
        median_record_quality = Union{Float64, Missing}[],
        record_length = Int[],
        record_sequence = String[],
    )

    file_type = :unknown
    if occursin(Mycelia.FASTA_REGEX, fastx)
        @info "Processing FASTA file"
        file_type = :fasta
    elseif occursin(Mycelia.FASTQ_REGEX, fastx)
        @info "Processing FASTQ file"
        file_type = :fastq
    else
        error("File is not FASTA or FASTQ")
    end
    
    for record in Mycelia.open_fastx(fastx)
        record_sequence = FASTX.sequence(record)
        if file_type == :fasta
            record_quality = missing
        else
            record_quality = collect(FASTX.quality_scores(record))
        end
        push!(normalized_table, (
            record_identifier = FASTX.identifier(record),
            record_description = FASTX.description(record),
            record_sha256 = Mycelia.seq2sha256(record_sequence),
            record_quality = record_quality,
            record_alphabet = join(sort(collect(Set(uppercase(record_sequence))))),
            record_type = Mycelia.detect_alphabet(record_sequence),
            mean_record_quality = file_type == :fastq ? Statistics.mean(record_quality) : missing,
            median_record_quality = file_type == :fastq ? Statistics.median(record_quality) : missing,
            record_length = length(record_sequence),
            record_sequence = record_sequence,
        ))
    end
    current_columns = names(normalized_table)
    normalized_table[!, "fastx_path"] .= basename(fastx)
    normalized_table[!, "fastx_sha256"] .= Mycelia.metasha256(normalized_table[!, "record_sha256"])
    return normalized_table[!, ["fastx_path", "fastx_sha256", current_columns...]]
end