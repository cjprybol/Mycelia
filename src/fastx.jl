"""
    join_fastqs_with_uuid(
        fastq_files::Vector{String};
        fastq_out::String
        tsv_out::String
    )

Note: does not keep track of paired-end data - assumes single end reads

Designed primarily to allow joint mapping of many long-read samples

Given a collection of fastq files, creates:
- A gzipped TSV mapping original file and read_id to a new UUID per read
- A gzipped joint fastq file with the new UUID as read header

Returns: Tuple of output file paths (tsv_out, fastq_out)
"""
function join_fastqs_with_uuid(
    fastq_files::Vector{String};
    fastq_out::String = Mycelia.normalized_current_datetime() * ".joint_reads.fq.gz",
    tsv_out::String = replace(fastq_out, Mycelia.FASTQ_REGEX => ".tsv.gz")
)
    # Build mapping as a DataFrame in memory
    mapping = DataFrames.DataFrame(
        input_file = String[],
        original_read_id = String[],
        new_uuid = String[]
    )

    if (!isfile(fastq_out) || filesize(fastq_out) == 0) || (!isfile(tsv_out) || filesize(tsv_out) == 0)
        # Prepare gzipped FASTQ writer using CodecZlib
        gz_out = CodecZlib.GzipCompressorStream(open(fastq_out, "w"))
        writer = FASTX.FASTQ.Writer(gz_out)

        for fastq_file in fastq_files
            for record in Mycelia.open_fastx(fastq_file)
                original_id = String(FASTX.identifier(record))
                uuid = string(UUIDs.uuid4())
                # Save mapping to DataFrame
                DataFrames.push!(mapping, (fastq_file, original_id, uuid))
                # Write record with new UUID as id
                new_record = FASTX.FASTQ.Record(uuid, FASTX.sequence(record), FASTX.quality(record))
                FASTX.write(writer, new_record)
            end
        end

        FASTX.close(writer)
        CodecZlib.close(gz_out)

        # Write mapping table as gzipped TSV using CodecZlib
        tsv_io = CodecZlib.GzipCompressorStream(open(tsv_out, "w"))
        CSV.write(tsv_io, mapping; delim='\t')
        CodecZlib.close(tsv_io)
    else
        @warn "Output files already exist and are not empty: $fastq_out, $tsv_out"
    end
    return (;tsv_out, fastq_out)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function run_fastqc(;fastq, outdir = replace(fastq, Mycelia.FASTQ_REGEX => "_fastqc"))
    Mycelia.add_bioconda_env("fastqc")
    if !isdir(outdir)
        mkpath(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n fastqc fastqc --outdir $(outdir) $(fastq)`)
    else
        @warn "$outdir already exists"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function run_fastqc(;forward, reverse, outdir = Mycelia.find_matching_prefix(forward, reverse) * "_fastqc")
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
    return results_table
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
function qc_filter_short_reads_fastp(;
        forward_reads::String,
        reverse_reads::String,
        out_forward::String=replace(forward_reads, Mycelia.FASTQ_REGEX => ".fastp.1.fq.gz"),
        out_reverse::String=replace(reverse_reads, Mycelia.FASTQ_REGEX => ".fastp.2.fq.gz"),
        report_title::String="$(forward_reads) $(reverse_reads) fastp report",
        html::String=Mycelia.find_matching_prefix(out_forward, out_reverse) * ".fastp_report.html",
        json::String=Mycelia.find_matching_prefix(out_forward, out_reverse) * ".fastp_report.json"
    )
    # usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
    # options:
    #   # I/O options
    #   -i, --in1                          read1 input file name (string)
    #   -o, --out1                         read1 output file name (string [=])
    #   -I, --in2                          read2 input file name (string [=])
    #   -O, --out2                           read2 output file name (string [=])
    #       --unpaired1                      for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. (string [=])
    #       --unpaired2                      for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
    #       --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
    #       --overlapped_out                 for each read pair, output the overlapped region if it has no any mismatched base. (string [=])
    #   -m, --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
    #       --merged_out                     in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output (string [=])
    #       --include_unmerged               in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. Disabled by default.
    #   -6, --phred64                      indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
    #   -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
    #       --stdin                          input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
    #       --stdout                         output passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end input. Disabled by default.
    #       --interleaved_in                 indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
    #       --reads_to_process             specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
    #       --dont_overwrite               don't overwrite existing files. Overwritting is allowed by default.
    #       --fix_mgi_id                     the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.
    
    #   # adapter trimming options
    #   -A, --disable_adapter_trimming     adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
    #   -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
    #       --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=])
    #       --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
    #       --detect_adapter_for_pe          by default, the adapter sequence auto-detection is enabled for SE data only, turn on this option to enable it for PE data.
    
    #   # global trimming options
    #   -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
    #   -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
    #   -b, --max_len1                       if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
    #   -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
    #   -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
    #   -B, --max_len2                       if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
    
    #   # duplication evaluation and deduplication
    #   -D, --dedup                          enable deduplication to drop the duplicated reads/pairs
    #       --dup_calc_accuracy              accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
    #       --dont_eval_duplication          don't evaluate duplication rate to save time and use less memory.
    
    #   # polyG tail trimming, useful for NextSeq/NovaSeq data
    #   -g, --trim_poly_g                  force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
    #       --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
    #   -G, --disable_trim_poly_g          disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
    
    #   # polyX tail trimming
    #   -x, --trim_poly_x                    enable polyX trimming in 3' ends.
    #       --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
    
    #   # per read cutting by quality options
    #   -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
    #   -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
    #   -r, --cut_right                      move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
    #   -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
    #   -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
    #       --cut_front_window_size          the window size option of cut_front, default to cut_window_size if not specified (int [=4])
    #       --cut_front_mean_quality         the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
    #       --cut_tail_window_size           the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
    #       --cut_tail_mean_quality          the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
    #       --cut_right_window_size          the window size option of cut_right, default to cut_window_size if not specified (int [=4])
    #       --cut_right_mean_quality         the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
    
    #   # quality filtering options
    #   -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is specified, quality filtering is disabled
    #   -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
    #   -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
    #   -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
    #   -e, --average_qual                 if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
    
    
    #   # length filtering options
    #   -L, --disable_length_filtering     length filtering is enabled by default. If this option is specified, length filtering is disabled
    #   -l, --length_required              reads shorter than length_required will be discarded, default is 15. (int [=15])
    #       --length_limit                 reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
    
    #   # low complexity filtering
    #   -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
    #   -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
    
    #   # filter reads with unwanted indexes (to remove possible contamination)
    #       --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
    #       --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
    #       --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
    
    #   # base correction by overlap analysis options
    #   -c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled
    #       --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
    #       --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
    #       --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
    
    #   # UMI processing
    #   -U, --umi                          enable unique molecular identifier (UMI) preprocessing
    #       --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
    #       --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])
    #       --umi_prefix                   if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
    #       --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])
    
    #   # overrepresented sequence analysis
    #   -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
    #   -P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])
    
    #   # reporting options
    #   -j, --json                         the json format report file name (string [=fastp.json])
    #   -h, --html                         the html format report file name (string [=fastp.html])
    #   -R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])
    
    #   # threading options
    #   -w, --thread                       worker thread number, default is 3 (int [=3])
    
    #   # output splitting options
    #   -s, --split                        split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
    #   -S, --split_by_lines               split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
    #   -d, --split_prefix_digits          the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
    
    #   # help
    #   -?, --help                         print this message
    if !isfile(out_forward) || !isfile(out_reverse) || !isfile(json) || !isfile(html)
        Mycelia.add_bioconda_env("fastp")
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastp fastp
                --in1 $(forward_reads)
                --in2 $(reverse_reads)
                --out1 $(out_forward)
                --out2 $(out_reverse)
                --json $(json)
                --html $(html)
                --report_title $(report_title)
                --dedup`
        run(cmd)
    else
        @show isfile(out_forward)
        @show isfile(out_reverse)
        @show isfile(json)
        @show isfile(html)
    end
    return (;out_forward, out_reverse, json, html)
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
    fastx2normalized_table(fastx::AbstractString) -> DataFrames.DataFrame

$(DocStringExtensions.TYPEDSIGNATURES)

Read a FASTA or FASTQ file and convert its records into a normalized `DataFrames.DataFrame` where each row represents a sequence record and columns provide standardized metadata and sequence statistics.

# Arguments

- `fastx::AbstractString`: Path to a FASTA or FASTQ file. The file must exist and be non-empty. The file type is inferred from the filename using `Mycelia.FASTA_REGEX` and `Mycelia.FASTQ_REGEX`.

# Returns

- `DataFrames.DataFrame`: A data frame where each row contains information for a record from the input file, and columns include:
    - `fastx_path`: Basename of the input file.
    - `fastx_sha256`: Aggregated SHA256 hash of all record SHA256s in the file.
    - `record_identifier`: Identifier from the record header.
    - `record_description`: Description from the record header.
    - `record_sha256`: SHA256 hash of the record sequence.
    - `record_quality`: Vector of quality scores (`Vector{Float64}`) for FASTQ, or `missing` for FASTA.
    - `record_alphabet`: Sorted, joined string of unique, uppercase characters in the record sequence.
    - `record_type`: Alphabet type detected by `Mycelia.detect_alphabet` (e.g., `:DNA`, `:RNA`, etc.).
    - `mean_record_quality`: Mean quality score (for FASTQ), or `missing` (for FASTA).
    - `median_record_quality`: Median quality score (for FASTQ), or `missing` (for FASTA).
    - `record_length`: Length of the sequence.
    - `record_sequence`: The sequence string itself.

# Notes

- The function asserts that the file exists and is not empty.
- File type is determined by regex matching on the filename.
- For FASTA files, quality-related columns are set to `missing`.
- For FASTQ files, quality scores are extracted and statistics are computed.
- Record SHA256 hashes are aggregated to compute a file-level SHA256 via `Mycelia.metasha256`.
- Requires the following namespaces: `DataFrames`, `Statistics`, `Mycelia`, `FASTX`, and `Base.basename`.
- The function returns the columns in the order: `fastx_path`, `fastx_sha256`, followed by all other record columns.

# Example

```julia
import DataFrames
import Mycelia
import FASTX

table = fastx2normalized_table("example.fasta")
DataFrames.first(table, 3)
```
"""
function fastx2normalized_table(fastx)
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
        # @info "Processing FASTA file"
        file_type = :fasta
    elseif occursin(Mycelia.FASTQ_REGEX, fastx)
        # @info "Processing FASTQ file"
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

"""
    write_fastq(;records, filename, gzip=false)

$(DocStringExtensions.TYPEDSIGNATURES)

Write FASTQ records to file using FASTX.jl.
Validates extension: .fastq, .fq, .fastq.gz, or .fq.gz.
If `gzip` is true or filename endswith .gz, output is gzipped.
`records` must be an iterable of FASTX.FASTQ.Record.
"""
function write_fastq(;records, filename, gzip=false)
    function is_valid_fastq_ext(fname)
        any(endswith.(fname, [".fastq", ".fq", ".fastq.gz", ".fq.gz"]))
    end

    if !is_valid_fastq_ext(filename)
        error("File extension must be .fastq, .fq, .fastq.gz, or .fq.gz")
    end

    gzip_out = gzip || endswith(filename, ".gz")
    if gzip_out
        io = CodecZlib.GzipCompressorStream(open(filename, "w"))
    else
        io = open(filename, "w")
    end

    try
        writer = FASTX.FASTQ.Writer(io)
        for record in records
            FASTX.write(writer, record)
        end
        FASTX.close(writer)
    finally
        close(io)
    end
    return filename
end

"""
    write_fastas_from_normalized_fastx_tables(
        table_paths::Vector{String};
        output_dir::String = pwd(),
        show_progress::Bool = true,
        overwrite::Bool = false,
        error_handler = (e, table_path)->display((e, table_path))
    ) -> NamedTuple

Given a vector of normalized fastx table paths, writes out gzipped FASTA files in parallel.
Each table must have columns: "fastx_sha256", "record_sha256", "record_sequence".
Automatically decompresses input files if they end with ".gz".
Returns a summary NamedTuple with successes, failures, failed tables, and output files.

# Keyword Arguments
- `output_dir`: Directory to write .fna.gz files to.
- `show_progress`: Show a progress bar (default: true).
- `overwrite`: Overwrite existing files (default: false).
- `error_handler`: Function called with (exception, table_path) on error.

"""
function write_fastas_from_normalized_fastx_tables(
    table_paths::Vector{String};
    output_dir::String = pwd(),
    show_progress::Bool = true,
    overwrite::Bool = false,
    error_handler = (e, table_path)->display((e, table_path))
)
    n = length(table_paths)
    successes = Base.Threads.Atomic{Int}(0)
    failures = Base.Threads.Atomic{Int}(0)
    failed_tables = String[]
    failed_tables_lock = Base.Threads.SpinLock()
    output_files = String[]
    output_files_lock = Base.Threads.SpinLock()

    progress = show_progress ? ProgressMeter.Progress(n, 1) : nothing

    Base.Threads.@threads for i in 1:n
        table_path = table_paths[i]
        try
            # Determine if the table should be decompressed based on extension
            decompress = endswith(table_path, ".gz")
            io = open(table_path)
            stream = decompress ? CodecZlib.GzipDecompressorStream(io) : io
            loaded_table = CSV.read(stream, DataFrames.DataFrame, delim='\t')
            close(io)
            # Remove .gz if present in the fastx_sha256 (basename)
            fastx_basename = loaded_table[1, "fastx_sha256"]
            if occursin(r"\.gz$", fastx_basename)
                fastx_basename = replace(fastx_basename, r"\.gz$" => "")
            end
            fastx_filename = fastx_basename * ".fna.gz"
            fastx_file = joinpath(output_dir, fastx_filename)
            # Write if missing or overwrite specified
            write_file = overwrite || !isfile(fastx_file)
            if write_file
                fastx_records = [
                    FASTX.FASTA.Record(row["record_sha256"], row["record_sequence"])
                    for row in DataFrames.eachrow(loaded_table)
                ]
                Mycelia.write_fasta(outfile = fastx_file, records = fastx_records)
            end
            # Record output file path
            Base.Threads.lock(output_files_lock) do
                push!(output_files, fastx_file)
            end
            Base.Threads.atomic_add!(successes, 1)
        catch e
            error_handler(e, table_path)
            Base.Threads.atomic_add!(failures, 1)
            Base.Threads.lock(failed_tables_lock) do
                push!(failed_tables, table_path)
            end
        end
        if show_progress
            ProgressMeter.next!(progress)
        end
    end

    if show_progress
        ProgressMeter.finish!(progress)
    end

    return (
        total = n,
        succeeded = successes[],
        failed = failures[],
        failed_tables = failed_tables,
        output_files = output_files
    )
end

# convert all genomes into normalized tables
function fastxs2normalized_tables(;fastxs, outdir, force=false)
    mkpath(outdir)
    n = length(fastxs)
    normalized_table_paths = Vector{Union{String, Nothing}}(undef, n)
    errors = Vector{Union{Nothing, Tuple{String, Exception}}}(undef, n)

    prog = ProgressMeter.Progress(n, desc = "Processing files")

    Threads.@threads for i in 1:n
        fastx = fastxs[i]
        outfile_base = basename(fastx)
        if occursin(r"\.gz$", outfile_base)
            outfile_base = replace(outfile_base, r"\.gz$" => "")
        end
        outfile = joinpath(outdir, outfile_base * ".tsv.gz")
        try
            if !isfile(outfile) || (filesize(outfile) == 0) || force
                normalized_table = Mycelia.fastx2normalized_table(fastx)
                open(outfile, "w") do file
                    io = CodecZlib.GzipCompressorStream(file)
                    CSV.write(io, normalized_table; delim='\t', bufsize=64*1024*1024)
                    close(io)
                end
            end
            normalized_table_paths[i] = outfile
            errors[i] = nothing
        catch e
            normalized_table_paths[i] = nothing
            errors[i] = (fastx, e)
            @warn "Failed to process $fastx: $e"
        end
        ProgressMeter.next!(prog)
    end

    ProgressMeter.finish!(prog)
    return (;normalized_table_paths, errors)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Detect sequence type from input and suggest appropriate file extension.

# Arguments
- `record`: A FASTA/FASTQ record
- `sequence`: A string or BioSequence containing sequence data

# Returns
- `String`: Suggested file extension:
  - ".fna" for DNA
  - ".frn" for RNA
  - ".faa" for protein
  - ".fa" for unrecognized sequences
"""
function detect_sequence_extension(record::Union{FASTX.FASTA.Record, FASTX.FASTQ.Record})
    return detect_sequence_extension(FASTX.sequence(record))
end
function detect_sequence_extension(sequence::AbstractString)
    sequence_type = detect_alphabet(sequence)
    return _detect_sequence_extension(sequence_type::Symbol)
end
function detect_sequence_extension(sequence::BioSequences.LongSequence)
    sequence_type = detect_alphabet(sequence)
    return _detect_sequence_extension(sequence_type::Symbol)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Internal helper function to convert sequence type to file extension.

Arguments
- sequence_type: Symbol representing sequence type (:DNA, :RNA, or :AA)

Returns
- String: Appropriate file extensions
"""
function _detect_sequence_extension(sequence_type::Symbol)
    @assert sequence_type in [:DNA, :RNA, :AA]
    if sequence_type == :DNA
        return ".fna"
    elseif sequence_type == :RNA
        return ".frn"
    elseif sequence_type == :AA
        return ".faa"
    else
        @warn "unrecognized sequence type: $(seq_type)"
        return ".fa"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compare two FASTA files to determine if they contain the same set of sequences,
regardless of sequence order.

# Arguments
- `fasta_1::String`: Path to first FASTA file
- `fasta_2::String`: Path to second FASTA file

# Returns
- `Bool`: `true` if both files contain exactly the same sequences, `false` otherwise

# Details
Performs a set-based comparison of DNA sequences by hashing each sequence.
Sequence order differences between files do not affect the result.
"""
function equivalent_fasta_sequences(fasta_1, fasta_2)
    fasta_1_hashes = Set(hash(BioSequences.LongDNA{2}(FASTX.sequence(record))) for record in Mycelia.open_fastx(fasta_1))
    fasta_2_hashes = Set(hash(BioSequences.LongDNA{2}(FASTX.sequence(record))) for record in Mycelia.open_fastx(fasta_2))
    @show setdiff(fasta_1_hashes, fasta_2_hashes)
    @show setdiff(fasta_2_hashes, fasta_1_hashes)
    return fasta_1_hashes == fasta_2_hashes
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Join fasta files while adding origin prefixes to the identifiers.

Does not guarantee uniqueness but will warn if conflicts arise
"""
function merge_fasta_files(;fasta_files, fasta_file)
    @info "merging $(length(fasta_files)) files..."
    identifiers = Set{String}()
    open(fasta_file, "w") do io
        fastx_io = FASTX.FASTA.Writer(io)
        ProgressMeter.@showprogress for f in fasta_files
            f_id = replace(basename(f), Mycelia.FASTA_REGEX => "")
            for record in Mycelia.open_fastx(f)
                new_record_id = f_id * "__" * FASTX.identifier(record)
                if new_record_id in identifiers
                    @warn "new identifier $(new_record_id) already in identifiers!!!"
                end
                push!(identifiers, new_record_id)
                new_record = FASTX.FASTA.Record(new_record_id, FASTX.sequence(record))
                write(fastx_io, new_record)
            end
        end
    end
    @info "$(length(identifiers)) records merged..."
    return fasta_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Remove duplicate sequences from a FASTA file while preserving headers.

# Arguments
- `in_fasta`: Path to input FASTA file
- `out_fasta`: Path where deduplicated FASTA will be written

# Returns
Path to the output FASTA file (same as `out_fasta` parameter)

# Details
- Sequences are considered identical if they match exactly (case-sensitive)
- For duplicate sequences, keeps the first header encountered
- Input sequences are sorted by identifier before deduplication
- Preserves the original sequence formatting
"""
function deduplicate_fasta_file(in_fasta, out_fasta)
    fasta_df = fasta_to_table(collect(open_fastx(in_fasta)))
    sort!(fasta_df, "identifier")
    unique_sequences = DataFrames.combine(DataFrames.groupby(fasta_df, "sequence"), first)
    fasta = fasta_table_to_fasta(unique_sequences)
    open(out_fasta, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for record in fasta
            write(writer, record)
        end
        close(writer)
    end
    return out_fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Translates nucleic acid sequences from a FASTA file into amino acid sequences.

# Arguments
- `fasta_nucleic_acid_file::String`: Path to input FASTA file containing nucleic acid sequences
- `fasta_amino_acid_file::String`: Path where the translated amino acid sequences will be written

# Returns
- `String`: Path to the output FASTA file containing translated amino acid sequences
"""
function translate_nucleic_acid_fasta(fasta_nucleic_acid_file, fasta_amino_acid_file)
    open(fasta_amino_acid_file, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for record in FASTX.FASTA.Reader(open(fasta_nucleic_acid_file))
            try
                raw_seq = FASTX.sequence(record)
                pruned_seq_length = Int(floor(length(raw_seq)/3)) * 3
                truncated_seq = raw_seq[1:pruned_seq_length]
                amino_acid_seq = BioSequences.translate(truncated_seq)
                amino_acid_record = FASTX.FASTA.Record(FASTX.identifier(record), FASTX.description(record), amino_acid_seq)
                write(writer, amino_acid_record)
            catch
                @warn "unable to translate record", record
            end
        end
        close(writer)
    end
    return fasta_amino_acid_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a FASTA file/record iterator to a DataFrame.

# Arguments
- `fasta`: FASTA record iterator from FASTX.jl

# Returns
- `DataFrame` with columns:
  - `identifier`: Sequence identifiers
  - `description`: Full sequence descriptions 
  - `sequence`: Biological sequences as strings
"""
function fasta_to_table(fasta)
    collected_fasta = collect(fasta)
    fasta_df = DataFrames.DataFrame(
        identifier = FASTX.identifier.(collected_fasta),
        description = FASTX.description.(collected_fasta),
        sequence = FASTX.sequence.(collected_fasta)
    )
    return fasta_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a DataFrame containing FASTA sequence information into a vector of FASTA records.

# Arguments
- `fasta_df::DataFrame`: DataFrame with columns "identifier", "description", and "sequence"

# Returns
- `Vector{FASTX.FASTA.Record}`: Vector of FASTA records
"""
function fasta_table_to_fasta(fasta_df)
    records = Vector{FASTX.FASTA.Record}(undef, DataFrames.nrow(fasta_df))
    for (i, row) in enumerate(DataFrames.eachrow(fasta_df))
        record = FASTX.FASTA.Record(row["identifier"], row["description"], row["sequence"])
        records[i] = record
    end
    return records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

This turns a 4-line FASTQ entry into a single tab separated line,
adds a column with the length of each read, passes it to Unix sort,
removes the length column, and converts it back into a FASTQ file.

sorts longest to shortest!!

http://thegenomefactory.blogspot.com/2012/11/sorting-fastq-files-by-sequence-length.html
"""
function sort_fastq(input_fastq, output_fastq="")
    
    if endswith(input_fastq, ".gz")
        p = pipeline(
                `gzip -dc $input_fastq`,
                `paste - - - -`,
                `perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'`,
                `sort -nr`,
                `cut -f2-`,
                `tr "\t" "\n"`,
                `gzip`
                )
    else
        p = pipeline(
                `cat $input_fastq`,
                `paste - - - -`,
                `perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'`,
                `sort -nr`,
                `cut -f2-`,
                `tr "\t" "\n"`
                )
    end
    run(pipeline(p, output_fastq))
    return output_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Counts the total number of records in a FASTA/FASTQ file.

# Arguments
- `fastx`: Path to a FASTA or FASTQ file (can be gzipped)

# Returns
- Number of records (sequences) in the file
"""
function count_records(fastx)
    n_records = 0
    for record in open_fastx(fastx)
        n_records += 1
    end
    return n_records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate sequence lengths for reads in a FASTQ file.

# Arguments
- `fastq_file::String`: Path to input FASTQ file
- `total_reads::Integer=Inf`: Number of reads to process (defaults to all reads)

# Returns
- `Vector{Int}`: Array containing the length of each sequence read
"""
function determine_read_lengths(fastq_file; total_reads = Inf)
    if total_reads == Inf
        total_reads = count_records(fastq_file)
    end
    read_lengths = zeros(Int, total_reads)
    @info "determining read lengths"
    p = ProgressMeter.Progress(total_reads, 1)
    for (i, record) in enumerate(open_fastx(fastq_file))
#         push!(read_lengths, length(FASTX.sequence(record)))
        read_lengths[i] = length(FASTX.sequence(record))
        ProgressMeter.next!(p)
    end
    return read_lengths
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Phred quality score (Q-value) to a probability of error.

# Arguments
- `q_value`: Phred quality score, typically ranging from 0 to 40

# Returns
- Error probability in range [0,1], where 0 indicates highest confidence

A Q-value of 10 corresponds to an error rate of 0.1 (10%), while a Q-value of 
30 corresponds to an error rate of 0.001 (0.1%).
"""
function q_value_to_error_rate(q_value)
    error_rate = 10^(q_value/(-10))
    return error_rate
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a sequencing error probability to a Phred quality score (Q-value).

The calculation uses the standard Phred formula: Q = -10 * log₁₀(error_rate)

# Arguments
- `error_rate::Float64`: Probability of error (between 0 and 1)

# Returns
- `q_value::Float64`: Phred quality score
"""
function error_rate_to_q_value(error_rate)
    q_value = -10 * log10(error_rate)
    return q_value
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the total size (in bases) of all sequences in a FASTA file.

# Arguments
- `fasta_file::AbstractString`: Path to the FASTA file

# Returns
- `Int`: Sum of lengths of all sequences in the FASTA file
"""
function fasta_genome_size(fasta_file)
    return reduce(sum, map(record -> length(FASTX.sequence(record)), Mycelia.open_fastx(fasta_file)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Subsample reads from a FASTQ file using seqkit.

# Arguments
- `in_fastq::String`: Path to input FASTQ file
- `out_fastq::String=""`: Path to output FASTQ file. If empty, auto-generated based on input filename
- `n_reads::Union{Missing,Int}=missing`: Number of reads to sample
- `proportion_reads::Union{Missing,Float64}=missing`: Proportion of reads to sample (0.0-1.0)

# Returns
- `String`: Path to the output FASTQ file
"""
function subsample_reads_seqkit(;in_fastq::String, out_fastq::String="", n_reads::Union{Missing, Int}=missing, proportion_reads::Union{Missing, Float64}=missing)
    Mycelia.add_bioconda_env("seqkit")
    if ismissing(n_reads) && ismissing(proportion_reads)
        error("please specify the number or proportion of reads")
    elseif !ismissing(n_reads)
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit sample --two-pass --number $(n_reads) $(in_fastq)`, `gzip`)
        if isempty(out_fastq)
            out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqkit.N$(n_reads).fq.gz")
        end
    elseif !ismissing(proportion_reads)
        p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit sample --proportion $(proportion_reads) $(in_fastq)`, `gzip`)
        if isempty(out_fastq)
            out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqkit.P$(proportion_reads).fq.gz")
        end
    end
    @assert !isempty(out_fastq)
    if !isfile(out_fastq)
        run(pipeline(p, out_fastq))
    else
        @info "$(out_fastq) already present"
    end
    return out_fastq
end

# function subsample_reads_seqtk(;in_fastq::String, out_fastq="", n_reads::Union{Missing, Int}=missing, proportion_reads::Union{Missing, Float64}=missing)
#     Mycelia.add_bioconda_env("seqtk")
#     if ismissing(n_reads) && ismissing(proportion_reads)
#         error("please specify the number or proportion of reads")
#     elseif !ismissing(n_reads)
#         p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqtk seqtk sample -2 $(n_reads) $(in_fastq)`, `gzip`)
#         if isempty(out_fastq)
#             out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqtk.N$(n_reads).fq.gz")
#         end
#     elseif !ismissing(proportion_reads)
#         p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqtk seqtk sample -2 $(proportion_reads) $(in_fastq)`, `gzip`)
#         if isempty(out_fastq)
#             out_fastq = replace(in_fastq, Mycelia.FASTQ_REGEX => ".seqtk.P$(proportion_reads).fq.gz")
#         end
#     end
#     @assert !isempty(out_fastq)
#     run(pipeline(p, out_fastq))
#     return out_fastq
# end

# subsample_reads_seqtk(in_fastq = fastq, n_reads=10)

# function filter_short_reads()
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates an index file (.fai) for a FASTA reference sequence using samtools.

The FASTA index allows efficient random access to the reference sequence. This is 
required by many bioinformatics tools that need to quickly fetch subsequences 
from the reference.

# Arguments
- `fasta`: Path to the input FASTA file

# Side Effects
- Creates a `{fasta}.fai` index file in the same directory as input
- Installs samtools via conda if not already present
"""
function samtools_index_fasta(;fasta)
    Mycelia.add_bioconda_env("samtools")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools faidx $(fasta)`)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Construct a FASTX FASTQ record from its components.

# Arguments
- `identifier::String`: The sequence identifier without the '@' prefix
- `sequence::String`: The nucleotide sequence
- `quality_scores::Vector{Int}`: Quality scores (0-93) as raw integers

# Returns
- `FASTX.FASTQRecord`: A parsed FASTQ record

# Notes
- Quality scores are automatically capped at 93 to ensure FASTQ compatibility
- Quality scores are converted to ASCII by adding 33 (Phred+33 encoding)
- The record is constructed in standard FASTQ format with four lines:
  1. Header line (@ + identifier)
  2. Sequence
  3. Plus line
  4. Quality scores (ASCII encoded)
"""
function fastq_record(;identifier, sequence, quality_scores)
    # Fastx wont parse anything higher than 93
    quality_scores = min.(quality_scores, 93)
    record_string = join(["@" * identifier, sequence, "+", join([Char(x+33) for x in quality_scores])], "\n")
    return FASTX.parse(FASTX.FASTQRecord, record_string)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate detailed mapping statistics for each reference sequence/contig in a XAM (SAM/BAM/CRAM) file.

# Arguments
- `xam`: Path to XAM file or XAM object

# Returns
A DataFrame with per-contig statistics including:
- `n_aligned_reads`: Number of aligned reads
- `total_aligned_bases`: Sum of alignment lengths
- `total_alignment_score`: Sum of alignment scores
- Mapping quality statistics (mean, std, median)
- Alignment length statistics (mean, std, median)
- Alignment score statistics (mean, std, median)
- Percent mismatches statistics (mean, std, median)

Note: Only primary alignments (isprimary=true) and mapped reads (ismapped=true) are considered.
"""
function fastx_to_contig_lengths(fastx)
    OrderedCollections.OrderedDict(String(FASTX.identifier(record)) => length(FASTX.sequence(record)) for record in Mycelia.open_fastx(fastx))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Open and return a reader for FASTA or FASTQ format files.

# Arguments
- `path::AbstractString`: Path to input file. Can be:
    - Local file path
    - HTTP/FTP URL
    - Gzip compressed (.gz extension)

# Supported formats
- FASTA (.fasta, .fna, .faa, .fa)
- FASTQ (.fastq, .fq)

# Returns
- `FASTX.FASTA.Reader` for FASTA files
- `FASTX.FASTQ.Reader` for FASTQ files
"""
function open_fastx(path::AbstractString)
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
        path_base = replace(path_base, ".gz" => "")
    end
    if occursin(r"\.(fasta|fna|faa|fa|frn)$", path_base)
        fastx_io = FASTX.FASTA.Reader(io)
    elseif occursin(r"\.(fastq|fq)$", path_base)
        fastx_io = FASTX.FASTQ.Reader(io)
    else
        error("attempting to open a FASTX file with an unsupported extension: $(path_base)")
    end
    return fastx_io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Writes FASTA records to a file, optionally gzipped.

# Arguments
- `outfile::AbstractString`: Path to the output FASTA file.  Will append ".gz" if `gzip` is true and ".gz" isn't already the extension.
- `records::Vector{FASTX.FASTA.Record}`: A vector of FASTA records.
- `gzip::Bool`: Optionally force compression of the output with gzip. By default will use the file name to infer.

# Returns
- `outfile::String`: The path to the output FASTA file (including ".gz" if applicable).
"""
function write_fasta(;outfile::AbstractString=tempname()*".fna", records::Vector{FASTX.FASTA.Record}, gzip::Bool=false)
    # Determine if gzip compression should be used based on both the filename and the gzip argument
    gzip = occursin(r"\.gz$", outfile) || gzip

    # Append ".gz" to the filename if gzip is true and the extension isn't already present
    outfile = gzip && !occursin(r"\.gz$", outfile) ? outfile * ".gz" : outfile  # More concise way to handle filename modification

    # Use open with do block for automatic resource management (closing the file)
    open(outfile, "w") do io
        if gzip
            io = CodecZlib.GzipCompressorStream(io)  # Wrap the io stream for gzip compression
        end

        FASTX.FASTA.Writer(io) do fastx_io  # Use do block for automatic closing of the FASTA writer
            for record in records
                write(fastx_io, record)
            end
        end # fastx_io automatically closed here
    end # io automatically closed here

    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GenBank format file to FASTA format using EMBOSS seqret.

# Arguments
- `genbank`: Path to input GenBank format file
- `fasta`: Optional output FASTA file path (defaults to input path with .fna extension)
- `force`: If true, overwrites existing output file (defaults to false)

# Returns
Path to the output FASTA file

# Notes
- Requires EMBOSS suite (installed automatically via Conda)
- Will not regenerate output if it already exists unless force=true
"""
function genbank_to_fasta(;genbank, fasta=genbank * ".fna", force=false)
    add_bioconda_env("emboss")
    if !isfile(fasta) || force
        run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret $(genbank) fasta:$(fasta)`)
    end
    return fasta
end