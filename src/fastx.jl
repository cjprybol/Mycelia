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