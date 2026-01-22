# function normalize_vcf(;reference_fasta, vcf, normalized_vcf=)
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize a VCF file using bcftools norm, with automated handling of compression and indexing.

# Arguments
- `reference_fasta::String`: Path to the reference FASTA file used for normalization
- `vcf_file::String`: Path to input VCF file (can be gzipped or uncompressed)

# Returns
- `String`: Path to the normalized, sorted, and compressed output VCF file (*.sorted.normalized.vcf.gz)

# Notes
- Requires bioconda packages: htslib, tabix, bcftools
- Creates intermediate files with extensions .tbi for indices
- Skips processing if output file already exists
- Performs left-alignment and normalization of variants
"""
function normalize_vcf(;reference_fasta, vcf_file)
    add_bioconda_env("htslib")
    add_bioconda_env("tabix")
    add_bioconda_env("bcftools")
    normalized_vcf=replace(vcf_file, Mycelia.VCF_REGEX => ".sorted.normalized.vcf")
    out_vcf = normalized_vcf * ".gz"
    if !isfile(out_vcf)
        if occursin(r"\.gz$", vcf_file)
            gzipped_vcf = vcf_file
        else
            gzipped_vcf = "$(vcf_file).gz"
            run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n htslib bgzip -c $(vcf_file)`, gzipped_vcf))
        end
        tabix_index = "$(gzipped_vcf).tbi"
        if !isfile(tabix_index)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n tabix tabix -f -p vcf $(gzipped_vcf)`)
        end
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bcftools bcftools norm -cs --fasta-ref $(reference_fasta) $(gzipped_vcf)`, normalized_vcf))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n htslib bgzip $(normalized_vcf)`)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n tabix tabix -p vcf $(out_vcf)`)
    end
    return out_vcf
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Apply variants from a VCF file to a reference FASTA sequence.

# Arguments
- `in_fasta`: Path to input reference FASTA file
- `vcf_file`: Path to input VCF file containing variants
- `out_fasta`: Optional output path for modified FASTA. Defaults to replacing '.vcf' with '.normalized.vcf.fna'

# Details
1. Normalizes indels in the VCF using bcftools norm
2. Applies variants to the reference sequence using bcftools consensus
3. Handles temporary files and compression with bgzip/tabix

# Requirements
Requires bioconda packages: htslib, tabix, bcftools

# Returns
Path to the output FASTA file containing the modified sequence
"""
function update_fasta_with_vcf(;in_fasta, vcf_file, out_fasta=replace(vcf_file, ".vcf" => ".normalized.vcf.fna"))
    add_bioconda_env("htslib")
    add_bioconda_env("tabix")
    add_bioconda_env("bcftools")
    isfile("$(vcf_file).gz") && rm("$(vcf_file).gz")
    isfile("$(vcf_file).gz.tbi") && rm("$(vcf_file).gz.tbi")
    normalized_vcf_file = replace(vcf_file, ".vcf" => ".normalized.vcf")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n htslib bgzip $(vcf_file)`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n tabix tabix -f -p vcf $(vcf_file).gz`)
    run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bcftools bcftools norm -cs --fasta-ref $(in_fasta) $(vcf_file).gz`, normalized_vcf_file))
    rm("$(vcf_file).gz")
    rm("$(vcf_file).gz.tbi")
    isfile("$(normalized_vcf_file).gz") && rm("$(normalized_vcf_file).gz")
    isfile("$(normalized_vcf_file).gz.tbi") && rm("$(normalized_vcf_file).gz.tbi")
    isfile("$(normalized_vcf_file).fna") && rm("$(normalized_vcf_file).fna")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n htslib bgzip $(normalized_vcf_file)`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n tabix tabix -p vcf $(normalized_vcf_file).gz`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bcftools bcftools consensus -f $(in_fasta) $(normalized_vcf_file).gz -o $(out_fasta)`)
    return out_fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write variant data to a VCF v4.3 format file.

# Arguments
- `vcf_file::String`: Output path for the VCF file
- `vcf_table::DataFrame`: Table containing variant data with standard VCF columns
- `fasta_file::String`: Path to the reference genome FASTA file

# Details
Automatically filters out equivalent variants where REF == ALT.
Includes standard VCF headers for substitutions, insertions, deletions, and inversions.
Adds GT (Genotype) and GQ (Genotype Quality) format fields.
"""
function write_vcf_table(;vcf_file, vcf_table, fasta_file)
    true_variant = vcf_table[!, "REF"] .!= vcf_table[!, "ALT"]
    if !all(true_variant)
        @warn "filtering equivalent variants"
        vcf_table = vcf_table[true_variant, :]
    end
    open(vcf_file, "w") do io
        VCF_HEADER = 
        """
        ##fileformat=VCFv4.3
        ##fileDate=$(Dates.today())
        ##source=simulated-variants
        ##reference=$(fasta_file)
        ##FILTER=<ID=substitution,Description="substitution variant">
        ##FILTER=<ID=insertion,Description="insertion variant">
        ##FILTER=<ID=deletion,Description="deletion variant">
        ##FILTER=<ID=inversion,Description="inversion variant">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        """
        print(io, VCF_HEADER)
        println(io, join(names(vcf_table), '\t'))
        for row in DataFrames.eachrow(vcf_table)
            println(io, join([row[col] for col in names(vcf_table)], '\t'))
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse RTG evaluation output from a gzipped tab-separated file.

# Arguments
- `f`: Path to a gzipped TSV file containing RTG evaluation output

# Format
Expected file format:
- Header line starting with '#' and tab-separated column names
- Data rows in tab-separated format
- Empty files return a DataFrame with empty columns matching header

# Returns
A DataFrame where:
- Column names are taken from the header line (stripped of '#')
- Data is parsed as Float64 values
- Empty files result in empty columns preserving header structure
"""
function parse_rtg_eval_output(f)
    # import CodecZlib
    flines = readlines(CodecZlib.GzipDecompressorStream(open(f)))
    header_line = last(filter(fline -> occursin(r"^#", fline), flines))
    header = lstrip.(split(header_line, "\t"), '#')
    data_lines = filter(fline -> !occursin(r"^#", fline), flines)
    if isempty(data_lines)
        data = [Float64[] for i in 1:length(header)]
    else
        data, h = uCSV.read(IOBuffer(join(data_lines, '\n')), delim='\t')
    end
    # data = [[parse(Float64, x)] for x in split(last(flines), '\t')]
    # @show data, header
    DataFrames.DataFrame(data, header)
end

# vcf-eval Integration Functions

"""
    run_vcfeval(baseline_vcf::String, calls_vcf::String, reference_fasta::String, output_dir::String;
                threads::Int=1, memory_gb::Int=8, squash_ploidy::Bool=true, all_records::Bool=true,
                score_field::String="QUAL")

Run RTG vcfeval to evaluate variant calling accuracy against a baseline truth set.

Uses the RTG Tools vcfeval utility to compare called variants against a baseline
truth set, generating precision/recall metrics and ROC curves.

# Arguments
- `baseline_vcf`: Path to baseline/truth VCF file
- `calls_vcf`: Path to variant calls VCF file to evaluate
- `reference_fasta`: Path to reference genome FASTA file
- `output_dir`: Directory for vcfeval output files
- `threads`: Number of threads for parallel processing (default: 1)
- `memory_gb`: Memory allocation in GB (default: 8)
- `squash_ploidy`: Squash samples to be haploid (default: true)
- `all_records`: Include all records in evaluation (default: true)
- `score_field`: VCF field to use for scoring (default: "QUAL")

# Returns
- Dictionary with paths to key output files (summary, ROC curves, etc.)

# Example
```julia
results = Mycelia.run_vcfeval(
    "truth.vcf.gz", "calls.vcf.gz", "reference.fasta", "eval_output",
    memory_gb=16, threads=4
)
```
"""
function run_vcfeval(baseline_vcf::String, calls_vcf::String, reference_fasta::String, output_dir::String;
                    threads::Int=1, memory_gb::Int=8, squash_ploidy::Bool=true, all_records::Bool=true,
                    score_field::String="QUAL")
    
    # Validate input files
    for file in [baseline_vcf, calls_vcf, reference_fasta]
        if !isfile(file)
            error("Input file does not exist: $(file)")
        end
    end
    
    # Ensure RTG Tools environment exists
    add_bioconda_env("rtg-tools")
    
    # Prepare reference template if not already formatted
    template_dir = reference_fasta * "_RTG"
    if !isdir(template_dir)
        println("Formatting reference genome for RTG Tools...")
        cmd = `$(CONDA_RUNNER) run --live-stream -n rtg-tools rtg format -o $(template_dir) $(reference_fasta)`
        try
            run(cmd)
        catch e
            error("RTG format failed: $(e)")
        end
    end
    
    # Prepare VCF files - ensure they are sorted, compressed, and indexed
    processed_baseline = prepare_vcf_for_rtg(baseline_vcf)
    processed_calls = prepare_vcf_for_rtg(calls_vcf)
    
    # Create output directory
    mkpath(output_dir)
    
    # Build vcfeval command
    vcfeval_args = [
        "RTG_MEM=$(memory_gb)G", "vcfeval",
        "--template", template_dir,
        "--baseline", processed_baseline,
        "--calls", processed_calls,
        "--output", output_dir,
        "--threads", string(threads),
        "--vcf-score-field", score_field
    ]
    
    # Add optional flags
    if all_records
        push!(vcfeval_args, "--all-records")
    end
    if squash_ploidy
        push!(vcfeval_args, "--squash-ploidy")
    end
    
    println("Running RTG vcfeval...")
    println("  Baseline: $(basename(baseline_vcf))")
    println("  Calls: $(basename(calls_vcf))")
    println("  Output: $(output_dir)")
    
    cmd = `$(CONDA_RUNNER) run --live-stream -n rtg-tools rtg $(vcfeval_args)`
    
    try
        run(cmd)
    catch e
        error("RTG vcfeval failed: $(e)")
    end
    
    # Collect output files
    output_files = Dict{String, String}()
    
    # Standard vcfeval outputs
    standard_outputs = [
        "summary.txt",
        "snp_roc.tsv.gz", 
        "non_snp_roc.tsv.gz",
        "weighted_roc.tsv.gz",
        "tp.vcf.gz",      # True positives
        "fp.vcf.gz",      # False positives  
        "fn.vcf.gz"       # False negatives
    ]
    
    for output in standard_outputs
        file_path = joinpath(output_dir, output)
        if isfile(file_path)
            key = replace(output, r"\.(txt|tsv\.gz|vcf\.gz)$" => "")
            output_files[key] = file_path
        end
    end
    
    println("vcfeval completed successfully.")
    println("Output files: $(length(output_files))")
    
    return output_files
end

"""
    prepare_vcf_for_rtg(vcf_file::String)

Prepare a VCF file for RTG Tools analysis by ensuring it is sorted, compressed, and indexed.

Performs necessary preprocessing steps to ensure VCF files are compatible with RTG Tools,
including sorting, compression with bgzip, and indexing with tabix.

# Arguments
- `vcf_file`: Path to input VCF file (can be compressed or uncompressed)

# Returns
- Path to the processed VCF file (sorted, compressed, indexed)

# Example
```julia
processed_vcf = Mycelia.prepare_vcf_for_rtg("variants.vcf")
```
"""
function prepare_vcf_for_rtg(vcf_file::String)
    
    if !isfile(vcf_file)
        error("VCF file does not exist: $(vcf_file)")
    end
    
    # Ensure required tools are available
    add_bioconda_env("bcftools")
    add_bioconda_env("htslib")
    add_bioconda_env("rtg-tools")
    
    # Determine output filename
    base_name = replace(vcf_file, r"\.vcf(\.gz)?$" => "")
    sorted_vcf = base_name * ".sorted.vcf"
    output_vcf = sorted_vcf * ".gz"
    
    # Skip if already processed
    if isfile(output_vcf) && isfile(output_vcf * ".tbi")
        return output_vcf
    end
    
    # Sort VCF file
    if !isfile(sorted_vcf)
        println("Sorting VCF file: $(vcf_file)")
        cmd = `$(CONDA_RUNNER) run --live-stream -n bcftools bcftools sort $(vcf_file) --output $(sorted_vcf)`
        run(cmd)
    end
    
    # Compress with RTG's bgzip
    if !isfile(output_vcf)
        println("Compressing VCF file: $(sorted_vcf)")
        cmd = `$(CONDA_RUNNER) run --live-stream -n rtg-tools rtg bgzip $(sorted_vcf)`
        run(cmd)
    end
    
    # Index with RTG's index command
    if !isfile(output_vcf * ".tbi")
        println("Indexing VCF file: $(output_vcf)")
        cmd = `$(CONDA_RUNNER) run --live-stream -n rtg-tools rtg index $(output_vcf)`
        run(cmd)
    end
    
    return output_vcf
end

"""
    generate_roc_plots(vcfeval_output_dir::String; output_formats::Vector{String}=["png", "svg"])

Generate ROC curve plots from vcfeval output using RTG Tools rocplot.

Creates publication-quality ROC curve plots comparing SNP and non-SNP performance
using the standard output files from RTG vcfeval.

# Arguments
- `vcfeval_output_dir`: Directory containing vcfeval output files
- `output_formats`: Image formats to generate (default: ["png", "svg"])

# Returns
- Dictionary mapping formats to output file paths

# Example
```julia
plots = Mycelia.generate_roc_plots("vcfeval_results", output_formats=["png", "svg", "pdf"])
```
"""
function generate_roc_plots(vcfeval_output_dir::String; output_formats::Vector{String}=["png", "svg"])
    
    if !isdir(vcfeval_output_dir)
        error("vcfeval output directory does not exist: $(vcfeval_output_dir)")
    end
    
    # Check for required ROC files
    roc_files = ["snp_roc.tsv.gz", "non_snp_roc.tsv.gz", "weighted_roc.tsv.gz"]
    missing_files = String[]
    
    for roc_file in roc_files
        file_path = joinpath(vcfeval_output_dir, roc_file)
        if !isfile(file_path)
            push!(missing_files, roc_file)
        end
    end
    
    if !isempty(missing_files)
        error("Missing ROC files: $(join(missing_files, ", "))")
    end
    
    # Ensure RTG Tools environment exists
    add_bioconda_env("rtg-tools")
    
    # Generate plots for each requested format
    plot_files = Dict{String, String}()
    
    for format in output_formats
        if !(format in ["png", "svg", "pdf"])
            @warn "Unsupported format: $(format). Skipping."
            continue
        end
        
        output_file = joinpath(vcfeval_output_dir, "roc.$(format)")
        
        # Build rocplot command
        cmd_args = [
            "RTG_MEM=4G", "rocplot",
            "--$(format)", output_file,
            "--curve", "$(joinpath(vcfeval_output_dir, "non_snp_roc.tsv.gz"))=Non-SNP",
            "--curve", "$(joinpath(vcfeval_output_dir, "snp_roc.tsv.gz"))=SNP", 
            "--curve", "$(joinpath(vcfeval_output_dir, "weighted_roc.tsv.gz"))=Weighted"
        ]
        
        println("Generating $(format) ROC plot...")
        cmd = `$(CONDA_RUNNER) run --live-stream -n rtg-tools rtg $(cmd_args)`
        
        try
            run(cmd)
            plot_files[format] = output_file
            println("ROC plot generated: $(output_file)")
        catch e
            @warn "Failed to generate $(format) plot: $(e)"
        end
    end
    
    return plot_files
end

"""
    evaluate_variant_calling_accuracy(baseline_vcf::String, calls_vcf::String, 
                                     reference_fasta::String, output_dir::String;
                                     generate_plots::Bool=true, 
                                     plot_formats::Vector{String}=["png", "svg"],
                                     kwargs...)

Complete variant calling accuracy evaluation workflow using RTG vcfeval.

Performs end-to-end evaluation of variant calling accuracy, including vcfeval
analysis, ROC curve generation, and results parsing.

# Arguments
- `baseline_vcf`: Path to baseline truth VCF file
- `calls_vcf`: Path to variant calls VCF file to evaluate  
- `reference_fasta`: Path to reference genome FASTA file
- `output_dir`: Directory for all output files
- `generate_plots`: Whether to generate ROC plots (default: true)
- `plot_formats`: Image formats for plots (default: ["png", "svg"])
- `kwargs...`: Additional arguments passed to run_vcfeval

# Returns
- Named tuple with vcfeval output files, parsed results, and plot files

# Example
```julia
results = Mycelia.evaluate_variant_calling_accuracy(
    "truth.vcf.gz", "calls.vcf.gz", "reference.fasta", "evaluation",
    threads=4, memory_gb=16
)

# Access parsed ROC data
snp_roc = results.parsed_results["snp_roc"]
println("SNP precision: \$(maximum(snp_roc.Precision))")
```
"""
function evaluate_variant_calling_accuracy(baseline_vcf::String, calls_vcf::String, 
                                         reference_fasta::String, output_dir::String;
                                         generate_plots::Bool=true, 
                                         plot_formats::Vector{String}=["png", "svg"],
                                         kwargs...)
    
    # Run vcfeval
    println("Step 1: Running RTG vcfeval...")
    vcfeval_files = run_vcfeval(baseline_vcf, calls_vcf, reference_fasta, output_dir; kwargs...)
    
    # Parse results
    println("Step 2: Parsing vcfeval results...")
    parsed_results = Dict{String, DataFrames.DataFrame}()
    
    for (key, file_path) in vcfeval_files
        if endswith(file_path, ".tsv.gz")
            try
                parsed_results[key] = parse_rtg_eval_output(file_path)
                println("  Parsed $(key): $(DataFrames.nrow(parsed_results[key])) rows")
            catch e
                @warn "Failed to parse $(key): $(e)"
            end
        end
    end
    
    # Generate plots if requested
    plot_files = Dict{String, String}()
    if generate_plots && haskey(vcfeval_files, "snp_roc")
        println("Step 3: Generating ROC plots...")
        try
            plot_files = generate_roc_plots(output_dir, output_formats=plot_formats)
        catch e
            @warn "Failed to generate plots: $(e)"
        end
    end
    
    # Calculate summary statistics
    summary_stats = calculate_evaluation_summary(parsed_results)
    
    println("Variant calling evaluation completed successfully!")
    println("Output directory: $(output_dir)")
    
    return (
        vcfeval_files = vcfeval_files,
        parsed_results = parsed_results,
        plot_files = plot_files,
        summary_stats = summary_stats,
        output_dir = output_dir
    )
end

"""
    calculate_evaluation_summary(parsed_results::Dict{String, DataFrames.DataFrame})

Calculate summary statistics from parsed vcfeval results.

Extracts key performance metrics from vcfeval ROC curve data, including
maximum F1 scores, area under curve estimates, and threshold values.

# Arguments
- `parsed_results`: Dictionary of parsed DataFrames from vcfeval output

# Returns
- Named tuple with summary statistics for SNPs, non-SNPs, and weighted results

# Example
```julia
summary = Mycelia.calculate_evaluation_summary(parsed_results)
println("Best SNP F1: \$(summary.snp_max_f1)")
```
"""
function calculate_evaluation_summary(parsed_results::Dict{String, DataFrames.DataFrame})
    
    summary_stats = Dict{String, Any}()
    
    for (roc_type, data) in parsed_results
        if endswith(roc_type, "_roc") && DataFrames.nrow(data) > 0
            # Calculate key metrics if columns exist
            metrics = Dict{String, Float64}()
            
            if "Precision" in names(data) && "Recall" in names(data)
                # Calculate F1 scores
                precision = data[!, "Precision"]
                recall = data[!, "Recall"]
                
                # Avoid division by zero
                valid_indices = (precision .+ recall) .> 0
                if any(valid_indices)
                    f1_scores = similar(precision, Float64)
                    f1_scores[valid_indices] = 2 .* (precision[valid_indices] .* recall[valid_indices]) ./ 
                                              (precision[valid_indices] .+ recall[valid_indices])
                    f1_scores[.!valid_indices] .= 0.0
                    
                    metrics["max_f1"] = maximum(f1_scores)
                    metrics["max_precision"] = maximum(precision)
                    metrics["max_recall"] = maximum(recall)
                    
                    # Find threshold at max F1
                    if "Score" in names(data)
                        max_f1_idx = argmax(f1_scores)
                        metrics["optimal_threshold"] = data[max_f1_idx, "Score"]
                    end
                end
            end
            
            # Store area under curve estimate (simple trapezoidal)
            if "Precision" in names(data) && "Recall" in names(data) && DataFrames.nrow(data) > 1
                recall_vals = sort(data[!, "Recall"])
                precision_vals = data[sortperm(data[!, "Recall"]), "Precision"]
                
                if length(recall_vals) > 1
                    auc = 0.0
                    for i in 2:length(recall_vals)
                        width = recall_vals[i] - recall_vals[i-1]
                        height = (precision_vals[i] + precision_vals[i-1]) / 2
                        auc += width * height
                    end
                    metrics["auc_estimate"] = auc
                end
            end
            
            summary_stats[roc_type] = metrics
        end
    end
    
    return NamedTuple(Symbol(k) => NamedTuple(Symbol(k2) => v2 for (k2, v2) in v) for (k, v) in summary_stats)
end

# Variant Caller Integration Functions

"""
    run_gatk_haplotypecaller(bam_file::String, reference_fasta::String, output_vcf::String;
                            ploidy::Int=1, memory_gb::Int=8, threads::Int=1,
                            additional_args::Vector{String}=String[])

Run GATK HaplotypeCaller for variant calling from BAM alignments.

Uses GATK4 HaplotypeCaller to identify variants from aligned reads, with support
for both haploid and diploid calling modes.

# Arguments
- `bam_file`: Path to sorted BAM file with alignments
- `reference_fasta`: Path to reference genome FASTA file
- `output_vcf`: Path for output VCF file
- `ploidy`: Sample ploidy level (default: 1 for haploid)
- `memory_gb`: Memory allocation in GB (default: 8)
- `threads`: Number of threads (default: 1)
- `additional_args`: Additional GATK arguments

# Returns
- Path to output VCF file

# Example
```julia
vcf_file = Mycelia.run_gatk_haplotypecaller(
    "aligned.bam", "reference.fasta", "variants.gatk.vcf",
    ploidy=2, memory_gb=16
)
```
"""
function run_gatk_haplotypecaller(bam_file::String, reference_fasta::String, output_vcf::String;
                                 ploidy::Int=1, memory_gb::Int=8, threads::Int=1,
                                 additional_args::Vector{String}=String[])
    
    # Validate input files
    for file in [bam_file, reference_fasta]
        if !isfile(file)
            error("Input file does not exist: $(file)")
        end
    end
    
    # Ensure required tools are available
    add_bioconda_env("gatk4")
    add_bioconda_env("picard")
    
    # Create sequence dictionary if it doesn't exist
    dict_file = replace(reference_fasta, r"\.(fasta|fa|fna)$" => ".dict")
    if !isfile(dict_file)
        println("Creating sequence dictionary for GATK...")
        cmd = `$(CONDA_RUNNER) run --live-stream -n picard picard CreateSequenceDictionary REFERENCE=$(reference_fasta) OUTPUT=$(dict_file)`
        run(cmd)
    end
    
    # Ensure BAM is indexed
    bam_index = bam_file * ".bai"
    if !isfile(bam_index)
        println("Indexing BAM file...")
        add_bioconda_env("samtools")
        run(`$(CONDA_RUNNER) run --live-stream -n samtools samtools index $(bam_file)`)
    end
    
    # Create output directory
    mkpath(dirname(output_vcf))
    
    # Build GATK command
    gatk_args = [
        "--java-options", "-Xmx$(memory_gb)G",
        "HaplotypeCaller",
        "--sample-ploidy", string(ploidy),
        "--reference", reference_fasta,
        "--input", bam_file,
        "--output", output_vcf
    ]
    
    # Add additional arguments
    append!(gatk_args, additional_args)
    
    println("Running GATK HaplotypeCaller...")
    println("  BAM: $(basename(bam_file))")
    println("  Reference: $(basename(reference_fasta))")
    println("  Ploidy: $(ploidy)")
    println("  Output: $(output_vcf)")
    
    cmd = `$(CONDA_RUNNER) run --live-stream -n gatk4 gatk $(gatk_args)`
    
    try
        run(cmd)
    catch e
        error("GATK HaplotypeCaller failed: $(e)")
    end
    
    println("GATK variant calling completed: $(output_vcf)")
    
    return output_vcf
end

"""
    run_freebayes(bam_file::String, reference_fasta::String, output_vcf::String;
                  ploidy::Int=1, best_n_alleles::Int=1, 
                  additional_args::Vector{String}=String[])

Run Freebayes for variant calling from BAM alignments.

Uses Freebayes to call variants from aligned reads with support for 
different ploidy levels and allele selection strategies.

# Arguments
- `bam_file`: Path to sorted BAM file with alignments
- `reference_fasta`: Path to reference genome FASTA file  
- `output_vcf`: Path for output VCF file
- `ploidy`: Sample ploidy level (default: 1)
- `best_n_alleles`: Number of best alleles to use (default: 1)
- `additional_args`: Additional Freebayes arguments

# Returns
- Path to output VCF file

# Example
```julia
vcf_file = Mycelia.run_freebayes(
    "aligned.bam", "reference.fasta", "variants.freebayes.vcf",
    ploidy=2
)
```
"""
function run_freebayes(bam_file::String, reference_fasta::String, output_vcf::String;
                      ploidy::Int=1, best_n_alleles::Int=1, 
                      additional_args::Vector{String}=String[])
    
    # Validate input files
    for file in [bam_file, reference_fasta]
        if !isfile(file)
            error("Input file does not exist: $(file)")
        end
    end
    
    # Ensure Freebayes environment exists
    add_bioconda_env("freebayes")
    
    # Create output directory
    mkpath(dirname(output_vcf))
    
    # Build Freebayes command
    freebayes_args = [
        "--ploidy", string(ploidy),
        "--use-best-n-alleles", string(best_n_alleles),
        "--fasta-reference", reference_fasta,
        "--bam", bam_file,
        "--vcf", output_vcf
    ]
    
    # Add additional arguments
    append!(freebayes_args, additional_args)
    
    println("Running Freebayes...")
    println("  BAM: $(basename(bam_file))")
    println("  Reference: $(basename(reference_fasta))")
    println("  Ploidy: $(ploidy)")
    println("  Output: $(output_vcf)")
    
    cmd = `$(CONDA_RUNNER) run --live-stream -n freebayes freebayes $(freebayes_args)`
    
    try
        run(cmd)
    catch e
        error("Freebayes failed: $(e)")
    end
    
    println("Freebayes variant calling completed: $(output_vcf)")
    
    return output_vcf
end

"""
    run_clair3(bam_file::String, reference_fasta::String, output_dir::String;
               platform::String="ilmn", threads::Int=4, haploid_precise::Bool=true,
               model_path::String="", additional_args::Vector{String}=String[])

Run Clair3 for variant calling from BAM alignments.

Uses Clair3 deep learning-based variant caller with platform-specific models
for Illumina or Oxford Nanopore data.

# Arguments
- `bam_file`: Path to sorted and indexed BAM file
- `reference_fasta`: Path to reference genome FASTA file
- `output_dir`: Directory for Clair3 output files
- `platform`: Sequencing platform ("ilmn" or "ont", default: "ilmn")
- `threads`: Number of threads (default: 4)
- `haploid_precise`: Use haploid precise mode (default: true)
- `model_path`: Custom model path (auto-detected if empty)
- `additional_args`: Additional Clair3 arguments

# Returns
- Path to main output VCF file (merge_output.vcf.gz)

# Example
```julia
vcf_file = Mycelia.run_clair3(
    "aligned.bam", "reference.fasta", "clair3_output",
    platform="ont", threads=8
)
```
"""
function run_clair3(bam_file::String, reference_fasta::String, output_dir::String;
                   platform::String="ilmn", threads::Int=4, haploid_precise::Bool=true,
                   model_path::String="", additional_args::Vector{String}=String[])
    
    # Validate input files
    for file in [bam_file, reference_fasta]
        if !isfile(file)
            error("Input file does not exist: $(file)")
        end
    end
    
    # Validate platform
    if !(platform in ["ilmn", "ont"])
        error("Platform must be 'ilmn' or 'ont', got: $(platform)")
    end
    
    # Ensure required tools are available
    if platform == "ilmn"
        add_bioconda_env("clair3-illumina")
        env_name = "clair3-illumina"
    else
        add_bioconda_env("clair3")
        env_name = "clair3"
    end
    add_bioconda_env("samtools")
    
    # Ensure BAM is indexed
    bam_index = bam_file * ".bai"
    if !isfile(bam_index)
        println("Indexing BAM file...")
        run(`$(CONDA_RUNNER) run --live-stream -n samtools samtools index $(bam_file)`)
    end
    
    # Set model path if not provided
    if isempty(model_path)
        if platform == "ilmn"
            # Get conda environment prefix and set model path
            prefix_cmd = `$(CONDA_RUNNER) run -n $(env_name) python -c "import sys; print(sys.prefix)"`
            conda_prefix = strip(read(prefix_cmd, String))
            model_path = joinpath(conda_prefix, "bin", "models", "ilmn")
        else
            # ONT model path
            prefix_cmd = `$(CONDA_RUNNER) run -n $(env_name) python -c "import sys; print(sys.prefix)"`
            conda_prefix = strip(read(prefix_cmd, String))
            model_path = joinpath(conda_prefix, "bin", "models", "r941_prom_sup_g5014")
        end
    end
    
    # Create output directory
    mkpath(output_dir)
    
    # Build Clair3 command
    clair3_args = [
        "--bam_fn", bam_file,
        "--ref_fn", reference_fasta,
        "--threads", string(threads),
        "--platform", platform,
        "--model_path", model_path,
        "--output", output_dir,
        "--include_all_ctgs",
        "--no_phasing_for_fa"
    ]
    
    # Add haploid mode if requested
    if haploid_precise
        push!(clair3_args, "--haploid_precise")
    end
    
    # Add additional arguments
    append!(clair3_args, additional_args)
    
    println("Running Clair3...")
    println("  BAM: $(basename(bam_file))")
    println("  Reference: $(basename(reference_fasta))")
    println("  Platform: $(platform)")
    println("  Threads: $(threads)")
    println("  Output: $(output_dir)")
    
    cmd = `$(CONDA_RUNNER) run --live-stream -n $(env_name) run_clair3.sh $(clair3_args)`
    
    try
        run(cmd)
    catch e
        error("Clair3 failed: $(e)")
    end
    
    # Return path to main output VCF
    output_vcf = joinpath(output_dir, "merge_output.vcf.gz")
    
    if !isfile(output_vcf)
        error("Clair3 output VCF not found: $(output_vcf)")
    end
    
    println("Clair3 variant calling completed: $(output_vcf)")
    
    return output_vcf
end

"""
    run_bcftools_call(bam_file::String, reference_fasta::String, output_vcf::String;
                      threads::Int=4, multiallelic::Bool=true, variants_only::Bool=true,
                      additional_mpileup_args::Vector{String}=String[],
                      additional_call_args::Vector{String}=String[])

Run BCFtools mpileup and call pipeline for variant calling.

Uses the BCFtools two-step process (mpileup followed by call) to identify
variants from aligned reads.

# Arguments
- `bam_file`: Path to sorted BAM file with alignments
- `reference_fasta`: Path to reference genome FASTA file
- `output_vcf`: Path for output VCF file
- `threads`: Number of threads (default: 4)
- `multiallelic`: Enable multiallelic calling (default: true)
- `variants_only`: Output variants only, skip reference sites (default: true)
- `additional_mpileup_args`: Additional arguments for mpileup step
- `additional_call_args`: Additional arguments for call step

# Returns
- Path to output VCF file

# Example
```julia
vcf_file = Mycelia.run_bcftools_call(
    "aligned.bam", "reference.fasta", "variants.bcftools.vcf",
    threads=8, multiallelic=true
)
```
"""
function run_bcftools_call(bam_file::String, reference_fasta::String, output_vcf::String;
                          threads::Int=4, multiallelic::Bool=true, variants_only::Bool=true,
                          additional_mpileup_args::Vector{String}=String[],
                          additional_call_args::Vector{String}=String[])
    
    # Validate input files
    for file in [bam_file, reference_fasta]
        if !isfile(file)
            error("Input file does not exist: $(file)")
        end
    end
    
    # Ensure BCFtools environment exists
    add_bioconda_env("bcftools")
    
    # Create output directory
    mkpath(dirname(output_vcf))
    
    # Build mpileup arguments
    mpileup_args = [
        "mpileup",
        "--threads", string(threads),
        "--output-type", "u",  # Uncompressed BCF
        "--fasta-ref", reference_fasta
    ]
    append!(mpileup_args, additional_mpileup_args)
    push!(mpileup_args, bam_file)
    
    # Build call arguments
    call_args = [
        "call",
        "--threads", string(threads),
        "--output-type", "v",  # VCF format
        "--output", output_vcf
    ]
    
    if multiallelic
        push!(call_args, "--multiallelic-caller")
    end
    
    if variants_only
        push!(call_args, "--variants-only")
    end
    
    append!(call_args, additional_call_args)
    
    println("Running BCFtools mpileup | call pipeline...")
    println("  BAM: $(basename(bam_file))")
    println("  Reference: $(basename(reference_fasta))")
    println("  Threads: $(threads)")
    println("  Output: $(output_vcf)")
    
    # Build pipeline command
    mpileup_cmd = `$(CONDA_RUNNER) run --live-stream -n bcftools bcftools $(mpileup_args)`
    call_cmd = `$(CONDA_RUNNER) run --live-stream -n bcftools bcftools $(call_args)`
    
    try
        run(pipeline(mpileup_cmd, call_cmd))
    catch e
        error("BCFtools pipeline failed: $(e)")
    end
    
    println("BCFtools variant calling completed: $(output_vcf)")
    
    return output_vcf
end

"""
    run_variant_calling_comparison(bam_file::String, reference_fasta::String, output_dir::String;
                                  callers::Vector{String}=["gatk", "freebayes", "clair3", "bcftools"],
                                  platform::String="ilmn", threads::Int=4, 
                                  run_evaluation::Bool=true, baseline_vcf::String="")

Run multiple variant callers and optionally compare their performance.

Executes multiple variant calling algorithms on the same BAM file and reference,
with optional evaluation against a baseline truth set using vcfeval.

# Arguments
- `bam_file`: Path to sorted BAM file with alignments
- `reference_fasta`: Path to reference genome FASTA file
- `output_dir`: Directory for all output files
- `callers`: List of variant callers to run (default: all four)
- `platform`: Sequencing platform for Clair3 ("ilmn" or "ont", default: "ilmn")
- `threads`: Number of threads for parallel processing
- `run_evaluation`: Whether to run vcfeval comparison (default: true)
- `baseline_vcf`: Path to baseline truth VCF for evaluation (required if run_evaluation=true)

# Returns
- Dictionary with paths to VCF files from each caller and evaluation results

# Example
```julia
results = Mycelia.run_variant_calling_comparison(
    "aligned.bam", "reference.fasta", "comparison_output",
    callers=["gatk", "freebayes"], 
    baseline_vcf="truth.vcf.gz"
)
```
"""
function run_variant_calling_comparison(bam_file::String, reference_fasta::String, output_dir::String;
                                      callers::Vector{String}=["gatk", "freebayes", "clair3", "bcftools"],
                                      platform::String="ilmn", threads::Int=4, 
                                      run_evaluation::Bool=true, baseline_vcf::String="")
    
    # Validate inputs
    for file in [bam_file, reference_fasta]
        if !isfile(file)
            error("Input file does not exist: $(file)")
        end
    end
    
    if run_evaluation && isempty(baseline_vcf)
        error("baseline_vcf must be provided when run_evaluation=true")
    end
    
    if run_evaluation && !isfile(baseline_vcf)
        error("Baseline VCF file does not exist: $(baseline_vcf)")
    end
    
    # Create output directory
    mkpath(output_dir)
    
    # Run each variant caller
    caller_outputs = Dict{String, String}()
    
    println("Running variant calling comparison with $(length(callers)) callers...")
    
    for caller in callers
        caller_dir = joinpath(output_dir, caller)
        mkpath(caller_dir)
        
        try
            if caller == "gatk"
                output_vcf = joinpath(caller_dir, "variants.gatk.vcf")
                caller_outputs[caller] = run_gatk_haplotypecaller(bam_file, reference_fasta, output_vcf, threads=threads)
                
            elseif caller == "freebayes"
                output_vcf = joinpath(caller_dir, "variants.freebayes.vcf")
                caller_outputs[caller] = run_freebayes(bam_file, reference_fasta, output_vcf)
                
            elseif caller == "clair3"
                caller_outputs[caller] = run_clair3(bam_file, reference_fasta, caller_dir, platform=platform, threads=threads)
                
            elseif caller == "bcftools"
                output_vcf = joinpath(caller_dir, "variants.bcftools.vcf")
                caller_outputs[caller] = run_bcftools_call(bam_file, reference_fasta, output_vcf, threads=threads)
                
            else
                @warn "Unknown variant caller: $(caller). Skipping."
                continue
            end
            
            println("Completed $(caller): $(caller_outputs[caller])")
            
        catch e
            @warn "$(caller) failed: $(e)"
        end
    end
    
    # Run evaluation if requested
    evaluation_results = Dict{String, Any}()
    
    if run_evaluation && !isempty(caller_outputs)
        println("Running vcfeval evaluation against baseline...")
        
        for (caller, vcf_file) in caller_outputs
            eval_dir = joinpath(output_dir, "evaluation", caller)
            
            try
                # Normalize VCF before evaluation
                normalized_vcf = normalize_vcf(reference_fasta=reference_fasta, vcf_file=vcf_file)
                
                # Run evaluation
                eval_results = evaluate_variant_calling_accuracy(
                    baseline_vcf, normalized_vcf, reference_fasta, eval_dir,
                    threads=threads
                )
                
                evaluation_results[caller] = eval_results
                println("Evaluation completed for $(caller)")
                
            catch e
                @warn "Evaluation failed for $(caller): $(e)"
            end
        end
    end
    
    # Generate summary report
    summary_file = joinpath(output_dir, "comparison_summary.txt")
    generate_comparison_summary(caller_outputs, evaluation_results, summary_file)
    
    println("Variant calling comparison completed!")
    println("Output directory: $(output_dir)")
    println("Summary: $(summary_file)")
    
    return (
        caller_outputs = caller_outputs,
        evaluation_results = evaluation_results,
        summary_file = summary_file,
        output_dir = output_dir
    )
end

"""
    generate_comparison_summary(caller_outputs::Dict{String, String}, 
                               evaluation_results::Dict{String, Any}, 
                               summary_file::String)

Generate a text summary comparing variant caller performance.

Creates a readable summary report comparing the number of variants called
and evaluation metrics across different variant callers.

# Arguments
- `caller_outputs`: Dictionary mapping caller names to VCF file paths
- `evaluation_results`: Dictionary mapping caller names to evaluation results  
- `summary_file`: Path for output summary text file

# Example
```julia
Mycelia.generate_comparison_summary(caller_outputs, evaluation_results, "summary.txt")
```
"""
function generate_comparison_summary(caller_outputs::Dict{String, String}, 
                                   evaluation_results::Dict{String, Any}, 
                                   summary_file::String)
    
    open(summary_file, "w") do io
        println(io, "# Variant Calling Comparison Summary")
        println(io, "Generated: $(Dates.now())")
        println(io, "")
        
        # Caller outputs section
        println(io, "## Variant Callers Run")
        for (caller, vcf_file) in caller_outputs
            println(io, "- $(caller): $(basename(vcf_file))")
        end
        println(io, "")
        
        # Evaluation results section
        if !isempty(evaluation_results)
            println(io, "## Evaluation Results")
            println(io, "")
            
            # Create comparison table
            println(io, "| Caller | SNP Max F1 | Non-SNP Max F1 | Weighted Max F1 |")
            println(io, "|--------|------------|-----------------|-----------------|")
            
            for (caller, results) in evaluation_results
                snp_f1 = haskey(results.summary_stats, :snp_roc) ? 
                        get(results.summary_stats.snp_roc, :max_f1, "N/A") : "N/A"
                nonsnp_f1 = haskey(results.summary_stats, :non_snp_roc) ? 
                           get(results.summary_stats.non_snp_roc, :max_f1, "N/A") : "N/A"
                weighted_f1 = haskey(results.summary_stats, :weighted_roc) ? 
                             get(results.summary_stats.weighted_roc, :max_f1, "N/A") : "N/A"
                
                println(io, "| $(caller) | $(snp_f1) | $(nonsnp_f1) | $(weighted_f1) |")
            end
            
            println(io, "")
            
            # Best performer
            best_caller = ""
            best_f1 = 0.0
            
            for (caller, results) in evaluation_results
                if haskey(results.summary_stats, :weighted_roc) && 
                   haskey(results.summary_stats.weighted_roc, :max_f1)
                    f1 = results.summary_stats.weighted_roc.max_f1
                    if f1 > best_f1
                        best_f1 = f1
                        best_caller = caller
                    end
                end
            end
            
            if !isempty(best_caller)
                println(io, "**Best Overall Performance**: $(best_caller) (Weighted F1: $(round(best_f1, digits=4)))")
            end
        else
            println(io, "## Evaluation Results")
            println(io, "No evaluation performed (no baseline VCF provided)")
        end
        
        println(io, "")
        println(io, "## Files Generated")
        for (caller, vcf_file) in caller_outputs
            println(io, "- $(caller) VCF: $(vcf_file)")
            if haskey(evaluation_results, caller)
                eval_dir = evaluation_results[caller].output_dir
                println(io, "  - Evaluation: $(eval_dir)/")
            end
        end
    end
    
    println("Summary report written to: $(summary_file)")
end
