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