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