"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the hifiasm genome assembler on PacBio HiFi reads.

# Arguments
- `fastq::String`: Path to input FASTQ file containing HiFi reads
- `outdir::String`: Output directory path (default: "\${basename(fastq)}_hifiasm")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `hifiasm_outprefix::String`: Prefix used for hifiasm output files

# Details
- Automatically creates and uses a conda environment with hifiasm
- Uses primary assembly mode (--primary) optimized for inbred samples
- Skips assembly if output files already exist at the specified prefix
- Utilizes all available CPU threads
"""
function run_hifiasm(;fastq, outdir=basename(fastq) * "_hifiasm")
    Mycelia.add_bioconda_env("hifiasm")
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm")
    hifiasm_outputs = filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join=true))
    # https://hifiasm.readthedocs.io/en/latest/faq.html#are-inbred-homozygous-genomes-supported
    if isempty(hifiasm_outputs)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm hifiasm --primary -l0 -o $(hifiasm_outprefix) -t $(Sys.CPU_THREADS) $(fastq)`)
    end
    return (;outdir, hifiasm_outprefix)
end