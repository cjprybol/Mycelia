"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run MEGAHIT assembler for metagenomic short read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "megahit_output")
- `min_contig_len::Int`: Minimum contig length (default: 200)
- `k_list::String`: k-mer sizes to use (default: "21,29,39,59,79,99,119,141")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to final contigs file

# Details
- Automatically creates and uses a conda environment with megahit
- Optimized for metagenomic assemblies with varying coverage
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_megahit(;fastq1, fastq2=nothing, outdir="megahit_output", min_contig_len=200, k_list="21,29,39,59,79,99,119,141")
    Mycelia.add_bioconda_env("megahit")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "final.contigs.fa"))
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit -r $(fastq1) -o $(outdir) --min-contig-len $(min_contig_len) --k-list $(k_list) -t $(Sys.CPU_THREADS)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit -1 $(fastq1) -2 $(fastq2) -o $(outdir) --min-contig-len $(min_contig_len) --k-list $(k_list) -t $(Sys.CPU_THREADS)`)
        end
    end
    return (;outdir, contigs=joinpath(outdir, "final.contigs.fa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaSPAdes assembler for metagenomic short read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "metaspades_output")
- `k_list::String`: k-mer sizes to use (default: "21,33,55,77")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `scaffolds::String`: Path to scaffolds file

# Details
- Automatically creates and uses a conda environment with spades
- Designed for metagenomic datasets with uneven coverage
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_metaspades(;fastq1, fastq2=nothing, outdir="metaspades_output", k_list="21,33,55,77")
    Mycelia.add_bioconda_env("spades")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "contigs.fasta"))
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades metaspades.py -s $(fastq1) -o $(outdir) -k $(k_list) -t $(Sys.CPU_THREADS)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades metaspades.py -1 $(fastq1) -2 $(fastq2) -o $(outdir) -k $(k_list) -t $(Sys.CPU_THREADS)`)
        end
    end
    return (;outdir, contigs=joinpath(outdir, "contigs.fasta"), scaffolds=joinpath(outdir, "scaffolds.fasta"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Flye assembler for long read assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "flye_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g")
- `read_type::String`: Type of reads ("pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Automatically creates and uses a conda environment with flye
- Supports various long read technologies and quality levels
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_flye(;fastq, outdir="flye_output", genome_size, read_type="pacbio-hifi")
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "assembly.fasta"))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye flye --$(read_type) $(fastq) --out-dir $(outdir) --genome-size $(genome_size) --threads $(Sys.CPU_THREADS)`)
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Canu assembler for long read assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "canu_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g")
- `read_type::String`: Type of reads ("pacbio", "nanopore")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Automatically creates and uses a conda environment with canu
- Includes error correction, trimming, and assembly stages
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_canu(;fastq, outdir="canu_output", genome_size, read_type="pacbio")
    Mycelia.add_bioconda_env("canu")
    mkpath(outdir)
    
    prefix = basename(fastq, ".fastq")
    if !isfile(joinpath(outdir, "$(prefix).contigs.fasta"))
        if read_type == "pacbio"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -pacbio $(fastq) maxThreads=$(Sys.CPU_THREADS)`)
        elseif read_type == "nanopore"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -nanopore $(fastq) maxThreads=$(Sys.CPU_THREADS)`)
        else
            error("Unsupported read type: $(read_type). Use 'pacbio' or 'nanopore'.")
        end
    end
    return (;outdir, assembly=joinpath(outdir, "$(prefix).contigs.fasta"))
end

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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run hybrid assembly combining short and long reads using Unicycler.

# Arguments
- `short_1::String`: Path to first short read FASTQ file
- `short_2::String`: Path to second short read FASTQ file (optional)
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "unicycler_output")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file

# Details
- Automatically creates and uses a conda environment with unicycler
- Combines short read accuracy with long read scaffolding
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_unicycler(;short_1, short_2=nothing, long_reads, outdir="unicycler_output")
    Mycelia.add_bioconda_env("unicycler")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "assembly.fasta"))
        if isnothing(short_2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -s $(short_1) -l $(long_reads) -o $(outdir) -t $(Sys.CPU_THREADS)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -1 $(short_1) -2 $(short_2) -l $(long_reads) -o $(outdir) -t $(Sys.CPU_THREADS)`)
        end
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"))
end