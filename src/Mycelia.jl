module Mycelia

__precompile__(false)

import AlgebraOfGraphics
import ArgParse
import Arrow
import BioAlignments
import BioSequences
import BioSymbols
import CairoMakie
import Clustering
import CodecBase
import CodecBzip2
import CodecZlib
import Colors
import ColorSchemes
import Conda
import CSV
import DataFrames
import DataStructures
import Dates
import DelimitedFiles
import Distances
import Distributions
import DocStringExtensions
import Downloads
import FASTX
import FileIO
import GenomicAnnotations
import GeoMakie
import GFF3
import GLM
import GraphMakie
import Graphs
import HDF5
import HTTP
import JLD2
import JSON
import Karnak
import Kmers
import LsqFit
import Luxor
import Makie
import MetaGraphs
import Mmap
import OrderedCollections
import Plots
import Primes
import ProgressMeter
import Random
import SHA
import SparseArrays
import Statistics
import StatsBase
import StatsPlots
import Tar
import TopoPlots
import TranscodingStreams
import uCSV
import XAM
import XMLDict
import UUIDs
import StableRNGs

import Pkg

# preserve definitions between code jldoctest code blocks
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Preserving-Definitions-Between-Blocks
# use this to build up a story as we go, where outputs of earlier defined functions feed into
# tests for further downstream functions

# if things fall out of date but look correct, update them automatically
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Fixing-Outdated-Doctests

const METADATA = joinpath(dirname(dirname(pathof(Mycelia))), "docs", "metadata")
const DNA_ALPHABET = BioSymbols.ACGT
const RNA_ALPHABET = BioSymbols.ACGU

# fix new error
# ENV["MAMBA_ROOT_PREFIX"] = joinpath(DEPOT_PATH[1], "conda", "3", "x86_64")

# Mycelia.NERSC_MEM * .95
# const NERSC_MEM=512
# const NERSC_MEM=480
const NERSC_MEM=460
# const NERSC_CPU=240
const NERSC_CPU=240
# const AA_ALPHABET = filter(
#     x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x) || BioSymbols.isterm(x)),
#     BioSymbols.alphabet(BioSymbols.AminoAcid))
const AA_ALPHABET = filter(
    x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x)),
    BioSymbols.alphabet(BioSymbols.AminoAcid))

# can add support for conda too if needed
# const CONDA_RUNNER = joinpath(Conda.BINDIR, "mamba")
const CONDA_RUNNER = joinpath(Conda.BINDIR, "conda")
const FASTQ_REGEX = r"\.(fq\.gz|fastq\.gz|fastq|fq)$"
const FASTA_REGEX = r"\.(fa\.gz|fasta\.gz|fna\.gz|fasta|fa|fna)$"
const VCF_REGEX = r"\.(vcf|vcf\.gz)$"
# none of this code currently supports CRAM
const XAM_REGEX = r"\.(sam|sam\.gz|bam)$"

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a new Conda environment with a specified Bioconda package.

# Arguments
- `pkg::String`: Package name to install. Can include channel specification using 
  the format "channel::package"

# Keywords
- `force::Bool=false`: If true, recreates the environment even if it already exists

# Details
The function creates a new Conda environment named after the package and installs
the package into it. It uses channel priority: conda-forge > bioconda > defaults.
If CONDA_RUNNER is set to 'mamba', it will ensure mamba is installed first.

# Examples
```julia
# Install basic package
add_bioconda_env("blast")

# Install from specific channel
add_bioconda_env("bioconda::blast")

# Force reinstallation
add_bioconda_env("blast", force=true)
```
# Notes
- Requires Conda.jl to be installed and configured
- Uses CONDA_RUNNER global variable to determine whether to use conda or mamba
- Cleans conda cache after installation
"""
function add_bioconda_env(pkg; force=false)
    # ensure conda environment is available
    if !isfile(CONDA_RUNNER)
        if (basename(CONDA_RUNNER) == "mamba")
            Conda.add("mamba")
        elseif (basename(CONDA_RUNNER) == "conda")
            Conda.update()
        end
    end
    # try
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
    channel = nothing
    if occursin("::", pkg)
        println("splitting $(pkg)")
        channel, pkg = split(pkg, "::")
        println("into channel:$(channel) pkg:$(pkg)")
    end
    if !(pkg in current_environments) || force
        @info "installing conda environment $(pkg)"
        if isnothing(channel)
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
        else
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y`)
        end
        run(`$(CONDA_RUNNER) clean --all -y`)
    # else
    #     # @info "conda environment $(pkg) already present; set force=true to update/re-install"
    end
    # catch
    #     add_bioconda_envs()
    #     current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
    #     if !(pkg in current_environments) || force
    #         @info "installing conda environment $(pkg)"
    #         run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
    #         run(`$(CONDA_RUNNER) clean --all -y`)
    #     else
    #         # @info "conda environment $(pkg) already present; set force=true to update/re-install"
    #     end
    # end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Update a package and its dependencies in its dedicated Conda environment.

# Arguments
- `pkg::String`: Name of the package/environment to update
"""
function update_bioconda_env(pkg)
    run(`$(CONDA_RUNNER) update -n $(pkg) $(pkg) -y`)
    # conda update --all -n <env_name>
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function add_bioconda_envs(;all=false, force=false)
#     if !isfile(CONDA_RUNNER) && (basename(CONDA_RUNNER) == "mamba")
#         Conda.add("mamba")
#     end
#     if !isfile(joinpath(Conda.BINDIR, "pigz"))
#         run(`$(CONDA_RUNNER) install pigz -y`)
#     end
#     current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
#     # https://github.com/JuliaPy/Conda.jl/issues/185#issuecomment-1145149905
#     if all
#         for pkg in [
#             "art",
#             # "bioconvert",
#             "badread",
#             "bcftools",
#             "bedtools",
#             "blast",
#             "clair3-illumina",
#             "clair3",    
#             # "bwa",
#             # "bwa-mem2",
#             # "deepvariant",
#             "emboss",
#             "filtlong",
#             # "freebayes",
#             "flye",
#             "gatk4",
#             # "gffread",
#             "htslib",
#             "megahit",
#             "medaka",
#             "minimap2",
#             "mmseqs2",
#             "nanocaller",
#             "nanovar",
#             # "nanoq",
#             # "nanosim",
#             # "nanosim-h",
#             "ncbi-datasets-cli",
#             "pggb",
#             "picard",
#             # "polypolish",
#             "prodigal",
#             "raven-assembler",
#             "rtg-tools",
#             "samtools",
#             "sniffles",
#             "sourmash",
#             "spades",
#             "tabix",
#             "transtermhp",
#             "trim-galore",
#             "vcftools",
#             "vg"
#             ]
#             if !(pkg in current_environments) || force
#                 @info "installing conda environment $(pkg)"
#                 add_bioconda_env(pkg)
#             else
#                 @info "conda environment $(pkg) already present; set force=true to update/re-install"
#             end
#         end
#     end
#     run(`$(CONDA_RUNNER) clean --all -y`)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to SLURM scheduler on Lawrence Berkeley Lab's Lawrencium cluster.

# Arguments
- `job_name`: Name identifier for the SLURM job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type ("ALL", "BEGIN", "END", "FAIL", or "NONE")
- `logdir`: Directory for SLURM output and error logs
- `partition`: Lawrencium compute partition
- `qos`: Quality of Service level
- `account`: Project account for billing
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to spawn
- `time`: Wall time limit in format "days-hours:minutes:seconds"
- `cpus_per_task`: CPU cores per task
- `mem_gb`: Memory per node in GB
- `cmd`: Shell command to execute

# Returns
- `true` if submission was successful

# Note
Function includes 5-second delays before and after submission for queue stability.
"""
function lawrencium_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        partition::String="lr3",
        qos::String="lr_normal",
        account::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="3-00:00:00",
        cpus_per_task::Int=16,
        mem_gb::Int=64,
        cmd::String
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --partition=$(partition)
    --qos=$(qos)
    --account=$(account)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to SLURM using sbatch with specified parameters.

# Arguments
- `job_name::String`: Name identifier for the SLURM job
- `mail_user::String`: Email address for job notifications
- `mail_type::String`: Type of mail notifications (default: "ALL")
- `logdir::String`: Directory for error and output logs (default: "~/workspace/slurmlogs")
- `partition::String`: SLURM partition to submit job to
- `account::String`: Account to charge for compute resources
- `nodes::Int`: Number of nodes to allocate (default: 1)
- `ntasks::Int`: Number of tasks to run (default: 1)
- `time::String`: Maximum wall time in format "days-hours:minutes:seconds" (default: "1-00:00:00")
- `cpus_per_task::Int`: CPUs per task (default: 1)
- `mem_gb::Int`: Memory in GB, defaults to 32GB per CPU
- `cmd::String`: Command to execute

# Returns
- `Bool`: Returns true if submission succeeds

# Notes
- Function includes 5-second delays before and after submission
- Memory is automatically scaled with CPU count
- Log files are named with job ID (%j) and job name (%x)
"""
function scg_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath("$(homedir())/workspace/slurmlogs"),
        partition::String,
        account::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 32,
        cmd::String
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --partition=$(partition)
    --account=$(account)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a job to NERSC's SLURM scheduler using the shared QOS (Quality of Service).

# Arguments
- `job_name`: Identifier for the job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type ("ALL", "BEGIN", "END", "FAIL", "REQUEUE", "STAGE_OUT")
- `logdir`: Directory for storing job output and error logs
- `qos`: Quality of Service level ("shared", "regular", "preempt", "premium")
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to run
- `time`: Maximum wall time in format "days-hours:minutes:seconds"
- `cpus_per_task`: Number of CPUs per task
- `mem_gb`: Memory per node in GB (default: 2GB per CPU)
- `cmd`: Command to execute
- `constraint`: Node type constraint ("cpu" or "gpu")

# Resource Limits
- Maximum memory per node: 512GB
- Maximum cores per node: 128
- Default memory allocation: 2GB per CPU requested

# QOS Options
- shared: Default QOS for shared node usage
- regular: Standard priority
- preempt: Reduced credit usage but preemptible
- premium: 5x throughput priority (limited usage)

# Returns
`true` if job submission succeeds


https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priority runs limited to 5x throughput)

max request is 512Gb memory and 128 cores per node

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
"""
function nersc_sbatch_shared(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath("$(homedir())/workspace/slurmlogs"),
        qos::String="shared",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="2-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 2,
        cmd::String,
        constraint::String="cpu"
    )
    submission = 
    `sbatch
    --job-name=$(job_name)
    --mail-user=$(mail_user)
    --mail-type=$(mail_type)
    --error=$(logdir)/%j.%x.err
    --output=$(logdir)/%j.%x.out
    --qos=$(qos)
    --nodes=$(nodes)
    --ntasks=$(ntasks)
    --time=$(time)   
    --cpus-per-task=$(cpus_per_task)
    --mem=$(mem_gb)G
    --constraint=cpu
    --wrap $(cmd)
    `
    sleep(5)
    run(submission)
    sleep(5)
    return true
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Submit a batch job to NERSC's SLURM workload manager.

# Arguments
- `job_name`: Identifier for the SLURM job
- `mail_user`: Email address for job notifications
- `mail_type`: Notification type ("ALL", "BEGIN", "END", "FAIL", or "NONE")
- `logdir`: Directory for storing job output/error logs
- `scriptdir`: Directory for storing generated SLURM scripts
- `qos`: Quality of Service level ("regular", "premium", or "preempt")
- `nodes`: Number of nodes to allocate
- `ntasks`: Number of tasks to run
- `time`: Maximum wall time in format "days-HH:MM:SS"
- `cpus_per_task`: CPU cores per task
- `mem_gb`: Memory per node in GB
- `cmd`: Command(s) to execute (String or Vector{String})
- `constraint`: Node type constraint ("cpu" or "gpu")

# Returns
- `true` if job submission succeeds
- `false` if submission fails

# QoS Options
- regular: Standard priority queue
- premium: High priority queue (5x throughput limit)
- preempt: Reduced credit usage but jobs may be interrupted


https://docs.nersc.gov/jobs/policy/
https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes

default is to use shared qos

use
- regular
- preempt (reduced credit usage but not guaranteed to finish)
- premium (priorty runs limited to 5x throughput)

https://docs.nersc.gov/systems/perlmutter/running-jobs/#tips-and-tricks
"""
function nersc_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=mkpath(joinpath(homedir(), "workspace/slurmlogs")),
        scriptdir::String=mkpath(joinpath(homedir(), "workspace/slurm")),
        qos::String="regular",
        nodes::Int=1,
        ntasks::Int=1,
        time::String="2-00:00:00",
        cpus_per_task::Int=Mycelia.NERSC_CPU,
        mem_gb::Int=Mycelia.NERSC_MEM,
        cmd::Union{String, Vector{String}},
        constraint::String="cpu"
    )
        
    # Generate timestamp
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
    
    # Create script filename
    script_name = "$(timestamp)-$(job_name).sh"
    script_path = joinpath(scriptdir, script_name)
    
    # Process commands
    cmd_block = if isa(cmd, String)
        cmd  # Single command as is
    else
        join(cmd, "\n")  # Multiple commands joined with newlines
    end
    
    # Create script content
    script_content = """
    #!/bin/bash
    #SBATCH --job-name=$(job_name)
    #SBATCH --mail-user=$(mail_user)
    #SBATCH --mail-type=$(mail_type)
    #SBATCH --error=$(logdir)/%j.%x.err
    #SBATCH --output=$(logdir)/%j.%x.out
    #SBATCH --qos=$(qos)
    #SBATCH --nodes=$(nodes)
    #SBATCH --ntasks=$(ntasks)
    #SBATCH --time=$(time)
    #SBATCH --cpus-per-task=$(cpus_per_task)
    #SBATCH --mem=$(mem_gb)G
    #SBATCH --constraint=$(constraint)

    $cmd_block
    """
    
    # Write script to file
    write(script_path, script_content)
    
    # Make script executable
    chmod(script_path, 0o755)
    
    # Submit the job
    sleep(5)
    try
        run(`sbatch $script_path`)
    catch e
        @error "Failed to submit job with sbatch: $e"
        return false
    end
    sleep(5)
    
    return true
end
# function nersc_sbatch_regular(;
#         job_name::String,
#         mail_user::String,
#         mail_type::String="ALL",
#         logdir::String=mkpath("$(homedir())/workspace/slurmlogs"),
#         qos::String="regular",
#         nodes::Int=1,
#         ntasks::Int=1,
#         time::String="2-00:00:00",
#         cpus_per_task::Int=Mycelia.NERSC_CPU,
#         mem_gb::Int=Mycelia.NERSC_MEM,
#         cmd::String,
#         constraint::String="cpu"
#     )
#     submission = 
#     `sbatch
#     --job-name=$(job_name)
#     --mail-user=$(mail_user)
#     --mail-type=$(mail_type)
#     --error=$(logdir)/%j.%x.err
#     --output=$(logdir)/%j.%x.out
#     --qos=$(qos)
#     --nodes=$(nodes)
#     --ntasks=$(ntasks)
#     --time=$(time)   
#     --cpus-per-task=$(cpus_per_task)
#     --mem=$(mem_gb)G
#     --constraint=cpu
#     --wrap $(cmd)
#     `
#     sleep(5)
#     run(submission)
#     sleep(5)
#     return true
# end

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

Convert a GFA (Graphical Fragment Assembly) file to FASTA format.

# Arguments
- `gfa::String`: Path to input GFA file
- `fasta::String=gfa * ".fna"`: Path for output FASTA file. Defaults to input filename with ".fna" extension

# Returns
- `String`: Path to the generated FASTA file

# Details
Uses gfatools (via Conda) to perform the conversion. The function will:
1. Ensure gfatools is available in the Conda environment
2. Execute the conversion using gfatools gfa2fa
3. Write sequences to the specified FASTA file
"""
function gfa_to_fasta(;gfa, fasta=gfa * ".fna")
    Mycelia.add_bioconda_env("gfatools")
    p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n gfatools gfatools gfa2fa $(gfa)`, fasta)
    run(p)
    return fasta
    # collect(GraphicalFragmentAssembly.Reader(open(primary_contig_gfa)))
    # open(fasta, "w") do io
    #     fastx_io = FASTX.FASTA.Writer(io)
    #     gfa_graph = Mycelia.parse_gfa(gfa)
    #     for v in Graphs.vertices(gfa_graph)
    #         record = FASTX.FASTA.Record(gfa_graph.vprops[v][:identifier], gfa_graph.vprops[v][:sequence])
    #         write(fastx_io, record)
    #     end
    #     close(fastx_io)
    # end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate per-base genomic coverage from a BAM file using bedtools.

# Arguments
- `bam::String`: Path to input BAM file

# Returns
- `String`: Path to the generated coverage file (`.coverage.txt`)

# Details
Uses bedtools genomecov to compute per-base coverage. Creates a coverage file 
with the format: <chromosome> <position> <coverage_depth>. 
If the coverage file already exists, returns the existing file path.

# Dependencies
Requires bedtools (automatically installed in conda environment)
"""
function determine_fasta_coverage(bam)
    Mycelia.add_bioconda_env("bedtools")
    genome_coverage_file = bam * ".coverage.txt"
    if !isfile(genome_coverage_file)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bedtools bedtools genomecov -d -ibam $(bam)`, genome_coverage_file))
    end
    return genome_coverage_file
end

#outdir="$(homedir())/software/bandage"
# I don't think that this is very portable - assumes sudo and linux
# can make a bandage_jll to fix this longer term
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and installs Bandage, a bioinformatics visualization tool for genome assembly graphs.

# Arguments
- `outdir="/usr/local/bin"`: Target installation directory for the Bandage executable

# Returns
- Path to the installed Bandage executable

# Details
- Downloads Bandage v0.8.1 for Ubuntu
- Installs required system dependencies (libxcb-glx0, libx11-xcb-dev, libfontconfig, libgl1-mesa-glx)
- Attempts installation with sudo, falls back to root if sudo fails
- Skips download if Bandage is already installed at target location

# Dependencies
Requires system commands: wget, unzip, apt
"""
function download_bandage(outdir="/usr/local/bin")
    bandage_executable = joinpath(outdir, "Bandage")
    if !isfile(bandage_executable)
        run(`wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip`)
        run(`unzip Bandage_Ubuntu_static_v0_8_1.zip`)
        isfile("sample_LastGraph") && rm("sample_LastGraph")
        isfile("Bandage_Ubuntu_static_v0_8_1.zip") && rm("Bandage_Ubuntu_static_v0_8_1.zip")
        try # not root
            run(`sudo mv Bandage $(outdir)`)
            run(`sudo apt install libxcb-glx0 libx11-xcb-dev libfontconfig libgl1-mesa-glx -y`)
        catch # root
            run(`mv Bandage $(outdir)`)
            run(`apt install libxcb-glx0 libx11-xcb-dev libfontconfig libgl1-mesa-glx -y`)
        end
    end
    return bandage_executable
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform comprehensive annotation of a FASTA file including gene prediction, protein homology search,
and terminator prediction.

# Arguments
- `fasta::String`: Path to input FASTA file
- `identifier::String`: Unique identifier for output directory (default: FASTA filename without extension)
- `basedir::String`: Base directory for output (default: current working directory)
- `mmseqsdb::String`: Path to MMseqs2 UniRef50 database (default: "~/workspace/mmseqs/UniRef50")
- `threads::Int`: Number of CPU threads to use (default: all available)

# Processing Steps
1. Creates output directory and copies input FASTA
2. Runs Prodigal for gene prediction (nucleotide, amino acid, and GFF output)
3. Performs MMseqs2 homology search against UniRef50
4. Predicts terminators using TransTerm
5. Combines annotations into a unified GFF file
6. Generates GenBank format output

# Returns
- `String`: Path to the output directory containing all generated files

# Files Generated
- `.prodigal.fna`: Predicted genes (nucleotide)
- `.prodigal.faa`: Predicted proteins
- `.prodigal.gff`: Prodigal GFF annotations
- `.gff`: Combined annotations
- `.gff.genbank`: Final GenBank format
"""
function annotate_fasta(;
        fasta,
        identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
        basedir = pwd(),        
        mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
        threads=Sys.CPU_THREADS
    )
    # @show basedir
    outdir = joinpath(basedir, identifier)
    @assert outdir != fasta
    # if !isdir(outdir)
    #     @show isdir(outdir)
    mkpath(outdir)
    f = joinpath(outdir, basename(fasta))
    # make this an rclone copy for portability
    !isfile(f) && cp(fasta, f)
    nucleic_acid_fasta = f * ".prodigal.fna"
    amino_acid_fasta = f * ".prodigal.faa"
    gff_file = f * ".prodigal.gff"
    if !isfile(nucleic_acid_fasta) || !isfile(amino_acid_fasta) || !isfile(gff_file)
        Mycelia.run_prodigal(fasta_file=f)
    end
    mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
    mmseqs_gff_file = Mycelia.write_gff(gff = Mycelia.update_gff_with_mmseqs(gff_file, mmseqs_outfile), outfile = mmseqs_outfile * ".gff")
    transterm_gff_file = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
    joint_gff = Mycelia.write_gff(
        gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file)), ["#seqid", "start", "end"]),
        outfile=f * ".gff")
    annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff, genbank = joint_gff * ".genbank")

    transterm_gff_file_raw_fasta = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
    joint_gff_raw_fasta = Mycelia.write_gff(
        gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file_raw_fasta)), ["#seqid", "start", "end"]),
        outfile=f * ".transterm_raw.gff")
    annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff_raw_fasta, genbank = joint_gff_raw_fasta * ".genbank")
    # else
    #     @info "$(outdir) already present, skipping..."
    # end
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Annotate amino acid sequences in a FASTA file using MMseqs2 search against UniRef50 database.

# Arguments
- `fasta`: Path to input FASTA file containing amino acid sequences
- `identifier`: Name for the output directory (defaults to FASTA filename without extension)
- `basedir`: Base directory for output (defaults to current directory)
- `mmseqsdb`: Path to MMseqs2 formatted UniRef50 database (defaults to ~/workspace/mmseqs/UniRef50)
- `threads`: Number of CPU threads to use (defaults to system thread count)

# Returns
- Path to the output directory containing MMseqs2 search results

The function creates a new directory named by `identifier` under `basedir`, copies the input FASTA file,
and runs MMseqs2 easy-search against the specified database. If the output directory already exists,
the function skips processing and returns the directory path.
"""
function annotate_aa_fasta(;
        fasta,
        identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
        basedir = pwd(),
        mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
        threads=Sys.CPU_THREADS
    )
    # @show basedir
    outdir = joinpath(basedir, identifier)
    @assert outdir != fasta
    if !isdir(outdir)
        @show isdir(outdir)
        mkpath(outdir)
        f = joinpath(outdir, basename(fasta))
        # make this an rclone copy for portability
        cp(fasta, f, force=true)
        amino_acid_fasta = f

        mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
    else
        @info "$(outdir) already present, skipping..."
    end
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

List all directories at the specified rclone path.

# Arguments
- `path::String`: Remote path to list directories from (e.g. "remote:/path/to/dir")

# Returns
- `Vector{String}`: Full paths to all directories found at the specified location
"""
function rclone_list_directories(path)
    directories = [join(split(line)[5:end], " ") for line in eachline(open(`rclone lsd $(path)`))]
    directories = joinpath.(path, directories)
    return directories
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform all-vs-all sequence search using MMseqs2's easy-search command.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequences to compare
- `output::String`: Output directory path (default: input filename + ".mmseqs_easy_search_pairwise")

# Returns
- `String`: Path to the output directory

# Details
Executes MMseqs2 with sensitive search parameters (7 sensitivity steps) and outputs results in 
tabular format with the following columns:
- query, qheader: Query sequence ID and header
- target, theader: Target sequence ID and header  
- pident: Percentage sequence identity
- fident: Fraction of identical matches
- nident: Number of identical matches
- alnlen: Alignment length
- mismatch: Number of mismatches
- gapopen: Number of gap openings
- qstart, qend, qlen: Query sequence coordinates and length
- tstart, tend, tlen: Target sequence coordinates and length
- evalue: Expected value
- bits: Bit score

Requires MMseqs2 to be available through Bioconda.
"""
function mmseqs_pairwise_search(;fasta, output=fasta*".mmseqs_easy_search_pairwise")
    Mycelia.add_bioconda_env("mmseqs2")
    mkpath(output)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-search
        $(fasta)
        $(fasta)
        $(output)/$(basename(fasta)).mmseqs_pairwise_search.txt $(tempdir())
        --format-mode 4
        --format-output query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
        --start-sens 1 -s 7 --sens-steps 7 --sort-results 1 --remove-tmp-files 1 --search-type 3`)
    return output
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function mmseqs_easy_linclust(;fasta, output=fasta*".mmseqs_easy_linclust", tmp=mktempdir())
#     Mycelia.add_bioconda_env("mmseqs2")   
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createdb $(fasta) $(fasta)_DB`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createindex --search-type 3 $(fasta)_DB $(tempdir())`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-linclust $(fasta)_DB $(output) $(tmp)`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta)_DB $(fasta)_DB $(output) $(output).tsv`)
#     return "$(output).tsv"
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Cluster protein or nucleotide sequences using MMseqs2 easy-cluster workflow.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequences to cluster
- `output::String`: Base path for output files (default: input path + ".mmseqs_easy_cluster")
- `tmp::String`: Path to temporary directory (default: auto-generated temp dir)

# Returns
- `String`: Path to the output cluster TSV file containing cluster assignments

# Details
Uses MMseqs2 with minimum sequence identity threshold of 50% (-min-seq-id 0.5) and 
minimum coverage threshold of 80% (-c 0.8). The output TSV file format contains 
tab-separated cluster representative and member sequences.
"""
# --cov-mode: coverage mode (0: coverage of query and target, 1: coverage of target, 2: coverage of query)
function mmseqs_easy_cluster(;fasta, output=fasta*".mmseqs_easy_cluster", tmp=mktempdir())
    outfile = "$(output)_cluster.tsv"
    if !isfile(outfile)
        Mycelia.add_bioconda_env("mmseqs2")
        # at least 50% equivalent
        # --min-seq-id 0.5 -c 0.8
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-cluster $(fasta) $(output) $(tmp)`)
    end
    rm(tmp, recursive=true)
    return "$(output)_cluster.tsv"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Export sequences from a BLAST database to a gzipped FASTA file.

# Arguments
- `path_to_db`: Path to the BLAST database
- `fasta`: Output path for the gzipped FASTA file (default: `path_to_db * ".fna.gz"`)

# Details
Uses conda's BLAST environment to extract sequences using `blastdbcmd`.
The output is automatically compressed using `pigz`.
If the output file already exists, the function will skip extraction.

"""
function export_blast_db(;path_to_db, fasta = path_to_db * ".fna.gz")
    Mycelia.add_bioconda_env("blast")
    if !isfile(fasta)
        # -long_seqids adds GI identifiers - these are cross-referenceable through other means so I'm dropping
        @time run(pipeline(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd  -entry all -outfmt '%f' -db $(path_to_db)`, `pigz`), fasta))
    else
        @info "$(fasta) already present"
    end
end

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
Prints statistics to standard output in tabular format
"""
function fastx_stats(fastq)
    Mycelia.add_bioconda_env("seqkit")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n seqkit seqkit stats $(fastq)`)
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a copy of a file in a temporary directory while preserving the original filename.

# Arguments
- `file_path::String`: Path to the source file to be copied

# Returns
- `String`: Path to the newly created temporary file
"""
function copy_to_tempdir(file_path::String)
    # Create a temporary directory
    temp_dir = mktempdir()

    # Get the file name from the original file path
    file_name = basename(file_path)

    # Create the new file path within the temporary directory
    temp_file_path = joinpath(temp_dir, file_name)

    # Copy the original file to the new path
    cp(file_path, temp_file_path)

    # Return the path of the temporary file
    return temp_file_path
end

# conda install -c bioconda kmer-jellyfish
# count, bc, info, stats, histo, dump, merge, query, cite, mem, jf
# cap at 4 threads, 8Gb per thread by default - this should be plenty fast enough for base usage, but open it up for higher performance!
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers in a FASTA/FASTQ file using Jellyfish.

# Arguments
- `fastx::String`: Path to input FASTA/FASTQ file (can be gzipped)
- `k::Integer`: k-mer length
- `threads::Integer=Sys.CPU_THREADS`: Number of threads to use
- `max_mem::Integer=Int(Sys.free_memory())`: Maximum memory in bytes (defaults to system free memory)
- `canonical::Bool=false`: Whether to count canonical k-mers (both strands combined)
- `outfile::String=auto`: Output filename (auto-generated based on input and parameters)
- `conda_check::Bool=true`: Whether to verify Jellyfish conda installation

# Returns
- `String`: Path to gzipped TSV file containing k-mer counts
"""
function jellyfish_count(;fastx, k, threads=Sys.CPU_THREADS, max_mem=Int(Sys.free_memory()), canonical=false, outfile = ifelse(canonical, "$(fastx).k$(k).canonical.jf", "$(fastx).k$(k).jf"), conda_check=true)
    if conda_check
        Mycelia.add_bioconda_env("kmer-jellyfish")
    end
    mem = Int(floor(max_mem * 0.8))
    
    jellyfish_mem_output = read(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish mem --mer-len $(k) --mem $(mem)`, String)
    # sample output = "68719476736 (68G)\n"
    # this version grabs the exact number at the beginning
    jellyfish_buffer_size = first(split(strip(jellyfish_mem_output)))
    
    # this grabs the human readable version in parentheses
    jellyfish_buffer_size = replace(split(strip(jellyfish_mem_output))[2], r"[\(\)]" => "")
    @show jellyfish_buffer_size

    @info "making a temporary copy of the input fastx"
    temp_fastx = copy_to_tempdir(fastx)
    if occursin(r"\.gz$", temp_fastx)
        run(`gzip -d $(temp_fastx)`)
        temp_fastx = replace(temp_fastx, r"\.gz$" => "")
    end
    @assert isfile(temp_fastx)
    if canonical
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish count --canonical --size $(jellyfish_buffer_size) --threads $(threads) --mer-len $(k) --output $(outfile) $(temp_fastx)`
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish count --size $(jellyfish_buffer_size) --threads $(threads) --mer-len $(k) --output $(outfile) $(temp_fastx)`
    end

    temp_tab = outfile * ".tsv"
    tabular_counts = temp_tab * ".gz"

    if !isfile(tabular_counts)
        if !isfile(outfile)
            @info "running kmer counts"
            @show cmd
            run(cmd)
            @info "done counting kmers"
        end
        if !isfile(temp_tab)
            @info "dumping counts to tab-delimited file"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kmer-jellyfish jellyfish dump --column --tab --output $(temp_tab) $(outfile)`)
            @info "done dumping counts to tab-delimited file"
        end
        run(`gzip $(temp_tab)`)
    end

    isfile(outfile) && rm(outfile)
    isfile(temp_fastx) && rm(temp_fastx)
    return tabular_counts
end


# 7
# 0.169580 seconds (587.93 k allocations: 39.458 MiB, 76.55% compilation time)
# 11
# 4.420040 seconds (194.49 k allocations: 11.566 MiB)
# 13
# 20.690993 seconds (521.75 k allocations: 37.253 MiB, 0.58% compilation time)
# 17
# 412.670050 seconds (529.86 k allocations: 37.007 MiB, 0.03% compilation time)
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Jellyfish k-mer count file into a frequency histogram.

# Arguments
- `jellyfish_counts_file::String`: Path to the gzipped TSV file containing Jellyfish k-mer counts
- `outfile::String=replace(jellyfish_counts_file, r"\\.tsv\\.gz\$" => ".count_histogram.tsv")`: Optional output file path

# Returns
- `String`: Path to the generated histogram file

# Description
Processes a Jellyfish k-mer count file to create a frequency histogram where:
- Column 1: Number of k-mers that share the same count
- Column 2: The count they share

Uses system sorting with LC_ALL=C for optimal performance on large files.

Notes
- Requires gzip, sort, uniq, and sed command line tools
- Uses intermediate disk storage for sorting large files
- Skips processing if output file already exists
"""
function jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts_file, outfile=replace(jellyfish_counts_file, r"\.tsv\.gz$" => ".count_histogram.tsv"))
    # sorting with LC_ALL=C is the biggest speed up here of anything I've found
    if !isfile(outfile)
        io = open(pipeline(
                `gzip -dc $(jellyfish_counts_file)`,
                `cut -f2`,
                Cmd(`sort --temporary-directory . --compress-program gzip --numeric --stable`, env=Dict("LC_ALL" => "C")),
                `uniq --count`,
                `sed 's/^ *//'`,
                `sed 's/ /\t/'`
                ))
        frequency_histogram_table = CSV.read(io, DataFrames.DataFrame, header=["number of kmers", "number of observations"], delim='\t')
        CSV.write(outfile, frequency_histogram_table, delim='\t')
    else
        @info "$(outfile) already exists"
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load k-mer counts from a Jellyfish output file into a DataFrame.

# Arguments
- `jellyfish_counts::String`: Path to a gzipped TSV file (*.jf.tsv.gz) containing Jellyfish k-mer counts

# Returns
- `DataFrame`: Table with columns:
  - `kmer`: Biologically encoded k-mers as `DNAKmer{k}` objects
  - `count`: Integer count of each k-mer's occurrences

# Notes
- Input file must be a gzipped TSV with exactly two columns (k-mer sequences and counts)
- K-mer length is automatically detected from the first entry
- Filename must end with '.jf.tsv.gz'
"""
function load_jellyfish_counts(jellyfish_counts)
    @assert occursin(r"\.jf\.tsv\.gz$", jellyfish_counts)
    open(jellyfish_counts) do io
        table = CSV.read(CodecZlib.GzipDecompressorStream(io), DataFrames.DataFrame, delim="\t", header=["kmer", "count"])
        k = length(table[1, "kmer"])
        table[!, "kmer"] = map(x -> Kmers.DNAKmer{k}(String(x)), collect(table[!, "kmer"]))
        return table
    end
end

# Usage: jellyfish merge [options] input:string+

# Merge jellyfish databases

# Options (default value in (), *required):
#  -o, --output=string                      Output file (mer_counts_merged.jf)
#  -m, --min                                Compute min count instead of sum (false)
#  -M, --max                                Compute max count instead of sum (false)
#  -j, --jaccard                            Compute the jaccard and weighted jaccard similarities (false)
#  -L, --lower-count=uint64                 Don't output k-mer with count < lower-count
#  -U, --upper-count=uint64                 Don't output k-mer with count > upper-count
#      --usage                              Usage
#  -h, --help                               This message
#  -V, --version                            Version

# Usage: jellyfish query [options] file:path mers:string+

# Query a Jellyfish database

# Options (default value in (), *required):
#  -s, --sequence=path                      Output counts for all mers in sequence
#  -o, --output=path                        Output file (stdout)
#  -i, --interactive                        Interactive, queries from stdin (false)
#  -l, --load                               Force pre-loading of database file into memory (false)
#  -L, --no-load                            Disable pre-loading of database file into memory (false)
#  -U, --usage                              Usage
#  -h, --help                               This message
#  -V, --version                            Version


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add random noise to create a vector of jittered values.

Generates `n` values by adding random noise to the input value `x`. 
The noise is uniformly distributed between -1/3 and 1/3.

# Arguments
- `x`: Base value to add jitter to
- `n`: Number of jittered values to generate

# Returns
- Vector of length `n` containing jittered values around `x`
"""
function jitter(x, n)
    return [x + rand() / 3 * (ifelse(rand(Bool), 1, -1)) for i in 1:n]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Counts k-mer occurrences in a FASTA file, considering both forward and reverse complement sequences.

# Arguments
- `kmer_type`: Type specification for k-mers (e.g., `DNAKmer{21}`)
- `fasta`: Path to FASTA file containing reference sequences

# Returns
- `Dict{kmer_type, Int}`: Dictionary mapping each k-mer to its total count across all sequences
"""
function fasta_to_reference_kmer_counts(;kmer_type, fasta)
    kmer_counts = Dict{kmer_type, Int}()
    for record in Mycelia.open_fastx(fasta)
        record_sequence = BioSequences.LongDNA{2}(FASTX.sequence(record))
        forward_counts = StatsBase.countmap(kmer for (i, kmer) in Kmers.EveryKmer{kmer_type}(record_sequence))
        reverse_counts = StatsBase.countmap(kmer for (i, kmer) in Kmers.EveryKmer{kmer_type}(BioSequences.reverse_complement(record_sequence)))
        record_counts = merge(+, forward_counts, reverse_counts)
        merge!(+, kmer_counts, record_counts)
    end
    return kmer_counts
end



# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# My standard pacbio aligning and sorting. No filtering done in this step.

# Use shell_only=true to get string command to submit to SLURM
# """
# function map_pacbio_reads(;
#         fastq,
#         reference_fasta,
#         temp_sam_outfile = fastq * "." * basename(reference_fasta) * "." * "minimap2.sam",
#         outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz"),
#         threads = Sys.CPU_THREADS,
#         memory = Sys.total_memory(),
#         shell_only = false
#     )
#     # 4G is the default
#     # smaller, higher diversity databases do better with 5+ as the denominator - w/ <=4 they run out of memory
#     index_chunk_size = "$(Int(floor(memory/5e9)))G"
#     @show index_chunk_size
#     @show threads
#     Mycelia.add_bioconda_env("minimap2")
#     # Mycelia.add_bioconda_env("samtools")
#     Mycelia.add_bioconda_env("pigz")
#     if shell_only
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
#         """
#         return cmd
#     else
#         if !isfile(outfile)
#             map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#             run(map)
#             run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`)
#             @assert isfile(outfile)
#         else
#             @info "$(outfile) already present"
#         end
#     end
# end

# """
# My standard pacbio aligning and sorting. No filtering done in this step.

# Use shell_only=true to get string command to submit to SLURM
# """
# function minimap_index_pacbio(;
#         reference_fasta,
#         outfile = replace(reference_fasta, Mycelia.FASTA_REGEX => ".pacbio.mmi"),
#         threads = Sys.CPU_THREADS,
#         shell_only = false
#     )
#     Mycelia.add_bioconda_env("minimap2")
#     Mycelia.add_bioconda_env("samtools")
#     if shell_only
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-pb $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile) \\
#         | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile) \\
#         && rm $(temp_sam_outfile)
#         """
#         return cmd
#     else
#         if !isfile(outfile)
#             map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-pb $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#             run(map)
#             p = pipeline(
#                 `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile)`,
#                 `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile)`
#             )
#             run(p)
#             rm(temp_sam_outfile)
#         else
#             @info "$(outfile) already present"
#         end
#     end
# end

# function filter_short_reads()
# end

# function map_short_reads()
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the current date and time as a normalized string with all non-word characters removed.

The output format is based on ISO datetime (YYYYMMDDThhmmss) but strips any special characters
like hyphens, colons or dots.
"""
function normalized_current_datetime()
    return replace(Dates.format(Dates.now(), Dates.ISODateTimeFormat), r"[^\w]" => "")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the current git commit hash of the repository.

# Arguments
- `short::Bool=false`: If true, returns abbreviated 8-character hash

# Returns
A string containing the git commit hash (full 40 characters by default)
"""
function githash(;short=false)
    git_hash = rstrip(read(`git rev-parse HEAD`, String))
    if short
        git_hash = git_hash[1:8]
    end
    return git_hash
end

# CSV is too memory inefficient, the others too slow :(
# # using uCSV
# # k=11
# # 3.444974 seconds (24.58 M allocations: 1.374 GiB, 34.65% gc time, 16.90% compilation time)
# # k=13
# # 362.285866 seconds (357.11 M allocations: 20.550 GiB, 91.60% gc time)

# # using DelimitedFiles.readdlm
# # k=11
# # 2.386620 seconds (16.11 M allocations: 632.732 MiB, 34.16% gc time, 24.25% compilation time)
# # k=13
# # 82.888552 seconds (227.49 M allocations: 8.766 GiB, 82.01% gc time)

# # CSV
# # k=11
# # 12.328422 seconds (7.62 M allocations: 732.639 MiB, 19091.67% compilation time: <1% of which was recompilation)
# # k=13
# # 37.098948 seconds (89.38 k allocations: 2.354 GiB, 93.56% gc time)

# function parse_jellyfish_counts(tabular_counts)
#     # load in the data
#     @assert occursin(r"\.gz$", tabular_counts) "this expects gzipped jellyfish tabular counts"
#     io = CodecZlib.GzipDecompressorStream(open(tabular_counts))
#     canonical_kmer_counts_table = DataFrames.DataFrame(CSV.File(io; delim='\t', header=false))
#     DataFrames.rename!(canonical_kmer_counts_table, [:Column1 => :kmer, :Column2 => :count])
    
#     # recode the kmers from strings to fixed sized kmer types
#     unique_kmer_lengths = unique(length.(canonical_kmer_counts_table[!, "kmer"]))
#     @assert length(unique_kmer_lengths) == 1
#     k = first(unique_kmer_lengths)
#     canonical_kmer_counts_table[!, "kmer"] = Kmers.DNAKmer{k}.(canonical_kmer_counts_table[!, "kmer"])
    
#     return canonical_kmer_counts_table
# end


# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/taxonomy/
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Retrieve taxonomic information for a given NCBI taxonomy ID.

# Arguments
- `taxa_id`: NCBI taxonomy identifier (integer)

# Returns
- `DataFrame`: Taxonomy summary containing fields like tax_id, rank, species, etc.
"""
function ncbi_taxon_summary(taxa_id)
    Mycelia.add_bioconda_env("ncbi-datasets")
    p = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets datasets summary taxonomy taxon $(taxa_id) --as-json-lines`,
        `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets dataformat tsv taxonomy --template tax-summary`
        )
    return DataFrames.DataFrame(uCSV.read(open(p), delim='\t', header=1))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find the closest prime number to the given integer `n`.

Returns the nearest prime number to `n`. If two prime numbers are equally distant 
from `n`, returns the smaller one.

# Arguments
- `n::Int`: The input integer to find the nearest prime for

# Returns
- `Int`: The closest prime number to `n`
"""
function nearest_prime(n::Int)
    if n < 2
        return 2
    end
    next_p = Primes.nextprime(n)
    prev_p = Primes.prevprime(n)
    if n - prev_p <= next_p - n
        return prev_p
    else
        return next_p
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a sequence of Fibonacci numbers strictly less than the input value.

# Arguments
- `n::Int`: Upper bound (exclusive) for the Fibonacci sequence

# Returns
- `Vector{Int}`: Array containing Fibonacci numbers less than n
"""
function fibonacci_numbers_less_than(n::Int)
    if n <= 0
        return []
    elseif n == 1
        return [0]
    else
        fib = [0, 1]
        next_fib = fib[end] + fib[end-1]
        while next_fib < n
            push!(fib, next_fib)
            next_fib = fib[end] + fib[end-1]
        end
        return fib
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generates a specialized sequence of prime numbers combining:
- Odd primes up to 23 (flip_point)
- Primes nearest to Fibonacci numbers above 23 up to max

# Arguments
- `min::Int=0`: Lower bound for the sequence
- `max::Int=10_000`: Upper bound for the sequence

# Returns
Vector of Int containing the specialized prime sequence
"""
function ks(;min=0, max=10_000)
    # flip from all odd primes to only nearest to fibonnaci primes
    flip_point = 23
    vcat(
        filter(isodd, Primes.primes(min, flip_point)),
        filter(x -> x > flip_point, nearest_prime.(fibonacci_numbers_less_than(max)))
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Copy files between local and remote storage using rclone with automated retry logic.

# Arguments
- `source::String`: Source path or remote (e.g. "local/path" or "gdrive:folder")
- `dest::String`: Destination path or remote (e.g. "gdrive:folder" or "local/path")

# Keywords
- `config::String=""`: Optional path to rclone config file
- `max_attempts::Int=3`: Maximum number of retry attempts
- `sleep_timer::Int=60`: Initial sleep duration between retries in seconds (doubles after each attempt)

# Details
Uses optimized rclone settings for large files:
- 2GB chunk size
- 1TB upload cutoff
- Rate limited to 1 transaction per second
"""
function rclone_copy(source, dest; config="", max_attempts=3, sleep_timer=60)
    done = false
    attempts = 0
    while !done && attempts < max_attempts
        attempts += 1
        try
            # https://forum.rclone.org/t/google-drive-uploads-failing-http-429/34147/9
            # --tpslimit                                       Limit HTTP transactions per second to this
            # --drive-chunk-size SizeSuffix                    Upload chunk size (default 8Mi)
            # --drive-upload-cutoff SizeSuffix                 Cutoff for switching to chunked upload (default 8Mi)
            # not currently using these but they may become helpful
            # --drive-pacer-burst int                          Number of API calls to allow without sleeping (default 100)
            # --drive-pacer-min-sleep Duration                 Minimum time to sleep between API calls (default 100ms)
            if isempty(config)
                cmd = `rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $(source) $(dest)`
            else
                cmd = `rclone --config $(config) copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $(source) $(dest)`
            end
            @info "copying $(source) to $(dest) with command: $(cmd)"
            run(cmd)
            done = true
        catch
            @info "copying incomplete, sleeping $(sleep_timer) seconds and trying again..."
            sleep(sleep_timer)
            sleep_timer *= 2
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify all columns that have only missing or empty values

Returns as a bit array

See also: drop_empty_columns, drop_empty_columns!
"""
function find_nonempty_columns(df)
    non_empty_columns = [eltype(col) != Missing || !all(v -> isnothing(v) || ismissing(v) || (!isa(v, Date) && isempty(v)), col) for col in DataFrames.eachcol(df)]
    return non_empty_columns
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify all columns that have only missing or empty values, and remove those columns from the dataframe.

Returns a modified copy of the dataframe.

See also: drop_empty_columns!
"""
function drop_empty_columns(df::DataFrames.AbstractDataFrame)
    # Filter the DataFrame columns by checking if not all values in the column are missing or empty
    non_empty_columns = find_nonempty_columns(df)
    filtered_df = df[:, non_empty_columns]
    return filtered_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identify all columns that have only missing or empty values, and remove those columns from the dataframe *in-place*.

Returns a modified version of the original dataframe. 

See also: drop_empty_columns
"""
function drop_empty_columns!(df::DataFrames.AbstractDataFrame)
    # Filter the DataFrame columns by checking if not all values in the column are missing or empty
    non_empty_columns = find_nonempty_columns(df)
    # df = df[!, non_empty_columns]
    DataFrames.select!(df, non_empty_columns)
    return df
end

# # Need to add hashdeep & logging
# function tarchive(;directory, tarchive=directory * ".tar.gz")
#     run(`tar --create --gzip --verbose --file=$(tarchive) $(directory)`)
#     return tarchive
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensures the hashdeep utility is installed on the system.

Checks if hashdeep is available in PATH and attempts to install it via apt package manager
if not found. Will try with sudo privileges first, then without sudo if that fails.

# Details
- Checks PATH for existing hashdeep executable
- Attempts installation using apt package manager
- Requires a Debian-based Linux distribution

# Returns
- Nothing, but prints status messages during execution
"""
function install_hashdeep()
    if Sys.which("hashdeep") !== nothing
        println("hashdeep executable found")
    else
        println("hashdeep executable not found in PATH, installing")
        try
            run(`sudo apt install hashdeep -y`)
        catch
            run(`apt install hashdeep -y`)
        end
    end
end

# Need to add hashdeep & logging
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a gzipped tar archive of the specified directory along with verification files.

# Arguments
- `directory`: Source directory path to archive
- `tarchive`: Optional output archive path (defaults to directory name with .tar.gz extension)

# Generated Files
- `{tarchive}`: The compressed tar archive
- `{tarchive}.log`: Contents listing of the archive
- `{tarchive}.hashdeep.dfxml`: Cryptographic hashes (MD5, SHA1, SHA256) of the archive

# Returns
- Path to the created tar archive file
"""
function create_tarchive(;directory, tarchive=directory * ".tar.gz")
    directory = normpath(directory)
    tarchive = normpath(tarchive)
    output_dir = mkpath(dirname(tarchive))
    
    install_hashdeep()


    working_dir, source = splitdir(directory)
    if isempty(source)
        source = working_dir
        working_dir = pwd()
    end
    if isempty(working_dir)
        working_dir = pwd()
    end
    println("output_dir: $output_dir\n" *
            "working_dir: $working_dir\n" *
            "source: $source\n" *
            "target_file: $tarchive")
    log_file = tarchive * ".log"
    hashdeep_file = tarchive * ".hashdeep.dfxml"
    
    if !isfile(tarchive)
        run(`tar --create --gzip --verbose --file=$(tarchive) $(directory)`)
    end
    if !isfile(log_file)
        run(pipeline(`tar -tvf $(tarchive)`,log_file))
    end
    if !isfile(hashdeep_file)
        run(pipeline(`hashdeep -c md5,sha1,sha256 -b -d $(tarchive)`,hashdeep_file))
    end
    return tarchive
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract contents of a gzipped tar archive file to a specified directory.

# Arguments
- `tarchive::AbstractString`: Path to the .tar.gz file to extract
- `directory::AbstractString=dirname(tarchive)`: Target directory for extraction (defaults to the archive's directory)

# Returns
- `AbstractString`: Path to the directory where contents were extracted
"""
function tar_extract(;tarchive, directory=dirname(tarchive))
    run(`tar --extract --gzip --verbose --file=$(tarchive) --directory=$(directory)`)
    return directory
end

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

Identifies sequence regions that require resampling based on kmer solidity patterns.

# Arguments
- `record_kmer_solidity::BitVector`: Boolean array where `true` indicates solid kmers
- `solid_branching_kmer_indices::Vector{Int}`: Indices of solid branching kmers

# Returns
- `Vector{UnitRange{Int64}}`: Array of ranges (start:stop) indicating stretches that need resampling

# Details
Finds continuous stretches of non-solid kmers and extends them to the nearest solid branching
kmers on either side. These stretches represent regions that need resampling.

If a stretch doesn't have solid branching kmers on both sides, it is excluded from the result.
Duplicate ranges are removed from the final output.
"""
function find_resampling_stretches(;record_kmer_solidity, solid_branching_kmer_indices)
    indices = findall(.!record_kmer_solidity)  # Find the indices of false values
    if isempty(indices)
        return UnitRange{Int64}[]
    end
    
    diffs = diff(indices)  # Calculate the differences between consecutive indices
    # @show diffs
    range_starts = [indices[1]]  # Start with the first false index
    range_ends = Int[]
    
    for (i, d) in enumerate(diffs)
        if d > 1
            push!(range_ends, indices[i])
            push!(range_starts, indices[i+1])
        end
    end
    
    push!(range_ends, indices[end])  # Add the last false index as a range end
    
    low_quality_runs = [(start, stop) for (start, stop) in zip(range_starts, range_ends)]
    
    resampling_stretches = UnitRange{Int64}[]
    
    for low_quality_run in low_quality_runs
        unders = filter(solid_branching_kmer -> solid_branching_kmer < first(low_quality_run), solid_branching_kmer_indices)
        overs = filter(solid_branching_kmer -> solid_branching_kmer > last(low_quality_run), solid_branching_kmer_indices)
        if isempty(overs) || isempty(unders)
            continue
        else
            nearest_under = maximum(unders)
            nearest_over = minimum(overs)
            push!(resampling_stretches, nearest_under:nearest_over)
        end
    end
    if !allunique(resampling_stretches)
        resampling_stretches = unique!(resampling_stretches)
    end
    return resampling_stretches
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

Process and error-correct a FASTQ sequence record using a kmer graph and path resampling.

# Arguments
- `record`: FASTQ record containing the sequence to process
- `kmer_graph`: MetaGraph containing the kmer network and associated properties
- `yen_k_shortest_paths_and_weights`: Cache of pre-computed k-shortest paths between nodes
- `yen_k`: Number of alternative paths to consider during resampling (default: 3)

# Description
Performs error correction by:
1. Trimming low-quality sequence ends
2. Identifying stretches requiring resampling between solid branching kmers
3. Selecting alternative paths through the kmer graph based on:
   - Path quality scores
   - Transition likelihoods
   - Path length similarity to original sequence

# Returns
- Modified FASTQ record with error-corrected sequence and updated quality scores
- Original record if no error correction was needed

# Required Graph Properties
The kmer_graph must contain the following properties:
- :ordered_kmers
- :likely_valid_kmer_indices  
- :kmer_indices
- :branching_nodes
- :assembly_k
- :transition_likelihoods
- :kmer_mean_quality
- :kmer_total_quality
"""
function process_fastq_record(;record, kmer_graph, yen_k_shortest_paths_and_weights, yen_k=3)
    ordered_kmers = MetaGraphs.get_prop(kmer_graph, :ordered_kmers)
    likely_valid_kmers = Set(ordered_kmers[MetaGraphs.get_prop(kmer_graph, :likely_valid_kmer_indices)])
    kmer_to_index_map = MetaGraphs.get_prop(kmer_graph, :kmer_indices)
    branching_nodes_set = MetaGraphs.get_prop(kmer_graph, :branching_nodes)
    assembly_k = MetaGraphs.get_prop(kmer_graph, :assembly_k)
    transition_likelihoods = MetaGraphs.get_prop(kmer_graph, :transition_likelihoods)
    kmer_mean_quality = MetaGraphs.get_prop(kmer_graph, :kmer_mean_quality)
    kmer_total_quality = MetaGraphs.get_prop(kmer_graph, :kmer_total_quality)
    
    new_record_identifier = FASTX.identifier(record) * ".k$(assembly_k)"
    record_sequence = BioSequences.LongDNA{4}(FASTX.sequence(record))

    kmer_type = Kmers.DNAKmer{assembly_k}
    record_kmers = last.(collect(Kmers.EveryKmer{kmer_type}(record_sequence)))
    record_quality_scores = collect(FASTX.quality_scores(record))
    record_kmer_quality_scores = [record_quality_scores[i:i+assembly_k-1] for i in 1:length(record_quality_scores)-assembly_k+1]
    
    record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
    record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
    record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    
    # trim beginning of fastq
    initial_solid_kmer = findfirst(record_kmer_solidity)
    if isnothing(initial_solid_kmer)
        return record
    elseif initial_solid_kmer > 1
        record_kmers = record_kmers[initial_solid_kmer:end]
        record_kmer_quality_scores = record_kmer_quality_scores[initial_solid_kmer:end]
        record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
        record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
        record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    end
    initial_solid_kmer = 1
    
    # trim end of fastq
    last_solid_kmer = findlast(record_kmer_solidity)
    if last_solid_kmer != length(record_kmer_solidity)
        record_kmers = record_kmers[1:last_solid_kmer]
        record_kmer_quality_scores = record_kmer_quality_scores[1:last_solid_kmer]
        record_kmer_solidity = map(kmer -> kmer in likely_valid_kmers, record_kmers)
        record_branching_kmers = [kmer_to_index_map[kmer] in branching_nodes_set for kmer in record_kmers]
        record_solid_branching_kmers = record_kmer_solidity .& record_branching_kmers
    end
    
    # identify low quality runs and the solid branchpoints we will use for resampling
    solid_branching_kmer_indices = findall(record_solid_branching_kmers)
    resampling_stretches = find_resampling_stretches(;record_kmer_solidity, solid_branching_kmer_indices)

    # nothing to do
    if isempty(resampling_stretches)
        return record
    end
    trusted_range = 1:max(first(first(resampling_stretches))-1, 1)
    
    new_record_kmers = record_kmers[trusted_range]
    new_record_kmer_qualities = record_kmer_quality_scores[trusted_range]
    
    
    for (i, resampling_stretch) in enumerate(resampling_stretches)
        starting_solid_kmer = record_kmers[first(resampling_stretch)]
        ending_solid_kmer = record_kmers[last(resampling_stretch)]
        
        current_quality_scores = record_quality_scores[resampling_stretch]
        u = kmer_to_index_map[starting_solid_kmer]
        v = kmer_to_index_map[ending_solid_kmer]
        if !haskey(yen_k_shortest_paths_and_weights, u => v)
            yen_k_result = Graphs.yen_k_shortest_paths(kmer_graph, u, v, Graphs.weights(kmer_graph), yen_k)
            yen_k_shortest_paths_and_weights[u => v] = Vector{Pair{Vector{Int}, Float64}}()
            for path in yen_k_result.paths
                path_weight = Statistics.mean([kmer_total_quality[ordered_kmers[node]] for node in path])
                path_transition_likelihoods = 1.0
                for (a, b) in zip(path[1:end-1], path[2:end])
                    path_transition_likelihoods *= transition_likelihoods[a, b]
                end
                joint_weight = path_weight * path_transition_likelihoods
                push!(yen_k_shortest_paths_and_weights[u => v], path => joint_weight)
            end
        end
        yen_k_path_weights = yen_k_shortest_paths_and_weights[u => v]      
        if length(yen_k_path_weights) > 1
            current_distance = length(resampling_stretch)
            initial_weights = last.(yen_k_path_weights)
            path_lengths = length.(first.(yen_k_path_weights))
            deltas = map(l -> abs(l-current_distance), path_lengths)
            adjusted_weights = initial_weights .* map(d -> exp(-d * log(2)), deltas)
            # make it more severe?
            # adjusted_weights = adjusted_weights.^2
            
            # and a bonus for usually being correct
            
            selected_path_index = StatsBase.sample(StatsBase.weights(adjusted_weights))
            selected_path, selected_path_weights = yen_k_path_weights[selected_path_index]
            selected_path_kmers = [ordered_kmers[kmer_index] for kmer_index in selected_path]
            
            if last(new_record_kmers) == first(selected_path_kmers)
                selected_path_kmers = selected_path_kmers[2:end]
            end
            append!(new_record_kmers, selected_path_kmers)
            selected_kmer_qualities = [Int8.(min.(typemax(Int8), floor.(kmer_mean_quality[kmer]))) for kmer in selected_path_kmers]
            append!(new_record_kmer_qualities, selected_kmer_qualities)
        else
            selected_path_kmers = record_kmers[resampling_stretch]
            if last(new_record_kmers) == first(selected_path_kmers)
                selected_path_kmers = selected_path_kmers[2:end]
            end
            append!(new_record_kmers, selected_path_kmers)
            selected_kmer_qualities = [Int8.(min.(typemax(Int8), floor.(kmer_mean_quality[kmer]))) for kmer in selected_path_kmers]
            append!(new_record_kmer_qualities, selected_kmer_qualities)
        end
        if i < length(resampling_stretches) # append high quality gap
            next_solid_start = last(resampling_stretch)+1
            next_resampling_stretch = resampling_stretches[i+1]
            next_solid_stop = first(next_resampling_stretch)-1
            if !isempty(next_solid_start:next_solid_stop)
                selected_path_kmers = record_kmers[next_solid_start:next_solid_stop]
                append!(new_record_kmers, selected_path_kmers)
                selected_kmer_qualities = record_kmer_quality_scores[next_solid_start:next_solid_stop]
                append!(new_record_kmer_qualities, selected_kmer_qualities)
            end
        else # append remainder of sequence
            @assert i == length(resampling_stretches)
            next_solid_start = last(resampling_stretch)+1
            if next_solid_start < length(record_kmers)
                selected_path_kmers = record_kmers[next_solid_start:end]
                append!(new_record_kmers, selected_path_kmers)
                selected_kmer_qualities = record_kmer_quality_scores[next_solid_start:end]
                append!(new_record_kmer_qualities, selected_kmer_qualities)
            end
        end
    end
    
    for (a, b) in zip(new_record_kmers[1:end-1], new_record_kmers[2:end])
        @assert a != b
    end
    new_record_sequence = Mycelia.kmer_path_to_sequence(new_record_kmers)
    new_record_quality_scores = new_record_kmer_qualities[1]
    for new_record_kmer_quality in new_record_kmer_qualities[2:end]
        push!(new_record_quality_scores, last(new_record_kmer_quality))
    end
    new_record = fastq_record(identifier=new_record_identifier, sequence=new_record_sequence, quality_scores=new_record_quality_scores)
    return new_record
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Polish FASTQ reads using a k-mer graph-based approach to correct potential sequencing errors.

# Arguments
- `fastq::String`: Path to input FASTQ file
- `k::Int=1`: Initial k-mer size parameter. Final assembly k-mer size may differ.

# Process
1. Builds a directed k-mer graph from input reads
2. Processes each read through the graph to find optimal paths
3. Writes corrected reads to a new FASTQ file
4. Automatically compresses output with gzip

# Returns
Named tuple with:
- `fastq::String`: Path to output gzipped FASTQ file
- `k::Int`: Final assembly k-mer size used
"""
function polish_fastq(;fastq, k=1)
    kmer_graph = build_directed_kmer_graph(fastq=fastq, k=k)
    assembly_k = MetaGraphs.get_prop(kmer_graph, :assembly_k)
    @info "polishing with k = $(assembly_k)"
    revised_records = []
    yen_k_shortest_paths_and_weights = Dict{Pair{Int, Int}, Vector{Pair{Vector{Int}, Float64}}}()
    ProgressMeter.@showprogress for record in collect(Mycelia.open_fastx(fastq))
        revised_record = process_fastq_record(;record, kmer_graph, yen_k_shortest_paths_and_weights)
        push!(revised_records, revised_record)
    end
    
    fastq_out = replace(fastq, Mycelia.FASTQ_REGEX => ".k$(assembly_k).fq")
    open(fastq_out, "w") do io
        fastx_io = FASTX.FASTQ.Writer(io)
        for record in revised_records
            write(fastx_io, record)
        end
        close(fastx_io)
    end
    run(`gzip --force $(fastq_out)`)
    return (fastq = fastq_out * ".gz", k=assembly_k)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a directed graph representation of k-mer transitions from FASTQ sequencing data.

# Arguments
- `fastq`: Path to input FASTQ file
- `k`: K-mer size (default: 1). Must be odd and prime. If k=1, optimal size is auto-determined
- `plot`: Boolean to display quality distribution plot (default: false)

# Returns
MetaDiGraph with properties:
- assembly_k: k-mer size used
- kmer_counts: frequency of each k-mer
- transition_likelihoods: edge weights between k-mers
- kmer_mean_quality, kmer_total_quality: quality metrics
- branching_nodes, unbranching_nodes: topological classification
- likely_valid_kmer_indices: k-mers above mean quality threshold
- likely_sequencing_artifact_indices: potential erroneous k-mers

# Note
For DNA assembly, quality scores are normalized across both strands.
"""
function build_directed_kmer_graph(;fastq, k=1, plot=false)
    if k == 1
        assembly_k = Mycelia.assess_dnamer_saturation([fastq])
    else
        @assert isodd(k)
        @assert Primes.isprime(k)
        assembly_k = k
    end
    kmer_type = Kmers.DNAKmer{assembly_k}

    # initializing the graph with kmer counts
    kmer_counts = Mycelia.count_kmers(kmer_type, fastq)
    ordered_kmers = collect(keys(kmer_counts))
    total_states = length(ordered_kmers)
    graph = MetaGraphs.MetaDiGraph(total_states)
    MetaGraphs.set_prop!(graph, :assembly_k, assembly_k)
    MetaGraphs.set_prop!(graph, :kmer_counts, kmer_counts)
    MetaGraphs.set_prop!(graph, :total_states, total_states)
    MetaGraphs.set_prop!(graph, :ordered_kmers, ordered_kmers)
    kmer_indices = sort(Dict(kmer => i for (i, kmer) in enumerate(keys(kmer_counts))))
    MetaGraphs.set_prop!(graph, :kmer_indices, kmer_indices)
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, fastq)
    MetaGraphs.set_prop!(graph, :canonical_kmer_counts, canonical_kmer_counts)
    canonical_kmer_indices = sort(Dict(kmer => i for (i, kmer) in enumerate(keys(canonical_kmer_counts))))
    MetaGraphs.set_prop!(graph, :canonical_kmer_indices, canonical_kmer_indices)
    
    
    # kmer quality and likelihoods
    
    records = collect(Mycelia.open_fastx(fastq))
    read_quality_scores = [collect(FASTX.quality_scores(record)) for record in records]
    all_kmer_quality_support = Dict{kmer_type, Vector{Float64}}()
    for record in records
        record_quality_scores = collect(FASTX.quality_scores(record))
        record_quality_score_slices = [record_quality_scores[i:i+assembly_k-1] for i in 1:length(record_quality_scores)-assembly_k+1]
        sequence = BioSequences.LongDNA{2}(FASTX.sequence(record))
        for ((i, kmer), kmer_base_qualities) in zip(Kmers.EveryKmer{kmer_type}(sequence), record_quality_score_slices)
            if haskey(all_kmer_quality_support, kmer)
                all_kmer_quality_support[kmer] = all_kmer_quality_support[kmer] .+ kmer_base_qualities
            else
                all_kmer_quality_support[kmer] = kmer_base_qualities
            end
        end
    end
    
    # strand normalization shares observational quality across strands - only relevant for non-stranded DNA genome assembly
    strand_normalized_quality_support = Dict{kmer_type, Vector{Float64}}()
    for (kmer, support) in all_kmer_quality_support
        strand_normalized_quality_support[kmer] = support
        if haskey(all_kmer_quality_support, BioSequences.reverse_complement(kmer))
            strand_normalized_quality_support[kmer] .+= all_kmer_quality_support[BioSequences.reverse_complement(kmer)]
        end
    end
    strand_normalized_quality_support
    kmer_mean_quality = sort(Dict(kmer => strand_normalized_quality_support[kmer] ./ canonical_kmer_counts[BioSequences.canonical(kmer)] for kmer in ordered_kmers))
    MetaGraphs.set_prop!(graph, :kmer_mean_quality, kmer_mean_quality)
    kmer_total_quality = sort(Dict(kmer => sum(quality_values) for (kmer, quality_values) in strand_normalized_quality_support))
    MetaGraphs.set_prop!(graph, :kmer_total_quality, kmer_total_quality)
    state_likelihoods = sort(Dict(kmer => total_quality / sum(values(kmer_total_quality)) for (kmer, total_quality) in kmer_total_quality))
    MetaGraphs.set_prop!(graph, :state_likelihoods, state_likelihoods)


    # all transition likelihood calculation
    transition_likelihoods = SparseArrays.spzeros(total_states, total_states)
    for record in records
        sequence = BioSequences.LongDNA{4}(FASTX.sequence(record))
        sources = Kmers.EveryKmer{kmer_type}(sequence[1:end-1])
        destinations = Kmers.EveryKmer{kmer_type}(sequence[2:end])
        for ((source_i, source), (destination_i, destination)) in zip(sources, destinations)
            source_index = kmer_indices[source]
            destination_index = kmer_indices[destination]
            transition_likelihoods[source_index, destination_index] += 1
        end
    end
    for source in 1:total_states
        outgoing_transition_counts = transition_likelihoods[source, :]
        if sum(outgoing_transition_counts) > 0
            transition_likelihoods[source, :] .= transition_likelihoods[source, :] ./ sum(transition_likelihoods[source, :]) 
        end
    end
    row_indices, column_indices, cell_values = SparseArrays.findnz(transition_likelihoods)
    for (row, col, value) in zip(row_indices, column_indices, cell_values)
        Graphs.add_edge!(graph, row, col)
        MetaGraphs.set_prop!(graph, row, col, :transition_likelihood, value)
    end
    MetaGraphs.set_prop!(graph, :transition_likelihoods, transition_likelihoods)

    # helpful for downstream processing
    unbranching_nodes = Set(Int[])
    for node in Graphs.vertices(graph)
        if (Graphs.indegree(graph, node) <= 1) && (Graphs.outdegree(graph, node) <= 1)
            push!(unbranching_nodes, node)
        end
    end
    branching_nodes = Set(setdiff(Graphs.vertices(graph), unbranching_nodes))
    MetaGraphs.set_prop!(graph, :unbranching_nodes, unbranching_nodes)
    MetaGraphs.set_prop!(graph, :branching_nodes, branching_nodes)
    
    
    # total_strand_normalized_quality_support = sum.(collect(values(strand_normalized_quality_support)))
    mean_total_support = Statistics.mean(collect(values(kmer_total_quality)))
    sorted_kmer_total_quality_values = collect(values(kmer_total_quality))
    mean_quality_value = Statistics.mean(sorted_kmer_total_quality_values)
    threshold = mean_quality_value

    xs = [
        [i for (i, y) in enumerate(sorted_kmer_total_quality_values) if y > threshold],
        [i for (i, y) in enumerate(sorted_kmer_total_quality_values) if y <= threshold]
        ]
    
    likely_valid_kmer_indices = xs[1]
    MetaGraphs.set_prop!(graph, :likely_valid_kmer_indices, likely_valid_kmer_indices)
    likely_sequencing_artifact_indices = xs[2]
    MetaGraphs.set_prop!(graph, :likely_sequencing_artifact_indices, likely_sequencing_artifact_indices)
    # likely_sequencing_artifact_kmers = Set(ordered_kmers[likely_sequencing_artifact_indices])
    # likely_valid_kmers = Set(ordered_kmers[likely_valid_kmer_indices])
    # kmer_to_index_map = Dict(kmer => i for (i, kmer) in enumerate(ordered_kmers))
    
    
    if plot
        ys = [
            [y for y in sorted_kmer_total_quality_values if y > threshold],
            [y for y in sorted_kmer_total_quality_values if y <= threshold]
        ]

        p = StatsPlots.scatter(
            xs,
            ys,
            title = "kmer qualities",
            ylabel = "canonical kmer cumulative QUAL value",
            label = ["above" "below"],
            legend = :outertopright,
            # size = (900, 500),
            margins=10StatsPlots.Plots.PlotMeasures.mm,
            xticks = false
        )
        p = StatsPlots.hline!(p, [mean_quality_value], label="mean")
        display(p)
    end
    return graph
end

# selected after trialing previous and next ks and finding those to be too unstable
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Performs iterative error correction on FASTQ sequences using progressively larger k-mer sizes.

Starting with the default k-mer size, this function repeatedly applies polishing steps,
incrementing the k-mer size until either reaching max_k or encountering instability.

# Arguments
- `fastq`: Path to input FASTQ file or FastqRecord object
- `max_k`: Maximum k-mer size to attempt (default: 89)
- `plot`: Whether to generate diagnostic plots (default: false)

# Returns
Vector of polishing results, where each element contains:
- k: k-mer size used
- fastq: resulting polished sequences
"""
function iterative_polishing(fastq, max_k = 89, plot=false)
    # initial polishing
    polishing_results = [polish_fastq(fastq=fastq)]
    while (!ismissing(last(polishing_results).k)) && (last(polishing_results).k < max_k)
        next_k = first(filter(k -> k > last(polishing_results).k, Mycelia.ks()))
        # @show next_k
        push!(polishing_results, polish_fastq(fastq=last(polishing_results).fastq, k=next_k))
    end
    return polishing_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Exports a taxonomy mapping table from a BLAST database in seqid2taxid format.

# Arguments
- `path_to_db::String`: Path to the BLAST database
- `outfile::String`: Output file path (defaults to input path + ".seqid2taxid.txt.gz")

# Returns
- `String`: Path to the created output file

# Details
Creates a compressed tab-delimited file mapping sequence IDs to taxonomy IDs.
Uses blastdbcmd without GI identifiers for better cross-referencing compatibility.
If the output file already exists, returns the path without regenerating.

# Dependencies
Requires BLAST+ tools installed via Bioconda.
"""
function export_blast_db_taxonomy_table(;path_to_db, outfile = path_to_db * ".seqid2taxid.txt.gz")
    Mycelia.add_bioconda_env("blast")
    if !isfile(outfile)
        # -long_seqids adds GI identifiers - these are cross-referenceable through other means so I'm dropping
        @time run(pipeline(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastdbcmd  -entry all -outfmt "%a %T" -db $(path_to_db)`, `gzip`), outfile))
    else
        @info "$(outfile) already present"
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Loads a BLAST database taxonomy mapping table from a gzipped file into a DataFrame.

# Arguments
- `compressed_blast_db_taxonomy_table_file::String`: Path to a gzipped file containing BLAST taxonomy mappings

# Returns
- `DataFrame`: A DataFrame with columns `:sequence_id` and `:taxid` containing the sequence-to-taxonomy mappings

# Format
Input file should be a space-delimited text file (gzipped) with two columns:
1. sequence identifier
2. taxonomy identifier (taxid)
"""
function load_blast_db_taxonomy_table(compressed_blast_db_taxonomy_table_file)
    return CSV.read(CodecZlib.GzipDecompressorStream(open(compressed_blast_db_taxonomy_table_file)), delim=' ', header=["sequence_id", "taxid"], DataFrames.DataFrame)
    # data, header = uCSV.read(CodecZlib.GzipDecompressorStream(open(compressed_blast_db_taxonomy_table_file)), delim=' ')
    # header = ["sequence_id", "taxid"]
    # DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse SAM/BAM files into a DataFrame containing mapped read alignments.

# Arguments
- `xam::String`: Path to input SAM/BAM file (supports .sam, .bam, or .sam.gz)
- `primary_only::Bool=false`: Flag to filter for primary alignments only (currently unused)

# Returns
DataFrame with columns:
- `template`: Read template name
- `flag`: SAM flag
- `reference`: Reference sequence name
- `position`: Alignment position range
- `mappingquality`: Mapping quality score
- `tlen`: Template length
- `alignlength`: Alignment length
- `ismapped`: Boolean indicating if read is mapped
- `isprimary`: Boolean indicating if alignment is primary
- `alignment_score`: Alignment score (AS tag)
- `mismatches`: Number of mismatches (NM tag)
"""
function parse_xam_to_mapped_records_table(xam, primary_only=false)
# merge name conflicts, leaving breadcrumb for reference
# function xam_records_to_dataframe(records)
    record_table = DataFrames.DataFrame(
        template = String[],
        flag = UInt16[],
        reference = String[],
        position = UnitRange{Int}[],
        mappingquality = UInt8[],
        tlen = Int[],
        alignlength = Int[],
        ismapped = Bool[],
        isprimary = Bool[],
        alignment_score = Int[],
        mismatches = Int[]
    )
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
        io = open(xam)
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
        io = open(xam)
    elseif occursin(r"\.sam.gz$", xam)
        MODULE = XAM.SAM
        io = CodecZlib.GzipDecompressorStream(open(xam))
    else
        error("unrecognized file extension in file: $xam")
    end
    # filter out header lines
    reader = MODULE.Reader(IOBuffer(join(Iterators.filter(line -> !startswith(line, '@'), eachline(io)), '\n')))
    # reader = MODULE.Reader(io)
    for record in reader
        if XAM.SAM.ismapped(record)
            row = (
                template = XAM.SAM.tempname(record),
                flag = XAM.flag(record),
                reference = XAM.SAM.refname(record),
                position = XAM.SAM.position(record):XAM.SAM.rightposition(record),
                mappingquality = XAM.SAM.mappingquality(record),
                tlen = XAM.SAM.templength(record),
                alignlength = XAM.SAM.alignlength(record),
                ismapped = XAM.SAM.ismapped(record),
                isprimary = XAM.SAM.isprimary(record),
                alignment_score = record["AS"],
                mismatches = record["NM"]
                )
            push!(record_table, row, promote=true)
        end
    end
    close(io)
    return record_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a SAM/BAM alignment file and extract template-to-reference mapping information.

# Arguments
- `xam::String`: Path to input alignment file (.sam, .sam.gz, or .bam format)

# Returns
- `DataFrame`: Table with columns:
  - `template`: Read template names
  - `reference`: Reference sequence names

Only includes primary alignments (not secondary/supplementary) that are mapped to references.
Skips header lines starting with '@'.
"""
function parse_xam_to_primary_mapping_table(xam)
    record_table = DataFrames.DataFrame(
        template = String[],
        reference = String[]
    )
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
        io = open(xam)
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
        io = open(xam)
    elseif occursin(r"\.sam.gz$", xam)
        MODULE = XAM.SAM
        io = CodecZlib.GzipDecompressorStream(open(xam))
    else
        error("unrecognized file extension in file: $xam")
    end
    # filter out header lines
    reader = MODULE.Reader(IOBuffer(join(Iterators.filter(line -> !startswith(line, '@'), eachline(io)), '\n')))
    # reader = MODULE.Reader(io)
    for record in reader
        if XAM.SAM.ismapped(record) && XAM.SAM.isprimary(record)
            row = (
                template = XAM.SAM.tempname(record),
                reference = XAM.SAM.refname(record)
                )
            push!(record_table, row, promote=true)
        end
    end
    close(io)
    return record_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse alignment data from SAM/BAM files into a structured DataFrame containing mapping quality
and taxonomic information.

# Arguments
- `xam::String`: Path to input file (.sam, .bam, or .sam.gz)
- `primary_only::Bool=false`: When true, include only primary alignments (currently not implemented)

# Returns
DataFrames.DataFrame with columns:
- `template`: Read name
- `flag`: SAM flag
- `reference`: Reference sequence name  
- `position`: Alignment position range
- `mappingquality`: Mapping quality score
- `tlen`: Template length
- `alignlength`: Alignment length
- `ismapped`: Boolean indicating if read is mapped
- `isprimary`: Boolean indicating if alignment is primary
- `alignment_score`: Alignment score (AS tag)
- `mismatches`: Number of mismatches (NM tag)

# Notes
- Skips unmapped reads
- Automatically detects and handles SAM/BAM/compressed SAM formats
"""
function parse_xam_to_taxonomic_mapping_quality(xam, primary_only=false)
# merge name conflicts, leaving breadcrumb for reference
# function xam_records_to_dataframe(records)
    record_table = DataFrames.DataFrame(
        template = String[],
        flag = UInt16[],
        reference = String[],
        position = UnitRange{Int}[],
        mappingquality = UInt8[],
        tlen = Int[],
        alignlength = Int[],
        ismapped = Bool[],
        isprimary = Bool[],
        alignment_score = Int[],
        mismatches = Int[]
    )
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
        io = open(xam)
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
        io = open(xam)
    elseif occursin(r"\.sam.gz$", xam)
        MODULE = XAM.SAM
        io = CodecZlib.GzipDecompressorStream(open(xam))
    else
        error("unrecognized file extension in file: $xam")
    end
    # filter out header lines
    reader = MODULE.Reader(IOBuffer(join(Iterators.filter(line -> !startswith(line, '@'), eachline(io)), '\n')))
    # reader = MODULE.Reader(io)
    for record in reader
        if XAM.SAM.ismapped(record)
            row = (
                template = XAM.SAM.tempname(record),
                flag = XAM.flag(record),
                reference = XAM.SAM.refname(record),
                position = XAM.SAM.position(record):XAM.SAM.rightposition(record),
                mappingquality = XAM.SAM.mappingquality(record),
                tlen = XAM.SAM.templength(record),
                alignlength = XAM.SAM.alignlength(record),
                ismapped = XAM.SAM.ismapped(record),
                isprimary = XAM.SAM.isprimary(record),
                alignment_score = record["AS"],
                mismatches = record["NM"]
                )
            push!(record_table, row, promote=true)
        end
    end
    close(io)
    return record_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Evaluate genome assembly quality by comparing k-mer distributions between assembled sequences and raw observations.

# Arguments
- `assembly`: Input assembled sequences to evaluate
- `observations`: Raw sequencing data for comparison
- `ks::Vector{Int}`: Vector of k-mer sizes to analyze (default: k=17 to 23)

# Returns
DataFrame containing quality metrics for each k-mer size:
- `k`: K-mer length used
- `cosine_distance`: Cosine similarity between k-mer distributions
- `js_divergence`: Jensen-Shannon divergence between distributions  
- `qv`: MerQury-style quality value score
"""
function assess_assembly_quality(;assembly, observations, ks::Vector{Int}=filter(x -> 17 <= x <= 23, Mycelia.ks()))
    results = DataFrames.DataFrame()
    @show ks
    ProgressMeter.@showprogress for k in ks
        @show k
        @info "counting assembly kmers..."
        assembled_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, assembly)
        # assembled_canonical_kmer_counts_file = Mycelia.jellyfish_count(fastx=assembly, k=k, canonical=true)
        # @info "loading assembly kmer counts..."
        # assembled_canonical_kmer_counts_table = Mycelia.load_jellyfish_counts(assembled_canonical_kmer_counts_file)
        # sort!(assembled_canonical_kmer_counts_table, "kmer")
        # assembled_canonical_kmer_counts = OrderedCollections.OrderedDict(row["kmer"] => row["count"] for row in DataFrames.eachrow(assembled_canonical_kmer_counts_table))
        
        @info "counting observation kmers..."
        observed_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, observations)
        # observed_canonical_kmer_counts_file = Mycelia.jellyfish_count(fastx=observations, k=k, canonical=true)
        # @info "loading observation kmer counts..."
        # observed_canonical_kmer_counts_table = Mycelia.load_jellyfish_counts(observed_canonical_kmer_counts_file)
        # sort!(observed_canonical_kmer_counts_table, "kmer")
        # observed_canonical_kmer_counts = OrderedCollections.OrderedDict(row["kmer"] => row["count"] for row in DataFrames.eachrow(observed_canonical_kmer_counts_table))
        cosine_distance = kmer_counts_to_cosine_similarity(observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        js_divergence = kmer_counts_to_js_divergence(observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        qv = kmer_counts_to_merqury_qv(raw_data_counts=observed_canonical_kmer_counts, assembly_counts=assembled_canonical_kmer_counts)
        push!(results, (;k, cosine_distance, js_divergence, qv))
    end
    return results
end

# function assess_assembly_quality(;assembled_sequence::BioSequences.LongDNA{2}, fastq::String, k::Int)
#     assess_assembly_quality(assembled_sequence=assembled_sequence, fastq=fastq, ks=[k])
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the cosine similarity between two k-mer count dictionaries.

# Arguments
- `kmer_counts_1::Dict{String,Int}`: First dictionary mapping k-mer sequences to their counts
- `kmer_counts_2::Dict{String,Int}`: Second dictionary mapping k-mer sequences to their counts

# Returns
- `Float64`: Cosine distance between the two k-mer count vectors, in range [0,1]
  where 0 indicates identical distributions and 1 indicates maximum dissimilarity

# Details
Converts k-mer count dictionaries into vectors using a unified set of keys,
then computes cosine distance. Missing k-mers are treated as count 0.
Result is invariant to input order and total counts (normalized internally).
"""
function kmer_counts_to_cosine_similarity(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    # Distances.cosine_dist(a, b) == Distances.cosine_dist(b, a) == Distances.cosine_dist(a ./ sum(a), b ./ sum(b)) == Distances.cosine_dist(b ./ sum(b), a ./ sum(a))
    return Distances.cosine_dist(a, b)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Jensen-Shannon divergence between two k-mer frequency distributions.

# Arguments
- `kmer_counts_1`: Dictionary mapping k-mers to their counts in first sequence
- `kmer_counts_2`: Dictionary mapping k-mers to their counts in second sequence

# Returns
- Normalized Jensen-Shannon divergence score between 0 and 1, where:
  - 0 indicates identical distributions
  - 1 indicates maximally different distributions

# Notes
- The measure is symmetric: JS(P||Q) = JS(Q||P)
- Counts are automatically normalized to probability distributions
"""
function kmer_counts_to_js_divergence(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    a_norm = a ./ sum(a)
    b_norm = b ./ sum(b)
    # Distances.js_divergence(a ./ sum(a), b ./ sum(b)) == Distances.js_divergence(b ./ sum(b), a ./ sum(a))
    # Distances.js_divergence(a, b) != Distances.js_divergence(a ./ sum(a), b ./ sum(b))
    return Distances.js_divergence(a_norm, b_norm)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Jaccard similarity coefficient between two sets.

The Jaccard similarity is defined as the size of the intersection divided by the size
of the union of two sets:

    J(A,B) = |A ∩ B| / |A ∪ B|

# Arguments
- `set1`: First set for comparison
- `set2`: Second set for comparison

# Returns
- `Float64`: A value between 0.0 and 1.0, where:
  * 1.0 indicates identical sets
  * 0.0 indicates completely disjoint sets
"""
function jaccard_similarity(set1, set2)
    return length(intersect(set1, set2)) / length(union(set1, set2))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Jaccard distance between two sets, which is the complement of the Jaccard similarity.

The Jaccard distance is defined as:
``J_d(A,B) = 1 - J_s(A,B) = 1 - \\frac{|A ∩ B|}{|A ∪ B|}``

# Arguments
- `set1`: First set to compare
- `set2`: Second set to compare

# Returns
- `Float64`: A value in [0,1] where 0 indicates identical sets and 1 indicates disjoint sets
"""
function jaccard_distance(set1, set2)
    return 1.0 - jaccard_similarity(set1, set2)
end

# function kmer_counts_to_jaccard(kmer_counts_1::AbstractDict{Kmers.DNAKmer{k}, Int64}, kmer_counts_2::AbstractDict{Kmers.DNAKmer{k}, Int64}) where k
#     # sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
#     # a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # a_indices = findall(a .> 0)
#     # b_indices = findall(b .> 0)
#     # return Distances.jaccard(a_indices, b_indices)
#     return jaccard(collect(keys(kmer_counts_1)), collect(keys(kmer_counts_2)))
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate assembly Quality Value (QV) score using the Merqury method.

Estimates base-level accuracy by comparing k-mer distributions between raw sequencing
data and assembly. Higher QV scores indicate better assembly quality.

# Arguments
- `raw_data_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}`: K-mer counts from raw sequencing data
- `assembly_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}`: K-mer counts from assembly

# Returns
- `Float64`: Quality Value score in Phred scale (-10log₁₀(error rate))

# Method
QV is calculated using:
1. Ktotal = number of unique kmers in assembly
2. Kshared = number of kmers shared between raw data and assembly
3. P = (Kshared/Ktotal)^(1/k) = estimated base-level accuracy
4. QV = -10log₁₀(1-P)

# Reference
Rhie et al. "Merqury: reference-free quality, completeness, and phasing assessment
for genome assemblies" Genome Biology (2020)
"""
function kmer_counts_to_merqury_qv(;raw_data_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}, assembly_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}) where {k,N}
    # Ktotal = # of kmers found in assembly
    Ktotal = length(keys(assembly_counts))
    # Kshared = # of shared kmers between assembly and readset
    Kshared = length(intersect(keys(raw_data_counts), keys(assembly_counts)))
    # probabilitiy_base_in_assembly_correct
    P = (Kshared/Ktotal)^(1/k)
    # # Error rate
    E = 1-P
    QV = -10log10(E)
    # return (;P, E, QV)
    return QV
end

# smaller, higher diversity databases do better with >=5 as the denominator - w/ <=4 they run out of memory
# denominator = 5 # produced OOM for NT on NERSC
# denominator = 8 # produced OOM for NT on Lawrencium
# denominator = 10 was only 56% efficient for NT on NERSC
const DEFAULT_MINIMAP_DENOMINATOR=10

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate appropriate minimap index size based on available system memory.

Converts total system memory to a recommended minimap index size by dividing the available
memory by a denominator factor. Returns the size as a string with 'G' suffix.

# Arguments
- `system_mem_gb::Number`: Total system memory in gigabytes
- `denominator::Number=DEFAULT_MINIMAP_DENOMINATOR`: Divisor to scale down memory allocation

# Returns
- `String`: Formatted memory size string (e.g. "4G")
"""
function system_mem_to_minimap_index_size(;system_mem_gb, denominator=DEFAULT_MINIMAP_DENOMINATOR)

    value = Int(floor(system_mem_gb/denominator))
    # 4G is the default
    # this value should be larger for larger memory machines, and smaller for smaller ones
    # it seems related to the total size of the sequences stored in memory, rather than the total size of the in-memory database
    return "$(value)G"
end
# system_mem_to_minimap_index_size(Mycelia.NERSC_MEM)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a minimap2 index command and output file path for mapping DNA sequencing reads.

Run this on the machine you intend to use to map the reads to confirm the index will fit.

# Arguments
- `fasta`: Path to the reference FASTA file to be indexed
- `mem_gb`: Available system memory in gigabytes for indexing
- `mapping_type`: Sequencing technology preset. One of:
  * "map-hifi": PacBio HiFi reads
  * "map-ont": Oxford Nanopore reads
  * "map-pb": PacBio CLR reads
  * "sr": Short reads
  * "lr:hq": High-quality long reads
- `threads`: Number of threads to use for indexing
- `as_string`: Return command as String instead of Cmd (default: false)
- `denominator`: Divisor for calculating index size (default: $(DEFAULT_MINIMAP_DENOMINATOR))

# Returns
Named tuple containing:
- `cmd`: The minimap2 indexing command (as String or Cmd)
- `outfile`: Path to the output index file (.mmi)
"""
function minimap_index(;fasta, mem_gb, mapping_type, threads, as_string=false, denominator=DEFAULT_MINIMAP_DENOMINATOR)
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    if as_string
        cmd = "$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -d $(index_file) $(fasta)"
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -d $(index_file) $(fasta)`
    end
    outfile = index_file
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate minimap2 alignment commands for sequence mapping.

aligning and compressing. No sorting or filtering.

Use shell_only=true to get string command to submit to SLURM

Creates a command to align reads in FASTQ format to a reference FASTA using minimap2, 
followed by SAM compression with pigz. Handles resource allocation and conda environment setup.

# Arguments
- `fasta`: Path to reference FASTA file
- `fastq`: Path to query FASTQ file
- `mapping_type`: Alignment preset ("map-hifi", "map-ont", "map-pb", "sr", or "lr:hq")
- `as_string`: If true, returns shell command as string; if false, returns command array
- `mem_gb`: Available memory in GB for indexing (defaults to system free memory)
- `threads`: Number of CPU threads to use (defaults to system threads)
- `denominator`: Divisor for calculating minimap2 index size

# Returns
Named tuple containing:
- `cmd`: Shell command (as string or array)
- `outfile`: Path to compressed output SAM file
"""
function minimap_map(;
        fasta,
        fastq,
        mapping_type,
        as_string=false,
        mem_gb=(Int(Sys.free_memory()) / 1e9),
        threads=Sys.CPU_THREADS,
        denominator=DEFAULT_MINIMAP_DENOMINATOR
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    temp_sam_outfile = fastq * "." * basename(fasta) * "." * "minimap2.sam"
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = [map, compress]
    end
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate minimap2 mapping commands with a pre-built index file.

# Arguments
- `fasta`: Path to the reference FASTA file
- `mem_gb`: Available system memory in gigabytes for index generation
- `mapping_type`: Mapping preset. Must be one of: "map-hifi", "map-ont", "map-pb", "sr", or "lr:hq"
- `threads`: Number of threads to use for mapping and compression
- `fastq`: Path to input FASTQ file
- `as_string`: If true, returns command as a string; if false, returns as command array
- `denominator`: Divisor for index size calculation (default: $(DEFAULT_MINIMAP_DENOMINATOR))

# Returns
Named tuple containing:
- `cmd`: The minimap2 mapping command (string or array)
- `outfile`: Path to the output compressed SAM file

# Notes
- Requires pre-built index file with pattern: `\${fasta}.x\${mapping_type}.I\${index_size}.mmi`
- Automatically installs required conda environments (minimap2, samtools, pigz)
- Output is automatically compressed with pigz
"""
function minimap_map_with_index(;
        fasta,
        mem_gb,
        mapping_type,
        threads,
        fastq,
        as_string=false,
        denominator=DEFAULT_MINIMAP_DENOMINATOR
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    @show index_file
    @assert isfile(index_file)
    temp_sam_outfile = fastq * "." * basename(index_file) * "." * "minimap2.sam"
    # outfile = temp_sam_outfile
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = [map, compress]
    end
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Find the longest common prefix between two filenames.

# Arguments
- `filename1::String`: First filename to compare
- `filename2::String`: Second filename to compare

# Keywords
- `strip_trailing_delimiters::Bool=true`: If true, removes trailing dots, hyphens, and underscores from the result

# Returns
- `String`: The longest common prefix found between the filenames
"""
function find_matching_prefix(filename1::String, filename2::String; strip_trailing_delimiters=true)
    min_length = min(length(filename1), length(filename2))
    matching_prefix = ""
    
    for i in 1:min_length
        if filename1[i] == filename2[i]
            matching_prefix *= filename1[i]
        else
            break
        end
    end
    if strip_trailing_delimiters
        matching_prefix = replace(matching_prefix, r"[\.\-_]+$" => "")
    end
    
    return matching_prefix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Map paired-end reads to a reference sequence using minimap2.

# Arguments
- `fasta::String`: Path to reference FASTA file
- `forward::String`: Path to forward reads FASTQ file
- `reverse::String`: Path to reverse reads FASTQ file
- `mem_gb::Integer`: Available system memory in GB
- `threads::Integer`: Number of threads to use
- `outdir::String`: Output directory (defaults to forward reads directory)
- `as_string::Bool=false`: Return command as string instead of Cmd array
- `mapping_type::String="sr"`: Minimap2 preset ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
- `denominator::Float64`: Memory scaling factor for index size

# Returns
Named tuple containing:
- `cmd`: Command(s) to execute (String or Array{Cmd})
- `outfile`: Path to compressed output SAM file (*.sam.gz)

# Notes
- Requires minimap2, samtools, and pigz conda environments
- Automatically compresses output using pigz
- Index file must exist at `\$(fasta).x\$(mapping_type).I\$(index_size).mmi`
"""
function minimap_map_paired_end_with_index(;
        fasta,
        forward,
        reverse,
        mem_gb,
        threads,
        outdir = dirname(forward),
        as_string=false,
        mapping_type="sr",
        denominator=DEFAULT_MINIMAP_DENOMINATOR
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    # @show index_file
    @assert isfile(index_file) "$(index_file) not found!!"
    @assert isfile(forward) "$(forward) not found!!"
    @assert isfile(reverse) "$(reverse) not found!!"
    fastq_prefix = find_matching_prefix(basename(forward), basename(reverse))
    temp_sam_outfile = joinpath(outdir, fastq_prefix) * "." * basename(index_file) * "." * "minimap2.sam"
    # outfile = temp_sam_outfile
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = [map, compress]
    end
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Maps paired-end reads to a reference genome using minimap2 and compresses the output.

# Arguments
- `fasta::String`: Path to reference genome FASTA file
- `forward::String`: Path to forward reads FASTQ file
- `reverse::String`: Path to reverse reads FASTQ file  
- `mem_gb::Integer`: Available system memory in GB
- `threads::Integer`: Number of threads to use
- `outdir::String`: Output directory (defaults to forward reads directory)
- `as_string::Bool`: Return command as string instead of Cmd array
- `mapping_type::String`: Mapping preset, e.g. "sr" for short reads (default)
- `denominator::Float64`: Memory scaling factor for minimap2 index

# Returns
Named tuple containing:
- `cmd`: Command(s) to execute (String or Vector{Cmd})
- `outfile`: Path to compressed output SAM file (*.sam.gz)

# Dependencies
Requires bioconda packages: minimap2, samtools, pigz
"""
function minimap_map_paired_end(;
        fasta,
        forward,
        reverse,
        mem_gb,
        threads,
        outdir = dirname(forward),
        as_string=false,
        mapping_type="sr",
        denominator=Mycelia.DEFAULT_MINIMAP_DENOMINATOR
    )
    index_size = Mycelia.system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    @assert isfile(forward) "$(forward) not found!!"
    @assert isfile(reverse) "$(reverse) not found!!"
    fastq_prefix = Mycelia.find_matching_prefix(basename(forward), basename(reverse))
    temp_sam_outfile = joinpath(outdir, fastq_prefix) * "." * "minimap2.sam"
    # outfile = temp_sam_outfile
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        fasta,
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_size) -ax $(mapping_type) $(fasta) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_size) -ax $(mapping_type) $(fasta) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = [map, compress]
    end
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a BAM file to FASTQ format with gzip compression.

# Arguments
- `bam`: Path to input BAM file
- `fastq`: Optional output path. Defaults to input path with ".fq.gz" extension

# Returns
- Path to the generated FASTQ file

# Details
- Uses samtools through conda environment
- Automatically skips if output file exists
- Output is gzip compressed
- Requires samtools to be available via conda

"""
function bam_to_fastq(;bam, fastq=bam * ".fq.gz")
    Mycelia.add_bioconda_env("samtools")
    bam_to_fastq_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools fastq $(bam)`
    gzip_cmd = `gzip`
    p = pipeline(bam_to_fastq_cmd, gzip_cmd)
    if !isfile(fastq)
        @time run(pipeline(p, fastq))
    else
        @info "$(fastq) already exists"
    end
    return fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize a dictionary of counts into a probability distribution where values sum to 1.0.

# Arguments
- `countmap::Dict`: Dictionary mapping keys to count values

# Returns
- `Dict`: New dictionary with same keys but values normalized by total sum
"""
function normalize_countmap(countmap)
    sum_total = sum(values(countmap))
    return Dict(k => v/sum_total for (k, v) in countmap)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract biosample and barcode information from a PacBio XML metadata file.

# Arguments
- `xml`: Path to PacBio XML metadata file

# Returns
DataFrame with two columns:
- `BioSampleName`: Name of the biological sample
- `BarcodeName`: Associated DNA barcode identifier
"""
function extract_pacbiosample_information(xml)
    xml_dict = XMLDict.parse_xml(read(xml, String))
    wellsample = xml_dict["ExperimentContainer"]["Runs"]["Run"]["Outputs"]["SubreadSets"]["SubreadSet"]["DataSetMetadata"]["Collections"]["CollectionMetadata"]["WellSample"]

    # Initialize empty arrays to store the data
    biosample_names = []
    barcode_names = []
    
    if haskey(wellsample, "BioSamples")
        # display(wellsample)
        for bs in wellsample["BioSamples"]["BioSample"]
            push!(biosample_names, bs[:Name])
            push!(barcode_names, bs["DNABarcodes"]["DNABarcode"][:Name])
        end
    end

    # Create the DataFrame
    df = DataFrames.DataFrame(BioSampleName=biosample_names, BarcodeName=barcode_names)
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a visualization of a genome assembly graph using Bandage.

# Arguments
- `gfa`: Path to input GFA (Graphical Fragment Assembly) file
- `img`: Optional output image path. Defaults to GFA filename with .png extension

# Returns
- Path to the generated image file
"""
function bandage_visualize(;gfa, img=gfa*".png")
    # run(`$(bandage) image --helpall`)
    bandage = Mycelia.download_bandage()
    if !isfile(img)
        run(`$(bandage) image $(gfa) $(img)`)
    end
    return img
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
function xam_to_contig_mapping_stats(xam)
    xam_results = Mycelia.parse_xam_to_summary_table(xam)
    xam_results = xam_results[xam_results[!, "isprimary"] .& xam_results[!, "ismapped"], :]
    # Calculate the percentage of mismatches
    xam_results.percent_mismatches = xam_results.mismatches ./ xam_results.alignlength * 100
    
    # Group by the 'reference' column and calculate the summary statistics
    contig_mapping_stats = DataFrames.combine(DataFrames.groupby(xam_results, :reference)) do subdf
        mappingquality_stats = StatsBase.summarystats(subdf.mappingquality)
        alignlength_stats = StatsBase.summarystats(subdf.alignlength)
        alignment_score_stats = StatsBase.summarystats(subdf.alignment_score)
        # mismatches_stats = StatsBase.summarystats(subdf.mismatches)
        percent_mismatches_stats = StatsBase.summarystats(subdf.percent_mismatches)

        (n_aligned_reads = length(subdf[!, "alignlength"]),
         total_aligned_bases = sum(subdf[!, "alignlength"]),
         total_alignment_score = sum(subdf[!, "alignment_score"]),
         mappingquality_mean = mappingquality_stats.mean,
         mappingquality_std = mappingquality_stats.sd,
         # mappingquality_min = mappingquality_stats.min,
         mappingquality_median = mappingquality_stats.median,
         # mappingquality_max = mappingquality_stats.max,

         alignlength_mean = alignlength_stats.mean,
         alignlength_std = alignlength_stats.sd,
         # alignlength_min = alignlength_stats.min,
         alignlength_median = alignlength_stats.median,
         # alignlength_max = alignlength_stats.max,

         alignment_score_mean = alignment_score_stats.mean,
         alignment_score_std = alignment_score_stats.sd,
         # alignment_score_min = alignment_score_stats.min,
         alignment_score_median = alignment_score_stats.median,
         # alignment_score_max = alignment_score_stats.max,

         # mismatches_mean = mismatches_stats.mean,
         # mismatches_std = mismatches_stats.sd,
         # mismatches_min = mismatches_stats.min,
         # mismatches_median = mismatches_stats.median,
         # mismatches_max = mismatches_stats.max,

         percent_mismatches_mean = percent_mismatches_stats.mean,
         percent_mismatches_std = percent_mismatches_stats.sd,
         # percent_mismatches_min = percent_mismatches_stats.min,
         percent_mismatches_median = percent_mismatches_stats.median,
         # percent_mismatches_max = percent_mismatches_stats.max)
        )
    end
    return contig_mapping_stats
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

Calculate mapping statistics by comparing sequence alignments (BAM/SAM) to a reference FASTA.

# Arguments
- `fasta::String`: Path to reference FASTA file
- `xam::String`: Path to alignment file (BAM or SAM format)

# Returns
DataFrame with columns:
- `contig`: Reference sequence name
- `contig_length`: Length of reference sequence
- `total_aligned_bases`: Total number of bases aligned to reference
- `mean_depth`: Average depth of coverage (total_aligned_bases/contig_length)
"""
function fasta_xam_mapping_stats(;fasta, xam)
    fastx_contig_lengths = fastx_to_contig_lengths(fasta)
    xam_stats = xam_to_contig_mapping_stats(xam)
    fastx_contig_lengths = fastx_to_contig_lengths(fasta)
    fastx_contig_lengths_table = DataFrames.DataFrame(contig = collect(keys(fastx_contig_lengths)), contig_length = collect(values(fastx_contig_lengths)))
    fastx_contig_mapping_stats_table = DataFrames.innerjoin(fastx_contig_lengths_table, xam_stats, on="contig" => "reference")
    mean_depth = fastx_contig_mapping_stats_table[!, "total_aligned_bases"] ./ fastx_contig_mapping_stats_table[!, "contig_length"]
    DataFrames.insertcols!(fastx_contig_mapping_stats_table, 4, :mean_depth => mean_depth)
    return fastx_contig_mapping_stats_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate alignment statistics for a SAM/BAM/CRAM file using samtools flagstat.

# Arguments
- `xam::AbstractString`: Path to input SAM/BAM/CRAM alignment file
- `samtools_flagstat::AbstractString`: Output path for flagstat results (default: input_path.samtools-flagstat.txt)

# Returns
- `String`: Path to the generated flagstat output file

# Details
Runs samtools flagstat to calculate statistics on the alignment file, including:
- Total reads
- Secondary alignments
- Supplementary alignments  
- Duplicates
- Mapped/unmapped reads
- Proper pairs
- Read 1/2 counts

# Requirements
- Requires samtools to be available via Bioconda
- Input file must be in SAM, BAM or CRAM format
"""
function run_samtools_flagstat(xam, samtools_flagstat=xam * ".samtools-flagstat.txt")
    Mycelia.add_bioconda_env("samtools")
    if !isfile(samtools_flagstat)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools flagstat $(xam)`, samtools_flagstat))
    end
    return samtools_flagstat
end

# not a very good function yet, but good enough for the pinches I need it for
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a GFA (Graphical Fragment Assembly) file into a MetaGraph representation.

# Arguments
- `gfa`: Path to GFA format file

# Returns
A `MetaGraph` where:
- Vertices represent segments (contigs)
- Edges represent links between segments
- Vertex properties include `:id` with segment identifiers
- Graph property `:records` contains the original FASTA records

# Format Support
Handles standard GFA v1 lines:
- `H`: Header lines (skipped)
- `S`: Segments (stored as nodes with FASTA records)
- `L`: Links (stored as edges)
- `P`: Paths (stored in paths dictionary)
- `A`: HiFiAsm specific lines (skipped)
"""
function parse_gfa(gfa)
    segments = Vector{FASTX.FASTA.Record}()
    links = Vector{Pair{String, String}}()
    paths = Dict{String, Vector{String}}()
    for l in eachline(open(gfa))
        s = split(l, '\t')
        if first(s) == "H"
            # header line
            continue
        elseif first(s) == "S"
            # segment
            # push!(segments, string(s[2]))
            identifier = string(s[2])
            description = string(s[4])
            sequence = string(s[3])
            push!(segments, FASTX.FASTA.Record("$(identifier) $(description)", sequence))
        elseif first(s) == "L"
            # link
            push!(links, string(s[2]) => string(s[4]))
        elseif first(s) == "P"
            # path
            paths[string(s[2])] = string.(split(replace(s[3], r"[+-]" => ""), ','))
        elseif first(s) == "A" # hifiasm https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#output-file-formats
            continue
        else
            println(l)
            error("unexpected line encountered while parsing GFA")
        end
    end

    g = MetaGraphs.MetaGraph(length(segments))

    for link in links
        (u, v) = link
        ui = findfirst(FASTX.identifier.(segments) .== u)
        vi = findfirst(FASTX.identifier.(segments) .== v)
        Graphs.add_edge!(g, ui => vi)
    end
    for (i, segment) in enumerate(segments)
        MetaGraphs.set_prop!(g, i, :id, FASTX.identifier(segment))
    end
    MetaGraphs.set_prop!(g, :records, segments)
    return g
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GFA (Graphical Fragment Assembly) file into a structured representation.

# Arguments
- `gfa`: Path to GFA file or GFA content as string

# Returns
Named tuple containing:
- `contig_table`: DataFrame with columns:
  - `connected_component`: Integer ID for each component
  - `contigs`: Comma-separated list of contig IDs
  - `is_circular`: Boolean indicating if component forms a cycle
  - `is_closed`: Boolean indicating if single contig forms a cycle
  - `lengths`: Comma-separated list of contig lengths
- `records`: FASTA records from the GFA
"""
function gfa_to_structure_table(gfa)
    gfa_metagraph = parse_gfa(gfa)
    contig_table = DataFrames.DataFrame()
    records = MetaGraphs.get_prop(gfa_metagraph, :records)
    contig_lengths = Dict(FASTX.identifier(record) => length(FASTX.sequence(record)) for record in records)
    # @show String.(FASTX.description.(records))
    # try
    # contig_depths = Dict(FASTX.identifier(record) => first(match(r"^.*?dp:i:(\d+).*$", String(FASTX.description(record))).captures) for record in records)
    for (i, connected_component) in enumerate(Graphs.connected_components(gfa_metagraph))
        subgraph, node_map = Graphs.induced_subgraph(gfa_metagraph, connected_component)
        # display(subgraph)
        contigs = [MetaGraphs.get_prop(subgraph, v, :id) for v in Graphs.vertices(subgraph)]
        row = (
            connected_component = i,
            contigs = join(contigs, ","),
            is_circular = Graphs.is_cyclic(subgraph),
            is_closed = (length(contigs) == 1) && Graphs.is_cyclic(subgraph),
            lengths = join([contig_lengths[contig] for contig in contigs], ","),
            # depths = join([contig_depths[contig] for contig in contigs], ","),
            )
        push!(contig_table, row)
    end
    
    return (;contig_table, records)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Copy a file to a new location with a unique identifier prepended to the filename.

# Arguments
- `infile::AbstractString`: Path to the source file to copy
- `out_directory::AbstractString`: Destination directory for the copied file
- `unique_identifier::AbstractString`: String to prepend to the filename
- `force::Bool=true`: If true, overwrite existing files

# Returns
- `String`: Path to the newly created file
"""
function copy_with_unique_identifier(infile, out_directory, unique_identifier; force=true)
    outfile = joinpath(out_directory, unique_identifier * "." * basename(infile))
    cp(infile, outfile, force=force)
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Clustal Omega multiple sequence alignment on a FASTA file.

# Arguments
- `fasta::String`: Path to input FASTA file
- `outfmt::String="clustal"`: Output format for the alignment

# Returns
- `String`: Path to the output alignment file

# Supported Output Formats
- `"fasta"`: FASTA format
- `"clustal"`: Clustal format
- `"msf"`: MSF format  
- `"phylip"`: PHYLIP format
- `"selex"`: SELEX format
- `"stockholm"`: Stockholm format
- `"vienna"`: Vienna format

# Notes
- Uses Bioconda to manage the Clustal Omega installation
- Caches results - will return existing output file if already generated
- Handles single sequence files gracefully by returning output path without error
"""
function run_clustal_omega(;fasta, outfmt="clustal")
    Mycelia.add_bioconda_env("clustalo")
    outfile = "$(fasta).clustal-omega.$(outfmt)"
    if !isfile(outfile)
        try
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n clustalo clustalo -i $(fasta) --outfmt $(outfmt) -o $(outfile)`)
        catch e
            # FATAL: File '...' contains 1 sequence, nothing to align
            return outfile
        end
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Gets the size of a file and returns it in a human-readable format.

# Arguments
- `f`: The path to the file, either as a `String` or an `AbstractString`.

# Returns
A string representing the file size in a human-readable format (e.g., "3.40 MB").

# Details
This function internally uses `filesize(f)` to get the file size in bytes, then leverages `Base.format_bytes` to convert it into a human-readable format with appropriate units (KB, MB, GB, etc.).

# Examples
```julia
julia> filesize_human_readable("my_image.jpg")
"2.15 MB"
```
See Also
* filesize: Gets the size of a file in bytes.
* Base.format_bytes: Converts a byte count into a human-readable string. 
"""
function filesize_human_readable(f)
    return Base.format_bytes(filesize(f))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the 'padloc' tool from the 'padlocbio' conda environment on a given FASTA file.

https://doi.org/10.1093/nar/gkab883
https://github.com/padlocbio/padloc

This function first ensures that the 'padloc' environment is available via Bioconda. 
It then attempts to update the 'padloc' database. 
If a 'padloc' output file (with a '_padloc.csv' suffix) does not already exist for the input FASTA file, 
it runs 'padloc' with the specified FASTA file as input.
"""
function run_padloc(fasta_file)
    Mycelia.add_bioconda_env("padlocbio::padloc")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --db-update`)
    padloc_outfile = replace(fasta_file, ".fna" => "") * "_padloc.csv"
    if !isfile(padloc_outfile)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --fna $(fasta_file)`)
    else
        @info "$(padloc_outfile) already present"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Multi-Locus Sequence Typing (MLST) analysis on a genome assembly.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing the genome assembly

# Returns
- Path to the output file containing MLST results (`<input>.mlst.out`)

# Details
Uses the `mlst` tool from PubMLST to identify sequence types by comparing allelic 
profiles of housekeeping genes against curated MLST schemes.

# Dependencies
- Requires Bioconda and the `mlst` package
- Automatically sets up conda environment if not present
"""
function run_mlst(fasta_file)
    Mycelia.add_bioconda_env("mlst")    
    mlst_outfile = "$(fasta_file).mlst.out"
    @show mlst_outfile
    if !isfile(mlst_outfile)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mlst mlst $(fasta_file)`, mlst_outfile))
    else
        @info "$(mlst_outfile) already present"
    end
    return mlst_outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run ECTyper for serotyping E. coli genome assemblies.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing assembled genome(s)

# Returns
- `String`: Path to output directory containing ECTyper results
"""
function run_ectyper(fasta_file)
    Mycelia.add_bioconda_env("ectyper")
    outdir = fasta_file * "_ectyper"
    if !isdir(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ectyper ectyper -i $(fasta_file) -o $(outdir)`)
    else
        @info "$(outdir) already present"
    end
    return outdir
end

 # -i ecoliA.fasta -o output_dir

# function run_mlst(ID, OUT_DIR, normalized_fasta_file)
#     mlst_dir="$(OUT_DIR)/mlst"
#     if !isdir(mlst_dir)
#         mkdir(mlst_dir)
#     end
#     p = pipeline(
#             `mlst $(normalized_fasta_file)`,
#             stdout="$(mlst_dir)/$(ID).mlst.out")
#     run(p)
#     return mlst_dir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a JSONL (JSON Lines) file into a vector of dictionaries.

# Arguments
- `filepath::String`: Path to the JSONL file to parse

# Returns
- `Vector{Dict{String, Any}}`: Vector containing parsed JSON objects, one per line

# Description
Reads a JSONL file line by line, parsing each line as a separate JSON object.
Uses pre-allocation and progress tracking for efficient processing of large files.
"""
function parse_jsonl(filepath::String)
    # Count the lines first
    num_lines = countlines(filepath)

    # Pre-allocate the array
    json_objects = Vector{Dict{String, Any}}(undef, num_lines)

    # Progress meter setup
    p = ProgressMeter.Progress(num_lines; desc="Parsing JSONL: ", showspeed=true)
    
    open(filepath, "r") do file
        for (i, line) in enumerate(eachline(file))
            json_objects[i] = JSON.parse(line)
            ProgressMeter.next!(p) # Update the progress meter
        end
    end
    return json_objects
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a JSONL (JSON Lines) file into a vector of dictionaries.

# Arguments
- `filepath::String`: Path to the JSONL file to parse

# Returns
- `Vector{Dict{String, Any}}`: Vector containing parsed JSON objects, one per line

# Description
Reads a JSONL file line by line, parsing each line as a separate JSON object.
Uses pre-allocation and progress tracking for efficient processing of large files.
"""
function system_overview(;path=pwd())
    total_memory = Base.format_bytes(Sys.total_memory())
    available_memory = Base.format_bytes(Sys.free_memory())
    occupied_memory = Base.format_bytes(Sys.total_memory() - Sys.free_memory())
    system_threads = Sys.CPU_THREADS
    julia_threads = Threads.nthreads()
    available_storage = Base.format_bytes(Base.diskstat(path).available)
    total_storage = Base.format_bytes(Base.diskstat(path).total)
    occupied_storage = Base.format_bytes(Base.diskstat(path).used)
    return (;
            system_threads,
            julia_threads,
            total_memory,
            available_memory,
            occupied_memory,
            total_storage,
            available_storage,
            occupied_storage)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Vertically concatenate DataFrames with different column structures by automatically handling missing values.

# Arguments
- `dfs`: Variable number of DataFrames to concatenate vertically

# Returns
- `DataFrame`: Combined DataFrame containing all rows and columns from input DataFrames, 
  with `missing` values where columns didn't exist in original DataFrames
"""
function vcat_with_missing(dfs::Vararg{DataFrames.AbstractDataFrame})
    # Get all unique column names across all DataFrames
    all_columns = unique(reduce(vcat, [names(df) for df in dfs]))

    # Add missing columns to each DataFrame and fill with missing
    for df in dfs
        for col in all_columns
            if !(col in names(df))
                df[!, col] .= missing
            end
        end
    end

    # Now you can safely vcat the DataFrames
    return vcat(dfs...)
end

# always interpret as strings to ensure changes in underlying biosequence representation don't change results
# results in 64 character string
# a = "TTANC"
# b = "ttANc"
# c = "ttanc"
# dna_a = BioSequences.LongDNA{4}(a)
# dna_b = BioSequences.LongDNA{4}(b)
# dna_c = BioSequences.LongDNA{4}(c)
# seq2sha256(a) == seq2sha256(dna_a)
# seq2sha256(b) == seq2sha256(dna_b)
# seq2sha256(c) == seq2sha256(dna_c)
# seq2sha256(BioSequences.LongDNA{2}("AAA")) == seq2sha256(BioSequences.LongDNA{4}("AAA"))

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the SHA-256 hash of a sequence string.

# Arguments
- `seq::AbstractString`: Input sequence to be hashed

# Returns
- `String`: Hexadecimal representation of the SHA-256 hash

# Details
The input sequence is converted to uppercase before hashing.
"""
function seq2sha256(seq::AbstractString)
    return SHA.bytes2hex(SHA.sha256(uppercase(seq)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a biological sequence to its SHA256 hash value.

Calculates a cryptographic hash of the sequence by first converting it to a string representation.
This method dispatches to the string version of `seq2sha256`.

# Arguments
- `seq::BioSequences.BioSequence`: The biological sequence to hash

# Returns
- `String`: A 64-character hexadecimal string representing the SHA256 hash
"""
function seq2sha256(seq::BioSequences.BioSequence)
    return seq2sha256(string(seq))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute a single SHA256 hash from multiple SHA256 hashes.

Takes a vector of hex-encoded SHA256 hashes and produces a new SHA256 hash by:
1. Sorting the input hashes lexicographically
2. Concatenating them in sorted order
3. Computing a new SHA256 hash over the concatenated data

# Arguments
- `vector_of_sha256s`: Vector of hex-encoded SHA256 hash strings

# Returns
- A hex-encoded string representing the computed meta-hash
"""
function metasha256(vector_of_sha256s::Vector{<:AbstractString})
    ctx = SHA.SHA2_256_CTX()
    for sha_hash in sort(vector_of_sha256s)
        SHA.update!(ctx, collect(codeunits(sha_hash)))
    end
    return SHA.bytes2hex(SHA.digest!(ctx))
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function fasta2normalized_table(fasta_file)
#     @assert isfile(fasta_file) && filesize(fasta_file) > 0
#     n_records = Mycelia.count_records(fasta_file)

#     # Pre-allocate arrays
#     sequence_sha256s = Vector{String}(undef, n_records)
#     sequence_identifiers = Vector{String}(undef, n_records)
#     sequence_descriptions = Vector{String}(undef, n_records)
#     sequences = Vector{String}(undef, n_records)

#     # Progress meter
#     p = ProgressMeter.Progress(n_records, desc="Processing FASTA records: ", color=:blue)

#     for (i, record) in enumerate(Mycelia.open_fastx(fasta_file))
#         sequence_identifiers[i] = FASTX.identifier(record)
#         sequence_descriptions[i] = FASTX.description(record)
#         sequences[i] = FASTX.sequence(record)
#         sequence_sha256s[i] = Mycelia.seq2sha256(FASTX.sequence(record))
#         ProgressMeter.next!(p)
#     end

#     normalized_fasta_table = DataFrames.DataFrame(
#         fasta_sha256 = metasha256(sequence_sha256s),
#         fasta_identifier = basename(fasta_file),
#         sequence_sha256 = sequence_sha256s,
#         sequence_identifier = sequence_identifiers,
#         sequence_description = sequence_descriptions,
#         sequence = sequences
#     )
#     return normalized_fasta_table
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract the base file extension from a filename, handling compressed files.

For regular files, returns the last extension. For gzipped files, returns the extension
before .gz.
"""
function get_base_extension(filename::String)
  parts = split(basename(filename), "."; limit=3)  # Limit to 3 to handle 2-part extensions
  extension = parts[end]  # Get the last part
  
  if extension == "gz" && length(parts) > 2  # Check for .gz and more parts
    extension = parts[end - 1]  # Get the part before .gz
  end
  
  return "." * extension
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a FASTA file into a normalized tab-separated table format with standardized sequence identifiers.

# Arguments
- `fasta_file::String`: Path to input FASTA file
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
function fasta2normalized_table(fasta_file;
    normalize_name=true,
    outfile="",
    force=false)
    @assert isfile(fasta_file) && filesize(fasta_file) > 0

    if !isempty(outfile) && normalize_name
        @warn "Both 'outfile' and 'normalize_name' are set. Ignoring 'outfile'."
    end
    if isempty(outfile) && !normalize_name
        @warn "Both 'outfile' and 'normalize_name' are unset. Setting 'outfile' to '$(fasta_file).tsv.gz'."
    end

    if isempty(outfile)
        if !normalize_name
            outfile = fasta_file * ".tsv.gz"
        else
            filename = Mycelia.normalized_current_datetime() * "." * basename(fasta_file) * ".tsv.gz"
            outfile = joinpath(dirname(fasta_file), filename)
        end
    end

    if isfile(outfile) && (filesize(outfile) > 0) && !force
        @warn "$(outfile) already present and non-empty"
        return outfile
    end

    # Progress meter
    n_records = Mycelia.count_records(fasta_file)
    p = ProgressMeter.Progress(n_records, desc="Processing FASTA records: ")

    # Open the output file for writing
    outfile_io = CodecZlib.GzipCompressorStream(open(outfile, "w"))

    # Write the header
    header = "fasta_identifier\tsequence_sha256\tsequence_identifier\tsequence_description\tsequence"
    println(outfile_io, header)

    # Vector to store SHAs
    sequence_sha256s = Vector{String}(undef, n_records)

    # Streaming processing
    for (i, record) in enumerate(Mycelia.open_fastx(fasta_file))
        seq_id = FASTX.identifier(record)
        seq_desc = FASTX.description(record)
        seq = String(FASTX.sequence(record))
        sha256 = Mycelia.seq2sha256(seq)

        # Write to the temporary file
        line = join([basename(fasta_file), string(sha256), String(seq_id), String(seq_desc), String(seq)], '\t')
        println(outfile_io, line)

        sequence_sha256s[i] = sha256
        ProgressMeter.next!(p)
    end

    # Close the output file
    close(outfile_io)

    if normalize_name
        full_extension = Mycelia.get_base_extension(fasta_file) * ".tsv.gz"
        final_outfile = joinpath(dirname(outfile), Mycelia.metasha256(sequence_sha256s) * full_extension)
        mv(outfile, final_outfile; force=true)
        outfile = final_outfile
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Select one random row from each group in a grouped DataFrame.

# Arguments
- `gdf::GroupedDataFrame`: A grouped DataFrame created using `groupby`

# Returns
- `DataFrame`: A new DataFrame containing exactly one randomly sampled row from each group
"""
function rand_of_each_group(gdf::DataFrames.GroupedDataFrame{DataFrames.DataFrame})
    result = DataFrames.combine(gdf) do sdf
        sdf[StatsBase.sample(1:DataFrames.nrow(sdf), 1), :]
    end
    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a random symmetric distance matrix of size n×n with zeros on the diagonal.

# Arguments
- `n`: Positive integer specifying the matrix dimensions

# Returns
- A symmetric n×n matrix with random values in [0,1), zeros on the diagonal

# Details
- The matrix is symmetric, meaning M[i,j] = M[j,i]
- Diagonal elements M[i,i] are set to 0.0
- Off-diagonal elements are uniformly distributed random values
"""
function random_symmetric_distance_matrix(n)
  # Generate a random matrix
  matrix = rand(n, n)

  # Make the matrix symmetric
  matrix = (matrix + matrix') / 2

  # Ensure the diagonal is zero
  for i in 1:n
    matrix[i, i] = 0.0
  end

  return matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return a DataFrame containing the first row from each group in a GroupedDataFrame.

# Arguments
- `gdf::GroupedDataFrame`: A grouped DataFrame created using `groupby`

# Returns
- `DataFrame`: A new DataFrame containing first row from each group
"""
function first_of_each_group(gdf::DataFrames.GroupedDataFrame{DataFrames.DataFrame})
    return DataFrames.combine(gdf, first)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate `n` colors that are maximally distinguishable from each other.

# Arguments
- `n::Integer`: The number of distinct colors to generate

# Returns
A vector of `n` RGB colors that are optimized for maximum perceptual distinction,
using white (RGB(1,1,1)) and black (RGB(0,0,0)) as anchor colors.
"""
function n_maximally_distinguishable_colors(n)
    return Colors.distinguishable_colors(n, [Colors.RGB(1,1,1), Colors.RGB(0,0,0)], dropseed=true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a hierarchical clustering tree into a directed metagraph representation.

# Arguments
- `hcl::Clustering.Hclust`: Hierarchical clustering result object

# Returns
- `MetaGraphs.MetaDiGraph`: Directed graph with metadata representing the clustering hierarchy

# Graph Properties
The resulting graph contains the following vertex properties:
- `:hclust_id`: String identifier for each node
- `:height`: Height/distance at each merge point (0.0 for leaves)
- `:x`: Horizontal position for visualization (0-1 range)
- `:y`: Vertical position based on normalized height
- `:hcl`: Original clustering object (stored as graph property)
"""
function hclust_to_metagraph(hcl::Clustering.Hclust)
    total_nodes = length(hcl.order) + size(hcl.merges, 1)
    mg = MetaGraphs.MetaDiGraph(total_nodes)
    MetaGraphs.set_prop!(mg, :hcl, hcl)
    for leaf_node in hcl.labels
        MetaGraphs.set_prop!(mg, leaf_node, :hclust_id, string(-leaf_node))
        MetaGraphs.set_prop!(mg, leaf_node, :height, 0.0)
    end
    for (i, ordered_leaf_node) in enumerate(hcl.order)
        x = i / (length(hcl.order) + 1)
        MetaGraphs.set_prop!(mg, ordered_leaf_node, :x, x)
    end
    for (i, (left, right)) in enumerate(eachrow(hcl.merges))
        graph_i = length(hcl.order) + i
        hclust_id = string(i)
        MetaGraphs.set_prop!(mg, graph_i, :hclust_id, hclust_id)
    end
    MetaGraphs.set_indexing_prop!(mg, :hclust_id)

    for (i, ((left, right), height)) in enumerate(zip(eachrow(hcl.merges), hcl.heights))
        parent_vertex = mg[string(i), :hclust_id]
        MetaGraphs.set_prop!(mg, parent_vertex, :height, height)
        
        left_child_vertex = mg[string(left), :hclust_id]
        right_child_vertex = mg[string(right), :hclust_id]

        x = Statistics.middle(mg.vprops[left_child_vertex][:x], mg.vprops[right_child_vertex][:x])
        MetaGraphs.set_prop!(mg, parent_vertex, :x, x)

        e1 = Graphs.Edge(parent_vertex, left_child_vertex)
        e2 = Graphs.Edge(parent_vertex, right_child_vertex)
        Graphs.add_edge!(mg, e1)
        Graphs.add_edge!(mg, e2)
    end

    max_height = maximum(get.(values(mg.vprops), :height, missing))
    for v in Graphs.vertices(mg)
        relative_height = (max_height - mg.vprops[v][:height]) / max_height
        # y = 1 - relative_height
        MetaGraphs.set_prop!(mg, v, :y, relative_height)
    end

    return mg
end

# https://www.giantfocal.com/toolkit/font-size-converter
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert pixel measurements to point measurements using the standard 4:3 ratio.

Points are the standard unit for typography (1 point = 1/72 inch), while pixels are 
used for screen measurements. This conversion uses the conventional 4:3 ratio where 
3 points equal 4 pixels.

# Arguments
- `pixels`: The number of pixels to convert

# Returns
- The equivalent measurement in points
"""
function pixels_to_points(pixels)
    return pixels / 4 * 3
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert typographic points to pixels using a 4:3 ratio (1 point = 4/3 pixels).

# Arguments
- `points`: Size in typographic points (pt)

# Returns
- Size in pixels (px)
"""
function points_to_pixels(points)
    return points / 3 * 4
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Draw a dendrogram visualization of hierarchical clustering results stored in a MetaDiGraph.

# Arguments
- `mg::MetaGraphs.MetaDiGraph`: Graph containing hierarchical clustering results. Must have
  `:hcl` in graph properties with clustering data and vertex properties containing `:x`, `:y` coordinates.

# Keywords
- `width::Integer=500`: Width of output image in pixels
- `height::Integer=500`: Height of output image in pixels 
- `fontsize::Integer=12`: Font size for node labels in points
- `margins::Float64`: Margin size in pixels, defaults to min(width,height)/20
- `mergenodesize::Float64=1`: Size of circular nodes at merge points
- `lineweight::Float64=1`: Thickness of dendrogram lines
- `filename::String`: Output filename, defaults to timestamp with .dendrogram.png extension

# Returns
Nothing, but saves dendrogram image to disk and displays preview.
"""
function draw_dendrogram_tree(
        mg::MetaGraphs.MetaDiGraph;
        width=500,
        height=500,
        fontsize=12,
        # margins=min(width, height)/25,
        margins=min(width, height)/20,
        # margins=20,
        mergenodesize=1,
        lineweight=1,
        filename=Dates.format(Dates.now(), "yyyymmddTHHMMSS") * ".dendrogram.png"
    )

    fontsizebuffer = points_to_pixels(fontsize) / 3
    available_width = width - 2 * margins
    available_height = height - 2 * margins

    leaf_nodes = mg.gprops[:hcl].labels

    # Create a new drawing
    Luxor.Drawing(width, height, filename)
    # Luxor.origin()  # Set the origin to the center of the drawing
    # Luxor.background("white")
    Luxor.sethue("black")
    Luxor.fontsize(fontsize)
    Luxor.setline(lineweight)

    for ordered_leaf_node in mg.gprops[:hcl].order
        x = mg.vprops[ordered_leaf_node][:x] * available_width + margins
        y = mg.vprops[ordered_leaf_node][:y] * available_height + margins
        Luxor.circle(x, y, mergenodesize, :fill)
        # Luxor.text(string(ordered_leaf_node), Luxor.Point(x, y), halign=:center, valign=:middle)
        Luxor.text(string(ordered_leaf_node), Luxor.Point(x, y + fontsizebuffer), halign=:center, valign=:top)
    end
    for (i, (left, right)) in enumerate(eachrow(mg.gprops[:hcl].merges))
        parent_vertex = mg[string(i), :hclust_id]
                                        
        x = mg.vprops[parent_vertex][:x] * available_width + margins
        y = mg.vprops[parent_vertex][:y] * available_height + margins
        # Luxor.text(string(ordered_leaf_node), Luxor.Point(x, y), halign=:center, valign=:middle)
        Luxor.circle(x, y, mergenodesize, :fill)

        left_child_vertex = mg[string(left), :hclust_id]
        left_child_x = mg.vprops[left_child_vertex][:x] * available_width + margins
        left_child_y = mg.vprops[left_child_vertex][:y] * available_height + margins
        # draw horizontal bar
        Luxor.line(Luxor.Point(left_child_x, y), Luxor.Point(x, y), action=:stroke)
        # draw vertical bar
        Luxor.line(Luxor.Point(left_child_x, left_child_y), Luxor.Point(left_child_x, y), action=:stroke)

        right_child_vertex = mg[string(right), :hclust_id]
        right_child_x = mg.vprops[right_child_vertex][:x] * available_width + margins
        right_child_y = mg.vprops[right_child_vertex][:y] * available_height + margins
        # draw horizontal bar
        Luxor.line(Luxor.Point(x, y), Luxor.Point(right_child_x, y), action=:stroke)
        # draw vertical bar
        Luxor.line(Luxor.Point(right_child_x, right_child_y), Luxor.Point(right_child_x, y), action=:stroke)
    end

    # Finish the drawing and save the file
    Luxor.finish()
    Luxor.preview()
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Draw a radial hierarchical clustering tree visualization and save it as an image file.

# Arguments
- `mg::MetaGraphs.MetaDiGraph`: A meta directed graph containing hierarchical clustering data
  with required graph properties `:hcl` containing clustering information.

# Keywords
- `width::Int=500`: Width of the output image in pixels
- `height::Int=500`: Height of the output image in pixels
- `fontsize::Int=12`: Font size for node labels
- `margins::Float64`: Margin size (automatically calculated as min(width,height)/20)
- `mergenodesize::Float64=1`: Size of the merge point nodes
- `lineweight::Float64=1`: Thickness of the connecting lines
- `filename::String`: Output filename (defaults to timestamp with ".radial.png" suffix)

# Details
The function creates a radial visualization of hierarchical clustering results where:
- Leaf nodes are arranged in a circle with labels
- Internal nodes represent merge points
- Connections show the hierarchical structure through arcs and lines

The visualization is saved as a PNG file and automatically previewed.

# Required Graph Properties
The input graph must have:
- `mg.gprops[:hcl].labels`: Vector of leaf node labels
- `mg.gprops[:hcl].order`: Vector of ordered leaf nodes
- `mg.gprops[:hcl].merges`: Matrix of merge operations
- `mg.vprops[v][:x]`: X coordinate for each vertex
- `mg.vprops[v][:y]`: Y coordinate for each vertex
"""
function draw_radial_tree(
        mg::MetaGraphs.MetaDiGraph;
        width=500,
        height=500,
        fontsize=12,
        # margins=min(width, height)/25,
        margins=min(width, height)/20,
        # margins=20,
        mergenodesize=1,
        lineweight=1,
        filename=Dates.format(Dates.now(), "yyyymmddTHHMMSS") * ".radial.png"
    )

    # check max label size against margin and update as necessary
    max_radius = (min(width, height) - (margins * 2)) / 2

    # fontsizebuffer = points_to_pixels(fontsize) / 3
    available_width = width - 2 * margins
    available_height = height - 2 * margins
    leaf_nodes = mg.gprops[:hcl].labels

    # Create a new drawing
    Luxor.Drawing(width, height, filename)
    Luxor.origin()  # Set the origin to the center of the drawing
    # Luxor.background("white")
    Luxor.sethue("black")
    Luxor.fontsize(fontsize)
    Luxor.setline(lineweight)

    for ordered_leaf_node in mg.gprops[:hcl].order
        original_x = mg.vprops[ordered_leaf_node][:x]
        original_y = mg.vprops[ordered_leaf_node][:y]
        polar_radian = original_x * 2 * MathConstants.pi
        adjusted_radius = original_y * max_radius
        x = cos(polar_radian) * adjusted_radius
        y = sin(polar_radian) * adjusted_radius
        Luxor.circle(x, y, mergenodesize, :fill)
        text_x = cos(polar_radian) * (max_radius + 10)
        text_y = sin(polar_radian) * (max_radius + 10)
        Luxor.text(string(ordered_leaf_node), Luxor.Point(text_x, text_y), angle=polar_radian, halign=:left, valign=:middle)
    end
    for (i, (left, right)) in enumerate(eachrow(mg.gprops[:hcl].merges))
        parent_vertex = mg[string(i), :hclust_id]

        original_x = mg.vprops[parent_vertex][:x]
        original_y = mg.vprops[parent_vertex][:y]

        polar_radian = original_x * 2 * MathConstants.pi
        adjusted_radius = original_y * max_radius
        x = cos(polar_radian) * adjusted_radius
        y = sin(polar_radian) * adjusted_radius
        Luxor.circle(x, y, mergenodesize, :fill)

        # left child
        left_child_vertex = mg[string(left), :hclust_id]
        left_child_original_x = mg.vprops[left_child_vertex][:x]
        left_child_original_y = mg.vprops[left_child_vertex][:y]
        left_child_polar_radian = left_child_original_x * 2 * MathConstants.pi
        left_child_adjusted_radius = left_child_original_y * max_radius
        left_child_x = cos(left_child_polar_radian) * left_child_adjusted_radius
        left_child_y = sin(left_child_polar_radian) * left_child_adjusted_radius
        # # draw horizontal bar
        Luxor.arc(Luxor.Point(0, 0), adjusted_radius, left_child_polar_radian, polar_radian, action=:stroke)
        # draw vertical bar
        join_x = cos(left_child_polar_radian) * adjusted_radius
        join_y = sin(left_child_polar_radian) * adjusted_radius
        Luxor.line(Luxor.Point(left_child_x, left_child_y), Luxor.Point(join_x, join_y), action=:stroke)

        # right child
        right_child_vertex = mg[string(right), :hclust_id]
        right_child_original_x = mg.vprops[right_child_vertex][:x]
        right_child_original_y = mg.vprops[right_child_vertex][:y]
        right_child_polar_radian = right_child_original_x * 2 * MathConstants.pi
        right_child_adjusted_radius = right_child_original_y * max_radius
        right_child_x = cos(right_child_polar_radian) * right_child_adjusted_radius
        right_child_y = sin(right_child_polar_radian) * right_child_adjusted_radius
        # # draw horizontal bar
        Luxor.arc(Luxor.Point(0, 0), adjusted_radius, polar_radian, right_child_polar_radian, action=:stroke)
        # draw vertical bar
        join_x = cos(right_child_polar_radian) * adjusted_radius
        join_y = sin(right_child_polar_radian) * adjusted_radius
        Luxor.line(Luxor.Point(right_child_x, right_child_y), Luxor.Point(join_x, join_y), action=:stroke)
    end

    # Finish the drawing and save the file
    Luxor.finish()
    Luxor.preview()
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Performs hierarchical clustering on a distance matrix using Ward's linkage method.

# Arguments
- `distance_matrix::Matrix{<:Real}`: A symmetric distance/dissimilarity matrix

# Returns
- `HierarchicalCluster`: A hierarchical clustering object from Clustering.jl

# Details
Uses Ward's method (minimum variance) for clustering, which:
- Minimizes total within-cluster variance
- Produces compact, spherical clusters
- Works well for visualization in radial layouts
"""
function heirarchically_cluster_distance_matrix(distance_matrix)
    # Perform hierarchical clustering using Ward's method.
    # Ward's method minimizes the total within-cluster variance, which tends to
    # produce compact, spherical clusters. This can be beneficial for visualization,
    # especially in radial layouts, as it can prevent branches from being overly long
    # and overlapping.  Other linkage methods like 'average' and 'complete' are also
    # common choices, however 'single' linkage is prone to chaining effects.
    hcl = Clustering.hclust(distance_matrix, linkage=:ward, branchorder=:optimal)
    # this is equivalent to UPGMA
    # hcl = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)
    return hcl
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determine the optimal number of clusters for a dataset using hierarchical clustering
and silhouette score analysis.

# Arguments
- `distance_matrix`: A symmetric matrix of pairwise distances between data points

# Returns
A tuple containing:
- `hcl`: The fitted hierarchical clustering object
- `optimal_number_of_clusters`: Integer indicating the optimal number of clusters

# Details
The function:
1. Performs hierarchical clustering on the distance matrix
2. Tests cluster counts from 2 to √n (where n is dataset size)
3. Evaluates each clustering using silhouette scores
4. Generates a diagnostic plot showing silhouette scores vs cluster counts
5. Selects the number of clusters that minimizes the silhouette score

The diagnostic plot is displayed automatically during execution.

# See Also
- `heirarchically_cluster_distance_matrix`
- `Clustering.silhouettes`
"""
function identify_optimal_number_of_clusters(distance_matrix)
    hcl = heirarchically_cluster_distance_matrix(distance_matrix)
    ks = 2:Int(ceil(sqrt(length(hcl.labels))))
    silhouette_scores = Float64[]
    ProgressMeter.@showprogress for k in ks
        v = sum(Clustering.silhouettes(Clustering.cutree(hcl, k=k), distance_matrix))
        push!(silhouette_scores, v)
    end
    optimal_number_of_clusters = ks[last(findmin(silhouette_scores))]
    p = StatsPlots.scatter(
        ks,
        silhouette_scores,
        title = "Clustering Performance vs. Number of Clusters\n(lower is better)",
        xlabel = "Number of Clusters",
        ylabel = "Silhouette Score",
        label = "",
        ylims=(0, maximum(silhouette_scores) * 1.1)
    )
    p = StatsPlots.vline!([optimal_number_of_clusters], color=:red, label="inferred optimum = $(optimal_number_of_clusters)")
    return (hcl, optimal_number_of_clusters)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Merge two colors by calculating their minimal color difference vector.

# Arguments
- `c1::Color`: First color input
- `c2::Color`: Second color input 

# Returns
- If colors are equal, returns the input color
- Otherwise returns the color difference vector (c1-c2 or c2-c1) with minimal RGB sum

# Details
Calculates two difference vectors:
- mix_a = c1 - c2 
- mix_b = c2 - c1
Returns the difference vector with the smallest sum of RGB components.
"""
function merge_colors(c1, c2)
    if c1 == c2
        return c1
    else
        mix_a = c1 - c2
        mix_b = c2 - c1
        mix_a_sum = mix_a.r + mix_a.g + mix_a.b
        mix_b_sum = mix_b.r + mix_b.g + mix_b.b
        min_value, min_index = findmin([mix_a_sum, mix_b_sum])
        mixed_color = [mix_a, mix_b][min_index]
        # return Colors.color_names["black"]
        return mixed_color
    end
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

Finds contiguous ranges of `true` values in a boolean vector.

# Arguments
- `bool_vec::AbstractVector{Bool}`: Input boolean vector to analyze
- `min_length=1`: Minimum length requirement for a range to be included

# Returns
Vector of tuples `(start, end)` where each tuple represents the indices of a
contiguous range of `true` values meeting the minimum length requirement.
"""
function find_true_ranges(bool_vec::AbstractVector{Bool}; min_length=1)
    indices = findall(bool_vec)  # Get indices of true values
    if isempty(indices)
    return []  # Handle the case of no true values
    end
    diffs = diff(indices)
    breakpoints = findall(>(1), diffs)  # Find where the difference is greater than 1
    starts = [first(indices); indices[breakpoints .+ 1]]
    ends = [indices[breakpoints]; last(indices)]
    # true_ranges_table = DataFrames.DataFrame(starts = starts, ends = ends, lengths = ends .- starts)
    return collect(Iterators.filter(x -> (x[2] - x[1]) >= min_length, zip(starts, ends)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Sample `n` equally spaced elements from `vector`.

# Arguments
- `vector`: Input vector to sample from
- `n`: Number of samples to return (must be positive)

# Returns
A vector containing `n` equally spaced elements from the input vector.
"""
function equally_spaced_samples(vector, n)
    indices = round.(Int, range(1, length(vector), length=n))
    return vector[indices]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a visualization of chromosome coverage data with statistical thresholds.

# Arguments
- `cdf::DataFrame`: Coverage data frame containing columns:
  - `index`: Chromosome position indices
  - `depth`: Coverage depth values
  - `chromosome`: Chromosome identifier
  - `mean_coverage`: Mean coverage value
  - `std_coverage`: Standard deviation of coverage
  - `3σ`: Boolean vector indicating +3 sigma regions
  - `-3σ`: Boolean vector indicating -3 sigma regions

# Returns
- A StatsPlots plot object showing:
  - Raw coverage data (black line)
  - Mean coverage and ±1,2,3σ thresholds (rainbow colors)
  - Highlighted regions exceeding ±3σ thresholds (red vertical lines)
"""
function chromosome_coverage_table_to_plot(cdf)
    p = StatsPlots.plot(
        xlims = extrema(cdf[!, "index"]),
        ylims=(1, maximum(cdf[!, "depth"]) * 1.1),
        title = cdf[1, "chromosome"],
        xlabel = "chromosome index",
        ylabel = "depth"
    )
    for r in find_true_ranges(cdf[!, "3σ"]; min_length=1000)
        range_mean = Statistics.mean(cdf[r[1]:r[2], :depth])
        StatsPlots.vline!(p,
            [r[1], r[2]],
            # [range_mean, range_mean],
            seriestype = :path,
            label="",
            c=:red,
            linewidth=3,
            alpha=0.1
        )
    end
    for r in find_true_ranges(cdf[!, "-3σ"]; min_length=1000)
        range_mean = Statistics.mean(cdf[r[1]:r[2], :depth])
        StatsPlots.vline!(p,
            [r[1], r[2]],
            # [range_mean, range_mean],
            seriestype = :path,
            label="",
            c=:red,
            linewidth=3,
            alpha=1/3
        )
    end

    color_vec = StatsPlots.cgrad(ColorSchemes.rainbow, 7, categorical = true)
    
    
    
    StatsPlots.plot!(p, equally_spaced_samples(cdf[!, "index"], 10_000), equally_spaced_samples(cdf[!, "depth"], 10_000), label="coverage", c=:black)
    # StatsPlots.plot!(p, equally_spaced_samples(cdf[!, "index"], 10_000), equally_spaced_samples(rolling_centered_avg(cdf[!, "depth"], window_size=1001), 10_000), label="101bp sliding window mean", c=:gray)
    mean_coverage = first(unique(cdf[!, "mean_coverage"]))
    stddev_coverage = first(unique(cdf[!, "std_coverage"]))
    StatsPlots.hline!(p, [mean_coverage + 3 * stddev_coverage], label="+3σ", c=color_vec[7])
    StatsPlots.hline!(p, [mean_coverage + 2 * stddev_coverage], label="+2σ", c=color_vec[6])
    StatsPlots.hline!(p, [mean_coverage + 1 * stddev_coverage], label="+σ", c=color_vec[5])
    StatsPlots.hline!(p, unique(cdf[!, "mean_coverage"]), label="mean_coverage", c=color_vec[4])
    StatsPlots.hline!(p, [mean_coverage + -1 * stddev_coverage], label="-σ", c=color_vec[3])
    StatsPlots.hline!(p, [mean_coverage + -2 * stddev_coverage], label="-2σ", c=color_vec[2])
    StatsPlots.hline!(p, [mean_coverage + -3 * stddev_coverage], label="-3σ", c=color_vec[1])
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute a centered moving average over a vector using a sliding window.

# Arguments
- `data::AbstractVector{T}`: Input vector to be averaged
- `window_size::Int`: Size of the sliding window (odd number recommended)

# Returns
- `Vector{Float64}`: Vector of same length as input containing moving averages

# Details
- For points near the edges, the window is truncated to available data
- Window is centered on each point, using floor(window_size/2) points on each side
- Result type is always Float64 regardless of input type T
"""
function rolling_centered_avg(data::AbstractVector{T}; window_size::Int) where T
    half_window = Int(floor(window_size / 2))
    result = Vector{Float64}(undef, length(data))
    for i in eachindex(data)
        start_idx = max(1, i - half_window)
        end_idx = min(length(data), i + half_window)
        result[i] = Statistics.mean(data[start_idx:end_idx])
    end
    return result
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a multi-panel visualization of genome coverage across chromosomes.

# Arguments
- `coverage_table`: DataFrame containing columns "chromosome" and "coverage" with genomic coverage data

# Returns
- `Plots.Figure`: A composite figure with coverage plots for each chromosome

# Details
Generates one subplot per chromosome, arranged vertically. Each subplot shows the coverage 
distribution across genomic positions for that chromosome.
"""
function visualize_genome_coverage(coverage_table)
    num_plots = length(unique(coverage_table[!, "chromosome"]))
    meta_figure = StatsPlots.plot(
        [chromosome_coverage_table_to_plot(cdf) for cdf in DataFrames.groupby(coverage_table, "chromosome")]...,
        layout = (num_plots, 1),
        size = (800, 600 * num_plots)) # Adjust size as needed
    return meta_figure
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run TransTermHP to predict rho-independent transcription terminators in DNA sequences.

# Arguments
- `fasta`: Path to input FASTA file containing DNA sequences
- `gff`: Optional path to GFF annotation file. If provided, improves prediction accuracy

# Returns
- `String`: Path to output file containing TransTermHP predictions

# Details
- Uses Conda environment 'transtermhp' for execution
- Automatically generates coordinate file from FASTA or GFF input
- Removes temporary coordinate file after completion
- Requires Mycelia's Conda setup
"""
function run_transterm(;fasta, gff="")
    # note in my one test with phage genomes, calling without gff yeilds more hits but average confidence is a bit lower
    if isempty(gff)
        coordinates = generate_transterm_coordinates_from_fasta(fasta)
    else
        coordinates = generate_transterm_coordinates_from_gff(gff)
    end
    transterm_calls_file = replace(coordinates, ".coords" => ".transterm.txt")
    Mycelia.add_bioconda_env("transtermhp")

    conda_base = dirname(dirname(Mycelia.CONDA_RUNNER))
    dat_file = "$(conda_base)/envs/transtermhp/data/expterm.dat"
    # @show isfile(dat_file)
    # @assert = first(readlines(`find $(conda_base) -name "expterm.dat"`))
    # @show dat_file
    @assert isfile(dat_file)
    # conda_base = chomp(read(`$(Mycelia.CONDA_RUNNER) info --base`, String))
    # dat_file = "$(conda_base)/envs/transtermhp/data/expterm.dat"
    # @assert isfile(dat_file)
    run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n transtermhp transterm -p $(dat_file) $(fasta) $(coordinates)`, transterm_calls_file))    
    rm(coordinates)
    return transterm_calls_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate minimal coordinate files required for TransTermHP analysis from FASTA sequences.

Creates artificial gene annotations at sequence boundaries to enable TransTermHP to run
without real gene annotations. For each sequence in the FASTA file, generates two
single-base-pair "genes" at positions 1-2 and (L-1)-L, where L is sequence length.

# Arguments
- `fasta`: Path to input FASTA file containing sequences to analyze

# Returns
- Path to generated coordinate file (original path with ".coords" extension)

# Format
Generated coordinate file follows TransTermHP format:
gene_id start stop chromosome

where chromosome matches FASTA sequence identifiers.

See also: [`run_transterm`](@ref)
"""
function generate_transterm_coordinates_from_fasta(fasta)
    # 10. USING TRANSTERM WITHOUT GENOME ANNOTATIONS

    # TransTermHP uses known gene information for only 3 things: (1) tagging the
    # putative terminators as either "inside genes" or "intergenic," (2) choosing the
    # background GC-content percentage to compute the scores, because genes often
    # have different GC content than the intergenic regions, and (3) producing
    # slightly more readable output. Items (1) and (3) are not really necessary, and
    # (2) has no effect if your genes have about the same GC-content as your
    # intergenic regions.

    # Unfortunately, TransTermHP doesn't yet have a simple option to run without an
    # annotation file (either .ptt or .coords), and requires at least 2 genes to be
    # present. The solution is to create fake, small genes that flank each
    # chromosome. To do this, make a fake.coords file that contains only these two
    # lines:

    # 	fakegene1	1 2	chome_id
    # 	fakegene2	L-1 L	chrom_id

    # where L is the length of the input sequence and L-1 is 1 less than the length
    # of the input sequence. "chrom_id" should be the word directly following the ">"
    # in the .fasta file containing your sequence. (If, for example, your .fasta file
    # began with ">seq1", then chrom_id = seq1).

    # This creates a "fake" annotation with two 1-base-long genes flanking the
    # sequence in a tail-to-tail arrangement: --> <--. TransTermHP can then be run
    # with:

    # 	transterm -p expterm.dat sequence.fasta fake.coords

    # If the G/C content of your intergenic regions is about the same as your genes,
    # then this won't have too much of an effect on the scores terminators receive.
    # On the other hand, this use of TransTermHP hasn't been tested much at all, so
    # it's hard to vouch for its accuracy.

    coords_table = DataFrames.DataFrame(
        gene_id = String[],
        start = Int[],
        stop = Int[],
        chromosome = String[]
    )
    for record in Mycelia.open_fastx(fasta)
        row = (
            gene_id = FASTX.identifier(record) * "_start",
            start = 1,
            stop = 2,
            chromosome = FASTX.identifier(record)
        )
        push!(coords_table, row)

        row = (
            gene_id = FASTX.identifier(record) * "_stop",
            start = length(FASTX.sequence(record))-1,
            stop = length(FASTX.sequence(record)),
            chromosome = FASTX.identifier(record)
        )
        push!(coords_table, row)
    end
    transterm_coordinates_file = fasta * ".coords"
    uCSV.write(transterm_coordinates_file, data = collect(DataFrames.eachcol(coords_table)), delim="  ")
    return transterm_coordinates_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GFF file to a coordinates file compatible with TransTermHP format.

# Arguments
- `gff_file::String`: Path to input GFF file

# Processing
- Converts 1-based to 0-based coordinates
- Extracts gene IDs from the attributes field
- Retains columns: gene_id, start, end, seqid

# Returns
- Path to the generated coordinates file (original filename with '.coords' suffix)

# Output Format
Space-delimited file with columns: gene_id, start, end, seqid
"""
function generate_transterm_coordinates_from_gff(gff_file)
    raw_gff = Mycelia.read_gff(gff_file)
    # switch start to be zero index by subtracting one
    raw_gff[!, "start"] = raw_gff[!, "start"] .- 1
    raw_gff[!, "gene_id"] = last.(split.(first.(split.(raw_gff[!, "attributes"], ";")), '='))
    raw_gff = raw_gff[!, ["gene_id", "start", "end", "#seqid"]]
    transterm_coordinates_file = gff_file * ".coords"
    uCSV.write(transterm_coordinates_file, data = collect(DataFrames.eachcol(raw_gff)), delim="  ")
    return transterm_coordinates_file
end

# function normalize_fasta(fasta_file, outdir)
#     mkpath("$(outdir)/normalized_fasta")
#     normalized_fasta_file = "$(outdir)/normalized_fasta/$(basename(fasta_file))"
#     if !isfile(normalized_fasta_file)
#         fasta_in = FASTX.FASTA.Reader(open(fasta_file))
#         fasta_out = FASTX.FASTA.Writer(open(normalized_fasta_file, "w"))
#         for (i, record) in enumerate(fasta_in)
#             updated_record = FASTX.FASTA.Record("$(i)", FASTX.FASTA.sequence(record))
#             write(fasta_out, updated_record)
#         end
#         close(fasta_in)
#         close(fasta_out)
#     end
#     return normalized_fasta_file
# end

# function run_prokka(ID, OUT_DIR, normalized_fasta_file)
#     prokka_dir="$(OUT_DIR)/prokka"
#     if !isdir(prokka_dir)
#         mkdir(prokka_dir)
#     end
#     prokka_cmd = `prokka --force --cpus 1 --outdir $(prokka_dir) --prefix $(ID) $(normalized_fasta_file)`
#     run(pipeline(prokka_cmd, stdout="$(prokka_dir)/prokka.out"))
#     return prokka_dir
# end

# function run_mlst(ID, OUT_DIR, normalized_fasta_file)
#     mlst_dir="$(OUT_DIR)/mlst"
#     if !isdir(mlst_dir)
#         mkdir(mlst_dir)
#     end
#     p = pipeline(
#             `mlst $(normalized_fasta_file)`,
#             stdout="$(mlst_dir)/$(ID).mlst.out")
#     run(p)
#     return mlst_dir
# end

# function run_phispy(ID, OUT_DIR, prokka_dir)
#     # 1 	prophage_coordinates.tsv
#     # 2 	GenBank format output
#     # 4 	prophage and bacterial sequences
#     # 8 	prophage_information.tsv
#     # 16 	prophage.tsv
#     # 32 	GFF3 format
#     # 64 	prophage.tbl
#     # 128 	test data used in the random forest
#     # 255   for all of them

#     phispy_dir="$(OUT_DIR)/phispy"
#     if !isdir(phispy_dir)
#         mkdir(phispy_dir)
#     end
#     if isempty(readdir(phispy_dir))
#         phisphy_cmd = `PhiSpy.py $(prokka_dir)/$(ID).gbk --output_dir $(phispy_dir) --file_prefix $(ID)-phispy --output_choice 255`
#         try
#             run(pipeline(phisphy_cmd, stdout="$(phispy_dir)/phisphy.out"))
#         catch
#             if isfile("$(phispy_dir)/$(ID)-phispy_prophage.gff3")
#                 @warn "phispy errored out after prophage gff3 was written"
#             else
#                 @error "phispy prophage gff3 not written"
#             end
#         end
#     end
#     return phispy_dir
# end

# function run_trnascan(ID, out_dir, normalized_fasta_file)
#     trnascan_dir = "$(out_dir)/trnascan"
#     # trnascan doesn't like to overwrite existing things
#     if !isdir(trnascan_dir)
#         mkdir(trnascan_dir)
#     end
#     if isempty(readdir(trnascan_dir))

#         #     -B for using Bacterial
#         trnascan_cmd = 
#         `tRNAscan-SE 
#             -B 
#             --output $(trnascan_dir)/$(ID).trnascan.out 
#             --bed $(trnascan_dir)/$(ID).trnascan.bed 
#             --fasta $(trnascan_dir)/$(ID).trnascan.fasta 
#             --struct $(trnascan_dir)/$(ID).trnascan.struct
#             --stats $(trnascan_dir)/$(ID).trnascan.stats 
#             --log $(trnascan_dir)/$(ID).trnascan.log
#             $(normalized_fasta_file)`
#         run(pipeline(trnascan_cmd, stdout="$(trnascan_dir)/trnascan.out", stderr="$(trnascan_dir)/trnascan.out"))
#     end
#     return trnascan_dir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run tRNAscan-SE to identify and annotate transfer RNA genes in the provided sequence file.

# Arguments
- `fna_file::String`: Path to input FASTA nucleotide file
- `outdir::String`: Output directory path (default: input_file_path + "_trnascan")

# Returns
- `String`: Path to the output directory containing tRNAscan-SE results

# Output Files
Creates the following files in `outdir`:
- `*.trnascan.out`: Main output with tRNA predictions
- `*.trnascan.bed`: BED format coordinates
- `*.trnascan.fasta`: FASTA sequences of predicted tRNAs
- `*.trnascan.struct`: Secondary structure predictions
- `*.trnascan.stats`: Summary statistics
- `*.trnascan.log`: Program execution log

# Notes
- Uses the general tRNA model (-G flag) suitable for all domains of life
- Automatically sets up tRNAscan-SE via Bioconda
- Skips processing if output directory contains files
"""
function run_trnascan(;fna_file, outdir=fna_file * "_trnascan")
    Mycelia.add_bioconda_env("trnascan-se")
    # trnascan doesn't like to overwrite existing things
    if !isdir(outdir)
        mkdir(outdir)
    else
        @info "$(outdir) already exists"
    end
    ID = basename(fna_file)
    if isempty(readdir(outdir))
        # -G : use general tRNA model (cytoslic tRNAs from all 3 domains included)
        #     -B for using Bacterial
        trnascan_cmd =
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n trnascan-se tRNAscan-SE
            -G
            --output $(outdir)/$(ID).trnascan.out
            --bed $(outdir)/$(ID).trnascan.bed
            --fasta $(outdir)/$(ID).trnascan.fasta
            --struct $(outdir)/$(ID).trnascan.struct
            --stats $(outdir)/$(ID).trnascan.stats
            --log $(outdir)/$(ID).trnascan.log
            --prefix
            --progress
            $(fna_file)`
        # run(pipeline(trnascan_cmd, stdout="$(outdir)/$(ID).trnascan.out", stderr="$(outdir)/$(ID).trnascan.out"))
        run(trnascan_cmd)
    end
    return outdir
end

# function run_counterselection_spacer_detection(strain, out_dir, normalized_fasta_file)
#     counter_selection_dir = "$(out_dir)/counter-selection"
#     if !isdir(counter_selection_dir)
#         mkdir(counter_selection_dir)
#     end

#     if isempty(readdir(counter_selection_dir))
#         regex = BioSequences.biore"TTT[CG][ACGT]{25}"dna
#         k = 29
#         KMER_TYPE = BioSequences.BigDNAMer{k}

#         spacer_table = DataFrames.DataFrame(
#             ID = [],
#             contig = [],
#             PAM_and_spacer = [],
#             spacer = [],
#             strand = [],
#             start = [],
#             stop = [],
# #             free_energy = [],
# #             visualization_url = []
#         )

#         ProgressMeter.@showprogress for record in collect(FASTX.FASTA.Reader(open(normalized_fasta_file)))
#             for (i, kmer, reverse_complement_kmer) in BioSequences.each(KMER_TYPE, FASTX.FASTA.sequence(record))
#                 strand = missing
#                 if occursin(regex, kmer)
#                     strand = "+"
#                 elseif occursin(regex, reverse_complement_kmer)
#                     strand = "-"
#                     kmer = reverse_complement_kmer
#                 end
#                 if !ismissing(strand)
#                     spacer = BioSequences.DNAMer(kmer[i] for i in 5:length(kmer)) 
# #                     RNAfold_output = read(pipeline(`echo "$(string(spacer))"`, `RNAfold --noLP`), String)

# #                     rna_sequence, structure, free_energy = match(r"([ACGU]{25})\n([.()]{25})\s\(\s*(.*?)\)", RNAfold_output).captures
# #                     url = "http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence=$(rna_sequence)&structure=$(structure)"

#                     kmer_range = i:i+k-1

#                     t = DataFrames.DataFrame(
#                         ID = ID,
#                         contig = FASTX.FASTA.identifier(record),
#                         PAM_and_spacer = kmer,
#                         spacer = spacer,
#                         strand = strand,
#                         start = i,
#                         stop = i+k-1,
# #                         free_energy = free_energy,
# #                         visualization_url = url
#                     )
#                     spacer_table = vcat(spacer_table, t)
#                 end
#             end
#         end


#         uCSV.write(
#             "$(counter_selection_dir)/$(ID)-cpf1-spacers.tsv",
#             delim='\t',
#             data = collect(DataFrames.eachcol(spacer_table)),
#             header = DataFrames.names(spacer_table)
#         )
#         if isfile("$(out_dir)/rna.ps")
#             rm("$(out_dir)/rna.ps")
#         end
#     end
#     return counter_selection_dir
# end



# function run_amrfinderplus(ID, out_dir, protein_fasta)
#     amrfinderplus_dir = "$(out_dir)/amrfinderplus"
#     if !isdir(amrfinderplus_dir)
#         mkdir(amrfinderplus_dir)
#     end

#     if isempty(readdir(amrfinderplus_dir))
# #         run(`amrfinder -u`)
        
#         # because the pipeline is set up with some hacky CONDA path manipulation, 
#         # explictly setting amrfinder directory path to the location in the docker host
#         amrfinder_db_path = get(ENV, "AMRFINDER_DB", "none")

#         if amrfinder_db_path != "none"
#             cmd = 
#             `amrfinder
#             -p $(protein_fasta)
#             --plus
#             --output $(amrfinderplus_dir)/$(ID).amrfinderplus.tsv
#             -d $(amrfinder_db_path)
#             `
#         else
#             cmd = 
#             `amrfinder
#             -p $(protein_fasta)
#             --plus
#             --output $(amrfinderplus_dir)/$(ID).amrfinderplus.tsv
#             `
#         end
            

#         p = pipeline(cmd, 
#                 stdout="$(amrfinderplus_dir)/$(ID).amrfinderplus.out",
#                 stderr="$(amrfinderplus_dir)/$(ID).amrfinderplus.err")
#         run(p)
#     end
#     return amrfinderplus_dir
# end

# function make_diamond_db(fasta_file, db_file=fasta_file)
#     @time run(`diamond makedb --in $(fasta_file) -d $(db_file)`)
# end

# # in order to change this to be a standard blast where we don't need all pairwise hits
# # just drop the parameters id, min-score, max-target-seqs
# function pairwise_diamond(joint_fasta_file)
#     if !isfile("$(joint_fasta_file).dmnd")
#         make_diamond_db(joint_fasta_file)
#     end
#     n_records = count_records(joint_fasta_file)
#     # max_target_seqs = Int(ceil(sqrt(n_records)))
#     # @show "here!"
#     sensitivity = "--iterate"
#     # --block-size/-b
#     # https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options
#     # set block size to total memory / 8
#     available_gigabytes = floor(Sys.free_memory() / 1e9)
#     block_size = floor(available_gigabytes / 8)
    
#     @time run(`diamond blastp $(sensitivity) --block-size $(block_size) --id 0 --min-score 0 --max-target-seqs $(n_records) --unal 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore -d $(joint_fasta_file).dmnd -q $(joint_fasta_file) -o $(joint_fasta_file).dmnd.tsv`)
#     # # pairwise output is all of the alignments, super helpful!
#     # # @time run(`diamond blastp $(sensitivity) --id 0 --min-score 0 --max-target-seqs $(N_RECORDS) --unal 1 --outfmt 0  -d $(joint_fasta_outfile).dmnd -q $(joint_fasta_outfile) -o $(joint_fasta_outfile).diamond.pairwise.txt`)
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run diamond search, returns path to diamond results.

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_diamond(;
#         identifier,
#         out_dir,
#         protein_fasta,
#         diamond_db,
#         force=false,
#         outfile="$(identifier).prodigal.faa.diamond.txt"
#     )
#     diamond_dir = mkpath("$(out_dir)/diamond")

#     # http://www.diamondsearch.org/index.php?pages/command_line_options/
#     # --block-size/-b #Block size in billions of sequence letters to be processed at a time.  
#     #     This is the main pa-rameter for controlling the program’s memory usage.  
#     #     Bigger numbers will increase the useof memory and temporary disk space, but also improve performance.  
#     #     The program can beexpected to use roughly six times this number of memory (in GB). So for the default value of-b2.0, 
#     #     the memory usage will be about 12 GB
#     system_memory_in_gigabytes = Int(Sys.total_memory()) / 1e9
#     # reference says 6 but let's round upwards towards 8
#     gb_per_block = 8
#     block_size = system_memory_in_gigabytes / gb_per_block
    
#     outfile = "$(diamond_dir)/$(outfile)"
    
#     if force || !isfile(outfile)
#         cmd = 
#         `diamond blastp
#         --threads $(Sys.CPU_THREADS)
#         --block-size $(block_size)
#         --db $(diamond_db)
#         --query $(protein_fasta)
#         --out $(outfile)
#         --evalue 0.001
#         --iterate
#         --outfmt 6 qseqid qtitle qlen sseqid sallseqid stitle salltitles slen qstart qend sstart send evalue bitscore length pident nident mismatch staxids
#         `

#         # --un                     file for unaligned queries
#         # --al                     file or aligned queries
#         # --unfmt                  format of unaligned query file (fasta/fastq)
#         # --alfmt                  format of aligned query file (fasta/fastq)
#         # --unal                   report unaligned queries (0=no, 1=yes)

# #         Value 6 may be followed by a space-separated list of these keywords:

# #         qseqid means Query Seq - id
# #         qtitle means Query title
# #         qlen means Query sequence length
# #         sseqid means Subject Seq - id
# #         sallseqid means All subject Seq - id(s), separated by a ';'
# #         stitle means Subject Title
# #         salltitles means All Subject Title(s), separated by a '<>'
# #         slen means Subject sequence length
# #         qstart means Start of alignment in query
# #         qend means End of alignment in query
# #         sstart means Start of alignment in subject
# #         send means End of alignment in subject
# #         evalue means Expect value
# #         bitscore means Bit score
# #         length means Alignment length
# #         pident means Percentage of identical matches
# #         nident means Number of identical matches
# #         mismatch means Number of mismatches
# #         staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
        
#         @time run(pipeline(cmd))
#     end
#     return outfile
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_mmseqs_easy_taxonomy(;out_dir, query_fasta, target_database, outfile, force=false)
#     add_bioconda_env("mmseqs2")
#     out_dir = mkpath(joinpath(out_dir, "mmseqs_easy_taxonomy"))
#     outfile = joinpath(out_dir, outfile * ".mmseqs_easy_taxonomy." * basename(target_database) * ".txt")
#     # note I tried adjusting all of the following, and none of them improved overall runtime
#     # in any meaningful way
#     # -s FLOAT                         Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
#     # https://github.com/soedinglab/MMseqs2/issues/577#issuecomment-1191584081
#     # apparently orf-filter 1 speeds up by 50%!
#     # --orf-filter INT                 Prefilter query ORFs with non-selective search
#     #                               Only used during nucleotide-vs-protein classification
#     #                               NOTE: Consider disabling when classifying short reads [0]
#     # --lca-mode INT                   LCA Mode 1: single search LCA , 2/3: approximate 2bLCA, 4: top hit [3]
#     # --lca-search BOOL                Efficient search for LCA candidates [0]
#     # ^ this looks like it actually runs @ 1 with s=1.0 & --orf-filter=1
    
#     # 112 days to process 600 samples at this rate....
#     # 278 minutes or 4.5 hours for a single sample classification!!
#     # 16688.050696 seconds (1.43 M allocations: 80.966 MiB, 0.02% gc time, 0.00% compilation time)
#     # this is for default parameters
#     # lowering sensitivity and taking LCA
#     # 16590.725343 seconds (1.06 M allocations: 53.487 MiB, 0.00% compilation time)
#     # took just as long!
#     # difference was only 10 minutes
#     # 15903.218456 seconds (969.92 k allocations: 48.624 MiB, 0.01% gc time)
#     # use default parameters
    
#     if force || (!force && !isfile(outfile))
#         cmd = 
#         `$(CONDA_RUNNER) run --no-capture-output -n mmseqs2 mmseqs
#          easy-taxonomy
#          $(query_fasta)
#          $(target_database)
#          $(outfile)
#          $(joinpath(out_dir, "tmp"))
#         `
#         @time run(pipeline(cmd))
#     else
#         @info "target outfile $(outfile) already exists, remove it or set force=true to re-generate"
#     end
#     return outfile
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Runs the MMseqs2 easy-search command on the given query FASTA file against the target database.

# Arguments
- `query_fasta::String`: Path to the query FASTA file.
- `target_database::String`: Path to the target database.
- `out_dir::String`: Directory to store the output file. Defaults to the directory of the query FASTA file.
- `outfile::String`: Name of the output file. Defaults to a combination of the query FASTA and target database filenames.
- `format_output::String`: Format of the output. Defaults to a predefined set of fields.
- `threads::Int`: Number of CPU threads to use. Defaults to the number of CPU threads available.
- `force::Bool`: If true, forces the re-generation of the output file even if it already exists. Defaults to false.

# Returns
- `outfile_path::String`: Path to the generated output file.

# Notes
- Adds the `mmseqs2` environment using Bioconda if not already present.
- Removes temporary files created during the process.
"""
function run_mmseqs_easy_search(;
        query_fasta,
        target_database,
        out_dir=dirname(query_fasta),
        outfile=basename(query_fasta) * ".mmseqs_easy_search." * basename(target_database) * ".txt",
        format_output = "query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,taxid,taxname",
        threads = Sys.CPU_THREADS,
        force=false)
    
    add_bioconda_env("mmseqs2")
    outfile_path = joinpath(out_dir, outfile)
    tmp_dir = joinpath(out_dir, "tmp")
    if force || (!force && !isfile(outfile_path))
        cmd = 
        `$(CONDA_RUNNER) run --no-capture-output -n mmseqs2 mmseqs
            easy-search
            $(query_fasta)
            $(target_database)
            $(outfile_path)
            $(tmp_dir)
            --threads $(threads)
            --format-mode 4
            --format-output $(format_output)
            --start-sens 1 -s 7 --sens-steps 7
            --sort-results 1
            --remove-tmp-files 1
        `
        @time run(pipeline(cmd))
    else
        @info "target outfile $(outfile_path) already exists, remove it or set force=true to re-generate"
    end
    # we set remote tmp files = 1 above, but it still doesn't seem to work?
    if isdir(tmp_dir)
        rm(tmp_dir, recursive=true)
    end
    return outfile_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)


Run the BLASTN (Basic Local Alignment Search Tool for Nucleotides) command with specified parameters.

# Arguments
- `out_dir::String`: The output directory where the BLASTN results will be saved.
- `fasta::String`: The path to the input FASTA file containing the query sequences.
- `blast_db::String`: The path to the BLAST database to search against.
- `task::String`: The BLASTN task to perform. Default is "megablast".
- `force::Bool`: If true, forces the BLASTN command to run even if the output file already exists. Default is false.
- `remote::Bool`: If true, runs the BLASTN command remotely. Default is false.
- `wait::Bool`: If true, waits for the BLASTN command to complete before returning. Default is true.

# Returns
- `outfile::String`: The path to the output file containing the BLASTN results.

# Description
This function constructs and runs a BLASTN command based on the provided parameters.
It creates an output directory if it doesn't exist, constructs the output file path, and checks if the BLASTN command needs to be run based on the existence and size of the output file.
The function supports running the BLASTN command locally or remotely, with options to force re-running and to wait for completion.
"""
function run_blastn(;out_dir, fasta, blast_db, task="megablast", force=false, remote=false, wait=true)
    blast_dir = mkpath(joinpath(out_dir, "blastn"))
    outfile = "$(blast_dir)/$(basename(fasta)).blastn.$(basename(blast_db)).$(task).txt"
    # if remote
        # outfile = replace(outfile, ".txt" => ".remote.txt")
    # end
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    
    # I want to speed this up more but don't know how
    # 
    # num_alignments
    # https://www.ncbi.nlm.nih.gov/books/NBK569845/
    # Windowmasker masks the over-represented sequence data and it can also mask the low complexity sequence data using the built-in dust algorithm (through the -dust option). To mask low-complexity sequences only, we will need to use dustmasker.
    # http://ftp.ncbi.nlm.nih.gov/pub/agarwala/dustmasker/README.dustmasker
    # http://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/README.windowmasker
    # ./windowmasker -ustat ustat.15 -in chr1.fa -out chr1.wm -dust true
    
    # windowmasker -in hs_chr -infmt blastdb -mk_counts -parse_seqids -out hs_chr_mask.counts
    # windowmasker -in hs_chr -infmt blastdb -ustat hs_chr_mask.count -outfmt maskinfo_asn1_bin -parse_seqids -out hs_chr_mask.asnb
    
    # dustmasker -in hs_chr -infmt blastdb -parse_seqids -outfmt maskinfo_asn1_bin -out hs_chr_dust.asnb
    
    # makeblastdb -in hs_chr –input_type blastdb -dbtype nucl -parse_seqids -mask_data hs_chr_mask.asnb -out hs_chr -title "Human Chromosome, Ref B37.1"
    # blastdbcmd -db hs_chr -info
    
    # https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/winmasker/README

    # windowmasker -mk_counts -checkdup true -infmt blastdb -in nt -out nt.wm.1
    # windowmasker -ustat nt.wm.1 -dust true -in nt -infmt blastdb -out nt.wm.2
    # makeblastdb -in nt –input_type blastdb -dbtype nucl -parse_seqids -mask_data nt.wm.2 -out nt_masked_deduped -title "nt masked and deduped"
    
    # windowmasker -mk_counts -infmt blastdb -in nt -out nt.wm.no_dup_check.1
    # windowmasker -ustat nt.wm.no_dup_check.1 -dust true -in nt -infmt blastdb -out nt.wm.no_dup_check.2
    # makeblastdb -in nt –input_type blastdb -dbtype nucl -parse_seqids -mask_data nt.wm.no_dup_check.2 -out nt_masked -title "nt masked"

    # windowmasker -convert -in input_file_name -out output_file_name [-sformat output_format] [-smem available_memory]
    
    if force || need_to_run
        # if remote
        #     cmd = 
        #         `
        #         blastn
        #         -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        #         -query $(fasta)
        #         -db $(basename(blast_db))
        #         -out $(outfile)
        #         -max_target_seqs 10
        #         -evalue 0.001
        #         -task $(task)
        #         -soft_masking true
        #         -subject_besthit
        #         -dust
        #         -remote
        #         `
        # else
        # https://www.ncbi.nlm.nih.gov/books/NBK571452/
        # cap @ 8 and also use -mt_mode = 1 based on empirical evidence from
        # above blog post
        cmd = 
        `
        blastn
        -num_threads $(min(Sys.CPU_THREADS, 8))
        -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        -query $(fasta)
        -db $(blast_db)
        -out $(outfile)
        -max_target_seqs 10
        -subject_besthit
        -task $(task)
        -evalue 0.001
        `
        # end
#         p = pipeline(cmd, 
#                 stdout="$(blastn_dir)/$(ID).blastn.out",
#                 stderr="$(blastn_dir)/$(ID).blastn.err")
        @info "running cmd $(cmd)"
        @time run(pipeline(cmd), wait=wait)
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run a BLAST (Basic Local Alignment Search Tool) command with the specified parameters.

# Arguments
- `out_dir::String`: The output directory where the BLAST results will be stored.
- `fasta::String`: The path to the input FASTA file.
- `blast_db::String`: The path to the BLAST database.
- `blast_command::String`: The BLAST command to be executed (e.g., `blastn`, `blastp`).
- `force::Bool`: If `true`, forces the BLAST command to run even if the output file already exists. Default is `false`.
- `remote::Bool`: If `true`, runs the BLAST command remotely. Default is `false`.
- `wait::Bool`: If `true`, waits for the BLAST command to complete before returning. Default is `true`.

# Returns
- `outfile::String`: The path to the output file containing the BLAST results.

# Description
This function constructs and runs a BLAST command based on the provided parameters. It creates the necessary output directory, constructs the output file name, and determines whether the BLAST command needs to be run based on the existence and size of the output file. The function supports both local and remote execution of the BLAST command.

If `force` is set to `true` or the output file does not exist or is empty, the BLAST command is executed. The function logs the command being run and measures the time taken for execution. The output file path is returned upon completion.
"""
function run_blast(;out_dir, fasta, blast_db, blast_command, force=false, remote=false, wait=true)
    blast_dir = mkpath(joinpath(out_dir, blast_command))
    outfile = "$(blast_dir)/$(basename(fasta)).$(blast_command).$(basename(blast_db)).txt"
    if remote
        outfile = replace(outfile, ".txt" => ".remote.txt")
    end
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    if force || need_to_run
        if remote
            cmd = 
                `
                $(blast_command)
                -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
                -query $(fasta)
                -db $(basename(blast_db))
                -out $(outfile)
                -max_target_seqs 10
                -evalue 0.001
                -remote
                `
        else
            cmd = 
            `
            $(blast_command)
            -num_threads $(Sys.CPU_THREADS)
            -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
            -query $(fasta)
            -db $(blast_db)
            -out $(outfile)
            -max_target_seqs 10
            -evalue 0.001
            `
        end
#         p = pipeline(cmd, 
#                 stdout="$(blastn_dir)/$(ID).blastn.out",
#                 stderr="$(blastn_dir)/$(ID).blastn.err")
        @info "running cmd $(cmd)"
        @time run(pipeline(cmd), wait=wait)
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Prodigal gene prediction software on input FASTA file to identify protein-coding genes
in metagenomes or single genomes.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing genomic sequences
- `out_dir::String=dirname(fasta_file)`: Directory for output files. Defaults to input file's directory

# Returns
Named tuple containing paths to all output files:
- `fasta_file`: Input FASTA file path
- `out_dir`: Output directory path  
- `gff`: Path to GFF format gene predictions
- `gene_scores`: Path to all potential genes and their scores
- `fna`: Path to nucleotide sequences of predicted genes
- `faa`: Path to protein translations of predicted genes
- `std_out`: Path to captured stdout
- `std_err`: Path to captured stderr
"""
function run_prodigal(;fasta_file, out_dir=dirname(fasta_file))
    
    # if isempty(out_dir)
    #     prodigal_dir = mkpath("$(fasta_file)_prodigal")
    # else
    #     prodigal_dir = mkpath(out_dir)
    # end

    # $ prodigal
    # -------------------------------------
    # PRODIGAL v2.6.3 [February, 2016]         
    # Univ of Tenn / Oak Ridge National Lab
    # Doug Hyatt, Loren Hauser, et al.     
    # -------------------------------------

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    gff = "$(out_dir)/$(basename(fasta_file)).prodigal.gff"
    faa = "$(out_dir)/$(basename(fasta_file)).prodigal.faa"
    fna = "$(out_dir)/$(basename(fasta_file)).prodigal.fna"
    gene_scores = "$(out_dir)/$(basename(fasta_file)).prodigal.all_potential_gene_scores.txt"
    std_out = "$(out_dir)/$(basename(fasta_file)).prodigal.out"
    std_err = "$(out_dir)/$(basename(fasta_file)).prodigal.err"
    
    # I usually delete the rest, so don't reprocess if outputs of interest are present
    if (!isfile(gff) && !isfile(faa))
        add_bioconda_env("prodigal")
        cmd = 
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n prodigal prodigal
        -f gff
        -m
        -p meta
        -o $(gff)
        -i $(fasta_file)
        -a $(faa)
        -d $(fna)
        -s $(gene_scores)
        `
        p = pipeline(cmd, stdout=std_out, stderr=std_err)
        run(p)
    end
    return (;fasta_file, out_dir, gff, gene_scores, fna, faa, std_out, std_err)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Pyrodigal gene prediction on a FASTA file using the meta procedure optimized for metagenomic sequences.

Pyrodigal is a reimplementation of the Prodigal gene finder, which identifies protein-coding sequences in bacterial and archaeal genomes.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing genomic sequences
- `out_dir::String`: Output directory path (default: input filename + "_pyrodigal")

# Returns
Named tuple containing:
- `fasta_file`: Input FASTA file path
- `out_dir`: Output directory path
- `gff`: Path to GFF output file with gene predictions
- `faa`: Path to FASTA file with predicted protein sequences 
- `fna`: Path to FASTA file with nucleotide sequences

# Notes
- Uses metagenomic mode (`-p meta`) optimized for mixed communities
- Masks runs of N nucleotides (`-m` flag)
- Minimum gene length set to 33bp
- Maximum overlap between genes set to 31bp
- Requires Pyrodigal to be available in a Conda environment
- Skips processing if output files already exist
"""
function run_pyrodigal(;fasta_file, out_dir=fasta_file * "_pyrodigal")
    # https://pyrodigal.readthedocs.io/en/stable/guide/cli.html#command-line-interface

    # -a trans_file         Write protein translations to the selected file.
    # -c                    Closed ends. Do not allow genes to run off edges.
    # -d nuc_file           Write nucleotide sequences of genes to the selected file.
    # -f output_type        Select output format.
    # -g tr_table           Specify a translation table to use.
    # -i input_file         Specify FASTA input file.
    # -m                    Treat runs of N as masked sequence and don't build genes across them.
    # -n                    Bypass Shine-Dalgarno trainer and force a full motif scan.
    # -o output_file        Specify output file.
    # -p mode               Select procedure.
    # -s start_file         Write all potential genes (with scores) to the selected file.
    # -t training_file      Write a training file (if none exists); otherwise, read and use the specified training file.
    # -j jobs, --jobs jobs           The number of threads to use if input contains multiple sequences.
    # --min-gene MIN_GENE            The minimum gene length.
    # --min-edge-gene MIN_EDGE_GENE  The minimum edge gene length.
    # --max-overlap MAX_OVERLAP      The maximum number of nucleotides that can overlap between two genes on the same strand.
    #                             This must be lower or equal to the minimum gene length.
    # --no-stop-codon                Disables translation of stop codons into star characters (*) for complete genes.
    # --pool {thread,process}        The sort of pool to use to process genomes in parallel. Processes may be faster than
    #                             threads on some machines, refer to documentation. (default: thread)

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    gff = "$(out_dir)/$(basename(fasta_file)).pyrodigal.gff"
    faa = "$(out_dir)/$(basename(fasta_file)).pyrodigal.faa"
    fna = "$(out_dir)/$(basename(fasta_file)).pyrodigal.fna"
    # -s $(gene_scores)
    # gene_scores = "$(out_dir)/$(basename(fasta_file)).pyrodigal.all_potential_gene_scores.txt"
    std_out = "$(out_dir)/$(basename(fasta_file)).pyrodigal.out"
    std_err = "$(out_dir)/$(basename(fasta_file)).pyrodigal.err"
    mkpath(out_dir)
    
    # I usually delete the rest, so don't reprocess if outputs of interest are present
    # `max_overlap` must be lower than `min_gene`
    if (!isfile(gff) || !isfile(faa))
        Mycelia.add_bioconda_env("pyrodigal")
        cmd = 
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n pyrodigal pyrodigal
        -f gff
        -m
        -p meta
        -o $(gff)
        -i $(fasta_file)
        -a $(faa)
        -d $(fna)
        --min-gene 33
        --max-overlap 31
        `
        p = pipeline(cmd, stdout=std_out, stderr=std_err)
        run(p)
        cp(fasta_file, "$(out_dir)/$(basename(fasta_file))", force=true)
        (filesize(std_out) == 0) && rm(std_out)
        (filesize(std_err) == 0) && rm(std_err)
    end
    # return (;fasta_file, out_dir, gff, gene_scores, fna, faa, std_out, std_err)
    return (;fasta_file, out_dir, gff, faa, fna)
end

# ```
# https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/

# GFF3 has 9 required fields, though not all are utilized (either blank or a default value of ‘.’).

#     Sequence ID
#     Source
#         Describes the algorithm or the procedure that generated this feature. Typically Genescane or Genebank, respectively.
#     Feature Type
#         Describes what the feature is (mRNA, domain, exon, etc.).
#         These terms are constrained to the [Sequence Ontology terms](http://www.sequenceontology.org/).
#     Feature Start
#     Feature End
#     Score
#         Typically E-values for sequence similarity and P-values for predictions.
#     Strand
#     Phase
#         Indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
#     Atributes
#         A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent . You can see the full list [here](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

# ```

# """
#     create_chromosome_genedata_table(chromosome)

# Take a chromosome from GenomicAnnotations.jl in GFF (and possibly genbank)
# and return the formed dataframe.
# """
# function create_chromosome_genedata_table(chromosome)
#     # genedata is already provided as a dataframe with all of the Attributes as columns
#     table = copy(chromosome.genedata)
    
#     # the rest of these need to be created
#     # I'm inserting them directly into their correct locations
#     DataFrames.insertcols!(table, 1, "sequence-id" => fill(chromosome.name, DataFrames.nrow(table)))
    
#     DataFrames.insertcols!(table, 3, "feature" => GenomicAnnotations.feature.(chromosome.genes))

#     loci = getproperty.(GenomicAnnotations.locus.(chromosome.genes), :position)
#     DataFrames.insertcols!(table, 4, "start" => first.(loci))
#     DataFrames.insertcols!(table, 5, "stop" => last.(loci))
    
#     DataFrames.insertcols!(table, 7, "strand" => .!GenomicAnnotations.iscomplement.(chromosome.genes))        
    
#     return table
# end

# function create_chromosome_genedata_table(chromosome)
#     table = chromosome.genedata
#     table[!, "chromosome"] .= chromosome.name
#     table[!, "feature"] = GenomicAnnotations.feature.(chromosome.genes)
#     table[!, "strand"] = .!GenomicAnnotations.iscomplement.(chromosome.genes)
#     table[!, "locus"] = getproperty.(GenomicAnnotations.locus.(chromosome.genes), :position)
#     table[!, "start"] = first.(chromosome.genedata[!, "locus"])
#     table[!, "stop"] = last.(chromosome.genedata[!, "locus"])
#     return table
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the optimal number of clusters for k-means clustering by maximizing the silhouette score.

# Algorithm
- Starts by evaluating the first 5 k values
- Continues evaluation if optimal k is at the edge of evaluated range
- Refines search by evaluating midpoints between k values around the current optimum
- Iterates until convergence (optimal k remains stable)

# Arguments
- `distance_matrix::Matrix`: Square matrix of pairwise distances between points
- `ks_to_try::Vector{Int}`: Vector of k values to evaluate. Defaults to [1, 2, ...] up to matrix size

# Returns
Named tuple containing:
- `optimal_number_of_clusters::Int`: The k value giving highest silhouette score
- `ks_assessed::Vector{Int}`: All k values that were evaluated
- `within_cluster_sum_of_squares::Vector{Float64}`: WCSS for each k assessed
- `silhouette_scores::Vector{Float64}`: Silhouette scores for each k assessed
"""
function fit_optimal_number_of_clusters(distance_matrix, ks_to_try=[1, 2, Mycelia.ks(max=size(distance_matrix, 1))...])
    # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
    # N = size(distance_matrix, 1)
    # ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
    # Int(round(size(distance_matrix, 1)/2))
    # insert!(ks_to_try, 2, Int(round(size(distance_matrix, 1)))))
    # ks_to_try = vcat([2^i for i in 0:Int(floor(log2(size(distance_matrix, 1))))], size(distance_matrix, 1))
    # @info "ks = $(ks_to_try)"
    @show ks_to_try
    
    # can calculate this for k >= 1
    # within_cluster_sum_of_squares = Union{Float64, Missing}[]
    within_cluster_sum_of_squares = Float64[]
    # these are only valid for k >= 2 so set initial value to missing
    # between_cluster_sum_of_squares = [missing, zeros(length(ks_to_try)-1)...]
    # silhouette_scores = Union{Float64, Missing}[]
    silhouette_scores = Float64[]

    for k in ks_to_try[1:5]
        @show k
        @time this_clustering = Clustering.kmeans(distance_matrix, k)
        push!(within_cluster_sum_of_squares, wcss(this_clustering))
        if k == 1
            this_silhouette_score = 0
        else
            this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
            
        end
        push!(silhouette_scores, this_silhouette_score)
        @show this_silhouette_score
    end      
    optimal_silhouette, optimal_index = findmax(silhouette_scores)
    while (optimal_index == length(silhouette_scores)) && (optimal_index != length(ks_to_try))
        k = ks_to_try[optimal_index+1]
        @show k
        @time this_clustering = Clustering.kmeans(distance_matrix, k)
        push!(within_cluster_sum_of_squares, wcss(this_clustering))
        this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
        push!(silhouette_scores, this_silhouette_score)
        @show this_silhouette_score
        optimal_silhouette, optimal_index = findmax(silhouette_scores)
    end
    previous_optimal_number_of_clusters = 0
    optimal_number_of_clusters = ks_to_try[optimal_index]
    done = false
    while optimal_number_of_clusters != previous_optimal_number_of_clusters
        @show optimal_number_of_clusters
        if optimal_index == 1
            window_of_focus = ks_to_try[optimal_index:optimal_index+1]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        elseif optimal_index == length(ks_to_try)
            window_of_focus = ks_to_try[optimal_index-1:optimal_index]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        else
            window_of_focus = ks_to_try[optimal_index-1:optimal_index+1]
        end
        # @show window_of_focus
        midpoints = [
            Int(round(Statistics.mean(window_of_focus[1:2]))),
            Int(round(Statistics.mean(window_of_focus[2:3])))
            ]
        # @show sort(vcat(midpoints, window_of_focus))
        @show midpoints
        
        for k in midpoints
            insertion_index = first(searchsorted(ks_to_try, k))
            if ks_to_try[insertion_index] != k
                insert!(ks_to_try, insertion_index, k)
                @show k
                @time this_clustering = Clustering.kmeans(distance_matrix, k)
                insert!(within_cluster_sum_of_squares, insertion_index, wcss(this_clustering))
                this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
                insert!(silhouette_scores, insertion_index, this_silhouette_score)
                @show this_silhouette_score
            end
        end
        
        previous_optimal_number_of_clusters = optimal_number_of_clusters
        optimal_silhouette, optimal_index = findmax(silhouette_scores)
        optimal_number_of_clusters = ks_to_try[optimal_index]
    end
    @assert length(within_cluster_sum_of_squares) == length(silhouette_scores)
    ks_assessed = ks_to_try[1:length(within_cluster_sum_of_squares)]
    return (;optimal_number_of_clusters, ks_assessed, within_cluster_sum_of_squares, silhouette_scores)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Visualizes cluster assessment metrics and saves the resulting plots.

# Arguments
- `clustering_results`: A named tuple containing:
    * `ks_assessed`: Vector of k values tested
    * `within_cluster_sum_of_squares`: Vector of WCSS scores
    * `silhouette_scores`: Vector of silhouette scores
    * `optimal_number_of_clusters`: Integer indicating optimal k

# Details
Creates two plots:
1. Within-cluster sum of squares (WCSS) vs number of clusters
2. Silhouette scores vs number of clusters

Both plots include a vertical line indicating the optimal number of clusters.

# Outputs
Saves two SVG files in the project directory:
- `wcss.svg`: WCSS plot
- `silhouette.svg`: Silhouette scores plot
"""
function plot_optimal_cluster_assessment_results(clustering_results)
    p1 = StatsPlots.plot(
        ks_assessed[1:length(within_cluster_sum_of_squares)],
        within_cluster_sum_of_squares,
        ylabel = "within cluster sum of squares\n(lower is better)",
        xlabel = "n clusters",
        legend=false
    )
    StatsPlots.vline!(p1, [optimal_number_of_clusters])
    p2 = StatsPlots.plot(
        ks_assessed[1:length(silhouette_scores)],
        silhouette_scores,
        ylabel = "silhouette scores\n(higher is better)",
        xlabel = "n clusters",
        title = "Optimal n clusters = $(optimal_number_of_clusters)",
        legend=false
    )
    StatsPlots.vline!(p2, [optimal_number_of_clusters])
    # TODO write me out
    display(p2)
    StatsPlots.savefig(p1, "$DIR/wcss.svg")
    display(p1)
    StatsPlots.savefig(p2, "$DIR/silhouette.svg")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the document frequency of tokens across a collection of documents.

# Arguments
- `documents`: Collection of text documents where each document is a string

# Returns
- Dictionary mapping each unique token to the number of documents it appears in

# Description
Computes how many documents contain each unique token. Each document is tokenized 
by splitting on whitespace. Tokens are counted only once per document, regardless 
of how many times they appear within that document.
"""
function document_frequency(documents)
    document_tokens = Set(split(strip(first(documents))))
    countmap = StatsBase.countmap(document_tokens)
    for document in documents[2:end]
        document_tokens = Set(split(strip(document)))
        this_countmap = StatsBase.countmap(document_tokens)
        merge!(+, countmap, this_countmap)
    end
    return countmap
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Within-Cluster Sum of Squares (WCSS) for a clustering result.

# Arguments
- `clustering_result`: A clustering result object containing:
  - `counts`: Vector with number of points in each cluster
  - `assignments`: Vector of cluster assignments for each point
  - `costs`: Vector of distances/costs from each point to its cluster center

# Returns
- `Float64`: The total within-cluster sum of squared distances

# Description
WCSS measures the compactness of clusters by summing the squared distances 
between each data point and its assigned cluster center.
"""
function wcss(clustering_result)
    n_clusters = length(clustering_result.counts)
    total_squared_cost = 0.0
    for cluster_id in 1:n_clusters
        cluster_indices = clustering_result.assignments .== cluster_id
        total_squared_cost += sum(clustering_result.costs[cluster_indices] .^ 2)
    end
    return total_squared_cost
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determine the optimal number of clusters using hierarchical clustering with iterative refinement.

# Arguments
- `distance_matrix::Matrix`: Square matrix of pairwise distances between observations
- `ks_to_try::Vector{Int}`: Vector of cluster counts to evaluate (default: 1, 2, and sequence from `Mycelia.ks()`)

# Returns
Named tuple containing:
- `optimal_number_of_clusters::Int`: Best number of clusters found
- `ks_assessed::Vector{Int}`: All cluster counts that were evaluated
- `silhouette_scores::Vector{Float64}`: Silhouette scores for each k assessed
- `hclust_result`: Hierarchical clustering result object

# Details
Uses Ward's linkage method and silhouette scores to evaluate cluster quality. 
Implements an adaptive search that focuses on promising regions between initially tested k values.
For k=1, silhouette score is set to 0 as a special case.
"""
function fit_optimal_number_of_clusters_hclust(distance_matrix, ks_to_try=[1, 2, Mycelia.ks(max=size(distance_matrix, 1))...])
    # ks_to_try = [1, Int(round(size(distance_matrix, 1)/2)), size(distance_matrix, 1)]
    # N = size(distance_matrix, 1)
    # ks_to_try = push!([2^i for i in 0:Int(floor(log2(N)))], N)
    
    # Int(round(size(distance_matrix, 1)/2))
    # insert!(ks_to_try, 2, Int(round(size(distance_matrix, 1)))))
    # ks_to_try = vcat([2^i for i in 0:Int(floor(log2(size(distance_matrix, 1))))], size(distance_matrix, 1))
    # @info "ks = $(ks_to_try)"
    @show ks_to_try
    
    # can calculate this for k >= 1
    # within_cluster_sum_of_squares = Union{Float64, Missing}[]
    # within_cluster_sum_of_squares = Float64[]
    # these are only valid for k >= 2 so set initial value to missing
    # between_cluster_sum_of_squares = [missing, zeros(length(ks_to_try)-1)...]
    # silhouette_scores = Union{Float64, Missing}[]
    silhouette_scores = Float64[]
    # wcss_scores = Float64[]



    @show "initial heirarchical clustering"
    @time hclust_result = Clustering.hclust(distance_matrix, linkage=:ward, branchorder=:optimal)
    # for k in ks_to_try[1:3]
    for k in ks_to_try
        @show k
        this_clustering = Clustering.cutree(hclust_result, k=k)
        # push!(within_cluster_sum_of_squares, wcss(this_clustering))
        if k == 1
            this_silhouette_score = 0
            # this_wcss = Inf
        else
            this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
            # this_wcss = wcss(this_clustering)
        end
        push!(silhouette_scores, this_silhouette_score)
        # push!(within_cluster_sum_of_squares, this_wcss)
        @show this_silhouette_score
    end      
    optimal_silhouette, optimal_index = findmax(silhouette_scores)
    # while (optimal_index == length(silhouette_scores)) && (optimal_index != length(ks_to_try))
    #     k = ks_to_try[optimal_index+1]
    #     @show k
    #     this_clustering = Clustering.cutree(hclust_result, k=k)
    #     # push!(within_cluster_sum_of_squares, wcss(this_clustering))
    #     this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
    #     push!(silhouette_scores, this_silhouette_score)
    #     @show this_silhouette_score
    #     optimal_silhouette, optimal_index = findmax(silhouette_scores)
    # end
    previous_optimal_number_of_clusters = 0
    optimal_number_of_clusters = ks_to_try[optimal_index]
    done = false
    while optimal_number_of_clusters != previous_optimal_number_of_clusters
        @show optimal_number_of_clusters
        if optimal_index == 1
            window_of_focus = ks_to_try[optimal_index:optimal_index+1]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        elseif optimal_index == length(ks_to_try)
            window_of_focus = ks_to_try[optimal_index-1:optimal_index]
            insert!(window_of_focus, 2, Int(round(Statistics.mean(window_of_focus))))
        else
            window_of_focus = ks_to_try[optimal_index-1:optimal_index+1]
        end
        # @show window_of_focus
        midpoints = [
            Int(round(Statistics.mean(window_of_focus[1:2]))),
            Int(round(Statistics.mean(window_of_focus[2:3])))
            ]
        # @show sort(vcat(midpoints, window_of_focus))
        @show midpoints
        
        for k in midpoints
            insertion_index = first(searchsorted(ks_to_try, k))
            if ks_to_try[insertion_index] != k
                insert!(ks_to_try, insertion_index, k)
                @show k
                this_clustering = Clustering.cutree(hclust_result, k=k)
                this_silhouette_score = Statistics.mean(Clustering.silhouettes(this_clustering, distance_matrix))
                insert!(silhouette_scores, insertion_index, this_silhouette_score)
                # this_wcss = wcss(this_clustering)
                # insert!(within_cluster_sum_of_squares, insertion_index, this_wcss)
                @show this_silhouette_score
            end
        end
        
        previous_optimal_number_of_clusters = optimal_number_of_clusters
        optimal_silhouette, optimal_index = findmax(silhouette_scores)
        optimal_number_of_clusters = ks_to_try[optimal_index]
    end
    # @assert length(within_cluster_sum_of_squares) == length(silhouette_scores)
    ks_assessed = ks_to_try[1:length(silhouette_scores)]
    return (;optimal_number_of_clusters, ks_assessed, silhouette_scores, hclust_result)
end



# taxon=694009
# $(CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome taxon 694009 --dehydrated --filename /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009.zip
# unzip /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009.zip -d /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009
# $(CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets rehydrate --directory /global/cfs/cdirs/m4269/cjprybol/coronaviridae/data/ncbi-datasets.694009



# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ncbi_datasets_genome(; kwargs...)

# Download and rehydrate a data package from [NCBI datasets genome tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)

# Specify the download using either

# - [taxon](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/)

# or

# - [accession(https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_accession/)

# # Arguments
# - `annotated::Bool=false`: Limit to annotated genomes.
# - `api_key::String=""`: Specify an NCBI API key.
# - `assembly_level::Array{String}=[]`: Limit to genomes at specific assembly levels (e.g., "chromosome", "complete", "contig", "scaffold"). Default is empty (no specific level).
# - `assembly_source::String="all"`: Limit to 'RefSeq' (GCF_) or 'GenBank' (GCA_) genomes. Default is "all".
# - `assembly_version::String=""`: Limit to 'latest' assembly accession version or include 'all' (latest + previous versions).
# - `chromosomes::Array{String}=[]`: Limit to a specified, comma-delimited list of chromosomes, or 'all' for all chromosomes.
# - `debug::Bool=false`: Emit debugging info.
# - `dehydrated::Bool=false`: Download a dehydrated zip archive including the data report and locations of data files (use the rehydrate command to retrieve data files).
# - `exclude_atypical::Bool=false`: Exclude atypical assemblies.
# - `filename::String=""`: Specify a custom file name for the downloaded data package. Default is "taxon.zip" or "accession.zip" if left blank.
# - `include::Array{String}=["genome"]`: Specify the data files to include (e.g., "genome", "rna", "protein"). Default includes genomic sequence files only.
# - `mag::String="all"`: Limit to metagenome assembled genomes (only) or remove them from the results (exclude). Default is "all".
# - `no_progressbar::Bool=false`: Hide the progress bar.
# - `preview::Bool=false`: Show information about the requested data package without downloading.
# - `reference::Bool=false`: Limit to reference genomes.
# - `released_after::String=""`: Limit to genomes released on or after a specified date (MM/DD/YYYY).
# - `released_before::String=""`: Limit to genomes released on or before a specified date (MM/DD/YYYY).
# - `search::Array{String}=[]`: Limit results to genomes with specified text in the searchable fields (e.g., species, assembly name).

# # Returns
# - The result of the API call.
# """

# function ncbi_datasets_genome(;
#         taxon=missing,
#         accession=missing,
#         annotated=false,
#         api_key="",
#         assembly_level=[],
#         assembly_source="all",
#         assembly_version="",
#         chromosomes=[],
#         debug=false,
#         dehydrated=false,
#         exclude_atypical=false,
#         filename="",
#         outdir="",
#         include=["genome"],
#         mag="all",
#         no_progressbar=false,
#         preview=false,
#         reference=false,
#         released_after="",
#         released_before="",
#         search=[])

#     # Base command
#     command = "$(CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome "
    
#     if !ismissing(taxon) && !ismissing(accession)
#         @error "can only provide taxon or accession, not both" taxon accession
#     elseif ismissing(taxon) && ismissing(accession)
#         @error "must provide either taxon or accession"
#     elseif !ismissing(taxon) && ismissing(accession)
#         command *= "taxon $(taxon) "
#         if isempty(filename)
#             filename = string(taxon) * ".zip"
#         end
#     elseif !ismissing(accession) && ismissing(taxon)
#         command *= "accession $(accession) "
#         if isempty(filename)
#             filename = string(accession) * ".zip"
#         end
#     end
    
#     @assert occursin(r"\.zip$", filename)
    
#     if !isempty(outdir)
#         filename = joinpath(outdir, filename)
#     end

#     annotated && (command *= "--annotated ")
#     !isempty(api_key) && (command *= "--api-key $api_key ")
#     !isempty(assembly_level) && (command *= "--assembly-level $(join(assembly_level, ',')) ")
#     command *= "--assembly-source $assembly_source "
#     !isempty(assembly_version) && (command *= "--assembly-version $assembly_version ")
#     !isempty(chromosomes) && (command *= "--chromosomes $(join(chromosomes, ',')) ")
#     debug && (command *= "--debug ")
#     dehydrated && (command *= "--dehydrated ")
#     exclude_atypical && (command *= "--exclude-atypical ")
#     command *= "--filename $filename "
#     !isempty(include) && (command *= "--include $(join(include, ',')) ")
#     command *= "--mag $mag "
#     no_progressbar && (command *= "--no-progressbar ")
#     preview && (command *= "--preview ")
#     reference && (command *= "--reference ")
#     !isempty(released_after) && (command *= "--released-after $released_after ")
#     !isempty(released_before) && (command *= "--released-before $released_before ")
#     for s in search
#         command *= "--search $s "
#     end

#     # Execute the command
#     println("Executing command: $command")
#     run(`$command`)
#     if dehydrated
#         @info "add code to rehydrate here"
#     end
#     return true
# end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Downloads and sets up MMseqs2 reference databases for sequence searching and analysis.

# Arguments
- `db::String`: Name of database to download (see table below)
- `dbdir::String`: Directory to store the downloaded database (default: "~/workspace/mmseqs")
- `force::Bool`: If true, force re-download even if database exists (default: false)
- `wait::Bool`: If true, wait for download to complete (default: true)

# Returns 
- Path to the downloaded database as a String

# Available Databases

| Database           | Type       | Taxonomy | Description                               |
|-------------------|------------|----------|-------------------------------------------|
| UniRef100         | Aminoacid  | Yes      | UniProt Reference Clusters - 100% identity|
| UniRef90          | Aminoacid  | Yes      | UniProt Reference Clusters - 90% identity |
| UniRef50          | Aminoacid  | Yes      | UniProt Reference Clusters - 50% identity |
| UniProtKB         | Aminoacid  | Yes      | Universal Protein Knowledge Base          |
| NR               | Aminoacid  | Yes      | NCBI Non-redundant proteins              |
| NT               | Nucleotide | No       | NCBI Nucleotide collection               |
| GTDB             | Aminoacid  | Yes      | Genome Taxonomy Database                  |
| PDB              | Aminoacid  | No       | Protein Data Bank structures             |
| Pfam-A.full      | Profile    | No       | Protein family alignments                |
| SILVA            | Nucleotide | Yes      | Ribosomal RNA database                   |

```
  Name                  Type            Taxonomy        Url                                                           
- UniRef100             Aminoacid            yes        https://www.uniprot.org/help/uniref
- UniRef90              Aminoacid            yes        https://www.uniprot.org/help/uniref
- UniRef50              Aminoacid            yes        https://www.uniprot.org/help/uniref
- UniProtKB             Aminoacid            yes        https://www.uniprot.org/help/uniprotkb
- UniProtKB/TrEMBL      Aminoacid            yes        https://www.uniprot.org/help/uniprotkb
- UniProtKB/Swiss-Prot  Aminoacid            yes        https://uniprot.org
- NR                    Aminoacid            yes        https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- NT                    Nucleotide             -        https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- GTDB                  Aminoacid            yes        https://gtdb.ecogenomic.org
- PDB                   Aminoacid              -        https://www.rcsb.org
- PDB70                 Profile                -        https://github.com/soedinglab/hh-suite
- Pfam-A.full           Profile                -        https://pfam.xfam.org
- Pfam-A.seed           Profile                -        https://pfam.xfam.org
- Pfam-B                Profile                -        https://xfam.wordpress.com/2020/06/30/a-new-pfam-b-is-released
- CDD                   Profile                -        https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
- eggNOG                Profile                -        http://eggnog5.embl.de
- VOGDB                 Profile                -        https://vogdb.org
- dbCAN2                Profile                -        http://bcb.unl.edu/dbCAN2
- SILVA                 Nucleotide           yes        https://www.arb-silva.de
- Resfinder             Nucleotide             -        https://cge.cbs.dtu.dk/services/ResFinder
- Kalamari              Nucleotide           yes        https://github.com/lskatz/Kalamari
```
"""
function download_mmseqs_db(;db, dbdir="$(homedir())/workspace/mmseqs", force=false, wait=true)
    Mycelia.add_bioconda_env("mmseqs2")
    mkpath(dbdir)
    # sanitized_db = replace(db, "/" => "_")
    db_path = joinpath(dbdir, db)
    mkpath(dirname(db_path))
    if !isfile(db_path) || force
        cmd = `$(CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs databases --compressed 1 $(db) --remove-tmp-files 1 $(dbdir)/$(db) $(dbdir)/tmp`
        @time run(cmd, wait=wait)
    else
        @info "db $db @ $(db_path) already exists, set force=true to overwrite"
    end
    return db_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists available BLAST databases from the specified source.

# Arguments
- `source::String="ncbi"`: Source of BLAST databases, defaults to "ncbi"

# Returns
- `Vector{String}`: Array of available BLAST database names

# Notes
- Requires sudo/root privileges to install BLAST+ tools if not already present
- Internet connection required to fetch database listings
"""
function list_blastdbs(;source="ncbi")
    # Mycelia.add_bioconda_env("blast")
    # readlines(`$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl --showall`)
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    readlines(`update_blastdb --source $(source) --showall`)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Display available BLAST databases from specified source.

# Arguments
- `source::String="ncbi"`: Database source (default: "ncbi")

# Returns
- `DataFrame`: Table of available BLAST databases with columns:
  - NAME: Database name
  - SIZE (GB): Database size
  - LAST_UPDATED: Update date
  - and other metadata fields
"""
function showall_blastdbs(;source="ncbi")
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    blast_table_header = filter(!isempty, split(readlines(`update_blastdb --source $(source) --showall pretty`)[2], "  "))
    data, header = uCSV.read(IOBuffer(join(readlines(`update_blastdb --source $(source) --showall tsv`)[2:end], "\n")), delim="\t")
    df = sort(DataFrames.DataFrame(data, blast_table_header), "SIZE (GB)", rev=true)
    df[!, "LAST_UPDATED"] = map(dt -> Dates.Date(Dates.DateTime(first(split(dt, '.')), "yyyy-mm-ddTHH:MM:SS")), df[!, "LAST_UPDATED"])
    return df
end

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
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
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

# function remove_non_ascii(input::String)
#     return String(filter(x -> x <= '\x7f', input))
# end

# function sanitize_string(input::String)
#     return String(filter(x -> isvalid(Char, x), input))
# end

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

# Dependencies
Requires NCBI BLAST+ and perl-doc packages. May need sudo privileges for installation.
"""
function blastdb2table(;blastdb, outfile="", force=false)
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
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
    io = `blastdbcmd -entry 'all' -db $(blastdb) -outfmt $(outfmt_string)`
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
    # Mycelia.add_bioconda_env("blast")
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    @assert source in ["", "aws", "gcp", "ncbi"]
    mkpath(dbdir)
    current_directory = pwd()
    cd(dbdir)
    if isempty(source)
        @info "source not provided, letting blast auto-detect fastest download option"
        # cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress`
        cmd = `update_blastdb --decompress $(db)`
    else
        @info "downloading from source $(source)"
        if source == "ncbi"
            # --timeout 360 --passive no 
            # cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source) --timeout 360 --passive no`
            cmd = `update_blastdb --timeout 360 --passive no --decompress --source $(source) $(db)`
        else
            # cmd = `$(CONDA_RUNNER) run --live-stream -n blast update_blastdb.pl $(db) --decompress --source $(source)`
            cmd = `update_blastdb --decompress --source $(source) $(db)`
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
    try
        run(`sudo apt-get install ncbi-blast+ perl-doc -y`)
    catch
        run(`apt-get install ncbi-blast+ perl-doc -y`)
    end
    # p = pipeline(`$(CONDA_RUNNER) run --live-stream -n blast blastdbcmd -db $(dbdir)/$(db) -entry all -outfmt %f`)
    p = pipeline(`blastdbcmd -db $(dbdir)/$(db) -entry all -outfmt %f`)
    if compressed
        p = pipeline(p, `gzip`)
    end
    run(pipeline(p, outfile))
    return outfile
end

# neo_import_dir = "/Users/cameronprybol/Library/Application Support/Neo4j Desktop/Application/relate-data/dbmss/dbms-8ab8baac-5dea-4137-bb24-e0b426447940/import"


# uploading over API is slow for remote and local connections
# Progress:   0%|▏                                        |  ETA: 8:13:39
# upload_nodes_over_api(graph, ADDRESS=local_neo4j_bolt_address, PASSWORD=local_neo4j_password)
# Progress:   0%|         

# # push to Neo4J Aura
# # run(`sudo touch /etc/neo4j/neo4j.conf`)
# run(`sudo neo4j stop`)
# # remote database needs to be running
# # needs to be big enough
# # leave off port from address
# # run(`neo4j-admin push-to-cloud --overwrite --verbose --bolt-uri=$(ADDRESS) --username=$(USERNAME) --password=$(PASSWORD)`)
# # run(`sudo neo4j-admin push-to-cloud --overwrite --verbose --dump-to "$(DIR)/test.db.dump" --bolt-uri=$(a) --username=$(USERNAME) --password=$(PASSWORD)`)
# run(`sudo neo4j-admin push-to-cloud --overwrite --verbose --bolt-uri=$(a) --username=$(USERNAME) --password=$(PASSWORD)`)


# cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.tax_id IS UNIQUE"
# @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.`scientific name` IS UNIQUE"
# @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.identifier IS UNIQUE"
# @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# parameters = ["$(n): row.$(n)" for n in filter(x -> x != "TYPE", names(node_table))]
# parameters = "{" * join(parameters, ", ") * "}"

# window_size = 10000
# V = DataFrames.nrow(node_table)
# windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
# ProgressMeter.@showprogress for (i, w) in enumerate(windows)
#     df_sub = node_table[w, :]
#     f = "node$i.tsv"
#     local_f_path = "$(temp_upload_dir)/$(f)"
#     uCSV.write(local_f_path, df_sub, delim='\t')
#     run(`chmod 777 $(local_f_path)`)
#     f_url = "file:///$(local_f_path)"
#     cmd =
#     """
#     LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
#     CREATE (node:$(NODE_TYPE) $(parameters))
#     """
#     cmd = rstrip(replace(cmd, '\n' => ' '))
#     cypher_cmd = Mycelia.cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd)
#     run(cypher_cmd) 
# end


# ProgressMeter.@showprogress for (i, w) in enumerate(windows)
#     df_sub = edge_table[w, :]
#     f = "edge$i.tsv"
#     local_f_path = "$(temp_upload_dir)/$(f)"
#     uCSV.write(local_f_path, df_sub, delim='\t')
#     run(`chmod 777 $(local_f_path)`)
#     f_url = "file:///$(local_f_path)"
#     cmd = 
#     """
#     LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
#     MATCH (src:$(src_type) {identifier: row.src})
#     MATCH (dst:$(dst_type) {identifier: row.dst})
#     MERGE (src)-[p:$(edge_type)]->(dst)
#     """
#     cmd = rstrip(replace(cmd, '\n' => ' '))
#     cypher_cmd = Mycelia.cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd)
#     run(cypher_cmd) 
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Upload all nodes from a MetaGraph to a Neo4j database, processing each unique node type separately.

# Arguments
- `graph`: A MetaGraph containing nodes to be uploaded
- `address`: Neo4j server address (e.g., "bolt://localhost:7687")
- `username`: Neo4j authentication username (default: "neo4j")
- `password`: Neo4j authentication password
- `format`: Data format for upload (default: "auto")
- `database`: Target Neo4j database name (default: "neo4j")
- `neo4j_import_directory`: Path to Neo4j's import directory for bulk loading
"""
function upload_nodes_to_neo4j(;graph, address, username="neo4j", password, format="auto", database="neo4j", neo4j_import_directory)
    
    node_types = unique(MetaGraphs.props(graph, v)[:TYPE] for v in Graphs.vertices(graph))
    # node_type_strings = Mycelia.type_to_string.(node_types)
    
    for node_type in node_types
        @info "uploading node_type => $(Mycelia.type_to_string(node_type))..."
        node_type_table = node_type_to_dataframe(node_type=node_type, graph=graph)
        try
            upload_node_table(table=node_type_table, address=address, password=password, neo4j_import_dir=neo4j_import_directory)
        catch e
            showerror(stdout, e)
        end
    end
    
    @info "done!"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert all nodes of a specific type in a MetaGraph to a DataFrame representation.

# Arguments
- `node_type`: The type of nodes to extract from the graph
- `graph`: A MetaGraph containing the nodes

# Returns
A DataFrame where:
- Each row represents a node of the specified type
- Columns correspond to all unique properties found across nodes
- Values are JSON-serialized strings for consistency

# Notes
- All values are normalized through JSON serialization
- Dictionary values receive double JSON encoding
- The TYPE column is converted using `type_to_string`
"""
function node_type_to_dataframe(;node_type, graph)
    node_type_indices = filter(v -> MetaGraphs.props(graph, v)[:TYPE] == node_type, Graphs.vertices(graph))
    node_type_parameters = unique(reduce(vcat, map(v -> collect(keys(MetaGraphs.props(graph, v))), node_type_indices)))
    node_type_table = DataFrames.DataFrame(Dict(p => [] for p in node_type_parameters))
    for node_index in node_type_indices      
        push!(node_type_table, MetaGraphs.props(graph, node_index))
    end
    # normalize
    node_type_table[!, "TYPE"] = Mycelia.type_to_string.(node_type_table[!, "TYPE"])
    for column in names(node_type_table)
        T = Union{unique(typeof.(node_type_table[!, column]))...}
        if T <: AbstractDict
            node_type_table[!, column] = map(d -> JSON.json(string(JSON.json(d))), node_type_table[!, column])
        else
            node_type_table[!, column] = JSON.json.(string.(node_type_table[!, column]))
        end
    end  
    return node_type_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Upload a DataFrame to Neo4j as nodes in batched windows.

# Arguments
- `table::DataFrame`: Input DataFrame where each row becomes a node. Must contain a "TYPE" column.
- `address::String`: Neo4j server address (e.g. "bolt://localhost:7687")
- `password::String`: Neo4j database password
- `neo4j_import_dir::String`: Directory path accessible to Neo4j for importing files
- `window_size::Int=1000`: Number of rows to process in each batch
- `username::String="neo4j"`: Neo4j database username
- `database::String="neo4j"`: Target Neo4j database name

# Notes
- All rows must have the same node type (specified in "TYPE" column)
- Column names become node properties
- Requires write permissions on neo4j_import_dir
- Large tables are processed in batches of size window_size
"""
function upload_node_table(;table, window_size=1000, address, password, username="neo4j", database="neo4j", neo4j_import_dir)
    nrows = DataFrames.nrow(table)
    windows = (i:min(i+window_size-1,nrows) for i in 1:window_size:nrows)
    
    node_types = unique(table[!, "TYPE"])
    @assert length(node_types) == 1
    NODE_TYPE = Mycelia.type_to_string(first(node_types))
    parameters = ["$(n): row.$(n)" for n in filter(x -> !(x in ["TYPE"]), names(table))]
    parameters = "{" * join(parameters, ", ") * "}"

    ProgressMeter.@showprogress for (i, window) in enumerate(windows)
        df_sub = table[window, :]
        f = "node$i.tsv"
        local_f_path = "$(neo4j_import_dir)/$(f)"
        uCSV.write(local_f_path, df_sub, delim='\t')
        run(`chmod 777 $(local_f_path)`)
        f_url = "file:///$(f)"
        cmd =
        """
        LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
        CREATE (:`$(NODE_TYPE)` $(parameters))
        """
        cmd = rstrip(replace(cmd, '\n' => ' '))
        cypher_cmd = Mycelia.cypher(cmd, address = address, username = username, password = password, database = database)
        run(cypher_cmd) 
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Upload a single node from a MetaGraph to a Neo4j database using the HTTP API.

# Arguments
- `graph`: MetaGraph containing the node to be uploaded
- `v`: Vertex identifier in the graph
- `ADDRESS`: Neo4j server address (e.g. "http://localhost:7474")
- `USERNAME`: Neo4j authentication username (default: "neo4j")
- `PASSWORD`: Neo4j authentication password
- `DATABASE`: Target Neo4j database name (default: "neo4j")

# Details
Generates and executes a Cypher MERGE command using the node's properties. The node's :TYPE 
and :identifier properties are used for node labeling, while other non-empty properties 
are added as node properties.
"""
function upload_node_over_api(graph, v; ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j")
    node_type = MetaGraphs.props(graph, v)[:TYPE]
    node_identifier = MetaGraphs.props(graph, v)[:identifier]
    node_parameters = filter(x -> 
            !(x[1] in (:TYPE, :identifier)) && 
            !(ismissing(x[2]) || isempty(x[2])), 
        MetaGraphs.props(graph, v))
    params_string = join(["$(string(key)): \"$(string(value))\"" for (key, value) in node_parameters], ", ")
    node_type_string = Mycelia.type_to_string(node_type)
    node_identifier_string = string(node_identifier)
    cmd = 
    """
    MERGE (`$(node_identifier_string)`:`$(node_type_string)` {$(params_string)})
    """
    cmd = strip(cmd)
    cypher_cmd = Mycelia.cypher(cmd, address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE)
    run(cypher_cmd)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Uploads all nodes from the given graph to a specified API endpoint.

# Arguments
- `graph`: The graph containing the nodes to be uploaded.
- `ADDRESS`: The API endpoint address.
- `USERNAME`: The username for authentication (default: "neo4j").
- `PASSWORD`: The password for authentication.
- `DATABASE`: The database name (default: "neo4j").
"""
function upload_nodes_over_api(graph; ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j")
    ProgressMeter.@showprogress for v in Graphs.vertices(graph)
        upload_node_over_api(graph, v, ADDRESS=ADDRESS, USERNAME=USERNAME, PASSWORD=PASSWORD, DATABASE=DATABASE)
    end    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates unique identifier constraints for each node type in a Neo4j database.

# Arguments
- `graph`: A MetaGraph containing nodes with TYPE properties
- `address`: Neo4j server address
- `username`: Neo4j username (default: "neo4j")
- `password`: Neo4j password
- `database`: Neo4j database name (default: "neo4j")

# Details
Extracts unique node types from the graph and creates Neo4j constraints ensuring
each node of a given type has a unique identifier property.

Failed constraint creation attempts are silently skipped.
"""
function create_node_constraints(graph; address, username="neo4j", password, database="neo4j")
    node_types = unique(MetaGraphs.props(graph, v)[:TYPE] for v in Graphs.vertices(graph))
    node_type_strings = map(t -> Mycelia.type_to_string(t), node_types)
    for t in node_type_strings
        cmd = "CREATE CONSTRAINT ON (t:`$(t)`) ASSERT t.identifier IS UNIQUE"
        try
            cypher = Mycelia.cypher(cmd, address = address, username = username, password = password, database = database)
            @show cypher
            @time run(cypher)
        catch
            continue
        end
    end
end

# function type_to_string(T::KMER_TYPE) where {KMER_TYPE <: Kmers.Kmer}
#     @show "here"
#     return 
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a type to its string representation, with special handling for Kmer types.

# Arguments
- `T`: The type to convert to string

# Returns
- String representation of the type
  - For Kmer types: Returns "Kmers.DNAKmer{K}" where K is the kmer length
  - For other types: Returns the standard string representation
"""
function type_to_string(T)
    if T <: Kmers.Kmer
        return "Kmers.DNAKmer{$(T.parameters[2])}"
    else
        return string(T)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Converts an AbstractString type to its string representation.

# Arguments
- `T::AbstractString`: The string type to convert

# Returns
A string representation of the input type
"""
function type_to_string(T::AbstractString)
    return string(T)
end

# function batch_upload_edge_type_over_url_from_graph(src_type, dst_type, edge_type, graph, ADDRESS, USERNAME, PASSWORD, DATABASE; window_size=100)    
#     src_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == src_type, Graphs.vertices(graph))
#     dst_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == dst_type, Graphs.vertices(graph))
#     edges_to_upload = []
#     for src_node in src_nodes
#         outneighbors = Graphs.outneighbors(graph, src_node)
#         outneighbors = filter(outneighbor -> outneighbor in dst_nodes, outneighbors)
#         for outneighbor in outneighbors
#             this_edge = Graphs.Edge(src_node, outneighbor)
#             @assert MetaGraphs.get_prop(graph, this_edge, :TYPE) == edge_type
#             push!(edges_to_upload, this_edge)
#         end
#     end
    
#     N = length(edges_to_upload)
#     windows = [i:min(i+window_size-1,N) for i in 1:window_size:N]
    
#     ProgressMeter.@showprogress for window in windows
#         cmds = []
#         for (i, e) in enumerate(edges_to_upload[window])
#             edge_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, e)))
#             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, e, param))'" for param in edge_params]
#             joint_params = join(params, ", ")
#             node_cmds = 
#             """
#             MERGE (src$(i):$(MetaGraphs.props(graph, e.src)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.src)[:identifier])'})
#             MERGE (dst$(i):$(MetaGraphs.props(graph, e.dst)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.dst)[:identifier])'})
#             """
#             if !isempty(joint_params)
#                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE]) {$(joint_params)}]->(dst$(i))"
#             else
#                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE])]->(dst$(i))"
#             end
#             cmd = node_cmds * relationship_cmd
#             cmd = replace(cmd, '\n' => ' ')
#             push!(cmds, cmd)
#         end
#         cmd = join(cmds, ' ')
#         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
#         run(cypher_cmd)
#     end    
# end

# function batch_upload_node_type_over_url_from_graph(node_type, graph, ADDRESS, USERNAME, PASSWORD, DATABASE, window_size=100)
#     node_type_params = Set{Symbol}()
#     vertices_of_type = [v for v in Graphs.vertices(graph) if (graph.vprops[v][:TYPE] == node_type)]
    
#     V = length(vertices_of_type)
#     windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
    
#     ProgressMeter.@showprogress for window in windows
#         cmds = []
#         for (i, v) in enumerate(vertices_of_type[window])
#             node_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, v)))
#             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, v, param))'" for param in node_params]
#             # params = ["$(string(param)):'$(escape_string(MetaGraphs.get_prop(graph, v, param)))'" for param in node_params]
#             joint_params = join(params, ", ")
#             cmd = "MERGE (node$(i):$(node_type) {$(joint_params)})"
#             push!(cmds, cmd)
#         end
#         cmd = join(cmds, ' ')
#         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
#         run(cypher_cmd)
#     end    
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Upload nodes of a specific type from a graph to a Neo4j database using MERGE operations.

# Arguments
- `node_type`: The type label for the nodes to upload
- `graph`: Source MetaGraph containing the nodes
- `ADDRESS`: Neo4j server address (e.g. "bolt://localhost:7687")
- `PASSWORD`: Neo4j database password
- `USERNAME="neo4j"`: Neo4j username (defaults to "neo4j")
- `DATABASE="neo4j"`: Target Neo4j database name (defaults to "neo4j")
- `window_size=100`: Number of nodes to upload in each batch (defaults to 100)

# Details
Performs batched uploads of nodes using Neo4j MERGE operations. Node properties are 
automatically extracted from the graph vertex properties, excluding the 'TYPE' property.
"""
function upload_node_type_over_url_from_graph(;node_type, graph, ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j", window_size=100)
    node_type_params = Set{Symbol}()
    vertices_of_type = [v for v in Graphs.vertices(graph) if (graph.vprops[v][:TYPE] == node_type)]
    
    V = length(vertices_of_type)
    windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
    
    ProgressMeter.@showprogress for window in windows
        cmds = []
        for (i, v) in enumerate(vertices_of_type[window])
            node_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, v)))
            params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, v, param))'" for param in node_params]
            # params = ["$(string(param)):'$(escape_string(MetaGraphs.get_prop(graph, v, param)))'" for param in node_params]
            joint_params = join(params, ", ")
            cmd = "MERGE (node$(i):$(node_type) {$(joint_params)})"
            push!(cmds, cmd)
        end
        cmd = join(cmds, ' ')
        cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
        run(cypher_cmd)
    end    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Upload edges of a specific type from a MetaGraph to a Neo4j database, batching uploads in windows.

# Arguments
- `src_type`: Type of source nodes to filter
- `dst_type`: Type of destination nodes to filter  
- `edge_type`: Type of edges to upload
- `graph`: MetaGraph containing the nodes and edges
- `ADDRESS`: Neo4j server URL
- `USERNAME`: Neo4j username (default: "neo4j")
- `PASSWORD`: Neo4j password
- `DATABASE`: Neo4j database name (default: "neo4j")
- `window_size`: Number of edges to upload in each batch (default: 100)

# Details
- Filters edges based on source, destination and edge types
- Preserves all edge properties except :TYPE when uploading
- Uses MERGE operations to avoid duplicate nodes/relationships
- Uploads are performed in batches for better performance
- Progress is shown via ProgressMeter

# Returns
Nothing
"""
function upload_edge_type_over_url_from_graph(;src_type, dst_type, edge_type, graph, ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j", window_size=100)    
    src_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == src_type, Graphs.vertices(graph))
    dst_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == dst_type, Graphs.vertices(graph))
    edges_to_upload = []
    for src_node in src_nodes
        outneighbors = Graphs.outneighbors(graph, src_node)
        outneighbors = filter(outneighbor -> outneighbor in dst_nodes, outneighbors)
        for outneighbor in outneighbors
            this_edge = Graphs.Edge(src_node, outneighbor)
            @assert MetaGraphs.get_prop(graph, this_edge, :TYPE) == edge_type
            push!(edges_to_upload, this_edge)
        end
    end
    
    N = length(edges_to_upload)
    windows = (i:min(i+window_size-1,N) for i in 1:window_size:N)
    
    ProgressMeter.@showprogress for window in windows
        cmds = []
        for (i, e) in enumerate(edges_to_upload[window])
            edge_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, e)))
            params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, e, param))'" for param in edge_params]
            joint_params = join(params, ", ")
            node_cmds = 
            """
            MERGE (src$(i):$(MetaGraphs.props(graph, e.src)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.src)[:identifier])'})
            MERGE (dst$(i):$(MetaGraphs.props(graph, e.dst)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.dst)[:identifier])'})
            """
            if !isempty(joint_params)
                relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE]) {$(joint_params)}]->(dst$(i))"
            else
                relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE])]->(dst$(i))"
            end
            cmd = node_cmds * relationship_cmd
            cmd = replace(cmd, '\n' => ' ')
            push!(cmds, cmd)
        end
        cmd = join(cmds, ' ')
        cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
        run(cypher_cmd)
    end    
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a command to execute Neo4j Cypher queries via cypher-shell.

# Arguments
- `cmd`: The Cypher query command to execute
- `address::String="neo4j://localhost:7687"`: Neo4j server address
- `username::String="neo4j"`: Neo4j authentication username
- `password::String="password"`: Neo4j authentication password 
- `format::String="auto"`: Output format (auto, verbose, or plain)
- `database::String="neo4j"`: Target Neo4j database name

# Returns
- `Cmd`: A command object ready for execution
"""
function cypher(cmd;
    address="neo4j://localhost:7687",
    username="neo4j",
    password="password",
    format="auto",
    database="neo4j"
    )
    # cmd = `cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
    cmd = `/home/jovyan/.local/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
#    cmd = `/home/jupyter-cjprybol/software/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
    return cmd
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists all available Neo4j databases on the specified server.

# Arguments
- `address::String`: Neo4j server address (e.g. "neo4j://localhost:7687")
- `username::String="neo4j"`: Neo4j authentication username
- `password::String`: Neo4j authentication password

# Returns
- `DataFrame`: Contains database information with columns typically including:
  - name: Database name
  - address: Database address
  - role: Database role (e.g., primary, secondary)
  - status: Current status (e.g., online, offline)
  - default: Boolean indicating if it's the default database
"""
function list_databases(;address, username="neo4j", password)
    cmd = "show databases"
    database = "system"
    cmd = cypher(cmd, address=address, username=username, password=password, database=database)
    return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, quotes='"', encodings=Dict("FALSE" => false, "TRUE" => true))...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a new Neo4j database instance if it doesn't already exist.

# Arguments
- `database::String`: Name of the database to create
- `address::String`: Neo4j server address (e.g. "neo4j://localhost:7687")
- `username::String`: Neo4j authentication username (defaults to "neo4j")
- `password::String`: Neo4j authentication password

# Notes
- Requires system database privileges to execute
- Silently returns if database already exists
- Temporarily switches to system database to perform creation
"""
function create_database(;database, address, username="neo4j", password)
    current_databases = list_databases(;address, username, password)
    if database in current_databases[!, "name"]
        return
    else
        f = run
        cmd = "create database $(database)"
        # switch database to system, so that we can create the user-specific database in the system
        database = "system"
        run(cypher(;address, username, password, database, cmd, f))
    end
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function cypher(;address, username, password, database, cmd)
#     return `cypher-shell --address $address --username $username --password $password --database $(database) --format auto $(cmd)`
# end

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

Query Neo4j database to find all descendant taxonomic IDs for a given taxonomic ID.

# Arguments
- `tax_id`: Source taxonomic ID to find children for
- `DATABASE_ID`: Neo4j database identifier (required)
- `USERNAME`: Neo4j database username (default: "neo4j")
- `PASSWORD`: Neo4j database password (required)

# Returns
`Vector{Int}`: Sorted array of unique child taxonomic IDs
"""
function taxonomic_id_to_children(tax_id; DATABASE_ID, USERNAME="neo4j", PASSWORD)
    DATABASE = "neo4j"
    ADDRESS="neo4j+s://$(database_id).databases.neo4j.io:7687"
    
    # NOTE! *, or 0 distance (e.g. [*0..2]) step range will include source node!!!!
    cmd = "MATCH (n)<-[*]-(n2) WHERE n.tax_id IS NOT NULL AND n.tax_id = \"$(tax_id)\" RETURN DISTINCT n2.tax_id AS tax_id"
    println(cmd)
    
    cypher = cypher(cmd, address=ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE)
    tax_ids = readlines(open(cypher))[2:end]
    tax_ids = strip.(tax_ids, '"')
    tax_ids = parse.(Int, tax_ids)
    return unique(tax_ids)
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

Generate and display frequency counts for all columns in a DataFrame.

# Arguments
- `table::DataFrame`: Input DataFrame to analyze

# Details
Iterates through each column in the DataFrame and displays:
1. The column name
2. A Dict mapping unique values to their frequencies using StatsBase.countmap
"""
function countmap_columns(table)
    for n in names(refseq_metadata)
        display(n)
        display(StatsBase.countmap(refseq_metadata[!, n]))
    end
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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function parse_xam(xam; filter_unmapped=false, primary_only=false, min_mapping_quality=0, min_align_length=1)
#     if occursin(r"\.bam$", xam)
#         MODULE = XAM.BAM
#         io = open(xam)
#     elseif occursin(r"\.sam$", xam)
#         MODULE = XAM.SAM
#         io = open(xam)
#     elseif occursin(r"\.sam.gz$", xam)
#         MODULE = XAM.SAM
#         io = CodecZlib.GzipDecompressorStream(open(xam))
#     else
#         error("unrecognized file extension in file: $xam")
#     end
#     # reader = open(MODULE.Reader, io)
#     reader = MODULE.Reader(io)
#     header = reader.header
#     record_iterator = Iterators.filter(record -> true, reader)
#     if filter_unmapped
#         record_iterator = Iterators.filter(record -> MODULE.ismapped(record), record_iterator)
#     end
#     if primary_only
#         record_iterator = Iterators.filter(record -> MODULE.isprimary(record), record_iterator)
#     end
#     record_iterator = Iterators.filter(record -> MODULE.mappingquality(record) >= min_mapping_quality, record_iterator)
#     record_iterator = Iterators.filter(record -> MODULE.alignlength(record) >= min_align_length, record_iterator)
#     records = sort(collect(record_iterator), by=x->[MODULE.refname(x), MODULE.position(x)])
#     # reset header to specify sorted
#     header.metainfo[1] = MODULE.MetaInfo("HD", ["VN" => 1.6, "SO" => "coordinate"])
#     close(io)
#     return (;records, header)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a SAM/BAM file into a summary DataFrame containing alignment metadata.

# Arguments
- `xam::AbstractString`: Path to input SAM (.sam), BAM (.bam), or gzipped SAM (.sam.gz) file

# Returns
DataFrame with columns:
- `template`: Read name
- `flag`: SAM flag
- `reference`: Reference sequence name
- `position`: Alignment position range (start:end)
- `mappingquality`: Mapping quality score
- `alignment_score`: Alignment score (AS tag)
- `isprimary`: Whether alignment is primary
- `alignlength`: Length of the alignment
- `ismapped`: Whether read is mapped
- `mismatches`: Number of mismatches (NM tag)

Note: Only mapped reads are included in the output DataFrame.
"""
function parse_xam_to_summary_table(xam)
    record_table = DataFrames.DataFrame(
        template = String[],
        flag = UInt16[],
        reference = String[],
        position = UnitRange{Int}[],
        mappingquality = UInt8[],
        alignment_score = Int[],
        isprimary = Bool[],
        # cigar = String[],
        # rnext = String[],
        # pnext = Int[],
        # tlen = Int[],
        # sequence = BioSequences.LongDNA{4}[],
        # quality = UInt8[],
        alignlength = Int[],
        ismapped = Bool[],
        # alignment = BioAlignments.Alignment[],
        mismatches = Int[]
    )
    if occursin(r"\.bam$", xam)
        MODULE = XAM.BAM
        io = open(xam)
    elseif occursin(r"\.sam$", xam)
        MODULE = XAM.SAM
        io = open(xam)
    elseif occursin(r"\.sam.gz$", xam)
        MODULE = XAM.SAM
        io = CodecZlib.GzipDecompressorStream(open(xam))
    else
        error("unrecognized file extension in file: $xam")
    end
    # reader = open(MODULE.Reader, io)
    reader = MODULE.Reader(io)
    header = reader.header
    for record in reader
        if XAM.SAM.ismapped(record)
            # @assert !ismissing()
            row = (
                template = XAM.SAM.tempname(record),
                flag = XAM.flag(record),
                reference = XAM.SAM.refname(record),
                position = XAM.SAM.position(record):XAM.SAM.rightposition(record),
                mappingquality = XAM.SAM.mappingquality(record),
                # cigar = XAM.SAM.cigar(record),
                # rnext = XAM.SAM.nextrefname(record),
                # pnext = XAM.SAM.nextposition(record),
                # tlen = XAM.SAM.templength(record),
                # sequence = XAM.SAM.sequence(record),
                # quality = XAM.SAM.quality(record),
                alignlength = XAM.SAM.alignlength(record),
                ismapped = XAM.SAM.ismapped(record),
                isprimary = XAM.SAM.isprimary(record),
                # alignment = XAM.SAM.alignment(record),
                alignment_score = record["AS"],
                mismatches = record["NM"]
                )
            push!(record_table, row, promote=true)
        end
    end
    # records = sort(collect(record_iterator), by=x->[MODULE.refname(x), MODULE.position(x)])
    # reset header to specify sorted
    # header.metainfo[1] = MODULE.MetaInfo("HD", ["VN" => 1.6, "SO" => "coordinate"])
    close(io)
    # return (;records, header)
    return record_table
end

# # map reads to the assembly and run qualimap QC
# bwt_index = "$(assembled_fasta).bwt"
# if !isfile(bwt_index)
#     run(`bwa index $(assembled_fasta)`)
# end

# mapped_reads_bam = "$(assembled_fasta).bwa.bam"
# if !isfile(mapped_reads_bam)
#     run(pipeline(
#         `bwa mem -R "@RG\tID:$(config["sample identifier"])\tSM:bar" -t $(Sys.CPU_THREADS) $(assembled_fasta) $(TRIMMED_FORWARD) $(TRIMMED_REVERSE)`,
#         `samtools collate -O - -`,
#         `samtools fixmate -m - -`,
#         `samtools sort`,
#         `samtools markdup - -`,
#         `samtools view -buh`,
#         mapped_reads_bam))
# end

# if !isfile("$(mapped_reads_bam).bai")
#     run(`samtools index $(mapped_reads_bam)`)
# end

# qualimap_report_pdf = "$(assembly_dir)/qualimap/report.pdf"
# qualimap_report_txt = "$(assembly_dir)/qualimap/genome_results.txt"

# if !isfile(qualimap_report_pdf) || !isfile(qualimap_report_txt)
#     run(`
#         qualimap bamqc
#         -nt $(Sys.CPU_THREADS)
#         -bam $(mapped_reads_bam)
#         -outdir $(assembly_dir)/qualimap
#         -outformat PDF:HTML
#         --output-genome-coverage $(mapped_reads_bam).genome_coverage.txt
#         `)
# end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run fastani with a query and reference list

Calculate Average Nucleotide Identity (ANI) between genome sequences using FastANI.

# Arguments
- `query_list::String`: Path to file containing list of query genome paths (one per line)
- `reference_list::String`: Path to file containing list of reference genome paths (one per line)
- `outfile::String`: Path to output file that will contain ANI results
- `threads::Int=Sys.CPU_THREADS`: Number of parallel threads to use
- `force::Bool=false`: If true, rerun analysis even if output file exists

# Output
Generates a tab-delimited file with columns:
- Query genome
- Reference genome  
- ANI value (%)
- Count of bidirectional fragment mappings
- Total query fragments

# Notes
- Requires FastANI to be available via Bioconda
- Automatically sets up required conda environment
"""
function fastani_list(;query_list="", reference_list="", outfile="", threads=Sys.CPU_THREADS, force=false)
    Mycelia.add_bioconda_env("fastani")
    if !isfile(outfile) || force
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n fastani fastANI --ql $(query_list) --rl $(reference_list) --threads $(threads) -o $(outfile)`)
        # run(
        # pipeline(
        #     `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastani fastANI --ql $(query_list) --rl $(reference_list) --threads $(threads) -o $(outfile)`,
        #     stdout=outfile * "fastani.stdout.txt",
        #     stderr=outfile * "fastani.stderr.txt"
        #     )
        # )
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Average Nucleotide Identity (ANI) between two genomes using FastANI.

# Arguments
- `query::String`: Path to query genome FASTA file
- `reference::String`: Path to reference genome FASTA file  
- `outfile::String`: Path to save FastANI results
- `force::Bool=false`: If true, overwrite existing output file

# Notes
- Requires FastANI to be available via Bioconda
- Stdout and stderr are captured in separate files with '.stdout.txt' and '.stderr.txt' suffixes
- ANI results are written to the specified outfile
"""
# ./fastANI -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
function fastani_pair(;query="", reference="", outfile="", force=false)
    Mycelia.add_bioconda_env("fastani")
    if !isfile(outfile) || force
        run(
        pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n fastani fastANI -q $(query) -r $(reference) -o $(outfile)`,
            stdout=outfile * "fastani.stdout.txt",
            stderr=outfile * "fastani.stderr.txt"
            )
        )
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the contig with the greatest number of total bases mapping to it

Identify the primary contig based on mapping coverage from Qualimap results.

# Arguments
- `qualimap_results::DataFrame`: DataFrame containing Qualimap alignment statistics with 
  columns "Contig" and "Mapped bases"

# Returns
- `String`: Name of the contig with the highest number of mapped bases

# Description
Takes Qualimap alignment results and determines which contig has the most total bases 
mapped to it, which often indicates the main chromosomal assembly.
"""
function determine_primary_contig(qualimap_results)
    primary_contig_index = last(findmax(qualimap_results[!, "Mapped bases"]))
    primary_contig = qualimap_results[primary_contig_index, "Contig"]
    return primary_contig
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Primary contig is defined as the contig with the most bases mapped to it

In the context of picking out phage from metagenomic assemblies
the longest contig is often bacteria whereas the highest coverage contigs are often primer-dimers or other PCR amplification artifacts.

Taking the contig that has the most bases mapped to it as a product of length * depth is cherry picked as our phage

Isolates and exports the primary contig from an assembly based on coverage depth × length.

The primary contig is defined as the contig with the highest total mapped bases 
(coverage depth × length). This method helps identify potential phage contigs in 
metagenomic assemblies, avoiding both long bacterial contigs and short high-coverage 
PCR artifacts.

# Arguments
- `assembled_fasta`: Path to the assembled contigs in FASTA format
- `assembled_gfa`: Path to the assembly graph in GFA format
- `qualimap_report_txt`: Path to Qualimap coverage report
- `identifier`: String identifier for the output file
- `k`: Integer representing k-mer size used in assembly
- `primary_contig_fasta`: Optional output filepath (default: "\${identifier}.primary_contig.fna")

# Returns
- Path to the output FASTA file containing the primary contig

# Notes
- For circular contigs, removes the k-mer closure scar if detected
- Trims k bases from the end if they match the first k bases
- Uses coverage × length to avoid both long bacterial contigs and short PCR artifacts
"""
function isolate_normalized_primary_contig(assembled_fasta, assembled_gfa, qualimap_report_txt, identifier, k::Int; primary_contig_fasta = "$(identifier).primary_contig.fna")
    
    qualimap_results = parse_qualimap_contig_coverage(qualimap_report_txt)
    primary_contig = determine_primary_contig(qualimap_results)

    # Find primary contig from scaffolds, then export as primary_contig.fasta
    for record in FASTX.FASTA.Reader(open(assembled_fasta))
        record_id = FASTX.identifier(record)
        if record_id == primary_contig
            primary_contig_sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)

            # If the primary contig is circular, need to trim to remove closure scar
            if Mycelia.contig_is_circular(assembled_gfa, primary_contig)
                # trim k-length from end before writing if it matches the first k of the contig
                if primary_contig_sequence[1:k] == primary_contig_sequence[end-k+1:end]
                    for i in 1:k pop!(primary_contig_sequence) end
                end
            end

        w = FASTX.FASTA.Writer(open(primary_contig_fasta, "w")) 
            write(w, FASTX.FASTA.Record(identifier, primary_contig_sequence))
            close(w)
        end
    end
    return primary_contig_fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns bool indicating whether the contig is cleanly assembled.

By cleanly assembled we mean that the contig does not have other contigs attached in the same connected component.

graph_file = path to assembly graph.gfa file
contig_name = name of the contig

Check if a contig exists in isolation within its connected component in an assembly graph.

# Arguments
- `graph_file::String`: Path to the assembly graph file in GFA format
- `contig_name::String`: Name/identifier of the contig to check

# Returns
- `Bool`: `true` if the contig exists alone in its connected component, `false` otherwise

# Details
A contig is considered "cleanly assembled" if it appears as a single entry in the 
assembly graph's connected components. This function parses the GFA file and checks
the contig's isolation status using the graph structure.
"""
function contig_is_cleanly_assembled(graph_file::String, contig_name::String)
    contig_table, records = Mycelia.gfa_to_structure_table(graph_file)
    matching_connected_components = findfirst(contig_table[!, "contigs"] .== contig_name)
    # contigs should be comma-seperated identifiers as strings, so if we get a hit then
    # it must be cleanly assembled.
    return !isnothing(matching_connected_components)
    # gfa_graph = Mycelia.parse_gfa(graph_file)
    # if !isempty(gfa_graph.gprops[:paths])
    #     # probably spades
    #     # gfa segment identifiers have _1 appended to fasta sequence identifiers for spades
    #     segment_identifier = contig_name * "_1"
    #     segment_node_string = gfa_graph.gprops[:paths][segment_identifier]["segments"]
    #     nodes_in_segment = replace.(split(segment_node_string, ','), r"[^\d]" => "")
    #     node_indices = [gfa_graph[n, :identifier] for n in nodes_in_segment]
    # else
    #     # megahit
    #     node_indices = [gfa_graph[contig_name, :identifier]]
    # end
    # connected_components = Graphs.connected_components(gfa_graph)
    # component_of_interest = first(filter(cc -> all(n -> n in cc, node_indices), connected_components))
    # if (1 <= length(component_of_interest) <= 2)
    #     return true
    # else
    #     return false
    # end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns bool indicating whether the contig is a circle

graph_file = path to assembly graph.gfa file
contig_name = name of the contig

Determine if a contig represents a circular structure in the assembly graph.

A circular contig is one where the sequence forms a complete loop in the assembly graph,
typically representing structures like plasmids, circular chromosomes, or other circular DNA elements.

# Arguments
- `graph_file::String`: Path to the assembly graph in GFA format
- `contig_name::String`: Name/identifier of the contig to check

# Returns
- `Bool`: `true` if the contig forms a circular structure, `false` otherwise
"""
function contig_is_circular(graph_file::String, contig_name::String)
    contig_table, records = Mycelia.gfa_to_structure_table(graph_file)
    # display(contig_name)
    # display(contig_table)
    matching_connected_components = findfirst(contig_table[!, "contigs"] .== contig_name)
    if !isnothing(matching_connected_components)
        return contig_table[matching_connected_components, "is_circular"] && contig_table[matching_connected_components, "is_closed"]
    else
        return false
    end
    # gfa_graph = Mycelia.parse_gfa(graph_file)
    # if !isempty(gfa_graph.gprops[:paths])
    #     # probably spades
    #     # gfa segment identifiers have _1 appended to fasta sequence identifiers for spades
    #     segment_identifier = contig_name * "_1"
    #     segment_node_string = gfa_graph.gprops[:paths][segment_identifier]["segments"]
    #     nodes_in_segment = replace.(split(segment_node_string, ','), r"[^\d]" => "")
    #     node_indices = [gfa_graph[n, :identifier] for n in nodes_in_segment]
    # else
    #     # megahit
    #     node_indices = [gfa_graph[contig_name, :identifier]]
    # end
    # connected_components = Graphs.connected_components(gfa_graph)
    # component_of_interest = first(filter(cc -> all(n -> n in cc, node_indices), connected_components))
    # subgraph, vertex_map = Graphs.induced_subgraph(gfa_graph, component_of_interest)
    # if Graphs.is_cyclic(subgraph)
    #     return true
    # else
    #     return false
    # end
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
function determine_percent_identity(reference_fasta, query_fasta)
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyzes k-mer frequency spectra from sequencing reads to detect sequencing errors and 
coverage patterns.

# Arguments
- `out_directory`: Directory where output files will be saved
- `forward_reads`: Path to forward reads FASTQ file
- `reverse_reads`: Path to reverse reads FASTQ file (optional)
- `k`: Length of k-mers to analyze (default: 17)
- `target_coverage`: Desired coverage level. If non-zero, returns downsampling rate (default: 0)
- `plot_size`: Dimensions of output plot in pixels (default: (600,400))

# Returns
- If `target_coverage > 0`: Returns downsampling rate as Float64
- If `target_coverage = 0`: Returns nothing, saves plots only

# Outputs
- `peak-detected.png`: K-mer spectra plot with detected boundaries
- `peak-detected.svg`: Vector format of the same plot
- `downsampling-rate.txt`: Created only if target_coverage is specified

# Algorithm
1. Counts canonical k-mers from input reads
2. Generates log-scale histogram of k-mer frequencies
3. Detects error threshold using regression or local minimum
4. Identifies coverage peak after error threshold
5. Calculates downsampling rate if target coverage specified
"""
# https://github.com/cjprybol/Mycelia/blob/e7fe50ffe2d18406fb70e0e24ebcfa45e0937596/notebooks/exploratory/2021-08-25-k-medoids-error-cluster-detection-multi-entity-graph-aligner-test.ipynb
function analyze_kmer_spectra(;out_directory, forward_reads="", reverse_reads="", k=17, target_coverage=0, plot_size=(600,400))
    @info "counting $k-mers"
    user_provided_reads = filter(x -> !isempty(x), [forward_reads, reverse_reads])
    canonical_kmer_counts = count_canonical_kmers(Kmers.Kmer{BioSequences.DNAAlphabet{4},k}, user_provided_reads)

    @info "determining max count"
    max_count = maximum(values(canonical_kmer_counts))
    @info "max count = $max_count"

    @info "generating histogram"
    kmer_counts_histogram = sort(collect(StatsBase.countmap(values(canonical_kmer_counts))), by=x->x[1])

    X = log2.(first.(kmer_counts_histogram))
    Y = log2.(last.(kmer_counts_histogram))
    
    @info "plotting kmer spectra"
    p = StatsPlots.scatter(
        X,
        Y,
        xlabel="log2(kmer_frequency)",
        ylabel="log2(# of kmers @ frequency)",
        label="",
        size=plot_size
    )

    earliest_y_min_index = last(findmin(Y))
    lower_boundary = X[earliest_y_min_index]
    lower_boundary_source = "first minimum"

    try
        # take the first 1/denominator datapoints in the set
        # to capture the error line on the left side of the graph
        @info "fitting error curve"
        denominators = [2^i for i in 1:5]
        coeficient_matrix = zeros(length(denominators), 2)
        for (i, denominator) in enumerate(denominators)
            prefix_index = Int(floor(length(X)/denominator))
            _x = X[1:prefix_index]
            _y = Y[1:prefix_index]
            model = GLM.lm(GLM.@formula(Y ~ X), DataFrames.DataFrame(X = _x, Y = _y))
            coeficient_matrix[i, :] = GLM.coef(model)
        end
        median_intercept = Statistics.median(coeficient_matrix[:, 1])
        median_slope = Statistics.median(coeficient_matrix[:, 2])

        X_intercept = (0 - median_intercept) / median_slope

        # some libraries detect the x_intercept being AFTER the end of the data
        # in these instances detect the earliest x-minimum
        if X_intercept < lower_boundary
            lower_boundary = X_intercept
            lower_boundary_source = "detected x-intercept"
        end
    catch
        @info "unable to fit regression"
    end

    p = StatsPlots.vline!(p,
        [lower_boundary],
        label="lower boundary ($(lower_boundary_source))"
    );
    
    is_above_lower_bounds = X .>= lower_boundary
    max_Y_post_error_intercept = first(findmax(Y[is_above_lower_bounds]))
    peak_indices = findall(is_above_lower_bounds .& (Y .== max_Y_post_error_intercept))
    peak_index = Int(round(Statistics.median(peak_indices)))

    p = StatsPlots.vline!([X[peak_index]], label="inferred sample coverage)")
    if isinteractive()
        display(p)
    end
    StatsPlots.savefig(p, "$out_directory/peak-detected.png")
    StatsPlots.savefig(p, "$out_directory/peak-detected.svg")
    
    if target_coverage != 0
        detected_coverage = 2^(X[peak_index])
        downsampling_rate = round(target_coverage/detected_coverage, sigdigits=3)
        downsampling_rate = min(downsampling_rate, 1)
        @info "downsampling rate = $downsampling_rate"

        outfile = "$out_directory/downsampling-rate.txt"
        open(outfile, "w") do io
            @info "writing downsampling rate to $outfile"
            println(io, downsampling_rate)
        end
        return downsampling_rate
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assess k-mer saturation in DNA sequences from FASTX files.

# Arguments
- `fastxs::AbstractVector{<:AbstractString}`: Vector of paths to FASTA/FASTQ files
- `kmer_type`: Type of k-mer to analyze (e.g., DNAKmer{21})
- `kmers_to_assess=Inf`: Maximum number of k-mers to process
- `power=10`: Base for exponential sampling intervals
- `min_count=1`: Minimum count threshold for considering a k-mer

# Returns
Named tuple containing:
- `sampling_points::Vector{Int}`: K-mer counts at which samples were taken
- `unique_kmer_counts::Vector{Int}`: Number of unique canonical k-mers at each sampling point
- `eof::Bool`: Whether the entire input was processed

# Details
Analyzes k-mer saturation by counting unique canonical k-mers at exponentially spaced 
intervals (powers of `power`). Useful for assessing sequence complexity and coverage.
Returns early if all possible k-mers are observed.
"""
function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}, kmer_type; kmers_to_assess=Inf, power=10, min_count = 1)
    # canonical_kmers = Set{kmer_type}()
    canonical_kmer_counts = Dict{kmer_type, Int}()
    
    @show kmer_type
    k = Kmers.ksize(Kmers.kmertype(kmer_type))
    
    max_possible_kmers = determine_max_canonical_kmers(k, DNA_ALPHABET)
    
    if kmers_to_assess == Inf
        # want to read the whole file and predict how long that will take
        # n_records = reduce(sum, map(f -> Mycelia.count_records(f), fastxs))
        kmers_to_assess = max_possible_kmers
        p = ProgressMeter.Progress(kmers_to_assess, 1)
    else
        p = ProgressMeter.Progress(kmers_to_assess, 1)
    end
    
    sampling_points = Int[0]
    i = 0
    while power^i <= kmers_to_assess
        push!(sampling_points, power^i)
        i += 1
    end
    
    unique_kmer_counts = zeros(Int, length(sampling_points))
    
    if length(sampling_points) < 3
        @info "increase the # of reads analyzed or decrease the power to acquire more data points"
        return (;sampling_points, unique_kmer_counts)
    end
    
    kmers_assessed = 0
    for fastx in fastxs
        for record in open_fastx(fastx)      
            record_sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
            for (index, kmer) in Kmers.EveryKmer{kmer_type}(record_sequence)
                canonical_kmer = BioSequences.canonical(kmer)
                if haskey(canonical_kmer_counts, canonical_kmer)
                    canonical_kmer_counts[canonical_kmer] += 1
                else
                    canonical_kmer_counts[canonical_kmer] = 1
                end
                kmers_assessed += 1
                if (length(canonical_kmer_counts) == max_possible_kmers)                 
                    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
                    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], length(canonical_kmer_counts))
                    return (;sampling_points, unique_kmer_counts, eof = false)
                elseif kmers_assessed in sampling_points
                    i = findfirst(sampling_points .== kmers_assessed)
                    unique_kmer_counts[i] = length(filter(x -> x[2] >= min_count, canonical_kmer_counts))
                    if i == length(sampling_points)
                        return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = false)
                    end
                end
                ProgressMeter.next!(p)
            end
        end
    end
    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], [length(canonical_kmer_counts)])    
    return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyze k-mer saturation in DNA sequences to determine optimal k value.

# Arguments
- `fastxs`: Vector of paths to FASTA/FASTQ files to analyze
- `power`: Base of logarithmic sampling points (default: 10)
- `outdir`: Optional output directory for plots and results
- `min_k`: Minimum k-mer size to test (default: 7)
- `max_k`: Maximum k-mer size to test (default: 17)
- `threshold`: Saturation threshold to determine optimal k (default: 0.1)
- `kmers_to_assess`: Maximum number of k-mers to sample (default: 10M)
- `plot`: Whether to generate saturation curves (default: true)

# Returns
Integer representing the first k value that achieves saturation below threshold.
If no k value meets the threshold, returns the k with minimum saturation.

# Details
- Tests only prime k values between min_k and max_k
- Generates saturation curves using logarithmic sampling
- Fits curves to estimate maximum unique k-mers
- If outdir is provided, saves plots as SVG and chosen k value to text file
"""
function assess_dnamer_saturation(fastxs::AbstractVector{<:AbstractString}; power=10, outdir::Union{Missing, String}=missing, min_k=7, max_k=17, threshold=0.1, kmers_to_assess=10_000_000, plot=true)
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        kmer_type = Kmers.kmertype(Kmers.Kmer{BioSequences.DNAAlphabet{4},k})
        sampling_points, kmer_counts, hit_eof = assess_dnamer_saturation(fastxs, kmer_type, kmers_to_assess=kmers_to_assess, power=power)
        @show sampling_points, kmer_counts, hit_eof
        observed_midpoint_index = findfirst(i -> kmer_counts[i] > last(kmer_counts)/2, 1:length(sampling_points))
        observed_midpoint = sampling_points[observed_midpoint_index]
        initial_parameters = Float64[maximum(kmer_counts), observed_midpoint]
        @time fit = LsqFit.curve_fit(calculate_v, sampling_points, kmer_counts, initial_parameters)
        max_canonical_kmers = determine_max_canonical_kmers(k, DNA_ALPHABET)
        if hit_eof
            inferred_maximum = last(kmer_counts)
        else
            inferred_maximum = max(Int(ceil(fit.param[1])), last(kmer_counts))
            if inferred_maximum > max_canonical_kmers
                inferred_maximum = max_canonical_kmers
            end
        end

        inferred_midpoint = Int(ceil(fit.param[2]))
        predicted_saturation = inferred_maximum / max_canonical_kmers
        @show k, predicted_saturation

        if plot
            scale = 300
            fontsize = 14
            p = StatsPlots.scatter(
                sampling_points,
                kmer_counts,
                label="observed counts",
                ylabel="# unique canonical kmers",
                xlabel="# kmers assessed",
                title = "sequencing saturation @ k = $k",
                # legend=:outertopright,
                # size=(2*scale, 1*scale),
                margins=5Plots.PlotMeasures.mm,
                titlefontsize=fontsize,
                xguidefontsize=fontsize,
                yguidefontsize=fontsize,
                legendfontsize=fontsize-2,
                xtickfontsize=fontsize-6,
                ytickfontsize=fontsize-4,
                # xrotation=45
                )
            StatsPlots.hline!(p, [max_canonical_kmers], label="absolute maximum", line = :solid, linewidth = 2)
            StatsPlots.hline!(p, [inferred_maximum], label="inferred maximum", line = :dash, linewidth = 2)
            StatsPlots.vline!(p, [inferred_midpoint], label="inferred midpoint", line = :dot, linewidth = 2)
            # xs = vcat(sampling_points, [last(sampling_points) * 2^i for i in 1:2])
            xs = sort([sampling_points..., inferred_midpoint])
            ys = calculate_v(xs, fit.param)
            StatsPlots.plot!(
                p,
                xs,
                ys,
                label="fit trendline",
                line=:dashdot,
                linewidth = 2)
            if k != first(ks)
                StatsPlots.plot!(p, legend=false)
            end
            display(p)
            if !ismissing(outdir)
                # StatsPlots.savefig(p, joinpath(outdir, "$k.png"))
                StatsPlots.savefig(p, joinpath(outdir, "$k.svg"))
            end
        end
            

        if predicted_saturation < minimum_saturation
            minimum_saturation = predicted_saturation
            min_k = k
            midpoint = inferred_midpoint 
        end
        if predicted_saturation < threshold
            if !ismissing(outdir)
                chosen_k_file = joinpath(outdir, "chosen_k.txt")
                println("chosen k = $k")
                open(chosen_k_file, "w") do io
                    println(io, k)
                end
            end
            return k
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyzes k-mer saturation in a FASTA/FASTQ file to determine optimal k-mer size.

# Arguments
- `fastx::AbstractString`: Path to input FASTA/FASTQ file
- `power::Int=10`: Exponent for downsampling k-mers (2^power)
- `outdir::String=""`: Output directory for results. Uses current directory if empty
- `min_k::Int=3`: Minimum k-mer size to evaluate
- `max_k::Int=17`: Maximum k-mer size to evaluate
- `threshold::Float64=0.1`: Saturation threshold for k-mer assessment
- `kmers_to_assess::Int=10_000_000`: Maximum number of k-mers to sample

# Returns
`Dict{Int,Float64}`: Dictionary mapping k-mer sizes to their saturation scores
"""
function assess_dnamer_saturation(fastx::AbstractString; power=10, outdir="", min_k=3, max_k=17, threshold=0.1, kmers_to_assess=10_000_000)
    assess_dnamer_saturation([fastx], outdir=outdir, min_k=min_k, max_k=max_k, threshold=threshold, power=power, kmers_to_assess=kmers_to_assess)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create an in-memory kmer-graph that records:
- all kmers
- counts
- all *observed* edges between kmers
- edge orientations
- edge counts

Construct a kmer-graph from one or more FASTX files (FASTA/FASTQ).

# Arguments
- `KMER_TYPE`: Type for kmer representation (e.g., `DNAKmer{K}`)
- `fastxs`: Vector of paths to FASTX files

# Returns
A `MetaGraph` where:
- Vertices represent unique kmers with properties:
  - `:kmer` => The kmer sequence
  - `:count` => Number of occurrences
- Edges represent observed kmer adjacencies with properties:
  - `:orientation` => Relative orientation of connected kmers
  - `:count` => Number of observed transitions
"""
function fastx_to_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString})
    
    @info "counting kmers"
    @time kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
    K = length(keys(kmer_counts))
    k = length(first(keys(kmer_counts)))
    
    @info "initializing graph with $(K) $(k)mer states"
    graph = MetaGraphs.MetaGraph(K)
    MetaGraphs.set_prop!(graph, :k, k)
    
    @info "adding kmers and kmer counts"
    ProgressMeter.@showprogress for (i, (kmer, count)) in enumerate(kmer_counts)
    #     @show i, kmer, count
        MetaGraphs.set_prop!(graph, i, :kmer, kmer)
        MetaGraphs.set_prop!(graph, i, :count, count)
    end
    
    @info "indexing kmers"
    # allow graph[kmer, :kmer] to dict-lookup the index of a kmer
    MetaGraphs.set_indexing_prop!(graph, :kmer)
    
    @info "adding records to graph"
    graph = add_fastx_records_to_graph!(graph, fastxs)

    @info "adding edges and edge counts"
    graph = add_record_edgemers_to_graph!(graph)
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a k-mer graph from a single FASTX format string.

# Arguments
- `KMER_TYPE`: The k-mer type specification (e.g., DNAKmer{K} where K is k-mer length)
- `fastx::AbstractString`: Input sequence in FASTX format (FASTA or FASTQ)

# Returns
- `KmerGraph`: A directed graph where vertices are k-mers and edges represent overlaps
"""
function fastx_to_kmer_graph(KMER_TYPE, fastx::AbstractString)
    return fastx_to_kmer_graph(KMER_TYPE, [fastx])
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add FASTX records from multiple files as a graph property.

# Arguments
- `graph`: A MetaGraph that will store the FASTX records
- `fastxs`: Collection of FASTA/FASTQ file paths to process

# Details
Creates a dictionary mapping sequence descriptions to their corresponding FASTX records,
then stores this dictionary as a graph property under the key `:records`.
Multiple input files are merged, with later files overwriting records with duplicate descriptions.

# Returns
The modified graph with added records property.
"""
function add_fastx_records_to_graph!(graph, fastxs)
    record_dict = Dict(String(FASTX.description(record)) => record for record in Mycelia.open_fastx(first(fastxs)))
    for fastx in fastxs[2:end]
        for record in Mycelia.open_fastx(fastx)
            record_dict[String(FASTX.description(record))] = record
        end
    end
    MetaGraphs.set_prop!(graph, :records, record_dict)
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Processes DNA sequence records stored in the graph and adds their edgemers (k+1 length subsequences) 
to build the graph structure.

# Arguments
- `graph`: A Mycelia graph object containing DNA sequence records and graph properties

# Details
- Uses the k-mer size specified in `graph.gprops[:k]` to generate k+1 length edgemers
- Iterates through each record in `graph.gprops[:records]`
- For each record, generates all possible overlapping edgemers
- Adds each edgemer to the graph with its position and record information

# Returns
- The modified graph object with added edgemer information
"""
function add_record_edgemers_to_graph!(graph)
    edgemer_size = graph.gprops[:k] + 1
    ProgressMeter.@showprogress for (record_id, record) in graph.gprops[:records]
        for (index, observed_edgemer) in Kmers.EveryKmer{Kmers.DNAKmer{edgemer_size}}(BioSequences.LongDNA{2}(FASTX.sequence(record)))
            add_edgemer_to_graph!(graph, record_id, index, observed_edgemer)
        end
    end
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert an edgemer to two vertex kmers.

This function takes an edgemer (a sequence of DNA nucleotides) and converts it into two vertex kmers. 
A kmer is a substring of length k from a DNA sequence. The first kmer is created from the first 
n-1 elements of the edgemer, and the second kmer is created from the last n-1 elements of the edgemer.

# Arguments
- `edgemer::AbstractVector{T}`: A vector of DNA nucleotides where T is a subtype of `BioSequences.DNAAlphabet{2}`.

# Returns
- `Tuple{Kmers.Kmer{BioSequences.DNAAlphabet{2}}, Kmers.Kmer{BioSequences.DNAAlphabet{2}}}`: A tuple containing two kmers.
"""
function edgemer_to_vertex_kmers(edgemer)
    a = Kmers.Kmer{BioSequences.DNAAlphabet{2}}(collect(edgemer[i] for i in 1:length(edgemer)-1))
    b = Kmers.Kmer{BioSequences.DNAAlphabet{2}}(collect(edgemer[i] for i in 2:length(edgemer)))
    return a, b
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add an observed edgemer to a graph with its associated metadata.

# Arguments
- `graph::MetaGraph`: The graph to modify
- `record_identifier`: Identifier for the source record
- `index`: Position where edgemer was observed
- `observed_edgemer`: The biological sequence representing the edgemer

# Details
Processes the edgemer by:
1. Splitting it into source and destination kmers
2. Converting kmers to their canonical forms
3. Creating or updating an edge with orientation metadata
4. Storing observation details (record, position, orientation)

# Returns
Modified graph with the new edge and metadata

# Note
If the edge already exists, the observation is added to the existing metadata.
"""
function add_edgemer_to_graph!(graph, record_identifier, index, observed_edgemer)
    # observed_orientation = BioSequences.iscanonical(observed_edgemer)
    # canonical_edgemer = BioSequences.canonical(observed_edgemer)
    observed_source_kmer, observed_destination_kmer = Mycelia.edgemer_to_vertex_kmers(observed_edgemer)
    
    canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
    canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
    
    source_kmer_orientation = BioSequences.iscanonical(observed_source_kmer)
    destination_kmer_orientation = BioSequences.iscanonical(observed_destination_kmer)
    
    source_kmer_index = graph[canonical_source_kmer, :kmer]
    destination_kmer_index = graph[canonical_destination_kmer, :kmer]
    
    edgemer = Graphs.Edge(source_kmer_index, destination_kmer_index)
    orientation = (source_kmer_orientation => destination_kmer_orientation)
    observation = (;record_identifier, index, orientation)
    
    if !Graphs.has_edge(graph, edgemer)
        Graphs.add_edge!(graph, edgemer)
        MetaGraphs.set_prop!(graph, edgemer, :observations, Set([observation]))
    else
        observations = push!(MetaGraphs.get_prop(graph, edgemer, :observations), observation)
        MetaGraphs.set_prop!(graph, edgemer, :observations, observations)
    end
    return graph
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_fastx_to_graph!(graph, fastx_file::AbstractString)
#     for record in Mycelia.open_fastx(fastx_file)
#         add_fastx_record_to_graph!(graph, record)
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_fastx_record_to_graph!(graph, record::FASTX.FASTA.Record)
#     try
#         graph[FASTX.identifier(record), :identifier]
#         @info "node $(FASTX.identifier(record)) already present"
#     catch
#         Graphs.add_vertex!(graph)
#         vertex_id = Graphs.nv(graph)

#         MetaGraphs.set_prop!(graph, vertex_id, :TYPE, typeof(record))

#         MetaGraphs.set_prop!(graph, vertex_id, :identifier, FASTX.identifier(record))

#         MetaGraphs.set_prop!(graph, vertex_id, :description, FASTX.description(record))

#         sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
#         MetaGraphs.set_prop!(graph, vertex_id, :sequence, sequence)
#     end
#     return graph
# end
    
# function add_fastx_record_to_graph!(graph, record::FASTX.FASTQ.Record)
#     Graphs.add_vertex!(graph)
#     vertex_id = Graphs.nv(graph)
    
#     MetaGraphs.set_prop!(graph, vertex_id, :TYPE, typeof(record))
    
#     MetaGraphs.set_prop!(graph, vertex_id, :identifier, FASTX.identifier(record))
    
#     MetaGraphs.set_prop!(graph, vertex_id, :description, FASTX.description(record))
    
#     sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
#     MetaGraphs.set_prop!(graph, vertex_id, :sequence, sequence)
    
#     MetaGraphs.set_prop!(graph, vertex_id, :quality, FASTX.quality_scores(record))
    
#     return graph
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# # for kmer_size in kmer_sizes
# function add_fasta_record_kmers_to_graph!(graph, kmer_size)
#     record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
#     for vertex in record_vertices
#         record_identifier = graph.vprops[vertex][:identifier]
#         record_sequence = graph.vprops[vertex][:sequence]
#         kmer_counts = Mycelia.count_canonical_kmers(Kmers.Kmer{BioSequences.DNAAlphabet{4},kmer_size}, record_sequence)
#         Mycelia.add_kmers_to_graph!(graph, keys(kmer_counts))
#         Mycelia.add_record_kmer_counts_to_graph!(graph, kmer_counts, record_identifier)
#         Mycelia.add_record_edgemers_to_graph!(graph, record_identifier, kmer_size)
#     end
#     return graph
# end



# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_record_kmer_counts_to_graph!(graph, kmer_counts, record_identifier)
#     for (kmer, count) in kmer_counts
#         kmer_vertex = graph[kmer, :identifier]
#         record_vertex = graph[record_identifier, :identifier]
#         edge = Graphs.Edge(kmer_vertex, record_vertex)
#         if !Graphs.has_edge(graph, edge)
#             Graphs.add_edge!(graph, edge)
#             MetaGraphs.set_prop!(graph, edge, :TYPE, "RECORD_KMER_COUNT")
#             MetaGraphs.set_prop!(graph, edge, :count, count)
#         else
#             graph_count = MetaGraphs.get_prop(graph, edge, :count)
#             if graph_count != count
#                 @warn "edge found but this count $(count) != current count $(graph_count)"
#             # else
#                 # @info "edge exists and matches current data"
#             end
#         end
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_kmers_to_graph!(graph, kmers)
#     for kmer in kmers
#         if !haskey(graph.metaindex[:identifier], kmer)
#             Graphs.add_vertex!(graph)
#             v = Graphs.nv(graph)
#             MetaGraphs.set_prop!(graph, v, :identifier, kmer)
#             MetaGraphs.set_prop!(graph, v, :sequence, kmer)
#             MetaGraphs.set_prop!(graph, v, :TYPE, typeof(kmer))
#         end
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_from_table!(
#         graph::MetaGraphs.AbstractMetaGraph,
#         table::DataFrames.AbstractDataFrame;
#         identifier_column::Union{Symbol, AbstractString} = :identifier)
#     for row in DataFrames.eachrow(table)
#         add_metadata_from_table_row!(graph, row, identifier_column)
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_from_table_row!(graph, row, identifier_column)
#     other_columns = filter(n -> n != :identifier, Symbol.(names(row)))
#     row_metadata_dict = Dict(column => row[column] for column in other_columns)
#     metadata_dict = Dict(row[identifier_column] => row_metadata_dict)
#     Mycelia.add_metadata_to_graph!(graph, metadata_dict::AbstractDict)
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_key_value_pair_to_node!(graph, identifier, key, value)
#     node_id = graph[identifier, :identifier]
#     MetaGraphs.set_prop!(graph, node_id, Symbol(key), value)
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_to_node!(graph, identifier, metadata::AbstractVector{<:Pair})
#     for (key, value) in metadata
#         add_key_value_pair_to_node!(graph, identifier, key, value)
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_to_graph!(graph, metadata_dict::AbstractDict)
#     for (identifier, metadata) in metadata_dict
#         add_metadata_to_node!(graph, identifier, collect(metadata))
#     end
#     return graph
# end

# need to get identifier column and then all non-identifier columns and then pass that to above
# function add_metadata_to_graph!(graph, metadata_table::DataFrames.AbstractDataFrame)
#     for (identifier, metadata) in metadata_dict
#         add_metadata_to_node!(graph, identifier, collect(metadata))
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function initialize_graph()
#     graph = MetaGraphs.MetaDiGraph()
#     MetaGraphs.set_indexing_prop!(graph, :identifier)
#     return graph
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function construct(KMER_TYPE, fastx, out)
#     mkpath(dirname(out))
#     if !occursin(r"\.jld2$", out)
#         out *= ".jld2"
#     end
#     if !isfile(out)
#         graph = fastx_to_kmer_graph(KMER_TYPE, fastx)
#         @info "saving graph"
#         FileIO.save(out, Dict("graph" => graph))
#         return graph
#     else
#         @info "graph $out already exists, loading existing"
#         return load_graph(out)
#     end
# end

# function construct(args)
#     @show args
#     @assert (0 < args["k"] < 64) && isodd(args["k"]) 
#     KMER_TYPE = BioSequences.BigDNAMer{args["k"]}
#     construct(KMER_TYPE, args["fastx"], args["out"])
# end


# # @inline function add_edge_to_simple_kmer_graph!(simple_kmer_graph, kmers, sequence_edge)
# @inline function add_edge_to_simple_kmer_graph!(simple_kmer_graph, sequence_edge)
#     observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.fw)
#     canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
#     canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
#     source_kmer_index = simple_kmer_graph[canonical_source_kmer, :kmer]
#     desination_kmer_index = simple_kmer_graph[canonical_destination_kmer, :kmer]
    
#     if source_kmer_index > desination_kmer_index
#         observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.bw)
#         canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
#         canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
#         source_kmer_index = simple_kmer_graph[canonical_source_kmer, :kmer]
#         desination_kmer_index = simple_kmer_graph[canonical_destination_kmer, :kmer]
#     end
    
#     @assert source_kmer_index <= desination_kmer_index

#     oriented_source_kmer = 
#         (canonical_kmer = canonical_source_kmer,
#          orientation = BioSequences.iscanonical(observed_source_kmer))

#     oriented_destination_kmer = 
#         (canonical_kmer = canonical_destination_kmer,
#          orientation = BioSequences.iscanonical(observed_destination_kmer))

#     oriented_source_vertex = 
# #         (vertex = searchsortedfirst(kmers, oriented_source_kmer.canonical_kmer),
#         (vertex = simple_kmer_graph[oriented_source_kmer.canonical_kmer, :kmer],
#          orientation = oriented_source_kmer.orientation)

#     oriented_destination_vertex = 
# #         (vertex = searchsortedfirst(kmers, oriented_destination_kmer.canonical_kmer),
#         (vertex = simple_kmer_graph[oriented_destination_kmer.canonical_kmer, :kmer],
#          orientation = oriented_destination_kmer.orientation)

# #     forward_edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)
# #     forward_edge_orientations = 
# #         (source_orientation = oriented_source_vertex.orientation,
# #          destination_orientation = oriented_destination_vertex.orientation)
    
# #     reverse_edge = Graphs.Edge(oriented_destination_vertex.vertex, oriented_source_vertex.vertex)
# #     reverse_edge_orientations = 
# #         (source_orientation = !oriented_destination_vertex.orientation,
# #          destination_orientation = !oriented_source_vertex.orientation)
    
#     edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)
#     edge_orientations = 
#         (source_orientation = oriented_source_vertex.orientation,
#          destination_orientation = oriented_destination_vertex.orientation)
    
# #     orientations = Set([forward_edge_orientations, reverse_edge_orientations])
#     orientations = Set([edge_orientations])
#     if Graphs.has_edge(simple_kmer_graph, edge)
#         edge_weight = MetaGraphs.get_prop(simple_kmer_graph, edge, :weight) + 1
#         orientations = union(MetaGraphs.get_prop(simple_kmer_graph, edge, :orientations), orientations)
#     else
#         Graphs.add_edge!(simple_kmer_graph, edge)
#         edge_weight = 1
# #         @show Graphs.ne(simple_kmer_graph)
# #         @show forward_edge
#     end
#     MetaGraphs.set_prop!(simple_kmer_graph, edge, :weight, edge_weight)
#     MetaGraphs.set_prop!(simple_kmer_graph, edge, :orientations, orientations)
# end

# function fastx_to_simple_kmer_graph(KMER_TYPE, fastx::AbstractString; minimum_coverage::Int=1)
#     fastx_to_simple_kmer_graph(KMER_TYPE, [fastx], minimum_coverage=minimum_coverage)
# end

# function fastx_to_simple_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString}; minimum_coverage::Int=1)
#     @info "counting kmers"
#     canonical_kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
#     # hard filter any nodes that are less frequent than minimum coverage threshold
#     canonical_kmer_counts = filter(canonical_kmer_count -> last(canonical_kmer_count) >= minimum_coverage, canonical_kmer_counts)
#     simple_kmer_graph = MetaGraphs.MetaDiGraph(length(canonical_kmer_counts))
    
#     k = length(first(keys(canonical_kmer_counts)))

#     MetaGraphs.set_prop!(simple_kmer_graph, :k, k)

#     @info "setting metadata on vertices"
#     ProgressMeter.@showprogress for (vertex, (kmer, count)) in enumerate(canonical_kmer_counts)
#         MetaGraphs.set_prop!(simple_kmer_graph, vertex, :kmer, kmer)
#         MetaGraphs.set_prop!(simple_kmer_graph, vertex, :weight, count)
#     end

#     kmers = collect(keys(canonical_kmer_counts))

#     EDGE_MER = BioSequences.BigDNAMer{k+1}
#     @info "loading fastx files into graph"
#     ProgressMeter.@showprogress for fastx in fastxs
#         n_records = 0
#         for record in (open_fastx(fastx))
#             n_records += 1
#         end
#         p = ProgressMeter.Progress(n_records, 1)   # minimum update interval: 1 second
#         for record in (open_fastx(fastx))
#             sequence = FASTX.sequence(record)
#             edge_iterator = BioSequences.each(EDGE_MER, sequence)
#             for sequence_edge in edge_iterator
#                 add_edge_to_simple_kmer_graph!(simple_kmer_graph, kmers, sequence_edge)
#             end
#             ProgressMeter.next!(p)
#         end
#     end
#     return simple_kmer_graph
# end



# @inline function add_edge_to_kmer_graph!(kmer_graph, kmers, sequence_edge, record_identifier)
#     observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.fw)

#     oriented_source_kmer = 
#         (canonical_kmer = BioSequences.canonical(observed_source_kmer),
#          orientation = BioSequences.iscanonical(observed_source_kmer))

#     oriented_destination_kmer = 
#         (canonical_kmer = BioSequences.canonical(observed_destination_kmer),
#          orientation = BioSequences.iscanonical(observed_destination_kmer))

#     oriented_source_vertex = 
#         (vertex = searchsortedfirst(kmers, oriented_source_kmer.canonical_kmer),
#          orientation = oriented_source_kmer.orientation)

#     oriented_destination_vertex = 
#         (vertex = searchsortedfirst(kmers, oriented_destination_kmer.canonical_kmer),
#          orientation = oriented_destination_kmer.orientation)

#     source_evidence = 
#         (record = record_identifier,
#          index = sequence_edge.position,
#          orientation = oriented_source_vertex.orientation)

#     destination_evidence = 
#         (record = record_identifier,
#          index = sequence_edge.position + 1,
#          orientation = oriented_destination_vertex.orientation)

#     set_metadata!(kmer_graph, oriented_source_vertex.vertex, :evidence, source_evidence)
#     new_weight = length(kmer_graph.vprops[oriented_source_vertex.vertex][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, oriented_source_vertex.vertex, :weight, new_weight)

#     set_metadata!(kmer_graph, oriented_destination_vertex.vertex, :evidence, destination_evidence)
#     new_weight = length(kmer_graph.vprops[oriented_destination_vertex.vertex][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, oriented_destination_vertex.vertex, :weight, new_weight)
    

#     forward_edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)

#     Graphs.add_edge!(kmer_graph, forward_edge)

#     forward_edge_orientations = 
#         (source_orientation = oriented_source_vertex.orientation,
#          destination_orientation = oriented_destination_vertex.orientation)

#     set_metadata!(kmer_graph, forward_edge, :orientations, forward_edge_orientations)

#     forward_edge_evidence = (
#         record = record_identifier,
#         index = sequence_edge.position,
#         orientation = true
#     )

#     set_metadata!(kmer_graph, forward_edge, :evidence, forward_edge_evidence)
#     new_weight = length(kmer_graph.eprops[forward_edge][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, forward_edge, :weight, new_weight)

#     reverse_edge = Graphs.Edge(oriented_destination_vertex.vertex, oriented_source_vertex.vertex)

#     Graphs.add_edge!(kmer_graph, reverse_edge)

#     reverse_edge_orientations = 
#         (source_orientation = !oriented_destination_vertex.orientation,
#          destination_orientation = !oriented_source_vertex.orientation)

#     set_metadata!(kmer_graph, reverse_edge, :orientations, reverse_edge_orientations)

#     reverse_edge_evidence = (
#         record = record_identifier,
#         index = sequence_edge.position,
#         orientation = false
#     )

#     set_metadata!(kmer_graph, reverse_edge, :evidence, reverse_edge_evidence)
#     new_weight = length(kmer_graph.eprops[reverse_edge][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, reverse_edge, :weight, new_weight)
# end

# # function fastx_to_kmer_graph(::Type{KMER_TYPE}, fastxs) where {KMER_TYPE <: BioSequences.AbstractMer{A, K}} where {A, K}
# function fastx_to_kmer_graph(KMER_TYPE, fastx::AbstractString)
#     fastx_to_kmer_graph(KMER_TYPE, [fastx])
# end

# function fastx_to_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString})
#     @info "assessing kmers"
#     kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
# #     kmer_set = Set{KMER_TYPE}()
# #     for fastxs in fastxs
# #         kmer_set = union!(kmer_set, collect(keys(count_canonical_kmers(KMER_TYPE, fastxs))))
# #     end
# #     kmers = unique(sort(collect(kmer_set)))
    
#     kmer_graph = MetaGraphs.MetaDiGraph(length(kmer_counts))
#     k = length(first(keys(kmer_counts)))
#     kmers = collect(keys(kmer_counts))
#     MetaGraphs.set_prop!(kmer_graph, :k, k)
#     # don't set this since when we filter an induced subgraph, these don't update
# #     MetaGraphs.set_prop!(kmer_graph, :kmers, kmers)
#     for (vertex, (kmer, count)) in enumerate(kmer_counts)
#         MetaGraphs.set_prop!(kmer_graph, vertex, :kmer, kmer)
#         MetaGraphs.set_prop!(kmer_graph, vertex, :weight, count)
#     end
#     EDGE_MER = BioSequences.BigDNAMer{k+1}
#     @info "creating graph"
#     ProgressMeter.@showprogress for fastx in fastxs
#         fastx_io = open_fastx(fastx)
#         for record in fastx_io
#             sequence = FASTX.sequence(record)
#             record_identifier = FASTX.identifier(record) 
#             edge_iterator = BioSequences.each(EDGE_MER, sequence)
#             for sequence_edge in edge_iterator
#                 add_edge_to_kmer_graph!(kmer_graph, kmers, sequence_edge, record_identifier)
#             end
#         end
#         close(fastx_io)
#     end
#     return kmer_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_edge_to_graph(graph, edge_mer, kmers)
#     edge = BioSequences.LongDNASeq(edge_mer.fw)
#     k = length(first(kmers))
# #     canonical_src = BioSequences.DNAMer{k}(BioSequences.canonical!(edge[1:end-1]))
# #     canonical_dst = BioSequences.DNAMer{k}(BioSequences.canonical!(edge[2:end]))

#     canonical_src = BioSequences.canonical(BioSequences.DNAMer{k}(edge[1:end-1]))
    
#     src_index_range = searchsorted(kmers, canonical_src)
#     if isempty(src_index_range)
#         return
#     else
#         @assert length(src_index_range) == 1
#     end
#     src_index = first(src_index_range)

#     canonical_dst = BioSequences.canonical(BioSequences.DNAMer{k}(edge[2:end]))
#     dst_index_range = searchsorted(kmers, canonical_dst)
#     if isempty(dst_index_range)
#         return
#     else
#         @assert length(dst_index_range) == 1
#     end
#     dst_index = first(dst_index_range)
#     graph_edge = Graphs.Edge(src_index, dst_index)
#     Graphs.add_edge!(graph, graph_edge)
# end



################################################################################
# defunct bcalm usage
# run(`conda install -c bioconda bcalm`)

# fasta_file_list = joinpath(working_directory, repr(hash(fastx_files)) * ".fasta_list.txt")
# open(fasta_file_list, "w") do io
#     for f in fastx_files
#         @show f
#         println(io, f)
#     end
# end

# k = 3
# outfile = fasta_file_list * ".bcalm.$(k).fna"
# cmd = `bcalm -in $(fastx_files[1]) -abundance-min 1 -kmer-size 11 -all-abundance-counts -out $(outfile)`
# run(cmd)

# cmds = [
#     `bcalm`,
#     `-in $(fasta_list_file)`,
#     `-abundance-min 1`,
#     `-kmer-size 3`,
#     `-all-abundance-counts`,
#     `-abundance-max $(typemax(UInt64))`
# ]
# run(cmds)

# ls -1 *.fastq > list_reads
# ./bcalm -in list_reads [..]
# ./bcalm -in [reads.fa] -kmer-size [kmer_size] -abundance-min [abundance_threshold]

# scripts/convertToGFA.py
##################################################################################

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function graph_to_kmers(g)
#     kmers = [g.vprops[v][:kmer] for v in Graphs.vertices(g)]
#     return kmers
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function graph_to_edge_sequences(g)
#     edges = Set{BioSequences.BigDNAMer{g.gprops[:k]+1}}()
#     for edge in Graphs.edges(g)
#         src_kmer = g.vprops[edge.src][:kmer]
#         dst_kmer = g.vprops[edge.dst][:kmer]
#         for orientation in g.eprops[edge][:orientations]
#             if orientation.source_orientation
#                 oriented_src_kmer = src_kmer
#             else
#                 oriented_src_kmer = BioSequences.reverse_complement(src_kmer)
#             end
#             if orientation.destination_orientation
#                 oriented_dst_kmer = dst_kmer
#             else
#                 oriented_dst_kmer = BioSequences.reverse_complement(dst_kmer)
#             end
#             for i in 1:g.gprops[:k]-1
#                 @assert oriented_src_kmer[i+1] == oriented_dst_kmer[i]
#             end
#             edge_mer = BioSequences.BigDNAMer((nuc for nuc in oriented_src_kmer)..., last(oriented_dst_kmer))
#             push!(edges, BioSequences.canonical(edge_mer))
#         end
#     end
#     return edges
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function kmer_graph_distances(g1, g2)
#     g1_kmers = Set(graph_to_kmers(g1))
#     g1_edges = graph_to_edge_sequences(g1)
    
#     g2_kmers = Set(graph_to_kmers(g2))
#     g2_edges = graph_to_edge_sequences(g2)
    
#     kmer_distance = 1 - LSHFunctions.jaccard(g1_kmers, g2_kmers)
#     edge_distance = 1 - LSHFunctions.jaccard(g1_edges, g2_edges)
    
#     result = (
#         kmer_distance = kmer_distance,
#         edge_distance = edge_distance
#     )
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function set_metadata!(kmer_graph, vertex::V, key, value) where V <: Integer
#     if MetaGraphs.has_prop(kmer_graph, vertex, key)
#         push!(kmer_graph.vprops[vertex][key], value)
#     else
#         MetaGraphs.set_prop!(kmer_graph, vertex, key, Set([value]))
#     end
#     return true
# end

# function set_metadata!(kmer_graph, edge::E, key, value) where E <: Graphs.Edge
#     if MetaGraphs.has_prop(kmer_graph, edge, key)
#         current_value = MetaGraphs.get_prop(kmer_graph, edge, key)
#         updated_value = push!(current_value, value)
#         MetaGraphs.set_prop!(kmer_graph, edge, key, updated_value)
#     else
#         MetaGraphs.set_prop!(kmer_graph, edge, key, Set([value]))
#     end
#     return true
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simple_polish_fastq(simple_kmer_graph, fastq_file; min_depth=3)
#     solid_vertices = filter(v -> simple_kmer_graph.vprops[v][:weight] >= min_depth, Graphs.vertices(simple_kmer_graph))
#     filtered_simple_kmer_graph, vertex_map = Graphs.induced_subgraph(simple_kmer_graph, solid_vertices)
# #     display(simple_kmer_graph)
# #     display(filtered_simple_kmer_graph)
#     kmers = sort!(graph_to_kmers(filtered_simple_kmer_graph))
    
#     old_kmers = graph_to_kmers(simple_kmer_graph)
    
#     k = filtered_simple_kmer_graph.gprops[:k]
    
#     polished_fastq_file = replace(fastq_file, ".fastq" => ".k$k.d$(min_depth).fastq")

#     transition_probabilities = initialize_transition_probabilities(filtered_simple_kmer_graph)
#     state_likelihoods = [Float64(filtered_simple_kmer_graph.vprops[v][:weight]) for v in Graphs.vertices(filtered_simple_kmer_graph)]
#     state_likelihoods ./= sum(state_likelihoods)
    
# #     @info "counting the number of records to establish runtime estimate"
#     number_of_records = 0
#     for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
#         number_of_records += 1
#     end
#     progress_bar = ProgressMeter.Progress(number_of_records, 1)

#     fastq_reader = FASTX.FASTQ.Reader(open(fastq_file))
#     fastq_writer = FASTX.FASTQ.Writer(open(polished_fastq_file, "w"))

#     for fastx_record in fastq_reader
#         ProgressMeter.next!(progress_bar)
#         bubble_start = 0
#         updated_path = Vector{Pair{Int, Bool}}()
#     #     @show FASTX.sequence(fastx_record)
#         for (i, kmer) in enumerate(BioSequences.each(BioSequences.BigDNAMer{k}, FASTX.sequence(fastx_record)))
#             canonical_kmer = min(kmer.fw, kmer.bw)
#             orientation = canonical_kmer == kmer.fw
#             kmer_index_range = searchsorted(kmers, canonical_kmer)
            
#             old_kmer_index = searchsortedfirst(old_kmers, canonical_kmer)
            
# #             FASTX.identifier(fastx_record) == "4" && @show kmer_index_range
# #             FASTX.identifier(fastx_record) == "4" && @show kmers[kmer_index_range]
# #             FASTX.identifier(fastx_record) == "4" && @show old_kmers[old_kmer_index]
#             kmer_is_solid = !isempty(kmer_index_range)
# #             FASTX.identifier(fastx_record) == "4" && @show canonical_kmer
# #             FASTX.identifier(fastx_record) == "4" &&  @show kmer_is_solid

#             if kmer_is_solid
#                 kmer_index = first(kmer_index_range)
# #                 FASTX.identifier(fastx_record) == "4" && @show filtered_simple_kmer_graph.vprops[kmer_index], simple_kmer_graph.vprops[old_kmer_index]
#             else
#                 kmer_index = 0
#             end

#             in_bubble = bubble_start > 0

#             if !kmer_is_solid
#                 if !in_bubble
# #                     FASTX.identifier(fastx_record) == "4" && @show "starting a bubble"
#                     bubble_start = i
#                 else
# #                     FASTX.identifier(fastx_record) == "4" && @show "continuing in a bubble"
#                 end
#             else
#                 if !in_bubble
# #                     FASTX.identifier(fastx_record) == "4" && @show "pushing solid kmer to updated path"
#                     push!(updated_path, kmer_index => orientation)
#                 else
#                     if bubble_start == 1
# #                         FASTX.identifier(fastx_record) == "4" && @show "ending an opening bubble"
#                         # we're in a bubble that began at the beginning of the read
#                         # we'll do nothing and just remove this
#                         # equivalent to tip clipping
# #                         FASTX.identifier(fastx_record) == "4" && @show "pushing solid kmer to updated path"
#                         push!(updated_path, kmer_index => orientation)
#                         bubble_start = 0
#                     else
# #                         FASTX.identifier(fastx_record) == "4" && @show "found end of an anchored bubble -- correcting"
#                         source_vertex, source_orientation = last(updated_path)
#                         destination_vertex, destination_orientation = kmer_index, orientation                

#                         shortest_paths = Graphs.yen_k_shortest_paths(
#                             filtered_simple_kmer_graph,
#                             source_vertex,
#                             destination_vertex,
#                             Graphs.weights(filtered_simple_kmer_graph),
#                             10).paths

#                         if isempty(shortest_paths)
#                             @show source_vertex, destination_vertex
#                             @warn "no valid alternate paths found for $(FASTX.identifier(fastx_record)), continuing"
#                             break
# #                             @show fastx_record
# #                             error("try increasing min_depth")
#                         end
#                         candidate_path_probabilities = ones(length(shortest_paths))
#                         oriented_candidate_paths = [
#                             [last(updated_path)] for i in 1:length(shortest_paths)
#                         ]

#                         for (i, candidate_path) in enumerate(shortest_paths)
#                             for dest_vertex in candidate_path[2:end]
#                                 source_vertex, source_orientation = last(oriented_candidate_paths[i])
#                                 candidate_path_probabilities[i] *= transition_probabilities[source_orientation][source_vertex, dest_vertex]
#                                 candidate_path_probabilities[i] *= state_likelihoods[dest_vertex]
#                                 if candidate_path_probabilities[i] > 0
#                                     edge = Graphs.Edge(source_vertex, dest_vertex)
#                                     destination_orientation = 
#                                     first(
#                                         filter(o -> o.source_orientation == source_orientation,
#                                             filtered_simple_kmer_graph.eprops[edge][:orientations])).destination_orientation
#                                     push!(oriented_candidate_paths[i], (dest_vertex => destination_orientation))
#                                 else
#                                     break # this path is no good, evaluate the next
#                                 end
#                             end
#                         end
#                         non_zero_indices = findall(p -> p > 0, candidate_path_probabilities)
#                         if isempty(non_zero_indices)
# #                             @show candidate_path_probabilities
#                             @warn "no resolution found for read $(FASTX.identifier(fastx_record))"
#                             break
# #                             error("no valid alternate path probabilities")
# #                             error("try increasing min_depth?")
#                         end

#                         candidate_path_probabilities = candidate_path_probabilities[non_zero_indices]
#                         oriented_candidate_paths = oriented_candidate_paths[non_zero_indices]

#                         # offset is for debugging
#                         # make sure that anchors on both sides are the same
#                         offset = 0
#                         observed_sequence = FASTX.sequence(fastx_record)[bubble_start+k-1-offset:i-1+offset]                    
#                         for (i, oriented_candidate_path) in enumerate(oriented_candidate_paths)
#                             candidate_sequence = oriented_path_to_sequence(
#                                 filtered_simple_kmer_graph, 
#                                 oriented_candidate_path)
#                             candidate_sequence = candidate_sequence[k+1-offset:end-k+offset]
#                             alignment_result = BioAlignments.pairalign(
#                                 BioAlignments.LevenshteinDistance(),
#                                 candidate_sequence,
#                                 observed_sequence)
# #                             @show alignment_result
#     #                         @show alignment_result
#                             average_error_rate = Statistics.mean(q_value_to_error_rate.(FASTX.quality(fastx_record)))
#                             for error in 1:alignment_result.value
#                                 candidate_path_probabilities[i] *= average_error_rate
#                             end
#                             for match in 1:BioAlignments.count_matches(alignment_result.aln)
#                                 candidate_path_probabilities[i] *= (1 - average_error_rate)
#                             end
#                         end

#                         chosen_replacement = StatsBase.sample(oriented_candidate_paths, StatsBase.weights(candidate_path_probabilities))

#                         for i in 2:length(chosen_replacement)
#                             oriented_state = chosen_replacement[i]
#                             push!(updated_path, oriented_state)
#                         end
#                         bubble_start = 0
#                     end
#                 end
#             end
#         end
#     #     @show updated_path
#         sequence = oriented_path_to_sequence(filtered_simple_kmer_graph, updated_path)
#         alignment_result = BioAlignments.pairalign(
#             BioAlignments.LevenshteinDistance(),
#             sequence,
#             FASTX.sequence(fastx_record))
# #         if alignment_result.value > 0
# # #             @show alignment_result
# #         end
#         quality = StatsBase.sample(FASTX.quality(fastx_record), length(sequence), replace=true, ordered=true)
#         description =  join(filter(!isempty, (FASTX.description(fastx_record), "k$k.d$(min_depth)")), '.')
#         identifier = FASTX.identifier(fastx_record)
#         new_record = FASTX.FASTQ.Record(identifier, description, sequence, quality)
#         write(fastq_writer, new_record)
#     end
#     close(fastq_reader)
#     close(fastq_writer)
#     return polished_fastq_file
# end



# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function assess_observations(graph::KmerGraph{KMER_TYPE}, observations, error_rate; verbose = isinteractive()) where {KMER_TYPE}
#     k = last(KMER_TYPE.parameters)
#     total_edits_accepted = 0
#     total_bases_evaluated = 0
#     reads_processed = 0
#     maximum_likelihood_observations = Vector{BioSequences.LongDNASeq}(undef, length(observations))
#     for (observation_index, observation) in enumerate(observations)
#         if length(observation) >= k
#             optimal_path, edit_distance, relative_likelihood = viterbi_maximum_likelihood_path(graph, observation, error_rate)
#             maximum_likelihood_observation = oriented_path_to_sequence(optimal_path, graph.kmers)
#             maximum_likelihood_observations[observation_index] = maximum_likelihood_observation
#             reads_processed += 1
#             total_bases_evaluated += length(observation)
#             total_edits_accepted += edit_distance
#         else
#             maximum_likelihood_observations[observation_index] = observation
#         end
#     end
#     inferred_error_rate = round(total_edits_accepted / total_bases_evaluated, digits = 3)
#     if verbose
#         display("reads_processed = $(reads_processed)")
#         display("total_edits_accepted = $(total_edits_accepted)")
#         display("inferred_error_rate = $(inferred_error_rate)")
#     end
#     if total_edits_accepted == 0
#         has_converged = true
#     else
#         has_converged = false
#     end
#     return maximum_likelihood_observations, has_converged
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function iterate_until_convergence(ks, observations, error_rate)
#     graph = nothing
#     for k in ks
#         graph = KmerGraph(BioSequences.DNAMer{k}, observations)
#         observations, has_converged = assess_observations(graph, observations, error_rate; verbose = verbose)
#     end
#     return graph, observations
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function clip_low_coverage_tips(graph, observations)
#     connected_components = Graphs.connected_components(graph.graph)
#     vertices_to_keep = Int[]
#     for connected_component in connected_components
        
#         component_coverage = graph.counts[connected_component]
#         median = Statistics.median(component_coverage)
#         standard_deviation = Statistics.std(component_coverage)
        
#         for vertex in connected_component
#             keep = true
#             if Graphs.degree(graph.graph, vertex) == 1
#                 this_coverage = graph.counts[vertex]
#                 is_low_coverage = (graph.counts[vertex] == 1) || 
#                                     (median-this_coverage) > (3*standard_deviation)
#                 if is_low_coverage
#                     keep = false
#                 end
#             end
#             if keep
#                 push!(vertices_to_keep, vertex)
#             end
#         end
#     end
    
#     KmerType = first(typeof(graph).parameters)
#     pruned_graph = KmerGraph(KmerType, observations, graph.kmers[vertices_to_keep], graph.counts[vertices_to_keep])
    
#     return pruned_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simplify_kmer_graph(kmer_graph)
#     @info "simplifying kmer graph"
#     @info "resolving untigs..."
#     @time untigs = resolve_untigs(kmer_graph)
#     @info "determining untig orientations..."
#     oriented_untigs = determine_oriented_untigs(kmer_graph, untigs)
#     simplified_graph = MetaGraphs.MetaDiGraph(length(oriented_untigs))
#     MetaGraphs.set_prop!(simplified_graph, :k, kmer_graph.gprops[:k])
#     @info "initializing graph node metadata"
#     for (vertex, untig) in enumerate(oriented_untigs)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :sequence, untig.sequence)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :path, untig.path)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :orientations, untig.orientations)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :evidence, untig.evidence)
#     end
    
#     # determine oriented edges of simplified graph
#     simplified_untigs = []
#     @info "creating simplified untigs to resolve connections"
#     for vertex in Graphs.vertices(simplified_graph)
#         in_kmer = simplified_graph.vprops[vertex][:path][1] => simplified_graph.vprops[vertex][:orientations][1]
#         out_kmer = simplified_graph.vprops[vertex][:path][end] => simplified_graph.vprops[vertex][:orientations][end]
#     #     @show vertex, in_kmer, out_kmer
#         push!(simplified_untigs, in_kmer => out_kmer)
#     end

#     @info "resolving connections"
#     ProgressMeter.@showprogress for (ui, u) in enumerate(simplified_untigs)
#         for (vi, v) in enumerate(simplified_untigs)
#     #         + => +
#             source_kmer_index, source_orientation = last(u)
#             destination_kmer_index, destination_orientation = first(v)
#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! + +"

#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         + => -
#             source_kmer_index, source_orientation = last(u)
#             destination_kmer_index, destination_orientation = last(v)
#             destination_orientation = !destination_orientation

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! + -"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         - => +
#             source_kmer_index, source_orientation = first(u)
#             source_orientation = !source_orientation
#             destination_kmer_index, destination_orientation = first(v)

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! - +"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         - => -
#             source_kmer_index, source_orientation = first(u)
#             source_orientation = !source_orientation
#             destination_kmer_index, destination_orientation = last(v)
#             destination_orientation = !destination_orientation

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! - -"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#         end
#     end
#     return simplified_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simplify_kmer_graph(kmer_graph)
#     @info "simplifying kmer graph"
#     @info "resolving untigs..."
#     @time untigs = resolve_untigs(kmer_graph)

#     @info "determining untig orientations..."
#     @time oriented_untigs = determine_oriented_untigs(kmer_graph, untigs)

#     simplified_graph = MetaGraphs.MetaDiGraph(length(oriented_untigs))
#     MetaGraphs.set_prop!(simplified_graph, :k, kmer_graph.gprops[:k])
#     @info "initializing graph node metadata"
#     ProgressMeter.@showprogress for (vertex, untig) in enumerate(oriented_untigs)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :sequence, untig.sequence)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :path, untig.path)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :orientations, untig.orientations)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :weight, untig.weight)
#     end

#     # determine oriented edges of simplified graph
#     simplified_untigs = Vector{Pair{Pair{Int64,Bool},Pair{Int64,Bool}}}(undef, length(Graphs.vertices(simplified_graph)))
#     @info "creating simplified unitgs to help resolve connections"
#     # use a pre-allocated array here to speed up
#     ProgressMeter.@showprogress for vertex in Graphs.vertices(simplified_graph)
#         in_kmer = simplified_graph.vprops[vertex][:path][1] => simplified_graph.vprops[vertex][:orientations][1]
#         out_kmer = simplified_graph.vprops[vertex][:path][end] => simplified_graph.vprops[vertex][:orientations][end]
#     #     @show vertex, in_kmer, out_kmer
#         simplified_untigs[vertex] = in_kmer => out_kmer
#     #     push!(simplified_untigs, )
#     end

#     # make a dictionary mapping endcap to oriented_untig index

#     end_mer_map = Dict()
#     ProgressMeter.@showprogress for (i, oriented_untig) in enumerate(oriented_untigs)
#         end_mer_map[first(oriented_untig.path)] = i
#         end_mer_map[last(oriented_untig.path)] = i
#     end

#     ProgressMeter.@showprogress for (untig_index, oriented_untig) in enumerate(oriented_untigs)
#     #     @show untig_index
#         true_in_overlap = oriented_untig.sequence[1:simplified_graph.gprops[:k]-1]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[1])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[2])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_out_overlap = neighboring_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_true_out_overlap == true_in_overlap
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = true => true
#                 o = (source_orientation = true, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_out_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_false_out_overlap == true_in_overlap        
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = false => true
#                 o = (source_orientation = false, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         true_out_overlap = oriented_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[end])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[end-1])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_in_overlap = neighboring_untig.sequence[1:simplified_graph.gprops[:k]-1]
#             if true_out_overlap == neighbor_true_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = true => true
#                 o = (source_orientation = true, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_in_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[1:simplified_graph.gprops[:k]-1]
#             if true_out_overlap == neighbor_false_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = true => false
#                 o = (source_orientation = true, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         false_in_overlap = BioSequences.reverse_complement(oriented_untig.sequence)[1:simplified_graph.gprops[:k]-1]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[end])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[end-1])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_out_overlap = neighboring_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_true_out_overlap == false_in_overlap
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = true => false
#                 o = (source_orientation = true, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_out_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_false_out_overlap == false_in_overlap        
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = false => false
#                 o = (source_orientation = false, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         false_out_overlap = BioSequences.reverse_complement(oriented_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[1])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[2])
#         end

#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_in_overlap = neighboring_untig.sequence[1:simplified_graph.gprops[:k]-1]
#             if false_out_overlap == neighbor_true_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = false => true
#                 o = (source_orientation = false, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_in_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[1:simplified_graph.gprops[:k]-1]
#             if false_out_overlap == neighbor_false_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = false => false
#                 o = (source_orientation = false, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end
#     end
#     return simplified_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_oriented_untigs(kmer_graph, untigs)
#     oriented_untigs = []
#     for path in untigs
#         sequence = BioSequences.LongDNASeq(kmer_graph.vprops[first(path)][:kmer])
#         if length(path) == 1
#             orientations = [true]
#         elseif length(path) > 1
#             initial_edge = Graphs.Edge(path[1], path[2])
#             initial_orientation = first(kmer_graph.eprops[initial_edge][:orientations]).source_orientation
#             orientations = [initial_orientation]
#             if !initial_orientation
#                 sequence = BioSequences.reverse_complement(sequence)
#             end

#             for (src, dst) in zip(path[1:end-1], path[2:end])
#                 edge = Graphs.Edge(src, dst)
#                 destination = BioSequences.LongDNASeq(kmer_graph.vprops[edge.dst][:kmer])
#                 destination_orientation = first(kmer_graph.eprops[edge][:orientations]).destination_orientation
#                 push!(orientations, destination_orientation)
#                 if !destination_orientation
#                     destination = BioSequences.reverse_complement(destination)
#                 end
#                 sequence_suffix = sequence[end-length(destination)+2:end]
#                 destination_prefix = destination[1:end-1]
#                 @assert sequence_suffix == destination_prefix
#                 push!(sequence, destination[end])
#             end
#         end

#         oriented_untig = 
#         (
#             sequence = BioSequences.canonical(sequence),
#             path = BioSequences.iscanonical(sequence) ? path : reverse(path),
#             orientations = BioSequences.iscanonical(sequence) ? orientations : reverse(.!orientations),
#             weight = Statistics.median([kmer_graph.vprops[v][:weight] for v in path])
#         )

#         push!(oriented_untigs, oriented_untig)
#     end
#     return oriented_untigs
# end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # Description

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# function resolve_untigs(kmer_graph)
#     untigs = Vector{Int}[]
#     visited = falses(Graphs.nv(kmer_graph))
#     first_unvisited = findfirst(!, visited)
#     while first_unvisited != nothing
#         forward_walk = oriented_unbranching_walk(kmer_graph, first_unvisited, true)
#         reverse_walk = oriented_unbranching_walk(kmer_graph, first_unvisited, false)
#         inverted_reverse_walk = [Graphs.Edge(e.dst, e.src) for e in reverse(reverse_walk)]
#         edges = vcat(inverted_reverse_walk, forward_walk)
#         if isempty(edges)
#             untig = [first_unvisited]
#         else
#             untig = vcat([first(edges).src], [edge.dst for edge in edges])
#         end
#         push!(untigs, untig)
#         for vertex in untig
#             visited[vertex] = true
#         end
#         first_unvisited = findfirst(!, visited)
#     end
#     return untigs
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function oriented_unbranching_walk(kmer_graph, vertex, orientation)
#     walk = []
#     viable_neighbors = find_unbranched_neighbors(kmer_graph, vertex, orientation)
#     while length(viable_neighbors) == 1
# #         @show "found a viable neighbor!!"
#         viable_neighbor = first(viable_neighbors)
#         edge = Graphs.Edge(vertex, viable_neighbor)
#         push!(walk, edge)
#         vertex = edge.dst
#         viable_neighbors = Set{Int}()
#         destination_orientations = [o.destination_orientation for o in kmer_graph.eprops[edge][:orientations]]
#         for destination_orientation in destination_orientations
#             union!(viable_neighbors, find_unbranched_neighbors(kmer_graph, vertex, destination_orientation))
#         end
#     end
#     return walk
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
# or get fasta directly from FTP site

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function find_unbranched_neighbors(kmer_graph, vertex, orientation)
#     downstream_vertices = find_downstream_vertices(kmer_graph, vertex, orientation)
# #     backtrack_vertices
# #     @show downstream_vertices
#     if length(downstream_vertices) == 1
#         downstream_vertex = first(downstream_vertices)
# #         @show downstream_vertex
#         edge = Graphs.Edge(vertex, downstream_vertex)
# #         @show edge
#         destination_orientations = [o.destination_orientation for o in kmer_graph.eprops[edge][:orientations]]
# #         @show destination_orientations
#         for destination_orientation in destination_orientations
#             backtrack_vertices = find_downstream_vertices(kmer_graph, downstream_vertex, !destination_orientation)
# #             @show backtrack_vertices
#             # if the only backtrack is the vertex we're on, then we can simplify
#             if backtrack_vertices == Set([vertex])
#                 return downstream_vertices
#             end
#         end
#     end
#     return Int[]
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
# or get fasta directly from FTP site

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function find_downstream_vertices(kmer_graph, vertex, orientation)
#     viable_neighbors = Set{Int}()
#     for neighbor in Graphs.neighbors(kmer_graph, vertex)
#         not_same_vertex = vertex != neighbor
#         candidate_edge = Graphs.Edge(vertex, neighbor)
#         # palindromes can have multiple viable orientations
#         # check each viable orientation individually
#         edge_src_orientations = [e.source_orientation for e in kmer_graph.eprops[candidate_edge][:orientations]]
#         for edge_src_orientation in edge_src_orientations
# #             edge_src_orientation = kmer_graph.eprops[candidate_edge][:orientations].source_orientation
#             viable_orientation = edge_src_orientation == orientation
#             if not_same_vertex && viable_orientation
#                 push!(viable_neighbors, neighbor)
#             end
#         end
#     end
#     return viable_neighbors
# end


# function apply_kmedoids_treshold(graph)
#     kmer_counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]

#     kmer_counts_histogram = sort(collect(StatsBase.countmap(values(kmer_counts))), by=x->x[1])

# #     scale = 250
# #     p = plot_kmer_frequency_spectra(values(kmer_counts), size=(2scale,scale), log_scale=log2, title="kmer frequencies")
# #     display(p)

# #     p = StatsPlots.scatter(log2.(first.(kmer_counts_histogram)))
# #     display(p)

#     kmer_depth_of_coverage_bins = log2.(first.(kmer_counts_histogram))

#     distance_matrix = zeros((length(kmer_depth_of_coverage_bins), length(kmer_depth_of_coverage_bins)))
#     for (row, depth_of_coverage_bin_1) in enumerate(kmer_depth_of_coverage_bins)
#         for (col, depth_of_coverage_bin_2) in enumerate(kmer_depth_of_coverage_bins)
#             distance = abs(depth_of_coverage_bin_1 - depth_of_coverage_bin_2)
#             distance_matrix[row, col] = distance
#         end
#     end
#     distance_matrix

#     # max out k at the same max k we use for DNAMers
#     max_k = min(length(kmer_depth_of_coverage_bins), 63)
#     ks = Primes.primes(2, max_k)
#     ys = map(k ->
#                 Statistics.mean(Statistics.mean(Clustering.silhouettes(Clustering.kmedoids(distance_matrix, k), distance_matrix)) for i in 1:100),
#                 ks)

#     p = StatsPlots.plot(ks, ys, label="silhouette score", ylabel = "silhouette score", xlabel = "number of clusters")
#     display(p)

#     ymax, ymax_index = findmax(ys)
#     optimal_k = ks[ymax_index]
#     clusterings = [Clustering.kmedoids(distance_matrix, optimal_k) for i in 1:10]
#     max_value, max_value_index = findmax(clustering -> Statistics.mean(Clustering.silhouettes(clustering, distance_matrix)), clusterings)
#     optimal_clustering = clusterings[max_value_index]
#     # optimal_clustering.assignments
#     min_medoid_value, min_medoid_index = findmin(optimal_clustering.medoids)
#     indices_to_include = map(assignment -> assignment .!= min_medoid_index, optimal_clustering.assignments)
#     # kmer_depth_of_coverage_bins
#     threshold = Int(ceil(2^maximum(kmer_depth_of_coverage_bins[.!indices_to_include]))) + 1

#     scale = 250
#     p = plot_kmer_frequency_spectra(values(kmer_counts), log_scale = log2, size=(2scale,scale), title="kmer frequencies")
#     StatsPlots.vline!(p, log2.([threshold]))
#     display(p)

#     # find all vertices with count > threshold
#     vertices_to_keep = [v for v in Graphs.vertices(graph) if (MetaGraphs.get_prop(graph, v, :count) > threshold)]
#     # induce subgraph
#     induced_subgraph, vertex_map = Graphs.induced_subgraph(graph, vertices_to_keep)

#     # set kmer as indexing prop
#     MetaGraphs.set_indexing_prop!(induced_subgraph, :kmer)
#     return induced_subgraph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_edge_probabilities(graph, strand)
#     kmers = graph.gprops[:kmers]
#     outgoing_edge_probabilities = SparseArrays.spzeros(length(kmers), length(kmers))
    
#     for (kmer_index, kmer) in enumerate(kmers)
#         if !strand
#             kmer = BioSequences.reverse_complement(kmer)
#         end
        
#         downstream_neighbor_indices = Int[]
#         for neighbor in BioSequences.neighbors(kmer)
#             index = get_kmer_index(kmers, BioSequences.canonical(neighbor))
#             # kmer must be in our dataset and there must be a connecting edge
#             if !isnothing(index) && Graphs.has_edge(graph, ordered_edge(kmer_index, index))
#                 push!(downstream_neighbor_indices, index)
#             end
#         end
#         sort!(unique!(downstream_neighbor_indices))
        
#         downstream_edge_weights = Int[
#             length(get(graph.edge_evidence, ordered_edge(kmer_index, neighbor_index), EdgeEvidence[])) for neighbor_index in downstream_neighbor_indices
#         ]
        
#         non_zero_indices = downstream_edge_weights .> 0
#         downstream_neighbor_indices = downstream_neighbor_indices[non_zero_indices]
#         downstream_edge_weights = downstream_edge_weights[non_zero_indices]
        
#         downstream_edge_likelihoods = downstream_edge_weights ./ sum(downstream_edge_weights)
        
#         for (neighbor_index, likelihood) in zip(downstream_neighbor_indices, downstream_edge_likelihoods)
#             outgoing_edge_probabilities[kmer_index, neighbor_index] = likelihood
#         end
#     end
#     return outgoing_edge_probabilities
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_edge_probabilities(graph)
#     outgoing_edge_probabilities = determine_edge_probabilities(graph, true)
#     incoming_edge_probabilities = determine_edge_probabilities(graph, false)
#     return outgoing_edge_probabilities, incoming_edge_probabilities
# end





# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function reverse_oriented_path(oriented_path)
# #     reversed_path = copy(oriented_path)
# #     for (index, state) in enumerate(oriented_path)
# #         reversed_path[index] = OrientedKmer(index = state.index, orientation = !state.orientation)
# #     end
# #     return reverse!(reversed_path)
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function maximum_likelihood_walk(graph, connected_component)
# #     max_count = maximum(graph.counts[connected_component])
# #     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
# #     initial_node_index = rand(max_count_indices)
# #     initial_node = connected_component[initial_node_index]
# #     outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
# #     forward_walk = maximum_likelihood_walk(graph, [OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
# #     reverse_walk = maximum_likelihood_walk(graph, [OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
# #     reversed_reverse_walk = reverse!(
# #         [
# #             OrientedKmer(index = oriented_kmer.index, orientation = oriented_kmer.orientation)
# #             for oriented_kmer in reverse_walk[2:end]
# #         ]
# #         )
# #     full_path = [reversed_reverse_walk..., forward_walk...]
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function maximum_likelihood_walk(graph, path::Vector{OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
# #     done = false
# #     while !done
# #         maximum_path_likelihood = 0.0
# #         maximum_likelihood_path = Vector{OrientedKmer}()
# #         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
# #             this_path = [last(path).index, neighbor]
# #             this_oriented_path, this_path_likelihood = 
# #                 assess_path(this_path,
# #                     graph.kmers,
# #                     graph.counts,
# #                     last(path).orientation,
# #                     outgoing_edge_probabilities,
# #                     incoming_edge_probabilities)
# #             if this_path_likelihood > maximum_path_likelihood
# #                 maximum_path_likelihood = this_path_likelihood
# #                 maximum_likelihood_path = this_oriented_path
# #             end
# #         end
# #         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
# #             done = true
# #         else
# #             append!(path, maximum_likelihood_path[2:end])
# #         end
# #     end
# #     return path
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function take_a_walk(graph, connected_component)
# #     max_count = maximum(graph.counts[connected_component])
# #     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
# #     initial_node_index = rand(max_count_indices)
# #     initial_node = connected_component[initial_node_index]
# #     outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
    
# #     # walk forwards from the initial starting node
# #     forward_walk = take_a_walk(graph, [OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
# #     # walk backwards from the initial starting node
# #     reverse_walk = take_a_walk(graph, [OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
# #     # we need to reverse everything to re-orient against the forward walk
# #     reverse_walk = reverse_oriented_path(reverse_walk)
    
# #     # also need to drop the last node, which is equivalent to the first node of the 
# #     @assert last(reverse_walk) == first(forward_walk)
# #     full_path = [reverse_walk[1:end-1]..., forward_walk...]
# # #     @show full_path
# # end


# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function take_a_walk(graph, path::Vector{OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
# #     done = false
# #     while !done
# #         maximum_path_likelihood = 0.0
# #         maximum_likelihood_path = Vector{OrientedKmer}()
# #         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
# #             this_path = [last(path).index, neighbor]
# #             this_oriented_path, this_path_likelihood = 
# #                 assess_path(this_path,
# #                     graph.kmers,
# #                     graph.counts,
# #                     last(path).orientation,
# #                     outgoing_edge_probabilities,
# #                     incoming_edge_probabilities)
# #             if this_path_likelihood > maximum_path_likelihood
# #                 maximum_path_likelihood = this_path_likelihood
# #                 maximum_likelihood_path = this_oriented_path
# #             end
# #         end
# #         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
# #             done = true
# #         else
# #             append!(path, maximum_likelihood_path[2:end])
# #         end
# #     end
# #     return path
# # end


# # use default cost 1 for breadth first seach
# # use default cost 0 for depth first search
# function dijkstra_step!(graph, distances, arrival_paths, queue; default_cost = 1)
#     current_kmer, cost = DataStructures.dequeue_pair!(queue)
#     current_orientation = BioSequences.iscanonical(current_kmer)
#     current_canonical_kmer = BioSequences.canonical(current_kmer)
#     current_index = graph[current_canonical_kmer, :kmer]

#     present_neighbors = Vector{typeof(current_kmer)}()
#     for neighbor in BioSequences.neighbors(current_kmer)
#         try
#             graph[BioSequences.canonical(neighbor), :kmer]
#             push!(present_neighbors, neighbor)
#         catch
#             continue
#         end
#     end
    
#     true_neighbors =
#         filter(neighbor -> 
#             Graphs.has_edge(
#                 graph,
#                 graph[BioSequences.canonical(current_kmer), :kmer],
#                 graph[BioSequences.canonical(neighbor), :kmer]), 
#             present_neighbors)
    
#     if !isempty(true_neighbors)
#         canonical_true_neighbors = BioSequences.canonical.(true_neighbors)
#         true_neighbor_is_canonical = BioSequences.iscanonical.(true_neighbors)
#         neighbor_indices = map(canonical_neighbor -> graph[canonical_neighbor, :kmer], canonical_true_neighbors)
#         neighbor_weights = map(neighbor_i -> MetaGraphs.get_prop(graph, current_index, neighbor_i, :weight), neighbor_indices)


#         # inverse probability!!
#         neighbor_costs = 1 .- (neighbor_weights ./ sum(neighbor_weights))
#         neighbor_costs .+= default_cost
#         # afterwards, 100% likelihood edges cost 1
#         # 50% likelihood costs 1.5
#         # 0% likelihood is impossible and not included

#         for (neighbor, cost) in zip(true_neighbors, neighbor_costs)    
#             candidate_path = vcat(arrival_paths[current_kmer], current_kmer)
#             candidate_path_cost = distances[current_kmer] + cost  
#             if distances[neighbor] > candidate_path_cost
#                 arrival_paths[neighbor] = candidate_path
#                 DataStructures.enqueue!(queue, neighbor => cost)
#                 distances[neighbor] = candidate_path_cost
#             end  
#         end
#     end
#     return distances, arrival_paths, queue
# end

# function evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)

#     lowest_cost = Inf
#     optimal_path = Int[]

#     for hit in hits
#         reverse_complement_hit = BioSequences.reverse_complement(hit)
#         total_cost = forward_distances[hit] + reverse_distances[reverse_complement_hit]
#         if total_cost < lowest_cost    
#             forward_path = vcat(forward_arrival_paths[hit], hit)   
#             reverse_path = vcat(reverse_arrival_paths[reverse_complement_hit], reverse_complement_hit)
#             reverse_path = reverse!(BioSequences.reverse_complement.(reverse_path))
#             full_path = vcat(forward_path, reverse_path[2:end])

#             lowest_cost = total_cost
#             optimal_path = full_path
#         end
#     end

#     return lowest_cost, optimal_path
# end

# function bidirectional_dijkstra(graph, a::T, b::T) where T <: Kmers.Kmer
#     forward_distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     forward_distances[a] = 0

#     forward_arrival_paths = Dict{T, Vector{T}}()
#     forward_arrival_paths[a] = Int[]

#     forward_queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(forward_queue, a, 0)

#     reverse_b = BioSequences.reverse_complement(b)
#     reverse_distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     reverse_distances[reverse_b] = 0

#     reverse_arrival_paths = Dict{T, Vector{T}}()
#     reverse_arrival_paths[reverse_b] = Int[]

#     reverse_queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(reverse_queue, reverse_b, 0)
    
#     hits = Set{T}()
    
#     while isempty(hits)
#         dijkstra_step!(graph, forward_distances, forward_arrival_paths, forward_queue)
#         dijkstra_step!(graph, reverse_distances, reverse_arrival_paths, reverse_queue)
#         hits = intersect(
#                     keys(forward_arrival_paths), 
#                     BioSequences.reverse_complement.(keys(reverse_arrival_paths)))
#     end
#     lowest_cost, optimal_path = evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)
#     return lowest_cost, optimal_path
# end

# # # simple point to point
# # function dijkstra(graph, a::T, b::T) where {T <: Kmers.Kmer}
# #     distances = DataStructures.DefaultDict{T, Float64}(Inf)
# #     distances[a] = 0

# #     arrival_paths = Dict{T, Vector{T}}()
# #     arrival_paths[a] = Int[]

# #     queue = DataStructures.PriorityQueue{T, Float64}()
# #     DataStructures.enqueue!(queue, a, 0)

# #     current_kmer, cost = a, Inf
# #     while current_kmer != b
# #         distances, arrival_paths, queue = dijkstra_step!(graph, distances, arrival_paths, queue)
# #     end
# #     shortest_path = vcat(arrival_paths[b], b)
# #     return shortest_path, distances[b]    
# # end

# function dijkstra(graph, a::T, targets::Set{T}; search_strategy::Union{Symbol, Missing}=missing) where {T <: Kmers.Kmer}
#     if ismissing(search_strategy)
#         # local breadth first search to find the single target
#         # a-b shortest path
#         if length(targets) == 1
#             search_strategy = :BFS
#         # a - no targets - go find the ends
#         # a - series of any target
#         else
#             search_strategy = :DFS
#         end
#     end 
#     if search_strategy == :BFS
#         default_cost = 1
#     elseif search_strategy == :DFS
#         default_cost = 0
#     else
#         error("please specify either :BFS (breadth-first search) or :DFS (depth-first search)")
#     end
    
#     distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     distances[a] = 0

#     arrival_paths = Dict{T, Vector{T}}()
#     arrival_paths[a] = Int[]

#     queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(queue, a, 0)
    
#     while !isempty(queue) && !(first(DataStructures.peek(queue)) in targets)
#         distances, arrival_paths, queue = 
#             dijkstra_step!(graph, distances, arrival_paths, queue, default_cost = default_cost)
#     end
#     if !isempty(queue)
#         discovered_destination = DataStructures.dequeue!(queue)
#         @assert discovered_destination in targets
#         shortest_path = vcat(arrival_paths[discovered_destination], discovered_destination)
#         return shortest_path, distances[discovered_destination]
#     else
#         longest_path = Int[]
#         distance = Inf
#         for (destination, arrival_path) in arrival_paths
#             this_distance = distances[destination]
#             this_path = vcat(arrival_path, destination)
#             if (length(this_path) > length(longest_path)) && (this_distance <= distance)
#                 longest_path = this_path
#                 distance = this_distance
#             end
            
#         end
#         return longest_path, distance
#     end
# end

# function update_remaining_targets(current_walk::AbstractVector{T}, remaining_targets::AbstractSet{T}) where T <: Kmers.Kmer
#     # assess whether targets have been hit in the canonical space
#     remaining_targets = setdiff(BioSequences.canonical.(remaining_targets), BioSequences.canonical.(current_walk))
#     # blow back out into forward and reverse_complement space
#     remaining_targets = Set{T}(vcat(remaining_targets, BioSequences.reverse_complement.(remaining_targets)))
#     return remaining_targets
# end

# function assess_downstream_weight(graph, kmer)
#     # here we look to see if walking forward or backward from the initial node gets us to heavier weight options
#     score = 0
#     for neighbor in BioSequences.neighbors(kmer)
#         try
#             score += MetaGraphs.get_prop(graph, graph[BioSequences.canonical(neighbor), :kmer], :count)
#         catch
#             continue
#         end
#     end
#     return score
# end

# # vertices should either be entire graph (by default) or a connected component
# # if people want to work on just the connected component, let them induce a subgraph
# function find_graph_core(graph; seed=rand(Int))
    
#     Random.seed!(seed)
    
#     T = typeof(MetaGraphs.get_prop(graph, 1, :kmer))
    
#     # targets = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
#     # sortperm
    
#     # targets = [
#     #         MetaGraphs.get_prop(graph, v_index, :count) for v_index in 
#     #         sortperm([MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)], rev=true)[1:Int(ceil(Graphs.nv(graph)/10))]]
    
#     targets = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph) if MetaGraphs.get_prop(graph, v, :count) > 1]
#             # sortperm([MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)], rev=true)[1:Int(ceil(Graphs.nv(graph)/10))]]
    
#     @show length(targets)
#     starting_kmer = first(targets)
#     max_degree = 0
#     for node in targets
#         node_degree = Graphs.degree(graph, graph[node, :kmer])
#         if node_degree > max_degree
#             max_degree = node_degree
#             starting_kmer = node
#         end
#     end
        
#     current_walk = [starting_kmer]
#     prior_walk_length = length(current_walk)
#     remaining_targets = update_remaining_targets(current_walk, Set(targets))
#     done = isempty(remaining_targets)
    
#     while !done
#         # here we look to see if walking forward or backward from the current ends gets us to heavier weight options
#         # we want to prioritize walks toward higher coverage nodes
#         forward_score = assess_downstream_weight(graph, last(current_walk))
#         reverse_score = assess_downstream_weight(graph, BioSequences.reverse_complement(first(current_walk)))
#         if reverse_score > forward_score
#             current_walk = reverse(BioSequences.reverse_complement.(current_walk))
#         end
        
#         forward_source = last(current_walk)
#         forward_walk, forward_distance = dijkstra(graph, forward_source, remaining_targets, search_strategy=:DFS)
#         current_walk = vcat(current_walk, forward_walk[2:end])
#         remaining_targets = update_remaining_targets(current_walk, remaining_targets)
#         if isempty(remaining_targets)
#             done = true
#         else
#             reverse_source = BioSequences.reverse_complement(first(current_walk))
#             reverse_walk, reverse_distance = dijkstra(graph, reverse_source, remaining_targets, search_strategy=:DFS)
#             current_walk = vcat(reverse(BioSequences.reverse_complement.(reverse_walk))[1:end-1], current_walk)
#             remaining_targets = update_remaining_targets(current_walk, remaining_targets)
#         end
#         @show length(current_walk)
#         failed_this_expansion = length(current_walk) == prior_walk_length
#         prior_walk_length = length(current_walk)
#         if isempty(remaining_targets)
#             done = true
#         elseif failed_this_expansion
#             done = true
#         end
#     end
#     return current_walk
# end

# Save the distance matrix to a JLD2 file
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Saves a matrix to a JLD2 file format.

# Arguments
- `matrix`: The matrix to be saved
- `filename`: String path where the file should be saved

# Returns
- The filename string that was used to save the matrix
"""
function save_matrix_jld2(;matrix, filename)
    if !isfile(filename) || (filesize(filename) == 0)
        JLD2.@save filename matrix
    else
        @warn "$(filename) already exists and is non-empty, skipping..."
    end
    return filename
end

# Load the distance matrix from a JLD2 file
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Loads a matrix from a JLD2 file.

# Arguments
- `filename::String`: Path to the JLD2 file containing the matrix under the key "matrix"

# Returns
- `Matrix`: The loaded matrix data
"""
function load_matrix_jld2(filename)
    return JLD2.load(filename, "matrix")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load data stored in a JLD2 file format.

# Arguments
- `filename::String`: Path to the JLD2 file to load

# Returns
- `Dict`: Dictionary containing the loaded data structures
"""
function load_jld2(filename)
    return JLD2.load(filename)
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

Load a graph structure from a serialized file.

# Arguments
- `file::AbstractString`: Path to the file containing the serialized graph data

# Returns
- The deserialized graph object
"""
function load_graph(file)
    return FileIO.load(file)["graph"]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse TransTerm terminator prediction output into a structured DataFrame.

Takes a TransTerm output file path and returns a DataFrame containing parsed terminator predictions.
Each row represents one predicted terminator with the following columns:

- `chromosome`: Identifier of the sequence being analyzed
- `term_id`: Unique terminator identifier (e.g. "TERM 19")
- `start`: Start position of the terminator
- `stop`: End position of the terminator
- `strand`: Strand orientation ("+" or "-")
- `location`: Context type, where:
    * G/g = in gene interior (≥50bp from ends)
    * F/f = between two +strand genes
    * R/r = between two -strand genes
    * T = between ends of +strand and -strand genes
    * H = between starts of +strand and -strand genes
    * N = none of the above
    Lowercase indicates opposite strand from region
- `confidence`: Overall confidence score (0-100)
- `hairpin_score`: Hairpin structure score
- `tail_score`: Tail sequence score  
- `notes`: Additional annotations (e.g. "bidir")

# Arguments
- `transterm_output::AbstractString`: Path to TransTerm output file

# Returns
- `DataFrame`: Parsed terminator predictions with columns as described above

See TransTerm HP documentation for details on scoring and location codes.
"""
function parse_transterm_output(transterm_output)
    
   #     3. FORMAT OF THE TRANSTERM OUTPUT

#     The organism's genes are listed sorted by their end coordinate and terminators
#     are output between them. A terminator entry looks like this:

#         TERM 19  15310 - 15327  -      F     99      -12.7 -4.0 |bidir
#         (name)   (start - end)  (sense)(loc) (conf) (hp) (tail) (notes)

#     where 'conf' is the overall confidence score, 'hp' is the hairpin score, and
#     'tail' is the tail score. 'Conf' (which ranges from 0 to 100) is what you
#     probably want to use to assess the quality of a terminator. Higher is better.
#     The confidence, hp score, and tail scores are described in the paper cited
#     above.  'Loc' gives type of region the terminator is in:

#         'G' = in the interior of a gene (at least 50bp from an end),
#         'F' = between two +strand genes,
#         'R' = between two -strand genes,
#         'T' = between the ends of a +strand gene and a -strand gene,
#         'H' = between the starts of a +strand gene and a -strand gene,
#         'N' = none of the above (for the start and end of the DNA)

#     Because of how overlapping genes are handled, these designations are not
#     exclusive. 'G', 'F', or 'R' can also be given in lowercase, indicating that
#     the terminator is on the opposite strand as the region.  Unless the
#     --all-context option is given, only candidate terminators that appear to be in
#     an appropriate genome context (e.g. T, F, R) are output. 

#     Following the TERM line is the sequence of the hairpin and the 5' and 3'
#     tails, always written 5' to 3'.
    
    transterm_table = DataFrames.DataFrame()
    chromosome = ""
    for line in Iterators.filter(x -> occursin(r"^\s*(SEQUENCE|TERM)", x), eachline(transterm_output))
        line = strip(line)
        if occursin(r"^SEQUENCE", line)
            chromosome = split(line)[2]
        else
            transterm_regex = r"(TERM \d+)\s+(\d+) - (\d+)\s+(\S)\s+(\w+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\|(.*)"
            term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes = match(transterm_regex, line).captures
            notes = strip(notes)
            row = (;chromosome, term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes)
            push!(transterm_table, row)
        end
    end
    return transterm_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert TransTerm terminator predictions output to GFF3 format.

Parses TransTerm output and generates a standardized GFF3 file with the following transformations:
- Sets source field to "transterm"
- Sets feature type to "terminator"  
- Converts terminator IDs to GFF attributes
- Renames fields to match GFF3 spec

# Arguments
- `transterm_output::String`: Path to the TransTerm output file

# Returns
- `String`: Path to the generated GFF3 file (original filename with .gff extension)
"""
function transterm_output_to_gff(transterm_output)
    transterm_table = parse_transterm_output(transterm_output)
    transterm_table[!, "source"] .= "transterm"
    transterm_table[!, "type"] .= "terminator"
    transterm_table[!, "phase"] .= "."
    transterm_table[!, "attributes"] = map(x -> "label=" * replace(x, " " => "_"), transterm_table[!, "term_id"])
    DataFrames.rename!(transterm_table,
        ["chromosome" => "#seqid",
            "stop" => "end",
            "confidence" => "score",
        ]
    )
    transterm_table = transterm_table[!, ["#seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    transterm_gff = write_gff(gff=transterm_table, outfile=transterm_output * ".gff")
    uCSV.write(transterm_gff, transterm_table, delim='\t')
    return transterm_gff
end

# TODO: switch to using GenomicAnnotations if GFF3 package isn't updated
# PR -> https://github.com/BioJulia/GFF3.jl/pull/12
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site

Retrieve GFF3 formatted genomic feature data from NCBI or direct FTP source.

# Arguments
- `db::String`: NCBI database to query ("nuccore" for DNA or "protein" for protein sequences)
- `accession::String`: NCBI accession number
- `ftp::String`: Direct FTP URL to GFF3 file (typically gzipped)

# Returns
- `IO`: IOBuffer containing uncompressed GFF3 data
"""
function get_gff(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=gff3&id=$(accession)"
        return IOBuffer(HTTP.get(url).body)
    elseif !isempty(ftp)
        return CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body))
    else
        @error "invalid call"
    end
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert FASTA sequence and GFF annotation files to GenBank format using EMBOSS seqret.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequence data
- `gff::String`: Path to input GFF file containing genomic features
- `genbank::String`: Path for output GenBank file

# Details
Requires EMBOSS toolkit (installed via Bioconda). The function will:
1. Create necessary output directories
2. Run seqret to combine sequence and features
3. Generate a GenBank format file at the specified location
"""
function fasta_and_gff_to_genbank(;fasta, gff, genbank)
    add_bioconda_env("emboss")
    # https://www.insdc.org/submitting-standards/feature-table/
    genbank_directory = dirname(genbank)
    genbank_basename = basename(genbank)
    genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
    # genbank_extension = "genbank"
    # https://bioinformatics.stackexchange.com/a/11140
    # https://www.biostars.org/p/72220/#72272
    # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
#     -osname
    # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
    run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
    # return genbank
end

# function gff_to_genbank(gff, genbank)
#     # https://www.insdc.org/submitting-standards/feature-table/
#     genbank_directory = dirname(genbank)
#     genbank_basename = basename(genbank)
#     genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
#     # genbank_extension = "genbank"
#     # https://bioinformatics.stackexchange.com/a/11140
#     # https://www.biostars.org/p/72220/#72272
#     # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
# #     -osname
#     # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
#     run(`seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
#     # return genbank
# end

# function genbank_to_fasta_and_gff(genbank_file)

# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Opens and parses a GenBank format file containing genomic sequence and annotation data.

# Arguments
- `genbank_file::AbstractString`: Path to the GenBank (.gb or .gbk) file

# Returns
- `Vector{GenomicAnnotations.Chromosome}`: Vector containing parsed chromosome data
"""
function open_genbank(genbank_file)
    return GenomicAnnotations.readgbk(genbank_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a VirSorter score TSV file and return a DataFrame.

# Arguments
- `virsorter_score_tsv::String`: The file path to the VirSorter score TSV file.

# Returns
- `DataFrame`: A DataFrame containing the parsed data from the TSV file. If the file is empty, returns a DataFrame with the appropriate headers but no data.
"""
function parse_virsorter_score_tsv(virsorter_score_tsv)
    data, header = uCSV.read(virsorter_score_tsv, delim='\t', header=1)
    if length(data) == 0
        data = [[] for i in 1:length(header)]
    end
    return DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse MMseqs2 tophit alignment output file into a structured DataFrame.

# Arguments
- `tophit_aln::AbstractString`: Path to tab-delimited MMseqs2 alignment output file

# Returns
DataFrame with columns:
- `query`: Query sequence/profile identifier
- `target`: Target sequence/profile identifier  
- `percent identity`: Sequence identity percentage
- `alignment length`: Length of alignment
- `number of mismatches`: Count of mismatched positions
- `number of gaps`: Count of gap openings
- `query start`: Start position in query sequence
- `query end`: End position in query sequence
- `target start`: Start position in target sequence
- `target end`: End position in target sequence
- `evalue`: E-value of alignment
- `bit score`: Bit score of alignment
"""
function parse_mmseqs_tophit_aln(tophit_aln)
    data, header = uCSV.read(tophit_aln, delim='\t')
    # (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.
    # query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
    header = [
        "query",
        "target",
        "percent identity",
        "alignment length",
        "number of mismatches",
        "number of gaps",
        "query start",
        "query end",
        "target start",
        "target end",
        "evalue",
        "bit score"
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse an MMseqs2 easy-taxonomy tophit report into a structured DataFrame.

# Arguments
- `tophit_report::String`: Path to the MMseqs2 easy-taxonomy tophit report file (tab-delimited)

# Returns
- `DataFrame`: A DataFrame with columns:
  - `target_id`: Target sequence identifier
  - `number of sequences aligning to target`: Count of aligned sequences
  - `unique coverage of target`: Ratio of uniqueAlignedResidues to targetLength
  - `Target coverage`: Ratio of alignedResidues to targetLength
  - `Average sequence identity`: Mean sequence identity
  - `taxon_id`: Taxonomic identifier
  - `taxon_rank`: Taxonomic rank
  - `taxon_name`: Species name and lineage
"""
function parse_mmseqs_easy_taxonomy_tophit_report(tophit_report)
    data, header = uCSV.read(tophit_report, delim='\t')
    # tophit_report
    # (1) Target identifier 
    # (2) Number of sequences aligning to target
    # (3) Unique coverage of target uniqueAlignedResidues / targetLength
    # (4) Target coverage alignedResidues / targetLength
    # (5) Average sequence identity
    # (6) Taxonomical information identifier, species, lineage
    header = [
        "target_id",
        "number of sequences aligning to target",
        "unique coverage of target (uniqueAlignedResidues / targetLength)",
        "Target coverage (alignedResidues / targetLength)",
        "Average sequence identity",
        "taxon_id",
        "taxon_rank",
        "taxon_name"
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse the taxonomic Last Common Ancestor (LCA) TSV output from MMseqs2's easy-taxonomy workflow.

# Arguments
- `lca_tsv`: Path to the TSV file containing MMseqs2 taxonomy results

# Returns
DataFrame with columns:
- `contig_id`: Sequence identifier
- `taxon_id`: NCBI taxonomy identifier 
- `taxon_rank`: Taxonomic rank (e.g. species, genus)
- `taxon_name`: Scientific name
- `fragments_retained`: Number of fragments kept
- `fragments_taxonomically_assigned`: Number of fragments with taxonomy
- `fragments_in_agreement_with_assignment`: Fragments matching contig taxonomy
- `support -log(E-value)`: Statistical support score
"""
function parse_mmseqs_easy_taxonomy_lca_tsv(lca_tsv)
    data, header = uCSV.read(lca_tsv, delim='\t')
    # contig
    # (1) a single taxonomy numeric identifier
    # (2) a taxonomic rank column
    # (3) taxonomic name column
    # fragments retained
    # fragments taxonomically assigned
    # fragments in agreement with the contig label (i.e. same taxid or have it as an ancestor)
    # the support received -log(E-value)
    header = [
        "contig_id",
        "taxon_id",
        "taxon_rank",
        "taxon_name",
        "fragments_retained",
        "fragments_taxonomically_assigned",
        "fragments_in_agreement_with_assignment",
        "support -log(E-value)"
    ]
    return DataFrames.DataFrame(data, header)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function read_mmseqs_easy_search(mmseqs_file; top_hit_only=false)
#     mmseqs_results = DataFrames.DataFrame(uCSV.read(mmseqs_file, header=1, delim='\t')...)
#     if top_hit_only
#         gdf = DataFrames.groupby(mmseqs_results, "query")
#         for g in gdf
#             @assert issorted(g[!, "evalue"])
#         end
#         top_hits = DataFrames.combine(gdf, first)
#         mmseqs_results = top_hits
#     end
#     return mmseqs_results
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Read results from MMSeqs2 easy-search output file into a DataFrame.

# Arguments
- `mmseqs_file::String`: Path to the tab-delimited output file from MMSeqs2 easy-search

# Returns
- `DataFrame`: Contains search results with columns:
  - `query`: Query sequence identifier
  - `target`: Target sequence identifier
  - `seqIdentity`: Sequence identity (0.0-1.0)
  - `alnLen`: Alignment length
  - `mismatch`: Number of mismatches
  - `gapOpen`: Number of gap openings
  - `qStart`: Query start position
  - `qEnd`: Query end position
  - `tStart`: Target start position 
  - `tEnd`: Target end position
  - `evalue`: Expected value
  - `bits`: Bit score
"""
function read_mmseqs_easy_search(mmseqs_file)
    mmseqs_results = CSV.read(mmseqs_file, DataFrames.DataFrame, header=1, delim='\t')
    return mmseqs_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Update GFF annotations with protein descriptions from MMseqs2 search results.

# Arguments
- `gff_file::String`: Path to input GFF3 format file
- `mmseqs_file::String`: Path to MMseqs2 easy-search output file

# Returns
- `DataFrame`: Modified GFF table with updated attribute columns containing protein descriptions

# Details
Takes sequence matches from MMseqs2 and adds their descriptions as 'label' and 'product' 
attributes in the GFF file. Only considers top hits from MMseqs2 results. Preserves existing 
GFF attributes while prepending new annotations.
"""
function update_gff_with_mmseqs(gff_file, mmseqs_file)
    top_hits = read_mmseqs_easy_search(mmseqs_file, top_hit_only=true)

    id_to_product = Dict{String, String}()
    for row in DataFrames.eachrow(top_hits)
        id = Dict(a => b for (a, b) in split.(split(last(split(row["qheader"], " # ")), ';'), '='))["ID"]
        product = replace(row["theader"], " " => "__")
        id_to_product[id] = product
    end

    gff_table = Mycelia.read_gff(gff_file)
    for (i, row) in enumerate(DataFrames.eachrow(gff_table))
        id = Dict(a => b for (a,b) in split.(filter(!isempty, split(row["attributes"], ';')), '='))["ID"]
        product = get(id_to_product, id, "")
        gff_table[i, "attributes"] = "label=\"$(product)\";product=\"$(product)\";" * row["attributes"]
    end
    return gff_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Opens a GFF (General Feature Format) file for reading.

# Arguments
- `path::String`: Path to GFF file. Can be:
    - Local file path
    - HTTP/FTP URL (FTP URLs are automatically converted to HTTP)
    - Gzipped file (automatically decompressed)

# Returns
- `IO`: An IO stream ready to read the GFF content
"""
function open_gff(path::String)
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
    end
    return io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Reads a GFF (General Feature Format) file and parses it into a DataFrame.

# Arguments
- `gff::AbstractString`: Path to the GFF file

# Returns
- `DataFrame`: A DataFrame containing the parsed GFF data with standard columns:
  seqid, source, type, start, end, score, strand, phase, and attributes
"""
function read_gff(gff::AbstractString)
    return read_gff(open_gff(gff))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Read a GFF (General Feature Format) file into a DataFrame.

# Arguments
- `gff_io`: An IO stream containing GFF formatted data

# Returns
- `DataFrame`: A DataFrame with standard GFF columns:
  - seqid: sequence identifier
  - source: feature source
  - type: feature type
  - start: start position (1-based)
  - end: end position
  - score: numeric score
  - strand: strand (+, -, or .)
  - phase: phase (0, 1, 2 or .)
  - attributes: semicolon-separated key-value pairs
"""
function read_gff(gff_io)
    data, header = uCSV.read(gff_io, delim='\t', header=0, comment='#')
    header = [
        "#seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
    types=[
        String,
        String,
        String,
        Int,
        Int,
        Int,
        String,
        Int,
        String
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Takes a GFF (General Feature Format) DataFrame and expands the attributes column into separate columns.

# Arguments
- `gff_df::DataFrame`: A DataFrame containing GFF data with an 'attributes' column formatted as key-value pairs
  separated by semicolons (e.g., "ID=gene1;Name=BRCA1;Type=gene")

# Returns
- `DataFrame`: The input DataFrame with additional columns for each attribute key found in the 'attributes' column
"""
function split_gff_attributes_into_columns(gff_df)
    # 1. Extract unique keys from the attribute column
    all_keys = Set{String}()
    for row in DataFrames.eachrow(gff_df)
        attributes = split(row["attributes"], ';')
        for attribute in attributes
            if !isempty(attribute)
                key = split(attribute, '=')[1]
                push!(all_keys, key)
            end
        end
    end

  # 2. Create new columns for each key
    for key in all_keys
        gff_df[!, key] = Vector{Union{Missing, String}}(missing, size(gff_df, 1))
    end

    # 3. Populate the new columns with values
    for (i, row) in enumerate(DataFrames.eachrow(gff_df))
        attributes = split(row["attributes"], ';')
        for attribute in attributes
            if !isempty(attribute)
            key_value = split(attribute, '=')
                if length(key_value) == 2
                    key, value = key_value
                    gff_df[i, key] = value
                end
            end
        end
    end
    return gff_df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write GFF (General Feature Format) data to a tab-delimited file.

# Arguments
- `gff`: DataFrame/Table containing GFF formatted data
- `outfile`: String path where the output file should be written

# Returns
- `String`: Path to the written output file
"""
function write_gff(;gff, outfile)
    uCSV.write(outfile, gff, delim='\t')
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a Kraken taxonomic classification report into a structured DataFrame.

# Arguments
- `kraken_report::AbstractString`: Path to a tab-delimited Kraken report file

# Returns
- `DataFrame`: A DataFrame with the following columns:
  - `percentage_of_fragments_at_or_below_taxon`: Percentage of fragments covered
  - `number_of_fragments_at_or_below_taxon`: Count of fragments at/below taxon
  - `number_of_fragments_assigned_directly_to_taxon`: Direct fragment assignments
  - `rank`: Taxonomic rank
  - `ncbi_taxonid`: NCBI taxonomy identifier
  - `scientific_name`: Scientific name (whitespace-trimmed)

# Notes
- Scientific names are automatically stripped of leading/trailing whitespace
- Input file must be tab-delimited
"""
function read_kraken_report(kraken_report)
    kraken_report_header = [
        "percentage_of_fragments_at_or_below_taxon",
        "number_of_fragments_at_or_below_taxon",
        "number_of_fragments_assigned_directly_to_taxon",
        "rank",
        "ncbi_taxonid",
        "scientific_name"
    ]

    data, header = uCSV.read(kraken_report, delim='\t')
    kraken_report_table = DataFrames.DataFrame(data, kraken_report_header)
    kraken_report_table[!, "scientific_name"] = string.(strip.(kraken_report_table[!, "scientific_name"]))
    return kraken_report_table
end

# function diamond_line_to_named_tuple(diamond_line)
#     sline = split(line)
#     values_named_tuple = (
#         qseqid = sline[1],
#         sseqid = sline[2],
#         pident = parse(Float64, sline[3]),
#         length = parse(Int, sline[4]),
#         mismatch = parse(Int, sline[5]),
#         gapopen = parse(Int, sline[6]),
#         qlen = parse(Int, sline[7]),
#         qstart = parse(Int, sline[8]),
#         qend = parse(Int, sline[9]),
#         slen = parse(Int, sline[10]),
#         sstart = parse(Int, sline[11]),
#         send = parse(Int, sline[12]),
#         evalue = parse(Float64, sline[13]),
#         bitscore = parse(Float64, sline[14])
#         )
#     return values_named_tuple
# end

# function read_diamond_alignments_file(diamond_file)
#     column_names_to_types = [
#         "qseqid" => String,
#         "sseqid" => String,
#         "pident" => Float64,
#         "length" => Int,
#         "mismatch" => Int,
#         "gapopen" => Int,
#         "qlen" => Int,
#         "qstart" => Int,
#         "qend" => Int,
#         "slen" => Int,
#         "sstart" => Int,
#         "send" => Int,
#         "evalue" => Float64,
#         "bitscore" => Float64,
#     ]
#     types = Dict(i => t for (i, t) in enumerate(last.(column_names_to_types)))
    
#     data, header = uCSV.read(diamond_file, header=0, delim='\t', types = types)
#     header = first.(column_names_to_types)    
    
#     # data, header = uCSV.read(diamond_file, header=1, delim='\t', types = types)
#     # @assert header == first.(column_names_to_types)
    
#     table = DataFrames.DataFrame(data, header)
#     return table
# end

# function add_header_to_diamond_file(infile, outfile=replace(infile, ".tsv" => ".with-header.tsv"))
#     column_names = [
#         "qseqid",
#         "sseqid",
#         "pident",
#         "length",
#         "mismatch",
#         "gapopen",
#         "qlen",
#         "qstart",
#         "qend",
#         "slen",
#         "sstart",
#         "send",
#         "evalue",
#         "bitscore"
#     ]
#     # dangerous but fast
#     # try
#     #     inserted_text = join(columns_names, '\t') * '\n'
#     #     sed_cmd = "1s/^/$(inserted_text)/"
#     #     full_cmd = `sed -i $sed_cmd $infile`
#     # catch
#     open(outfile, "w") do io
#         println(io, join(column_names, "\t"))
#         for line in eachline(infile)
#             println(io, line)
#         end
#     end
#     return outfile
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse contig coverage statistics from a Qualimap BAM QC report file.

# Arguments
- `qualimap_report_txt::String`: Path to Qualimap bamqc report text file

# Returns
- `DataFrame`: Coverage statistics with columns:
  - `Contig`: Contig identifier
  - `Length`: Contig length in bases
  - `Mapped bases`: Number of bases mapped to contig
  - `Mean coverage`: Average coverage depth
  - `Standard Deviation`: Standard deviation of coverage
  - `% Mapped bases`: Percentage of total mapped bases on this contig

# Supported Assemblers
Handles output from both SPAdes and MEGAHIT assemblers:
- SPAdes format: NODE_X_length_Y_cov_Z
- MEGAHIT format: kXX_Y 

Parse the contig coverage information from qualimap bamqc text report, which looks like the following:

```
# this is spades
>>>>>>> Coverage per contig

	NODE_1_length_107478_cov_9.051896	107478	21606903	201.0355886786133	60.39424208607496
	NODE_2_length_5444_cov_1.351945	5444	153263	28.152645113886848	5.954250612823136
	NODE_3_length_1062_cov_0.154390	1062	4294	4.043314500941619	1.6655384692688975
	NODE_4_length_776_cov_0.191489	776	3210	4.13659793814433	2.252009588980858

# below is megahit
>>>>>>> Coverage per contig

	k79_175	235	3862	16.43404255319149	8.437436249612457
	k79_89	303	3803	12.551155115511552	5.709975376279777
	k79_262	394	6671	16.931472081218274	7.579217802849293
	k79_90	379	1539	4.060686015831134	1.2929729111266581
	k79_91	211	3749	17.767772511848342	11.899185693011933
	k79_0	2042	90867	44.49902056807052	18.356525483516613
```

To make this more robust, consider reading in the names of the contigs from the assembled fasta
"""
function parse_qualimap_contig_coverage(qualimap_report_txt)
    coverage_line_regex = r"\t.*?\t\d+\t\d+\t[\d\.]+\t[\d\.]+$"
    lines = filter(x -> occursin(coverage_line_regex, x), readlines("$(qualimap_report_txt)"))
    io = IOBuffer(join(map(x -> join(split(x, '\t')[2:end], '\t'), lines), '\n'))
    header = ["Contig", "Length", "Mapped bases", "Mean coverage", "Standard Deviation"]
    types = [String, Int, Int, Float64, Float64]
    data, _ = uCSV.read(io, delim='\t', types=types)
    qualimap_results = DataFrames.DataFrame(data, header)
    qualimap_results[!, "% Mapped bases"] = qualimap_results[!, "Mapped bases"] ./ sum(qualimap_results[!, "Mapped bases"]) .* 100
    return qualimap_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Imports results of fastani

Reads and processes FastANI output results from a tab-delimited file.

# Arguments
- `path::String`: Path to the FastANI output file

# Returns
DataFrame with columns:
- `query`: Original query filepath
- `query_identifier`: Extracted filename without extension
- `reference`: Original reference filepath
- `reference_identifier`: Extracted filename without extension
- `%_identity`: ANI percentage identity
- `fragments_mapped`: Number of fragments mapped
- `total_query_fragments`: Total number of query fragments

# Notes
- Expects tab-delimited input file from FastANI
- Automatically strips .fasta, .fna, or .fa extensions from filenames
- Column order is preserved as listed above
"""
function read_fastani(path::String)
    data, header = uCSV.read(path, delim='\t', typedetectrows=100)
    header = [
        "query",
        "reference",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]
    ani_table = DataFrames.DataFrame(data, header)
    ani_table[!, "query_identifier"] = replace.(basename.(ani_table[!, "query"]), r"\.(fasta|fna|fa)$" => "")
    ani_table[!, "reference_identifier"] = replace.(basename.(ani_table[!, "reference"]), r"\.(fasta|fna|fa)$" => "")
    columns = [
        "query",
        "query_identifier",
        "reference",
        "reference_identifier",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]
    ani_table = ani_table[!, columns]    
    return ani_table
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Imports results of Diamond (or blast) in outfmt 6 as a DataFrame
# """
# function read_diamond(path::String)
#   diamond_colnames = [ "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
  
#   diamond_coltypes = Dict(
#      1 => String, 
#      2 => String, 
#      3 => Float64,
#      4 => Int,
#      5 => Int,
#      6 => Int,
#      7 => Int,
#      8 => Int,
#      9 => Int,
#      10 => Int,
#      11 => String,
#      12 => String
#   )
#     return DataFrames.DataFrame(uCSV.read(open(path), delim = '\t', header = diamond_colnames, types = diamond_coltypes)...)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Expects output type 7 from BLAST, default output type 6 doesn't have the header comments and won't auto-parse

Parse a BLAST output file into a structured DataFrame.

# Arguments
- `blast_report::AbstractString`: Path to a BLAST output file in format 7 (tabular with comments)

# Returns
- `DataFrame`: Table containing BLAST results with columns matching the header fields.
  Returns empty DataFrame if no hits found.

# Details
- Requires BLAST output format 7 (`-outfmt 7`), which includes header comments
- Handles missing values (encoded as "N/A") automatically
- Infers column types based on BLAST field names
- Supports standard BLAST tabular fields including sequence IDs, scores, alignments and taxonomic information
"""
function parse_blast_report(blast_report)
    # example header line 
    # "# Fields: query id, subject id, subject acc., subject acc.ver, subject title, query length, subject length, q. start, q. end, s. start, s. end, evalue, bit score, score, alignment length, % identity, identical, mismatches, subject tax id"
    header_lines = collect(Iterators.filter(x -> occursin(r"# Fields:", x), eachline(blast_report)))
    if isempty(header_lines)
        @info "not hits found, returning empty table"
        return DataFrames.DataFrame()
    end
    header_line = first(header_lines)
    header = split(last(split(header_line, ": ")), ", ")
    blast_col_types = Dict(
        "query id" => String,
        "query title" => String,
        "subject id" => String,
        "subject gi" => String,
        "subject acc." => String,
        "subject acc.ver" => String,
        "subject title" => String,
        "query length" => Int,
        "subject length" => Int,
        "q. start" => Int,
        "q. end" => Int,
        "s. start" => Int,
        "s. end" => Int,
        "evalue" => Float64,
        "bit score" => Float64,
        "score" => Float64,
        "alignment length" => Int,
        "% identity" => Float64,
        "identical" => Int,
        "mismatches" => Int,
        "subject tax id" => Int,
        "subject sci name" => String,
        "subject com names" => String,
        "subject blast name" => String,
        "subject super kingdom" => String,
        "subject tax ids" => String,
        "subject sci names" => String,
        "subject com names" => String,
        "subject blast names" => String,
        "subject super kingdoms" => String,
        "subject title" => String,
        "subject titles" => String
    )
    data, _ = uCSV.read(
        blast_report,
        delim='\t',
        comment='#',
        skipmalformed=true,
        allowmissing=true,
        encodings=Dict("N/A" => missing),
        types=[blast_col_types[h] for h in header])
    return DataFrames.DataFrame(data, header, makeunique=true)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Parse a GFA file into a genome graph - need to finish implementation and assert contig normalization (i.e. is canonical) before using with my code
# """
# function parse_gfa(gfa)
    
#     # collect(GraphicalFragmentAssembly.Reader(open(primary_contig_gfa)))
    
#     gfa_record_types = Dict(
#         '#' => "Comment",
#         'H' => "Header",
#         'S' => "Segment",
#         'L' => "Link",
#         'J' => "Jump",
#         'C' => "Containment",
#         'P' => "Path",
#         'W' => "Walk"
#     )

#     gfa_graph = MetaGraphs.MetaDiGraph()
#     MetaGraphs.set_prop!(gfa_graph, :paths, Dict{String, Any}())
#     for line in eachline(gfa)
#         record_type = gfa_record_types[line[1]]
#         if record_type == "Header"
#             # metadata
#             sline = split(line)
#             # add me later
#         elseif record_type == "Comment"
#             # metadata
#             # add me later
#         elseif record_type == "Segment"
#             # node
#             record_type, record_name, sequence = split(line, '\t')
#             Graphs.add_vertex!(gfa_graph)
#             node_index = Graphs.nv(gfa_graph)
#             MetaGraphs.set_prop!(gfa_graph, node_index, :identifier, record_name)
#             MetaGraphs.set_indexing_prop!(gfa_graph, :identifier)
#             MetaGraphs.set_prop!(gfa_graph, node_index, :sequence, sequence)
#         elseif record_type == "Link"
#             record_type, source_identifier, source_orientation, destination_identifier, destination_orientation, overlap_CIGAR = split(line, '\t')
#             source_index = gfa_graph[source_identifier, :identifier]
#             destination_index = gfa_graph[destination_identifier, :identifier]
#             edge = Graphs.Edge(source_index, destination_index)
#             Graphs.add_edge!(gfa_graph, edge)
#             MetaGraphs.set_prop!(gfa_graph, edge, :source_identifier, source_identifier)
#             MetaGraphs.set_prop!(gfa_graph, edge, :source_orientation, source_orientation)
#             MetaGraphs.set_prop!(gfa_graph, edge, :destination_identifier, destination_identifier)
#             MetaGraphs.set_prop!(gfa_graph, edge, :destination_orientation, destination_orientation)
#             MetaGraphs.set_prop!(gfa_graph, edge, :overlap_CIGAR, overlap_CIGAR)
#         elseif record_type == "Path"
#             record_type, path_identifier, segments, overlaps = split(line, '\t')
#             gfa_graph.gprops[:paths][path_identifier] = Dict("segments" => segments, "overlaps" => overlaps)
#         else
#             @warn "GFA line type $(record_type) not currently handled by the import - please add"
#         end
#     end
#     return gfa_graph
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Mycelia graph to GFA (Graphical Fragment Assembly) format.

Writes a graph to GFA format, including:
- Header (H) line with GFA version
- Segment (S) lines for each vertex with sequence and depth
- Link (L) lines for edges with overlap size and orientations

# Arguments
- `graph`: MetaGraph containing sequence vertices and their relationships
- `outfile`: Path where the GFA file should be written

# Returns
- Path to the written GFA file
"""
function graph_to_gfa(;graph, outfile)
    # kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
    # add fastq here too???
    # record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
    # edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0")
        for vertex in Graphs.vertices(graph)
            if haskey(graph.vprops[vertex], :kmer)
                sequence = graph.vprops[vertex][:kmer]
            else
                sequence = graph.vprops[vertex][:sequence]
            end
            depth = graph.vprops[vertex][:count]
            # depth = length(graph.vprops[vertex][:evidence])
            # total_count = 0
            # vertex_outneighbors = Graphs.outneighbors(graph, vertex)
            # for connected_record in vertex_outneighbors
            #     edge = Graphs.Edge(vertex, connected_record)
            #     total_count += MetaGraphs.get_prop(graph, edge, :count)
            # end

            fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
            line = join(fields, '\t')
            println(io, line)
        end
        overlap = MetaGraphs.get_prop(graph, :k) - 1
        for edge in collect(Graphs.edges(graph))
            # @show edge
            # @show MetaGraphs.props(graph, edge)
            unique_orientations = unique(observation.orientation for observation in MetaGraphs.get_prop(graph, edge, :observations))
            # @show unique_orientations
            source_kmer = edge.src
            destination_kmer = edge.dst
            for unique_orientation in unique_orientations
                # source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
                link = ["L",
                            edge.src,
                            first(unique_orientation) ? '+' : '-',
                            edge.dst,
                            last(unique_orientation) ? '+' : '-',
                            "$(overlap)M"]
                line = join(link, '\t')
                println(io, line)
            end
                
            # source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
            # link = ["L",
            #             edgemer_edge.src,
            #             BioSequences.iscanonical(source_kmer) ? '+' : '-',
            #             edgemer_edge.dst,
            #             BioSequences.iscanonical(dest_kmer) ? '+' : '-',
            #             "$(overlap)M"]
            # line = join(link, '\t')
            # println(io, line)
        end
    end
    return outfile
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function graph_to_gfa(graph, kmer_size, outfile="$(kmer_size).gfa")
#     kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
#     # add fastq here too???
#     record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
#     edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0")
#         for vertex in kmer_vertices
#             if haskey(graph.vprops[vertex], :kmer)
#                 sequence = graph.vprops[vertex][:kmer]
#             else
#                 sequence = graph.vprops[vertex][:sequence]
#             end
#             # depth = graph.vprops[vertex][:count]
#             # depth = length(graph.vprops[vertex][:evidence])
#             total_count = 0
#             vertex_outneighbors = Graphs.outneighbors(graph, vertex)
#             for connected_record in intersect(vertex_outneighbors, record_vertices)
#                 edge = Graphs.Edge(vertex, connected_record)
#                 total_count += MetaGraphs.get_prop(graph, edge, :count)
#             end

#             fields = ["S", "$vertex", sequence, "DP:f:$(total_count)"]
#             line = join(fields, '\t')
#             println(io, line)
#         end
#         overlap = kmer_size - 1
#         for edgemer_edge in edgemer_edges
#             source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
#             link = ["L",
#                         edgemer_edge.src,
#                         BioSequences.iscanonical(source_kmer) ? '+' : '-',
#                         edgemer_edge.dst,
#                         BioSequences.iscanonical(dest_kmer) ? '+' : '-',
#                         "$(overlap)M"]
#             line = join(link, '\t')
#             println(io, line)
#         end
#     end
#     return
# end

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
- `outfile::AbstractString`: Path to the output FASTA file.  Will append ".gz" if `gzip` is true.
- `records::Vector{FASTX.FASTA.Record}`: A vector of FASTA records.
- `gzip::Bool=false`: Whether to compress the output with gzip.

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

Plots a histogram of kmer counts against # of kmers with those counts

Returns the plot object for adding additional layers and saving

Creates a scatter plot visualizing the k-mer frequency spectrum - the relationship
between k-mer frequencies and how many k-mers occur at each frequency.

# Arguments
- `counts::AbstractVector{<:Integer}`: Vector of k-mer counts/frequencies
- `log_scale::Union{Function,Nothing} = log2`: Function to apply logarithmic scaling to both axes.
  Set to `nothing` to use linear scaling.
- `kwargs...`: Additional keyword arguments passed to `StatsPlots.plot()`

# Returns
- `Plots.Plot`: A scatter plot object that can be further modified or saved

# Details
The x-axis shows k-mer frequencies (how many times each k-mer appears),
while the y-axis shows how many distinct k-mers appear at each frequency.
Both axes are log-scaled by default using log2.
"""
function plot_kmer_frequency_spectra(counts; log_scale = log2, kwargs...)
    kmer_counts_hist = StatsBase.countmap(c for c in counts)
    xs = collect(keys(kmer_counts_hist))
    ys = collect(values(kmer_counts_hist))
    if isa(log_scale, Function)
        xs = log_scale.(xs)
        ys = log_scale.(ys)
    end
    
    p = StatsPlots.plot(
        xs,
        ys,
        xlims = (0, maximum(xs) + max(1, ceil(0.1 * maximum(xs)))),
        ylims = (0, maximum(ys) + max(1, ceil(0.1 * maximum(ys)))),
        seriestype = :scatter,
        legend = false,
        xlabel = isa(log_scale, Function) ? "$(log_scale)(observed frequency)" : "observed frequency",
        ylabel = isa(log_scale, Function) ? "$(log_scale)(# of kmers)" : "observed frequency",
        ;kwargs...
    )
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert between different graph file formats.

# Arguments
- `args`: Dictionary with required keys:
    - `"in"`: Input filepath (supported: .jld2, .gfa, .neo4j)
    - `"out"`: Output filepath (supported: .jld2, .gfa, .neo4j)

# Details
Performs format conversion based on file extensions. For non-JLD2 to non-JLD2
conversions, uses JLD2 as an intermediate format.
"""
function convert(args)
    @show args
    in_type = missing
    out_type = missing
    if occursin(r"\.jld2$", args["in"])
        in_type = :jld2
    elseif occursin(r"\.gfa$", args["in"])
        in_type = :gfa
    elseif occursin(r"\.neo4j$", args["in"])
        in_type = :neo4j
    end

    if occursin(r"\.jld2$", args["out"])
        out_type = :jld2
    elseif occursin(r"\.gfa$", args["out"])
        out_type = :gfa
    elseif occursin(r"\.neo4j$", args["out"])
        out_type = :neo4j
    end
    
    if ismissing(in_type) || ismissing(out_type)
        error("unable to determine in and out types")
    end
    
    if (in_type == :jld2) && (out_type == :jld2)
        # done
    elseif (in_type == :jld2) && (out_type != :jld2)
        # convert out of jld2
        if out_type == :gfa
            loaded = FileIO.load(args["in"])
            @assert haskey(loaded, "graph")
            graph = loaded["graph"]
            graph_to_gfa(graph, args["out"])
        end
    elseif (in_type != :jld2) && (out_type == :jld2)
        # convert into jld2
    else
        # need to convert from input to jld2 as an intermediate first
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a visualization of a kmer graph where nodes represent kmers and their sizes reflect counts.

# Arguments
- `graph`: A MetaGraph where vertices have `:kmer` and `:count` properties

# Returns
- A Plots.jl plot object showing the graph visualization

# Details
- Node sizes are scaled based on kmer counts
- Plot dimensions scale logarithmically with number of vertices
- Each node is labeled with its kmer sequence
"""
function plot_graph(graph)
    
#     kmer_counts = MetaGraphs.get_prop(graph, :kmer_counts)
    kmers = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
    counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]
    scale = 150
    
    n = Graphs.nv(graph)
    p = GraphRecipes.graphplot(
        graph,
#         markersize = 1/log2(n),
        markersize = 1/2^2,
        size = (2 * scale * log(n), scale * log(n)),
        node_weights = counts,
        names = kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Saves the given graph to a file in JLD2 format.

# Arguments
- `graph::Graphs.AbstractGraph`: The graph to be saved.
- `outfile::String`: The name of the output file. If the file extension is not `.jld2`, it will be appended automatically.

# Returns
- `String`: The name of the output file with the `.jld2` extension.
"""
function save_graph(graph::Graphs.AbstractGraph, outfile::String)
    if !occursin(r"\.jld2$", outfile)
        outfile *= ".jld2"
    end
    FileIO.save(outfile, Dict("graph" => graph))
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Loads a graph object from a serialized file.

# Arguments
- `file::String`: Path to the file containing the serialized graph data. The file should have been created using `save_graph`.

# Returns
- The deserialized graph object stored under the "graph" key.
"""
function load_graph(file::String)
    return FileIO.load(file)["graph"]
end

# OLD FOR SIMPLE KMER GRAPHS
# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # Description

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# function graph_to_gfa(graph, outfile)
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0")
#         for vertex in Graphs.vertices(graph)
#             if haskey(graph.vprops[vertex], :kmer)
#                 sequence = graph.vprops[vertex][:kmer]
#             else
#                 sequence = graph.vprops[vertex][:sequence]
#             end
# #             depth = graph.vprops[vertex][:weight]
#             depth = graph.vprops[vertex][:count]
#             fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
#             line = join(fields, '\t')
#             println(io, line)
#         end
#         for edge in Graphs.edges(graph)
#             overlap = graph.gprops[:k] - 1
#             for o in graph.eprops[edge][:orientations]
# #                 if !(!o.source_orientation && !o.destination_orientation)
#                 link = ["L",
#                             edge.src,
#                             o.source_orientation ? '+' : '-',
#                             edge.dst,
#                             o.destination_orientation ? '+' : '-',
#                             "$(overlap)M"]
#                 line = join(link, '\t')
#                 println(io, line)
# #                 end
#             end
#         end
#     end
#     return outfile
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function my_plot(graph::KmerGraph)
#     graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))
#     filename = "/assets/images/$(graph_hash).svg"
#     p = plot_graph(graph)
#     Plots.savefig(p, dirname(pwd()) * filename)
#     display(p)
#     display("text/markdown", "![]($filename)")
# end

# function observed_kmer_frequencies(seq::BioSequences.BioSequence{A}, k::Int) where A<:BioSequences.Alphabet
#     kmer_count_dict = BioSequences.kmercounts(seq, k)
#     total_kmers = sum(values(kmer_count_dict))
#     return Dict(kmer => count / total_kmers for (kmer, count) in kmer_count_dict)
# end

# function expected_kmer_frequencies(seq::BioSequences.BioSequence{A}, k::Int) where A<:BioSequences.Alphabet
#     k_minus_1_mer_counts = Dict{BioSequences.BioSequence{A}, Int}()
#     for i in 1:(length(seq) - k + 1)
#         k_minus_1_mer = seq[i:(i + k - 2)]
#         k_minus_1_mer_counts[k_minus_1_mer] = get(k_minus_1_mer_counts, k_minus_1_mer, 0) + 1
#     end

#     nucleotide_counts = Dict{BioSequences.BioSequence{A}, Int}()
#     for nucleotide in seq
#         nucleotide_seq = BioSequences.BioSequence{A}([nucleotide])
#         nucleotide_counts[nucleotide_seq] = get(nucleotide_counts, nucleotide_seq, 0) + 1
#     end

#     total_k_minus_1_mers = sum(values(k_minus_1_mer_counts))
#     expected_frequencies = Dict{BioSequences.BioSequence{A}, Float64}()
#     for (k_minus_1_mer, count) in k_minus_1_mer_counts
#         last_nucleotide = k_minus_1_mer[end]
#         last_nucleotide_seq = BioSequences.BioSequence{A}([last_nucleotide])
#         expected_frequencies[k_minus_1_mer] = (count / total_k_minus_1_mers) * (nucleotide_counts[last_nucleotide_seq] / length(seq))
#     end

#     return expected_frequencies
# end

# function precompute_expected_frequencies(sequences::Vector{BioSequences.BioSequence{A}}, k::Int) where A<:BioSequences.Alphabet
#     expected_freqs = Vector{Dict{BioSequences.BioSequence{A}, Float64}}(undef, length(sequences))
#     for (i, seq) in enumerate(sequences)
#         expected_freqs[i] = expected_kmer_frequencies(seq, k)
#     end
#     return expected_freqs
# end

# function D2_star(obs_freq1::Dict{BioSequences.BioSequence{A}, Float64},
#                  obs_freq2::Dict{BioSequences.BioSequence{A}, Float64},
#                  exp_freq::Dict{BioSequences.BioSequence{A}, Float64}) where A<:BioSequences.Alphabet
#     numerator = sum((obs_freq1[kmer] - exp_freq[kmer]) * (obs_freq2[kmer] - exp_freq[kmer]) for kmer in keys(exp_freq))
#     denominator = sqrt(sum((obs_freq1[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq)) *
#                        sum((obs_freq2[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq)))
#     return numerator / denominator
# end

# function d2_star_normalized(obs_freq1::Dict{BioSequences.BioSequence{A}, Float64},
#                             obs_freq2::Dict{BioSequences.BioSequence{A}, Float64},
#                             exp_freq::Dict{BioSequences.BioSequence{A}, Float64}) where A<:BioSequences.Alphabet
#     numerator = sum((obs_freq1[kmer] - exp_freq[kmer]) * (obs_freq2[kmer] - exp_freq[kmer]) for kmer in keys(exp_freq))
#     var1 = sum((obs_freq1[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq))
#     var2 = sum((obs_freq2[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq))
#     denominator = sqrt(var1 * var2)
#     return numerator / denominator
# end

# function compute_observed_frequencies(sequences::Vector{BioSequences.BioSequence{A}}, k::Int) where A<:BioSequences.Alphabet
#     observed_freqs = Vector{Dict{BioSequences.BioSequence{A}, Float64}}(undef, length(sequences))
#     for (i, seq) in enumerate(sequences)
#         observed_freqs[i] = observed_kmer_frequencies(seq, k)
#     end
#     return observed_freqs
# end

# function compute_pairwise_statistics(observed_freqs::Vector{Dict{BioSequences.BioSequence{A}, Float64}},
#                                      expected_freqs::Vector{Dict{BioSequences.BioSequence{A}, Float64}}) where A<:BioSequences.Alphabet
#     n = length(observed_freqs)
#     d2_star_matrix = Matrix{Float64}(undef, n, n)
#     d2_star_norm_matrix = Matrix{Float64}(undef, n, n)
#     for i in 1:n
#         for j in i:n
#             d2_star_value = D2_star(observed_freqs[i], observed_freqs[j], expected_freqs[i])
#             d2_star_norm_value = d2_star_normalized(observed_freqs[i], observed_freqs[j], expected_freqs[i])
#             d2_star_matrix[i, j] = d2_star_value
#             d2_star_matrix[j, i] = d2_star_value
#             d2_star_norm_matrix[i, j] = d2_star_norm_value
#             d2_star_norm_matrix[j, i] = d2_star_norm_value
#         end
#     end
#     return d2_star_matrix, d2_star_norm_matrix
# end


# # Example sequences
# sequences = [
#     BioSequences.LongSequence{BioSymbols.DNAAlphabet{4}}("ACGTACGTGACG"),
#     BioSequences.LongSequence{BioSymbols.DNAAlphabet{4}}("TGCATGCATGCA"),
#     # Add more sequences as needed
# ]

# # k-mer length
# k = 3

# # Precompute expected k-mer frequencies
# expected_freqs = precompute_expected_frequencies(sequences, k)

# # Compute observed k-mer frequencies
# observed_freqs = compute_observed_frequencies(sequences, k)

# # Compute pairwise D2* and d2* statistics
# d2_star_matrix, d2_star_norm_matrix = compute_pairwise_statistics(observed_freqs, expected_freqs)

# println("D2* Matrix:")
# println(d2_star_matrix)

# println("Normalized d2* Matrix:")
# println(d2_star_norm_matrix)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a mapping from amino acids to representative DNA codons using the standard genetic code.

# Returns
- Dictionary mapping each amino acid (including stop codon `AA_Term`) to a valid DNA codon that encodes it
"""
function amino_acids_to_codons()
    amino_acid_to_codon_map = Dict(a => Kmers.DNACodon for a in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term]))
    for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
        amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
        amino_acid_to_codon_map[amino_acid] = codon
    end   
    return amino_acid_to_codon_map
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a mapping between DNA codons and their corresponding amino acids using the standard genetic code.

Returns a dictionary where:
- Keys are 3-letter DNA codons (e.g., "ATG")
- Values are the corresponding amino acids from BioSequences.jl
"""
function codons_to_amino_acids()
    codons = Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
    codon_to_amino_acid_map = Dict(codon => BioSequences.translate(BioSequences.LongDNA{2}(codon)))
    return codon_to_amino_acid_map
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyze codon usage frequencies from genes in a GenBank file.

# Arguments
- `genbank`: Path to GenBank format file containing genomic sequences and annotations
- `allow_all`: If true, initializes frequencies for all possible codons with count=1 (default: true)

# Returns
Nested dictionary mapping amino acids to their corresponding codon usage counts:
- Outer key: AminoAcid (including stop codon)
- Inner key: DNACodon
- Value: Count of codon occurrences

# Details
- Only processes genes marked as ':misc_feature' in the GenBank file
- Analyzes both forward and reverse complement sequences
- Determines coding strand based on presence of stop codons and start codons
- Skips ambiguous sequences that cannot be confidently oriented
"""
function genbank_to_codon_frequencies(genbank; allow_all=true)
    # create an initial codon frequency table, where we initialize all possible codons with equal probability
    # this way, if we don't see the amino acid in the observed proteins we're optimizing, we can still produce an codon profile
    codon_frequencies = Dict(a => Dict{Kmers.DNACodon, Int}() for a in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term]))
    if allow_all
        for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
            amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
            codon_frequencies[amino_acid][codon] = get(codon_frequencies[amino_acid], codon, 0) + 1
        end    
    end
    genome_genbank_data = GenomicAnnotations.readgbk(genbank)
    for chromosome in genome_genbank_data
        for gene in chromosome.genes
            gene_range = GenomicAnnotations.locus(gene).position
            gene_type = GenomicAnnotations.feature(gene)
            # is_terminator = occursin(r"^TERM", gene.label)
            if (length(gene_range) % 3 == 0) && (gene_type == :misc_feature)

                fw_dnaseq = GenomicAnnotations.sequence(gene)
                revcom_dnaseq = BioSequences.reverse_complement(fw_dnaseq)

                fw_aaseq = BioSequences.translate(fw_dnaseq)
                revcom_aaseq = BioSequences.translate(revcom_dnaseq)

                if last(fw_aaseq) == BioSequences.AA_Term
                    dnaseq = fw_dnaseq
                    aaseq = fw_aaseq
                elseif last(revcom_aaseq) == BioSequences.AA_Term
                    dnaseq = revcom_dnaseq
                    aaseq = revcom_aaseq
                elseif first(fw_aaseq) == BioSequences.AA_M
                    dnaseq = fw_dnaseq
                    aaseq = fw_aaseq
                elseif first(revcom_aaseq) == BioSequences.AA_M
                    dnaseq = revcom_dnaseq
                    aaseq = revcom_aaseq
                else
                    # @show "ambiguous"
                    continue
                end
                
                # for (mer, amino_acid) in zip(BioSequences.each(BioSequences.DNAMer{3}, dnaseq, 3), aaseq)
                for ((i, codon), amino_acid) in zip(Kmers.SpacedKmers{Kmers.DNACodon}(BioSequences.LongDNA{4}(dnaseq), 3), aaseq)
                    @assert amino_acid == first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
                    codon_frequencies[amino_acid][codon] = get(codon_frequencies[amino_acid], codon, 0) + 1
                end
            end
        end
    end
    return codon_frequencies
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalizes codon frequencies for each amino acid such that frequencies sum to 1.0.

# Arguments
- `codon_frequencies`: Nested dictionary mapping amino acids to their codon frequency distributions

# Returns
- Normalized codon frequencies where values for each amino acid sum to 1.0
"""
function normalize_codon_frequencies(codon_frequencies)
    normalized_codon_frequencies = Dict{BioSymbols.AminoAcid, Dict{Kmers.DNACodon, Float64}}()
    for (amino_acid, amino_acid_codon_frequencies) in codon_frequencies
        total_count = sum(values(amino_acid_codon_frequencies))
        normalized_codon_frequencies[amino_acid] = Dict(
            amino_acid_codon => amino_acid_codon_frequency/total_count for (amino_acid_codon, amino_acid_codon_frequency) in amino_acid_codon_frequencies
        )
        if !isempty(normalized_codon_frequencies[amino_acid])
            @assert abs(1-sum(values(normalized_codon_frequencies[amino_acid]))) <= eps()
        end
    end
    return normalized_codon_frequencies
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a protein sequence back to a possible DNA coding sequence using weighted random codon selection.

# Arguments
- `protein_sequence::BioSequences.LongAA`: The amino acid sequence to reverse translate

# Returns
- `BioSequences.LongDNA{2}`: A DNA sequence that would translate to the input protein sequence

# Details
Uses codon usage frequencies to randomly select codons for each amino acid, weighted by their 
natural occurrence. Each selected codon is guaranteed to translate back to the original amino acid.
"""
function reverse_translate(protein_sequence::BioSequences.LongAA)
    this_sequence = BioSequences.LongDNA{2}()
    codon_frequencies = Dict(a => Dict{Kmers.DNACodon, Int}() for a in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term]))
    for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
        amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
        codon_frequencies[amino_acid][codon] = get(codon_frequencies[amino_acid], codon, 0) + 1
    end    
    for amino_acid in protein_sequence
        # I'm collecting first because I'm worried about the keys and values not being sorted the same between queries, but that feels like it's not a viable worry
        collected = collect(codon_frequencies[amino_acid])
        codons = first.(collected)
        frequencies = last.(collected)
        chosen_codon_index = StatsBase.sample(1:length(codons), StatsBase.weights(frequencies))
        chosen_codon = codons[chosen_codon_index]
        @assert first(BioSequences.translate(BioSequences.LongDNA{2}(chosen_codon))) == amino_acid
        this_sequence *= chosen_codon
    end
    @assert BioSequences.translate(this_sequence) == protein_sequence
    return this_sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Optimizes the DNA sequence encoding for a given protein sequence using codon usage frequencies.

# Arguments
- `normalized_codon_frequencies`: Dictionary mapping amino acids to their codon frequencies
- `protein_sequence::BioSequences.LongAA`: Target protein sequence to optimize
- `n_iterations::Integer`: Number of optimization iterations to perform

# Algorithm
1. Creates initial DNA sequence through reverse translation
2. Iteratively generates new sequences by sampling codons based on their frequencies
3. Keeps track of the sequence with highest codon usage likelihood

# Returns
- `BioSequences.LongDNA{2}`: Optimized DNA sequence encoding the input protein
"""
function codon_optimize(;normalized_codon_frequencies, protein_sequence::BioSequences.LongAA, n_iterations)
    best_sequence = reverse_translate(protein_sequence)
    codons = last.(collect(Kmers.SpacedKmers{Kmers.DNACodon}(BioSequences.LongDNA{4}(best_sequence), 3)))
    initial_log_likelihood = -log10(1.0)
    for (codon, amino_acid) in collect(zip(codons, protein_sequence))
        this_codon_likelihood = normalized_codon_frequencies[amino_acid][codon]
        initial_log_likelihood -= log10(this_codon_likelihood)
    end
    best_likelihood = initial_log_likelihood

    ProgressMeter.@showprogress for i in 1:n_iterations
    # for iteration in 1:n_iterations
        this_sequence = BioSequences.LongDNA{2}()
        this_log_likelihood = -log10(1.0)
        for amino_acid in protein_sequence
            # I'm collecting first because I'm worried about the keys and values not being sorted the same between queries, but that feels like it's not a viable worry
            collected = collect(normalized_codon_frequencies[amino_acid])
            codons = first.(collected)
            frequencies = last.(collected)
            chosen_codon_index = StatsBase.sample(1:length(codons), StatsBase.weights(frequencies))
            chosen_codon = codons[chosen_codon_index]
            chosen_codon_frequency = frequencies[chosen_codon_index]
            this_log_likelihood += -log10(chosen_codon_frequency)
            this_sequence *= chosen_codon
        end
        if this_log_likelihood < best_likelihood
            best_likelihood = this_log_likelihood
            best_sequence = this_sequence
        end
        @assert BioSequences.translate(this_sequence) == protein_sequence
    end
    @show (best_likelihood)^-10 / (initial_log_likelihood)^-10
    return best_sequence
end

# function codon_optimize(;normalized_codon_frequencies, optimization_sequence::BioSequences.LongDNA, n_iterations)
#     protein_sequence = BioSequences.translate(optimization_sequence)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a dictionary of k-mer counts to a fixed-length numeric vector based on a predefined mapping.

# Arguments
- `kmer_to_index_map`: Dictionary mapping k-mer sequences to their corresponding vector indices
- `kmer_counts`: Dictionary containing k-mer sequences and their occurrence counts

# Returns
- A vector where each position corresponds to a k-mer count, with zeros for absent k-mers
"""
function kmer_counts_dict_to_vector(kmer_to_index_map, kmer_counts)
    kmer_counts_vector = zeros(length(kmer_to_index_map))
    for (kmer, count) in kmer_counts
        kmer_index = kmer_to_index_map[kmer]
        kmer_counts_vector[kmer_index] = count
    end
    return kmer_counts_vector
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Generate a sorted list of all possible k-mers for a given alphabet.

# Arguments
- `k::Integer`: Length of k-mers to generate
- `alphabet`: Collection of symbols (DNA, RNA, or amino acids) from BioSymbols

# Returns
- Sorted Vector of Kmers of the appropriate type (DNA, RNA, or amino acid)
"""
function generate_all_possible_kmers(k, alphabet)
    kmer_iterator = Iterators.product([alphabet for i in 1:k]...)
    kmer_vectors = collect.(vec(collect(kmer_iterator)))
    if eltype(alphabet) == BioSymbols.AminoAcid
        kmers = [Kmers.Kmer{BioSequences.AminoAcidAlphabet}(BioSequences.LongAA(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.DNA
        kmers = [Kmers.Kmer{BioSequences.DNAAlphabet{2}}(BioSequences.LongDNA{2}(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.RNA
        kmers = [Kmers.Kmer{BioSequences.RNAAlphabet{2}}(BioSequences.LongRNA{2}(kv)) for kv in kmer_vectors]
    else
        error()
    end
    return sort!(kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Generate all possible canonical k-mers of length `k` from the given `alphabet`.

For DNA/RNA sequences, returns unique canonical k-mers where each k-mer is represented by
the lexicographically smaller of itself and its reverse complement.
For amino acid sequences, returns all possible k-mers without canonicalization.

# Arguments
- `k`: Length of k-mers to generate
- `alphabet`: Vector of BioSymbols (DNA, RNA or AminoAcid)

# Returns
- Vector of k-mers, canonicalized for DNA/RNA alphabets
"""
function generate_all_possible_canonical_kmers(k, alphabet)
    kmers = generate_all_possible_kmers(k, alphabet)
    if eltype(alphabet) == BioSymbols.AminoAcid
        return kmers
    elseif eltype(alphabet) in (BioSymbols.DNA, BioSymbols.RNA)
        return unique!(BioSequences.canonical.(kmers))
    else
        error()
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the index position of a given k-mer in a sorted list of k-mers.

# Arguments
- `kmers`: A sorted vector of k-mers to search within
- `kmer`: The k-mer sequence to find

# Returns
Integer index position where `kmer` is found in `kmers`

# Throws
- `AssertionError`: If the k-mer is not found in the list
"""
function get_kmer_index(kmers, kmer)
    index = searchsortedfirst(kmers, kmer)
    @assert kmers[index] == kmer "$kmer not found in kmer list"
    return index
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a dense kmer counts table (canonical for DNA, stranded for RNA & AA) for each fasta provided in a list.
Scales very well for large numbers of organisms/fasta files, but not for k.
Recommended for k <= 13, although 17 may still be possible

Generate a dense k-mer frequency matrix from multiple FASTA files.

# Arguments
- `fasta_list`: Vector of paths to FASTA files
- `k`: Length of k-mers to count (recommended k ≤ 13)
- `alphabet`: Symbol specifying sequence type (`:DNA`, `:RNA`, or `:AA`)

# Returns 
Named tuple containing:
- `sorted_kmers`: Vector of sorted k-mers
- `kmer_counts_matrix`: Dense matrix where rows are k-mers and columns are sequences

# Details
- For DNA: Uses canonical k-mers (strand-neutral)
- For RNA/AA: Uses stranded k-mers
- Parallelized using Julia's multi-threading

# Performance
- Efficient for large numbers of sequences
- Memory usage grows exponentially with k
- For k > 13, use `fasta_list_to_sparse_counts_table` instead
"""
function fasta_list_to_dense_counts_table(;fasta_list, k, alphabet)
    k > 13 && error("use fasta_list_to_sparse_counts_table")
    if alphabet == :AA
        KMER_TYPE = BioSequences.AminoAcidAlphabet
        sorted_kmers = sort(generate_all_possible_kmers(k, AA_ALPHABET))
        COUNT = count_kmers
    elseif alphabet == :DNA
        KMER_TYPE = BioSequences.DNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_canonical_kmers(k, DNA_ALPHABET))
        COUNT = count_canonical_kmers
    elseif alphabet == :RNA
        KMER_TYPE = BioSequences.RNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_kmers(k, RNA_ALPHABET))
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    kmer_counts_matrix = zeros(length(sorted_kmers), length(fasta_list))
    progress = ProgressMeter.Progress(length(fasta_list))
    reenrantlock = ReentrantLock()
    Threads.@threads for (entity_index, fasta_file) in collect(enumerate(fasta_list))
        # Acquire the lock before updating the progress
        lock(reenrantlock) do
            # Update the progress meter
            ProgressMeter.next!(progress)
        end
        entity_mer_counts = COUNT(Kmers.Kmer{KMER_TYPE, k}, fasta_file)
        for (i, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[i, entity_index] = get(entity_mer_counts, kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a collection of biological sequences into a dense k-mer count matrix.

# Arguments
- `biosequences`: Collection of DNA, RNA, or amino acid sequences (BioSequence types)
- `k::Integer`: Length of k-mers to count (must be ≤ 13)

# Returns
Named tuple containing:
- `sorted_kmers`: Vector of all possible k-mers in sorted order
- `kmer_counts_matrix`: Dense matrix where rows are k-mers and columns are sequences

# Details
- For DNA sequences, counts canonical k-mers (both strands)
- For RNA and protein sequences, counts exact k-mers
- Uses parallel processing with threads
- For k > 13, use `fasta_list_to_sparse_counts_table` instead
"""
function biosequences_to_dense_counts_table(;biosequences, k)
    k > 13 && error("use fasta_list_to_sparse_counts_table")    
    if eltype(first(biosequences)) == BioSymbols.AminoAcid
        KMER_TYPE = BioSequences.AminoAcidAlphabet
        sorted_kmers = sort(generate_all_possible_kmers(k, AA_ALPHABET))
        COUNT = count_kmers
    elseif eltype(first(biosequences)) == BioSymbols.DNA
        KMER_TYPE = BioSequences.DNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_canonical_kmers(k, DNA_ALPHABET))
        COUNT = count_canonical_kmers
    elseif eltype(first(biosequences)) == BioSymbols.RNA
        KMER_TYPE = BioSequences.RNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_kmers(k, RNA_ALPHABET))
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    kmer_counts_matrix = zeros(length(sorted_kmers), length(biosequences))
    progress = ProgressMeter.Progress(length(biosequences))
    reenrantlock = ReentrantLock()
    Threads.@threads for (entity_index, biosequence) in collect(enumerate(biosequences))
        # Acquire the lock before updating the progress
        lock(reenrantlock) do
            # Update the progress meter
            ProgressMeter.next!(progress)
        end
        entity_mer_counts = COUNT(Kmers.Kmer{KMER_TYPE, k}, biosequence)
        for (i, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[i, entity_index] = get(entity_mer_counts, kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Create a sparse kmer counts table in memory for each fasta provided in a list
# """
# function fasta_list_to_counts_table(;fasta_list::AbstractVector{<:AbstractString}, k, alphabet)
#     if alphabet == :AA
#         KMER_TYPE = Kmers.AAKmer{k}
#         COUNT = count_kmers
#     elseif alphabet == :DNA
#         KMER_TYPE = Kmers.DNAKmer{k}
#         COUNT = count_canonical_kmers
#     elseif alphabet == :RNA
#         KMER_TYPE = Kmers.RNAKmer{k}
#         COUNT = count_kmers
#     else
#         error("invalid alphabet, please choose from :AA, :DNA, :RNA")
#     end
    
#     fasta_kmer_counts_dict = Dict()
#     progress = ProgressMeter.Progress(length(fasta_list))
#     reenrantlock = ReentrantLock()
#     Threads.@threads for fasta_file in fasta_list
#         # Acquire the lock before updating the progress
#         these_kmer_counts = COUNT(KMER_TYPE, fasta_file)
#         lock(reenrantlock) do
#             # Update the progress meter
#             ProgressMeter.next!(progress)
#             fasta_kmer_counts_dict[fasta_file] = these_kmer_counts
#         end
#     end
#     sorted_kmers = sort(collect(union(keys(x) for x in fasta_kmer_counts_dict)))
#     kmer_counts_matrix = zeros(length(sorted_kmers), length(fasta_list))
#     @info "populating sparse counts matrix..."
#     for (col, fasta) in enumerate(fasta_list)
#         for (row, kmer) in enumerate(sorted_kmers)
#             kmer_counts_matrix[row, col] = get(fasta_kmer_counts_dict[fasta], kmer, 0)
#         end
#     end
#     return (;sorted_kmers, kmer_counts_matrix)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a collection of biological sequences into a k-mer count matrix.

# Arguments
- `biosequences`: Vector of biological sequences (DNA, RNA, or Amino Acids)
- `k`: Length of k-mers to count

# Returns
Named tuple with:
- `sorted_kmers`: Vector of all unique k-mers found, lexicographically sorted
- `kmer_counts_matrix`: Sparse matrix where rows are k-mers and columns are sequences

# Details
- For DNA sequences, counts canonical k-mers (both strands)
- Uses parallel processing with Thread-safe progress tracking
- Memory efficient sparse matrix representation
- Supports DNA, RNA and Amino Acid sequences
"""
function biosequences_to_counts_table(;biosequences, k)
    if eltype(first(biosequences)) == BioSymbols.AminoAcid
        KMER_TYPE = Kmers.AAKmer{k}
        COUNT = count_kmers
    elseif eltype(first(biosequences)) == BioSymbols.DNA
        KMER_TYPE = Kmers.DNAKmer{k}
        COUNT = count_canonical_kmers
    elseif eltype(first(biosequences)) == BioSymbols.RNA
        KMER_TYPE = Kmers.RNAKmer{k}
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    
    kmer_counts = Vector{OrderedCollections.OrderedDict{KMER_TYPE, Int}}(undef, length(biosequences))
    progress = ProgressMeter.Progress(length(biosequences))
    reenrantlock = ReentrantLock()
    Threads.@threads for i in eachindex(biosequences)
        lock(reenrantlock) do
            ProgressMeter.next!(progress)
        end
        kmer_counts[i] = COUNT(KMER_TYPE, biosequences[i])
    end
    sorted_kmers = sort(collect(reduce(union, keys.(kmer_counts))))
    kmer_counts_matrix = SparseArrays.spzeros(Int, length(sorted_kmers), length(biosequences))
    @info "populating sparse counts matrix..."
    for (col, biosequence) in enumerate(biosequences)
        for (row, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[row, col] = get(kmer_counts[col], kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Normalize a distance matrix by dividing all elements by the maximum non-NaN value.

# Arguments
- `distance_matrix`: A matrix of distance values that may contain `NaN`, `nothing`, or `missing` values

# Returns
- Normalized distance matrix with values scaled to [0, 1] range

# Details
- Filters out `NaN`, `nothing`, and `missing` values when finding the maximum
- All elements are divided by the same maximum value to preserve relative distances
- If all values are NaN/nothing/missing, may return NaN values
"""
function normalize_distance_matrix(distance_matrix)
    max_non_nan_value = maximum(filter(x -> !isnan(x) && !isnothing(x) && !ismissing(x), vec(distance_matrix)))
    return distance_matrix ./ max_non_nan_value
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)
# """
# function count_matrix_to_probability_matrix(
#         counts_matrix,
#         probability_matrix_file = replace(counts_matrix_file, ".bin" => ".probability_matrix.bin")
#     )
#     probability_matrix = Mmap.mmap(probability_matrix_file, Array{Float64, 2}, size(counts_matrix))
#     if !isfile(probability_matrix_file)
#         println("creating new probability matrix $probability_matrix_file")
#         # probability_matrix .= count_matrix_to_probability_matrix(counts_matrix)
#         for (i, col) in enumerate(eachcol(counts_matrix))
#             probability_matrix[:, i] .= col ./ sum(col)
#         end
#     else
#         println("probability matrix found $probability_matrix_file")
#     end
#     return probability_matrix, probability_matrix_file
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a matrix of counts into a probability matrix by normalizing each column to sum to 1.0.

# Arguments
- `counts_matrix::Matrix{<:Number}`: Input matrix where each column represents counts/frequencies

# Returns
- `Matrix{Float64}`: Probability matrix where each column sums to 1.0
"""
function count_matrix_to_probability_matrix(counts_matrix)
    probability_matrix = zeros(size(counts_matrix))
    for (i, col) in enumerate(eachcol(counts_matrix))
        probability_matrix[:, i] .= col ./ sum(col)
    end
    return probability_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Convert a distance matrix into a Newick tree format using UPGMA hierarchical clustering.

# Arguments
- `distance_matrix`: Square matrix of pairwise distances between entities
- `labels`: Vector of labels corresponding to the entities in the distance matrix
- `outfile`: Path where the Newick tree file will be written

# Returns
- Path to the generated Newick tree file

# Details
Performs hierarchical clustering using the UPGMA (average linkage) method and 
converts the resulting dendrogram into Newick tree format. The branch lengths 
in the tree represent the heights from the clustering.
"""
function distance_matrix_to_newick(;distance_matrix, labels, outfile)
    # phage_names = phage_host_table[indices, :name]
    # this is equivalent to UPGMA
    tree = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)
    # reference_phage_indices = findall(x -> x in reference_phages, phage_names)
    newick = Dict()
    for row in 1:size(tree.merges, 1)
        left, right = tree.merges[row, :]
        if left < 0
            l = string(labels[abs(left)])
        else
            l = newick[left]
        end
        if right < 0
            r = string(labels[abs(right)])
        else
            r = newick[right]
        end
        height = tree.heights[row]
        newick[row] = "($l:$height, $r:$height)"
    end
    open(outfile, "w") do io
        println(io, newick[size(tree.merges, 1)] * ";")
    end
    return outfile
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# DEPRECATED: THIS WAS THE MEASURE WITH THE LEAST AGREEMENT TO EXISTING MEASURES LIKE BLAST AND % AVERAGE NUCLEOTIDE IDENTITY
# Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)
# """
# function counts_matrix_to_size_normalized_cosine_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     for entity_1_index in 1:n_entities
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             sa = sum(a)
#             sb = sum(b)
#             size_dist = 1-(min(sa, sb)/max(sa, sb))
#             cosine_dist = Distances.cosine_dist(a, b)
#             distances = filter(x -> x > 0, (size_dist, cosine_dist))
#             if isempty(distances)
#                 dist = 0.0
#             else
#                 dist = reduce(*, distances)
#             end
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = dist
#         end
#     end
#     return distance_matrix
# end

# """
# Create euclidean distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to total feature count magnitude (size)
# """
# function frequency_matrix_to_euclidean_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     ProgressMeter.@showprogress for entity_1_index in 1:n_entities
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = 
#                 Distances.euclidean(a, b)
#         end
#     end
#     return distance_matrix
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a Euclidean distance matrix from a column-major counts matrix
(features as rows and entities as columns), where distance is proportional
to total feature count magnitude (size).

Compute pairwise Euclidean distances between entity profiles in a counts matrix.

# Arguments
- `counts_table`: A matrix where rows represent features and columns represent entities (column-major format).
  Each column contains the feature counts/frequencies for one entity.

# Returns
- A symmetric N×N matrix of Euclidean distances between each pair of entities, where N is the number of entities.

# Details
- Parallelized computation using multi-threading
- Progress tracking via ProgressMeter
- Memory efficient: only upper triangle is computed, then mirrored
- Distance between entities increases with total feature magnitude differences
"""
function frequency_matrix_to_euclidean_distance_matrix(counts_table)
    n_entities = size(counts_table, 2)
    distance_matrix = zeros(n_entities, n_entities)

    # Initialize a thread-safe progress meter
    progress = ProgressMeter.Progress(n_entities * (n_entities - 1) ÷ 2, desc = "Computing distances", dt = 0.1)

    Threads.@threads for entity_1_index in 1:n_entities
        for entity_2_index in entity_1_index+1:n_entities
            a = counts_table[:, entity_1_index]
            b = counts_table[:, entity_2_index]
            dist = Distances.euclidean(a, b)
            distance_matrix[entity_1_index, entity_2_index] = dist
            distance_matrix[entity_2_index, entity_1_index] = dist
            ProgressMeter.next!(progress)  # Update the progress meter
        end
    end

    ProgressMeter.finish!(progress)  # Ensure the progress meter completes
    return distance_matrix
end

# didn't work
# function frequency_matrix_to_euclidean_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     progress = ProgressMeter.Progress(n_entities)
#     reenrantlock = ReentrantLock()
#     Threads.@threads for entity_1_index in 1:n_entities
#         lock(reenrantlock) do
#             ProgressMeter.next!(progress)
#         end
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = 
#                 Distances.euclidean(a, b)
#         end
#     end
#     return distance_matrix
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create cosine distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to cosine similarity (relative frequency)

Compute pairwise cosine distances between entities based on their feature distributions.

# Arguments
- `probability_matrix`: Column-major matrix where rows represent features and columns represent entities.
  Each column should contain frequency/probability values for one entity.

# Returns
- Symmetric matrix of size (n_entities × n_entities) containing pairwise cosine distances.
  Distance values range from 0 (identical distributions) to 1 (orthogonal distributions).

# Details
- Computes upper triangle and mirrors to lower triangle for efficiency
- Uses `Distances.cosine_dist` for the core computation
- Time complexity is O(n²) where n is the number of entities
"""
function frequency_matrix_to_cosine_distance_matrix(probability_matrix)
    n_entities = size(probability_matrix, 2)
    distance_matrix = zeros(n_entities, n_entities)
    for entity_1_index in 1:n_entities
        for entity_2_index in entity_1_index+1:n_entities
            a = probability_matrix[:, entity_1_index]
            b = probability_matrix[:, entity_2_index]
            distance_matrix[entity_1_index, entity_2_index] = 
                distance_matrix[entity_2_index, entity_1_index] = 
                Distances.cosine_dist(a, b)
        end
    end
    return distance_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a path of overlapping k-mers into a single DNA sequence.

# Arguments
- `kmer_path`: Vector of k-mers (DNA sequences) where each consecutive pair overlaps by k-1 bases

# Returns
- `BioSequences.LongDNA{2}`: Assembled DNA sequence from the k-mer path

# Description
Reconstructs the original DNA sequence by joining k-mers, validating that consecutive k-mers 
overlap correctly. The first k-mer is used in full, then each subsequent k-mer contributes 
its last base.
"""
function kmer_path_to_sequence(kmer_path)
    sequence = BioSequences.LongDNA{2}(first(kmer_path))
    for kmer in kmer_path[2:end]
        for i in 1:length(kmer)-1
            a = kmer[i]
            b = sequence[end-(length(kmer)-1)+i]
            @assert a == b
        end
        push!(sequence, kmer[end])
    end
    return sequence
end
    
# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function observe(records::AbstractVector{R}; outfile = "", error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
#     if isempty(outfile)
#         error("no file name supplied")
#     end
#     io = open(outfile, 'w')
#     fastx_io = FASTX.FASTQ.Writer(io)
#     for record in records
#         new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
#         new_seq_id = string(hash(new_seq)) * "-" * Random.randstring(32)
#     new_seq_description = FASTX.identifier(record)
#     quality = fill(UInt8(60), length(new_seq))
#     return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function random_fasta_record(;seed=rand(0:typemax(Int)), L = rand(0:Int(typemax(UInt16))))
#     id = Random.randstring(Int(ceil(log(L + 1))))
#     seq = BioSequences.randdnaseq(Random.seed!(seed), L)
#     return FASTX.FASTA.Record(id, seq)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generates a random FASTA record with a specified molecular type and sequence length.

# Arguments
- `moltype::Symbol=:DNA`: The type of molecule to generate (`:DNA`, `:RNA`, or `:AA` for amino acids).
- `seed`: The random seed used for sequence generation (default: a random integer).
- `L`: The length of the sequence (default: a random integer up to `typemax(UInt16)`).

# Returns
- A `FASTX.FASTA.Record` containing:
  - A randomly generated UUID identifier.
  - A randomly generated sequence of the specified type.

# Errors
- Throws an error if `moltype` is not one of `:DNA`, `:RNA`, or `:AA`.
"""
function random_fasta_record(;moltype::Symbol=:DNA, seed=rand(0:typemax(Int)), L = rand(0:Int(typemax(UInt16))))
    Random.seed!(seed)
    if moltype == :DNA
        seq = BioSequences.randdnaseq(StableRNGs.StableRNG(seed), L)
    elseif moltype == :RNA
        seq = BioSequences.randrnaseq(StableRNGs.StableRNG(seed), L)
    elseif moltype == :AA
        seq = BioSequences.randaaseq(StableRNGs.StableRNG(seed), L)
    else
        error("unrecognized molecule type: $(moltype) ! found in [:DNA, :RNA, :AA]")
    end
    # id = Mycelia.seq2sha256(seq)
    id = string(UUIDs.uuid4())
    return FASTX.FASTA.Record(id, seq)
end

"""
    observe(sequence::BioSequences.LongSequence{T}; error_rate=nothing, tech::Symbol=:illumina) where T

Simulates the “observation” of a biological polymer (DNA, RNA, or protein) by introducing realistic errors along with base‐quality scores.
The simulation takes into account both random and systematic error components. In particular, for technologies:
  
- **illumina**: (mostly substitution errors) the per‐base quality decays along the read (from ~Q40 at the start to ~Q20 at the end);
- **nanopore**: errors are more frequent and include both substitutions and indels (with overall lower quality scores, and an extra “homopolymer” penalty);
- **pacbio**: errors are dominated by indels (with quality scores typical of raw reads);
- **ultima**: (UG 100/ppmSeq™) correct bases are assigned very high quality (~Q60) while errors are extremely rare and, if they occur, are given a modest quality.

An error is introduced at each position with a (possibly position‐dependent) probability. For Illumina, the error probability increases along the read; additionally, if a base is part of a homopolymer run (length ≥ 3) and the chosen technology is one that struggles with homopolymers (nanopore, pacbio, ultima), then the local error probability is multiplied by a constant factor.

Returns a tuple `(new_seq, quality_scores)` where:
- `new_seq` is a `BioSequences.LongSequence{T}` containing the “observed” sequence (which may be longer or shorter than the input if insertions or deletions occur), and 
- `quality_scores` is a vector of integers representing the Phred quality scores (using the Sanger convention) for each base in the output sequence.
"""
function observe(sequence::BioSequences.LongSequence{T}; error_rate=nothing, tech::Symbol=:illumina) where T
    # Determine the appropriate alphabet based on the type T.
    if T <: BioSequences.DNAAlphabet
        alphabet = Mycelia.DNA_ALPHABET
    elseif T <: BioSequences.RNAAlphabet
        alphabet = Mycelia.RNA_ALPHABET
    else
        @assert T <: BioSequences.AminoAcidAlphabet "For amino acid sequences, T must be a subtype of BioSequences.AminoAcidAlphabet."
        alphabet = Mycelia.AA_ALPHABET
    end

    # Set a default baseline error rate if not provided.
    base_error_rate = isnothing(error_rate) ?
         (tech == :illumina  ? 0.005 :
          tech == :nanopore  ? 0.10  :
          tech == :pacbio    ? 0.11  :
          tech == :ultima    ? 1e-6  : error("Unknown technology")) : error_rate

    # Define error type probabilities (mismatch, insertion, deletion) for each technology.
    error_probs = Dict{Symbol, Float64}()
    if tech == :illumina
        error_probs[:mismatch]  = 0.90
        error_probs[:insertion] = 0.05
        error_probs[:deletion]  = 0.05
    elseif tech == :nanopore
        error_probs[:mismatch]  = 0.40
        error_probs[:insertion] = 0.30
        error_probs[:deletion]  = 0.30
    elseif tech == :pacbio
        error_probs[:mismatch]  = 0.20
        error_probs[:insertion] = 0.40
        error_probs[:deletion]  = 0.40
    elseif tech == :ultima
        error_probs[:mismatch]  = 0.95
        error_probs[:insertion] = 0.025
        error_probs[:deletion]  = 0.025
    end

    # Parameters for boosting error probability in homopolymer regions.
    homopolymer_threshold = 3    # If the homopolymer run length is ≥ 3, boost error probability.
    homopolymer_factor    = 2.0

    new_seq = BioSequences.LongSequence{T}()
    quality_scores = Int[]
    n = length(sequence)
    pos = 1

    while pos <= n
        base = sequence[pos]
        # For Illumina, increase error probability along the read; otherwise use a constant baseline.
        pos_error_rate = (tech == :illumina) ?
            base_error_rate * (1 + (pos - 1) / (n - 1)) : base_error_rate

        # Check for a homopolymer run by looking backwards.
        run_length = 1
        if pos > 1 && sequence[pos] == sequence[pos - 1]
            run_length = 2
            j = pos - 2
            while j ≥ 1 && sequence[j] == base
                run_length += 1
                j -= 1
            end
        end
        if run_length ≥ homopolymer_threshold && (tech in (:nanopore, :pacbio, :ultima))
            pos_error_rate *= homopolymer_factor
        end

        # Decide whether to observe the base correctly or introduce an error.
        if rand() > pos_error_rate
            # No error: add the correct base.
            push!(new_seq, base)
            push!(quality_scores, get_correct_quality(tech, pos, n))
            pos += 1
        else
            # An error occurs; choose the error type by sampling.
            r = rand()
            if r < error_probs[:mismatch]
                # Mismatch: choose a random base different from the true base.
                new_base = rand(filter(x -> x != base, alphabet))
                push!(new_seq, new_base)
                push!(quality_scores, get_error_quality(tech))
                pos += 1
            elseif r < (error_probs[:mismatch] + error_probs[:insertion])
                # Insertion: insert one or more random bases (simulate an extra insertion error),
                # then add the correct base.
                num_insertions = 1 + rand(Distributions.Poisson(pos_error_rate))
                for _ in 1:num_insertions
                    ins_base = rand(alphabet)
                    push!(new_seq, ins_base)
                    push!(quality_scores, get_error_quality(tech))
                end
                # Append the original base as a correct base.
                push!(new_seq, base)
                push!(quality_scores, get_correct_quality(tech, pos, n))
                pos += 1
            else
                # Deletion: skip adding the base (simulate a deletion).
                pos += 1
            end
        end
    end
    return new_seq, quality_scores
end

"""
    get_correct_quality(tech::Symbol, pos::Int, read_length::Int) -> Int

Simulates a Phred quality score (using the Sanger convention) for a correctly observed base.
For Illumina, the quality score is modeled to decay linearly from ~40 at the start to ~20 at the end of the read.
For other technologies, the score is sampled from a normal distribution with parameters typical for that platform.

Returns an integer quality score.
"""
function get_correct_quality(tech::Symbol, pos::Int, read_length::Int)
    if tech == :illumina
        q = 40 - 20 * (pos - 1) / (read_length - 1)
        q = clamp(round(Int, rand(Distributions.Normal(q, 2))), 20, 40)
        return q
    elseif tech == :nanopore
        q = clamp(round(Int, rand(Distributions.Normal(12, 2))), 10, 15)
        return q
    elseif tech == :pacbio
        q = clamp(round(Int, rand(Distributions.Normal(15, 2))), 12, 18)
        return q
    elseif tech == :ultima
        q = clamp(round(Int, rand(Distributions.Normal(60, 3))), 55, 65)
        return q
    else
        return 30
    end
end

"""
    get_error_quality(tech::Symbol) -> Int

Simulates a Phred quality score (using the Sanger convention) for a base observed with an error.
Error bases are assigned lower quality scores than correctly observed bases.
For Illumina, scores typically range between 5 and 15; for nanopore and pacbio, slightly lower values are used;
and for ultima, a modest quality score is assigned.

Returns an integer quality score.
"""
function get_error_quality(tech::Symbol)
    if tech == :illumina
        q = clamp(round(Int, rand(Distributions.Normal(10, 2))), 5, 15)
        return q
    elseif tech == :nanopore
        q = clamp(round(Int, rand(Distributions.Normal(7, 2))), 5, 10)
        return q
    elseif tech == :pacbio
        q = clamp(round(Int, rand(Distributions.Normal(7, 2))), 5, 10)
        return q
    elseif tech == :ultima
        q = clamp(round(Int, rand(Distributions.Normal(20, 3))), 15, 25)
        return q
    else
        return 10
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulate sequencing of a DNA/RNA record by introducing random errors at the specified rate.

# Arguments
- `record`: A FASTA or FASTQ record containing the sequence to be "observed"
- `error_rate`: Probability of error at each position (default: 0.0)

# Returns
A new FASTQ.Record with:
- Random UUID as identifier
- Original record's description 
- Modified sequence with introduced errors
- Generated quality scores
"""
function observe(record::R; error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    converted_sequence = convert_sequence(FASTX.sequence(record))
    new_seq, quality_scores = observe(converted_sequence, error_rate=error_rate)
    new_seq_id = UUIDs.uuid4()
    new_seq_description = FASTX.identifier(record)
    return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality_scores)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Determines the alphabet of a sequence. The function scans through `seq` only once:
- If a 'T' or 't' is found (and no 'U/u'), the sequence is classified as DNA.
- If a 'U' or 'u' is found (and no 'T/t'), it is classified as RNA.
- If both T and U occur, an error is thrown.
- If a character outside the canonical nucleotide and ambiguity codes is encountered,
  the sequence is assumed to be protein.
- If neither T nor U are found, the sequence is assumed to be DNA.
"""
function detect_alphabet(seq::AbstractString)::Symbol
    hasT = false
    hasU = false
    # Define allowed nucleotide characters (both for DNA and RNA, including common ambiguity codes)
    # TODO: define this by merging the alphabets from BioSymbols
    valid_nucleotides = "ACGTacgtACGUacguNRYSWKMBDHnryswkmbdh"
    for c in seq
        if c == 'T' || c == 't'
            hasT = true
        elseif c == 'U' || c == 'u'
            hasU = true
        elseif !(c in valid_nucleotides)
            # If an unexpected character is encountered, assume it's a protein sequence.
            return :AA
        end
        if hasT && hasU
            throw(ArgumentError("Sequence contains both T and U, ambiguous alphabet"))
        end
    end
    if hasT
        return :DNA
    elseif hasU
        return :RNA
    else
        # In the absence of explicit T or U, default to DNA.
        return :DNA
    end
end

# Specific dispatch for DNA sequences
function detect_alphabet(sequence::BioSequences.LongDNA)
    return :DNA
end

# Specific dispatch for RNA sequences
function detect_alphabet(sequence::BioSequences.LongRNA)
    return :RNA
end

# Specific dispatch for protein/amino acid sequences
function detect_alphabet(sequence::BioSequences.LongAA)
    return :AA
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Converts the given sequence (output from FASTX.sequence) into the appropriate BioSequence type:
- DNA sequences are converted using `BioSequences.LongDNA`
- RNA sequences are converted using `BioSequences.LongRNA`
- AA sequences are converted using `BioSequences.LongAA`
"""
function convert_sequence(seq::AbstractString)
    alphabet = detect_alphabet(seq)
    if alphabet == :DNA
        return BioSequences.LongDNA{4}(seq)
    elseif alphabet == :RNA
        return BioSequences.LongRNA{4}(seq)
    elseif alphabet == :AA
        return BioSequences.LongAA(seq)
    else
        throw(ArgumentError("Unrecognized alphabet type"))
    end
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function observe(records::AbstractVector{R};
#                 weights=ones(length(records)),
#                 N = length(records),
#                 outfile = "",
#                 error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
#     if isempty(outfile)
#         error("no file name supplied")
#     end
#     io = open(outfile, "w")
#     fastx_io = FASTX.FASTA.Writer(io)
#     for i in 1:N
#         record = StatsBase.sample(records, StatsBase.weights(weights))
#         new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
#         new_seq_id = Random.randstring(Int(ceil(log(length(new_seq) + 1))))
#         new_seq_description = FASTX.identifier(record)
#         observed_record = FASTX.FASTA.Record(new_seq_id, new_seq_description, new_seq)
#         write(fastx_io, observed_record)
#     end
#     close(fastx_io)
#     close(io)
#     return outfile
# end

# currently this is only for amino acid sequences, expand to include DNA and RNA via multiple dispatch
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a single random mutation in an amino acid sequence.

# Arguments
- `reference_sequence`: Input amino acid sequence to be mutated

# Returns
- `mutant_sequence`: The sequence after applying the mutation
- `haplotype`: A `SequenceVariation.Haplotype` object containing the mutation details

# Details
Performs one of three possible mutation types:
- Substitution: Replace one amino acid with another
- Insertion: Insert 1+ random amino acids at a position
- Deletion: Remove 1+ amino acids from a position

Insertion and deletion sizes follow a truncated Poisson distribution (λ=1, min=1).
"""
function mutate_sequence(reference_sequence)
    i = rand(1:length(reference_sequence))
    mutation_type = rand([SequenceVariation.Substitution, SequenceVariation.Insertion, SequenceVariation.Deletion])
    if mutation_type == SequenceVariation.Substitution
        # new_amino_acid = BioSequences.AA_A
        new_amino_acid = rand(amino_acids)
        mutation_string = "$(reference_sequence[i])$(i)$(new_amino_acid)"
    else
        indel_size = rand(Distributions.truncated(Distributions.Poisson(1), lower=1))
        if mutation_type == SequenceVariation.Insertion
            inserted_sequence = join([rand(amino_acids) for i in 1:indel_size], "")
            mutation_string = "$(i)$(inserted_sequence)"
        else
            @assert mutation_type == SequenceVariation.Deletion
            i_stop = min(i+indel_size, length(reference_sequence))
            mutation_string = "Δ$(i)-$(i_stop)"
        end
    end
    # println(mutation_string)
    mutation = SequenceVariation.Variation(reference_sequence, mutation_string)
    haplotype = SequenceVariation.Haplotype(reference_sequence, [mutation])
    mutant_sequence = SequenceVariation.reconstruct(haplotype)
    return mutant_sequence, haplotype
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return proportion of matched bases in alignment to total matches + edits.

Calculate the accuracy of a sequence alignment by computing the ratio of matched bases 
to total alignment operations (matches + edits).

# Arguments
- `alignment_result`: Alignment result object containing `total_matches` and `total_edits` fields

# Returns
Float64 between 0.0 and 1.0 representing alignment accuracy, where:
- 1.0 indicates perfect alignment (all matches)
- 0.0 indicates no matches
"""
function assess_alignment_accuracy(alignment_result)
    return alignment_result.total_matches / (alignment_result.total_matches + alignment_result.total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Used to determine which orientation provides an optimal alignment for initiating path likelihood analyses in viterbi analysis

Compare alignment scores between a query k-mer and an observed k-mer in both forward and
reverse complement orientations to determine optimal alignment.

# Arguments
- `kmer`: Query k-mer sequence to align
- `observed_kmer`: Target k-mer sequence to align against

# Returns
A tuple containing:
- `alignment_result`: The alignment result object for the optimal orientation
- `orientation`: Boolean indicating orientation (`true` = forward, `false` = reverse complement, `missing` = tied scores)

# Details
- Performs pairwise alignment in both orientations using `assess_alignment()`
- Calculates accuracy scores using `assess_alignment_accuracy()`
- For tied alignment scores, randomly selects one orientation
- Uses BioSequences.reverse_complement for reverse orientation comparison
"""
function assess_optimal_kmer_alignment(kmer, observed_kmer)

    forward_alignment_result = assess_alignment(kmer, observed_kmer)
    forward_alignment_accuracy = assess_alignment_accuracy(forward_alignment_result)

    reverse_alignment_result = assess_alignment(kmer, BioSequences.reverse_complement(observed_kmer))
    reverse_alignment_accuracy = assess_alignment_accuracy(reverse_alignment_result)

    if forward_alignment_accuracy > reverse_alignment_accuracy
        alignment_result = forward_alignment_result
        orientation = true
    elseif forward_alignment_accuracy < reverse_alignment_accuracy
        alignment_result = reverse_alignment_result
        orientation = false
    elseif forward_alignment_accuracy == reverse_alignment_accuracy
        alignment_result, orientation = rand(((forward_alignment_result, missing), (reverse_alignment_result, missing)))
    end

    return (alignment_result, orientation)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Aligns two sequences using the Levenshtein distance and returns the total number of matches and edits.

# Arguments
- `a::AbstractString`: The first sequence to be aligned.
- `b::AbstractString`: The second sequence to be aligned.

# Returns
- `NamedTuple{(:total_matches, :total_edits), Tuple{Int, Int}}`: A named tuple containing:
    - `total_matches::Int`: The total number of matching bases in the alignment.
    - `total_edits::Int`: The total number of edits (insertions, deletions, substitutions) in the alignment.
"""
function assess_alignment(a, b)
    pairwise_alignment = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), a, b)
    alignment_result = BioAlignments.alignment(pairwise_alignment)
    total_aligned_bases = BioAlignments.count_aligned(alignment_result)
    total_matches = Int(BioAlignments.count_matches(alignment_result))
    total_edits = Int(total_aligned_bases - total_matches)
    return (total_matches = total_matches, total_edits = total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Canonicalizes the k-mer counts in the given dictionary.

This function iterates over the provided dictionary `kmer_counts`, which maps k-mers to their respective counts. For each k-mer that is not in its canonical form, it converts the k-mer to its canonical form and updates the count in the dictionary accordingly. If the canonical form of the k-mer already exists in the dictionary, their counts are summed. The original non-canonical k-mer is then removed from the dictionary.

# Arguments
- `kmer_counts::Dict{BioSequences.Kmer, Int}`: A dictionary where keys are k-mers and values are their counts.

# Returns
- The input dictionary `kmer_counts` with all k-mers in their canonical form, sorted by k-mers.
"""
function canonicalize_kmer_counts!(kmer_counts)
    for (kmer, count) in kmer_counts
        if !BioSequences.iscanonical(kmer)
            canonical_kmer = BioSequences.canonical(kmer)
            if haskey(kmer_counts, canonical_kmer)
                kmer_counts[canonical_kmer] += count
            else
                kmer_counts[canonical_kmer] = count
            end
            delete!(kmer_counts, kmer)
        end
    end
    return sort!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalize k-mer counts into a canonical form by creating a non-mutating copy.

# Arguments
- `kmer_counts`: Dictionary or collection of k-mer count data

# Returns
- A new normalized k-mer count collection
"""
function canonicalize_kmer_counts(kmer_counts)
    return canonicalize_kmer_counts!(copy(kmer_counts))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count canonical k-mers in biological sequences. A canonical k-mer is the lexicographically 
smaller of a DNA sequence and its reverse complement, ensuring strand-independent counting.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer size and structure
- `sequences`: Iterator of biological sequences to analyze

# Returns
- `Dict{KMER_TYPE,Int}`: Dictionary mapping canonical k-mers to their counts
"""
function count_canonical_kmers(::Type{KMER_TYPE}, sequences) where KMER_TYPE
    kmer_counts = count_kmers(KMER_TYPE, sequences)
    return canonicalize_kmer_counts!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of each k-mer in a DNA sequence.

# Arguments
- `::Type{Kmers.Kmer{A,K}}`: K-mer type with alphabet A and length K
- `sequence::BioSequences.LongSequence`: Input DNA sequence to analyze

# Returns
A sorted dictionary mapping each k-mer to its frequency count in the sequence.

# Type Parameters
- `A <: BioSequences.DNAAlphabet`: DNA alphabet type
- `K`: Length of k-mers
"""
function count_kmers(::Type{Kmers.Kmer{A, K}}, sequence::BioSequences.LongSequence) where {A <: BioSequences.DNAAlphabet, K}
    # return sort(StatsBase.countmap(Kmers.FwDNAMers{K}(sequence)))
    return sort(StatsBase.countmap([kmer for (kmer, index) in Kmers.UnambiguousDNAMers{K}(sequence)]))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of each k-mer in an RNA sequence.

# Arguments
- `Kmer`: Type parameter specifying the k-mer length K and RNA alphabet
- `sequence`: Input RNA sequence to analyze

# Returns
- `Dict{Kmers.Kmer, Int}`: Sorted dictionary mapping each k-mer to its frequency count
"""
function count_kmers(::Type{Kmers.Kmer{A, K}}, sequence::BioSequences.LongSequence) where {A <: BioSequences.RNAAlphabet, K}
    # return sort(StatsBase.countmap(Kmers.FwRNAMers{K}(sequence)))
    return sort(StatsBase.countmap([kmer for (kmer, index) in Kmers.UnambiguousRNAMers{K}(sequence)]))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of amino acid k-mers in a biological sequence.

# Arguments
- `Kmers.Kmer{A,K}`: Type parameter specifying amino acid alphabet (A) and k-mer length (K)
- `sequence`: Input biological sequence to analyze

# Returns
A sorted dictionary mapping each k-mer to its frequency count in the sequence.
"""
function count_kmers(::Type{Kmers.Kmer{A, K}}, sequence::BioSequences.LongSequence) where {A <: BioSequences.AminoAcidAlphabet, K}
    return sort(StatsBase.countmap(Kmers.FwAAMers{K}(sequence)))
end

# TODO add a way to handle ambiguity or not
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count the frequency of amino acid k-mers in a biological sequence.

# Arguments
- `Kmers.Kmer{A,K}`: Type parameter specifying amino acid alphabet (A) and k-mer length (K)
- `sequence`: Input biological sequence to analyze

# Returns
A sorted dictionary mapping each k-mer to its frequency count in the sequence.
"""
function count_kmers(::Type{KMER_TYPE}, record::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    # TODO: need to figure out how to infer the sequence type
    if eltype(KMER_TYPE) == BioSymbols.DNA
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongDNA{4}, record))
    elseif eltype(KMER_TYPE) == BioSymbols.RNA
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongRNA{4}, record))
    elseif eltype(KMER_TYPE) == BioSymbols.AminoAcid
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongAA, record))
    else
        @error KMER_TYPE
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers across multiple sequence records and return a sorted frequency table.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer length (e.g., `DNAKmer{3}` for 3-mers)
- `records`: Vector of FASTA/FASTQ records to analyze

# Returns
- `Dict{KMER_TYPE, Int}`: Sorted dictionary mapping k-mers to their frequencies
"""
function count_kmers(::Type{KMER_TYPE}, records::AbstractVector{T}) where {KMER_TYPE, T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    kmer_counts = count_kmers(KMER_TYPE, first(records))
    for record in records[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, record)
        merge!(+, kmer_counts, _kmer_counts)
    end
    sort!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Counts k-mer occurrences in biological sequences from a FASTA/FASTQ reader.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer length and encoding (e.g., `DNAKmer{4}` for 4-mers)
- `sequences`: A FASTA or FASTQ reader containing the biological sequences to analyze

# Returns
A dictionary mapping k-mers to their counts in the input sequences
"""
function count_kmers(::Type{KMER_TYPE}, sequences::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Reader, FASTX.FASTQ.Reader}}
    return count_kmers(KMER_TYPE, collect(sequences))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers across multiple FASTA/FASTQ files and merge the results.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer length (e.g., `DNAKmer{4}` for 4-mers)
- `fastx_files`: Vector of paths to FASTA/FASTQ files

# Returns
- `Dict{KMER_TYPE, Int}`: Dictionary mapping k-mers to their total counts across all files
"""
function count_kmers(::Type{KMER_TYPE}, fastx_files::AbstractVector{T}) where {KMER_TYPE, T <: AbstractString}
    kmer_counts = count_kmers(KMER_TYPE, first(fastx_files))
    for file in fastx_files[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, file)
        kmer_counts = merge!(+, kmer_counts, _kmer_counts)
    end
    return kmer_counts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Count k-mers in a FASTA/FASTQ file and return their frequencies.

# Arguments
- `KMER_TYPE`: Type parameter specifying the k-mer type (e.g., `DNAKmer{K}`)
- `fastx_file`: Path to input FASTA/FASTQ file

# Returns
- `Dict{KMER_TYPE, Int}`: Dictionary mapping each k-mer to its frequency
"""
function count_kmers(::Type{KMER_TYPE}, fastx_file::AbstractString) where {KMER_TYPE}
    fastx_io = open_fastx(fastx_file)
    kmer_counts = count_kmers(KMER_TYPE, fastx_io)
    close(fastx_io)
    return kmer_counts
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

Check if two biological sequences are equivalent, considering both direct and reverse complement matches.

# Arguments
- `a`: First biological sequence (BioSequence or compatible type)
- `b`: Second biological sequence (BioSequence or compatible type)

# Returns
- `Bool`: `true` if sequences are identical or if one is the reverse complement of the other, `false` otherwise
"""
function is_equivalent(a, b)
    a == b || a == BioSequences.reverse_complement(b)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Converts a path of edges in a kmer graph into a DNA sequence by concatenating overlapping kmers.

# Arguments
- `kmer_graph`: A directed graph where vertices represent kmers and edges represent overlaps
- `edge_path`: Vector of edges representing a path through the graph

# Returns
A `BioSequences.LongDNASeq` containing the merged sequence represented by the path

# Details
The function:
1. Takes the first kmer from the source vertex of first edge
2. For each edge, handles orientation (forward/reverse complement)
3. Verifies overlaps between consecutive kmers
4. Concatenates unique bases to build final sequence
"""
function edge_path_to_sequence(kmer_graph, edge_path)
    edge = first(edge_path)
    sequence = BioSequences.LongDNASeq(kmers[edge.src])
    if !kmer_graph.eprops[edge][:orientations].source_orientation
        sequence = BioSequences.reverse_complement(sequence)
    end
    for edge in edge_path
        destination = BioSequences.LongDNASeq(kmers[edge.dst])
        if !kmer_graph.eprops[edge][:orientations].destination_orientation
            destination = BioSequences.reverse_complement(destination)
        end
        sequence_suffix = sequence[end-length(destination)+2:end]
        destination_prefix = destination[1:end-1]
        @assert sequence_suffix == destination_prefix
        push!(sequence, destination[end])
    end
    sequence
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

Calculate the maximum number of possible canonical k-mers for a given alphabet.

# Arguments
- `k::Integer`: Length of k-mer
- `ALPHABET::Vector{Char}`: Character set (nucleotides or amino acids)

# Returns
- `Int`: Maximum number of possible canonical k-mers

# Details
- For amino acids (AA_ALPHABET): returns total possible k-mers
- For nucleotides: returns half of total possible k-mers (canonical form)
- Requires odd k-mer length for nucleotide alphabets
"""
function determine_max_canonical_kmers(k, ALPHABET)
    max_possible_kmers = determine_max_possible_kmers(k, ALPHABET)
    if ALPHABET == AA_ALPHABET
        return max_possible_kmers
    else
        @assert isodd(k) "this calculation is not valid for even length kmers"
        return Int(max_possible_kmers / 2)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the total number of possible unique k-mers that can be generated from a given alphabet.

# Arguments
- `k`: Length of k-mers to consider
- `ALPHABET`: Vector containing the allowed characters/symbols

# Returns
- Integer representing the maximum number of possible unique k-mers (|Σ|ᵏ)
"""
function determine_max_possible_kmers(k, ALPHABET)
    return length(ALPHABET)^k
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate reaction velocity using Michaelis-Menten enzyme kinetics equation.

# Arguments
- `s::Vector{Float64}`: Substrate concentration(s) [mol/L]
- `p::Vector{Float64}`: Parameters vector where:
    - `p[1]`: vmax - Maximum reaction velocity [mol/(L⋅s)]
    - `p[2]`: km - Michaelis constant [mol/L]

# Returns
- `v::Vector{Float64}`: Reaction velocity [mol/(L⋅s)]

# Description
Implements the Michaelis-Menten equation: v = (vmax * s)/(km + s)

# Reference
Michaelis L., Menten M.L. (1913). Die Kinetik der Invertinwirkung. 
Biochem Z 49:333-369.
"""
# Michaelis–Menten
function calculate_v(s,p)
    vmax = p[1]
    km = p[2]
    v = (vmax .* s) ./ (km .+ s)
    return v
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

# seqkit concat $(cat fasta_files.txt) > merged.fasta
# seqtk seq -L $(cat fasta_files.txt) > merged.fasta
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Join fasta files without any regard to record uniqueness.

A cross-platform version of `cat *.fasta > joint.fasta`

See merge_fasta_files

Concatenate multiple FASTA files into a single output file by simple appending.

# Arguments
- `files`: Vector of paths to input FASTA files
- `file`: Path where the concatenated output will be written

# Returns
- Path to the output concatenated file

# Details
Platform-independent implementation of `cat *.fasta > combined.fasta`.
Files are processed sequentially with a progress indicator.
"""
function concatenate_files(;files, file)
    close(open(file, "w"))
    ProgressMeter.@showprogress for f in files
        # stderr=file_path
        run(pipeline(`cat $(f)`, stdout=file, append=true))
    end
    return file
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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function path_to_sequence(kmer_graph, path)
#     sequence = BioSequences.LongDNASeq(oriented_kmer_to_sequence(kmer_graph, first(path)))
#     for oriented_kmer in path[2:end]
#         nucleotide = last(oriented_kmer_to_sequence(kmer_graph, oriented_kmer))
#         push!(sequence, nucleotide)
#     end
#     return sequence
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function oriented_kmer_to_sequence(kmer_graph, oriented_kmer)
#     kmer_sequence = kmer_graph.kmers[oriented_kmer.index]
#     if !oriented_kmer.orientation
#         kmer_sequence = BioSequences.reverse_complement(kmer_sequence)
#     end
#     return kmer_sequence
# end

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

# this is faster than NCBI version
# run(pipeline(
#         `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon 1 --children --as-json-lines`,
#         `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli dataformat tsv taxonomy --template tax-summary`
#     )
# )
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Retrieves and formats the complete NCBI taxonomy hierarchy into a structured DataFrame.

# Details
- Automatically sets up taxonkit environment and downloads taxonomy database if needed
- Starts from root taxid (1) and includes all descendant taxa
- Reformats lineage information into separate columns for each taxonomic rank

# Returns
DataFrame with columns:
- `taxid`: Taxonomy identifier
- `lineage`: Full taxonomic lineage string
- `taxid_lineage`: Lineage with taxonomy IDs
- Individual rank columns:
  - superkingdom, kingdom, phylum, class, order, family, genus, species
  - corresponding taxid columns (e.g., superkingdom_taxid)

# Dependencies
Requires taxonkit (installed automatically via Bioconda)
"""
function list_full_taxonomy()
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    p = pipeline(
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit list --ids 1`,
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --add-prefix --fill-miss-rank --show-lineage-taxids --format '{k};{K};{p};{c};{o};{f};{g};{s}'`
    )
    data, header = uCSV.read(open(p), delim='\t')
    header = ["taxid", "lineage", "taxid_lineage"]
    table = DataFrames.DataFrame(data, header)
    ranks = [
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]
    for rank in ranks
        table[!, rank] .= ""
        table[!, "$(rank)_taxid"] .= ""
    end

    for (i, row) in enumerate(DataFrames.eachrow(table))
        for (rank, x) in zip(ranks, split(row["lineage"], ';'))
            table[i, rank] = x
        end
        for (rank, x) in zip(ranks, split(row["taxid_lineage"], ';'))
            table[i, "$(rank)_taxid"] = x
        end
    end
    return table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return an ordered list of taxonomic ranks from highest (top) to lowest (species).

# Arguments
- `synonyms::Bool=false`: If true, includes alternative names for certain ranks (e.g. "domain" for "superkingdom")

# Returns
- `Vector{String}`: An array of taxonomic rank names in hierarchical order
"""
function list_ranks(;synonyms=false)
    if !synonyms
        return [
            "top",
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        ]
    else
        return [
            "top",
            "superkingdom/domain",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        ]
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a DataFrame containing the top-level taxonomic nodes.

The DataFrame has two fixed rows representing the most basic taxonomic classifications:
- taxid=0: "unclassified"
- taxid=1: "root"

Returns
-------
DataFrame
    Columns:
    - taxid::Int : Taxonomic identifier
    - name::String : Node name
"""
function list_toplevel()
    return DataFrames.DataFrame(taxid=[0, 1], name=["unclassified", "root"])
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

List all taxonomic entries at the specified rank level.

# Arguments
- `rank::String`: Taxonomic rank to query. Must be one of:
  - "top" (top level)
  - "superkingdom"/"domain"  
  - "kingdom"
  - "phylum" 
  - "class"
  - "order"
  - "family"
  - "genus"
  - "species"

# Returns
DataFrame with columns:
- `taxid`: NCBI taxonomy ID
- `name`: Scientific name at the specified rank
"""
function list_rank(rank)
    if rank == "top"
        return list_toplevel()
    else
        ranks_to_shorthand = Dict(
            "superkingdom" => "k",
            "kingdom" => "K",
            "phylum" => "p",
            "class" => "c",
            "order" => "o",
            "family" => "f",
            "genus" => "g",
            "species" => "s"
        )
        shorthand = ranks_to_shorthand[rank]
        p = pipeline(
            `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit list --ids 1`,
            `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit filter --equal-to "$(rank)"`,
            `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit reformat --taxid-field 1 --format "{$(shorthand)}"`
        )
        data, header = uCSV.read(open(p), delim='\t')
        header = ["taxid", "name"]
        return DataFrames.DataFrame(data, header)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns an array of all taxonomic superkingdoms (e.g., Bacteria, Archaea, Eukaryota).

# Returns
- `Vector{String}`: Array containing names of all superkingdoms in the taxonomy database
"""
function list_superkingdoms()
    return list_rank("superkingdom")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists all taxonomic kingdoms in the database.

Returns a vector of kingdom names as strings. Kingdoms represent the highest
major taxonomic rank in biological classification.

# Returns
- `Vector{String}`: Array of kingdom names
"""
function list_kingdoms()
    return list_rank("kingdom")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted list of all unique phyla in the database.
"""
function list_phylums()
    return list_rank("phylum")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns an array of all taxonomic classes in the database.

Classes represent a major taxonomic rank between phylum and order in biological classification.

# Returns
- `Vector{String}`: Array of class names sorted alphabetically
"""
function list_classes()
    return list_rank("class")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Lists all orders in the taxonomic database.

Returns a vector of strings containing valid order names according to current mycological taxonomy.
Uses the underlying `list_rank()` function with rank="order".

# Returns
- `Vector{String}`: Alphabetically sorted list of order names
"""
function list_orders()
    return list_rank("order")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted vector of all family names present in the database.
"""
function list_families()
    return list_rank("family")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted vector of all genera names present in the database.
"""
function list_genera()
    return list_rank("genus")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns a sorted vector of all species names present in the database.
"""
function list_species()
    return list_rank("species")
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function list_subtaxa(taxid)
#     return parse.(Int, filter(!isempty, strip.(readlines(`conda run --no-capture-output -n taxonkit taxonkit list --ids 10239`))))
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns an array of Integer taxon IDs representing all sub-taxa under the specified taxonomic ID.

# Arguments
- `taxid`: NCBI taxonomy identifier for the parent taxon

# Returns
Vector{Int} containing all descendant taxon IDs

# Details
- Requires taxonkit to be installed via Bioconda
- Automatically sets up taxonkit database if not present
- Uses local taxonomy database in ~/.taxonkit/
"""
function list_subtaxa(taxid)
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    return parse.(Int, filter(!isempty, strip.(readlines(`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit list --ids $(taxid)`))))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert scientific name(s) to NCBI taxonomy ID(s) using taxonkit.

# Arguments
- `name::AbstractString`: Scientific name(s) to query. Can be a single name or multiple names separated by newlines.

# Returns
- `DataFrame` with columns:
  - `name`: Input scientific name
  - `taxid`: NCBI taxonomy ID
  - `rank`: Taxonomic rank (e.g., "species", "genus")

# Dependencies
Requires taxonkit package (installed automatically via Bioconda)
"""
function name2taxid(name)
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    p = pipeline(`echo $(name)`, `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n taxonkit taxonkit name2taxid --show-rank`)
    data, header = uCSV.read(open(p), delim='\t')
    header = ["name", "taxid", "rank"]
    return DataFrames.DataFrame(data, header)
end

# other ncbi-datasets reports that I didn't find as useful initially
# run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon 10114 --report names`))
# x = JSON.parse(open(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon "rattus norvegicus"`)))
# run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download taxonomy taxon 33554 --children`))

# more useful
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a vector of NCBI taxonomy IDs into a detailed taxonomy table using NCBI Datasets CLI.

# Arguments
- `taxids::AbstractVector{Int}`: Vector of NCBI taxonomy IDs to query

# Returns
- `DataFrame`: Table containing taxonomy information with columns including:
  - tax_id
  - species
  - genus
  - family
  - order
  - class
  - phylum
  - kingdom

# Dependencies
Requires ncbi-datasets-cli Conda package (automatically installed if missing)
"""
function taxids2ncbi_taxonomy_table(taxids::AbstractVector{Int})
    Mycelia.add_bioconda_env("ncbi-datasets-cli")
    joint_table = DataFrames.DataFrame()
    ProgressMeter.@showprogress for taxid in taxids
        cmd1 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets summary taxonomy taxon $(taxid) --as-json-lines`
        cmd2 = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli dataformat tsv taxonomy --template tax-summary`
        io = open(pipeline(cmd1, cmd2))
        append!(joint_table, DataFrames.DataFrame(uCSV.read(io, delim='\t', header=1)), promote=true)
    end
    return joint_table
end

# more complete
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert NCBI taxonomic IDs to their complete taxonomic lineage information using taxonkit.

# Arguments
- `taxids::AbstractVector{Int}`: Vector of NCBI taxonomy IDs

# Returns
A DataFrame with columns:
- `taxid`: Original query taxonomy ID
- `lineage`: Full taxonomic lineage as semicolon-separated string
- `lineage-taxids`: Corresponding taxonomy IDs for each rank in lineage
- `lineage-ranks`: Taxonomic ranks for each level in lineage
"""
# function taxids2taxonkit_lineage_table(taxids::AbstractVector{Int})
function taxids2taxonkit_full_lineage_table(taxids::AbstractVector{Int})
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    f = tempname()
    open(f, "w") do io
        for taxid in taxids
            println(io, taxid)
        end
    end
    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n taxonkit taxonkit lineage --show-lineage-taxids --show-lineage-ranks $(f)`
    data, header = uCSV.read(open(pipeline(cmd)), delim='\t', header=false)
    rm(f)
    header = ["taxid", "lineage", "lineage-taxids", "lineage-ranks"]
    return DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert taxonomic IDs to a structured lineage rank mapping.

Takes a vector of taxonomic IDs and returns a nested dictionary mapping each input taxid 
to its complete taxonomic lineage information. For each taxid, creates a dictionary where:
- Keys are taxonomic ranks (e.g., "species", "genus", "family")
- Values are NamedTuples containing:
  - `lineage::String`: The taxonomic name at that rank
  - `taxid::Union{Int, Missing}`: The corresponding taxonomic ID (if available)

Excludes "no rank" entries from the final output.

Returns:
    Dict{Int, Dict{String, NamedTuple{(:lineage, :taxid), Tuple{String, Union{Int, Missing}}}}}
"""
function taxids2taxonkit_taxid2lineage_ranks(taxids::AbstractVector{Int})
    table = taxids2taxonkit_full_lineage_table(taxids)
    # table = taxids2taxonkit_lineage_table(taxids)
    taxid_to_lineage_ranks = Dict{Int, Dict{String, @NamedTuple{lineage::String, taxid::Union{Int, Missing}}}}()
    for row in DataFrames.eachrow(table)
        lineage_ranks = String.(split(row["lineage-ranks"], ';'))
        lineage_taxids = [something(tryparse(Int, x), missing) for x in split(row["lineage-taxids"], ';')]
        # lineage_taxids = something.(tryparse.(Int, split(row["lineage-taxids"], ';')), missing)
        lineage = String.(split(row["lineage"], ';'))
        row_dict = Dict(rank => (;lineage, taxid) for (lineage, rank, taxid) in zip(lineage, lineage_ranks, lineage_taxids))
        delete!(row_dict, "no rank")
        taxid_to_lineage_ranks[row["taxid"]] = row_dict
    end
    return taxid_to_lineage_ranks
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a vector of taxonomy IDs to a summarized lineage table using taxonkit.

# Arguments
- `taxids::AbstractVector{Int}`: Vector of NCBI taxonomy IDs

# Returns
DataFrame with the following columns:
- `taxid`: Original input taxonomy ID
- `species_taxid`, `species`: Species level taxonomy ID and name
- `genus_taxid`, `genus`: Genus level taxonomy ID and name  
- `family_taxid`, `family`: Family level taxonomy ID and name
- `superkingdom_taxid`, `superkingdom`: Superkingdom level taxonomy ID and name

Missing values are used when a taxonomic rank is not available.
"""
function taxids2taxonkit_summarized_lineage_table(taxids::AbstractVector{Int})
    taxid_to_lineage_ranks = taxids2taxonkit_taxid2lineage_ranks(taxids)
    taxids_to_lineage_table = DataFrames.DataFrame()
    for (taxid, lineage_ranks) in taxid_to_lineage_ranks
        # 
        row = (
            taxid = taxid,
            species_taxid = haskey(lineage_ranks, "species") ? lineage_ranks["species"].taxid : missing,
            species = haskey(lineage_ranks, "species") ? lineage_ranks["species"].lineage : missing,
            genus_taxid = haskey(lineage_ranks, "genus") ? lineage_ranks["genus"].taxid : missing,
            genus = haskey(lineage_ranks, "genus") ? lineage_ranks["genus"].lineage : missing,
            family_taxid = haskey(lineage_ranks, "family") ? lineage_ranks["family"].taxid : missing,
            family = haskey(lineage_ranks, "family") ? lineage_ranks["family"].lineage : missing,
            superkingdom_taxid = haskey(lineage_ranks, "superkingdom") ? lineage_ranks["superkingdom"].taxid : missing,
            superkingdom = haskey(lineage_ranks, "superkingdom") ? lineage_ranks["superkingdom"].lineage : missing,
        )
        push!(taxids_to_lineage_table, row, promote=true)
    end
    return taxids_to_lineage_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Lowest Common Ancestor (LCA) taxonomic ID for a set of input taxonomic IDs.

# Arguments
- `ids::Vector{Int}`: Vector of NCBI taxonomic IDs

# Returns
- `Int`: The taxonomic ID of the lowest common ancestor

# Details
Uses taxonkit to compute the LCA. Automatically sets up the required taxonomy database 
if not already present in `~/.taxonkit/`.

# Dependencies
- Requires taxonkit (installed via Bioconda)
- Requires taxonomy database (downloaded automatically if missing)
"""
function taxids2lca(ids::Vector{Int})
    Mycelia.add_bioconda_env("taxonkit")
    if !isdir("$(homedir())/.taxonkit") || isempty(readdir("$(homedir())/.taxonkit"))
        setup_taxonkit_taxonomy()
    end
    # Convert the list of integers to a space-separated string
    input_str = join(ids, " ")

    # Pass the input string to the `taxonkit lca` command and capture the output
    output = read(pipeline(`echo $(input_str)`, `$(Mycelia.CONDA_RUNNER) run --live-stream -n taxonkit taxonkit lca`), String)

    # Split the output string and select the last item
    lca_id = split(chomp(output), "\t")[end]

    # Convert the LCA identifier to an integer and return it
    return parse(Int, lca_id)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a vector of species/taxon names to their corresponding NCBI taxonomy IDs.

# Arguments
- `names::AbstractVector{<:AbstractString}`: Vector of scientific names or common names

# Returns
- `Vector{Int}`: Vector of NCBI taxonomy IDs corresponding to the input names

Progress is displayed using ProgressMeter.
"""
function names2taxids(names::AbstractVector{<:AbstractString})
    results = []
    ProgressMeter.@showprogress for name in names
        push!(results, Mycelia.name2taxid(name))
    end
    return reduce(vcat, results)
end

# """
# Downloads and unpacks the desired .tar.gz prebuilt kraken index

# Go to https://benlangmead.github.io/aws-indexes/k2 and identify the appropriate ".tar.gz" url download
# """
# function download_kraken_index(;url, directory="$(homedir())/workspace/kraken")
#     @assert occursin(r"\.tar\.gz$", url)
#     filename = last(split(url, '/'))
#     output_path = joinpath(directory, filename)
#     # @show url
#     # @show filename
#     # @show output_path
#     mkpath(directory)
#     extracted_directory = replace(basename(output_path), ".tar.gz" => "")
#     if !isdir(extracted_directory)
#         mkpath(extracted_directory)
#     end
#     if isempty(readdir(extracted_directory))
#         # download the file only if needed
#         if !isfile(output_path)
#             download(url, output_path)
#             # run(`wget -P $(directory) $(url)`)
#         end
#         run(`tar -xvzf $(output_path) -C $(extracted_directory)`)
#     end
#     return extracted_directory
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulates genetic variants (substitutions, insertions, deletions, inversions) in a DNA sequence.

# Arguments
- `fasta_record`: Input DNA sequence in FASTA format

# Keywords
- `n_variants=√(sequence_length)`: Number of variants to generate
- `window_size=sequence_length/n_variants`: Size of windows for variant placement
- `variant_size_disbribution=Geometric(1/√window_size)`: Distribution for variant sizes
- `variant_type_likelihoods`: Vector of pairs mapping variant types to probabilities
    - `:substitution => 10⁻¹`
    - `:insertion => 10⁻²` 
    - `:deletion => 10⁻²`
    - `:inversion => 10⁻²`

# Returns
DataFrame in VCF format containing simulated variants with columns:
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE

# Notes
- Variants are distributed across sequence windows to ensure spread
- Variant sizes are capped by window size
- Equivalent variants are filtered out
- FILTER column indicates variant type
"""
function simulate_variants(fasta_record::FASTX.FASTA.Record;        
        n_variants = Int(floor(sqrt(length(FASTX.sequence(fasta_record))))),
        window_size = Int(ceil(length(FASTX.sequence(fasta_record)) / n_variants)),
        variant_size_disbribution = Distributions.Geometric(1/sqrt(window_size)),
        variant_type_likelihoods = [
            :substitution => 10^-1,
            :insertion => 10^-2,
            :deletion => 10^-2,
            :inversion => 10^-2,
            # special case insertion/deletions, skipping
            # :translocations => 10^-3,
            # :duplication => 10^-3,
        ]
    )
    vcf_table = DataFrames.DataFrame(
        "#CHROM" => String[],
        "POS" => Int[],
        "ID" => String[],
        "REF" => String[],
        "ALT" => String[],
        "QUAL" => Int[],
        "FILTER" => String[],
        "INFO" => String[],
        "FORMAT" => String[],
        "SAMPLE" => String[]
    )
    @assert join(names(vcf_table), '\t') == "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE"
    
    original_sequence = BioSequences.LongDNA{4}(FASTX.sequence(fasta_record))
    original_sequence_id = string(first(split(FASTX.description(fasta_record))))
    modified_sequence_id = "modified__" * original_sequence_id
    
    variant_sizes = rand(variant_size_disbribution, n_variants) .+ 1
    @assert all(variant_sizes .>= 1)
    any(variant_sizes .>= window_size) && @warn "variant size distribution is too large, truncating variant sizes to fit in windows"
    variant_sizes = map(x -> min(x, window_size - 1), variant_sizes)
    # display(StatsPlots.plot(variant_size_disbribution, ylabel = "probability density", xlabel = "variant size", title = "Variant Size Distribution", legend=false))
    # display(StatsPlots.histogram(variant_sizes, ylabel="# of variants", xlabel="variant size", title="Actual samples drawn", nbins=length(unique(variant_sizes)), legend=false))
    variant_types = StatsBase.sample(first.(variant_type_likelihoods), StatsBase.weights(last.(variant_type_likelihoods)), n_variants)
    window_starts = 1:window_size:length(original_sequence)
    window_ends = window_size:window_size:length(original_sequence)
    windows = zip(window_starts, window_ends)
    variant_type_sizes = zip(variant_types, variant_sizes)
        
    for ((variant_type, variant_size), (start, stop)) in collect(zip(variant_type_sizes, windows))
        selected_start = rand(start:stop-variant_size)
        @assert selected_start <= stop
        original_subsequence = original_sequence[selected_start:selected_start+variant_size-1]
        @assert length(original_subsequence) == variant_size
        if variant_type == :substitution
            ref = original_subsequence
            alt = BioSequences.randdnaseq(variant_size)
            while alt == ref
                @info "substitution collision, resampling..."
                alt = BioSequences.randdnaseq(variant_size)
            end
        elseif variant_type == :insertion
            ref = original_subsequence
            alt = original_subsequence * BioSequences.randdnaseq(variant_size)
        elseif variant_type == :deletion
            ref = original_sequence[selected_start:selected_start+variant_size]
            alt = original_sequence[selected_start:selected_start]
        elseif variant_type == :inversion
            ref = original_subsequence
            alt = BioSequences.reverse_complement(original_subsequence)
        end 
        # @show selected_start
        row = Dict(
            "#CHROM" => original_sequence_id,
            "POS" => selected_start,
            "ID" => ".",
            "REF" => string(ref),
            "ALT" => string(alt),
            "QUAL" => 60,
            "FILTER" => string(variant_type),
            "INFO" => ".",
            "FORMAT" => "GT:GQ",
            "SAMPLE" => "1:60"
        )
        push!(vcf_table, row)
        original_prefix = original_sequence[start:selected_start-1]
        if variant_type == :deletion
            original_suffix = original_sequence[selected_start+variant_size+1:stop]
        else
            original_suffix = original_sequence[selected_start+variant_size:stop]
        end
        reconstructed_window = original_prefix * ref * original_suffix
        original_window = original_sequence[start:stop]
        @assert reconstructed_window == original_window
    end
    true_variant = vcf_table[!, "REF"] .!= vcf_table[!, "ALT"]
    if !all(true_variant)
        @warn "filtering equivalent variants"
        vcf_table = vcf_table[true_variant, :]
    end
    return vcf_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Simulates genetic variants from sequences in a FASTA file and generates corresponding VCF records.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing sequences to analyze

# Details
1. Processes each record in the input FASTA file
2. Generates simulated variants for each sequence
3. Creates a VCF file with the same base name as input file (.vcf extension)
4. Updates sequences with simulated variants in a new FASTA file (.vcf.fna extension)

# Returns
Path to the modified FASTA file containing sequences with simulated variants
"""
function simulate_variants(fasta_file::String)
    vcf_file = fasta_file * ".vcf"
    modified_fasta_file = vcf_file * ".fna"
    vcf_table = DataFrames.DataFrame()
    for fasta_record in open_fastx(fasta_file)
        append!(vcf_table, simulate_variants(fasta_record))
    end
    Mycelia.write_vcf_table(vcf_file=vcf_file, vcf_table=vcf_table, fasta_file=fasta_file)
    return Mycelia.update_fasta_with_vcf(in_fasta = fasta_file, vcf_file = vcf_file)
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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function initialize_transition_probabilities(kmer_graph)
    
#     total_kmers = Graphs.nv(kmer_graph)
#     transition_likelihoods = Dict(
#         true => SparseArrays.spzeros(total_kmers, total_kmers),
#         false => SparseArrays.spzeros(total_kmers, total_kmers)
#     )

#     for edge in collect(Graphs.edges(kmer_graph))
# #         weight = length(kmer_graph.eprops[edge][:evidence])
#         weight = kmer_graph.eprops[edge][:weight]
#         for o in kmer_graph.eprops[edge][:orientations]
#             transition_likelihoods[o.source_orientation][edge.src, edge.dst] = weight
#         end
#     end

#     for source_orientation in (true, false)
#         for src in 1:total_kmers
#             transition_weights = transition_likelihoods[source_orientation][src, :]
#             total_weight = sum(transition_weights)
#             dsts, vals = SparseArrays.findnz(transition_weights)
#             for (dst, val) in zip(dsts, vals) 
#                 transition_likelihoods[source_orientation][src, dst] = val / total_weight
#             end
#             normalized_probability = sum(transition_likelihoods[source_orientation][src, :])
#             @assert isapprox(normalized_probability, 0) || isapprox(normalized_probability, 1)
#         end
#     end
#     return transition_likelihoods
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function set_initial_state_likelihoods!(
#         kmer_graph,
#         initial_state,
#         kmer_likelihoods,
#         error_rate,
#         state_likelihoods,
#         arrival_paths
#     )
#     for vertex in collect(Graphs.vertices(kmer_graph))
#         hidden_kmer = kmer_graph.vprops[vertex][:kmer]

#         fw_alignment = 
#             BioAlignments.pairalign(
#                 BioAlignments.LevenshteinDistance(), 
#                 initial_state.fw, 
#                 hidden_kmer)

#         fw_probability = kmer_likelihoods[vertex]

#         for match in 1:BioAlignments.count_matches(BioAlignments.alignment(fw_alignment))
#             fw_probability *= 1 - error_rate
#         end

#         for edit in 1:fw_alignment.value
#             fw_probability *= error_rate
#         end

#         bw_alignment = 
#             BioAlignments.pairalign(
#                 BioAlignments.LevenshteinDistance(),
#                 initial_state.bw,
#                 hidden_kmer)

#         bw_probability = kmer_likelihoods[vertex]

#         for match in 1:BioAlignments.count_matches(BioAlignments.alignment(bw_alignment))
#             bw_probability *= 1 - error_rate
#         end

#         for edit in 1:bw_alignment.value
#             bw_probability *= error_rate
#         end

#         if fw_probability > bw_probability
#             state_probability = fw_probability
#             state_orientation = true
#         elseif fw_probability < bw_probability
#             state_probability = bw_probability
#             state_orientation = false
#         else fw_probability == bw_probability
#             state_probability = fw_probability
#             state_orientation = missing
#         end
#         state_likelihoods[vertex, 1] = state_probability
#         arrival_paths[vertex, 1] = [vertex => state_orientation]
#     end
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_viterbi!(
#         current_state,
#         prior_state,
#         observed_nucleotide,
#         observed_quality_score,
#         observed_error_rate,
#         current_vertex,
#         prior_vertex,
#         state_likelihoods,
#         transition_likelihoods,
#         shortest_paths,
#         arrival_paths,
#         kmer_graph,
#         kmer_likelihoods
#         )
#     # if probability of prior state is lower than current probability, skip

# #     @show current_state
# #     @show prior_state
# #     @show current_vertex
# #     @show prior_vertex
    
    
#     current_state_likelihood = state_likelihoods[current_vertex, current_state]
#     prior_state_likelihood = state_likelihoods[prior_vertex, prior_state]

#     # if we already have a better possible path, skip calculating anything
#     if prior_state_likelihood < current_state_likelihood
# #         @show prior_state_likelihood < current_state_likelihood
#         return
#     end

#     # take shortest path and assume it's the maximum likelihood path
#     # this assumption seems fair because in an ideal situation
#     # we're just moving to an adjacent kmer
#     # and the shortest path and most likely path should be the same
#     shortest_path = shortest_paths[prior_vertex][current_vertex]
    
# #     no path & not considering insertion
#     if isempty(shortest_path) && (prior_vertex != current_vertex)
# #         @show "no path, skipping"
#         return
#     end
    
#     # if shortest path isn't viable, exit
#     if !isempty(shortest_path)
# #         @show "checking if path is viable"

#         terminal_orientation_prior_state = last(last(arrival_paths[prior_vertex, prior_state]))
# #         @show arrival_paths[prior_vertex, prior_state]
# #         @show "we were at vertex $(prior_vertex) in orientation $(terminal_orientation_prior_state)"
#         candidate_edge = Graphs.Edge(shortest_path[1], shortest_path[2])
                
#         if !ismissing(terminal_orientation_prior_state) && 
#             !any(o -> o.source_orientation == terminal_orientation_prior_state, kmer_graph.eprops[candidate_edge][:orientations])
            
# #             @show "no viable orientation matching edges detected between $(candidate_edge)"
# #             @show "full candidate path was $(shortest_path)"
# #             @show "orientation options were:"
# #             @show kmer_graph.eprops[candidate_edge][:orientations]
#             return
#         end
#     end
    
#     # zero step path - insertion in observed sequence relative to kmer graph
#     is_same_vertex = (current_vertex == prior_vertex)
#     has_edge = Graphs.has_edge(kmer_graph, Graphs.Edge(prior_vertex, current_vertex))
#     if is_same_vertex && has_edge
#         shortest_path = [prior_vertex, current_vertex]
#     end
    
#     if is_same_vertex
# #         @show "same vertex, considering insertion potential"
#         emission_likelihood = observed_error_rate
#         transition_likelihood = observed_error_rate
#         state_likelihood = kmer_likelihoods[current_vertex]
#         path_likelihood = prior_state_likelihood * emission_likelihood * transition_likelihood * state_likelihood
#         path = [last(arrival_paths[prior_vertex, prior_state])]

#         if current_state_likelihood > state_likelihoods[current_vertex, current_state]
# #             @show "selecting path"
# #             @show path
# #             @show path_likelihood
#             state_likelihoods[current_vertex, current_state] = path_likelihood
#             arrival_paths[current_vertex, current_state] = path
#         end
#     # one or more step path - match, mismatch, or deletion in observed sequence relative to kmer graph
#     elseif !isempty(shortest_path)
# #         @show "path is viable!"
# #         @show "considering shortest path: $(shortest_path)"

#         initial_path_state = last(arrival_paths[prior_vertex, prior_state])

#         path = Vector{typeof(initial_path_state)}(undef, length(shortest_path))
#         path[1] = initial_path_state

#         path_likelihood::Float64 = state_likelihoods[prior_vertex, prior_state]

#         for i in 2:length(shortest_path)

#             this_vertex = shortest_path[i]
#             prior_vertex, prior_orientation = path[i-1]
#             edge = Graphs.Edge(prior_vertex, this_vertex)

#             possible_edge_orientations::Set{NamedTuple{(:source_orientation, :destination_orientation), Tuple{Bool, Bool}}} = kmer_graph.eprops[edge][:orientations]
            
# #             @show possible_edge_orientations
            
#             if !ismissing(prior_orientation)
#                 possible_edge_orientations = filter(o -> o.source_orientation == prior_orientation, possible_edge_orientations)
#             end
            
# #             @show possible_edge_orientations
            
#             if isempty(possible_edge_orientations)
#                 path_likelihood *= 0.0
#                 path = Vector{eltype(path)}()
# #                 @show "no possible orientations, bailing early"
#                 break
#             end

# #             @show prior_orientation
#             if ismissing(prior_orientation)
#                 if transition_likelihoods[true][prior_vertex, this_vertex] > transition_likelihoods[false][prior_vertex, this_vertex]
#                     prior_orientation = true
#                     transition_likelihood = transition_likelihoods[true][prior_vertex, this_vertex]::Float64
#                 elseif transition_likelihoods[true][prior_vertex, this_vertex] < transition_likelihoods[false][prior_vertex, this_vertex]
#                     prior_orientation = false
#                     transition_likelihood = transition_likelihoods[false][prior_vertex, this_vertex]::Float64
#                 else transition_likelihoods[true][prior_vertex, this_vertex] == transition_likelihoods[false][prior_vertex, this_vertex]
#                     prior_orientation = missing
#                     transition_likelihood = transition_likelihoods[true][prior_vertex, this_vertex]::Float64
#                 end
#             else
#                 transition_likelihood = transition_likelihoods[prior_orientation][prior_vertex, this_vertex]::Float64
#             end
#             state_likelihood::Float64 = kmer_likelihoods[this_vertex]
#             path_likelihood *= transition_likelihood * state_likelihood
            
#             if length(possible_edge_orientations) == 1
#                 orientation = first(possible_edge_orientations).destination_orientation
#                 path[i] = this_vertex => orientation
#             else
#                 path[i] = this_vertex => missing
#             end
#         end

#         # see if new nucleotide is a match or mismatch to terminal kmer in path
#         if !isempty(path) && path_likelihood > 0
#             terminal_kmer_index, terminal_kmer_orientation = last(path)
#             terminal_kmer = BioSequences.LongDNASeq(kmer_graph.vprops[terminal_kmer_index][:kmer])::BioSequences.LongDNASeq
#             if ismissing(terminal_kmer_orientation)
#                 fw_is_match = observed_nucleotide == last(terminal_kmer)
#                 bw_is_match = observed_nucleotide == last(BioSequences.reverse_complement!(terminal_kmer))
#                 if fw_ismatch && !bw_is_match
#                     path[end] = terminal_kmer_index => true
#                     path_likelihood *= 1 - observed_error_rate
#                 elseif !fw_ismatch && bw_is_match
#                     path[end] = terminal_kmer_index => false
#                     path_likelihood *= 1 - observed_error_rate
#                 elseif fw_ismatch && bw_is_match
#                     path_likelihood *= 1 - observed_error_rate
#                 elseif !fw_ismatch && !bw_is_match
#                     path_likelihood *= observed_error_rate
#                 end
#             elseif terminal_kmer_orientation
#                 is_match = observed_nucleotide == last(terminal_kmer)
#                 if is_match
#                     path_likelihood *= 1 - observed_error_rate
#                 else
#                     path_likelihood *= observed_error_rate
#                 end
#             else
#                 terminal_kmer = BioSequences.reverse_complement!(terminal_kmer)
#                 is_match = observed_nucleotide == last(terminal_kmer)
#                 if is_match
#                     path_likelihood *= 1 - observed_error_rate
#                 else
#                     path_likelihood *= observed_error_rate
#                 end
#             end
#         end

#         if path_likelihood > state_likelihoods[current_vertex, current_state]
# #             @show "selecting path"
# #             @show path
# #             @show path_likelihood
#             state_likelihoods[current_vertex, current_state] = path_likelihood
#             arrival_paths[current_vertex, current_state] = path
#         end
#     end
#     return
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function polish_fastq(kmer_graph, fastq_file)

# #     @info "Assessing kmer likelihoods"
#     kmers = [kmer_graph.vprops[v][:kmer] for v in Graphs.vertices(kmer_graph)]
# #     kmer_counts = [length(kmer_graph.vprops[v][:evidence]) for v in Graphs.vertices(kmer_graph)]
#     kmer_counts = [kmer_graph.vprops[v][:weight] for v in Graphs.vertices(kmer_graph)]
#     kmer_likelihoods = kmer_counts ./ sum(kmer_counts)
#     k = kmer_graph.gprops[:k]
#     kmer_type = BioSequences.BigDNAMer{k}
#     total_kmers = length(kmers)
    
# #     @info "determining shortest paths between kmers"
#     shortest_paths = Graphs.enumerate_paths(Graphs.floyd_warshall_shortest_paths(kmer_graph));
    
#     @info "counting the number of records to establish runtime estimate"
#     number_of_records = 0
#     for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
#         number_of_records += 1
#     end
#     progress_bar = ProgressMeter.Progress(number_of_records, 1)
    
#     output_fastq_file = replace(fastq_file, ".fastq" => ".k$(kmer_graph.gprops[:k]).fastq")
#     fastq_writer = FASTX.FASTQ.Writer(open(output_fastq_file, "w"))
#     for fastq_record in FASTX.FASTQ.Reader(open(fastq_file))
#         ProgressMeter.next!(progress_bar)
        
# #         @info "Initializing matrices"
#         total_states = length(FASTX.sequence(fastq_record))-k+1
#         transition_likelihoods = initialize_transition_probabilities(kmer_graph)
#         state_likelihoods = zeros(total_kmers, total_states)
#         arrival_paths = fill(Pair{Int, Union{Bool, Missing}}[], total_kmers, total_states)

# #         @info "Determining Likelihoods of initial states"
#         initial_state = first(BioSequences.each(kmer_type, FASTX.sequence(fastq_record)))
#         current_state = 1
#         # note this is a place for potential improvement, use the q value at each base to guide probability rather than median
#         median_q_value = Statistics.median(Int.(FASTX.quality(fastq_record)[1:k]))
#         current_error_rate = q_value_to_error_rate(median_q_value)
#         # canonical_kmer = BioSequences.canonical(initial_state.fw)
#         set_initial_state_likelihoods!(
#                 kmer_graph,
#                 initial_state,
#                 kmer_likelihoods,
#                 current_error_rate,
#                 state_likelihoods,
#                 arrival_paths
#             )

# #         @info "Determining likelihood of downstream states"

# #         non_singleton_states = findall(kmer_counts .> 1)

#         for current_state in 2:total_states
#             prior_state = current_state - 1

#         #     observed_kmer = BioSequences.BigDNAMer{k}(FASTX.sequence(fastq_record)[current_state:current_state+k-1])

#         #     @assert observed_kmer == collect(BioSequences.each(kmer_type, FASTX.sequence(fastq_record)))[current_state].fw

#         #     canonical_kmer = BioSequences.canonical(observed_kmer)

#             observed_nucleotide = FASTX.sequence(fastq_record)[k-1+current_state]
#         #     observed_nucleotide = last(observed_kmer)
#             observed_quality_score = FASTX.quality(fastq_record)[k-1+current_state]
#             observed_error_rate = q_value_to_error_rate(observed_quality_score)

#             # we'll assess prior states in order of decreasing likelihood
#             # such that we maximize how frequently we are able to utilize the
#             # current_state_likelihood > candidate prior state
#             # break that won't bother evaluating lower likelihood possibilities
#             prior_states_in_decreasing_likelihood = sortperm(state_likelihoods[:, prior_state], rev=true)

#             # and skip all prior states with zero probability

#             for current_vertex in total_states
#                 for prior_vertex in prior_states_in_decreasing_likelihood
#                     if state_likelihoods[prior_vertex, prior_state] > 0
#                         run_viterbi!(
#                                 current_state,
#                                 prior_state,
#                                 observed_nucleotide,
#                                 observed_quality_score,
#                                 observed_error_rate,
#                                 current_vertex,
#                                 prior_vertex,
#                                 state_likelihoods,
#                                 transition_likelihoods,
#                                 shortest_paths,
#                                 arrival_paths,
#                                 kmer_graph,
#                                 kmer_likelihoods
#                                 )
#                     end
#                 end
#             end
#         end

# #         try
#         maximum_likelihood_path, maximum_likelihood_value = 
#             determine_maximum_likelihood_path(
#                 state_likelihoods,
#                 arrival_paths
#                 )
# #         catch
# #             return state_likelihoods, arrival_paths
# #         end

#         sequence = oriented_path_to_sequence(kmer_graph, maximum_likelihood_path)

# #         @info "comparing to original path"
#         original_sequence_likelihood = oriented_path_to_likelihood(kmer_graph, kmers, kmer_likelihoods, transition_likelihoods, fastq_record)
#         relative_likelihood = maximum_likelihood_value / original_sequence_likelihood
# #         relative_likelihood_formatted = NumericIO.formatted(relative_likelihood, ndigits=1, charset=:ASCII)
# #         println("relative likelihood of new path to old path is $(relative_likelihood_formatted)")

# #         @info "writing updated record"
#         identifier = FASTX.identifier(fastq_record) * "_k$(k)"
#         description = string(relative_likelihood)
#         # because the sequences won't always be the same length, we take an ordered sampling with replacement
#         # which introduces some random error but preserves overall patterns and areas of high/low accuracy
#         quality_scores = StatsBase.sample(FASTX.quality(fastq_record), length(sequence), ordered=true)

#         new_fastq_record = FASTX.FASTQ.Record(
#             identifier,
#             description,
#             sequence,
#             quality_scores
#         )
#         write(fastq_writer, new_fastq_record)
#     end
#     close(fastq_writer)
#     return output_fastq_file
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_maximum_likelihood_path(
#     state_likelihoods,
#     arrival_paths
#     )
#     maximum_likelihood_value = maximum(state_likelihoods[:, end])

#     maximum_likelihood_path_indices = findall(state_likelihoods[:, end] .== maximum_likelihood_value)

#     # if multiple paths are tied, randomly choose one
#     maximum_likelihood_path_index = rand(maximum_likelihood_path_indices)

#     maximum_likelihood_path = arrival_paths[maximum_likelihood_path_index, end]

#     for state_index in size(state_likelihoods, 2)-1:-1:1
#         next_kmer, next_orientation = first(maximum_likelihood_path)
#         maximum_likelihood_arrival_path = arrival_paths[next_kmer, state_index]
        
#         is_match = last(maximum_likelihood_arrival_path) == (next_kmer => next_orientation)
#         if !ismissing(is_match) && !is_match
#             error("breaking")
#         end
#         maximum_likelihood_path = vcat(maximum_likelihood_arrival_path[1:end-1], maximum_likelihood_path)
#     end
#     return maximum_likelihood_path, maximum_likelihood_value
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a DNA sequence into a path through a collection of stranded k-mers.

# Arguments
- `stranded_kmers`: Collection of unique k-mers representing possible path vertices
- `sequence`: Input DNA sequence to convert to a path

# Returns
Vector of `Pair{Int,Bool}` where:
- First element (Int) is the index of the k-mer in `stranded_kmers`
- Second element (Bool) indicates orientation (true=forward, false=reverse)
"""
function sequence_to_stranded_path(stranded_kmers, sequence)
    KMER_TYPE = typeof(first(stranded_kmers))
    path = Vector{Pair{Int, Bool}}()
    for (i, kmer) in Kmers.EveryKmer{KMER_TYPE}(BioSequences.LongDNA{4}(sequence))
        kmer_index = findfirst(stranded_kmer -> kmer == stranded_kmer, stranded_kmers)
        orientation = true
        push!(path, kmer_index => orientation)
    end
    return path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a path through k-mers into a single DNA sequence.

Takes a vector of k-mers and a path representing the order to traverse them,
reconstructs the original sequence by joining the k-mers according to the path.
The first k-mer is used in full, then only the last nucleotide from each subsequent k-mer is added.

# Arguments
- `kmers`: Vector of DNA k-mers (as LongDNA{4})
- `path`: Vector of tuples representing the path through the k-mers

# Returns
- `LongDNA{4}`: The reconstructed DNA sequence
"""
function path_to_sequence(kmers, path)
    # @show path
    sequence = BioSequences.LongDNA{4}(kmers[first(first(path))])
    for i in 2:length(path)
        push!(sequence, kmers[first(path[i])][end])
    end
    return sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the probability of traversing a specific edge in a stranded k-mer graph.

The probability is computed as the ratio of this edge's coverage weight to the sum
of all outgoing edge weights from the source vertex.

$(DocStringExtensions.TYPEDSIGNATURES)

# Arguments
- `stranded_kmer_graph`: A directed graph where edges represent k-mer connections
- `edge`: The edge for which to calculate the probability

# Returns
- `Float64`: Probability in range [0,1] representing likelihood of traversing this edge
  Returns 0.0 if sum of all outgoing edge weights is zero

# Note
Probability is based on the :coverage property of edges, using their length as weights
"""
function edge_probability(stranded_kmer_graph, edge)
    neighbors = Graphs.outneighbors(stranded_kmer_graph, edge.src)
    neighbor_under_consideration = findfirst(neighbor -> neighbor == edge.dst, neighbors)
    edge_weights = [length(stranded_kmer_graph.eprops[Graphs.Edge(edge.src, neighbor)][:coverage]) for neighbor in neighbors]
    if sum(edge_weights) == 0
        p = 0.0
    else
        edge_probabilities = edge_weights ./ sum(edge_weights)
        p = edge_probabilities[neighbor_under_consideration]
    end
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Finds maximum likelihood paths through a stranded k-mer graph using the Viterbi algorithm
to correct sequencing errors.

# Arguments
- `stranded_kmer_graph`: A directed graph where vertices represent k-mers and edges represent overlaps
- `error_rate::Float64`: Expected per-base error rate (default: 1/(k+1)). Must be < 0.5
- `verbosity::String`: Output detail level ("debug", "reads", or "dataset")

# Returns
Vector of FASTX.FASTA.Record containing error-corrected sequences

# Details
- Uses dynamic programming to find most likely path through k-mer graph
- Accounts for matches, mismatches, insertions and deletions
- State likelihoods based on k-mer coverage counts
- Transition probabilities derived from error rate
- Progress tracking based on verbosity level

# Notes
- Error rate should be probability of error (e.g. 0.01 for 1%), not accuracy
- Higher verbosity levels ("debug", "reads") provide detailed path finding information
- "dataset" verbosity shows only summary statistics
"""
function viterbi_maximum_likelihood_traversals(stranded_kmer_graph;
                                               error_rate::Float64=1/(stranded_kmer_graph.gprops[:k] + 1),
                                               verbosity::String="dataset")
    @assert verbosity in ["debug", "reads", "dataset"]
    if error_rate >= .5
        error("Error rate >= 50%. Did you enter the accuracy by mistake?")
    end

    if verbosity in ["debug", "reads", "dataset"]
        println("computing kmer counts...")
    end
    stranded_kmer_counts = [length(stranded_kmer_graph.vprops[vertex][:coverage]) for vertex in Graphs.vertices(stranded_kmer_graph)]
    if verbosity in ["debug", "reads", "dataset"]
        println("computing kmer state likelihoods...")
    end
    stranded_kmer_likelihoods = stranded_kmer_counts ./ sum(stranded_kmer_counts)
    accuracy = 1 - error_rate

    if verbosity in ["debug"]
        println("STATE LIKELIHOODS:")
        println("\tkmer\tcount\tlikelihood")
        for vertex in Graphs.vertices(stranded_kmer_graph)
            kmer = stranded_kmer_graph.gprops[:stranded_kmers][vertex]
            count = stranded_kmer_counts[vertex]
            likelihood = stranded_kmer_likelihoods[vertex]
            println("\t$kmer\t$count\t$likelihood")
        end
    end
    if verbosity in ["debug", "reads", "dataset"]
        println("finding shortest paths between kmers...")
    end
    shortest_paths = Graphs.enumerate_paths(Graphs.floyd_warshall_shortest_paths(stranded_kmer_graph))
    K = stranded_kmer_graph.gprops[:K]
    for K1 in 1:K
        for K2 in 1:K
            if K1 != K2
                shortest_path = shortest_paths[K1][K2]
                path_likelihood = 1.0
                for ui in 1:length(shortest_path)-1
                    u = shortest_path[ui]
                    v = shortest_path[ui + 1]
                    # likelihood of the transition
                    path_likelihood *= edge_probability(stranded_kmer_graph, Graphs.Edge(u, v))
                end
                if path_likelihood == 0.0
                    shortest_paths[K1][K2] = Vector{Int}()
                end
            elseif K1 == K2
                # the shortest path from a kmer to itself is an insertion (no edge)
                # so need to manually check for self loops
                if Graphs.has_edge(stranded_kmer_graph, Graphs.Edge(K1, K2))
                    if edge_probability(stranded_kmer_graph, Graphs.Edge(K1, K2)) != 0.0
                        shortest_paths[K1][K2] = [K1, K2]
                    else
                        shortest_paths[K1][K2] = Vector{Int}()
                    end
                # otherwise, check to see if any outneighbors connect back to the kmer
                else
                    connected_outneighbors = filter(outneighbor -> Graphs.has_path(stranded_kmer_graph, outneighbor, K2), Graphs.outneighbors(stranded_kmer_graph, K1))
                    if !isempty(connected_outneighbors)
                        outneighbor_cycles = [[K1, shortest_paths[outneighbor][K2]...] for outneighbor in connected_outneighbors]
                        cycle_likelihoods = ones(length(outneighbor_cycles))
                        for (i, cycle) in enumerate(outneighbor_cycles)
                            for ui in 1:length(cycle)-1
                                u = cycle[ui]
                                v = cycle[ui + 1]
                                # likelihood of the transition
                                cycle_likelihoods[i] *= edge_probability(stranded_kmer_graph, Graphs.Edge(u, v))
                            end
                            # include likelihoods of states
                            for vertex in cycle[2:end-1]
                                cycle_likelihoods[i] *= stranded_kmer_likelihoods[vertex]
                            end
                        end
                        path_likelihood = maximum(cycle_likelihoods)
                        max_likelihood_cycle_indices = findall(cycle_likelihoods .== path_likelihood)
                        shortest_paths[K1][K2] = outneighbor_cycles[first(max_likelihood_cycle_indices)]
                    else
                        shortest_paths[K1][K2] = Vector{Int}()
                    end
                end
            end
            if length(shortest_paths[K1][K2]) == 1
                shortest_paths[K1][K2] = Vector{Int}()
            end
        end
    end

    if verbosity in ["debug"]
        for K1 in 1:K
            for K2 in 1:K
                println("\t$K1\t$K2\t$(shortest_paths[K1][K2])")
            end
        end
    end

    total_bases_observed = 0
    total_edits_accepted = 0

    corrected_observations = FASTX.FASTA.Record[]
    if verbosity in ["debug", "reads", "dataset"]
        println("finding viterbi maximum likelihood paths for observed sequences...")
    end
    # p = Progress(length(stranded_kmer_graph.gprops[:observed_paths]))
    for (observation_index, observed_path) in enumerate(stranded_kmer_graph.gprops[:observed_paths])
        if verbosity in ["debug", "reads"]
            println("\nevaluating sequence $observation_index of $(length(stranded_kmer_graph.gprops[:observed_paths]))")
        end
        # consider switching to log transform
        kmer_likelihoods = zeros(Graphs.nv(stranded_kmer_graph), length(observed_path))
        kmer_arrival_paths = Array{Vector{Int}}(undef, Graphs.nv(stranded_kmer_graph), length(observed_path))
        edit_distances = zeros(Int, Graphs.nv(stranded_kmer_graph), length(observed_path))
        # changed here!!
        observed_kmer_index, observed_kmer_orientation = observed_path[1]
        
        observed_kmer_sequence = stranded_kmer_graph.gprops[:stranded_kmers][observed_kmer_index]
        for hidden_kmer_index in Graphs.vertices(stranded_kmer_graph)
            hidden_kmer_sequence = stranded_kmer_graph.gprops[:stranded_kmers][hidden_kmer_index]
            alignment_result = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), observed_kmer_sequence, hidden_kmer_sequence)
            number_of_matches = BioAlignments.count_matches(BioAlignments.alignment(alignment_result))
            number_of_edits = stranded_kmer_graph.gprops[:k] - number_of_matches
            kmer_likelihoods[hidden_kmer_index, 1] = stranded_kmer_likelihoods[hidden_kmer_index]
            for match in 1:number_of_matches
                kmer_likelihoods[hidden_kmer_index, 1] *= accuracy
            end
            for edit in 1:number_of_edits
                kmer_likelihoods[hidden_kmer_index, 1] *= error_rate
            end
            kmer_arrival_paths[hidden_kmer_index, 1] = Vector{Int}()
            edit_distances[hidden_kmer_index, 1] = number_of_edits
        end
        kmer_likelihoods[:, 1] ./= sum(kmer_likelihoods[:, 1])
        # from here on, all probabilities are log transformed
        kmer_likelihoods[:, 1] .= log.(kmer_likelihoods[:, 1])
        if verbosity in ["debug"]
            println("\tconsidering path state 1")
            println("\t\tobserved kmer $observed_kmer_sequence")
            println("\t\tInitial state log likelihoods:")
            for line in split(repr(MIME("text/plain"), kmer_likelihoods[:, 1]), '\n')
                println("\t\t\t$line")
            end
        end
        for observed_path_index in 2:length(observed_path)
            # changed!!
            observed_kmer_index, observed_kmer_orientation = observed_path[observed_path_index]
            observed_base = stranded_kmer_graph.gprops[:stranded_kmers][observed_kmer_index][end]

            if verbosity in ["debug"]
                println("\tconsidering path state $observed_path_index")
                println("\t\tobserved base $observed_base")
            end

            MATCH = 1
            MISMATCH = 2
            DELETION = 3
            INSERTION = 4
            arrival_likelihoods = ones(K, 4)
            arrival_paths = fill(Vector{Int}(), K, 4)

            for K2 in 1:K
                kmer_base = stranded_kmer_graph.gprops[:stranded_kmers][K2][end]
                base_is_match = kmer_base == observed_base

                maximum_likelihood = log(0.0)
                maximum_likelihood_path = Vector{Int}()
                maximum_likelihood_edit_distance = 0

                for K1 in 1:K
                    shortest_path = shortest_paths[K1][K2]
                    if length(shortest_path) >= 2
                        edit_distance = Int(!base_is_match) + length(shortest_path) - 2
                        if edit_distance == 0
                            p = kmer_likelihoods[K1, observed_path_index-1] +
                                log(accuracy) + log(stranded_kmer_likelihoods[K2])
                        else
                            p = kmer_likelihoods[K1, observed_path_index-1] +
                                log(error_rate^edit_distance) + log(stranded_kmer_likelihoods[K2])
                        end
                        edit_distance += edit_distances[K1, observed_path_index-1]
                    else
                        p = log(0.0)
                    end
                    if K1 == K2 # consider insertion
                        # in theory, I don't think we should care if the base
                        # matches or not because it's an inserted & erroneous
                        # base, but in practice it's necessary to balance
                        # insertion probabilities with deletion probabilities
                        insertion_p = kmer_likelihoods[K1, observed_path_index-1] +
                                      log(error_rate^(1 + Int(!base_is_match))) + log(stranded_kmer_likelihoods[K2])
                        if insertion_p > p
                            p = insertion_p
                            edit_distance = edit_distances[K1, observed_path_index-1] + 1
                            shortest_path = [K2]
                        end
                    end
                    if p > maximum_likelihood
                        maximum_likelihood = p
                        maximum_likelihood_path = shortest_path
                        maximum_likelihood_edit_distance = edit_distance
                    end
                end
                kmer_likelihoods[K2, observed_path_index] = maximum_likelihood
                kmer_arrival_paths[K2, observed_path_index] = maximum_likelihood_path
                edit_distances[K2, observed_path_index] = maximum_likelihood_edit_distance
            end

            if verbosity in ["debug"]
                println("\t\tkmer log likelihoods")
                for line in split(repr(MIME("text/plain"), kmer_likelihoods), '\n')
                    println("\t\t\t$line")
                end
                println("\t\tarrival paths")
                for line in split(repr(MIME("text/plain"), kmer_arrival_paths), '\n')
                    println("\t\t\t$line")
                end
            end
        end

        if verbosity in ["debug"]
            println("\n\tInputs for viterbi maximum likelihood traversal evaluation:")
            println("\t\tkmer log likelihoods")
            for line in split(repr(MIME("text/plain"), kmer_likelihoods), '\n')
                println("\t\t\t$line")
            end
            println("\t\tkmer arrival paths")
            for line in split(repr(MIME("text/plain"), kmer_arrival_paths), '\n')
                println("\t\t\t$line")
            end
            println("\t\tedit distances")
            for line in split(repr(MIME("text/plain"), edit_distances), '\n')
                println("\t\t\t$line")
            end
        end

        ## backtrack
        maximum_likelihood_path_value = maximum(kmer_likelihoods[:, end])
        maximum_likelihood_path_indices = findall(kmer_likelihoods[:, end] .== maximum_likelihood_path_value)
        # if multiple paths are tied, randomly choose one
        maximum_likelihood_path_index = rand(maximum_likelihood_path_indices)
        maximum_likelihood_edit_distance = edit_distances[maximum_likelihood_path_index, end]

        if length(kmer_arrival_paths[maximum_likelihood_path_index, end]) > 0
            maximum_likelihood_path = last(kmer_arrival_paths[maximum_likelihood_path_index, end])
            for observed_path_index in length(observed_path):-1:1
                maximum_likelihood_arrival_path = kmer_arrival_paths[maximum_likelihood_path_index, observed_path_index]
                maximum_likelihood_path = vcat(maximum_likelihood_arrival_path[1:end-1], maximum_likelihood_path)
                maximum_likelihood_path_index = first(maximum_likelihood_path)
            end
        else
            maximum_likelihood_path = [maximum_likelihood_path_index]
        end
        observed_sequence = path_to_sequence(stranded_kmer_graph.gprops[:stranded_kmers], observed_path)
        maximum_likelihood_sequence = path_to_sequence(stranded_kmer_graph.gprops[:stranded_kmers], maximum_likelihood_path)
        if verbosity in ["debug", "reads"]
            println("\tobserved sequence                 $observed_sequence")
            println("\tmaximum likelihood sequence       $maximum_likelihood_sequence")
            println("\tmaximum likelihood edit distance  $maximum_likelihood_edit_distance")
        end
        total_bases_observed += length(observed_sequence)
        total_edits_accepted += maximum_likelihood_edit_distance
        id = stranded_kmer_graph.gprops[:observation_ids][observation_index]
        kmer_stamped_id = id * "_" * string(stranded_kmer_graph.gprops[:k])
        push!(corrected_observations, FASTX.FASTA.Record(kmer_stamped_id, maximum_likelihood_sequence))
        # progress meter
        # next!(p)
    end
    if verbosity in ["debug", "reads", "dataset"]
        println("\nDATASET STATISTICS:")
        println("\tassumed error rate    $(error_rate * 100)%")
        println("\ttotal bases observed  $total_bases_observed")
        println("\ttotal edits accepted  $total_edits_accepted")
        println("\tinferred error rate   $((total_edits_accepted/total_bases_observed) * 100)%")
    end
    return corrected_observations
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a weighted, strand-specific kmer (de bruijn) graph from a set of kmers
and a series of sequence observations in FASTA format.
"""
function build_stranded_kmer_graph(kmer_type, observations::AbstractVector{<:Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}})
    
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, observations)
    canonical_kmers = collect(keys(canonical_kmer_counts))
    # @show canonical_kmers
    
    # if isempty(canonical_kmers)
    #     @error "isempty(canonical_kmers) = $(isempty(canonical_kmers))"
    # elseif isempty(observations)
    #     @error "isempty(observations) = $(isempty(observations))"
    # end
    stranded_kmers = sort!(vcat(canonical_kmers, [BioSequences.reverse_complement(kmer) for kmer in canonical_kmers]))
    stranded_kmer_to_reverse_complement_map = [
        findfirst(stranded_kmer -> BioSequences.reverse_complement(stranded_kmer) == kmer, stranded_kmers) for kmer in stranded_kmers
    ]
    stranded_kmer_graph = MetaGraphs.MetaDiGraph(length(stranded_kmers))
    stranded_kmer_graph.gprops[:stranded_kmers] = stranded_kmers
    stranded_kmer_graph.gprops[:reverse_complement_map] = stranded_kmer_to_reverse_complement_map
    stranded_kmer_graph.gprops[:k] = length(first(stranded_kmers))
    stranded_kmer_graph.gprops[:K] = length(stranded_kmers)
    stranded_kmer_graph.gprops[:observation_color_map] = Vector{Int}()
    stranded_kmer_graph.gprops[:observation_ids] = Vector{String}()
    stranded_kmer_graph.gprops[:observed_paths] = Vector{Vector{Pair{Int, Bool}}}()
    for vertex in 1:Graphs.nv(stranded_kmer_graph)
        stranded_kmer_graph.vprops[vertex] = Dict(:coverage => Vector{Pair{Int, Pair{Int, Bool}}}())
    end
    for (observation_index, observation) in enumerate(observations)
        observation_id = FASTX.FASTA.identifier(observation)
        observed_sequence = FASTX.FASTA.sequence(observation)
        if length(observed_sequence) < stranded_kmer_graph.gprops[:k]
            @error "skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        else
            observed_path = sequence_to_stranded_path(stranded_kmer_graph.gprops[:stranded_kmers], observed_sequence)
            i = 1
            ui, ui_orientation = observed_path[i]
            ui_coverage = (observation_index => (i => ui_orientation ))
            push!(stranded_kmer_graph.vprops[ui][:coverage], ui_coverage)
            for i in 2:length(observed_path)
                vi, vi_orientation = observed_path[i]
                vi_coverage = (observation_index => (i => vi_orientation))
                push!(stranded_kmer_graph.vprops[vi][:coverage], vi_coverage)
                edge_coverage = ui_coverage => vi_coverage
                if Graphs.has_edge(stranded_kmer_graph, ui, vi)
                    push!(stranded_kmer_graph.eprops[Graphs.Edge(ui, vi)][:coverage], edge_coverage)
                else
                    Graphs.add_edge!(stranded_kmer_graph, ui, vi, Dict(:coverage => [edge_coverage]))
                end
                # not sure this is necessary
#                 ui′ = stranded_kmer_graph.gprops[:reverse_complement_map][ui]
#                 vi′ = stranded_kmer_graph.gprops[:reverse_complement_map][vi]
#                 if !LightGraphs.has_edge(stranded_kmer_graph, vi′, ui′)
#                     LightGraphs.add_edge!(stranded_kmer_graph, vi′, ui′, Dict(:coverage => Vector{typeof(edge_coverage)}()))
#                 end
                ui, ui_orientation = vi, vi_orientation
                ui_coverage = vi_coverage
            end
            push!(stranded_kmer_graph.gprops[:observed_paths], observed_path)
            push!(stranded_kmer_graph.gprops[:observation_ids], observation_id)
            push!(stranded_kmer_graph.gprops[:observation_color_map], observation_index)
        end
    end
    @info Graphs.nv(stranded_kmer_graph)
    @info Graphs.ne(stranded_kmer_graph)
    return stranded_kmer_graph
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

Compute the SHA-256 hash of the contents of a file.

# Arguments
- `file::AbstractString`: The path to the file for which the SHA-256 hash is to be computed.

# Returns
- `String`: The SHA-256 hash of the file contents, represented as a hexadecimal string.
"""
function sha256_file(file::AbstractString)
    @assert isfile(file)
    return SHA.bytes2hex(SHA.sha256(read(file)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run a function `f` in parallel over a collection of `items` with a progress meter.

# Arguments
- `f::Function`: The function to be applied to each item in the collection.
- `items::AbstractVector`: A collection of items to be processed.

# Description
This function creates a progress meter to track the progress of processing each item in the `items` collection. 
It uses multithreading to run the function `f` on each item in parallel, updating the progress meter after each item is processed.
"""
function run_parallel_progress(f::Function, items::AbstractVector)
    # Create a progress meter
    p = ProgressMeter.Progress(length(items))
    
    # Create a lock and vector to store errors
    lock = ReentrantLock()
    errors = Vector{Union{Nothing, Tuple{Any, Any}}}(nothing, length(items))
    
    Threads.@threads for (i, item) in enumerate(items)
        try
            f(item)
        catch e
            lock() do
                errors[i] = (e, item)
            end
        end
        lock() do
            ProgressMeter.next!(p)
        end
    end
    return errors
end

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module
