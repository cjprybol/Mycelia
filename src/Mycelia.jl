module Mycelia

__precompile__(false)

import AlgebraOfGraphics
import Arrow
import BioAlignments
import BioSequences
import BioSymbols
import Clustering
# import CodecBase
# import CodecBzip2
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
# import GeoMakie
import GFF3
import GLM
import GraphMakie
import Graphs
# import HDF5
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
import SankeyPlots
import SHA
import SparseArrays
import StableRNGs
import Statistics
import StatsBase
import StatsPlots
import Tar
# import TopoPlots
# import TranscodingStreams
import uCSV
import UUIDs
import XAM
import XMLDict

using CairoMakie

import Pkg

# preserve definitions between code jldoctest code blocks
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Preserving-Definitions-Between-Blocks
# use this to build up a story as we go, where outputs of earlier defined functions feed into
# tests for further downstream functions

# if things fall out of date but look correct, update them automatically
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Fixing-Outdated-Doctests

# const METADATA = joinpath(dirname(dirname(pathof(Mycelia))), "docs", "metadata")
const DNA_ALPHABET = BioSymbols.ACGT
const RNA_ALPHABET = BioSymbols.ACGU
const DEFAULT_BLASTDB_PATH = "$(homedir())/workspace/blastdb"

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
const FASTA_REGEX = r"\.(fa|fasta|fna|fas|fsa|ffn|faa|mpfa|frn)(\.gz)?$"
const FASTQ_REGEX = r"\.(fq|fastq)(\.gz)?$"
const XAM_REGEX = r"\.(sam|bam|cram|sam\.gz)$"
const VCF_REGEX = r"\.vcf(\.gz)?$"

ProgressMeter.ijulia_behavior(:clear)

function check_bioconda_env_is_installed(pkg)
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
    if pkg in current_environments
        return true
    else
        return false
    end
end

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
    channel = nothing
    if occursin("::", pkg)
        println("splitting $(pkg)")
        channel, pkg = split(pkg, "::")
        println("into channel:$(channel) pkg:$(pkg)")
    end
    already_installed = check_bioconda_env_is_installed(pkg)
    if !already_installed || force
        @info "installing conda environment $(pkg)"
        if isnothing(channel)
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
        else
            run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y`)
        end
        run(`$(CONDA_RUNNER) clean --all -y`)
    end
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

Returns the current date as a normalized string with all non-word characters removed.

The output format is based on ISO datetime (YYYYMMDD) but strips any special characters
like hyphens, colons or dots.
"""
function normalized_current_date()
    return replace(Dates.format(Dates.today(), Dates.ISODateFormat), r"[^\w]" => "")
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
    # skip 19 because it is so similar to 17
    results = vcat(
        filter(x -> x != 19, filter(isodd, Primes.primes(0, flip_point))),
        filter(x -> x > flip_point, nearest_prime.(fibonacci_numbers_less_than(max*10)))
    )
    return filter(x -> min <= x <= max, results)
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

Copy files between local and remote storage using rclone with automated retry logic.

# Arguments
- `source::String`: Source path or remote (e.g. "local/path" or "gdrive:folder")
- `dest::String`: Destination path or remote (e.g. "gdrive:folder" or "local/path")

# Keywords
- `config::String=""`: Optional path to rclone config file
- `max_attempts::Int=3`: Maximum number of retry attempts
- `sleep_timer::Int=60`: Initial sleep duration between retries in seconds (doubles after each attempt)
- `includes::Vector{String}=[]`: One or more include patterns (each will be passed using `--include`)
- `excludes::Vector{String}=[]`: One or more exclude patterns (each will be passed using `--exclude`)
- `recursive::Bool=false`: If true, adds the flag for recursive traversal
"""
function rclone_copy2(source, dest;
                     config = "",
                     max_attempts = 3, sleep_timer = 60,
                     includes = String[],
                     excludes = String[],
                     recursive = false)
    done = false
    attempts = 0
    while !done && attempts < max_attempts
        attempts += 1
        try
            # Define base flags optimized for large files
            flags = ["--drive-chunk-size", "2G",
                     "--drive-upload-cutoff", "1T",
                     "--tpslimit", "1",
                     "--verbose"]

            # Append each include pattern with its flag
            for pattern in includes
                push!(flags, "--include")
                push!(flags, pattern)
            end

            # Append each exclude pattern with its flag
            for pattern in excludes
                push!(flags, "--exclude")
                push!(flags, pattern)
            end

            # Optionally add the recursive flag
            if recursive
                push!(flags, "--recursive")
            end

            # Build the full argument list as an array of strings.
            args = String[]
            # Add base command and optional config
            push!(args, "rclone")
            if !isempty(config)
                push!(args, "--config")
                push!(args, config)
            end
            push!(args, "copy")
            # Insert all flags (each flag and its parameter are separate elements)
            append!(args, flags)
            # Add source and destination paths
            push!(args, source)
            push!(args, dest)

            # Convert the argument vector into a Cmd object
            cmd = Cmd(args)

            @info "copying $(source) to $(dest) with command: $(cmd)"
            run(cmd)
            done = true
        catch e
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

# smaller, higher diversity databases do better with >=5 as the denominator - w/ <=4 they run out of memory
# denominator = 5 # produced OOM for NT on NERSC
# denominator = 8 # produced OOM for NT on Lawrencium
# denominator = 10 was only 56% efficient for NT on NERSC
const DEFAULT_MINIMAP_DENOMINATOR=10

function system_mem_to_minimap_index_size(;system_mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85), denominator=DEFAULT_MINIMAP_DENOMINATOR)
    value = Int(floor(system_mem_gb/denominator))
    # 4G is the default
    # this value should be larger for larger memory machines, and smaller for smaller ones
    # it seems related to the total size of the sequences stored in memory, rather than the total size of the in-memory database
    return "$(value)G"
end

function minimap_index(;fasta, mapping_type, mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85), threads=Sys.CPU_THREADS, as_string=false, denominator=DEFAULT_MINIMAP_DENOMINATOR)
    Mycelia.add_bioconda_env("minimap2")
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


function minimap_map_with_index(;
        fasta,
        mapping_type,
        fastq,
        index_file="",
        mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85),
        threads=Sys.CPU_THREADS,
        as_string=false,
        denominator=DEFAULT_MINIMAP_DENOMINATOR
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    
    if !isempty(index_file) && !isfile(index_file)
        error("user-specific index file $index_file does not exist")
    else
        # index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
        index_file_result = Mycelia.minimap_index(fasta=fasta, mapping_type=mapping_type, mem_gb = mem_gb, threads=threads)
        index_file = index_file_result.outfile
    end
    @show index_file
    @assert isfile(index_file)
    outfile = fastq * "." * basename(index_file) * "." * "minimap2.bam"
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(outfile).tmp \\
        | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(outfile).tmp`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -`
        cmd = pipeline(map, compress)
    end
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
        cmd = pipeline(map, compress)
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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Map paired-end reads to a reference sequence using minimap2.

# # Arguments
# - `fasta::String`: Path to reference FASTA file
# - `forward::String`: Path to forward reads FASTQ file
# - `reverse::String`: Path to reverse reads FASTQ file
# - `mem_gb::Integer`: Available system memory in GB
# - `threads::Integer`: Number of threads to use
# - `outdir::String`: Output directory (defaults to forward reads directory)
# - `as_string::Bool=false`: Return command as string instead of Cmd array
# - `mapping_type::String="sr"`: Minimap2 preset ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
# - `denominator::Float64`: Memory scaling factor for index size

# # Returns
# Named tuple containing:
# - `cmd`: Command(s) to execute (String or Array{Cmd})
# - `outfile`: Path to compressed output SAM file (*.sam.gz)

# # Notes
# - Requires minimap2, samtools, and pigz conda environments
# - Automatically compresses output using pigz
# - Index file must exist at `\$(fasta).x\$(mapping_type).I\$(index_size).mmi`
# """
# function minimap_map_paired_end_with_index(;
#         fasta,
#         forward,
#         reverse,
#         mem_gb,
#         threads,
#         outdir = dirname(forward),
#         as_string=false,
#         mapping_type="sr",
#         denominator=DEFAULT_MINIMAP_DENOMINATOR
#     )
#     @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
#     index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
#     index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
#     # @show index_file
#     @assert isfile(index_file) "$(index_file) not found!!"
#     @assert isfile(forward) "$(forward) not found!!"
#     @assert isfile(reverse) "$(reverse) not found!!"
#     fastq_prefix = find_matching_prefix(basename(forward), basename(reverse))
#     temp_sam_outfile = joinpath(outdir, fastq_prefix) * "." * basename(index_file) * "." * "minimap2.sam"
#     # outfile = temp_sam_outfile
#     outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
#     Mycelia.add_bioconda_env("minimap2")
#     Mycelia.add_bioconda_env("samtools")
#     Mycelia.add_bioconda_env("pigz")
#     if as_string
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
#         """
#     else
#         map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#         compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
#         cmd = [map, compress]
#     end
#     return (;cmd, outfile)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Maps paired-end reads to a reference genome using minimap2 and compresses the output.

# # Arguments
# - `fasta::String`: Path to reference genome FASTA file
# - `forward::String`: Path to forward reads FASTQ file
# - `reverse::String`: Path to reverse reads FASTQ file  
# - `mem_gb::Integer`: Available system memory in GB
# - `threads::Integer`: Number of threads to use
# - `outdir::String`: Output directory (defaults to forward reads directory)
# - `as_string::Bool`: Return command as string instead of Cmd array
# - `mapping_type::String`: Mapping preset, e.g. "sr" for short reads (default)
# - `denominator::Float64`: Memory scaling factor for minimap2 index

# # Returns
# Named tuple containing:
# - `cmd`: Command(s) to execute (String or Vector{Cmd})
# - `outfile`: Path to compressed output SAM file (*.sam.gz)

# # Dependencies
# Requires bioconda packages: minimap2, samtools, pigz
# """
# function minimap_map_paired_end(;
#         fasta,
#         forward,
#         reverse,
#         mem_gb,
#         threads,
#         outdir = dirname(forward),
#         as_string=false,
#         mapping_type="sr",
#         denominator=Mycelia.DEFAULT_MINIMAP_DENOMINATOR
#     )
#     index_size = Mycelia.system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
#     @assert isfile(forward) "$(forward) not found!!"
#     @assert isfile(reverse) "$(reverse) not found!!"
#     fastq_prefix = Mycelia.find_matching_prefix(basename(forward), basename(reverse))
#     temp_sam_outfile = joinpath(outdir, fastq_prefix) * "." * "minimap2.sam"
#     # outfile = temp_sam_outfile
#     outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
#     Mycelia.add_bioconda_env("minimap2")
#     Mycelia.add_bioconda_env("samtools")
#     Mycelia.add_bioconda_env("pigz")
#     if as_string
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_size) -ax $(mapping_type) $(fasta) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
#         """
#     else
#         map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_size) -ax $(mapping_type) $(fasta) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#         compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
#         cmd = [map, compress]
#     end
#     return (;cmd, outfile)
# end

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
function fastx_to_contig_lengths(fastx)
    OrderedCollections.OrderedDict(String(FASTX.identifier(record)) => length(FASTX.sequence(record)) for record in Mycelia.open_fastx(fastx))
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

function setup_padloc()
    padloc_is_already_installed = check_bioconda_env_is_installed("padloc")
    if !padloc_is_already_installed
        Mycelia.add_bioconda_env("padlocbio::padloc")
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --db-update`)
    end
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
function run_padloc(;fasta_file, outdir=dirname(abspath(fasta_file)), threads=Sys.CPU_THREADS)
    padloc_outfile = joinpath(outdir, replace(basename(fasta_file), ".fna" => "") * "_padloc.csv")
    if !isfile(padloc_outfile)
        setup_padloc()
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --fna $(fasta_file) --outdir $(outdir) --cpu $(threads)`)
    else
        @info "$(padloc_outfile) already present"
    end
    padloc_faa = replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.faa")
    padloc_gff = replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.gff")
    padloc_domtblout = replace(basename(fasta_file), Mycelia.FASTA_REGEX => ".domtblout")
    if !isfile(padloc_outfile)
        padloc_outfile = missing
    end
    return (csv = padloc_outfile, faa = padloc_faa, gff = padloc_gff, domtblout = padloc_domtblout)
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

Identifies the optimal number of clusters using hierarchical clustering
and maximizing the average silhouette score, displaying progress.

Uses Clustering.clustering_quality for score calculation.

Args:
    distance_matrix: A square matrix of pairwise distances between items.

Returns:
    A tuple containing:
    - hcl: The hierarchical clustering result object.
    - optimal_number_of_clusters: The inferred optimal number of clusters (k).
"""
function identify_optimal_number_of_clusters(distance_matrix, min_k = Int(ceil(log2(size(distance_matrix, 2)))) * 2, max_k = Int(ceil(sqrt(size(distance_matrix, 2)))))
    # Ensure the input is a square matrix
    if size(distance_matrix, 1) != size(distance_matrix, 2)
         error("Input distance_matrix must be square.")
    end
    n_items = size(distance_matrix, 1)
    if n_items < 2
        error("Need at least 2 items to perform clustering.")
    end

    # Compute hierarchical clustering using your custom function.
    println("Performing hierarchical clustering...")
    hcl = heirarchically_cluster_distance_matrix(distance_matrix)
    println("Hierarchical clustering complete.")

    # # Determine N based on hcl object or matrix size
    # local N::Int
    # if hasfield(typeof(hcl), :labels) && !isempty(hcl.labels)
    N = length(hcl.labels)
    if N != n_items
         @warn "Number of labels in Hclust object ($(N)) does not match distance matrix size ($(n_items)). Using Hclust labels count."
    end
    # else
    #      @warn "Hierarchical clustering result type $(typeof(hcl)) might not have a `.labels` field or it's empty. Using distance matrix size for N."
    #      N = n_items
    # end
    
    # Define the range of k values to test: 2 to min(N-1, ceil(sqrt(N)))
    # Cannot have more clusters than N items, and silhouette requires at least 1 item per cluster,
    # implicitly k <= N. Silhouettes are ill-defined for k=1. Need k < N for non-trivial clustering.
    if max_k < 2
        println("Not enough items (N=$N) to test multiple cluster numbers (k >= 2). Cannot determine optimal k.")
        # Optionally, plot a message or return a specific value indicating failure/trivial case
        p_trivial = StatsPlots.plot(title="Cannot optimize k (N=$N, max tested k=$(max_k))", xlabel="Number of Clusters", ylabel="Avg Silhouette Score")
        display(p_trivial)
        # Returning N clusters (each item its own) might be misleading as optimal
        # Consider returning nothing or throwing error, or returning k=N with a warning
        @warn "Returning trivial clustering with k=$N as optimal k could not be determined."
        return (hcl, N)
    end
    ks = min_k:max_k
    n_ks = length(ks)
    println("Evaluating average silhouette scores for k from $(min_k) to $(max_k)...")

    # Preallocate an array for average silhouette scores for each value of k.
    avg_silhouette_scores = zeros(Float64, n_ks)

    # --- ProgressMeter Setup ---
    pm = ProgressMeter.Progress(n_ks, 1, "Calculating Avg Silhouette Scores: ", 50)
    # --------------------------

    # Parallel loop: each iteration computes the average silhouette score for a given k.
    Threads.@threads for i in 1:n_ks
        k = ks[i]
        avg_silhouette_scores[i] = Statistics.mean(Clustering.silhouettes(Clustering.cutree(hcl, k=k), distance_matrix))
        ProgressMeter.next!(pm)
    end

    # Check if all scores resulted in errors
    if all(s -> s == Float64(-Inf) || isnan(s), avg_silhouette_scores)
        error("Failed to calculate valid average silhouette scores for all values of k.")
    end

    # Determine the optimal number of clusters by finding the MAXIMUM average score
    # Filter out -Inf/NaN before finding the maximum.
    valid_indices = findall(s -> !isnan(s) && s != Float64(-Inf), avg_silhouette_scores)
    if isempty(valid_indices)
        error("No valid average silhouette scores were computed.")
    end

    # Find maximum among valid scores
    max_score_valid, local_max_idx_valid = findmax(avg_silhouette_scores[valid_indices])
    # Map back to the original index in avg_silhouette_scores and ks
    max_idx_original = valid_indices[local_max_idx_valid]
    optimal_number_of_clusters = ks[max_idx_original]

    println("\nOptimal number of clusters inferred (max avg silhouette score): ", optimal_number_of_clusters)
    println("Maximum average silhouette score: ", max_score_valid)


    # Plotting the average silhouette scores.
    p = StatsPlots.scatter(
        ks,
        avg_silhouette_scores,
        title = "Clustering Performance vs. Number of Clusters\n(higher is better)",
        xlabel = "Number of Clusters (k)",
        ylabel = "Average Silhouette Score", # Updated label
        label = "Avg Score",
        markersize = 5,
        markerstrokewidth = 0.5
    )

    # Adjust y-limits based on valid scores, keeping [-1, 1] range in mind
    valid_scores = avg_silhouette_scores[valid_indices]
    min_val = minimum(valid_scores)
    max_val = maximum(valid_scores) # This is max_score_valid
    range = max_val - min_val
    # Set sensible limits slightly beyond observed range, but clamp to [-1.05, 1.05]
    ylims_lower = max(-1.05, min_val - range * 0.1)
    ylims_upper = min(1.05, max_val + range * 0.1)
    # Ensure lower < upper, handle case where range is 0
    if isapprox(ylims_lower, ylims_upper)
         ylims_lower -= 0.1
         ylims_upper += 0.1
         # Clamp again if needed
         ylims_lower = max(-1.05, ylims_lower)
         ylims_upper = min(1.05, ylims_upper)
    end
    StatsPlots.ylims!(p, (ylims_lower, ylims_upper))

    # Add reference line at y=0 (separation between reasonable/poor clusters)
    StatsPlots.hline!(p, [0], color=:gray, linestyle=:dot, label="Score = 0")

    # Add vertical line for the optimum
    StatsPlots.vline!(p, [optimal_number_of_clusters], color = :red, linestyle = :dash, label = "Inferred optimum = $(optimal_number_of_clusters)")

    # Display the plot
    display(p)

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
- `outdir::String`: The output directory where the BLASTN results will be saved.
- `fasta::String`: The path to the input FASTA file containing the query sequences.
- `blastdb::String`: The path to the BLAST database to search against.
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
function run_blastn(;outdir=pwd(), fasta, blastdb, threads=min(Sys.CPU_THREADS, 8), task="megablast", force=false, remote=false, wait=true)
    Mycelia.add_bioconda_env("blast")
    outdir = mkpath(outdir)
    outfile = "$(outdir)/$(basename(fasta)).blastn.$(basename(blastdb)).$(task).txt"
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    
    if force || need_to_run
        cmd = 
        `
        $(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastn
        -num_threads $(threads)
        -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        -query $(fasta)
        -db $(blastdb)
        -out $(outfile)
        -max_target_seqs 10
        -subject_besthit
        -task $(task)
        -evalue 0.001
        `
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

# function remove_non_ascii(input::String)
#     return String(filter(x -> x <= '\x7f', input))
# end

# function sanitize_string(input::String)
#     return String(filter(x -> isvalid(Char, x), input))
# end

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
function fasta_and_gff_to_genbank(;fasta, gff, genbank=gff * ".genbank")
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
    return genbank
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

function normalize_kmer_counts(kmer_counts)
    total_kmer_counts = sum(values(kmer_counts))
    normalized_kmer_frequencies = DataStructures.OrderedDict(k => v/total_kmer_counts for (k,v) in kmer_counts)
    return normalized_kmer_frequencies
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
    # codons = last.(collect(Kmers.SpacedKmers{Kmers.DNACodon}(BioSequences.LongDNA{4}(best_sequence), 3)))
    codons = first.(collect(Kmers.UnambiguousDNAMers{3}(BioSequences.LongDNA{4}(best_sequence))))[1:3:end]
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
    # @show (best_likelihood)^-10 / (initial_log_likelihood)^-10
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
        kmers = [Kmers.AAKmer{k}(BioSequences.LongAA(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.DNA
        kmers = [Kmers.DNAKmer{k}(BioSequences.LongDNA{2}(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.RNA
        kmers = [Kmers.RNAKmer{k}(BioSequences.LongRNA{2}(kv)) for kv in kmer_vectors]
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
"""
function fasta_list_to_dense_counts_table(; fasta_list, k, alphabet)
    k >= 11 && error("use fasta_list_to_sparse_counts_table")
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

    progress = ProgressMeter.Progress(length(fasta_list))
    # Prepare thread-safe containers for successful results and error reporting.
    successful_results = Vector{Vector{Int}}()  # each will be a vector of counts for one file
    successful_indices = Vector{Int}()
    errors = Vector{Tuple{Int, String}}()  # (file index, error message)
    
    # Use a single lock for both progress updating and pushing into shared arrays.
    reentrant_lock = ReentrantLock()
    
    Threads.@threads for (entity_index, fasta_file) in collect(enumerate(fasta_list))
        # Update the progress meter within the lock so output remains synchronized.
        lock(reentrant_lock) do
            ProgressMeter.next!(progress)
        end
        try
            # Process the FASTA (or FASTQ) file to count kmers.
            entity_mer_counts = COUNT(Kmers.Kmer{KMER_TYPE, k}, fasta_file)
            # Build a vector of counts for each sorted kmer.
            local_counts = [ get(entity_mer_counts, kmer, 0) for kmer in sorted_kmers ]
            lock(reentrant_lock) do
                push!(successful_results, local_counts)
                push!(successful_indices, entity_index)
            end
        catch e
            # Report the issue and save the error details.
            lock(reentrant_lock) do
                @warn "Error processing file: $fasta_file. Error: $e"
                push!(errors, (entity_index, string(e)))
            end
        end
    end
    
    if isempty(successful_results)
        error("All files failed to process.")
    end
    # Assemble the final counts matrix from the successful results.
    kmer_counts_matrix = hcat(successful_results...)  # Each column is the count vector for one file
    successful_fasta_list = fasta_list[successful_indices]
    
    return (; sorted_kmers, kmer_counts_matrix, successful_fasta_list)
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
"""
function biosequences_to_dense_counts_table(;biosequences, k)
    k >= 11 && error("use sparse counts to count k >= 11")    
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the current time as a Unix timestamp (seconds since epoch).

# Returns
- `Int`: Current time as an integer Unix timestamp (seconds since January 1, 1970 UTC)

# Examples
```julia
unix_time = current_unix_datetime()
# => 1709071368 (example value, will differ based on current time)
```
"""
function current_unix_datetime()
    return Int(floor(Dates.datetime2unix(Dates.now())))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the optimal subsequence length based on error rate distribution.

# Arguments
- `error_rate`: Single error rate or array of error rates (between 0 and 1)
- `threshold`: Desired probability that a subsequence is error-free (default: 0.95)
- `sequence_length`: Maximum sequence length to consider for plotting
- `plot_result`: If true, returns a plot of probability vs. length

# Returns
- If `plot_result=false`: Integer representing optimal subsequence length
- If `plot_result=true`: Tuple of (optimal_length, plot)

# Examples
```julia
# Single error rate
optimal_subsequence_length(error_rate=0.01)

# Array of error rates
optimal_subsequence_length(error_rate=[0.01, 0.02, 0.01])

# With more stringent threshold
optimal_subsequence_length(error_rate=0.01, threshold=0.99)

# Generate plot
length, p = optimal_subsequence_length(error_rate=0.01, plot_result=true)
Plots.display(p)
```
"""
function optimal_subsequence_length(;error_rate::Union{Real, AbstractArray{<:Real}},
                                   threshold::Float64=0.95,
                                   sequence_length::Union{Nothing, Int}=nothing,
                                   plot_result::Bool=false)
    # Handle array input by calculating mean error rate
    avg_error_rate = isa(error_rate, AbstractArray) ? Statistics.mean(error_rate) : error_rate
    
    # Validate inputs
    if avg_error_rate <= 0
        optimal_length = typemax(Int)
    elseif avg_error_rate >= 1
        optimal_length = 1
    else
        # Calculate optimal length where P(error-free) >= threshold
        optimal_length = floor(Int, log(threshold) / log(1 - avg_error_rate))
        optimal_length = max(1, optimal_length)  # Ensure at least length 1
    end

    # Return early if no plot requested
    if !plot_result
        return optimal_length
    end
    
    # For plotting, determine sequence length to display
    max_length = isnothing(sequence_length) ? 2 * optimal_length : sequence_length
    
    # Calculate probabilities for different lengths
    lengths = 1:max_length
    probabilities = [(1 - avg_error_rate)^len for len in lengths]
    
    # Create DataFrame for plotting
    df = DataFrames.DataFrame(
        Length = collect(lengths),
        Probability = probabilities,
        Optimal = lengths .== optimal_length
    )

    quality_score = Mycelia.error_rate_to_q_value(error_rate)
    rounded_quality_score = Int(floor(quality_score))
    rounded_error_rate = round(avg_error_rate, digits=4)
    rounded_threshold_rate = round(threshold, digits=2)
    plot_title = 
    """
    Optimal Kmer Length Inference
    Error Rate: $(rounded_error_rate * 100)%≈Q$(rounded_quality_score)
    Threshold: $(rounded_threshold_rate * 100)% of kmers expected to be correct
    """
    
    # Create plot
    p = Plots.plot(
        df.Length, df.Probability,
        linewidth=2, 
        label="P(error-free)",
        xlabel="Subsequence Length",
        ylabel="Probability of Error-free Match",
        title=plot_title,
        grid=true,
        ylims=(0,1),
        alpha=0.8
    )
    
    # Add horizontal line for threshold
    Plots.hline!([threshold], linestyle=:dash, color=:red, label="Threshold")
    
    # Add vertical line for optimal length
    Plots.vline!([optimal_length], linestyle=:dash, color=:green, label="Optimal length: $optimal_length")
    
    # Highlight optimal point
    Plots.scatter!([optimal_length], [threshold], color=:orange, markersize=8, label="")
    
    return optimal_length, p
end

function get_biosequence_alphabet(s::T) where T<:BioSequences.BioSequence
    return first(T.parameters)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Write a DataFrame to a JLD2 file using a standardized internal name.
"""
function JLD2_write_table(;df::DataFrames.DataFrame, filename::String)
    JLD2.jldopen(filename, "w") do file
        file["dataframe"] = df  # Always use the same internal name
    end
    return filename
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Read a DataFrame from a JLD2 file without needing to know the internal name.
If the file contains multiple DataFrames, returns the first one found.
"""
function JLD2_read_table(filename::String)
    df = JLD2.jldopen(filename, "r") do file
        # Try standard name first
        if haskey(file, "dataframe")
            return file["dataframe"]
        end
        
        # Otherwise search for any DataFrame
        for key in keys(file)
            if typeof(file[key]) <: DataFrames.DataFrame
                return file[key]
            end
        end
        
        # No DataFrame found
        error("No DataFrame found in file: $filename")
    end
    
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert all InlineString columns in a DataFrame to standard Strings.
Modifies the dataframe in-place and returns it.
"""
function sanitize_inline_strings!(df::DataFrames.DataFrame)
    for col in names(df)
        if eltype(df[!, col]) <: InlineStrings.InlineString
            df[!, col] = String.(df[!, col])
        end
    end
    return df
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a column to standard Strings if it contains InlineStrings,
otherwise return the original column unchanged.
"""
function sanitize_inline_strings(v::AbstractVector)
    if eltype(v) <: InlineStrings.InlineString
        return String.(v)
    else
        return v
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function dataframe_convert_dicts_to_json(df)
    df_copy = DataFrames.copy(df)
    for col in DataFrames.names(df_copy)
        if eltype(df_copy[!, col]) <: AbstractDict || any(x -> isa(x, AbstractDict), df_copy[!, col])
            df_copy[!, col] = [isa(cell, AbstractDict) ? JSON.json(cell) : cell for cell in df_copy[!, col]]
        end
    end
    return df_copy
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return a string representation of the vector `v` with each element on a new line,
mimicking valid Julia syntax. The output encloses the elements in square brackets
and separates them with a comma followed by a newline.
"""
function repr_long(v)
    buf = IOBuffer()
    println(buf, "[")
    for x in v
        println(buf, "    \"$x\",")
    end
    println(buf, "]")
    return String(take!(buf))
end


function rclone_copy_list(;source::String, destination::String, relative_paths::Vector{String})
    # Create a temporary file for storing file paths
    temp_file = joinpath(tempdir(), "rclone_sources_$(Random.randstring(8)).txt")
    
    try
        # Write paths to temp file
        open(temp_file, "w") do file
            for path in relative_paths
                println(file, path)
            end
        end
        
        println("Starting transfer of $(length(relative_paths)) files from $source to $destination...")
        
        # Make sure the destination directory exists if it's local
        if !occursin(":", destination) && !isdir(destination)
            mkdir(destination)
        end
        
        # Download files using rclone with progress reporting
        # --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 
        run(`rclone copy $source $destination --verbose --files-from $temp_file`)
        
        println("Download completed successfully")
        return true
    catch e
        println("Error: $e")
        return false
    finally
        # Clean up temp file
        if isfile(temp_file)
            rm(temp_file)
            println("Temporary file removed")
        end
    end
end
        
"""
    dataframe_to_ndjson(df::DataFrame; outfile::Union{String,Nothing}=nothing)

Converts a DataFrame `df` into a newline-delimited JSON (NDJSON) string.
Each line in the returned string represents one DataFrame row in JSON format,
suitable for upload to Google BigQuery.

# Keyword Arguments
- `outfile::Union{String,Nothing}`: If provided, writes the resulting NDJSON to the file path given.

# Examples
```julia
using DataFrames, Dates

# Sample DataFrame
df = DataFrame(
    id = [1, 2, 3],
    name = ["Alice", "Bob", "Carol"],
    created = [DateTime(2025, 4, 8, 14, 30), DateTime(2025, 4, 8, 15, 0), missing]
)

ndjson_str = dataframe_to_ndjson(df)
println(ndjson_str)

# Optionally, write to a file
dataframe_to_ndjson(df; outfile="output.ndjson")
"""
function dataframe_to_ndjson(df::DataFrames.DataFrame; outfile::Union{String, Nothing}=nothing) ndjson_lines = String[]
    # Iterate over each row in the DataFrame
    for row in eachrow(df)
        row_dict = Dict{String, Any}()
    
        # Build a dictionary for the current row
        for (col, value) in pairs(row)
            # Convert column names to strings
            col_name = string(col)
            if value === missing
                # Convert missing values to `nothing` which JSON prints as null
                row_dict[col_name] = nothing
            elseif isa(value, Dates.DateTime)
                # Format DateTime value in ISO 8601 format with millisecond precision
                # Modify the format string if a different format is needed for BigQuery.
                formatted = Dates.format(value, Dates.dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
                row_dict[col_name] = formatted
            else
                row_dict[col_name] = value
            end
        end
    
        # Convert the dictionary to a JSON string and push into the list
        push!(ndjson_lines, JSON.json(row_dict))
    end
    
    # Join all JSON strings with newline delimiters (one JSON object per line)
    ndjson_str = join(ndjson_lines, "\n")
    
    # Optionally write the NDJSON content to a file if an outfile path is provided.
    if outfile !== nothing
        open(outfile, "w") do io
            write(io, ndjson_str)
        end
        return outfile
    else
        return ndjson_str
    end
end

function install_cloud_cli()
    # Get the user's home directory
    home_dir = ENV["HOME"]
    
    # Define the target installation directory ($HOME/google-cloud-sdk)
    sdk_install_dir = joinpath(home_dir, "google-cloud-sdk")
    sdk_bin_dir = joinpath(sdk_install_dir, "bin")
    
    # Function to update the PATH in current session and persist to .bashrc if needed.
    function update_path()
        # Update current session PATH
        if !occursin(sdk_bin_dir, ENV["PATH"])
            ENV["PATH"] = "$sdk_bin_dir:" * ENV["PATH"]
        end
        
        bashrc_file = joinpath(home_dir, ".bashrc")
        path_is_set = false
        if isfile(bashrc_file)
            for line in eachline(bashrc_file)
                if occursin("google-cloud-sdk/bin", line)
                    path_is_set = true
                    break
                end
            end
        end
        
        if !path_is_set
            println("Appending SDK bin path to $bashrc_file ...")
            open(bashrc_file, "a") do io
                println(io, "\n# Added by install_cloud_cli() on $(Dates.now())")
                println(io, "export PATH=\"$sdk_bin_dir:\$PATH\"")
            end
            println("Please reload your shell (e.g., run: source ~/.bashrc) to update your PATH.")
        else
            println("PATH already includes google-cloud-sdk/bin")
        end
    end
    
    # Check if the SDK is already installed.
    if isdir(sdk_install_dir)
        println("Google Cloud CLI is already installed.")
        update_path()
        println("No further installation required.")
        return
    end

    # Not installed: proceed with downloading and installing.
    tarball = "google-cloud-cli-linux-x86_64.tar.gz"
    sdk_url = "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/$tarball"
    
    println("Downloading Google Cloud CLI from: $sdk_url ...")
    run(`curl -O $sdk_url`)
    
    println("Extracting $tarball to $home_dir ...")
    run(`tar -xf $tarball -C $home_dir`)
    
    # Run the install script from within the extracted folder.
    install_script = joinpath(sdk_install_dir, "install.sh")
    println("Running the installation script ...")
    run(`bash $install_script --quiet --usage-reporting false --command-completion false --path-update false`)
    
    # Cleanup the downloaded tarball on successful install.
    println("Cleaning up the downloaded tarball...")
    rm(tarball, force=true)
    
    update_path()
    
    println("Installation complete. You now have gsutil and bq installed as part of the Google Cloud CLI.")
    println("Please run `gcloud init` or `gcloud auth login` and follow the interactive prompts to log in.")
end


function upload_dataframe_to_bigquery(;
    ndjson_file::String,
    project_id::String,
    dataset_id::String,
    table_id::String,
    gcs_bucket_name::String,
    bq_location::String, # e.g., "US"
    verbose=true
)
    full_table_id = "$(project_id):$(dataset_id).$(table_id)"
    if verbose
        println("Target BigQuery Table: ", full_table_id)
        println("Target GCS Bucket: ", gcs_bucket_name)
        println("BigQuery Location: ", bq_location)
    end
    
    try
        @assert isfile(ndjson_file) && filesize(ndjson_file) > 0
        verbose && println("NDJSON detected locally at: ", ndjson_file)

        # --- 3. Upload NDJSON to GCS ---
        verbose && println("Uploading to GCS bucket: ", gcs_bucket_name)
        gcs_uri = "gs://$(gcs_bucket_name)/$(basename(ndjson_file))"

        # Using gsutil command (simpler integration)
        upload_cmd = `gsutil cp $(ndjson_file) $(gcs_uri)`
        verbose && println("Running command: ", upload_cmd)
        run(upload_cmd)
        verbose && println("Upload to GCS complete: ", gcs_uri)

        # --- 4. Trigger BigQuery Load Job ---
        verbose && println("Triggering BigQuery load job for: ", full_table_id)

        # Using bq command (simpler integration)
        # Explicitly setting schema source_format and location
        load_cmd = `bq load --source_format=NEWLINE_DELIMITED_JSON --location=$(bq_location) $(full_table_id) $(gcs_uri)`
        # If you have a local schema file (e.g., schema.json):
        # load_cmd = `bq load --source_format=NEWLINE_DELIMITED_JSON --location=$(bq_location) $(full_table_id) $(gcs_uri) ./schema.json`
        # Or rely on auto-detect (use with caution):
        # load_cmd = `bq load --source_format=NEWLINE_DELIMITED_JSON --autodetect --location=$(bq_location) $(full_table_id) $(gcs_uri)`

        verbose && println("Running command: ", load_cmd)
        run(load_cmd)
        verbose && println("BigQuery load job initiated.")

        # --- 5. Cleanup GCS File (Optional) ---
        # Rely on bucket lifecycle rules or explicitly remove the temp file
        cleanup_cmd = `gsutil rm $(gcs_uri)`
        verbose && println("Running command: ", cleanup_cmd)
        run(cleanup_cmd)
        verbose && println("Cleaned up GCS file: ", gcs_uri)

        verbose && println("Process complete for table: ", full_table_id)

    catch e
        println("Error during process: ", e)
        # Add more robust error handling
    end
end

function fasta_to_kmer_and_codon_frequencies(fasta)
    kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{5}, fasta)
    pyrodigal_results = Mycelia.run_pyrodigal(fasta_file = fasta)
    if occursin(r"\.gz$", fasta)
        unzipped_fasta = replace(fasta, ".gz" => "")
        run(pipeline(`gzip -dc $(fasta)`, unzipped_fasta))
        @assert isfile(unzipped_fasta)
    end
    genbank = Mycelia.fasta_and_gff_to_genbank(fasta = unzipped_fasta, gff = pyrodigal_results.gff)
    codon_frequencies = Mycelia.genbank_to_codon_frequencies(genbank)
    return (;kmer_counts, codon_frequencies)
end

function calculate_sequence_likelihood_from_kmer_profile(sequence, normalized_kmer_counts)
    initial_likelihood = 1.0
    for (k, i) in Kmers.UnambiguousDNAMers{5}(sequence)
        canonical_k = BioSequences.canonical(k)
        # # display(k)
        # try
        #     @assert BioSequences.iscanonical(k)
        # catch
        #     display(k)
        # end
        initial_likelihood *= normalized_kmer_counts[canonical_k]
        # end
    end
    initial_likelihood
end

function load_bvbrc_genome_metadata(; 
    summary_url = "ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary",
    metadata_url = "ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata")
    
    # Create a unique temporary directory
    temp_dir = joinpath(tempdir(), "bvbrc_temp_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))")
    mkpath(temp_dir)
    
    try
        # Define temporary file paths
        summary_file = joinpath(temp_dir, "genome_summary.tsv")
        metadata_file = joinpath(temp_dir, "genome_metadata.tsv")
        
        # Download files to temporary location
        @info "Downloading genome summary from $(summary_url)"
        Downloads.download(summary_url, summary_file)
        
        @info "Downloading genome metadata from $(metadata_url)"
        Downloads.download(metadata_url, metadata_file)
        
        # Read files into DataFrames
        @info "Reading genome summary file"
        genome_summary = CSV.read(summary_file, DataFrames.DataFrame, delim='\t', header=1, 
                                 types=Dict("genome_id" => String))
        
        @info "Reading genome metadata file"
        genome_metadata = CSV.read(metadata_file, DataFrames.DataFrame, delim='\t', header=1, 
                                  types=Dict("genome_id" => String))
        
        # Join the DataFrames
        @info "Joining genome summary and metadata"
        bvbrc_genome_summary = DataFrames.innerjoin(genome_summary, genome_metadata, 
                                                  on="genome_id", makeunique=true)
        
        return bvbrc_genome_summary
    finally
        # Clean up temporary files regardless of success or failure
        @info "Cleaning up temporary files"
        rm(temp_dir, recursive=true, force=true)
    end
end

"""
    parallel_pyrodigal(normalized_fastas::Vector{String})

Runs Mycelia.run_pyrodigal on a list of FASTA files in parallel using Threads.

Args:
    normalized_fastas: A vector of strings, where each string is a path to a FASTA file.

Returns:
    A tuple containing two elements:
    1. successes (Vector{Tuple{String, Any}}): A vector of tuples, where each tuple contains the
       filename and the result returned by a successful Mycelia.run_pyrodigal call.
    2. failures (Vector{Tuple{String, String}}): A vector of tuples, where each tuple contains the
       filename and the error message string for a failed Mycelia.run_pyrodigal call.
"""
function parallel_pyrodigal(normalized_fastas::Vector{String})
    num_files = Base.length(normalized_fastas)
    Base.println("Processing $(num_files) FASTA files using $(Threads.nthreads()) threads...")

    # Create a Progress object for manual updates
    p = ProgressMeter.Progress(num_files, 1, "Running Pyrodigal: ", 50)

    # Use Channels to collect results and failures thread-safely
    # Channel{Tuple{Filename, ResultType}} - adjust ResultType if known
    successes = Base.Channel{Tuple{String, Any}}(num_files)
    failures = Base.Channel{Tuple{String, String}}(num_files)

    # Use Threads.@threads for parallel execution
    Threads.@threads for fasta_file in normalized_fastas
        result = nothing # Initialize result variable in the loop's scope
        try
            # --- Execute the function ---
            # Base.println("Thread $(Threads.threadid()) processing: $(fasta_file)") # Optional: for debugging
            result = Mycelia.run_pyrodigal(fasta_file = fasta_file) # Capture the result

            # --- Store success ---
            Base.put!(successes, (fasta_file, result))

        catch e
            # --- Store failure ---
            err_msg = Base.sprint(Base.showerror, e) # Get the error message as a string
            Base.println(Base.stderr, "ERROR processing $(fasta_file) on thread $(Threads.threadid()): $(err_msg)")
            Base.put!(failures, (fasta_file, err_msg))
        finally
            # --- Always update progress ---
            ProgressMeter.next!(p)
        end
    end

    # Close channels now that all threads are done writing
    Base.close(successes)
    Base.close(failures)

    # Collect results and failures from the channels
    successful_results = Base.collect(successes)
    failed_files = Base.collect(failures)

    # --- Report Summary ---
    Base.println("\n--- Pyrodigal Processing Summary ---")
    num_success = Base.length(successful_results)
    num_failed = Base.length(failed_files)
    Base.println("Successfully processed: $(num_success)")
    Base.println("Failed: $(num_failed)")

    if !Base.isempty(failed_files)
        Base.println("\nFailures:")
        for (file, err) in failed_files
            Base.println("- File: $(file)\n  Error: $(err)")
        end
    end
    Base.println("------------------------------------")

    return successful_results, failed_files # Return both successes and failures
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a sparse kmer counts table (SparseMatrixCSC) from a list of FASTA files using a 3-pass approach.
Pass 1 (Parallel): Counts kmers per file and writes to temporary JLD2 files.
Pass 2 (Serial): Aggregates unique kmers, max count, nnz per file, and rarefaction data from temp files.
                 Generates and saves a k-mer rarefaction plot.
Pass 3 (Parallel): Reads temporary counts again to construct the final sparse matrix.

# Arguments
- `fasta_list::AbstractVector{<:AbstractString}`: A list of paths to FASTA files.
- `k::Integer`: The length of the kmer.
- `alphabet::Symbol`: The alphabet type (:AA, :DNA, :RNA).
- `temp_dir_parent::AbstractString`: Parent directory for creating the temporary working directory. Defaults to `Base.tempdir()`.
- `count_element_type::Union{Type{<:Unsigned}, Nothing}`: Optional. Specifies the unsigned integer type for the counts. If `nothing` (default), the smallest `UInt` type capable of holding the maximum observed count is used.
- `rarefaction_data_filename::AbstractString`: Filename for the TSV output of rarefaction data. Defaults to "kmer_rarefaction_data_3pass.tsv".
- `rarefaction_plot_basename::AbstractString`: Basename for the output rarefaction plots. Defaults to "kmer_rarefaction_curve_3pass".
- `show_rarefaction_plot::Bool`: Whether to display the rarefaction plot after generation. Defaults to `true`.
- `rarefaction_plot_kwargs...`: Keyword arguments to pass to `plot_kmer_rarefaction` for plot customization.

# Returns
- `NamedTuple{(:kmers, :counts, :rarefaction_data_path)}`:
    - `kmers`: A sorted `Vector` of unique kmer objects.
    - `counts`: A `SparseArrays.SparseMatrixCSC{V, Int}` storing kmer counts.
    - `rarefaction_data_path`: Path to the saved TSV file with rarefaction data.

# Raises
- `ErrorException`: If input `fasta_list` is empty, alphabet is invalid, or required Kmer/counting functions are not found.
"""
function fasta_list_to_sparse_counts_table(;
    fasta_list::AbstractVector{<:AbstractString},
    k::Integer,
    alphabet::Symbol,
    temp_dir_parent::AbstractString = Base.tempdir(),
    count_element_type::Union{Type{<:Unsigned}, Nothing} = nothing,
    rarefaction_data_filename::AbstractString = "$(normalized_current_datetime()).$(lowercase(string(alphabet)))$(k)mer_rarefaction.tsv",
    rarefaction_plot_basename::AbstractString = replace(rarefaction_data_filename, ".tsv" => ""),
    show_rarefaction_plot::Bool = true,
    rarefaction_plot_kwargs... 
)
    # --- 0. Input Validation and Setup ---
    num_files = length(fasta_list)
    if num_files == 0
        error("Input fasta_list is empty.")
    end

    KMER_TYPE = if alphabet == :AA
        isdefined(Kmers, :AAKmer) ? Kmers.AAKmer{k} : error("Kmers.AAKmer not found or Kmers.jl not loaded correctly.")
    elseif alphabet == :DNA
        isdefined(Kmers, :DNAKmer) ? Kmers.DNAKmer{k} : error("Kmers.DNAKmer not found or Kmers.jl not loaded correctly.")
    elseif alphabet == :RNA
        isdefined(Kmers, :RNAKmer) ? Kmers.RNAKmer{k} : error("Kmers.RNAKmer not found or Kmers.jl not loaded correctly.")
    else
        error("Invalid alphabet: $alphabet. Choose from :AA, :DNA, :RNA")
    end

    COUNT_FUNCTION = if alphabet == :DNA
        func_name = :count_canonical_kmers
        isdefined(Main, func_name) ? getfield(Main, func_name) : (isdefined(@__MODULE__, func_name) ? getfield(@__MODULE__, func_name) : error("$func_name not found"))
    else
        func_name = :count_kmers
        isdefined(Main, func_name) ? getfield(Main, func_name) : (isdefined(@__MODULE__, func_name) ? getfield(@__MODULE__, func_name) : error("$func_name not found"))
    end

    temp_dir = Base.Filesystem.mktempdir(temp_dir_parent; prefix="kmer_counts_3pass_")
    Base.@info "Using temporary directory for intermediate counts: $temp_dir"
    temp_file_paths = [Base.Filesystem.joinpath(temp_dir, "counts_$(i).jld2") for i in 1:num_files]
    
    # --- Pass 1: Count Kmers per file and Write to Temp Files (Parallel) ---
    Base.@info "Pass 1: Counting $(KMER_TYPE) for $num_files files and writing temps using $(Base.Threads.nthreads()) threads..."
    progress_pass1 = ProgressMeter.Progress(num_files; desc="Pass 1 (Counting & Writing Temps): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:cyan)
    progress_pass1_lock = Base.ReentrantLock()

    Base.Threads.@threads for original_file_idx in 1:num_files
        fasta_file = fasta_list[original_file_idx]
        temp_filename = temp_file_paths[original_file_idx]
        counts_dict = Dict{KMER_TYPE, Int}()
        try
            counts_dict = COUNT_FUNCTION(KMER_TYPE, fasta_file)
            JLD2.save_object(temp_filename, counts_dict)
        catch e
            Base.println("Error processing file $fasta_file (idx $original_file_idx) during Pass 1: $e")
            # Save an empty dict if error occurred to ensure file exists, simplifying later checks
            try
                JLD2.save_object(temp_filename, Dict{KMER_TYPE, Int}())
            catch save_err
                 Base.@error "Failed to save empty placeholder for errored file $fasta_file to $temp_filename: $save_err"
            end
        end
        Base.lock(progress_pass1_lock) do
            ProgressMeter.next!(progress_pass1)
        end
    end
    ProgressMeter.finish!(progress_pass1)
    Base.@info "Pass 1 finished."

    # --- Pass 2: Aggregate Unique Kmers, Max Count, NNZ per file, Rarefaction Data (Serial) ---
    Base.@info "Pass 2: Aggregating stats and rarefaction data from temporary files (serially)..."
    all_kmers_set = Set{KMER_TYPE}()
    max_observed_count = 0
    total_non_zero_entries = 0
    nnz_per_file = Vector{Int}(undef, num_files)
    rarefaction_points = Vector{Tuple{Int, Int}}() # Grow as we process

    progress_pass2 = ProgressMeter.Progress(num_files; desc="Pass 2 (Aggregating Serially): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:yellow)

    for i in 1:num_files
        temp_filename = temp_file_paths[i]
        current_file_nnz = 0
        try
            if Base.Filesystem.isfile(temp_filename)
                counts_dict = JLD2.load_object(temp_filename)
                current_file_nnz = length(counts_dict)
                if current_file_nnz > 0
                    for k in keys(counts_dict)
                        push!(all_kmers_set, k)
                    end
                    current_file_max_val = maximum(values(counts_dict))
                    if current_file_max_val > max_observed_count
                        max_observed_count = current_file_max_val
                    end
                end
            else
                Base.@warn "Temporary file not found during Pass 2: $temp_filename. Assuming 0 counts for this file."
            end
        catch e
            Base.@error "Error reading or processing temporary file $temp_filename (for original file $i) during Pass 2: $e. Assuming 0 counts for this file."
        end
        nnz_per_file[i] = current_file_nnz
        total_non_zero_entries += current_file_nnz
        push!(rarefaction_points, (i, length(all_kmers_set))) # Record after processing file i
        ProgressMeter.next!(progress_pass2; showvalues = [(:unique_kmers, length(all_kmers_set))])
    end
    ProgressMeter.finish!(progress_pass2)
    Base.@info "Pass 2 aggregation finished."

    # --- Process and Save/Plot Rarefaction Data (after Pass 2) ---
    default_output_dir = pwd() 
    rarefaction_data_path = Base.Filesystem.joinpath(default_output_dir, rarefaction_data_filename)
    Base.@info "Saving rarefaction data to $rarefaction_data_path..."
    try
        data_to_write = [ [pt[1], pt[2]] for pt in rarefaction_points ]
        if !isempty(data_to_write)
             DelimitedFiles.writedlm(rarefaction_data_path, data_to_write, '\t')
        else
             Base.@warn "No rarefaction points recorded. TSV file will be empty or not created."
             DelimitedFiles.writedlm(rarefaction_data_path, Array{Int}(undef,0,2), '\t')
        end
    catch e
        Base.@error "Failed to write rarefaction data to $rarefaction_data_path: $e"
    end
    
    if isfile(rarefaction_data_path) && !isempty(rarefaction_points)
        Base.@info "Generating k-mer rarefaction plot..."
        try
            # Ensure plot_kmer_rarefaction is accessible.
            # This function must be defined in the current scope or imported.
            plot_kmer_rarefaction( 
                rarefaction_data_path;
                output_dir = default_output_dir,
                output_basename = rarefaction_plot_basename,
                display_plot = show_rarefaction_plot,
                rarefaction_plot_kwargs...
            )
        catch e
            Base.@error "Failed to generate rarefaction plot: $e. Ensure Makie and a backend are correctly set up. Also ensure 'plot_kmer_rarefaction' function is loaded."
        end
    else
        Base.@warn "Skipping rarefaction plot generation as data file is missing or empty."
    end

    # --- Determine Value Type, Prepare for Pass 3 ---
    ValType = if isnothing(count_element_type)
        if max_observed_count <= typemax(UInt8)
            UInt8
        elseif max_observed_count <= typemax(UInt16)
            UInt16
        elseif max_observed_count <= typemax(UInt32)
            UInt32
        else
            UInt64
        end
    else
        count_element_type
    end
    if !isnothing(count_element_type) && max_observed_count > typemax(count_element_type)
        Base.@warn "User-specified count_element_type ($count_element_type) may be too small for the maximum observed count ($max_observed_count)."
    end
    Base.@info "Using element type $ValType for kmer counts. Max observed count: $max_observed_count."

    if isempty(all_kmers_set) && total_non_zero_entries == 0
        Base.@warn "No kmers found across any files, or errors prevented aggregation."
        return (; kmers=Vector{KMER_TYPE}(), counts=SparseArrays.spzeros(ValType, Int, 0, num_files), rarefaction_data_path=rarefaction_data_path)
    end

    sorted_kmers = sort(collect(all_kmers_set))
    num_kmers = length(sorted_kmers)
    empty!(all_kmers_set); all_kmers_set = nothing; GC.gc() # Release memory

    kmer_to_row_map = Dict{KMER_TYPE, Int}(kmer => i for (i, kmer) in enumerate(sorted_kmers))
    Base.@info "Found $num_kmers unique kmers. Total non-zero entries: $total_non_zero_entries."

    # --- Pass 3: Prepare Sparse Matrix Data (In Parallel from Temp Files) ---
    Base.@info "Pass 3: Preparing data for sparse matrix construction using $(Base.Threads.nthreads()) threads..."
    
    # Ensure total_non_zero_entries is accurate based on nnz_per_file sum
    # This is important if any files failed in Pass 2 and nnz_per_file[i] was set to 0 for them.
    actual_total_nnz_from_files = sum(nnz_per_file)
    if total_non_zero_entries != actual_total_nnz_from_files
        Base.@warn "Sum of nnz_per_file ($actual_total_nnz_from_files) differs from serially accumulated total_non_zero_entries ($total_non_zero_entries). Using sum of nnz_per_file."
        total_non_zero_entries = actual_total_nnz_from_files # Use the more robust sum for allocation
    end
    
    if total_non_zero_entries == 0 && num_kmers > 0
         Base.@warn "Found $num_kmers unique kmers, but $total_non_zero_entries non-zero entries. Matrix will be empty of values."
    elseif total_non_zero_entries == 0 && num_kmers == 0
         Base.@warn "No k-mers and no non-zero entries. Resulting matrix will be empty."
    end

    row_indices = Vector{Int}(undef, total_non_zero_entries)
    col_indices = Vector{Int}(undef, total_non_zero_entries)
    values_vec = Vector{ValType}(undef, total_non_zero_entries)

    write_offsets = Vector{Int}(undef, num_files + 1)
    write_offsets[1] = 0
    for i in 1:num_files
        write_offsets[i+1] = write_offsets[i] + nnz_per_file[i] 
    end
    
    progress_pass3 = ProgressMeter.Progress(num_files; desc="Pass 3 (Filling Sparse Data): ", barglyphs=ProgressMeter.BarGlyphs("[=> ]"), color=:blue)
    progress_pass3_lock = Base.ReentrantLock()

    Base.Threads.@threads for original_file_idx in 1:num_files
        if nnz_per_file[original_file_idx] > 0 
            temp_filename = temp_file_paths[original_file_idx]
            file_start_offset = write_offsets[original_file_idx] 
            current_entry_in_file = 0 
            try
                if Base.Filesystem.isfile(temp_filename)
                    counts_dict_pass3 = JLD2.load_object(temp_filename) 
                    for (kmer, count_val) in counts_dict_pass3
                        if count_val > 0
                            row_idx_val = Base.get(kmer_to_row_map, kmer, 0)
                            if row_idx_val > 0 
                                global_idx = file_start_offset + current_entry_in_file + 1 
                                if global_idx <= total_non_zero_entries # total_non_zero_entries is now sum(nnz_per_file)
                                    row_indices[global_idx] = row_idx_val
                                    col_indices[global_idx] = original_file_idx
                                    values_vec[global_idx] = ValType(count_val)
                                    current_entry_in_file += 1
                                else
                                    Base.@error "Internal error: Exceeded sparse matrix capacity (total_non_zero_entries = $total_non_zero_entries) for file $original_file_idx. Global Idx: $global_idx. Kmer: $kmer."
                                    break 
                                end
                            end
                        end
                    end
                    if current_entry_in_file != nnz_per_file[original_file_idx]
                        Base.@warn "Mismatch in NNZ for file $original_file_idx ($(fasta_list[original_file_idx])) during Pass 3. Expected $(nnz_per_file[original_file_idx]), wrote $current_entry_in_file."
                    end
                else
                     Base.@warn "Temp file $temp_filename (for original file $original_file_idx) not found in Pass 3, though nnz_per_file was $(nnz_per_file[original_file_idx])."
                end
            catch e
                Base.@error "Error reading or processing temporary file $temp_filename (for original file $original_file_idx) in Pass 3: $e."
            end
        end
        Base.lock(progress_pass3_lock) do
            ProgressMeter.next!(progress_pass3)
        end
    end
    ProgressMeter.finish!(progress_pass3)
    Base.@info "Pass 3 finished."
    
    Base.@info "Constructing sparse matrix ($num_kmers rows, $num_files columns, $total_non_zero_entries non-zero entries)..."
    
    kmer_counts_sparse_matrix = if total_non_zero_entries > 0 && num_kmers > 0
        SparseArrays.sparse(
            row_indices[1:total_non_zero_entries], 
            col_indices[1:total_non_zero_entries], 
            values_vec[1:total_non_zero_entries], 
            num_kmers, 
            num_files
        )
    else
        SparseArrays.spzeros(ValType, Int, num_kmers, num_files)
    end

    Base.@info "Done. Returning sorted kmer list, sparse counts matrix, and rarefaction data path."
    final_result = (; kmers=sorted_kmers, counts=kmer_counts_sparse_matrix, rarefaction_data_path=rarefaction_data_path)
    
    # Consider adding cleanup for temp_dir here or instructing the user.
    # Base.@info "Temporary directory $temp_dir can be manually removed."
    # try
    #     Base.Filesystem.rm(temp_dir; recursive=true, force=true)
    #     Base.@info "Successfully removed temporary directory: $temp_dir"
    # catch e
    #     Base.@warn "Could not remove temporary directory $temp_dir: $e"
    # end

    return final_result
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Plots a k-mer rarefaction curve from data stored in a TSV file.
The TSV file should contain two columns:
1. Number of FASTA files processed.
2. Cumulative unique k-mers observed at that point.

The plot is displayed and saved in PNG, PDF, and SVG formats.

# Arguments
- `rarefaction_data_path::AbstractString`: Path to the TSV file containing rarefaction data.
- `output_dir::AbstractString`: Directory where the output plots will be saved. Defaults to the directory of `rarefaction_data_path`.
- `output_basename::AbstractString`: Basename for the output plot files (without extension). Defaults to the basename of `rarefaction_data_path` without its original extension.
- `display_plot::Bool`: Whether to display the plot interactively. Defaults to `true`.

# Keyword Arguments
- `fig_size::Tuple{Int, Int}`: Size of the output figure, e.g., `(1000, 750)`.
- `title::AbstractString`: Title of the plot.
- `xlabel::AbstractString`: Label for the x-axis.
- `ylabel::AbstractString`: Label for the y-axis.
- `line_color`: Color of the plotted line.
- `line_style`: Style of the plotted line (e.g. `:dash`, `:dot`).
- `marker`: Marker style for points (e.g. `:circle`, `:xcross`).
- `markersize::Number`: Size of the markers.
- Any other keyword arguments will be passed to `Makie.Axis`.
"""
function plot_kmer_rarefaction(
    rarefaction_data_path::AbstractString;
    output_dir::AbstractString = dirname(rarefaction_data_path),
    output_basename::AbstractString = first(Base.Filesystem.splitext(basename(rarefaction_data_path))),
    display_plot::Bool = true,
    # Makie specific customizations
    fig_size::Tuple{Int, Int} = (1000, 750),
    title::AbstractString = "K-mer Rarefaction Curve",
    xlabel::AbstractString = "Number of FASTA Files Processed",
    ylabel::AbstractString = "Cumulative Unique K-mers",
    line_color = :blue,
    line_style = nothing,
    marker = nothing,
    markersize::Number = 10,
    axis_kwargs... # Capture other axis properties
)
    if !isfile(rarefaction_data_path)
        Base.@error "Rarefaction data file not found: $rarefaction_data_path"
        return nothing
    end

    data = try
        DelimitedFiles.readdlm(rarefaction_data_path, '\t', Int, header=false)
    catch e
        Base.@error "Failed to read rarefaction data from $rarefaction_data_path: $e"
        return nothing
    end

    if size(data, 2) != 2
        Base.@error "Rarefaction data file $rarefaction_data_path must have exactly two columns."
        return nothing
    end

    files_processed = data[:, 1]
    unique_kmers = data[:, 2]

    # Sort data by files_processed for a proper line plot,
    # as the input might be from unordered parallel processing.
    sort_indices = sortperm(files_processed)
    files_processed = files_processed[sort_indices]
    unique_kmers = unique_kmers[sort_indices]
    
    Base.mkpath(output_dir) # Ensure output directory exists

    fig = Makie.Figure(size = fig_size)
    ax = Makie.Axis(
        fig[1, 1],
        title = title,
        xlabel = xlabel,
        ylabel = ylabel;
        axis_kwargs... # Pass through other axis settings
    )

    Makie.lines!(ax, files_processed, unique_kmers, color = line_color, linestyle = line_style)
    if !isnothing(marker)
        Makie.scatter!(ax, files_processed, unique_kmers, color = line_color, marker = marker, markersize = markersize)
    end


    if display_plot
        Base.@info "Displaying rarefaction plot..."
        Makie.display(fig)
    end

    output_path_png = Base.Filesystem.joinpath(output_dir, output_basename * ".png")
    output_path_pdf = Base.Filesystem.joinpath(output_dir, output_basename * ".pdf")
    output_path_svg = Base.Filesystem.joinpath(output_dir, output_basename * ".svg")

    try
        Base.@info "Saving plot to $output_path_png"
        Makie.save(output_path_png, fig)
        Base.@info "Saving plot to $output_path_pdf"
        Makie.save(output_path_pdf, fig)
        Base.@info "Saving plot to $output_path_svg"
        Makie.save(output_path_svg, fig)
    catch e
        Base.@error "Failed to save plot: $e. Make sure a Makie backend (e.g., CairoMakie for saving, GLMakie for display) is active and correctly configured in your environment."
    end
    
    return fig # Return the figure object
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
        outfile = joinpath(outdir, basename(fastx) * ".tsv.gz")
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

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module
