module Mycelia

__precompile__(false)

# import ArgParse
import BioAlignments
import BioSequences
import BioSymbols
# import CairoMakie
import Clustering
import CodecZlib
import Conda
import DataFrames
import DataStructures
import Dates
import Distances
import Distributions
import DocStringExtensions
import FASTX
import FileIO
import GenomicAnnotations
import GFF3
import GraphRecipes
import Graphs
import GLM
import HTTP
import JSON
import Kmers
import LsqFit
import MetaGraphs
import Mmap
import OrderedCollections
import Plots
import Primes
import ProgressMeter
import Random
import Statistics
import StatsBase
import StatsPlots
import XAM
import uCSV

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
# const AA_ALPHABET = filter(
#     x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x) || BioSymbols.isterm(x)),
#     BioSymbols.alphabet(BioSymbols.AminoAcid))
const AA_ALPHABET = filter(
    x -> !(BioSymbols.isambiguous(x) || BioSymbols.isgap(x)),
    BioSymbols.alphabet(BioSymbols.AminoAcid))

const MAMBA = joinpath(Conda.BINDIR, "mamba")

function add_bioconda_env(pkg; force=false)
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(Conda.conda) env list`))))))
    if !(pkg in current_environments) || force
        @info "installing conda environment $(pkg)"
        run(`$(MAMBA) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
        run(`$(MAMBA) clean --all -y`)
    else
        @info "conda environment $(pkg) already present; set force=true to update/re-install"
    end
end

function add_bioconda_envs(;all=false, force=false)
    if !isfile(MAMBA)
        Conda.add("mamba")
    end
    if !isfile(joinpath(Conda.BINDIR, "pigz"))
        run(`$(MAMBA) install pigz -y`)
    end
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(Conda.conda) env list`))))))
    # https://github.com/JuliaPy/Conda.jl/issues/185#issuecomment-1145149905
    if all
        for pkg in [
            "art",
            # "bioconvert",
            "badread",
            "bcftools",
            "bedtools",
            "blast",
            "clair3-illumina",
            "clair3",    
            # "bwa",
            # "bwa-mem2",
            # "deepvariant",
            "emboss",
            "filtlong",
            # "freebayes",
            "flye",
            "gatk4",
            # "gffread",
            "htslib",
            "megahit",
            "medaka",
            "minimap2",
            "mmseqs2",
            "nanocaller",
            "nanovar",
            # "nanoq",
            # "nanosim",
            # "nanosim-h",
            "ncbi-datasets-cli",
            "pggb",
            "picard",
            # "polypolish",
            "prodigal",
            "raven-assembler",
            "rtg-tools",
            "samtools",
            "sniffles",
            "sourmash",
            "spades",
            "tabix",
            "transtermhp",
            "trim-galore",
            "vcftools",
            "vg"
            ]
            if !(pkg in current_environments) || force
                @info "installing conda environment $(pkg)"
                add_bioconda_env(pkg)
            else
                @info "conda environment $(pkg) already present; set force=true to update/re-install"
            end
        end
    end
    run(`$(MAMBA) clean --all -y`)
end

"""
Submit a command to SLURM using sbatch
"""
function scg_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        partition::String,
        account::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 8,
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
    run(submission)
    sleep(1)
    return true
end

"""
Submit a command to SLURM using sbatch
"""
function nersc_sbatch(;
        job_name::String,
        mail_user::String,
        mail_type::String="ALL",
        logdir::String=pwd(),
        qos::String,
        nodes::Int=1,
        ntasks::Int=1,
        time::String="1-00:00:00",
        cpus_per_task::Int=1,
        mem_gb::Int=cpus_per_task * 2,
        cmd::String,
        constrain::String="cpu"
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
    run(submission)
    sleep(1)
    return true
end

function fasta_genome_size(fasta_file)
    return reduce(sum, map(record -> length(FASTX.sequence(record)), Mycelia.open_fastx(fasta_file)))
end

function gfa_to_fasta(;gfa, fasta=gfa * ".fna")
    open(fasta, "w") do io
        fastx_io = FASTX.FASTA.Writer(io)
        gfa_graph = Mycelia.parse_gfa(gfa)
        for v in Graphs.vertices(gfa_graph)
            record = FASTX.FASTA.Record(gfa_graph.vprops[v][:identifier], gfa_graph.vprops[v][:sequence])
            write(fastx_io, record)
        end
        close(fastx_io)
    end
end

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module