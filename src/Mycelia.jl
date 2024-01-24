module Mycelia

__precompile__(false)

# import ArgParse
import BioAlignments
import BioSequences
import BioSymbols
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

function add_bioconda_env(pkg)
    run(`$(Conda.conda) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
    # Conda.runconda(`create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
end

function add_bioconda_envs(;force=false)
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(Conda.conda) env list`))))))    
    for pkg in [
            "htslib",
            "tabix",
            "bcftools",
            "vcftools",
            "samtools",
            "flye",
            "prodigal",
            "mmseqs2",
            "minimap2"
        ]
        if !(pkg in current_environments) && !force
            @info "installing conda environment $(pkg)"
            add_bioconda_env(pkg)
        else
            @info "conda environment $(pkg) already present; set force=true to update/re-install"
        end
    end
end

function sbatch(;
            job_name::String,
            mail_user::String="ALL",
            mail_type::String,
            logdir::String=pwd(),
            partition::String,
            account::String,
            nodes::Int=1,
            ntasks::Int=1,
            time::String="1-00:00:00",
            cpus_per_task::Int,
            mem_gb::Int,
            cmd::String
        )
        submission = 
        `sbatch
        --job-name=$(job_name)
        --mail-user=$(mail_user)
        --mail-type=$(mail_type)
        --error=$(logdir)/$(job_name).%x-%j.err
        --output=$(logdir)/$(job_name).%x-%j.out
        --partition=$(batch)
        --account=$(account)
        --nodes=$(nodes)
        --ntasks=$(ntasks)
        --time=$(time)   
        --cpus-per-task=$(cpus_per_task)
        --mem=$(mem_gb)G
        --wrap $(cmd)
        `
        run(submission)
end

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't import yourself :)
all_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_julia_files
    include(f)
end

# for f in [
#         "database-interaction.jl",
#         "genome-annotations.jl",
#         "graph-construction.jl",
#         "graph-interaction.jl",
#         "graph-polishing.jl",
#         "graph-simplification.jl",
#         "graph-traversal.jl",
#         "io-transformations.jl",
#         "proteome-analysis.jl",
#         "sequence-interactions.jl",
#         # "viterbi.jl"
#     ]
#     include(f)
# end

end # module