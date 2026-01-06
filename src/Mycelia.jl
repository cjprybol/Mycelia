module Mycelia

# Enable precompilation for faster loading; set JULIA_PKG_PRECOMPILE=0 to skip
__precompile__()

import AlgebraOfGraphics
import Arrow
import Base58
import Base64
import BioAlignments
import BioSequences
import BioSymbols
import Blake3Hash
import CairoMakie
import Clustering
# import CodecBase
# import CodecBzip2
import CodecBGZF
import CodecZlib
import Colors
import ColorSchemes
import Conda
import CRC32c
import CSV
import DataFrames
import DataStructures
import Dates
import DelimitedFiles
import Distances
import Distributions
import DocStringExtensions
import Downloads
# import ExpFamilyPCA
import EzXML
import FASTX
import FileIO
import GenomicAnnotations
import GFF3
# import GeoMakie
import GLM
import GraphMakie
import GraphPlot
import Graphs
# import HDF5
import Hungarian
import HTTP
import JLD2
import JSON
import Karnak
import Kmers
import LinearAlgebra
import LsqFit
import Luxor
import Makie
import MD5
import MetaGraphsNext
import Mmap
import MultivariateStats
import OrderedCollections
import Pkg
import Plots
import POMDPs
import PooledArrays
import PrecompileTools
import Primes
import Printf
import ProgressMeter
import Random
import ReinforcementLearning
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
import UMAP
import UUIDs
import XAM
import XLSX
import XMLDict

# preserve definitions between code jldoctest code blocks
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Preserving-Definitions-Between-Blocks
# use this to build up a story as we go, where outputs of earlier defined functions feed into
# tests for further downstream functions

# if things fall out of date but look correct, update them automatically
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Fixing-Outdated-Doctests

# Dependency-ordered loading to ensure type-stable definitions
# Core types and utilities must be loaded first
include("utility-functions.jl")
include("alphabets.jl")
include("constants.jl")
include("fastx.jl")

# Topological data analysis utilities (depends on Graphs)
include("tda.jl")

# Graph type definitions and enums (foundation layer)
include("graph-core.jl")            # Shared graph enums (OLD - to be deprecated)

# Rhizomorph graph ecosystem (NEW - correct implementation)
include("rhizomorph/rhizomorph.jl")
include("kmer-graphs.jl")           # K-mer graph construction utilities (OLD - to be replaced)
include("sequence-graphs-next.jl")  # Higher-level graph algorithms (depends on enums & k-mer types)
include("string-graphs.jl")         # N-gram graphs
include("qualmer-analysis.jl")      # Qualmer types and functions
include("qualmer-graphs.jl")        # Qualmer graph construction utilities

# Variable-length graph implementations (depend on fixed-length types)
include("fasta-graphs.jl")          # BioSequence graphs (depends on sequence-graphs-next.jl)
include("fastq-graphs.jl")          # Quality-aware BioSequence graphs (depends on qualmer-analysis.jl)

# Assembly pipeline (depends on all graph types and GraphMode enum)
include("assembly.jl")

# Intelligent assembly algorithms (depends on assembly.jl and core graph types)
# include("intelligent-assembly.jl")

# Iterative assembly algorithms (depends on intelligent-assembly.jl)
include("iterative-assembly.jl")

# Advanced algorithms (depend on core graph types)
include("viterbi-next.jl")

# # Cross-validation pipeline (depends on both intelligent and iterative assembly)
# include("cross-validation.jl")

# to be added back as we finish their implementations
# Reinforcement learning framework (depends on intelligent-assembly.jl and iterative-assembly.jl)
# include("reinforcement-learning.jl")
# include("reinforcement-learning-rl-jl.jl")
# include("reinforcement-learning-pomdp.jl")
# include("reinforcement-learning-mcts.jl")
# include("reinforcement-learning-comparison.jl")

# Load remaining files in alphabetical order (no critical dependencies)
remaining_files = [
    "alignments-and-mapping.jl",
    "amino-acid-analysis.jl",
    "annotation.jl",
    "autocycler.jl",
    "bcalm.jl",
    "bioconda.jl",
    "binning.jl",
    "classification.jl",
    "metagraph.jl",
    "metagenomic-classification.jl",
    "clustering.jl",
    "codon-optimization.jl",
    "dimensionality-reduction.jl",
    "distance-metrics.jl",
    "foldseek.jl",
    "genome-features.jl",
    "ggcat.jl",
    "kmer-analysis.jl",
    "neo4jl.jl",
    "pangenome-analysis.jl",
    "pantools.jl",
    "performance-benchmarks.jl",
    "plotting-and-visualization.jl",
    "prokrustean.jl",
    "quality-control-and-benchmarking.jl",
    "read-quality-control.jl",
    "busco-datasets.jl",
    "ncbi-datasets-cli.jl",
    "rclone.jl",
    "reference-databases.jl",
    "sentencepiece.jl",
    "sequence-comparison.jl",
    "sequence-graphs.jl",
    "simulation.jl",
    "slurm-sbatch.jl",
    "taxonomy-and-trees.jl",
    "testing-utilities.jl",
    "variant-analysis.jl",
    "viterbi-polishing-and-error-correction.jl",
    "xam.jl"
]

for file in remaining_files
    file_path = joinpath(dirname(pathof(Mycelia)), file)
    if isfile(file_path)
        include(file)
    end
end

# PrecompileTools workload for faster startup
# This must be included last, after all other definitions are loaded
include("precompile_workload.jl")

end # module
