module Mycelia

# Enable precompilation for faster loading
__precompile__(false)  # Disabled for testing framework fix

import AlgebraOfGraphics
import Arrow
import BioAlignments
import BioSequences
import BioSymbols
import CairoMakie
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
import ExpFamilyPCA
import FASTX
import FileIO
import GenomicAnnotations
# import GeoMakie
import GFF3
import GLM
import GraphMakie
import Graphs
import GraphPlot
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
import MetaGraphs
import MetaGraphsNext
import Mmap
import MultivariateStats

# Auto-include all Julia files in src directory
# This ensures our next-generation modules are automatically loaded:
# - sequence-graphs-next.jl (strand-aware k-mer graphs)
# - gfa-io-next.jl (MetaGraphsNext GFA I/O)
# - probabilistic-algorithms-next.jl (probabilistic walks, shortest paths)
# - viterbi-next.jl (enhanced Viterbi with batch processing)
# - graph-algorithms-next.jl (Eulerian paths, bubble detection, repeat resolution)
import OrderedCollections
import Plots
import POMDPs
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
import XMLDict



import Pkg
# TODO remove this import and update function call usage
import Base.Filesystem: stat

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

# Graph type definitions and enums (foundation layer)
include("sequence-graphs-next.jl")  # Defines GraphMode, StrandOrientation enums
include("string-graphs.jl")         # N-gram graphs
include("qualmer-analysis.jl")      # Qualmer types and functions

# Variable-length graph implementations (depend on fixed-length types)
include("fasta-graphs.jl")          # BioSequence graphs (depends on sequence-graphs-next.jl)
include("fastq-graphs.jl")          # Quality-aware BioSequence graphs (depends on qualmer-analysis.jl)

# Assembly pipeline (depends on all graph types and GraphMode enum)
include("assembly.jl")

# Intelligent assembly algorithms (depends on assembly.jl and core graph types)
include("intelligent-assembly.jl")

# Iterative assembly algorithms (depends on intelligent-assembly.jl)
include("iterative-assembly.jl")

# Cross-validation pipeline (depends on both intelligent and iterative assembly)
include("cross-validation.jl")

# Reinforcement learning framework (depends on intelligent-assembly.jl and iterative-assembly.jl)
include("reinforcement-learning.jl")

# Advanced algorithms (depend on core graph types)
include("viterbi-next.jl")

# Load remaining files in alphabetical order (no critical dependencies)
remaining_files = [
    "alignments-and-mapping.jl",
    "annotation.jl",
    "bioconda.jl",
    "classification.jl",
    "clustering.jl",
    "codon-optimization.jl",
    "dimensionality-reduction.jl",
    "distance-metrics.jl",
    "genome-features.jl",
    "kmer-analysis.jl",
    "neo4jl.jl",
    "pangenome-analysis.jl",
    "performance-benchmarks.jl",
    "plotting-and-visualization.jl",
    "quality-control-and-benchmarking.jl",
    "rclone.jl",
    "reference-databases.jl",
    "sequence-comparison.jl",
    "sequence-graphs.jl",
    "simulation.jl",
    "slurm-sbatch.jl",
    "taxonomy-and-trees.jl",
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

end # module
