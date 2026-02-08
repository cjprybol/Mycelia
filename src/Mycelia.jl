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
import Logging
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

# Rhizomorph graph ecosystem (primary implementation)
include("rhizomorph/rhizomorph.jl")

# Legacy qualmer utilities (still used by iterative-assembly.jl)
include("qualmer-analysis.jl")

# Assembly pipeline (depends on Rhizomorph graph ecosystem)
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
# Using explicit include() statements enables full static analysis by ExplicitImports.jl
include("alignments-and-mapping.jl")
include("amino-acid-analysis.jl")
include("annotation.jl")
include("autocycler.jl")
include("bcalm.jl")
include("binning.jl")
include("bioconda.jl")
include("busco-datasets.jl")
include("checkpointing.jl")
include("classification.jl")
include("clustering.jl")
include("codon-optimization.jl")
include("coverage-clustering.jl")
include("dimensionality-reduction.jl")
include("distance-metrics.jl")
include("foldseek.jl")
include("genome-features.jl")
include("ggcat.jl")
include("graph-cleanup.jl")
include("kmer-analysis.jl")
include("kmer-saturation-analysis.jl")
include("metagraph.jl")
include("metagenomic-classification.jl")
include("ncbi-datasets-cli.jl")
include("pangenome-analysis.jl")
include("pantools.jl")
include("plotting-and-visualization.jl")
include("prokrustean.jl")
include("quality-control-and-benchmarking.jl")
include("rclone.jl")
include("read-quality-control.jl")
include("reference-databases.jl")
include("relational-matrices.jl")
include("sentencepiece.jl")
include("sequence-comparison.jl")
include("simulation.jl")
include("slurm-templates.jl")
include("slurm-sbatch.jl")
include("taxonomy-and-trees.jl")
include("testing-utilities.jl")
include("variant-analysis.jl")
include("viterbi-polishing-and-error-correction.jl")
include("xam.jl")

# PrecompileTools workload for faster startup
# This must be included last, after all other definitions are loaded
include("precompile_workload.jl")

end # module
