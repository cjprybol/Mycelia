module Mycelia

__precompile__(false)

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
# import HDF5
import Hungarian
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
import MultivariateStats
import OrderedCollections
import Plots
import Primes
import Printf
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

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module
