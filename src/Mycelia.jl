module Mycelia

import BioAlignments
import BioSequences
import BioSymbols
import Clustering
import Dates
import DataFrames
import DataStructures
import Distributions
import Distances
import DocStringExtensions
# import GFF3
import GenomicAnnotations
import Kmers
import GraphRecipes
import Graphs
import LsqFit
import MetaGraphs
import Plots
import PrettyTables
import Primes
import Random
import SparseArrays
import Statistics
import StatsBase
import StatsPlots
import FASTX
import CodecZlib
import HTTP
import ProgressMeter
# import LSHFunctions
import FileIO
import uCSV
import Mmap
import OrderedCollections

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

# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
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