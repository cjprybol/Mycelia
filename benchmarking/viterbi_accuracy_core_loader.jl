# Lightweight benchmark-only Mycelia entry point for the Viterbi accuracy harness.
#
# This file intentionally lives under benchmarking/ rather than src/ so the
# precompiled Mycelia package has a stable module layout. HPC smoke runs can
# include this local module to avoid unrelated plotting/ML imports while the
# package itself remains unaffected by MYCELIA_CORE_BENCHMARK.

module Mycelia

const _MYCELIA_SRC = normpath(joinpath(@__DIR__, "..", "src"))

import BioSequences
import BioSymbols
import Conda
import DataStructures
import DocStringExtensions
import FASTX
import Graphs
import Kmers
import LinearAlgebra
import MetaGraphs
import MetaGraphsNext
import ProgressMeter
import Random
import Statistics

include(joinpath(_MYCELIA_SRC, "alphabets.jl"))
include(joinpath(_MYCELIA_SRC, "constants.jl"))

module Rhizomorph

import ..Mycelia
import BioSequences
import DataStructures
import DocStringExtensions
import FASTX
import Graphs
import Kmers
import LinearAlgebra
import MetaGraphsNext
import Random
import Statistics

const _RHIZOMORPH_SRC = joinpath(Mycelia._MYCELIA_SRC, "rhizomorph")

include(joinpath(_RHIZOMORPH_SRC, "core", "enums.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "evidence-structures.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "vertex-data.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "edge-data.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "evidence-functions.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "quality-functions.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "graph-query.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "graph-construction.jl"))
include(joinpath(_RHIZOMORPH_SRC, "core", "graph-type-conversions.jl"))

include(joinpath(_RHIZOMORPH_SRC, "algorithms", "path-finding.jl"))

include(joinpath(_RHIZOMORPH_SRC, "fixed-length", "kmer-graphs.jl"))
include(joinpath(_RHIZOMORPH_SRC, "fixed-length", "qualmer-graphs.jl"))
include(joinpath(_RHIZOMORPH_SRC, "fixed-length", "ngram-graphs.jl"))

end

include(joinpath(_MYCELIA_SRC, "viterbi-next.jl"))

end
