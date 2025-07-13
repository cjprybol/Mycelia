module Mycelia

__precompile__(false)

# ============================================================================
# DEPENDENCY IMPORTS
# ============================================================================

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
import LsqFit
import Luxor
import Makie
import MetaGraphs
import MetaGraphsNext
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

# ============================================================================
# DEPRECATION SYSTEM
# ============================================================================

# Include deprecation macros
include("deprecation.jl")

# ============================================================================
# SUBMODULE SYSTEM
# ============================================================================

# Include submodules
include("submodules/IO.jl")
include("submodules/Sequences.jl")
include("submodules/Assembly.jl")
include("submodules/Alignment.jl")
include("submodules/Annotation.jl")
include("submodules/Taxonomy.jl")
include("submodules/SequenceClustering.jl")
include("submodules/Variants.jl")
include("submodules/Visualization.jl")
include("submodules/QualityControl.jl")
include("submodules/Utils.jl")

# Expose submodules (but do NOT export their functions)
# Note: We don't use 'using' to avoid namespace conflicts with imported packages
# Instead, the submodules are available as Mycelia.IO, Mycelia.Sequences, etc.

# ============================================================================
# LEGACY SYSTEM (TEMPORARY - for backward compatibility)
# ============================================================================

# Include legacy files temporarily to avoid breaking existing code
# These will be removed once functions are fully migrated to submodules
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# Don't recursively import these files
excluded_files = ["Mycelia.jl", "deprecation.jl", "api_migration_plan.jl"]
legacy_files = filter(x -> !(x in excluded_files) && !startswith(x, "submodules/"), all_julia_files)

for f in legacy_files
    include(f)
end

# ============================================================================
# NO TOP-LEVEL EXPORTS
# ============================================================================
# All functions must be called via Mycelia.Submodule.function_name
# This forces users to be explicit about which functionality they're using

# ============================================================================
# DEPRECATION SHIMS - Backward compatibility warnings
# ============================================================================

# Create deprecation warnings for functions currently called at top level
# These direct users to the new submodule structure

# Example: when someone calls Mycelia.open_fastx(), they get a warning
# directing them to use Mycelia.IO.open_fastx() instead
@deprecated_toplevel open_fastx "Mycelia.IO.open_fastx"
@deprecated_toplevel write_fasta "Mycelia.IO.write_fasta"
@deprecated_toplevel write_fastq "Mycelia.IO.write_fastq"
@deprecated_toplevel count_kmers "Mycelia.Sequences.count_kmers"
@deprecated_toplevel count_canonical_kmers "Mycelia.Sequences.count_canonical_kmers"
@deprecated_toplevel run_megahit "Mycelia.Assembly.run_megahit"
@deprecated_toplevel run_flye "Mycelia.Assembly.run_flye"

# More deprecation shims will be added as we identify all currently used functions

end # module
