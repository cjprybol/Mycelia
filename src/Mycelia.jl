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
import MultivariateStats
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

import JSON
import DataFrames
import ProgressMeter
import CodecZlib
import Base.Filesystem: stat

# preserve definitions between code jldoctest code blocks
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Preserving-Definitions-Between-Blocks
# use this to build up a story as we go, where outputs of earlier defined functions feed into
# tests for further downstream functions

# if things fall out of date but look correct, update them automatically
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Fixing-Outdated-Doctests



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




# dynamic import of files??
all_julia_files = filter(x -> occursin(r"\.jl$", x), readdir(dirname(pathof(Mycelia))))
# don't recusively import this file
all_other_julia_files = filter(x -> x != "Mycelia.jl", all_julia_files)
for f in all_other_julia_files
    include(f)
end

end # module
