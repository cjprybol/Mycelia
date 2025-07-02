import Pkg
if isinteractive()
    Pkg.activate("..")
else
    # this should be set via the `--project=` flag
end

using Revise
using Test
import Mycelia
import FASTX
import DocStringExtensions
import BioSequences
import Dates
import Random
import SHA
import CSV
import uCSV
import DataFrames
import Arrow

# Dynamically include all test files in test/sets/
for testfile in sort(readdir(joinpath(@__DIR__, "sets")))
    if endswith(testfile, ".jl")
        include(joinpath("sets", testfile))
    end
end
