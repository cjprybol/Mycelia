# run me from the root of the Mycelia package directory:
# julia --project=. -e "import Pkg; Pkg.test()"
# julia --project=. --color=yes -e "import Pkg; Pkg.test()"
# julia --project=. --code-coverage=user -e "import Pkg; Pkg.test()"
# julia --project=. -e "import Pkg; Pkg.instantiate(); Pkg.test()"
# julia --project=test -e 'using Pkg; Pkg.test()'

import Pkg
if isinteractive()
    Pkg.activate(".")
else
    # this should be set via the `--project=` flag
end

# using Revise
# using Test
# import Mycelia
# import FASTX
# import DocStringExtensions
# import BioSequences
# import Dates
# import Random
# import SHA
# import CSV
# import uCSV
# import DataFrames
# import Arrow
# import CodecZlib

# Recursively include all test files in the new nested structure, in logical order
function include_all_tests(dir)
    for (root, dirs, files) in walkdir(dir)
        for file in sort(files)
            if endswith(file, ".jl")
                include(joinpath(root, file))
            end
        end
    end
end

# Include in logical order
include_all_tests(joinpath(@__DIR__, "1_data_acquisition"))
include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"))
include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"))
include_all_tests(joinpath(@__DIR__, "4_assembly"))
include_all_tests(joinpath(@__DIR__, "5_validation"))
include_all_tests(joinpath(@__DIR__, "6_annotation"))
include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"))
include_all_tests(joinpath(@__DIR__, "8_tool_integration"))
