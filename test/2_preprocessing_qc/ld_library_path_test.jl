# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/ld_library_path_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/ld_library_path_test.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia

Test.@testset "LD_LIBRARY_PATH sanitization" begin
    ## __init__() should clear LD_LIBRARY_PATH on Linux to prevent system
    ## libraries from conflicting with Julia's bundled JLL artifacts.

    if Sys.islinux()
        Test.@testset "clears non-empty LD_LIBRARY_PATH" begin
            old = get(ENV, "LD_LIBRARY_PATH", nothing)
            try
                ENV["LD_LIBRARY_PATH"] = "/usr/lib/x86_64-linux-gnu"
                Mycelia.__init__()
                Test.@test ENV["LD_LIBRARY_PATH"] == ""
            finally
                ## Restore original state
                if old === nothing
                    delete!(ENV, "LD_LIBRARY_PATH")
                else
                    ENV["LD_LIBRARY_PATH"] = old
                end
            end
        end

        Test.@testset "no-op when LD_LIBRARY_PATH is empty" begin
            old = get(ENV, "LD_LIBRARY_PATH", nothing)
            try
                ENV["LD_LIBRARY_PATH"] = ""
                Mycelia.__init__()
                Test.@test ENV["LD_LIBRARY_PATH"] == ""
            finally
                if old === nothing
                    delete!(ENV, "LD_LIBRARY_PATH")
                else
                    ENV["LD_LIBRARY_PATH"] = old
                end
            end
        end

        Test.@testset "no-op when LD_LIBRARY_PATH is unset" begin
            old = get(ENV, "LD_LIBRARY_PATH", nothing)
            try
                delete!(ENV, "LD_LIBRARY_PATH")
                Mycelia.__init__()
                Test.@test !haskey(ENV, "LD_LIBRARY_PATH")
            finally
                if old === nothing
                    delete!(ENV, "LD_LIBRARY_PATH")
                else
                    ENV["LD_LIBRARY_PATH"] = old
                end
            end
        end
    else
        Test.@testset "no-op on non-Linux platforms" begin
            old = get(ENV, "LD_LIBRARY_PATH", nothing)
            try
                ENV["LD_LIBRARY_PATH"] = "/some/path"
                Mycelia.__init__()
                ## On non-Linux, __init__() should leave LD_LIBRARY_PATH untouched
                Test.@test ENV["LD_LIBRARY_PATH"] == "/some/path"
            finally
                if old === nothing
                    delete!(ENV, "LD_LIBRARY_PATH")
                else
                    ENV["LD_LIBRARY_PATH"] = old
                end
            end
        end
    end
end
