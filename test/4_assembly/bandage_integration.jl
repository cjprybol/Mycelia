# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/4_assembly/bandage_integration.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/bandage_integration.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia

Test.@testset "Bandage download and CLI (opt-in)" begin
    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    if !Sys.islinux() || !run_external
        @info "Skipping Bandage download/invoke test (requires Linux and MYCELIA_RUN_EXTERNAL=true)"
        return
    end

    mktempdir() do dir
        bandage_bin = try
            Mycelia.download_bandage(joinpath(dir, "bin"))
        catch e
            if e isa Mycelia.BandageCompatibilityError
                @info "Skipping Bandage download/invoke test (incompatible binary)" error=sprint(showerror, e)
                return
            end
            rethrow()
        end
        Test.@test isfile(bandage_bin)

        gfa = joinpath(dir, "graph.gfa")
        open(gfa, "w") do io
            write(io, "H\tVN:Z:1.0\nS\t1\tACGT\n")
        end

        fastg = joinpath(dir, "graph.fastg")
        open(fastg, "w") do io
            write(io, ">NODE_1_length_4\nACGT\n")
        end

        Base.withenv(
            "MYCELIA_BANDAGE_CMD" => bandage_bin,
            "QT_QPA_PLATFORM" => "offscreen",
            "XDG_RUNTIME_DIR" => dir
        ) do
            img_path = Mycelia.bandage_visualize(
                gfa = gfa,
                img = joinpath(dir, "graph.png"),
                extra_args = ["--dpi", "1"],
                force = true
            )
            Test.@test isfile(img_path)

            reduced = Mycelia.bandage_reduce(
                fastg = fastg,
                gfa = joinpath(dir, "reduced.gfa")
            )
            Test.@test isfile(reduced)
        end
    end
end
