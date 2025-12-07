import Test
import Mycelia

Test.@testset "Bandage download and CLI (opt-in)" begin
    if !Sys.islinux() || get(ENV, "MYCELIA_RUN_BANDAGE_DOWNLOAD", "false") != "true"
        @info "Skipping Bandage download/invoke test (requires Linux and MYCELIA_RUN_BANDAGE_DOWNLOAD=true)"
        return
    end

    mktempdir() do dir
        bandage_bin = Mycelia.download_bandage(joinpath(dir, "bin"))
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
                gfa=gfa,
                img=joinpath(dir, "graph.png"),
                extra_args=["--dpi", "1"],
                force=true
            )
            Test.@test isfile(img_path)

            reduced = Mycelia.bandage_reduce(
                fastg=fastg,
                gfa=joinpath(dir, "reduced.gfa")
            )
            Test.@test isfile(reduced)
        end
    end
end
