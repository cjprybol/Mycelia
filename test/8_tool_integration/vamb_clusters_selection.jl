# Pure-Julia coverage for selecting VAMB's contig-to-bin clusters table.
#
# VAMB writes `<prefix>_clusters_metadata.tsv` (per-cluster statistics) next to
# `<prefix>_clusters_unsplit.tsv` and, when binsplitting, `<prefix>_clusters_split.tsv`
# (both `clustername\tcontigname`). The metadata table sorts first by name, so a
# name-sorted pattern match returned it instead of an assignments table.
#
# Runs in default CI (no conda, no external tools):
#   julia --project=. -e 'include("test/8_tool_integration/vamb_clusters_selection.jl")'

import Test
import Mycelia

function _touch_all(dir::String, names::Vector{String})
    for name in names
        touch(joinpath(dir, name))
    end
    return dir
end

Test.@testset "VAMB clusters table selection" begin
    Test.@testset "prefers split over unsplit, never metadata" begin
        mktempdir() do dir
            _touch_all(dir,
                ["vae_clusters_metadata.tsv", "vae_clusters_split.tsv",
                    "vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_split.tsv")
        end
    end

    Test.@testset "falls back to unsplit when binsplitting is off" begin
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters_metadata.tsv", "vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
    end

    Test.@testset "metadata alone is not an assignments table" begin
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters_metadata.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") === nothing
        end
    end

    Test.@testset "legacy single-file names" begin
        mktempdir() do dir
            _touch_all(dir, ["vae_clusters.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters.tsv")
        end
        mktempdir() do dir
            _touch_all(dir, ["clusters.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "clusters.tsv")
        end
    end

    Test.@testset "prefix selects the model's own tables" begin
        mktempdir() do dir
            # TaxVAMB names contain "vae_clusters" as a substring; the prefix must
            # match exactly so a `vae` lookup does not pick up `vaevae` tables.
            _touch_all(dir,
                ["vaevae_clusters_metadata.tsv", "vaevae_clusters_split.tsv",
                    "vaevae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") === nothing
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vaevae") ==
                       joinpath(dir, "vaevae_clusters_split.tsv")
        end
        mktempdir() do dir
            # AVAMB writes both VAE and adversarial tables in one directory.
            _touch_all(dir, ["aae_z_clusters_split.tsv", "vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
    end

    Test.@testset "ignores directories and missing outdir" begin
        mktempdir() do dir
            mkdir(joinpath(dir, "vae_clusters_split.tsv"))
            _touch_all(dir, ["vae_clusters_unsplit.tsv"])
            Test.@test Mycelia._vamb_clusters_tsv(dir; prefix = "vae") ==
                       joinpath(dir, "vae_clusters_unsplit.tsv")
        end
        mktempdir() do dir
            Test.@test Mycelia._vamb_clusters_tsv(joinpath(dir, "absent"); prefix = "vae") ===
                       nothing
        end
    end
end
