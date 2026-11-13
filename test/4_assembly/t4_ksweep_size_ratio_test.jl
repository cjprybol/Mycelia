# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/t4_ksweep_size_ratio_test.jl")'
# ```
#
# td-ouix: t4_ksweep.jl should name size_ratio, pick the most contiguous k by n50/n_contigs,
# fix the stale real_genome_benchmark.jl header, and caveat the committed
# CSVs. This test reads the script (it does not run the benchmark).

import Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SCRIPT = joinpath(REPO_ROOT, "benchmarking", "t4_ksweep.jl")
const CAVEAT = joinpath(REPO_ROOT, "benchmarking", "results",
    "real_genome_benchmark_SIZE_RATIO_CAVEAT.md")

Test.@testset "t4_ksweep size_ratio and best-k (td-ouix)" begin
    src = read(SCRIPT, String)

    Test.@testset "file header names t4_ksweep.jl" begin
        Test.@test occursin("benchmarking/t4_ksweep.jl", src)
        Test.@test !occursin("benchmarking/real_genome_benchmark.jl", src)
    end

    Test.@testset "column and selection use size_ratio and most contiguous k" begin
        Test.@test occursin(":size_ratio", src) || occursin("size_ratio", src)
        Test.@test occursin("most contiguous", src)
        Test.@test occursin("most_contiguous_score", src)
        Test.@test occursin("most_contiguous_k_row", src)
        Test.@test occursin(r"n50\s*/\s*n_contigs", src)
        Test.@test !occursin("assembly_quality_score", src)
        Test.@test !occursin("function best_k_row", src)
        # The live column must not still be named genome_fraction.
        Test.@test !occursin(r"assembly_results\[!, :genome_fraction\]", src)
        Test.@test !occursin(r"argmax\(subset\.genome_fraction\)", src)
        Test.@test !occursin("by genome fraction", src)
    end

    Test.@testset "touched helpers are type-annotated" begin
        Test.@test occursin(r"function download_reference\([^)]*::", src)
        Test.@test occursin(r"function run_rhizomorph_assembly\([^)]*::", src)
    end

    Test.@testset "committed evidence CSVs carry a caveat" begin
        Test.@test isfile(CAVEAT)
        text = isfile(CAVEAT) ? read(CAVEAT, String) : ""
        Test.@test occursin("genome_fraction", text)
        Test.@test occursin("size_ratio", text)
        Test.@test occursin(r"(?i)total_length/ref_size|total_length / ref_size", text)
        Test.@test occursin(r"100%", text) || occursin(r"exceed", text)
        Test.@test occursin("most contiguous", text)
    end
end
