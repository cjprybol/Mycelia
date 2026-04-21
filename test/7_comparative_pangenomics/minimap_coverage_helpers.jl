# Pure-Julia unit tests for the internal helpers used by
# pairwise_minimap_fasta_comparison. These run in default CI because they do
# not invoke minimap2 or any bioconda tool.
#
# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/minimap_coverage_helpers.jl")'
# ```

import Test
import Mycelia
import FASTX
import DataFrames

Test.@testset "Minimap coverage helpers" begin
    Test.@testset "_merge_intervals — half-open [start, end) semantics" begin
        Test.@testset "empty input" begin
            Test.@test Mycelia._merge_intervals(Tuple{Int, Int}[]) == Tuple{Int, Int}[]
        end

        Test.@testset "single interval is returned unchanged" begin
            Test.@test Mycelia._merge_intervals([(10, 20)]) == [(10, 20)]
        end

        Test.@testset "disjoint intervals are preserved and sorted" begin
            Test.@test Mycelia._merge_intervals([(30, 40), (10, 20)]) ==
                       [(10, 20), (30, 40)]
        end

        Test.@testset "overlapping intervals merge" begin
            Test.@test Mycelia._merge_intervals([(10, 25), (20, 30)]) == [(10, 30)]
        end

        Test.@testset "adjacent intervals (touching endpoints) merge under half-open semantics" begin
            # [10,20) and [20,30) touch at 20; half-open treats them as mergeable
            Test.@test Mycelia._merge_intervals([(10, 20), (20, 30)]) == [(10, 30)]
        end

        Test.@testset "nested intervals collapse" begin
            Test.@test Mycelia._merge_intervals([(10, 50), (20, 30)]) == [(10, 50)]
        end

        Test.@testset "many overlapping collapse to one" begin
            intervals = [(1, 5), (3, 9), (8, 12), (11, 15)]
            Test.@test Mycelia._merge_intervals(intervals) == [(1, 15)]
        end

        Test.@testset "mixed overlapping + disjoint" begin
            intervals = [(1, 5), (10, 15), (3, 4), (20, 25), (14, 18)]
            Test.@test Mycelia._merge_intervals(intervals) == [(1, 5), (10, 18), (20, 25)]
        end

        Test.@testset "duplicates collapse" begin
            Test.@test Mycelia._merge_intervals([(10, 20), (10, 20), (10, 20)]) ==
                       [(10, 20)]
        end

        Test.@testset "zero-length intervals are dropped" begin
            Test.@test Mycelia._merge_intervals([(10, 10)]) == Tuple{Int, Int}[]
            Test.@test Mycelia._merge_intervals([(10, 20), (15, 15)]) == [(10, 20)]
        end
    end

    Test.@testset "_fasta_total_length" begin
        Test.@testset "single-contig FASTA" begin
            path = tempname() * ".fna"
            write(path, ">seq1\nACGTACGTAC\n")
            Test.@test Mycelia._fasta_total_length(path) == 10
        end

        Test.@testset "multi-contig FASTA sums across sequences" begin
            path = tempname() * ".fna"
            # 10 + 7 + 3 = 20
            write(path, ">contig_1\nACGTACGTAC\n>contig_2\nAAAAAAA\n>contig_3\nTTT\n")
            Test.@test Mycelia._fasta_total_length(path) == 20
        end

        Test.@testset "multi-line sequence is handled correctly" begin
            path = tempname() * ".fna"
            write(path, ">seq1\nACGT\nACGT\nAC\n")  # 4 + 4 + 2 = 10
            Test.@test Mycelia._fasta_total_length(path) == 10
        end

        Test.@testset "empty sequence FASTA" begin
            path = tempname() * ".fna"
            write(path, ">empty\n\n")
            Test.@test Mycelia._fasta_total_length(path) == 0
        end
    end

    Test.@testset "_minimap_merged_coverage per-sequence" begin
        Test.@testset "secondary-alignment duplication collapses to unique coverage" begin
            # Simulate a minimap output where a repeat on query `contig_1` was
            # hit 3 times over the same coordinate range. Sum-of-lengths would
            # be 300; merged (unique) coverage is 100.
            df = DataFrames.DataFrame(
                Query = ["contig_1", "contig_1", "contig_1"],
                var"Query start" = [0, 0, 0],
                var"Query end" = [100, 100, 100],
                Target = ["ref_1", "ref_1", "ref_1"],
                var"Target start" = [0, 500, 1000],
                var"Target end" = [100, 600, 1100]
            )
            query_unique = Mycelia._minimap_merged_coverage(df, :query)
            ref_unique = Mycelia._minimap_merged_coverage(df, :target)
            # Query has 3 hits all at [0, 100) → 100 unique bases
            Test.@test query_unique == 100
            # Reference hits are disjoint → 300 unique bases (3 * 100)
            Test.@test ref_unique == 300
        end

        Test.@testset "multi-sequence query is summed across sequences" begin
            df = DataFrames.DataFrame(
                Query = ["contig_1", "contig_2", "contig_1"],
                var"Query start" = [0, 0, 50],
                var"Query end" = [100, 200, 150],
                Target = ["ref", "ref", "ref"],
                var"Target start" = [0, 100, 200],
                var"Target end" = [100, 300, 300]
            )
            # contig_1: [0,100) ∪ [50,150) = [0,150) → 150 bases
            # contig_2: [0,200) → 200 bases
            # total: 350
            Test.@test Mycelia._minimap_merged_coverage(df, :query) == 350
        end

        Test.@testset "empty DataFrame returns 0" begin
            df = DataFrames.DataFrame(
                Query = String[],
                var"Query start" = Int[],
                var"Query end" = Int[],
                Target = String[],
                var"Target start" = Int[],
                var"Target end" = Int[]
            )
            Test.@test Mycelia._minimap_merged_coverage(df, :query) == 0
            Test.@test Mycelia._minimap_merged_coverage(df, :target) == 0
        end

        Test.@testset "invalid `which` symbol errors with informative message" begin
            df = DataFrames.DataFrame(
                Query = ["c"], var"Query start" = [0], var"Query end" = [10],
                Target = ["r"], var"Target start" = [0], var"Target end" = [10]
            )
            try
                Mycelia._minimap_merged_coverage(df, :nonsense)
                Test.@test false  # unreachable
            catch e
                Test.@test occursin("which", String(sprint(showerror, e)))
            end
        end
    end
end
