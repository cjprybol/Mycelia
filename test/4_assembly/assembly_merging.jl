# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/assembly_merging.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/assembly_merging.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# TODO: Add more tests for merging assemblies, including edge cases and larger graphs.

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import FASTX

Test.@testset "Assembly Merging Tests" begin


    Test.@testset "1. merge_fasta_files" begin
        fasta1_path = joinpath(tempdir(), "merge_test1.fa")
        fasta2_path = joinpath(tempdir(), "merge_test2.fa")
        merged_path = joinpath(tempdir(), "merged_actual.fa")

        records1 = [FASTX.FASTA.Record("seq1", "ATCG")]
        records2 = [FASTX.FASTA.Record("seq2", "GCTA")]

        Mycelia.write_fasta(records=records1, outfile=fasta1_path)
        Mycelia.write_fasta(records=records2, outfile=fasta2_path)

        Mycelia.merge_fasta_files(fasta_files=[fasta1_path, fasta2_path], fasta_file=merged_path)

        Test.@test isfile(merged_path)

        merged_records = [r for r in FASTX.FASTA.Reader(open(merged_path))]
        Test.@test length(merged_records) == 2
        Test.@test FASTX.FASTA.identifier(merged_records[1]) == "merge_test1__seq1"
        Test.@test FASTX.FASTA.sequence(merged_records[1]) == "ATCG"
        Test.@test FASTX.FASTA.identifier(merged_records[2]) == "merge_test2__seq2"
        Test.@test FASTX.FASTA.sequence(merged_records[2]) == "GCTA"

        println("✓ merge_fasta_files successfully concatenated FASTA files.")
    end

    # Test.@testset "2. run_quickmerge" begin
    #     # QuickMerge is an external tool, so we test that Mycelia can call it correctly.
    #     # We'll create two overlapping assemblies.
    #     assembly1_path = joinpath(tempdir(), "assembly1.fa")
    #     assembly2_path = joinpath(tempdir(), "assembly2.fa")
    #     outdir = joinpath(tempdir(), "quickmerge_out")
    #     mkpath(outdir)

    #     # Assembly 1 has a contig that extends into Assembly 2's contig
    #     records1 = [FASTX.FASTA.Record("contig1_part1", "ATCGATCGATCG" * "GATTACA")]
    #     # Assembly 2 has an overlapping contig
    #     records2 = [FASTX.FASTA.Record("contig2_part2", "GATTACA" * "CTAGCTAG")]

    #     Mycelia.write_fasta(records=records1, outfile=assembly1_path)
    #     Mycelia.write_fasta(records=records2, outfile=assembly2_path)

    #     # Try to run QuickMerge - the function will automatically install it if needed
    #     try
    #         result = Mycelia.run_quickmerge(assembly1_path, assembly2_path, outdir=outdir)

    #         Test.@test isfile(result.merged_assembly)
            
    #         # Check if the merged assembly contains a merged contig.
    #         # The exact output of QuickMerge can be complex, so we'll just check
    #         # that the output file is non-empty and seems reasonable.
    #         merged_records = [r for r in FASTX.FASTA.Reader(open(result.merged_assembly))]
    #         Test.@test !isempty(merged_records)
            
    #         # A simple check: the merged sequence should be longer than either of the inputs
    #         # that could have been merged.
    #         merged_seq_length = maximum(length(FASTX.FASTA.sequence(r)) for r in merged_records)
    #         Test.@test merged_seq_length > length("ATCGATCGATCGGATTACA")
    #         Test.@test merged_seq_length > length("GATTACACTAGCTAG")
    #         # The merged length should be the sum of the two parts minus the overlap
    #         Test.@test merged_seq_length == length("ATCGATCGATCGGATTACACTAGCTAG")

    #         println("✓ run_quickmerge successfully executed and produced a merged assembly.")
    #     catch e
    #         Test.@test_skip true  # Skip test if QuickMerge fails for any reason
    #         @warn "run_quickmerge test failed - skipping: $e"
    #         println("⚠ run_quickmerge test skipped due to error.")
    #     end
    # end
end
