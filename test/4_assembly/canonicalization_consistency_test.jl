# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/canonicalization_consistency_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/canonicalization_consistency_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "K-mer Canonicalization Consistency Tests" begin
    ## Test that validates the fix for canonical k-mer consistency between
    ## graph construction and path extraction (the main issue we resolved)

    Test.@testset "DoubleStrand DNA Canonicalization Consistency" begin
        # Create a DNA sequence that will have different canonical forms
        # if lexicographic vs BioSequences.canonical() are used
        reference_seq = BioSequences.dna"ATCGTTTT"  # Forward
        reverse_comp = BioSequences.reverse_complement(reference_seq)  # AAAACGAT

        # Test with both forward and reverse complement sequences
        reads = [
            FASTX.FASTA.Record("forward", reference_seq),
            FASTX.FASTA.Record("reverse", reverse_comp)
        ]

        kmer_type = Kmers.DNAKmer{5}
        graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.DoubleStrand)

        Test.@test !isempty(MetaGraphsNext.labels(graph))

        # Test that we have proper coverage and strand orientations
        # This validates the canonicalization consistency fix
        total_coverage_entries = 0
        has_both_orientations = false

        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
            total_coverage_entries += length(strand_orientations)

            # Check if we see both orientations (should happen with reverse complements)
            if Mycelia.Forward in strand_orientations && Mycelia.Reverse in strand_orientations
                has_both_orientations = true
            end
        end

        Test.@test total_coverage_entries > 0
        Test.@test has_both_orientations  # Should see both orientations with forward and reverse sequences
    end

    Test.@testset "SingleStrand RNA Non-Canonicalization Consistency" begin
        # Test that SingleStrand mode uses k-mers as-is without canonicalization
        reference_seq = BioSequences.rna"AUCGUUUU"
        reads = [FASTX.FASTA.Record("rna", reference_seq)]

        kmer_type = Kmers.RNAKmer{5}
        graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.SingleStrand)

        Test.@test !isempty(MetaGraphsNext.labels(graph))

        # Test that all coverage entries use Forward orientation in SingleStrand mode
        # This validates the SingleStrand non-canonicalization fix
        total_coverage_entries = 0
        all_forward = true

        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            for (obs_id, pos, strand) in vertex_data.coverage
                total_coverage_entries += 1
                if strand != Mycelia.Forward
                    all_forward = false
                end
            end
        end

        Test.@test total_coverage_entries > 0
        Test.@test all_forward  # All orientations should be Forward in SingleStrand mode
    end
end