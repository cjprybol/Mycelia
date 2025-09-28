# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/comprehensive_fixes_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/comprehensive_fixes_tests.jl", "test/4_assembly", execute=false)'
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

# Test comprehensive fixes for canonicalization consistency and FASTQ compatibility
Test.@testset "Canonicalization and FASTQ Compatibility Regression Tests" begin

    Test.@testset "DoubleStrand DNA Mode Fix" begin
        reference_seq = BioSequences.dna"ATCGATCGATCG"
        reads = [FASTX.FASTA.Record("read1", reference_seq)]

        kmer_type = Kmers.DNAKmer{5}
        graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.DoubleStrand)

        Test.@test !isempty(MetaGraphsNext.labels(graph))

        # Check that coverage exists (the original failing test)
        has_coverage = false
        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
            if !isempty(strand_orientations)
                has_coverage = true
            end
        end
        Test.@test has_coverage
    end

    Test.@testset "SingleStrand RNA Mode Fix" begin
        reference_seq = BioSequences.rna"AUCGAUCGAUCG"
        reads = [FASTX.FASTA.Record("read1", reference_seq)]

        kmer_type = Kmers.RNAKmer{5}
        graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.SingleStrand)

        Test.@test !isempty(MetaGraphsNext.labels(graph))

        # Check that coverage exists (was failing due to k-mer mismatch)
        has_coverage = false
        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            strand_orientations = [strand for (obs_id, pos, strand) in vertex_data.coverage]
            if !isempty(strand_orientations)
                has_coverage = true
                # In SingleStrand mode, all should be Forward
                Test.@test all(s == Mycelia.Forward for s in strand_orientations)
            end
        end
        Test.@test has_coverage
    end

    Test.@testset "Amino Acid FASTQ Fix" begin
        reference_seq = BioSequences.aa"ACDLMVFGHY"

        # Test that observe doesn't produce invalid characters
        observed_seq, quality_scores = Mycelia.observe(reference_seq, error_rate=0.1)
        observed_seq_str = string(observed_seq)

        # No termination characters
        Test.@test !occursin('*', observed_seq_str)

        # Can create FASTQ record (was failing before)
        quality_string = String([Char(q + 33) for q in quality_scores])
        record = FASTX.FASTQ.Record("test", observed_seq_str, quality_string)
        Test.@test record isa FASTX.FASTQ.Record

        # Can build k-mer graph with amino acids
        kmer_type = Kmers.AAKmer{3}
        reads = [FASTX.FASTA.Record("read1", reference_seq)]
        graph = Mycelia.build_kmer_graph_next(kmer_type, reads; graph_mode=Mycelia.SingleStrand)
        Test.@test !isempty(MetaGraphsNext.labels(graph))
    end
end