# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/dna_kmer_doublestrand_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/dna_kmer_doublestrand_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# DNA K-mer Doublestrand (Canonical) Graph Test
# Tests for canonical k-mer graph construction where forward and reverse
# complement k-mers are merged into canonical representation.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "DNA K-mer Doublestrand Graph" begin

    Test.@testset "Canonical K-mer Merging - Basic" begin
        # Read with ATG and its RC CAT
        seq1 = "ATGCAT"
        record = FASTX.FASTA.Record("read_001", seq1)

        # Use canonical builder for canonical behavior (merging forward and RC)
        graph = Mycelia.Rhizomorph.build_kmer_graph_canonical([record], 3; dataset_id="test")

        # Canonical form of ATG is ATG (ATG < CAT lexicographically)
        # Canonical form of CAT is ATG
        # So we should have fewer vertices than singlestrand

        # In "ATGCAT": ATG, TGC, GCA, CAT
        # Canonical: ATG (for ATG and CAT), GCA (for TGC and GCA)
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2
    end

    Test.@testset "Evidence Merging from Both Strands" begin
        # Forward read: ATGC
        record1 = FASTX.FASTA.Record("read_fwd", "ATGC")
        # Reverse read: GCAT (RC of ATGC)
        record2 = FASTX.FASTA.Record("read_rev", "GCAT")

        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record1, record2], 3; dataset_id="test")

        # Both reads should contribute to canonical k-mers
        # Check that observation count reflects both strands
        canon_atg = BioSequences.canonical(Kmers.DNAKmer{3}("ATG"))
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, canon_atg)

        # Should have 2 observations (one from each read)
        Test.@test count == 2
    end

    Test.@testset "Strand Evidence is Preserved" begin
        record = FASTX.FASTA.Record("read_001", "ATGC")
        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 3; dataset_id="test")

        canon_atg = BioSequences.canonical(Kmers.DNAKmer{3}("ATG"))
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, canon_atg)

        # Should have evidence with Forward strand
        evidence_set = vertex_data.evidence["test"]["read_001"]
        Test.@test length(evidence_set) >= 1

        # At least one evidence entry should be Forward
        has_forward = any(e -> e.strand == Mycelia.Rhizomorph.Forward, evidence_set)
        Test.@test has_forward
    end

    Test.@testset "Palindromic K-mers" begin
        # Palindromic k-mer: GCGC has RC = GCGC
        record = FASTX.FASTA.Record("read_001", "GCGCGC")
        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 4; dataset_id="test")

        kmer = Kmers.DNAKmer{4}("GCGC")
        canon = BioSequences.canonical(kmer)

        Test.@test Mycelia.Rhizomorph.has_vertex(graph, canon)
        # GCGCGC contains: GCGC (palindrome), CGCG (palindrome), GCGC (duplicate)
        # Both are palindromes so they are their own canonical forms
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2  # GCGC and CGCG
    end

    Test.@testset "Graph Query Functions Work on Doublestrand" begin
        record = FASTX.FASTA.Record("read_001", "ATGCGATCG")
        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 3; dataset_id="test")

        # Query functions should work
        sources = Mycelia.Rhizomorph.get_all_sources(graph)

        # In doublestrand graphs, RC edges can create cycles, so there may not be sinks
        # Just verify that query functions work without error
        Test.@test length(sources) > 0
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) > 0
    end

    Test.@testset "AA Sequences Rejected for Doublestrand" begin
        record = FASTX.FASTA.Record("protein_001", "MKVLW")

        # AA sequences don't have reverse complement
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 3)
    end
end

println("✓ DNA doublestrand k-mer graph tests completed")
