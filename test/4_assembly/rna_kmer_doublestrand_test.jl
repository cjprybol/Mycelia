# RNA K-mer Doublestrand (Canonical) Graph Test
#
# Tests for canonical RNA k-mer graph construction where forward and reverse
# complement k-mers are merged into canonical representation.
#
# Run with: julia --project=. test/4_assembly/rna_kmer_doublestrand_test.jl

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "RNA K-mer Doublestrand Graph" begin

    Test.@testset "Canonical K-mer Merging - Basic" begin
        # RNA sequence
        seq1 = "AUGCAU"
        record = FASTX.FASTA.Record("read_001", seq1)

        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 3; dataset_id="test")

        # In "AUGCAU": AUG, UGC, GCA, CAU
        # Canonical: AUG (for AUG and CAU), GCA (for UGC and GCA)
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2
    end

    Test.@testset "Evidence Merging from Both Strands" begin
        # Forward read: AUGC
        record1 = FASTX.FASTA.Record("read_fwd", "AUGC")
        # Reverse read: GCAU (RC of AUGC)
        record2 = FASTX.FASTA.Record("read_rev", "GCAU")

        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record1, record2], 3; dataset_id="test")

        # Both reads should contribute to canonical k-mers
        # Check that observation count reflects both strands
        canon_aug = BioSequences.canonical(Kmers.RNAKmer{3}("AUG"))
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, canon_aug)

        # Should have 2 observations (one from each read)
        Test.@test count == 2
    end

    Test.@testset "Strand Evidence is Preserved" begin
        record = FASTX.FASTA.Record("read_001", "AUGC")
        graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 3; dataset_id="test")

        canon_aug = BioSequences.canonical(Kmers.RNAKmer{3}("AUG"))
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, canon_aug)

        # Should have evidence with Forward strand
        evidence_set = vertex_data.evidence["test"]["read_001"]
        Test.@test length(evidence_set) >= 1

        # At least one evidence entry should be Forward
        has_forward = any(e -> e.strand == Mycelia.Rhizomorph.Forward, evidence_set)
        Test.@test has_forward
    end

    Test.@testset "Graph Query Functions Work on Doublestrand" begin
        record = FASTX.FASTA.Record("read_001", "AUGCGAUCG")
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

println("✓ RNA doublestrand k-mer graph tests completed")