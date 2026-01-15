# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rna_kmer_graph_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rna_kmer_graph_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# RNA K-mer Graph Construction Test
# Tests for RNA k-mer graph construction in singlestrand mode.
# Validates that RNA sequences are properly processed with UnambiguousRNAMers.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "RNA K-mer Graph - Singlestrand" begin

    Test.@testset "RNA K-mer Extraction" begin
        # RNA sequence with U instead of T
        rna_seq = BioSequences.LongRNA{4}("AUGCGAUCG")
        k = 3

        # Extract k-mers using UnambiguousRNAMers
        kmers = [kmer for (kmer, pos) in Kmers.UnambiguousRNAMers{k}(rna_seq)]

        Test.@test length(kmers) == 7  # 9 - 3 + 1 = 7
        Test.@test kmers[1] == Kmers.RNAKmer{3}("AUG")
        Test.@test kmers[2] == Kmers.RNAKmer{3}("UGC")
        Test.@test kmers[end] == Kmers.RNAKmer{3}("UCG")
    end

    Test.@testset "RNA Graph Construction - Single Read" begin
        records = [
            FASTX.FASTA.Record("read_001", "AUGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 3 unique RNA k-mers
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 3
        Test.@test Mycelia.Rhizomorph.edge_count(graph) == 2

        # Verify RNA k-mers are present
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.RNAKmer{3}("AUG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.RNAKmer{3}("UGC"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.RNAKmer{3}("GCG"))
    end

    Test.@testset "RNA Graph Construction - Multiple Reads" begin
        records = [
            FASTX.FASTA.Record("read_001", "AUGCG"),
            FASTX.FASTA.Record("read_002", "AUGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Same k-mers, but observed twice
        kmer = Kmers.RNAKmer{3}("AUG")
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, kmer)
        Test.@test count == 2
    end

    Test.@testset "RNA Graph - Evidence Tracking" begin
        records = [
            FASTX.FASTA.Record("read_001", "AUGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="rna_test")

        kmer = Kmers.RNAKmer{3}("AUG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        Test.@test !isnothing(vertex_data)
        Test.@test haskey(vertex_data.evidence, "rna_test")
        Test.@test haskey(vertex_data.evidence["rna_test"], "read_001")

        # Check evidence entry
        evidence_set = vertex_data.evidence["rna_test"]["read_001"]
        Test.@test length(evidence_set) == 1
        evidence = first(evidence_set)
        Test.@test evidence.position == 1
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "RNA Graph - Strand Specific" begin
        # Forward and reverse RNA k-mers should be separate
        records = [
            FASTX.FASTA.Record("read_fwd", "AUG"),
            FASTX.FASTA.Record("read_rev", "CAU")  # Reverse complement of AUG
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 2 separate vertices (strand-specific)
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.RNAKmer{3}("AUG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.RNAKmer{3}("CAU"))
    end

    Test.@testset "RNA Graph - Ambiguous Bases Skipped" begin
        # N bases should be skipped by UnambiguousRNAMers
        records = [
            FASTX.FASTA.Record("read_001", "AUGNNCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should only have AUG and NCG is invalid (contains N)
        # Actually, NCG won't be extracted at all
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.RNAKmer{3}("AUG"))
        # K-mers containing N should not exist
        vertices = Mycelia.Rhizomorph.get_all_vertices(graph)
        for v in vertices
            seq = string(v)
            Test.@test !occursin('N', seq)
        end
    end

    Test.@testset "RNA Graph - Edge Evidence" begin
        records = [
            FASTX.FASTA.Record("read_001", "AUGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.RNAKmer{3}("AUG")
        dst = Kmers.RNAKmer{3}("UGC")

        Test.@test Mycelia.Rhizomorph.has_edge(graph, src, dst)

        edge_data = Mycelia.Rhizomorph.get_edge_data(graph, src, dst)
        Test.@test !isnothing(edge_data)
        Test.@test haskey(edge_data.evidence, "test")

        # Check edge evidence
        edge_evidence = edge_data.evidence["test"]["read_001"]
        Test.@test length(edge_evidence) == 1
        evidence = first(edge_evidence)
        Test.@test evidence.from_position == 1
        Test.@test evidence.to_position == 2
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "RNA Graph - Path Assembly" begin
        records = [
            FASTX.FASTA.Record("read_001", "AUGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Get path
        path = [
            Kmers.RNAKmer{3}("AUG"),
            Kmers.RNAKmer{3}("UGC"),
            Kmers.RNAKmer{3}("GCG")
        ]

        # Assemble sequence from path
        sequence = Mycelia.Rhizomorph.assemble_path_sequence(path)

        Test.@test sequence == BioSequences.LongRNA{4}("AUGCG")
        Test.@test sequence isa BioSequences.LongRNA{4}
    end

    Test.@testset "RNA Graph - Query Functions" begin
        records = [
            FASTX.FASTA.Record("read_001", "AUGCGAUC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Test source/sink detection
        sources = Mycelia.Rhizomorph.get_all_sources(graph)
        sinks = Mycelia.Rhizomorph.get_all_sinks(graph)

        Test.@test Kmers.RNAKmer{3}("AUG") in sources
        Test.@test Kmers.RNAKmer{3}("AUC") in sinks

        # Test neighbor queries
        kmer = Kmers.RNAKmer{3}("UGC")
        incoming = Mycelia.Rhizomorph.get_incoming_neighbors(graph, kmer)
        outgoing = Mycelia.Rhizomorph.get_outgoing_neighbors(graph, kmer)

        Test.@test Kmers.RNAKmer{3}("AUG") in incoming
        Test.@test Kmers.RNAKmer{3}("GCG") in outgoing
    end
end

println("✓ RNA k-mer graph tests completed")
