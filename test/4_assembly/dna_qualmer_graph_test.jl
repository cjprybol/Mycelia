# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/dna_qualmer_graph_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/dna_qualmer_graph_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# DNA Qualmer Graph Construction Test
# Tests for DNA qualmer (quality-aware k-mer) graph construction in singlestrand mode.
# Validates that DNA sequences with quality scores are properly processed.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "DNA Qualmer Graph - Singlestrand" begin

    Test.@testset "DNA Qualmer Graph Construction - Single Read" begin
        # Create FASTQ record with quality scores
        seq = "ATGCG"
        qual = [30, 35, 32, 28, 40]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 3 unique DNA k-mers
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 3
        Test.@test Mycelia.Rhizomorph.edge_count(graph) == 2

        # Verify DNA k-mers are present
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("ATG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("TGC"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("GCG"))
    end

    Test.@testset "DNA Qualmer Graph - Quality Evidence Tracking" begin
        seq = "ATGCG"
        qual = [30, 35, 32, 28, 40]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        Test.@test !isnothing(vertex_data)
        Test.@test haskey(vertex_data.evidence, "test")
        Test.@test haskey(vertex_data.evidence["test"], "read_001")

        # Check that evidence includes quality scores
        evidence_set = vertex_data.evidence["test"]["read_001"]
        Test.@test length(evidence_set) == 1
        evidence = first(evidence_set)

        # QualityEvidenceEntry should have quality_scores field
        Test.@test evidence isa Mycelia.Rhizomorph.QualityEvidenceEntry
        Test.@test evidence.position == 1
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
        # Quality scores are stored as Phred+33 encoded (30+33=63, 35+33=68, 32+33=65)
        Test.@test evidence.quality_scores == UInt8[63, 68, 65]
    end

    Test.@testset "DNA Qualmer Graph - Multiple Reads Same Kmer" begin
        # Two reads with same k-mer but different quality scores
        seq1 = "ATGC"
        qual1 = [30, 30, 30, 30]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq1, qual_str1)

        seq2 = "ATGC"
        qual2 = [40, 40, 40, 40]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual_str2)

        records = [record1, record2]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Same k-mers, but observed with different qualities
        kmer = Kmers.DNAKmer{3}("ATG")
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, kmer)
        Test.@test count == 2

        # Check that both quality profiles are stored
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)
        Test.@test length(vertex_data.evidence["test"]) == 2  # Two observations
    end

    Test.@testset "DNA Qualmer Graph - Edge Quality Evidence" begin
        seq = "ATGC"
        qual = [30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.DNAKmer{3}("ATG")
        dst = Kmers.DNAKmer{3}("TGC")

        Test.@test Mycelia.Rhizomorph.has_edge(graph, src, dst)

        edge_data = Mycelia.Rhizomorph.get_edge_data(graph, src, dst)
        Test.@test !isnothing(edge_data)
        Test.@test haskey(edge_data.evidence, "test")

        # Check edge quality evidence
        edge_evidence = edge_data.evidence["test"]["read_001"]
        Test.@test length(edge_evidence) == 1
        evidence = first(edge_evidence)

        Test.@test evidence isa Mycelia.Rhizomorph.EdgeQualityEvidenceEntry
        Test.@test evidence.from_position == 1
        Test.@test evidence.to_position == 2
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
        # Quality scores are Phred+33 encoded
        Test.@test evidence.from_quality == UInt8[63, 68, 65]  # 30+33, 35+33, 32+33
        Test.@test evidence.to_quality == UInt8[68, 65, 61]    # 35+33, 32+33, 28+33
    end

    Test.@testset "DNA Qualmer Graph - Low Quality Bases" begin
        # Test that low quality bases are still included (not filtered at graph construction)
        seq = "ATGC"
        qual = [5, 5, 5, 5]  # Very low quality
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # K-mers should still be added (filtering happens later if needed)
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2

        # Quality scores should be preserved (Phred+33 encoded)
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)
        evidence = first(vertex_data.evidence["test"]["read_001"])
        Test.@test evidence.quality_scores == UInt8[38, 38, 38]  # 5+33
    end

    Test.@testset "DNA Qualmer Graph - Ambiguous Bases Skipped" begin
        # N bases should be skipped by UnambiguousDNAMers
        seq = "ATGNNCG"
        qual = [30, 30, 30, 30, 30, 30, 30]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should only have ATG (NNG, NCG are invalid)
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("ATG"))

        # K-mers containing N should not exist
        vertices = Mycelia.Rhizomorph.get_all_vertices(graph)
        for v in vertices
            seq_str = string(v)
            Test.@test !occursin('N', seq_str)
        end
    end

    Test.@testset "DNA Qualmer Graph - Strand Specific" begin
        # Forward and reverse DNA k-mers should be separate
        seq1 = "ATG"
        qual1 = [30, 30, 30]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_fwd", seq1, qual_str1)

        seq2 = "CAT"  # Reverse complement of ATG
        qual2 = [30, 30, 30]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_rev", seq2, qual_str2)

        records = [record1, record2]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 2 separate vertices (strand-specific)
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("ATG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("CAT"))
    end

    Test.@testset "DNA Qualmer Graph - Quality Score Range" begin
        # Test various quality score ranges
        seq = "ATGC"

        # Test with Phred+33 scores (typical Illumina)
        qual = [0, 20, 40, 60]  # Full range
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)
        evidence = first(vertex_data.evidence["test"]["read_001"])

        # Quality scores should be in UInt8 range (Phred+33 encoded)
        Test.@test all(q -> 0 <= q <= 255, evidence.quality_scores)
        Test.@test evidence.quality_scores == UInt8[33, 53, 73]  # 0+33, 20+33, 40+33
    end

    Test.@testset "DNA Qualmer Graph - Query Functions Work" begin
        seq = "ATGCGATCG"
        qual = [30 for _ in 1:length(seq)]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Test that graph query functions work on qualmer graphs
        sources = Mycelia.Rhizomorph.get_all_sources(graph)
        sinks = Mycelia.Rhizomorph.get_all_sinks(graph)

        Test.@test Kmers.DNAKmer{3}("ATG") in sources
        Test.@test length(sinks) > 0

        # Test neighbor queries
        kmer = Kmers.DNAKmer{3}("TGC")
        incoming = Mycelia.Rhizomorph.get_incoming_neighbors(graph, kmer)
        outgoing = Mycelia.Rhizomorph.get_outgoing_neighbors(graph, kmer)

        Test.@test length(incoming) > 0
        Test.@test length(outgoing) > 0
    end

    Test.@testset "DNA Qualmer Graph - Path Assembly" begin
        seq = "ATGCG"
        qual = [30 for _ in 1:length(seq)]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Get path
        path = [
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("TGC"),
            Kmers.DNAKmer{3}("GCG")
        ]

        # Assemble sequence from path (works same as k-mer graphs)
        sequence = Mycelia.Rhizomorph.assemble_path_sequence(path)

        Test.@test sequence == BioSequences.LongDNA{4}("ATGCG")
        Test.@test sequence isa BioSequences.LongDNA{4}
    end

    Test.@testset "DNA Qualmer Graph - Multiple Datasets" begin
        # Test evidence tracking across multiple datasets
        seq = "ATGC"
        qual = [30, 30, 30, 30]
        qual_str = String([Char(q + 33) for q in qual])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str)
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str)

        # Build graph with dataset 1
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1], 3; dataset_id="dataset_01")

        # Add observations from dataset 2
        Mycelia.Rhizomorph.add_observations_to_graph!(graph, [record2], 3; dataset_id="dataset_02")

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        # Should have evidence from both datasets
        Test.@test haskey(vertex_data.evidence, "dataset_01")
        Test.@test haskey(vertex_data.evidence, "dataset_02")

        # Count observations per dataset
        count_d1 = Mycelia.Rhizomorph.count_dataset_observations(vertex_data, "dataset_01")
        count_d2 = Mycelia.Rhizomorph.count_dataset_observations(vertex_data, "dataset_02")

        Test.@test count_d1 == 1
        Test.@test count_d2 == 1
    end
end

println("✓ DNA qualmer graph tests completed")
