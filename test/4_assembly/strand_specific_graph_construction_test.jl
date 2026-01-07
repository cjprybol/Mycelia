# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/strand_specific_graph_construction_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/strand_specific_graph_construction_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Strand-Specific K-mer Graph Construction Test
# Tests for strand-specific k-mer de Bruijn graph construction as specified
# in the rhizomorph graph ecosystem planning document section 1.4.
# CRITICAL REQUIREMENTS:
# 1. Store k-mers AS OBSERVED (not canonical)
# 2. Strand-specific construction by default
# 3. Evidence tracked with dataset_id and observation_id
# 4. Graph uses MetaGraphsNext with proper vertex/edge data

import Test
import Mycelia
import Kmers
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "Strand-Specific Graph Construction - Basic K-mer Extraction" begin

    Test.@testset "Extract k-mers from simple sequence" begin
        # Simple DNA sequence
        seq = BioSequences.LongDNA{4}("ATCGATCG")
        k = 3

        # Extract k-mers using UnambiguousDNAMers (strand-specific)
        kmers = [kmer for (kmer, pos) in Kmers.UnambiguousDNAMers{k}(seq)]

        Test.@test length(kmers) == length(seq) - k + 1
        Test.@test length(kmers) == 6  # ATCGATCG -> 6 3-mers

        # Verify actual k-mers
        Test.@test String(kmers[1]) == "ATC"
        Test.@test String(kmers[2]) == "TCG"
        Test.@test String(kmers[3]) == "CGA"
        Test.@test String(kmers[4]) == "GAT"
        Test.@test String(kmers[5]) == "ATC"  # Repeated k-mer
        Test.@test String(kmers[6]) == "TCG"  # Repeated k-mer
    end

    Test.@testset "K-mer extraction tracks positions" begin
        seq = BioSequences.LongDNA{4}("ATCGATCG")
        k = 3

        kmers_with_positions = [(kmer, pos) for (kmer, pos) in Kmers.UnambiguousDNAMers{k}(seq)]

        # Each k-mer should know its position (1-indexed)
        Test.@test kmers_with_positions[1][1] == Kmers.DNAKmer{3}("ATC")
        Test.@test kmers_with_positions[1][2] == 1

        Test.@test kmers_with_positions[2][1] == Kmers.DNAKmer{3}("TCG")
        Test.@test kmers_with_positions[2][2] == 2

        Test.@test kmers_with_positions[3][1] == Kmers.DNAKmer{3}("CGA")
        Test.@test kmers_with_positions[3][2] == 3
    end

    Test.@testset "Extract k-mers stores OBSERVED not canonical" begin
        # Sequence where some k-mers are their own reverse complement
        seq = BioSequences.LongDNA{4}("ATGCAT")
        k = 3

        kmers = [kmer for (kmer, pos) in Kmers.UnambiguousDNAMers{k}(seq)]

        # Should get ATG, TGC, GCA, CAT
        Test.@test String(kmers[1]) == "ATG"
        Test.@test String(kmers[2]) == "TGC"
        Test.@test String(kmers[3]) == "GCA"
        Test.@test String(kmers[4]) == "CAT"

        # ATG RC is CAT - but we store both as observed
        Test.@test Kmers.DNAKmer{3}("ATG") != BioSequences.reverse_complement(Kmers.DNAKmer{3}("ATG"))
        Test.@test String(BioSequences.reverse_complement(Kmers.DNAKmer{3}("ATG"))) == "CAT"
    end

    Test.@testset "UnambiguousDNAMers skips ambiguous bases" begin
        # Sequence with ambiguous base (N)
        seq = BioSequences.LongDNA{4}("ATNGCAT")  # N at position 3
        k = 3

        kmers_with_positions = [(kmer, pos) for (kmer, pos) in Kmers.UnambiguousDNAMers{k}(seq)]

        # Should skip: ATN (pos 1), TNG (pos 2), NGC (pos 3)
        # Should keep: GCA (pos 5), CAT (pos 6)
        # Note: Position 4 would be NCA which contains N

        # Only unambiguous k-mers are returned
        Test.@test length(kmers_with_positions) <= 4

        # None should contain ambiguous bases
        for (kmer, pos) in kmers_with_positions
            kmer_str = String(kmer)
            Test.@test !occursin('N', kmer_str)
        end
    end
end

Test.@testset "Strand-Specific Graph Construction - Single Read" begin

    Test.@testset "Build graph from single short read" begin
        # Create a simple FASTQ record
        seq = BioSequences.LongDNA{4}("ATCGATCG")
        qual = fill(UInt8(30), length(seq))
        record = FASTX.FASTQ.Record("read_001", seq, qual)

        k = 3
        dataset_id = "test_dataset"

        # Build strand-specific k-mer graph
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [record], k; dataset_id=dataset_id
        )

        # Graph should have vertices for observed k-mers
        # ATCGATCG gives: ATC(x2), TCG(x2), CGA, GAT
        # Unique k-mers: ATC, TCG, CGA, GAT (4 unique)
        Test.@test MetaGraphsNext.nv(graph) == 4

        # Should have edges connecting adjacent k-mers
        # ATC->TCG, TCG->CGA, CGA->GAT, GAT->ATC, ATC->TCG
        Test.@test MetaGraphsNext.ne(graph) >= 4
    end

    Test.@testset "Graph stores k-mers as observed" begin
        seq = BioSequences.LongDNA{4}("ATGCAT")
        qual = fill(UInt8(30), length(seq))
        record = FASTX.FASTQ.Record("read_001", seq, qual)

        k = 3
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], k)

        # Should have ATG, TGC, GCA, CAT as separate vertices
        # (not merged into canonical pairs)
        Test.@test MetaGraphsNext.nv(graph) == 4

        # Verify ATG and CAT are separate vertices
        atg_kmer = Kmers.DNAKmer{3}("ATG")
        cat_kmer = Kmers.DNAKmer{3}("CAT")

        Test.@test haskey(graph, atg_kmer)
        Test.@test haskey(graph, cat_kmer)
        Test.@test atg_kmer != cat_kmer
    end

    Test.@testset "Graph vertex data stores evidence" begin
        seq = BioSequences.LongDNA{4}("ATCGATCG")
        qual = fill(UInt8(30), length(seq))
        record = FASTX.FASTQ.Record("read_001", seq, qual)

        k = 3
        dataset_id = "test_dataset"

        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [record], k; dataset_id=dataset_id
        )

        # Get vertex data for first k-mer (ATC)
        atc_kmer = Kmers.DNAKmer{3}("ATC")
        vertex_data = graph[atc_kmer]

        # Vertex data should be KmerVertexData
        Test.@test isa(vertex_data, Mycelia.Rhizomorph.KmerVertexData)

        # Should have evidence from our read
        Test.@test haskey(vertex_data.evidence, dataset_id)
        Test.@test haskey(vertex_data.evidence[dataset_id], "read_001")

        # ATC appears at positions 1 and 5
        evidence_set = vertex_data.evidence[dataset_id]["read_001"]
        Test.@test length(evidence_set) == 2
    end
end

Test.@testset "Strand-Specific Graph Construction - Multiple Reads" begin

    Test.@testset "Build graph from multiple reads" begin
        # Two reads observing same k-mers
        seq1 = BioSequences.LongDNA{4}("ATCGAT")
        seq2 = BioSequences.LongDNA{4}("TCGATC")
        qual = fill(UInt8(30), 6)

        record1 = FASTX.FASTQ.Record("read_001", seq1, qual)
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual)

        k = 3
        dataset_id = "test_dataset"

        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [record1, record2], k; dataset_id=dataset_id
        )

        # Should have combined evidence from both reads
        tcg_kmer = Kmers.DNAKmer{3}("TCG")
        vertex_data = graph[tcg_kmer]

        # Evidence from both reads
        Test.@test haskey(vertex_data.evidence[dataset_id], "read_001")
        Test.@test haskey(vertex_data.evidence[dataset_id], "read_002")
    end

    Test.@testset "Different reads contribute different positions" begin
        # Same k-mer at different positions in different reads
        seq1 = BioSequences.LongDNA{4}("ATCGGG")  # ATC at pos 1
        seq2 = BioSequences.LongDNA{4}("GGATCG")  # ATC at pos 3
        qual = fill(UInt8(30), 6)

        record1 = FASTX.FASTQ.Record("read_001", seq1, qual)
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual)

        k = 3
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record1, record2], k)

        atc_kmer = Kmers.DNAKmer{3}("ATC")
        vertex_data = graph[atc_kmer]

        # read_001 should have evidence at position 1
        read1_evidence = Mycelia.Rhizomorph.get_observation_evidence(
            vertex_data, "dataset_01", "read_001"
        )
        Test.@test any(e -> e.position == 1, read1_evidence)

        # read_002 should have evidence at position 3
        read2_evidence = Mycelia.Rhizomorph.get_observation_evidence(
            vertex_data, "dataset_01", "read_002"
        )
        Test.@test any(e -> e.position == 3, read2_evidence)
    end
end

Test.@testset "Strand-Specific Graph Construction - Forward vs Reverse" begin

    Test.@testset "Forward and reverse complement k-mers are separate" begin
        # Read on forward strand
        seq_fwd = BioSequences.LongDNA{4}("ATGCCC")
        # Read on reverse strand (RC of forward)
        seq_rev = BioSequences.reverse_complement(seq_fwd)  # GGGCAT

        qual = fill(UInt8(30), 6)
        record_fwd = FASTX.FASTQ.Record("read_fwd", seq_fwd, qual)
        record_rev = FASTX.FASTQ.Record("read_rev", seq_rev, qual)

        k = 3
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [record_fwd, record_rev], k
        )

        # ATG from forward read
        atg_kmer = Kmers.DNAKmer{3}("ATG")
        # CAT from reverse read (RC of ATG)
        cat_kmer = Kmers.DNAKmer{3}("CAT")

        # Both should exist as separate vertices
        Test.@test haskey(graph, atg_kmer)
        Test.@test haskey(graph, cat_kmer)
        Test.@test atg_kmer != cat_kmer

        # They should have different evidence
        atg_data = graph[atg_kmer]
        cat_data = graph[cat_kmer]

        Test.@test haskey(atg_data.evidence["dataset_01"], "read_fwd")
        Test.@test haskey(cat_data.evidence["dataset_01"], "read_rev")
    end

    Test.@testset "Evidence tracks strand orientation" begin
        seq = BioSequences.LongDNA{4}("ATCGAT")
        qual = fill(UInt8(30), 6)
        record = FASTX.FASTQ.Record("read_001", seq, qual)

        k = 3
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], k)

        # Get evidence for ATC k-mer
        atc_kmer = Kmers.DNAKmer{3}("ATC")
        vertex_data = graph[atc_kmer]
        evidence_set = vertex_data.evidence["dataset_01"]["read_001"]

        # All evidence should be Forward strand (read as-is)
        Test.@test all(e -> e.strand == Mycelia.Rhizomorph.Forward, evidence_set)
    end
end

Test.@testset "Strand-Specific Graph Construction - Edge Creation" begin

    Test.@testset "Edges connect overlapping k-mers" begin
        seq = BioSequences.LongDNA{4}("ATCG")
        qual = fill(UInt8(30), 4)
        record = FASTX.FASTQ.Record("read_001", seq, qual)

        k = 3
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], k)

        # ATCG gives k-mers: ATC, TCG
        # Should have edge ATC -> TCG
        atc_kmer = Kmers.DNAKmer{3}("ATC")
        tcg_kmer = Kmers.DNAKmer{3}("TCG")

        Test.@test haskey(graph, atc_kmer, tcg_kmer)
    end

    Test.@testset "Edge data stores evidence" begin
        seq = BioSequences.LongDNA{4}("ATCG")
        qual = fill(UInt8(30), 4)
        record = FASTX.FASTQ.Record("read_001", seq, qual)

        k = 3
        dataset_id = "test_dataset"
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [record], k; dataset_id=dataset_id
        )

        atc_kmer = Kmers.DNAKmer{3}("ATC")
        tcg_kmer = Kmers.DNAKmer{3}("TCG")

        edge_data = graph[atc_kmer, tcg_kmer]

        # Edge data should be KmerEdgeData
        Test.@test isa(edge_data, Mycelia.Rhizomorph.KmerEdgeData)

        # Should have evidence from our read
        Test.@test haskey(edge_data.evidence, dataset_id)
        Test.@test haskey(edge_data.evidence[dataset_id], "read_001")

        # Edge evidence tracks transition positions
        edge_evidence = edge_data.evidence[dataset_id]["read_001"]
        Test.@test length(edge_evidence) >= 1

        # Evidence should record transition from position 1 to position 2
        Test.@test any(e -> e.from_position == 1 && e.to_position == 2, edge_evidence)
    end

    Test.@testset "Multiple reads strengthen edge evidence" begin
        seq1 = BioSequences.LongDNA{4}("ATCG")
        seq2 = BioSequences.LongDNA{4}("ATCG")
        qual = fill(UInt8(30), 4)

        record1 = FASTX.FASTQ.Record("read_001", seq1, qual)
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual)

        k = 3
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record1, record2], k)

        atc_kmer = Kmers.DNAKmer{3}("ATC")
        tcg_kmer = Kmers.DNAKmer{3}("TCG")

        edge_data = graph[atc_kmer, tcg_kmer]

        # Should have evidence from both reads
        Test.@test haskey(edge_data.evidence["dataset_01"], "read_001")
        Test.@test haskey(edge_data.evidence["dataset_01"], "read_002")
    end
end

Test.@testset "Strand-Specific Graph Construction - Multiple Datasets" begin

    Test.@testset "Build graph from multiple datasets" begin
        seq1 = BioSequences.LongDNA{4}("ATCGAT")
        seq2 = BioSequences.LongDNA{4}("ATCGAT")
        qual = fill(UInt8(30), 6)

        record1 = FASTX.FASTQ.Record("read_001", seq1, qual)
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual)

        k = 3

        # Build from dataset 1
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [record1], k; dataset_id="dataset_01"
        )

        # Add dataset 2
        Mycelia.Rhizomorph.add_observations_to_graph!(
            graph, [record2], k; dataset_id="dataset_02"
        )

        # Should have evidence from both datasets
        atc_kmer = Kmers.DNAKmer{3}("ATC")
        vertex_data = graph[atc_kmer]

        Test.@test haskey(vertex_data.evidence, "dataset_01")
        Test.@test haskey(vertex_data.evidence, "dataset_02")

        # Can query dataset-specific evidence
        dataset1_evidence = Mycelia.Rhizomorph.get_dataset_evidence(
            vertex_data, "dataset_01"
        )
        dataset2_evidence = Mycelia.Rhizomorph.get_dataset_evidence(
            vertex_data, "dataset_02"
        )

        Test.@test !isnothing(dataset1_evidence)
        Test.@test !isnothing(dataset2_evidence)
    end
end

println("✓ Strand-specific graph construction tests defined (will fail until implementation)")
