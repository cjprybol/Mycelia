# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/aa_kmer_graph_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/aa_kmer_graph_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Amino Acid K-mer Graph Construction Test
# Tests for amino acid k-mer graph construction in singlestrand mode.
# Validates that protein sequences are properly processed with FwAAMers
# (or UnambiguousAAMers if available).

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "AA K-mer Graph - Singlestrand" begin

    Test.@testset "AA K-mer Extraction" begin
        # Protein sequence
        aa_seq = BioSequences.LongAA("MKVLWAALL")
        k = 3

        # Extract k-mers using FwAAMers
        # Note: FwAAMers returns just k-mers, not (kmer, position) tuples
        kmers = collect(Kmers.FwAAMers{k}(aa_seq))

        Test.@test length(kmers) == 7  # 9 - 3 + 1 = 7
        # Check that k-mers are AAKmers
        Test.@test all(kmer -> kmer isa Kmers.AAKmer{3}, kmers)
    end

    Test.@testset "AA Graph Construction - Single Read" begin
        records = [
            FASTX.FASTA.Record("protein_001", "MKVLW")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 3 unique AA k-mers
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 3
        Test.@test Mycelia.Rhizomorph.edge_count(graph) == 2

        # Verify AA k-mers are present
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("MKV"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("KVL"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("VLW"))
    end

    Test.@testset "AA Graph Construction - Multiple Reads" begin
        records = [
            FASTX.FASTA.Record("protein_001", "MKVLW"),
            FASTX.FASTA.Record("protein_002", "MKVLW")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Same k-mers, but observed twice
        kmer = Kmers.AAKmer{3}("MKV")
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, kmer)
        Test.@test count == 2
    end

    Test.@testset "AA Graph - Evidence Tracking" begin
        records = [
            FASTX.FASTA.Record("protein_001", "MKVLW")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="aa_test")

        kmer = Kmers.AAKmer{3}("MKV")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        Test.@test !isnothing(vertex_data)
        Test.@test haskey(vertex_data.evidence, "aa_test")
        Test.@test haskey(vertex_data.evidence["aa_test"], "protein_001")

        # Check evidence entry
        evidence_set = vertex_data.evidence["aa_test"]["protein_001"]
        Test.@test length(evidence_set) == 1
        evidence = first(evidence_set)
        Test.@test evidence.position == 1
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "AA Graph - No Reverse Complement" begin
        # Protein sequences don't have reverse complement
        # Each sequence should be stored as-is
        records = [
            FASTX.FASTA.Record("protein_001", "MKV"),
            FASTX.FASTA.Record("protein_002", "VKM")  # Different sequence, not RC
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 2 separate vertices
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 2
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("MKV"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("VKM"))
    end

    Test.@testset "AA Graph - Ambiguous Amino Acids" begin
        # X represents any amino acid (ambiguous)
        # Test depends on whether FwAAMers skips X or includes it
        records = [
            FASTX.FASTA.Record("protein_001", "MKVXXLW")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have MKV at minimum
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("MKV"))

        # Check if any vertices contain X (iterator-dependent behavior)
        vertices = Mycelia.Rhizomorph.get_all_vertices(graph)
        # We're just testing that graph construction doesn't crash with X
        Test.@test length(vertices) >= 1
    end

    Test.@testset "AA Graph - Edge Evidence" begin
        records = [
            FASTX.FASTA.Record("protein_001", "MKVL")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.AAKmer{3}("MKV")
        dst = Kmers.AAKmer{3}("KVL")

        Test.@test Mycelia.Rhizomorph.has_edge(graph, src, dst)

        edge_data = Mycelia.Rhizomorph.get_edge_data(graph, src, dst)
        Test.@test !isnothing(edge_data)
        Test.@test haskey(edge_data.evidence, "test")

        # Check edge evidence
        edge_evidence = edge_data.evidence["test"]["protein_001"]
        Test.@test length(edge_evidence) == 1
        evidence = first(edge_evidence)
        Test.@test evidence.from_position == 1
        Test.@test evidence.to_position == 2
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "AA Graph - Path Assembly" begin
        records = [
            FASTX.FASTA.Record("protein_001", "MKVLW")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Get path
        path = [
            Kmers.AAKmer{3}("MKV"),
            Kmers.AAKmer{3}("KVL"),
            Kmers.AAKmer{3}("VLW")
        ]

        # Assemble sequence from path
        sequence = Mycelia.Rhizomorph.assemble_path_sequence(path)

        Test.@test sequence == BioSequences.LongAA("MKVLW")
        Test.@test sequence isa BioSequences.LongAA
    end

    Test.@testset "AA Graph - Query Functions" begin
        records = [
            FASTX.FASTA.Record("protein_001", "MKVLWALL")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Test source/sink detection
        sources = Mycelia.Rhizomorph.get_all_sources(graph)
        sinks = Mycelia.Rhizomorph.get_all_sinks(graph)

        Test.@test Kmers.AAKmer{3}("MKV") in sources
        Test.@test Kmers.AAKmer{3}("ALL") in sinks

        # Test neighbor queries
        kmer = Kmers.AAKmer{3}("KVL")
        incoming = Mycelia.Rhizomorph.get_incoming_neighbors(graph, kmer)
        outgoing = Mycelia.Rhizomorph.get_outgoing_neighbors(graph, kmer)

        Test.@test Kmers.AAKmer{3}("MKV") in incoming
        Test.@test Kmers.AAKmer{3}("VLW") in outgoing
    end

    Test.@testset "AA Graph - Common Protein Motifs" begin
        # Test with realistic protein sequence containing common motifs
        records = [
            FASTX.FASTA.Record("protein_001", "GPGPGPG")  # Collagen-like repeat
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should handle repeated k-mers
        kmer_gpg = Kmers.AAKmer{3}("GPG")
        kmer_pgp = Kmers.AAKmer{3}("PGP")

        Test.@test Mycelia.Rhizomorph.has_vertex(graph, kmer_gpg)
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, kmer_pgp)

        # Check observation count for repeated k-mer
        # GPG appears at positions 1, 3, 5
        count_gpg = Mycelia.Rhizomorph.get_vertex_observation_count(graph, kmer_gpg)
        Test.@test count_gpg == 1  # 1 observation (same read), multiple positions in evidence
    end

    Test.@testset "AA Graph - Reduced Alphabet Input (MURPHY2)" begin
        full_seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
        reduced_seq = Mycelia.reduce_amino_acid_alphabet(full_seq, :MURPHY2)
        records = [FASTX.FASTA.Record("protein_reduced", reduced_seq)]

        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="reduced")

        # Graph should build and contain only reduced alphabet symbols (I/E) in labels
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) > 0
        vertex_labels = collect(Mycelia.Rhizomorph.get_all_vertices(graph))
        Test.@test all(label -> all(c -> c in ('I','E'), String(label)), vertex_labels)
    end
end

println("✓ AA k-mer graph tests completed")
