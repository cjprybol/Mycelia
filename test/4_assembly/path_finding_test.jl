# Path Finding Tests - Eulerian Path Detection and Traversal
#
# Tests for find_eulerian_paths_next() and related path finding algorithms
#
# Run with: julia --project=. -e 'include("test/4_assembly/path_finding_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "Path Finding - Eulerian Paths" begin

    Test.@testset "Simple Linear Path - DNA K-mer" begin
        # Create simple linear sequence: ATCGATCG
        # k=3: ATC -> TCG -> CGA -> GAT -> ATC (cycle!)
        # Actually: ATCG: ATC -> TCG -> CGA
        sequence = BioSequences.dna"ATCG"
        records = [FASTX.FASTA.Record("seq1", sequence)]

        # Build k=3 DNA k-mer graph
        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{3},
            records,
            dataset_id="test_dataset",
            graph_mode=:singlestrand
        )

        # Find Eulerian paths
        paths = Mycelia.find_eulerian_paths_next(graph)

        # Should find at least one path
        Test.@test !isempty(paths)

        # Check first path
        path = first(paths)
        Test.@test !isempty(path)
        Test.@test length(path) >= 2

        # Verify path contains valid k-mers
        for kmer in path
            Test.@test kmer isa Kmers.DNAKmer{3}
        end

        # Verify path connectivity (each k-mer overlaps with next by k-1)
        for i in 1:(length(path)-1)
            kmer1 = path[i]
            kmer2 = path[i+1]
            # Check k-1 overlap
            suffix = string(kmer1)[2:end]
            prefix = string(kmer2)[1:end-1]
            Test.@test suffix == prefix
        end
    end

    Test.@testset "Multiple Sequences with Overlap" begin
        # Two sequences that overlap
        seq1 = BioSequences.dna"ATCGATCG"
        seq2 = BioSequences.dna"TCGATCGA"
        records = [
            FASTX.FASTA.Record("seq1", seq1),
            FASTX.FASTA.Record("seq2", seq2)
        ]

        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{4},
            records,
            dataset_id="multi_seq",
            graph_mode=:singlestrand
        )

        paths = Mycelia.find_eulerian_paths_next(graph)

        Test.@test !isempty(paths) "Should find paths in multi-sequence graph"

        # Check that paths are valid
        for path in paths
            if !isempty(path)
                # Verify type
                Test.@test all(k -> k isa Kmers.DNAKmer{4}, path) "All k-mers should be DNAKmer{4}"

                # Verify connectivity
                for i in 1:(length(path)-1)
                    suffix = string(path[i])[2:end]
                    prefix = string(path[i+1])[1:end-1]
                    Test.@test suffix == prefix "Path should be connected"
                end
            end
        end
    end

    Test.@testset "Circular Path (Cycle)" begin
        # Create sequence that forms a cycle: ATCATCATC
        # k=3: ATC -> TCA -> CAT -> ATC (back to start)
        sequence = BioSequences.dna"ATCATCATC"
        records = [FASTX.FASTA.Record("circular", sequence)]

        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{3},
            records,
            dataset_id="circular_test",
            graph_mode=:singlestrand
        )

        # Check graph has cycle
        vertex_count = length(MetaGraphsNext.labels(graph))
        edge_count = length(MetaGraphsNext.edge_labels(graph))

        Test.@test vertex_count > 0 "Graph should have vertices"
        Test.@test edge_count >= vertex_count "Circular graph should have edges >= vertices"

        # Find paths
        paths = Mycelia.find_eulerian_paths_next(graph)

        # May or may not find Eulerian path depending on graph structure
        # (Eulerian path exists only if at most 2 vertices have odd degree)
        # Document the behavior
        if isempty(paths)
            @info "No Eulerian path found in circular graph (expected if all vertices have even degree)"
        else
            @info "Found $(length(paths)) path(s) in circular graph"
            # Verify first path is valid
            path = first(paths)
            for i in 1:(length(path)-1)
                suffix = string(path[i])[2:end]
                prefix = string(path[i+1])[1:end-1]
                Test.@test suffix == prefix "Circular path should be connected"
            end
        end
    end

    Test.@testset "Branching Structure (Bubble)" begin
        # Create bubble: common prefix, two alternatives, common suffix
        # seq1: ATCGAT
        # seq2: ATCCAT (SNP: G->C)
        seq1 = BioSequences.dna"ATCGAT"
        seq2 = BioSequences.dna"ATCCAT"
        records = [
            FASTX.FASTA.Record("path1", seq1),
            FASTX.FASTA.Record("path2", seq2)
        ]

        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{3},
            records,
            dataset_id="bubble_test",
            graph_mode=:singlestrand
        )

        # Should have branching structure
        vertex_count = length(MetaGraphsNext.labels(graph))
        Test.@test vertex_count > 2 "Bubble should have multiple vertices"

        # Find paths
        paths = Mycelia.find_eulerian_paths_next(graph)

        # Should find multiple paths through bubble
        Test.@test length(paths) >= 1 "Should find at least one path through bubble"

        # Verify paths are valid
        for path in paths
            if !isempty(path)
                for i in 1:(length(path)-1)
                    suffix = string(path[i])[2:end]
                    prefix = string(path[i+1])[1:end-1]
                    Test.@test suffix == prefix "Bubble paths should be connected"
                end
            end
        end
    end

    Test.@testset "Disconnected Components" begin
        # Two unrelated sequences
        seq1 = BioSequences.dna"AAAA"
        seq2 = BioSequences.dna"TTTT"
        records = [
            FASTX.FASTA.Record("comp1", seq1),
            FASTX.FASTA.Record("comp2", seq2)
        ]

        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{3},
            records,
            dataset_id="disconnected",
            graph_mode=:singlestrand
        )

        # Find paths
        paths = Mycelia.find_eulerian_paths_next(graph)

        # Should find separate paths for each component
        # Or may find no paths if neither component has Eulerian path
        @info "Found $(length(paths)) path(s) in disconnected graph"

        # Verify any paths found are valid
        for path in paths
            if !isempty(path)
                for i in 1:(length(path)-1)
                    suffix = string(path[i])[2:end]
                    prefix = string(path[i+1])[1:end-1]
                    Test.@test suffix == prefix "Disconnected component paths should be valid"
                end
            end
        end
    end

    Test.@testset "RNA K-mer Paths" begin
        # Test with RNA sequence
        sequence = BioSequences.rna"AUCGAUCG"
        records = [FASTX.FASTA.Record("rna_seq", sequence)]

        graph = Mycelia.build_kmer_graph(
            Kmers.RNAKmer{3},
            records,
            dataset_id="rna_test",
            graph_mode=:singlestrand
        )

        paths = Mycelia.find_eulerian_paths_next(graph)

        Test.@test !isempty(paths) "Should find RNA k-mer paths"

        path = first(paths)
        if !isempty(path)
            Test.@test all(k -> k isa Kmers.RNAKmer{3}, path) "RNA paths should contain RNAKmer"

            # Verify connectivity
            for i in 1:(length(path)-1)
                suffix = string(path[i])[2:end]
                prefix = string(path[i+1])[1:end-1]
                Test.@test suffix == prefix "RNA path should be connected"
            end
        end
    end

    Test.@testset "AA K-mer Paths" begin
        # Test with amino acid sequence
        sequence = BioSequences.aa"MKKLAVAA"
        records = [FASTX.FASTA.Record("protein", sequence)]

        graph = Mycelia.build_kmer_graph(
            Kmers.AAKmer{3},
            records,
            dataset_id="aa_test",
            graph_mode=:singlestrand
        )

        paths = Mycelia.find_eulerian_paths_next(graph)

        Test.@test !isempty(paths) "Should find AA k-mer paths"

        path = first(paths)
        if !isempty(path)
            Test.@test all(k -> k isa Kmers.AAKmer{3}, path) "AA paths should contain AAKmer"

            # Verify connectivity
            for i in 1:(length(path)-1)
                suffix = string(path[i])[2:end]
                prefix = string(path[i+1])[1:end-1]
                Test.@test suffix == prefix "AA path should be connected"
            end
        end
    end

    Test.@testset "Large K-mer Size (k=31)" begin
        # Test with larger k
        sequence = BioSequences.dna"ATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        records = [FASTX.FASTA.Record("long_seq", sequence)]

        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{31},
            records,
            dataset_id="large_k",
            graph_mode=:singlestrand
        )

        paths = Mycelia.find_eulerian_paths_next(graph)

        # Should find path even with large k
        Test.@test !isempty(paths) "Should find paths with k=31"

        path = first(paths)
        if !isempty(path)
            Test.@test all(k -> k isa Kmers.DNAKmer{31}, path) "Paths should contain DNAKmer{31}"

            # Verify connectivity with k=31
            for i in 1:(length(path)-1)
                suffix = string(path[i])[2:end]
                prefix = string(path[i+1])[1:end-1]
                Test.@test suffix == prefix "Large k-mer path should be connected"
            end
        end
    end

    Test.@testset "DoubleStrand Mode" begin
        # Test with doublestrand graph
        sequence = BioSequences.dna"ATCGAT"
        records = [FASTX.FASTA.Record("ds_test", sequence)]

        graph = Mycelia.build_kmer_graph(
            Kmers.DNAKmer{3},
            records,
            dataset_id="doublestrand_test",
            graph_mode=:doublestrand
        )

        # DoubleStrand graph includes both forward and reverse complement
        paths = Mycelia.find_eulerian_paths_next(graph)

        # Should find paths (behavior may differ from singlestrand)
        @info "DoubleStrand mode found $(length(paths)) path(s)"

        # Verify any paths found are valid
        for path in paths
            if !isempty(path)
                Test.@test all(k -> k isa Kmers.DNAKmer{3}, path) "DoubleStrand paths should contain valid k-mers"

                # Verify connectivity
                for i in 1:(length(path)-1)
                    suffix = string(path[i])[2:end]
                    prefix = string(path[i+1])[1:end-1]
                    Test.@test suffix == prefix "DoubleStrand path should be connected"
                end
            end
        end
    end

    Test.@testset "Empty Graph Error Handling" begin
        # Test with empty input
        Test.@test_throws ArgumentError begin
            Mycelia.build_kmer_graph(
                Kmers.DNAKmer{3},
                FASTX.FASTA.Record[],
                dataset_id="empty",
                graph_mode=:singlestrand
            )
        end
    end

    @info "✅ Path finding tests completed"
end
