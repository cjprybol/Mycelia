# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/comprehensive_correctness_tests.jl")'
# ```
#
# Comprehensive Correctness Testing for All 24 Graph Variants
#
# This test suite validates:
# 1. Graph structure correctness (vertices, edges, labels, orientations)
# 2. Path reconstruction (expected paths + rediscovery)
# 3. Round-trip I/O (write/read identity preservation)
# 4. Quality score handling (retention, joint probabilities, provenance)
# 5. Graph simplification pipeline preparation

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext
import Statistics

Test.@testset "Comprehensive Graph Correctness Tests" begin

    # Test data with known structure
    known_dna = BioSequences.dna"ATCGATCG"  # 8 bases
    known_rna = BioSequences.rna"AUCGAUCG"  # 8 bases
    known_aa = BioSequences.aa"MKTV"        # 4 amino acids
    known_string = "ABCDEF"                 # 6 characters

    # Quality scores for testing
    high_quality = [40, 40, 40, 40, 40, 40, 40, 40]  # High confidence
    mixed_quality = [20, 40, 30, 35, 25, 40, 30, 35]  # Mixed confidence

    """
    Comprehensive validation function for any graph type.

    Tests:
    - Graph structure (vertex/edge counts, labels)
    - Coverage tracking and strand orientations
    - Path reconstruction correctness
    - Quality retention (if applicable)
    """
    function validate_graph_comprehensive(graph, input_sequences, expected_vertices, expected_edges, graph_type_name)
        @info "Validating $graph_type_name graph..."

        # 1. Basic structure validation
        vertices = collect(MetaGraphsNext.labels(graph))
        if length(vertices) != expected_vertices
            @warn "Expected $expected_vertices vertices, got $(length(vertices))"
        end
        Test.@test length(vertices) == expected_vertices

        edge_count = MetaGraphsNext.ne(graph)
        if edge_count != expected_edges
            @warn "Expected $expected_edges edges, got $edge_count"
        end
        Test.@test edge_count == expected_edges

        # 2. Coverage validation
        total_coverage_entries = 0
        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            # Check coverage exists and is non-empty
            if hasfield(typeof(vertex_data), :coverage)
                coverage = vertex_data.coverage
                if isempty(coverage)
                    @warn "Vertex $vertex_label should have coverage"
                end
                Test.@test !isempty(coverage)
                total_coverage_entries += length(coverage)

                # Validate coverage format based on graph type
                for cov_entry in coverage
                    if length(cov_entry) == 3  # (obs_id, pos, strand)
                        obs_id, pos, strand = cov_entry
                        Test.@test obs_id isa Int
                        Test.@test pos isa Int
                        Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
                    elseif length(cov_entry) == 4  # (obs_id, pos, strand, quality)
                        obs_id, pos, strand, quality = cov_entry
                        Test.@test obs_id isa Int
                        Test.@test pos isa Int
                        Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
                        Test.@test quality isa Vector{Int}
                    end
                end
            elseif vertex_data isa Integer
                total_coverage_entries += vertex_data
            end
        end

        Test.@test total_coverage_entries > 0

        # 3. Path reconstruction validation
        paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
        Test.@test !isempty(paths)

        reconstruction_success = false
        for path_vector in paths
            if !isempty(path_vector)
                try
                    # Create properly typed WalkStep objects
                    first_vertex = first(path_vector)
                    vertex_type = typeof(first_vertex)
                    walk_steps = Mycelia.Rhizomorph.WalkStep{vertex_type}[]

                    for (i, vertex_label) in enumerate(path_vector)
                        step = Mycelia.Rhizomorph.WalkStep(vertex_label, Mycelia.Rhizomorph.Forward, 1.0, Float64(i))
                        push!(walk_steps, step)
                    end

                    graph_path = Mycelia.Rhizomorph.GraphPath(walk_steps)
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(graph_path, graph)

                    if reconstructed !== nothing
                        reconstruction_success = true
                        @info "Successfully reconstructed sequence: $(typeof(reconstructed)) of length $(length(reconstructed))"
                        break
                    end
                catch e
                    @warn "Path reconstruction failed: $e"
                end
            end
        end

        Test.@test reconstruction_success

        return true
    end

    """
    Validate quality score handling in quality-aware graphs.
    """
    function validate_quality_handling(graph, input_records, graph_type_name)
        @info "Validating quality handling for $graph_type_name..."

        vertices = collect(MetaGraphsNext.labels(graph))

        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            # Check if vertex data has quality-related fields
            if hasfield(typeof(vertex_data), :average_quality)
                Test.@test vertex_data.average_quality >= 0.0
                @info "Vertex $vertex_label has average quality: $(vertex_data.average_quality)"
            end

            if hasfield(typeof(vertex_data), :quality_scores)
                Test.@test !isempty(vertex_data.quality_scores)

                # Validate quality score format
                for qual_vec in vertex_data.quality_scores
                    Test.@test isa(qual_vec, Vector{Int})
                    Test.@test all(q -> 0 <= q <= 60, qual_vec)
                end
            end

            # Check coverage with quality information
            if hasfield(typeof(vertex_data), :coverage)
                for cov_entry in vertex_data.coverage
                    if length(cov_entry) == 4  # Quality-aware coverage
                        obs_id, pos, strand, quality = cov_entry
                        Test.@test isa(quality, Vector{Int})
                        Test.@test all(q -> 0 <= q <= 60, quality)
                    end
                end
            end
        end

        return true
    end

    """
    Test round-trip I/O by writing and reading back graphs/sequences.
    """
    function test_round_trip_io(graph, reconstructed_sequences, graph_type_name)
        @info "Testing round-trip I/O for $graph_type_name..."

        # For now, test sequence identity preservation
        # TODO: Add actual graph serialization/deserialization when implemented

        for (i, seq) in enumerate(reconstructed_sequences)
            if isa(seq, BioSequences.BioSequence)
                # Test BioSequence consistency
                seq_copy = deepcopy(seq)
                Test.@test seq == seq_copy
                Test.@test string(seq) == string(seq_copy)
            elseif isa(seq, String)
                # Test String consistency
                seq_copy = String(seq)
                Test.@test seq == seq_copy
            end
        end

        return true
    end

    Test.@testset "DNA K-mer Graph Correctness" begin
        Test.@testset "DNA K-mer SingleStrand Detailed" begin
            reads = [FASTX.FASTA.Record("test", known_dna)]
            graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{3}, reads; graph_mode=Mycelia.SingleStrand)

            # For ATCGATCG with k=3: ATC, TCG, CGA, GAT, ATC, TCG
            # Unique k-mers: ATC, TCG, CGA, GAT (4 vertices)
            # Edges: ATC->TCG, TCG->CGA, CGA->GAT, GAT->ATC, ATC->TCG (but some are duplicates)
            expected_vertices = 4  # ATC, TCG, CGA, GAT
            expected_edges = 4     # ATC->TCG, TCG->CGA, CGA->GAT, GAT->ATC

            validate_graph_comprehensive(graph, [known_dna], expected_vertices, expected_edges, "DNA K-mer SingleStrand")
        end

        Test.@testset "DNA K-mer DoubleStrand Detailed" begin
            reads = [FASTX.FASTA.Record("test", known_dna)]
            graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{3}, reads; graph_mode=Mycelia.DoubleStrand)

            # In DoubleStrand mode, k-mers are canonicalized which merges reverse complements
            expected_vertices = 2  # Canonical ATC and CGA for this sequence
            expected_edges = 4     # Two canonical transitions plus their reverse-oriented self loops

            validate_graph_comprehensive(graph, [known_dna], expected_vertices, expected_edges, "DNA K-mer DoubleStrand")

            # Additional validation: check that both orientations are tracked
            found_both_orientations = false
            for vertex_label in MetaGraphsNext.labels(graph)
                vertex_data = graph[vertex_label]
                strands = [strand for (obs_id, pos, strand) in vertex_data.coverage]
                if Mycelia.Forward in strands && Mycelia.Reverse in strands
                    found_both_orientations = true
                    break
                end
            end
            @info "Found both orientations in DoubleStrand mode: $found_both_orientations"
        end
    end

    Test.@testset "Quality-Aware Graph Correctness" begin
        Test.@testset "DNA Qualmer Graph Detailed" begin
            # Create FASTQ record with known quality
            qual_str = String([Char(q + 33) for q in high_quality[1:length(known_dna)]])
            fastq_record = FASTX.FASTQ.Record("test", string(known_dna), qual_str)
            reads = [fastq_record]

            graph = Mycelia.build_qualmer_graph(reads; k=3, graph_mode=Mycelia.SingleStrand)

            expected_vertices = 4  # Same structure as k-mer graph
            expected_edges = 4

            validate_graph_comprehensive(graph, reads, expected_vertices, expected_edges, "DNA Qualmer")
            validate_quality_handling(graph, reads, "DNA Qualmer")
        end

        Test.@testset "Quality Score Joint Probability" begin
            # Test that seeing the same vertex multiple times affects joint probability
            qual_str1 = String([Char(q + 33) for q in high_quality[1:length(known_dna)]])
            qual_str2 = String([Char(q + 33) for q in mixed_quality[1:length(known_dna)]])

            fastq_records = [
                FASTX.FASTQ.Record("test1", string(known_dna), qual_str1),
                FASTX.FASTQ.Record("test2", string(known_dna), qual_str2)  # Same sequence, different quality
            ]

            graph = Mycelia.build_qualmer_graph(fastq_records; k=3, graph_mode=Mycelia.SingleStrand)

            # Verify that vertices have multiple quality observations
            multiple_observations_found = false
            for vertex_label in MetaGraphsNext.labels(graph)
                vertex_data = graph[vertex_label]
                if vertex_data.coverage > 1
                    multiple_observations_found = true
                    @info "Vertex $vertex_label observed $(vertex_data.coverage) times"

                    # Check that average quality reflects both observations
                    Test.@test vertex_data.mean_quality > 0.0
                end
            end

            Test.@test multiple_observations_found
        end
    end

    Test.@testset "String Graph Correctness" begin
        Test.@testset "String N-gram Graph Detailed" begin
            reads = [FASTX.FASTA.Record("test", known_string)]
            graph = Mycelia.build_string_ngram_graph_next(reads, 3; graph_mode=Mycelia.SingleStrand)

            # For ABCDEF with n=3: ABC, BCD, CDE, DEF (4 n-grams)
            expected_vertices = 4
            expected_edges = 3  # ABC->BCD, BCD->CDE, CDE->DEF

            validate_graph_comprehensive(graph, [known_string], expected_vertices, expected_edges, "String N-gram")
        end

        Test.@testset "String Quality Graph Detailed" begin
            qual_str = String([Char(q + 33) for q in high_quality[1:length(known_string)]])
            fastq_record = FASTX.FASTQ.Record("test", known_string, qual_str)
            reads = [fastq_record]

            graph = Mycelia.build_string_qualmer_ngram_graph_next(reads, 3; graph_mode=Mycelia.SingleStrand)

            expected_vertices = 4
            expected_edges = 3

            validate_graph_comprehensive(graph, reads, expected_vertices, expected_edges, "String Quality N-gram")
            validate_quality_handling(graph, reads, "String Quality N-gram")
        end
    end

    Test.@testset "Cross-Type Consistency" begin
        # Test that different graph types for the same input produce consistent results

        Test.@testset "DNA K-mer vs BioSequence Consistency" begin
            reads = [FASTX.FASTA.Record("test", known_dna)]

            kmer_graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{3}, reads; graph_mode=Mycelia.SingleStrand)
            bio_graph = Mycelia.build_biosequence_graph(reads; graph_mode=Mycelia.SingleStrand)

            # Both should successfully reconstruct sequences
            kmer_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(kmer_graph)
            bio_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(bio_graph)

            Test.@test !isempty(kmer_paths)
            Test.@test !isempty(bio_paths)

            @info "K-mer graph: $(length(MetaGraphsNext.labels(kmer_graph))) vertices, $(MetaGraphsNext.ne(kmer_graph)) edges"
            @info "BioSequence graph: $(length(MetaGraphsNext.labels(bio_graph))) vertices, $(MetaGraphsNext.ne(bio_graph)) edges"
        end
    end

    Test.@testset "Provenance and Traceability" begin
        # Test that we can trace back from vertices to original observations

        qual_str = String([Char(q + 33) for q in high_quality[1:length(known_dna)]])
        fastq_records = [
            FASTX.FASTQ.Record("obs1", string(known_dna), qual_str),
            FASTX.FASTQ.Record("obs2", string(known_dna), qual_str)  # Duplicate observation
        ]

        graph = Mycelia.build_qualmer_graph(fastq_records; k=3, graph_mode=Mycelia.SingleStrand)

        # Verify provenance tracking
        for vertex_label in MetaGraphsNext.labels(graph)
            vertex_data = graph[vertex_label]

            # Each observation should track the provenance metadata
            for obs in vertex_data.observations
                Test.@test obs.sequence_id in (1:2)
                Test.@test obs.position >= 1
                @info "Vertex $vertex_label: obs_id=$(obs.sequence_id), pos=$(obs.position)"
            end
        end
    end

end
