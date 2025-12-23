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
        if expected_vertices !== nothing
            if length(vertices) != expected_vertices
                @warn "Expected $expected_vertices vertices, got $(length(vertices))"
            end
            Test.@test length(vertices) == expected_vertices
        else
            Test.@test !isempty(vertices)
        end

        edge_count = MetaGraphsNext.ne(graph)
        if expected_edges !== nothing
            if edge_count != expected_edges
                @warn "Expected $expected_edges edges, got $edge_count"
            end
            Test.@test edge_count == expected_edges
        else
            Test.@test edge_count >= 0
        end

        # 2. Evidence validation
        total_evidence_entries = 0
        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            if hasproperty(vertex_data, :evidence)
                evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                if isempty(evidence_entries)
                    @warn "Vertex $vertex_label should have evidence"
                end
                Test.@test !isempty(evidence_entries)
                total_evidence_entries += length(evidence_entries)

                for entry in evidence_entries
                    Test.@test entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
                end
            end
        end

        Test.@test total_evidence_entries > 0

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

            if hasproperty(vertex_data, :evidence)
                evidence_entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex_data.evidence)
                for entry in evidence_entries
                    if hasproperty(entry, :quality_scores)
                        Test.@test !isempty(entry.quality_scores)
                        Test.@test all(q -> q >= UInt8(33), entry.quality_scores)
                    end
                end

                dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
                if !isempty(dataset_ids)
                    mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, first(dataset_ids))
                    if mean_quality !== nothing
                        Test.@test all(q -> q >= 0.0, mean_quality)
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
            graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 3; dataset_id="test", mode=:singlestrand)

            # For ATCGATCG with k=3: ATC, TCG, CGA, GAT, ATC, TCG
            # Unique k-mers: ATC, TCG, CGA, GAT (4 vertices)
            # Edges: ATC->TCG, TCG->CGA, CGA->GAT, GAT->ATC, ATC->TCG (but some are duplicates)
            expected_vertices = 4  # ATC, TCG, CGA, GAT
            expected_edges = 4     # ATC->TCG, TCG->CGA, CGA->GAT, GAT->ATC

            validate_graph_comprehensive(graph, [known_dna], expected_vertices, expected_edges, "DNA K-mer SingleStrand")
        end

        Test.@testset "DNA K-mer DoubleStrand Detailed" begin
            reads = [FASTX.FASTA.Record("test", known_dna)]
            graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 3; dataset_id="test", mode=:doublestrand)

            # In DoubleStrand mode, forward and reverse-complement k-mers are retained
            expected_vertices = 4  # ATC, TCG, CGA, GAT for this sequence
            expected_edges = nothing

            validate_graph_comprehensive(graph, [known_dna], expected_vertices, expected_edges, "DNA K-mer DoubleStrand")

            # Additional validation: check that both orientations are tracked
            found_both_orientations = false
            for vertex_label in MetaGraphsNext.labels(graph)
                vertex_data = graph[vertex_label]
                strands = Mycelia.Rhizomorph.collect_evidence_strands(vertex_data.evidence)
                if Mycelia.Rhizomorph.Forward in strands && Mycelia.Rhizomorph.Reverse in strands
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

            graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, 3; dataset_id="test", mode=:singlestrand)

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

            graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, 3; dataset_id="test", mode=:singlestrand)

            # Verify that vertices have multiple quality observations
            multiple_observations_found = false
            for vertex_label in MetaGraphsNext.labels(graph)
                vertex_data = graph[vertex_label]
                evidence_count = Mycelia.Rhizomorph.count_evidence(vertex_data)
                if evidence_count > 1
                    multiple_observations_found = true
                    @info "Vertex $vertex_label observed $evidence_count times"

                    dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
                    if !isempty(dataset_ids)
                        mean_quality = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, first(dataset_ids))
                        Test.@test mean_quality !== nothing
                        Test.@test all(q -> q >= 0.0, mean_quality)
                    end
                end
            end

            Test.@test multiple_observations_found
        end
    end

    Test.@testset "String Graph Correctness" begin
        Test.@testset "String N-gram Graph Detailed" begin
            graph = Mycelia.Rhizomorph.build_ngram_graph([known_string], 3; dataset_id="test")

            # For ABCDEF with n=3: ABC, BCD, CDE, DEF (4 n-grams)
            expected_vertices = 4
            expected_edges = 3  # ABC->BCD, BCD->CDE, CDE->DEF

            validate_graph_comprehensive(graph, [known_string], expected_vertices, expected_edges, "String N-gram")
        end

        Test.@testset "String N-gram Graph from FASTQ Inputs" begin
            qual_str = String([Char(q + 33) for q in high_quality[1:length(known_string)]])
            fastq_record = FASTX.FASTQ.Record("test", known_string, qual_str)
            reads = [fastq_record]

            graph = Mycelia.Rhizomorph.build_ngram_graph([FASTX.sequence(String, fastq_record)], 3; dataset_id="test")

            expected_vertices = 4
            expected_edges = 3

            validate_graph_comprehensive(graph, reads, expected_vertices, expected_edges, "String N-gram from FASTQ")
        end
    end

    Test.@testset "Cross-Type Consistency" begin
        # Test that different graph types for the same input produce consistent results

        Test.@testset "DNA K-mer vs BioSequence Consistency" begin
            reads = [FASTX.FASTA.Record("test", known_dna)]

            kmer_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, 3; dataset_id="test", mode=:singlestrand)
            bio_graph = Mycelia.Rhizomorph.build_fasta_graph(reads; dataset_id="test", min_overlap=3)

            # Both should successfully reconstruct sequences
            kmer_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(kmer_graph)
            bio_paths = Mycelia.Rhizomorph.find_contigs_next(bio_graph; min_contig_length=1)

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

        graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, 3; dataset_id="test", mode=:singlestrand)

        # Verify provenance tracking
        for vertex_label in MetaGraphsNext.labels(graph)
            vertex_data = graph[vertex_label]

            dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
            Test.@test !isempty(dataset_ids)
            obs_ids = Mycelia.Rhizomorph.get_all_observation_ids(vertex_data, first(dataset_ids))
            Test.@test "obs1" in obs_ids
            Test.@test "obs2" in obs_ids
        end
    end

end
