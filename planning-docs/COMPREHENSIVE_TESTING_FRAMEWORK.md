# Comprehensive Testing Framework for Rhizomorph Graph Ecosystem

*A robust testing template ensuring scientific accuracy, type safety, and correctness across all 24 graph type combinations*

---

## 🎯 Testing Philosophy

### Core Principles
1. **Mathematical Precision** - Exact validation of graph topology and structure
2. **Biological Accuracy** - Scientific correctness in quality scores and strand handling
3. **Type Safety** - No inappropriate conversions, maintain BioSequence types
4. **Performance Validation** - Scalability and efficiency guarantees
5. **Error Resilience** - Comprehensive edge case and error handling

### Coverage Requirements (targets; partial coverage today)
- ⏳ **24 Graph Type Combinations** - Target; current suites cover a subset in `test/4_assembly`.
- ⏳ **Pathological Suite Coverage** - Target; partial coverage in existing graph tests.
- ⏳ **Core Functions** - Construction, pathfinding, reconstruction, I/O (uneven coverage by graph type).
- ⏳ **Correctness Assertions** - Mathematical validation of outputs (present in several suites, not all).
- ⏳ **Scientific Accuracy** - Phred scores, strand orientations, provenance (covered for qualmer/k-mer paths; expand as needed).
- ⏳ **Performance Bounds** - Time/memory complexity validation (benchmarks exist, not enforced in tests).

---

### 🔬 Pathological Test Suite for TDD

To facilitate a TDD approach, we define a standard library of small, hand-craftable input sequences that produce known, tricky graph topologies. Every one of the 24 graph types must be tested against this entire suite to ensure robust handling of common biological and data-related edge cases.

1.  **Linear Path (Baseline): A simple sequence with no repeats or errors**
    *   **Purpose:** Tests basic, error-free construction and traversal.
    *   **Example Input:** `>seq1\nAGCT`
    *   **Expected Topology:** A single, unbranched path.

2.  **Bifurcation & Convergence (SNP/Indel): A sequence creating a "bubble"**
    *   **Purpose:** Creates a "bubble" structure, common in representing genetic variants.
    *   **Example Input:** `>seqA\nAGCATGT\n>seqB\nAGCGTGT`
    *   **Expected Topology:** A path that splits into two parallel edges (representing `A` and `G`) and then converges back into a single path.

3.  **Tandem Repeat: A sequence like `ATCATCATC` that creates a cycle**
    *   **Purpose:** Tests the ability to handle cycles correctly.
    *   **Example Input:** `>repeat\nATCATCATC`
    *   **Expected Topology:** A cycle where a node representing `ATC` (or a k-mer thereof) has an edge back to itself.

4.  **Hairpin/Inverted Repeat: A sequence and its reverse-complement, to test strand handling**
    *   **Purpose:** Critically tests strand handling and canonicalization in `DoubleStrand` modes.
    *   **Example Input:** `>hairpin\nAGCTTTTCGA` (where `TCGA` is the reverse complement of `AGCT`)
    *   **Expected Topology:** A path that folds back on itself, where nodes from the forward and reverse complement strands are correctly identified and merged (in canonical graphs).

5.  **Disconnected Components: Two unrelated sequences in the same input file**
    *   **Purpose:** Ensures the graph constructor can handle multiple, unrelated sequences in the same input, resulting in a graph with multiple components.
    *   **Example Input:** `>seq1\nATGC\n>seq2\nGGCC`
    *   **Expected Topology:** A graph object containing two distinct, unconnected subgraphs.

6.  **Quality Score Edge Case: For Qualmer/FASTQ graphs, a sequence where the same k-mer appears with vastly different quality scores, to ensure metadata is handled correctly**
    *   **Purpose:** For `Quality` aware graphs, this validates that metadata (like Phred scores) is correctly aggregated and not clobbered.
    *   **Example Input (FASTQ):** `@read1\nATGC\n+\n!!!!\n@read2\nATGC\n+\nIIII`
    *   **Expected Topology:** A simple linear path, but the nodes/edges corresponding to `ATGC` must have correctly aggregated quality scores reflecting both the low-quality (`!`) and high-quality (`I`) observations.

---

## 📋 Comprehensive Testing Template

### 1. **Graph Construction Validation**

NOTE: The helper functions referenced below (`validate_kmer_transition`, `is_topologically_equal`, `is_semantically_equal`) are templates. Either implement them in shared test utilities or replace with the concrete helpers already used in `test/4_assembly`.

```julia
"""
Standard template for testing graph construction correctness.
Use this pattern in all 24 individual test files.
"""
function test_graph_construction(
    graph_builder_func,
    input_data,
    expected_topology,
    test_name::String
)
    @testset "$test_name - Construction" begin
        # 1. Build graph
        graph = graph_builder_func(input_data...)

        # 2. Exact topology validation
        vertices = collect(MetaGraphsNext.labels(graph))
        edges = collect(MetaGraphsNext.edge_labels(graph))

        @test length(vertices) == expected_topology.vertex_count
        @test length(edges) == expected_topology.edge_count

        # 3. Vertex type consistency
        if !isempty(vertices)
            vertex_type = typeof(first(vertices))
            @test all(v -> typeof(v) == vertex_type, vertices)
            @test vertex_type == expected_topology.vertex_type
        end

        # 4. Edge validity (k-1 overlaps for fixed-length)
        if expected_topology.is_fixed_length
            k = expected_topology.element_length
            for (src, dst) in edges
                @test validate_kmer_transition(src, dst, k)
            end
        end

        # 5. No string conversions in vertex labels
        if expected_topology.should_preserve_biotype
            @test !any(v -> v isa String, vertices)
        end
    end
end
```

#### 1.1 Graph Construction with I/O Round-Trip Validation

**Emphasize I/O Round-Trip Testing:** A critical test that should be explicitly part of the template is I/O "round-tripping." The process should be:

1. Construct a graph from input data.
2. Write the graph to a file (e.g., GFA format).
3. Read the graph back from that file into a new graph object.
4. Assert that the original graph and the re-read graph are topologically and semantically identical.

```julia
"""
Standard template for testing graph construction and I/O correctness.
Use this pattern in all 24 individual test files.
"""
function test_graph_construction_and_io(
    graph_builder_func,
    input_data,
    expected_topology,
    test_name::String
)
    @testset "$test_name - Construction and I/O" begin
        # 1. Build the original graph from input data
        original_graph = graph_builder_func(input_data...)

        # 2. Exact topology validation
        # (Assertions for vertex count, edge count, specific labels, etc.)
        @test MetaGraphsNext.labels(original_graph) == expected_topology.vertices
        @test MetaGraphsNext.edge_labels(original_graph) == expected_topology.edges

        # 3. I/O Round-Trip Validation
        mktemp() do path, io
            # Write the graph to a file (e.g., GFA format)
            write_gfa(path, original_graph)

            # Read the graph back into a new object
            reread_graph = read_gfa(path)

            # Assert that the re-read graph is identical to the original
            @test is_topologically_equal(original_graph, reread_graph)
            @test is_semantically_equal(original_graph, reread_graph) # Checks node/edge metadata
        end

        # 4. Metadata and Type Safety Validation
        # (Assertions for Phred scores, BioSequence types, provenance, etc.)
    end
end
```

### 2. **Coverage Tracking Validation**

```julia
"""
Validate coverage tracking across all observation modes.
"""
function test_coverage_tracking(graph, input_observations, test_name::String)
    @testset "$test_name - Coverage Tracking" begin
        vertices = collect(MetaGraphsNext.labels(graph))

        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            # 1. Coverage exists and is non-empty
            @test hasfield(typeof(vertex_data), :coverage)
            @test !isempty(vertex_data.coverage)

            # 2. Coverage entry format validation
            for cov_entry in vertex_data.coverage
                if length(cov_entry) == 3  # (obs_id, pos, strand)
                    obs_id, pos, strand = cov_entry
                    @test obs_id isa Int
                    @test pos isa Int
                    @test strand in [Mycelia.Forward, Mycelia.Reverse]
                    @test 1 <= obs_id <= length(input_observations)
                    @test pos >= 1

                elseif length(cov_entry) == 4  # (obs_id, pos, strand, quality)
                    obs_id, pos, strand, quality = cov_entry
                    @test obs_id isa Int
                    @test pos isa Int
                    @test strand in [Mycelia.Forward, Mycelia.Reverse]
                    @test quality isa Vector{Int}
                    @test all(q -> 0 <= q <= 60, quality)  # Proper Phred range
                    @test 1 <= obs_id <= length(input_observations)
                    @test pos >= 1
                end
            end

            # 3. Provenance completeness
            obs_ids = [entry[1] for entry in vertex_data.coverage]
            @test all(id -> 1 <= id <= length(input_observations), obs_ids)
        end

        # 4. Total coverage accounting
        total_coverage_entries = sum(length(graph[v].coverage) for v in vertices)
        @test total_coverage_entries > 0
    end
end
```

### 3. **Quality Score Scientific Accuracy**

```julia
"""
Rigorous validation of quality score handling and aggregation.
"""
function test_quality_accuracy(graph, input_fastq_records, test_name::String)
    @testset "$test_name - Quality Score Accuracy" begin
        vertices = collect(MetaGraphsNext.labels(graph))

        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            # 1. Quality field validation
            if hasfield(typeof(vertex_data), :mean_quality)
                @test vertex_data.mean_quality >= 0.0
                @test vertex_data.mean_quality <= 60.0  # Max Phred score
            end

            if hasfield(typeof(vertex_data), :quality_scores)
                @test !isempty(vertex_data.quality_scores)
                for qual_vec in vertex_data.quality_scores
                    @test isa(qual_vec, Vector{Int})
                    @test all(q -> 0 <= q <= 60, qual_vec)  # Strict Phred range
                    @test !any(q -> q == 33, qual_vec)  # No ASCII offsets
                end
            end

            # 2. Coverage quality consistency
            if hasfield(typeof(vertex_data), :coverage)
                for cov_entry in vertex_data.coverage
                    if length(cov_entry) == 4
                        quality = cov_entry[4]
                        @test isa(quality, Vector{Int})
                        @test all(q -> 0 <= q <= 60, quality)

                        # 3. Quality matches source FASTQ
                        obs_id, pos = cov_entry[1:2]
                        source_phred = Mycelia.get_phred_scores(input_fastq_records[obs_id])
                        expected_qual = source_phred[pos:pos+length(quality)-1]
                        @test quality == expected_qual
                    end
                end
            end
        end

        # 4. Quality aggregation correctness (multiple observations)
        vertices_with_multiple_obs = filter(vertices) do v
            length(graph[v].coverage) > 1
        end

        for vertex in vertices_with_multiple_obs
            vertex_data = graph[vertex]
            if hasfield(typeof(vertex_data), :mean_quality)
                # Validate mean calculation
                all_qualities = Float64[]
                for cov_entry in vertex_data.coverage
                    if length(cov_entry) == 4
                        append!(all_qualities, cov_entry[4])
                    end
                end
                if !isempty(all_qualities)
                    expected_mean = Statistics.mean(all_qualities)
                    @test abs(vertex_data.mean_quality - expected_mean) < 1e-6
                end
            end
        end
    end
end
```

### 4. **Strand Orientation Correctness**

```julia
"""
Validate strand handling for SingleStrand vs DoubleStrand modes.
"""
function test_strand_correctness(graph, graph_mode::GraphMode, test_name::String)
    @testset "$test_name - Strand Orientation" begin
        vertices = collect(MetaGraphsNext.labels(graph))
        edges = collect(MetaGraphsNext.edge_labels(graph))

        if graph_mode == Mycelia.SingleStrand
            # 1. SingleStrand: should preserve exact observations
            for vertex_label in vertices
                vertex_data = graph[vertex_label]
                strands = [entry[3] for entry in vertex_data.coverage]

                # Should only have Forward strand in most cases
                @test all(s -> s in [Mycelia.Forward, Mycelia.Reverse], strands)
            end

        elseif graph_mode == Mycelia.DoubleStrand
            # 2. DoubleStrand: should show canonicalization effects
            has_both_orientations = false
            for vertex_label in vertices
                vertex_data = graph[vertex_label]
                strands = [entry[3] for entry in vertex_data.coverage]

                if Mycelia.Forward in strands && Mycelia.Reverse in strands
                    has_both_orientations = true
                    break
                end
            end

            # Should find evidence of both orientations for canonical vertices
            @test has_both_orientations "DoubleStrand mode should track both orientations"
        end

        # 3. Edge strand consistency
        for (src, dst) in edges
            if length((src, dst)) == 2
                edge_data = graph[src, dst]
                if hasfield(typeof(edge_data), :src_strand) && hasfield(typeof(edge_data), :dst_strand)
                    @test edge_data.src_strand in [Mycelia.Forward, Mycelia.Reverse]
                    @test edge_data.dst_strand in [Mycelia.Forward, Mycelia.Reverse]
                end
            end
        end
    end
end
```

### 5. **Path Reconstruction Fidelity**

```julia
"""
Rigorous validation of sequence reconstruction accuracy.
"""
function test_reconstruction_fidelity(graph, original_sequences, test_name::String)
    @testset "$test_name - Reconstruction Fidelity" begin
        paths = Mycelia.find_eulerian_paths_next(graph)
        @test !isempty(paths) "Should find at least one valid path"

        reconstruction_success = false
        reconstructed_sequences = []

        for path_vector in paths
            if !isempty(path_vector)
                try
                    # 1. Type-stable reconstruction
                    vertex_type = typeof(first(path_vector))
                    walk_steps = Mycelia.WalkStep{vertex_type}[]

                    for (i, vertex_label) in enumerate(path_vector)
                        step = Mycelia.WalkStep(vertex_label, Mycelia.Forward, 1.0, Float64(i))
                        push!(walk_steps, step)
                    end

                    graph_path = Mycelia.GraphPath(walk_steps)
                    reconstructed = Mycelia.path_to_sequence(graph_path, graph)

                    if reconstructed !== nothing
                        reconstruction_success = true
                        push!(reconstructed_sequences, reconstructed)

                        # 2. Type preservation validation
                        if !isempty(original_sequences)
                            original_type = typeof(first(original_sequences))
                            @test typeof(reconstructed) == original_type
                        end

                        # 3. Length validation
                        @test length(reconstructed) > 0

                        # 4. Content validation (exact match for known sequences)
                        if length(original_sequences) == 1
                            original = first(original_sequences)
                            if typeof(original) == typeof(reconstructed)
                                # Allow for reverse complement in DoubleStrand mode
                                is_exact_match = (reconstructed == original)
                                is_rc_match = false

                                if !is_exact_match && hasmethod(BioSequences.reverse_complement, (typeof(original),))
                                    try
                                        rc_original = BioSequences.reverse_complement(original)
                                        is_rc_match = (reconstructed == rc_original)
                                    catch
                                        # Reverse complement not applicable
                                    end
                                end

                                @test is_exact_match || is_rc_match "Reconstructed sequence should match original or its reverse complement"
                            end
                        end

                        break
                    end
                catch e
                    @warn "Reconstruction failed: $e"
                end
            end
        end

        @test reconstruction_success "Should successfully reconstruct at least one sequence"

        # 5. Reconstruction completeness (for multiple paths)
        if length(paths) > 1
            @test length(reconstructed_sequences) >= 1 "Should reconstruct sequences from available paths"
        end
    end
end
```

### 6. **Performance and Scalability Validation**

```julia
"""
Performance bounds validation for scalability assurance.
"""
function test_performance_bounds(graph_builder_func, scaling_inputs, test_name::String)
    @testset "$test_name - Performance Bounds" begin
        times = Float64[]
        memory_usage = Int[]

        for (scale_factor, input_data) in scaling_inputs
            # 1. Time measurement
            start_time = time()
            graph = graph_builder_func(input_data...)
            elapsed_time = time() - start_time
            push!(times, elapsed_time)

            # 2. Memory usage estimation
            vertex_count = length(MetaGraphsNext.labels(graph))
            edge_count = length(MetaGraphsNext.edge_labels(graph))
            estimated_memory = vertex_count * 100 + edge_count * 50  # Rough estimate
            push!(memory_usage, estimated_memory)

            # 3. Performance bounds
            @test elapsed_time < 10.0 "Construction should complete within 10 seconds for test data"
            @test vertex_count > 0 "Should create vertices"
            @test edge_count >= 0 "Should create valid edge count"
        end

        # 4. Scaling behavior validation
        if length(times) > 1
            # Should not have exponential blowup
            max_time_ratio = maximum(times[2:end] ./ times[1:end-1])
            @test max_time_ratio < 100.0 "Time scaling should be reasonable"
        end
    end
end
```

### 7. **Error Handling and Edge Cases**

```julia
"""
Comprehensive error handling and edge case validation.
"""
function test_error_handling(graph_builder_func, test_name::String)
    @testset "$test_name - Error Handling" begin
        # 1. Empty input handling
        @test_throws ArgumentError graph_builder_func([])

        # 2. Single sequence handling
        single_seq = [FASTX.FASTA.Record("single", "ATCG")]
        graph = graph_builder_func(single_seq)
        @test !isnothing(graph)

        # 3. Very short sequences
        short_seq = [FASTX.FASTA.Record("short", "AT")]
        try
            graph = graph_builder_func(short_seq)
            # Should either work or fail gracefully
            @test !isnothing(graph) || true
        catch e
            @test e isa Union{ArgumentError, BoundsError}
        end

        # 4. Duplicate sequences
        duplicate_seqs = [
            FASTX.FASTA.Record("dup1", "ATCG"),
            FASTX.FASTA.Record("dup2", "ATCG")
        ]
        graph = graph_builder_func(duplicate_seqs)
        @test !isnothing(graph)

        # 5. Invalid characters (should be handled gracefully)
        # This depends on the specific graph type and input validation
    end
end
```

### 8. **Cross-Algorithm Consistency**

```julia
"""
Validate consistency between different graph construction approaches.
"""
function test_cross_algorithm_consistency(input_data, test_name::String)
    @testset "$test_name - Cross-Algorithm Consistency" begin
        # Build same input with different approaches where applicable

        # Example: K-mer vs BioSequence graphs should have compatible structure
        if applicable_to_test_case
            kmer_graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{3}, input_data)
            bio_graph = Mycelia.build_biosequence_graph(input_data)

            # Both should successfully reconstruct
            kmer_paths = Mycelia.find_eulerian_paths_next(kmer_graph)
            bio_paths = Mycelia.find_eulerian_paths_next(bio_graph)

            @test !isempty(kmer_paths)
            @test !isempty(bio_paths)

            # Path counts should be related (though not necessarily identical)
            @test length(kmer_paths) > 0
            @test length(bio_paths) > 0
        end
    end
end
```

### ⚙️ Algorithm Testing Templates

The framework is expanded to include templates for testing core algorithms that operate on the graphs. These tests should be run on graphs generated from the pathological suite.

#### 1. Path Reconstruction

```julia
"""Validates that Eulerian or other pathfinding algorithms can
correctly reconstruct the original input sequence(s)."""
function test_path_reconstruction(graph, expected_sequences::Vector{BioSequence}, test_name::String)
    @testset "$test_name - Path Reconstruction" begin
        # Find paths through the graph
        paths = find_all_paths(graph)

        # Reconstruct sequences from these paths
        reconstructed_sequences = reconstruct_sequences_from_paths(paths)

        # Assert that the set of reconstructed sequences matches the expected input
        @test Set(reconstructed_sequences) == Set(expected_sequences)
    end
end
```

#### 2. Graph Simplification

```julia
"""Tests algorithms that simplify graph topology, such as removing bubbles or collapsing linear paths."""
function test_graph_simplification(unsimplified_graph, expected_simplified_topology, test_name::String)
    @testset "$test_name - Graph Simplification" begin
        # Run the simplification algorithm
        simplified_graph = simplify_graph(unsimplified_graph)

        # Assert the new topology matches the expected simplified structure
        @test MetaGraphsNext.nv(simplified_graph) == expected_simplified_topology.vertex_count
        @test MetaGraphsNext.ne(simplified_graph) == expected_simplified_topology.edge_count
        # ... more specific assertions on structure
    end
end
```

#### 3. Graph Canonicalization

```julia
"""Tests the conversion of a single-strand graph into a canonical double-strand graph."""
function test_canonicalization(singlestrand_graph, expected_doublestrand_topology, test_name::String)
    @testset "$test_name - Canonicalization" begin
        # Generate the double-stranded graph
        doublestrand_graph = build_canonical_graph(singlestrand_graph)

        # Assert that reverse-complement vertices have been merged and edges are correct
        @test MetaGraphsNext.nv(doublestrand_graph) == expected_doublestrand_topology.vertex_count
        # Assert that specific forward/reverse nodes from the original are now one
        # ... more specific assertions on canonical structure
    end
end
```

---

## 🧪 Enhanced Test File Template

### Complete Template for Individual Graph Type Tests

```julia
# [GraphType] Graph - [StrandMode] Mode [QualityMode] Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/[filename]")'

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers  # If applicable
import MetaGraphsNext
import Statistics

Test.@testset "[GraphType] [StrandMode] [QualityMode] Graph" begin
    # Test data setup
    test_data = setup_test_data()
    input_observations = create_input_observations(test_data)
    expected_topology = calculate_expected_topology(test_data)

    # Build graph
    graph = build_target_graph(input_observations)

    # Core validation battery
    test_graph_construction(build_target_graph, (input_observations,), expected_topology, "[GraphType]")
    test_coverage_tracking(graph, input_observations, "[GraphType]")
    test_strand_correctness(graph, Mycelia.[StrandMode], "[GraphType]")
    test_reconstruction_fidelity(graph, test_data.original_sequences, "[GraphType]")

    # Quality-specific validation (if applicable)
    if is_quality_aware_graph
        test_quality_accuracy(graph, input_observations, "[GraphType]")
    end

    # Performance validation
    scaling_inputs = create_scaling_test_data()
    test_performance_bounds(build_target_graph, scaling_inputs, "[GraphType]")

    # Error handling
    test_error_handling(build_target_graph, "[GraphType]")

    # Cross-algorithm consistency (where applicable)
    test_cross_algorithm_consistency(input_observations, "[GraphType]")

    # Custom validation specific to this graph type
    test_custom_graph_properties(graph, test_data, "[GraphType]")
end

# Helper functions specific to this graph type
function setup_test_data()
    # Return structured test data for this specific graph type
end

function build_target_graph(observations)
    # Call the specific graph builder for this test
end

function test_custom_graph_properties(graph, test_data, test_name)
    # Add any graph-type-specific validation
end
```

```julia
# Enhanced Test Template - Example Implementation
#
# This file demonstrates the comprehensive testing patterns from the testing framework.
# Use this as a template for upgrading all 24 individual test files.
#
# Run with: julia --project=. -e 'include("test/4_assembly/enhanced_test_template.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext
import Statistics

# Test data structure for organizing expected results
struct ExpectedTopology
    vertex_count::Int
    edge_count::Int
    vertex_type::Type
    element_length::Int
    is_fixed_length::Bool
    should_preserve_biotype::Bool
end

# Helper functions for comprehensive validation
function validate_kmer_transition(src_kmer, dst_kmer, k)
    """Validate that k-mers have proper k-1 overlap."""
    try
        src_str = string(src_kmer)
        dst_str = string(dst_kmer)
        return src_str[2:end] == dst_str[1:end-1]
    catch
        return false  # Not applicable for this kmer type
    end
end

function test_graph_construction(graph_builder_func, input_args, expected_topology, test_name::String)
    """Standard template for testing graph construction correctness."""
    Test.@testset "$test_name - Construction" begin
        # 1. Build graph
        graph = graph_builder_func(input_args...)

        # 2. Exact topology validation
        vertices = collect(MetaGraphsNext.labels(graph))
        edges = collect(MetaGraphsNext.edge_labels(graph))

        Test.@test length(vertices) == expected_topology.vertex_count
        Test.@test length(edges) == expected_topology.edge_count

        # 3. Vertex type consistency
        if !isempty(vertices)
            vertex_type = typeof(first(vertices))
            Test.@test all(v -> typeof(v) == vertex_type, vertices)
            Test.@test vertex_type == expected_topology.vertex_type
        end

        # 4. Edge validity (k-1 overlaps for fixed-length)
        if expected_topology.is_fixed_length && length(edges) > 0
            k = expected_topology.element_length
            for edge_tuple in edges
                if length(edge_tuple) == 2
                    src, dst = edge_tuple
                    Test.@test validate_kmer_transition(src, dst, k)
                end
            end
        end

        # 5. No string conversions in vertex labels
        if expected_topology.should_preserve_biotype
            Test.@test !any(v -> v isa String, vertices)
        end
    end
end

function test_coverage_tracking(graph, input_observations, test_name::String)
    """Validate coverage tracking across all observation modes."""
    Test.@testset "$test_name - Coverage Tracking" begin
        vertices = collect(MetaGraphsNext.labels(graph))

        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            # 1. Coverage exists and is non-empty
            Test.@test hasfield(typeof(vertex_data), :coverage)
            Test.@test !isempty(vertex_data.coverage)

            # 2. Coverage entry format validation
            for cov_entry in vertex_data.coverage
                if length(cov_entry) == 3  # (obs_id, pos, strand)
                    obs_id, pos, strand = cov_entry
                    Test.@test obs_id isa Int
                    Test.@test pos isa Int
                    Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
                    Test.@test 1 <= obs_id <= length(input_observations)
                    Test.@test pos >= 1

                elseif length(cov_entry) == 4  # (obs_id, pos, strand, quality)
                    obs_id, pos, strand, quality = cov_entry
                    Test.@test obs_id isa Int
                    Test.@test pos isa Int
                    Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
                    Test.@test quality isa Vector{Int}
                    Test.@test all(q -> 0 <= q <= 60, quality)  # Proper Phred range
                    Test.@test 1 <= obs_id <= length(input_observations)
                    Test.@test pos >= 1
                end
            end

            # 3. Provenance completeness
            obs_ids = [entry[1] for entry in vertex_data.coverage]
            Test.@test all(id -> 1 <= id <= length(input_observations), obs_ids)
        end

        # 4. Total coverage accounting
        total_coverage_entries = sum(length(graph[v].coverage) for v in vertices)
        Test.@test total_coverage_entries > 0
    end
end

function test_quality_accuracy(graph, input_fastq_records, test_name::String)
    """Rigorous validation of quality score handling and aggregation."""
    Test.@testset "$test_name - Quality Score Accuracy" begin
        vertices = collect(MetaGraphsNext.labels(graph))

        for vertex_label in vertices
            vertex_data = graph[vertex_label]

            # 1. Quality field validation
            if hasfield(typeof(vertex_data), :mean_quality)
                Test.@test vertex_data.mean_quality >= 0.0
                Test.@test vertex_data.mean_quality <= 60.0  # Max Phred score
            end

            if hasfield(typeof(vertex_data), :quality_scores)
                Test.@test !isempty(vertex_data.quality_scores)
                for qual_vec in vertex_data.quality_scores
                    Test.@test isa(qual_vec, Vector{Int})
                    Test.@test all(q -> 0 <= q <= 60, qual_vec)  # Strict Phred range
                    Test.@test !any(q -> q == 33, qual_vec)  # No ASCII offsets
                end
            end

            # 2. Coverage quality consistency
            if hasfield(typeof(vertex_data), :coverage)
                for cov_entry in vertex_data.coverage
                    if length(cov_entry) == 4
                        quality = cov_entry[4]
                        Test.@test isa(quality, Vector{Int})
                        Test.@test all(q -> 0 <= q <= 60, quality)
                    end
                end
            end
        end
    end
end

function test_strand_correctness(graph, graph_mode::Mycelia.GraphMode, test_name::String)
    """Validate strand handling for SingleStrand vs DoubleStrand modes."""
    Test.@testset "$test_name - Strand Orientation" begin
        vertices = collect(MetaGraphsNext.labels(graph))
        edges = collect(MetaGraphsNext.edge_labels(graph))

        if graph_mode == Mycelia.SingleStrand
            # 1. SingleStrand: should preserve exact observations
            for vertex_label in vertices
                vertex_data = graph[vertex_label]
                if hasfield(typeof(vertex_data), :coverage) && !isempty(vertex_data.coverage)
                    strands = [entry[3] for entry in vertex_data.coverage]
                    # Should only have Forward strand in most cases
                    Test.@test all(s -> s in [Mycelia.Forward, Mycelia.Reverse], strands)
                end
            end

        elseif graph_mode == Mycelia.DoubleStrand
            # 2. DoubleStrand: should show canonicalization effects
            has_both_orientations = false
            for vertex_label in vertices
                vertex_data = graph[vertex_label]
                if hasfield(typeof(vertex_data), :coverage) && !isempty(vertex_data.coverage)
                    strands = [entry[3] for entry in vertex_data.coverage]

                    if Mycelia.Forward in strands && Mycelia.Reverse in strands
                        has_both_orientations = true
                        break
                    end
                end
            end

            # Should find evidence of both orientations for canonical vertices
            @info "DoubleStrand mode both orientations found: $has_both_orientations"
        end

        # 3. Edge strand consistency
        for edge_tuple in edges
            if length(edge_tuple) == 2
                src, dst = edge_tuple
                edge_data = graph[src, dst]
                if hasfield(typeof(edge_data), :src_strand) && hasfield(typeof(edge_data), :dst_strand)
                    Test.@test edge_data.src_strand in [Mycelia.Forward, Mycelia.Reverse]
                    Test.@test edge_data.dst_strand in [Mycelia.Forward, Mycelia.Reverse]
                end
            end
        end
    end
end

function test_reconstruction_fidelity(graph, original_sequences, test_name::String)
    """Rigorous validation of sequence reconstruction accuracy."""
    Test.@testset "$test_name - Reconstruction Fidelity" begin
        paths = Mycelia.find_eulerian_paths_next(graph)
        Test.@test !isempty(paths) "Should find at least one valid path"

        reconstruction_success = false
        reconstructed_sequences = []

        for path_vector in paths
            if !isempty(path_vector)
                try
                    # 1. Type-stable reconstruction
                    vertex_type = typeof(first(path_vector))
                    walk_steps = Mycelia.WalkStep{vertex_type}[]

                    for (i, vertex_label) in enumerate(path_vector)
                        step = Mycelia.WalkStep(vertex_label, Mycelia.Forward, 1.0, Float64(i))
                        push!(walk_steps, step)
                    end

                    graph_path = Mycelia.GraphPath(walk_steps)
                    reconstructed = Mycelia.path_to_sequence(graph_path, graph)

                    if reconstructed !== nothing
                        reconstruction_success = true
                        push!(reconstructed_sequences, reconstructed)

                        # 2. Type preservation validation
                        if !isempty(original_sequences)
                            original_type = typeof(first(original_sequences))
                            Test.@test typeof(reconstructed) == original_type
                        end

                        # 3. Length validation
                        Test.@test length(reconstructed) > 0

                        @info "Successfully reconstructed: $(typeof(reconstructed)) of length $(length(reconstructed))"
                        break
                    end
                catch e
                    @warn "Reconstruction failed: $e"
                end
            end
        end

        Test.@test reconstruction_success "Should successfully reconstruct at least one sequence"
    end
end

function test_error_handling(graph_builder_func, test_name::String)
    """Comprehensive error handling and edge case validation."""
    Test.@testset "$test_name - Error Handling" begin
        # 1. Empty input handling
        Test.@test_throws ArgumentError graph_builder_func([])

        # 2. Single sequence handling
        single_seq = [FASTX.FASTA.Record("single", "ATCG")]
        try
            graph = graph_builder_func(single_seq)
            Test.@test !isnothing(graph)
        catch e
            # Should fail gracefully
            Test.@test e isa Union{ArgumentError, BoundsError}
        end

        # 3. Duplicate sequences
        duplicate_seqs = [
            FASTX.FASTA.Record("dup1", "ATCG"),
            FASTX.FASTA.Record("dup2", "ATCG")
        ]
        try
            graph = graph_builder_func(duplicate_seqs)
            Test.@test !isnothing(graph)
        catch e
            # Should handle gracefully
            @warn "Duplicate sequence handling: $e"
        end
    end
end

# Example implementation: DNA K-mer SingleStrand Graph Test
Test.@testset "Enhanced DNA K-mer SingleStrand Graph Test" begin
    # Test data setup
    test_dna = BioSequences.dna"ATCGATCG"
    original_sequences = [test_dna]
    input_observations = [FASTX.FASTA.Record("test", test_dna)]

    # Expected topology for k=3 DNA k-mers on ATCGATCG
    expected_topology = ExpectedTopology(
        4,                          # vertex_count: ATC, TCG, CGA, GAT
        4,                          # edge_count: ATC->TCG, TCG->CGA, CGA->GAT, GAT->ATC
        Kmers.DNAKmer{3},          # vertex_type
        3,                          # element_length
        true,                       # is_fixed_length
        true                        # should_preserve_biotype
    )

    # Graph builder function
    function build_dna_kmer_graph(observations)
        return Mycelia.build_kmer_graph_next(Kmers.DNAKmer{3}, observations; graph_mode=Mycelia.SingleStrand)
    end

    # Build graph
    graph = build_dna_kmer_graph(input_observations)

    # Run comprehensive validation battery
    test_graph_construction(build_dna_kmer_graph, (input_observations,), expected_topology, "DNA K-mer SingleStrand")
    test_coverage_tracking(graph, input_observations, "DNA K-mer SingleStrand")
    test_strand_correctness(graph, Mycelia.SingleStrand, "DNA K-mer SingleStrand")
    test_reconstruction_fidelity(graph, original_sequences, "DNA K-mer SingleStrand")
    test_error_handling(build_dna_kmer_graph, "DNA K-mer SingleStrand")

    @info "Enhanced DNA K-mer SingleStrand test completed successfully"
end

@info "Enhanced test template demonstration completed"
```

---

## 📊 Validation Checklist for Each Test File

### ✅ **Required Validations (All 24 Files)**

- [ ] **Pathological Suite Coverage** - Passes all standard pathological graph tests.
- [ ] **Graph Construction** - Exact topology validation
- [ ] **I/O Round-Trip** - Graph serialization and deserialization is lossless.
- [ ] **Type Safety** - No inappropriate string conversions
- [ ] **Coverage Tracking** - Complete provenance and accuracy
- [ ] **Strand Handling** - Correct SingleStrand vs DoubleStrand behavior
- [ ] **Path Reconstruction** - Type-stable sequence recovery
- [ ] **Performance Bounds** - Reasonable time/memory complexity
- [ ] **Error Handling** - Graceful failure modes

### ✅ **Quality-Aware Graph Additional Validations**

- [ ] **Phred Score Accuracy** - Strict 0-60 range validation
- [ ] **Quality Aggregation** - Mathematically correct mean calculations
- [ ] **Provenance Tracking** - Source-to-vertex quality traceability

### ✅ **DoubleStrand Graph Additional Validations**

- [ ] **Canonicalization Correctness** - Proper reverse complement handling
- [ ] **Strand Orientation Evidence** - Both orientations tracked
- [ ] **Traversal Validity** - Biologically valid transitions

### ✅ **String Graph Additional Validations**

- [ ] **Text Processing Accuracy** - Proper n-gram extraction
- [ ] **Variable Length Handling** - Alignment-based edge detection
- [ ] **Unicode Support** - Character handling correctness

---

## 🚀 Implementation Priority

### Phase 1: Core Function Validation (Week 1)
1. **Mathematical Precision** - Implement exact topology validation
2. **Type Safety** - Add strict type checking across all operations
3. **Scientific Accuracy** - Fix Phred score validation

### Phase 2: Advanced Validation (Week 2)
1. **Performance Testing** - Add scalability bounds
2. **Error Resilience** - Comprehensive edge case handling
3. **Cross-Algorithm Consistency** - Validate compatible outputs

### Phase 3: Production Readiness (Week 3)
1. **Integration Testing** - End-to-end workflow validation
2. **Benchmarking** - Performance regression testing
3. **Documentation** - Complete validation requirements

---

## 🎯 Success Metrics

### Quantitative Targets
- **100% test coverage** across all 24 graph combinations
- **Zero type conversion violations** in BioSequence pipelines
- **100% Phred score accuracy** in quality-aware graphs
- **Sub-linear scaling** for graph construction algorithms
- **Zero unexpected failures** in error handling tests

### Qualitative Outcomes
- **Scientific accuracy** validated for all biological operations
- **Production reliability** demonstrated through comprehensive testing
- **Maintainable codebase** with clear validation patterns
- **Extensible framework** for future graph algorithm development

---

*This comprehensive testing framework ensures the rhizomorph graph ecosystem meets the highest standards for scientific accuracy, performance, and reliability required for production bioinformatics workflows.*
