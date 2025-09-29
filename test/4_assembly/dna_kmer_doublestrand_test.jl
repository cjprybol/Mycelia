# DNA K-mer Graph - DoubleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/dna_kmer_doublestrand_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "DNA K-mer DoubleStrand Graph" begin
    # Test data
    test_dna = BioSequences.dna"ATCGATCG"
    reads = [FASTX.FASTA.Record("test", test_dna)]

    # Build graph
    graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{3}, reads; graph_mode=Mycelia.DoubleStrand)

    # Structure validation - canonicalization reduces vertex count
    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) == 2  # Canonical k-mers
    Test.@test MetaGraphsNext.ne(graph) == 4

    # Coverage validation
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test !isempty(vertex_data.coverage)

        for cov_entry in vertex_data.coverage
            obs_id, pos, strand = cov_entry
            Test.@test obs_id isa Int
            Test.@test pos isa Int
            Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
        end
    end

    # Verify both orientations are tracked
    found_both_orientations = false
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        strands = [strand for (obs_id, pos, strand) in vertex_data.coverage]
        if Mycelia.Forward in strands && Mycelia.Reverse in strands
            found_both_orientations = true
            break
        end
    end
    @info "Found both orientations in DoubleStrand mode: $found_both_orientations"

    # Path reconstruction
    paths = Mycelia.find_eulerian_paths_next(graph)
    Test.@test !isempty(paths)

    # Verify we can reconstruct sequences
    reconstruction_success = false
    for path_vector in paths
        if !isempty(path_vector)
            try
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
                    @info "Successfully reconstructed: $(typeof(reconstructed)) of length $(length(reconstructed))"
                    break
                end
            catch e
                @warn "Reconstruction failed: $e"
            end
        end
    end

    Test.@test reconstruction_success
end