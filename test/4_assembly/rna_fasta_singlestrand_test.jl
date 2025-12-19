# RNA BioSequence Graph - SingleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/rna_biosequence_singlestrand_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "RNA BioSequence SingleStrand Graph" begin
    # Test data
    test_rna = BioSequences.rna"AUCGAUCG"
    reads = [FASTX.FASTA.Record("test", test_rna)]

    # Build graph
    graph = Mycelia.build_biosequence_graph(reads; graph_mode=Mycelia.SingleStrand)

    # Structure validation
    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) >= 1
    Test.@test MetaGraphsNext.ne(graph) >= 0

    # Coverage validation
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        if hasfield(typeof(vertex_data), :coverage)
            Test.@test !isempty(vertex_data.coverage)

            for cov_entry in vertex_data.coverage
                obs_id, pos, strand = cov_entry
                Test.@test obs_id isa Int
                Test.@test pos isa Int
                Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
            end
        elseif vertex_data isa Integer
            Test.@test vertex_data > 0
        end
    end

    # Path reconstruction
    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
    Test.@test !isempty(paths)

    reconstruction_success = false
    for path_vector in paths
        if !isempty(path_vector)
            try
                vertex_type = typeof(first(path_vector))
                walk_steps = Mycelia.Rhizomorph.WalkStep{vertex_type}[]

                for (i, vertex_label) in enumerate(path_vector)
                    step = Mycelia.Rhizomorph.WalkStep(vertex_label, Mycelia.Rhizomorph.Forward, 1.0, Float64(i))
                    push!(walk_steps, step)
                end

                graph_path = Mycelia.Rhizomorph.GraphPath(walk_steps)
                reconstructed = Mycelia.Rhizomorph.path_to_sequence(graph_path, graph)

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