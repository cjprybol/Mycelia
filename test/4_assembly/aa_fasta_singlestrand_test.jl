# Amino Acid BioSequence Graph - SingleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/aa_fasta_singlestrand_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "AA BioSequence SingleStrand Graph" begin
    test_aa = BioSequences.aa"MKTV"
    reads = [FASTX.FASTA.Record("test", test_aa)]

    graph = Mycelia.Rhizomorph.build_fasta_graph(reads; dataset_id="aa_fasta_test", min_overlap=3)

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) == 1
    Test.@test MetaGraphsNext.ne(graph) == 0

    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test hasfield(typeof(vertex_data), :evidence)
        dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
        Test.@test "aa_fasta_test" in dataset_ids
        obs_ids = Mycelia.Rhizomorph.get_all_observation_ids(vertex_data, "aa_fasta_test")
        Test.@test "test" in obs_ids
        obs_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex_data, "aa_fasta_test", "test")
        Test.@test !isnothing(obs_evidence)
        if !isnothing(obs_evidence)
            for entry in obs_evidence
                Test.@test entry.position == 1
                Test.@test entry.strand == Mycelia.Rhizomorph.Forward
            end
        end
    end

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
