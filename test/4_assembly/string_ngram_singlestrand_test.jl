# String N-gram Graph - SingleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/string_ngram_singlestrand_test.jl")'

import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "String N-gram SingleStrand Graph" begin
    test_string = "ABCDEF"
    reads = [FASTX.FASTA.Record("test", test_string)]

    graph = Mycelia.build_string_ngram_graph_next(reads, 3; graph_mode=Mycelia.SingleStrand)

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) == 4  # ABC, BCD, CDE, DEF
    Test.@test MetaGraphsNext.ne(graph) == 3  # ABC->BCD, BCD->CDE, CDE->DEF

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

    paths = Mycelia.find_eulerian_paths_next(graph)
    Test.@test !isempty(paths)

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

    Test.@testset "Reduced Amino Acid Alphabet as Unicode Input" begin
        # Use reduced alphabet that includes non-standard symbols to ensure Unicode graph support
        full_seq = "ACDEFGHIKLMNPQRSTVWY"
        reduced = Mycelia.reduce_amino_acid_alphabet(full_seq, :CHEMICAL6)  # includes '-' and '+'
        records = [FASTX.FASTA.Record("reduced_seq", reduced)]

        graph_reduced = Mycelia.build_string_ngram_graph_next(records, 3; graph_mode=Mycelia.SingleStrand)
        vertices_reduced = collect(MetaGraphsNext.labels(graph_reduced))

        Test.@test !isempty(vertices_reduced)
        # Ensure characters from reduced alphabet appear in labels
        Test.@test any(label -> occursin("+", String(label)) || occursin("-", String(label)), vertices_reduced)
    end
end
