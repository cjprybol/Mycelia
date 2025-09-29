# String N-gram Graph - SingleStrand Mode with Quality Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/string_ngram_singlestrand_quality_test.jl")'

import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "String N-gram SingleStrand Quality Graph" begin
    test_string = "ABCDEF"
    high_quality = [40, 40, 40, 40, 40, 40]
    qual_str = String([Char(q + 33) for q in high_quality])
    fastq_record = FASTX.FASTQ.Record("test", test_string, qual_str)
    reads = [fastq_record]

    graph = Mycelia.build_string_qualmer_ngram_graph_next(reads, 3; graph_mode=Mycelia.SingleStrand)

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) == 4  # ABC, BCD, CDE, DEF
    Test.@test MetaGraphsNext.ne(graph) == 3  # ABC->BCD, BCD->CDE, CDE->DEF

    # Quality-aware coverage validation
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        if hasfield(typeof(vertex_data), :coverage)
            Test.@test !isempty(vertex_data.coverage)
            for cov_entry in vertex_data.coverage
                if length(cov_entry) == 4  # Quality-aware coverage
                    obs_id, pos, strand, quality = cov_entry
                    Test.@test obs_id isa Int
                    Test.@test pos isa Int
                    Test.@test strand in [Mycelia.Forward, Mycelia.Reverse]
                    Test.@test quality isa Vector{Int}
                    Test.@test all(q -> 0 <= q <= 60, quality)
                end
            end
        end

        # Check quality-related fields
        if hasfield(typeof(vertex_data), :average_quality)
            Test.@test vertex_data.average_quality >= 0.0
        end

        if hasfield(typeof(vertex_data), :quality_scores)
            Test.@test !isempty(vertex_data.quality_scores)
            for qual_vec in vertex_data.quality_scores
                Test.@test isa(qual_vec, Vector{Int})
                Test.@test all(q -> 0 <= q <= 60, qual_vec)
            end
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
end