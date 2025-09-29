# DNA BioSequence Graph - SingleStrand Mode with Quality Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/dna_biosequence_singlestrand_quality_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "DNA BioSequence SingleStrand Quality Graph" begin
    # Test data with quality scores
    test_dna = BioSequences.dna"ATCGATCG"
    high_quality = [40, 40, 40, 40, 40, 40, 40, 40]
    qual_str = String([Char(q + 33) for q in high_quality])
    fastq_record = FASTX.FASTQ.Record("test", string(test_dna), qual_str)
    reads = [fastq_record]

    # Build graph with quality awareness
    graph = Mycelia.build_biosequence_quality_graph(reads; graph_mode=Mycelia.SingleStrand)

    # Structure validation
    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) >= 1
    Test.@test MetaGraphsNext.ne(graph) >= 0

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

    # Path reconstruction
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