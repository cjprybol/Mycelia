# String Variable Length Graph - SingleStrand Mode with Quality Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/string_variable_singlestrand_quality_test.jl")'

import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "Quality BioSequence SingleStrand Graph" begin
    test_string = "ATCGAT"
    high_quality = [40, 40, 40, 40, 40, 40]
    qual_str = String([Char(q + 33) for q in high_quality])
    fastq_record = FASTX.FASTQ.Record("test", test_string, qual_str)
    reads = [fastq_record]

    graph = Mycelia.Rhizomorph.build_fastq_graph(reads; dataset_id="test", min_overlap=3)

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) >= 1
    Test.@test MetaGraphsNext.ne(graph) >= 0

    # Quality-aware coverage validation
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
        Test.@test !isempty(vertex_data.evidence)
        for evidence_map in values(vertex_data.evidence)
            for entries in values(evidence_map)
                for entry in entries
                    Test.@test entry isa Mycelia.Rhizomorph.QualityEvidenceEntry
                    Test.@test entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
                    Test.@test all(q -> 0 <= q - UInt8(33) <= 60, entry.quality_scores)
                end
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
