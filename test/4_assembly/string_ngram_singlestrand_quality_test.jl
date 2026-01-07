# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/string_ngram_singlestrand_quality_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/string_ngram_singlestrand_quality_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# String N-gram Graph - SingleStrand Mode with Quality Test

import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "String N-gram SingleStrand Graph from FASTQ" begin
    test_string = "ATCGAT"
    high_quality = [40, 40, 40, 40, 40, 40]
    qual_str = String([Char(q + 33) for q in high_quality])
    fastq_record = FASTX.FASTQ.Record("test", test_string, qual_str)
    reads = [fastq_record]

    graph = Mycelia.Rhizomorph.build_ngram_graph([String(FASTX.sequence(reads[1]))], 3; dataset_id="test")

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) == 4  # ABC, BCD, CDE, DEF
    Test.@test MetaGraphsNext.ne(graph) == 3  # ABC->BCD, BCD->CDE, CDE->DEF

    # Evidence validation
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.StringVertexData
        Test.@test !isempty(vertex_data.evidence)
        for evidence_map in values(vertex_data.evidence)
            for entries in values(evidence_map)
                for entry in entries
                    Test.@test entry isa Mycelia.Rhizomorph.EvidenceEntry
                    Test.@test entry.strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
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
