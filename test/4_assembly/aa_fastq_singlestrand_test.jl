# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/aa_fastq_singlestrand_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/aa_fastq_singlestrand_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Amino Acid BioSequence Graph - SingleStrand Mode with Quality Test

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "AA BioSequence SingleStrand Quality Graph" begin
    test_aa = BioSequences.aa"MKTV"
    # TODO: use Mycelia functions to do this for additional testing depth
    high_quality = [40, 40, 40, 40]
    qual_str = String([Char(q + 33) for q in high_quality])
    expected_quality = UInt8.(high_quality .+ 33)
    fastq_record = FASTX.FASTQ.Record("test", string(test_aa), qual_str)
    reads = [fastq_record]

    graph = Mycelia.Rhizomorph.build_fastq_graph(reads; dataset_id="aa_fastq_test", min_overlap=3)

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) == 1
    Test.@test MetaGraphsNext.ne(graph) == 0

    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test hasfield(typeof(vertex_data), :evidence)
        dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
        Test.@test "aa_fastq_test" in dataset_ids
        obs_ids = Mycelia.Rhizomorph.get_all_observation_ids(vertex_data, "aa_fastq_test")
        Test.@test "test" in obs_ids
        obs_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex_data, "aa_fastq_test", "test")
        Test.@test !isnothing(obs_evidence)
        if !isnothing(obs_evidence)
            for entry in obs_evidence
                Test.@test entry.position == 1
                Test.@test entry.strand == Mycelia.Rhizomorph.Forward
                Test.@test entry.quality_scores isa Vector{UInt8}
                Test.@test length(entry.quality_scores) == length(test_aa)
                Test.@test entry.quality_scores == expected_quality
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
