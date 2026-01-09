# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/dna_fasta_doublestrand_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/dna_fasta_doublestrand_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# DNA BioSequence Graph - DoubleStrand Mode Test

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "DNA BioSequence DoubleStrand Graph" begin
    # Test data
    test_dna = BioSequences.dna"ATCGATCG"
    reads = [FASTX.FASTA.Record("test", test_dna)]

    # Build graph
    singlestrand = Mycelia.Rhizomorph.build_fasta_graph(reads; dataset_id="test", min_overlap=3)
    graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(singlestrand)

    # Structure validation
    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test length(vertices) >= 1
    Test.@test MetaGraphsNext.ne(graph) >= 0

    # Coverage validation
    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.BioSequenceVertexData
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
