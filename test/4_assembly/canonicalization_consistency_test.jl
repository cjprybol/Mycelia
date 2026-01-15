# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/canonicalization_consistency_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/canonicalization_consistency_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "K-mer Canonicalization Consistency Tests (Rhizomorph)" begin
    Test.@testset "DoubleStrand DNA Canonicalization Consistency" begin
        reference_seq = BioSequences.dna"ATCGTTTT"  # Forward
        reverse_comp = BioSequences.reverse_complement(reference_seq)  # AAAACGAT

        reads = [
            FASTX.FASTA.Record("forward", reference_seq),
            FASTX.FASTA.Record("reverse", reverse_comp),
        ]

        graph = Mycelia.Rhizomorph.build_kmer_graph(
            reads,
            5;
            dataset_id="canon_ds",
            mode=:doublestrand,
        )

        Test.@test !isempty(MetaGraphsNext.labels(graph))

        total_evidence = 0
        has_both_orientations = false

        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
            dataset_evidence = values(get(vertex_data.evidence, "canon_ds", Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()))
            total_evidence += sum(length, dataset_evidence)
            strands = Set(obs.strand for obs in Iterators.flatten(dataset_evidence))
            if Mycelia.Rhizomorph.Forward in strands && Mycelia.Rhizomorph.Reverse in strands
                has_both_orientations = true
            end
        end

        Test.@test total_evidence > 0
        Test.@test has_both_orientations
    end

    Test.@testset "SingleStrand RNA Non-Canonicalization Consistency" begin
        reference_seq = BioSequences.rna"AUCGUUUU"
        reads = [FASTX.FASTA.Record("rna", reference_seq)]

        graph = Mycelia.Rhizomorph.build_kmer_graph(
            reads,
            5;
            dataset_id="canon_ss",
            mode=:singlestrand,
        )

        Test.@test !isempty(MetaGraphsNext.labels(graph))

        total_evidence = 0
        all_forward = true

        for label in MetaGraphsNext.labels(graph)
            vertex_data = graph[label]
            Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
            dataset_evidence = values(get(vertex_data.evidence, "canon_ss", Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()))
            total_evidence += sum(length, dataset_evidence)
            strands = Set(obs.strand for obs in Iterators.flatten(dataset_evidence))
            if isempty(strands)
                all_forward = false
            else
                all_forward &= all(strand == Mycelia.Rhizomorph.Forward for strand in strands)
            end
        end

        Test.@test total_evidence > 0
        Test.@test all_forward
    end
end
