# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/sequence_graphs_next.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/sequence_graphs_next.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers

Test.@testset "Rhizomorph Sequence Graphs" begin
    dna_records = [FASTX.FASTA.Record("dna", "ATCGATCG")]
    mixed_records = [
        FASTX.FASTA.Record("dna1", "ATCGATCG"),
        FASTX.FASTA.Record("dna2", "TCGATCGA"),
    ]

    Test.@testset "Enums and evidence structs" begin
        Test.@test Mycelia.Rhizomorph.Forward isa Mycelia.Rhizomorph.StrandOrientation
        Test.@test Mycelia.Rhizomorph.Reverse isa Mycelia.Rhizomorph.StrandOrientation
        Test.@test Bool(Mycelia.Rhizomorph.Forward)
        Test.@test !Bool(Mycelia.Rhizomorph.Reverse)
        Test.@test Mycelia.Rhizomorph.SingleStrand isa Mycelia.Rhizomorph.GraphMode
        Test.@test Mycelia.Rhizomorph.DoubleStrand isa Mycelia.Rhizomorph.GraphMode

        kmer = Kmers.DNAKmer{3}("ATC")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(vertex, "ds1", "obs1", Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
        Test.@test Mycelia.Rhizomorph.count_evidence(vertex) == 1

        edge = Mycelia.Rhizomorph.KmerEdgeData()
        Mycelia.Rhizomorph.add_evidence!(edge, "ds1", "obs1", Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Reverse))
        Test.@test Mycelia.Rhizomorph.compute_edge_weight(edge) == 1
    end

    Test.@testset "Evidence helpers" begin
        kmer = Kmers.DNAKmer{3}("ATC")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(vertex, "ds", "obsA", Mycelia.Rhizomorph.EvidenceEntry(2, Mycelia.Rhizomorph.Forward))
        ds_evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex, "ds")
        Test.@test length(ds_evidence) == 1
        obs_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex, "ds", "obsA")
        Test.@test !isnothing(obs_evidence)

        edge = Mycelia.Rhizomorph.KmerEdgeData()
        Mycelia.Rhizomorph.add_evidence!(edge, "ds", "obsA", Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward))
        Test.@test Mycelia.Rhizomorph.compute_edge_weight(edge) == 1
        Test.@test Mycelia.Rhizomorph.compute_edge_coverage(edge) == 1
    end

    Test.@testset "K-mer graph construction" begin
        ss_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, 3; dataset_id="rhizo_ss", mode=:singlestrand)
        ds_graph = Mycelia.Rhizomorph.build_kmer_graph(mixed_records, 3; dataset_id="rhizo_ds", mode=:doublestrand)

        for (graph, dsid) in ((ss_graph, "rhizo_ss"), (ds_graph, "rhizo_ds"))
            Test.@test graph isa MetaGraphsNext.MetaGraph
            Test.@test !isempty(MetaGraphsNext.labels(graph))

            first_label = first(MetaGraphsNext.labels(graph))
            vdata = graph[first_label]
            Test.@test vdata isa Mycelia.Rhizomorph.KmerVertexData
            Test.@test haskey(vdata.evidence, dsid)

            edges = collect(MetaGraphsNext.edge_labels(graph))
            if !isempty(edges)
                edge_data = graph[first(edges)...]
                Test.@test Mycelia.Rhizomorph.compute_edge_weight(edge_data) >= 1
            end
        end

        ds_labels = collect(MetaGraphsNext.labels(ds_graph))
        Test.@test all(label isa Kmers.DNAKmer{3} for label in ds_labels)
        rc_present = BioSequences.reverse_complement(first(ds_labels)) in ds_labels
        Test.@test rc_present  # doublestrand includes RC orientations
    end

    Test.@testset "Graph modes and strand awareness" begin
        k = 3
        seq = FASTX.FASTA.Record("test", "ATCG")
        ds_graph = Mycelia.Rhizomorph.build_kmer_graph([seq], k; dataset_id="ds_mode", mode=:doublestrand)
        ss_graph = Mycelia.Rhizomorph.build_kmer_graph([seq], k; dataset_id="ss_mode", mode=:singlestrand)

        Test.@test length(MetaGraphsNext.labels(ds_graph)) >= length(MetaGraphsNext.labels(ss_graph))
        Test.@test !haskey(ss_graph, BioSequences.reverse_complement(first(MetaGraphsNext.labels(ss_graph))))

        ds_label = first(MetaGraphsNext.labels(ds_graph))
        ds_data = ds_graph[ds_label]
        Test.@test Mycelia.Rhizomorph.count_evidence(ds_data) > 0
    end

    Test.@testset "Canonical conversion" begin
        records = [
            FASTX.FASTA.Record("f1", "ATCGA"),
            FASTX.FASTA.Record("f2", "TCGAT"),
        ]
        ss_graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; dataset_id="canon_ss", mode=:singlestrand)
        canonical = Mycelia.Rhizomorph.convert_to_canonical(ss_graph)
        Test.@test canonical isa MetaGraphsNext.MetaGraph
        Test.@test length(MetaGraphsNext.labels(canonical)) <= length(MetaGraphsNext.labels(ss_graph))
    end

    Test.@testset "Empty and short inputs" begin
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph(FASTX.FASTA.Record[], 3; dataset_id="empty", mode=:singlestrand)

        short_graph = Mycelia.Rhizomorph.build_kmer_graph([FASTX.FASTA.Record("short", "AT")], 3; dataset_id="short", mode=:singlestrand)
        Test.@test short_graph isa MetaGraphsNext.MetaGraph
        Test.@test isempty(MetaGraphsNext.labels(short_graph))
    end

    Test.@testset "Strand-aware coverage tracking" begin
        seq1 = FASTX.FASTA.Record("seq1", "ATCGATCG")
        seq2 = FASTX.FASTA.Record("seq2", "CGATCG")
        graph = Mycelia.Rhizomorph.build_kmer_graph([seq1, seq2], 3; dataset_id="cov", mode=:doublestrand)

        for label in MetaGraphsNext.labels(graph)
            vdata = graph[label]
            Test.@test vdata isa Mycelia.Rhizomorph.KmerVertexData
            Test.@test Mycelia.Rhizomorph.count_evidence(vdata) >= 1
        end

        for edge_label in MetaGraphsNext.edge_labels(graph)
            edge_data = graph[edge_label...]
            Test.@test edge_data isa Mycelia.Rhizomorph.KmerEdgeData
            Test.@test Mycelia.Rhizomorph.compute_edge_weight(edge_data) >= 1
        end
    end

    Test.@testset "N-gram graph integration" begin
        ngram = Mycelia.Rhizomorph.build_ngram_graph(["ABCABC"], 2; dataset_id="text")
        Test.@test ngram isa MetaGraphsNext.MetaGraph
        Test.@test !isempty(MetaGraphsNext.labels(ngram))
        vdata = ngram[first(MetaGraphsNext.labels(ngram))]
        Test.@test Mycelia.Rhizomorph.count_evidence(vdata) >= 1
    end
end
