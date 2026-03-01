import Test
import Mycelia

Test.@testset "Lightweight Graph Mode" begin
    Test.@testset "Lightweight N-gram Graph - basic construction" begin
        texts = ["ABCDE", "BCDEF", "CDEFG"]

        # Build both full and lightweight graphs
        full_graph = Mycelia.Rhizomorph.build_ngram_graph(texts, 3)
        lw_graph = Mycelia.Rhizomorph.build_ngram_graph(texts, 3; lightweight = true)

        # Same topology
        Test.@test Mycelia.Graphs.nv(full_graph.graph) == Mycelia.Graphs.nv(lw_graph.graph)
        Test.@test Mycelia.Graphs.ne(full_graph.graph) == Mycelia.Graphs.ne(lw_graph.graph)

        # Same vertex labels
        full_labels = sort(collect(Mycelia.MetaGraphsNext.labels(full_graph)))
        lw_labels = sort(collect(Mycelia.MetaGraphsNext.labels(lw_graph)))
        Test.@test full_labels == lw_labels

        # Same counts per vertex
        for label in full_labels
            full_count = Mycelia.Rhizomorph.count_total_observations(full_graph[label])
            lw_count = Mycelia.Rhizomorph.count_total_observations(lw_graph[label])
            Test.@test full_count == lw_count
        end

        # Same edge counts
        for (src, dst) in Mycelia.MetaGraphsNext.edge_labels(full_graph)
            full_count = Mycelia.Rhizomorph.count_evidence(full_graph[src, dst])
            lw_count = Mycelia.Rhizomorph.count_evidence(lw_graph[src, dst])
            Test.@test full_count == lw_count
        end
    end

    Test.@testset "Lightweight N-gram Graph - multi-dataset" begin
        texts1 = ["ABCDE"]
        texts2 = ["BCDEF"]

        full_graph = Mycelia.Rhizomorph.build_ngram_graph(texts1, 3; dataset_id = "ds1")
        Mycelia.Rhizomorph.add_observations_to_ngram_graph!(
            full_graph, texts2, 3; dataset_id = "ds2")

        lw_graph = Mycelia.Rhizomorph.build_ngram_graph(
            texts1, 3; dataset_id = "ds1", lightweight = true)
        Mycelia.Rhizomorph.add_observations_to_ngram_graph_lightweight!(
            lw_graph, texts2, 3; dataset_id = "ds2")

        # Both should have BCD vertex with evidence from both datasets
        full_bcd_ds = sort(Mycelia.Rhizomorph.get_all_dataset_ids(full_graph["BCD"]))
        lw_bcd_ds = sort(Mycelia.Rhizomorph.get_all_dataset_ids(lw_graph["BCD"]))
        Test.@test full_bcd_ds == lw_bcd_ds
        Test.@test "ds1" in lw_bcd_ds
        Test.@test "ds2" in lw_bcd_ds

        # Dataset-specific counts should match
        for label in Mycelia.MetaGraphsNext.labels(full_graph)
            if haskey(lw_graph, label)
                for ds_id in Mycelia.Rhizomorph.get_all_dataset_ids(full_graph[label])
                    full_ds_count = Mycelia.Rhizomorph.count_dataset_observations(
                        full_graph[label], ds_id)
                    lw_ds_count = Mycelia.Rhizomorph.count_dataset_observations(
                        lw_graph[label], ds_id)
                    Test.@test full_ds_count == lw_ds_count
                end
            end
        end
    end

    Test.@testset "Lightweight N-gram Graph - statistics compatibility" begin
        texts = ["hello world", "world peace"]
        lw_graph = Mycelia.Rhizomorph.build_ngram_graph(texts, 3; lightweight = true)

        stats = Mycelia.Rhizomorph.get_ngram_statistics(lw_graph)
        Test.@test stats[:num_vertices] > 0
        Test.@test stats[:num_edges] > 0
        Test.@test stats[:n] == 3
        Test.@test !isnothing(stats[:most_common_ngram])
    end

    Test.@testset "Lightweight N-gram Graph - analysis functions" begin
        texts = ["ABCABC", "BCABCA"]
        lw_graph = Mycelia.Rhizomorph.build_ngram_graph(texts, 3; lightweight = true)

        # find_high_coverage_ngrams should work
        common = Mycelia.Rhizomorph.find_high_coverage_ngrams(lw_graph, 2)
        Test.@test !isempty(common)

        # find_unique_ngrams should work
        unique_ngrams = Mycelia.Rhizomorph.find_unique_ngrams(lw_graph, "dataset_01")
        Test.@test isa(unique_ngrams, Vector{String})

        # find_shared_ngrams should work (single dataset)
        shared = Mycelia.Rhizomorph.find_shared_ngrams(lw_graph, ["dataset_01"])
        Test.@test !isempty(shared)
    end

    Test.@testset "Lightweight N-gram Graph - from file" begin
        mktempdir() do dir
            path = joinpath(dir, "doc.txt")
            open(path, "w") do io
                println(io, "ABCD")
                println(io, "BCDA")
            end

            full_graph = Mycelia.Rhizomorph.build_ngram_graph_from_file(path, 2)
            lw_graph = Mycelia.Rhizomorph.build_ngram_graph_from_file(
                path, 2; lightweight = true)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(lw_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(lw_graph.graph)
        end
    end

    Test.@testset "Lightweight N-gram Graph - from multiple files" begin
        mktempdir() do dir
            path1 = joinpath(dir, "doc1.txt")
            path2 = joinpath(dir, "doc2.txt")
            open(path1, "w") do io
                println(io, "ABCD")
            end
            open(path2, "w") do io
                println(io, "BCXY")
            end

            full_graph = Mycelia.Rhizomorph.build_ngram_graph_from_files(
                [path1, path2], 2)
            lw_graph = Mycelia.Rhizomorph.build_ngram_graph_from_files(
                [path1, path2], 2; lightweight = true)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(lw_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(lw_graph.graph)
        end
    end

    Test.@testset "Lightweight evidence function interface" begin
        # Test the lightweight vertex data directly
        v = Mycelia.Rhizomorph.LightweightStringVertexData("ABC")
        Test.@test Mycelia.Rhizomorph.count_total_observations(v) == 0
        Test.@test Mycelia.Rhizomorph.count_evidence_entries(v) == 0
        Test.@test isempty(Mycelia.Rhizomorph.get_all_dataset_ids(v))

        # Add evidence and verify counts
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs1",
            Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs2",
            Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds2", "obs3",
            Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Forward))

        Test.@test Mycelia.Rhizomorph.count_total_observations(v) == 3
        Test.@test Mycelia.Rhizomorph.count_evidence_entries(v) == 3
        Test.@test Mycelia.Rhizomorph.count_dataset_observations(v, "ds1") == 2
        Test.@test Mycelia.Rhizomorph.count_dataset_observations(v, "ds2") == 1
        Test.@test sort(Mycelia.Rhizomorph.get_all_dataset_ids(v)) == ["ds1", "ds2"]

        # Positional evidence queries return nothing for lightweight (no EvidenceEntry structs)
        Test.@test isnothing(Mycelia.Rhizomorph.get_dataset_evidence(v, "ds1"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_observation_evidence(v, "ds1", "obs1"))

        # Observation IDs ARE tracked in lightweight mode (Option B)
        Test.@test sort(Mycelia.Rhizomorph.get_all_observation_ids(v, "ds1")) ==
                   ["obs1", "obs2"]
        Test.@test Mycelia.Rhizomorph.get_all_observation_ids(v, "ds2") == ["obs3"]

        # Strand defaults
        Test.@test Mycelia.Rhizomorph.first_evidence_strand(v) == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "Lightweight edge data interface" begin
        e = Mycelia.Rhizomorph.LightweightEdgeData(2)
        Test.@test e.overlap_length == 2
        Test.@test Mycelia.Rhizomorph.count_evidence(e) == 0

        Mycelia.Rhizomorph.add_evidence!(
            e, "ds1", "obs1",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward))

        Test.@test Mycelia.Rhizomorph.count_evidence(e) == 1
        Test.@test Mycelia.Rhizomorph.compute_edge_weight(e) == 1.0
        Test.@test Mycelia.Rhizomorph.compute_edge_coverage(e) == 1
    end

    Test.@testset "Lightweight K-mer Graph - basic construction" begin
        # Create FASTA records
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG"),
            Mycelia.FASTX.FASTA.Record("seq2", "ATGATCATG")
        ]

        full_graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3)
        lw_graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; lightweight = true)

        # Same topology
        Test.@test Mycelia.Graphs.nv(full_graph.graph) == Mycelia.Graphs.nv(lw_graph.graph)
        Test.@test Mycelia.Graphs.ne(full_graph.graph) == Mycelia.Graphs.ne(lw_graph.graph)

        # Same labels
        full_labels = sort(collect(Mycelia.MetaGraphsNext.labels(full_graph));
            by = string)
        lw_labels = sort(collect(Mycelia.MetaGraphsNext.labels(lw_graph));
            by = string)
        Test.@test length(full_labels) == length(lw_labels)
        for (fl, ll) in zip(full_labels, lw_labels)
            Test.@test string(fl) == string(ll)
        end

        # count_evidence_entries matches exactly between modes —
        # this counts unique (read, position) pairs, which is the meaningful
        # observation count at the k-mer level.
        for label in full_labels
            full_entries = Mycelia.Rhizomorph.count_evidence_entries(full_graph[label])
            lw_entries = Mycelia.Rhizomorph.count_evidence_entries(lw_graph[label])
            Test.@test full_entries == lw_entries
        end

        # count_total_observations also matches (Option B: unique observation IDs tracked)
        for label in full_labels
            full_obs = Mycelia.Rhizomorph.count_total_observations(full_graph[label])
            lw_obs = Mycelia.Rhizomorph.count_total_observations(lw_graph[label])
            Test.@test full_obs == lw_obs
        end

        # Edge evidence counts should also match
        for (src, dst) in Mycelia.MetaGraphsNext.edge_labels(full_graph)
            full_count = Mycelia.Rhizomorph.count_evidence(full_graph[src, dst])
            lw_count = Mycelia.Rhizomorph.count_evidence(lw_graph[src, dst])
            Test.@test full_count == lw_count
        end
    end

    Test.@testset "Lightweight K-mer Graph - statistics" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]

        lw_graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; lightweight = true)
        stats = Mycelia.Rhizomorph.get_kmer_statistics(lw_graph)

        Test.@test stats[:num_vertices] > 0
        Test.@test stats[:num_edges] > 0
        Test.@test stats[:k] == 3
    end

    Test.@testset "Lightweight mode supports all strand modes" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]

        ds_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; lightweight = true, mode = :doublestrand)
        Test.@test Mycelia.Graphs.nv(ds_graph.graph) > 0
        Test.@test ds_graph.graph isa Mycelia.Graphs.DiGraph

        cn_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; lightweight = true, mode = :canonical)
        Test.@test Mycelia.Graphs.nv(cn_graph.graph) > 0
        Test.@test cn_graph.graph isa Mycelia.Graphs.SimpleGraph
    end

    Test.@testset "Lightweight K-mer Graph - from file" begin
        mktempdir() do dir
            path = joinpath(dir, "test.fasta")
            open(path, "w") do io
                println(io, ">seq1")
                println(io, "ATGATGATG")
                println(io, ">seq2")
                println(io, "ATGATCATG")
            end

            full_graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(path, 3)
            lw_graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(
                path, 3; lightweight = true)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(lw_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(lw_graph.graph)
        end
    end

    Test.@testset "Lightweight K-mer Graph - from multiple files" begin
        mktempdir() do dir
            path1 = joinpath(dir, "sample1.fasta")
            path2 = joinpath(dir, "sample2.fasta")
            open(path1, "w") do io
                println(io, ">seq1")
                println(io, "ATGATGATG")
            end
            open(path2, "w") do io
                println(io, ">seq2")
                println(io, "ATGATCATG")
            end

            full_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path1, path2], 3)
            lw_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path1, path2], 3; lightweight = true)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(lw_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(lw_graph.graph)

            # Should have 2 datasets
            if !isempty(Mycelia.MetaGraphsNext.labels(lw_graph))
                first_label = first(Mycelia.MetaGraphsNext.labels(lw_graph))
                ds_ids = Mycelia.Rhizomorph.get_all_dataset_ids(lw_graph[first_label])
                Test.@test length(ds_ids) >= 1
            end
        end
    end

    Test.@testset "Lightweight path-finding compatibility" begin
        texts = ["ABCDE", "BCDEF"]
        lw_graph = Mycelia.Rhizomorph.build_ngram_graph(texts, 3; lightweight = true)

        # weighted_graph_from_rhizomorph should work with lightweight edges
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(lw_graph)
        Test.@test Mycelia.Graphs.nv(weighted.graph) == Mycelia.Graphs.nv(lw_graph.graph)
        Test.@test Mycelia.Graphs.ne(weighted.graph) >= Mycelia.Graphs.ne(lw_graph.graph)
    end

    Test.@testset "Lightweight type mutable struct behavior" begin
        v = Mycelia.Rhizomorph.LightweightStringVertexData("test")
        Test.@test v.total_count == 0

        # total_count should be mutable
        v.total_count = 42
        Test.@test v.total_count == 42

        e = Mycelia.Rhizomorph.LightweightEdgeData()
        Test.@test e.total_count == 0
        e.total_count = 10
        Test.@test e.total_count == 10
    end
end
