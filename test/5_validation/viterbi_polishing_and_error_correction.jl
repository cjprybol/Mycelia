# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/5_validation/viterbi_polishing_and_error_correction.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/viterbi_polishing_and_error_correction.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "Rhizomorph Traversal and Polishing Tests" begin
    Test.@testset "Traversal Variants Across Graph Types" begin
        dna_records = [FASTX.FASTA.Record("dna1", "ATCGATCGATCG")]
        for mode in (:singlestrand, :doublestrand, :canonical)
            graph = Mycelia.Rhizomorph.build_kmer_graph(dna_records, 4; dataset_id="dna", mode=mode)
            weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
            start = first(MetaGraphsNext.labels(weighted))
            prob_path = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start, 4; seed=1)
            greedy_path = Mycelia.Rhizomorph.maximum_weight_walk_next(weighted, start, 4)
            Test.@test prob_path isa Mycelia.Rhizomorph.GraphPath
            Test.@test greedy_path isa Mycelia.Rhizomorph.GraphPath
        end

        rna_records = [FASTX.FASTA.Record("rna1", "AUCGAUCGAUCG")]
        for mode in (:singlestrand, :doublestrand, :canonical)
            graph = Mycelia.Rhizomorph.build_kmer_graph(rna_records, 4; dataset_id="rna", mode=mode)
            weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
            start = first(MetaGraphsNext.labels(weighted))
            prob_path = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start, 4; seed=2)
            greedy_path = Mycelia.Rhizomorph.maximum_weight_walk_next(weighted, start, 4)
            Test.@test prob_path isa Mycelia.Rhizomorph.GraphPath
            Test.@test greedy_path isa Mycelia.Rhizomorph.GraphPath
        end

        aa_records = [FASTX.FASTA.Record("aa1", "ACDEFGHIK")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(aa_records, 3; dataset_id="aa", mode=:singlestrand)
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        start = first(MetaGraphsNext.labels(weighted))
        prob_path = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start, 3; seed=3)
        greedy_path = Mycelia.Rhizomorph.maximum_weight_walk_next(weighted, start, 3)
        Test.@test prob_path isa Mycelia.Rhizomorph.GraphPath
        Test.@test greedy_path isa Mycelia.Rhizomorph.GraphPath

        strings = ["ABCD", "BCDA", "CDAB"]
        graph = Mycelia.Rhizomorph.build_ngram_graph(strings, 3; dataset_id="strings")
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        start = first(MetaGraphsNext.labels(weighted))
        prob_path = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start, 3; seed=4)
        greedy_path = Mycelia.Rhizomorph.maximum_weight_walk_next(weighted, start, 3)
        Test.@test prob_path isa Mycelia.Rhizomorph.GraphPath
        Test.@test greedy_path isa Mycelia.Rhizomorph.GraphPath
    end

    Test.@testset "Quality-weighted Traversal" begin
        fastq_records = [FASTX.FASTQ.Record("read1", "ATCGATCG", "IIIIIIII")]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, 4; dataset_id="qual", mode=:singlestrand)
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(
            graph;
            edge_weight=Mycelia.Rhizomorph.edge_quality_weight
        )
        start = first(MetaGraphsNext.labels(weighted))
        path = Mycelia.Rhizomorph.maximum_weight_walk_next(weighted, start, 3)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
    end

    Test.@testset "K-mer Graph Properties and Coverage" begin
        test_sequence = "ATCGATCGATCGATCG"
        for k in (3, 4, 5, 6)
            records = [FASTX.FASTA.Record("seq", test_sequence)]
            graph = Mycelia.Rhizomorph.build_kmer_graph(records, k; dataset_id="props", mode=:singlestrand)
            labels = collect(MetaGraphsNext.labels(graph))
            Test.@test !isempty(labels)
            Test.@test Mycelia.Kmers.ksize(typeof(first(labels))) == k
        end

        high_records = [FASTX.FASTA.Record("seq$(i)", "ATCGATCG") for i in 1:6]
        low_records = [FASTX.FASTA.Record("seq1", "ATCGATCG")]
        high_graph = Mycelia.Rhizomorph.build_kmer_graph(high_records, 3; dataset_id="high", mode=:singlestrand)
        low_graph = Mycelia.Rhizomorph.build_kmer_graph(low_records, 3; dataset_id="low", mode=:singlestrand)

        high_weights = Dict{Tuple{Any,Any}, Float64}()
        for (src, dst) in MetaGraphsNext.edge_labels(high_graph)
            high_weights[(src, dst)] = Mycelia.Rhizomorph.compute_edge_weight(high_graph[src, dst])
        end

        low_weights = Dict{Tuple{Any,Any}, Float64}()
        for (src, dst) in MetaGraphsNext.edge_labels(low_graph)
            low_weights[(src, dst)] = Mycelia.Rhizomorph.compute_edge_weight(low_graph[src, dst])
        end

        Test.@test !isempty(high_weights)
        for (edge, low_weight) in low_weights
            Test.@test haskey(high_weights, edge)
            Test.@test high_weights[edge] > low_weight
        end
    end

    Test.@testset "Traversal Parameter Effects" begin
        records = [FASTX.FASTA.Record("seq", "ATCGATCGATCG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 4; dataset_id="seeded", mode=:singlestrand)
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        start = first(MetaGraphsNext.labels(weighted))

        path_a = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start, 4; seed=42)
        path_b = Mycelia.Rhizomorph.probabilistic_walk_next(weighted, start, 4; seed=42)
        seq_a = Mycelia.Rhizomorph.path_to_sequence(path_a, weighted)
        seq_b = Mycelia.Rhizomorph.path_to_sequence(path_b, weighted)
        Test.@test seq_a == seq_b

        custom_weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(
            graph;
            edge_weight=edge_data -> 42.0
        )
        edge = first(MetaGraphsNext.edge_labels(custom_weighted))
        Test.@test custom_weighted[edge...].weight == 42.0
    end

    Test.@testset "Traversal Error Handling" begin
        records = [FASTX.FASTA.Record("seq", "ATCGATCGATCG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 4; dataset_id="errors", mode=:singlestrand)
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        Test.@test_throws ArgumentError Mycelia.Rhizomorph.probabilistic_walk_next(
            weighted,
            "missing",
            3;
            seed=1
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.maximum_weight_walk_next(
            weighted,
            "missing",
            3
        )
    end

    Test.@testset "Graph Construction Edge Cases" begin
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph(
            FASTX.FASTA.Record[],
            3;
            dataset_id="empty",
            mode=:singlestrand
        )

        short_records = [FASTX.FASTA.Record("short", "ATC")]
        short_graph = Mycelia.Rhizomorph.build_kmer_graph(short_records, 4; dataset_id="short", mode=:singlestrand)
        Test.@test isempty(MetaGraphsNext.labels(short_graph))

        weighted_short = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(short_graph)
        short_fastq = FASTX.FASTQ.Record("short_read", "ATC", "III")
        short_result = Mycelia.process_fastq_record(
            record=short_fastq,
            graph=short_graph,
            weighted_graph=weighted_short,
            traversal=:greedy
        )
        Test.@test short_result.status == :empty_graph
        Test.@test FASTX.sequence(short_result.polished_record) == FASTX.sequence(short_fastq)
    end

    Test.@testset "FASTQ Record Processing" begin
        test_record = FASTX.FASTQ.Record("test_read", "ATCGATCG", "IIIIIIII")
        fasta_records = [FASTX.FASTA.Record("ref", "ATCGATCGATCG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(fasta_records, 4; dataset_id="ref", mode=:singlestrand)
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        result = Mycelia.process_fastq_record(
            record=test_record,
            graph=graph,
            weighted_graph=weighted,
            traversal=:greedy
        )

        Test.@test result isa NamedTuple
        Test.@test haskey(result, :polished_record)
    end

    Test.@testset "FASTQ Record Edge Cases" begin
        fasta_records = [FASTX.FASTA.Record("ref", "ATCGATCGATCG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(fasta_records, 4; dataset_id="ref", mode=:singlestrand)
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        short_record = FASTX.FASTQ.Record("short_read", "ATC", "III")
        short_result = Mycelia.process_fastq_record(
            record=short_record,
            graph=graph,
            weighted_graph=weighted,
            traversal=:greedy
        )
        Test.@test short_result.status == :short_record
        Test.@test FASTX.sequence(short_result.polished_record) == FASTX.sequence(short_record)

        test_record = FASTX.FASTQ.Record("test_read", "ATCGATCG", "IIIIIIII")
        Test.@test_throws ArgumentError Mycelia.process_fastq_record(
            record=test_record,
            graph=graph,
            weighted_graph=weighted,
            traversal=:unknown
        )

        truncated = Mycelia.process_fastq_record(
            record=test_record,
            graph=graph,
            weighted_graph=weighted,
            traversal=:greedy,
            max_steps=0
        )
        Test.@test length(FASTX.sequence(truncated.polished_record)) == 4
        Test.@test length(collect(FASTX.quality_scores(truncated.polished_record))) == 4
    end

    Test.@testset "Quality-weighted Evidence" begin
        high_records = [FASTX.FASTQ.Record("hq", "ATCGATCG", "IIIIIIII")]
        low_records = [FASTX.FASTQ.Record("lq", "ATCGATCG", "!!!!!!!!")]

        high_graph = Mycelia.Rhizomorph.build_qualmer_graph(high_records, 4; dataset_id="hq", mode=:singlestrand)
        low_graph = Mycelia.Rhizomorph.build_qualmer_graph(low_records, 4; dataset_id="lq", mode=:singlestrand)

        high_edge = first(MetaGraphsNext.edge_labels(high_graph))
        low_edge = first(MetaGraphsNext.edge_labels(low_graph))

        high_weight = Mycelia.Rhizomorph.edge_quality_weight(high_graph[high_edge...])
        low_weight = Mycelia.Rhizomorph.edge_quality_weight(low_graph[low_edge...])

        Test.@test high_weight > low_weight
    end

    Test.@testset "Polish FASTQ Pipeline" begin
        temp_fastq = tempname() * ".fastq"
        fastq_content = """@read1
ATCGATCGATCG
+
IIIIIIIIIIII
@read2
GCTAGCTAGCTA
+
IIIIIIIIIIII
"""
        write(temp_fastq, fastq_content)

        result = Mycelia.polish_fastq(
            fastq=temp_fastq,
            k=4,
            traversal=:probabilistic,
            weighting=:quality
        )

        Test.@test result isa NamedTuple
        Test.@test haskey(result, :fastq)
        Test.@test isfile(result.fastq)

        rm(temp_fastq, force=true)
        rm(result.fastq, force=true)
    end

    Test.@testset "Polish FASTQ Error Handling" begin
        Test.@test_throws ErrorException Mycelia.polish_fastq(
            fastq="does_not_exist.fastq",
            k=3
        )

        empty_fastq = tempname() * ".fastq"
        write(empty_fastq, "")
        Test.@test_throws ErrorException Mycelia.polish_fastq(
            fastq=empty_fastq,
            k=3
        )

        temp_fasta = tempname() * ".fasta"
        write(temp_fasta, ">read1\nATCG\n")
        Test.@test_throws ErrorException Mycelia.polish_fastq(
            fastq=temp_fasta,
            k=3,
            weighting=:quality
        )

        rm(empty_fastq, force=true)
        rm(temp_fasta, force=true)
    end

    Test.@testset "Iterative Polishing" begin
        temp_fastq = tempname() * ".fastq"
        fastq_content = """@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
@read2
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
"""
        write(temp_fastq, fastq_content)

        result = Mycelia.iterative_polishing(
            temp_fastq,
            13;
            traversal=:greedy,
            weighting=:evidence
        )

        Test.@test result isa Vector
        Test.@test !isempty(result)

        for entry in result
            if entry.k isa Int && isfile(entry.fastq)
                rm(entry.fastq, force=true)
            end
        end
        rm(temp_fastq, force=true)
    end

    Test.@testset "Parameter Validation" begin
        temp_fastq = tempname() * ".fastq"
        fastq_content = """@read1
ATCGATCG
+
IIIIIIII
"""
        write(temp_fastq, fastq_content)

        Test.@test_throws ArgumentError Mycelia.polish_fastq(
            fastq=temp_fastq,
            k=3,
            traversal=:unknown
        )
        Test.@test_throws ArgumentError Mycelia.polish_fastq(
            fastq=temp_fastq,
            k=3,
            weighting=:unknown
        )
        Test.@test_throws ArgumentError Mycelia.polish_fastq(
            fastq=temp_fastq,
            k=3,
            graph_mode=:unknown
        )

        rm(temp_fastq, force=true)
    end
end
