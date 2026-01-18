import Test
import Mycelia
import FASTX
import MetaGraphsNext
import BioSequences
import Kmers

Test.@testset "Graph Type Conversions" begin
    Test.@testset "Fixed to variable conversion" begin
        records = [FASTX.FASTA.Record("read1", "ATCG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; dataset_id="test")

        variable_graph = Mycelia.Rhizomorph.convert_fixed_to_variable(graph)
        labels = collect(MetaGraphsNext.labels(variable_graph))

        Test.@test all(label -> label isa BioSequences.LongDNA{4}, labels)
        Test.@test MetaGraphsNext.ne(variable_graph) == MetaGraphsNext.ne(graph)

        edge = first(MetaGraphsNext.edge_labels(variable_graph))
        edge_data = variable_graph[edge[1], edge[2]]
        Test.@test edge_data.overlap_length == 2

        original_label = first(MetaGraphsNext.labels(graph))
        converted_label = BioSequences.LongDNA{4}(string(original_label))
        Test.@test Mycelia.Rhizomorph.count_evidence(variable_graph[converted_label]) ==
            Mycelia.Rhizomorph.count_evidence(graph[original_label])
    end

    Test.@testset "Variable to fixed conversion" begin
        fasta_records = [
            FASTX.FASTA.Record("read1", "ATG"),
            FASTX.FASTA.Record("read2", "TGC"),
        ]
        variable_graph = Mycelia.Rhizomorph.build_fasta_graph(fasta_records; min_overlap=2)
        fixed_graph = Mycelia.Rhizomorph.convert_variable_to_fixed(variable_graph, Kmers.DNAKmer{3})

        Test.@test all(label -> label isa Kmers.DNAKmer{3}, MetaGraphsNext.labels(fixed_graph))
        Test.@test MetaGraphsNext.ne(fixed_graph) == MetaGraphsNext.ne(variable_graph)

        bad_records = [FASTX.FASTA.Record("read3", "ATCG")]
        bad_graph = Mycelia.Rhizomorph.build_fasta_graph(bad_records; min_overlap=2)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.convert_variable_to_fixed(bad_graph, Kmers.DNAKmer{3})
    end

    Test.@testset "Drop quality scores" begin
        fastq_record = FASTX.FASTQ.Record("read1", "ATCG", "IIII")
        qual_graph = Mycelia.Rhizomorph.build_qualmer_graph([fastq_record], 3; dataset_id="test")
        dropped = Mycelia.Rhizomorph.drop_quality_scores(qual_graph)

        first_label = first(MetaGraphsNext.labels(dropped))
        Test.@test dropped[first_label] isa Mycelia.Rhizomorph.KmerVertexData

        fastq_graph = Mycelia.Rhizomorph.build_fastq_graph([fastq_record]; min_overlap=3)
        dropped_var = Mycelia.Rhizomorph.drop_quality_scores(fastq_graph)
        var_label = first(MetaGraphsNext.labels(dropped_var))
        Test.@test dropped_var[var_label] isa Mycelia.Rhizomorph.BioSequenceVertexData
    end
end
