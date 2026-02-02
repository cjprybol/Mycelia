import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "Alphabet Inference" begin
    Test.@testset "Ambiguous ACGT defaults to DNA" begin
        records = [FASTX.FASTA.Record("test", "ACGT")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3)
        labels = collect(MetaGraphsNext.labels(graph))
        Test.@test !isempty(labels)
        Test.@test all(label -> label isa Kmers.DNAKmer{3}, labels)
    end

    Test.@testset "Ambiguous ACGT error policy" begin
        records = [FASTX.FASTA.Record("test", "ACGT")]
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            ambiguous_action = :error
        )
    end

    Test.@testset "FAA extension hint selects AA" begin
        mktempdir() do temp_dir
            fasta_path = joinpath(temp_dir, "hint.faa")
            open(fasta_path, "w") do io
                write(io, ">seq1\nACGT\n")
            end

            graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta_path, 3)
            labels = collect(MetaGraphsNext.labels(graph))
            Test.@test !isempty(labels)
            Test.@test all(label -> label isa Kmers.AAKmer{3}, labels)
        end
    end

    Test.@testset "FASTA graph file hints" begin
        mktempdir() do temp_dir
            fasta_path = joinpath(temp_dir, "contigs.faa")
            open(fasta_path, "w") do io
                write(io, ">seq1\nACGT\n")
            end

            graph = Mycelia.Rhizomorph.build_fasta_graph_from_file(fasta_path)
            labels = collect(MetaGraphsNext.labels(graph))
            Test.@test !isempty(labels)
            Test.@test all(label -> label isa BioSequences.LongAA, labels)
        end
    end

    Test.@testset "Conflicting file hints error" begin
        mktempdir() do temp_dir
            dna_path = joinpath(temp_dir, "sample.fna")
            aa_path = joinpath(temp_dir, "sample.faa")
            open(dna_path, "w") do io
                write(io, ">seq1\nACGT\n")
            end
            open(aa_path, "w") do io
                write(io, ">seq1\nACGT\n")
            end

            Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path, aa_path],
                3
            )
            Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_fasta_graph_from_files(
                [dna_path, aa_path]
            )
        end
    end
end
