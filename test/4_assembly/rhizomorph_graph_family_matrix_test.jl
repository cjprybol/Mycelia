import Test
import Mycelia
import FASTX

Test.@testset "Rhizomorph Graph Family Matrix" begin
    fasta_reads = [
        FASTX.FASTA.Record("r1", "ATGCA"),
        FASTX.FASTA.Record("r2", "TGCAT")
    ]
    fastq_reads = [
        FASTX.FASTQ.Record("q1", "ATGCA", "IIIII"),
        FASTX.FASTQ.Record("q2", "TGCAT", "IIIII")
    ]

    Test.@testset "K-mer family memory profiles" begin
        result_fasta = Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :kmer,
            k = 3,
            memory_profile = :ultralight
        )
        Test.@test result_fasta.assembly_stats["effective_graph_family"] == "kmer"
        Test.@test result_fasta.assembly_stats["effective_memory_profile"] == "ultralight"

        result_fastq = Mycelia.Rhizomorph.assemble_genome(
            fastq_reads;
            graph_family = :kmer,
            k = 3,
            memory_profile = :ultralight
        )
        Test.@test result_fastq.assembly_stats["effective_graph_family"] == "kmer"
        Test.@test result_fastq.assembly_stats["effective_memory_profile"] == "ultralight_quality"
    end

    Test.@testset "N-gram family across FASTA inputs" begin
        result = Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :ngram,
            k = 3,
            memory_profile = :lightweight
        )
        Test.@test result.assembly_stats["effective_graph_family"] == "ngram"
        Test.@test result.assembly_stats["effective_memory_profile"] == "lightweight"
    end

    Test.@testset "OLC family" begin
        result = Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :olc,
            min_overlap = 3,
            memory_profile = :full
        )
        Test.@test result.assembly_stats["effective_graph_family"] == "olc"
        Test.@test result.assembly_stats["effective_memory_profile"] == "full"

        Test.@test_throws ErrorException Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :olc,
            min_overlap = 3,
            memory_profile = :lightweight
        )
    end

    Test.@testset "Token family fallback tokenization" begin
        result = Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :token,
            k = 3,
            tokenizer = :bpe,
            memory_profile = :full
        )
        Test.@test result.assembly_stats["effective_graph_family"] == "token"
        Test.@test result.assembly_stats["tokenizer_used"] == "char_fallback"

        Test.@test_throws ErrorException Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :token,
            k = 3,
            memory_profile = :lightweight
        )
    end

    Test.@testset "Invalid profile combinations" begin
        Test.@test_throws ErrorException Mycelia.Rhizomorph.assemble_genome(
            fasta_reads;
            graph_family = :ngram,
            k = 3,
            memory_profile = :lightweight_quality
        )
    end
end
