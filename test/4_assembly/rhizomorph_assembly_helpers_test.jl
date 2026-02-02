import Test
import Mycelia
import FASTX

Test.@testset "Rhizomorph Assembly Helpers" begin
    Test.@testset "Sequence type detection" begin
        fasta_reads = [FASTX.FASTA.Record("read1", "ATCG")]
        seq_type = Mycelia.Rhizomorph._detect_sequence_type(fasta_reads)
        Test.@test seq_type <: Mycelia.BioSequences.LongDNA

        fastq_reads = [FASTX.FASTQ.Record("read1", "ATCG", "IIII")]
        seq_type_fastq = Mycelia.Rhizomorph._detect_sequence_type(fastq_reads)
        Test.@test seq_type_fastq <: Mycelia.BioSequences.LongDNA

        string_reads = ["ATCG"]
        seq_type_string = Mycelia.Rhizomorph._detect_sequence_type(string_reads)
        Test.@test seq_type_string == String
    end

    Test.@testset "Observation preparation" begin
        fastq_reads = [FASTX.FASTQ.Record("read1", "ATCG", "IIII")]
        observations = Mycelia.Rhizomorph._prepare_observations(fastq_reads)
        Test.@test observations[1] isa FASTX.FASTA.Record
        Test.@test FASTX.FASTA.sequence(observations[1]) == "ATCG"

        mktempdir() do dir
            fasta_path = joinpath(dir, "reads.fasta")
            open(fasta_path, "w") do io
                println(io, ">read1")
                println(io, "ATCG")
            end
            file_obs = Mycelia.Rhizomorph._prepare_observations([fasta_path])
            Test.@test length(file_obs) == 1
            Test.@test file_obs[1] isa FASTX.FASTA.Record
        end
    end

    Test.@testset "Auto configuration" begin
        fasta_reads = [FASTX.FASTA.Record("read1", "ATCG")]
        config = Mycelia.Rhizomorph._auto_configure_assembly(fasta_reads; k = 3)
        Test.@test config.graph_mode == Mycelia.Rhizomorph.DoubleStrand
        Test.@test config.use_quality_scores == false

        string_reads = ["ABCD"]
        string_config = Mycelia.Rhizomorph._auto_configure_assembly(string_reads; k = 2)
        Test.@test string_config.sequence_type == String
        Test.@test string_config.graph_mode == Mycelia.Rhizomorph.SingleStrand
    end

    Test.@testset "Assembly statistics helpers" begin
        lengths = [10, 5, 5]
        Test.@test Mycelia.Rhizomorph._calculate_n_statistic(lengths, 0.5) == 10
        Test.@test Mycelia.Rhizomorph._calculate_l_statistic(lengths, 0.5) == 1
    end
end
