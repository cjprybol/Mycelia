import Test
import Mycelia
import FASTX
import BioSequences
import Random

Test.@testset "Simulation Sequence Helper Coverage" begin
    Test.@testset "record observation preserves sequence typing" begin
        fasta_record = FASTX.FASTA.Record("read0", BioSequences.LongDNA{4}("ACGTACGT"))

        observed_record = Mycelia.observe(fasta_record; error_rate = 0.0)
        Test.@test FASTX.sequence(BioSequences.LongDNA{4}, observed_record) ==
                   BioSequences.LongDNA{4}("ACGTACGT")
        Test.@test length(FASTX.quality(observed_record)) == 8
    end

    Test.@testset "random_fasta_record invalid moltype" begin
        Test.@test_throws ErrorException Mycelia.random_fasta_record(
            moltype = :protein_like,
            seed = 1,
            L = 4
        )
    end

    Test.@testset "paired read helpers" begin
        reference = BioSequences.LongDNA{4}("ACGTACGTACGTACGTACGT")

        Random.seed!(1)
        reads_1, reads_2 = Mycelia.generate_paired_end_reads(
            reference,
            2,
            4,
            10;
            error_rate = 0.0
        )
        Test.@test length(reads_1) == 5
        Test.@test length(reads_2) == 5
        Test.@test all(length.(reads_1) .== 4)
        Test.@test all(length.(reads_2) .== 4)

        Test.@test Mycelia.introduce_sequencing_errors(reference, 0.0) == reference

        Random.seed!(2)
        errored = Mycelia.introduce_sequencing_errors(reference, 1.0)
        Test.@test errored isa BioSequences.LongDNA{4}
        Test.@test length(errored) > 0
        Test.@test errored != reference || length(errored) != length(reference)

        Random.seed!(3)
        noisy_reads_1, noisy_reads_2 = Mycelia.generate_paired_end_reads(
            reference,
            1,
            5,
            12;
            error_rate = 0.5
        )
        Test.@test !isempty(noisy_reads_1)
        Test.@test !isempty(noisy_reads_2)
    end

    Test.@testset "observe default tech and error paths" begin
        dna = BioSequences.LongDNA{4}("ACGTACGTACGT")
        fastq_record = FASTX.FASTQ.Record("read1", "ACGTACGT", "!!!!!!!!")

        Random.seed!(4)
        observed_dna, dna_qualities = Mycelia.observe(dna)
        Test.@test length(observed_dna) > 0
        Test.@test length(dna_qualities) == length(observed_dna)

        Random.seed!(5)
        observed_pacbio, pacbio_qualities = Mycelia.observe(dna; tech = :pacbio)
        Test.@test length(observed_pacbio) > 0
        Test.@test length(pacbio_qualities) == length(observed_pacbio)

        Random.seed!(6)
        observed_fastq = Mycelia.observe(fastq_record; error_rate = 0.0)
        Test.@test FASTX.sequence(BioSequences.LongDNA{4}, observed_fastq) ==
                   BioSequences.LongDNA{4}("ACGTACGT")

        Test.@test_throws ErrorException Mycelia.observe(dna; tech = :unknown)
    end
end
