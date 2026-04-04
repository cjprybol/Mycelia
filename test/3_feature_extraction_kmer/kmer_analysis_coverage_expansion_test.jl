import Test
import Mycelia
import BioSequences
import CodecZlib
import CSV
import DataFrames
import FASTX
import JLD2
import Kmers
import SparseArrays

Test.@testset "K-mer Analysis Coverage Expansion" begin
    Test.@testset "k-mer results persistence edge cases" begin
        temp_dir = mktempdir()
        kmers = Mycelia.generate_all_possible_canonical_kmers(3, Mycelia.DNA_ALPHABET)[1:2]
        counts = reshape(UInt16[3, 5], 2, 1)
        fasta_list = ["example.fasta"]

        missing_path = joinpath(temp_dir, "missing.jld2")
        Test.@test isnothing(Mycelia.load_kmer_results(missing_path))

        no_extension_path = joinpath(temp_dir, "results")
        Mycelia.save_kmer_results(
            filename = no_extension_path,
            kmers = kmers,
            counts = counts,
            fasta_list = fasta_list,
            k = 3,
            alphabet = :DNA
        )
        appended_path = no_extension_path * ".jld2"
        Test.@test isfile(appended_path)

        loaded_without_extension = Mycelia.load_kmer_results(appended_path)
        Test.@test loaded_without_extension.kmers == kmers
        Test.@test loaded_without_extension.counts == counts
        Test.@test loaded_without_extension.fasta_list == fasta_list
        Test.@test loaded_without_extension.metadata["alphabet"] == :DNA

        missing_keys_path = joinpath(temp_dir, "missing_keys.jld2")
        JLD2.jldopen(missing_keys_path, "w") do file
            file["counts"] = counts
        end
        Test.@test isnothing(Mycelia.load_kmer_results(missing_keys_path))

        invalid_alphabet_path = joinpath(temp_dir, "invalid_alphabet.jld2")
        JLD2.jldopen(invalid_alphabet_path, "w") do file
            file["kmers"] = kmers
            file["counts"] = counts
            file["fasta_list"] = fasta_list
            file["metadata/k"] = 3
            file["metadata/alphabet"] = 17
            file["metadata/custom_field"] = "kept"
        end
        loaded_invalid_alphabet = Mycelia.load_kmer_results(invalid_alphabet_path)
        Test.@test loaded_invalid_alphabet.metadata["alphabet"] == 17
        Test.@test loaded_invalid_alphabet.metadata["k"] == 3
        Test.@test loaded_invalid_alphabet.metadata["custom_field"] == "kept"
    end

    Test.@testset "jellyfish tabular helpers" begin
        temp_dir = mktempdir()
        jellyfish_counts_path = joinpath(temp_dir, "counts.jf.tsv.gz")

        open(jellyfish_counts_path, "w") do io
            gzip_io = CodecZlib.GzipCompressorStream(io)
            try
                write(gzip_io, "AAA\t2\nAAC\t3\nAAG\t2\n")
            finally
                close(gzip_io)
            end
        end

        histogram_path = Mycelia.jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts_path)
        Test.@test isfile(histogram_path)

        histogram_table = CSV.read(histogram_path, DataFrames.DataFrame; delim = '\t')
        histogram_map = Dict(
            histogram_table[!, "number of observations"] .=> histogram_table[!, "number of kmers"]
        )
        Test.@test histogram_map == Dict(2 => 2, 3 => 1)

        Test.@test Mycelia.jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts_path) ==
                   histogram_path

        loaded_counts = Mycelia.load_jellyfish_counts(jellyfish_counts_path)
        Test.@test isa(loaded_counts, DataFrames.DataFrame)
        Test.@test Kmers.ksize(eltype(loaded_counts[!, "kmer"])) == 3
        Test.@test loaded_counts[!, "kmer"] == [
            Kmers.DNAKmer{3}("AAA"),
            Kmers.DNAKmer{3}("AAC"),
            Kmers.DNAKmer{3}("AAG")
        ]
        Test.@test loaded_counts[!, "count"] == [2, 3, 2]

        Test.@test_throws AssertionError Mycelia.load_jellyfish_counts(
            joinpath(temp_dir, "counts.tsv")
        )
    end

    Test.@testset "reference counting and sequence table helpers" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "reference.fasta")
        open(fasta_path, "w") do io
            FASTX.write(io, FASTX.FASTA.Record("ref", "ATGC"))
        end

        reference_counts = Mycelia.fasta_to_reference_kmer_counts(
            kmer_type = Kmers.DNAKmer{2},
            fasta = fasta_path
        )
        Test.@test reference_counts[Kmers.DNAKmer{2}("AT")] == 2
        Test.@test reference_counts[Kmers.DNAKmer{2}("TG")] == 1
        Test.@test reference_counts[Kmers.DNAKmer{2}("GC")] == 2
        Test.@test reference_counts[Kmers.DNAKmer{2}("CA")] == 1

        dna_sequences = [
            BioSequences.LongDNA{4}("ATGC"),
            BioSequences.LongDNA{4}("ATAT")
        ]
        dense_counts = Mycelia.biosequences_to_dense_counts_table(
            biosequences = dna_sequences,
            k = 3
        )
        Test.@test length(dense_counts.sorted_kmers) ==
                   Mycelia.determine_max_canonical_kmers(3, Mycelia.DNA_ALPHABET)
        Test.@test size(dense_counts.kmer_counts_matrix) == (32, 2)
        Test.@test vec(sum(dense_counts.kmer_counts_matrix; dims = 1)) == [2.0, 2.0]

        sparse_counts = Mycelia.biosequences_to_counts_table(
            biosequences = dna_sequences,
            k = 3
        )
        Test.@test SparseArrays.issparse(sparse_counts.kmer_counts_matrix)
        Test.@test size(sparse_counts.kmer_counts_matrix, 2) == 2
        Test.@test sum(sparse_counts.kmer_counts_matrix) == 4

        Test.@test_throws ErrorException Mycelia.biosequences_to_dense_counts_table(
            biosequences = dna_sequences,
            k = 11
        )
    end

    Test.@testset "reader overloads and genome estimation variants" begin
        temp_dir = mktempdir()
        fastq_path = joinpath(temp_dir, "reads.fastq")
        Mycelia.write_fastq(
            records = [
                FASTX.FASTQ.Record("read1", "ATGC", "IIII"),
                FASTX.FASTQ.Record("read2", "ATGA", "IIII")
            ],
            filename = fastq_path
        )

        fastq_record = FASTX.FASTQ.Record("record", "AUGC", "IIII")
        record_counts = Mycelia.count_kmers(Kmers.RNAKmer{2}, fastq_record)
        Test.@test sum(values(record_counts)) == 3

        reader = Mycelia.open_fastx(fastq_path)
        reader_counts = Mycelia.count_kmers(Kmers.DNAKmer{2}, reader)
        Test.@test reader_counts[Kmers.DNAKmer{2}("AT")] == 2
        Test.@test reader_counts[Kmers.DNAKmer{2}("TG")] == 2
        close(reader)

        rna_estimate = Mycelia.estimate_genome_size_from_kmers(
            BioSequences.LongRNA{4}("AUGCAU"),
            2
        )
        Test.@test rna_estimate["unique_kmers"] == 4
        Test.@test rna_estimate["actual_size"] == 6

        aa_estimate = Mycelia.estimate_genome_size_from_kmers(
            BioSequences.LongAA("ACDE"),
            2
        )
        Test.@test aa_estimate["total_kmers"] == 3
        Test.@test aa_estimate["estimated_genome_size"] == 2

        fastq_records = [
            FASTX.FASTQ.Record("read1", "ATGC", "IIII"),
            FASTX.FASTQ.Record("read2", "ATGA", "IIII")
        ]
        record_estimate = Mycelia.estimate_genome_size_from_kmers(fastq_records, 2)
        Test.@test record_estimate["num_records"] == 2
        Test.@test record_estimate["actual_size"] == 8
    end

    Test.@testset "saturation and adaptive selection helpers" begin
        temp_dir = mktempdir()
        fasta_path = joinpath(temp_dir, "reads.fasta")
        open(fasta_path, "w") do io
            FASTX.write(io, FASTX.FASTA.Record("seq1", "ATGCATGC"))
        end

        low_sample = Mycelia.assess_dnamer_saturation(
            [fasta_path],
            Kmers.DNAKmer{3};
            kmers_to_assess = 1,
            power = 10
        )
        Test.@test low_sample.sampling_points == [0, 1]
        Test.@test low_sample.unique_kmer_counts == [0, 0]

        saturation = Mycelia.assess_dnamer_saturation(
            [fasta_path],
            Kmers.DNAKmer{3};
            kmers_to_assess = 20,
            power = 2
        )
        Test.@test saturation.eof == true
        Test.@test last(saturation.sampling_points) == 6
        Test.@test last(saturation.unique_kmer_counts) == 2

        no_match_results = (; coverage_estimates = Dict(3 => 0.5, 5 => 0.7))
        Test.@test Mycelia.adaptive_kmer_selection(
            no_match_results;
            target_coverage = 10.0
        ) == [5]
    end
end
