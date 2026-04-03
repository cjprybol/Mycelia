import Test
import Mycelia
import FASTX
import BioSequences
import Kmers
import CodecZlib
import JLD2

Test.@testset "K-mer Analysis Additional Coverage" begin
    mktempdir() do temp_dir
        Test.@testset "Jellyfish counts loader" begin
            jellyfish_counts = joinpath(temp_dir, "counts.jf.tsv.gz")
            open(jellyfish_counts, "w") do io
                gzip_stream = CodecZlib.GzipCompressorStream(io)
                write(gzip_stream, "AAA\t3\nTTT\t5\n")
                close(gzip_stream)
            end

            table = Mycelia.load_jellyfish_counts(jellyfish_counts)
            Test.@test size(table, 1) == 2
            Test.@test table[1, 2] == 3
            Test.@test table[2, 2] == 5
            Test.@test table[1, 1] isa Kmers.DNAKmer{3}
            Test.@test table[2, 1] isa Kmers.DNAKmer{3}

            Test.@test_throws AssertionError Mycelia.load_jellyfish_counts(
                joinpath(temp_dir, "counts.tsv.gz"))
        end

        Test.@testset "Reference k-mer counting" begin
            fasta = joinpath(temp_dir, "reference.fasta")
            open(fasta, "w") do io
                FASTX.write(io, FASTX.FASTA.Record("ref", "ATGC"))
            end

            counts = Mycelia.fasta_to_reference_kmer_counts(
                kmer_type = Kmers.DNAKmer{2},
                fasta = fasta
            )

            Test.@test counts[Kmers.DNAKmer{2}("AT")] == 2
            Test.@test counts[Kmers.DNAKmer{2}("TG")] == 1
            Test.@test counts[Kmers.DNAKmer{2}("GC")] == 2
            Test.@test counts[Kmers.DNAKmer{2}("CA")] == 1
        end

        Test.@testset "Sequence count tables" begin
            dna_sequences = [
                BioSequences.LongDNA{4}("ATGC"),
                BioSequences.LongDNA{4}("ATAT")
            ]
            dense = Mycelia.biosequences_to_dense_counts_table(
                biosequences = dna_sequences,
                k = 2
            )
            Test.@test size(dense.kmer_counts_matrix, 2) == 2
            Test.@test size(dense.kmer_counts_matrix, 1) ==
                       length(Mycelia.generate_all_possible_canonical_kmers(2, Mycelia.DNA_ALPHABET))
            Test.@test sum(dense.kmer_counts_matrix[:, 1]) == 3
            Test.@test sum(dense.kmer_counts_matrix[:, 2]) == 3
            Test.@test_throws ErrorException Mycelia.biosequences_to_dense_counts_table(
                biosequences = dna_sequences,
                k = 11
            )

            rna_sequences = [
                BioSequences.LongRNA{4}("ACGU"),
                BioSequences.LongRNA{4}("AGGU")
            ]
            sparse_rna = Mycelia.biosequences_to_counts_table(
                biosequences = rna_sequences,
                k = 2
            )
            Test.@test size(sparse_rna.kmer_counts_matrix, 2) == 2
            Test.@test sum(sparse_rna.kmer_counts_matrix[:, 1]) == 3
            Test.@test sum(sparse_rna.kmer_counts_matrix[:, 2]) == 3
            Test.@test all(kmer -> kmer isa Kmers.RNAKmer{2}, sparse_rna.sorted_kmers)

            aa_sequences = [
                BioSequences.LongAA("ACDE"),
                BioSequences.LongAA("ACDF")
            ]
            sparse_aa = Mycelia.biosequences_to_counts_table(
                biosequences = aa_sequences,
                k = 2
            )
            Test.@test size(sparse_aa.kmer_counts_matrix, 2) == 2
            Test.@test sum(sparse_aa.kmer_counts_matrix[:, 1]) == 3
            Test.@test sum(sparse_aa.kmer_counts_matrix[:, 2]) == 3
            Test.@test all(kmer -> kmer isa Kmers.AAKmer{2}, sparse_aa.sorted_kmers)

            Test.@test_throws ErrorException Mycelia.biosequences_to_counts_table(
                biosequences = ["ATGC"],
                k = 2
            )
        end

        Test.@testset "Result loading edge cases" begin
            missing_path = joinpath(temp_dir, "missing_results.jld2")
            Test.@test isnothing(Mycelia.load_kmer_results(missing_path))

            incomplete_path = joinpath(temp_dir, "incomplete_results.jld2")
            JLD2.jldopen(incomplete_path, "w") do file
                file["kmers"] = [Kmers.DNAKmer{3}("AAA")]
            end
            Test.@test isnothing(Mycelia.load_kmer_results(incomplete_path))

            weird_metadata_path = joinpath(temp_dir, "weird_metadata.jld2")
            JLD2.jldopen(weird_metadata_path, "w") do file
                file["kmers"] = [Kmers.DNAKmer{3}("AAA")]
                file["counts"] = reshape(UInt16[4, 2], 1, 2)
                file["fasta_list"] = ["a.fasta", "b.fasta"]
                file["metadata/k"] = 3
                file["metadata/alphabet"] = 123
            end

            loaded = Mycelia.load_kmer_results(weird_metadata_path)
            Test.@test loaded !== nothing
            Test.@test loaded.metadata["k"] == 3
            Test.@test loaded.metadata["alphabet"] == Symbol("123")
        end

        Test.@testset "Additional helper coverage" begin
            dna_kmer = Kmers.DNAKmer{4}("ATGC")
            rna_kmer = Kmers.RNAKmer{4}("AUGC")
            aa_kmer = Kmers.AAKmer{3}("ACD")

            Test.@test Mycelia._is_nucleotide_sequence(dna_kmer)
            Test.@test Mycelia._is_nucleotide_sequence(rna_kmer)
            Test.@test !Mycelia._is_nucleotide_sequence(aa_kmer)

            Test.@test length(Mycelia._collect_kmers(dna_kmer, 2)) == 3
            Test.@test length(Mycelia._collect_kmers(rna_kmer, 2)) == 3
            Test.@test length(Mycelia._collect_kmers(aa_kmer, 2)) == 2

            rna_result = Mycelia.estimate_genome_size_from_kmers(
                BioSequences.LongRNA{4}("ACGUAC"),
                2
            )
            Test.@test rna_result["actual_size"] == 6
            Test.@test rna_result["total_kmers"] == 5
            Test.@test rna_result["estimated_genome_size"] == 4

            aa_result = Mycelia.estimate_genome_size_from_kmers(
                BioSequences.LongAA("ACDE"),
                2
            )
            Test.@test aa_result["actual_size"] == 4
            Test.@test aa_result["total_kmers"] == 3
            Test.@test aa_result["estimated_genome_size"] == 2

            fallback = Mycelia.adaptive_kmer_selection((
                coverage_estimates = Dict(3 => 1.0, 5 => 2.0, 7 => 3.0),
            ); target_coverage = 10.0)
            Test.@test fallback == [7]

            top3 = Mycelia.adaptive_kmer_selection((
                coverage_estimates = Dict(3 => 9.0, 5 => 10.5, 7 => 11.0, 11 => 9.8),
            ); target_coverage = 10.0)
            Test.@test length(top3) == 3
            Test.@test issorted(top3)
            Test.@test all(k -> k in (3, 5, 7, 11), top3)

            Test.@test_throws ErrorException Mycelia.generate_all_possible_kmers(2, [1, 2])
        end
    end
end
