import Test
import Mycelia

Test.@testset "Kmer Analysis Helpers" begin
    Test.@testset "Histogram and peaks" begin
        counts = Dict("AAA" => 2, "CCC" => 2, "GGG" => 3)
        hist = Mycelia.kmer_frequency_histogram(counts)
        Test.@test hist[2] == 2
        Test.@test hist[3] == 1

        vec_hist = Mycelia.kmer_frequency_histogram([1, 1, 2, 2, 2])
        Test.@test vec_hist[1] == 2
        Test.@test vec_hist[2] == 3

        peak = Mycelia.coverage_peak_from_hist(hist; min_coverage=2)
        Test.@test peak.coverage == 2
        Test.@test peak.kmers == 2

        analysis = Mycelia.analyze_kmer_frequency_spectrum(counts; min_coverage=2)
        Test.@test analysis.total_kmers == 7
        Test.@test analysis.unique_kmers == 3
    end

    Test.@testset "Minimizers and syncmers" begin
        sequence = Mycelia.BioSequences.LongDNA{4}("ATGCAT")
        minimizers = Mycelia.canonical_minimizers(sequence, 3, 2)
        Test.@test length(minimizers) == 3

        syncmers = Mycelia.open_syncmers(sequence, 3, 2, 1)
        Test.@test syncmers == [Mycelia.Kmers.DNAKmer{3}("ATG")]

        strobes = Mycelia.strobemers(sequence, 3, 1, 1)
        Test.@test length(strobes) == 3
        Test.@test strobes[1] == (Mycelia.Kmers.DNAKmer{3}("ATG"), Mycelia.Kmers.DNAKmer{3}("TGC"))
    end

    Test.@testset "Canonicalization and counting" begin
        kmer_counts = Dict(
            Mycelia.Kmers.DNAKmer{3}("CAT") => 2,
            Mycelia.Kmers.DNAKmer{3}("ATG") => 1
        )
        canonical = Mycelia.canonicalize_kmer_counts(kmer_counts)
        Test.@test length(canonical) == 1
        Test.@test canonical[Mycelia.Kmers.DNAKmer{3}("ATG")] == 3

        seq = Mycelia.BioSequences.LongDNA{4}("ATGC")
        seq_counts = Mycelia.count_kmers(Mycelia.Kmers.DNAKmer{2}, seq)
        Test.@test seq_counts[Mycelia.Kmers.DNAKmer{2}("AT")] == 1
        Test.@test seq_counts[Mycelia.Kmers.DNAKmer{2}("TG")] == 1
        Test.@test seq_counts[Mycelia.Kmers.DNAKmer{2}("GC")] == 1

        record = Mycelia.FASTX.FASTA.Record("id", "ATGC")
        record_counts = Mycelia.count_kmers(Mycelia.Kmers.DNAKmer{2}, record)
        Test.@test record_counts == seq_counts
    end

    Test.@testset "Kmer space helpers" begin
        max_canonical = Mycelia.determine_max_canonical_kmers(3, Mycelia.DNA_ALPHABET)
        Test.@test max_canonical == 32

        possible = Mycelia.generate_all_possible_kmers(1, Mycelia.DNA_ALPHABET)
        Test.@test length(possible) == length(Mycelia.DNA_ALPHABET)
        Test.@test Mycelia.Kmers.DNAKmer{1}("A") in possible

        canonical = Mycelia.generate_all_possible_canonical_kmers(1, Mycelia.DNA_ALPHABET)
        Test.@test length(canonical) == Mycelia.determine_max_canonical_kmers(1, Mycelia.DNA_ALPHABET)
    end

    Test.@testset "Genome size estimates" begin
        result = Mycelia.estimate_genome_size_from_kmers("ATGCAT", 3)
        Test.@test result["unique_kmers"] == 4
        Test.@test result["total_kmers"] == 4
        Test.@test result["estimated_genome_size"] == 2
        Test.@test result["actual_size"] == 6

        records = [
            Mycelia.FASTX.FASTA.Record("r1", "ATGC"),
            Mycelia.FASTX.FASTA.Record("r2", "ATGC")
        ]
        rec_result = Mycelia.estimate_genome_size_from_kmers(records, 2)
        Test.@test rec_result["unique_kmers"] == 3
        Test.@test rec_result["total_kmers"] == 6
        Test.@test rec_result["estimated_genome_size"] == 5
        Test.@test rec_result["actual_size"] == 8
    end

    Test.@testset "Counting across files" begin
        temp_dir = mktempdir()
        fasta_a = joinpath(temp_dir, "a.fna")
        fasta_b = joinpath(temp_dir, "b.fna")
        records = [Mycelia.FASTX.FASTA.Record("seq1", "ATGC")]
        Mycelia.write_fasta(outfile=fasta_a, records=records, gzip=false, show_progress=false)
        Mycelia.write_fasta(outfile=fasta_b, records=records, gzip=false, show_progress=false)

        file_counts = Mycelia.count_kmers(Mycelia.Kmers.DNAKmer{2}, [fasta_a, fasta_b])
        Test.@test file_counts[Mycelia.Kmers.DNAKmer{2}("AT")] == 2
        Test.@test file_counts[Mycelia.Kmers.DNAKmer{2}("TG")] == 2
        Test.@test file_counts[Mycelia.Kmers.DNAKmer{2}("GC")] == 2

        reader = Mycelia.open_fastx(fasta_a)
        reader_counts = Mycelia.count_kmers(Mycelia.Kmers.DNAKmer{2}, reader)
        close(reader)
        Test.@test reader_counts[Mycelia.Kmers.DNAKmer{2}("AT")] == 1
    end
end
