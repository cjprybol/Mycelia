import Test
import Mycelia
import FASTX
import BioSequences
import Kmers
import CodecZlib
import JLD2
import DelimitedFiles

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

        Test.@testset "Dense in-memory counts helper" begin
            fasta_a = joinpath(temp_dir, "dense_a.fasta")
            fasta_b = joinpath(temp_dir, "dense_b.fasta")
            records_a = [FASTX.FASTA.Record("a1", "ATGCAT")]
            records_b = [FASTX.FASTA.Record("b1", "ATATAT")]
            Mycelia.write_fasta(
                outfile = fasta_a,
                records = records_a,
                gzip = false,
                show_progress = false
            )
            Mycelia.write_fasta(
                outfile = fasta_b,
                records = records_b,
                gzip = false,
                show_progress = false
            )

            sorted_kmers = sort(Mycelia.generate_all_possible_canonical_kmers(3, Mycelia.DNA_ALPHABET))
            result_file = joinpath(temp_dir, "dense_counts.jld2")
            dense = Mycelia._dense_kmer_counts_in_memory(
                [fasta_a, fasta_b],
                sorted_kmers,
                Kmers.DNAKmer{3},
                Mycelia.count_canonical_kmers,
                false,
                false,
                UInt8,
                result_file
            )

            Test.@test size(dense.counts) == (length(sorted_kmers), 2)
            Test.@test eltype(dense.counts) == UInt8
            Test.@test sum(dense.counts[:, 1]) == 4
            Test.@test sum(dense.counts[:, 2]) == 4
            Test.@test isfile(result_file)
            Test.@test JLD2.load_object(result_file) == dense
        end

        Test.@testset "Jellyfish histogram conversion" begin
            jellyfish_counts = joinpath(temp_dir, "synthetic_counts.tsv.gz")
            open(jellyfish_counts, "w") do io
                gzip_stream = CodecZlib.GzipCompressorStream(io)
                write(gzip_stream, "AAA\t1\nCCC\t2\nGGG\t2\nTTT\t4\n")
                close(gzip_stream)
            end

            histogram_file = Mycelia.jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts)
            Test.@test isfile(histogram_file)

            histogram = DelimitedFiles.readdlm(histogram_file, '\t', Int, '\n';
                skipstart = 1)
            Test.@test size(histogram, 1) == 3
            Test.@test vec(histogram[1, :]) == [1, 1]
            Test.@test vec(histogram[2, :]) == [2, 2]
            Test.@test vec(histogram[3, :]) == [1, 4]

            Test.@test Mycelia.jellyfish_counts_to_kmer_frequency_histogram(jellyfish_counts) ==
                       histogram_file
        end

        Test.@testset "DNA saturation helpers" begin
            fasta = joinpath(temp_dir, "saturation.fasta")
            records = [
                FASTX.FASTA.Record("sat1", "ATGCATGCATGC"),
                FASTX.FASTA.Record("sat2", "ATATATATATAT")
            ]
            Mycelia.write_fasta(
                outfile = fasta,
                records = records,
                gzip = false,
                show_progress = false
            )

            short_probe = Mycelia.assess_dnamer_saturation(
                [fasta],
                Kmers.DNAKmer{3};
                kmers_to_assess = 1,
                power = 10
            )
            Test.@test short_probe.sampling_points == [0, 1]
            Test.@test short_probe.unique_kmer_counts == [0, 0]

            full_probe = Mycelia.assess_dnamer_saturation(
                [fasta],
                Kmers.DNAKmer{3};
                kmers_to_assess = 16,
                power = 2
            )
            Test.@test full_probe.eof == false
            Test.@test length(full_probe.sampling_points) >= 3
            Test.@test full_probe.unique_kmer_counts[1] == 0
            Test.@test issorted(full_probe.sampling_points)

            chosen_k = Mycelia.assess_dnamer_saturation(
                fasta;
                power = 2,
                min_k = 3,
                max_k = 5,
                threshold = 2.0,
                kmers_to_assess = 16
            )
            Test.@test chosen_k in (3, 5)
        end

        Test.@testset "K-mer spectra analysis" begin
            forward_reads = joinpath(temp_dir, "spectra_reads.fasta")
            records = [
                FASTX.FASTA.Record("sp1", "ATGCATGCATGCATGC"),
                FASTX.FASTA.Record("sp2", "ATGCATGCATTTATTT"),
                FASTX.FASTA.Record("sp3", "GGGGATGCATGCATGC")
            ]
            Mycelia.write_fasta(
                outfile = forward_reads,
                records = records,
                gzip = false,
                show_progress = false
            )

            output_dir = joinpath(temp_dir, "spectra_output")
            mkpath(output_dir)
            ENV["GKS_WSTYPE"] = "100"
            rate = Mycelia.analyze_kmer_spectra(
                out_directory = output_dir,
                forward_reads = forward_reads,
                k = 3,
                target_coverage = 2
            )

            Test.@test 0 < rate <= 1
            Test.@test isfile(joinpath(output_dir, "peak-detected.png"))
            Test.@test isfile(joinpath(output_dir, "peak-detected.svg"))
            Test.@test isfile(joinpath(output_dir, "downsampling-rate.txt"))
        end

        Test.@testset "Generate and cache k-mer count results" begin
            fasta = joinpath(temp_dir, "generator_input.fasta")
            Mycelia.write_fasta(
                outfile = fasta,
                records = [FASTX.FASTA.Record("gen", "ATGCATGCATGC")],
                gzip = false,
                show_progress = false
            )

            generated = Mycelia.generate_and_save_kmer_counts(
                alphabet = :DNA,
                fastas = [fasta],
                k = 3,
                output_dir = temp_dir,
                filename = "generated_dense.jld2"
            )
            Test.@test isfile(generated)
            loaded = Mycelia.load_kmer_results(generated)
            Test.@test loaded !== nothing
            Test.@test loaded.metadata["k"] == 3
            Test.@test loaded.metadata["alphabet"] == :DNA

            sentinel = Mycelia.load_kmer_results(generated)
            Test.@test sentinel !== nothing
            touched_before = stat(generated).mtime
            repeated = Mycelia.generate_and_save_kmer_counts(
                alphabet = :DNA,
                fastas = [fasta],
                k = 3,
                output_dir = temp_dir,
                filename = "generated_dense.jld2"
            )
            Test.@test repeated == generated
            Test.@test stat(generated).mtime == touched_before

            aa_generated = Mycelia.generate_and_save_kmer_counts(
                alphabet = :AA,
                fastas = [fasta],
                k = 3,
                output_dir = temp_dir,
                filename = "generated_aa.jld2"
            )
            Test.@test isfile(aa_generated)
            Test.@test_throws ErrorException Mycelia.generate_and_save_kmer_counts(
                alphabet = :INVALID,
                fastas = [fasta],
                k = 3,
                output_dir = temp_dir,
                filename = "bad.jld2"
            )
        end

        Test.@testset "Save/load normalization and spectrum edge cases" begin
            save_base = joinpath(temp_dir, "normalized_kmers")
            kmers = [Kmers.DNAKmer{3}("AAA"), Kmers.DNAKmer{3}("AAC")]
            counts = reshape(UInt8[3, 1, 2, 4], 2, 2)

            Mycelia.save_kmer_results(
                filename = save_base,
                kmers = kmers,
                counts = counts,
                fasta_list = ["a.fasta", "b.fasta"],
                k = 3,
                alphabet = :DNA
            )
            normalized_path = save_base * ".jld2"
            Test.@test isfile(normalized_path)

            loaded = Mycelia.load_kmer_results(normalized_path)
            Test.@test loaded !== nothing
            Test.@test loaded.metadata["alphabet"] == :DNA
            Test.@test loaded.counts == counts

            alternate_path = joinpath(temp_dir, "normalized_kmers.txt")
            cp(normalized_path, alternate_path; force = true)
            same_loaded = Mycelia.load_kmer_results(alternate_path)
            Test.@test same_loaded !== nothing
            Test.@test same_loaded.kmers == kmers

            no_peak = Mycelia.coverage_peak_from_hist(Dict(1 => 4); min_coverage = 2)
            Test.@test ismissing(no_peak.coverage)
            Test.@test no_peak.kmers == 0

            vector_analysis = Mycelia.analyze_kmer_frequency_spectrum([1, 1, 2, 4];
                min_coverage = 3)
            Test.@test vector_analysis.histogram[1] == 2
            Test.@test vector_analysis.total_kmers == 8
            Test.@test vector_analysis.unique_kmers == 4
            Test.@test vector_analysis.peak.coverage == 4
            Test.@test vector_analysis.peak.kmers == 1
        end

        Test.@testset "Sketching helper edge cases" begin
            dna_sequence = BioSequences.LongDNA{4}("ATGCAT")

            Test.@test_throws ErrorException Mycelia.canonical_minimizers(
                dna_sequence, 0, 2)
            Test.@test_throws ErrorException Mycelia.canonical_minimizers(
                dna_sequence, 3, 0)
            Test.@test isempty(Mycelia.canonical_minimizers(dna_sequence, 3, 10))

            canonical_syncmers = Mycelia.open_syncmers(
                dna_sequence, 3, 2, 1; canonical = true)
            Test.@test all(kmer -> kmer isa Kmers.DNAKmer{3}, canonical_syncmers)
            Test.@test length(canonical_syncmers) <=
                       length(Mycelia._collect_kmers(dna_sequence, 3))
            Test.@test_throws ErrorException Mycelia.open_syncmers(
                dna_sequence, 3, 0, 1)
            Test.@test_throws ErrorException Mycelia.open_syncmers(
                dna_sequence, 3, 2, 3)

            canonical_strobes = Mycelia.strobemers(
                dna_sequence, 3, 1, 2; canonical = true)
            Test.@test !isempty(canonical_strobes)
            Test.@test all(strobe ->
                               strobe isa Tuple{Kmers.DNAKmer{3}, Kmers.DNAKmer{3}},
                           canonical_strobes)
            Test.@test_throws ErrorException Mycelia.strobemers(dna_sequence, 3, 0, 1)
        end

        Test.@testset "K-mer math and ladder helpers" begin
            Test.@test Mycelia.determine_max_possible_kmers(2, Mycelia.DNA_ALPHABET) == 16
            Test.@test Mycelia.determine_max_canonical_kmers(2, Mycelia.AA_ALPHABET) ==
                       length(Mycelia.AA_ALPHABET)^2
            Test.@test_throws AssertionError Mycelia.determine_max_canonical_kmers(
                2, Mycelia.DNA_ALPHABET)

            possible_aa = Mycelia.generate_all_possible_canonical_kmers(
                2, Mycelia.AA_ALPHABET)
            Test.@test length(possible_aa) == length(Mycelia.AA_ALPHABET)^2

            kmer_index = Dict(
                Kmers.DNAKmer{2}("AT") => 1,
                Kmers.DNAKmer{2}("TG") => 2,
                Kmers.DNAKmer{2}("GC") => 3
            )
            kmer_counts = Dict(
                Kmers.DNAKmer{2}("AT") => 4,
                Kmers.DNAKmer{2}("GC") => 2
            )
            Test.@test Mycelia.kmer_counts_dict_to_vector(kmer_index, kmer_counts) ==
                       [4.0, 0.0, 2.0]

            prime_sequence = Mycelia.ks(min = 3, max = 35)
            Test.@test issorted(prime_sequence)
            Test.@test all(isodd, prime_sequence)
            Test.@test 19 ∉ prime_sequence

            ladder = Mycelia.k_ladder(
                max_k = 31,
                seed_primes = [2, 3, 5, 7],
                read_length = 30,
                read_margin = 20
            )
            Test.@test issorted(ladder)
            Test.@test all(isodd, ladder)
            Test.@test minimum(ladder) >= 3
            Test.@test maximum(ladder) <= 10
        end

        Test.@testset "Window quality and multi-scale analysis" begin
            short_window = BioSequences.LongDNA{4}("AT")
            Test.@test Mycelia.calculate_window_quality(short_window, 3) == 0.0

            balanced_window = BioSequences.LongDNA{4}("ATGCATGC")
            quality = Mycelia.calculate_window_quality(balanced_window, 3)
            Test.@test 0.0 <= quality <= 1.0

            multi_k = Mycelia.multi_scale_kmer_analysis(
                [
                    BioSequences.LongDNA{4}("ATGCATGC"),
                    BioSequences.LongDNA{4}("ATGCATTT")
                ];
                prime_ks = [3, 5],
                window_size = 4,
                step_size = 2
            )
            Test.@test keys(multi_k.multi_k_counts) == keys(multi_k.coverage_estimates)
            Test.@test Set(keys(multi_k.coverage_estimates)) == Set([3, 5])
            Test.@test all(k -> haskey(multi_k.consensus_kmers, k), [3, 5])
            Test.@test all(k -> !isempty(multi_k.quality_profiles[k]), [3, 5])
            Test.@test all(v -> v >= 0.0, values(multi_k.coverage_estimates))
        end

        Test.@testset "Dense and sparse file count helpers" begin
            fasta = joinpath(temp_dir, "matrix_input.fasta")
            Mycelia.write_fasta(
                outfile = fasta,
                records = [FASTX.FASTA.Record("matrix", "ATGCATGC")],
                gzip = false,
                show_progress = false
            )

            dense_cache = joinpath(temp_dir, "dense_cache.jld2")
            dense_result = Mycelia.fasta_list_to_dense_kmer_counts(
                fasta_list = [fasta],
                k = 3,
                alphabet = :DNA,
                result_file = dense_cache,
                force_threading = false,
                force_temp_files = false,
                force_progress_bars = false
            )
            Test.@test isfile(dense_cache)
            Test.@test size(dense_result.counts, 2) == 1

            dense_cached = Mycelia.fasta_list_to_dense_kmer_counts(
                fasta_list = [fasta],
                k = 3,
                alphabet = :DNA,
                result_file = dense_cache
            )
            Test.@test dense_cached == dense_result

            sparse_cache = joinpath(temp_dir, "sparse_cache.jld2")
            sparse_result = Mycelia.fasta_list_to_sparse_kmer_counts(
                fasta_list = [fasta],
                k = 3,
                alphabet = :DNA,
                result_file = sparse_cache,
                skip_rarefaction_plot = true,
                skip_rarefaction_data = true,
                force_threading = false,
                force_temp_files = true,
                force_progress_bars = false
            )
            Test.@test isfile(sparse_cache)
            Test.@test size(sparse_result.counts, 2) == 1

            sparse_cached = Mycelia.fasta_list_to_sparse_kmer_counts(
                fasta_list = [fasta],
                k = 3,
                alphabet = :DNA,
                result_file = sparse_cache,
                skip_rarefaction_plot = true,
                skip_rarefaction_data = true
            )
            Test.@test sparse_cached == sparse_result

            Test.@test_throws ErrorException Mycelia.fasta_list_to_dense_kmer_counts(
                fasta_list = String[],
                k = 3,
                alphabet = :DNA
            )
            Test.@test_throws ErrorException Mycelia.fasta_list_to_dense_kmer_counts(
                fasta_list = [joinpath(temp_dir, "missing.fasta")],
                k = 3,
                alphabet = :DNA
            )
            Test.@test_throws ErrorException Mycelia.fasta_list_to_sparse_kmer_counts(
                fasta_list = String[],
                k = 3,
                alphabet = :DNA
            )
            Test.@test_throws ErrorException Mycelia.fasta_list_to_sparse_kmer_counts(
                fasta_list = [fasta],
                k = 3,
                alphabet = :INVALID,
                skip_rarefaction_plot = true,
                skip_rarefaction_data = true
            )
        end

        Test.@testset "Count dispatch and empty sparse results" begin
            dna_fasta = joinpath(temp_dir, "dispatch_dna.fasta")
            dna_fasta_2 = joinpath(temp_dir, "dispatch_dna_2.fasta")
            rna_fasta = joinpath(temp_dir, "dispatch_rna.fasta")
            aa_fasta = joinpath(temp_dir, "dispatch_aa.fasta")

            Mycelia.write_fasta(
                outfile = dna_fasta,
                records = [FASTX.FASTA.Record("dna1", "ATGCAT")],
                gzip = false,
                show_progress = false
            )
            Mycelia.write_fasta(
                outfile = dna_fasta_2,
                records = [FASTX.FASTA.Record("dna2", "GGGAAA")],
                gzip = false,
                show_progress = false
            )
            Mycelia.write_fasta(
                outfile = rna_fasta,
                records = [FASTX.FASTA.Record("rna1", "AUGCAU")],
                gzip = false,
                show_progress = false
            )
            Mycelia.write_fasta(
                outfile = aa_fasta,
                records = [FASTX.FASTA.Record("aa1", "ACDEFG")],
                gzip = false,
                show_progress = false
            )

            dna_record_counts = Mycelia.count_kmers(
                Kmers.DNAKmer{3},
                FASTX.FASTA.Record("dna_record", "ATGCAT")
            )
            Test.@test sum(values(dna_record_counts)) == 4

            rna_record_counts = Mycelia.count_kmers(
                Kmers.RNAKmer{3},
                FASTX.FASTA.Record("rna_record", "AUGCAU")
            )
            Test.@test sum(values(rna_record_counts)) == 4

            aa_record_counts = Mycelia.count_kmers(
                Kmers.AAKmer{2},
                FASTX.FASTA.Record("aa_record", "ACDEFG")
            )
            Test.@test sum(values(aa_record_counts)) == 5

            open(dna_fasta) do io
                reader_counts = Mycelia.count_kmers(Kmers.DNAKmer{3}, FASTX.FASTA.Reader(io))
                Test.@test sum(values(reader_counts)) == 4
            end

            merged_counts = Mycelia.count_kmers(Kmers.DNAKmer{3}, [dna_fasta, dna_fasta_2])
            Test.@test sum(values(merged_counts)) == 8

            canonical_counts = Mycelia.count_canonical_kmers(
                Kmers.DNAKmer{3},
                [dna_fasta, dna_fasta_2]
            )
            Test.@test sum(values(canonical_counts)) == 8

            short_fasta = joinpath(temp_dir, "too_short.fasta")
            Mycelia.write_fasta(
                outfile = short_fasta,
                records = [FASTX.FASTA.Record("short", "AT")],
                gzip = false,
                show_progress = false
            )

            empty_sparse_cache = joinpath(temp_dir, "empty_sparse_cache.jld2")
            empty_sparse = Mycelia.fasta_list_to_sparse_kmer_counts(
                fasta_list = [short_fasta],
                k = 3,
                alphabet = :DNA,
                result_file = empty_sparse_cache,
                skip_rarefaction_plot = true,
                skip_rarefaction_data = true,
                force_threading = false,
                force_temp_files = true,
                force_progress_bars = false
            )
            Test.@test isempty(empty_sparse.kmers)
            Test.@test size(empty_sparse.counts) == (0, 1)
            Test.@test isfile(empty_sparse_cache)
        end
    end
end
