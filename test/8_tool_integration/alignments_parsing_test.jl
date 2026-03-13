import Test
import Mycelia

Test.@testset "Alignment Parsing Helpers" begin
    Test.@testset "dnadiff report parsing" begin
        mktempdir() do dir
            report_path = joinpath(dir, "sample.report")
            open(report_path, "w") do io
                println(io, "1-to-1")
                println(io, "TotalBases 1000 900")
                println(io, "AlignedBases 800 700 80.0 77.8")
                println(io, "UnalignedBases 200 200 20.0 22.2")
                println(io, "AvgIdentity 98.5")
                println(io, "AvgIdentity(Aligned) 99.0")
                println(io, "SNPs 5")
                println(io, "Indels 2 3")
            end

            parsed = Mycelia.parse_dnadiff_report(report_path)
            Test.@test parsed.summary.avg_identity == 98.5
            Test.@test parsed.summary.snps == 5
            Test.@test parsed.distance ≈ 0.015
            Test.@test parsed.distance_aligned ≈ 0.01
        end
    end

    Test.@testset "MUMmer coords parsing and summary" begin
        mktempdir() do dir
            coords_path = joinpath(dir, "sample.coords")
            open(coords_path, "w") do io
                println(io, "1 100 5 104 100 100 99.5 90.0 90.0 ref1 qry1")
            end

            coords = Mycelia.parse_mummer_coords_table(coords_path)
            Test.@test Mycelia.DataFrames.nrow(coords) == 1
            Test.@test coords[1, "ref_start"] == 1
            Test.@test coords[1, "identity"] == 99.5

            summary = Mycelia.summarize_mummer_coords(coords; reference_length = 1000, query_length = 900)
            Test.@test summary.num_alignments == 1
            Test.@test summary.aligned_bases_ref == 100
            Test.@test summary.distance ≈ 0.005
            Test.@test summary.aligned_pct_ref ≈ 10.0
        end
    end

    Test.@testset "Alignment scoring helpers" begin
        alignment = Mycelia.assess_alignment("ATCG", "ATCG")
        Test.@test alignment.total_matches == 4
        Test.@test alignment.total_edits == 0
        Test.@test Mycelia.assess_alignment_accuracy(alignment) == 1.0

        seq = Mycelia.BioSequences.LongDNA{4}("ATG")
        alignment_result, orientation = Mycelia.assess_optimal_kmer_alignment(seq, seq)
        Test.@test alignment_result.total_matches == 3
        Test.@test orientation == true

        rc_seq = Mycelia.BioSequences.reverse_complement(seq)
        Test.@test Mycelia.is_equivalent(seq, rc_seq)
    end

    Test.@testset "HMMER parsing and extraction" begin
        mktempdir() do dir
            domtbl_path = joinpath(dir, "sample.domtblout")
            open(domtbl_path, "w") do io
                println(io, "PF14859 PF14859.1 128 seqA - 300 1e-40 120.0 0.0 1 1 1e-42 1e-42 119.5 0.0 1 120 150 271 148 273 0.98 Colicin_M")
                println(io, "PF01024 PF01024.1 210 seqB - 520 1e-30 95.0 0.0 1 1 1e-31 1e-31 94.5 0.0 5 205 301 510 298 512 0.95 Colicin")
                println(io, "PF11504 PF11504.1 180 seqB - 520 1e-10 40.0 0.0 2 2 1e-11 1e-11 39.5 0.0 8 170 40 210 38 212 0.91 Colicin_Ia")
            end

            fasta_path = joinpath(dir, "sample.faa")
            records = [
                Mycelia.FASTX.FASTA.Record("seqA annotated sequence", repeat("A", 300)),
                Mycelia.FASTX.FASTA.Record("seqB second sequence", repeat("C", 520)),
            ]
            Mycelia.write_fasta(outfile = fasta_path, records = records, gzip = false)

            parsed = Mycelia.parse_hmmer_domtblout(domtbl_path)
            Test.@test Mycelia.DataFrames.nrow(parsed) == 3
            Test.@test parsed[1, :target_name] == "PF14859"
            Test.@test parsed[1, :independent_evalue] ≈ 1e-42

            per_query = Mycelia.summarize_hmmer_domain_hits(parsed)
            Test.@test Mycelia.DataFrames.nrow(per_query) == 2
            Test.@test per_query[per_query.query_name .== "seqB", :n_unique_domains][1] == 2

            boundaries = Mycelia.summarize_hmmer_domain_boundaries(parsed)
            Test.@test Mycelia.DataFrames.nrow(boundaries) == 3
            Test.@test boundaries[boundaries.target_name .== "PF01024", :start_min][1] == 301

            extracted = Mycelia.extract_hmmer_domain_sequences(parsed; fasta = fasta_path)
            Test.@test Mycelia.DataFrames.nrow(extracted) == 3
            Test.@test extracted[1, :length] == 122
            Test.@test extracted[1, :domain_sequence] == repeat("A", 122)

            outdir = joinpath(dir, "domains")
            domain_fastas = Mycelia.write_hmmer_domain_fastas(extracted; outdir = outdir)
            Test.@test Mycelia.DataFrames.nrow(domain_fastas) == 3
            Test.@test all(isfile, domain_fastas.fasta)
        end
    end

    Test.@testset "Pfam subset extraction" begin
        mktempdir() do dir
            pfam_a = joinpath(dir, "Pfam-A.hmm")
            open(pfam_a, "w") do io
                println(io, "HMMER3/f [3.4 | Aug 2023]")
                println(io, "NAME  DomainA")
                println(io, "ACC   PF00001.1")
                println(io, "//")
                println(io, "NAME  DomainB")
                println(io, "ACC   PF00002.4")
                println(io, "//")
                println(io, "NAME  DomainC")
                println(io, "ACC   PF00003.2")
                println(io, "//")
            end

            outdir = joinpath(dir, "subset")
            combined = joinpath(outdir, "subset.hmm")
            result = Mycelia.extract_pfam_hmm_subset(
                pfam_a = pfam_a,
                accessions = ["PF00003", "PF00001"],
                outdir = outdir,
                combined_hmm = combined,
                press = false,
            )

            Test.@test length(result.hmm_files) == 2
            Test.@test all(isfile, result.hmm_files)
            Test.@test isfile(combined)
            combined_text = read(combined, String)
            Test.@test occursin("PF00001.1", combined_text)
            Test.@test occursin("PF00003.2", combined_text)
            Test.@test !occursin("PF00002.4", combined_text)
        end
    end

    Test.@testset "Minimap and label helpers" begin
        size = Mycelia.system_mem_to_minimap_index_size(system_mem_gb = 8.0, denominator = 4.0)
        Test.@test size == "2G"

        tag = Mycelia.build_sample_tag("foo bar.fastq.gz")
        Test.@test tag == "foo_bar"

        label = Mycelia.build_output_label(["/path/a.fastq.gz", "/path/b.fastq.gz"])
        Test.@test label == "a__b"
    end
end
