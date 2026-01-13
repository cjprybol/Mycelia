# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/7_comparative_pangenomics/gold_standard_genome_comparison.jl")'
# ```

import Test
import Mycelia
import FASTX
import DataFrames
import StableRNGs

Test.@testset "Gold-Standard Genome Comparison Tests" begin
    Test.@testset "dnadiff Report Parsing" begin
        temp_dir = mktempdir()
        report_file = joinpath(temp_dir, "sample.report")
        sample_report = """
1-to-1
TotalBases  1000  900
AlignedBases  800  700  80.0  77.8
UnalignedBases  200  200  20.0  22.2
AvgIdentity  99.5
AvgIdentity(Aligned)  99.7
SNPs  10
Indels  5  20
"""
        write(report_file, sample_report)

        parsed = Mycelia.parse_dnadiff_report(report_file)
        summary = parsed.summary

        Test.@test summary.total_bases_ref == 1000
        Test.@test summary.total_bases_query == 900
        Test.@test summary.aligned_bases_ref == 800
        Test.@test summary.aligned_bases_query == 700
        Test.@test summary.avg_identity ≈ 99.5
        Test.@test summary.avg_identity_aligned ≈ 99.7
        Test.@test parsed.distance ≈ 0.005
        Test.@test parsed.distance_aligned ≈ 0.003

        rm(temp_dir; recursive=true, force=true)
    end

    Test.@testset "Genome Fragmentation" begin
        temp_dir = mktempdir()
        input_fasta = joinpath(temp_dir, "input.fasta")
        output_fasta = joinpath(temp_dir, "fragments.fasta")

        write(input_fasta, ">seq1\nACGTACGTAC\n>seq2\nTTGACCA\n")

        result_path = Mycelia.fragment_genome(input_fasta, 4; out_path=output_fasta)
        Test.@test result_path == output_fasta
        Test.@test isfile(output_fasta)

        records = open(output_fasta) do io
            collect(FASTX.FASTA.Reader(io))
        end

        Test.@test length(records) == 5
        Test.@test FASTX.identifier(records[1]) == "seq1_frag_1"
        Test.@test FASTX.identifier(records[2]) == "seq1_frag_2"
        Test.@test FASTX.identifier(records[3]) == "seq1_frag_3"
        Test.@test FASTX.identifier(records[4]) == "seq2_frag_1"
        Test.@test FASTX.identifier(records[5]) == "seq2_frag_2"
        Test.@test length(FASTX.sequence(records[1])) == 4
        Test.@test length(FASTX.sequence(records[2])) == 4
        Test.@test length(FASTX.sequence(records[3])) == 2
        Test.@test length(FASTX.sequence(records[4])) == 4
        Test.@test length(FASTX.sequence(records[5])) == 3

        rm(temp_dir; recursive=true, force=true)
    end

    run_all = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
    run_external = run_all || lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
    conda_available = haskey(ENV, "CONDA_PREFIX") || isfile(Mycelia.CONDA_RUNNER)

    Test.@testset "Gold-Standard ANI (external)" begin
        if !(run_external && conda_available)
            Test.@test_skip "Gold-standard ANI requires MYCELIA_RUN_EXTERNAL=true and a working conda runner."
        else
            rng = StableRNGs.StableRNG(42)
            dna_alphabet = ['A', 'C', 'G', 'T']

            random_dna = len -> String(rand(rng, dna_alphabet, len))
            function mutate_dna(seq::AbstractString, mutation_rate::Float64)
                seq_chars = collect(seq)
                for i in eachindex(seq_chars)
                    if rand(rng) < mutation_rate
                        original = seq_chars[i]
                        candidates = filter(base -> base != original, dna_alphabet)
                        seq_chars[i] = candidates[rand(rng, 1:length(candidates))]
                    end
                end
                return String(seq_chars)
            end

            temp_dir = mktempdir()
            query_fasta = joinpath(temp_dir, "query.fasta")
            reference_fasta = joinpath(temp_dir, "reference.fasta")

            ref_seq1 = random_dna(2048)
            ref_seq2 = random_dna(1536)
            query_seq1 = mutate_dna(ref_seq1, 0.01)
            query_seq2 = mutate_dna(ref_seq2, 0.01)

            write(reference_fasta, ">reference_1\n$(ref_seq1)\n>reference_2\n$(ref_seq2)\n")
            write(query_fasta, ">query_1\n$(query_seq1)\n>query_2\n$(query_seq2)\n")

            ani, df = Mycelia.calculate_gold_standard_ani(
                query=query_fasta,
                reference=reference_fasta,
                outdir=joinpath(temp_dir, "ani"),
                threads=1,
                fragment_size=512,
                min_coverage_pct=70.0,
                min_identity_pct=30.0,
                task="blastn",
                force=true
            )

            Test.@test ani >= 95.0
            Test.@test ani <= 100.0
            Test.@test DataFrames.nrow(df) >= 2

            rm(temp_dir; recursive=true, force=true)
        end
    end

    Test.@testset "Gold-Standard AAI (external)" begin
        if !(run_external && conda_available)
            Test.@test_skip "Gold-standard AAI requires MYCELIA_RUN_EXTERNAL=true and a working conda runner."
        else
            rng = StableRNGs.StableRNG(7)
            aa_alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

            random_protein = len -> String(rand(rng, aa_alphabet, len))
            function mutate_protein(seq::AbstractString, mutation_rate::Float64)
                seq_chars = collect(seq)
                for i in eachindex(seq_chars)
                    if rand(rng) < mutation_rate
                        original = seq_chars[i]
                        candidates = filter(residue -> residue != original, aa_alphabet)
                        seq_chars[i] = candidates[rand(rng, 1:length(candidates))]
                    end
                end
                return String(seq_chars)
            end

            temp_dir = mktempdir()
            query_fasta = joinpath(temp_dir, "query_proteins.fasta")
            reference_fasta = joinpath(temp_dir, "reference_proteins.fasta")

            ref_p1 = random_protein(220)
            ref_p2 = random_protein(180)
            query_p1 = mutate_protein(ref_p1, 0.02)
            query_p2 = mutate_protein(ref_p2, 0.02)

            write(reference_fasta, ">ref_p1\n$(ref_p1)\n>ref_p2\n$(ref_p2)\n")
            write(query_fasta, ">query_p1\n$(query_p1)\n>query_p2\n$(query_p2)\n")

            aai, df = Mycelia.calculate_gold_standard_aai(
                query_proteins=query_fasta,
                reference_proteins=reference_fasta,
                outdir=joinpath(temp_dir, "aai"),
                threads=1,
                min_coverage_pct=70.0,
                min_identity_pct=30.0,
                force=true
            )

            Test.@test aai >= 90.0
            Test.@test aai <= 100.0
            Test.@test DataFrames.nrow(df) >= 2

            rm(temp_dir; recursive=true, force=true)
        end
    end

    Test.@testset "MUMmer NUCmer Wrapper (external)" begin
        if !(run_external && conda_available)
            Test.@test_skip "NUCmer wrapper test requires MYCELIA_RUN_EXTERNAL=true and a working conda runner."
        else
            temp_dir = mktempdir()
            reference_fasta = joinpath(temp_dir, "reference.fasta")
            query_fasta = joinpath(temp_dir, "query.fasta")
            outdir = joinpath(temp_dir, "nucmer_out")

            rng = StableRNGs.StableRNG(99)
            dna_alphabet = ['A', 'C', 'G', 'T']
            seq = String(rand(rng, dna_alphabet, 2000))
            write(reference_fasta, ">reference\n$(seq)\n")
            write(query_fasta, ">query\n$(seq)\n")

            delta_path = Mycelia.run_nucmer(
                reference=reference_fasta,
                query=query_fasta,
                outdir=outdir,
                prefix="test",
                threads=1,
                min_match=10
            )

            Test.@test isfile(delta_path)
            Test.@test endswith(delta_path, ".delta")
            Test.@test filesize(delta_path) > 0

            rm(temp_dir; recursive=true, force=true)
        end
    end

    Test.@testset "MUMmer Genome Distance (external)" begin
        if !(run_external && conda_available)
            Test.@test_skip "Genome distance test requires MYCELIA_RUN_EXTERNAL=true and a working conda runner."
        else
            temp_dir = mktempdir()
            reference_fasta = joinpath(temp_dir, "reference.fasta")
            query_fasta = joinpath(temp_dir, "query.fasta")

            rng = StableRNGs.StableRNG(123)
            dna_alphabet = ['A', 'C', 'G', 'T']
            ref_seq = String(rand(rng, dna_alphabet, 2500))
            query_seq = ref_seq
            write(reference_fasta, ">reference\n$(ref_seq)\n")
            write(query_fasta, ">query\n$(query_seq)\n")

            results = Mycelia.calculate_mummer_genome_distance(
                reference=reference_fasta,
                query=query_fasta,
                outdir=joinpath(temp_dir, "distance"),
                prefix="distance",
                force=true
            )

            Test.@test isfile(results.dnadiff.report)
            Test.@test isfile(results.coords)
            Test.@test DataFrames.nrow(results.coords_table) > 0
            Test.@test results.report_metrics.distance isa Float64
            Test.@test results.report_metrics.distance <= 0.01

            rm(temp_dir; recursive=true, force=true)
        end
    end
end
