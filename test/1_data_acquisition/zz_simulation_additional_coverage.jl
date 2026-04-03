import Test
import Mycelia
import FASTX
import BioSequences
import Distributions
import Random

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_fake_simulation_conda_runner(path::AbstractString)
    script = """
#!/usr/bin/env bash
set -euo pipefail

if [[ "\${1:-}" == "run" ]]; then
    shift
fi

if [[ "\${1:-}" == "--live-stream" ]]; then
    shift
fi

if [[ "\${1:-}" == "-n" ]]; then
    shift 2
fi

tool="\${1:-}"
if [[ -n "\$tool" ]]; then
    shift
fi

if [[ "\$tool" == "art_illumina" ]]; then
    outbase=""

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            --out)
                outbase="\${2:-}"
                shift 2
                ;;
            --samout|--errfree|--noALN)
                shift
                ;;
            --fcov|--rcount|--seqSys|--len|--rndSeed|--id|-ir|-dr|--in)
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done

    mkdir -p "\$(dirname \"\$outbase\")"
    printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}.fq"
    printf '@SQ\\tSN:ref\\tLN:4\\n' > "\${outbase}.sam"
    printf '@SQ\\tSN:ref\\tLN:4\\n' > "\${outbase}_errFree.sam"
    exit 0
fi

if [[ "\$tool" == "badread" ]]; then
    printf '@read\\nACGT\\n+\\n!!!!\\n'
    exit 0
fi

printf 'unexpected tool: %s\\n' "\$tool" >&2
exit 1
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_fake_simulation_conda_runner(f::Function)
    mktempdir() do dir
        runner_path = Mycelia.CONDA_RUNNER
        backup_path = joinpath(dir, "conda-runner-backup")
        runner_existed = isfile(runner_path)
        if runner_existed
            cp(runner_path, backup_path; force = true)
        else
            mkpath(dirname(runner_path))
        end
        write_fake_simulation_conda_runner(runner_path)
        try
            return f()
        finally
            if runner_existed
                mv(backup_path, runner_path; force = true)
            else
                rm(runner_path; force = true)
            end
        end
    end
end

Test.@testset "Simulation Additional Coverage" begin
    Test.@testset "matrix generators" begin
        Random.seed!(11)
        binary = Mycelia.generate_binary_matrix(3, 4, 0.0)
        Test.@test size(binary) == (3, 4)
        Test.@test eltype(binary) == Bool
        Test.@test !any(binary)

        poisson = Mycelia.generate_poisson_matrix(2, 5, 0.0)
        Test.@test size(poisson) == (2, 5)
        Test.@test all(==(0), poisson)
    end

    Test.@testset "simulate_variants record helper" begin
        Random.seed!(7)
        record = FASTX.FASTA.Record("seq1", BioSequences.LongDNA{4}("ACGTACGT"))
        vcf_table = Mycelia.simulate_variants(
            record;
            n_variants = 1,
            window_size = 4,
            variant_size_disbribution = Distributions.Geometric(1.0),
            variant_type_likelihoods = [:substitution => 1.0]
        )

        Test.@test Mycelia.DataFrames.nrow(vcf_table) == 1
        Test.@test vcf_table[1, "#CHROM"] == "seq1"
        Test.@test vcf_table[1, "FILTER"] == "substitution"
        Test.@test vcf_table[1, "REF"] != vcf_table[1, "ALT"]
    end

    Test.@testset "observe helpers" begin
        dna = BioSequences.LongDNA{4}("ACGTACGT")
        rna = BioSequences.LongRNA{4}("ACGUACGU")
        amino = BioSequences.LongAA("MKWV")

        Random.seed!(5)
        observed_dna, dna_quals = Mycelia.observe(dna; error_rate = 0.0, tech = :illumina)
        Test.@test observed_dna == dna
        Test.@test length(dna_quals) == length(dna)
        Test.@test all(q -> 20 <= q <= 40, dna_quals)

        Random.seed!(6)
        observed_rna, rna_quals = Mycelia.observe(rna; error_rate = 0.0, tech = :nanopore)
        Test.@test observed_rna == rna
        Test.@test length(rna_quals) == length(rna)
        Test.@test all(q -> 10 <= q <= 15, rna_quals)

        Random.seed!(7)
        observed_amino, amino_quals = Mycelia.observe(amino; error_rate = 0.0, tech = :mystery)
        Test.@test observed_amino == amino
        Test.@test length(amino_quals) == length(amino)
        Test.@test all(==(30), amino_quals)

        Random.seed!(8)
        observed_record = Mycelia.observe(FASTX.FASTA.Record("read1", dna); error_rate = 0.0)
        Test.@test !isempty(FASTX.identifier(observed_record))
        Test.@test FASTX.sequence(BioSequences.LongDNA{4}, observed_record) == dna
        Test.@test length(FASTX.quality(observed_record)) == length(dna)

        Random.seed!(9)
        Test.@test 5 <= Mycelia.get_error_quality(:illumina) <= 15
        Random.seed!(10)
        Test.@test 55 <= Mycelia.get_correct_quality(:ultima, 2, 10) <= 65
        Test.@test Mycelia.get_error_quality(:unknown) == 10
        Test.@test Mycelia.get_correct_quality(:unknown, 1, 5) == 30
    end

    Test.@testset "stubbed simulation wrappers" begin
        with_fake_simulation_conda_runner() do
            mktempdir() do dir
                fasta = joinpath(dir, "reference.fna")
                write(fasta, ">ref\nACGTACGT\n")

                ultima_default = Mycelia.simulate_ultima_reads(
                    fasta = fasta,
                    coverage = 5,
                    quiet = true
                )
                Test.@test ultima_default.forward_reads == fasta * ".ultima_5x.art.fq.gz"
                Test.@test ultima_default.reverse_reads === nothing
                Test.@test isfile(ultima_default.forward_reads)
                Test.@test isfile(ultima_default.sam)
                Test.@test isfile(ultima_default.error_free_sam)

                ultima_count = Mycelia.simulate_ultima_reads(
                    fasta = fasta,
                    coverage = nothing,
                    read_count = 7,
                    outbase = joinpath(dir, "ultima_count"),
                    read_length = 200,
                    insertion_rate = 0.005,
                    deletion_rate = 0.003,
                    id = "ultima_case",
                    quiet = false,
                    rndSeed = 123
                )
                Test.@test ultima_count.forward_reads == joinpath(dir, "ultima_count.fq.gz")
                Test.@test ultima_count.reverse_reads === nothing
                Test.@test isfile(ultima_count.forward_reads)
                Test.@test isfile(ultima_count.sam)
                Test.@test isfile(ultima_count.error_free_sam)

                cached_ultima_base = joinpath(dir, "cached_ultima")
                for path in [
                        cached_ultima_base * ".fq.gz",
                        cached_ultima_base * ".sam.gz",
                        cached_ultima_base * "_errFree.sam.gz"
                    ]
                    write(path, "cached")
                end
                cached_ultima = Mycelia.simulate_ultima_reads(
                    fasta = fasta,
                    coverage = 9,
                    outbase = cached_ultima_base,
                    quiet = true
                )
                Test.@test cached_ultima.forward_reads == cached_ultima_base * ".fq.gz"
                Test.@test read(cached_ultima.forward_reads, String) == "cached"

                Test.@test_throws ErrorException Mycelia.simulate_ultima_reads(
                    fasta = fasta,
                    coverage = nothing,
                    read_count = nothing,
                    quiet = true
                )

                pacbio = Mycelia.simulate_pacbio_reads(
                    fasta = fasta,
                    quantity = "6x",
                    quiet = true,
                    seed = 22
                )
                Test.@test isfile(pacbio)
                Test.@test endswith(pacbio, ".badread.pacbio_hifi.6x.fq.gz")

                nanopore = Mycelia.simulate_nanopore_reads(
                    fasta = fasta,
                    quantity = "4x",
                    quiet = false,
                    seed = 33
                )
                Test.@test isfile(nanopore)
                Test.@test endswith(nanopore, ".badread.nanopore_r10.4x.fq.gz")

                nearly_perfect = Mycelia.simulate_nearly_perfect_long_reads(
                    fasta = fasta,
                    quantity = "3x",
                    length_mean = 12000,
                    length_sd = 3000,
                    quiet = true
                )
                Test.@test isfile(nearly_perfect)
                Test.@test endswith(nearly_perfect, ".badread.perfect.3x.fq.gz")

                nanopore_r941 = Mycelia.simulate_nanopore_r941_reads(
                    fasta = fasta,
                    quantity = "2x",
                    quiet = false
                )
                Test.@test isfile(nanopore_r941)
                Test.@test endswith(nanopore_r941, ".badread.nanopore_r941.2x.fq.gz")

                very_bad = Mycelia.simulate_very_bad_reads(
                    fasta = fasta,
                    quantity = "2x",
                    quiet = true
                )
                Test.@test isfile(very_bad)
                Test.@test endswith(very_bad, ".badread.very_bad.2x.fq.gz")

                pretty_good = Mycelia.simulate_pretty_good_reads(
                    fasta = fasta,
                    quantity = "2x",
                    quiet = false
                )
                Test.@test isfile(pretty_good)
                Test.@test endswith(pretty_good, ".badread.pretty_good.2x.fq.gz")

                custom_badread = Mycelia.simulate_badread_reads(
                    fasta = fasta,
                    quantity = "5x",
                    error_model = "pacbio2021",
                    qscore_model = "ideal",
                    seed = 44,
                    small_plasmid_bias = true,
                    outfile = joinpath(dir, "custom_badread.fq.gz"),
                    quiet = true
                )
                Test.@test custom_badread == joinpath(dir, "custom_badread.fq.gz")
                Test.@test isfile(custom_badread)

                default_badread = Mycelia.simulate_badread_reads(
                    fasta = fasta,
                    quantity = "1x",
                    quiet = false
                )
                Test.@test isfile(default_badread)
                Test.@test endswith(
                    default_badread,
                    ".badread.custom.nanopore2023.nanopore2023.1x.fq.gz"
                )

                cached_badread = joinpath(dir, "cached_badread.fq.gz")
                write(cached_badread, "cached")
                reused_badread = Mycelia.simulate_badread_reads(
                    fasta = fasta,
                    quantity = "8x",
                    outfile = cached_badread,
                    quiet = true
                )
                Test.@test reused_badread == cached_badread
                Test.@test read(reused_badread, String) == "cached"
            end
        end
    end
end
