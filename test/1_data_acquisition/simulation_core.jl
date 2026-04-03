import Test
import Mycelia
import BioSequences
import Random
import StableRNGs

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_fake_art_conda_runner(path::AbstractString)
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
    paired=0
    errfree=0
    samout=0

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            --out)
                outbase="\${2:-}"
                shift 2
                ;;
            --paired)
                paired=1
                shift
                ;;
            --errfree)
                errfree=1
                shift
                ;;
            --samout)
                samout=1
                shift
                ;;
            --fcov|--rcount|--seqSys|--len|--mflen|--sdev|--rndSeed|--in)
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done

    mkdir -p "\$(dirname \"\$outbase\")"
    if [[ \$paired -eq 1 ]]; then
        printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}1.fq"
        printf '@r2\\nTGCA\\n+\\n!!!!\\n' > "\${outbase}2.fq"
    else
        printf '@r1\\nACGT\\n+\\n!!!!\\n' > "\${outbase}.fq"
    fi

    if [[ \$samout -eq 1 ]]; then
        printf '@SQ\\tSN:ref\\tLN:4\\n' > "\${outbase}.sam"
    fi

    if [[ \$errfree -eq 1 ]]; then
        printf '@SQ\\tSN:ref\\tLN:4\\n' > "\${outbase}_errFree.sam"
    fi
    exit 0
fi

printf 'unexpected tool: %s\\n' "\$tool" >&2
exit 1
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_fake_art_conda_runner(f::Function)
    mktempdir() do dir
        runner_path = Mycelia.CONDA_RUNNER
        backup_path = joinpath(dir, "conda-runner-backup")
        runner_existed = isfile(runner_path)
        if runner_existed
            cp(runner_path, backup_path; force = true)
        else
            mkpath(dirname(runner_path))
        end
        write_fake_art_conda_runner(runner_path)
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

Test.@testset "Simulation Core Helpers" begin
    Test.@testset "mutate_string" begin
        Test.@test Mycelia.mutate_string("AB"; error_rate = 0.0) == "AB"

        Random.seed!(1)
        Test.@test Mycelia.mutate_string(
            "AB";
            alphabet = ['A', 'B', 'C'],
            error_rate = 1.0
        ) == "CA"

        Random.seed!(9)
        Test.@test Mycelia.mutate_string(
            "AB";
            alphabet = ['A', 'B', 'C'],
            error_rate = 1.0
        ) == "BCB"

        Random.seed!(2)
        Test.@test Mycelia.mutate_string(
            "AB";
            alphabet = ['A', 'B', 'C'],
            error_rate = 1.0
        ) == "B"

        Random.seed!(1)
        Test.@test Mycelia.mutate_string(
            "AA";
            alphabet = ['A'],
            error_rate = 1.0
        ) == "AA"
    end

    Test.@testset "mutate_dna_substitution_fraction" begin
        dna_string = "AAAAAAAA"
        mutated_string = Mycelia.mutate_dna_substitution_fraction(
            dna_string;
            fraction = 0.25,
            rng = StableRNGs.StableRNG(11)
        )
        Test.@test mutated_string isa String
        Test.@test length(mutated_string) == length(dna_string)
        Test.@test sum(dna_string[i] != mutated_string[i] for i in eachindex(dna_string)) == 2
        Test.@test all(base -> base in ('A', 'C', 'G', 'T'), mutated_string)

        dna_sequence = BioSequences.LongDNA{4}("ACGTACGT")
        mutated_sequence = Mycelia.mutate_dna_substitution_fraction(
            dna_sequence;
            fraction = 0.5,
            rng = StableRNGs.StableRNG(29)
        )
        mutated_sequence_string = String(mutated_sequence)
        original_sequence_string = String(dna_sequence)
        Test.@test mutated_sequence isa BioSequences.LongDNA{4}
        Test.@test length(mutated_sequence) == length(dna_sequence)
        Test.@test sum(original_sequence_string[i] != mutated_sequence_string[i]
                       for i in eachindex(mutated_sequence_string)) == 4
        Test.@test all(base -> base in ('A', 'C', 'G', 'T'), mutated_sequence_string)

        Test.@test Mycelia.mutate_dna_substitution_fraction(
            dna_string;
            fraction = 0.0,
            rng = StableRNGs.StableRNG(7)
        ) == dna_string

        Test.@test_throws ErrorException Mycelia.mutate_dna_substitution_fraction(
            dna_string;
            fraction = -0.1
        )
        Test.@test_throws ErrorException Mycelia.mutate_dna_substitution_fraction(
            dna_string;
            fraction = 1.1
        )
    end

    Test.@testset "printable string generators" begin
        ascii_greek = Mycelia.rand_ascii_greek_string(32)
        ascii_greek_alphabet = Set(vcat(Mycelia.PRINTABLE_ASCII_ALPHABET, Mycelia.PRINTABLE_GREEK_ALPHABET))
        Test.@test length(ascii_greek) == 32
        Test.@test all(ch -> ch in ascii_greek_alphabet, ascii_greek)

        latin1 = Mycelia.rand_latin1_string(32)
        latin1_alphabet = Set(Mycelia.PRINTABLE_LATIN1_ALPHABET)
        Test.@test length(latin1) == 32
        Test.@test all(ch -> ch in latin1_alphabet, latin1)

        unicode_string = Mycelia.rand_printable_unicode_string(16)
        Test.@test length(unicode_string) == 16
        Test.@test all(isprint, unicode_string)
        Test.@test all(ch -> !(0xD800 <= Int(ch) <= 0xDFFF), unicode_string)

        bmp_string = Mycelia.rand_bmp_printable_string(16)
        Test.@test length(bmp_string) == 16
        Test.@test all(isprint, bmp_string)
        Test.@test all(ch -> Int(ch) <= 0xFFFD, bmp_string)
    end

    Test.@testset "sample_abundance_weights" begin
        equal_weights = Mycelia.sample_abundance_weights(n_organisms = 4, balance = :equal)
        Test.@test equal_weights == fill(0.25, 4)

        random_weights = Mycelia.sample_abundance_weights(
            n_organisms = 5,
            balance = :random,
            rng = StableRNGs.StableRNG(5)
        )
        Test.@test length(random_weights) == 5
        Test.@test isapprox(sum(random_weights), 1.0; atol = 1e-12)
        Test.@test all(weight -> weight > 0, random_weights)

        log_normal_weights = Mycelia.sample_abundance_weights(
            n_organisms = 5,
            balance = :log_normal,
            rng = StableRNGs.StableRNG(17),
            lognorm_mu = 0.2,
            lognorm_sigma = 0.4
        )
        Test.@test length(log_normal_weights) == 5
        Test.@test isapprox(sum(log_normal_weights), 1.0; atol = 1e-12)
        Test.@test all(weight -> weight > 0, log_normal_weights)

        Test.@test_throws AssertionError Mycelia.sample_abundance_weights(
            n_organisms = 0,
            balance = :equal
        )
        Test.@test_throws AssertionError Mycelia.sample_abundance_weights(
            n_organisms = 2,
            balance = :unsupported
        )
    end

    Test.@testset "simulate_illumina_reads cached outputs" begin
        mktempdir() do dir
            outbase = joinpath(dir, "paired")
            for path in [
                    outbase * "1.fq.gz",
                    outbase * "2.fq.gz",
                    outbase * ".sam.gz",
                    outbase * "_errFree.sam.gz"
                ]
                write(path, "cached")
            end

            result = Mycelia.simulate_illumina_reads(
                fasta = joinpath(dir, "reference.fna"),
                coverage = 12,
                outbase = outbase,
                errfree = true,
                quiet = true
            )

            Test.@test result.forward_reads == outbase * "1.fq.gz"
            Test.@test result.reverse_reads == outbase * "2.fq.gz"
            Test.@test result.sam == outbase * ".sam.gz"
            Test.@test result.error_free_sam == outbase * "_errFree.sam.gz"
        end

        mktempdir() do dir
            outbase = joinpath(dir, "single")
            for path in [outbase * ".fq.gz", outbase * ".sam.gz"]
                write(path, "cached")
            end

            result = Mycelia.simulate_illumina_reads(
                fasta = joinpath(dir, "reference.fna"),
                read_count = 24,
                outbase = outbase,
                paired = false,
                quiet = true
            )

            Test.@test result.forward_reads == outbase * ".fq.gz"
            Test.@test result.reverse_reads === nothing
            Test.@test result.sam == outbase * ".sam.gz"
            Test.@test result.error_free_sam === nothing
        end

        Test.@test_throws ErrorException Mycelia.simulate_illumina_reads(
            fasta = "missing_reference.fna",
            quiet = true
        )
    end

    Test.@testset "simulate_illumina_reads command execution" begin
        with_fake_art_conda_runner() do
            mktempdir() do dir
                fasta = joinpath(dir, "reference.fna")
                write(fasta, ">ref\nACGT\n")

                paired = Mycelia.simulate_illumina_reads(
                    fasta = fasta,
                    coverage = 8,
                    outbase = joinpath(dir, "paired_run"),
                    errfree = true,
                    quiet = true
                )

                Test.@test isfile(paired.forward_reads)
                Test.@test paired.reverse_reads == joinpath(dir, "paired_run2.fq.gz")
                Test.@test isfile(paired.reverse_reads)
                Test.@test paired.sam == joinpath(dir, "paired_run.sam.gz")
                Test.@test isfile(paired.sam)
                Test.@test paired.error_free_sam == joinpath(dir, "paired_run_errFree.sam.gz")
                Test.@test isfile(paired.error_free_sam)

                single = Mycelia.simulate_illumina_reads(
                    fasta = fasta,
                    read_count = 12,
                    outbase = joinpath(dir, "single_run"),
                    paired = false,
                    quiet = false
                )

                Test.@test single.forward_reads == joinpath(dir, "single_run.fq.gz")
                Test.@test isfile(single.forward_reads)
                Test.@test single.reverse_reads === nothing
                Test.@test single.sam == joinpath(dir, "single_run.sam.gz")
                Test.@test isfile(single.sam)
                Test.@test single.error_free_sam === nothing
            end
        end
    end
end
