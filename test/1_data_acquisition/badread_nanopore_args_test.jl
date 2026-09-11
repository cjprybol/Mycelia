# Contract tests for the Badread nanopore argument builder.
#
# `simulate_nanopore_reads` used to pass Badread only --reference and
# --quantity, inheriting the error model, qscore model, identity distribution
# and read-length distribution from whichever Badread version happened to be
# installed. Its docstring nonetheless asserted a specific chemistry
# ("R10.4.1 / nanopore2023"), which was true of the installed binary but was
# documentation ABOUT Badread rather than anything the wrapper enforced.
#
# That gap is the kind that cannot fail loudly. A Badread release that changed a
# default would keep producing well-formed reads at a different error rate, and
# every simulated-read benchmark in the repo would silently move with it.
#
# These tests run with no conda environment and no Badread binary: the argument
# vector is built by a pure function, so the pinning is checkable directly.
#
# Run:
#   julia --project=. test/1_data_acquisition/badread_nanopore_args_test.jl

module BadreadNanoporeArgsTest

import Test
import Mycelia

# Pair each flag with the value it must carry, so a test failure names the
# parameter that came loose rather than just reporting a length mismatch.
const REQUIRED_DEFAULTS = [
    "--error_model" => "nanopore2023",
    "--qscore_model" => "nanopore2023",
    "--identity" => "95,99,2.5",
    "--length" => "15000,13000"
]

"""Value following `flag` in `args`, or `nothing` if the flag is absent."""
function flag_value(args, flag)
    i = findfirst(==(flag), args)
    return (i === nothing || i == length(args)) ? nothing : args[i + 1]
end

Test.@testset "Badread nanopore argument pinning" begin
    Test.@testset "every model parameter is passed explicitly" begin
        args = Mycelia._badread_nanopore_args(fasta = "ref.fna", quantity = "30x")
        for (flag, expected) in REQUIRED_DEFAULTS
            Test.@test flag_value(args, flag) == expected
        end
        Test.@test flag_value(args, "--reference") == "ref.fna"
        Test.@test flag_value(args, "--quantity") == "30x"
        Test.@test args[1:2] == ["badread", "simulate"]
    end

    Test.@testset "defaults match Badread 0.4.1/0.4.2 so pinning is behaviour-preserving" begin
        # Verified out-of-band: the pinned form and the bare form produce
        # byte-identical FASTQ for a fixed reference and seed under both
        # versions. If these literals are ever changed, that equivalence — and
        # comparability with every previously recorded ONT benchmark — breaks.
        Test.@test flag_value(
            Mycelia._badread_nanopore_args(fasta = "r.fna", quantity = "1x"),
            "--identity") == "95,99,2.5"
        Test.@test flag_value(
            Mycelia._badread_nanopore_args(fasta = "r.fna", quantity = "1x"),
            "--error_model") == "nanopore2023"
    end

    Test.@testset "seed is included only when supplied" begin
        without_seed = Mycelia._badread_nanopore_args(fasta = "r.fna", quantity = "1x")
        Test.@test !("--seed" in without_seed)

        with_seed = Mycelia._badread_nanopore_args(
            fasta = "r.fna", quantity = "1x", seed = 42)
        Test.@test flag_value(with_seed, "--seed") == "42"
    end

    Test.@testset "overrides reach the argument vector" begin
        # The R9.4.1-style settings must be expressible through the same entry
        # point, otherwise callers would go back to relying on binary defaults.
        args = Mycelia._badread_nanopore_args(
            fasta = "r.fna", quantity = "5x",
            # DISTINCT values on purpose: setting both to the same string
            # cannot detect the two being wired to each other's flag. A mutant
            # that swaps them passed all 30 assertions when both were
            # "nanopore2020".
            error_model = "nanopore2020", qscore_model = "nanopore2018",
            identity = "90,98,5", length_dist = "8000,6000")
        Test.@test flag_value(args, "--error_model") == "nanopore2020"
        Test.@test flag_value(args, "--qscore_model") == "nanopore2018"
        Test.@test flag_value(args, "--identity") == "90,98,5"
        Test.@test flag_value(args, "--length") == "8000,6000"
    end

    Test.@testset "no flag is ever emitted without a value" begin
        args = Mycelia._badread_nanopore_args(
            fasta = "r.fna", quantity = "1x", seed = 7)
        # A dangling final flag would make Badread consume the next token, or
        # error, depending on position — cheap to assert, hard to spot by eye.
        Test.@test iseven(length(args) - 2)
        for (i, token) in enumerate(args)
            if startswith(token, "--")
                Test.@test i < length(args)
                Test.@test !startswith(args[i + 1], "--")
            end
        end
    end
end

end  # module
