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

    Test.@testset "the default outfile is keyed to the whole profile" begin
        # `simulate_nanopore_reads` SKIPS Badread when its output path already
        # exists. The default path used to be derived from only the reference and
        # the quantity, so a call overriding the error model, the identity
        # distribution, or the seed resolved to the same path as a default-profile
        # call and silently returned that call's reads. Nothing errors: the FASTQ
        # is well-formed and the read count is right, and only the error process —
        # the variable being manipulated — is wrong.
        default_path = Mycelia._badread_nanopore_outfile("ref.fna", "30x";
            error_model = "nanopore2023", qscore_model = "nanopore2023",
            identity = "95,99,2.5", length_dist = "15000,13000", seed = nothing)

        # BACKWARD COMPATIBILITY: the pinned defaults must still produce exactly
        # the historical path, or every already-generated file silently becomes a
        # cache miss and every recorded ONT benchmark is regenerated under a new
        # name.
        Test.@test default_path == "ref.badread.nanopore_r10.30x.fq.gz"

        # Each profile parameter must move the path ON ITS OWN. Asserting only
        # that "some override differs" would pass a key that incorporated, say,
        # error_model and ignored identity.
        overrides = [
            (; error_model = "nanopore2020"),
            (; qscore_model = "nanopore2018"),
            (; identity = "90,98,5"),
            (; length_dist = "8000,6000"),
            (; seed = 42)
        ]
        paths = String[]
        for override in overrides
            settings = merge(
                (; error_model = "nanopore2023", qscore_model = "nanopore2023",
                    identity = "95,99,2.5", length_dist = "15000,13000",
                    seed = nothing),
                override)
            path = Mycelia._badread_nanopore_outfile("ref.fna", "30x"; settings...)
            Test.@test path != default_path
            push!(paths, path)
        end
        # ...and they must be distinct FROM EACH OTHER, not merely from the
        # default. A key that collapsed every non-default profile onto one
        # "non-default" path would satisfy the assertions above while still
        # serving one profile's reads for another's request.
        Test.@test length(unique(paths)) == length(paths)

        # Deterministic: the same profile must resolve to the same path across
        # processes, or the cache never hits and every call re-runs Badread.
        Test.@test Mycelia._badread_nanopore_outfile("ref.fna", "30x";
            error_model = "nanopore2020", qscore_model = "nanopore2023",
            identity = "95,99,2.5", length_dist = "15000,13000", seed = nothing) ==
                   paths[1]

        # Quantity and reference remain part of the key.
        Test.@test Mycelia._badread_nanopore_outfile("ref.fna", "50x";
            error_model = "nanopore2023", qscore_model = "nanopore2023",
            identity = "95,99,2.5", length_dist = "15000,13000",
            seed = nothing) != default_path

        # The suffix is empty for the defaults and a filesystem-safe token
        # otherwise — no path separators, no shell metacharacters.
        Test.@test Mycelia._badread_nanopore_profile_suffix(
            error_model = "nanopore2023", qscore_model = "nanopore2023",
            identity = "95,99,2.5", length_dist = "15000,13000",
            seed = nothing) == ""
        nondefault_suffix = Mycelia._badread_nanopore_profile_suffix(
            error_model = "nanopore2023", qscore_model = "nanopore2023",
            identity = "90,98,5", length_dist = "15000,13000", seed = nothing)
        Test.@test occursin(r"^\.profile-[0-9a-f]+$", nondefault_suffix)
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
