# Wiring test for the Rhizomorph correction-validation sweep's technology path.
#
# This is the PROOF that the corrector actually receives the technology the
# sweep simulated with. The sweep simulates long reads with the :nanopore error
# model, but before this fix it never told the CORRECTOR that: `assemble_genome`
# was called without `sequencing_tech`, so it fell back to the :illumina
# substitution-only profile and every "nanopore" long-read row was silently
# substitution-only.
#
# Two properties are asserted here that the fast, dependency-free
# `rhizomorph_correction_validation_sweep_test.jl` cannot cover:
#
#   T3 — `run_arm` actually forwards `sequencing_tech` to the assembler, and the
#        `corrector_sequencing_tech` it reports is READ BACK from the assembler's
#        own `assembly_stats` stamp rather than echoed from the input argument.
#        A stub assembler that stamps a DIFFERENT value than it received proves
#        the read-back (an echo would report the input and pass regardless).
#
#   T5 — passing `sequencing_tech` to the NAIVE arm (corrector=:none) is inert.
#        `run_arm` passes the tech to BOTH arms. `AssemblyConfig` consults it in
#        THREE places: unconditionally at construction to VALIDATE the symbol,
#        in the `:olc` layout branch, and in the iterative corrector's error
#        profile. Only the last two can change output and the sweep pins
#        `layout=:native`, so the naive arm now merely VALIDATES a tech it
#        previously never saw. That reasoning is VERIFIED here byte-for-byte on
#        the emitted contig FASTA, not assumed.
#
# This file loads Mycelia (unlike the fast test file, which must stay
# dependency-free). Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/rhizomorph_correction_validation_sweep_wiring_test.jl

import Test
import Random
import FASTX
import BioSequences

# Defines run_arm / simulate_regime_reads / regime_for_readlen without executing
# the sweep (the script's `abspath(PROGRAM_FILE) == @__FILE__` guard).
include(joinpath(@__DIR__, "rhizomorph_correction_validation_sweep.jl"))

Test.@testset "run_arm forwards sequencing_tech and reads the stamp back" begin
    mktempdir() do dir
        # --- Stub assembler: records the kwargs it was handed, and stamps the
        # tech it RECEIVED into assembly_stats (the same key assemble_genome
        # stamps at src/rhizomorph/assembly.jl).
        recorded = Ref{Any}(nothing)
        stub = (reads;
            kwargs...) -> begin
            recorded[] = Dict(kwargs)
            return (contigs = String[], contig_names = String[],
                assembly_stats = Dict{String, Any}(
                    "sequencing_tech" => String(kwargs[:sequencing_tech]),
                    "indel_engaged" => 7))
        end

        # --- Long-read cell: the regime maps to :nanopore, which must reach the
        # assembler AND come back out on the row.
        regime, tech = regime_for_readlen(5000)
        Test.@test tech === :nanopore

        res_iter = run_arm(FASTX.FASTQ.Record[], :iterative, 21, 1_000, dir, "iter";
            sequencing_tech = tech, assembler = stub)
        Test.@test recorded[] !== nothing
        Test.@test recorded[][:sequencing_tech] === :nanopore
        Test.@test recorded[][:corrector] === :iterative
        Test.@test recorded[][:k] == 21
        Test.@test recorded[][:graph_mode] === Mycelia.Rhizomorph.DoubleStrand
        Test.@test res_iter.ok
        Test.@test res_iter.corrector_sequencing_tech == "nanopore"
        # The runtime-engagement counter is read back from the same stats dict.
        Test.@test res_iter.corrector_indel_engaged == 7

        # --- Short-read cell: the naive arm gets the tech too (see T5 for the
        # inertness proof) and reports the assembler's stamp.
        regime_short, tech_short = regime_for_readlen(150)
        Test.@test tech_short === :illumina

        res_none = run_arm(FASTX.FASTQ.Record[], :none, 21, 1_000, dir, "none";
            sequencing_tech = tech_short, assembler = stub)
        Test.@test recorded[][:sequencing_tech] === :illumina
        Test.@test recorded[][:corrector] === :none
        Test.@test res_none.corrector_sequencing_tech == "illumina"
    end
end

Test.@testset "corrector_sequencing_tech is read back, never echoed" begin
    mktempdir() do dir
        # The column's whole purpose is to be structurally incapable of reading
        # green while the wiring is broken. Prove that by stamping a value that
        # DISAGREES with the input: an echoed column would report "nanopore"
        # (the argument); a read-back column reports the stamp.
        liar = (reads;
            kwargs...) -> (contigs = String[], contig_names = String[],
            assembly_stats = Dict{String, Any}("sequencing_tech" => "illumina"))
        res = run_arm(FASTX.FASTQ.Record[], :iterative, 21, 1_000, dir, "liar";
            sequencing_tech = :nanopore, assembler = liar)
        Test.@test res.corrector_sequencing_tech == "illumina"

        # An assembler that stamps NOTHING (the real naive route's behavior)
        # reports the "n/a" sentinel rather than inventing the input value.
        silent = (reads;
            kwargs...) -> (contigs = String[], contig_names = String[],
            assembly_stats = Dict{String, Any}())
        res_silent = run_arm(FASTX.FASTQ.Record[], :none, 21, 1_000, dir, "silent";
            sequencing_tech = :nanopore, assembler = silent)
        Test.@test res_silent.corrector_sequencing_tech == "n/a"
        # An UNSTAMPED route must stay `missing`, not collapse to 0 — 0 would
        # assert "the gap moves were available and never fired", which is a
        # different (and false) claim about a route that has no gap moves.
        Test.@test res_silent.corrector_indel_engaged === missing
    end
end

Test.@testset "a failed arm records the distinct \"error\" sentinel" begin
    # An arm that THREW and a healthy naive arm both used to report "n/a",
    # disambiguated only by ok=false. Give the exception path its own sentinel so
    # the CSV distinguishes "this route does not stamp the field" from "this arm
    # never got far enough to stamp anything".
    mktempdir() do dir
        boom = (reads; kwargs...) -> error("assembler exploded")
        res = Test.@test_logs (:warn,) run_arm(
            FASTX.FASTQ.Record[], :iterative, 21, 1_000, dir, "boom";
            sequencing_tech = :nanopore, assembler = boom)
        Test.@test !res.ok
        Test.@test res.corrector_sequencing_tech == "error"
        Test.@test res.corrector_indel_engaged === missing
    end
end

Test.@testset "naive arm is inert to sequencing_tech (T5)" begin
    # `run_arm` passes sequencing_tech to BOTH arms. On corrector=:none with the
    # sweep's layout=:native that should change nothing at all. Assert it on the
    # produced bytes rather than trusting the code reading.
    rng = Random.MersenneTwister(20260726)
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = 7, L = 600)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    reads, _ = simulate_regime_reads(refseq, 120, 20, 0.01, :illumina, rng)
    Test.@test !isempty(reads)

    mktempdir() do dir
        with_dir = joinpath(dir, "with_tech")
        without_dir = joinpath(dir, "without_tech")
        mkpath(with_dir)
        mkpath(without_dir)

        # Arm A: the production assembler, tech explicitly :nanopore.
        res_with = run_arm(reads, :none, 21, length(refseq), with_dir, "naive";
            sequencing_tech = :nanopore)

        # Arm B: literally OMIT the kwarg, so the assembler takes its own
        # :illumina default — the pre-fix call shape.
        omitting = (r; kwargs...) -> begin
            kw = Dict{Symbol, Any}(kwargs)
            delete!(kw, :sequencing_tech)
            return Mycelia.Rhizomorph.assemble_genome(r; kw...)
        end
        res_without = run_arm(reads, :none, 21, length(refseq), without_dir, "naive";
            sequencing_tech = :nanopore, assembler = omitting)

        Test.@test res_with.ok
        Test.@test res_without.ok
        Test.@test res_with.n_contigs > 0

        # The decisive assertion: identical contig FASTA bytes.
        Test.@test read(res_with.contigs_path) == read(res_without.contigs_path)

        # And the derived metrics agree, so nothing downstream shifted either.
        Test.@test res_with.n_contigs == res_without.n_contigs
        Test.@test res_with.total_length == res_without.total_length
        Test.@test res_with.n50 == res_without.n50

        # The naive route does not stamp sequencing_tech, so the column records
        # the "n/a" sentinel on both — NOT the tech that was passed in.
        Test.@test res_with.corrector_sequencing_tech == "n/a"
        Test.@test res_without.corrector_sequencing_tech == "n/a"
        # Likewise the runtime-engagement counter: the naive route stamps no
        # corrector telemetry at all, so both arms report `missing`.
        Test.@test res_with.corrector_indel_engaged === missing
        Test.@test res_without.corrector_indel_engaged === missing
    end
end
