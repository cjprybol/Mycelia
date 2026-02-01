# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/4_assembly/third_party_assemblers_hybrid.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_hybrid.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import StableRNGs
import BioSequences
import FASTX

Test.@testset "Hybrid Assembly Tools" begin
    Test.@testset "Input validation" begin
        outdir = mktempdir()
        Test.@test_throws ErrorException Mycelia.run_hybracter_hybrid_single(
            "missing_long.fq",
            "missing_r1.fq",
            "missing_r2.fq",
            5_000_000;
            sample_name = "sample",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_hybracter_long_single(
            "missing_long.fq",
            5_000_000;
            sample_name = "sample",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_hybracter_hybrid_single(
            "missing_long.fq",
            "missing_r1.fq",
            "missing_r2.fq",
            5_000_000;
            sample_name = "",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_hybracter_long_single(
            "missing_long.fq",
            5_000_000;
            sample_name = "",
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_plassembler(
            "missing_long.fq",
            "missing_r1.fq",
            "missing_r2.fq",
            5_000_000;
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_plassembler_long(
            "missing_long.fq",
            5_000_000;
            outdir = outdir
        )
        Test.@test_throws ErrorException Mycelia.run_dnaapler_all(
            "missing_input.fasta";
            outdir = outdir
        )
    end

    run_external = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true" ||
                   lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"

    if run_external
        threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

        Test.@testset "DNAapler (external)" begin
            mktempdir() do dir
                input_fasta = joinpath(dir, "contigs.fasta")
                open(input_fasta, "w") do io
                    write(io, ">contig1\nACGTACGTACGTACGT\n")
                end

                outdir = joinpath(dir, "dnaapler_out")
                try
                    result = Mycelia.run_dnaapler_all(
                        input_fasta; outdir = outdir, force = true, threads = threads)
                    Test.@test isfile(result.reoriented)
                catch e
                    @error "DNAapler test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            end
        end

        Test.@testset "Hybracter/Plassembler (external)" begin
            mktempdir() do dir
                ref_fasta = joinpath(dir, "hybrid_ref.fasta")
                rng = StableRNGs.StableRNG(101)
                genome = BioSequences.randdnaseq(rng, 30_000)
                Mycelia.write_fasta(outfile = ref_fasta,
                    records = [FASTX.FASTA.Record("hybrid_ref", genome)])

                chrom_size = length(genome)
                plassembler_chrom_size = max(10_000, Int(floor(chrom_size * 0.8)))

                illumina = Mycelia.simulate_illumina_reads(
                    fasta = ref_fasta,
                    coverage = 10,
                    outbase = joinpath(dir, "hybrid_short"),
                    read_length = 150,
                    mflen = 300,
                    seqSys = "HS25",
                    paired = true,
                    errfree = true,
                    rndSeed = 101,
                    quiet = true
                )

                long_reads_gz = Mycelia.simulate_nanopore_reads(
                    fasta = ref_fasta,
                    quantity = "12x",
                    quiet = true,
                    seed = 101
                )
                long_reads = joinpath(dir, "hybrid_long.fq")
                run(pipeline(`gunzip -c $(long_reads_gz)`, long_reads))

                Test.@testset "Hybracter hybrid-single" begin
                    outdir = joinpath(dir, "hybracter_out")
                    try
                        result = Mycelia.run_hybracter_hybrid_single(
                            long_reads,
                            illumina.forward_reads,
                            illumina.reverse_reads,
                            chrom_size;
                            sample_name = "hybracter_sample",
                            outdir = outdir,
                            threads = threads
                        )
                        Test.@test isfile(result.final_fasta)
                    catch e
                        @error "Hybracter test failed." exception=(e, catch_backtrace())
                        Test.@test false
                    finally
                        isdir(outdir) && rm(outdir; recursive = true, force = true)
                    end
                end

                Test.@testset "Plassembler run" begin
                    outdir = joinpath(dir, "plassembler_out")
                    try
                        result = Mycelia.run_plassembler(
                            long_reads,
                            illumina.forward_reads,
                            illumina.reverse_reads,
                            plassembler_chrom_size;
                            outdir = outdir,
                            threads = threads
                        )
                        Test.@test isfile(result.plasmids_fasta)
                    catch e
                        @error "Plassembler test failed." exception=(e, catch_backtrace())
                        Test.@test false
                    finally
                        isdir(outdir) && rm(outdir; recursive = true, force = true)
                    end
                end
            end
        end
    else
        @info "External tool runs disabled; set MYCELIA_RUN_EXTERNAL=true to enable."
    end
end

Test.@testset "Hybrid Metagenomic Assembly - HyLight" begin
    mktempdir() do outdir
        Test.@test_throws ErrorException Mycelia.run_hylight(
            "missing_r1.fq",
            "missing_r2.fq",
            "missing_long.fq";
            outdir = outdir
        )
    end
end
