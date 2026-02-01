# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/4_assembly/hpc_assembly_wrappers_test.jl")'
# ```

import Test
import Mycelia
import StableRNGs

const RUN_ALL = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
const RUN_EXTERNAL = RUN_ALL || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

# Helper to materialize gzipped FASTQ into a plain FASTQ file for tools that need it.
function materialize_fastq(input_fastq::AbstractString, output_fastq::AbstractString)
    reader = Mycelia.open_fastx(input_fastq)
    records = collect(reader)
    close(reader)
    Mycelia.write_fastq(records = records, filename = output_fastq)
    return output_fastq
end

Test.@testset "HPC Assembly Wrappers" begin
    if RUN_EXTERNAL
        mktempdir() do dir
            rng = StableRNGs.StableRNG(2025)
            record = Mycelia.random_fasta_record(moltype = :DNA, seed = rand(rng, 1:typemax(Int)), L = 5000)
            reference = joinpath(dir, "reference.fasta")
            Mycelia.write_fasta(outfile = reference, records = [record])

            illumina = Mycelia.simulate_illumina_reads(
                fasta = reference,
                coverage = 20,
                read_length = 100,
                rndSeed = 1234,
                quiet = true
            )
            pacbio_reads = Mycelia.simulate_pacbio_reads(
                fasta = reference,
                quantity = "10x",
                quiet = true,
                seed = 5678
            )

            short_1 = materialize_fastq(illumina.forward_reads, joinpath(dir, "short_1.fastq"))
            short_2 = materialize_fastq(illumina.reverse_reads, joinpath(dir, "short_2.fastq"))
            long_reads = materialize_fastq(pacbio_reads, joinpath(dir, "long.fastq"))

            Test.@testset "metaVelvet assembly" begin
                outdir = joinpath(dir, "metavelvet_assembly")
                try
                    result = Mycelia.run_metavelvet(
                        short_1; fastq2 = short_2, outdir = outdir,
                        k = 21, min_contig_lgth = 100)
                    Test.@test result.outdir == outdir
                    Test.@test isfile(result.contigs)
                catch e
                    @error "metaVelvet test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            end

            Test.@testset "Unicycler hybrid assembly" begin
                outdir = joinpath(dir, "unicycler_assembly")
                try
                    result = Mycelia.run_unicycler(short_1 = short_1, short_2 = short_2,
                        long_reads = long_reads, outdir = outdir, threads = 2,
                        spades_options = "-m 4", kmers = "21,33,55")
                    Test.@test result.outdir == outdir
                    Test.@test isfile(result.assembly)
                catch e
                    @error "Unicycler test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            end

            Test.@testset "Apollo polishing" begin
                outdir = joinpath(dir, "apollo_polish")
                try
                    result = Mycelia.run_apollo(reference, long_reads; outdir = outdir)
                    Test.@test result.outdir == outdir
                    Test.@test isfile(result.polished_assembly)
                catch e
                    @error "Apollo test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            end

            Test.@testset "Homopolish polishing" begin
                outdir = joinpath(dir, "homopolish_polish")
                try
                    result = Mycelia.run_homopolish(reference, reference; outdir = outdir, threads = 2)
                    Test.@test result.outdir == outdir
                    Test.@test isfile(result.polished_assembly)
                catch e
                    @error "Homopolish test failed." exception=(e, catch_backtrace())
                    Test.@test false
                end
            end
        end
    else
        @info "HPC assembly wrappers disabled; set MYCELIA_RUN_EXTERNAL=true to enable."
    end
end
