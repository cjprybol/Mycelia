# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/5_validation/coverm_integration_extended.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/coverm_integration_extended.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import DataFrames
import StableRNGs

# CoverM integration test on a real (or fallback simulated) genome.
# Runs only when `MYCELIA_RUN_EXTERNAL=true` or `MYCELIA_RUN_ALL=true`.
# Uses: ART read simulation -> minimap2 mapping -> CoverM contig/genome.
Test.@testset "CoverM integration (extended)" begin
    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if !run_external
        Test.@test_skip "CoverM integration skipped (set MYCELIA_RUN_EXTERNAL=true to enable)"
    else
        mktempdir() do tmp
            rng = StableRNGs.StableRNG(7)

            # Fetch a small NCBI genome (fallbacks to simulated if download fails)
            genome_info = Mycelia.get_test_genome_fasta(use_ncbi=true)
            assembly = genome_info.fasta

            try
                # Simulate reads and map to BAM
                reads = Mycelia.simulate_illumina_reads(
                    fasta=assembly,
                    read_count=1_000,
                    read_length=150,
                    rndSeed=rand(rng, 0:typemax(Int)),
                    quiet=true
                )
                map_result = Mycelia.minimap_map(
                    fasta=assembly,
                    fastq=reads.forward_reads,
                    mapping_type="sr",
                    threads=2
                )
                run(map_result.cmd)
                bam_path, _ = Mycelia.index_bam(map_result.outfile; threads=2)

                # Contig mode
                contig_df = Mycelia.run_coverm_contig(
                    bam_files=[bam_path],
                    methods=["mean", "covered_fraction"],
                    threads=2,
                    outdir=joinpath(tmp, "coverm_contig"),
                    quiet=true
                )
                Test.@test DataFrames.nrow(contig_df) >= 1
                mean_cols = filter(name -> endswith(String(name), "_Mean"), names(contig_df))
                Test.@test !isempty(mean_cols)
                for col in mean_cols
                    Test.@test all(contig_df[!, col] .>= 0)
                end
                covered_cols = filter(name -> endswith(String(name), "_Covered_Fraction"), names(contig_df))
                for col in covered_cols
                    Test.@test all((0 .<= contig_df[!, col]) .& (contig_df[!, col] .<= 1))
                end

                # Genome mode (single FASTA treated as one bin)
                genome_df = Mycelia.run_coverm_genome(
                    bam_files=[bam_path],
                    genome_fasta_files=[assembly],
                    methods=["relative_abundance", "mean"],
                    threads=2,
                    outdir=joinpath(tmp, "coverm_genome"),
                    quiet=true
                )
                Test.@test DataFrames.nrow(genome_df) >= 1
                mapped_df = genome_df[genome_df.Genome .!= "unmapped", :]
                Test.@test DataFrames.nrow(mapped_df) >= 1
                abundance_cols = filter(name -> occursin("_Relative_Abundance", String(name)), names(genome_df))
                for col in abundance_cols
                    Test.@test all(mapped_df[!, col] .> 0)
                end
                mean_cols = filter(name -> endswith(String(name), "_Mean"), names(genome_df))
                for col in mean_cols
                    values = mapped_df[!, col]
                    if eltype(values) <: AbstractString
                        parsed = parse.(Float64, values)
                        Test.@test all(parsed .> 0)
                    else
                        Test.@test all(values .> 0)
                    end
                end

            finally
                genome_info.cleanup()
            end
        end
    end
end
