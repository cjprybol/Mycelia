import Test
import Mycelia
import DataFrames
import StableRNGs

"""
CoverM integration test on a real (or fallback simulated) genome.
Runs only when `MYCELIA_RUN_COVERM=true` or `MYCELIA_RUN_EXTENDED=true`.
Uses: ART read simulation -> minimap2 mapping -> CoverM contig/genome.
"""
Test.@testset "CoverM integration (extended)" begin
    should_run = get(ENV, "MYCELIA_RUN_COVERM", "false") == "true" ||
                 get(ENV, "MYCELIA_RUN_EXTENDED", "false") == "true"

    if !should_run
        Test.@test_skip "CoverM integration skipped (set MYCELIA_RUN_COVERM=true or MYCELIA_RUN_EXTENDED=true)"
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
                    len=150,
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
                bam_path = Mycelia.sort_bam(map_result.outfile; threads=2)

                # Contig mode
                contig_df = Mycelia.run_coverm_contig(
                    bam_files=[bam_path],
                    reference_fasta=assembly,
                    methods=["mean", "covered_fraction"],
                    threads=2,
                    outdir=joinpath(tmp, "coverm_contig"),
                    quiet=true
                )
                Test.@test DataFrames.nrow(contig_df) >= 1
                Test.@test all(contig_df.mean .>= 0)
                if "covered_fraction" in Symbol.(names(contig_df))
                    Test.@test all((0 .<= contig_df.covered_fraction) .<= 1)
                end

                # Genome mode (single FASTA treated as one bin)
                genome_df = Mycelia.run_coverm_genome(
                    bam_files=[bam_path],
                    genome_fasta_files=[assembly],
                    methods=["relative_abundance", "mean_coverage"],
                    threads=2,
                    outdir=joinpath(tmp, "coverm_genome"),
                    quiet=true
                )
                Test.@test DataFrames.nrow(genome_df) == 1
                if "relative_abundance" in Symbol.(names(genome_df))
                    Test.@test genome_df.relative_abundance[1] > 0.0
                end
                if "mean_coverage" in Symbol.(names(genome_df))
                    Test.@test genome_df.mean_coverage[1] > 0.0
                end

            finally
                genome_info.cleanup()
            end
        end
    end
end
