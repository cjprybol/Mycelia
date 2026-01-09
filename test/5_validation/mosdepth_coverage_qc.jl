# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/5_validation/mosdepth_coverage_qc.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/mosdepth_coverage_qc.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import FASTX
import DataFrames
import StableRNGs

Test.@testset "Mosdepth Coverage QC" begin
    ## Setup: Create test directory
    test_dir = mktempdir()

    Test.@testset "Multi-Contig BAM Sorting and Coverage Analysis" begin
        ## Step 1: Create multi-contig FASTA reference
        rng = StableRNGs.StableRNG(42)

        ## Generate 3 different contigs/chromosomes with different lengths
        contig1_record = Mycelia.random_fasta_record(L=5000, seed=rand(rng, 0:typemax(Int)))
        contig2_record = Mycelia.random_fasta_record(L=3000, seed=rand(rng, 0:typemax(Int)))
        contig3_record = Mycelia.random_fasta_record(L=4000, seed=rand(rng, 0:typemax(Int)))

        ## Write multi-contig FASTA file
        multi_contig_fasta = joinpath(test_dir, "reference.fasta")
        Mycelia.write_fasta(outfile=multi_contig_fasta, records=[contig1_record, contig2_record, contig3_record])

        Test.@test isfile(multi_contig_fasta)

        ## Step 2: Simulate reads from multi-contig reference
        ## Use lower coverage for speed but enough to get meaningful results
        simulated_reads = Mycelia.simulate_illumina_reads(fasta=multi_contig_fasta, coverage=10, quiet=true)

        Test.@test isfile(simulated_reads.forward_reads)

        ## Step 3: Map reads to create unsorted BAM
        map_result = Mycelia.minimap_map(
            fasta=multi_contig_fasta,
            fastq=simulated_reads.forward_reads,
            mapping_type="sr",
            output_format="bam",
            sorted=false,
            quiet=true
        )
        run(map_result.cmd)
        unsorted_bam = map_result.outfile

        Test.@test isfile(unsorted_bam)
        Test.@test !occursin("sorted", unsorted_bam)

        ## Step 4: Test is_bam_coordinate_sorted on unsorted BAM
        Test.@testset "Unsorted BAM Detection" begin
            is_sorted = Mycelia.is_bam_coordinate_sorted(unsorted_bam)
            Test.@test is_sorted == false
        end

        ## Step 5: Sort the BAM file
        sorted_bam = Mycelia.sort_bam(unsorted_bam)

        Test.@test isfile(sorted_bam)
        Test.@test occursin("sorted", sorted_bam)

        ## Step 6: Test is_bam_coordinate_sorted on sorted BAM
        Test.@testset "Sorted BAM Detection" begin
            is_sorted = Mycelia.is_bam_coordinate_sorted(sorted_bam)
            Test.@test is_sorted == true
        end

        ## Step 7: Index the sorted BAM
        indexed_bam, bai_path = Mycelia.index_bam(sorted_bam)

        Test.@test isfile(indexed_bam)
        Test.@test isfile(bai_path)
        Test.@test indexed_bam == sorted_bam
        Test.@test bai_path == sorted_bam * ".bai"

        ## Step 8: Run mosdepth on the indexed BAM
        mosdepth_output = Mycelia.run_mosdepth(sorted_bam, no_per_base=true)

        Test.@test isfile(mosdepth_output.global_dist)
        Test.@test isfile(mosdepth_output.summary)

        ## Step 9: Parse mosdepth distribution file
        Test.@testset "Parse Mosdepth Distribution" begin
            dist_df = Mycelia.parse_mosdepth_distribution(mosdepth_output.global_dist)

            Test.@test dist_df isa DataFrames.DataFrame
            Test.@test "chromosome" in names(dist_df)
            Test.@test "coverage" in names(dist_df)
            Test.@test "proportion" in names(dist_df)

            ## Should have data for each contig plus "total"
            unique_chroms = unique(dist_df.chromosome)
            Test.@test "total" in unique_chroms
            Test.@test length(unique_chroms) >= 2  ## At least total + one contig
        end

        ## Step 10: Summarize QC metrics from distribution
        Test.@testset "Summarize Mosdepth QC Metrics" begin
            dist_df = Mycelia.parse_mosdepth_distribution(mosdepth_output.global_dist)
            qc_metrics = Mycelia.summarize_mosdepth_qc(dist_df)

            Test.@test qc_metrics isa DataFrames.DataFrame
            Test.@test "chromosome" in names(qc_metrics)

            ## Check for expected coverage threshold columns
            expected_thresholds = [1, 3, 5, 10, 30, 50, 100, 300, 500, 1000]
            for threshold in expected_thresholds
                col_name = Symbol("coverage_$(threshold)X")
                Test.@test String(col_name) in names(qc_metrics)
            end

            ## Should have metrics for each contig plus total
            unique_chroms = unique(qc_metrics.chromosome)
            Test.@test "total" in unique_chroms

            ## Coverage proportions should be between 0 and 1 (skip missing values)
            for threshold in expected_thresholds
                col_name = Symbol("coverage_$(threshold)X")
                Test.@test all(0 .<= x .<= 1 for x in skipmissing(qc_metrics[!, col_name]))
            end

            ## Higher coverage thresholds should have equal or lower proportions
            ## (proportions should increase as coverage threshold decreases)
            ## Check monotonicity within each chromosome separately
            for chrom_df in DataFrames.groupby(qc_metrics, :chromosome)
                for j in 1:(length(expected_thresholds)-1)
                    lower_threshold = expected_thresholds[j]
                    higher_threshold = expected_thresholds[j+1]
                    lower_col = Symbol("coverage_$(lower_threshold)X")
                    higher_col = Symbol("coverage_$(higher_threshold)X")
                    lower_val = chrom_df[1, lower_col]
                    higher_val = chrom_df[1, higher_col]
                    ## Proportion at lower threshold should be >= proportion at higher threshold
                    ## Skip if either value is missing (threshold not in raw data)
                    if !ismissing(lower_val) && !ismissing(higher_val)
                        Test.@test lower_val >= higher_val
                    end
                end
            end
        end

        ## Step 11: Parse mosdepth summary statistics
        Test.@testset "Parse Mosdepth Summary" begin
            summary_df = Mycelia.parse_mosdepth_summary(mosdepth_output.summary)

            Test.@test summary_df isa DataFrames.DataFrame
            Test.@test "chrom" in names(summary_df) || "chromosome" in names(summary_df)
            Test.@test "length" in names(summary_df)
            Test.@test "bases" in names(summary_df)
            Test.@test "mean" in names(summary_df)
            Test.@test "min" in names(summary_df)
            Test.@test "max" in names(summary_df)

            ## Should have summary for each contig plus total
            chrom_col = "chrom" in names(summary_df) ? "chrom" : "chromosome"
            unique_chroms = unique(summary_df[!, chrom_col])
            Test.@test "total" in unique_chroms
            Test.@test length(unique_chroms) >= 2  ## At least total + one contig

            ## Mean coverage should be positive for all contigs
            Test.@test all(summary_df[!, "mean"] .>= 0)
        end

        ## Step 12: Merge QC metrics and summary tables
        Test.@testset "Merge QC and Summary Tables" begin
            dist_df = Mycelia.parse_mosdepth_distribution(mosdepth_output.global_dist)
            qc_metrics = Mycelia.summarize_mosdepth_qc(dist_df)
            summary_df = Mycelia.parse_mosdepth_summary(mosdepth_output.summary)

            ## Standardize chromosome column names for merging
            chrom_col_summary = "chrom" in names(summary_df) ? "chrom" : "chromosome"
            if chrom_col_summary != "chromosome"
                DataFrames.rename!(summary_df, chrom_col_summary => "chromosome")
            end

            ## Merge QC metrics with summary statistics
            comprehensive_table = DataFrames.innerjoin(qc_metrics, summary_df, on="chromosome")

            Test.@test comprehensive_table isa DataFrames.DataFrame
            Test.@test DataFrames.nrow(comprehensive_table) >= 2  ## At least total + one contig

            ## Check that merged table has columns from both sources
            Test.@test "chromosome" in names(comprehensive_table)
            Test.@test "coverage_1X" in names(comprehensive_table)
            Test.@test "coverage_10X" in names(comprehensive_table)
            Test.@test "coverage_30X" in names(comprehensive_table)
            Test.@test "length" in names(comprehensive_table)
            Test.@test "mean" in names(comprehensive_table)

            ## Verify we have data for all expected chromosomes
            unique_chroms = unique(comprehensive_table.chromosome)
            Test.@test "total" in unique_chroms

            ## For debugging/validation, show the comprehensive table structure
            ## (in actual tests this would be commented out)
            @info "Comprehensive coverage table columns:" names(comprehensive_table)
            @info "Number of chromosomes/contigs:" DataFrames.nrow(comprehensive_table)
        end

        ## Cleanup temporary files
        rm(multi_contig_fasta, force=true)
        rm(simulated_reads.forward_reads, force=true)
        rm(simulated_reads.reverse_reads, force=true)
        if simulated_reads.sam !== nothing
            rm(simulated_reads.sam, force=true)
        end
        if simulated_reads.error_free_sam !== nothing
            rm(simulated_reads.error_free_sam, force=true)
        end
        rm(unsorted_bam, force=true)
        rm(sorted_bam, force=true)
        rm(bai_path, force=true)
        rm(mosdepth_output.global_dist, force=true)
        rm(mosdepth_output.summary, force=true)
        ## Also clean up mosdepth index files if they exist
        for ext in [".mosdepth.region.dist.txt", ".mosdepth.summary.txt.csi"]
            cleanup_file = sorted_bam * ext
            if isfile(cleanup_file)
                rm(cleanup_file, force=true)
            end
        end
    end

    ## Cleanup test directory
    rm(test_dir, recursive=true, force=true)
end
