import Test
import Mycelia
import FASTX
import BioSequences
import StableRNGs

include(joinpath(dirname(@__DIR__), "..", "benchmarking", "rhizomorph_benchmark_suite_common.jl"))

Test.@testset "Rhizomorph Benchmark Suite Utilities" begin
    Test.@testset "Isolate matrix counts (tiered)" begin
        genomes = [
            (genome_id = "g1", genome_tier = "viral", genome_length = 5000, seed = 1)
        ]
        q_values = [10, 20]
        coverage_levels = [10, 25]
        anchor_q_values = [20]
        anchor_coverage_levels = [25]
        replicates = 2

        matrix = build_isolate_benchmark_matrix(
            genomes = genomes,
            q_values = q_values,
            coverage_levels = coverage_levels,
            replicates = replicates,
            strategy = :tiered,
            include_advanced = false,
            anchor_q_values = anchor_q_values,
            anchor_coverage_levels = anchor_coverage_levels
        )

        specs = isolate_algorithm_specs(include_advanced = false)
        primary_algorithms = count(spec -> !spec.anchor_only, specs)
        anchor_algorithms = count(spec -> spec.anchor_only, specs)
        expected_count = length(genomes) * replicates * (
            length(q_values) * length(coverage_levels) * primary_algorithms +
            length(anchor_q_values) * length(anchor_coverage_levels) * anchor_algorithms
        )
        Test.@test length(matrix) == expected_count
    end

    Test.@testset "Isolate matrix counts (full factorial)" begin
        genomes = [
            (genome_id = "g1", genome_tier = "viral", genome_length = 5000, seed = 1)
        ]
        q_values = [10, 20]
        coverage_levels = [10, 25]
        replicates = 1

        matrix = build_isolate_benchmark_matrix(
            genomes = genomes,
            q_values = q_values,
            coverage_levels = coverage_levels,
            replicates = replicates,
            strategy = :full_factorial,
            include_advanced = false
        )

        specs = isolate_algorithm_specs(include_advanced = false)
        expected_count = length(genomes) * replicates * length(q_values) *
                         length(coverage_levels) * length(specs)
        Test.@test length(matrix) == expected_count
    end

    Test.@testset "Deterministic isolate read simulation" begin
        reference = BioSequences.randdnaseq(StableRNGs.StableRNG(42), 8000)
        bundle_a = generate_isolate_read_bundle(
            reference_sequence = reference,
            coverage = 10,
            q_value = 30,
            seed = 12345
        )
        bundle_b = generate_isolate_read_bundle(
            reference_sequence = reference,
            coverage = 10,
            q_value = 30,
            seed = 12345
        )

        Test.@test bundle_a.n_reads == bundle_b.n_reads
        Test.@test length(bundle_a.fastq) == length(bundle_b.fastq)
        if !isempty(bundle_a.fastq)
            Test.@test String(FASTX.sequence(first(bundle_a.fastq))) ==
                       String(FASTX.sequence(first(bundle_b.fastq)))
        end

        capped_bundle = generate_isolate_read_bundle(
            reference_sequence = reference,
            coverage = 100,
            q_value = 30,
            seed = 12345,
            max_records = 10
        )
        Test.@test capped_bundle.n_reads == 10
        Test.@test length(capped_bundle.fastq) == 10
    end

    Test.@testset "Metagenome matrix quality profile coverage" begin
        complexity_levels = [10]
        depth_levels = [10, 25]
        abundance_profiles = [:equal]
        quality_profiles_input = [:const_q10, :const_q20, :mixed_balanced]
        anchor_depth_levels = [25]
        anchor_quality_profiles = [:mixed_balanced]

        matrix = build_metagenome_benchmark_matrix(
            complexity_levels = [10],
            depth_levels = depth_levels,
            abundance_profiles = abundance_profiles,
            quality_profiles = quality_profiles_input,
            replicates = 1,
            strategy = :tiered,
            include_advanced = false,
            anchor_depth_levels = anchor_depth_levels,
            anchor_quality_profiles = anchor_quality_profiles
        )
        quality_profiles = unique(row.quality_profile for row in matrix)
        Test.@test :const_q10 in quality_profiles
        Test.@test :const_q20 in quality_profiles
        Test.@test :mixed_balanced in quality_profiles

        specs = metagenome_algorithm_specs(include_advanced = false)
        primary_algorithms = count(spec -> !spec.anchor_only, specs)
        anchor_algorithms = count(spec -> spec.anchor_only, specs)
        expected_count = length(complexity_levels) * length(abundance_profiles) * (
            primary_algorithms * length(depth_levels) * length(quality_profiles_input) +
            anchor_algorithms * length(anchor_depth_levels) * length(anchor_quality_profiles)
        )
        Test.@test length(matrix) == expected_count
    end

    Test.@testset "Assembly metric contract" begin
        truth = ["ACGTACGTACGT"]
        contigs = ["ACGTACGTACGT"]
        metrics = compute_assembly_metrics(truth, contigs; k = 3)

        Test.@test hasproperty(metrics, :num_contigs)
        Test.@test hasproperty(metrics, :total_contig_length)
        Test.@test hasproperty(metrics, :n50)
        Test.@test hasproperty(metrics, :l50)
        Test.@test hasproperty(metrics, :assembly_size_ratio)
        Test.@test hasproperty(metrics, :kmer_precision)
        Test.@test hasproperty(metrics, :kmer_recall)
        Test.@test hasproperty(metrics, :kmer_f1)
        Test.@test metrics.kmer_precision ≈ 1.0
        Test.@test metrics.kmer_recall ≈ 1.0
        Test.@test metrics.kmer_f1 ≈ 1.0
    end
end
