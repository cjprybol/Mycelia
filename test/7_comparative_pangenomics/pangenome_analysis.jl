# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/pangenome_analysis.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/pangenome_analysis.jl", "test/7_comparative_pangenomics", execute=false)'
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
import Kmers
import BioSequences
import Statistics

Test.@testset "Pangenome analysis - k-mer workflows" begin
    temp_dir = mktempdir()
    try
        genomes = [
            ("genome1.fasta", BioSequences.LongDNA{4}("ATGCGATGCA")),
            ("genome2.fasta", BioSequences.LongDNA{4}("ATGCGTTTAA")),
            ("genome3.fasta", BioSequences.LongDNA{4}("TTTATGCGAA"))
        ]

        genome_files = String[]
        for (name, seq) in genomes
            path = joinpath(temp_dir, name)
            open(path, "w") do io
                FASTX.write(io, FASTX.FASTA.Record(name, seq))
            end
            push!(genome_files, path)
        end

        kmer_type = Kmers.DNAKmer{3}
        result = Mycelia.analyze_pangenome_kmers(
            genome_files;
            kmer_type = kmer_type,
            distance_metric = :jaccard
        )

        genome_names = [basename(path) for path in genome_files]
        Test.@test result.genome_names == genome_names

        expected_counts = Dict{String, Dict}()
        for (name, path) in zip(genome_names, genome_files)
            expected_counts[name] = Mycelia.count_canonical_kmers(kmer_type, path)
            Test.@test result.kmer_counts_by_genome[name] == expected_counts[name]
        end

        kmer_sets = [Set(keys(expected_counts[name])) for name in genome_names]
        all_kmers = reduce(union, kmer_sets)
        expected_core = reduce(intersect, kmer_sets)

        presence_counts = Dict{Any, Int}()
        for kmer_set in kmer_sets
            for kmer in kmer_set
                presence_counts[kmer] = get(presence_counts, kmer, 0) + 1
            end
        end

        expected_accessory = Set(
            kmer
        for (kmer, count) in presence_counts
        if count > 1 && count < length(kmer_sets)
        )
        expected_shared = union(expected_core, expected_accessory)

        expected_unique = Dict{String, Set}()
        for (name, kmer_set) in zip(genome_names, kmer_sets)
            other_kmers = Set{Any}()
            for (other_name, other_set) in zip(genome_names, kmer_sets)
                if other_name != name
                    union!(other_kmers, other_set)
                end
            end
            expected_unique[name] = setdiff(kmer_set, other_kmers)
        end

        Test.@test Set(result.core_kmers) == expected_core
        Test.@test Set(result.accessory_kmers) == expected_accessory
        Test.@test Set(result.shared_kmers) == expected_shared

        for name in genome_names
            Test.@test Set(result.unique_kmers_by_genome[name]) == expected_unique[name]
        end

        Test.@test size(result.presence_absence_matrix, 1) == length(all_kmers)
        Test.@test size(result.presence_absence_matrix, 2) == length(genome_names)
        Test.@test result.distance_matrix ==
                   Mycelia.jaccard_distance(result.presence_absence_matrix)

        expected_mean = Statistics.mean(result.distance_matrix[result.distance_matrix .> 0])
        Test.@test result.similarity_stats.core_size == length(expected_core)
        Test.@test result.similarity_stats.accessory_size == length(expected_accessory)
        Test.@test result.similarity_stats.pangenome_size == length(all_kmers)
        Test.@test result.similarity_stats.unique_total ==
                   sum(length(v) for v in values(expected_unique))
        Test.@test isapprox(result.similarity_stats.mean_pairwise_distance, expected_mean; atol = 1e-12)

        bray_result = Mycelia.analyze_pangenome_kmers(
            genome_files;
            kmer_type = kmer_type,
            distance_metric = :bray_curtis
        )
        Test.@test size(bray_result.distance_matrix) ==
                   (length(genome_files), length(genome_files))
        Test.@test bray_result.distance_matrix == transpose(bray_result.distance_matrix)
        Test.@test all(bray_result.distance_matrix[i, i] == 0
        for i in 1:length(genome_files))

        jaccard_similarity = Mycelia.compare_genome_kmer_similarity(
            genome_files[1],
            genome_files[2];
            kmer_type = kmer_type,
            metric = :jaccard
        )

        kmer_counts_1 = expected_counts[genome_names[1]]
        kmer_counts_2 = expected_counts[genome_names[2]]
        shared_kmers = length(intersect(keys(kmer_counts_1), keys(kmer_counts_2)))
        total_kmers = length(union(keys(kmer_counts_1), keys(kmer_counts_2)))
        expected_distance = 1.0 - shared_kmers / total_kmers

        Test.@test isapprox(jaccard_similarity.distance, expected_distance; atol = 1e-12)
        Test.@test isapprox(
            jaccard_similarity.jaccard_similarity, shared_kmers /
                                                   total_kmers; atol = 1e-12)
        Test.@test jaccard_similarity.total_kmers == total_kmers

        js_divergence = Mycelia.compare_genome_kmer_similarity(
            genome_files[1],
            genome_files[2];
            kmer_type = kmer_type,
            metric = :js_divergence
        )
        js_divergence_same = Mycelia.compare_genome_kmer_similarity(
            genome_files[1],
            genome_files[1];
            kmer_type = kmer_type,
            metric = :js_divergence
        )
        Test.@test js_divergence.metric == :js_divergence
        Test.@test js_divergence.distance > 0
        Test.@test isapprox(js_divergence_same.distance, 0.0; atol = 1e-12)

        cosine_similarity = Mycelia.compare_genome_kmer_similarity(
            genome_files[1],
            genome_files[2];
            kmer_type = kmer_type,
            metric = :cosine
        )
        cosine_similarity_same = Mycelia.compare_genome_kmer_similarity(
            genome_files[1],
            genome_files[1];
            kmer_type = kmer_type,
            metric = :cosine
        )
        Test.@test cosine_similarity.metric == :cosine
        Test.@test cosine_similarity.distance > 0
        Test.@test isapprox(cosine_similarity_same.distance, 0.0; atol = 1e-12)

        distance_result = Mycelia.build_genome_distance_matrix(
            genome_files;
            kmer_type = kmer_type,
            metric = :jaccard
        )
        Test.@test distance_result.genome_names == genome_names
        Test.@test size(distance_result.distance_matrix) ==
                   (length(genome_files), length(genome_files))
        Test.@test all(distance_result.distance_matrix[i, i] == 0
        for i in 1:length(genome_files))
        Test.@test distance_result.distance_matrix[1, 2] == jaccard_similarity.distance
        Test.@test distance_result.distance_matrix[2, 1] == jaccard_similarity.distance
    finally
        isdir(temp_dir) && rm(temp_dir, recursive = true, force = true)
    end
end

Test.@testset "Pangenome analysis - input validation" begin
    Test.@test_throws ErrorException Mycelia.analyze_pangenome_kmers(String[])

    temp_dir = mktempdir()
    try
        missing_file = joinpath(temp_dir, "missing.fasta")
        Test.@test_throws ErrorException Mycelia.analyze_pangenome_kmers([missing_file])

        valid_file = joinpath(temp_dir, "valid.fasta")
        open(valid_file, "w") do io
            FASTX.write(io, FASTX.FASTA.Record("seq1", BioSequences.LongDNA{4}("ATGCGT")))
        end
        Test.@test_throws ErrorException Mycelia.analyze_pangenome_kmers(
            [valid_file]; distance_metric = :cosine)

        other_file = joinpath(temp_dir, "other.fasta")
        open(other_file, "w") do io
            FASTX.write(io, FASTX.FASTA.Record("seq2", BioSequences.LongDNA{4}("ATGCAT")))
        end
        Test.@test_throws ErrorException Mycelia.compare_genome_kmer_similarity(
            valid_file,
            other_file;
            metric = :euclidean
        )
    finally
        isdir(temp_dir) && rm(temp_dir, recursive = true, force = true)
    end
end
