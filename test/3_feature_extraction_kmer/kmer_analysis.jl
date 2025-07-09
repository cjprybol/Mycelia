import Pkg
Pkg.activate("..")
using Test
import Mycelia
import FASTX
import Random
import ProgressMeter

# unique sequences
Random.seed!(42)
fasta_list = String[]
ProgressMeter.@showprogress for i in 1:100
    fasta_file = Random.randstring() * ".fna.gz"
    isfile(fasta_file) && rm(fasta_file)
    Mycelia.write_fasta(outfile = fasta_file, records = [Mycelia.random_fasta_record(moltype=:DNA, L=100_000)])
    push!(fasta_list, fasta_file)
end

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=1)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=3)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=5)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=7)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=9)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=11)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=1)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=3)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=5)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=7)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=9)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=11)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=13)

for f in fasta_list
    if isfile(f)
        rm(f)
    end
end



# equivalent sequences
Random.seed!(42)
fasta_list = String[]
records = [Mycelia.random_fasta_record(moltype=:DNA, L=100_000)]
ProgressMeter.@showprogress for i in 1:100
    fasta_file = Random.randstring() * ".fna.gz"
    isfile(fasta_file) && rm(fasta_file)
    Mycelia.write_fasta(outfile = fasta_file, records = records)
    push!(fasta_list, fasta_file)
end

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=1)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=3)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=5)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=7)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=9)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=11)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=1)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=3)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=5)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=7)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=9)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=11)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:DNA, k=13)

for f in fasta_list
    if isfile(f)
        rm(f)
    end
end







# unique sequences
Random.seed!(42)
fasta_list = String[]
ProgressMeter.@showprogress for i in 1:100
    fasta_file = Random.randstring() * ".fna.gz"
    isfile(fasta_file) && rm(fasta_file)
    Mycelia.write_fasta(outfile = fasta_file, records = [Mycelia.random_fasta_record(moltype=:AA, L=100_000)])
    push!(fasta_list, fasta_file)
end

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=1)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=3)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=5)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=1)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=3)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=5)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=7)

for f in fasta_list
    if isfile(f)
        rm(f)
    end
end



# equivalent sequences
Random.seed!(42)
fasta_list = String[]
records = [Mycelia.random_fasta_record(moltype=:AA, L=100_000)]
ProgressMeter.@showprogress for i in 1:100
    fasta_file = Random.randstring() * ".fna.gz"
    isfile(fasta_file) && rm(fasta_file)
    Mycelia.write_fasta(outfile = fasta_file, records = records)
    push!(fasta_list, fasta_file)
end

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=1)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=3)

Mycelia.fasta_list_to_dense_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=5)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=1)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=3)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=5)

Mycelia.fasta_list_to_sparse_kmer_counts(fasta_list = fasta_list, alphabet=:AA, k=7)

for f in fasta_list
    if isfile(f)
        rm(f)
    end
end


# Reference Graph and K-mer Analysis tests
@testset "Reference Graph and K-mer Analysis" begin
    @testset "Pangenome Construction" begin
        @test true  # placeholder
    end
    @testset "Optimal K-mer Selection" begin
        @test true  # placeholder
    end
    @testset "statistical kmer analyses" begin
        @test Mycelia.optimal_subsequence_length(error_rate=0.001, sequence_length=100, threshold=0.99) == 10
    end
end
