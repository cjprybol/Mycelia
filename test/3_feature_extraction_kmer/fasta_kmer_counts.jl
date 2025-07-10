using Test
import Mycelia
import FASTX

@testset "fasta_list_to_*_kmer_counts" begin
    sequences = ["ACGTACGT", "AAAACCCC", "GGGGAAAA"]
    fasta_files = String[]
    for (i, seq) in enumerate(sequences)
        fname = "test_seq$(i).fasta"
        open(fname, "w") do io
            FASTX.write(io, FASTX.FASTA.Record("seq$(i)", seq))
        end
        push!(fasta_files, fname)
    end

    k = 3
    dense_res = Mycelia.fasta_list_to_dense_kmer_counts(
        fasta_list=fasta_files,
        k=k,
        alphabet=:DNA,
    )
    sparse_res = Mycelia.fasta_list_to_sparse_kmer_counts(
        fasta_list=fasta_files,
        k=k,
        alphabet=:DNA,
    )

    expected_rows = length(Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET))
    nfiles = length(fasta_files)
    @test size(dense_res.counts, 1) == expected_rows
    @test size(dense_res.counts, 2) == nfiles
    @test size(sparse_res.counts, 1) == expected_rows
    @test size(sparse_res.counts, 2) == nfiles

    for (i, seq) in enumerate(sequences)
        expected_total = length(seq) - k + 1
        @test sum(dense_res.counts[:, i]) == expected_total
        @test sum(sparse_res.counts[:, i]) == expected_total
    end

    for f in fasta_files
        isfile(f) && rm(f)
    end
end
