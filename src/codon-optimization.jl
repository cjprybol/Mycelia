function fasta_to_kmer_and_codon_frequencies(fasta)
    kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{5}, fasta)
    pyrodigal_results = Mycelia.run_pyrodigal(fasta_file = fasta)
    if occursin(r"\.gz$", fasta)
        unzipped_fasta = replace(fasta, ".gz" => "")
        run(pipeline(`gzip -dc $(fasta)`, unzipped_fasta))
        @assert isfile(unzipped_fasta)
    end
    genbank = Mycelia.fasta_and_gff_to_genbank(fasta = unzipped_fasta, gff = pyrodigal_results.gff)
    codon_frequencies = Mycelia.genbank_to_codon_frequencies(genbank)
    return (;kmer_counts, codon_frequencies)
end

function calculate_sequence_likelihood_from_kmer_profile(sequence, normalized_kmer_counts)
    initial_likelihood = 1.0
    for (k, i) in Kmers.UnambiguousDNAMers{5}(sequence)
        canonical_k = BioSequences.canonical(k)
        # # display(k)
        # try
        #     @assert BioSequences.iscanonical(k)
        # catch
        #     display(k)
        # end
        initial_likelihood *= normalized_kmer_counts[canonical_k]
        # end
    end
    initial_likelihood
end