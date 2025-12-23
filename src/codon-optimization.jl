# function fasta_to_kmer_and_codon_frequencies(fasta)
#     kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{5}, fasta)
#     pyrodigal_results = Mycelia.run_pyrodigal(fasta_file = fasta)
#     if occursin(r"\.gz$", fasta)
#         unzipped_fasta = replace(fasta, ".gz" => "")
#         run(pipeline(`gzip -dc $(fasta)`, unzipped_fasta))
#         @assert isfile(unzipped_fasta)
#     end
#     genbank = Mycelia.fasta_and_gff_to_genbank(fasta = unzipped_fasta, gff = pyrodigal_results.gff)
#     codon_frequencies = Mycelia.genbank_to_codon_frequencies(genbank)
#     return (;kmer_counts, codon_frequencies)
# end

# function calculate_sequence_likelihood_from_kmer_profile(sequence, normalized_kmer_counts)
#     initial_likelihood = 1.0
#     for (k, i) in Kmers.UnambiguousDNAMers{5}(sequence)
#         canonical_k = BioSequences.canonical(k)
#         # # display(k)
#         # try
#         #     @assert BioSequences.iscanonical(k)
#         # catch
#         #     display(k)
#         # end
#         initial_likelihood *= normalized_kmer_counts[canonical_k]
#         # end
#     end
#     initial_likelihood
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a protein sequence back to a possible DNA coding sequence using weighted random codon selection.

# Arguments
- `protein_sequence::BioSequences.LongAA`: The amino acid sequence to reverse translate

# Returns
- `BioSequences.LongDNA{2}`: A DNA sequence that would translate to the input protein sequence

# Details
Uses codon usage frequencies to randomly select codons for each amino acid, weighted by their 
natural occurrence. Each selected codon is guaranteed to translate back to the original amino acid.
"""
function reverse_translate(protein_sequence::BioSequences.LongAA)
    this_sequence = BioSequences.LongDNA{2}()
    codon_frequencies = Dict(a => Dict{Kmers.DNACodon, Int}() for a in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term]))
    for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
        amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
        codon_frequencies[amino_acid][codon] = get(codon_frequencies[amino_acid], codon, 0) + 1
    end    
    for amino_acid in protein_sequence
        # I'm collecting first because I'm worried about the keys and values not being sorted the same between queries, but that feels like it's not a viable worry
        collected = collect(codon_frequencies[amino_acid])
        codons = first.(collected)
        frequencies = last.(collected)
        chosen_codon_index = StatsBase.sample(1:length(codons), StatsBase.weights(frequencies))
        chosen_codon = codons[chosen_codon_index]
        @assert first(BioSequences.translate(BioSequences.LongDNA{2}(chosen_codon))) == amino_acid
        this_sequence *= chosen_codon
    end
    @assert BioSequences.translate(this_sequence) == protein_sequence
    return this_sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Optimizes the DNA sequence encoding for a given protein sequence using codon usage frequencies.

# Arguments
- `normalized_codon_frequencies`: Dictionary mapping amino acids to their codon frequencies
- `protein_sequence::BioSequences.LongAA`: Target protein sequence to optimize
- `n_iterations::Integer`: Number of optimization iterations to perform

# Algorithm
1. Creates initial DNA sequence through reverse translation
2. Iteratively generates new sequences by sampling codons based on their frequencies
3. Keeps track of the sequence with highest codon usage likelihood

# Returns
- `BioSequences.LongDNA{2}`: Optimized DNA sequence encoding the input protein
"""
function codon_optimize(;normalized_codon_frequencies, protein_sequence::BioSequences.LongAA, n_iterations)
    best_sequence = reverse_translate(protein_sequence)
    # codons = last.(collect(Kmers.SpacedKmers{Kmers.DNACodon}(BioSequences.LongDNA{4}(best_sequence), 3)))
    codons = first.(collect(Kmers.UnambiguousDNAMers{3}(BioSequences.LongDNA{4}(best_sequence))))[1:3:end]
    initial_log_likelihood = -log10(1.0)
    for (codon, amino_acid) in collect(zip(codons, protein_sequence))
        this_codon_likelihood = normalized_codon_frequencies[amino_acid][codon]
        initial_log_likelihood -= log10(this_codon_likelihood)
    end
    best_likelihood = initial_log_likelihood

    ProgressMeter.@showprogress for i in 1:n_iterations
    # for iteration in 1:n_iterations
        this_sequence = BioSequences.LongDNA{2}()
        this_log_likelihood = -log10(1.0)
        for amino_acid in protein_sequence
            # I'm collecting first because I'm worried about the keys and values not being sorted the same between queries, but that feels like it's not a viable worry
            collected = collect(normalized_codon_frequencies[amino_acid])
            codons = first.(collected)
            frequencies = last.(collected)
            chosen_codon_index = StatsBase.sample(1:length(codons), StatsBase.weights(frequencies))
            chosen_codon = codons[chosen_codon_index]
            chosen_codon_frequency = frequencies[chosen_codon_index]
            this_log_likelihood += -log10(chosen_codon_frequency)
            this_sequence *= chosen_codon
        end
        if this_log_likelihood < best_likelihood
            best_likelihood = this_log_likelihood
            best_sequence = this_sequence
        end
        @assert BioSequences.translate(this_sequence) == protein_sequence
    end
    # @show (best_likelihood)^-10 / (initial_log_likelihood)^-10
    return best_sequence
end

# function codon_optimize(;normalized_codon_frequencies, optimization_sequence::BioSequences.LongDNA, n_iterations)
#     protein_sequence = BioSequences.translate(optimization_sequence)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a mapping from amino acids to representative DNA codons using the standard genetic code.

# Returns
- Dictionary mapping each amino acid (including stop codon `AA_Term`) to a valid DNA codon that encodes it
"""
function amino_acids_to_codons()
    amino_acid_to_codon_map = Dict{BioSymbols.AminoAcid, Kmers.DNACodon}()
    for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
        amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
        amino_acid_to_codon_map[amino_acid] = codon
    end   
    return amino_acid_to_codon_map
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Creates a mapping between DNA codons and their corresponding amino acids using the standard genetic code.

Returns a dictionary where:
- Keys are 3-letter DNA codons (e.g., "ATG")
- Values are the corresponding amino acids from BioSequences.jl
"""
function codons_to_amino_acids()
    codon_to_amino_acid_map = Dict{Kmers.DNACodon, BioSequences.LongAA}()
    for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
        codon_to_amino_acid_map[codon] = BioSequences.translate(BioSequences.LongDNA{2}(codon))
    end
    return codon_to_amino_acid_map
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Analyze codon usage frequencies from genes in a GenBank file.

# Arguments
- `genbank`: Path to GenBank format file containing genomic sequences and annotations
- `allow_all`: If true, initializes frequencies for all possible codons with count=1 (default: true)

# Returns
Nested dictionary mapping amino acids to their corresponding codon usage counts:
- Outer key: AminoAcid (including stop codon)
- Inner key: DNACodon
- Value: Count of codon occurrences

# Details
- Only processes genes marked as ':misc_feature' in the GenBank file
- Analyzes both forward and reverse complement sequences
- Determines coding strand based on presence of stop codons and start codons
- Skips ambiguous sequences that cannot be confidently oriented
"""
function genbank_to_codon_frequencies(genbank; allow_all=true)
    # create an initial codon frequency table, where we initialize all possible codons with equal probability
    # this way, if we don't see the amino acid in the observed proteins we're optimizing, we can still produce an codon profile
    codon_frequencies = Dict(a => Dict{Kmers.DNACodon, Int}() for a in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term]))
    if allow_all
        for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
            amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
            codon_frequencies[amino_acid][codon] = get(codon_frequencies[amino_acid], codon, 0) + 1
        end    
    end
    genome_genbank_data = GenomicAnnotations.readgbk(genbank)
    for chromosome in genome_genbank_data
        for gene in chromosome.genes
            gene_range = GenomicAnnotations.locus(gene).position
            gene_type = GenomicAnnotations.feature(gene)
            # is_terminator = occursin(r"^TERM", gene.label)
            if (length(gene_range) % 3 == 0) && (gene_type == :misc_feature)

                fw_dnaseq = GenomicAnnotations.sequence(gene)
                revcom_dnaseq = BioSequences.reverse_complement(fw_dnaseq)

                fw_aaseq = BioSequences.translate(fw_dnaseq)
                revcom_aaseq = BioSequences.translate(revcom_dnaseq)

                if last(fw_aaseq) == BioSequences.AA_Term
                    dnaseq = fw_dnaseq
                    aaseq = fw_aaseq
                elseif last(revcom_aaseq) == BioSequences.AA_Term
                    dnaseq = revcom_dnaseq
                    aaseq = revcom_aaseq
                elseif first(fw_aaseq) == BioSequences.AA_M
                    dnaseq = fw_dnaseq
                    aaseq = fw_aaseq
                elseif first(revcom_aaseq) == BioSequences.AA_M
                    dnaseq = revcom_dnaseq
                    aaseq = revcom_aaseq
                else
                    # @show "ambiguous"
                    continue
                end
                
                # for (mer, amino_acid) in zip(BioSequences.each(BioSequences.DNAMer{3}, dnaseq, 3), aaseq)
                for ((i, codon), amino_acid) in zip(Kmers.SpacedKmers{Kmers.DNACodon}(BioSequences.LongDNA{4}(dnaseq), 3), aaseq)
                    @assert amino_acid == first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
                    codon_frequencies[amino_acid][codon] = get(codon_frequencies[amino_acid], codon, 0) + 1
                end
            end
        end
    end
    return codon_frequencies
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Normalizes codon frequencies for each amino acid such that frequencies sum to 1.0.

# Arguments
- `codon_frequencies`: Nested dictionary mapping amino acids to their codon frequency distributions

# Returns
- Normalized codon frequencies where values for each amino acid sum to 1.0
"""
function normalize_codon_frequencies(codon_frequencies)
    normalized_codon_frequencies = Dict{BioSymbols.AminoAcid, Dict{Kmers.DNACodon, Float64}}()
    for (amino_acid, amino_acid_codon_frequencies) in codon_frequencies
        total_count = sum(values(amino_acid_codon_frequencies))
        normalized_codon_frequencies[amino_acid] = Dict(
            amino_acid_codon => amino_acid_codon_frequency/total_count for (amino_acid_codon, amino_acid_codon_frequency) in amino_acid_codon_frequencies
        )
        if !isempty(normalized_codon_frequencies[amino_acid])
            @assert abs(1-sum(values(normalized_codon_frequencies[amino_acid]))) <= eps()
        end
    end
    return normalized_codon_frequencies
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert raw k‑mer counts into normalized frequencies.

# Arguments
- `kmer_counts::Dict`: Mapping of k-mers to counts.

# Returns
`OrderedDict` with values scaled so the sum equals 1.
"""
function normalize_kmer_counts(kmer_counts)
    total_kmer_counts = sum(values(kmer_counts))
    if total_kmer_counts == 0
        return DataStructures.OrderedDict(k => 0.0 for (k, _) in kmer_counts)
    end
    normalized_kmer_frequencies = DataStructures.OrderedDict(k => v/total_kmer_counts for (k,v) in kmer_counts)
    return normalized_kmer_frequencies
end
