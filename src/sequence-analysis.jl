function trim_galore(;outdir="", identifier="")
    
    trim_galore_dir = joinpath(outdir, "trim_galore")
    
    forward_reads = joinpath(outdir, "$(identifier)_1.fastq.gz")
    reverse_reads = joinpath(outdir, "$(identifier)_2.fastq.gz")
    
    trimmed_forward_reads = joinpath(trim_galore_dir, "$(identifier)_1_val_1.fq.gz")
    trimmed_reverse_reads = joinpath(trim_galore_dir, "$(identifier)_2_val_2.fq.gz")
    
    # mamba create -n trim_galore -c bioconda trim_galore
    if !isfile(trimmed_forward_reads) && !isfile(trimmed_reverse_reads)
        cmd = `conda run -n trim_galore trim_galore --suppress_warn --cores $(min(Sys.CPU_THREADS, 4)) --output_dir $(trim_galore_dir) --paired $(forward_reads) $(reverse_reads)`
        run(cmd)
    else
        @info "$(trimmed_forward_reads) & $(trimmed_reverse_reads) already present"
    end
end

function fasterq_dump(;outdir="", srr_identifier="")
    
    forward_reads = joinpath(outdir, "$(srr_identifier)_1.fastq")
    reverse_reads = joinpath(outdir, "$(srr_identifier)_2.fastq")
    
    forward_reads_gz = forward_reads * ".gz"
    reverse_reads_gz = reverse_reads * ".gz"
    
    if !isfile(forward_reads_gz) && !isfile(reverse_reads_gz)
        # --progress doesn't work well for jupyter output
        fasterq_dump_cmd = `
            fasterq-dump
                --outdir $(outdir)
                --mem 1G
                --split-3
                --threads $(min(Sys.CPU_THREADS, 4))
                --skip-technical
                $(srr_identifier)`
        @time run(fasterq_dump_cmd)
        run(`pigz $(forward_reads)`)
        run(`pigz $(reverse_reads)`)
    else
        @info "$(forward_reads_gz) & $(reverse_reads_gz) already present"
    end
end

function download_and_filter_sra_reads(;outdir="", srr_identifier="")
    forward_reads = joinpath(outdir, "$(srr_identifier)_1.fastq")
    reverse_reads = joinpath(outdir, "$(srr_identifier)_2.fastq")
    forward_reads_gz = forward_reads * ".gz"
    reverse_reads_gz = reverse_reads * ".gz"
    trimmed_forward_reads = joinpath(outdir, "trim_galore", "$(srr_identifier)_1_val_1.fq.gz")
    trimmed_reverse_reads = joinpath(outdir, "trim_galore", "$(srr_identifier)_2_val_2.fq.gz")

    if !(isfile(trimmed_forward_reads) && isfile(trimmed_reverse_reads))
        @info "processing $(srr_identifier)"
        fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        trim_galore(outdir=outdir, identifier=srr_identifier)
    # else
        # @info "$(srr_identifier) already processed..."
    end
    isfile(forward_reads_gz) && rm(forward_reads_gz)
    isfile(reverse_reads_gz) && rm(reverse_reads_gz)
end

function amino_acids_to_codons()
    amino_acid_to_codon_map = Dict(a => Kmers.DNACodon for a in vcat(Mycelia.AA_ALPHABET..., [BioSequences.AA_Term]))
    for codon in Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
        amino_acid = first(BioSequences.translate(BioSequences.LongDNA{2}(codon)))
        amino_acid_to_codon_map[amino_acid] = codon
    end   
    return amino_acid_to_codon_map
end

function codons_to_amino_acids()
    codons = Mycelia.generate_all_possible_kmers(3, Mycelia.DNA_ALPHABET)
    codon_to_amino_acid_map = Dict(codon => BioSequences.translate(BioSequences.LongDNA{2}(codon)))
    return codon_to_amino_acid_map
end

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

function codon_optimize(;normalized_codon_frequencies, protein_sequence::BioSequences.LongAA, n_iterations)
    best_sequence = reverse_translate(protein_sequence)
    codons = last.(collect(Kmers.SpacedKmers{Kmers.DNACodon}(BioSequences.LongDNA{4}(best_sequence), 3)))
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
    @show (best_likelihood)^-10 / (initial_log_likelihood)^-10
    return best_sequence
    
end


# function codon_optimize(;normalized_codon_frequencies, optimization_sequence::BioSequences.LongDNA, n_iterations)
#     protein_sequence = BioSequences.translate(optimization_sequence)

# end

function kmer_counts_dict_to_vector(kmer_to_index_map, kmer_counts)
    kmer_counts_vector = zeros(length(kmer_to_index_map))
    for (kmer, count) in kmer_counts
        kmer_index = kmer_to_index_map[kmer]
        kmer_counts_vector[kmer_index] = count
    end
    return kmer_counts_vector
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function generate_all_possible_kmers(k, alphabet)
    kmer_iterator = Iterators.product([alphabet for i in 1:k]...)
    kmer_vectors = collect.(vec(collect(kmer_iterator)))
    if eltype(alphabet) == BioSymbols.AminoAcid
        kmers = [Kmers.Kmer{BioSequences.AminoAcidAlphabet}(BioSequences.LongAA(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.DNA
        kmers = [Kmers.Kmer{BioSequences.DNAAlphabet{2}}(BioSequences.LongDNA{2}(kv)) for kv in kmer_vectors]
    elseif eltype(alphabet) == BioSymbols.RNA
        kmers = [Kmers.Kmer{BioSequences.RNAAlphabet{2}}(BioSequences.LongRNA{2}(kv)) for kv in kmer_vectors]
    else
        error()
    end
    return sort!(kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function generate_all_possible_canonical_kmers(k, alphabet)
    kmers = generate_all_possible_kmers(k, alphabet)
    if eltype(alphabet) == BioSymbols.AminoAcid
        return kmers
    elseif eltype(alphabet) in (BioSymbols.DNA, BioSymbols.RNA)
        return unique!(BioSequences.canonical.(kmers))
    else
        error()
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function get_kmer_index(kmers, kmer)
    index = searchsortedfirst(kmers, kmer)
    @assert kmers[index] == kmer "$kmer not found in kmer list"
    return index
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a dense kmer counts table (canonical for DNA, stranded for RNA & AA) for each fasta provided in a list.
Scales very well for large numbers of organisms/fasta files, but not for k.
Recommended for k <= 13, although 17 may still be possible

```jldoctest
julia> 1 + 1
2
```
"""
function fasta_list_to_dense_counts_table(;fasta_list, k, alphabet)
    k > 13 && error("use fasta_list_to_sparse_counts_table")
    if alphabet == :AA
        KMER_TYPE = BioSequences.AminoAcidAlphabet
        sorted_kmers = sort(generate_all_possible_kmers(k, AA_ALPHABET))
        COUNT = count_kmers
    elseif alphabet == :DNA
        KMER_TYPE = BioSequences.DNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_canonical_kmers(k, DNA_ALPHABET))
        COUNT = count_canonical_kmers
    elseif alphabet == :RNA
        KMER_TYPE = BioSequences.RNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_kmers(k, RNA_ALPHABET))
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    kmer_counts_matrix = zeros(length(sorted_kmers), length(fasta_list))
    progress = ProgressMeter.Progress(length(fasta_list))
    reenrantlock = ReentrantLock()
    Threads.@threads for (entity_index, fasta_file) in collect(enumerate(fasta_list))
        # Acquire the lock before updating the progress
        lock(reenrantlock) do
            # Update the progress meter
            ProgressMeter.next!(progress)
        end
        entity_mer_counts = COUNT(Kmers.Kmer{KMER_TYPE, k}, fasta_file)
        for (i, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[i, entity_index] = get(entity_mer_counts, kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

function biosequences_to_dense_counts_table(;biosequences, k)
    k > 13 && error("use fasta_list_to_sparse_counts_table")    
    if eltype(first(biosequences)) == BioSymbols.AminoAcid
        KMER_TYPE = BioSequences.AminoAcidAlphabet
        sorted_kmers = sort(generate_all_possible_kmers(k, AA_ALPHABET))
        COUNT = count_kmers
    elseif eltype(first(biosequences)) == BioSymbols.DNA
        KMER_TYPE = BioSequences.DNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_canonical_kmers(k, DNA_ALPHABET))
        COUNT = count_canonical_kmers
    elseif eltype(first(biosequences)) == BioSymbols.RNA
        KMER_TYPE = BioSequences.RNAAlphabet{2}
        sorted_kmers = sort(generate_all_possible_kmers(k, RNA_ALPHABET))
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    kmer_counts_matrix = zeros(length(sorted_kmers), length(biosequences))
    progress = ProgressMeter.Progress(length(biosequences))
    reenrantlock = ReentrantLock()
    Threads.@threads for (entity_index, biosequence) in collect(enumerate(biosequences))
        # Acquire the lock before updating the progress
        lock(reenrantlock) do
            # Update the progress meter
            ProgressMeter.next!(progress)
        end
        entity_mer_counts = COUNT(Kmers.Kmer{KMER_TYPE, k}, biosequence)
        for (i, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[i, entity_index] = get(entity_mer_counts, kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a sparse kmer counts table in memory for each fasta provided in a list

```jldoctest
julia> 1 + 1
2
```
"""
# function fasta_list_to_counts_table(;fasta_list::AbstractVector{<:AbstractString}, k, alphabet)
#     if alphabet == :AA
#         KMER_TYPE = Kmers.AAKmer{k}
#         COUNT = count_kmers
#     elseif alphabet == :DNA
#         KMER_TYPE = Kmers.DNAKmer{k}
#         COUNT = count_canonical_kmers
#     elseif alphabet == :RNA
#         KMER_TYPE = Kmers.RNAKmer{k}
#         COUNT = count_kmers
#     else
#         error("invalid alphabet, please choose from :AA, :DNA, :RNA")
#     end
    
#     fasta_kmer_counts_dict = Dict()
#     progress = ProgressMeter.Progress(length(fasta_list))
#     reenrantlock = ReentrantLock()
#     Threads.@threads for fasta_file in fasta_list
#         # Acquire the lock before updating the progress
#         these_kmer_counts = COUNT(KMER_TYPE, fasta_file)
#         lock(reenrantlock) do
#             # Update the progress meter
#             ProgressMeter.next!(progress)
#             fasta_kmer_counts_dict[fasta_file] = these_kmer_counts
#         end
#     end
#     sorted_kmers = sort(collect(union(keys(x) for x in fasta_kmer_counts_dict)))
#     kmer_counts_matrix = zeros(length(sorted_kmers), length(fasta_list))
#     @info "populating sparse counts matrix..."
#     for (col, fasta) in enumerate(fasta_list)
#         for (row, kmer) in enumerate(sorted_kmers)
#             kmer_counts_matrix[row, col] = get(fasta_kmer_counts_dict[fasta], kmer, 0)
#         end
#     end
#     return (;sorted_kmers, kmer_counts_matrix)
# end

function biosequences_to_counts_table(;biosequences, k)
    if eltype(first(biosequences)) == BioSymbols.AminoAcid
        KMER_TYPE = Kmers.AAKmer{k}
        COUNT = count_kmers
    elseif eltype(first(biosequences)) == BioSymbols.DNA
        KMER_TYPE = Kmers.DNAKmer{k}
        COUNT = count_canonical_kmers
    elseif eltype(first(biosequences)) == BioSymbols.RNA
        KMER_TYPE = Kmers.RNAKmer{k}
        COUNT = count_kmers
    else
        error("invalid alphabet, please choose from :AA, :DNA, :RNA")
    end
    
    kmer_counts = Vector{OrderedCollections.OrderedDict{KMER_TYPE, Int}}(undef, length(biosequences))
    progress = ProgressMeter.Progress(length(biosequences))
    reenrantlock = ReentrantLock()
    Threads.@threads for i in eachindex(biosequences)
        lock(reenrantlock) do
            ProgressMeter.next!(progress)
        end
        kmer_counts[i] = COUNT(KMER_TYPE, biosequences[i])
    end
    sorted_kmers = sort(collect(reduce(union, keys.(kmer_counts))))
    kmer_counts_matrix = SparseArrays.spzeros(Int, length(mers), length(biosequences))
    @info "populating sparse counts matrix..."
    for (col, biosequence) in enumerate(biosequences)
        for (row, kmer) in enumerate(sorted_kmers)
            kmer_counts_matrix[row, col] = get(kmer_counts[col], kmer, 0)
        end
    end
    return (;sorted_kmers, kmer_counts_matrix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function normalize_distance_matrix(distance_matrix)
    max_non_nan_value = maximum(filter(x -> !isnan(x) && !isnothing(x) && !ismissing(x), vec(distance_matrix)))
    return distance_matrix ./ max_non_nan_value
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
# function count_matrix_to_probability_matrix(
#         counts_matrix,
#         probability_matrix_file = replace(counts_matrix_file, ".bin" => ".probability_matrix.bin")
#     )
#     probability_matrix = Mmap.mmap(probability_matrix_file, Array{Float64, 2}, size(counts_matrix))
#     if !isfile(probability_matrix_file)
#         println("creating new probability matrix $probability_matrix_file")
#         # probability_matrix .= count_matrix_to_probability_matrix(counts_matrix)
#         for (i, col) in enumerate(eachcol(counts_matrix))
#             probability_matrix[:, i] .= col ./ sum(col)
#         end
#     else
#         println("probability matrix found $probability_matrix_file")
#     end
#     return probability_matrix, probability_matrix_file
# end

function count_matrix_to_probability_matrix(counts_matrix)
    probability_matrix = zeros(size(counts_matrix))
    for (i, col) in enumerate(eachcol(counts_matrix))
        probability_matrix[:, i] .= col ./ sum(col)
    end
    return probability_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function distance_matrix_to_newick(;distance_matrix, labels, outfile)
    # phage_names = phage_host_table[indices, :name]
    # this is equivalent to UPGMA
    tree = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)
    # reference_phage_indices = findall(x -> x in reference_phages, phage_names)
    newick = Dict()
    for row in 1:size(tree.merges, 1)
        left, right = tree.merges[row, :]
        if left < 0
            l = string(labels[abs(left)])
        else
            l = newick[left]
        end
        if right < 0
            r = string(labels[abs(right)])
        else
            r = newick[right]
        end
        height = tree.heights[row]
        newick[row] = "($l:$height, $r:$height)"
    end
    open(outfile, "w") do io
        println(io, newick[size(tree.merges, 1)] * ";")
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

DEPRECATED: THIS WAS THE MEASURE WITH THE LEAST AGREEMENT TO EXISTING MEASURES LIKE BLAST AND % AVERAGE NUCLEOTIDE IDENTITY
Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function counts_matrix_to_size_normalized_cosine_distance_matrix(counts_table)
    n_entities = size(counts_table, 2)
    distance_matrix = zeros(n_entities, n_entities)
    for entity_1_index in 1:n_entities
        for entity_2_index in entity_1_index+1:n_entities
            a = counts_table[:, entity_1_index]
            b = counts_table[:, entity_2_index]
            sa = sum(a)
            sb = sum(b)
            size_dist = 1-(min(sa, sb)/max(sa, sb))
            cosine_dist = Distances.cosine_dist(a, b)
            distances = filter(x -> x > 0, (size_dist, cosine_dist))
            if isempty(distances)
                dist = 0.0
            else
                dist = reduce(*, distances)
            end
            distance_matrix[entity_1_index, entity_2_index] = 
                distance_matrix[entity_2_index, entity_1_index] = dist
        end
    end
    return distance_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create euclidean distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size)

```jldoctest
julia> 1 + 1
2
```
"""
function frequency_matrix_to_euclidean_distance_matrix(counts_table)
    n_entities = size(counts_table, 2)
    distance_matrix = zeros(n_entities, n_entities)
    for entity_1_index in 1:n_entities
        for entity_2_index in entity_1_index+1:n_entities
            a = counts_table[:, entity_1_index]
            b = counts_table[:, entity_2_index]
            distance_matrix[entity_1_index, entity_2_index] = 
                distance_matrix[entity_2_index, entity_1_index] = 
                Distances.euclidean(a, b)
        end
    end
    return distance_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create cosine distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function frequency_matrix_to_cosine_distance_matrix(probability_matrix)
    n_entities = size(probability_matrix, 2)
    distance_matrix = zeros(n_entities, n_entities)
    for entity_1_index in 1:n_entities
        for entity_2_index in entity_1_index+1:n_entities
            a = probability_matrix[:, entity_1_index]
            b = probability_matrix[:, entity_2_index]
            distance_matrix[entity_1_index, entity_2_index] = 
                distance_matrix[entity_2_index, entity_1_index] = 
                Distances.cosine_dist(a, b)
        end
    end
    return distance_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function kmer_path_to_sequence(kmer_path)
    sequence = BioSequences.LongDNA{2}(first(kmer_path))
    for kmer in kmer_path[2:end]
        for i in 1:length(kmer)-1
            a = kmer[i]
            b = sequence[end-(length(kmer)-1)+i]
            @assert a == b
        end
        push!(sequence, kmer[end])
    end
    return sequence
end
    
# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function observe(records::AbstractVector{R}; outfile = "", error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
#     if isempty(outfile)
#         error("no file name supplied")
#     end
#     io = open(outfile, 'w')
#     fastx_io = FASTX.FASTQ.Writer(io)
#     for record in records
#         new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
#         new_seq_id = string(hash(new_seq)) * "-" * Random.randstring(32)
#     new_seq_description = FASTX.identifier(record)
#     quality = fill(UInt8(60), length(new_seq))
#     return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function random_fasta_record(;seed=rand(0:typemax(Int)), L = rand(0:Int(typemax(UInt16))))
    id = Random.randstring(Int(ceil(log(L + 1))))
    seq = BioSequences.randdnaseq(Random.seed!(seed), L)
    return FASTX.FASTA.Record(id, seq)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function observe(record::R; error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    
    new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
    new_seq_id = string(hash(new_seq)) * "-" * Random.randstring(32)
    new_seq_description = FASTX.identifier(record)
    quality = fill(UInt8(60), length(new_seq))
    return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality)
end

function observe(sequence::BioSequences.LongSequence{T}; error_rate = 0.0) where T
    
    if T <: BioSequences.DNAAlphabet
        alphabet = DNA_ALPHABET
    elseif T <: BioSequences.RNAAlphabet
        alphabet = RNA_ALPHABET
    else
        @assert T <: BioSequences.AminoAcidAlphabet
        alphabet = AA_ALPHABET
    end
    
    new_seq = BioSequences.LongSequence{T}()
    for character in sequence
        if rand() > error_rate
            # match
            push!(new_seq, character)
        else
            error_type = rand(1:3)
            if error_type == 1
                # mismatch
                new_character = rand(alphabet)
                while new_character == character
                    new_character = rand(alphabet)
                end
                push!(new_seq, new_character)
            elseif error_type == 2
                # insertion
                total_insertions = 1 + rand(Distributions.Poisson(error_rate))
                for i in 1:total_insertions
                    push!(new_seq, rand(alphabet))
                end
                push!(new_seq, character)
            else
                # deletion
                continue
            end
        end
    end
    if (T <: BioSequences.DNAAlphabet || T <: BioSequences.RNAAlphabet) && rand(Bool)
        BioSequences.reverse_complement!(new_seq)
    end
    return new_seq
end

function observe(records::AbstractVector{R};
                weights=ones(length(records)),
                N = length(records),
                outfile = "",
                error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    if isempty(outfile)
        error("no file name supplied")
    end
    io = open(outfile, "w")
    fastx_io = FASTX.FASTA.Writer(io)
    for i in 1:N
        record = StatsBase.sample(records, StatsBase.weights(weights))
        new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
        new_seq_id = Random.randstring(Int(ceil(log(length(new_seq) + 1))))
        new_seq_description = FASTX.identifier(record)
        observed_record = FASTX.FASTA.Record(new_seq_id, new_seq_description, new_seq)
        write(fastx_io, observed_record)
    end
    close(fastx_io)
    close(io)
    return outfile
end

# currently this is only for amino acid sequences, expand to include DNA and RNA via multiple dispatch
function mutate_sequence(reference_sequence)
    i = rand(1:length(reference_sequence))
    mutation_type = rand([SequenceVariation.Substitution, SequenceVariation.Insertion, SequenceVariation.Deletion])
    if mutation_type == SequenceVariation.Substitution
        # new_amino_acid = BioSequences.AA_A
        new_amino_acid = rand(amino_acids)
        mutation_string = "$(reference_sequence[i])$(i)$(new_amino_acid)"
    else
        indel_size = rand(Distributions.truncated(Distributions.Poisson(1), lower=1))
        if mutation_type == SequenceVariation.Insertion
            inserted_sequence = join([rand(amino_acids) for i in 1:indel_size], "")
            mutation_string = "$(i)$(inserted_sequence)"
        else
            @assert mutation_type == SequenceVariation.Deletion
            i_stop = min(i+indel_size, length(reference_sequence))
            mutation_string = "Δ$(i)-$(i_stop)"
        end
    end
    # println(mutation_string)
    mutation = SequenceVariation.Variation(reference_sequence, mutation_string)
    haplotype = SequenceVariation.Haplotype(reference_sequence, [mutation])
    mutant_sequence = SequenceVariation.reconstruct(haplotype)
    return mutant_sequence, haplotype
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return proportion of matched bases in alignment to total matches + edits.

0-1, not %

```jldoctest
julia> 1 + 1
2
```
"""
function assess_alignment_accuracy(alignment_result)
    return alignment_result.total_matches / (alignment_result.total_matches + alignment_result.total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Used to determine which orientation provides an optimal alignment for initiating path likelihood analyses in viterbi analysis

```jldoctest
julia> 1 + 1
2
```
"""
function assess_optimal_kmer_alignment(kmer, observed_kmer)

    forward_alignment_result = assess_alignment(kmer, observed_kmer)
    forward_alignment_accuracy = assess_alignment_accuracy(forward_alignment_result)

    reverse_alignment_result = assess_alignment(kmer, BioSequences.reverse_complement(observed_kmer))
    reverse_alignment_accuracy = assess_alignment_accuracy(reverse_alignment_result)

    if forward_alignment_accuracy > reverse_alignment_accuracy
        alignment_result = forward_alignment_result
        orientation = true
    elseif forward_alignment_accuracy < reverse_alignment_accuracy
        alignment_result = reverse_alignment_result
        orientation = false
    elseif forward_alignment_accuracy == reverse_alignment_accuracy
        alignment_result, orientation = rand(((forward_alignment_result, missing), (reverse_alignment_result, missing)))
    end

    return (alignment_result, orientation)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function assess_alignment(a, b)
    pairwise_alignment = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), a, b)
    alignment_result = BioAlignments.alignment(pairwise_alignment)
    total_aligned_bases = BioAlignments.count_aligned(alignment_result)
    total_matches = Int(BioAlignments.count_matches(alignment_result))
    total_edits = Int(total_aligned_bases - total_matches)
    return (total_matches = total_matches, total_edits = total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function canonicalize_kmer_counts!(kmer_counts)
    for (kmer, count) in kmer_counts
        if !BioSequences.iscanonical(kmer)
            canonical_kmer = BioSequences.canonical(kmer)
            if haskey(kmer_counts, canonical_kmer)
                kmer_counts[canonical_kmer] += count
            else
                kmer_counts[canonical_kmer] = count
            end
            delete!(kmer_counts, kmer)
        end
    end
    return sort!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function canonicalize_kmer_counts(kmer_counts)
    return canonicalize_kmer_counts!(copy(kmer_counts))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function count_canonical_kmers(::Type{KMER_TYPE}, sequences) where KMER_TYPE
    kmer_counts = count_kmers(KMER_TYPE, sequences)
    return canonicalize_kmer_counts!(kmer_counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function count_kmers(::Type{KMER_TYPE}, sequence::BioSequences.LongSequence) where KMER_TYPE
    return sort(StatsBase.countmap([kmer for (index, kmer) in Kmers.EveryKmer{KMER_TYPE}(sequence)]))
end

function count_kmers(::Type{KMER_TYPE}, record::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    # TODO: need to figure out how to infer the sequence type
    if eltype(KMER_TYPE) == BioSymbols.DNA
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongDNA{4}, record))
    elseif eltype(KMER_TYPE) == BioSymbols.RNA
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongRNA{4}, record))
    elseif eltype(KMER_TYPE) == BioSymbols.AminoAcid
        return count_kmers(KMER_TYPE, FASTX.sequence(BioSequences.LongAA, record))
    else
        @error KMER_TYPE
    end
end

function count_kmers(::Type{KMER_TYPE}, records::AbstractVector{T}) where {KMER_TYPE, T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    kmer_counts = count_kmers(KMER_TYPE, first(records))
    for record in records[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, record)
        merge!(+, kmer_counts, _kmer_counts)
    end
    sort!(kmer_counts)
end

function count_kmers(::Type{KMER_TYPE}, sequences::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Reader, FASTX.FASTQ.Reader}}
    return count_kmers(KMER_TYPE, collect(sequences))
end

function count_kmers(::Type{KMER_TYPE}, fastx_files::AbstractVector{T}) where {KMER_TYPE, T <: AbstractString}
    kmer_counts = count_kmers(KMER_TYPE, first(fastx_files))
    for file in fastx_files[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, file)
        kmer_counts = merge!(+, kmer_counts, _kmer_counts)
    end
    return kmer_counts
end

function count_kmers(::Type{KMER_TYPE}, fastx_file::AbstractString) where {KMER_TYPE}
    fastx_io = open_fastx(fastx_file)
    kmer_counts = count_kmers(KMER_TYPE, fastx_io)
    close(fastx_io)
    return kmer_counts
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function q_value_to_error_rate(q_value)
    error_rate = 10^(q_value/(-10))
    return error_rate
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function error_rate_to_q_value(error_rate)
    q_value = -10 * log10(error_rate)
    return q_value
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function is_equivalent(a, b)
    a == b || a == BioSequences.reverse_complement(b)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function edge_path_to_sequence(kmer_graph, edge_path)
    edge = first(edge_path)
    sequence = BioSequences.LongDNASeq(kmers[edge.src])
    if !kmer_graph.eprops[edge][:orientations].source_orientation
        sequence = BioSequences.reverse_complement(sequence)
    end
    for edge in edge_path
        destination = BioSequences.LongDNASeq(kmers[edge.dst])
        if !kmer_graph.eprops[edge][:orientations].destination_orientation
            destination = BioSequences.reverse_complement(destination)
        end
        sequence_suffix = sequence[end-length(destination)+2:end]
        destination_prefix = destination[1:end-1]
        @assert sequence_suffix == destination_prefix
        push!(sequence, destination[end])
    end
    sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

This turns a 4-line FASTQ entry into a single tab separated line,
adds a column with the length of each read, passes it to Unix sort,
removes the length column, and converts it back into a FASTQ file.

sorts longest to shortest!!

http://thegenomefactory.blogspot.com/2012/11/sorting-fastq-files-by-sequence-length.html
```jldoctest
julia> 1 + 1
2
```
"""
function sort_fastq(input_fastq, output_fastq="")
    
    if endswith(input_fastq, ".gz")
        p = pipeline(
                `gzip -dc $input_fastq`,
                `paste - - - -`,
                `perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'`,
                `sort -nr`,
                `cut -f2-`,
                `tr "\t" "\n"`,
                `gzip`
                )
    else
        p = pipeline(
                `cat $input_fastq`,
                `paste - - - -`,
                `perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'`,
                `sort -nr`,
                `cut -f2-`,
                `tr "\t" "\n"`
                )
    end
    run(pipeline(p, output_fastq))
    return output_fastq
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function count_records(fastx)
    n_records = 0
    for record in open_fastx(fastx)
        n_records += 1
    end
    return n_records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_read_lengths(fastq_file; total_reads = Inf)
    if total_reads == Inf
        total_reads = count_records(fastq_file)
    end
    read_lengths = zeros(Int, total_reads)
    @info "determining read lengths"
    p = ProgressMeter.Progress(total_reads, 1)
    for (i, record) in enumerate(open_fastx(fastq_file))
#         push!(read_lengths, length(FASTX.sequence(record)))
        read_lengths[i] = length(FASTX.sequence(record))
        ProgressMeter.next!(p)
    end
    return read_lengths
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_max_canonical_kmers(k, ALPHABET)
    max_possible_kmers = determine_max_possible_kmers(k, ALPHABET)
    if ALPHABET == AA_ALPHABET
        return max_possible_kmers
    else
        @assert isodd(k) "this calculation is not valid for even length kmers"
        return Int(max_possible_kmers / 2)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function determine_max_possible_kmers(k, ALPHABET)
    return length(ALPHABET)^k
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
# Michaelis–Menten
function calculate_v(s,p)
    vmax = p[1]
    km = p[2]
    v = (vmax .* s) ./ (km .+ s)
    return v
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function translate_nucleic_acid_fasta(fasta_nucleic_acid_file, fasta_amino_acid_file)
    open(fasta_amino_acid_file, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for record in FASTX.FASTA.Reader(open(fasta_nucleic_acid_file))
            try
                raw_seq = FASTX.sequence(record)
                pruned_seq_length = Int(floor(length(raw_seq)/3)) * 3
                truncated_seq = raw_seq[1:pruned_seq_length]
                amino_acid_seq = BioSequences.translate(truncated_seq)
                amino_acid_record = FASTX.FASTA.Record(FASTX.identifier(record), FASTX.description(record), amino_acid_seq)
                write(writer, amino_acid_record)
            catch
                @warn "unable to translate record", record
            end
        end
        close(writer)
    end
    return fasta_amino_acid_file
end

function fasta_to_table(fasta)
    collected_fasta = collect(fasta)
    fasta_df = DataFrames.DataFrame(
        identifier = FASTX.identifier.(collected_fasta),
        description = FASTX.description.(collected_fasta),
        sequence = FASTX.sequence.(collected_fasta)
    )
    return fasta_df
end

function fasta_table_to_fasta(fasta_df)
    records = Vector{FASTX.FASTA.Record}(undef, DataFrames.nrow(fasta_df))
    for (i, row) in enumerate(DataFrames.eachrow(fasta_df))
        record = FASTX.FASTA.Record(row["identifier"], row["description"], row["sequence"])
        records[i] = record
    end
    return records
end

function deduplicate_fasta_file(in_fasta, out_fasta)
    fasta_df = fasta_to_table(collect(open_fastx(in_fasta)))
    sort!(fasta_df, "identifier")
    unique_sequences = DataFrames.combine(DataFrames.groupby(fasta_df, "sequence"), first)
    fasta = fasta_table_to_fasta(unique_sequences)
    open(out_fasta, "w") do io
        writer = FASTX.FASTA.Writer(io)
        for record in fasta
            write(writer, record)
        end
        close(writer)
    end
    return out_fasta
end


# seqkit concat $(cat fasta_files.txt) > merged.fasta
# seqtk seq -L $(cat fasta_files.txt) > merged.fasta
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Join fasta files without any regard to record uniqueness.

A cross-platform version of `cat *.fasta > joint.fasta`

See merge_fasta_files
"""
function concatenate_files(;files, file)
    close(open(file, "w"))
    ProgressMeter.@showprogress for f in files
        # stderr=file_path
        run(pipeline(`cat $(f)`, stdout=file, append=true))
    end
    return file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Join fasta files while adding origin prefixes to the identifiers.

Does not guarantee uniqueness but will warn if conflicts arise
"""
function merge_fasta_files(;fasta_files, fasta_file)
    @info "merging $(length(fasta_files)) files..."
    identifiers = Set{String}()
    open(fasta_file, "w") do io
        fastx_io = FASTX.FASTA.Writer(io)
        ProgressMeter.@showprogress for f in fasta_files
            f_id = replace(basename(f), Mycelia.FASTA_REGEX => "")
            for record in Mycelia.open_fastx(f)
                new_record_id = f_id * "__" * FASTX.identifier(record)
                if new_record_id in identifiers
                    @warn "new identifier $(new_record_id) already in identifiers!!!"
                end
                push!(identifiers, new_record_id)
                new_record = FASTX.FASTA.Record(new_record_id, FASTX.sequence(record))
                write(fastx_io, new_record)
            end
        end
    end
    @info "$(length(identifiers)) records merged..."
    return fasta_file
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function path_to_sequence(kmer_graph, path)
#     sequence = BioSequences.LongDNASeq(oriented_kmer_to_sequence(kmer_graph, first(path)))
#     for oriented_kmer in path[2:end]
#         nucleotide = last(oriented_kmer_to_sequence(kmer_graph, oriented_kmer))
#         push!(sequence, nucleotide)
#     end
#     return sequence
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function oriented_kmer_to_sequence(kmer_graph, oriented_kmer)
#     kmer_sequence = kmer_graph.kmers[oriented_kmer.index]
#     if !oriented_kmer.orientation
#         kmer_sequence = BioSequences.reverse_complement(kmer_sequence)
#     end
#     return kmer_sequence
# end