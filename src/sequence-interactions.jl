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
        kmers = BioSequences.LongAminoAcidSeq.(kmer_vectors)
    elseif eltype(alphabet) == BioSymbols.DNA
        kmers = BioSequences.LongDNASeq.(kmer_vectors)
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
    elseif eltype(alphabet) == BioSymbols.DNA
        return BioSequences.DNAMer.(unique!(BioSequences.canonical.(kmers)))
    else
        error()
    end
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
function count_aamers_by_file(k, fastx_file)
    kmer_counts = StatsBase.countmap(FASTX.sequence(record)[i:i+k-1] for record in Mycelia.open_fastx(fastx_file) for i in 1:length(FASTX.sequence(record))-k+1)
    kmer_counts = sort(kmer_counts)
    kmer_counts_table = 
    DataFrames.DataFrame(
        kmer = collect(keys(kmer_counts)),
        count = collect(values(kmer_counts))
    )
    return kmer_counts_table
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
function count_aamers_by_record(k, fastx_file)
    kmer_counts_table = 
    DataFrames.DataFrame(
        record_identifier = String[],
        kmer = BioSequences.LongAminoAcidSeq[],
        count = Int[]
        )
    for record in Mycelia.open_fastx(fastx_file)
        kmer_counts = StatsBase.countmap(FASTX.sequence(record)[i:i+k-1] for i in 1:length(FASTX.sequence(record))-k+1)
        kmer_counts = sort(kmer_counts)
        for (kmer, count) in sort(kmer_counts)
            row = (
                record_identifier = FASTX.identifier(record),
                kmer = kmer,
                count = count
                )
            push!(kmer_counts_table, row)
        end
    end
    return kmer_counts_table
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
function count_aamers(k, fasta_proteins::FP) where {FP <: Union{AbstractVector{FASTX.FASTA.Record}, FASTX.FASTA.Reader}}
    aamer_counts = OrderedCollections.OrderedDict{BioSequences.LongAminoAcidSeq, Int64}()
    for protein in fasta_proteins
        if !(FASTX.sequence(protein) isa BioSequences.LongAminoAcidSeq)
            # @warn "record $(protein) is not encoded as a protein sequence, skipping..."
            continue
        end
        these_counts = count_aamers(k, protein)
        merge!(+, aamer_counts, these_counts)
    end
    return sort(aamer_counts)
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
function count_aamers(k, fasta_protein::FASTX.FASTA.Record)
    s = FASTX.sequence(fasta_protein)
    these_counts = sort(StatsBase.countmap([s[i:i+k-1] for i in 1:length(s)-k-1]))
    return these_counts    
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
function update_counts_matrix!(matrix, sample_index, countmap, sorted_kmers)
    for (i, kmer) in enumerate(sorted_kmers)
        matrix[i, sample_index] = get(countmap, kmer, 0)
    end
    return matrix
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
function fasta_list_to_counts_table(;fasta_list, k, alphabet, outfile="")
    if alphabet == :AA
        canonical_mers = generate_all_possible_canonical_kmers(k, Mycelia.AA_ALPHABET)
    elseif alphabet == :DNA
        canonical_mers = generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET)
    else
        error("invalid alphabet")
    end
    if isempty(outfile)
        outfile = joinpath(pwd(), "$(hash(fasta_list)).$(alphabet).k$(k).bin")
    end
    if isfile(outfile)
        println("$outfile found, loading into memory")
        mer_counts_matrix = Mmap.mmap(open(outfile), Array{Int, 2}, (length(canonical_mers), length(fasta_list)))
    else
        println("creating new counts matrix $outfile")
        mer_counts_matrix = Mmap.mmap(open(outfile, "w+"), Array{Int, 2}, (length(canonical_mers), length(fasta_list)))
        mer_counts_matrix .= 0
        ProgressMeter.@showprogress for (entity_index, fasta_file) in enumerate(fasta_list)
            if alphabet == :DNA
                entity_mer_counts = Mycelia.count_canonical_kmers(BioSequences.DNAMer{dna_k}, fasta_file)
            elseif alphabet == :AA
                # faa_file = "$(accession).fna.faa"
                # if !isfile(faa_file)
                #     run(pipeline(`prodigal -i $(fna_file) -o $(fna_file).genes -a $(faa_file) -p meta`, stderr="$(fna_file).prodigal.stderr"))
                # end
                try
                    entity_mer_counts = count_aamers(k, Mycelia.open_fastx(fasta_file))
                catch e
                    @warn "investigate possible malformed fasta file $(fasta_file)"
                    # println(e.msg)
                    showerror(stdout, e)
                    continue
                end
            end
            update_counts_matrix!(mer_counts_matrix, entity_index, entity_mer_counts, canonical_mers)            
        end
    end
    return mer_counts_matrix, outfile
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
function count_matrix_to_probability_matrix(counts_matrix, counts_matrix_file)
    probability_matrix_file = replace(counts_matrix_file, ".bin" => ".probability_matrix.bin")
    already_there = isfile(probability_matrix_file)
    probability_matrix = Mmap.mmap(probability_matrix_file, Array{Float64, 2}, size(counts_matrix))
    if !already_there
        println("creating new probability matrix $probability_matrix_file")
        for (i, col) in enumerate(eachcol(counts_matrix))
            probability_matrix[:, i] .= col ./ sum(col)
        end
    else
        println("probability matrix found $probability_matrix_file")
    end
    return probability_matrix, probability_matrix_file
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
function distance_matrix_to_newick(distance_matrix, labels, outfile)
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

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

```jldoctest
julia> 1 + 1
2
```
"""
function counts_matrix_to_distance_matrix(counts_table)
    # TODO, if a file path is provided, make this a mmap table and return that
    distance_matrix = zeros(size(counts_table, 2), size(counts_table, 2))
    for i1 in 1:size(counts_table, 2)
        for i2 in i1+1:size(counts_table, 2)
            a = counts_table[:, i1]
            b = counts_table[:, i2]
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
            distance_matrix[i1, i2] = distance_matrix[i2, i1] = dist
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
    sequence = BioSequences.LongDNASeq(first(kmer_path))
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function generate_all_possible_kmers(k)
    product_iterator = Iterators.product([Mycelia.DNA_ALPHABET for i in 1:k]...)
    return sort(vec(BioSequences.BigDNAMer.(product_iterator)))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function random_fasta_record(;seed=rand(Int), L = rand(0:Int(typemax(UInt16))))
    id = Random.randstring(Int(ceil(log(L + 1))))
    seq = BioSequences.randdnaseq(Random.seed!(seed), L)
    return FASTX.FASTA.Record(id, seq)
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
function observe(record::R; error_rate = 0.0) where {R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    
    new_seq = observe(FASTX.sequence(record), error_rate=error_rate)
    new_seq_id = string(hash(new_seq)) * "-" * Random.randstring(32)
    new_seq_description = FASTX.identifier(record)
    quality = fill(UInt8(60), length(new_seq))
    return FASTX.FASTQ.Record(new_seq_id, new_seq_description, new_seq, quality)
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
                push!(new_seq, rand(setdiff(alphabet, character)))
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
        new_seq = Mycelia.observe(FASTX.sequence(record), error_rate=error_rate)
        new_seq_id = Random.randstring(Int(ceil(log(length(new_seq) + 1))))
        new_seq_description = FASTX.identifier(record)
        observed_record = FASTX.FASTA.Record(new_seq_id, new_seq_description, new_seq)
        write(fastx_io, observed_record)
    end
    close(fastx_io)
    close(io)
    return outfile
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
function count_canonical_kmers(::Type{KMER_TYPE}, sequence::BioSequences.LongSequence) where KMER_TYPE
    canonical_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    canonical_kmer_iterator = (BioSequences.canonical(kmer.fw) for kmer in BioSequences.each(KMER_TYPE, sequence))
    for canonical_kmer in canonical_kmer_iterator
        canonical_kmer_counts[canonical_kmer] = get(canonical_kmer_counts, canonical_kmer, 0) + 1
    end
    return canonical_kmer_counts
end

function count_canonical_kmers(::Type{KMER_TYPE}, record::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    return count_canonical_kmers(KMER_TYPE, FASTX.sequence(record))    
end

function count_canonical_kmers(::Type{KMER_TYPE}, sequences::AbstractVector{T}) where {KMER_TYPE, T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    joint_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    for sequence in sequences
        sequence_kmer_counts = count_canonical_kmers(KMER_TYPE, sequence)
        merge!(+, joint_kmer_counts, sequence_kmer_counts)
    end
    sort!(joint_kmer_counts)
end

function count_canonical_kmers(::Type{KMER_TYPE}, sequences::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Reader, FASTX.FASTQ.Reader}}
    joint_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    for sequence in sequences
        sequence_kmer_counts = count_canonical_kmers(KMER_TYPE, sequence)
        merge!(+, joint_kmer_counts, sequence_kmer_counts)
    end
    sort!(joint_kmer_counts)
end

function count_canonical_kmers(::Type{KMER_TYPE}, fastx_files::AbstractVector{S}) where {KMER_TYPE, S <: AbstractString}
    kmer_counts = count_canonical_kmers(KMER_TYPE, first(fastx_files))
    for fastx_file in fastx_files[2:end]
        _kmer_counts = count_canonical_kmers(KMER_TYPE, fastx_file)
        kmer_counts = merge!(+, kmer_counts, _kmer_counts)
    end
    return kmer_counts
end

function count_canonical_kmers(::Type{KMER_TYPE}, fastx_file::S) where {KMER_TYPE, S <: AbstractString}
    fastx_io = open_fastx(fastx_file)
    kmer_counts = Mycelia.count_canonical_kmers(KMER_TYPE, fastx_io)
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
function count_kmers(::Type{KMER_TYPE}, sequence::BioSequences.LongSequence) where KMER_TYPE
    kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    kmer_iterator = (kmer.fw for kmer in BioSequences.each(KMER_TYPE, sequence))
    for kmer in kmer_iterator
        kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
    end
    return kmer_counts
end

function count_kmers(::Type{KMER_TYPE}, record::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    return count_kmers(KMER_TYPE, FASTX.sequence(record))    
end

function count_kmers(::Type{KMER_TYPE}, sequences::AbstractVector{T}) where {KMER_TYPE, T <: Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}}
    joint_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    for sequence in sequences
        sequence_kmer_counts = count_kmers(KMER_TYPE, sequence)
        merge!(+, joint_kmer_counts, sequence_kmer_counts)
    end
    sort!(joint_kmer_counts)
end

function count_kmers(::Type{KMER_TYPE}, sequences::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Reader, FASTX.FASTQ.Reader}}
    joint_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    for sequence in sequences
        sequence_kmer_counts = count_kmers(KMER_TYPE, sequence)
        merge!(+, joint_kmer_counts, sequence_kmer_counts)
    end
    sort!(joint_kmer_counts)
end

function count_kmers(KMER_TYPE, fastx_files::AbstractVector{AbstractString})
    kmer_counts = count_kmers(KMER_TYPE, first(files))
    for file in files[2:end]
        _kmer_counts = count_kmers(KMER_TYPE, file)
        kmer_counts = merge!(+, kmer_counts, _kmer_counts)
    end
    return kmer_counts
end

function count_kmers(KMER_TYPE, fastx_file::AbstractString)
    fastx_io = open_fastx(fastx_file)
    kmer_counts = Mycelia.count_kmers(KMER_TYPE, fastx_io)
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
function assess_dnamer_saturation(fastxs, kmer_type; kmers_to_assess=Inf, power=10)
    canonical_kmers = Set{kmer_type}()
    
    max_possible_kmers = max_canonical_kmers(kmer_type)
    
    if kmers_to_assess == Inf
        kmers_to_assess = max_possible_kmers
    end
    
    sampling_points = Int[0]
    i = 0
    while power^i <= kmers_to_assess
        push!(sampling_points, power^i)
        i += 1
    end
    
    unique_kmer_counts = zeros(Int, length(sampling_points))
    
    if length(sampling_points) < 3
        @info "increase the # of reads analyzed or decrease the power to acquire more data points"
        return (;sampling_points, unique_kmer_counts)
    end
    
    p = ProgressMeter.Progress(kmers_to_assess, 1)
    
    kmers_assessed = 0
    for fastx in fastxs
        for record in Mycelia.open_fastx(fastx)
            for kmer in BioSequences.each(kmer_type, FASTX.sequence(record))
                canonical_kmer = kmer.fw < kmer.bw ? kmer.fw : kmer.bw
                push!(canonical_kmers, canonical_kmer)
                kmers_assessed += 1
                if (length(canonical_kmers) == max_possible_kmers)                 
                    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
                    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], length(canonical_kmers))
                    return (;sampling_points, unique_kmer_counts, eof = false)
                elseif kmers_assessed in sampling_points
                    i = findfirst(sampling_points .== kmers_assessed)
                    unique_kmer_counts[i] = length(canonical_kmers)
                    if i == length(sampling_points)
                        return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = false)
                    end
                end
                ProgressMeter.next!(p)
            end
        end
    end
    sampling_points = vcat(filter(s -> s < kmers_assessed, sampling_points), [kmers_assessed])
    unique_kmer_counts = vcat(unique_kmer_counts[1:length(sampling_points)-1], [length(canonical_kmers)])    
    return (sampling_points = sampling_points, unique_kmer_counts = unique_kmer_counts, eof = true)
end

function assess_dnamer_saturation(fastxs; outdir="", min_k=3, max_k=31, threshold=0.1)
    
    if isempty(outdir)
        outdir = joinpath(pwd(), "kmer-saturation")
    end
    mkpath(outdir)
    
    ks = Primes.primes(min_k, max_k)
    minimum_saturation = Inf
    midpoint = Inf
    for k in ks
        kmer_type = BioSequences.BigDNAMer{k}
        kmers_to_assess = 10_000_000
        sampling_points, kmer_counts, hit_eof = assess_dnamer_saturation(fastxs, kmer_type, kmers_to_assess=kmers_to_assess)
        @show sampling_points, kmer_counts, hit_eof
        observed_midpoint_index = findfirst(i -> kmer_counts[i] > last(kmer_counts)/2, 1:length(sampling_points))
        observed_midpoint = sampling_points[observed_midpoint_index]
        initial_parameters = Float64[maximum(kmer_counts), observed_midpoint]
        @time fit = LsqFit.curve_fit(calculate_v, sampling_points, kmer_counts, initial_parameters)
        if hit_eof
            inferred_maximum = last(kmer_counts)
        else
            inferred_maximum = max(Int(ceil(fit.param[1])), last(kmer_counts))
        end

        inferred_midpoint = Int(ceil(fit.param[2]))
        predicted_saturation = inferred_maximum / max_canonical_kmers(kmer_type)
        @show k, predicted_saturation

        p = StatsPlots.scatter(
            sampling_points,
            kmer_counts,
            label="observed kmer counts",
            ylabel="# unique kmers",
            xlabel="# kmers assessed",
            title = "sequencing saturation @ k = $k",
            legend=:outertopright,
            size=(800, 400),
            margins=3Plots.PlotMeasures.mm
            )
        StatsPlots.hline!(p, [max_canonical_kmers(kmer_type)], label="absolute maximum")
        StatsPlots.hline!(p, [inferred_maximum], label="inferred maximum")
        StatsPlots.vline!(p, [inferred_midpoint], label="inferred midpoint")
        # xs = vcat(sampling_points, [last(sampling_points) * 2^i for i in 1:2])
        xs = sort([sampling_points..., inferred_midpoint])
        ys = calculate_v(xs, fit.param)
        StatsPlots.plot!(
            p,
            xs,
            ys,
            label="fit trendline")
        display(p)
        StatsPlots.savefig(p, joinpath(outdir, "$k.png"))
        StatsPlots.savefig(p, joinpath(outdir, "$k.svg"))

        if predicted_saturation < minimum_saturation
            minimum_saturation = predicted_saturation
            min_k = k
            midpoint = inferred_midpoint 
        end
        if predicted_saturation < threshold
            chosen_k_file = joinpath(outdir, "chosen_k.txt")
            println("chosen k = $k")
            open(chosen_k_file, "w") do io
                println(io, k)
            end
            return k
        end
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
    for record in Mycelia.open_fastx(fastx)
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
        total_reads = count_reads(fastq_file)
    end
    read_lengths = zeros(Int, total_reads)
    @info "determining read lengths"
    p = ProgressMeter.Progress(total_reads, 1)
    for (i, record) in enumerate(Mycelia.open_fastx(fastq_file))
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
function max_canonical_kmers(kmer_type)
    k_size = last(kmer_type.parameters)
    # we only consider canonical kmers so cut in 1/2
    max_canonical_kmers = Int(4^k_size / 2)
    return max_canonical_kmers
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



"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""