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
        kmers = BioSequences.LongAA.(kmer_vectors)
        # filter out any kmers where the stop codon is not the last codon
        # may want to add a similar filter step for if the start codon is not the first codon
        # filter!(kmer -> any(BioSymbols.isterm.(kmer[1:end-1])), kmers)     
    elseif eltype(alphabet) == BioSymbols.DNA
        kmers = Kmers.DNAKmer.(BioSequences.LongDNA{2}.(kmer_vectors))
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
        return Kmers.DNAKmer.(unique!(BioSequences.canonical.(kmers)))
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
        canonical_mers = generate_all_possible_canonical_kmers(k, AA_ALPHABET)
    elseif alphabet == :DNA
        canonical_mers = generate_all_possible_canonical_kmers(k, DNA_ALPHABET)
    elseif alphabet == :RNA
        canonical_mers = generate_all_possible_canonical_kmers(k, RNA_ALPHABET)
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
                entity_mer_counts = count_canonical_kmers(BioSequences.DNAMer{k}, fasta_file)
            elseif alphabet == :RNA
                entity_mer_counts = count_canonical_kmers(BioSequences.RNAMer{k}, fasta_file)
            elseif alphabet == :AA
                # faa_file = "$(accession).fna.faa"
                # if !isfile(faa_file)
                #     run(pipeline(`prodigal -i $(fna_file) -o $(fna_file).genes -a $(faa_file) -p meta`, stderr="$(fna_file).prodigal.stderr"))
                # end
                try
                    entity_mer_counts = count_aamers(k, open_fastx(fasta_file))
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
function count_matrix_to_probability_matrix(
        counts_matrix,
        probability_matrix_file = replace(counts_matrix_file, ".bin" => ".probability_matrix.bin")
    )
    probability_matrix = Mmap.mmap(probability_matrix_file, Array{Float64, 2}, size(counts_matrix))
    if !isfile(probability_matrix_file)
        println("creating new probability matrix $probability_matrix_file")
        # probability_matrix .= count_matrix_to_probability_matrix(counts_matrix)
        for (i, col) in enumerate(eachcol(counts_matrix))
            probability_matrix[:, i] .= col ./ sum(col)
        end
    else
        println("probability matrix found $probability_matrix_file")
    end
    return probability_matrix, probability_matrix_file
end

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
function counts_matrix_to_cosine_distance_matrix(counts_table)
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

function counts_matrix_to_euclidean_distance_matrix(counts_table)
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
function random_fasta_record(;seed=rand(Int), L = rand(0:Int(typemax(UInt16))))
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
    canonical_kmer_iterator = Kmers.EveryCanonicalKmer{KMER_TYPE}(sequence)
    for (index, canonical_kmer) in canonical_kmer_iterator
        canonical_kmer_counts[canonical_kmer] = get(canonical_kmer_counts, canonical_kmer, 0) + 1
    end
    return sort!(canonical_kmer_counts)
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
    return sort!(joint_kmer_counts)
end

function count_canonical_kmers(::Type{KMER_TYPE}, sequences::R) where {KMER_TYPE, R <: Union{FASTX.FASTA.Reader, FASTX.FASTQ.Reader}}
    joint_kmer_counts = DataStructures.OrderedDict{KMER_TYPE, Int}()
    for sequence in sequences
        sequence_kmer_counts = count_canonical_kmers(KMER_TYPE, sequence)
        merge!(+, joint_kmer_counts, sequence_kmer_counts)
    end
    return sort!(joint_kmer_counts)
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
    kmer_counts = count_canonical_kmers(KMER_TYPE, fastx_io)
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
        total_reads = count_reads(fastq_file)
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

function merge_fasta_files(list_of_fasta_files, joint_fasta_file)
    open(joint_fasta_file, "w") do io
        for f in list_of_fasta_files
            write(io, read(f))
        end
    end
    return joint_fasta_file
end