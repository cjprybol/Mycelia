"""
    per_base_quality_scores(records::AbstractVector{<:FASTX.FASTQ.Record}) -> Vector{Float64}

Return the mean Phred quality score at each base position across reads.
"""
function per_base_quality_scores(records::AbstractVector{<:FASTX.FASTQ.Record})
    isempty(records) && return Float64[]
    max_len = maximum(length(FASTX.quality_scores(r)) for r in records)
    totals = zeros(Float64, max_len)
    counts = zeros(Int, max_len)
    for record in records
        quals = FASTX.quality_scores(record)
        for (i, q) in enumerate(quals)
            totals[i] += q
            counts[i] += 1
        end
    end
    return [counts[i] == 0 ? 0.0 : totals[i] / counts[i] for i in 1:max_len]
end

"""
    per_read_quality_scores(records::AbstractVector{<:FASTX.FASTQ.Record}) -> Vector{Float64}

Return the mean Phred quality score for each read.
"""
function per_read_quality_scores(records::AbstractVector{<:FASTX.FASTQ.Record})
    return [Statistics.mean(FASTX.quality_scores(r)) for r in records]
end

"""
    gc_content_per_read(records::AbstractVector{<:FASTX.FASTQ.Record}) -> Vector{Float64}

Return GC fraction per read (0-1).
"""
function gc_content_per_read(records::AbstractVector{<:FASTX.FASTQ.Record})
    gc_vals = Float64[]
    for record in records
        seq = FASTX.sequence(record)
        gc_count = 0
        for base in seq
            if base == BioSequences.DNA_G || base == BioSequences.DNA_C ||
               base == BioSequences.RNA_G || base == BioSequences.RNA_C
                gc_count += 1
            end
        end
        total = length(seq)
        push!(gc_vals, total == 0 ? 0.0 : gc_count / total)
    end
    return gc_vals
end

"""
    read_length_distribution(records::AbstractVector{<:FASTX.FASTQ.Record}) -> Vector{Int}

Return read lengths for each record.
"""
function read_length_distribution(records::AbstractVector{<:FASTX.FASTQ.Record})
    return [length(FASTX.sequence(r)) for r in records]
end

"""
    duplication_stats(records::AbstractVector{<:FASTX.FASTQ.Record};
                      min_fraction::Float64=0.05,
                      top_n::Int=5)

Summarize sequence duplication and overrepresentation.

Returns a named tuple with:
- `unique_sequences`
- `duplicate_reads`
- `duplicate_fraction`
- `overrepresented` (vector of named tuples with sequence, count, fraction)
"""
function duplication_stats(records::AbstractVector{<:FASTX.FASTQ.Record};
        min_fraction::Float64 = 0.05,
        top_n::Int = 5)
    counts = Dict{String, Int}()
    for record in records
        seq = String(FASTX.sequence(record))
        counts[seq] = get(counts, seq, 0) + 1
    end
    total_reads = length(records)
    duplicate_reads = sum(v - 1 for v in values(counts))
    duplicate_fraction = total_reads == 0 ? 0.0 : duplicate_reads / total_reads
    overrepresented = NamedTuple{
        (:sequence, :count, :fraction), Tuple{String, Int, Float64}}[]
    for (seq, count) in counts
        if total_reads > 0 && count / total_reads >= min_fraction
            push!(overrepresented, (
                sequence = seq, count = count, fraction = count / total_reads))
        end
    end
    sort!(overrepresented, by = x -> x.count, rev = true)
    return (
        unique_sequences = length(counts),
        duplicate_reads = duplicate_reads,
        duplicate_fraction = duplicate_fraction,
        overrepresented = overrepresented[1:min(top_n, length(overrepresented))]
    )
end

"""
    summarize_fastq(label, fastq_file, genome_size)

Print basic FASTQ summary metrics and return the underlying statistics.
"""
function summarize_fastq(label, fastq_file, genome_size)
    if fastq_file === nothing || !isfile(fastq_file)
        println("  $(label): missing file, skipping")
        return nothing
    end

    quality_stats = Mycelia.analyze_fastq_quality(fastq_file)
    read_lengths = Mycelia.determine_read_lengths(fastq_file)

    min_len = isempty(read_lengths) ? 0 : minimum(read_lengths)
    max_len = isempty(read_lengths) ? 0 : maximum(read_lengths)
    median_len = isempty(read_lengths) ? 0 : Statistics.median(read_lengths)
    coverage = isempty(read_lengths) ? 0.0 : sum(read_lengths) / genome_size
    error_rate_est = Mycelia.q_value_to_error_rate(quality_stats.mean_quality)

    println("  $(label):")
    println("    reads: $(quality_stats.n_reads)")
    println("    length (min/median/max): $(min_len)/$(median_len)/$(max_len)")
    println("    mean quality: $(round(quality_stats.mean_quality, digits=2))")
    println("    GC%: $(round(quality_stats.gc_content, digits=2))")
    println("    Q30+ reads: $(round(quality_stats.quality_distribution.q30_percent, digits=2))%")
    println("    estimated coverage: $(round(coverage, digits=2))x")
    println("    estimated error rate: $(round(error_rate_est * 100, digits=3))%")

    return (; quality_stats, read_lengths, coverage, error_rate_est)
end
