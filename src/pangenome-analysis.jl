# Pangenome Analysis Module
# Functions for comparative genomics using existing k-mer infrastructure
# Includes integration with PGGB, Cactus, and vg toolkit

"""
    PangenomeAnalysisResult

Results of k-mer based pangenome analysis.
"""
struct PangenomeAnalysisResult
    genome_names::Vector{String}
    kmer_counts_by_genome::Dict{String, Dict}
    shared_kmers::Vector
    core_kmers::Vector  # Present in ALL genomes
    accessory_kmers::Vector  # Present in SOME genomes
    unique_kmers_by_genome::Dict{String, Vector}
    presence_absence_matrix::BitMatrix
    distance_matrix::Matrix{Float64}
    similarity_stats::NamedTuple
end

"""
    analyze_pangenome_kmers(genome_files::Vector{String}; kmer_type=Kmers.DNAKmer{21}, distance_metric=:jaccard)

Perform comprehensive k-mer based pangenome analysis using existing Mycelia infrastructure.

Leverages existing `count_canonical_kmers` and distance metric functions to analyze
genomic content across multiple genomes, identifying core, accessory, and unique regions.

# Arguments
- `genome_files`: Vector of FASTA file paths containing genome sequences
- `kmer_type`: K-mer type from Kmers.jl (default: `Kmers.DNAKmer{21}`)
- `distance_metric`: Distance metric (`:jaccard`, `:bray_curtis`, `:cosine`, `:js_divergence`)

# Returns
- `PangenomeAnalysisResult` with comprehensive pangenome statistics

# Example
```julia
genome_files = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
result = Mycelia.analyze_pangenome_kmers(genome_files, kmer_type=Kmers.DNAKmer{31})
println("Core k-mers: \$(length(result.core_kmers))")
println("Total pangenome size: \$(size(result.presence_absence_matrix, 1)) k-mers")
```
"""
function analyze_pangenome_kmers(genome_files::Vector{String};
        kmer_type = Kmers.DNAKmer{21}, distance_metric = :jaccard)
    if isempty(genome_files)
        error("No genome files provided")
    end

    # Validate files exist
    for file in genome_files
        if !isfile(file)
            error("Genome file does not exist: $(file)")
        end
    end

    genome_names = [basename(file) for file in genome_files]
    n_genomes = length(genome_files)

    # Count k-mers for each genome using existing infrastructure
    println("Counting k-mers for $(n_genomes) genomes...")
    kmer_counts_by_genome = Dict{String, Dict}()

    for (i, file) in enumerate(genome_files)
        println("  Processing $(genome_names[i])...")
        # Use existing count_canonical_kmers function
        kmer_counts = count_canonical_kmers(kmer_type, file)
        kmer_counts_by_genome[genome_names[i]] = kmer_counts
    end

    # Get all unique k-mers across all genomes
    all_kmers = Set()
    for kmer_counts in values(kmer_counts_by_genome)
        union!(all_kmers, keys(kmer_counts))
    end
    all_kmers = collect(all_kmers)
    n_kmers = length(all_kmers)

    println("Found $(n_kmers) unique k-mers across all genomes")

    # Build presence/absence matrix
    presence_absence_matrix = BitMatrix(undef, n_kmers, n_genomes)

    for (i, kmer) in enumerate(all_kmers)
        for (j, genome_name) in enumerate(genome_names)
            presence_absence_matrix[i, j] = haskey(kmer_counts_by_genome[genome_name], kmer)
        end
    end

    # Classify k-mers into categories
    core_kmers = []
    accessory_kmers = []
    shared_kmers = []
    unique_kmers_by_genome = Dict{String, Vector}()

    for genome_name in genome_names
        unique_kmers_by_genome[genome_name] = []
    end

    for (i, kmer) in enumerate(all_kmers)
        presence_count = sum(presence_absence_matrix[i, :])

        if presence_count == n_genomes
            # Present in all genomes - core
            push!(core_kmers, kmer)
            push!(shared_kmers, kmer)
        elseif presence_count == 1
            # Present in only one genome - unique
            genome_idx = findfirst(presence_absence_matrix[i, :])
            genome_name = genome_names[genome_idx]
            push!(unique_kmers_by_genome[genome_name], kmer)
        elseif presence_count > 1
            # Present in some genomes - accessory
            push!(accessory_kmers, kmer)
            push!(shared_kmers, kmer)
        end
    end

    # Calculate distance matrix using existing distance functions
    println("Calculating distance matrix...")
    distance_matrix = if distance_metric == :jaccard
        jaccard_distance(presence_absence_matrix)
    elseif distance_metric == :bray_curtis
        # Convert to count matrix for Bray-Curtis
        count_matrix = zeros(Int, n_kmers, n_genomes)
        for (i, kmer) in enumerate(all_kmers)
            for (j, genome_name) in enumerate(genome_names)
                count_matrix[i, j] = get(kmer_counts_by_genome[genome_name], kmer, 0)
            end
        end
        bray_curtis_distance(count_matrix)
    else
        error("Unsupported distance metric: $(distance_metric)")
    end

    # Calculate summary statistics
    core_size = length(core_kmers)
    accessory_size = length(accessory_kmers)
    pangenome_size = n_kmers
    unique_total = sum(length(kmers) for kmers in values(unique_kmers_by_genome))

    positive_distances = distance_matrix[distance_matrix .> 0]
    mean_pairwise_distance = isempty(positive_distances) ? 0.0 :
                             Statistics.mean(positive_distances)

    similarity_stats = (
        n_genomes = n_genomes,
        pangenome_size = pangenome_size,
        core_size = core_size,
        accessory_size = accessory_size,
        unique_total = unique_total,
        core_percentage = (core_size / pangenome_size) * 100,
        accessory_percentage = (accessory_size / pangenome_size) * 100,
        unique_percentage = (unique_total / pangenome_size) * 100,
        mean_pairwise_distance = mean_pairwise_distance
    )

    println("Pangenome Analysis Results:")
    println("  Total k-mers: $(pangenome_size)")
    println("  Core k-mers: $(core_size) ($(round(similarity_stats.core_percentage, digits=1))%)")
    println("  Accessory k-mers: $(accessory_size) ($(round(similarity_stats.accessory_percentage, digits=1))%)")
    println("  Unique k-mers: $(unique_total) ($(round(similarity_stats.unique_percentage, digits=1))%)")
    println("  Mean pairwise distance: $(round(similarity_stats.mean_pairwise_distance, digits=3))")

    return PangenomeAnalysisResult(
        genome_names,
        kmer_counts_by_genome,
        shared_kmers,
        core_kmers,
        accessory_kmers,
        unique_kmers_by_genome,
        presence_absence_matrix,
        distance_matrix,
        similarity_stats
    )
end

# =============================================================================
# Streaming pangenome k-mer FDR
#
# `analyze_pangenome_kmers` materializes a `BitMatrix(n_observed, n_genomes)`
# to compute presence/absence after sparse counting. That step's RAM peak
# scales with the observed-kmer set, which can balloon at high k where most
# kmers are near-singleton across distinct strains. For the common
# "what fraction of each subject genome's kmers also appear in this
# reference panel?" question the full presence/absence matrix is unnecessary
# — we only need set intersections against a single reference union.
#
# `pangenome_kmer_fdr` answers that narrower question with bounded memory:
# peak ≈ |reference_union| + max(|single subject k-mer set|).
# The dense path remains for callers that genuinely need the BitMatrix or
# unique-per-genome enumeration.
# =============================================================================

"""
    KmerMembership{K}

Abstract membership-test container for canonical k-mers.
Implementations must support:
  Base.in(kmer, m)        -> Bool
  Base.push!(m, kmer)     -> m
  Base.union!(m, kmers)   -> m
  Base.length(m)          -> Int   (cardinality; estimate for probabilistic impls)

Two concrete implementations ship today:
  - `ExactKmerMembership` (Set-backed; deterministic; default)
  - `BloomKmerMembership`  (probabilistic; opt-in for very large pangenomes;
    documented FPR adds a small upper-bound bias to reported FDRs)
"""
abstract type KmerMembership{K} end

struct ExactKmerMembership{K} <: KmerMembership{K}
    data::Set{K}
end
ExactKmerMembership{K}() where {K} = ExactKmerMembership{K}(Set{K}())
Base.in(kmer, m::ExactKmerMembership) = kmer in m.data
Base.push!(m::ExactKmerMembership, kmer) = (push!(m.data, kmer); m)
Base.union!(m::ExactKmerMembership, kmers) = (union!(m.data, kmers); m)
Base.length(m::ExactKmerMembership) = length(m.data)

"""
Probabilistic membership using k independent hash slices over a BitVector.
False-positive rate is configurable; false negatives never occur. The
caller is responsible for supplying a capacity estimate so we can size the
bit array correctly.
"""
mutable struct BloomKmerMembership{K} <: KmerMembership{K}
    bits::BitVector
    n_hashes::Int
    n_inserted::Int
    capacity::Int
    target_fpr::Float64
end
function BloomKmerMembership{K}(;
        capacity::Int, target_fpr::Float64 = 0.01) where {K}
    # Optimal Bloom sizing: m = -n*ln(p)/(ln 2)^2; k = (m/n)*ln 2.
    m = max(64, ceil(Int, -capacity * log(target_fpr) / (log(2)^2)))
    k = max(1, round(Int, (m / capacity) * log(2)))
    return BloomKmerMembership{K}(falses(m), k, 0, capacity, target_fpr)
end
function _bloom_indices(m::BloomKmerMembership, kmer)
    h1 = hash(kmer, UInt(0x9E3779B97F4A7C15))
    h2 = hash(kmer, UInt(0xBF58476D1CE4E5B9))
    n = length(m.bits)
    return ((mod(h1 + UInt(i) * h2, n) + 1) for i in 0:(m.n_hashes - 1))
end
function Base.in(kmer, m::BloomKmerMembership)
    @inbounds for idx in _bloom_indices(m, kmer)
        m.bits[idx] || return false
    end
    return true
end
function Base.push!(m::BloomKmerMembership, kmer)
    novel = false
    @inbounds for idx in _bloom_indices(m, kmer)
        if !m.bits[idx]
            ;
            m.bits[idx] = true;
            novel = true;
        end
    end
    novel && (m.n_inserted += 1)
    return m
end
Base.union!(m::BloomKmerMembership, kmers) = (foreach(k -> push!(m, k), kmers); m)
Base.length(m::BloomKmerMembership) = m.n_inserted

"""
    PangenomeKmerFDRResult

Per-subject (or per-group, when the `groups` kwarg is set) k-mer FDR vs
a reference union. `subject_names` is either genome basenames or group
labels, depending on the call.
"""
struct PangenomeKmerFDRResult
    subject_names::Vector{String}
    n_kmers::Vector{Int}
    n_shared::Vector{Int}
    fdr::Vector{Float64}
    reference_union_size::Int
    kmer_type::Type
    membership::Symbol      # :exact or :bloom (resolved value)
    cache_dir::Union{Nothing, String}
    elapsed_seconds::Float64
end

"""
    _per_genome_kmer_set(file, kmer_type, cache_dir)

Sparse per-genome canonical k-mer set extraction with optional JLD2 cache.
Cache key includes `K`, so different k values do not collide.

Shared by `analyze_pangenome_kmers` (dense path; refactor opportunity) and
`pangenome_kmer_fdr` (streaming path).
"""
function _per_genome_kmer_set(
        file::AbstractString, kmer_type::Type,
        cache_dir::Union{Nothing, AbstractString})
    K = _kmer_length(kmer_type)
    if !isnothing(cache_dir)
        Base.Filesystem.mkpath(cache_dir)
        cache_path = Base.Filesystem.joinpath(
            cache_dir, "$(Base.basename(file)).k$K.kmerset.jld2")
        if Base.Filesystem.isfile(cache_path)
            # JLD2 deserializes kmers with their concrete storage-tuple
            # parameter (`Kmer{Alphabet,K,N}`), which is a strict subtype of
            # the UnionAll `kmer_type` (`Kmer{Alphabet,K}`). `Set` is
            # invariant, so an `::Set{kmer_type}` assertion would falsely
            # reject the cached value. Return whatever JLD2 gives back —
            # downstream only iterates / membership-tests it.
            return JLD2.load(cache_path, "kmerset")
        end
        s = Set(keys(count_canonical_kmers(kmer_type, file)))
        JLD2.jldsave(cache_path; kmerset = s)
        return s
    end
    return Set(keys(count_canonical_kmers(kmer_type, file)))
end

# Extract K from `Kmers.DNAKmer{K}` (alias for `Kmers.Kmer{Alphabet,K,_}`).
# Note: `Kmers.DNAKmer{5}` is a UnionAll because the third internal type
# parameter (packed-storage tuple) is left free; unwrap to the inner body
# before scanning parameters for the Integer K.
function _kmer_length(::Type{T}) where {T}
    cur = T
    while cur isa UnionAll
        cur = cur.body
    end
    for p in cur.parameters
        p isa Integer && return Int(p)
    end
    error("cannot infer K from $T")
end

"""
    _resolve_cache(cache, subject_files, reference_files, kmer_type)

Resolve the `cache` kwarg to a concrete directory (or `nothing` for RAM).
`:auto` heuristic: if estimated peak in-RAM kmer-set size is < 50% of free
memory, use RAM; else create `mktempdir()`. A user-supplied path is used
verbatim. `:memory` forces RAM.
"""
function _resolve_cache(cache, subject_files, reference_files, kmer_type)
    cache === :memory && return nothing
    cache isa AbstractString && return String(cache)
    cache === :auto ||
        error("`cache` must be :auto, :memory, or a directory path; got $(cache)")
    # `sizeof(kmer_type)` fails on UnionAlls like `Kmers.DNAKmer{K}` (the
    # storage tuple parameter is left free), so derive from K instead:
    # 2-bit alphabet packs 4 bases/byte; +16 for Set entry bookkeeping.
    K = _kmer_length(kmer_type)
    bytes_per_kmer = max(8, cld(K, 4) + 16)
    # Estimate kmer count as ~1× FASTA bytes (canonical kmer ≈ unique within genome).
    total_bytes = sum(Base.Filesystem.filesize, vcat(subject_files, reference_files))
    estimated_peak = bytes_per_kmer * total_bytes  # very rough upper bound
    free_mem = Sys.free_memory()
    if estimated_peak < (free_mem ÷ 2)
        return nothing
    end
    return Base.Filesystem.mktempdir(prefix = "pangenome_kmer_fdr_")
end

"""
    _resolve_membership(membership_type, capacity_estimate, kmer_type, bloom_fpr)

Resolve `membership_type` to a concrete `KmerMembership` constructor.
`:auto`: exact when reference union is < 4 GB; else bloom.
`:exact`: forced exact.
`:bloom`: forced bloom (requires `capacity_estimate`).
"""
function _resolve_membership(
        membership_type, capacity_estimate, kmer_type, bloom_fpr)
    # Approximate per-kmer Set entry size in bytes. 2-bit alphabets (DNA/RNA)
    # pack 4 bases/byte; we add Set bookkeeping overhead. `sizeof(kmer_type)`
    # cannot be called on a UnionAll like `Kmers.DNAKmer{K}` (the storage
    # tuple parameter is left free), so derive from K instead.
    K = _kmer_length(kmer_type)
    bytes_per_kmer = max(8, cld(K, 4) + 16)
    if membership_type === :exact ||
       (membership_type === :auto &&
        capacity_estimate * bytes_per_kmer < 4 * 2^30)
        return (:exact, ExactKmerMembership{kmer_type}())
    elseif membership_type === :bloom || membership_type === :auto
        return (:bloom,
            BloomKmerMembership{kmer_type}(
                capacity = max(1, capacity_estimate), target_fpr = bloom_fpr))
    else
        error("`membership_type` must be :auto, :exact, or :bloom; got $(membership_type)")
    end
end

"""
    pangenome_kmer_fdr(subject_files, reference_files; kwargs...)

Streaming per-subject (or per-group) FDR vs a reference union of k-mers.
Memory bounded by `|reference union| + max(|single subject set|)` —
does NOT build the dense `(n_kmers × n_genomes)` presence/absence matrix
that `analyze_pangenome_kmers` uses.

# Arguments
- `subject_files`: vector of FASTA paths whose k-mers we are testing.
- `reference_files`: vector of FASTA paths whose union forms the
  reference set the subjects are compared against.

# Keyword arguments
- `kmer_type = Kmers.DNAKmer{21}`: canonical k-mer type.
- `groups = nothing`: when a `Dict{filename, label}` is supplied, subjects
  sharing a label are unioned and a single FDR row is emitted per label
  (per-cluster pangenome FDR). When `nothing`, one row per subject.
- `cache = :auto`: `:auto` picks RAM vs disk based on free memory;
  `:memory` forces RAM; a path is used as a per-genome JLD2 cache dir
  (cache keys include K so multiple K runs share extracted sets).
- `membership_type = :auto`: `:auto` picks exact for < 4 GB unions, else
  bloom; `:exact` forces deterministic Set-backed; `:bloom` opts into
  probabilistic membership for very large pangenomes.
- `bloom_fpr = 0.01`: Bloom target false-positive rate (only used when
  bloom is selected). Reported FDR is then a slight upper bound by
  approximately `bloom_fpr` for kmers genuinely absent from the reference.
- `progress = true`: print per-genome progress lines.
- `reproducibility_check = false`: when `true` and ≤ 10 subjects, runs
  both this function and `analyze_pangenome_kmers` and asserts agreement
  to 1e-12. Off by default to keep production runs fast.

# Example
```julia
result = Mycelia.pangenome_kmer_fdr(
    cstrain_fastas, public_ecoli_fastas;
    kmer_type = Kmers.DNAKmer{31},
    cache = mktempdir(),
)
println("Per-strain FDR median: ", Statistics.median(result.fdr))
```
"""
function pangenome_kmer_fdr(
        subject_files::AbstractVector{<:AbstractString},
        reference_files::AbstractVector{<:AbstractString};
        kmer_type::Type = Kmers.DNAKmer{21},
        groups::Union{Nothing, AbstractDict} = nothing,
        cache = :auto,
        membership_type::Symbol = :auto,
        bloom_fpr::Float64 = 0.01,
        progress::Bool = true,
        reproducibility_check::Bool = false
)::PangenomeKmerFDRResult
    isempty(subject_files) && error("no subject_files supplied")
    isempty(reference_files) && error("no reference_files supplied")
    for f in vcat(subject_files, reference_files)
        Base.Filesystem.isfile(f) || error("file not found: $f")
    end

    t0 = time()
    cache_dir = _resolve_cache(cache, subject_files, reference_files, kmer_type)
    progress &&
        println("pangenome_kmer_fdr: cache_dir = ", isnothing(cache_dir) ? ":memory" :
                                                    cache_dir)

    # 1) Build reference union (always exact during build; switch storage
    #    type only after we know its capacity).
    progress &&
        println("Building reference union over $(length(reference_files)) genome(s)...")
    ref_union_set = Set{kmer_type}()
    for (i, ref) in enumerate(reference_files)
        s = _per_genome_kmer_set(ref, kmer_type, cache_dir)
        union!(ref_union_set, s)
        progress && (i % 50 == 0 || i == length(reference_files)) &&
            println("  ref $i/$(length(reference_files))  |union| = $(length(ref_union_set))")
    end
    reference_union_size = length(ref_union_set)

    membership_kind,
    membership = _resolve_membership(
        membership_type, reference_union_size, kmer_type, bloom_fpr)
    if membership_kind === :exact
        # Reuse the set we already built; no extra allocation.
        membership = ExactKmerMembership{kmer_type}(ref_union_set)
    else
        union!(membership, ref_union_set)
        ref_union_set = Set{kmer_type}()  # release backing memory
    end
    progress &&
        println("Reference membership: $(membership_kind), |U| = $(length(membership))")

    # 2) Stream subjects (or grouped subjects) and tally intersections.
    if groups === nothing
        subject_names = [Base.basename(f) for f in subject_files]
        n_kmers = Vector{Int}(undef, length(subject_files))
        n_shared = Vector{Int}(undef, length(subject_files))
        for (i, file) in enumerate(subject_files)
            s = _per_genome_kmer_set(file, kmer_type, cache_dir)
            n_kmers[i] = length(s)
            n_shared[i] = count(in(membership), s)
            progress && (i % 50 == 0 || i == length(subject_files)) &&
                println("  subj $i/$(length(subject_files))  |s|=$(n_kmers[i])  shared=$(n_shared[i])")
        end
    else
        # Group → union of member kmer sets → FDR vs reference union.
        # Subjects keyed by basename so callers can build groups against
        # full paths or basenames interchangeably.
        labels = sort(unique(values(groups)))
        subject_names = String.(labels)
        n_kmers = zeros(Int, length(labels))
        n_shared = zeros(Int, length(labels))
        label_to_idx = Dict(lab => i for (i, lab) in enumerate(labels))
        # Build per-label kmer unions one label at a time to bound memory.
        for (li, lab) in enumerate(labels)
            members = filter(
                f -> get(groups, f, get(groups, Base.basename(f), nothing)) == lab,
                subject_files)
            grp_set = Set{kmer_type}()
            for f in members
                union!(grp_set, _per_genome_kmer_set(f, kmer_type, cache_dir))
            end
            n_kmers[li] = length(grp_set)
            n_shared[li] = count(in(membership), grp_set)
            progress &&
                println("  group $li/$(length(labels)) ($(lab)): |U|=$(n_kmers[li])  shared=$(n_shared[li])")
        end
    end

    fdr = [n_kmers[i] == 0 ? 0.0 : n_shared[i] / n_kmers[i]
           for i in eachindex(n_kmers)]

    elapsed = time() - t0
    result = PangenomeKmerFDRResult(
        subject_names, n_kmers, n_shared, fdr,
        reference_union_size, kmer_type, membership_kind,
        cache_dir, elapsed
    )

    if reproducibility_check && length(subject_files) <= 10 && groups === nothing
        progress &&
            println("reproducibility_check: comparing against analyze_pangenome_kmers...")
        dense = analyze_pangenome_kmers(
            vcat(subject_files, reference_files); kmer_type = kmer_type)
        ref_set = Set{kmer_type}()
        for ref in reference_files
            union!(ref_set, keys(dense.kmer_counts_by_genome[Base.basename(ref)]))
        end
        for (i, file) in enumerate(subject_files)
            kset = Set(keys(dense.kmer_counts_by_genome[Base.basename(file)]))
            expected = isempty(kset) ? 0.0 :
                       count(k -> k in ref_set, kset) / length(kset)
            @assert isapprox(result.fdr[i], expected; atol = 1e-12) (
                "reproducibility check failed for $(Base.basename(file)): " *
                "streaming=$(result.fdr[i]) dense=$expected")
        end
        progress && println("reproducibility_check: PASS")
    end

    return result
end

"""
    compare_genome_kmer_similarity(genome1_file::String, genome2_file::String; kmer_type=Kmers.DNAKmer{21}, metric=:js_divergence)

Compare two genomes using existing k-mer distance metrics.

Leverages existing distance metric functions to compare genomic similarity
between pairs of genomes using various distance measures.

# Arguments
- `genome1_file`: Path to first genome FASTA file
- `genome2_file`: Path to second genome FASTA file  
- `kmer_type`: K-mer type from Kmers.jl (default: `Kmers.DNAKmer{21}`)
- `metric`: Distance metric (`:js_divergence`, `:cosine`, `:jaccard`)

# Returns
- Named tuple with similarity/distance metrics and k-mer statistics

# Example
```julia
similarity = Mycelia.compare_genome_kmer_similarity(
    "genome1.fasta", "genome2.fasta", 
    kmer_type=Kmers.DNAKmer{31}, 
    metric=:js_divergence
)
println("JS divergence: \$(similarity.distance)")
println("Shared k-mers: \$(similarity.shared_kmers)")
```
"""
function compare_genome_kmer_similarity(genome1_file::String, genome2_file::String;
        kmer_type = Kmers.DNAKmer{21}, metric = :js_divergence)
    # Count k-mers using existing infrastructure
    kmer_counts1 = count_canonical_kmers(kmer_type, genome1_file)
    kmer_counts2 = count_canonical_kmers(kmer_type, genome2_file)

    # Calculate distance using existing functions
    distance = if metric == :js_divergence
        kmer_counts_to_js_divergence(kmer_counts1, kmer_counts2)
    elseif metric == :cosine
        kmer_counts_to_cosine_similarity(kmer_counts1, kmer_counts2)
    elseif metric == :jaccard
        # Calculate Jaccard manually for k-mer dicts
        shared_kmers = intersect(keys(kmer_counts1), keys(kmer_counts2))
        union_kmers = union(keys(kmer_counts1), keys(kmer_counts2))
        1.0 - length(shared_kmers) / length(union_kmers)
    else
        error("Unsupported metric: $(metric)")
    end

    # Calculate additional statistics
    shared_kmers = length(intersect(keys(kmer_counts1), keys(kmer_counts2)))
    total_kmers = length(union(keys(kmer_counts1), keys(kmer_counts2)))
    genome1_unique = length(setdiff(keys(kmer_counts1), keys(kmer_counts2)))
    genome2_unique = length(setdiff(keys(kmer_counts2), keys(kmer_counts1)))

    return (
        distance = distance,
        metric = metric,
        shared_kmers = shared_kmers,
        total_kmers = total_kmers,
        genome1_unique = genome1_unique,
        genome2_unique = genome2_unique,
        jaccard_similarity = shared_kmers / total_kmers
    )
end

"""
    build_genome_distance_matrix(genome_files::Vector{String}; kmer_type=Kmers.DNAKmer{21}, metric=:js_divergence)

Build a distance matrix between all genome pairs using existing distance metrics.

Creates a comprehensive pairwise distance matrix using established k-mer
distance functions, suitable for phylogenetic analysis and clustering.

# Arguments
- `genome_files`: Vector of genome FASTA file paths
- `kmer_type`: K-mer type from Kmers.jl (default: `Kmers.DNAKmer{21}`)
- `metric`: Distance metric (`:js_divergence`, `:cosine`, `:jaccard`)

# Returns
- Named tuple with distance matrix and genome names

# Example
```julia
genomes = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
result = Mycelia.build_genome_distance_matrix(genomes, kmer_type=Kmers.DNAKmer{31})
println("Distance matrix: \$(result.distance_matrix)")
```
"""
function build_genome_distance_matrix(genome_files::Vector{String};
        kmer_type = Kmers.DNAKmer{21}, metric = :js_divergence)
    n_genomes = length(genome_files)
    distance_matrix = zeros(Float64, n_genomes, n_genomes)
    genome_names = [basename(file) for file in genome_files]

    println("Building $(n_genomes)x$(n_genomes) distance matrix using $(metric)...")

    # Calculate pairwise distances
    for i in 1:n_genomes
        for j in (i + 1):n_genomes
            println("  Comparing $(genome_names[i]) vs $(genome_names[j])...")
            similarity = compare_genome_kmer_similarity(
                genome_files[i], genome_files[j];
                kmer_type = kmer_type, metric = metric
            )
            distance_matrix[i, j] = similarity.distance
            distance_matrix[j, i] = similarity.distance  # Symmetric
        end
    end

    return (
        distance_matrix = distance_matrix,
        genome_names = genome_names,
        metric = metric
    )
end

# PGGB Integration Functions

"""
    construct_pangenome_pggb(genome_files::Vector{String}, output_dir::String; 
                            threads::Int=2, segment_length::Int=5000, 
                            block_length::Int=3*segment_length, 
                            mash_kmer::Int=16, min_match_length::Int=19,
                            transclose_batch::Int=10000000, 
                            additional_args::Vector{String}=String[])

Construct a pangenome using PGGB (PanGenome Graph Builder).

Uses the PGGB tool to build pangenome graphs from multiple genome assemblies,
following the workflows established in the Mycelia-Dev benchmarking notebooks.

# Arguments
- `genome_files`: Vector of FASTA file paths to include in pangenome
- `output_dir`: Directory for PGGB output files
- `threads`: Number of threads for parallel processing (default: 2)
- `segment_length`: Segment length for mapping (default: 5000)
- `block_length`: Block length for mapping (default: 3*segment_length)
- `mash_kmer`: Kmer size for mash sketching (default: 16)
- `min_match_length`: Minimum match length (default: 19)
- `transclose_batch`: Batch size for transitive closure (default: 10000000)
- `additional_args`: Additional command line arguments

# Returns
- Path to the main GFA output file

# Example
```julia
genomes = ["reference.fasta", "assembly1.fasta", "assembly2.fasta"]
gfa_file = Mycelia.construct_pangenome_pggb(genomes, "pangenome_output")
```
"""
function construct_pangenome_pggb(genome_files::Vector{String}, output_dir::String;
        threads::Int = 2, segment_length::Int = 5000,
        block_length::Int = 3*segment_length,
        mash_kmer::Int = 16, min_match_length::Int = 19,
        transclose_batch::Int = 10000000,
        additional_args::Vector{String} = String[],
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "pggb",
        time_limit::String = "4-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)

    # Validate input files
    for file in genome_files
        if !isfile(file)
            error("Genome file does not exist: $(file)")
        end
    end

    # Create joint FASTA file for PGGB input
    joint_fasta = joinpath(output_dir, "joint_genomes.fasta")
    mkpath(dirname(joint_fasta))

    # Build PGGB command
    pggb_args = [
        "-i", joint_fasta,
        "-o", output_dir,
        "-t", string(threads),
        "-s", string(segment_length),
        "-l", string(block_length),
        "-k", string(mash_kmer),
        "-G", string(min_match_length),
        "-B", string(transclose_batch),
        "-n", string(length(genome_files))
    ]

    # Add any additional arguments
    append!(pggb_args, additional_args)
    cmd = `$(CONDA_RUNNER) run --live-stream -n pggb pggb $(pggb_args)`

    if executor !== nothing
        add_bioconda_env("pggb")
        quoted_inputs = join(["\"$(file)\"" for file in genome_files], " ")
        script = join(
            [
                "set -euo pipefail",
                "mkdir -p \"$(output_dir)\"",
                "cat $(quoted_inputs) > \"$(joint_fasta)\"",
                "if [ ! -f \"$(joint_fasta).fai\" ]; then",
                "  $(CONDA_RUNNER) run --live-stream -n pggb samtools faidx \"$(joint_fasta)\"",
                "fi",
                Mycelia.command_string(cmd)
            ],
            "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return joinpath(output_dir, "pangenome.gfa")
    end

    # Concatenate genome files (don't merge to preserve identifiers)
    concatenate_files(files = genome_files, file = joint_fasta)

    # Index the joint FASTA if not already indexed
    if !isfile(joint_fasta * ".fai")
        samtools_index_fasta(fasta = joint_fasta)
    end

    # Ensure PGGB conda environment exists
    add_bioconda_env("pggb")

    # Run PGGB
    println("Running PGGB pangenome construction with $(length(genome_files)) genomes...")

    try
        run(cmd)
    catch e
        error("PGGB failed: $(e)")
    end

    # Find and return the main GFA output file
    gfa_files = filter(x -> endswith(x, ".gfa"), readdir(output_dir, join = true))
    if isempty(gfa_files)
        error("No GFA output file found in $(output_dir)")
    end

    main_gfa = first(sort(gfa_files, by = filesize, rev = true))  # Return largest GFA file
    println("PGGB completed. Main GFA file: $(main_gfa)")

    return main_gfa
end

"""
    call_variants_from_pggb_graph(gfa_file::String, reference_prefix::String; 
                                 threads::Int=2, ploidy::Int=1, output_file::String="")

Call variants from a PGGB-generated pangenome graph using vg deconstruct.

Uses the vg toolkit to extract variants from pangenome graphs, following
the methodology established in Mycelia-Dev benchmarking workflows.

# Arguments
- `gfa_file`: Path to GFA pangenome graph file from PGGB
- `reference_prefix`: Prefix for reference paths in the graph
- `threads`: Number of threads (default: 2)
- `ploidy`: Ploidy level (default: 1)
- `output_file`: Output VCF file path (default: gfa_file + ".vcf")

# Returns
- Path to output VCF file

# Example
```julia
gfa_file = "pangenome.gfa"
vcf_file = Mycelia.call_variants_from_pggb_graph(gfa_file, "reference")
```
"""
function call_variants_from_pggb_graph(gfa_file::String, reference_prefix::String;
        threads::Int = 2, ploidy::Int = 1, output_file::String = "")
    if !isfile(gfa_file)
        error("GFA file does not exist: $(gfa_file)")
    end

    # Set default output file name
    if isempty(output_file)
        output_file = gfa_file * ".vcf"
    end

    # Ensure vg conda environment exists
    add_bioconda_env("vg")

    # Build vg deconstruct command
    vg_args = [
        "deconstruct",
        "--path-prefix", reference_prefix,
        "--ploidy", string(ploidy),
        "--path-traversals",
        "--all-snarls",
        "--threads", string(threads),
        gfa_file
    ]

    println("Calling variants from pangenome graph: $(gfa_file)")
    cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg $(vg_args)`

    try
        run(pipeline(cmd, output_file))
    catch e
        error("vg deconstruct failed: $(e)")
    end

    println("Variant calling completed. VCF file: $(output_file)")

    return output_file
end

# Cactus Integration Functions

"""
    construct_pangenome_cactus(genome_files::Vector{String}, genome_names::Vector{String}, 
                              output_dir::String, reference_name::String;
                              max_cores::Int=8, max_memory_gb::Int=32,
                              output_formats::Vector{String}=["gbz", "gfa", "vcf", "odgi"])

Construct a pangenome using Cactus alignment-based approach.

Uses the Cactus pangenome pipeline via containerized execution to build
pangenome graphs from multiple genome assemblies using progressive alignment.

# Arguments
- `genome_files`: Vector of FASTA file paths
- `genome_names`: Vector of sample names corresponding to genome files
- `output_dir`: Directory for Cactus output
- `reference_name`: Name of the reference sample for pangenome construction
- `max_cores`: Maximum number of CPU cores (default: 8)
- `max_memory_gb`: Maximum memory in GB (default: 32)
- `output_formats`: Output formats to generate (default: ["gbz", "gfa", "vcf", "odgi"])

# Returns
- Dictionary with paths to output files for each format

# Example
```julia
genomes = ["ref.fasta", "asm1.fasta", "asm2.fasta"]
names = ["REFERENCE", "SAMPLE1", "SAMPLE2"]
outputs = Mycelia.construct_pangenome_cactus(genomes, names, "cactus_out", "REFERENCE")
```
"""
function construct_pangenome_cactus(
        genome_files::Vector{String}, genome_names::Vector{String},
        output_dir::String, reference_name::String;
        max_cores::Int = 8, max_memory_gb::Int = 32,
        output_formats::Vector{String} = ["gbz", "gfa", "vcf", "odgi"],
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "cactus_pangenome",
        time_limit::String = "8-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)

    # Validate inputs
    if length(genome_files) != length(genome_names)
        error("Number of genome files must match number of genome names")
    end

    for file in genome_files
        if !isfile(file)
            error("Genome file does not exist: $(file)")
        end
    end

    if !(reference_name in genome_names)
        error("Reference name '$(reference_name)' not found in genome names")
    end

    # Create output directory
    mkpath(output_dir)

    # Create Cactus configuration file
    config_table = DataFrames.DataFrame(
        samples = genome_names,
        file_paths = genome_files
    )

    config_file = joinpath(output_dir, "cactus_config.txt")
    uCSV.write(config_file, data = collect(DataFrames.eachcol(config_table)),
        header = missing, delim = '\t')

    # Set up job store and output names
    jobstore = joinpath(output_dir, "cactus-job-store")
    output_name = "pangenome"

    # Prepare output format flags
    format_flags = String[]
    for format in output_formats
        if format in ["gbz", "gfa", "vcf", "odgi"]
            push!(format_flags, "--$(format)")
        else
            @warn "Unknown output format: $(format)"
        end
    end

    # Build Cactus command using podman-hpc
    cactus_args = [
        "run", "-it",
        "-v", "$(output_dir):/app",
        "-w", "/app",
        "quay.io/comparative-genomics-toolkit/cactus:v2.8.1",
        "cactus-pangenome",
        "./$(basename(jobstore))",
        "./$(basename(config_file))",
        "--maxCores", string(max_cores),
        "--maxMemory", "$(max_memory_gb)Gb",
        "--outDir", ".",
        "--outName", output_name,
        "--reference", reference_name
    ]

    # Add format flags
    append!(cactus_args, format_flags)

    println("Running Cactus pangenome construction...")
    println("  Genomes: $(length(genome_files))")
    println("  Reference: $(reference_name)")
    println("  Output formats: $(join(output_formats, ", "))")

    cmd = `podman-hpc $(cactus_args)`

    # Set up logging
    log_file = joinpath(output_dir, "cactus.log")

    if executor !== nothing
        script = join(
            [
                "set -euo pipefail",
                "mkdir -p \"$(output_dir)\"",
                "$(Mycelia.command_string(cmd)) > \"$(log_file)\" 2>&1"
            ],
            "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = max_cores,
            mem_gb = isnothing(mem_gb) ? max_memory_gb : mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        output_files = Dict{String, String}()
        for format in output_formats
            output_files[format] = joinpath(output_dir, "$(output_name).$(format)")
        end
        return output_files
    end

    try
        run(pipeline(cmd, stdout = log_file, stderr = log_file))
    catch e
        @warn "Cactus may have failed. Check log file: $(log_file)"
        error("Cactus pangenome construction failed: $(e)")
    end

    # Find output files
    output_files = Dict{String, String}()
    for format in output_formats
        pattern = "$(output_name).$(format)"
        files = filter(x -> occursin(pattern, x), readdir(output_dir, join = true))
        if !isempty(files)
            output_files[format] = first(files)
        else
            @warn "Output file for format $(format) not found"
        end
    end

    println("Cactus pangenome construction completed.")
    println("Output files: $(length(output_files))")

    return output_files
end

# vg Toolkit Integration Functions

"""
    convert_gfa_to_vg_format(gfa_file::String; output_file::String="")

Convert a GFA pangenome graph to vg native format.

Converts GFA format pangenome graphs to vg's optimized binary format
for improved performance in downstream analysis.

# Arguments
- `gfa_file`: Input GFA file path
- `output_file`: Output vg file path (default: gfa_file with .vg extension)

# Returns
- Path to output vg file

# Example
```julia
vg_file = Mycelia.convert_gfa_to_vg_format("pangenome.gfa")
```
"""
function convert_gfa_to_vg_format(gfa_file::String; output_file::String = "")
    if !isfile(gfa_file)
        error("GFA file does not exist: $(gfa_file)")
    end

    if isempty(output_file)
        output_file = replace(gfa_file, r"\.gfa$" => ".vg")
    end

    # Ensure vg conda environment exists
    add_bioconda_env("vg")

    println("Converting GFA to vg format: $(gfa_file)")
    cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg convert -g $(gfa_file) -v`

    try
        run(pipeline(cmd, output_file))
    catch e
        error("GFA to vg conversion failed: $(e)")
    end

    println("Conversion completed. vg file: $(output_file)")

    return output_file
end

"""
    index_pangenome_graph(graph_file::String; index_types::Vector{String}=["xg", "gcsa"])

Create indexes for a pangenome graph to enable efficient querying and mapping.

Builds various index types for pangenome graphs to support different
analysis workflows including read mapping and path queries.

# Arguments
- `graph_file`: Input graph file (GFA or vg format)
- `index_types`: Types of indexes to build (default: ["xg", "gcsa"])

# Returns
- Dictionary mapping index types to their file paths

# Example
```julia
indexes = Mycelia.index_pangenome_graph("pangenome.vg", index_types=["xg", "gcsa", "snarls"])
```
"""
function index_pangenome_graph(graph_file::String; index_types::Vector{String} = [
        "xg", "gcsa"])
    if !isfile(graph_file)
        error("Graph file does not exist: $(graph_file)")
    end

    # Ensure vg conda environment exists
    add_bioconda_env("vg")

    index_files = Dict{String, String}()
    base_name = replace(graph_file, r"\.(gfa|vg)$" => "")

    for index_type in index_types
        index_file = "$(base_name).$(index_type)"

        if index_type == "xg"
            cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg index -x $(index_file) $(graph_file)`
        elseif index_type == "gcsa"
            # GCSA indexing requires pruning first
            pruned_file = "$(base_name).pruned.vg"
            prune_cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg prune $(graph_file)`
            run(pipeline(prune_cmd, pruned_file))

            cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg index -g $(index_file) $(pruned_file)`
        elseif index_type == "snarls"
            cmd = `$(CONDA_RUNNER) run --live-stream -n vg vg snarls $(graph_file)`
        else
            @warn "Unknown index type: $(index_type)"
            continue
        end

        println("Building $(index_type) index for $(graph_file)")

        try
            if index_type == "snarls"
                run(pipeline(cmd, index_file))
            else
                run(cmd)
            end
            index_files[index_type] = index_file
            println("Index completed: $(index_file)")
        catch e
            @warn "Failed to build $(index_type) index: $(e)"
        end
    end

    return index_files
end
