# Tests for `Mycelia.pangenome_kmer_fdr` — the streaming pangenome FDR
# primitive. Run with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/pangenome_kmer_fdr.jl")'
# ```
#
# Or as part of the full suite:
#
# ```bash
# julia --project=. -e 'import Pkg; Pkg.test()'
# ```

import Test
import Mycelia
import FASTX
import Kmers
import BioSequences
import Statistics
import StableRNGs
import JLD2

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""Write a single-record FASTA with the given DNA sequence and return its path."""
function _write_fasta(dir, name, seq)
    path = joinpath(dir, name)
    open(path, "w") do io
        writer = FASTX.FASTA.Writer(io)
        write(writer, FASTX.FASTA.Record(name, BioSequences.LongDNA{4}(seq)))
        close(writer)
    end
    return path
end

"""Hand-compute expected FDR for a subject's canonical kmers vs a reference union."""
function _expected_fdr(subject_path, reference_paths, kmer_type)
    subj_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, subject_path)))
    ref_union = Set{eltype(subj_set)}()
    for ref in reference_paths
        union!(ref_union, keys(Mycelia.count_canonical_kmers(kmer_type, ref)))
    end
    n_total = length(subj_set)
    n_total == 0 && return (n_total = 0, n_shared = 0, fdr = 0.0)
    n_shared = count(in(ref_union), subj_set)
    return (n_total = n_total, n_shared = n_shared, fdr = n_shared / n_total)
end

# ---------------------------------------------------------------------------
# Per-genome FDR — golden-value test against hand-computed expected
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — per-genome FDR" begin
    temp = mktempdir()
    try
        rng = StableRNGs.StableRNG(0x4d59c3)
        # Three subjects + one reference. Construct so each subject has a
        # known, distinguishable FDR vs the reference.
        subj_a = _write_fasta(temp, "subj_a.fasta", "ACGTACGTACGTACGTACGT")
        subj_b = _write_fasta(temp, "subj_b.fasta", "ACGTACGTGGGGGGGGGGGG")
        subj_c = _write_fasta(temp, "subj_c.fasta", "TTTTTTTTTTTTTTTTTTTT")
        ref_x = _write_fasta(temp, "ref_x.fasta", "ACGTACGTACGTAAAAAAAA")

        kmer_type = Kmers.DNAKmer{5}
        result = Mycelia.pangenome_kmer_fdr(
            [subj_a, subj_b, subj_c], [ref_x];
            kmer_type = kmer_type,
            cache = :memory,
            progress = false
        )

        Test.@test result.subject_names ==
                   [basename(subj_a), basename(subj_b), basename(subj_c)]
        Test.@test result.kmer_type == kmer_type
        Test.@test result.membership === :exact
        Test.@test result.cache_dir === nothing
        Test.@test length(result.fdr) == 3

        for (i, subj) in enumerate([subj_a, subj_b, subj_c])
            exp = _expected_fdr(subj, [ref_x], kmer_type)
            Test.@test result.n_kmers[i] == exp.n_total
            Test.@test result.n_shared[i] == exp.n_shared
            Test.@test isapprox(result.fdr[i], exp.fdr; atol = 1e-12)
        end

        # Reference union size matches hand count.
        ref_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, ref_x)))
        Test.@test result.reference_union_size == length(ref_set)
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# `reproducibility_check` — streaming must match analyze_pangenome_kmers
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — reproducibility vs dense path" begin
    temp = mktempdir()
    try
        # 4 subject + 3 reference panel; small enough for the dense path.
        seqs = [
            ("g1.fasta", "ACGTACGTACGTACGTACGT"),
            ("g2.fasta", "ACGTAAGGCCCCCCCCCCCC"),
            ("g3.fasta", "TGCATGCATGCATGCATGCA"),
            ("g4.fasta", "ATATATATATATATATATAT"),
            ("r1.fasta", "ACGTACGTACGTGGGGGGGG"),
            ("r2.fasta", "TTTTTTTTTTTTTTTTTTTT"),
            ("r3.fasta", "GCGCGCGCGCGCGCGCGCGC")
        ]
        paths = [_write_fasta(temp, name, seq) for (name, seq) in seqs]
        subjects = paths[1:4]
        refs = paths[5:7]

        result = Mycelia.pangenome_kmer_fdr(
            subjects, refs;
            kmer_type = Kmers.DNAKmer{4},
            cache = :memory,
            progress = false,
            reproducibility_check = true  # asserts agreement to 1e-12 internally
        )
        Test.@test length(result.fdr) == 4
        # If reproducibility_check passed, the @assert inside the function
        # didn't fire. Re-check explicitly here as belt + suspenders.
        for (i, subj) in enumerate(subjects)
            exp = _expected_fdr(subj, refs, Kmers.DNAKmer{4})
            Test.@test isapprox(result.fdr[i], exp.fdr; atol = 1e-12)
        end
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# `groups` kwarg — per-group union FDR
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — groups (per-cluster union)" begin
    temp = mktempdir()
    try
        s1 = _write_fasta(temp, "s1.fasta", "ACGTACGTACGTACGT")
        s2 = _write_fasta(temp, "s2.fasta", "GGGGGGGGGGGGGGGG")
        s3 = _write_fasta(temp, "s3.fasta", "TTTTTTTTTTTTTTTT")
        r1 = _write_fasta(temp, "r1.fasta", "ACGTACGTAAAAAAAA")

        kmer_type = Kmers.DNAKmer{4}
        groups = Dict(s1 => "g1", s2 => "g1", s3 => "g2")
        result = Mycelia.pangenome_kmer_fdr(
            [s1, s2, s3], [r1];
            kmer_type = kmer_type,
            groups = groups,
            cache = :memory,
            progress = false
        )

        Test.@test sort(result.subject_names) == ["g1", "g2"]

        # g1's combined kmer set must equal union of s1 and s2 sets.
        s1_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, s1)))
        s2_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, s2)))
        s3_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, s3)))
        ref_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, r1)))

        idx_g1 = findfirst(==("g1"), result.subject_names)
        idx_g2 = findfirst(==("g2"), result.subject_names)
        g1_union = union(s1_set, s2_set)
        g2_union = s3_set

        Test.@test result.n_kmers[idx_g1] == length(g1_union)
        Test.@test result.n_kmers[idx_g2] == length(g2_union)
        Test.@test result.n_shared[idx_g1] == count(in(ref_set), g1_union)
        Test.@test result.n_shared[idx_g2] == count(in(ref_set), g2_union)
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# Cache modes — :memory, explicit path, and round-trip across calls
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — cache modes" begin
    temp = mktempdir()
    try
        s1 = _write_fasta(temp, "s1.fasta", "ACGTACGTACGTACGT")
        r1 = _write_fasta(temp, "r1.fasta", "ACGTAAAAAAAAAAAA")
        kmer_type = Kmers.DNAKmer{5}

        # :memory
        rmem = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type, cache = :memory, progress = false)
        Test.@test rmem.cache_dir === nothing

        # explicit path → caches written, K-keyed
        cache_dir = joinpath(temp, "kmer_cache")
        rdisk = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type, cache = cache_dir, progress = false)
        Test.@test rdisk.cache_dir == cache_dir
        Test.@test isfile(joinpath(cache_dir, "s1.fasta.k5.kmerset.jld2"))
        Test.@test isfile(joinpath(cache_dir, "r1.fasta.k5.kmerset.jld2"))

        # second call against same cache must produce identical result.
        rdisk2 = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type, cache = cache_dir, progress = false)
        Test.@test rdisk2.fdr == rdisk.fdr
        Test.@test rdisk2.n_kmers == rdisk.n_kmers

        # different K must use a different cache key (no collision).
        Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = Kmers.DNAKmer{4}, cache = cache_dir,
            progress = false)
        Test.@test isfile(joinpath(cache_dir, "s1.fasta.k4.kmerset.jld2"))
        Test.@test isfile(joinpath(cache_dir, "s1.fasta.k5.kmerset.jld2"))

        # :auto resolves to either :memory or a tempdir without erroring.
        rauto = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type, cache = :auto, progress = false)
        Test.@test rauto.fdr == rmem.fdr
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# Membership types — :exact, :bloom, :auto
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — membership types" begin
    temp = mktempdir()
    try
        s1 = _write_fasta(temp, "s1.fasta",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
        r1 = _write_fasta(temp, "r1.fasta",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
        kmer_type = Kmers.DNAKmer{8}

        # Exact: deterministic; for these disjoint inputs FDR must be 0.
        rex = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type,
            membership_type = :exact, cache = :memory, progress = false)
        Test.@test rex.membership === :exact
        Test.@test rex.fdr[1] == 0.0

        # Auto: same panel, must pick exact (small union) and match.
        rauto = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type,
            membership_type = :auto, cache = :memory, progress = false)
        Test.@test rauto.membership === :exact
        Test.@test rauto.fdr == rex.fdr

        # Bloom: forced. Reported FDR is an upper bound (false positives
        # never reduce shared count). With true FDR=0, bloom FDR is
        # ≤ bloom_fpr × |subject| / |subject| = bloom_fpr in expectation.
        # No-false-negatives invariant: bloom_fdr >= exact_fdr always.
        rbloom = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = kmer_type,
            membership_type = :bloom, bloom_fpr = 0.01,
            cache = :memory, progress = false)
        Test.@test rbloom.membership === :bloom
        Test.@test rbloom.fdr[1] >= rex.fdr[1] - 1e-12
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# Edge cases — input validation, identical genomes, empty groups
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — edge cases" begin
    temp = mktempdir()
    try
        s1 = _write_fasta(temp, "s1.fasta", "ACGTACGTACGT")
        r1 = _write_fasta(temp, "r1.fasta", "ACGTACGTACGT")  # identical to s1

        # Identical genomes: every subject kmer is in the reference union.
        result = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = Kmers.DNAKmer{4},
            cache = :memory, progress = false)
        Test.@test result.fdr[1] == 1.0
        Test.@test result.n_shared[1] == result.n_kmers[1]

        # Empty subject_files / reference_files must error.
        Test.@test_throws ErrorException Mycelia.pangenome_kmer_fdr(
            String[], [r1]; kmer_type = Kmers.DNAKmer{4}, progress = false)
        Test.@test_throws ErrorException Mycelia.pangenome_kmer_fdr(
            [s1], String[]; kmer_type = Kmers.DNAKmer{4}, progress = false)

        # Missing file must error.
        Test.@test_throws ErrorException Mycelia.pangenome_kmer_fdr(
            [s1], [joinpath(temp, "nonexistent.fasta")];
            kmer_type = Kmers.DNAKmer{4}, progress = false)

        # Bad membership_type must error.
        Test.@test_throws ErrorException Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = Kmers.DNAKmer{4},
            membership_type = :bogus, progress = false)

        # Bad cache value must error.
        Test.@test_throws ErrorException Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = Kmers.DNAKmer{4},
            cache = 42, progress = false)

        # Single-member group is allowed.
        single_group = Dict(s1 => "only_one")
        rg = Mycelia.pangenome_kmer_fdr(
            [s1], [r1]; kmer_type = Kmers.DNAKmer{4},
            groups = single_group, cache = :memory, progress = false)
        Test.@test rg.subject_names == ["only_one"]
        Test.@test rg.fdr[1] == 1.0
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# K extraction — _kmer_length must work across K values
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — K extraction across kmer types" begin
    temp = mktempdir()
    try
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        s1 = _write_fasta(temp, "s1.fasta", seq)
        r1 = _write_fasta(temp, "r1.fasta", seq)

        for K in (3, 5, 11, 21)
            kmer_type = Kmers.DNAKmer{K}
            result = Mycelia.pangenome_kmer_fdr(
                [s1], [r1]; kmer_type = kmer_type,
                cache = :memory, progress = false)
            Test.@test result.kmer_type == kmer_type
            Test.@test result.fdr[1] == 1.0  # identical inputs
        end
    finally
        rm(temp; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# KmerMembership interface — Exact + Bloom satisfy the same contract
# ---------------------------------------------------------------------------

Test.@testset "KmerMembership interface contract" begin
    K = Kmers.DNAKmer{4}
    # Get a kmer instance via count_canonical_kmers (the canonical Mycelia
    # entry point) so the test doesn't depend on Kmers.jl internal API.
    temp = mktempdir()
    fpath = _write_fasta(temp, "tiny.fasta", "ACGTACGT")
    counts = Mycelia.count_canonical_kmers(K, fpath)
    kmer_a = first(keys(counts))

    # Exact
    ex = Mycelia.ExactKmerMembership{K}()
    Test.@test length(ex) == 0
    push!(ex, kmer_a)
    Test.@test kmer_a in ex
    Test.@test length(ex) == 1
    push!(ex, kmer_a)  # idempotent
    Test.@test length(ex) == 1

    # Bloom
    bl = Mycelia.BloomKmerMembership{K}(capacity = 100, target_fpr = 0.001)
    Test.@test length(bl) == 0
    push!(bl, kmer_a)
    Test.@test kmer_a in bl
    Test.@test length(bl) == 1
    # No-false-negative: same kmer must still be present after re-insertion.
    push!(bl, kmer_a)
    Test.@test kmer_a in bl

    # union! with a Set
    ex2 = Mycelia.ExactKmerMembership{K}()
    union!(ex2, [kmer_a])
    Test.@test kmer_a in ex2

    bl2 = Mycelia.BloomKmerMembership{K}(capacity = 10, target_fpr = 0.01)
    union!(bl2, [kmer_a])
    Test.@test kmer_a in bl2

    rm(temp; recursive = true, force = true)
end

# ---------------------------------------------------------------------------
# Reference scaling — multi-reference union must agree with hand union
# ---------------------------------------------------------------------------

Test.@testset "pangenome_kmer_fdr — multi-reference union" begin
    temp = mktempdir()
    try
        s1 = _write_fasta(temp, "s1.fasta", "ACGTACGTACGTACGT")
        r1 = _write_fasta(temp, "r1.fasta", "ACGTACGTAAAAAAAA")
        r2 = _write_fasta(temp, "r2.fasta", "GCGCGCGCGCGCGCGC")
        r3 = _write_fasta(temp, "r3.fasta", "TTTTTTTTTTTTTTTT")

        kmer_type = Kmers.DNAKmer{4}
        result = Mycelia.pangenome_kmer_fdr(
            [s1], [r1, r2, r3]; kmer_type = kmer_type,
            cache = :memory, progress = false)

        ref_union = Set{Kmers.DNAKmer{4}}()
        for r in (r1, r2, r3)
            union!(ref_union, keys(Mycelia.count_canonical_kmers(kmer_type, r)))
        end
        Test.@test result.reference_union_size == length(ref_union)

        s1_set = Set(keys(Mycelia.count_canonical_kmers(kmer_type, s1)))
        Test.@test result.n_shared[1] == count(in(ref_union), s1_set)
    finally
        rm(temp; recursive = true, force = true)
    end
end

println("pangenome_kmer_fdr tests: all assertions passed")

# ---------------------------------------------------------------------------
# Bloom FPR characterization — empirical convergence to configured target
# ---------------------------------------------------------------------------

Test.@testset "BloomKmerMembership FPR convergence" begin
    rng = StableRNGs.StableRNG(0xb100f)
    kmer_type = Kmers.DNAKmer{21}
    kmer_length = 21
    n_inserted = 10_000
    n_queries = 100_000
    alphabet = ('A', 'C', 'G', 'T')

    random_kmers = function (n; excluded = nothing)
        kmers = Set{kmer_type}()
        while length(kmers) < n
            kmer = kmer_type(join(rand(rng, alphabet, kmer_length)))
            if excluded === nothing || !(kmer in excluded)
                push!(kmers, kmer)
            end
        end
        return kmers
    end

    println("BloomKmerMembership empirical FPR:")
    println("target_fpr\tobserved_fpr\tratio")

    for target_fpr in (0.01, 0.001, 0.0001)
        bloom = Mycelia.BloomKmerMembership{kmer_type}(
            capacity = n_inserted,
            target_fpr = target_fpr
        )

        inserted_kmers = random_kmers(n_inserted)
        foreach(kmer -> push!(bloom, kmer), inserted_kmers)

        query_kmers = random_kmers(n_queries; excluded = inserted_kmers)
        Test.@test isempty(intersect(inserted_kmers, query_kmers))
        Test.@test all(kmer -> kmer in bloom, inserted_kmers)

        false_positives = count(kmer -> kmer in bloom, query_kmers)
        observed_fpr = false_positives / n_queries
        ratio = observed_fpr / target_fpr

        println(
            "$(target_fpr)\t$(observed_fpr)\t$(round(ratio; sigdigits = 4))")
        Test.@test observed_fpr <= 1.5 * target_fpr
        Test.@test observed_fpr >= 0.3 * target_fpr
    end
end
