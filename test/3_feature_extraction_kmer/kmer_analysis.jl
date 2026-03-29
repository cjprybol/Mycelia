# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/3_feature_extraction_kmer/kmer_analysis.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/3_feature_extraction_kmer/kmer_analysis.jl", "test/3_feature_extraction_kmer", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import FASTX
import Random
import BioSequences
import Kmers
import OrderedCollections

Test.@testset "fasta_list_to_*_kmer_counts" begin
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
        fasta_list = fasta_files,
        k = k,
        alphabet = :DNA
    )
    sparse_res = Mycelia.fasta_list_to_sparse_kmer_counts(
        fasta_list = fasta_files,
        k = k,
        alphabet = :DNA,
        skip_rarefaction_plot = true,
        skip_rarefaction_data = true
    )

    expected_rows_dense = length(Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET))
    nfiles = length(fasta_files)

    ## Test dense results (should have all possible canonical k-mers)
    Test.@test size(dense_res.counts, 1) == expected_rows_dense
    Test.@test size(dense_res.counts, 2) == nfiles

    ## Test sparse results (should only have observed k-mers, which is fewer)
    Test.@test size(sparse_res.counts, 1) <= expected_rows_dense  ## Sparse has only observed k-mers
    Test.@test size(sparse_res.counts, 2) == nfiles
    Test.@test length(sparse_res.kmers) == size(sparse_res.counts, 1)  ## Kmers and rows should match

    for (i, seq) in enumerate(sequences)
        expected_total = length(seq) - k + 1
        Test.@test sum(dense_res.counts[:, i]) == expected_total
        Test.@test sum(sparse_res.counts[:, i]) == expected_total
    end

    ## Clean up test FASTA files
    for f in fasta_files
        isfile(f) && rm(f)
    end

    ## Clean up generated rarefaction files (these are created by sparse function)
    rarefaction_files = filter(f -> occursin("rarefaction", f), readdir("."))
    for f in rarefaction_files
        if any(endswith(f, ext) for ext in [".png", ".pdf", ".svg", ".tsv"])
            Base.@info "Cleaning up generated file: $f"
            rm(f)
        end
    end
end

Test.@testset "K-mer Analysis Tests" begin

    ## Test data setup
    test_dna_sequence = "ACGTACGTACGT"
    test_rna_sequence = "ACGUACGUACGU"
    test_aa_sequence = "ACDEFGHIKLMN"

    ## Create test files
    test_files = String[]

    # DNA FASTA file
    dna_file = "test_dna.fasta"
    open(dna_file, "w") do io
        FASTX.write(io, FASTX.FASTA.Record("dna_seq", test_dna_sequence))
    end
    push!(test_files, dna_file)

    # RNA FASTA file  
    rna_file = "test_rna.fasta"
    open(rna_file, "w") do io
        FASTX.write(io, FASTX.FASTA.Record("rna_seq", test_rna_sequence))
    end
    push!(test_files, rna_file)

    # AA FASTA file
    aa_file = "test_aa.fasta"
    open(aa_file, "w") do io
        FASTX.write(io, FASTX.FASTA.Record("aa_seq", test_aa_sequence))
    end
    push!(test_files, aa_file)

    try
        ## Test 1: Core K-mer Generation Functions
        Test.@testset "K-mer Generation Functions" begin
            k = 3

            # Test generate_all_possible_kmers for different alphabets
            dna_kmers = Mycelia.generate_all_possible_kmers(k, Mycelia.DNA_ALPHABET)
            Test.@test length(dna_kmers) == 4^k  ## 4^3 = 64 for DNA
            Test.@test isa(dna_kmers, Vector)  ## Should be vector of k-mer objects

            rna_kmers = Mycelia.generate_all_possible_kmers(k, Mycelia.RNA_ALPHABET)
            Test.@test length(rna_kmers) == 4^k  ## 4^3 = 64 for RNA
            Test.@test isa(rna_kmers, Vector)  ## Should be vector of k-mer objects

            # Test AA alphabet (22 amino acids: 20 standard + O, U, no ambiguous symbols)
            aa_kmers = Mycelia.generate_all_possible_kmers(2, Mycelia.AA_ALPHABET)
            Test.@test length(aa_kmers) == length(Mycelia.AA_ALPHABET)^2  ## 22^2 = 484 (20 standard + O + U)
            Test.@test length(aa_kmers) == 22^2  ## 22^2 = 484 (20 standard + O + U)
            Test.@test isa(aa_kmers, Vector)  ## Should be vector of k-mer objects

            # Test canonical k-mer generation
            canonical_dna_kmers = Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET)
            Test.@test length(canonical_dna_kmers) <= length(dna_kmers)  ## Should be <= due to canonicalization
            Test.@test length(canonical_dna_kmers) == 32  ## Expected canonical 3-mers for DNA
        end

        ## Test 2: Individual K-mer Counting Functions  
        Test.@testset "Individual K-mer Counting" begin
            k = 3

            # Test counting on DNA sequence
            dna_seq = BioSequences.LongDNA{4}(test_dna_sequence)
            dna_kmer_counts = Mycelia.count_kmers(Kmers.DNAKmer{k}, dna_seq)
            Test.@test isa(dna_kmer_counts, OrderedCollections.OrderedDict)
            Test.@test sum(values(dna_kmer_counts)) == length(test_dna_sequence) - k + 1

            # Test canonical counting
            canonical_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, dna_file)
            Test.@test isa(canonical_counts, OrderedCollections.OrderedDict)
            Test.@test sum(values(canonical_counts)) == length(test_dna_sequence) - k + 1

            # Test counting from file
            file_counts = Mycelia.count_kmers(Kmers.DNAKmer{k}, dna_file)
            Test.@test isa(file_counts, OrderedCollections.OrderedDict)
            Test.@test sum(values(file_counts)) == length(test_dna_sequence) - k + 1
        end

        ## Test 3: Different Alphabet Support
        Test.@testset "Different Alphabet Support" begin
            k = 3

            # RNA k-mer counting
            rna_counts = Mycelia.count_kmers(Kmers.RNAKmer{k}, rna_file)
            Test.@test isa(rna_counts, OrderedCollections.OrderedDict)
            Test.@test sum(values(rna_counts)) == length(test_rna_sequence) - k + 1

            # Amino acid k-mer counting (use k=2 for the 23-symbol alphabet)
            k_aa = 2
            aa_counts = Mycelia.count_kmers(Kmers.AAKmer{k_aa}, aa_file)
            Test.@test isa(aa_counts, OrderedCollections.OrderedDict)
            Test.@test sum(values(aa_counts)) == length(test_aa_sequence) - k_aa + 1
        end

        ## Test 4: K-mer Canonicalization
        Test.@testset "K-mer Canonicalization" begin
            # Create a test case with reverse complements
            test_seq = "ACGTACGT"
            k = 4
            kmer_counts = Mycelia.count_kmers(Kmers.DNAKmer{k}, BioSequences.LongDNA{4}(test_seq))

            # Test canonicalize_kmer_counts (non-mutating)
            canonical_counts = Mycelia.canonicalize_kmer_counts(kmer_counts)
            Test.@test isa(canonical_counts, OrderedCollections.OrderedDict)
            Test.@test sum(values(canonical_counts)) == sum(values(kmer_counts))

            # Test canonicalize_kmer_counts! (mutating)
            original_counts = copy(kmer_counts)
            Mycelia.canonicalize_kmer_counts!(kmer_counts)
            Test.@test sum(values(kmer_counts)) == sum(values(original_counts))
        end

        ## Test 5: Genome Size Estimation
        Test.@testset "Genome Size Estimation" begin
            k = 3

            # Test with BioSequence - function returns a Dict with detailed info
            dna_bioseq = BioSequences.LongDNA{4}(test_dna_sequence)
            estimated_result = Mycelia.estimate_genome_size_from_kmers(dna_bioseq, k)
            Test.@test isa(estimated_result, Dict)
            Test.@test haskey(estimated_result, "estimated_genome_size")
            Test.@test estimated_result["estimated_genome_size"] > 0

            # Test with string input (should be converted to BioSequence internally)
            estimated_result_string = Mycelia.estimate_genome_size_from_kmers(test_dna_sequence, k)
            Test.@test isa(estimated_result_string, Dict)
            Test.@test haskey(estimated_result_string, "estimated_genome_size")
            Test.@test estimated_result_string["estimated_genome_size"] > 0
            Test.@test estimated_result_string == estimated_result  ## Should be same result

            # Test with FASTA records
            records = FASTX.FASTA.Reader(open(dna_file)) |> collect
            try
                estimated_result_records = Mycelia.estimate_genome_size_from_kmers(records, k)
                Test.@test isa(estimated_result_records, Dict)
                Test.@test haskey(estimated_result_records, "estimated_genome_size")
                Test.@test estimated_result_records["estimated_genome_size"] > 0
            finally
                # File should already be closed from collect operation
            end
        end

        ## Test 6: File I/O Operations
        Test.@testset "K-mer File I/O" begin
            k = 3

            # Create some k-mer results to save
            test_results = (
                kmers = Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET),
                counts = rand(UInt16, 32, 2),  ## 32 canonical k-mers, 2 files
                metadata = Dict("k" => k, "alphabet" => "DNA")
            )

            # Test saving
            save_file = "test_kmer_results.jld2"
            try
                Mycelia.save_kmer_results(
                    filename = save_file,
                    kmers = test_results.kmers,
                    counts = test_results.counts,
                    fasta_list = test_files[1:2],  ## Use first 2 test files
                    k = k,
                    alphabet = :DNA
                )
                Test.@test isfile(save_file)

                # Test loading
                loaded_results = Mycelia.load_kmer_results(save_file)
                Test.@test haskey(loaded_results, :kmers)
                Test.@test haskey(loaded_results, :counts)
                Test.@test haskey(loaded_results, :metadata)
                Test.@test size(loaded_results.counts) == size(test_results.counts)
                Test.@test loaded_results.metadata["k"] == k
            finally
                isfile(save_file) && rm(save_file)
            end
        end

        ## Test 7: Multi-scale Analysis
        Test.@testset "Multi-scale K-mer Analysis" begin
            # Create test sequences
            sequences = [
                BioSequences.LongDNA{4}("ACGTACGTACGTACGT"),
                BioSequences.LongDNA{4}("GGGGCCCCAAAATTTT"),
                BioSequences.LongDNA{4}("ATCGATCGATCGATCG")
            ]

            # Test multi-scale analysis with smaller k values for faster testing
            prime_ks = [3, 5, 7]
            results = Mycelia.multi_scale_kmer_analysis(
                sequences,
                prime_ks = prime_ks,
                window_size = 8,
                step_size = 4
            )

            Test.@test isa(results, NamedTuple)  ## Returns NamedTuple, not Dict
            Test.@test haskey(results, :multi_k_counts)  ## Check for expected fields
            Test.@test haskey(results, :quality_profiles)
            Test.@test haskey(results, :coverage_estimates)

            # Check that the multi_k_counts contains our k values
            for k in prime_ks
                Test.@test haskey(results.multi_k_counts, k)
                Test.@test isa(results.multi_k_counts[k], Dict)  ## Should be Dict of k-mer counts
            end

            # Test adaptive k-mer selection
            adaptive_results = Mycelia.adaptive_kmer_selection(results, target_coverage = 5.0)
            Test.@test isa(adaptive_results, Vector)  ## Returns vector of recommended k values
            Test.@test all(k -> k in prime_ks, adaptive_results)  ## Should be subset of input ks
        end

    catch e
        Base.@error "Test failed with error: $e"
        rethrow(e)
    finally
        ## Cleanup test files
        for file in test_files
            isfile(file) && rm(file)
        end

        ## Clean up any generated output files
        for file in readdir(".")
            if occursin(r"\.(jld2|tsv|png|pdf|svg)$", file) &&
               (occursin("test", file) || occursin("kmer", file))
                Base.@info "Cleaning up generated file: $file"
                rm(file)
            end
        end
    end
end

# Reference Graph and K-mer Analysis tests
Test.@testset "Reference Graph and K-mer Analysis" begin
    Test.@testset "Pangenome Construction" begin
        Test.@test true  # placeholder
    end
    Test.@testset "Optimal K-mer Selection" begin
        Test.@test true  # placeholder
    end
    Test.@testset "statistical kmer analyses" begin
        Test.@test Mycelia.optimal_subsequence_length(
            error_rate = 0.001, sequence_length = 100, threshold = 0.99) == 10
    end
end

import StatsBase

Test.@testset "kmer_frequency_histogram" begin
    # Dict overload
    kmer_counts = OrderedCollections.OrderedDict(
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("ACG")) => 5,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("CGT")) => 5,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("GTA")) => 2,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("TAC")) => 1
    )
    hist = Mycelia.kmer_frequency_histogram(kmer_counts)
    Test.@test hist[5] == 2  # two kmers with count 5
    Test.@test hist[2] == 1  # one kmer with count 2
    Test.@test hist[1] == 1  # one kmer with count 1

    # Vector overload
    counts_vec = [5, 5, 2, 1]
    hist_vec = Mycelia.kmer_frequency_histogram(counts_vec)
    Test.@test hist_vec == hist
end

Test.@testset "coverage_peak_from_hist" begin
    hist = Dict(1 => 100, 5 => 50, 10 => 200, 20 => 30)

    # Default min_coverage=2 should skip coverage=1
    peak = Mycelia.coverage_peak_from_hist(hist)
    Test.@test peak.coverage == 10
    Test.@test peak.kmers == 200

    # With min_coverage=1 the coverage=1 entry has most kmers
    peak_low = Mycelia.coverage_peak_from_hist(hist; min_coverage = 1)
    Test.@test peak_low.coverage == 10  # 200 > 100

    # Empty histogram
    empty_peak = Mycelia.coverage_peak_from_hist(Dict{Int, Int}())
    Test.@test ismissing(empty_peak.coverage)
    Test.@test empty_peak.kmers == 0
end

Test.@testset "analyze_kmer_frequency_spectrum" begin
    kmer_counts = OrderedCollections.OrderedDict(
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("ACG")) => 5,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("CGT")) => 5,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("GTA")) => 2,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{2}("TAC")) => 1
    )
    result = Mycelia.analyze_kmer_frequency_spectrum(kmer_counts)
    Test.@test haskey(result, :histogram)
    Test.@test haskey(result, :peak)
    Test.@test haskey(result, :total_kmers)
    Test.@test haskey(result, :unique_kmers)
    Test.@test result.total_kmers == 13  # 5+5+2+1
    Test.@test result.unique_kmers == 4

    # Vector overload
    result_vec = Mycelia.analyze_kmer_frequency_spectrum([5, 5, 2, 1])
    Test.@test result_vec.total_kmers == 13
    Test.@test result_vec.unique_kmers == 4
end

Test.@testset "canonical_minimizers" begin
    seq = BioSequences.LongDNA{4}("ACGTACGTACGT")
    k = 3
    window = 3

    minimizers = Mycelia.canonical_minimizers(seq, k, window)
    Test.@test isa(minimizers, Vector)
    Test.@test length(minimizers) > 0
    # Number of minimizers = num_kmers - window + 1
    num_kmers = length(seq) - k + 1  # 10
    Test.@test length(minimizers) == num_kmers - window + 1  # 8

    # Edge case: window larger than number of kmers returns empty
    empty_min = Mycelia.canonical_minimizers(seq, k, 100)
    Test.@test isempty(empty_min)

    # Error on invalid k
    Test.@test_throws ErrorException Mycelia.canonical_minimizers(seq, 0, 3)
    # Error on invalid window
    Test.@test_throws ErrorException Mycelia.canonical_minimizers(seq, 3, 0)

    # RNA minimizers (non-canonical path for nucleotides still works)
    rna_seq = BioSequences.LongRNA{4}("ACGUACGUACGU")
    rna_min = Mycelia.canonical_minimizers(rna_seq, 3, 3)
    Test.@test length(rna_min) > 0
end

Test.@testset "open_syncmers" begin
    seq = BioSequences.LongDNA{4}("ACGTACGTACGTACGT")
    k = 5
    s = 3
    t = 1

    syncmers = Mycelia.open_syncmers(seq, k, s, t)
    Test.@test isa(syncmers, Vector)
    # Each syncmer should be a k-mer extracted from the sequence
    for sm in syncmers
        Test.@test length(sm) == k
    end

    # With canonical=true
    syncmers_canon = Mycelia.open_syncmers(seq, k, s, t; canonical = true)
    Test.@test isa(syncmers_canon, Vector)

    # Error on invalid s
    Test.@test_throws ErrorException Mycelia.open_syncmers(seq, 5, 0, 1)
    Test.@test_throws ErrorException Mycelia.open_syncmers(seq, 5, 6, 1)

    # Error on invalid t
    Test.@test_throws ErrorException Mycelia.open_syncmers(seq, 5, 3, 0)
    Test.@test_throws ErrorException Mycelia.open_syncmers(seq, 5, 3, 4)  # max_position = 5-3+1 = 3

    # AA syncmers (SingleStrand, no canonical)
    aa_seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
    aa_syncmers = Mycelia.open_syncmers(aa_seq, 4, 2, 1)
    Test.@test isa(aa_syncmers, Vector)
end

Test.@testset "strobemers" begin
    seq = BioSequences.LongDNA{4}("ACGTACGTACGTACGT")
    k = 3
    w_min = 2
    w_max = 4

    strobes = Mycelia.strobemers(seq, k, w_min, w_max)
    Test.@test isa(strobes, Vector{<:Tuple})
    Test.@test length(strobes) > 0
    # Each element is a pair of k-mers
    for (s1, s2) in strobes
        Test.@test length(s1) == k
        Test.@test length(s2) == k
    end

    # canonical=true
    strobes_canon = Mycelia.strobemers(seq, k, w_min, w_max; canonical = true)
    Test.@test length(strobes_canon) > 0

    # Error on invalid w_min/w_max
    Test.@test_throws ErrorException Mycelia.strobemers(seq, 3, 0, 4)
    Test.@test_throws ErrorException Mycelia.strobemers(seq, 3, 5, 3)
end

Test.@testset "determine_max_possible_kmers and determine_max_canonical_kmers" begin
    # max possible kmers
    Test.@test Mycelia.determine_max_possible_kmers(3, Mycelia.DNA_ALPHABET) == 4^3
    Test.@test Mycelia.determine_max_possible_kmers(2, Mycelia.AA_ALPHABET) ==
               length(Mycelia.AA_ALPHABET)^2
    Test.@test Mycelia.determine_max_possible_kmers(3, Mycelia.RNA_ALPHABET) == 4^3

    # max canonical kmers - for AA returns all, for DNA/RNA returns half (odd k)
    Test.@test Mycelia.determine_max_canonical_kmers(3, Mycelia.DNA_ALPHABET) == 4^3 ÷ 2
    Test.@test Mycelia.determine_max_canonical_kmers(2, Mycelia.AA_ALPHABET) ==
               length(Mycelia.AA_ALPHABET)^2
    Test.@test Mycelia.determine_max_canonical_kmers(3, Mycelia.RNA_ALPHABET) == 4^3 ÷ 2

    # Even k should error for nucleotide alphabets
    Test.@test_throws AssertionError Mycelia.determine_max_canonical_kmers(4, Mycelia.DNA_ALPHABET)
end

Test.@testset "kmer_counts_dict_to_vector" begin
    k = 3
    all_kmers = Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET)
    kmer_to_index = Dict(kmer => i for (i, kmer) in enumerate(all_kmers))

    # Build a small counts dict
    seq = BioSequences.LongDNA{4}("ACGTACGT")
    kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, seq)

    vec = Mycelia.kmer_counts_dict_to_vector(kmer_to_index, kmer_counts)
    Test.@test length(vec) == length(all_kmers)
    Test.@test sum(vec) == sum(values(kmer_counts))
    # Zeros for absent kmers
    Test.@test count(==(0), vec) == length(all_kmers) - length(kmer_counts)
end

Test.@testset "canonicalize_kmer_counts! with AA kmers (no-op)" begin
    aa_seq = BioSequences.LongAA("ACDEFGHIKLMN")
    aa_counts = Mycelia.count_kmers(Kmers.AAKmer{2}, aa_seq)
    original_counts = copy(aa_counts)

    # AA kmers should not be modified by canonicalization
    Mycelia.canonicalize_kmer_counts!(aa_counts)
    Test.@test aa_counts == original_counts
end

Test.@testset "_is_nucleotide_sequence helper" begin
    Test.@test Mycelia._is_nucleotide_sequence(BioSequences.LongDNA{4}("ACGT")) == true
    Test.@test Mycelia._is_nucleotide_sequence(BioSequences.LongRNA{4}("ACGU")) == true
    Test.@test Mycelia._is_nucleotide_sequence(BioSequences.LongAA("ACD")) == false
end

Test.@testset "_sequence_kmer_iterator and _collect_kmers" begin
    dna = BioSequences.LongDNA{4}("ACGTACGT")
    rna = BioSequences.LongRNA{4}("ACGUACGU")
    aa = BioSequences.LongAA("ACDEFGHI")

    dna_kmers = Mycelia._collect_kmers(dna, 3)
    Test.@test length(dna_kmers) == length(dna) - 3 + 1

    rna_kmers = Mycelia._collect_kmers(rna, 3)
    Test.@test length(rna_kmers) == length(rna) - 3 + 1

    aa_kmers = Mycelia._collect_kmers(aa, 3)
    Test.@test length(aa_kmers) == length(aa) - 3 + 1

    # Unsupported type should error
    Test.@test_throws ErrorException Mycelia._sequence_kmer_iterator("ACGT", 3)
end

Test.@testset "calculate_window_quality" begin
    window = BioSequences.LongDNA{4}("ACGTACGTACGT")
    quality = Mycelia.calculate_window_quality(window, 3)
    Test.@test quality isa Float64
    Test.@test 0.0 <= quality <= 1.0

    # Short window (shorter than k) returns 0.0
    short = BioSequences.LongDNA{4}("AC")
    Test.@test Mycelia.calculate_window_quality(short, 5) == 0.0

    # Repetitive sequence should have lower diversity
    repetitive = BioSequences.LongDNA{4}("AAAAAAAAAAAA")
    rep_quality = Mycelia.calculate_window_quality(repetitive, 3)
    Test.@test rep_quality < quality  # less diverse = lower quality
end

Test.@testset "ks function" begin
    result = Mycelia.ks()
    Test.@test isa(result, Vector{Int})
    Test.@test !isempty(result)
    Test.@test issorted(result)
    # Should start with small odd primes (3, 5, 7, ...)
    Test.@test 3 in result
    Test.@test 5 in result
    Test.@test 7 in result
    # 19 should be skipped
    Test.@test !(19 in result)

    # With min/max filters
    filtered = Mycelia.ks(min = 5, max = 50)
    Test.@test all(x -> 5 <= x <= 50, filtered)
end

Test.@testset "k_ladder function" begin
    result = Mycelia.k_ladder()
    Test.@test isa(result, Vector{Int})
    Test.@test !isempty(result)
    Test.@test issorted(result)
    # Default seeds include 3, 5, 7
    Test.@test 3 in result
    Test.@test 5 in result
    Test.@test 7 in result
    # All should be odd primes (check manually for small primes)
    Test.@test all(x -> isodd(x) && x > 1 && all(d -> x % d != 0, 2:isqrt(x)), result)

    # With read_length cap
    short_result = Mycelia.k_ladder(read_length = 50, read_margin = 10)
    Test.@test all(x -> x <= 40, short_result)

    # Custom seed primes
    custom = Mycelia.k_ladder(seed_primes = [11, 13, 17], max_k = 100)
    Test.@test 11 in custom
    Test.@test 13 in custom
    Test.@test 17 in custom
end

Test.@testset "count_kmers with RNA and AA sequences directly" begin
    # RNA from sequence
    rna_seq = BioSequences.LongRNA{4}("ACGUACGUACGU")
    k = 3
    rna_counts = Mycelia.count_kmers(Kmers.RNAKmer{k}, rna_seq)
    Test.@test isa(rna_counts, OrderedCollections.OrderedDict)
    Test.@test sum(values(rna_counts)) == length(rna_seq) - k + 1

    # AA from sequence
    aa_seq = BioSequences.LongAA("ACDEFGHIKLMN")
    k_aa = 2
    aa_counts = Mycelia.count_kmers(Kmers.AAKmer{k_aa}, aa_seq)
    Test.@test isa(aa_counts, OrderedCollections.OrderedDict)
    Test.@test sum(values(aa_counts)) == length(aa_seq) - k_aa + 1
end

Test.@testset "count_kmers with FASTA records vector" begin
    records = FASTX.FASTA.Record[]
    push!(records, FASTX.FASTA.Record("s1", "ACGTACGT"))
    push!(records, FASTX.FASTA.Record("s2", "GGGGAAAA"))
    k = 3
    counts = Mycelia.count_kmers(Kmers.DNAKmer{k}, records)
    Test.@test isa(counts, OrderedCollections.OrderedDict)
    expected_total = (8 - k + 1) * 2  # 6 kmers per 8-base seq, 2 seqs
    Test.@test sum(values(counts)) == expected_total
end

Test.@testset "generate_all_possible_kmers RNA" begin
    k = 2
    rna_kmers = Mycelia.generate_all_possible_kmers(k, Mycelia.RNA_ALPHABET)
    Test.@test length(rna_kmers) == 4^k
    Test.@test issorted(rna_kmers)
end

Test.@testset "generate_all_possible_canonical_kmers RNA" begin
    k = 3
    canonical_rna = Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.RNA_ALPHABET)
    all_rna = Mycelia.generate_all_possible_kmers(k, Mycelia.RNA_ALPHABET)
    Test.@test length(canonical_rna) <= length(all_rna)
end

Test.@testset "generate_all_possible_canonical_kmers AA (returns all)" begin
    k = 2
    canonical_aa = Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.AA_ALPHABET)
    all_aa = Mycelia.generate_all_possible_kmers(k, Mycelia.AA_ALPHABET)
    # AA has no reverse complement, so canonical == all
    Test.@test length(canonical_aa) == length(all_aa)
end

Test.@testset "save_kmer_results with non-.jld2 filename appends extension" begin
    k = 3
    kmers = Mycelia.generate_all_possible_canonical_kmers(k, Mycelia.DNA_ALPHABET)
    counts = zeros(UInt16, length(kmers), 1)
    save_file = "test_no_ext"
    expected_file = "test_no_ext.jld2"
    try
        Mycelia.save_kmer_results(
            filename = save_file,
            kmers = kmers,
            counts = counts,
            fasta_list = ["dummy.fasta"],
            k = k,
            alphabet = :DNA
        )
        Test.@test isfile(expected_file)
    finally
        isfile(expected_file) && rm(expected_file)
        isfile(save_file) && rm(save_file)
    end
end
