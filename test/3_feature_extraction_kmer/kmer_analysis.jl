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
##     Pkg.activate("../..")
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
        fasta_list=fasta_files,
        k=k,
        alphabet=:DNA,
    )
    sparse_res = Mycelia.fasta_list_to_sparse_kmer_counts(
        fasta_list=fasta_files,
        k=k,
        alphabet=:DNA,
        skip_rarefaction_plot=true,
        skip_rarefaction_data=true,
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
                    filename=save_file,
                    kmers=test_results.kmers,
                    counts=test_results.counts,
                    fasta_list=test_files[1:2],  ## Use first 2 test files
                    k=k,
                    alphabet=:DNA
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
                prime_ks=prime_ks,
                window_size=8,
                step_size=4
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
            adaptive_results = Mycelia.adaptive_kmer_selection(results, target_coverage=5.0)
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
        Test.@test Mycelia.optimal_subsequence_length(error_rate=0.001, sequence_length=100, threshold=0.99) == 10
    end
end
