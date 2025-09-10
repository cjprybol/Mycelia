# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/utility_functions.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/utility_functions.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import Random

Test.@testset "scientific notation" begin
    Test.@test Mycelia.scientific_notation(100) == "1.00e+02"
    Test.@test Mycelia.scientific_notation(1000, precision=3) == "1.000e+03"
    Test.@test Mycelia.scientific_notation(0) == "0.00e+00"
    Test.@test_throws ErrorException Mycelia.scientific_notation(1; precision=-1)
end

Test.@testset "byte formatting" begin
    Test.@test Mycelia.bytes_human_readable(0) == "0 bytes"
    Test.@test Mycelia.bytes_human_readable(1) == "1 byte"
    Test.@test Mycelia.bytes_human_readable(1024^0 + 1024^0) == "2 bytes"
    Test.@test Mycelia.bytes_human_readable(1024^0 + 1024^1) == "1.001 KiB"
    Test.@test Mycelia.bytes_human_readable(1024^1 + 1024^2) == "1.001 MiB"
    Test.@test Mycelia.bytes_human_readable(1024^2 + 1024^3) == "1.001 GiB"
    Test.@test Mycelia.bytes_human_readable(1024^3 + 1024^4) == "1.001 TiB"
    Test.@test Mycelia.bytes_human_readable(1024^4 + 1024^5) == "1.001 PiB"
    Test.@test Mycelia.bytes_human_readable(1024^5 + 1024^6) == "1025.000 PiB"
end

Test.@testset "Estimate memory utilization of dense and sparse Matrices by datatype" begin
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float64, 10_000, 10_000)) == "762.939 MiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float64, 100_000, 100_000)) == "74.506 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float64, 10^5, 10^5, density = 0.5)) == "74.507 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float16, 100_000, 100_000)) == "18.626 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float16, 10^5, 10^5, density = 0.2)) == "18.627 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float16, 10^6, 10^6)) == "1.819 TiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float16, 10^6, 10^6, density = 0.1)) == "931.330 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float16, 10^6, 10^6, density = 0.2)) == "1.819 TiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float64, 10^5, 10^5, density = 0.2)) == "29.803 GiB"
    # uses default byte size for Float64
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(10^5, 10^5, density = 0.2)) == "29.803 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Int64, 10^5, 10^5, density = 0.2)) == "29.803 GiB"
    Test.@test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(UInt64, 10^5, 10^5, density = 0.2)) == "29.803 GiB"

    Test.@test Mycelia.estimate_dense_matrix_memory(Int32, 2, 2) == 16
    Test.@test_throws ArgumentError Mycelia.estimate_dense_matrix_memory(10)

    Test.@test Mycelia.estimate_sparse_matrix_memory(Float32, 2, 2, nnz=2) == 48
    Test.@test Mycelia.estimate_sparse_matrix_memory(Float32, 2, 2, density=0.5) == 48
    Test.@test Mycelia.estimate_sparse_matrix_memory(2, 2, density=0.25) == 40
    Test.@test_throws ArgumentError Mycelia.estimate_sparse_matrix_memory(Float32, 2, 2)
end

Test.@testset "check matrix fits in memory" begin
    Mycelia.system_overview()
    assessed_memory_needs = Mycelia.check_matrix_fits_in_memory(
        Mycelia.estimate_dense_matrix_memory(10^3, 10^3)
    )
    Test.@test assessed_memory_needs.will_fit_available == (assessed_memory_needs.bytes_needed <= assessed_memory_needs.free_memory)
    Test.@test assessed_memory_needs.will_fit_total == (assessed_memory_needs.bytes_needed <= assessed_memory_needs.total_memory)

    too_big = Sys.total_memory() * 2
    Test.@test_throws ErrorException Mycelia.check_matrix_fits_in_memory(too_big; severity=:error)
end

Test.@testset "base hash functions" begin
    Test.@testset "case normalization" begin
        test_seq_upper = "ATCGATCGATCG"
        test_seq_lower = "atcgatcgatcg"
        
        # Test case normalization enabled (default)
        hash_upper_norm = Mycelia.create_blake3_hash(test_seq_upper, encoded_length=16)
        hash_lower_norm = Mycelia.create_blake3_hash(test_seq_lower, encoded_length=16)
        Test.@test hash_upper_norm == hash_lower_norm
        
        # Test case normalization disabled
        hash_upper_no_norm = Mycelia.create_blake3_hash(test_seq_upper, encoded_length=16, normalize_case=false)
        hash_lower_no_norm = Mycelia.create_blake3_hash(test_seq_lower, encoded_length=16, normalize_case=false)
        Test.@test hash_upper_no_norm != hash_lower_no_norm
    end
    
    Test.@testset "multiple encodings" begin
        test_seq = "ATCGATCGATCG"
        
        # Test different encodings
        blake3_hex = Mycelia.create_blake3_hash(test_seq, encoding=:hex, encoded_length=32)
        blake3_base58 = Mycelia.create_blake3_hash(test_seq, encoding=:base58, encoded_length=32)
        blake3_base64 = Mycelia.create_blake3_hash(test_seq, encoding=:base64, encoded_length=32)
        
        Test.@test length(blake3_hex) == 32
        Test.@test length(blake3_base58) == 32
        Test.@test length(blake3_base64) == 32
        
        # Different encodings should produce different representations
        Test.@test blake3_hex != blake3_base58
        Test.@test blake3_base58 != blake3_base64
        Test.@test blake3_hex != blake3_base64
        
        # But all should represent the same underlying hash
        Test.@test all(c -> c in "0123456789abcdef", blake3_hex)  # Hex characters only
    end
    
    Test.@testset "multiple hash algorithms" begin
        test_seq = "ATCGATCGATCG"
        
        # Test different hash algorithms with truncation allowed
        md5_hash = Mycelia.create_md5_hash(test_seq, encoded_length=16, allow_truncation=true)
        sha1_hash = Mycelia.create_sha1_hash(test_seq, encoded_length=16, allow_truncation=true)
        sha256_hash = Mycelia.create_sha256_hash(test_seq, encoded_length=16, allow_truncation=true)
        sha512_hash = Mycelia.create_sha512_hash(test_seq, encoded_length=16, allow_truncation=true)
        blake3_hash = Mycelia.create_blake3_hash(test_seq, encoded_length=16)
        
        Test.@test length(md5_hash) == 16
        Test.@test length(sha1_hash) == 16
        Test.@test length(sha256_hash) == 16
        Test.@test length(sha512_hash) == 16
        Test.@test length(blake3_hash) == 16
        
        # All algorithms should produce different hashes
        hashes = [md5_hash, sha1_hash, sha256_hash, sha512_hash, blake3_hash]
        Test.@test length(Set(hashes)) == 5
    end
    
    Test.@testset "CRC32 functions" begin
        test_seq = "ATCGATCGATCG"
        
        # Test raw CRC32 checksum
        crc32_raw = Mycelia.crc32_checksum(test_seq)
        Test.@test isa(crc32_raw, UInt32)
        
        # Test formatted CRC32 hash
        crc32_hex = Mycelia.create_crc32_hash(test_seq)
        Test.@test isa(crc32_hex, String)
        Test.@test length(crc32_hex) == 8  # 4 bytes * 2 hex chars = 8
        
        # Case normalization should affect CRC32
        crc32_upper = Mycelia.crc32_checksum("ATCG")
        crc32_lower = Mycelia.crc32_checksum("atcg")
        Test.@test crc32_upper == crc32_lower  # Default case normalization
        
        crc32_upper_no_norm = Mycelia.crc32_checksum("ATCG", normalize_case=false)
        crc32_lower_no_norm = Mycelia.crc32_checksum("atcg", normalize_case=false)
        Test.@test crc32_upper_no_norm != crc32_lower_no_norm
    end
    
    Test.@testset "hash length validation and truncation" begin
        test_seq = "ATCGATCGATCG"
        
        # Test that truncation is required for shorter lengths
        Test.@test_throws ErrorException Mycelia.create_sha256_hash(test_seq, encoded_length=16, allow_truncation=false)
        
        # Test that truncation works when allowed
        sha256_truncated = Mycelia.create_sha256_hash(test_seq, encoded_length=16, allow_truncation=true)
        Test.@test length(sha256_truncated) == 16
        
        # Test natural lengths work without truncation
        sha256_natural = Mycelia.create_sha256_hash(test_seq)
        Test.@test length(sha256_natural) == 64  # 32 bytes * 2 hex chars
        
        md5_natural = Mycelia.create_md5_hash(test_seq)
        Test.@test length(md5_natural) == 32  # 16 bytes * 2 hex chars
    end
    
    Test.@testset "Blake3 dynamic byte calculation" begin
        test_seq = "ATCGATCGATCG"
        
        # Test that Blake3 calculates required bytes correctly for different encodings
        hex_32 = Mycelia.create_blake3_hash(test_seq, encoding=:hex, encoded_length=32)
        base58_32 = Mycelia.create_blake3_hash(test_seq, encoding=:base58, encoded_length=32)
        
        Test.@test length(hex_32) == 32
        Test.@test length(base58_32) == 32
        
        # Test default Blake3 length (64 characters)
        blake3_default = Mycelia.create_blake3_hash(test_seq)
        Test.@test length(blake3_default) == 64  # Default is 64 characters for tree-of-life scale
    end
end
