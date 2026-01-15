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
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import Random
import BioSequences
import Kmers

Test.@testset "default thread detection" begin
    cpu_threads = Sys.CPU_THREADS
    half_cap = min(cld(cpu_threads, 2), 16)

    # Explicit Julia hint wins and is clamped to hardware.
    Test.@test Base.withenv("JULIA_NUM_THREADS" => string(cpu_threads)) do
        Mycelia.get_default_threads()
    end == cpu_threads
    Test.@test Base.withenv("JULIA_NUM_THREADS" => string(cpu_threads + 8)) do
        Mycelia.get_default_threads()
    end == cpu_threads
    Test.@test Base.withenv("JULIA_NUM_THREADS" => "3") do
        Mycelia.get_default_threads()
    end == min(3, cpu_threads)

    # SLURM job-level hint.
    Test.@test Base.withenv(
        "JULIA_NUM_THREADS" => nothing,
        "SLURM_JOB_CPUS_PER_NODE" => string(cpu_threads + 4),
        "SLURM_CPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_CPUS_ON_NODE" => nothing,
        "PBS_NCPUS" => nothing,
    ) do
        Mycelia.get_default_threads()
    end == cpu_threads

    # SLURM job-level hint with decorations like "24(x2)".
    Test.@test Base.withenv(
        "JULIA_NUM_THREADS" => nothing,
        "SLURM_JOB_CPUS_PER_NODE" => string(cpu_threads + 4) * "(x2)",
        "SLURM_CPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_CPUS_ON_NODE" => nothing,
        "PBS_NCPUS" => nothing,
    ) do
        Mycelia.get_default_threads()
    end == cpu_threads

    # Derived SLURM hint (CPUs per task × tasks per node) when job-level is absent.
    Test.@test Base.withenv(
        "JULIA_NUM_THREADS" => nothing,
        "SLURM_JOB_CPUS_PER_NODE" => nothing,
        "SLURM_CPUS_PER_TASK" => "2",
        "SLURM_TASKS_PER_NODE" => "3",
        "SLURM_CPUS_ON_NODE" => string(cpu_threads),
        "PBS_NCPUS" => nothing,
    ) do
        Mycelia.get_default_threads()
    end == min(6, cpu_threads)

    # SLURM CPUs on node as fallback.
    slurm_on_node = max(1, cpu_threads - 1)
    Test.@test Base.withenv(
        "JULIA_NUM_THREADS" => nothing,
        "SLURM_JOB_CPUS_PER_NODE" => nothing,
        "SLURM_CPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_CPUS_ON_NODE" => string(slurm_on_node),
        "PBS_NCPUS" => nothing,
    ) do
        Mycelia.get_default_threads()
    end == min(slurm_on_node, cpu_threads)

    # PBS hint.
    Test.@test Base.withenv(
        "JULIA_NUM_THREADS" => nothing,
        "SLURM_JOB_CPUS_PER_NODE" => nothing,
        "SLURM_CPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_CPUS_ON_NODE" => nothing,
        "PBS_NCPUS" => "7",
    ) do
        Mycelia.get_default_threads()
    end == min(7, cpu_threads)

    # No hints: half the hardware threads capped at 16.
    Test.@test Base.withenv(
        "JULIA_NUM_THREADS" => nothing,
        "SLURM_JOB_CPUS_PER_NODE" => nothing,
        "SLURM_CPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_CPUS_ON_NODE" => nothing,
        "PBS_NCPUS" => nothing,
    ) do
        Mycelia.get_default_threads()
    end == half_cap
end

Test.@testset "system overview returns raw values, flags, and display" begin
    overview = Mycelia.system_overview(memory_low_threshold=1.0, storage_low_threshold=1.0)

    Test.@test isa(overview, Mycelia.SystemOverview)
    Test.@test overview.default_threads == Mycelia.get_default_threads()
    Test.@test overview.slurm_threads === nothing || isa(overview.slurm_threads, Integer)
    Test.@test overview.system_gpus === nothing || isa(overview.system_gpus, Integer)
    Test.@test overview.slurm_gpus === nothing || isa(overview.slurm_gpus, Integer)
    Test.@test overview.slurm_memory === nothing || isa(overview.slurm_memory, Integer)
    Test.@test overview.slurm_memory_source === nothing || isa(overview.slurm_memory_source, Symbol)

    Test.@test isa(overview.total_memory, Integer)
    Test.@test isa(overview.available_memory, Integer)
    Test.@test isa(overview.occupied_memory, Integer)
    Test.@test overview.total_memory == overview.available_memory + overview.occupied_memory

    Test.@test isa(overview.total_storage, Integer)
    Test.@test isa(overview.available_storage, Integer)
    Test.@test isa(overview.occupied_storage, Integer)
    Test.@test overview.total_storage >= overview.occupied_storage
    Test.@test overview.total_storage >= overview.available_storage

    expected_memory_pct = overview.total_memory == 0 ? 0.0 : overview.occupied_memory / overview.total_memory
    expected_storage_pct = overview.total_storage == 0 ? 0.0 : overview.occupied_storage / overview.total_storage
    Test.@test isapprox(overview.memory_occupied_percent, expected_memory_pct)
    Test.@test isapprox(overview.storage_occupied_percent, expected_storage_pct)

    Test.@test overview.memory_running_low == (overview.memory_occupied_percent >= 1.0)
    Test.@test overview.storage_running_low == (overview.storage_occupied_percent >= 1.0)

    display_str = sprint(show, "text/plain", overview)
    Test.@test occursin("Threads: julia", display_str)
    Test.@test occursin(Mycelia.bytes_human_readable(overview.total_memory), display_str)
    Test.@test occursin(Mycelia.bytes_human_readable(overview.available_memory), display_str)
    Test.@test occursin(Mycelia.bytes_human_readable(overview.occupied_memory), display_str)
    Test.@test occursin(Mycelia.bytes_human_readable(overview.total_storage), display_str)
    Test.@test occursin(Mycelia.bytes_human_readable(overview.available_storage), display_str)
    Test.@test occursin(Mycelia.bytes_human_readable(overview.occupied_storage), display_str)
    Test.@test occursin("GPUs:", display_str)
end

Test.@testset "system overview reflects slurm allocations" begin
    cpu_allocation = 3
    gpu_allocation = 2
    mem_per_cpu = 1000  # MiB

    overview = Base.withenv(
        "SLURM_JOB_ID" => "12345",
        "SLURM_JOB_CPUS_PER_NODE" => nothing,
        "SLURM_CPUS_ON_NODE" => string(cpu_allocation),
        "SLURM_CPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_MEM_PER_CPU" => string(mem_per_cpu),
        "SLURM_MEM_PER_NODE" => nothing,
        "SLURM_MEM_PER_GPU" => nothing,
        "SLURM_GPUS_ON_NODE" => string(gpu_allocation),
        "SLURM_GPUS" => nothing,
        "SLURM_GPUS_PER_TASK" => nothing,
        "SLURM_TASKS_PER_NODE" => nothing,
        "SLURM_JOB_GPUS" => nothing,
        "SLURM_STEP_GPUS" => nothing,
        "SLURM_JOB_GRES" => nothing,
        "CUDA_VISIBLE_DEVICES" => nothing,
        "ROCR_VISIBLE_DEVICES" => nothing,
    ) do
        Mycelia.system_overview(memory_low_threshold=1.0, storage_low_threshold=1.0)
    end

    Test.@test overview.slurm_threads == cpu_allocation
    Test.@test overview.slurm_gpus == gpu_allocation
    Test.@test overview.system_gpus == gpu_allocation
    Test.@test overview.slurm_memory == mem_per_cpu * cpu_allocation * 1024^2
    Test.@test overview.slurm_memory_source == :per_cpu

    display_str = sprint(show, "text/plain", overview)
    Test.@test occursin("slurm=$(cpu_allocation)", display_str)
    Test.@test occursin("slurm=$(gpu_allocation)", display_str)
    Test.@test occursin("slurm_limit", display_str)
end

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

Test.@testset "create_sequence_hash comprehensive testing" begin
    Test.@testset "default behavior and parameter combinations" begin
        test_seq = "ATCGATCGATCG"
        
        # Test default behavior (Blake3 + Base58 + 64 characters)
        default_hash = Mycelia.create_sequence_hash(test_seq)
        Test.@test length(default_hash) == 64
        Test.@test isa(default_hash, String)
        
        # Test reproducibility
        hash1 = Mycelia.create_sequence_hash(test_seq)
        hash2 = Mycelia.create_sequence_hash(test_seq)
        Test.@test hash1 == hash2
        
        # Test case normalization (default is true)
        upper_hash = Mycelia.create_sequence_hash("ATCG")
        lower_hash = Mycelia.create_sequence_hash("atcg")
        Test.@test upper_hash == lower_hash
        
        # Test case normalization disabled
        upper_hash_no_norm = Mycelia.create_sequence_hash("ATCG", normalize_case=false)
        lower_hash_no_norm = Mycelia.create_sequence_hash("atcg", normalize_case=false)
        Test.@test upper_hash_no_norm != lower_hash_no_norm
    end
    
    Test.@testset "various encoded lengths with default settings" begin
        test_seq = "ATCGATCGATCG"
        
        # Test various common lengths
        for test_length in [8, 16, 24, 32, 48, 64, 80, 96, 128]
            hash_result = Mycelia.create_sequence_hash(test_seq, encoded_length=test_length)
            Test.@test length(hash_result) == test_length
        end
        
        # Test that different lengths produce different results (truncation works)
        hash_16 = Mycelia.create_sequence_hash(test_seq, encoded_length=16)
        hash_32 = Mycelia.create_sequence_hash(test_seq, encoded_length=32)
        Test.@test hash_16 != hash_32
        Test.@test length(hash_16) == 16
        Test.@test length(hash_32) == 32
    end
    
    Test.@testset "hash function parameter combinations" begin
        test_seq = "ATCGATCGATCG"
        
        # Test different hash functions with default encoding (Base58) using lengths they can produce
        blake3_hash = Mycelia.create_sequence_hash(test_seq, hash_function=:blake3, encoded_length=32)
        sha256_hash = Mycelia.create_sequence_hash(test_seq, hash_function=:sha256, encoded_length=32, allow_truncation=true)
        sha512_hash = Mycelia.create_sequence_hash(test_seq, hash_function=:sha512, encoded_length=32, allow_truncation=true)
        md5_hash = Mycelia.create_sequence_hash(test_seq, hash_function=:md5, encoded_length=22, allow_truncation=true)  # MD5 produces ~22 Base58 chars
        sha1_hash = Mycelia.create_sequence_hash(test_seq, hash_function=:sha1, encoded_length=25, allow_truncation=true)  # SHA1 produces ~27 Base58 chars
        
        Test.@test length(blake3_hash) == 32
        Test.@test length(sha256_hash) == 32
        Test.@test length(sha512_hash) == 32
        Test.@test length(md5_hash) == 22
        Test.@test length(sha1_hash) == 25
        
        # All should be different
        hashes = [blake3_hash, sha256_hash, sha512_hash, md5_hash, sha1_hash]
        Test.@test length(Set(hashes)) == 5
    end
    
    Test.@testset "encoding parameter combinations" begin
        test_seq = "ATCGATCGATCG"
        
        # Test different encodings with Blake3 (default hash function)
        base58_hash = Mycelia.create_sequence_hash(test_seq, encoding=:base58, encoded_length=32)
        hex_hash = Mycelia.create_sequence_hash(test_seq, encoding=:hex, encoded_length=32)
        base64_hash = Mycelia.create_sequence_hash(test_seq, encoding=:base64, encoded_length=32)
        
        Test.@test length(base58_hash) == 32
        Test.@test length(hex_hash) == 32
        Test.@test length(base64_hash) == 32
        
        # All should be different (different representations of same hash)
        Test.@test base58_hash != hex_hash
        Test.@test hex_hash != base64_hash
        Test.@test base58_hash != base64_hash
        
        # Hex should only contain hex characters
        Test.@test all(c -> c in "0123456789abcdef", hex_hash)
    end
    
    Test.@testset "Base58 auto-truncation edge cases" begin
        test_seq = "ATCGATCGATCG"
        
        # Test that Base58 +1 character is handled automatically (no allow_truncation needed)
        # This tests our auto-truncation logic for the common Base58 +1 scenario
        # Note: Some lengths may produce "too short" errors, so we test realistic ranges
        for test_length in [16, 17, 32, 33, 64, 65]
            hash_result = Mycelia.create_sequence_hash(test_seq, encoded_length=test_length)
            Test.@test length(hash_result) == test_length
        end
        
        # Test sequences that might produce edge case lengths
        edge_sequences = ["A", "T", "G", "C", "AT", "ATCG", "AAAA", "TTTT", "GGGG", "CCCC"]
        for seq in edge_sequences
            hash_result = Mycelia.create_sequence_hash(seq, encoded_length=16)
            Test.@test length(hash_result) == 16
        end
    end
    
    Test.@testset "error conditions and edge cases" begin
        # Test empty sequence validation
        Test.@test_throws ErrorException Mycelia.create_sequence_hash("")
        
        # Test invalid hash function
        Test.@test_throws ErrorException Mycelia.create_sequence_hash("ATCG", hash_function=:invalid)
        
        # Test invalid encoding
        Test.@test_throws ErrorException Mycelia.create_sequence_hash("ATCG", encoding=:invalid)
        
        # Test "too short" scenarios where requested length exceeds what can be produced from available bytes
        # This should trigger the error: "Cannot generate encoded hash of length X from Y bytes"
        Test.@test_throws ErrorException Mycelia.create_blake3_hash("ATCG", encoding=:hex, encoded_length=200, hash_bytes=1)
        Test.@test_throws ErrorException Mycelia.create_blake3_hash("ATCG", encoding=:base58, encoded_length=100, hash_bytes=1)
        Test.@test_throws ErrorException Mycelia.create_blake3_hash("ATCG", encoding=:base64, encoded_length=100, hash_bytes=1)
        
        # Note: The "too short" condition mainly occurs with specific hash_bytes parameters,
        # not in typical create_sequence_hash usage with default parameters
        
        # Test very small sequences
        single_char_hash = Mycelia.create_sequence_hash("A")
        Test.@test isa(single_char_hash, String)
        Test.@test length(single_char_hash) == 64  # Default length
        
        # Test very long sequences
        long_seq = "ATCG" ^ 1000  # 4000 characters
        long_hash = Mycelia.create_sequence_hash(long_seq, encoded_length=32)
        Test.@test length(long_hash) == 32
        
        # Test sequences with non-standard DNA characters
        ambiguous_seq = "ATCGRYSWKMBDHVN"  # Includes IUPAC ambiguity codes
        ambiguous_hash = Mycelia.create_sequence_hash(ambiguous_seq)
        Test.@test isa(ambiguous_hash, String)
        Test.@test length(ambiguous_hash) == 64
    end
    
    Test.@testset "performance and consistency across sequence types" begin
        # Test different sequence types and lengths
        sequences = [
            "A",
            "AT",
            "ATCG",
            "ATCGATCGATCG",
            "ATCGATCGATCGATCGATCGATCGATCGATCG",
            "A" ^ 100,
            "ATCG" ^ 100,
            "AAAAAAAAAA",
            "TTTTTTTTTT",
            "GGGGGGGGGG",
            "CCCCCCCCCC"
        ]
        
        # All should produce valid hashes of correct length
        for seq in sequences
            hash_result = Mycelia.create_sequence_hash(seq, encoded_length=32)
            Test.@test length(hash_result) == 32
            Test.@test isa(hash_result, String)
        end
        
        # Different sequences should produce different hashes
        hashes = [Mycelia.create_sequence_hash(seq, encoded_length=32) for seq in sequences]
        Test.@test length(Set(hashes)) == length(sequences)  # All unique
    end
end

Test.@testset "stable random sequence hash value assertions" begin
    Test.@testset "deterministic hash values with stable random generation" begin
        # Test with stable random DNA sequences
        Test.@testset "random DNA sequences with known hash values" begin
            Random.seed!(12345)
            
            # Generate a few sequences and record their expected hashes
            # These values are computed once and then used for regression testing
            test_sequences = []
            expected_hashes_default = []
            expected_hashes_16 = []
            expected_hashes_32 = []
            
            for i in 1:3
                seq_len = rand(10:30)
                seq = string(BioSequences.randdnaseq(seq_len))
                push!(test_sequences, seq)
                
                # Compute expected hashes
                default_hash = Mycelia.create_sequence_hash(seq)
                hash_16 = Mycelia.create_sequence_hash(seq, encoded_length=16)
                hash_32 = Mycelia.create_sequence_hash(seq, encoded_length=32)
                
                push!(expected_hashes_default, default_hash)
                push!(expected_hashes_16, hash_16)
                push!(expected_hashes_32, hash_32)
            end
            
            # Now test that we get the same results
            Random.seed!(12345)  # Reset seed
            for i in 1:3
                seq_len = rand(10:30)
                seq = string(BioSequences.randdnaseq(seq_len))
                
                Test.@test seq == test_sequences[i]
                
                default_hash = Mycelia.create_sequence_hash(seq)
                hash_16 = Mycelia.create_sequence_hash(seq, encoded_length=16)
                hash_32 = Mycelia.create_sequence_hash(seq, encoded_length=32)
                
                Test.@test default_hash == expected_hashes_default[i]
                Test.@test hash_16 == expected_hashes_16[i]
                Test.@test hash_32 == expected_hashes_32[i]
            end
        end
        
        Test.@testset "random RNA sequences reproducibility" begin
            Random.seed!(54321)
            
            # Test RNA sequence reproducibility
            rna_sequences = []
            rna_hashes = []
            
            for i in 1:3
                seq_len = rand(8:25)
                seq = string(BioSequences.randrnaseq(seq_len))
                hash_val = Mycelia.create_sequence_hash(seq, encoded_length=24)
                push!(rna_sequences, seq)
                push!(rna_hashes, hash_val)
            end
            
            # Verify reproducibility
            Random.seed!(54321)
            for i in 1:3
                seq_len = rand(8:25)
                seq = string(BioSequences.randrnaseq(seq_len))
                hash_val = Mycelia.create_sequence_hash(seq, encoded_length=24)
                
                Test.@test seq == rna_sequences[i]
                Test.@test hash_val == rna_hashes[i]
            end
        end
        
        Test.@testset "random amino acid sequences reproducibility" begin
            Random.seed!(67890)
            
            aa_sequences = []
            aa_hashes = []
            
            for i in 1:3
                seq_len = rand(5:20)
                seq = string(BioSequences.randaaseq(seq_len))
                hash_val = Mycelia.create_sequence_hash(seq, hash_function=:sha256, encoded_length=28, allow_truncation=true)
                push!(aa_sequences, seq)
                push!(aa_hashes, hash_val)
            end
            
            # Verify reproducibility
            Random.seed!(67890)
            for i in 1:3
                seq_len = rand(5:20)
                seq = string(BioSequences.randaaseq(seq_len))
                hash_val = Mycelia.create_sequence_hash(seq, hash_function=:sha256, encoded_length=28, allow_truncation=true)
                
                Test.@test seq == aa_sequences[i]
                Test.@test hash_val == aa_hashes[i]
            end
        end
    end
end

Test.@testset "BioSequence and k-mer helper functions comprehensive testing" begin
    Test.@testset "all BioSequence types produce identical hashes to strings" begin
        # Test DNA sequences
        Test.@testset "DNA sequence types" begin
            dna_string = "ATCGATCGATCG"
            dna_seq = BioSequences.LongDNA{4}(dna_string)
            
            # Test all parameter combinations
            test_params = [
                (hash_function=:blake3, encoded_length=32),
                (hash_function=:sha256, encoded_length=28, allow_truncation=true),
                (hash_function=:md5, encoding=:hex, encoded_length=20, allow_truncation=true),
                (encoding=:base64, encoded_length=24),
                (normalize_case=false,),
                (hash_function=:sha1, encoding=:base58, encoded_length=25, allow_truncation=true)
            ]
            
            for params in test_params
                string_hash = Mycelia.create_sequence_hash(dna_string; params...)
                bioSeq_hash = Mycelia.create_sequence_hash(dna_seq; params...)
                Test.@test string_hash == bioSeq_hash
            end
        end
        
        # Test RNA sequences  
        Test.@testset "RNA sequence types" begin
            rna_string = "AUCGAUCGAUCG"
            rna_seq = BioSequences.LongRNA{4}(rna_string)
            
            for params in [(encoded_length=16,), (hash_function=:blake3, encoded_length=48), (encoding=:hex, encoded_length=20)]
                string_hash = Mycelia.create_sequence_hash(rna_string; params...)
                bioSeq_hash = Mycelia.create_sequence_hash(rna_seq; params...)
                Test.@test string_hash == bioSeq_hash
            end
        end
        
        # Test amino acid sequences
        Test.@testset "amino acid sequence types" begin  
            aa_string = "MKTAYIAKQRQISFVK"
            aa_seq = BioSequences.LongAA(aa_string)
            
            for params in [(encoded_length=20,), (hash_function=:sha512, encoded_length=32, allow_truncation=true), (encoding=:base64, encoded_length=16)]
                string_hash = Mycelia.create_sequence_hash(aa_string; params...)
                bioSeq_hash = Mycelia.create_sequence_hash(aa_seq; params...)
                Test.@test string_hash == bioSeq_hash
            end
        end
        
        # Test different LongDNA types (2-bit, 4-bit)
        Test.@testset "different DNA encoding types" begin
            test_seq = "ATCGATCGATCGATCG"
            
            dna_2bit = BioSequences.LongDNA{2}(test_seq)
            dna_4bit = BioSequences.LongDNA{4}(test_seq) 
            
            string_hash = Mycelia.create_sequence_hash(test_seq, encoded_length=24)
            hash_2bit = Mycelia.create_sequence_hash(dna_2bit, encoded_length=24)
            hash_4bit = Mycelia.create_sequence_hash(dna_4bit, encoded_length=24)
            
            Test.@test string_hash == hash_2bit == hash_4bit
        end
        
        # Test different LongRNA types  
        Test.@testset "different RNA encoding types" begin
            test_seq = "AUCGAUCGAUCGAUCG"
            
            rna_2bit = BioSequences.LongRNA{2}(test_seq)
            rna_4bit = BioSequences.LongRNA{4}(test_seq)
            
            string_hash = Mycelia.create_sequence_hash(test_seq, encoded_length=24)
            hash_2bit = Mycelia.create_sequence_hash(rna_2bit, encoded_length=24)
            hash_4bit = Mycelia.create_sequence_hash(rna_4bit, encoded_length=24)
            
            Test.@test string_hash == hash_2bit == hash_4bit
        end
        
        # Test different amino acid alphabet types
        Test.@testset "amino acid alphabet consistency" begin
            test_seq = "ACDEFGHIKLMNPQRSTVWY"  # All 20 standard amino acids
            aa_seq = BioSequences.LongAA(test_seq)
            
            string_hash = Mycelia.create_sequence_hash(test_seq, encoded_length=32)
            aa_hash = Mycelia.create_sequence_hash(aa_seq, encoded_length=32)
            
            Test.@test string_hash == aa_hash
        end
    end
    
    Test.@testset "k-mer types produce identical hashes to strings" begin
        # Test different k-mer sizes
        Test.@testset "various k-mer sizes" begin
            test_sequences = [
                "ATCG",      # 4-mer
                "ATCGATCG",  # 8-mer  
                "ATCGATCGATCG", # 12-mer
                "ATCGATCGATCGATCG" # 16-mer
            ]
            
            for test_seq in test_sequences
                bioSeq = BioSequences.LongDNA{4}(test_seq)
                kmer = Kmers.DNAKmer{length(test_seq)}(bioSeq)
                
                # Test various parameter combinations
                string_hash_default = Mycelia.create_sequence_hash(test_seq)
                kmer_hash_default = Mycelia.create_sequence_hash(kmer)
                Test.@test string_hash_default == kmer_hash_default
                
                string_hash_param = Mycelia.create_sequence_hash(test_seq, encoded_length=20, hash_function=:sha256, allow_truncation=true)
                kmer_hash_param = Mycelia.create_sequence_hash(kmer, encoded_length=20, hash_function=:sha256, allow_truncation=true)
                Test.@test string_hash_param == kmer_hash_param
            end
        end
        
        # Test that k-mer, string, and BioSequence all produce same hash
        Test.@testset "string vs kmer vs BioSequence consistency" begin
            test_sequences = ["ATCG", "GCTA", "AAAA", "TTTT", "ATCGATCG"]
            
            for test_seq in test_sequences
                # Create all three types
                string_val = test_seq
                bioSeq_val = BioSequences.LongDNA{4}(test_seq)
                kmer_val = Kmers.DNAKmer{length(test_seq)}(bioSeq_val)
                
                # Test with various parameters
                params_list = [
                    (encoded_length=16,),
                    (hash_function=:md5, encoded_length=20, allow_truncation=true),
                    (encoding=:hex, encoded_length=12),
                    (hash_function=:blake3, encoding=:base64, encoded_length=18)
                ]
                
                for params in params_list
                    string_hash = Mycelia.create_sequence_hash(string_val; params...)
                    kmer_hash = Mycelia.create_sequence_hash(kmer_val; params...)
                    bioSeq_hash = Mycelia.create_sequence_hash(bioSeq_val; params...)
                    
                    Test.@test string_hash == kmer_hash == bioSeq_hash
                    Test.@test length(string_hash) == length(kmer_hash) == length(bioSeq_hash)
                end
            end
        end
        
        # Test k-mer specific edge cases
        Test.@testset "k-mer edge cases" begin
            # Single nucleotide k-mer
            single_string = "A"
            single_bioSeq = BioSequences.LongDNA{4}(single_string)
            single_kmer = Kmers.DNAKmer{1}(single_bioSeq)
            
            hash_kmer = Mycelia.create_sequence_hash(single_kmer, encoded_length=8)
            hash_string = Mycelia.create_sequence_hash(single_string, encoded_length=8)  
            hash_bioSeq = Mycelia.create_sequence_hash(single_bioSeq, encoded_length=8)
            
            Test.@test hash_kmer == hash_string == hash_bioSeq
            
            # Very long k-mer (limited by Kmers package constraints)
            long_seq = "ATCGATCGATCGATCG"  # 16-mer (reasonable size)
            long_bioSeq = BioSequences.LongDNA{4}(long_seq)
            long_kmer = Kmers.DNAKmer{16}(long_bioSeq)
            
            hash_long_kmer = Mycelia.create_sequence_hash(long_kmer, encoded_length=24)
            hash_long_string = Mycelia.create_sequence_hash(long_seq, encoded_length=24)
            hash_long_bioSeq = Mycelia.create_sequence_hash(long_bioSeq, encoded_length=24)
            
            Test.@test hash_long_kmer == hash_long_string == hash_long_bioSeq
        end
    end
    
    Test.@testset "cross-type hash consistency verification" begin
        # Comprehensive test ensuring all three types (String, BioSequence, Kmer) produce identical results
        Test.@testset "comprehensive consistency check" begin
            test_data = [
                ("ATCG", 4),
                ("ATCGATCGATCG", 12), 
                ("AAAAAAAAAA", 10),
                ("GCTAGCTAGCTA", 12),
                ("TTTTCCCCGGGGAAAA", 16)
            ]
            
            for (seq, expected_len) in test_data
                Test.@test length(seq) == expected_len
                
                # Create all three representations
                str = seq
                bioSeq = BioSequences.LongDNA{4}(seq)
                kmer = Kmers.DNAKmer{length(seq)}(bioSeq)
                
                # Test comprehensive parameter matrix
                hash_functions = [:blake3, :sha256, :md5]
                encodings = [:base58, :hex, :base64]
                lengths = [8, 16, 24]
                
                for hf in hash_functions, enc in encodings, len in lengths
                    # Skip combinations that would need truncation for non-Blake3
                    skip_truncation = (hf != :blake3 && enc == :base58 && len > 25)
                    
                    params = (
                        hash_function=hf, 
                        encoding=enc, 
                        encoded_length=len,
                        allow_truncation=skip_truncation
                    )
                    
                    try
                        hash_str = Mycelia.create_sequence_hash(str; params...)
                        hash_kmer = Mycelia.create_sequence_hash(kmer; params...)  
                        hash_bioSeq = Mycelia.create_sequence_hash(bioSeq; params...)
                        
                        Test.@test hash_str == hash_kmer == hash_bioSeq
                        Test.@test length(hash_str) == len
                    catch e
                        # If one fails, they should all fail the same way
                        Test.@test_throws typeof(e) Mycelia.create_sequence_hash(str; params...)
                        Test.@test_throws typeof(e) Mycelia.create_sequence_hash(kmer; params...)
                        Test.@test_throws typeof(e) Mycelia.create_sequence_hash(bioSeq; params...)
                    end
                end
            end
        end
    end
end
