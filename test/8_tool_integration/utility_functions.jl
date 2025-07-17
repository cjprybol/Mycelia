import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
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
