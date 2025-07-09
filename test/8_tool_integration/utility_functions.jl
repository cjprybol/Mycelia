import Pkg
Pkg.activate("..")
using Revise
using Test
import Mycelia
import Random

@testset "scientific notation" begin
    @test Mycelia.scientific_notation(100) == "1.00e+02"
    @test Mycelia.scientific_notation(1000, precision=3) == "1.000e+03"
end

@testset "byte formatting" begin
    @test Mycelia.bytes_human_readable(1) == "1 byte"
    @test Mycelia.bytes_human_readable(1024^0 + 1024^0) == "2 bytes"
    @test Mycelia.bytes_human_readable(1024^0 + 1024^1) == "1.001 KiB"
    @test Mycelia.bytes_human_readable(1024^1 + 1024^2) == "1.001 MiB"
    @test Mycelia.bytes_human_readable(1024^2 + 1024^3) == "1.001 GiB"
    @test Mycelia.bytes_human_readable(1024^3 + 1024^4) == "1.001 TiB"
    @test Mycelia.bytes_human_readable(1024^4 + 1024^5) == "1.001 PiB"
    @test Mycelia.bytes_human_readable(1024^5 + 1024^6) == "1025.000 PiB"
end

@testset "Estimate memory utilization of dense and sparse Matrices by datatype" begin
    @test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float64, 10_000, 10_000)) == "762.939 MiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float64, 100_000, 100_000)) == "74.506 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float64, 10^5, 10^5, density = 0.5)) == "74.507 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float16, 100_000, 100_000)) == "18.626 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float16, 10^5, 10^5, density = 0.2)) == "18.627 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_dense_matrix_memory(Float16, 10^6, 10^6)) == "1.819 TiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float16, 10^6, 10^6, density = 0.1)) == "931.330 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float16, 10^6, 10^6, density = 0.2)) == "1.819 TiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Float64, 10^5, 10^5, density = 0.2)) == "29.803 GiB"
    # uses default byte size for Float64
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(10^5, 10^5, density = 0.2)) == "29.803 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(Int64, 10^5, 10^5, density = 0.2)) == "29.803 GiB"
    @test Mycelia.bytes_human_readable(Mycelia.estimate_sparse_matrix_memory(UInt64, 10^5, 10^5, density = 0.2)) == "29.803 GiB"
end

@testset "check matrix fits in memory" begin
    Mycelia.system_overview()
    assessed_memory_needs = Mycelia.check_matrix_fits_in_memory(
        Mycelia.estimate_dense_matrix_memory(10^3, 10^3)
    )
    @test assessed_memory_needs.will_fit_available == (assessed_memory_needs.bytes_needed <= assessed_memory_needs.free_memory)
    @test assessed_memory_needs.will_fit_total == (assessed_memory_needs.bytes_needed <= assessed_memory_needs.total_memory)
end