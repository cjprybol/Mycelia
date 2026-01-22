import Test
import Mycelia

Test.@testset "Plotting Utilities" begin
    Test.@testset "Optimal subsequence length" begin
        Test.@test Mycelia.optimal_subsequence_length(error_rate=0.0) == typemax(Int)
        Test.@test Mycelia.optimal_subsequence_length(error_rate=1.0) == 1
        Test.@test Mycelia.optimal_subsequence_length(error_rate=0.01, threshold=0.99) == 1
    end

    Test.@testset "Jitter" begin
        values = Mycelia.jitter(5.0, 20)
        Test.@test length(values) == 20
        Test.@test maximum(abs.(values .- 5.0)) <= (1.0 / 3.0)
    end

    Test.@testset "Color helpers" begin
        colors = Mycelia.n_maximally_distinguishable_colors(4)
        Test.@test length(colors) == 4

        red = Mycelia.Colors.RGB(1, 0, 0)
        Test.@test Mycelia.merge_colors(red, red) == red
    end

    Test.@testset "Unit conversion helpers" begin
        Test.@test Mycelia.pixels_to_points(12) == 9
        Test.@test Mycelia.points_to_pixels(9) == 12
    end
end
