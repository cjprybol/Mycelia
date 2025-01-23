using Test

@testset "My Project Tests" begin
    @testset "Test Case 1" begin
        @test 1 + 1 == 2
    end

    @testset "Test Case 2" begin
        @test "hello" * " world" == "hello world"
    end
end