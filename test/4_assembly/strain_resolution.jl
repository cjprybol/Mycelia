# Strain Resolution tests
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia

Test.@testset "Strain Resolution" begin
    Test.@testset "Strain-aware Reassembly" begin
        Test.@test Mycelia.ks(min=5, max=7) == [5, 7]
    end
end
