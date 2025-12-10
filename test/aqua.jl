# Aqua.jl - Auto QUality Assurance for Julia packages
# https://juliatesting.github.io/Aqua.jl/stable/

import Aqua
import Mycelia
import Test

Test.@testset "Aqua.jl" begin
    Aqua.test_all(
        Mycelia;
        ambiguities=(broken=true),
        deps_compat=false,
        # persistent_tasks test fails due to background tasks spawned by dependencies
        # (e.g., HTTP.jl, Makie.jl, etc.) during package loading - not a Mycelia issue
        persistent_tasks=false,
    )
end
