# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/aqua.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/aqua.jl", "test", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, ".."))
## end
## using Revise

# Aqua.jl - Auto QUality Assurance for Julia packages
# https://juliatesting.github.io/Aqua.jl/stable/

import Aqua
import Mycelia
import Test
import TOML

Test.@testset "Aqua.jl" begin
    Aqua.test_all(
        Mycelia;
        ambiguities = (broken=true),
        deps_compat = false,
        # BenchmarkTools is used by repo-level benchmarking scripts outside
        # the Mycelia module, so it is intentionally not loaded by the package.
        #
        # XML is a direct dep ONLY so that [compat] can pin the transitive
        # XLSX -> XML resolution away from 0.4.5, which breaks XLSX's
        # precompilation on Julia 1.10. Mycelia never imports it, so Aqua would
        # otherwise correctly flag it as stale. The testset below couples this
        # entry to the pin so it cannot outlive it; see Project.toml's [compat]
        # comment for the full account and the removal condition.
        stale_deps = (ignore = [:BenchmarkTools, :XML],),
        # persistent_tasks test fails due to background tasks spawned by dependencies
        # (e.g., HTTP.jl, Makie.jl, etc.) during package loading - not a Mycelia issue
        persistent_tasks = false
    )
end

# The stale_deps ignore above suppresses a warning that is REAL and CORRECT:
# XML genuinely is a declared-but-unloaded dependency. That suppression is
# justified only while the [compat] pin it exists to serve is still there.
#
# Nothing otherwise couples the two, so the ignore would silently outlive the
# pin -- and a suppression whose justification has evaporated is exactly how a
# detector goes quiet without anyone deciding that it should. This is the
# failing test that pairs with the suppression.
Test.@testset "XML is a transitive pin only" begin
    project = TOML.parsefile(joinpath(pkgdir(Mycelia), "Project.toml"))

    # If the pin is removed, this fails and points at the ignore to remove with
    # it. Matched loosely on the blocked version rather than on the exact bound,
    # so that widening the self-healing half (currently "=0.4.4, 0.4.6 - 0.4")
    # does not spuriously fail -- what matters is that 0.4.5 is still excluded.
    xml_compat = get(project["compat"], "XML", nothing)
    Test.@test xml_compat !== nothing
    Test.@test occursin("0.4.4", something(xml_compat, ""))

    # And the ignore is justified only while Mycelia really does not load XML.
    # If XML ever becomes a genuine dependency, the ignore must go and the dep
    # should be declared for its own sake.
    Test.@test !(:XML in names(Mycelia; all = true, imported = true))
end
