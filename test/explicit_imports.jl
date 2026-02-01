# ExplicitImports.jl - Ensure explicit imports are used consistently
# https://github.com/ericphanson/ExplicitImports.jl
#
# This test enforces Cameron's Julia style preference:
# - Always use `import`, never `using`
# - Always use fully-qualified function calls (e.g., DataFrames.DataFrame())

import ExplicitImports
import Mycelia
import Test

Test.@testset "ExplicitImports.jl" begin
    # Check that all exports are explicitly qualified
    # This enforces the "import, never using" convention
    # Note: allow_unanalyzable=(Mycelia,) is needed because Mycelia uses dynamic includes
    # (loading files from a list), which ExplicitImports cannot statically analyze
    Test.@test ExplicitImports.check_no_implicit_imports(Mycelia; allow_unanalyzable = (Mycelia,)) ===
               nothing

    # Check for stale explicit imports (imported but not used)
    Test.@test ExplicitImports.check_no_stale_explicit_imports(Mycelia; allow_unanalyzable = (Mycelia,)) ===
               nothing
end
