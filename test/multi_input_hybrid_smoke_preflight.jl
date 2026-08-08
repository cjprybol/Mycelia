# Pre-discovery driver for the multi-input hybrid private-fixture smoke gates.
#
# `test/runtests.jl` includes this file before filtering external test files.
# Discovery regression tests also execute this file directly in a subprocess.
# The Autocycler-polished gate parses `MYCELIA_AUTOCYCLER_TEST_JOBS` (default 1,
# range 1:4) alongside the shared thread and exact read-chemistry controls.

if !isdefined(@__MODULE__, :_multi_input_hybrid_smoke_prerequisites)
    Base.include(
        @__MODULE__,
        joinpath(@__DIR__, "multi_input_hybrid_smoke_support.jl"),
    )
end

_multi_input_hybrid_smoke_prerequisites(ENV)
_autocycler_smoke_prerequisites(ENV)
