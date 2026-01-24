# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers_legacy.jl")'
# ```
#
# To turn one of the legacy group files into a Jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_legacy_short_read.jl", "test/4_assembly", execute=false)'
# ```
#
# NOTE: Legacy notebook-style tests are split into group files:
# - `test/4_assembly/third_party_assemblers_legacy_short_read.jl`
# - `test/4_assembly/third_party_assemblers_legacy_metagenome.jl`
# - `test/4_assembly/third_party_assemblers_legacy_long_read.jl`
# - `test/4_assembly/third_party_assemblers_legacy_hylight.jl`
#
# The lightweight test-suite entrypoint is `test/4_assembly/third_party_assemblers.jl`.

include(joinpath(@__DIR__, "third_party_assemblers_legacy_short_read.jl"))
include(joinpath(@__DIR__, "third_party_assemblers_legacy_metagenome.jl"))
include(joinpath(@__DIR__, "third_party_assemblers_legacy_long_read.jl"))
# NOTE: HyLight legacy tests are currently disabled pending fixture refinement.
# include(joinpath(@__DIR__, "third_party_assemblers_legacy_hylight.jl"))
