# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/panproteome.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/panproteome.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Panproteome construction tests on toy protein sequences
# import Revise

import Test
import Mycelia
import FASTX
import Kmers
import Graphs
import FASTX
import Kmers

Test.@testset "panproteome construction" begin
    recs = [
        Mycelia.random_fasta_record(moltype = :AA, seed = 1, L = 10),
        Mycelia.random_fasta_record(moltype = :AA, seed = 2, L = 10),
    ]

    counts = Mycelia.count_canonical_kmers(Kmers.AAKmer{2}, recs)
    expected = sum(length(FASTX.sequence(r)) - 1 for r in recs)

    Test.@test sum(values(counts)) == expected
    Test.@test length(counts) <= 400
end
