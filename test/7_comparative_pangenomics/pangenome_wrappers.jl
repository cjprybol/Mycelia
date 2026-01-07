# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/pangenome_wrappers.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/pangenome_wrappers.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia

Test.@testset "Pangenome wrapper validation" begin
    # PGGB should fail fast on missing genome files
    genomes = ["missing1.fasta", "missing2.fasta"]
    outdir = mktempdir()
    Test.@test_throws ErrorException Mycelia.construct_pangenome_pggb(genomes, outdir)

    # Cactus should catch mismatched names before touching the filesystem
    genome_files = ["g1.fasta", "g2.fasta"]
    genome_names = ["sample1"]
    Test.@test_throws ErrorException Mycelia.construct_pangenome_cactus(
        genome_files, genome_names, outdir, "sample1")

    # vg deconstruct validation
    Test.@test_throws ErrorException Mycelia.call_variants_from_pggb_graph("missing.gfa", "ref")
end
