# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_conversion_errors_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_conversion_errors_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Conversion error handling tests for Rhizomorph graphs

import Test
import Mycelia
import FASTX
import Kmers

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

Test.@testset "Rhizomorph conversion errors" begin
    # AA k-mer graph should not support doublestrand/canonical
    aa_records = [FASTX.FASTA.Record("aa1", "MKVLW")]
    aa_graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(aa_records, 3; dataset_id = "aa")
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        Mycelia.Rhizomorph.convert_to_doublestrand(aa_graph)
    end
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        Mycelia.Rhizomorph.convert_to_canonical(aa_graph)
    end

    # N-gram (String) graph should not support doublestrand/canonical
    strings = ["hello world"]
    ngram_graph = Mycelia.Rhizomorph.build_ngram_graph(strings, 3; dataset_id = "txt")
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        Mycelia.Rhizomorph.convert_to_doublestrand(ngram_graph)
    end
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        Mycelia.Rhizomorph.convert_to_canonical(ngram_graph)
    end
end

println("✓ Rhizomorph conversion error tests completed")
