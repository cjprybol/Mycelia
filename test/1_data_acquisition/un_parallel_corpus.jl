# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/1_data_acquisition/un_parallel_corpus.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/1_data_acquisition/un_parallel_corpus.jl", "test/1_data_acquisition", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia

Test.@testset "UN Parallel Corpus archive resolution" begin
    archives = Mycelia._resolve_un_archives(["all"], ["txt"])
    expected_count = length(Mycelia.UN_PAIRS) + length(Mycelia.UN_LANGUAGES) + 1
    Test.@test length(archives) == expected_count
    Test.@test "UNv1.0.6way.tar.gz" in archives
    Test.@test "UNv1.0.en-fr.tar.gz" in archives
    Test.@test "UNv1.0-TEI.en.tar.gz" in archives

    pair_archives = Mycelia._resolve_un_archives(["fr-en"], ["txt"])
    Test.@test pair_archives == ["UNv1.0.en-fr.tar.gz"]
end

Test.@testset "UN Parallel Corpus validation" begin
    Test.@test_throws ErrorException Mycelia._resolve_un_archives(String[], ["txt"])
    Test.@test_throws ErrorException Mycelia._resolve_un_archives(["bitext"], ["xml"])
    Test.@test_throws ErrorException Mycelia._resolve_un_archives(["tei"], ["txt"])
    Test.@test_throws ErrorException Mycelia._resolve_un_archives(["zz-yy"], ["txt"])
end
