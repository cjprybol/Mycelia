# NCBI genetic code lookup tests
#
# julia --project=. --color=yes -e 'include("test/6_annotation/genetic_code.jl")'

import Test
import Mycelia
import MetaGraphs
import Graphs

Test.@testset "NCBI genetic code lookup" begin
    ncbi_taxonomy = MetaGraphs.MetaDiGraph(3)

    MetaGraphs.set_prop!(ncbi_taxonomy, 1, :tax_id, 1)
    MetaGraphs.set_prop!(ncbi_taxonomy, 2, :tax_id, 2)
    MetaGraphs.set_prop!(ncbi_taxonomy, 3, :tax_id, 3)

    MetaGraphs.set_prop!(ncbi_taxonomy, 1, :genetic_code_id, 11)
    MetaGraphs.set_prop!(ncbi_taxonomy, 1, :mitochondrial_genetic_code_id, 5)
    MetaGraphs.set_prop!(ncbi_taxonomy, 2, :genetic_code_id, 0)
    MetaGraphs.set_prop!(ncbi_taxonomy, 2, :mitochondrial_genetic_code_id, 0)

    MetaGraphs.set_prop!(ncbi_taxonomy, :node_2_taxid_map, [1, 2, 3])

    Graphs.add_edge!(ncbi_taxonomy, 2, 1)
    Graphs.add_edge!(ncbi_taxonomy, 3, 2)

    Test.@test Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 1) == 11
    Test.@test Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 2) == 11
    Test.@test Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 2; inherit=false) == 0
    Test.@test Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 3) == 11
    Test.@test Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 3; type=:mitochondrial) == 5
    Test.@test Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 42) === missing
    Test.@test_throws ArgumentError Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, 1; type=:nuclear)
end
