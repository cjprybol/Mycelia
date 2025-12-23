# NCBI genetic code lookup tests
#
# julia --project=. --color=yes -e 'include("test/6_annotation/genetic_code.jl")'

import Test
import Mycelia
import MetaGraphsNext
import Graphs

Test.@testset "NCBI genetic code lookup" begin
    ncbi_taxonomy = MetaGraphsNext.MetaGraph(
        MetaGraphsNext.DiGraph();
        label_type=Int,
        vertex_data_type=Dict{Symbol, Any},
        edge_data_type=Nothing,
        graph_data=Dict{Symbol, Any}(),
    )

    ncbi_taxonomy[1] = Dict{Symbol, Any}(
        :tax_id => 1,
        :genetic_code_id => 11,
        :mitochondrial_genetic_code_id => 5,
    )
    ncbi_taxonomy[2] = Dict{Symbol, Any}(
        :tax_id => 2,
        :genetic_code_id => 0,
        :mitochondrial_genetic_code_id => 0,
    )
    ncbi_taxonomy[3] = Dict{Symbol, Any}(:tax_id => 3)

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
