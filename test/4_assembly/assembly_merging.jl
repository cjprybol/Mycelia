# TODO: Add more tests for merging assemblies, including edge cases and larger graphs.

# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/assembly_merging.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/assembly_merging.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import MetaGraphs
import MetaGraphsNext
import Graphs

Test.@testset "Assembly Merging" begin
    Test.@testset "Contig Merging" begin
        lines = [
            "H\tVN:Z:1.0",
            "S\t1\tACGT\t*\tDP:i:1",
            "S\t2\tCGTA\t*\tDP:i:1",
            "L\t1\t+\t2\t+\t3M",
        ]
        gfa = joinpath(tempdir(), "simple.gfa")
        open(gfa, "w") do io
            for l in lines
                println(io, l)
            end
        end

        g = Mycelia.parse_gfa(gfa)
        Test.@test Graphs.nv(g) == 2
        Test.@test Graphs.ne(g) == 1
        Test.@test length(MetaGraphs.get_prop(g, :records)) == 2

        result = Mycelia.gfa_to_structure_table(gfa)
        Test.@test size(result.contig_table, 1) == 1
        row = result.contig_table[1, :]
        Test.@test row.connected_component == 1
        Test.@test row.contigs == "1,2"
        Test.@test row.is_circular == false
        Test.@test row.is_closed == false
        Test.@test row.lengths == "4,4"
        Test.@test length(result.records) == 2
    end
end
