# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/5_validation/validation.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/validation.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Validation & Quality Assessment tests

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import Kmers
Test.@testset "Assembly Validation" begin
    Test.@testset "Reference-Free Validation" begin
        raw_counts = Dict(
            Kmers.DNAKmer{3}("AAA") => 10,
            Kmers.DNAKmer{3}("CCC") => 20,
            Kmers.DNAKmer{3}("TTT") => 30
        )
        asm_counts = Dict(
            Kmers.DNAKmer{3}("AAA") => 8,
            Kmers.DNAKmer{3}("CCC") => 15,
            Kmers.DNAKmer{3}("GGG") => 5
        )
        qv = Mycelia.kmer_counts_to_merqury_qv(
            raw_data_counts = raw_counts,
            assembly_counts = asm_counts
        )
        expected_P = (2 / 3) ^ (1 / 3)
        expected_QV = -10 * log10(1 - expected_P)
        Test.@test isapprox(qv, expected_QV; atol = 1e-6)
    end

    Test.@testset "Reference-Based Validation" begin
        fastani_text = "query.fasta\treference.fasta\t99.0\t50\t50"
        fastani_file = tempname()
        open(fastani_file, "w") do io
            write(io, fastani_text)
        end
        ani = Mycelia.read_fastani(fastani_file)
        Test.@test size(ani, 1) == 1
        Test.@test ani[1, "query_identifier"] == "query"
        Test.@test ani[1, "reference_identifier"] == "reference"
        Test.@test ani[1, "%_identity"] == 99.0
        Test.@test ani[1, "fragments_mapped"] == 50
        Test.@test ani[1, "total_query_fragments"] == 50
        rm(fastani_file, force = true)

        qualimap_text = """>>>>>>> Coverage per contig
\tcontig1\t100\t500\t5.0\t1.0
\tcontig2\t50\t250\t5.0\t1.0
"""
        qualimap_file = tempname()
        open(qualimap_file, "w") do io
            write(io, qualimap_text)
        end
        cov = Mycelia.parse_qualimap_contig_coverage(qualimap_file)
        Test.@test size(cov, 1) == 2
        Test.@test cov[1, "Contig"] == "contig1"
        Test.@test cov[1, "Length"] == 100
        Test.@test isapprox(cov[1, "% Mapped bases"], 66.6667; atol = 1e-3)
        Test.@test isapprox(cov[2, "% Mapped bases"], 33.3333; atol = 1e-3)
        rm(qualimap_file, force = true)
    end

    Test.@testset "Marker Gene Completeness" begin
        Test.@test true  # placeholder
    end
end
