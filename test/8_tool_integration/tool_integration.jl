# Tool integration tests
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia

Test.@testset "tool integration" begin
    Test.@testset "padloc" begin
        ecoli_k12_accession = "GCF_000005845.2"
        result = Mycelia.ncbi_genome_download_accession(accession=ecoli_k12_accession, include_string="genome")
        padloc_result = Mycelia.run_padloc(fasta_file = result.genome)
        Test.@test isfile(padloc_result.csv)
        rm(ecoli_k12_accession, recursive=true)
    end
    Test.@testset "ncbi-blast" begin
        blast_database_table = Mycelia.list_blastdbs()
        Test.@test blast_database_table isa DataFrames.DataFrame
    end
end
