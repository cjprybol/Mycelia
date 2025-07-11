# BLAST Database Integration and Metadata Tests
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia

Test.@testset "BLAST Database Integration" begin
    Test.@testset "BLAST DB Search Paths" begin
        search_paths = Mycelia.blastdb_search_paths()
        Test.@test isa(search_paths, Vector)
        Test.@test all(isa.(search_paths, String))
    end
    Test.@testset "Download and Metadata Extraction" begin
        db_name = "ref_viroids_rep_genomes"
        db_dir = mkpath("test-blastdb")
        db_result = Mycelia.download_blastdb(db_name, outdir=db_dir)
        Test.@test isfile(db_result.archive)
        Test.@test isfile(db_result.metadata)
        metadata = Mycelia.read_blastdb_metadata(db_result.metadata)
        Test.@test haskey(metadata, :title)
        Test.@test haskey(metadata, :taxid)
        rm(db_dir, recursive=true)
    end
    Test.@testset "BLAST DB to Arrow and FASTA" begin
        db_name = "ref_viroids_rep_genomes"
        db_dir = mkpath("test-blastdb")
        db_result = Mycelia.download_blastdb(db_name, outdir=db_dir)
        arrow_table = Mycelia.blastdb_to_arrow(db_result.archive)
        Test.@test arrow_table isa Arrow.Table
        fasta_file = Mycelia.blastdb_to_fasta(db_result.archive, outdir=db_dir)
        Test.@test isfile(fasta_file)
        rm(db_dir, recursive=true)
    end
    Test.@testset "Filter by Taxid/Entry" begin
        db_name = "ref_viroids_rep_genomes"
        db_dir = mkpath("test-blastdb")
        db_result = Mycelia.download_blastdb(db_name, outdir=db_dir)
        filtered = Mycelia.filter_blastdb_by_taxid(db_result.archive, taxid=12884)
        Test.@test !isempty(filtered)
        rm(db_dir, recursive=true)
    end
end
