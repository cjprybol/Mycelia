# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/blastdb_integration.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/blastdb_integration.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import Arrow
import DataFrames

if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    Test.@testset "BLAST Database Integration" begin
        Test.@testset "BLAST DB Search Paths" begin
            env_paths = split(get(ENV, "BLASTDB", ""), ":")
            search_paths = filter(isdir, unique(vcat(env_paths, [Mycelia.DEFAULT_BLASTDB_PATH])))
            Test.@test isa(search_paths, Vector)
            Test.@test all(isa.(search_paths, String))
        end
        Test.@testset "Download and Metadata Extraction" begin
            db_name = "ref_viroids_rep_genomes"
            db_dir = mkpath("test-blastdb")
            db_path = Mycelia.download_blast_db(db=db_name, dbdir=db_dir)
            Test.@test isa(db_path, String)
            Test.@test isdir(db_dir)
            metadata = Mycelia.get_blastdb_metadata(blastdb=db_path)
            Test.@test haskey(metadata, "dbtype")
            Test.@test haskey(metadata, "last-updated")
            rm(db_dir, recursive=true)
        end
        Test.@testset "BLAST DB to Arrow and FASTA" begin
            db_name = "ref_viroids_rep_genomes"
            db_dir = mkpath("test-blastdb")
            db_path = Mycelia.download_blast_db(db=db_name, dbdir=db_dir)
            table = Mycelia.blastdb2table(
                blastdb=db_path,
                ALL_FIELDS=false,
                accession=true,
                taxid=true,
                sequence=false,
                sequence_sha256=false
            )
            Test.@test table isa DataFrames.DataFrame
            Test.@test DataFrames.nrow(table) > 0
            arrow_table = Arrow.Table(table)
            Test.@test arrow_table isa Arrow.Table
            fasta_file = Mycelia.blastdb_to_fasta(blastdb=db_path, outfile=joinpath(db_dir, "$(db_name).fna.gz"), force=true, max_cores=1)
            Test.@test isfile(fasta_file)
            rm(db_dir, recursive=true)
        end
        Test.@testset "Filter by Taxid/Entry" begin
            db_name = "ref_viroids_rep_genomes"
            db_dir = mkpath("test-blastdb")
            db_path = Mycelia.download_blast_db(db=db_name, dbdir=db_dir)
            filtered = Mycelia.blastdb2table(
                blastdb=db_path,
                ALL_FIELDS=false,
                taxid=true,
                sequence=false,
                sequence_sha256=false
            )
            Test.@test !isempty(filtered)
            if "taxid" in DataFrames.names(filtered)
                Test.@test any(filtered.taxid .== 12884)
            end
            rm(db_dir, recursive=true)
        end
    end
else
    @info "Skipping BLAST Database Integration tests; external BLAST downloads are opt-in via MYCELIA_RUN_EXTERNAL=true"
end
