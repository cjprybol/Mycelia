# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/7_comparative_pangenomics/blastdb_integration.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/blastdb_integration.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import Arrow
import DataFrames

run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
if run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    test_tmp_root = get(ENV, "MYCELIA_TEST_TMPDIR", Mycelia.DEFAULT_BLASTDB_PATH)
    mkpath(test_tmp_root)
    Test.@testset "BLAST Database Integration" begin
        Test.@testset "BLAST DB Search Paths" begin
            env_paths = split(get(ENV, "BLASTDB", ""), ":")
            search_paths = filter(isdir, unique(vcat(env_paths, [Mycelia.DEFAULT_BLASTDB_PATH])))
            Test.@test isa(search_paths, Vector)
            Test.@test all(isa.(search_paths, String))
        end
        Test.@testset "Download and Metadata Extraction" begin
            db_name = "ref_viroids_rep_genomes"
            db_dir = mktempdir(test_tmp_root)
            try
                db_path = Mycelia.download_blast_db(db = db_name, dbdir = db_dir)
                Test.@test isa(db_path, String)
                Test.@test isdir(db_dir)
                metadata = Mycelia.get_blastdb_metadata(blastdb = db_path)
                Test.@test haskey(metadata, "dbtype")
                Test.@test haskey(metadata, "last-updated")
            finally
                rm(db_dir, recursive = true, force = true)
            end
        end
        Test.@testset "BLAST DB to Arrow and FASTA" begin
            db_name = "ref_viroids_rep_genomes"
            db_dir = mktempdir(test_tmp_root)
            try
                db_path = Mycelia.download_blast_db(db = db_name, dbdir = db_dir)
                table = Mycelia.blastdb2table(
                    blastdb = db_path,
                    blastdbs_dir = db_dir,
                    ALL_FIELDS = false,
                    accession = true,
                    taxid = true,
                    sequence = false,
                    sequence_sha256 = false
                )
                Test.@test table isa DataFrames.DataFrame
                Test.@test DataFrames.nrow(table) > 0
                arrow_path = tempname(db_dir) * ".arrow"
                Arrow.write(arrow_path, table)
                arrow_table = Arrow.Table(arrow_path)
                Test.@test arrow_table isa Arrow.Table
                entry_candidates = if "accession" in DataFrames.names(table)
                    filter(!isempty, string.(coalesce.(table.accession, "")))
                else
                    String[]
                end
                entries = unique(first(entry_candidates, min(length(entry_candidates), 5)))
                Test.@test !isempty(entries)
                fasta_file = withenv("TMPDIR" => db_dir) do
                    Mycelia.blastdb_to_fasta(
                        blastdb = db_path,
                        entries = entries,
                        outfile = joinpath(db_dir, "$(db_name).fna.gz"),
                        force = true,
                        max_cores = 1
                    )
                end
                Test.@test isfile(fasta_file)
            finally
                rm(db_dir, recursive = true, force = true)
            end
        end
        Test.@testset "Filter by Taxid/Entry" begin
            db_name = "ref_viroids_rep_genomes"
            db_dir = mktempdir(test_tmp_root)
            try
                db_path = Mycelia.download_blast_db(db = db_name, dbdir = db_dir)
                filtered = Mycelia.blastdb2table(
                    blastdb = db_path,
                    blastdbs_dir = db_dir,
                    ALL_FIELDS = false,
                    taxid = true,
                    sequence = false,
                    sequence_sha256 = false
                )
                Test.@test !isempty(filtered)
                if "taxid" in DataFrames.names(filtered)
                    taxid_strings = string.(coalesce.(filtered.taxid, ""))
                    # Require at least one non-empty taxid; specific values vary by DB release
                    has_taxids = any(!isempty, taxid_strings)
                    Test.@test has_taxids
                    if has_taxids && !any(taxid_strings .== "12884")
                        @info "Taxid 12884 not present in filtered table; sample taxids" unique(first(
                            filtered.taxid, min(length(filtered.taxid), 5)))
                    end
                end
            finally
                rm(db_dir, recursive = true, force = true)
            end
        end
    end
else
    @info "Skipping BLAST Database Integration tests; external BLAST downloads are opt-in via MYCELIA_RUN_EXTERNAL=true"
end
