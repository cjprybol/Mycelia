# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/blast_taxonkit.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/blast_taxonkit.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import DataFrames
import FASTX
import BioSequences

run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

if run_external
    Test.@testset "BLAST DB helpers" begin
        mktempdir() do dir
            fasta = joinpath(dir, "sample.fna")
            record = FASTX.FASTA.Record("seq1", BioSequences.LongDNA{4}("ATGCGT"))
            Mycelia.write_fasta(outfile = fasta, records = [record])
            db_prefix = Mycelia.ensure_blast_db(
                fasta = fasta,
                dbtype = "nucl",
                output_dir = dir,
                db_name = "sample_db",
                title = "sample_db",
                parse_seqids = true,
                force = true
            )
            Test.@test isfile(db_prefix * ".nsq") || isfile(db_prefix * ".00.nsq")
        end
    end

    Test.@testset "BLAST hit taxonkit lineage annotation" begin
        df = DataFrames.DataFrame("subject tax id" => [562, missing, 9606])
        annotated = Mycelia.annotate_blast_hits_with_taxonkit(df)
        Test.@test "lineage" in DataFrames.names(annotated)
        Test.@test annotated[2, "lineage"] === missing
        Test.@test !ismissing(annotated[1, "lineage"])
        Test.@test !isempty(annotated[1, "lineage"])
        Test.@test !ismissing(annotated[3, "lineage"])
        Test.@test !isempty(annotated[3, "lineage"])
    end
else
    @info "Skipping BLAST/taxonkit integration tests; set MYCELIA_RUN_EXTERNAL=true to enable"
end
