# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/mash.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/mash.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import DataFrames
import FASTX

Test.@testset "Mash Tool Integration" begin
    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    Test.@testset "Mash Output Parsing" begin
        dist_file = tempname() * ".tsv"
        screen_file = tempname() * ".tsv"

        write(dist_file, "ref.fasta\tquery.fasta\t0.05\t1e-05\t123/1000\n")
        write(screen_file, "#mash screen output\n0.98\t150/1000\t2\t1e-05\treads.fq\treference.fasta\n")

        dist_table = Mycelia.parse_mash_dist_output(dist_file)
        Test.@test dist_table isa DataFrames.DataFrame
        Test.@test names(dist_table) == ["reference", "query", "distance", "p_value", "matching_hashes"]
        Test.@test dist_table[1, "reference"] == "ref.fasta"

        screen_table = Mycelia.parse_mash_screen_output(screen_file)
        Test.@test screen_table isa DataFrames.DataFrame
        Test.@test names(screen_table) == ["identity", "shared_hashes", "median_multiplicity", "p_value", "query", "reference"]
        Test.@test screen_table[1, "reference"] == "reference.fasta"
    end

    Test.@testset "Mash External Tool Integration" begin
        if run_external
            workdir = mktempdir()

            ref_record = Mycelia.random_fasta_record(moltype=:DNA, seed=7, L=400)
            ref_fasta = joinpath(workdir, "reference.fasta")
            Mycelia.write_fasta(outfile=ref_fasta, records=[ref_record])

            ref_seq = FASTX.sequence(ref_record)
            reads = Mycelia.create_test_reads(String(ref_seq), 10, 0.01)
            reads_fastq = joinpath(workdir, "reads.fastq")
            Mycelia.write_fastq(records=reads, filename=reads_fastq)

            ref_sketch = Mycelia.run_mash_sketch(
                input_files=[ref_fasta],
                outdir=joinpath(workdir, "sketches"),
                k=21,
                s=1000
            ).sketches[1]

            read_sketch = Mycelia.run_mash_sketch(
                input_files=[reads_fastq],
                outdir=joinpath(workdir, "sketches"),
                k=21,
                s=1000,
                r=true,
                min_copies=2
            ).sketches[1]

            dist_result = Mycelia.run_mash_dist(
                reference=ref_sketch,
                query=read_sketch,
                outdir=joinpath(workdir, "dist")
            )
            Test.@test isfile(dist_result.results_tsv)

            screen_result = Mycelia.run_mash_screen(
                reference=ref_sketch,
                query=reads_fastq,
                outdir=joinpath(workdir, "screen"),
                winner_takes_all=true
            )
            Test.@test isfile(screen_result.results_tsv)

            dist_table = Mycelia.parse_mash_dist_output(dist_result.results_tsv)
            Test.@test dist_table isa DataFrames.DataFrame

            screen_table = Mycelia.parse_mash_screen_output(screen_result.results_tsv)
            Test.@test screen_table isa DataFrames.DataFrame
        else
            @info "Skipping mash integration test; external tool execution is opt-in via MYCELIA_RUN_EXTERNAL=true"
        end
    end
end
