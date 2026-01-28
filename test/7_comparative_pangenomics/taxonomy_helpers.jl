# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/taxonomy_helpers.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/taxonomy_helpers.jl", "test/7_comparative_pangenomics", execute=false)'
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

Test.@testset "Taxonomy helpers - column utilities" begin
    df = DataFrames.DataFrame("subject tax id" => [1, 2], "other" => [missing, 3])

    Test.@test Mycelia.first_nonmissing([missing, 3, 4]) == 3
    Test.@test Mycelia.first_nonmissing([missing, missing]) === missing

    Test.@test Mycelia.hascol(df, "subject tax id")
    Test.@test Mycelia.hascol(df, Symbol("subject tax id"))
    Test.@test !Mycelia.hascol(df, "missing_col")

    Test.@test Mycelia.actual_name(df, "subject tax id") == "subject tax id"
    Test.@test Mycelia.actual_name(df, Symbol("subject tax id")) == "subject tax id"
    Test.@test_throws ErrorException Mycelia.actual_name(df, "does_not_exist")

    returned = Mycelia.ensure_as_string!(df, "subject tax id")
    Test.@test returned === df
    Test.@test "subject tax id" in DataFrames.names(df)
end

Test.@testset "Taxonomy helpers - top call summary" begin
    conf_df = DataFrames.DataFrame(
        "query id" => ["q1", "q1", "q2"],
        "rank" => ["species", "species", "genus"],
        "taxid" => [11, 22, 33],
        "taxon" => ["A", "B", "C"],
        "rel_conf" => [0.7, 0.3, 1.0],
    )

    summary = Mycelia.summarize_top_calls(conf_df)
    Test.@test DataFrames.nrow(summary) == 2

    species_row = summary[summary[!, "rank"] .== "species", :]
    Test.@test DataFrames.nrow(species_row) == 1
    Test.@test species_row[1, "taxid"] == 11
    Test.@test species_row[1, "taxon"] == "A"
    Test.@test isapprox(species_row[1, "top_rel_conf"], 0.7; atol=1e-12)
    Test.@test isapprox(species_row[1, "delta_to_second"], 0.4; atol=1e-12)
    Test.@test species_row[1, "n_taxa_considered"] == 2

    missing_df = DataFrames.DataFrame(
        "query id" => ["q1", "q1"],
        "rank" => ["species", "species"],
        "taxid" => [11, 22],
        "taxon" => ["A", "B"],
        "rel_conf" => [missing, missing],
    )
    missing_summary = Mycelia.summarize_top_calls(missing_df)
    Test.@test DataFrames.nrow(missing_summary) == 1
    Test.@test ismissing(missing_summary[1, "taxid"])
    Test.@test ismissing(missing_summary[1, "taxon"])
    Test.@test missing_summary[1, "top_rel_conf"] == 0.0
    Test.@test missing_summary[1, "delta_to_second"] == 0.0
    Test.@test missing_summary[1, "n_taxa_considered"] == 0

    bad_df = DataFrames.DataFrame("query id" => ["q1"], "rank" => ["species"])
    Test.@test_throws ErrorException Mycelia.summarize_top_calls(bad_df)
end

Test.@testset "Taxonomy helpers - top2 indices" begin
    Test.@test Mycelia._top2_indices_and_values(Float64[]) == (0, -Inf, 0, -Inf)

    vals = [missing, 2.0, 1.0]
    i1, v1, i2, v2 = Mycelia._top2_indices_and_values(vals)
    Test.@test (i1, v1) == (2, 2.0)
    Test.@test (i2, v2) == (3, 1.0)
end

Test.@testset "Taxonomy helpers - aggregation and streamlining" begin
    df = DataFrames.DataFrame(
        "template" => ["t1", "t1", "t2"],
        "domain_taxid" => [1, missing, 2],
        "species_taxid" => [missing, 10, 20],
        "relative_alignment_proportion" => [0.5, 0.3, 0.2],
        "percent_identity" => [95.0, 90.0, 88.0],
        "mappingquality" => [60, 55, 50],
    )

    aggregated = Mycelia.aggregate_by_rank_nonmissing(df, ["domain", "species"])
    Test.@test DataFrames.nrow(aggregated) == 4
    Test.@test "rank" in DataFrames.names(aggregated)
    Test.@test "rank_taxid" in DataFrames.names(aggregated)
    Test.@test "total_relative_alignment_proportion" in DataFrames.names(aggregated)
    Test.@test "max_percent_identity" in DataFrames.names(aggregated)
    Test.@test "max_mappingquality" in DataFrames.names(aggregated)

    ranks = Set(aggregated[!, "rank"])
    Test.@test ranks == Set(["1_domain", "2_species"])

    counts = ["A" => 10, "B" => 1, "C" => 1]
    relative = Mycelia.streamline_counts(counts; threshold=0.2, min_len=0)
    Test.@test length(relative) == 2
    Test.@test relative[1] == ("A" => 10)
    Test.@test relative[2] == ("Other" => 2)

    absolute = Mycelia.streamline_counts(counts; threshold=5, min_len=0)
    Test.@test length(absolute) == 2
    Test.@test absolute[1] == ("A" => 10)
    Test.@test absolute[2] == ("Other" => 2)

    untouched = Mycelia.streamline_counts(counts; min_len=10)
    Test.@test untouched == counts
end
