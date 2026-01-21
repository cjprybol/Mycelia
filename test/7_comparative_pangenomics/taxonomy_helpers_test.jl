import Test
import Mycelia

Test.@testset "Taxonomy Helpers" begin
    Test.@testset "Column helpers" begin
        df = Mycelia.DataFrames.DataFrame(
            "subject tax id" => [1, 2],
            "foo" => ["a", "b"]
        )

        Test.@test Mycelia.actual_name(df, "subject tax id") == "subject tax id"
        Test.@test Mycelia.actual_name(df, Symbol("subject tax id")) == "subject tax id"

        Mycelia.ensure_as_string!(df, "foo")
        Test.@test "foo" in Mycelia.DataFrames.names(df)
    end

    Test.@testset "Rank defaults" begin
        ranks = Mycelia.default_ranks()
        Test.@test length(ranks) >= 5
        Test.@test ranks[1].rank == "species"
    end

    Test.@testset "Top-call summarization" begin
        conf_df = Mycelia.DataFrames.DataFrame(
            "query id" => ["q1", "q1", "q1", "q2"],
            "rank" => ["species", "species", "species", "genus"],
            "taxid" => [1, 2, 3, 4],
            "taxon" => ["A", "B", "C", "D"],
            "rel_conf" => [0.6, 0.3, 0.1, 0.9]
        )

        summary = Mycelia.summarize_top_calls(conf_df)
        Test.@test Mycelia.DataFrames.nrow(summary) == 2
        Test.@test summary[1, "top_rel_conf"] >= 0.6
    end

    Test.@testset "Streamline counts" begin
        counts = ["A" => 10, "B" => 1, "C" => 1]
        streamlined = Mycelia.streamline_counts(counts; threshold=0.2, min_len=1)
        Test.@test length(streamlined) == 2
        Test.@test streamlined[2][1] == "Other"
    end

    Test.@testset "Top-2 index helper" begin
        i1, v1, i2, v2 = Mycelia._top2_indices_and_values([missing, 0.2, 0.5])
        Test.@test i1 == 3
        Test.@test v1 == 0.5
        Test.@test i2 == 2
        Test.@test v2 == 0.2
    end
end
