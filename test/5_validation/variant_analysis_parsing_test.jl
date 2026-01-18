import Test
import Mycelia

Test.@testset "Variant Analysis Parsing" begin
    Test.@testset "RTG eval parsing" begin
        temp_dir = mktempdir()
        eval_path = joinpath(temp_dir, "eval.tsv.gz")
        content = "#Precision\tRecall\tScore\n0.9\t0.6\t10\n0.5\t0.8\t20\n"
        open(eval_path, "w") do io
            gz = Mycelia.CodecZlib.GzipCompressorStream(io)
            write(gz, content)
            close(gz)
        end

        df = Mycelia.parse_rtg_eval_output(eval_path)
        Test.@test Mycelia.DataFrames.nrow(df) == 2
        Test.@test isapprox(df[1, "Precision"], 0.9; atol=1e-6)
        Test.@test isapprox(df[2, "Recall"], 0.8; atol=1e-6)

        empty_path = joinpath(temp_dir, "empty.tsv.gz")
        empty_content = "#Precision\tRecall\tScore\n"
        open(empty_path, "w") do io
            gz = Mycelia.CodecZlib.GzipCompressorStream(io)
            write(gz, empty_content)
            close(gz)
        end

        empty_df = Mycelia.parse_rtg_eval_output(empty_path)
        Test.@test Mycelia.DataFrames.nrow(empty_df) == 0
        Test.@test Mycelia.DataFrames.names(empty_df) == ["Precision", "Recall", "Score"]
    end

    Test.@testset "Evaluation summary" begin
        df = Mycelia.DataFrames.DataFrame(
            Precision=[0.9, 0.5],
            Recall=[0.6, 0.8],
            Score=[10.0, 20.0]
        )
        summary = Mycelia.calculate_evaluation_summary(Dict("snp_roc" => df))
        Test.@test isapprox(summary.snp_roc.max_f1, 0.72; atol=1e-6)
        Test.@test summary.snp_roc.max_precision == 0.9
        Test.@test summary.snp_roc.max_recall == 0.8
        Test.@test summary.snp_roc.optimal_threshold == 10.0
        Test.@test isapprox(summary.snp_roc.auc_estimate, 0.14; atol=1e-6)
    end
end
