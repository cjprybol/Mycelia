import Test
import Mycelia

Test.@testset "QC Metrics" begin
    Test.@testset "Robust statistics" begin
        values = [1.0, 2.0, 3.0, 4.0]
        cv = Mycelia.robust_cv(values)
        Test.@test cv > 0.0

        threshold = Mycelia.robust_threshold(values; k=2.0)
        Test.@test threshold > 0.0

        Test.@test Mycelia.robust_cv(Float64[]) == 0.0
        Test.@test Mycelia.robust_threshold(Float64[]) == 0.0
    end

    Test.@testset "Classification metrics" begin
        true_labels = ["a", "b", "a", "b"]
        pred_labels = ["a", "b", "b", "b"]

        cm_out = Mycelia.confusion_matrix(true_labels, pred_labels)
        Test.@test size(cm_out.cm) == (2, 2)
        Test.@test length(cm_out.labels) == 2

        prf = Mycelia.precision_recall_f1(true_labels, pred_labels)
        Test.@test prf.macro_precision <= 1.0
        Test.@test prf.macro_recall <= 1.0
        Test.@test prf.macro_f1 <= 1.0

        acc = Mycelia.accuracy(true_labels, pred_labels)
        Test.@test acc == 0.75

        eval_out = Mycelia.evaluate_classification(true_labels, pred_labels; verbose=false)
        Test.@test eval_out.accuracy == 0.75

        mapped, mapping = Mycelia.best_label_mapping(true_labels, pred_labels)
        Test.@test length(mapped) == length(pred_labels)
        Test.@test !isempty(mapping)
    end

    Test.@testset "Presence metrics" begin
        truth = ["A", "B", "C"]
        pred = ["B", "C", "D"]
        scores = Mycelia.presence_precision_recall_f1(truth, pred)
        Test.@test scores.tp == 2
        Test.@test scores.fp == 1
        Test.@test scores.fn == 1
    end
end
