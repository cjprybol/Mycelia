"""
    confusion_matrix(true_labels, pred_labels)

Returns the confusion matrix as a Matrix{Int}, row = true, col = predicted.
Also returns the list of unique labels in sorted order and a heatmap plot.
"""
function confusion_matrix(true_labels, pred_labels)
    labels = sort(collect(Set(true_labels) ∪ Set(pred_labels)))
    label_to_index = Dict(l => i for (i, l) in enumerate(labels))
    cm = zeros(Int, length(labels), length(labels))
    for (t, p) in zip(true_labels, pred_labels)
        cm[label_to_index[t], label_to_index[p]] += 1
    end
    # Visualization: heatmap
    cm_plot = Plots.heatmap(
        string.(labels), string.(labels), cm;
        xlabel="Predicted", ylabel="True", title="Confusion Matrix", colorbar_title="Count", c=:blues
    )
    return (cm=cm, labels=labels, plot=cm_plot)
end

"""
    precision_recall_f1(true_labels, pred_labels)

Returns dictionaries mapping each label to its precision, recall, and F1 score.
Also returns macro-averaged (unweighted mean) precision, recall, F1, and a grouped bar plot.
"""
function precision_recall_f1(true_labels, pred_labels)
    cm_out = confusion_matrix(true_labels, pred_labels)
    cm, labels = cm_out.cm, cm_out.labels
    n_labels = length(labels)

    precisions = Dict{Any, Float64}()
    recalls = Dict{Any, Float64}()
    f1s = Dict{Any, Float64}()

    for (i, label) in enumerate(labels)
        tp = cm[i, i]
        fp = Statistics.sum(cm[:, i]) - tp
        fn = Statistics.sum(cm[i, :]) - tp
        precision = tp == 0 && fp == 0 ? 1.0 : tp / (tp + fp + 1e-10)
        recall    = tp == 0 && fn == 0 ? 1.0 : tp / (tp + fn + 1e-10)
        f1        = (precision + recall) == 0 ? 0.0 : 2 * precision * recall / (precision + recall)
        precisions[label] = precision
        recalls[label] = recall
        f1s[label] = f1
    end

    macro_precision = Statistics.mean(collect(values(precisions)))
    macro_recall = Statistics.mean(collect(values(recalls)))
    macro_f1 = Statistics.mean(collect(values(f1s)))

    # Visualization: one bar plot per metric, all labels on x-axis
    label_strs = string.(labels)
    precision_vals = [precisions[label] for label in labels]
    recall_vals = [recalls[label] for label in labels]
    f1_vals = [f1s[label] for label in labels]

    precision_plot = Plots.bar(
        label_strs, precision_vals;
        xlabel = "Label",
        ylabel = "Precision",
        legend = false,
        title = "Precision per Label",
        ylim = (0, 1),
        lw = 0.5,
        framestyle = :box
    )

    recall_plot = Plots.bar(
        label_strs, recall_vals;
        xlabel = "Label",
        ylabel = "Recall",
        legend = false,
        title = "Recall per Label",
        ylim = (0, 1),
        lw = 0.5,
        framestyle = :box
    )

    f1_plot = Plots.bar(
        label_strs, f1_vals;
        xlabel = "Label",
        ylabel = "F1 Score",
        legend = false,
        title = "F1 Score per Label",
        ylim = (0, 1),
        lw = 0.5,
        framestyle = :box
    )

    return (
        precisions=precisions,
        recalls=recalls,
        f1s=f1s,
        macro_precision=macro_precision,
        macro_recall=macro_recall,
        macro_f1=macro_f1,
        precision_plot=precision_plot,
        recall_plot=recall_plot,
        f1_plot=f1_plot
    )
end

"""
    accuracy(true_labels, pred_labels)

Returns the overall accuracy.
"""
function accuracy(true_labels, pred_labels)
    return Statistics.mean(true_labels .== pred_labels)
end

"""
    evaluate_classification(true_labels, pred_labels)

Runs confusion_matrix, precision_recall_f1, and accuracy.
Pretty-prints macro metrics and accuracy.
Returns a named tuple with all results and plots.
"""
function evaluate_classification(true_labels, pred_labels)
    cm_out = confusion_matrix(true_labels, pred_labels)
    prf_out = precision_recall_f1(true_labels, pred_labels)
    acc = accuracy(true_labels, pred_labels)

    macro_f1 = prf_out.macro_f1
    macro_precision = prf_out.macro_precision
    macro_recall = prf_out.macro_recall

    println("==== Evaluation Results ====")
    Printf.@printf("%-18s: %.4f\n", "Macro F1", macro_f1)
    Printf.@printf("%-18s: %.4f\n", "Macro Precision", macro_precision)
    Printf.@printf("%-18s: %.4f\n", "Macro Recall", macro_recall)
    Printf.@printf("%-18s: %.4f\n", "Accuracy", acc)

    return (
        confusion_matrix = cm_out.cm,
        labels = cm_out.labels,
        confusion_matrix_plot = cm_out.plot,
        precisions = prf_out.precisions,
        recalls = prf_out.recalls,
        f1s = prf_out.f1s,
        macro_precision = macro_precision,
        macro_recall = macro_recall,
        macro_f1 = macro_f1,
        precision_plot = prf_out.precision_plot,
        recall_plot = prf_out.recall_plot,
        f1_plot = prf_out.f1_plot,
        accuracy = acc
    )
end

"""
    best_label_mapping(true_labels, pred_labels)

Finds the optimal mapping from predicted labels to true labels using the Hungarian algorithm,
so that the total overlap (confusion matrix diagonal) is maximized.
Returns the remapped predicted labels and the mapping as a Dict.
"""
function best_label_mapping(true_labels, pred_labels)
    # Get sorted unique labels
    labels_true = sort(collect(Set(true_labels)))
    labels_pred = sort(collect(Set(pred_labels)))
    n = max(length(labels_true), length(labels_pred))
    # Build confusion matrix (rows: true, cols: pred)
    cm = zeros(Int, n, n)
    label_to_index_true = Dict(l => i for (i, l) in enumerate(labels_true))
    label_to_index_pred = Dict(l => i for (i, l) in enumerate(labels_pred))
    for (t, p) in zip(true_labels, pred_labels)
        i = label_to_index_true[t]
        j = label_to_index_pred[p]
        cm[i, j] += 1
    end
    # Hungarian algorithm minimizes cost, so use -cm
    assign, _ = Hungarian.hungarian(-cm)
    # According to Hungarian.jl, assign is a 1-based vector with zeros for unassigned
    mapping = Dict()
    for (i, j) in enumerate(assign)
        # Only map if assigned (j > 0) and indices are valid
        if j > 0 && i <= length(labels_true) && j <= length(labels_pred)
            mapping[labels_pred[j]] = labels_true[i]
        end
    end
    # Remap predicted labels using the mapping
    remapped_pred_labels = [haskey(mapping, p) ? mapping[p] : p for p in pred_labels]
    return remapped_pred_labels, mapping
end

#