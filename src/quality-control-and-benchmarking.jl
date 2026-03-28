"""
    parse_fastplong_json(json_file::String)

Robustly parse a fastplong JSON report, handling various nesting structures.
Returns `nothing` if parsing fails, preventing crashes in batch operations.
"""
function parse_fastplong_json(json_file::String)
    if !isfile(json_file)
        @warn "File not found: $json_file"
        return nothing
    end

    try
        data = JSON.parsefile(json_file)
        sample_id = replace(basename(json_file),
            r"\.json$" => "",
            r"\.fastplong.*$" => "", # Cleanup common suffixes
            r"\.hifi_reads.*$" => ""
        )

        # 1. Locate the summary section (handle nesting variations)
        summary_node = nothing
        if haskey(data, "summary")
            summary_node = data["summary"]
        elseif haskey(data, "before_filtering")
            summary_node = data # Data is at root
        end

        if isnothing(summary_node)
            @warn "Could not find 'summary' or 'before_filtering' in $json_file"
            return nothing
        end

        # 2. Extract Metrics
        rows = []
        for stage in ["before_filtering", "after_filtering"]
            if haskey(summary_node, stage)
                section = summary_node[stage]
                push!(rows,
                    (
                        sample_id = sample_id,
                        stage = stage,
                        total_reads = get(section, "total_reads", 0),
                        total_bases = get(section, "total_bases", 0),
                        q20_bases = get(section, "q20_bases", 0),
                        q30_bases = get(section, "q30_bases", 0),
                        # Handle naming variations (fastp uses mean_length, fastplong uses read_mean_length)
                        mean_length = get(section, "read_mean_length", get(section, "mean_length", 0.0)),
                        n50 = get(section, "read_n50", get(section, "n50", 0.0)),
                        gc_content = get(section, "gc_content", 0.0)
                    ))
            end
        end

        if isempty(rows)
            @warn "No stats found in $json_file"
            return nothing
        end

        df = DataFrames.DataFrame(rows)

        # 3. Safe Calculation of Derived Metrics
        # Use safe division to avoid errors if total_bases is 0
        df.q20_percent = [r.total_bases > 0 ? r.q20_bases / r.total_bases * 100 : 0.0
                          for r in eachrow(df)]
        df.q30_percent = [r.total_bases > 0 ? r.q30_bases / r.total_bases * 100 : 0.0
                          for r in eachrow(df)]
        df.yield_gb = df.total_bases ./ 1e9

        # 4. Extract Histograms (if available) for advanced plotting
        # fastplong often puts these at the root level, e.g., "read_length_histogram"
        histograms = Dict()
        if haskey(data, "read_length_histogram")
            histograms[:read_length] = data["read_length_histogram"]
        end

        return (; sample_id, summary = df, json_data = data, histograms)

    catch e
        @warn "Error parsing $json_file: $e"
        return nothing
    end
end

"""
    aggregate_fastplong_reports(parsed_results::Vector)

Aggregate a list of pre-parsed fastplong results into a single master DataFrame.
Useful if you have already run `parse_fastplong_json` (e.g., via pmap or a loop).
"""
function aggregate_fastplong_reports(parsed_results::Vector)
    # Filter out any 'nothing' entries from failed parses
    valid_results = filter(!isnothing, parsed_results)

    if isempty(valid_results)
        return DataFrames.DataFrame()
    end

    # Combine all individual summary DataFrames
    all_summaries = vcat([r.summary for r in valid_results]...)

    # By default, we usually want the 'after_filtering' stats for batch QC
    # But returning the full table allows the user to filter later if they want 'before_filtering'
    return all_summaries
end

"""
    aggregate_fastplong_reports(input_paths::Vector{String})

Aggregate QC metrics from a list of file paths. Parses them on the fly.
"""
function aggregate_fastplong_reports(input_paths::Vector{String})
    # Parse all files
    parsed = map(parse_fastplong_json, input_paths)
    return aggregate_fastplong_reports(parsed)
end

"""
    aggregate_fastplong_reports(input_dir::String; recursive=false)

Aggregate QC metrics from all .json files in a directory.
"""
function aggregate_fastplong_reports(input_dir::String; recursive = false)
    if recursive
        files = [joinpath(root, f) for (root, dirs, files) in walkdir(input_dir)
                 for f in files if endswith(f, ".json")]
    else
        files = filter(f -> endswith(f, ".json"), readdir(input_dir, join = true))
    end
    return aggregate_fastplong_reports(files)
end

"""
    robust_cv(x::AbstractVector{<:Real}; epsilon::Float64 = eps(Float64))::Float64

Compute a robust coefficient of variation for a vector of distances x
using MAD as the robust scale and the median as the robust location.

Returns 0.0 if the robust location is too close to zero to avoid numerical instability.
"""
function robust_cv(x::AbstractVector{<:Real}; epsilon::Float64 = eps(Float64))::Float64
    xx = collect(skipmissing(x))
    isempty(xx) && return 0.0
    loc = Statistics.median(xx)
    sca = StatsBase.mad(xx; normalize = true)
    denom = max(abs(loc), epsilon)
    return sca / denom
end

"""
    robust_threshold(values::AbstractVector{<:Real}; k::Float64 = 3.0)::Float64

Compute a robust high-side threshold using the median plus k times the normalized MAD.
"""
function robust_threshold(values::AbstractVector{<:Real}; k::Float64 = 3.0)::Float64
    vv = collect(skipmissing(values))
    isempty(vv) && return 0.0
    m = Statistics.median(vv)
    s = StatsBase.mad(vv; normalize = true)
    return m + k * s
end

"""
    filter_genome_outliers(
        df_in::DataFrames.AbstractDataFrame;
        group_col::AbstractString = "reference_entity",
        path_col::AbstractString = "normalized_fastx_path",
        k_local_median::Float64 = 3.0,
        k_local_rcv::Float64 = 3.0,
        k_global_median::Float64 = 3.0,
        k_global_rcv::Float64 = 3.0,
        combine_mode::Symbol = :or,
    )

Compute pairwise distance matrices per group, derive per-genome robust statistics from
row-wise distances, estimate global and local robust thresholds, and flag/filter genomes
that violate the automated checks.

Returns a named tuple with filtered, flagged, metrics, and thresholds.
"""
function filter_genome_outliers(
        df_in::DataFrames.AbstractDataFrame;
        group_col::AbstractString = "reference_entity",
        path_col::AbstractString = "normalized_fastx_path",
        k_local_median::Float64 = 3.0,
        k_local_rcv::Float64 = 3.0,
        k_global_median::Float64 = 3.0,
        k_global_rcv::Float64 = 3.0,
        combine_mode::Symbol = :and        # combine_mode::Symbol = :or,
)
    df = DataFrames.DataFrame(df_in)
    n = DataFrames.nrow(df)
    if n == 0
        return (
            filtered = df,
            flagged = DataFrames.DataFrame(),
            metrics = DataFrames.DataFrame(),
            thresholds = (global_median = 0.0, global_rcv = 0.0)
        )
    end

    # Track original row order to map per-genome metrics back to df
    DataFrames.insertcols!(df, 1, :__rowid => collect(1:n))

    groups = DataFrames.groupby(df, group_col)

    # Storage for global threshold estimation
    global_medians = Float64[]
    global_rcvs = Float64[]

    # Per-group cache to compute local thresholds later
    group_results = Vector{NamedTuple{(:rowids, :group_label, :medians, :rcvs),
        Tuple{Vector{Int}, Any, Vector{Float64}, Vector{Float64}}}}()

    # Compute per-genome metrics per group
    for sub in groups
        m = DataFrames.nrow(sub)
        if m <= 1
            continue
        end

        paths = Vector{String}(sub[!, path_col])
        D = Mycelia.pairwise_mash_distance_matrix(input_fasta_list = paths)
        M = Array{Float64}(D)

        medians = Vector{Float64}(undef, m)
        rcvs = Vector{Float64}(undef, m)

        for i in 1:m
            # distances from genome i to others in the same group (exclude self-distance)
            di = Vector{Float64}(undef, m - 1)
            idx = 1
            for j in 1:m
                if j == i
                    continue
                end
                di[idx] = M[i, j]
                idx += 1
            end
            medians[i] = Statistics.median(di)
            rcvs[i] = robust_cv(di)
        end

        append!(global_medians, medians)
        append!(global_rcvs, rcvs)

        push!(group_results,
            (
                rowids = Vector{Int}(sub[!, :__rowid]),
                group_label = sub[1, group_col],
                medians = medians,
                rcvs = rcvs
            ))
    end

    # Compute global thresholds
    T_global_median = robust_threshold(global_medians; k = k_global_median)
    T_global_rcv = robust_threshold(global_rcvs; k = k_global_rcv)

    # Prepare column vectors for metrics_df
    rowids_v = Int[]
    group_vals_v = Any[]
    median_v = Float64[]
    rcv_v = Float64[]
    flag_l_m_v = Bool[]
    flag_l_r_v = Bool[]
    flag_g_m_v = Bool[]
    flag_g_r_v = Bool[]
    flag_out_v = Bool[]
    T_local_m_v = Float64[]
    T_local_r_v = Float64[]

    flagged_rowids = Set{Int}()

    # Evaluate local thresholds per group and flag rows
    for gr in group_results
        # Local thresholds for this group
        T_local_median = robust_threshold(gr.medians; k = k_local_median)
        T_local_rcv = robust_threshold(gr.rcvs; k = k_local_rcv)

        for (idx, rowid) in enumerate(gr.rowids)
            m_i = gr.medians[idx]
            rcv_i = gr.rcvs[idx]

            flag_l_m = m_i > T_local_median
            flag_l_r = rcv_i > T_local_rcv
            flag_g_m = m_i > T_global_median
            flag_g_r = rcv_i > T_global_rcv

            local_trigger = (flag_l_m || flag_l_r)
            global_trigger = (flag_g_m || flag_g_r)

            flag = (combine_mode === :or) ? (local_trigger || global_trigger) :
                   (combine_mode === :and) ? (local_trigger && global_trigger) :
                   (local_trigger || global_trigger)

            if flag
                push!(flagged_rowids, rowid)
            end

            push!(rowids_v, rowid)
            push!(group_vals_v, gr.group_label)
            push!(median_v, m_i)
            push!(rcv_v, rcv_i)
            push!(flag_l_m_v, flag_l_m)
            push!(flag_l_r_v, flag_l_r)
            push!(flag_g_m_v, flag_g_m)
            push!(flag_g_r_v, flag_g_r)
            push!(flag_out_v, flag)
            push!(T_local_m_v, T_local_median)
            push!(T_local_r_v, T_local_rcv)
        end
    end

    # Build metrics_df from proper column vectors
    metrics_df = DataFrames.DataFrame(
        :__rowid => rowids_v,
        Symbol(group_col) => group_vals_v,
        :median_distance => median_v,
        :robust_cv => rcv_v,
        :flag_local_median => flag_l_m_v,
        :flag_local_rcv => flag_l_r_v,
        :flag_global_median => flag_g_m_v,
        :flag_global_rcv => flag_g_r_v,
        :flag_outlier => flag_out_v,
        :T_local_median => T_local_m_v,
        :T_local_rcv => T_local_r_v
    )

    # Drop the group column on the right side to avoid duplicate names during join
    metrics_join = DataFrames.select(metrics_df, DataFrames.Not(Symbol(group_col)))

    df_joined = DataFrames.leftjoin(df, metrics_join, on = :__rowid => :__rowid)

    is_flagged = hasproperty(df_joined, :flag_outlier) ?
                 coalesce.(df_joined[!, :flag_outlier], false) :
                 falses(DataFrames.nrow(df_joined))
    flagged_df = df_joined[is_flagged, :]
    filtered_df = df_joined[.!is_flagged, :]

    # Clean helper column
    if :__rowid in DataFrames.names(flagged_df)
        DataFrames.select!(flagged_df, DataFrames.Not(:__rowid))
    end
    if :__rowid in DataFrames.names(filtered_df)
        DataFrames.select!(filtered_df, DataFrames.Not(:__rowid))
    end
    if :__rowid in DataFrames.names(metrics_df)
        DataFrames.select!(metrics_df, DataFrames.Not(:__rowid))
    end

    thresholds = (global_median = T_global_median, global_rcv = T_global_rcv)

    return (filtered = filtered_df, flagged = flagged_df,
        metrics = metrics_df, thresholds = thresholds)
end

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
        xlabel = "Predicted", ylabel = "True", title = "Confusion Matrix", colorbar_title = "Count", c = :blues
    )
    return (cm = cm, labels = labels, plot = cm_plot)
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
        recall = tp == 0 && fn == 0 ? 1.0 : tp / (tp + fn + 1e-10)
        f1 = (precision + recall) == 0 ? 0.0 : 2 * precision * recall / (precision + recall)
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
        precisions = precisions,
        recalls = recalls,
        f1s = f1s,
        macro_precision = macro_precision,
        macro_recall = macro_recall,
        macro_f1 = macro_f1,
        precision_plot = precision_plot,
        recall_plot = recall_plot,
        f1_plot = f1_plot
    )
end

"""
    presence_precision_recall_f1(truth_ids, predicted_ids)

Compute presence/absence precision, recall, and F1 for two sets of identifiers.
"""
function presence_precision_recall_f1(truth_ids::AbstractVector, predicted_ids::AbstractVector)
    truth_set = Set(truth_ids)
    pred_set = Set(predicted_ids)
    tp = length(intersect(truth_set, pred_set))
    fp = length(setdiff(pred_set, truth_set))
    fn = length(setdiff(truth_set, pred_set))
    precision = tp == 0 ? 0.0 : tp / (tp + fp)
    recall = tp == 0 ? 0.0 : tp / (tp + fn)
    f1 = (precision + recall) == 0 ? 0.0 : 2 * precision * recall / (precision + recall)
    return (tp = tp, fp = fp, fn = fn, precision = precision, recall = recall, f1 = f1)
end

"""
    evaluate_taxonomy_presence_metrics(reference_table, truth_ids, predicted_ids; ranks=String[])

Compute presence/absence precision/recall/F1 for taxonomy ranks derived from a reference table.
Returns a DataFrame with one row per rank.
"""
function evaluate_taxonomy_presence_metrics(reference_table::DataFrames.DataFrame,
        truth_ids::AbstractVector,
        predicted_ids::AbstractVector;
        ranks::Vector{String} = String[])
    if isempty(ranks)
        ranks = ["species", "genus", "family", "order", "class",
            "phylum", "kingdom", "realm", "domain"]
    end
    ranks = filter(col -> col in DataFrames.names(reference_table), ranks)
    metrics = DataFrames.DataFrame(
        rank = String[], precision = Float64[], recall = Float64[],
        f1 = Float64[], tp = Int[], fp = Int[], fn = Int[])
    truth_set = Set(String.(truth_ids))
    pred_set = Set(String.(predicted_ids))
    for rank in ranks
        truth_taxa = unique(skipmissing(reference_table[reference_table.sequence_id .∈ Ref(truth_set), rank]))
        pred_taxa = unique(skipmissing(reference_table[reference_table.sequence_id .∈ Ref(pred_set), rank]))
        scores = presence_precision_recall_f1(String.(truth_taxa), String.(pred_taxa))
        DataFrames.push!(metrics,
            (rank = rank, precision = scores.precision, recall = scores.recall,
                f1 = scores.f1, tp = scores.tp, fp = scores.fp, fn = scores.fn))
    end
    return metrics
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
function evaluate_classification(true_labels, pred_labels; verbose::Bool = true)
    cm_out = confusion_matrix(true_labels, pred_labels)
    prf_out = precision_recall_f1(true_labels, pred_labels)
    acc = accuracy(true_labels, pred_labels)

    macro_f1 = prf_out.macro_f1
    macro_precision = prf_out.macro_precision
    macro_recall = prf_out.macro_recall

    if verbose
        println("==== Evaluation Results ====")
        Printf.@printf("%-18s: %.4f\n", "Macro F1", macro_f1)
        Printf.@printf("%-18s: %.4f\n", "Macro Precision", macro_precision)
        Printf.@printf("%-18s: %.4f\n", "Macro Recall", macro_recall)
        Printf.@printf("%-18s: %.4f\n", "Accuracy", acc)
    end

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

"""
    assess_assembly_quality(contigs_file)

Assess basic assembly quality metrics from a FASTA file.

Calculates standard assembly quality metrics including contig count, total length,
and N50 statistic for assembly evaluation.

# Arguments
- `contigs_file`: Path to FASTA file containing assembly contigs

# Returns
- Tuple of (n_contigs, total_length, n50, l50)
  - `n_contigs`: Number of contigs in the assembly
  - `total_length`: Total length of all contigs in base pairs
  - `n50`: N50 statistic (length of shortest contig in the set covering 50% of assembly)
  - `l50`: L50 statistic (number of contigs needed to reach 50% of assembly length)

# Example
```julia
n_contigs, total_length, n50, l50 = assess_assembly_quality("assembly.fasta")
println("Assembly has \$n_contigs contigs, \$total_length bp total, N50=\$n50, L50=\$l50")
```

# See Also
- `assess_assembly_kmer_quality`: For k-mer based assembly quality assessment
"""
function assess_assembly_quality(contigs_file)
    contigs = []

    open(FASTX.FASTA.Reader, contigs_file) do reader
        for record in reader
            push!(contigs, length(FASTX.FASTA.sequence(record)))
        end
    end

    n_contigs = length(contigs)
    total_length = sum(contigs)

    # Calculate N50 and L50
    sorted_lengths = sort(contigs, rev = true)
    cumsum_lengths = cumsum(sorted_lengths)
    target_length = total_length / 2
    n50_idx = findfirst(x -> x >= target_length, cumsum_lengths)
    n50 = n50_idx !== nothing ? sorted_lengths[n50_idx] : 0
    l50 = n50_idx !== nothing ? n50_idx : length(contigs)

    return n_contigs, total_length, n50, l50
end

"""
    compute_nl_stats(sorted_lengths, fraction) -> Tuple{Int,Int}

Compute Nx/Lx statistics for a sorted list of contig lengths.
"""
function compute_nl_stats(sorted_lengths, fraction)
    total = sum(sorted_lengths)
    target = total * fraction
    cumulative = 0
    for (idx, length_val) in enumerate(sorted_lengths)
        cumulative += length_val
        if cumulative >= target
            return length_val, idx
        end
    end
    return 0, length(sorted_lengths)
end

"""
    assembly_metrics(contigs_file) -> Union{NamedTuple,Nothing}

Return standard assembly metrics including N50/N90 and GC content.
"""
function assembly_metrics(contigs_file)
    if contigs_file === nothing || !isfile(contigs_file)
        return nothing
    end

    n_contigs, total_length, n50, l50 = Mycelia.assess_assembly_quality(contigs_file)
    records = open(FASTX.FASTA.Reader, contigs_file) do reader
        collect(reader)
    end
    lengths = [length(FASTX.sequence(record)) for record in records]
    sorted_lengths = sort(lengths, rev = true)
    n90, l90 = compute_nl_stats(sorted_lengths, 0.9)

    gc_count = 0
    total_bases = 0
    for record in records
        seq = FASTX.sequence(record)
        total_bases += length(seq)
        for base in seq
            if base == BioSequences.DNA_G || base == BioSequences.DNA_C
                gc_count += 1
            end
        end
    end
    gc_content = total_bases == 0 ? 0.0 : (gc_count / total_bases)
    largest_contig = isempty(sorted_lengths) ? 0 : first(sorted_lengths)

    return (
        n_contigs = n_contigs,
        total_length = total_length,
        n50 = n50,
        l50 = l50,
        n90 = n90,
        l90 = l90,
        largest_contig = largest_contig,
        gc_content = gc_content
    )
end

"""
    contig_gc_outliers(contigs_file; min_length=500)

Detect contigs with GC fraction more than two standard deviations from the mean.
"""
function contig_gc_outliers(contigs_file; min_length = 500)
    outliers = []
    if contigs_file === nothing || !isfile(contigs_file)
        return outliers
    end

    gc_values = Float64[]
    records = open(FASTX.FASTA.Reader, contigs_file) do reader
        collect(reader)
    end
    for record in records
        seq = FASTX.sequence(record)
        if length(seq) < min_length
            continue
        end
        gc_count = count(base -> base == BioSequences.DNA_G || base == BioSequences.DNA_C, seq)
        push!(gc_values, gc_count / length(seq))
    end

    if isempty(gc_values)
        return outliers
    end

    mean_gc = Statistics.mean(gc_values)
    std_gc = Statistics.std(gc_values)

    for record in records
        seq = FASTX.sequence(record)
        if length(seq) < min_length
            continue
        end
        gc_count = count(base -> base == BioSequences.DNA_G || base == BioSequences.DNA_C, seq)
        gc = gc_count / length(seq)
        if abs(gc - mean_gc) > 2 * std_gc
            push!(outliers, (id = FASTX.identifier(record), length = length(seq), gc = gc))
        end
    end

    return outliers
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Evaluate genome assembly quality by comparing k-mer distributions between assembled sequences and raw observations.

# Arguments
- `assembly`: Input assembled sequences to evaluate
- `observations`: Raw sequencing data for comparison
- `ks::Vector{Int}`: Vector of k-mer sizes to analyze (default: k=17 to 23)

# Returns
DataFrame containing quality metrics for each k-mer size:
- `k`: K-mer length used
- `cosine_distance`: Cosine similarity between k-mer distributions
- `js_divergence`: Jensen-Shannon divergence between distributions  
- `qv`: MerQury-style quality value score
"""
function assess_assembly_kmer_quality(;
        assembly, observations, ks::Vector{Int} = filter(x -> 17 <= x <= 23, Mycelia.ks()))
    results = DataFrames.DataFrame()
    @show ks
    ProgressMeter.@showprogress for k in ks
        @show k
        @info "counting assembly kmers..."
        assembled_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, assembly)
        # assembled_canonical_kmer_counts_file = Mycelia.jellyfish_count(fastx=assembly, k=k, canonical=true)
        # @info "loading assembly kmer counts..."
        # assembled_canonical_kmer_counts_table = Mycelia.load_jellyfish_counts(assembled_canonical_kmer_counts_file)
        # sort!(assembled_canonical_kmer_counts_table, "kmer")
        # assembled_canonical_kmer_counts = OrderedCollections.OrderedDict(row["kmer"] => row["count"] for row in DataFrames.eachrow(assembled_canonical_kmer_counts_table))

        @info "counting observation kmers..."
        observed_canonical_kmer_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, observations)
        # observed_canonical_kmer_counts_file = Mycelia.jellyfish_count(fastx=observations, k=k, canonical=true)
        # @info "loading observation kmer counts..."
        # observed_canonical_kmer_counts_table = Mycelia.load_jellyfish_counts(observed_canonical_kmer_counts_file)
        # sort!(observed_canonical_kmer_counts_table, "kmer")
        # observed_canonical_kmer_counts = OrderedCollections.OrderedDict(row["kmer"] => row["count"] for row in DataFrames.eachrow(observed_canonical_kmer_counts_table))
        cosine_distance = kmer_counts_to_cosine_similarity(
            observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        js_divergence = kmer_counts_to_js_divergence(
            observed_canonical_kmer_counts, assembled_canonical_kmer_counts)
        qv = kmer_counts_to_merqury_qv(raw_data_counts = observed_canonical_kmer_counts,
            assembly_counts = assembled_canonical_kmer_counts)
        push!(results, (; k, cosine_distance, js_divergence, qv))
    end
    return results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate assembly Quality Value (QV) score using the Merqury method.

Estimates base-level accuracy by comparing k-mer distributions between raw sequencing
data and assembly. Higher QV scores indicate better assembly quality.

# Arguments
- `raw_data_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}`: K-mer counts from raw sequencing data
- `assembly_counts::AbstractDict{Kmers.DNAKmer{k,N}, Int}`: K-mer counts from assembly

# Returns
- `Float64`: Quality Value score in Phred scale (-10log₁₀(error rate))

# Method
QV is calculated using:
1. Ktotal = number of unique kmers in assembly
2. Kshared = number of kmers shared between raw data and assembly
3. P = (Kshared/Ktotal)^(1/k) = estimated base-level accuracy
4. QV = -10log₁₀(1-P)

# Reference
Rhie et al. "Merqury: reference-free quality, completeness, and phasing assessment
for genome assemblies" Genome Biology (2020)
"""
function kmer_counts_to_merqury_qv(;
        raw_data_counts::AbstractDict{Kmers.DNAKmer{k, N}, Int},
        assembly_counts::AbstractDict{Kmers.DNAKmer{k, N}, Int}) where {k, N}
    # Ktotal = # of kmers found in assembly
    Ktotal = length(keys(assembly_counts))
    # Kshared = # of shared kmers between assembly and readset
    Kshared = length(intersect(keys(raw_data_counts), keys(assembly_counts)))
    # probability_base_in_assembly_correct
    P = (Kshared/Ktotal)^(1/k)
    # # Error rate
    E = 1-P
    QV = -10log10(E)
    # return (;P, E, QV)
    return QV
end

"""
    compute_merqury_qv(assembly_file, reads_file; k=21)

Compute Merqury-style QV for an assembly using raw read data.
"""
function compute_merqury_qv(assembly_file, reads_file; k = 21)
    if assembly_file === nothing || reads_file === nothing
        return missing
    end
    if !isfile(assembly_file) || !isfile(reads_file)
        return missing
    end

    assembly_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, assembly_file)
    read_counts = Mycelia.count_canonical_kmers(Kmers.DNAKmer{k}, reads_file)
    return Mycelia.kmer_counts_to_merqury_qv(raw_data_counts = read_counts, assembly_counts = assembly_counts)
end

# ---------- CheckM, CheckM2, CheckV Functions ----------

"""
    setup_checkv(; db_dir::String=joinpath(homedir(), "workspace", ".checkv"))

Install CheckV via Bioconda and set up its database.

CheckV is a tool for assessing the quality and completeness of viral genomes.

# Arguments
- `db_dir`: Directory to store CheckV database (default: ~/.checkv)

# Example
```julia
setup_checkv()
```
"""
function setup_checkv(; db_dir::String = joinpath(homedir(), "workspace", ".checkv"))
    add_bioconda_env("checkv")
    if !isdir(db_dir)
        mkpath(db_dir)
        run(`$(CONDA_RUNNER) run -n checkv checkv download_database $(db_dir)`)
    end
    return db_dir
end

"""
    setup_checkm(; db_dir::String=joinpath(homedir(), "workspace", ".checkm"))

Install CheckM via Bioconda and set up its database.

CheckM is a tool for assessing the quality and completeness of bacterial genomes.

# Arguments
- `db_dir`: Directory to store CheckM database (default: ~/.checkm)

# Example
```julia
setup_checkm()
```
"""
function setup_checkm(; db_dir::String = joinpath(homedir(), "workspace", ".checkm"))
    add_bioconda_env("checkm-genome")
    if !isdir(db_dir)
        mkpath(db_dir)
        run(`$(CONDA_RUNNER) run -n checkm-genome checkm data setRoot $(db_dir)`)
    end
    return db_dir
end

"""
    setup_checkm2(; db_dir::String=joinpath(homedir(), "workspace", ".checkm2"))

Install CheckM2 via Bioconda and set up its database.

CheckM2 is a rapid tool for assessing the quality and completeness of bacterial genomes.

# Arguments
- `db_dir`: Directory to store CheckM2 database (default: ~/.checkm2)

# Example
```julia
setup_checkm2()
```
"""
function setup_checkm2(; db_dir::String = joinpath(homedir(), "workspace", ".checkm2"))
    add_bioconda_env("checkm2")
    if !isdir(db_dir)
        mkpath(db_dir)
        try
            run(`$(CONDA_RUNNER) run -n checkm2 checkm2 database --download --path $(db_dir)`)
        catch e
            @warn "CheckM2 database download failed; CheckM2 will not work until the database is installed." exception=e
        end
    end
    return db_dir
end

"""
    setup_ncbi_fcs_gx(; db_dir::Union{String,Nothing}=joinpath(homedir(), "workspace", ".ncbi-fcs-gx", "gxdb"), manifest_url::Union{String,Nothing}=nothing, download_db::Bool=false)

Install NCBI FCS-GX via Bioconda and optionally download its GX database.

The Bioconda package is `ncbi-fcs-gx`. The GX database is large (~470 GB) and may
require a high-memory filesystem for optimal performance (/dev/shm on systems with
>= 512 GiB RAM).

# Arguments
- `db_dir`: Directory to store the GX database (only used when `download_db=true`)
- `manifest_url`: Manifest URL for `sync_files.py` (required when `download_db=true`)
- `download_db`: If true, download the GX DB using `sync_files.py`
"""
function setup_ncbi_fcs_gx(;
        db_dir::Union{String, Nothing} = joinpath(homedir(), "workspace", ".ncbi-fcs-gx", "gxdb"),
        manifest_url::Union{String, Nothing} = nothing,
        download_db::Bool = false)
    add_bioconda_env("ncbi-fcs-gx")
    if isnothing(db_dir)
        return db_dir
    end

    if download_db
        if manifest_url === nothing
            error("manifest_url is required when download_db=true")
        end
        mkpath(db_dir)
        run(`$(CONDA_RUNNER) run --live-stream -n ncbi-fcs-gx sync_files.py get --mft $(manifest_url) --dir $(db_dir)`)
    end
    return db_dir
end

function resolve_ncbi_fcs_gx_db_path(gx_db::String)
    if isfile(gx_db)
        if endswith(gx_db, ".gxi")
            return gx_db
        end
        error("GX database file must end with .gxi: $(gx_db)")
    end

    if isdir(gx_db)
        gxi_files = filter(path -> endswith(path, ".gxi"), readdir(gx_db; join = true))
        if isempty(gxi_files)
            error("No .gxi file found in GX database directory: $(gx_db)")
        elseif length(gxi_files) > 1
            error("Multiple .gxi files found in GX database directory: $(gx_db)")
        end
        return gxi_files[1]
    end

    parent_dir = dirname(gx_db)
    if isdir(parent_dir)
        prefix = basename(gx_db)
        gxi_files = filter(
            path -> startswith(basename(path), prefix) && endswith(path, ".gxi"),
            readdir(parent_dir; join = true)
        )
        if isempty(gxi_files)
            error("GX database not found for prefix: $(gx_db)")
        elseif length(gxi_files) > 1
            error("Multiple GX database files match prefix: $(gx_db)")
        end
        return gxi_files[1]
    end

    error("GX database path does not exist: $(gx_db)")
end

"""
    run_checkv(fasta_file::String; outdir::String=fasta_file * "_checkv", db_dir::String=joinpath(homedir(), "workspace", ".checkv"))

Run CheckV on a single genome FASTA file.

CheckV assesses the quality and completeness of viral genomes.

# Arguments
- `fasta_file`: Path to single FASTA file (can be gzipped)
- `outdir`: Output directory for CheckV results (default: fasta_file * "_checkv")
- `threads`: Number of threads to use (default: all available CPU threads)
- `db_dir`: CheckV database directory (default: ~/.checkv)

# Returns
- Named tuple with fields:
  - `outdir`: Output directory path
  - `complete_genomes`: Path to complete_genomes.tsv
  - `completeness`: Path to completeness.tsv
  - `contamination`: Path to contamination.tsv
  - `proviruses`: Path to proviruses.fna
  - `quality_summary`: Path to quality_summary.tsv
  - `viruses`: Path to viruses.fna

# Example
```julia
result = run_checkv("genome.fasta")
println("Quality summary: ", result.quality_summary)
```
"""
function run_checkv(fasta_file::String; outdir::String = fasta_file * "_checkv",
        db_dir::String = joinpath(homedir(), "workspace", ".checkv"),
        threads::Int = get_default_threads())
    if !isfile(fasta_file)
        error("Input file does not exist: $(fasta_file)")
    end

    if !occursin(FASTA_REGEX, fasta_file)
        error("Input file does not match FASTA format: $(fasta_file)")
    end

    # Define expected output files
    output_files = (
        complete_genomes = joinpath(outdir, "complete_genomes.tsv"),
        completeness = joinpath(outdir, "completeness.tsv"),
        contamination = joinpath(outdir, "contamination.tsv"),
        proviruses = joinpath(outdir, "proviruses.fna"),
        quality_summary = joinpath(outdir, "quality_summary.tsv"),
        viruses = joinpath(outdir, "viruses.fna")
    )

    # Check if output directory exists and all expected files are present
    if !(isdir(outdir) && all(isfile, output_files))
        # @warn "CheckV output directory already exists with all expected files. Skipping analysis: $(outdir)"
        # else
        latest_db = last(readdir(Mycelia.setup_checkv(), join = true))
        run(`$(CONDA_RUNNER) run -n checkv checkv end_to_end $(fasta_file) $(outdir) -d $(latest_db) -t $(threads)`)
    end

    # Clean up temporary directory if it exists
    tmp_dir = joinpath(outdir, "tmp")
    Mycelia.cleanup_directory(tmp_dir, verbose = false, force = true)

    return (outdir = outdir, output_files...)
end

"""
    run_ncbi_fcs(
        fasta_file::String;
        tax_id::Int,
        outdir::String=fasta_file * "_fcs_gx",
        out_basename::Union{String,Nothing}=nothing,
        gx_db::Union{String,Nothing}=nothing,
        manifest_url::Union{String,Nothing}=nothing,
        download_db::Bool=false,
        species::Union{String,Nothing}=nothing,
        split_fasta::Bool=true,
        div::Union{String,Nothing}=nothing,
        mask_transposons::Union{Bool,Nothing}=nothing,
        allow_same_species::Union{Bool,Nothing}=nothing,
        ignore_same_kingdom::Bool=false,
        action_report::Bool=true,
        save_hits::Bool=false,
        generate_logfile::Bool=false,
        phone_home_label::Union{String,Nothing}=nothing,
        debug::Bool=false,
    )

Run NCBI FCS-GX contamination screening on a single FASTA file.

Notes:
- Requires the GX database (`.gxi` + companion files). Pass `download_db=true`
  with `manifest_url`, or prepare the DB separately (e.g., via s5cmd download).
- The tool may query NCBI for taxonomy division information if `--div` is not set.
"""
function run_ncbi_fcs(
        fasta_file::String;
        tax_id::Int,
        outdir::String = fasta_file * "_fcs_gx",
        out_basename::Union{String, Nothing} = nothing,
        gx_db::Union{String, Nothing} = nothing,
        manifest_url::Union{String, Nothing} = nothing,
        download_db::Bool = false,
        species::Union{String, Nothing} = nothing,
        split_fasta::Bool = true,
        div::Union{String, Nothing} = nothing,
        mask_transposons::Union{Bool, Nothing} = nothing,
        allow_same_species::Union{Bool, Nothing} = nothing,
        ignore_same_kingdom::Bool = false,
        action_report::Bool = true,
        save_hits::Bool = false,
        generate_logfile::Bool = false,
        phone_home_label::Union{String, Nothing} = nothing,
        debug::Bool = false
)
    if !isfile(fasta_file)
        error("Input file does not exist: $(fasta_file)")
    end

    if !occursin(FASTA_REGEX, fasta_file)
        error("Input file does not match FASTA format: $(fasta_file)")
    end

    default_db_dir = joinpath(homedir(), "workspace", ".ncbi-fcs-gx", "gxdb")
    gx_db_path = isnothing(gx_db) ? default_db_dir : gx_db
    download_dir = isdir(gx_db_path) ? gx_db_path : dirname(gx_db_path)

    setup_ncbi_fcs_gx(
        db_dir = download_db ? download_dir : nothing,
        manifest_url = manifest_url,
        download_db = download_db
    )

    resolved_gx_db = resolve_ncbi_fcs_gx_db_path(gx_db_path)

    mkpath(outdir)

    cmd_parts = [
        "run_gx.py",
        "--fasta", fasta_file,
        "--tax-id", string(tax_id),
        "--out-dir", outdir,
        "--gx-db", resolved_gx_db,
        "--split-fasta", string(split_fasta),
        "--action-report", string(action_report),
        "--save-hits", string(save_hits),
        "--generate-logfile", string(generate_logfile)
    ]

    if !isnothing(out_basename)
        push!(cmd_parts, "--out-basename", out_basename)
    end

    if !isnothing(species)
        push!(cmd_parts, "--species", species)
    end

    if !isnothing(div)
        push!(cmd_parts, "--div", div)
    end

    if !isnothing(mask_transposons)
        push!(cmd_parts, "--mask-transposons", string(mask_transposons))
    end

    if !isnothing(allow_same_species)
        push!(cmd_parts, "--allow-same-species", string(allow_same_species))
    end

    if ignore_same_kingdom
        push!(cmd_parts, "--ignore-same-kingdom")
    end

    if !isnothing(phone_home_label)
        push!(cmd_parts, "--phone-home-label", phone_home_label)
    end

    if debug
        push!(cmd_parts, "--debug")
    end

    run(`$(CONDA_RUNNER) run --live-stream -n ncbi-fcs-gx $cmd_parts`)

    base = out_basename === nothing ?
           splitext(basename(fasta_file))[1] * "." * string(tax_id) : out_basename
    out_prefix = joinpath(outdir, base)
    taxonomy_report = out_prefix * ".taxonomy.rpt"
    report = action_report ? out_prefix * ".fcs_gx_report.txt" : nothing
    summary = generate_logfile ? out_prefix * ".summary.txt" : nothing
    hits = save_hits ? out_prefix * ".hits.tsv.gz" : nothing

    if !isfile(taxonomy_report)
        error("FCS-GX taxonomy report missing: $(taxonomy_report)")
    end

    if action_report && !isfile(report)
        error("FCS-GX action report missing: $(report)")
    end

    return (outdir = outdir, out_basename = base, taxonomy_report = taxonomy_report,
        report = report, summary = summary, hits = hits)
end

"""
    run_ncbi_fcs_clean_genome(
        fasta_file::String;
        action_report::String,
        outdir::String=fasta_file * "_fcs_gx_clean",
        output_fasta::Union{String,Nothing}=nothing,
        contam_fasta::Union{String,Nothing}=nothing,
        min_seq_len::Int=200,
    )

Apply an FCS-GX action report to produce a cleaned FASTA and optional contamination FASTA.
"""
function run_ncbi_fcs_clean_genome(
        fasta_file::String;
        action_report::String,
        outdir::String = fasta_file * "_fcs_gx_clean",
        output_fasta::Union{String, Nothing} = nothing,
        contam_fasta::Union{String, Nothing} = nothing,
        min_seq_len::Int = 200
)
    if !isfile(fasta_file)
        error("Input file does not exist: $(fasta_file)")
    end

    if !isfile(action_report)
        error("Action report does not exist: $(action_report)")
    end

    setup_ncbi_fcs_gx(db_dir = nothing)
    mkpath(outdir)

    output_fasta = output_fasta === nothing ? joinpath(outdir, "clean.fasta") : output_fasta
    contam_fasta = contam_fasta === nothing ? joinpath(outdir, "contam.fasta") :
                   contam_fasta

    cmd_parts = [
        "gx",
        "clean-genome",
        "--input", fasta_file,
        "--action-report", action_report,
        "--output", output_fasta,
        "--contam-fasta-out", contam_fasta,
        "--min-seq-len", string(min_seq_len)
    ]

    run(`$(CONDA_RUNNER) run --live-stream -n ncbi-fcs-gx $cmd_parts`)

    if !isfile(output_fasta)
        error("Cleaned FASTA not generated: $(output_fasta)")
    end

    return (
        outdir = outdir, cleaned_fasta = output_fasta, contamination_fasta = contam_fasta)
end

"""
    run_checkm(input_path::String; outdir::String=input_path * "_checkm", db_dir::String=joinpath(homedir(), "workspace", ".checkm"), extension::String="fasta")

Run CheckM on directory containing FASTA files.

CheckM requires a directory of genome files as input.

# Arguments
- `input_path`: Path to directory containing FASTA files
- `outdir`: Output directory for CheckM results (default: input_path * "_checkm")
- `db_dir`: CheckM database directory (default: ~/.checkm)
- `extension`: File extension for genomes (default: "fasta")
- `threads`: Number of threads to use (default: all available CPU threads)

# Example
```julia
run_checkm("./genomes/")
```
"""
function run_checkm(input_path::String; outdir::String = input_path * "_checkm",
        db_dir::String = joinpath(homedir(), "workspace", ".checkm"),
        extension::String = "fasta", threads::Int = get_default_threads())
    setup_checkm(db_dir = db_dir)

    if !isdir(input_path)
        error("CheckM requires a directory of genome files as input: $(input_path)")
    end

    fasta_files = find_fasta_files(input_path)
    if isempty(fasta_files)
        error("No FASTA files found in directory: $(input_path)")
    end

    # Check for potential performance issues
    num_genomes = length(fasta_files)
    if num_genomes > 1000
        @warn "Processing $(num_genomes) genomes. CheckM documentation suggests splitting large datasets (>1000 genomes) into smaller groups for better performance and memory management."
    end

    # Check available memory (convert from bytes to GB)
    available_memory_gb = Sys.total_memory() / (1024^3)
    if available_memory_gb < 64
        @warn "Available system memory: $(round(available_memory_gb, digits=1)) GB. CheckM may require 64+ GB of RAM for optimal performance with large datasets. Consider reducing the number of genomes or running on a system with more memory."
    end

    run(`$(CONDA_RUNNER) run -n checkm-genome checkm lineage_wf $(input_path) $(outdir) -x $(extension) --data $(db_dir) -t $(threads)`)

    return outdir
end

"""
    run_checkm2(input_path::String; outdir::String=input_path * "_checkm2", db_dir::String=joinpath(homedir(), "workspace", ".checkm2"))

Run CheckM2 on FASTA file(s) or directory containing FASTA files.

# Arguments
- `input_path`: Path to FASTA file or directory containing FASTA files
- `outdir`: Output directory for CheckM2 results (default: input_path * "_checkm2")
- `threads`: Number of threads to use (default: all available CPU threads)
- `db_dir`: CheckM2 database directory (default: ~/.checkm2)

# Returns
A named tuple with the following fields:
- `outdir`: Output directory path
- `quality_report`: Path to quality_report.tsv
- `log_file`: Path to checkm2.log
- `diamond_results`: Path to DIAMOND_RESULTS.tsv
- `protein_file`: Path to the single .faa file
"""
function run_checkm2(input_path::String; outdir::String = input_path * "_checkm2",
        db_dir::String = joinpath(homedir(), "workspace", ".checkm2"),
        threads::Int = get_default_threads())
    # /home/jupyter/workspace/.checkm2/CheckM2_database/uniref100.KO.1.dmnd
    # latest_db = first(readdir(last(readdir(setup_checkm2(db_dir=db_dir), join=true)), join=true))
    db_path = joinpath(db_dir, "CheckM2_database", "uniref100.KO.1.dmnd")
    @assert isfile(db_path)

    fasta_files = find_fasta_files(input_path)

    if isempty(fasta_files)
        error("No FASTA files found in: $(input_path)")
    end

    # Define expected output file paths
    quality_report = joinpath(outdir, "quality_report.tsv")
    log_file = joinpath(outdir, "checkm2.log")
    diamond_results = joinpath(outdir, "diamond_output", "DIAMOND_RESULTS.tsv")
    protein_files_dir = joinpath(outdir, "protein_files")

    # Check if all expected files already exist
    files_exist = isfile(quality_report) && isfile(log_file) && isfile(diamond_results) &&
                  isdir(protein_files_dir)

    if files_exist
        # Check for exactly one .faa file
        faa_files = filter(f -> endswith(f, ".faa"), readdir(protein_files_dir, join = true))

        if length(faa_files) == 1
            protein_file = first(faa_files)
            # @warn "All expected CheckM2 output files already exist. Skipping execution and returning existing results." outdir

            return (
                outdir = outdir,
                quality_report = quality_report,
                log_file = log_file,
                diamond_results = diamond_results,
                protein_file = protein_file
            )
        end
    end

    # CheckM2 can take comma-separated fasta files as input
    input_str = join(fasta_files, ",")

    run(`$(CONDA_RUNNER) run -n checkm2 checkm2 predict -i $(input_str) -o $(outdir) --database_path $(db_path) --threads $(threads) --force`)

    # Verify output files were created
    if !isdir(protein_files_dir)
        error("Protein files directory not found: $(protein_files_dir)")
    end

    faa_files = filter(f -> endswith(f, ".faa"), readdir(protein_files_dir, join = true))

    if length(faa_files) != 1
        error("Expected exactly 1 .faa file, but found $(length(faa_files)) in $(protein_files_dir)")
    end

    protein_file = first(faa_files)

    return (
        outdir = outdir,
        quality_report = quality_report,
        log_file = log_file,
        diamond_results = diamond_results,
        protein_file = protein_file
    )
end

"""
    run_checkm2_list(fasta_files::Vector{String}; outdir::String=normalized_current_datetime() * "_checkm2", db_dir::String=joinpath(homedir(), "workspace", ".checkm2"))

Run CheckM2 on a list of FASTA files.

CheckM2 can automatically handle mixed lists of gzipped and non-gzipped files when given a list.

# Arguments
- `fasta_files`: Vector of FASTA file paths (can be mixed gzipped and non-gzipped)
- `outdir`: Output directory for CheckM2 results (default: normalized_current_datetime() * "_checkm2")
- `threads`: Number of threads to use (default: all available CPU threads)
- `db_dir`: CheckM2 database directory (default: ~/.checkm2)

# Example
```julia
files = ["genome1.fasta.gz", "genome2.fasta", "genome3.fasta.gz"]
run_checkm2_list(files)
```
"""
function run_checkm2_list(fasta_files::Vector{String};
        outdir::String = normalized_current_datetime() * "_checkm2",
        db_dir::String = joinpath(homedir(), "workspace", ".checkm2"),
        threads::Int = get_default_threads())
    # /home/jupyter/workspace/.checkm2/CheckM2_database/uniref100.KO.1.dmnd
    # latest_db = first(readdir(last(readdir(setup_checkm2(db_dir=db_dir), join=true)), join=true))
    db_path = joinpath(db_dir, "CheckM2_database", "uniref100.KO.1.dmnd")
    @assert isfile(db_path)

    if isempty(fasta_files)
        error("No FASTA files provided")
    end

    # Validate that all files exist and match FASTA format
    for file in fasta_files
        if !isfile(file)
            error("File does not exist: $(file)")
        end
        if !occursin(FASTA_REGEX, file)
            error("File does not match FASTA format: $(file)")
        end
    end

    # CheckM2 can take comma-separated fasta files as input and handles mixed gz/non-gz automatically
    input_str = join(fasta_files, ",")

    run(`$(CONDA_RUNNER) run -n checkm2 checkm2 predict -i $(input_str) -o $(outdir) --database_path $(db_path) --threads $(threads)`)

    return outdir
end

# FASTQ Quality Analysis
#

"""
Stores quality distribution statistics for FASTQ analysis.
"""
struct QualityDistribution
    q20_percent::Float64
    q30_percent::Float64
    q40_percent::Float64
end

"""
Stores comprehensive FASTQ quality analysis results.
"""
struct FastqQualityResults
    n_reads::Int
    mean_quality::Float64
    mean_length::Float64
    gc_content::Float64
    quality_distribution::QualityDistribution
end

"""
    analyze_fastq_quality(fastq_file::String)

Analyzes quality metrics for a FASTQ file.

Calculates comprehensive quality statistics including read count, quality scores,
length distribution, GC content, and quality threshold percentages.

# Arguments
- `fastq_file`: Path to FASTQ file (can be gzipped)

# Returns
FastqQualityResults with the following fields:
- `n_reads`: Total number of reads
- `mean_quality`: Average Phred quality score across all reads
- `mean_length`: Average read length
- `gc_content`: GC content percentage
- `quality_distribution`: QualityDistribution with Q20+, Q30+, Q40+ percentages

# Example
```julia
quality_stats = Mycelia.analyze_fastq_quality("reads.fastq")
println("Total reads: \$(quality_stats.n_reads)")
println("Mean quality: \$(quality_stats.mean_quality)")
println("Q30+ reads: \$(quality_stats.quality_distribution.q30_percent)%")
```
"""
function analyze_fastq_quality(fastq_file::String)
    if !isfile(fastq_file)
        error("FASTQ file does not exist: $(fastq_file)")
    end

    n_reads = 0
    total_quality = 0.0
    total_length = 0
    gc_count = 0
    total_bases = 0
    q20_reads = 0
    q30_reads = 0
    q40_reads = 0

    reader = FASTX.FASTQ.Reader(open(fastq_file, "r"))

    try
        for record in reader
            n_reads += 1

            # Get sequence and quality
            sequence = FASTX.sequence(record)
            quality_scores = FASTX.quality_scores(record)

            # Calculate length
            read_length = length(sequence)
            total_length += read_length
            total_bases += read_length

            # Calculate quality metrics
            mean_read_quality = Statistics.mean(quality_scores)
            total_quality += mean_read_quality

            # Count quality thresholds
            if mean_read_quality >= 20
                q20_reads += 1
            end
            if mean_read_quality >= 30
                q30_reads += 1
            end
            if mean_read_quality >= 40
                q40_reads += 1
            end

            # Count GC content
            for base in sequence
                if base == BioSequences.DNA_G || base == BioSequences.DNA_C
                    gc_count += 1
                end
            end
        end
    finally
        close(reader)
    end

    if n_reads == 0
        error("No reads found in FASTQ file: $(fastq_file)")
    end

    # Calculate final metrics
    mean_quality = total_quality / n_reads
    mean_length = total_length / n_reads
    gc_content = (gc_count / total_bases) * 100.0

    quality_dist = QualityDistribution(
        (q20_reads / n_reads) * 100.0,
        (q30_reads / n_reads) * 100.0,
        (q40_reads / n_reads) * 100.0
    )

    return FastqQualityResults(
        n_reads,
        mean_quality,
        mean_length,
        gc_content,
        quality_dist
    )
end

"""
    run_quast(assembly_files::Vector{String}; outdir::Union{String,Nothing}=nothing, reference::Union{String,Nothing}=nothing, threads::Int=get_default_threads(), min_contig::Int=500, gene_finding::Bool=false)

Run QUAST (Quality Assessment Tool for Genome Assemblies) to evaluate assembly quality.

# Arguments
- `assembly_files::Vector{String}`: Vector of paths to assembly FASTA files to evaluate
- `outdir::Union{String,Nothing}=nothing`: Output directory for QUAST results. If `nothing`, defaults to `dirname(first_assembly)/basename_without_fasta_quast`
- `reference::Union{String,Nothing}=nothing`: Optional reference genome for reference-based metrics
- `threads::Int=get_default_threads()`: Number of threads to use
- `min_contig::Int=500`: Minimum contig length to consider
- `gene_finding::Bool=false`: Whether to run gene finding (requires GeneMark-ES/ET)

# Returns
- `String`: Path to the output directory containing QUAST results

# Output Files
- `report.html`: Interactive HTML report
- `report.txt`: Text summary report
- `report.tsv`: Tab-separated values report for programmatic access
- `transposed_report.tsv`: Transposed TSV format
- `icarus.html`: Icarus contig browser (if reference provided)

# Examples
```julia
# Basic assembly evaluation
assemblies = ["assembly1.fasta", "assembly2.fasta"]
quast_dir = Mycelia.run_quast(assemblies)

# With reference genome
ref_genome = "reference.fasta"
quast_dir = Mycelia.run_quast(assemblies, reference=ref_genome)

# Custom parameters
quast_dir = Mycelia.run_quast(assemblies, 
                             outdir="my_quast_results",
                             min_contig=1000,
                             threads=8)
```

# Notes
- Requires QUAST to be installed via Bioconda
- Without reference: provides basic metrics (N50, total length, # contigs, etc.)
- With reference: adds reference-based metrics (genome fraction, misassemblies, etc.)
- Gene finding requires additional dependencies and is disabled by default
"""
function run_quast(assembly_files::Vector{String};
        outdir::Union{String, Nothing} = nothing,
        reference::Union{String, Nothing} = nothing,
        threads::Int = get_default_threads(),
        min_contig::Int = 500,
        gene_finding::Bool = false)

    # Install QUAST via Bioconda if needed
    add_bioconda_env("quast")

    # Validate input files
    for assembly_file in assembly_files
        if !isfile(assembly_file)
            error("Assembly file does not exist: $(assembly_file)")
        end
    end

    # Validate reference if provided
    if reference !== nothing && !isfile(reference)
        error("Reference file does not exist: $(reference)")
    end

    # Derive default outdir from first assembly if not provided
    if outdir === nothing
        first_asm = assembly_files[1]
        base = replace(basename(first_asm), Mycelia.FASTA_REGEX => "")
        outdir = joinpath(dirname(first_asm), base * "_quast")
    end

    # Derive default outdir from first assembly if not provided
    if outdir === nothing
        first_asm = assembly_files[1]
        base = replace(basename(first_asm), Mycelia.FASTA_REGEX => "")
        outdir = joinpath(dirname(first_asm), base * "_busco")
    end

    # Create output directory
    mkpath(outdir)

    # Build QUAST command
    cmd_parts = [
        "quast.py",
        "--output-dir", outdir,
        "--threads", string(threads),
        "--min-contig", string(min_contig)
    ]

    # Add reference if provided
    if reference !== nothing
        push!(cmd_parts, "--reference", reference)
    end

    # Add gene finding if requested
    if gene_finding
        push!(cmd_parts, "--gene-finding")
    end

    # Add assembly files
    append!(cmd_parts, assembly_files)

    # Run QUAST
    try
        println("Running QUAST on $(length(assembly_files)) assemblies...")
        if reference !== nothing
            println("Using reference: $(reference)")
        end

        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n quast $cmd_parts`)

        # Verify output files were created
        report_txt = joinpath(outdir, "report.txt")
        report_tsv = joinpath(outdir, "report.tsv")

        if !isfile(report_txt) || !isfile(report_tsv)
            error("QUAST failed to generate expected output files")
        end

        println("QUAST analysis completed successfully")
        println("Results saved in: $(outdir)")
        println("View HTML report: $(joinpath(outdir, "report.html"))")

        return outdir

    catch e
        error("QUAST execution failed: $(e)")
    end
end

"""
    run_quast(assembly_file::String; kwargs...)

Run QUAST on a single assembly file. See `run_quast(::Vector{String})` for details.
"""
function run_quast(assembly_file::String; kwargs...)
    return run_quast([assembly_file]; kwargs...)
end

"""
    run_busco(assembly_files::Vector{String}; outdir::Union{String,Nothing}=nothing, lineage::Union{String,Nothing}=nothing, mode::String="genome", threads::Int=get_default_threads(), force::Bool=false, auto_lineage::Bool=true, auto_lineage_euk::Bool=false, auto_lineage_prok::Bool=false)

Run BUSCO (Benchmarking Universal Single-Copy Orthologs) to assess genome assembly completeness.

# Arguments
- `assembly_files::Vector{String}`: Vector of paths to assembly FASTA files to evaluate
- `outdir::Union{String,Nothing}=nothing`: Output directory for BUSCO results. If `nothing`, defaults to `dirname(first_assembly)/basename_without_fasta_busco`
- `lineage::Union{String,Nothing}=nothing`: BUSCO lineage dataset to use (e.g., "bacteria_odb10", "eukaryota_odb10"); if `nothing`, auto-lineage is used
- `mode::String="genome"`: BUSCO mode ("genome", "transcriptome", "proteins")
- `threads::Int=get_default_threads()`: Number of threads to use (`--cpu`)
- `force::Bool=false`: Force overwrite existing results
- `auto_lineage::Bool=true`: Use BUSCO auto-lineage placement (default when lineage is not provided)
- `auto_lineage_euk::Bool=false`: Auto-lineage limited to eukaryotic tree
- `auto_lineage_prok::Bool=false`: Auto-lineage limited to prokaryotic tree

# Returns
- `String`: Path to the output directory containing BUSCO results

# Output Files
- `short_summary.specific.lineage.txt`: Summary statistics
- `full_table.tsv`: Complete BUSCO results table
- `missing_busco_list.tsv`: List of missing BUSCOs
- `run_lineage/`: Detailed results directory

# Examples
```julia
# Basic completeness assessment
assemblies = ["assembly1.fasta", "assembly2.fasta"]
busco_dir = Mycelia.run_busco(assemblies)

# Specific lineage
busco_dir = Mycelia.run_busco(assemblies, lineage="bacteria_odb10")

# Custom parameters
busco_dir = Mycelia.run_busco(assemblies,
                             outdir="my_busco_results",
                             lineage="enterobacterales_odb10",
                             threads=8)
```

# Notes
- Requires BUSCO to be installed via Bioconda
- Auto lineage detection requires internet connection for first run
- Available lineages: bacteria_odb10, archaea_odb10, eukaryota_odb10, fungi_odb10, etc.
- Results provide Complete, Fragmented, and Missing BUSCO counts
"""
function run_busco(assembly_files::Vector{String};
        outdir::Union{String, Nothing} = nothing,
        lineage::Union{String, Nothing} = nothing,
        mode::String = "genome",
        threads::Int = get_default_threads(),
        force::Bool = false,
        auto_lineage::Bool = true,
        auto_lineage_euk::Bool = false,
        auto_lineage_prok::Bool = false)

    # Install BUSCO via Bioconda if needed
    add_bioconda_env("busco")

    # Validate input files
    for assembly_file in assembly_files
        if !isfile(assembly_file)
            error("Assembly file does not exist: $(assembly_file)")
        end
    end

    # Validate mode
    valid_modes = ["genome", "transcriptome", "proteins"]
    if !(mode in valid_modes)
        error("Invalid mode: $(mode). Must be one of: $(join(valid_modes, ", "))")
    end

    # Validate lineage / auto-lineage flags
    if auto_lineage + auto_lineage_euk + auto_lineage_prok > 1
        error("Specify only one of auto_lineage, auto_lineage_euk, or auto_lineage_prok.")
    end
    if lineage !== nothing && (auto_lineage || auto_lineage_euk || auto_lineage_prok)
        error("Provide either lineage or an auto-lineage flag, not both.")
    end
    # Default to auto-lineage if nothing provided
    if lineage === nothing && !(auto_lineage || auto_lineage_euk || auto_lineage_prok)
        auto_lineage = true
    end

    # Derive default outdir from first assembly if not provided
    if outdir === nothing
        first_asm = assembly_files[1]
        base = replace(basename(first_asm), Mycelia.FASTA_REGEX => "")
        outdir = joinpath(dirname(first_asm), base * "_busco")
    end

    # Create output directory
    mkpath(outdir)

    # Process each assembly file
    results_dirs = String[]

    for (i, assembly_file) in enumerate(assembly_files)
        assembly_name = splitext(basename(assembly_file))[1]
        assembly_outdir = joinpath(outdir, assembly_name)

        # Build BUSCO command
        cmd_parts = [
            "busco",
            "--in", assembly_file,
            "--out", assembly_name,
            "--out_path", outdir,
            "--mode", mode,
            "--cpu", string(threads)
        ]

        # Lineage / auto-lineage selection
        if lineage !== nothing
            push!(cmd_parts, "--lineage_dataset", lineage)
        elseif auto_lineage_euk
            push!(cmd_parts, "--auto-lineage-euk")
        elseif auto_lineage_prok
            push!(cmd_parts, "--auto-lineage-prok")
        else
            push!(cmd_parts, "--auto-lineage")
        end

        # Add force flag if requested
        if force
            push!(cmd_parts, "--force")
        end

        # Add quiet flag to reduce output
        push!(cmd_parts, "--quiet")

        try
            println("Running BUSCO on $(assembly_file)...")
            println("Lineage: $(lineage), Mode: $(mode)")

            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n busco $cmd_parts`)

            # Verify output files were created
            summary_files = String[]
            if isdir(assembly_outdir)
                summary_files = filter(
                    name -> occursin(r"^short_summary.*\.txt$", name), readdir(assembly_outdir))
                if isempty(summary_files)
                    for (root, _, files) in walkdir(assembly_outdir)
                        for name in files
                            if occursin(r"^short_summary.*\.txt$", name)
                                push!(summary_files, joinpath(root, name))
                            end
                        end
                    end
                end
            end

            if isempty(summary_files)
                error("BUSCO failed to generate summary file for $(assembly_file)")
            end

            println("BUSCO analysis completed for $(assembly_file)")
            println("Results saved in: $(assembly_outdir)")

            push!(results_dirs, assembly_outdir)

        catch e
            error("BUSCO execution failed for $(assembly_file): $(e)")
        end
    end

    # Print summary information
    println("\nBUSCO analysis completed for all assemblies")
    println("Results directories:")
    for dir in results_dirs
        println("  - $(dir)")
    end

    return outdir
end

"""
    run_busco(assembly_file::String; kwargs...)

Run BUSCO on a single assembly file. See `run_busco(::Vector{String})` for details.
"""
function run_busco(assembly_file::String; kwargs...)
    return run_busco([assembly_file]; kwargs...)
end

"""
    run_mummer(reference::String, query::String; outdir::String="mummer_results", prefix::String="out", mincluster::Int=65, minmatch::Int=20, threads::Int=1)

Run MUMmer for genome comparison and alignment between reference and query sequences.

# Arguments
- `reference::String`: Path to reference genome FASTA file
- `query::String`: Path to query genome FASTA file  
- `outdir::String="mummer_results"`: Output directory for MUMmer results
- `prefix::String="out"`: Prefix for output files
- `mincluster::Int=65`: Minimum cluster length for nucmer
- `minmatch::Int=20`: Minimum match length for nucmer
- `threads::Int=1`: Number of threads (note: MUMmer has limited multithreading)

# Returns
- `String`: Path to the output directory containing MUMmer results

# Output Files
- `prefix.delta`: Delta alignment file (main output)
- `prefix.coords`: Human-readable coordinates file
- `prefix.snps`: SNP/indel report (if show-snps is run)
- `prefix.plot.png`: Dot plot visualization (if mummerplot is run)

# Examples
```julia
# Basic genome comparison
ref_genome = "reference.fasta"
query_genome = "assembly.fasta"
mummer_dir = Mycelia.run_mummer(ref_genome, query_genome)

# Custom parameters
mummer_dir = Mycelia.run_mummer(ref_genome, query_genome,
                               outdir="comparison_results",
                               prefix="comparison",
                               mincluster=100,
                               minmatch=30)
```

# Notes
- Requires MUMmer to be installed via Bioconda
- nucmer is used for DNA sequence alignment
- show-coords generates human-readable coordinate output
- Results include alignment coordinates, percent identity, and coverage
- For visualization, use mummerplot (requires gnuplot)
"""
function run_mummer(reference::String, query::String;
        outdir::String = "mummer_results",
        prefix::String = "out",
        mincluster::Int = 65,
        minmatch::Int = 20,
        threads::Int = 1)

    # Install MUMmer via Bioconda if needed
    add_bioconda_env("mummer")

    # Validate input files
    if !isfile(reference)
        error("Reference file does not exist: $(reference)")
    end

    if !isfile(query)
        error("Query file does not exist: $(query)")
    end

    # Create output directory
    mkpath(outdir)

    # Define output files
    delta_file = joinpath(outdir, "$(prefix).delta")
    coords_file = joinpath(outdir, "$(prefix).coords")

    try
        println("Running MUMmer comparison...")
        println("Reference: $(reference)")
        println("Query: $(query)")

        # Run nucmer for alignment
        nucmer_cmd = [
            "nucmer",
            "--prefix", joinpath(outdir, prefix),
            "--mincluster", string(mincluster),
            "--minmatch", string(minmatch),
            reference,
            query
        ]

        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mummer $nucmer_cmd`)

        # Generate coordinate output
        coords_cmd = [
            "show-coords",
            "-r", "-c", "-l",  # -r: sort by reference, -c: show percent coverage, -l: show length
            delta_file
        ]

        coords_output = read(`$(Mycelia.CONDA_RUNNER) run -n mummer $coords_cmd`, String)

        # Write coordinates to file
        open(coords_file, "w") do io
            write(io, coords_output)
        end

        # Verify output files were created
        if !isfile(delta_file)
            error("MUMmer failed to generate delta file")
        end

        if !isfile(coords_file)
            error("Failed to generate coordinates file")
        end

        println("MUMmer analysis completed successfully")
        println("Results saved in: $(outdir)")
        println("Delta file: $(delta_file)")
        println("Coordinates file: $(coords_file)")

        return outdir

    catch e
        error("MUMmer execution failed: $(e)")
    end
end

"""
    parse_mummer_coords(coords_file) -> Vector{NamedTuple}

Parse a MUMmer `show-coords` output file into alignment records.
"""
function parse_mummer_coords(coords_file)
    alignments = NamedTuple[]
    if !isfile(coords_file)
        return alignments
    end

    open(coords_file) do io
        for raw_line in eachline(io)
            line = strip(raw_line)
            isempty(line) && continue
            fields = split(line)
            first_val = tryparse(Int, fields[1])
            if first_val === nothing || length(fields) < 7
                continue
            end
            start_ref = parse(Int, fields[1])
            end_ref = parse(Int, fields[2])
            start_query = parse(Int, fields[3])
            end_query = parse(Int, fields[4])
            len_ref = parse(Int, fields[5])
            identity = parse(Float64, fields[7])
            push!(alignments, (;
                start_ref, end_ref, start_query, end_query, len_ref, identity))
        end
    end
    return alignments
end

"""
    run_mummer_plot(delta_file::String; outdir::String="", prefix::String="plot", plot_type::String="png")

Generate dot plot visualization from MUMmer delta file using mummerplot.

# Arguments
- `delta_file::String`: Path to MUMmer delta file
- `outdir::String=""`: Output directory (defaults to same as delta file)
- `prefix::String="plot"`: Prefix for plot files
- `plot_type::String="png"`: Plot format ("png", "ps", "x11")

# Returns
- `String`: Path to the generated plot file

# Notes
- Requires gnuplot to be installed
- Useful for visualizing genome alignments and rearrangements
"""
function run_mummer_plot(delta_file::String;
        outdir::String = "",
        prefix::String = "plot",
        plot_type::String = "png")

    # Install MUMmer via Bioconda if needed
    add_bioconda_env("mummer")

    if !isfile(delta_file)
        error("Delta file does not exist: $(delta_file)")
    end

    # Set output directory
    if outdir == ""
        outdir = dirname(delta_file)
    else
        mkpath(outdir)
    end

    plot_file = joinpath(outdir, "$(prefix).$(plot_type)")

    try
        # Run mummerplot
        plot_cmd = [
            "mummerplot",
            "--$(plot_type)",
            "--prefix", joinpath(outdir, prefix),
            delta_file
        ]

        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mummer $plot_cmd`)

        if !isfile(plot_file)
            error("Failed to generate plot file")
        end

        println("MUMmer plot generated: $(plot_file)")
        return plot_file

    catch e
        error("MUMmer plot generation failed: $(e)")
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assess assembly quality using k-mer based metrics including QV scores.

This function provides a comprehensive quality assessment of genome assemblies by
comparing k-mer distributions between the assembly and source sequencing data.
It calculates multiple quality metrics across different k-mer sizes, with the 
QV (Quality Value) score being the primary measure of assembly accuracy.

# Arguments
- `assembly`: Path to assembly file (FASTA format) or assembly sequences
- `observations`: Path to sequencing data file(s) (FASTQ format) or sequence observations
- `ks::Vector{Int}`: K-mer sizes to analyze (default: 11,13,17,19,23,31,53 for comprehensive evaluation)

# Returns
DataFrame containing quality metrics for each k-mer size:
- `k::Int`: K-mer length used
- `cosine_distance::Float64`: Cosine similarity between k-mer distributions  
- `js_divergence::Float64`: Jensen-Shannon divergence between distributions
- `qv::Float64`: Merqury-style Quality Value score (primary assembly accuracy metric)

# QV Score Interpretation
- **QV ≥ 40**: High quality assembly (>99.99% base accuracy)
- **QV 30-40**: Good quality assembly (99.9-99.99% base accuracy) 
- **QV 20-30**: Moderate quality assembly (99-99.9% base accuracy)
- **QV < 20**: Lower quality assembly (<99% base accuracy)

# Implementation Details
The QV score follows the Merqury methodology:
1. Count k-mers in both assembly and sequencing data
2. Calculate shared k-mers between datasets
3. Estimate base-level accuracy: P = (Kshared/Ktotal)^(1/k)
4. Convert to Phred scale: QV = -10log₁₀(1-P)

# Performance Notes
- Uses multiple k-mer sizes for robust quality assessment
- Larger k-mer sizes provide more discriminative power
- Parallel k-mer counting for computational efficiency
- Progress tracking for long-running analyses

# References
Rhie, A. et al. "Merqury: reference-free quality, completeness, and phasing 
assessment for genome assemblies." Genome Biology 21, 245 (2020).
"""
function assess_assembly_quality(; assembly, observations,
        ks::Vector{Int} = filter(x -> x in [11, 13, 17, 19, 23, 31, 53], Mycelia.ks()))
    return assess_assembly_kmer_quality(assembly = assembly, observations = observations, ks = ks)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate QV score heatmap visualization for comparing assembler performance.

Creates a heatmap showing QV scores across different k-mer sizes and assemblers,
enabling easy comparison of assembly quality across multiple tools and parameters.

# Arguments
- `qv_results`: DataFrame containing QV results (from assess_assembly_quality)
- `assembler_column::String`: Column name identifying different assemblers (default: "assembler")
- `k_column::String`: Column name with k-mer sizes (default: "k") 
- `qv_column::String`: Column name with QV scores (default: "qv")
- `title::String`: Plot title (default: "Assembly Quality (QV Scores)")
- `clims::Tuple`: Color scale limits (default: (0, 60))

# Returns
Plots.jl heatmap object showing QV scores by assembler and k-mer size

# Visualization Features
- Color-coded heatmap with QV scores as values
- K-mer sizes on y-axis, assemblers on x-axis
- Color scale optimized for typical QV ranges
- Missing data handled gracefully
- Customizable styling and limits

# Usage Examples
```julia
# Basic heatmap
results = assess_assembly_quality(assembly="contigs.fa", observations=["reads.fq"])
plot = generate_qv_heatmap(results)

# Custom styling
plot = generate_qv_heatmap(results, 
                          title="HiFi Assembly Comparison",
                          clims=(20, 50))
```
"""
function generate_qv_heatmap(qv_results::DataFrames.DataFrame;
        assembler_column::String = "assembler",
        k_column::String = "k",
        qv_column::String = "qv",
        title::String = "Assembly Quality (QV Scores)",
        clims::Tuple = (0, 60))

    # Get unique assemblers and k-values
    assemblers = sort(unique(qv_results[!, assembler_column]))
    ks = sort(unique(qv_results[!, k_column]))

    # Create QV matrix
    qv_matrix = zeros(length(ks), length(assemblers))

    for (i, k) in enumerate(ks)
        for (j, assembler) in enumerate(assemblers)
            matching_rows = (qv_results[!, k_column] .== k) .&
                            (qv_results[!, assembler_column] .== assembler)
            if any(matching_rows)
                qv_matrix[i, j] = first(qv_results[matching_rows, qv_column])
            else
                qv_matrix[i, j] = NaN  # Missing data
            end
        end
    end

    # Create heatmap
    heatmap_plot = Plots.heatmap(
        qv_matrix,
        xlabel = "Assembler",
        ylabel = "K-mer Size",
        title = title,
        xticks = (1:length(assemblers), assemblers),
        yticks = (1:length(ks), ks),
        colorbar_title = "QV Score",
        clims = clims,
        size = (600, 400),
        margin = 5Plots.PlotMeasures.mm
    )

    return heatmap_plot
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run ALE (Assembly Likelihood Evaluation) for reference-free quality assessment.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `reads_file::String`: Path to reads FASTQ file (can be comma-separated for paired reads)
- `outdir::String`: Output directory path (default: "\${assembly_file}_ale")
- `threads::Int`: Number of threads to use for read mapping (default: `get_default_threads()`)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `ale_score::String`: Path to ALE score file
- `ale_plot::String`: Path to ALE plot file

# Details
- Uses ALE's per-base likelihood scoring for reference-free quality assessment
- Evaluates assembly quality based on read mapping consistency
- Identifies potential misassemblies and low-quality regions
- Automatically creates and uses a conda environment with ale
- Utilizes requested CPU threads for BWA mapping
- Skips analysis if output files already exist
"""
function run_ale(assembly_file::String, reads_file::String;
        outdir::String = assembly_file * "_ale", threads::Int = get_default_threads())
    Mycelia.add_bioconda_env("ale")
    mkpath(outdir)

    basename_assembly = basename(assembly_file, ".fasta")
    ale_score_file = joinpath(outdir, basename_assembly * ".ale")
    ale_plot_file = joinpath(outdir, basename_assembly * "_ALE_plot.png")

    if !isfile(ale_score_file)
        # Map reads to assembly first
        sam_file = joinpath(outdir, basename_assembly * ".sam")
        if !isfile(sam_file)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ale bwa index $(assembly_file)`)
            if occursin(",", reads_file)
                # Paired-end reads
                reads = split(reads_file, ",")
                run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ale bwa mem -t $(threads) $(assembly_file) $(reads[1]) $(reads[2]) -o $(sam_file)`)
            else
                # Single-end reads
                run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ale bwa mem -t $(threads) $(assembly_file) $(reads_file) -o $(sam_file)`)
            end
        end

        # Run ALE evaluation
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ale ALE $(sam_file) $(assembly_file) $(ale_score_file)`)

        # Generate plot if available
        if isfile(ale_score_file)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ale ALEplot $(ale_score_file) $(assembly_file) $(ale_plot_file)`)
        end
    end

    return (; outdir, ale_score = ale_score_file, ale_plot = ale_plot_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run FRCbam for reference-free misassembly detection using Feature Response Curves.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `reads_file::String`: Path to reads FASTQ file (can be comma-separated for paired reads)
- `outdir::String`: Output directory path (default: "\${assembly_file}_frcbam")
- `threads::Int`: Number of threads to use for read mapping (default: `get_default_threads()`)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `frc_txt::String`: Path to FRC results text file
- `frc_plot::String`: Path to FRC plot file

# Details
- Uses Feature Response Curves for reference-free misassembly detection
- Evaluates paired-end consistency and insert size distributions
- Identifies structural variations and assembly errors
- Automatically creates and uses a conda environment with frcbam
- Utilizes requested CPU threads for BWA mapping
- Skips analysis if output files already exist
"""
function run_frcbam(assembly_file::String, reads_file::String;
        outdir::String = assembly_file * "_frcbam", threads::Int = get_default_threads())
    Mycelia.add_bioconda_env("frcbam")
    mkpath(outdir)

    basename_assembly = basename(assembly_file, ".fasta")
    frc_txt_file = joinpath(outdir, basename_assembly * "_FRC.txt")
    frc_plot_file = joinpath(outdir, basename_assembly * "_FRC.png")

    if !isfile(frc_txt_file)
        # Map reads to assembly first
        bam_file = joinpath(outdir, basename_assembly * ".bam")
        if !isfile(bam_file)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n frcbam bwa index $(assembly_file)`)
            if occursin(",", reads_file)
                # Paired-end reads
                reads = split(reads_file, ",")
                run(
                    pipeline(
                    `$(Mycelia.CONDA_RUNNER) run --live-stream -n frcbam bwa mem -t $(threads) $(assembly_file) $(reads[1]) $(reads[2])`,
                    `samtools sort -@ $(get_default_threads()) -o $(bam_file)`
                )
                )
            else
                # Single-end reads (FRCbam works best with paired-end)
                run(
                    pipeline(
                    `$(Mycelia.CONDA_RUNNER) run --live-stream -n frcbam bwa mem -t $(threads) $(assembly_file) $(reads_file)`,
                    `samtools sort -@ $(get_default_threads()) -o $(bam_file)`
                )
                )
            end
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n frcbam samtools index $(bam_file)`)
        end

        # Run FRCbam analysis
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n frcbam FRC --pe-bam $(bam_file) --output $(frc_txt_file)`)

        # Generate plot if available
        if isfile(frc_txt_file)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n frcbam FRCplot --FRC-file $(frc_txt_file) --output $(frc_plot_file)`)
        end
    end

    return (; outdir, frc_txt = frc_txt_file, frc_plot = frc_plot_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run 4CAC for contig classification (virus/plasmid/prokaryote/eukaryote).

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `outdir::String`: Output directory path (default: "\${assembly_file}_4cac")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `predictions::String`: Path to classification predictions file

# Details
- Uses machine learning to classify contigs into virus/plasmid/prokaryote/eukaryote
- Based on tetranucleotide frequency and other sequence features
- Useful for binning and filtering metagenomic assemblies
- Automatically creates and uses a conda environment with 4cac
- Skips analysis if output files already exist
"""
function run_4cac(assembly_file::String; outdir::String = assembly_file * "_4cac")
    Mycelia.add_bioconda_env("4cac")
    mkpath(outdir)

    basename_assembly = basename(assembly_file, ".fasta")
    predictions_file = joinpath(outdir, basename_assembly * "_predictions.txt")

    if !isfile(predictions_file)
        # Run 4CAC classification
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n 4cac 4CAC -i $(assembly_file) -o $(predictions_file)`)
    end

    return (; outdir, predictions = predictions_file)
end

#
