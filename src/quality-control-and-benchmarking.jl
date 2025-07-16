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

"""
    assess_assembly_quality(contigs_file)

Assess basic assembly quality metrics from a FASTA file.

Calculates standard assembly quality metrics including contig count, total length,
and N50 statistic for assembly evaluation.

# Arguments
- `contigs_file`: Path to FASTA file containing assembly contigs

# Returns
- Tuple of (n_contigs, total_length, n50)
  - `n_contigs`: Number of contigs in the assembly
  - `total_length`: Total length of all contigs in base pairs
  - `n50`: N50 statistic (length of shortest contig in the set covering 50% of assembly)

# Example
```julia
n_contigs, total_length, n50 = assess_assembly_quality("assembly.fasta")
println("Assembly has \$n_contigs contigs, \$total_length bp total, N50=\$n50")
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
    
    # Calculate N50
    sorted_lengths = sort(contigs, rev=true)
    cumsum_lengths = cumsum(sorted_lengths)
    target_length = total_length / 2
    n50_idx = findfirst(x -> x >= target_length, cumsum_lengths)
    n50 = n50_idx !== nothing ? sorted_lengths[n50_idx] : 0
    
    return n_contigs, total_length, n50
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
function setup_checkv(; db_dir::String=joinpath(homedir(), "workspace", ".checkv"))
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
function setup_checkm(; db_dir::String=joinpath(homedir(), "workspace", ".checkm"))
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
function setup_checkm2(; db_dir::String=joinpath(homedir(), "workspace", ".checkm2"))
    add_bioconda_env("checkm2")
    if !isdir(db_dir)
        mkpath(db_dir)
        run(`$(CONDA_RUNNER) run -n checkm2 checkm2 database --download --path $(db_dir)`)
    end
    return db_dir
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
function run_checkv(fasta_file::String; outdir::String=fasta_file * "_checkv", db_dir::String=joinpath(homedir(), "workspace", ".checkv"), threads::Int=Sys.CPU_THREADS)
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
    if isdir(outdir) && all(isfile, output_files)
        @warn "CheckV output directory already exists with all expected files. Skipping analysis: $(outdir)"
    else
        latest_db = last(readdir(Mycelia.setup_checkv(), join=true))
        run(`$(CONDA_RUNNER) run -n checkv checkv end_to_end $(fasta_file) $(outdir) -d $(latest_db) -t $(threads)`)
    end
    
    # Clean up temporary directory if it exists
    tmp_dir = joinpath(outdir, "tmp")
    Mycelia.cleanup_directory(tmp_dir, verbose=true, force=true)
    
    return (outdir=outdir, output_files...)
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
function run_checkm(input_path::String; outdir::String=input_path * "_checkm", db_dir::String=joinpath(homedir(), "workspace", ".checkm"), extension::String="fasta", threads::Int=Sys.CPU_THREADS)
    setup_checkm(db_dir=db_dir)
    
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

# Example
```julia
run_checkm2("genome.fasta")
run_checkm2("./genomes/")
```
"""
function run_checkm2(input_path::String; outdir::String=input_path * "_checkm2", db_dir::String=joinpath(homedir(), "workspace", ".checkm2"), threads::Int=Sys.CPU_THREADS)
    setup_checkm2(db_dir=db_dir)
    
    fasta_files = find_fasta_files(input_path)
    
    if isempty(fasta_files)
        error("No FASTA files found in: $(input_path)")
    end
    
    # CheckM2 can take comma-separated fasta files as input
    input_str = join(fasta_files, ",")
    
    run(`$(CONDA_RUNNER) run -n checkm2 checkm2 predict -i $(input_str) -o $(outdir) --database_path $(db_dir) --threads $(threads)`)
    
    return outdir
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
function run_checkm2_list(fasta_files::Vector{String}; outdir::String=normalized_current_datetime() * "_checkm2", db_dir::String=joinpath(homedir(), "workspace", ".checkm2"), threads::Int=Sys.CPU_THREADS)
    setup_checkm2(db_dir=db_dir)
    
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
    
    run(`$(CONDA_RUNNER) run -n checkm2 checkm2 predict -i $(input_str) -o $(outdir) --database_path $(db_dir) --threads $(threads)`)
    
    return outdir
end

#