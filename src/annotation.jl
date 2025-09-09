function run_phage_annotation(;
    fasta::AbstractString,
    pharroka_dbdir = joinpath(homedir(), "workspace", "pharroka"),
    phold_dbdir = joinpath(homedir(), "workspace", "phold"),
    phynteny_dbdir = joinpath(homedir(), "workspace", "phynteny"),
    prefix = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
    pharroka_outdir = replace(fasta, Mycelia.FASTA_REGEX => "_pharroka"),
    phold_outdir = joinpath(pharroka_outdir, "phold"),
    phynteny_outdir = joinpath(phold_outdir, "phynteny"),
    threads = Sys.CPU_THREADS
)
    Mycelia.add_bioconda_env("pharokka")
    Mycelia.add_bioconda_env("phold")
    Mycelia.add_bioconda_env("phynteny_transformer")
    if !isdir(pharroka_dbdir) || isempty(readdir(pharroka_dbdir))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pharroka install_databases.py -o $(pharroka_dbdir)`)
    end
    if !isdir(pharroka_outdir) || isempty(readdir(pharroka_outdir))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pharroka  pharokka.py -i $(fasta) -o $(pharroka_outdir) -d $(pharroka_dbdir) -t $(threads)  -m  -g prodigal-gv --force`)
    end
    if !isdir(phold_dbdir) || isempty(readdir(phold_dbdir))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phold phold install --database $(phold_dbdir)`)
    end
    pharroko_gbk = joinpath(pharroka_outdir, prefix * ".gbk")
    @assert isfile(pharroko_gbk)
    if !isdir(phold_outdir) || isempty(readdir(phold_outdir))
        # skip GPU integration at the moment
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phold phold run --input $(pharroko_gbk) --prefix $(prefix) --output $(phold_outdir) --database $(phold_dbdir) --threads $(threads) --cpu --force`)
    end
    if !isdir(phynteny_dbdir) || isempty(readddir(phynteny_dbdir))
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phynteny_transformer install_models -o $(phynteny_dbdir)`)
    end
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n phynteny_transformer  phynteny_transformer --help --advanced`)
end

"""
    run_pgap(; fasta, organism, pgap_dir, outdir, as_string, force, threads, prefix)

Run the PGAP (Prokaryotic Genome Annotation Pipeline) tool on a FASTA file.

This function automatically handles compressed FASTA files (gzip, bzip2, xz) by creating 
temporary uncompressed versions when needed. The temporary files are automatically cleaned 
up after processing, regardless of whether the PGAP run succeeds or fails.

# Arguments
- `fasta::AbstractString`: Path to the input FASTA file. Can be compressed (.gz, .bz2, .xz) or uncompressed.
- `organism::AbstractString`: Organism name for taxonomic classification and annotation.
- `pgap_dir`: Directory containing PGAP installation. Defaults to `~/workspace/pgap`.
- `outdir`: Output directory for results. Defaults to input filename with "_pgap" suffix.
- `as_string::Bool=false`: If true, returns the command string instead of executing it.
- `force::Bool=false`: If true, overwrites existing output directory.
- `threads::Int`: Number of CPU threads to use. Defaults to system CPU count.
- `prefix`: Prefix for output files. Defaults to input filename without extension.

# Returns
- `String`: Path to output directory when `as_string=false`, or bash command string when `as_string=true`.

# Notes
The function will automatically download and set up PGAP if not already present in the specified directory.
For compressed input files, temporary uncompressed versions are created and automatically cleaned up.
The PGAP tool requires uncompressed FASTA files, so this function handles the decompression transparently.
"""
function run_pgap(;
        fasta::AbstractString,
        organism::AbstractString,
        pgap_dir = joinpath(homedir(), "workspace", "pgap"),
        outdir = replace(fasta, Mycelia.FASTA_REGEX => "_pgap"),
        as_string = false,
        force = false,
        threads = Sys.CPU_THREADS,
        prefix = replace(basename(fasta), Mycelia.FASTA_REGEX => "")
        # --memory
    )
    pgap_py_script = joinpath(pgap_dir, "pgap.py")
    if !isfile(pgap_py_script)
        mkpath(dirname(pgap_py_script))
        run(`wget -O $(pgap_py_script) https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py`)
        # run(`chmod +x $(pgap_py_script)`)
        @assert isfile(pgap_py_script)
    end
    if isempty(filter(x -> occursin(r"input\-", x), readdir(pgap_dir)))
        @time run(setenv(`python3 $(pgap_py_script) --update`, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
    end 
    if force && isdir(outdir)
        rm(outdir, recursive=true)
    elseif !force && isdir(outdir)
        @warn "$outdir already exists, use `force=true` to overwrite existing results"
        return outdir
    end    
    @assert !isdir(outdir)
    
    # Detect if the input FASTA is compressed and handle accordingly
    is_compressed = endswith(lowercase(fasta), ".gz") || endswith(lowercase(fasta), ".bz2") || endswith(lowercase(fasta), ".xz")
    
    if is_compressed
        # Create temporary uncompressed file
        temp_fasta = tempname() * ".fasta"
        
        try
            # Decompress based on file extension
            if endswith(lowercase(fasta), ".gz")
                run(pipeline(`gunzip -c $(fasta)`, temp_fasta))
            elseif endswith(lowercase(fasta), ".bz2")
                run(pipeline(`bunzip2 -c $(fasta)`, temp_fasta))
            elseif endswith(lowercase(fasta), ".xz")
                run(pipeline(`xz -dc $(fasta)`, temp_fasta))
            end
            
            # Use the temporary uncompressed file for PGAP
            fasta_to_use = temp_fasta
            
            if as_string
                # Generate the fully expanded bash command
                bash_command = "PGAP_INPUT_DIR=$pgap_dir python3 $pgap_py_script --output $outdir --genome $fasta_to_use --report-usage-true --taxcheck --auto-correct-tax --organism '$(organism)' --cpu $(threads) --prefix $(prefix)"
                return bash_command
            else
                @time run(setenv(`python3 $(pgap_py_script) --output $(outdir) --genome $(fasta_to_use) --report-usage-true --taxcheck --auto-correct-tax --organism "$(organism)" --cpu $(threads)  --prefix $(prefix)`, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
                return outdir
            end
            
        finally
            # Clean up temporary file regardless of success or failure
            if isfile(temp_fasta)
                rm(temp_fasta)
            end
        end
    else
        # File is not compressed, use original logic
        if as_string
            # Generate the fully expanded bash command
            bash_command = "PGAP_INPUT_DIR=$pgap_dir python3 $pgap_py_script --output $outdir --genome $fasta --report-usage-true --taxcheck --auto-correct-tax --organism '$(organism)' --cpu $(threads) --prefix $(prefix)"
            return bash_command
        else
            @time run(setenv(`python3 $(pgap_py_script) --output $(outdir) --genome $(fasta) --report-usage-true --taxcheck --auto-correct-tax --organism "$(organism)" --cpu $(threads)  --prefix $(prefix)`, merge(ENV, Dict("PGAP_INPUT_DIR" => pgap_dir))))
            return outdir
        end
    end
end

function run_bakta(;
        fasta::AbstractString,
        outdir = replace(fasta, Mycelia.FASTA_REGEX => "_bakta"),
        baktadb_dir=joinpath(homedir(), "workspace", "bakta"),
        mode="full",
        threads=Sys.CPU_THREADS,
        as_string = false,
        force = false
    )
    @assert mode in ["full", "light"]
    @assert isfile(fasta)
    Mycelia.add_bioconda_env("bakta")
    bakta_db_path = joinpath(baktadb_dir, mode)
    current_db_path_contents = filter(x -> !occursin(r"^\.", x), readdir(bakta_db_path))
    if !isdir(bakta_db_path) || isempty(current_db_path_contents)
        @info "downloading bakta db"
        mkpath(bakta_db_path)
        @time run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bakta bakta_db download --output $(bakta_db_path) --type $(mode)`)
        @info "bakta db downloaded!"
    end
    bakta_db_path = joinpath(bakta_db_path, "db")
    @assert isdir(bakta_db_path)
    # --prefix ecoli123
    # --locus-tag eco634
    # --prodigal-tf eco.tf
    # --replicons replicon.tsv
    if !isdir(outdir) || isempty(filter(x -> !occursin(r"^\.", x), readdir(outdir))) || force
        if as_string
            bash_command = "$(Mycelia.CONDA_RUNNER) run --live-stream -n bakta bakta --verbose --db $(bakta_db_path) --output $(outdir) --threads $(threads) $(fasta) --force"
            return bash_command
        else
            @info "running bakta"
            @time run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n bakta bakta --verbose --db $(bakta_db_path) --output $(outdir) --threads $(threads) $(fasta) --force`)
            return outdir
        end
    else
        @warn "$(outdir) exists and is non-empty - manually remove the directory before running or use `force=true` to overwrite existing outputs"
        return outdir
    end
end

"""
Run VirSorter2 viral sequence identification tool.

VirSorter2 identifies viral sequences in genomic and metagenomic data using
machine learning models and database comparisons.

# Arguments
- `input_fasta`: Path to input FASTA file
- `output_directory`: Output directory path
- `database_path`: Path to VirSorter2 database directory
- `include_groups`: Comma-separated viral groups to include
    - full set = `dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae`
    - tools original default set = `dsDNAphage,ssDNA`
    - Lavidaviridae = A family of small double‑stranded DNA "virophages" that parasitize the replication machinery of certain NCLDVs
    - NCLDV = Nucleocytoplasmic Large DNA Viruses = An informal clade of large double‑stranded DNA viruses that replicate (at least in part) in the cytoplasm of eukaryotic cells.
- `min_score`: Minimum score threshold for viral sequences
- `min_length`: Minimum sequence length threshold
- `threads`: Number of CPU threads to use
- `provirus_off`: Disable provirus detection
- `max_orf_per_seq`: Maximum ORFs per sequence
- `prep_for_dramv`: Prepare output for DRAMv annotation
- `label`: Label for output files
- `forceall`: Force rerun all steps
- `force`: Force rerun even if output files already exist

# Returns
NamedTuple containing paths to all generated output files and directories
"""
function run_virsorter2(;
    input_fasta,
    output_directory = input_fasta * "_virsorter2",
    database_path = mkpath(joinpath(homedir(), "workspace", "virsorter2_db")),
    include_groups = "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae", # full set is dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae
    min_score = 0.5,
    min_length = 1500,
    threads = Sys.CPU_THREADS,
    provirus_off = false,
    max_orf_per_seq = nothing,
    prep_for_dramv = false,
    label = nothing,
    forceall = false,
    force = false
)
    
    # Define output file paths based on VirSorter2 structure
    final_viral_combined = joinpath(output_directory, "final-viral-combined.fa")
    final_viral_score = joinpath(output_directory, "final-viral-score.tsv")
    final_viral_boundary = joinpath(output_directory, "final-viral-boundary.tsv")
    
    # Label-specific files if label was provided
    if label !== nothing
        labeled_viral_combined = joinpath(output_directory, "$(label)-final-viral-combined.fa")
        labeled_viral_score = joinpath(output_directory, "$(label)-final-viral-score.tsv")
        labeled_viral_boundary = joinpath(output_directory, "$(label)-final-viral-boundary.tsv")
    else
        labeled_viral_combined = nothing
        labeled_viral_score = nothing
        labeled_viral_boundary = nothing
    end
    
    # DRAMv output files if prep_for_dramv was enabled
    dramv_dir = joinpath(output_directory, "for-dramv")
    dramv_viral_combined = prep_for_dramv ? joinpath(dramv_dir, "final-viral-combined-for-dramv.fa") : nothing
    dramv_affi_contigs = prep_for_dramv ? joinpath(dramv_dir, "affi-contigs.tab") : nothing
    
    # Check if output files already exist (unless force is true)
    if !force
        # Build list of expected output files to check
        expected_files = [final_viral_combined, final_viral_score, final_viral_boundary]
        
        # Add label-specific files if label was provided
        if label !== nothing
            push!(expected_files, labeled_viral_combined, labeled_viral_score, labeled_viral_boundary)
        end
        
        # Add DRAMv files if prep_for_dramv is enabled
        if prep_for_dramv
            push!(expected_files, dramv_viral_combined, dramv_affi_contigs)
        end
        
        # Filter out nothing values and check if all files exist
        files_to_check = filter(x -> x !== nothing, expected_files)
        all_files_exist = all(Base.Filesystem.isfile, files_to_check)
        
        if all_files_exist
            @warn "All VirSorter2 output files already exist in $(output_directory). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            # Intermediate directories and files that might be useful
            iter_dir = joinpath(output_directory, "iter-0")
            checkpoints_dir = joinpath(output_directory, "checkpoints")
            logs_dir = joinpath(output_directory, "logs")
            
            # Configuration and temporary files
            config_yaml = joinpath(output_directory, "config.yaml")
            
            # Return comprehensive NamedTuple with all output paths (same as below)
            return (
                # Main output directory
                output_directory = output_directory,
                database_path = database_path,
                
                # Primary output files (most commonly used)
                final_viral_combined = final_viral_combined,
                final_viral_score = final_viral_score,
                final_viral_boundary = final_viral_boundary,
                
                # Label-specific outputs (if label was provided)
                labeled_viral_combined = labeled_viral_combined,
                labeled_viral_score = labeled_viral_score,
                labeled_viral_boundary = labeled_viral_boundary,
                
                # DRAMv compatibility outputs (if enabled)
                dramv_directory = prep_for_dramv ? dramv_dir : nothing,
                dramv_viral_combined = dramv_viral_combined,
                dramv_affi_contigs = dramv_affi_contigs,
                
                # Intermediate directories
                iter_directory = iter_dir,
                checkpoints_directory = checkpoints_dir,
                logs_directory = logs_dir,
                
                # Configuration
                config_file = config_yaml,
                
                # Environment information
                conda_environment = "vs2"
            )
        end
    end
    
    # Install VirSorter2 environment with specific version convention
    env_name = "vs2"
    
    # Check if environment exists, if not create it
    env_exists = try
        run(pipeline(`$(Mycelia.CONDA_RUNNER) env list`, `grep -q "^$(env_name) "`))
        true
    catch
        false
    end
    
    if !env_exists
        println("Creating VirSorter2 conda environment...")
        run(`$(Mycelia.CONDA_RUNNER) create -n $(env_name) -c conda-forge -c bioconda virsorter=2 -y`)
        println("Setting up VirSorter2 database (this may take a while)...")
        
        # Remove incomplete database directory if it exists
        if Base.Filesystem.isdir(database_path)
            Base.Filesystem.rm(database_path, recursive=true, force=true)
        end
        
        # Setup database
        setup_cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n $(env_name) virsorter setup -d $(database_path) -j $(threads)`
        
        # Run database setup
        process = run(setup_cmd)
        
        if !success(process)
            error("VirSorter2 database setup failed. You may need to download manually from https://osf.io/v46sc/download")
        end
    end
    
    # Setup database if it doesn't exist or is incomplete
    if !Base.Filesystem.isdir(database_path) || isempty(Base.Filesystem.readdir(database_path))
        @warn "VirSorter2 database not detected"
        @warn "to install, delete any existing directory at the output location, and then run"
        @warn "`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n $(env_name) virsorter setup -d $(database_path) -j $(threads)`"
        # # Remove incomplete database directory if it exists
        # if Base.Filesystem.isdir(database_path)
        #     Base.Filesystem.rm(database_path, recursive=true, force=true)
        # end
        
        # # Setup database
        # setup_cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n $(env_name) virsorter setup -d $(database_path) -j $(threads)`
        
        # # Run database setup
        # process = run(setup_cmd)
        
        # if !success(process)
        #     error("VirSorter2 database setup failed. You may need to download manually from https://osf.io/v46sc/download")
        # end
    end
    
    # Build VirSorter2 command arguments
    cmd_args = ["$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", env_name, "virsorter", "run"]
    
    # Add required arguments
    push!(cmd_args, "-w", output_directory)
    push!(cmd_args, "-i", input_fasta)
    push!(cmd_args, "-j", string(threads))
    
    # Add optional arguments
    push!(cmd_args, "--include-groups", include_groups)
    push!(cmd_args, "--min-score", string(min_score))
    push!(cmd_args, "--min-length", string(min_length))
    
    if provirus_off
        push!(cmd_args, "--provirus-off")
    end
    
    if max_orf_per_seq !== nothing
        push!(cmd_args, "--max-orf-per-seq", string(max_orf_per_seq))
    end
    
    if prep_for_dramv
        push!(cmd_args, "--prep-for-dramv")
    end
    
    if label !== nothing
        push!(cmd_args, "--label", label)
    end
    
    if forceall
        push!(cmd_args, "--forceall")
    end
    
    # Add positional argument (default is 'all')
    push!(cmd_args, "all")
    
    # Run VirSorter2
    println("Running VirSorter2...")
    run(Cmd(cmd_args))
    
    # Intermediate directories and files that might be useful
    iter_dir = joinpath(output_directory, "iter-0")
    checkpoints_dir = joinpath(output_directory, "checkpoints")
    logs_dir = joinpath(output_directory, "logs")
    
    # Configuration and temporary files
    config_yaml = joinpath(output_directory, "config.yaml")
    
    # Return comprehensive NamedTuple with all output paths
    return (
        # Main output directory
        output_directory = output_directory,
        database_path = database_path,
        
        # Primary output files (most commonly used)
        final_viral_combined = final_viral_combined,
        final_viral_score = final_viral_score,
        final_viral_boundary = final_viral_boundary,
        
        # Label-specific outputs (if label was provided)
        labeled_viral_combined = labeled_viral_combined,
        labeled_viral_score = labeled_viral_score,
        labeled_viral_boundary = labeled_viral_boundary,
        
        # DRAMv compatibility outputs (if enabled)
        dramv_directory = prep_for_dramv ? dramv_dir : nothing,
        dramv_viral_combined = dramv_viral_combined,
        dramv_affi_contigs = dramv_affi_contigs,
        
        # Intermediate directories
        iter_directory = iter_dir,
        checkpoints_directory = checkpoints_dir,
        logs_directory = logs_dir,
        
        # Configuration
        config_file = config_yaml,
        
        # Environment information
        conda_environment = env_name
    )
end

"""
Run geNomad mobile genetic element identification tool.

geNomad identifies viruses and plasmids in genomic and metagenomic data
using machine learning and database comparisons.

# Arguments
- `input_fasta`: Path to input FASTA file
- `output_directory`: Output directory path
- `genomad_dbpath`: Path to geNomad database directory
- `threads`: Number of CPU threads to use
- `cleanup`: Remove intermediate files after completion
- `splits`: Number of splits for memory management
- `force`: Force rerun even if output files already exist

# Returns
NamedTuple containing paths to all generated output files and directories
"""
function run_genomad(;
    input_fasta,
    output_directory=input_fasta * "_genomad",
    genomad_dbpath = mkpath(joinpath(homedir(), "workspace", "genomad")),
    threads = Sys.CPU_THREADS,
    cleanup = true,
    splits = nothing,
    force = false
)
    # Get the base name for output files (remove path and extensions)
    input_basename = splitext(basename(input_fasta))[1]
    if endswith(input_basename, ".fna") || endswith(input_basename, ".fa") || endswith(input_basename, ".fasta")
        input_basename = splitext(input_basename)[1]
    end
    
    # Define expected output paths
    summary_dir = joinpath(output_directory, "$(input_basename)_summary")
    
    # Core output files
    virus_summary = joinpath(summary_dir, "$(input_basename)_virus_summary.tsv")
    plasmid_summary = joinpath(summary_dir, "$(input_basename)_plasmid_summary.tsv")
    virus_fasta = joinpath(summary_dir, "$(input_basename)_virus.fna")
    plasmid_fasta = joinpath(summary_dir, "$(input_basename)_plasmid.fna")
    virus_proteins = joinpath(summary_dir, "$(input_basename)_virus_proteins.faa")
    plasmid_proteins = joinpath(summary_dir, "$(input_basename)_plasmid_proteins.faa")
    virus_genes = joinpath(summary_dir, "$(input_basename)_virus_genes.tsv")
    plasmid_genes = joinpath(summary_dir, "$(input_basename)_plasmid_genes.tsv")
    summary_json = joinpath(summary_dir, "$(input_basename)_summary.json")
    
    # Module directories
    marker_classification_dir = joinpath(output_directory, "$(input_basename)_marker_classification")
    nn_classification_dir = joinpath(output_directory, "$(input_basename)_nn_classification")
    find_proviruses_dir = joinpath(output_directory, "$(input_basename)_find_proviruses")
    annotate_dir = joinpath(output_directory, "$(input_basename)_annotate")
    aggregated_classification_dir = joinpath(output_directory, "$(input_basename)_aggregated_classification")
    
    # Log files
    marker_classification_log = joinpath(output_directory, "$(input_basename)_marker_classification.log")
    nn_classification_log = joinpath(output_directory, "$(input_basename)_nn_classification.log")
    find_proviruses_log = joinpath(output_directory, "$(input_basename)_find_proviruses.log")
    annotate_log = joinpath(output_directory, "$(input_basename)_annotate.log")
    aggregated_classification_log = joinpath(output_directory, "$(input_basename)_aggregated_classification.log")
    summary_log = joinpath(output_directory, "$(input_basename)_summary.log")
    
    # Check if output files already exist (unless force is true)
    if !force
        # Build list of expected core output files to check
        expected_files = [
            virus_summary,
            plasmid_summary,
            virus_fasta,
            plasmid_fasta,
            virus_proteins,
            plasmid_proteins,
            virus_genes,
            plasmid_genes,
            summary_json
        ]
        
        # Check if all core files exist
        all_files_exist = all(Base.Filesystem.isfile, expected_files)
        
        if all_files_exist
            @warn "All geNomad output files already exist in $(output_directory). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            # Return NamedTuple with all paths (same as below)
            return (
                # Main output directory
                output_directory = output_directory,
                summary_directory = summary_dir,
                
                # Summary files (most commonly used)
                virus_summary = virus_summary,
                plasmid_summary = plasmid_summary,
                virus_fasta = virus_fasta,
                plasmid_fasta = plasmid_fasta,
                virus_proteins = virus_proteins,
                plasmid_proteins = plasmid_proteins,
                virus_genes = virus_genes,
                plasmid_genes = plasmid_genes,
                summary_json = summary_json,
                
                # Module directories
                marker_classification_dir = marker_classification_dir,
                nn_classification_dir = nn_classification_dir,
                find_proviruses_dir = find_proviruses_dir,
                annotate_dir = annotate_dir,
                aggregated_classification_dir = aggregated_classification_dir,
                
                # Log files
                marker_classification_log = marker_classification_log,
                nn_classification_log = nn_classification_log,
                find_proviruses_log = find_proviruses_log,
                annotate_log = annotate_log,
                aggregated_classification_log = aggregated_classification_log,
                summary_log = summary_log
            )
        end
    end
    
    # Add genomad environment
    Mycelia.add_bioconda_env("genomad")
    
    # Download database if it doesn't exist or is empty
    final_genomad_dbpath = joinpath(genomad_dbpath, "genomad_db")
    if !Base.Filesystem.isdir(final_genomad_dbpath) || isempty(Base.Filesystem.readdir(final_genomad_dbpath))
        run(`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n genomad genomad download-database $(genomad_dbpath)`)
    end
    
    # Build command arguments
    cmd_args = ["$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "genomad", "genomad", "end-to-end"]
    
    if cleanup
        push!(cmd_args, "--cleanup")
    end
    
    push!(cmd_args, "--threads", string(threads))
    
    if splits !== nothing
        push!(cmd_args, "--splits", string(splits))
    end
    
    push!(cmd_args, input_fasta, output_directory, final_genomad_dbpath)
    
    # Run genomad
    run(Cmd(cmd_args))
    
    # Return NamedTuple with all paths
    return (
        # Main output directory
        output_directory = output_directory,
        summary_directory = summary_dir,
        
        # Summary files (most commonly used)
        virus_summary = virus_summary,
        plasmid_summary = plasmid_summary,
        virus_fasta = virus_fasta,
        plasmid_fasta = plasmid_fasta,
        virus_proteins = virus_proteins,
        plasmid_proteins = plasmid_proteins,
        virus_genes = virus_genes,
        plasmid_genes = plasmid_genes,
        summary_json = summary_json,
        
        # Module directories
        marker_classification_dir = marker_classification_dir,
        nn_classification_dir = nn_classification_dir,
        find_proviruses_dir = find_proviruses_dir,
        annotate_dir = annotate_dir,
        aggregated_classification_dir = aggregated_classification_dir,
        
        # Log files
        marker_classification_log = marker_classification_log,
        nn_classification_log = nn_classification_log,
        find_proviruses_log = find_proviruses_log,
        annotate_log = annotate_log,
        aggregated_classification_log = aggregated_classification_log,
        summary_log = summary_log
    )
end

"""
    run_phispy(input_file::String; output_dir::String="", 
           phage_genes::Int=2, color::Bool=false, prefix::String="",
           phmms::String="", threads::Int=1, metrics::Vector{String}=String[],
           expand_slope::Bool=false, window_size::Int=30, 
           min_contig_size::Int=5000, skip_search::Bool=false,
           output_choice::Int=3, training_set::String="", 
           prokka_args::NamedTuple=NamedTuple(), force::Bool=false)

Run PhiSpy to identify prophages in bacterial genomes.

PhiSpy identifies prophage regions in bacterial (and archaeal) genomes using
multiple approaches including gene composition, AT/GC skew, and optional HMM searches.

# Arguments
- `input_file::String`: Path to input file (FASTA or GenBank format)
- `output_dir::String`: Output directory (default: input_file * "_phispy")  
- `phage_genes::Int`: Minimum phage genes required per prophage region (default: 2, set to 0 for mobile elements)
- `color::Bool`: Add color annotations for CDS based on function (default: false)
- `prefix::String`: Prefix for output filenames (default: basename of input)
- `phmms::String`: Path to HMM database for additional phage gene detection
- `threads::Int`: Number of threads for HMM searches (default: 1)
- `metrics::Vector{String}`: Metrics to use for prediction (default: all standard metrics)
- `expand_slope::Bool`: Expand Shannon slope calculations (default: false)
- `window_size::Int`: Window size for calculations (default: 30)
- `min_contig_size::Int`: Minimum contig size to analyze (default: 5000)
- `skip_search::Bool`: Skip HMM search if already done (default: false)
- `output_choice::Int`: Bitmask for output files (default: 3 for coordinates + GenBank)
- `training_set::String`: Path to custom training set
- `prokka_args::NamedTuple`: Additional arguments to pass to Prokka if FASTA input is provided
- `force::Bool`: Force rerun even if output files already exist (default: false)

# Output Choice Codes (add values for multiple outputs)
- 1: prophage_coordinates.tsv
- 2: GenBank format output  
- 4: prophage and bacterial sequences
- 8: prophage_information.tsv
- 16: prophage.tsv
- 32: GFF3 format (prophages only)
- 64: prophage.tbl
- 128: test data used in random forest
- 256: GFF3 format (full genome)
- 512: all output files

# Returns
A NamedTuple with paths to generated output files (contents depend on output_choice):
- `prophage_coordinates`: prophage_coordinates.tsv file path
- `genbank_output`: Updated GenBank file with prophage annotations
- `prophage_sequences`: Prophage and bacterial sequence files
- `prophage_information`: prophage_information.tsv file path
- `prophage_simple`: prophage.tsv file path
- `gff3_prophages`: GFF3 file with prophage regions only
- `prophage_table`: prophage.tbl file path
- `test_data`: Random forest test data file path
- `gff3_genome`: GFF3 file with full genome annotations
- `output_dir`: Path to output directory
- `input_genbank`: Path to GenBank file used (original or generated by Prokka)
"""
function run_phispy(input_file::String; 
                output_dir::String="",
                phage_genes::Int=2,
                color::Bool=false,
                prefix::String="",
                phmms::String="",
                threads::Int=1,
                metrics::Vector{String}=String[],
                expand_slope::Bool=false,
                window_size::Int=30,
                min_contig_size::Int=5000,
                skip_search::Bool=false,
                output_choice::Int=512,
                training_set::String="",
                prokka_args::NamedTuple=NamedTuple(),
                force::Bool=false)
    
    # Validate input file exists
    if !Base.Filesystem.isfile(input_file)
        throw(ArgumentError("Input file does not exist: $input_file"))
    end
    
    # Set default output directory
    if Base.isempty(output_dir)
        output_dir = input_file * "_phispy"
    end
    
    # Set default prefix from input filename if not provided
    if Base.isempty(prefix)
        prefix = Base.Filesystem.splitext(Base.Filesystem.basename(input_file))[1]
    end
    
    # Build expected output file paths based on output_choice
    output_files = Dict{Symbol, String}()
    
    # Map output choice bits to file paths
    output_mapping = [
        (1, :prophage_coordinates, "prophage_coordinates.tsv"),
        (2, :genbank_output, prefix * ".gbk"),
        (4, :prophage_sequences, prefix * "_phage.fasta"), # Multiple files, returning one representative
        (8, :prophage_information, "prophage_information.tsv"),
        (16, :prophage_simple, "prophage.tsv"),
        (32, :gff3_prophages, prefix * "_prophage.gff3"),
        (64, :prophage_table, "prophage.tbl"),
        (128, :test_data, "test_data.tsv"),
        (256, :gff3_genome, prefix * ".gff3")
    ]
    
    # Check which outputs were requested and add paths
    expected_files = String[]
    for (bit, key, filename) in output_mapping
        if (output_choice & bit) != 0
            file_path = Base.Filesystem.joinpath(output_dir, filename)
            output_files[key] = file_path
            Base.push!(expected_files, file_path)
        end
    end
    
    # Check if output files already exist (unless force is true)
    if !force && Base.Filesystem.isdir(output_dir)
        # Check if all expected files exist
        # all_files_exist = Base.all(Base.Filesystem.isfile, expected_files)
        
        # if all_files_exist && !Base.isempty(expected_files)
        if !isempty(readdir(output_dir))
            @warn "PhiSpy output files already exist in $(output_dir). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            # Add standard output files
            output_files[:output_dir] = Base.Filesystem.abspath(output_dir)
            
            # For existing runs, we need to determine the input genbank path
            # Check if there's a Prokka temp directory
            # prokka_temp_dir = output_dir * "_prokka_temp"
            # if Base.Filesystem.isdir(prokka_temp_dir)
            #     prokka_prefix = prefix * "_prokka"
            #     potential_genbank = Base.Filesystem.joinpath(prokka_temp_dir, prokka_prefix * ".gbk")
            #     if Base.Filesystem.isfile(potential_genbank)
            #         output_files[:input_genbank] = Base.Filesystem.abspath(potential_genbank)
            #     else
            #         output_files[:input_genbank] = Base.Filesystem.abspath(input_file)
            #     end
            # else
            #     output_files[:input_genbank] = Base.Filesystem.abspath(input_file)
            # end
            
            return NamedTuple(output_files)
        end
    end
    
    # Determine input file type and prepare GenBank file
    input_genbank = ""
    file_extension = Base.lowercase(Base.Filesystem.splitext(input_file)[2])
    
    if file_extension in [".fasta", ".fa", ".fna", ".fas"]
        # Input is FASTA - need to run Prokka first
        @info "FASTA input detected. Running Prokka for annotation..."
        
        # Ensure bioconda environment is set up
        Mycelia.add_bioconda_env("phispy")
        
        prokka_output_dir = output_dir * "_prokka_temp"
        prokka_prefix = prefix * "_prokka"
        
        # Set up Prokka arguments
        prokka_kwargs = Dict{Symbol, Any}(
            :output_dir => prokka_output_dir,
            :prefix => prokka_prefix
        )
        
        # Add any additional Prokka arguments provided
        for (key, value) in Base.pairs(prokka_args)
            prokka_kwargs[key] = value
        end
        
        # Run Prokka
        prokka_results = Mycelia.run_prokka(input_file; prokka_kwargs...)
        input_genbank = prokka_results.gbk
        
        @info "Prokka completed. Using generated GenBank file: $input_genbank"
        
    elseif file_extension in [".gb", ".gbk", ".genbank", ".gbf"] || 
           (file_extension == ".gz" && Base.any(x -> Base.occursin(x, Base.lowercase(input_file)), [".gb", ".gbk", ".genbank", ".gbf"]))
        # Input is already GenBank format
        input_genbank = input_file
        @info "GenBank input detected: $input_genbank"
    else
        throw(ArgumentError("Unsupported file format. Please provide FASTA (.fasta, .fa, .fna, .fas) or GenBank (.gb, .gbk, .genbank, .gbf) files."))
    end
    
    # Ensure bioconda environment is set up
    Mycelia.add_bioconda_env("phispy")
    
    # Create output directory if it doesn't exist
    if !Base.Filesystem.isdir(output_dir)
        Base.Filesystem.mkdir(output_dir)
    end
    
    # Build PhiSpy command
    cmd_args = String[
        "$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "phispy", "PhiSpy.py",
        input_genbank,
        "-o", output_dir
    ]
    
    # Add optional arguments
    if phage_genes != 2
        Base.push!(cmd_args, "--phage_genes", Base.string(phage_genes))
    end
    
    if color
        Base.push!(cmd_args, "--color")
    end
    
    # if !Base.isempty(prefix)
    #     Base.push!(cmd_args, "--prefix", prefix)
    # end
    
    if !Base.isempty(phmms)
        Base.push!(cmd_args, "--phmms", phmms)
    end
    
    if threads != 1
        Base.push!(cmd_args, "--threads", Base.string(threads))
    end
    
    if !Base.isempty(metrics)
        Base.push!(cmd_args, "--metrics")
        Base.append!(cmd_args, metrics)
    end
    
    if expand_slope
        Base.push!(cmd_args, "--expand_slope")
    end
    
    if window_size != 30
        Base.push!(cmd_args, "--window_size", Base.string(window_size))
    end
    
    if min_contig_size != 5000
        Base.push!(cmd_args, "--min_contig_size", Base.string(min_contig_size))
    end
    
    if skip_search
        Base.push!(cmd_args, "--skip_search")
    end
    
    if output_choice != 3
        Base.push!(cmd_args, "--output_choice", Base.string(output_choice))
    end
    
    if !Base.isempty(training_set)
        Base.push!(cmd_args, "-t", training_set)
    end
    
    # Run PhiSpy
    try
        @info "Running PhiSpy with command: $(Base.join(cmd_args, " "))"
        Base.run(Base.Cmd(cmd_args))
    catch e
        throw(ErrorException("PhiSpy failed to run: $e"))
    end
    
    # Verify output directory exists
    if !Base.Filesystem.isdir(output_dir)
        throw(ErrorException("PhiSpy output directory was not created: $output_dir"))
    end
    
    # Add standard output files
    output_files[:output_dir] = Base.Filesystem.abspath(output_dir)
    # output_files[:input_genbank] = Base.Filesystem.abspath(input_genbank)
    
    # # Check if key output files exist and warn if missing
    # key_files = [
    #     Base.get(output_files, :prophage_coordinates, ""),
    #     Base.get(output_files, :genbank_output, ""),
    #     Base.get(output_files, :prophage_information, "")
    # ]
    
    # existing_files = Base.filter(f -> !Base.isempty(f) && Base.Filesystem.isfile(f), key_files)
    # missing_files = Base.filter(f -> !Base.isempty(f) && !Base.Filesystem.isfile(f), key_files)
    
    # if !Base.isempty(missing_files)
    #     @warn "Some expected output files were not created: $(Base.join(missing_files, ", "))"
    # end
    
    # if Base.isempty(existing_files)
    #     @warn "No standard PhiSpy output files were found. Check PhiSpy logs for errors."
    # else
    #     @info "PhiSpy completed successfully. Generated $(Base.length(existing_files)) output files."
    # end

    @assert isdir(output_files[:output_dir]) && !isempty(readdir(output_files[:output_dir]))
    
    return NamedTuple(output_files)
end

"""
    run_prokka(input_fasta::String; output_dir::String="", prefix::String="", 
           cpus::Int=0, kingdom::String="Bacteria", genus::String="", 
           species::String="", strain::String="", force_overwrite::Bool=false,
           addgenes::Bool=false, compliant::Bool=false, fast::Bool=false,
           evalue::Float64=1e-06, mincontiglen::Int=1, force::Bool=false)

Run Prokka for rapid prokaryotic genome annotation.

Prokka annotates bacterial, archaeal and viral genomes quickly and produces 
standards-compliant output files including GFF3, GenBank, and FASTA formats.

# Arguments
- `input_fasta::String`: Path to input FASTA file containing contigs
- `output_dir::String`: Output directory (default: input_fasta * "_prokka")
- `prefix::String`: Output file prefix (default: basename of input file)
- `cpus::Int`: Number of CPUs to use, 0 for all available (default: 0)
- `kingdom::String`: Annotation mode - "Bacteria", "Archaea", "Viruses", or "Mitochondria" (default: "Bacteria")
- `genus::String`: Genus name for annotation
- `species::String`: Species name for annotation  
- `strain::String`: Strain name for annotation
- `force_overwrite::Bool`: Force overwrite existing output directory (default: false)
- `addgenes::Bool`: Add 'gene' features for each 'CDS' feature (default: false)
- `compliant::Bool`: Force GenBank/ENA/DDJB compliance (default: false)
- `fast::Bool`: Fast mode - skip CDS product searching (default: false)
- `evalue::Float64`: Similarity e-value cut-off (default: 1e-06)
- `mincontiglen::Int`: Minimum contig size (default: 1, NCBI needs 200)
- `force::Bool`: Force rerun even if output files already exist (default: false)

# Returns
A NamedTuple with paths to all generated output files:
- `gff`: Master annotation in GFF3 format
- `gbk`: Standard GenBank file  
- `fna`: Nucleotide FASTA of input contigs
- `faa`: Protein FASTA of translated CDS sequences
- `ffn`: Nucleotide FASTA of all transcripts
- `sqn`: ASN1 Sequin file for GenBank submission
- `fsa`: Nucleotide FASTA for tbl2asn
- `tbl`: Feature table file
- `err`: NCBI discrepancy report
- `log`: Complete run log
- `txt`: Annotation statistics
- `tsv`: Tab-separated feature table
- `output_dir`: Path to output directory
"""
function run_prokka(input_fasta::String; 
                output_dir::String="",
                prefix::String="",
                cpus::Int=0,
                kingdom::String="Bacteria",
                genus::String="",
                species::String="", 
                strain::String="",
                force_overwrite::Bool=false,
                addgenes::Bool=false,
                compliant::Bool=false,
                fast::Bool=false,
                evalue::Float64=1e-06,
                mincontiglen::Int=1,
                force::Bool=false)
    
    # Validate input file exists
    if !Base.Filesystem.isfile(input_fasta)
        Base.throw(ArgumentError("Input FASTA file does not exist: $input_fasta"))
    end
    
    # Set default output directory
    if Base.isempty(output_dir)
        output_dir = input_fasta * "_prokka"
    end
    
    # Set default prefix from input filename if not provided
    if Base.isempty(prefix)
        prefix = Base.Filesystem.splitext(Base.Filesystem.basename(input_fasta))[1]
    end
    
    # Build output file paths
    output_files = (
        gff = Base.Filesystem.joinpath(output_dir, prefix * ".gff"),
        gbk = Base.Filesystem.joinpath(output_dir, prefix * ".gbk"), 
        fna = Base.Filesystem.joinpath(output_dir, prefix * ".fna"),
        faa = Base.Filesystem.joinpath(output_dir, prefix * ".faa"),
        ffn = Base.Filesystem.joinpath(output_dir, prefix * ".ffn"),
        sqn = Base.Filesystem.joinpath(output_dir, prefix * ".sqn"),
        fsa = Base.Filesystem.joinpath(output_dir, prefix * ".fsa"),
        tbl = Base.Filesystem.joinpath(output_dir, prefix * ".tbl"),
        err = Base.Filesystem.joinpath(output_dir, prefix * ".err"),
        log = Base.Filesystem.joinpath(output_dir, prefix * ".log"), 
        txt = Base.Filesystem.joinpath(output_dir, prefix * ".txt"),
        tsv = Base.Filesystem.joinpath(output_dir, prefix * ".tsv"),
        output_dir = Base.Filesystem.abspath(output_dir)
    )
    
    # Check if output files already exist (unless force is true)
    if !force
        # Build list of expected core output files to check
        expected_files = [
            output_files.gff,
            output_files.gbk,
            output_files.fna,
            output_files.faa,
            output_files.ffn,
            output_files.txt,
            output_files.tsv
        ]
        
        # Check if all core files exist
        all_files_exist = Base.all(Base.Filesystem.isfile, expected_files)
        
        if all_files_exist
            @warn "All Prokka output files already exist in $(output_dir). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
            
            return output_files
        end
    end
    
    # Ensure bioconda environment is set up
    Mycelia.add_bioconda_env("prokka")
    
    # Build prokka command
    cmd_args = String[
        "$(Mycelia.CONDA_RUNNER)", "run", "--no-capture-output", "-n", "prokka", "prokka",
        "--outdir", output_dir,
        "--prefix", prefix,
        "--kingdom", kingdom
    ]
    
    # Add optional arguments
    if cpus > 0
        Base.push!(cmd_args, "--cpus", Base.string(cpus))
    elseif cpus == 0
        Base.push!(cmd_args, "--cpus", "0")  # Use all available CPUs
    end
    
    if !Base.isempty(genus)
        Base.push!(cmd_args, "--genus", genus)
    end
    
    if !Base.isempty(species) 
        Base.push!(cmd_args, "--species", species)
    end
    
    if !Base.isempty(strain)
        Base.push!(cmd_args, "--strain", strain)
    end
    
    if force_overwrite
        Base.push!(cmd_args, "--force")
    end
    
    if addgenes
        Base.push!(cmd_args, "--addgenes")  
    end
    
    if compliant
        Base.push!(cmd_args, "--compliant")
    end
    
    if fast
        Base.push!(cmd_args, "--fast")
    end
    
    if evalue != 1e-06
        Base.push!(cmd_args, "--evalue", Base.string(evalue))
    end
    
    if mincontiglen != 1
        Base.push!(cmd_args, "--mincontiglen", Base.string(mincontiglen))
    end
    
    # Add input file
    Base.push!(cmd_args, Base.Filesystem.abspath(input_fasta))
    
    # Run prokka
    try
        Base.run(Base.Cmd(cmd_args))
    catch e
        Base.throw(ErrorException("Prokka failed to run: $e"))
    end
    
    # Verify output directory was created
    if !Base.Filesystem.isdir(output_dir)
        Base.throw(ErrorException("Prokka output directory was not created: $output_dir"))
    end
    
    # Verify key output files exist
    key_files = [output_files.gff, output_files.gbk, output_files.fna]
    missing_files = Base.filter(f -> !Base.Filesystem.isfile(f), key_files)
    
    if !Base.isempty(missing_files)
        @warn "Some expected output files were not created: $(Base.join(missing_files, ", "))"
    end
    
    return output_files
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run AMRFinderPlus on FASTA input to identify antimicrobial resistance genes.

# Arguments
- `fasta::String`: Path to input FASTA file (must match Mycelia.FASTA_REGEX pattern)
- `output_dir::String`: Output directory path (default: input filename + "_amrfinderplus")
- `force::Bool`: Force rerun even if output files already exist (default: false)

# Returns
Path to the output directory containing AMRFinderPlus results

# Details
- For nucleotide FASTA files, automatically runs Mycelia.run_pyrodigal to generate protein sequences
- For protein FASTA files, runs AMRFinderPlus directly  
- Validates input file extension against Mycelia.FASTA_REGEX
- Creates output directory if it doesn't exist
- Skips processing if results already exist in output directory unless force=true
- Uses --plus flag for enhanced detection capabilities

# Files Generated
- `<basename>.amrfinderplus.tsv`: AMRFinderPlus results table
- For nucleotide inputs: intermediate pyrodigal outputs in subdirectory
"""
function run_amrfinderplus(;
        fasta::String,
        output_dir::String = fasta * "_amrfinderplus",
        force::Bool = false
    )
    
    # Validate input file extension
    if !Base.occursin(Mycelia.FASTA_REGEX, fasta)
        Base.error("Input file does not match FASTA format: $(fasta)")
    end
    
    if !Base.Filesystem.isfile(fasta)
        Base.error("Input FASTA file not found: $(fasta)")
    end
    
    # Get base filename for outputs
    base_name = Base.replace(Base.Filesystem.basename(fasta), Mycelia.FASTA_REGEX => "")
    amrfinder_output = Base.Filesystem.joinpath(output_dir, "$(base_name).amrfinderplus.tsv")
    
    # Check if output files already exist (unless force is true)
    if !force
        if Base.Filesystem.isfile(amrfinder_output)
            @warn "AMRFinderPlus output file already exists: $(amrfinder_output). Skipping analysis."
            @warn "Use `force=true` to rerun anyway."
                return (;output_dir, amrfinder_output)
        end
    end
    
    # Create output directory
    if !Base.Filesystem.isdir(output_dir)
        Base.Filesystem.mkpath(output_dir)
    end
    
    # Determine sequence type by checking file extension first, then by detection
    sequence_type = :unknown
    
    # Check common file extensions first for efficiency
    fasta_lower = Base.lowercase(fasta)
    if Base.endswith(fasta_lower, ".faa") || Base.endswith(fasta_lower, ".faa.gz")
        sequence_type = :protein
    elseif Base.endswith(fasta_lower, ".fna") || Base.endswith(fasta_lower, ".fna.gz")
        sequence_type = :nucleotide
    else
        # For generic extensions (.fa, .fasta, etc.), detect from sequence content
        @info "Generic FASTA extension detected, analyzing sequence content to determine type"
        for (i, record) in Base.enumerate(Mycelia.open_fastx(fasta))
            if i > 3 Base.break end  # Sample first 3 sequences
            seq_ext = Mycelia.detect_sequence_extension(record)
            if seq_ext == ".faa"
                sequence_type = :protein
                Base.break
            elseif seq_ext == ".fna"
                sequence_type = :nucleotide
                Base.break
            end
        end
        
        if sequence_type == :unknown
            Base.error("Could not determine sequence type for: $(fasta)")
        end
    end
    
    # Determine protein FASTA file to use based on sequence type
    protein_fasta = ""
    
    if sequence_type == :protein
        # Input is already protein - use directly
        protein_fasta = fasta
        @info "Using protein FASTA directly: $(fasta)"
    elseif sequence_type == :nucleotide
        # Input is nucleotide - need to run pyrodigal first
        @info "Nucleotide FASTA detected, running pyrodigal to generate protein sequences"
        pyrodigal_dir = Base.Filesystem.joinpath(output_dir, "pyrodigal")
        pyrodigal_results = Mycelia.run_pyrodigal(fasta_file=fasta, out_dir=pyrodigal_dir)
        protein_fasta = pyrodigal_results.faa
    else
        Base.error("Unsupported sequence type: $(sequence_type)")
    end
    
    Base.@assert Base.Filesystem.isfile(protein_fasta) "Protein FASTA file not found: $(protein_fasta)"
    
    # Run AMRFinderPlus
    @info "Running AMRFinderPlus on protein sequences: $(protein_fasta)"
    
    # Ensure AMRFinderPlus is available
    Mycelia.add_bioconda_env("ncbi-amrfinderplus")
    Base.run(`$(Mycelia.CONDA_RUNNER) run --no-capture-output -n ncbi-amrfinderplus amrfinder -u`)
    
    cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n ncbi-amrfinderplus amrfinder
           -p $(protein_fasta)
           --plus
           --output $(amrfinder_output)`
    
    Base.run(cmd)
    
    Base.@assert Base.Filesystem.isfile(amrfinder_output) "AMRFinderPlus output not generated: $(amrfinder_output)"
    @info "AMRFinderPlus completed successfully. Results: $(amrfinder_output)"
    
    return (;output_dir, amrfinder_output)
end

# VIBRANT
function run_vibrant(;input_fasta, output_dir=input_fasta * "_vibrant", threads=Sys.CPU_THREADS)
    Mycelia._install_vibrant()
    # mkpath(output_dir)
    # run(`bash -lc "VIBRANT_run.py -i $input_fasta -folder $output_dir"`)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n vibrant VIBRANT_run.py -i $input_fasta -folder $output_dir -t $threads`)
    # phages_circular.fna
    # integrated_prophage_coordinates.tsv
    # figure_PCA.pdf
    # figure_PCA.tsv
    # summary_normalized.tsv
    # for more information
    # https://github.com/AnantharamanLab/VIBRANT?tab=readme-ov-file#output-explanations--
    return output_dir
end


"""
    run_phageboost(input_fasta::AbstractString, output_dir::AbstractString; force_reinstall::Bool=false)

Run PhageBoost on the provided FASTA file, automatically handling conda environment setup.

This function will:
1. Check if the phageboost_env conda environment exists
2. Create and set up the environment if it doesn't exist
3. Validate that PhageBoost is properly installed
4. Run PhageBoost on the input FASTA file
5. Return the output directory path and list of generated files

# Arguments
- `input_fasta::AbstractString`: Path to the input FASTA file
- `output_dir::AbstractString`: Directory where PhageBoost outputs will be saved
- `force_reinstall::Bool=false`: If true, recreate the environment even if it exists

# Returns
- `NamedTuple` with fields:
  - `output_dir::String`: Path to the output directory
  - `files::Vector{String}`: List of files generated in the output directory
"""
function run_phageboost(;input_fasta::AbstractString, output_dir::AbstractString=input_fasta * "_phageboost", force_reinstall::Bool=false)
    # Check if input file exists
    if !Base.isfile(input_fasta)
        throw(ArgumentError("Input FASTA file does not exist: $input_fasta"))
    end
    
    # Setup environment if needed
    _setup_phageboost_environment(force_reinstall)
    
    # # Validate PhageBoost installation
    # _validate_phageboost_installation()
    
    # Ensure output directory exists
    Base.mkpath(output_dir)
    
    # Run PhageBoost
    println("Running PhageBoost on $input_fasta...")
    try
        # Base.run(`bash -lc "conda activate phageboost_env && PhageBoost -f $input_fasta -o $output_dir"`)
        cmd = `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n phageboost_env PhageBoost -f $input_fasta -o $output_dir`
        display(cmd)
        run(cmd)
        println("PhageBoost completed successfully")
    catch e
        throw(ErrorException("PhageBoost execution failed: $e"))
    end
    
    # Get list of output files
    output_files = _get_output_files(output_dir)
    
    return (output_dir=output_dir, files=output_files)
end


"""
    parallel_pyrodigal(normalized_fastas::Vector{String})

Runs Mycelia.run_pyrodigal on a list of FASTA files in parallel using Threads.

Args:
    normalized_fastas: A vector of strings, where each string is a path to a FASTA file.

Returns:
    A tuple containing two elements:
    1. successes (Vector{Tuple{String, Any}}): A vector of tuples, where each tuple contains the
       filename and the result returned by a successful Mycelia.run_pyrodigal call.
    2. failures (Vector{Tuple{String, String}}): A vector of tuples, where each tuple contains the
       filename and the error message string for a failed Mycelia.run_pyrodigal call.
"""
function parallel_pyrodigal(normalized_fastas::Vector{String})
    num_files = Base.length(normalized_fastas)
    Base.println("Processing $(num_files) FASTA files using $(Threads.nthreads()) threads...")

    # Create a Progress object for manual updates
    p = ProgressMeter.Progress(num_files, 1, "Running Pyrodigal: ", 50)

    # Use Channels to collect results and failures thread-safely
    # Channel{Tuple{Filename, ResultType}} - adjust ResultType if known
    successes = Base.Channel{Tuple{String, Any}}(num_files)
    failures = Base.Channel{Tuple{String, String}}(num_files)

    # Use Threads.@threads for parallel execution
    Threads.@threads for fasta_file in normalized_fastas
        result = nothing # Initialize result variable in the loop's scope
        try
            # --- Execute the function ---
            # Base.println("Thread $(Threads.threadid()) processing: $(fasta_file)") # Optional: for debugging
            result = Mycelia.run_pyrodigal(fasta_file = fasta_file) # Capture the result

            # --- Store success ---
            Base.put!(successes, (fasta_file, result))

        catch e
            # --- Store failure ---
            err_msg = Base.sprint(Base.showerror, e) # Get the error message as a string
            Base.println(Base.stderr, "ERROR processing $(fasta_file) on thread $(Threads.threadid()): $(err_msg)")
            Base.put!(failures, (fasta_file, err_msg))
        finally
            # --- Always update progress ---
            ProgressMeter.next!(p)
        end
    end

    # Close channels now that all threads are done writing
    Base.close(successes)
    Base.close(failures)

    # Collect results and failures from the channels
    successful_results = Base.collect(successes)
    failed_files = Base.collect(failures)

    # --- Report Summary ---
    Base.println("\n--- Pyrodigal Processing Summary ---")
    num_success = Base.length(successful_results)
    num_failed = Base.length(failed_files)
    Base.println("Successfully processed: $(num_success)")
    Base.println("Failed: $(num_failed)")

    if !Base.isempty(failed_files)
        Base.println("\nFailures:")
        for (file, err) in failed_files
            Base.println("- File: $(file)\n  Error: $(err)")
        end
    end
    Base.println("------------------------------------")

    return successful_results, failed_files # Return both successes and failures
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Perform comprehensive annotation of a FASTA file including gene prediction, protein homology search,
# and terminator prediction.

# # Arguments
# - `fasta::String`: Path to input FASTA file
# - `identifier::String`: Unique identifier for output directory (default: FASTA filename without extension)
# - `basedir::String`: Base directory for output (default: current working directory)
# - `mmseqsdb::String`: Path to MMseqs2 UniRef50 database (default: "~/workspace/mmseqs/UniRef50")
# - `threads::Int`: Number of CPU threads to use (default: all available)

# # Processing Steps
# 1. Creates output directory and copies input FASTA
# 2. Runs Prodigal for gene prediction (nucleotide, amino acid, and GFF output)
# 3. Performs MMseqs2 homology search against UniRef50
# 4. Predicts terminators using TransTerm
# 5. Combines annotations into a unified GFF file
# 6. Generates GenBank format output

# # Returns
# - `String`: Path to the output directory containing all generated files

# # Files Generated
# - `.prodigal.fna`: Predicted genes (nucleotide)
# - `.prodigal.faa`: Predicted proteins
# - `.prodigal.gff`: Prodigal GFF annotations
# - `.gff`: Combined annotations
# - `.gff.genbank`: Final GenBank format
# """
# function annotate_fasta(;
#         fasta,
#         identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
#         basedir = pwd(),        
#         mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
#         threads=Sys.CPU_THREADS
#     )
#     # @show basedir
#     outdir = joinpath(basedir, identifier)
#     @assert outdir != fasta
#     # if !isdir(outdir)
#     #     @show isdir(outdir)
#     mkpath(outdir)
#     f = joinpath(outdir, basename(fasta))
#     # make this an rclone copy for portability
#     !isfile(f) && cp(fasta, f)
#     nucleic_acid_fasta = f * ".prodigal.fna"
#     amino_acid_fasta = f * ".prodigal.faa"
#     gff_file = f * ".prodigal.gff"
#     if !isfile(nucleic_acid_fasta) || !isfile(amino_acid_fasta) || !isfile(gff_file)
#         Mycelia.run_prodigal(fasta_file=f)
#     end
#     mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
#     mmseqs_gff_file = Mycelia.write_gff(gff = Mycelia.update_gff_with_mmseqs(gff_file, mmseqs_outfile), outfile = mmseqs_outfile * ".gff")
#     transterm_gff_file = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
#     joint_gff = Mycelia.write_gff(
#         gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file)), ["#seqid", "start", "end"]),
#         outfile=f * ".gff")
#     annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff, genbank = joint_gff * ".genbank")

#     transterm_gff_file_raw_fasta = Mycelia.transterm_output_to_gff(Mycelia.run_transterm(fasta=f))
#     joint_gff_raw_fasta = Mycelia.write_gff(
#         gff=sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file_raw_fasta)), ["#seqid", "start", "end"]),
#         outfile=f * ".transterm_raw.gff")
#     annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f, gff=joint_gff_raw_fasta, genbank = joint_gff_raw_fasta * ".genbank")
#     # else
#     #     @info "$(outdir) already present, skipping..."
#     # end
#     return outdir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform basic annotation of a FASTA file including gene prediction, protein homology search,
and terminator prediction applicable for phage and bacteria.

# Arguments
- `fasta::String`: Path to input FASTA file
- `identifier::String`: Unique identifier for output directory (default: FASTA filename without extension)
- `basedir::String`: Base directory for output (default: current working directory)
- `mmseqsdb::String`: Path to MMseqs2 UniRef50 database (default: `joinpath(homedir(), "workspace/mmseqs/UniRef50")`)
- `threads::Int`: Number of CPU threads to use (default: all available). Note: This argument is not explicitly used by Pyrodigal or MMseqs2 in this version of the function, they might use their own defaults or require modifications to `run_pyrodigal` or `run_mmseqs_easy_search` to respect it.

# Processing Steps
1. Creates output directory and copies input FASTA.
2. Runs Pyrodigal for gene prediction (nucleotide, amino acid, and GFF output).
3. Performs MMseqs2 homology search against UniRef50.
4. Predicts terminators using TransTerm.
5. Combines annotations into a unified GFF file.
6. Generates GenBank format output.

# Returns
- `String`: Path to the output directory containing all generated files.

# Files Generated (within the output directory specified by `identifier`)
- `(basename(fasta)).pyrodigal.fna`: Predicted genes (nucleotide) from Pyrodigal.
- `(basename(fasta)).pyrodigal.faa`: Predicted proteins from Pyrodigal.
- `(basename(fasta)).pyrodigal.gff`: Pyrodigal GFF annotations.
- `(basename(fasta)).gff`: Combined GFF annotations (MMseqs2 and TransTerm).
- `(basename(fasta)).gff.genbank`: Final GenBank format from the first combined GFF.
- `(basename(fasta)).transterm_raw.gff`: Combined GFF (MMseqs2 and a second TransTerm run).
- `(basename(fasta)).transterm_raw.gff.genbank`: Final GenBank format from the second combined GFF.
"""
function annotate_fasta(;
        fasta::String,
        # identifier::String = replace(basename(fasta), Mycelia.FASTA_REGEX => ""), # Assuming Mycelia.FASTA_REGEX is defined
        # basedir::String = pwd(),      
        mmseqsdb::String = joinpath(homedir(), "workspace/mmseqs/UniRef50"),
        threads::Int = Sys.CPU_THREADS,
        outdir::AbstractString = replace(fasta, Mycelia.FASTA_REGEX => "") * "_annotation"
    )
    
    # outdir = joinpath(basedir, identifier)
    # @assert outdir != fasta "Output directory cannot be the same as the input FASTA file path."

    if isdir(outdir) && !isempty(readdir(outdir))
        @warn "$outdir already exists, remove to regenerate"
        return outdir
    end
    mkpath(outdir)
    
    # Path to the FASTA file copied into the output directory
    f_in_outdir = joinpath(outdir, basename(fasta))
    # Copy input FASTA to output directory if it's not already there or needs update
    if !isfile(f_in_outdir) || mtime(fasta) > mtime(f_in_outdir)
        cp(fasta, f_in_outdir, force=true)
    end

    # --- Gene Prediction using Pyrodigal ---
    # Call run_pyrodigal, assuming it's part of the Mycelia module.
    # run_pyrodigal handles its own output file naming and existence checks.
    # We direct its output to be within our main `outdir`.
    pyrodigal_outputs = Mycelia.run_pyrodigal(fasta_file=f_in_outdir, out_dir=outdir)
    
    nucleic_acid_fasta = pyrodigal_outputs.fna
    amino_acid_fasta = pyrodigal_outputs.faa
    gff_file_pyrodigal = pyrodigal_outputs.gff # Renamed to avoid confusion with later gff_file variables
    # --- End of Pyrodigal section ---

    mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
    # Update GFF with MMseqs results, using Pyrodigal's GFF as base
    mmseqs_gff_file = Mycelia.write_gff(gff = Mycelia.update_gff_with_mmseqs(gff_file_pyrodigal, mmseqs_outfile), outfile = mmseqs_outfile * ".gff")
    
    # Predict terminators using TransTerm on the original sequence copy
    if occursin(r"\.gz$", f_in_outdir)
        @warn "transterm doesn't seem to work with gzip compressed fasta files"
    end
    transterm_results = Mycelia.run_transterm(fasta=f_in_outdir) # Assuming run_transterm returns path or object usable by transterm_output_to_gff
    transterm_gff_file = Mycelia.transterm_output_to_gff(transterm_results)
    
    # Combine MMseqs and TransTerm GFFs
    combined_gff_path = joinpath(outdir, basename(f_in_outdir) * ".gff")
    joint_gff_df = DataFrames.sort!(DataFrames.vcat(Mycelia.read_gff(mmseqs_gff_file), Mycelia.read_gff(transterm_gff_file)), ["#seqid", "start", "end"])
    Mycelia.write_gff(gff=joint_gff_df, outfile=combined_gff_path)
    
    genbank_path = combined_gff_path * ".genbank"
    annotated_genbank = Mycelia.fasta_and_gff_to_genbank(fasta=f_in_outdir, gff=combined_gff_path, genbank=genbank_path)
    
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Annotate amino acid sequences in a FASTA file using MMseqs2 search against UniRef50 database.

# Arguments
- `fasta`: Path to input FASTA file containing amino acid sequences
- `identifier`: Name for the output directory (defaults to FASTA filename without extension)
- `basedir`: Base directory for output (defaults to current directory)
- `mmseqsdb`: Path to MMseqs2 formatted UniRef50 database (defaults to ~/workspace/mmseqs/UniRef50)
- `threads`: Number of CPU threads to use (defaults to system thread count)

# Returns
- Path to the output directory containing MMseqs2 search results

The function creates a new directory named by `identifier` under `basedir`, copies the input FASTA file,
and runs MMseqs2 easy-search against the specified database. If the output directory already exists,
the function skips processing and returns the directory path.
"""
function annotate_aa_fasta(;
        fasta,
        identifier = replace(basename(fasta), Mycelia.FASTA_REGEX => ""),
        basedir = pwd(),
        mmseqsdb = "$(homedir())/workspace/mmseqs/UniRef50",
        threads=Sys.CPU_THREADS
    )
    # @show basedir
    outdir = joinpath(basedir, identifier)
    @assert outdir != fasta
    if !isdir(outdir)
        @show isdir(outdir)
        mkpath(outdir)
        f = joinpath(outdir, basename(fasta))
        # make this an rclone copy for portability
        cp(fasta, f, force=true)
        amino_acid_fasta = f

        mmseqs_outfile = Mycelia.run_mmseqs_easy_search(query_fasta=amino_acid_fasta, target_database=mmseqsdb)
    else
        @info "$(outdir) already present, skipping..."
    end
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensure the `padloc` environment and database are installed.

Downloads the environment if missing and updates the padloc database.
"""
function setup_padloc()
    padloc_is_already_installed = check_bioconda_env_is_installed("padloc")
    if !padloc_is_already_installed
        Mycelia.add_bioconda_env("padlocbio::padloc")
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --db-update`)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the 'padloc' tool from the 'padlocbio' conda environment on a given FASTA file.

https://doi.org/10.1093/nar/gkab883

https://github.com/padlocbio/padloc

This function first ensures that the 'padloc' environment is available via Bioconda. 
It then attempts to update the 'padloc' database. 
If a 'padloc' output file (with a '_padloc.csv' suffix) does not already exist for the input FASTA file, 
it runs 'padloc' with the specified FASTA file as input.

If the input FASTA is compressed (recognized by extensions .gz, .bz2, .xz, .zip) this function
creates a temporary uncompressed copy, runs padloc on that temporary file, and then removes the
temporary files and directory after the run completes.
"""
function run_padloc(;fasta_file, outdir=replace(fasta_file, Mycelia.FASTA_REGEX => "") * "_padloc", threads=Sys.CPU_THREADS)
    padloc_outfile = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "") * "_padloc.csv")
    if !isfile(padloc_outfile)
        setup_padloc()
        if !isdir(outdir)
            mkpath(outdir)
        end

        # Detect common compressed FASTA extensions
        compressed_ext_regex = r"\.(gz|bz2|xz|zip)$"i
        is_compressed = occursin(compressed_ext_regex, fasta_file)

        # Prepare variables for optional temporary uncompressed FASTA
        temp_dir = ""
        temp_fasta = ""
        input_fasta = fasta_file

        if is_compressed
            # Create a temporary directory to hold the uncompressed FASTA
            temp_dir = mktempdir()
            # Strip only the compression extension for the temp filename
            temp_fname = replace(basename(fasta_file), compressed_ext_regex => "")
            temp_fasta = joinpath(temp_dir, temp_fname)
            @info "Detected compressed FASTA; creating temporary uncompressed copy at $(temp_fasta)"

            # Decompress into the temporary file using the system's decompressors
            open(temp_fasta, "w") do io
                lf = lowercase(fasta_file)
                if endswith(lf, ".gz")
                    run(pipeline(`gzip -dc $fasta_file`, stdout=io))
                elseif endswith(lf, ".bz2")
                    run(pipeline(`bzip2 -dc $fasta_file`, stdout=io))
                elseif endswith(lf, ".xz")
                    run(pipeline(`xz -dc $fasta_file`, stdout=io))
                elseif endswith(lf, ".zip")
                    run(pipeline(`unzip -p $fasta_file`, stdout=io))
                else
                    # Fallback: try gzip -dc (may fail if not actually gzipped)
                    run(pipeline(`gzip -dc $fasta_file`, stdout=io))
                end
            end

            input_fasta = temp_fasta
        end

        # Run padloc against the chosen input_fasta; ensure temporary files get cleaned up
        try
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --fna $(input_fasta) --outdir $(outdir) --cpu $(threads)`)
        finally
            if is_compressed && temp_dir != ""
                try
                    # Remove the temporary directory and its contents
                    rm(temp_dir; force=true, recursive=true)
                    @info "Removed temporary directory $(temp_dir)"
                catch err
                    @warn "Failed to remove temporary files at $(temp_dir): $err"
                end
            end
        end
    else
        @info "$(padloc_outfile) already present"
    end

    padloc_faa = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.faa"))
    padloc_gff = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.gff"))
    padloc_domtblout = joinpath(outdir, replace(basename(fasta_file), Mycelia.FASTA_REGEX => ".domtblout"))
    if !isfile(padloc_outfile)
        padloc_outfile = missing
    end
    return (csv = padloc_outfile, faa = padloc_faa, gff = padloc_gff, domtblout = padloc_domtblout)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Multi-Locus Sequence Typing (MLST) analysis on a genome assembly.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing the genome assembly

# Returns
- Path to the output file containing MLST results (`<input>.mlst.out`)

# Details
Uses the `mlst` tool from PubMLST to identify sequence types by comparing allelic 
profiles of housekeeping genes against curated MLST schemes.

# Dependencies
- Requires Bioconda and the `mlst` package
- Automatically sets up conda environment if not present
"""
function run_mlst(fasta_file)
    Mycelia.add_bioconda_env("mlst")    
    mlst_outfile = "$(fasta_file).mlst.out"
    @show mlst_outfile
    if !isfile(mlst_outfile)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mlst mlst $(fasta_file)`, mlst_outfile))
    else
        @info "$(mlst_outfile) already present"
    end
    return mlst_outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run ECTyper for serotyping E. coli genome assemblies.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing assembled genome(s)

# Returns
- `String`: Path to output directory containing ECTyper results
"""
function run_ectyper(fasta_file)
    Mycelia.add_bioconda_env("ectyper")
    outdir = fasta_file * "_ectyper"
    if !isdir(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ectyper ectyper -i $(fasta_file) -o $(outdir)`)
    else
        @info "$(outdir) already present"
    end
    return outdir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run TransTermHP to predict rho-independent transcription terminators in DNA sequences.

# Arguments
- `fasta`: Path to input FASTA file containing DNA sequences
- `gff`: Optional path to GFF annotation file. If provided, improves prediction accuracy

# Returns
- `String`: Path to output file containing TransTermHP predictions

# Details
- Uses Conda environment 'transtermhp' for execution
- Automatically generates coordinate file from FASTA or GFF input
- Removes temporary coordinate file after completion
- Requires Mycelia's Conda setup
"""
function run_transterm(;fasta, gff="")
    # note in my one test with phage genomes, calling without gff yeilds more hits but average confidence is a bit lower
    if isempty(gff)
        coordinates = generate_transterm_coordinates_from_fasta(fasta)
    else
        coordinates = generate_transterm_coordinates_from_gff(gff)
    end
    transterm_calls_file = replace(coordinates, ".coords" => ".transterm.txt")
    Mycelia.add_bioconda_env("transtermhp")

    conda_base = dirname(dirname(Mycelia.CONDA_RUNNER))
    dat_file = "$(conda_base)/envs/transtermhp/data/expterm.dat"
    # @show isfile(dat_file)
    # @assert = first(readlines(`find $(conda_base) -name "expterm.dat"`))
    # @show dat_file
    @assert isfile(dat_file)
    # conda_base = chomp(read(`$(Mycelia.CONDA_RUNNER) info --base`, String))
    # dat_file = "$(conda_base)/envs/transtermhp/data/expterm.dat"
    # @assert isfile(dat_file)
    run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n transtermhp transterm -p $(dat_file) $(fasta) $(coordinates)`, transterm_calls_file))    
    rm(coordinates)
    return transterm_calls_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate minimal coordinate files required for TransTermHP analysis from FASTA sequences.

Creates artificial gene annotations at sequence boundaries to enable TransTermHP to run
without real gene annotations. For each sequence in the FASTA file, generates two
single-base-pair "genes" at positions 1-2 and (L-1)-L, where L is sequence length.

# Arguments
- `fasta`: Path to input FASTA file containing sequences to analyze

# Returns
- Path to generated coordinate file (original path with ".coords" extension)

# Format
Generated coordinate file follows TransTermHP format:
gene_id start stop chromosome

where chromosome matches FASTA sequence identifiers.

See also: [`run_transterm`](@ref)
"""
function generate_transterm_coordinates_from_fasta(fasta)
    # 10. USING TRANSTERM WITHOUT GENOME ANNOTATIONS

    # TransTermHP uses known gene information for only 3 things: (1) tagging the
    # putative terminators as either "inside genes" or "intergenic," (2) choosing the
    # background GC-content percentage to compute the scores, because genes often
    # have different GC content than the intergenic regions, and (3) producing
    # slightly more readable output. Items (1) and (3) are not really necessary, and
    # (2) has no effect if your genes have about the same GC-content as your
    # intergenic regions.

    # Unfortunately, TransTermHP doesn't yet have a simple option to run without an
    # annotation file (either .ptt or .coords), and requires at least 2 genes to be
    # present. The solution is to create fake, small genes that flank each
    # chromosome. To do this, make a fake.coords file that contains only these two
    # lines:

    # 	fakegene1	1 2	chome_id
    # 	fakegene2	L-1 L	chrom_id

    # where L is the length of the input sequence and L-1 is 1 less than the length
    # of the input sequence. "chrom_id" should be the word directly following the ">"
    # in the .fasta file containing your sequence. (If, for example, your .fasta file
    # began with ">seq1", then chrom_id = seq1).

    # This creates a "fake" annotation with two 1-base-long genes flanking the
    # sequence in a tail-to-tail arrangement: --> <--. TransTermHP can then be run
    # with:

    # 	transterm -p expterm.dat sequence.fasta fake.coords

    # If the G/C content of your intergenic regions is about the same as your genes,
    # then this won't have too much of an effect on the scores terminators receive.
    # On the other hand, this use of TransTermHP hasn't been tested much at all, so
    # it's hard to vouch for its accuracy.

    coords_table = DataFrames.DataFrame(
        gene_id = String[],
        start = Int[],
        stop = Int[],
        chromosome = String[]
    )
    for record in Mycelia.open_fastx(fasta)
        row = (
            gene_id = FASTX.identifier(record) * "_start",
            start = 1,
            stop = 2,
            chromosome = FASTX.identifier(record)
        )
        push!(coords_table, row)

        row = (
            gene_id = FASTX.identifier(record) * "_stop",
            start = length(FASTX.sequence(record))-1,
            stop = length(FASTX.sequence(record)),
            chromosome = FASTX.identifier(record)
        )
        push!(coords_table, row)
    end
    transterm_coordinates_file = fasta * ".coords"
    uCSV.write(transterm_coordinates_file, data = collect(DataFrames.eachcol(coords_table)), delim="  ")
    return transterm_coordinates_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GFF file to a coordinates file compatible with TransTermHP format.

# Arguments
- `gff_file::String`: Path to input GFF file

# Processing
- Converts 1-based to 0-based coordinates
- Extracts gene IDs from the attributes field
- Retains columns: gene_id, start, end, seqid

# Returns
- Path to the generated coordinates file (original filename with '.coords' suffix)

# Output Format
Space-delimited file with columns: gene_id, start, end, seqid
"""
function generate_transterm_coordinates_from_gff(gff_file)
    raw_gff = Mycelia.read_gff(gff_file)
    # switch start to be zero index by subtracting one
    raw_gff[!, "start"] = raw_gff[!, "start"] .- 1
    raw_gff[!, "gene_id"] = last.(split.(first.(split.(raw_gff[!, "attributes"], ";")), '='))
    raw_gff = raw_gff[!, ["gene_id", "start", "end", "#seqid"]]
    transterm_coordinates_file = gff_file * ".coords"
    uCSV.write(transterm_coordinates_file, data = collect(DataFrames.eachcol(raw_gff)), delim="  ")
    return transterm_coordinates_file
end

# function run_prokka(ID, OUT_DIR, normalized_fasta_file)
#     prokka_dir="$(OUT_DIR)/prokka"
#     if !isdir(prokka_dir)
#         mkdir(prokka_dir)
#     end
#     prokka_cmd = `prokka --force --cpus 1 --outdir $(prokka_dir) --prefix $(ID) $(normalized_fasta_file)`
#     run(pipeline(prokka_cmd, stdout="$(prokka_dir)/prokka.out"))
#     return prokka_dir
# end

# function run_mlst(ID, OUT_DIR, normalized_fasta_file)
#     mlst_dir="$(OUT_DIR)/mlst"
#     if !isdir(mlst_dir)
#         mkdir(mlst_dir)
#     end
#     p = pipeline(
#             `mlst $(normalized_fasta_file)`,
#             stdout="$(mlst_dir)/$(ID).mlst.out")
#     run(p)
#     return mlst_dir
# end

# function run_phispy(ID, OUT_DIR, prokka_dir)
#     # 1 	prophage_coordinates.tsv
#     # 2 	GenBank format output
#     # 4 	prophage and bacterial sequences
#     # 8 	prophage_information.tsv
#     # 16 	prophage.tsv
#     # 32 	GFF3 format
#     # 64 	prophage.tbl
#     # 128 	test data used in the random forest
#     # 255   for all of them

#     phispy_dir="$(OUT_DIR)/phispy"
#     if !isdir(phispy_dir)
#         mkdir(phispy_dir)
#     end
#     if isempty(readdir(phispy_dir))
#         phisphy_cmd = `PhiSpy.py $(prokka_dir)/$(ID).gbk --output_dir $(phispy_dir) --file_prefix $(ID)-phispy --output_choice 255`
#         try
#             run(pipeline(phisphy_cmd, stdout="$(phispy_dir)/phisphy.out"))
#         catch
#             if isfile("$(phispy_dir)/$(ID)-phispy_prophage.gff3")
#                 @warn "phispy errored out after prophage gff3 was written"
#             else
#                 @error "phispy prophage gff3 not written"
#             end
#         end
#     end
#     return phispy_dir
# end

# function run_trnascan(ID, out_dir, normalized_fasta_file)
#     trnascan_dir = "$(out_dir)/trnascan"
#     # trnascan doesn't like to overwrite existing things
#     if !isdir(trnascan_dir)
#         mkdir(trnascan_dir)
#     end
#     if isempty(readdir(trnascan_dir))

#         #     -B for using Bacterial
#         trnascan_cmd = 
#         `tRNAscan-SE 
#             -B 
#             --output $(trnascan_dir)/$(ID).trnascan.out 
#             --bed $(trnascan_dir)/$(ID).trnascan.bed 
#             --fasta $(trnascan_dir)/$(ID).trnascan.fasta 
#             --struct $(trnascan_dir)/$(ID).trnascan.struct
#             --stats $(trnascan_dir)/$(ID).trnascan.stats 
#             --log $(trnascan_dir)/$(ID).trnascan.log
#             $(normalized_fasta_file)`
#         run(pipeline(trnascan_cmd, stdout="$(trnascan_dir)/trnascan.out", stderr="$(trnascan_dir)/trnascan.out"))
#     end
#     return trnascan_dir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run tRNAscan-SE to identify and annotate transfer RNA genes in the provided sequence file.

# Arguments
- `fna_file::String`: Path to input FASTA nucleotide file
- `outdir::String`: Output directory path (default: input_file_path + "_trnascan")

# Returns
- `String`: Path to the output directory containing tRNAscan-SE results

# Output Files
Creates the following files in `outdir`:
- `*.trnascan.out`: Main output with tRNA predictions
- `*.trnascan.bed`: BED format coordinates
- `*.trnascan.fasta`: FASTA sequences of predicted tRNAs
- `*.trnascan.struct`: Secondary structure predictions
- `*.trnascan.stats`: Summary statistics
- `*.trnascan.log`: Program execution log

# Notes
- Uses the general tRNA model (-G flag) suitable for all domains of life
- Automatically sets up tRNAscan-SE via Bioconda
- Skips processing if output directory contains files
"""
function run_trnascan(;fna_file, outdir=fna_file * "_trnascan")
    Mycelia.add_bioconda_env("trnascan-se")
    # trnascan doesn't like to overwrite existing things
    if !isdir(outdir)
        mkdir(outdir)
    else
        @info "$(outdir) already exists"
    end
    ID = basename(fna_file)
    if isempty(readdir(outdir))
        # -G : use general tRNA model (cytoslic tRNAs from all 3 domains included)
        #     -B for using Bacterial
        trnascan_cmd =
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n trnascan-se tRNAscan-SE
            -G
            --output $(outdir)/$(ID).trnascan.out
            --bed $(outdir)/$(ID).trnascan.bed
            --fasta $(outdir)/$(ID).trnascan.fasta
            --struct $(outdir)/$(ID).trnascan.struct
            --stats $(outdir)/$(ID).trnascan.stats
            --log $(outdir)/$(ID).trnascan.log
            --prefix
            --progress
            $(fna_file)`
        # run(pipeline(trnascan_cmd, stdout="$(outdir)/$(ID).trnascan.out", stderr="$(outdir)/$(ID).trnascan.out"))
        run(trnascan_cmd)
    end
    return outdir
end

# function run_counterselection_spacer_detection(strain, out_dir, normalized_fasta_file)
#     counter_selection_dir = "$(out_dir)/counter-selection"
#     if !isdir(counter_selection_dir)
#         mkdir(counter_selection_dir)
#     end

#     if isempty(readdir(counter_selection_dir))
#         regex = BioSequences.biore"TTT[CG][ACGT]{25}"dna
#         k = 29
#         KMER_TYPE = BioSequences.BigDNAMer{k}

#         spacer_table = DataFrames.DataFrame(
#             ID = [],
#             contig = [],
#             PAM_and_spacer = [],
#             spacer = [],
#             strand = [],
#             start = [],
#             stop = [],
# #             free_energy = [],
# #             visualization_url = []
#         )

#         ProgressMeter.@showprogress for record in collect(FASTX.FASTA.Reader(open(normalized_fasta_file)))
#             for (i, kmer, reverse_complement_kmer) in BioSequences.each(KMER_TYPE, FASTX.FASTA.sequence(record))
#                 strand = missing
#                 if occursin(regex, kmer)
#                     strand = "+"
#                 elseif occursin(regex, reverse_complement_kmer)
#                     strand = "-"
#                     kmer = reverse_complement_kmer
#                 end
#                 if !ismissing(strand)
#                     spacer = BioSequences.DNAMer(kmer[i] for i in 5:length(kmer)) 
# #                     RNAfold_output = read(pipeline(`echo "$(string(spacer))"`, `RNAfold --noLP`), String)

# #                     rna_sequence, structure, free_energy = match(r"([ACGU]{25})\n([.()]{25})\s\(\s*(.*?)\)", RNAfold_output).captures
# #                     url = "http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence=$(rna_sequence)&structure=$(structure)"

#                     kmer_range = i:i+k-1

#                     t = DataFrames.DataFrame(
#                         ID = ID,
#                         contig = FASTX.FASTA.identifier(record),
#                         PAM_and_spacer = kmer,
#                         spacer = spacer,
#                         strand = strand,
#                         start = i,
#                         stop = i+k-1,
# #                         free_energy = free_energy,
# #                         visualization_url = url
#                     )
#                     spacer_table = vcat(spacer_table, t)
#                 end
#             end
#         end


#         uCSV.write(
#             "$(counter_selection_dir)/$(ID)-cpf1-spacers.tsv",
#             delim='\t',
#             data = collect(DataFrames.eachcol(spacer_table)),
#             header = DataFrames.names(spacer_table)
#         )
#         if isfile("$(out_dir)/rna.ps")
#             rm("$(out_dir)/rna.ps")
#         end
#     end
#     return counter_selection_dir
# end



# function run_amrfinderplus(ID, out_dir, protein_fasta)
#     amrfinderplus_dir = "$(out_dir)/amrfinderplus"
#     if !isdir(amrfinderplus_dir)
#         mkdir(amrfinderplus_dir)
#     end

#     if isempty(readdir(amrfinderplus_dir))
# #         run(`amrfinder -u`)
        
#         # because the pipeline is set up with some hacky CONDA path manipulation, 
#         # explictly setting amrfinder directory path to the location in the docker host
#         amrfinder_db_path = get(ENV, "AMRFINDER_DB", "none")

#         if amrfinder_db_path != "none"
#             cmd = 
#             `amrfinder
#             -p $(protein_fasta)
#             --plus
#             --output $(amrfinderplus_dir)/$(ID).amrfinderplus.tsv
#             -d $(amrfinder_db_path)
#             `
#         else
#             cmd = 
#             `amrfinder
#             -p $(protein_fasta)
#             --plus
#             --output $(amrfinderplus_dir)/$(ID).amrfinderplus.tsv
#             `
#         end
            

#         p = pipeline(cmd, 
#                 stdout="$(amrfinderplus_dir)/$(ID).amrfinderplus.out",
#                 stderr="$(amrfinderplus_dir)/$(ID).amrfinderplus.err")
#         run(p)
#     end
#     return amrfinderplus_dir
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Prodigal gene prediction software on input FASTA file to identify protein-coding genes
in metagenomes or single genomes.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing genomic sequences
- `out_dir::String=dirname(fasta_file)`: Directory for output files. Defaults to input file's directory

# Returns
Named tuple containing paths to all output files:
- `fasta_file`: Input FASTA file path
- `out_dir`: Output directory path  
- `gff`: Path to GFF format gene predictions
- `gene_scores`: Path to all potential genes and their scores
- `fna`: Path to nucleotide sequences of predicted genes
- `faa`: Path to protein translations of predicted genes
- `std_out`: Path to captured stdout
- `std_err`: Path to captured stderr
"""
function run_prodigal(;fasta_file, out_dir=dirname(fasta_file))
    
    # if isempty(out_dir)
    #     prodigal_dir = mkpath("$(fasta_file)_prodigal")
    # else
    #     prodigal_dir = mkpath(out_dir)
    # end

    # $ prodigal
    # -------------------------------------
    # PRODIGAL v2.6.3 [February, 2016]         
    # Univ of Tenn / Oak Ridge National Lab
    # Doug Hyatt, Loren Hauser, et al.     
    # -------------------------------------

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    gff = "$(out_dir)/$(basename(fasta_file)).prodigal.gff"
    faa = "$(out_dir)/$(basename(fasta_file)).prodigal.faa"
    fna = "$(out_dir)/$(basename(fasta_file)).prodigal.fna"
    gene_scores = "$(out_dir)/$(basename(fasta_file)).prodigal.all_potential_gene_scores.txt"
    std_out = "$(out_dir)/$(basename(fasta_file)).prodigal.out"
    std_err = "$(out_dir)/$(basename(fasta_file)).prodigal.err"
    
    # I usually delete the rest, so don't reprocess if outputs of interest are present
    if (!isfile(gff) && !isfile(faa))
        add_bioconda_env("prodigal")
        cmd = 
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n prodigal prodigal
        -f gff
        -m
        -p meta
        -o $(gff)
        -i $(fasta_file)
        -a $(faa)
        -d $(fna)
        -s $(gene_scores)
        `
        p = pipeline(cmd, stdout=std_out, stderr=std_err)
        run(p)
    end
    return (;fasta_file, out_dir, gff, gene_scores, fna, faa, std_out, std_err)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Pyrodigal gene prediction on a FASTA file using the meta procedure optimized for metagenomic sequences.

Pyrodigal is a reimplementation of the Prodigal gene finder, which identifies protein-coding sequences in bacterial and archaeal genomes.

# Arguments
- `fasta_file::String`: Path to input FASTA file containing genomic sequences
- `out_dir::String`: Output directory path (default: input filename + "_pyrodigal")

# Returns
Named tuple containing:
- `fasta_file`: Input FASTA file path
- `out_dir`: Output directory path
- `gff`: Path to GFF output file with gene predictions
- `faa`: Path to FASTA file with predicted protein sequences 
- `fna`: Path to FASTA file with nucleotide sequences

# Notes
- Uses metagenomic mode (`-p meta`) optimized for mixed communities
- Masks runs of N nucleotides (`-m` flag)
- Minimum gene length set to 33bp
- Maximum overlap between genes set to 31bp
- Requires Pyrodigal to be available in a Conda environment
- Skips processing if output files already exist
"""
function run_pyrodigal(;fasta_file, out_dir=fasta_file * "_pyrodigal")
    # https://pyrodigal.readthedocs.io/en/stable/guide/cli.html#command-line-interface

    # -a trans_file         Write protein translations to the selected file.
    # -c                    Closed ends. Do not allow genes to run off edges.
    # -d nuc_file           Write nucleotide sequences of genes to the selected file.
    # -f output_type        Select output format.
    # -g tr_table           Specify a translation table to use.
    # -i input_file         Specify FASTA input file.
    # -m                    Treat runs of N as masked sequence and don't build genes across them.
    # -n                    Bypass Shine-Dalgarno trainer and force a full motif scan.
    # -o output_file        Specify output file.
    # -p mode               Select procedure.
    # -s start_file         Write all potential genes (with scores) to the selected file.
    # -t training_file      Write a training file (if none exists); otherwise, read and use the specified training file.
    # -j jobs, --jobs jobs           The number of threads to use if input contains multiple sequences.
    # --min-gene MIN_GENE            The minimum gene length.
    # --min-edge-gene MIN_EDGE_GENE  The minimum edge gene length.
    # --max-overlap MAX_OVERLAP      The maximum number of nucleotides that can overlap between two genes on the same strand.
    #                             This must be lower or equal to the minimum gene length.
    # --no-stop-codon                Disables translation of stop codons into star characters (*) for complete genes.
    # --pool {thread,process}        The sort of pool to use to process genomes in parallel. Processes may be faster than
    #                             threads on some machines, refer to documentation. (default: thread)

    # Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
    #                  [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
    #                  [-p mode] [-q] [-s start_file] [-t training_file] [-v]

    #          -a:  Write protein translations to the selected file.
    #          -c:  Closed ends.  Do not allow genes to run off edges.
    #          -d:  Write nucleotide sequences of genes to the selected file.
    #          -f:  Select output format (gbk, gff, or sco).  Default is gbk.
    #          -g:  Specify a translation table to use (default 11).
    #          -h:  Print help menu and exit.
    #          -i:  Specify FASTA/Genbank input file (default reads from stdin).
    #          -m:  Treat runs of N as masked sequence; don't build genes across them.
    #          -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
    #          -o:  Specify output file (default writes to stdout).
    #          -p:  Select procedure (single or meta).  Default is single.
    #          -q:  Run quietly (suppress normal stderr output).
    #          -s:  Write all potential genes (with scores) to the selected file.
    #          -t:  Write a training file (if none exists); otherwise, read and use
    #               the specified training file.
    #          -v:  Print version number and exit.
    gff = "$(out_dir)/$(basename(fasta_file)).pyrodigal.gff"
    faa = "$(out_dir)/$(basename(fasta_file)).pyrodigal.faa"
    fna = "$(out_dir)/$(basename(fasta_file)).pyrodigal.fna"
    gene_scores = "$(out_dir)/$(basename(fasta_file)).pyrodigal.gene_scores.txt"
    std_out = "$(out_dir)/$(basename(fasta_file)).pyrodigal.out"
    std_err = "$(out_dir)/$(basename(fasta_file)).pyrodigal.err"
    mkpath(out_dir)
    
    # I usually delete the rest, so don't reprocess if outputs of interest are present
    # `max_overlap` must be lower than `min_gene`
    if (!isfile(gff) || !isfile(faa))
        Mycelia.add_bioconda_env("pyrodigal")
        cmd = 
        `$(Mycelia.CONDA_RUNNER) run --no-capture-output -n pyrodigal pyrodigal
        -f gff
        -m
        -p meta
        -o $(gff)
        -i $(fasta_file)
        -a $(faa)
        -d $(fna)
        -s $(gene_scores)
        --min-gene 33
        --max-overlap 31
        `
        p = pipeline(cmd, stdout=std_out, stderr=std_err)
        run(p)
        if fasta_file != "$(out_dir)/$(basename(fasta_file))"
            cp(fasta_file, "$(out_dir)/$(basename(fasta_file))", force=true)
        end
        if isfile(std_out) && (filesize(std_out) == 0)
            rm(std_out)
        end
        if isfile(std_err) && (filesize(std_err) == 0)
            rm(std_err)
        end
    end
    # return (;fasta_file, out_dir, gff, gene_scores, fna, faa, std_out, std_err)
    return (;fasta_file, out_dir, gff, faa, fna)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse TransTerm terminator prediction output into a structured DataFrame.

Takes a TransTerm output file path and returns a DataFrame containing parsed terminator predictions.
Each row represents one predicted terminator with the following columns:

- `chromosome`: Identifier of the sequence being analyzed
- `term_id`: Unique terminator identifier (e.g. "TERM 19")
- `start`: Start position of the terminator
- `stop`: End position of the terminator
- `strand`: Strand orientation ("+" or "-")
- `location`: Context type, where:
    * G/g = in gene interior (≥50bp from ends)
    * F/f = between two +strand genes
    * R/r = between two -strand genes
    * T = between ends of +strand and -strand genes
    * H = between starts of +strand and -strand genes
    * N = none of the above
    Lowercase indicates opposite strand from region
- `confidence`: Overall confidence score (0-100)
- `hairpin_score`: Hairpin structure score
- `tail_score`: Tail sequence score  
- `notes`: Additional annotations (e.g. "bidir")

# Arguments
- `transterm_output::AbstractString`: Path to TransTerm output file

# Returns
- `DataFrame`: Parsed terminator predictions with columns as described above

See TransTerm HP documentation for details on scoring and location codes.
"""
function parse_transterm_output(transterm_output)
    
   #     3. FORMAT OF THE TRANSTERM OUTPUT

#     The organism's genes are listed sorted by their end coordinate and terminators
#     are output between them. A terminator entry looks like this:

#         TERM 19  15310 - 15327  -      F     99      -12.7 -4.0 |bidir
#         (name)   (start - end)  (sense)(loc) (conf) (hp) (tail) (notes)

#     where 'conf' is the overall confidence score, 'hp' is the hairpin score, and
#     'tail' is the tail score. 'Conf' (which ranges from 0 to 100) is what you
#     probably want to use to assess the quality of a terminator. Higher is better.
#     The confidence, hp score, and tail scores are described in the paper cited
#     above.  'Loc' gives type of region the terminator is in:

#         'G' = in the interior of a gene (at least 50bp from an end),
#         'F' = between two +strand genes,
#         'R' = between two -strand genes,
#         'T' = between the ends of a +strand gene and a -strand gene,
#         'H' = between the starts of a +strand gene and a -strand gene,
#         'N' = none of the above (for the start and end of the DNA)

#     Because of how overlapping genes are handled, these designations are not
#     exclusive. 'G', 'F', or 'R' can also be given in lowercase, indicating that
#     the terminator is on the opposite strand as the region.  Unless the
#     --all-context option is given, only candidate terminators that appear to be in
#     an appropriate genome context (e.g. T, F, R) are output. 

#     Following the TERM line is the sequence of the hairpin and the 5' and 3'
#     tails, always written 5' to 3'.
    
    transterm_table = DataFrames.DataFrame()
    chromosome = ""
    for line in Iterators.filter(x -> occursin(r"^\s*(SEQUENCE|TERM)", x), eachline(transterm_output))
        line = strip(line)
        if occursin(r"^SEQUENCE", line)
            chromosome = split(line)[2]
        else
            transterm_regex = r"(TERM \d+)\s+(\d+) - (\d+)\s+(\S)\s+(\w+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\|(.*)"
            term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes = match(transterm_regex, line).captures
            notes = strip(notes)
            row = (;chromosome, term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes)
            push!(transterm_table, row)
        end
    end
    return transterm_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert TransTerm terminator predictions output to GFF3 format.

Parses TransTerm output and generates a standardized GFF3 file with the following transformations:
- Sets source field to "transterm"
- Sets feature type to "terminator"  
- Converts terminator IDs to GFF attributes
- Renames fields to match GFF3 spec

# Arguments
- `transterm_output::String`: Path to the TransTerm output file

# Returns
- `String`: Path to the generated GFF3 file (original filename with .gff extension)
"""
function transterm_output_to_gff(transterm_output)
    transterm_table = parse_transterm_output(transterm_output)
    transterm_table[!, "source"] .= "transterm"
    transterm_table[!, "type"] .= "terminator"
    transterm_table[!, "phase"] .= "."
    transterm_table[!, "attributes"] = map(x -> "label=" * replace(x, " " => "_"), transterm_table[!, "term_id"])
    DataFrames.rename!(transterm_table,
        ["chromosome" => "#seqid",
            "stop" => "end",
            "confidence" => "score",
        ]
    )
    transterm_table = transterm_table[!, ["#seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    transterm_gff = write_gff(gff=transterm_table, outfile=transterm_output * ".gff")
    uCSV.write(transterm_gff, transterm_table, delim='\t')
    return transterm_gff
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a VirSorter score TSV file and return a DataFrame.

# Arguments
- `virsorter_score_tsv::String`: The file path to the VirSorter score TSV file.

# Returns
- `DataFrame`: A DataFrame containing the parsed data from the TSV file. If the file is empty, returns a DataFrame with the appropriate headers but no data.
"""
function parse_virsorter_score_tsv(virsorter_score_tsv)
    data, header = uCSV.read(virsorter_score_tsv, delim='\t', header=1)
    if length(data) == 0
        data = [[] for i in 1:length(header)]
    end
    return DataFrames.DataFrame(data, header)
end

"""
    count_predicted_genes(gff_file)

Count the number of predicted genes from a GFF file.

Parses a GFF/GTF file and counts the number of CDS (coding sequence) features,
which correspond to predicted genes.

# Arguments
- `gff_file`: Path to GFF/GTF file

# Returns
- Integer count of predicted genes (CDS features)

# See Also
- `run_pyrodigal`: For gene prediction that generates GFF files
- `parse_transterm_output`: For parsing other annotation tool outputs
"""
function count_predicted_genes(gff_file)
    gene_count = 0
    
    open(gff_file, "r") do f
        for line in eachline(f)
            # Skip comment lines
            if startswith(line, "#")
                continue
            end
            
            # Parse GFF line
            fields = split(line, "\t")
            if length(fields) >= 3 && fields[3] == "CDS"
                gene_count += 1
            end
        end
    end
    
    return gene_count
end