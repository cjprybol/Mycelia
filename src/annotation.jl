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

Perform comprehensive annotation of a FASTA file including gene prediction, protein homology search,
and terminator prediction.

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
        identifier::String = replace(basename(fasta), Mycelia.FASTA_REGEX => ""), # Assuming Mycelia.FASTA_REGEX is defined
        basedir::String = pwd(),        
        mmseqsdb::String = joinpath(homedir(), "workspace/mmseqs/UniRef50"),
        threads::Int = Sys.CPU_THREADS
    )
    
    outdir = joinpath(basedir, identifier)
    @assert outdir != fasta "Output directory cannot be the same as the input FASTA file path."
    
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

Perform all-vs-all sequence search using MMseqs2's easy-search command.

# Arguments
- `fasta::String`: Path to input FASTA file containing sequences to compare
- `output::String`: Output directory path (default: input filename + ".mmseqs_easy_search_pairwise")

# Returns
- `String`: Path to the output directory

# Details
Executes MMseqs2 with sensitive search parameters (7 sensitivity steps) and outputs results in 
tabular format with the following columns:
- query, qheader: Query sequence ID and header
- target, theader: Target sequence ID and header  
- pident: Percentage sequence identity
- fident: Fraction of identical matches
- nident: Number of identical matches
- alnlen: Alignment length
- mismatch: Number of mismatches
- gapopen: Number of gap openings
- qstart, qend, qlen: Query sequence coordinates and length
- tstart, tend, tlen: Target sequence coordinates and length
- evalue: Expected value
- bits: Bit score

Requires MMseqs2 to be available through Bioconda.
"""
function mmseqs_pairwise_search(;fasta, output=fasta*".mmseqs_easy_search_pairwise")
    Mycelia.add_bioconda_env("mmseqs2")
    mkpath(output)
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-search
        $(fasta)
        $(fasta)
        $(output)/$(basename(fasta)).mmseqs_pairwise_search.txt $(tempdir())
        --format-mode 4
        --format-output query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
        --start-sens 1 -s 7 --sens-steps 7 --sort-results 1 --remove-tmp-files 1 --search-type 3`)
    return output
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function mmseqs_easy_linclust(;fasta, output=fasta*".mmseqs_easy_linclust", tmp=mktempdir())
#     Mycelia.add_bioconda_env("mmseqs2")   
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createdb $(fasta) $(fasta)_DB`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createindex --search-type 3 $(fasta)_DB $(tempdir())`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-linclust $(fasta)_DB $(output) $(tmp)`)
#     run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createtsv $(fasta)_DB $(fasta)_DB $(output) $(output).tsv`)
#     return "$(output).tsv"
# end