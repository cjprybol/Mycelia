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
"""
function run_padloc(;fasta_file, outdir=dirname(abspath(fasta_file)), threads=Sys.CPU_THREADS)
    padloc_outfile = joinpath(outdir, replace(basename(fasta_file), ".fna" => "") * "_padloc.csv")
    if !isfile(padloc_outfile)
        setup_padloc()
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n padloc padloc --fna $(fasta_file) --outdir $(outdir) --cpu $(threads)`)
    else
        @info "$(padloc_outfile) already present"
    end
    padloc_faa = replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.faa")
    padloc_gff = replace(basename(fasta_file), Mycelia.FASTA_REGEX => "_prodigal.gff")
    padloc_domtblout = replace(basename(fasta_file), Mycelia.FASTA_REGEX => ".domtblout")
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
        (filesize(std_out) == 0) && rm(std_out)
        (filesize(std_err) == 0) && rm(std_err)
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