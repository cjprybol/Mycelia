"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform DIAMOND BLASTP search between query and reference protein FASTA files.

# Arguments
- `query_fasta::String`: Path to query protein FASTA file
- `reference_fasta::String`: Path to reference protein FASTA file
- `output_dir::String`: Output directory (defaults to query filename + "_diamond")
- `threads::Int`: Number of threads (defaults to system CPU count)
- `evalue::Float64`: E-value threshold (default: 1e-3)
- `block_size::Float64`: Block size in GB (default: auto-calculated from system memory)
- `sensitivity::String`: Sensitivity mode (default: "--iterate")

# Returns
- `String`: Path to the DIAMOND results file (.tsv format)

# Throws
- `AssertionError`: If input files don't exist or are invalid
- `SystemError`: If DIAMOND execution fails
"""
function run_diamond_search(;
    query_fasta::String,
    reference_fasta::String,
    output_dir::String = replace(basename(query_fasta), Mycelia.FASTA_REGEX => "") * "_diamond",
    threads::Int = get_default_threads(),
    evalue::Float64 = 1e-3,
    block_size::Float64 = floor(Sys.total_memory() / 1e9 / 8), # Auto-calculate from memory
    sensitivity::String = "--iterate"
)
    # Input validation and assertions
    @assert isfile(query_fasta) "Query FASTA file does not exist: $(query_fasta)"
    @assert isfile(reference_fasta) "Reference FASTA file does not exist: $(reference_fasta)"
    @assert threads > 0 "Thread count must be positive: $(threads)"
    @assert evalue > 0 "E-value must be positive: $(evalue)"
    @assert block_size > 0 "Block size must be positive: $(block_size)"
    
    # Validate FASTA files have content
    @assert filesize(query_fasta) > 0 "Query FASTA file is empty: $(query_fasta)"
    @assert filesize(reference_fasta) > 0 "Reference FASTA file is empty: $(reference_fasta)"
    
    # Setup output directory and files
    mkpath(output_dir)
    diamond_db = joinpath(output_dir, "diamond_db.dmnd")
    results_file = joinpath(output_dir, replace(basename(query_fasta), Mycelia.FASTA_REGEX => "") * "__" * replace(basename(reference_fasta), Mycelia.FASTA_REGEX => "") * "_diamond_results.tsv")

    if isfile(results_file) && (filesize(results_file) > 0)
        return results_file
    end
    
    # Ensure DIAMOND environment exists
    Mycelia.add_bioconda_env("diamond")

    # Map for column header expansion
    # outfmt_fields = [
    #     "qseqid", "qtitle", "qlen", "sseqid", "sallseqid", "stitle",
    #     "salltitles", "slen", "qstart", "qend", "sstart", "send",
    #     "evalue", "bitscore", "length", "pident", "nident", "mismatch", "gapopen"
    # ]
    outfmt_headers = [
        "Query Seq - id", "Query title", "Query sequence length", "Subject Seq - id", "All subject Seq - id(s)",
        "Subject Title", "All Subject Title(s)", "Subject sequence length", "Start of alignment in query",
        "End of alignment in query", "Start of alignment in subject", "End of alignment in subject",
        "Expect value", "Bit score", "Alignment length", "Percentage of identical matches",
        "Number of identical matches", "Number of mismatches", "Number of gap openings"
    ]

    try
        # Create DIAMOND database
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n diamond diamond makedb --in $(reference_fasta) --db $(diamond_db)`)
        
        @assert isfile(diamond_db) "DIAMOND database creation failed: $(diamond_db)"
        
        # Run DIAMOND search
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n diamond diamond blastp --query $(query_fasta) --db $(diamond_db) --out $(results_file) --evalue $(evalue) --threads $(threads) --block-size $(block_size) $(sensitivity) --outfmt 6 qseqid qtitle qlen sseqid sallseqid stitle salltitles slen qstart qend sstart send evalue bitscore length pident nident mismatch gapopen`)
        
        @assert isfile(results_file) "DIAMOND results file was not created: $(results_file)"
        @assert filesize(results_file) > 0 "DIAMOND results file is empty"

        # Insert the header row
        results_tmp = results_file * ".tmp"
        open(results_tmp, "w") do out_io
            # Write header
            println(out_io, join(outfmt_headers, '\t'))
            # Write results
            open(results_file, "r") do in_io
                for line in eachline(in_io)
                    println(out_io, line)
                end
            end
        end
        mv(results_tmp, results_file; force=true)
        
        return results_file
        
    catch e
        @error "DIAMOND execution failed" exception=e
        rethrow(e)
    finally
        # Cleanup database file
        rm(diamond_db, force=true)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform BLASTP search between query and reference protein FASTA files.

# Arguments
- `query_fasta::String`: Path to query protein FASTA file
- `reference_fasta::String`: Path to reference protein FASTA file
- `output_dir::String`: Output directory (defaults to query filename + "_blastp")
- `threads::Int`: Number of threads (defaults to system CPU count)
- `evalue::Float64`: E-value threshold (default: 1e-3)
- `max_target_seqs::Int`: Maximum target sequences (default: 500)

# Returns
- `String`: Path to the BLASTP results file (.tsv format)

# Throws
- `AssertionError`: If input files don't exist or are invalid
- `SystemError`: If BLAST execution fails
"""
function run_blastp_search(;
    query_fasta::String,
    reference_fasta::String,
    output_dir::String = replace(basename(query_fasta), Mycelia.FASTA_REGEX => "") * "_blastp",
    threads::Int = get_default_threads(),
    evalue::Float64 = 1e-3,
    max_target_seqs::Int = 500
)
    # Input validation and assertions
    @assert isfile(query_fasta) "Query FASTA file does not exist: $(query_fasta)"
    @assert isfile(reference_fasta) "Reference FASTA file does not exist: $(reference_fasta)"
    @assert threads > 0 "Thread count must be positive: $(threads)"
    @assert evalue > 0 "E-value must be positive: $(evalue)"
    @assert max_target_seqs > 0 "Max target sequences must be positive: $(max_target_seqs)"
    
    # Validate FASTA files have content
    @assert filesize(query_fasta) > 0 "Query FASTA file is empty: $(query_fasta)"
    @assert filesize(reference_fasta) > 0 "Reference FASTA file is empty: $(reference_fasta)"
    
    # Setup output directory and files
    mkpath(output_dir)
    blast_db = joinpath(output_dir, "blast_db")
    # results_file = joinpath(output_dir, "blastp_results.tsv")
    results_file = joinpath(output_dir, replace(basename(query_fasta), Mycelia.FASTA_REGEX => "") * "__" * replace(basename(reference_fasta), Mycelia.FASTA_REGEX => "") * "_blastp_results.tsv")

    if isfile(results_file) && (filesize(results_file) > 0)
        return results_file
    end
    
    # Ensure BLAST environment exists
    Mycelia.add_bioconda_env("blast")
    
    try
        # Create BLAST database
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast makeblastdb -in $(reference_fasta) -dbtype prot -out $(blast_db)`)
        
        # Verify database was created
        @assert isfile("$(blast_db).phr") "BLAST database creation failed"
        
        # Run BLASTP search
        outfmt = "7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid"
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastp -query $(query_fasta) -db $(blast_db) -out $(results_file) -evalue $(evalue) -max_target_seqs $(max_target_seqs) -num_threads $(threads) -outfmt "$(outfmt)"`)
        
        @assert isfile(results_file) "BLASTP results file was not created: $(results_file)"
        
        return results_file
        
    catch e
        @error "BLASTP execution failed" exception=e
        rethrow(e)
    finally
        # Cleanup database files
        for ext in [".phr", ".pin", ".psq"]
            rm("$(blast_db)$(ext)", force=true)
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Perform MMseqs2 easy-search between query and reference FASTA files.

# Arguments
- `query_fasta::String`: Path to query FASTA file
- `reference_fasta::String`: Path to reference FASTA file  
- `output_dir::String`: Output directory (defaults to query filename + "_mmseqs")
- `threads::Int`: Number of threads (defaults to system CPU count)
- `evalue::Float64`: E-value threshold (default: 1e-3)
- `sensitivity::Float64`: Sensitivity parameter (default: 4.0)

# Returns
- `String`: Path to the MMseqs2 results file (.tsv format)

# Throws
- `AssertionError`: If input files don't exist or are invalid
- `SystemError`: If MMseqs2 execution fails
"""
function run_mmseqs_search(;
    query_fasta::String,
    reference_fasta::String,
    output_dir::String = replace(basename(query_fasta), Mycelia.FASTA_REGEX => "") * "_mmseqs",
    threads::Int = get_default_threads(),
    evalue::Float64 = 1e-3,
    sensitivity::Float64 = 4.0
)
    # Input validation and assertions
    @assert isfile(query_fasta) "Query FASTA file does not exist: $(query_fasta)"
    @assert isfile(reference_fasta) "Reference FASTA file does not exist: $(reference_fasta)"
    @assert threads > 0 "Thread count must be positive: $(threads)"
    @assert evalue > 0 "E-value must be positive: $(evalue)"
    @assert sensitivity > 0 "Sensitivity must be positive: $(sensitivity)"
    
    # Validate FASTA files have content
    @assert filesize(query_fasta) > 0 "Query FASTA file is empty: $(query_fasta)"
    @assert filesize(reference_fasta) > 0 "Reference FASTA file is empty: $(reference_fasta)"
    
    # Setup output directory and files
    mkpath(output_dir)
    query_db = joinpath(output_dir, "query_db")
    ref_db = joinpath(output_dir, "ref_db") 
    results_file = joinpath(output_dir, replace(basename(query_fasta), Mycelia.FASTA_REGEX => "") * "__" * replace(basename(reference_fasta), Mycelia.FASTA_REGEX => "") * "_mmseqs-easy-search.tsv")
    tmp_dir = joinpath(output_dir, "tmp")
    mkpath(tmp_dir)

    if isfile(results_file) && (filesize(results_file) > 0)
        return results_file
    end
    
    # Ensure MMseqs2 environment exists
    Mycelia.add_bioconda_env("mmseqs2")

    try
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs easy-search
            $(query_fasta)
            $(reference_fasta)
            $(results_file) $(tempdir())
            --format-mode 4
            --format-output query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
            --start-sens 1 -s 7 --sens-steps 7 --sort-results 1 --remove-tmp-files 1 --search-type 3`)
    
        # # Create databases
        # run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createdb $(query_fasta) $(query_db)`)
        # run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs createdb $(reference_fasta) $(ref_db)`)
        
        # # Run search
        # search_db = joinpath(output_dir, "search_results")
        # run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs search $(query_db) $(ref_db) $(search_db) $(tmp_dir) --threads $(threads) -e $(evalue) -s $(sensitivity)`)
        
        # # Convert to readable format
        # run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mmseqs2 mmseqs convertalis $(query_db) $(ref_db) $(search_db) $(results_file) --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"`)
        
        @assert isfile(results_file) "MMseqs2 results file was not created: $(results_file)"
        @assert filesize(results_file) > 0 "MMseqs2 results file is empty"

        # Cleanup temporary files
        rm(tmp_dir, recursive=true, force=true)
        
        return results_file
        
    catch e
        @error "MMseqs2 execution failed" exception=e
        rethrow(e)
    finally
        # Cleanup temporary files
        rm(tmp_dir, recursive=true, force=true)
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Qualimap BAM QC on an alignment file.

# Arguments
- `bam::String`: Path to input BAM file
- `outdir::String`: Output directory (default: sibling `qualimap` folder)
- `threads::Int`: Number of threads to use (default: all available)
- `outformat::String`: Output formats (default: "PDF:HTML")
- `java_mem::String`: Java memory string (default: "4G")

# Returns
Named tuple with `report_pdf`, `report_txt`, and `coverage` file paths.
"""
function run_qualimap_bamqc(;
    bam::String,
    outdir::String = joinpath(dirname(bam), "qualimap"),
    threads::Int = get_default_threads(),
    outformat::String = "PDF:HTML",
    java_mem::String = "4G"
)
    @assert isfile(bam) "BAM file does not exist: $(bam)"
    Mycelia.add_bioconda_env("qualimap")
    mkpath(outdir)

    report_pdf = joinpath(outdir, "report.pdf")
    report_txt = joinpath(outdir, "genome_results.txt")
    coverage_txt = bam * ".genome_coverage.txt"

    if !isfile(report_pdf) || !isfile(report_txt)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n qualimap qualimap bamqc -nt $(threads) -bam $(bam) -outdir $(outdir) -outformat $(outformat) --output-genome-coverage $(coverage_txt) --java-mem-size=$(java_mem)`)
    end

    return (;report_pdf, report_txt, coverage=coverage_txt)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Return proportion of matched bases in alignment to total matches + edits.

Calculate the accuracy of a sequence alignment by computing the ratio of matched bases 
to total alignment operations (matches + edits).

# Arguments
- `alignment_result`: Alignment result object containing `total_matches` and `total_edits` fields

# Returns
Float64 between 0.0 and 1.0 representing alignment accuracy, where:
- 1.0 indicates perfect alignment (all matches)
- 0.0 indicates no matches
"""
function assess_alignment_accuracy(alignment_result)
    return alignment_result.total_matches / (alignment_result.total_matches + alignment_result.total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Used to determine which orientation provides an optimal alignment for initiating path likelihood analyses in viterbi analysis

Compare alignment scores between a query k-mer and an observed k-mer in both forward and
reverse complement orientations to determine optimal alignment.

# Arguments
- `kmer`: Query k-mer sequence to align
- `observed_kmer`: Target k-mer sequence to align against

# Returns
A tuple containing:
- `alignment_result`: The alignment result object for the optimal orientation
- `orientation`: Boolean indicating orientation (`true` = forward, `false` = reverse complement, `missing` = tied scores)

# Details
- Performs pairwise alignment in both orientations using `assess_alignment()`
- Calculates accuracy scores using `assess_alignment_accuracy()`
- For tied alignment scores, randomly selects one orientation
- Uses BioSequences.reverse_complement for reverse orientation comparison
"""
function assess_optimal_kmer_alignment(kmer, observed_kmer)

    forward_alignment_result = assess_alignment(kmer, observed_kmer)
    forward_alignment_accuracy = assess_alignment_accuracy(forward_alignment_result)

    reverse_alignment_result = assess_alignment(kmer, BioSequences.reverse_complement(observed_kmer))
    reverse_alignment_accuracy = assess_alignment_accuracy(reverse_alignment_result)

    if forward_alignment_accuracy > reverse_alignment_accuracy
        alignment_result = forward_alignment_result
        orientation = true
    elseif forward_alignment_accuracy < reverse_alignment_accuracy
        alignment_result = reverse_alignment_result
        orientation = false
    elseif forward_alignment_accuracy == reverse_alignment_accuracy
        alignment_result, orientation = rand(((forward_alignment_result, missing), (reverse_alignment_result, missing)))
    end

    return (alignment_result, orientation)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Aligns two sequences using the Levenshtein distance and returns the total number of matches and edits.

# Arguments
- `a::AbstractString`: The first sequence to be aligned.
- `b::AbstractString`: The second sequence to be aligned.

# Returns
- `NamedTuple{(:total_matches, :total_edits), Tuple{Int, Int}}`: A named tuple containing:
    - `total_matches::Int`: The total number of matching bases in the alignment.
    - `total_edits::Int`: The total number of edits (insertions, deletions, substitutions) in the alignment.
"""
function assess_alignment(a, b)
    pairwise_alignment = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), a, b)
    alignment_result = BioAlignments.alignment(pairwise_alignment)
    total_aligned_bases = BioAlignments.count_aligned(alignment_result)
    total_matches = Int(BioAlignments.count_matches(alignment_result))
    total_edits = Int(total_aligned_bases - total_matches)
    return (total_matches = total_matches, total_edits = total_edits)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Check if two biological sequences are equivalent, considering both direct and reverse complement matches.

# Arguments
- `a`: First biological sequence (BioSequence or compatible type)
- `b`: Second biological sequence (BioSequence or compatible type)

# Returns
- `Bool`: `true` if sequences are identical or if one is the reverse complement of the other, `false` otherwise
"""
function is_equivalent(a, b)
    a == b || a == BioSequences.reverse_complement(b)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# My standard pacbio aligning and sorting. No filtering done in this step.

# Use shell_only=true to get string command to submit to SLURM
# """
# function map_pacbio_reads(;
#         fastq,
#         reference_fasta,
#         temp_sam_outfile = fastq * "." * basename(reference_fasta) * "." * "minimap2.sam",
#         outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz"),
#         threads = get_default_threads(),
#         memory = Sys.total_memory(),
#         shell_only = false
#     )
#     # 4G is the default
#     # smaller, higher diversity databases do better with 5+ as the denominator - w/ <=4 they run out of memory
#     index_chunk_size = "$(Int(floor(memory/5e9)))G"
#     @show index_chunk_size
#     @show threads
#     Mycelia.add_bioconda_env("minimap2")
#     # Mycelia.add_bioconda_env("samtools")
#     Mycelia.add_bioconda_env("pigz")
#     if shell_only
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
#         """
#         return cmd
#     else
#         if !isfile(outfile)
#             map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_chunk_size) -ax map-hifi $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#             run(map)
#             run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`)
#             @assert isfile(outfile)
#         else
#             @info "$(outfile) already present"
#         end
#     end
# end

# """
# My standard pacbio aligning and sorting. No filtering done in this step.

# Use shell_only=true to get string command to submit to SLURM
# """
# function minimap_index_pacbio(;
#         reference_fasta,
#         outfile = replace(reference_fasta, Mycelia.FASTA_REGEX => ".pacbio.mmi"),
#         threads = get_default_threads(),
#         shell_only = false
#     )
#     Mycelia.add_bioconda_env("minimap2")
#     Mycelia.add_bioconda_env("samtools")
#     if shell_only
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-pb $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile) \\
#         | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile) \\
#         && rm $(temp_sam_outfile)
#         """
#         return cmd
#     else
#         if !isfile(outfile)
#             map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -ax map-pb $(reference_fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#             run(map)
#             p = pipeline(
#                 `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort --threads $(threads) $(temp_sam_outfile)`,
#                 `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -bh -o $(outfile)`
#             )
#             run(p)
#             rm(temp_sam_outfile)
#         else
#             @info "$(outfile) already present"
#         end
#     end
# end

# function map_short_reads()
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the minimap2 index size string based on available system memory.

# Arguments
- `system_mem_gb::Real`: Amount of memory (GB) to allocate for indexing.
- `denominator::Real`: Factor to scale the memory usage (default `Mycelia.DEFAULT_MINIMAP_DENOMINATOR`).

# Returns
- `String`: Value such as `"4G"` suitable for the minimap2 `-I` option.
"""
function system_mem_to_minimap_index_size(;system_mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85), denominator=DEFAULT_MINIMAP_DENOMINATOR)
    value = Int(floor(system_mem_gb/denominator))
    # 4G is the default
    # this value should be larger for larger memory machines, and smaller for smaller ones
    # it seems related to the total size of the sequences stored in memory, rather than the total size of the in-memory database
    return "$(value)G"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a minimap2 index for the provided reference sequence.

# Arguments
- `fasta::String`: Path to the reference FASTA.
- `mapping_type::String`: Preset (e.g. `"map-hifi"`).
- `mem_gb::Real`: Memory available in GB.
- `threads::Integer`: Number of threads.
- `as_string::Bool=false`: If true, return the command string instead of `Cmd`.
- `denominator::Real`: Scaling factor passed to `system_mem_to_minimap_index_size`.

# Returns
Named tuple `(cmd, outfile)` where `outfile` is the generated `.mmi` index path.
"""
function minimap_index(;fasta, mapping_type, mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85), threads=get_default_threads(), as_string=false, denominator=DEFAULT_MINIMAP_DENOMINATOR)
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    # if lr:hq, deal with : in the name
    index_file = "$(fasta).x" * replace(mapping_type, ":" => "-") * ".I$(index_size).mmi"
    if !isfile(index_file)
        Mycelia.add_bioconda_env("minimap2")
    end
    if as_string
        cmd = "$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -d $(index_file) $(fasta)"
    else
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -d $(index_file) $(fasta)`
    end
    outfile = index_file
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Map reads using an existing minimap2 index file.

# Arguments
- `fasta`: Path to the reference FASTA (used only if an index must be created).
- `mapping_type`: Minimap2 preset.
- `fastq`: Input reads.
- `index_file::String=""`: Optional prebuilt index path. If empty, one is created.
- `mem_gb`, `threads`, `as_string`, `denominator`: Parameters forwarded to `minimap_index`.

# Returns
Named tuple `(cmd, outfile)` producing a BAM file from the mapping.
"""
function minimap_map_with_index(;
        fasta,
        mapping_type,
        fastq,
        index_file="",
        mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85),
        threads=get_default_threads(),
        as_string=false,
        denominator=DEFAULT_MINIMAP_DENOMINATOR
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    
    if !isempty(index_file) && !isfile(index_file)
        error("user-specific index file $index_file does not exist")
    else
        # index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
        index_file_result = Mycelia.minimap_index(fasta=fasta, mapping_type=mapping_type, mem_gb = mem_gb, threads=threads)
        index_file = index_file_result.outfile
    end
    @show index_file
    @assert isfile(index_file)
    outfile = fastq * "." * basename(index_file) * "." * "minimap2.bam"
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(outfile).tmp \\
        | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(fastq) --split-prefix=$(outfile).tmp`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -`
        cmd = pipeline(map, compress)
    end
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate minimap2 alignment commands for sequence mapping.

aligning and compressing. No sorting or filtering.

Use shell_only=true to get string command to submit to SLURM

Creates a command to align reads in FASTQ format to a reference FASTA using minimap2, 
followed by SAM compression with pigz. Handles resource allocation and conda environment setup.

# Arguments
- `fasta`: Path to reference FASTA file
- `fastq`: Path to query FASTQ file
- `mapping_type`: Alignment preset ("map-hifi", "map-ont", "map-pb", "sr", or "lr:hq")
- `as_string`: If true, returns shell command as string; if false, returns command array
- `mem_gb`: Available memory in GB for indexing (defaults to system free memory)
- `threads`: Number of CPU threads to use (defaults to system threads)
- `denominator`: Divisor for calculating minimap2 index size

# Returns
Named tuple containing:
- `cmd`: Shell command (as string or array)
- `outfile`: Path to compressed output SAM file
"""
function minimap_map(;
        fasta,
        fastq,
        mapping_type,
        as_string=false,
        mem_gb=(Int(Sys.free_memory()) / 1e9),
        threads=get_default_threads(),
        denominator=DEFAULT_MINIMAP_DENOMINATOR,
        output_format="bam",
        sorted=true,
        quiet=true
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    @assert output_format in ["sam", "sam.gz", "bam"] "output_format must be 'sam', 'sam.gz', or 'bam'"

    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)

    # Construct output filename based on format
    base_name = fastq * "." * basename(fasta) * ".minimap2"
    if sorted
        base_name *= ".sorted"
    end
    outfile = base_name * "." * output_format

    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")

    # Build command based on output format
    if output_format == "sam"
        # Direct SAM output
        if as_string
            if sorted
                cmd = """
                $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp \\
                | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp -o $(outfile) -
                """
            else
                cmd = """
                $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp -o $(outfile)
                """
            end
        else
            if sorted
                map_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp`
                sort_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp -o $(outfile) -`
                cmd = pipeline(map_cmd, sort_cmd)
                if quiet
                    cmd = pipeline(cmd, stderr=devnull)
                end
            else
                if quiet
                    cmd = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp -o $(outfile)`, stdout=devnull, stderr=devnull)
                else
                    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp -o $(outfile)`
                end
            end
        end

    elseif output_format == "sam.gz"
        # Gzipped SAM output using samtools for proper BGZF compression
        if as_string
            if sorted
                cmd = """
                $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp \\
                | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp -O sam,level=6 -o $(outfile) -
                """
            else
                cmd = """
                $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp \\
                | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -O sam,level=6 -o $(outfile) -
                """
            end
        else
            if sorted
                map_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp`
                sort_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp -O sam,level=6 -o $(outfile) -`
                cmd = pipeline(map_cmd, sort_cmd)
                if quiet
                    cmd = pipeline(cmd, stderr=devnull)
                end
            else
                map_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp`
                view_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -O sam,level=6 -o $(outfile) -`
                cmd = pipeline(map_cmd, view_cmd)
                if quiet
                    cmd = pipeline(cmd, stderr=devnull)
                end
            end
        end

    else  # output_format == "bam"
        # BAM output (default behavior)
        if as_string
            if sorted
                cmd = """
                $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp \\
                | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp - \\
                | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -
                """
            else
                cmd = """
                $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp \\
                | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -
                """
            end
        else
            if sorted
                map_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp`
                sort_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp -`
                compress_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -`
                cmd = pipeline(map_cmd, sort_cmd, compress_cmd)
                if quiet
                    cmd = pipeline(cmd, stderr=devnull)
                end
            else
                map_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(outfile).tmp`
                compress_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -`
                cmd = pipeline(map_cmd, compress_cmd)
                if quiet
                    cmd = pipeline(cmd, stderr=devnull)
                end
            end
        end
    end

    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Map paired-end reads to a reference sequence using minimap2.

# Arguments
- `fasta::String`: Path to reference FASTA file
- `forward::String`: Path to forward reads FASTQ file
- `reverse::String`: Path to reverse reads FASTQ file
- `mem_gb::Integer`: Available system memory in GB
- `threads::Integer`: Number of threads to use
- `outdir::String`: Output directory (defaults to forward reads directory)
- `as_string::Bool=false`: Return command as string instead of Cmd array
- `mapping_type::String="sr"`: Minimap2 preset ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
- `denominator::Float64`: Memory scaling factor for index size

# Returns
Named tuple containing:
- `cmd`: Command(s) to execute (String or Array{Cmd})
- `outfile`: Path to output BAM file

# Notes
- Requires minimap2, and samtools conda environments
- Index file must exist at `\$(fasta).x\$(mapping_type).I\$(index_size).mmi`
"""
function minimap_map_paired_end_with_index(;
        forward,
        reverse,
        mem_gb=(Int(Sys.free_memory()) / 1e9),
        threads=get_default_threads(),
        outdir = dirname(forward),
        as_string=false,
        denominator=DEFAULT_MINIMAP_DENOMINATOR,
        fasta="",
        index_file = ""
    )
    # determine index_file and index_size
    if isempty(index_file)
        @assert !isempty(fasta) "must supply index file or fasta + mem_gb + denominator values to infer index file"
        index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
        index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    else
        # parse basename for ".I<index_size>.mmi"
        filename = basename(index_file)
        m = match(r"\.I(\d+G)\.mmi$", filename)
        if m === nothing
            error("Could not parse index size from index_file basename: $filename. Expected pattern '.I<index_size>.mmi' (e.g. genome.xsr.I42.mmi).")
        end
        parsed_index_size = m.captures[1]
        expected_index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
        @assert parsed_index_size == expected_index_size "Index size in index_file ($parsed_index_size) does not match expected index size ($expected_index_size)."
        index_size = parsed_index_size
    end
    # @show index_file
    @assert isfile(index_file) "$(index_file) not found!!"
    @assert isfile(forward) "$(forward) not found!!"
    @assert isfile(reverse) "$(reverse) not found!!"
    fastq_prefix = find_matching_prefix(basename(forward), basename(reverse))
    outfile = joinpath(outdir, fastq_prefix) * "." * basename(index_file) * "." * "minimap2.bam"
    # only run if we will need to do work
    if !isfile(outfile)
        Mycelia.add_bioconda_env("minimap2")
        Mycelia.add_bioconda_env("samtools")
    end
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x sr -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(outfile).tmp \\
        | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp - \\
        | $(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x sr -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(outfile).tmp`
        sort_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -T $(outfile).sort.tmp -`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools view -@ $(threads) -bS --no-header -o $(outfile) -`
        cmd = pipeline(map, sort_cmd, compress)
    end
    
    return (;cmd, outfile)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Merge FASTQs, map with minimap2, and split the sorted BAM back into per-sample files.

Read identifiers are rewritten with a sample-specific prefix to guarantee uniqueness and
enable efficient splitting. Supports single-end reads and paired collections (forward,
reverse). Additional minimap2 flags can be supplied via `minimap_extra_args`.

# Arguments
- `reference_fasta::Union{Nothing,AbstractString}`: Reference FASTA (optionally gzipped). Required unless `minimap_index` is supplied.
- `mapping_type::AbstractString`: Minimap2 preset (`"map-hifi"`, `"map-ont"`, `"map-pb"`, `"sr"`, `"lr:hq"`).

# Keywords
- `single_end_fastqs::Vector{<:AbstractString}`: FASTQs treated as single-end.
- `paired_end_fastqs::Vector{Tuple{<:AbstractString,<:AbstractString}}`: Forward/reverse pairs.
- `minimap_index::AbstractString=""`: Optional prebuilt minimap2 `.mmi` index. When supplied, it is used as the minimap2 target.
- `build_index::Bool=false`: When `true` and `minimap_index==""`, build a persistent `.mmi` via `Mycelia.minimap_index` and use it for mapping.
- `outdir::Union{Nothing,String}=nothing`: Destination for per-sample BAMs (defaults to each sample's FASTQ dir).
- `tmpdir::Union{Nothing,String}=nothing`: Location for prefixed FASTQs and merged BAM (defaults to a tempdir).
- `minimap_extra_args::Vector{<:AbstractString}=String[]`: Additional minimap2 flags (e.g., `["-N", "10"]`).
- `threads::Integer=get_default_threads()`, `mem_gb`, `denominator`: Control minimap2 index sizing.
- `id_delimiter::AbstractString="::"`: Separator between sample tag and original id.
- `write_read_map::Bool=false`: Optionally emit read-id map files.
- `read_map_format::Symbol=:arrow`: `:arrow` (buffered) or `:tsv` (streaming) mapping output.
- `gzip_prefixed_fastqs::Bool=true`: Gzip the prefixed FASTQs using external gzip/pigz.
- `gzip_read_map_tsv::Bool=true`: If `read_map_format==:tsv`, gzip the mapping TSVs using external gzip/pigz.
- `show_progress::Union{Bool,Nothing}=nothing`: If `true`, show progress meters for long-running local stages (prefixing and splitting). If `nothing`, auto-enables for large batches.
- `run_mapping::Bool=true`, `run_splitting::Bool=true`: Execute minimap2 stage and per-sample splits.
- `keep_prefixed_fastqs::Bool=false`: Retain prefixed intermediates if true.
- `force::Bool=false`: If `true`, regenerate intermediates even if present.
- `merged_bam::Union{Nothing,String}=nothing`: Override merged BAM path.
- `as_string::Bool=false`: Return command strings instead of `Cmd`/pipelines.
#
# Notes
# - This function passes all (prefixed) FASTQ paths directly to minimap2; extremely large numbers of inputs
#   can exceed the OS `ARG_MAX` limit (“argument list too long”). In that case, use a shorter `tmpdir` or
#   run in batches and merge BAMs with `samtools merge`.
# - For very large references, minimap2 may use a multi-part index; we pass `--split-prefix` to ensure SAM
#   output contains @SQ header lines required by samtools.
# - Per-sample BAMs are named as `<fastq_basename_without_gz>.<minimap_target_basename>.sorted.bam`.
#   Paired inputs join both FASTQ basenames (with `.gz` removed) using `__` before the index basename.

# Returns
Named tuple with commands, paths, and per-sample output metadata.
"""
function minimap_merge_map_and_split(;
    reference_fasta::Union{Nothing,AbstractString}=nothing,
    mapping_type::AbstractString,
    single_end_fastqs::Vector{<:AbstractString}=String[],
    paired_end_fastqs::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}}=Tuple{String,String}[],
    minimap_index::AbstractString="",
    build_index::Bool=false,
    outdir::Union{Nothing,String}=nothing,
    tmpdir::Union{Nothing,String}=nothing,
    minimap_extra_args::Vector{<:AbstractString}=String[],
    threads::Integer=get_default_threads(),
    mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85),
    denominator::Real=DEFAULT_MINIMAP_DENOMINATOR,
    id_delimiter::AbstractString="::",
    write_read_map::Bool=false,
    read_map_format::Symbol=:arrow,
    gzip_prefixed_fastqs::Bool=true,
    gzip_read_map_tsv::Bool=true,
    show_progress::Union{Bool,Nothing}=nothing,
    run_mapping::Bool=true,
    run_splitting::Bool=true,
    keep_prefixed_fastqs::Bool=false,
    force::Bool=false,
    merged_bam::Union{Nothing,String}=nothing,
    as_string::Bool=false
)
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    if isempty(single_end_fastqs) && isempty(paired_end_fastqs)
        error("Provide at least one FASTQ via single_end_fastqs or paired_end_fastqs")
    end
    isempty(minimap_index) && isnothing(reference_fasta) && error("Provide reference_fasta or minimap_index")
    if !isnothing(reference_fasta)
        @assert isfile(reference_fasta) "Reference FASTA not found: $reference_fasta"
    end
    if !isempty(minimap_index)
        @assert isfile(minimap_index) "minimap_index supplied but not found: $minimap_index"
    end
    tmpdir = isnothing(tmpdir) ? mktempdir() : tmpdir
    mkpath(tmpdir)
    nonempty_file(path::AbstractString) = isfile(path) && filesize(path) > 0

    # Helper to build sample tag from filenames
    function build_sample_tag(name::AbstractString)
        base = replace(basename(name), Mycelia.FASTQ_REGEX => "")
        return replace(base, r"[^\w\.\-]+" => "_")
    end

    strip_gz_extension(name::AbstractString) = replace(name, r"\.gz$" => "")
    function build_output_label(paths::AbstractVector{<:AbstractString})
        labels = strip_gz_extension.(basename.(paths))
        return length(labels) == 1 ? labels[1] : join(labels, "__")
    end

    prefix_jobs = NamedTuple[]
    sample_infos = NamedTuple[]
    mapping_suffix = if read_map_format == :arrow
        "arrow"
    else
        gzip_read_map_tsv ? "tsv.gz" : "tsv"
    end

    # Single-end inputs
    for fq in single_end_fastqs
        sample_tag = build_sample_tag(fq)
        output_label = build_output_label([fq])
        outfq = joinpath(tmpdir, sample_tag * (gzip_prefixed_fastqs ? ".prefixed.fq.gz" : ".prefixed.fq"))
        map_path = write_read_map ? joinpath(tmpdir, "$(sample_tag).read_map.$mapping_suffix") : nothing
        push!(prefix_jobs, (fastq=fq, sample_tag=sample_tag, out_fastq=outfq, mapping_out=map_path))
        push!(sample_infos, (
            sample_tag=sample_tag,
            output_label=output_label,
            source_fastqs=[fq],
            prefixed_fastqs=[outfq],
            mapping_files=map_path === nothing ? String[] : [map_path],
            paired=false,
            output_bam=nothing
        ))
    end

    # Paired-end inputs
    for (fwd, rev) in paired_end_fastqs
        sample_tag = find_matching_prefix(basename(fwd), basename(rev))
        if isempty(sample_tag)
            sample_tag = build_sample_tag(fwd)
        end
        sample_tag = replace(sample_tag, r"[^\w\.\-]+" => "_")
        output_label = build_output_label([fwd, rev])
        out1 = joinpath(tmpdir, "$(sample_tag).R1" * (gzip_prefixed_fastqs ? ".prefixed.fq.gz" : ".prefixed.fq"))
        out2 = joinpath(tmpdir, "$(sample_tag).R2" * (gzip_prefixed_fastqs ? ".prefixed.fq.gz" : ".prefixed.fq"))
        map1 = write_read_map ? joinpath(tmpdir, "$(sample_tag).R1.read_map.$mapping_suffix") : nothing
        map2 = write_read_map ? joinpath(tmpdir, "$(sample_tag).R2.read_map.$mapping_suffix") : nothing
        push!(prefix_jobs, (fastq=fwd, sample_tag=sample_tag, out_fastq=out1, mapping_out=map1))
        push!(prefix_jobs, (fastq=rev, sample_tag=sample_tag, out_fastq=out2, mapping_out=map2))
        push!(sample_infos, (
            sample_tag=sample_tag,
            output_label=output_label,
            source_fastqs=[fwd, rev],
            prefixed_fastqs=[out1, out2],
            mapping_files=String[m for m in (map1, map2) if m !== nothing],
            paired=true,
            output_bam=nothing
        ))
    end

    prefixed_fastqs = String[j.out_fastq for j in prefix_jobs]

    index_cmd = nothing
    index_file = nothing
    target = minimap_index
    if isempty(target)
        @assert !isnothing(reference_fasta)
        if build_index
            index_result = Mycelia.minimap_index(
                fasta=reference_fasta,
                mapping_type=mapping_type,
                mem_gb=mem_gb,
                threads=threads,
                as_string=false,
                denominator=denominator
            )
            index_file = index_result.outfile
            index_cmd = as_string ? string(index_result.cmd) : index_result.cmd
            if run_mapping && !isfile(index_file)
                run(index_result.cmd)
            end
            target = index_file
        else
            target = reference_fasta
        end
    end

    target_label = basename(target)
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    if isnothing(merged_bam)
        paired_str = join(sort(string.(paired_end_fastqs)), "\n")
        fingerprint = Mycelia.create_sha1_hash(
            string(
                "mapping_type=", mapping_type, "\n",
                "target=", target, "\n",
                "id_delimiter=", id_delimiter, "\n",
                "single_end_fastqs=", join(sort(single_end_fastqs), "\n"), "\n",
                "paired_end_fastqs=", paired_str, "\n",
                "minimap_extra_args=", join(minimap_extra_args, " ")
            );
            encoding=:hex,
            encoded_length=12,
            normalize_case=false,
            allow_truncation=true
        )
        merged_bam = joinpath(tmpdir, "$(target_label).$(mapping_type).$(fingerprint).joint.minimap2.sorted.bam")
    end

    # Prefixing (and optional gzip) in parallel across samples.
    merged_ready = nonempty_file(merged_bam) && !force
    need_prefixing = force || write_read_map || keep_prefixed_fastqs || (run_mapping && !merged_ready)
    if need_prefixing && !isempty(prefix_jobs)
        prefix_workers = min(length(prefix_jobs), max(1, Threads.nthreads()))
        show_prefix_progress = if show_progress === nothing
            length(prefix_jobs) >= 10
        else
            show_progress
        end
        total_prefix_bytes = show_prefix_progress ? sum(filesize(j.fastq) for j in prefix_jobs) : 0
        progress_chan = show_prefix_progress ? Channel{NamedTuple}(prefix_workers) : nothing
        progress_task = nothing
        if show_prefix_progress
            progress_task = @async begin
                bytes_done = 0
                files_done = 0
                p = ProgressMeter.Progress(total_prefix_bytes; desc="Prefixing FASTQs: ")
                for msg in progress_chan
                    bytes_done += msg.bytes
                    files_done += 1
                    ProgressMeter.update!(p, min(bytes_done, total_prefix_bytes); showvalues=[(:files, files_done), (:total, length(prefix_jobs))])
                end
                ProgressMeter.finish!(p)
            end
        end
        compress_threads_per_job = max(1, fld(max(1, threads), prefix_workers))
        sem = Base.Semaphore(prefix_workers)
        errs = Channel{Any}(length(prefix_jobs))
        @sync for job in prefix_jobs
            Threads.@spawn begin
                Base.acquire(sem)
                try
                    Mycelia.prefix_fastq_reads(
                        job.fastq;
                        sample_tag=job.sample_tag,
                        out_fastq=job.out_fastq,
                        mapping_out=job.mapping_out,
                        mapping_format=read_map_format,
                        id_delimiter=id_delimiter,
                        force=force,
                        compress_threads=compress_threads_per_job
                    )
                    if progress_chan !== nothing
                        try
                            put!(progress_chan, (bytes=filesize(job.fastq),))
                        catch
                        end
                    end
                catch err
                    put!(errs, err)
                finally
                    Base.release(sem)
                end
            end
        end
        if progress_chan !== nothing
            close(progress_chan)
        end
        if progress_task !== nothing
            wait(progress_task)
        end
        if isready(errs)
            err = take!(errs)
            throw(err)
        end
    end

    minimap_parts = [
        Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", "minimap2", "minimap2",
        "-t", string(threads), "-x", mapping_type, "-I$(index_size)", "-a"
    ]
    # Required for multi-part indices to ensure @SQ lines are present in SAM output.
    # Also safe for single-part indices and helps avoid samtools header/truncation errors.
    split_prefix = merged_bam * ".minimap2.split"
    push!(minimap_parts, "--split-prefix=$(split_prefix)")
    append!(minimap_parts, minimap_extra_args)
    push!(minimap_parts, target)
    append!(minimap_parts, prefixed_fastqs)

    # Preflight argv size for many-input runs to avoid opaque "argument list too long" errors.
    if run_mapping
        arg_max = Mycelia.get_arg_max()
        if arg_max !== nothing
            argv_bytes = Mycelia.estimate_argv_bytes(minimap_parts)
            headroom = 256_000 # env + runtime variance
            if argv_bytes > (arg_max - headroom)
                error(
                    "minimap2 argv likely exceeds ARG_MAX (estimated argv bytes=$(argv_bytes), ARG_MAX=$(arg_max)).\n" *
                    "Inputs: $(length(prefixed_fastqs)) FASTQ(s). Consider:\n" *
                    "  - Using a shorter `tmpdir` (and/or shorter input paths)\n" *
                    "  - Running samples in batches and merging BAMs with `samtools merge`\n"
                )
            end
        end
    end

    map_cmd = Cmd(minimap_parts)
    sort_cmd = Cmd([Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", "samtools", "samtools", "sort", "-@", string(threads), "-o", merged_bam, "-"])
    minimap_cmd = pipeline(map_cmd, sort_cmd)

    if run_mapping
        Mycelia.add_bioconda_env("minimap2")
        Mycelia.add_bioconda_env("samtools")
        if nonempty_file(merged_bam) && !force
            # resume/caching: keep existing merged BAM
        else
            try
                run(minimap_cmd)
                # Best-effort cleanup of minimap2 split-prefix temp files.
                try
                    split_base = basename(split_prefix)
                    for path in readdir(tmpdir; join=true)
                        startswith(basename(path), split_base) && rm(path; force=true)
                    end
                catch
                end
            catch err
                msg = sprint(showerror, err)
                if occursin("E2BIG", msg) || occursin("argument list too long", lowercase(msg))
                    error(
                        "minimap2 failed with an argument-length error (too many/too-long FASTQ paths).\n" *
                        "Mitigations:\n" *
                        "  - Use a shorter `tmpdir` (and/or shorten input paths)\n" *
                        "  - Run samples in batches and merge BAMs with `samtools merge`\n" *
                        "Original error: $(msg)"
                    )
                end
                if occursin("no SQ lines present in the header", msg) || occursin("samtools sort: truncated file", msg)
                    # Leave a clearer error and avoid caching a possibly corrupt BAM.
                    try
                        isfile(merged_bam) && rm(merged_bam; force=true)
                    catch
                    end
                    error(
                        "Mapping/sorting failed due to a missing SAM @SQ header (common with minimap2 multi-part indices without --split-prefix).\n" *
                        "This run now uses `--split-prefix=$(split_prefix)`; rerun with `force=true` if needed.\n" *
                        "Original error: $(msg)"
                    )
                end
                rethrow()
            end
        end
    end

    split_cmds = Dict{String,Cmd}()
    if run_splitting
        @assert isfile(merged_bam) "Merged BAM not found: $(merged_bam). Run mapping or supply existing BAM."
        Mycelia.add_bioconda_env("samtools")
        show_split_progress = if show_progress === nothing
            length(sample_infos) >= 10
        else
            show_progress
        end
        split_progress = show_split_progress ? ProgressMeter.Progress(length(sample_infos); desc="Splitting BAMs: ") : nothing
        split_done = 0
        for info in sample_infos
            sample_tag = info.sample_tag
            sample_outdir = isnothing(outdir) ? dirname(first(info.source_fastqs)) : outdir
            mkpath(sample_outdir)
            sample_bam = joinpath(sample_outdir, string(info.output_label, ".", target_label, ".sorted.bam"))
            awk_script = "BEGIN{FS=\"\t\"; OFS=\"\t\"} /^@/ {print; next} { split(" * "\\\$1" * ", parts, \"" * id_delimiter * "\"); if (parts[1]==\"" * sample_tag * "\") print }"
            view_cmd = Cmd([Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", "samtools", "samtools", "view", "-@", string(threads), "-h", merged_bam])
            filter_cmd = Cmd(["awk", awk_script])
            sort_split_cmd = Cmd([Mycelia.CONDA_RUNNER, "run", "--live-stream", "-n", "samtools", "samtools", "sort", "-@", string(threads), "-o", sample_bam, "-"])
            split_cmd = pipeline(view_cmd, filter_cmd, sort_split_cmd)
            split_cmds[sample_tag] = split_cmd
            if nonempty_file(sample_bam) && !force
                # resume/caching: keep existing per-sample BAM
            else
                run(split_cmd)
            end
            split_done += 1
            if split_progress !== nothing
                ProgressMeter.update!(split_progress, split_done)
            end
        end
        if split_progress !== nothing
            ProgressMeter.finish!(split_progress)
        end
    end

    sample_outputs = map(sample_infos) do info
        outdir_for_sample = isnothing(outdir) ? dirname(first(info.source_fastqs)) : outdir
        sample_bam = joinpath(outdir_for_sample, string(info.output_label, ".", target_label, ".sorted.bam"))
        merge(info, (;output_bam=sample_bam))
    end

    if !keep_prefixed_fastqs && run_mapping && run_splitting
        try
            for fq in prefixed_fastqs
                rm(fq; force=true)
            end
        catch
            # best-effort cleanup
        end
    end

    return (
        index_cmd = index_cmd,
        index_file = index_file,
        minimap_cmd = as_string ? string(minimap_cmd) : minimap_cmd,
        merged_bam = merged_bam,
        split_cmds = as_string ? Dict(k => string(v) for (k, v) in split_cmds) : split_cmds,
        prefixed_fastqs = prefixed_fastqs,
        sample_outputs = sample_outputs,
        tmpdir = tmpdir
    )
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Maps paired-end reads to a reference genome using minimap2 and compresses the output.

# # Arguments
# - `fasta::String`: Path to reference genome FASTA file
# - `forward::String`: Path to forward reads FASTQ file
# - `reverse::String`: Path to reverse reads FASTQ file  
# - `mem_gb::Integer`: Available system memory in GB
# - `threads::Integer`: Number of threads to use
# - `outdir::String`: Output directory (defaults to forward reads directory)
# - `as_string::Bool`: Return command as string instead of Cmd array
# - `mapping_type::String`: Mapping preset, e.g. "sr" for short reads (default)
# - `denominator::Float64`: Memory scaling factor for minimap2 index

# # Returns
# Named tuple containing:
# - `cmd`: Command(s) to execute (String or Vector{Cmd})
# - `outfile`: Path to compressed output SAM file (*.sam.gz)

# # Dependencies
# Requires bioconda packages: minimap2, samtools, pigz
# """
# function minimap_map_paired_end(;
#         fasta,
#         forward,
#         reverse,
#         mem_gb,
#         threads,
#         outdir = dirname(forward),
#         as_string=false,
#         mapping_type="sr",
#         denominator=Mycelia.DEFAULT_MINIMAP_DENOMINATOR
#     )
#     index_size = Mycelia.system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
#     @assert isfile(forward) "$(forward) not found!!"
#     @assert isfile(reverse) "$(reverse) not found!!"
#     fastq_prefix = Mycelia.find_matching_prefix(basename(forward), basename(reverse))
#     temp_sam_outfile = joinpath(outdir, fastq_prefix) * "." * "minimap2.sam"
#     # outfile = temp_sam_outfile
#     outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
#     Mycelia.add_bioconda_env("minimap2")
#     Mycelia.add_bioconda_env("samtools")
#     Mycelia.add_bioconda_env("pigz")
#     if as_string
#         cmd =
#         """
#         $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_size) -ax $(mapping_type) $(fasta) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
#         && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
#         """
#     else
#         map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -I$(index_size) -ax $(mapping_type) $(fasta) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
#         compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
#         cmd = [map, compress]
#     end
#     return (;cmd, outfile)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Clustal Omega multiple sequence alignment on a FASTA file.

# Arguments
- `fasta::String`: Path to input FASTA file
- `outfmt::String="clustal"`: Output format for the alignment

# Returns
- `String`: Path to the output alignment file

# Supported Output Formats
- `"fasta"`: FASTA format
- `"clustal"`: Clustal format
- `"msf"`: MSF format  
- `"phylip"`: PHYLIP format
- `"selex"`: SELEX format
- `"stockholm"`: Stockholm format
- `"vienna"`: Vienna format

# Notes
- Uses Bioconda to manage the Clustal Omega installation
- Caches results - will return existing output file if already generated
- Handles single sequence files gracefully by returning output path without error
"""
function run_clustal_omega(;fasta, outfmt="clustal")
    Mycelia.add_bioconda_env("clustalo")
    outfile = "$(fasta).clustal-omega.$(outfmt)"
    if !isfile(outfile)
        try
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n clustalo clustalo -i $(fasta) --outfmt $(outfmt) -o $(outfile)`)
        catch e
            # FATAL: File '...' contains 1 sequence, nothing to align
            return outfile
        end
    end
    return outfile
end

# function make_diamond_db(fasta_file, db_file=fasta_file)
#     @time run(`diamond makedb --in $(fasta_file) -d $(db_file)`)
# end

# # in order to change this to be a standard blast where we don't need all pairwise hits
# # just drop the parameters id, min-score, max-target-seqs
# function pairwise_diamond(joint_fasta_file)
#     if !isfile("$(joint_fasta_file).dmnd")
#         make_diamond_db(joint_fasta_file)
#     end
#     n_records = count_records(joint_fasta_file)
#     # max_target_seqs = Int(ceil(sqrt(n_records)))
#     # @show "here!"
#     sensitivity = "--iterate"
#     # --block-size/-b
#     # https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options
#     # set block size to total memory / 8
#     available_gigabytes = floor(Sys.free_memory() / 1e9)
#     block_size = floor(available_gigabytes / 8)
    
#     @time run(`diamond blastp $(sensitivity) --block-size $(block_size) --id 0 --min-score 0 --max-target-seqs $(n_records) --unal 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore -d $(joint_fasta_file).dmnd -q $(joint_fasta_file) -o $(joint_fasta_file).dmnd.tsv`)
#     # # pairwise output is all of the alignments, super helpful!
#     # # @time run(`diamond blastp $(sensitivity) --id 0 --min-score 0 --max-target-seqs $(N_RECORDS) --unal 1 --outfmt 0  -d $(joint_fasta_outfile).dmnd -q $(joint_fasta_outfile) -o $(joint_fasta_outfile).diamond.pairwise.txt`)
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run diamond search, returns path to diamond results.

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_diamond(;
#         identifier,
#         out_dir,
#         protein_fasta,
#         diamond_db,
#         force=false,
#         outfile="$(identifier).prodigal.faa.diamond.txt"
#     )
#     diamond_dir = mkpath("$(out_dir)/diamond")

#     # http://www.diamondsearch.org/index.php?pages/command_line_options/
#     # --block-size/-b #Block size in billions of sequence letters to be processed at a time.  
#     #     This is the main pa-rameter for controlling the program’s memory usage.  
#     #     Bigger numbers will increase the useof memory and temporary disk space, but also improve performance.  
#     #     The program can beexpected to use roughly six times this number of memory (in GB). So for the default value of-b2.0, 
#     #     the memory usage will be about 12 GB
#     system_memory_in_gigabytes = Int(Sys.total_memory()) / 1e9
#     # reference says 6 but let's round upwards towards 8
#     gb_per_block = 8
#     block_size = system_memory_in_gigabytes / gb_per_block
    
#     outfile = "$(diamond_dir)/$(outfile)"
    
#     if force || !isfile(outfile)
#         cmd = 
#         `diamond blastp
#         --threads $(get_default_threads())
#         --block-size $(block_size)
#         --db $(diamond_db)
#         --query $(protein_fasta)
#         --out $(outfile)
#         --evalue 0.001
#         --iterate
#         --outfmt 6 qseqid qtitle qlen sseqid sallseqid stitle salltitles slen qstart qend sstart send evalue bitscore length pident nident mismatch staxids
#         `

#         # --un                     file for unaligned queries
#         # --al                     file or aligned queries
#         # --unfmt                  format of unaligned query file (fasta/fastq)
#         # --alfmt                  format of aligned query file (fasta/fastq)
#         # --unal                   report unaligned queries (0=no, 1=yes)

# #         Value 6 may be followed by a space-separated list of these keywords:

# #         qseqid means Query Seq - id
# #         qtitle means Query title
# #         qlen means Query sequence length
# #         sseqid means Subject Seq - id
# #         sallseqid means All subject Seq - id(s), separated by a ';'
# #         stitle means Subject Title
# #         salltitles means All Subject Title(s), separated by a '<>'
# #         slen means Subject sequence length
# #         qstart means Start of alignment in query
# #         qend means End of alignment in query
# #         sstart means Start of alignment in subject
# #         send means End of alignment in subject
# #         evalue means Expect value
# #         bitscore means Bit score
# #         length means Alignment length
# #         pident means Percentage of identical matches
# #         nident means Number of identical matches
# #         mismatch means Number of mismatches
# #         staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
        
#         @time run(pipeline(cmd))
#     end
#     return outfile
# end

# """
# Downloads and unpacks the desired .tar.gz prebuilt kraken index

# Go to https://benlangmead.github.io/aws-indexes/k2 and identify the appropriate ".tar.gz" url download
# """
# function download_kraken_index(;url, directory="$(homedir())/workspace/kraken")
#     @assert occursin(r"\.tar\.gz$", url)
#     filename = last(split(url, '/'))
#     output_path = joinpath(directory, filename)
#     # @show url
#     # @show filename
#     # @show output_path
#     mkpath(directory)
#     extracted_directory = replace(basename(output_path), ".tar.gz" => "")
#     if !isdir(extracted_directory)
#         mkpath(extracted_directory)
#     end
#     if isempty(readdir(extracted_directory))
#         # download the file only if needed
#         if !isfile(output_path)
#             download(url, output_path)
#             # run(`wget -P $(directory) $(url)`)
#         end
#         run(`tar -xvzf $(output_path) -C $(extracted_directory)`)
#     end
#     return extracted_directory
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function run_mmseqs_easy_taxonomy(;out_dir, query_fasta, target_database, outfile, force=false)
#     add_bioconda_env("mmseqs2")
#     out_dir = mkpath(joinpath(out_dir, "mmseqs_easy_taxonomy"))
#     outfile = joinpath(out_dir, outfile * ".mmseqs_easy_taxonomy." * basename(target_database) * ".txt")
#     # note I tried adjusting all of the following, and none of them improved overall runtime
#     # in any meaningful way
#     # -s FLOAT                         Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]
#     # https://github.com/soedinglab/MMseqs2/issues/577#issuecomment-1191584081
#     # apparently orf-filter 1 speeds up by 50%!
#     # --orf-filter INT                 Prefilter query ORFs with non-selective search
#     #                               Only used during nucleotide-vs-protein classification
#     #                               NOTE: Consider disabling when classifying short reads [0]
#     # --lca-mode INT                   LCA Mode 1: single search LCA , 2/3: approximate 2bLCA, 4: top hit [3]
#     # --lca-search BOOL                Efficient search for LCA candidates [0]
#     # ^ this looks like it actually runs @ 1 with s=1.0 & --orf-filter=1
    
#     # 112 days to process 600 samples at this rate....
#     # 278 minutes or 4.5 hours for a single sample classification!!
#     # 16688.050696 seconds (1.43 M allocations: 80.966 MiB, 0.02% gc time, 0.00% compilation time)
#     # this is for default parameters
#     # lowering sensitivity and taking LCA
#     # 16590.725343 seconds (1.06 M allocations: 53.487 MiB, 0.00% compilation time)
#     # took just as long!
#     # difference was only 10 minutes
#     # 15903.218456 seconds (969.92 k allocations: 48.624 MiB, 0.01% gc time)
#     # use default parameters
    
#     if force || (!force && !isfile(outfile))
#         cmd = 
#         `$(CONDA_RUNNER) run --no-capture-output -n mmseqs2 mmseqs
#          easy-taxonomy
#          $(query_fasta)
#          $(target_database)
#          $(outfile)
#          $(joinpath(out_dir, "tmp"))
#         `
#         @time run(pipeline(cmd))
#     else
#         @info "target outfile $(outfile) already exists, remove it or set force=true to re-generate"
#     end
#     return outfile
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Runs the MMseqs2 easy-search command on the given query FASTA file against the target database.

Arguments
- query_fasta::String: Path to the query FASTA file.
- target_database::String: Path to the target database.
- out_dir::String: Directory to store the output file. Defaults to the directory of the query FASTA file.
- outfile::String: Base name of the output file (may include .gz). Defaults to a combination of the query FASTA and target database filenames with a .txt extension (and .gz automatically added if gzip=true).
- format_output::String: Format of the output. Defaults to a predefined set of fields.
- threads::Int: Number of CPU threads to use. Defaults to the number of CPU threads available.
- force::Bool: If true, forces re-generation of the (uncompressed) output even if it or a compressed variant already exists.
- gzip::Bool: If true (default), compress the final output with pigz. If true and the given outfile lacks a .gz suffix, it is appended. If false and the provided outfile ends with .gz, an error is raised.
- validate_compression::Bool: If true (default) and gzip=true, performs a pigz -t integrity test before deleting the uncompressed file. Validation reads the entire compressed file and verifies CRC32 and size modulo 2^32.
- keep_uncompressed::Bool: If true, retains the uncompressed file even after successful compression/validation. If validate_compression=true and keep_uncompressed=false, the uncompressed file is removed only after a successful pigz -t test.

Behavior Notes
- If gzip=true and the compressed file is missing but the corresponding uncompressed file exists (from a prior run or failed compression) and force=false, the function reuses the uncompressed file and (re)attempts compression instead of re-running MMseqs2.
- If force=true, MMseqs2 is always re-run and existing output artifacts are removed first.
- If a zero-byte compressed file exists, it is treated as corrupt, removed, and regeneration proceeds.
- Validation (pigz -t) is only run during the compression step (i.e., not re-run on an already existing compressed file being returned early).
- Compression uses pigz with the same thread count as MMseqs2.
- No staleness check is performed to confirm the uncompressed file matches current inputs.

Returns
- outfile_path::String: Path to the final (compressed or uncompressed) output file.

Integrity Caveats
- pigz -t validates internal integrity (CRC32 and length modulo 2^32) of the compressed file; it does not compare against a preserved checksum of the original unless you add such a mechanism externally.
"""
function run_mmseqs_easy_search(;
    query_fasta,
    target_database,
    out_dir = dirname(query_fasta),
    outfile = replace(basename(query_fasta), r"\.gz$"i => "") * ".mmseqs_easy_search." * basename(target_database) * ".txt",
    format_output = "query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,taxid,taxname",
    threads = get_default_threads(),
    force::Bool = false,
    gzip::Bool = true,
    validate_compression::Bool = true,
    keep_uncompressed::Bool = false,
)
    # Validate interaction of flags
    if !gzip && validate_compression
        @warn "validate_compression requested but gzip=false; validation will be skipped."
    end
    if !gzip && keep_uncompressed
        @warn "keep_uncompressed is irrelevant because gzip=false."
    end

    # Normalize outfile names and semantics
    if gzip
        if endswith(outfile, ".gz")
            compressed_name = outfile
            uncompressed_name = replace(outfile, r"\.gz$" => "")
        else
            uncompressed_name = outfile
            compressed_name = outfile * ".gz"
        end
    else
        if endswith(outfile, ".gz")
            error("gzip was set to false but the provided outfile ends with .gz: $(outfile)")
        else
            uncompressed_name = outfile
            compressed_name = nothing
        end
    end

    final_outfile_name = gzip ? compressed_name :: String : uncompressed_name
    final_outfile_path = joinpath(out_dir, final_outfile_name)
    work_outfile_path  = gzip ? joinpath(out_dir, uncompressed_name) : final_outfile_path
    tmp_dir = joinpath(out_dir, "tmp")

    # Ensure required environments (pigz only if needed)
    add_bioconda_env("mmseqs2")
    if gzip
        add_bioconda_env("pigz")
    end

    # Helper: compress (optionally validate) with pigz
    function _compress_with_pigz(uncompressed_path::String, final_path::String, threads::Int;
                                 validate::Bool,
                                 keep_uncompressed::Bool)
        if !isfile(uncompressed_path)
            error("Cannot compress: uncompressed file $(uncompressed_path) does not exist.")
        end

        # Remove zero-byte or clearly invalid existing compressed file artifact if present
        if isfile(final_path) && filesize(final_path) == 0
            @warn "Existing compressed file $(final_path) is zero bytes; removing before retry."
            rm(final_path; force=true)
        end

        # Skip if final already exists (rare race case)
        if !isfile(final_path)
            # Construct pigz command flags
            # -k if we must keep original either for validation or explicit retention
            keep_flag = (validate || keep_uncompressed) ? "-k" : ""
            cmd_pigz = `$(CONDA_RUNNER) run --no-capture-output -n pigz pigz $keep_flag -p $(threads) -f $(uncompressed_path)`
            @info "Compressing with pigz (threads=$(threads), validate=$(validate), keep_uncompressed=$(keep_uncompressed))"
            @time run(cmd_pigz)
        end

        if !isfile(final_path)
            error("Compression failed: expected compressed file $(final_path) was not created.")
        end

        if validate
            # Integrity test
            cmd_test = `$(CONDA_RUNNER) run --no-capture-output -n pigz pigz -t $(final_path)`
            @info "Validating compressed file integrity with pigz -t: $(final_path)"
            try
                @time run(cmd_test)
            catch err
                @error "pigz -t integrity test failed; retaining uncompressed file $(uncompressed_path)" error=err
                # Remove corrupted compressed file to avoid confusion
                try
                    rm(final_path; force=true)
                catch
                end
                error("Validation failed for compressed file $(final_path): $(err)")
            end

            # Remove uncompressed only if we do not want to keep it
            if !keep_uncompressed && isfile(uncompressed_path)
                try
                    rm(uncompressed_path; force=true)
                catch err
                    @warn "Failed to remove uncompressed file after successful validation: $(uncompressed_path)" error=err
                end
            end
        else
            # If not validating and not keeping, pigz (without -k) already removed uncompressed.
            # If keeping (keep_uncompressed=true), pigz was invoked with -k above.
            nothing
        end

        return final_path
    end

    # Handle zero-byte final compressed artifact early (force regeneration path)
    if gzip && isfile(final_outfile_path) && filesize(final_outfile_path) == 0
        @warn "Existing compressed file $(final_outfile_path) is zero bytes; removing for regeneration."
        rm(final_outfile_path; force=true)
    end

    # Early return if final artifact already exists and not forcing
    if !force && isfile(final_outfile_path)
        @info "Target outfile $(final_outfile_path) already exists; returning existing file."
        return final_outfile_path
    end

    # Determine if we can reuse an existing uncompressed file (avoid re-running MMseqs2)
    reused_uncompressed = false
    if gzip && !force && !isfile(final_outfile_path) && isfile(work_outfile_path)
        @info "Found existing uncompressed file $(work_outfile_path); will skip MMseqs2 run and proceed to compression."
        reused_uncompressed = true
    end

    # Force cleanup if required
    if force
        if isfile(final_outfile_path)
            rm(final_outfile_path; force=true)
        end
        if gzip && isfile(work_outfile_path)
            rm(work_outfile_path; force=true)
        end
    end

    # Run MMseqs2 only if we do not reuse an existing uncompressed result
    if !(gzip && reused_uncompressed)
        cmd_search =
            `$(CONDA_RUNNER) run --no-capture-output -n mmseqs2 mmseqs
                easy-search
                $(query_fasta)
                $(target_database)
                $(work_outfile_path)
                $(tmp_dir)
                --threads $(threads)
                --format-mode 4
                --format-output $(format_output)
                --start-sens 1 -s 7 --sens-steps 7
                --sort-results 1
                --remove-tmp-files 1
            `
        @info "Running MMseqs2 easy-search -> $(work_outfile_path)"
        @time run(cmd_search)
    end

    # Compression (and optional validation)
    if gzip
        if !isfile(final_outfile_path)
            _compress_with_pigz(work_outfile_path, final_outfile_path, threads;
                                validate = validate_compression,
                                keep_uncompressed = keep_uncompressed)
        else
            @info "Compressed file already present after generation: $(final_outfile_path)"
            # If user requested validation but file already existed, we skip re-validation to avoid a full read.
            # A custom flag could be added later to force re-validation of existing compressed outputs.
        end
    end

    # Cleanup tmp directory if present
    if isdir(tmp_dir)
        rm(tmp_dir; recursive=true)
    end

    return final_outfile_path
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)


Run the BLASTN (Basic Local Alignment Search Tool for Nucleotides) command with specified parameters.

# Arguments
- `outdir::String`: The output directory where the BLASTN results will be saved.
- `fasta::String`: The path to the input FASTA file containing the query sequences.
- `blastdb::String`: The path to the BLAST database to search against.
- `task::String`: The BLASTN task to perform. Default is "megablast".
- `force::Bool`: If true, forces the BLASTN command to run even if the output file already exists. Default is false.
- `remote::Bool`: If true, runs the BLASTN command remotely. Default is false.
- `wait::Bool`: If true, waits for the BLASTN command to complete before returning. Default is true.

# Returns
- `outfile::String`: The path to the output file containing the BLASTN results.

# Description
This function constructs and runs a BLASTN command based on the provided parameters.
It creates an output directory if it doesn't exist, constructs the output file path, and checks if the BLASTN command needs to be run based on the existence and size of the output file.
The function supports running the BLASTN command locally or remotely, with options to force re-running and to wait for completion.
"""
function run_blastn(;outdir=pwd(), fasta, blastdb, threads=min(get_default_threads(), 8), task="megablast", force=false, remote=false, wait=true)
    Mycelia.add_bioconda_env("blast")
    outdir = mkpath(outdir)
    outfile = "$(outdir)/$(basename(fasta)).blastn.$(basename(blastdb)).$(task).txt"
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    
    if force || need_to_run
        cmd = 
        `
        $(Mycelia.CONDA_RUNNER) run --live-stream -n blast blastn
        -num_threads $(threads)
        -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
        -query $(fasta)
        -db $(blastdb)
        -out $(outfile)
        -max_target_seqs 10
        -subject_besthit
        -task $(task)
        -evalue 0.001
        `
        @info "running cmd $(cmd)"
        @time run(pipeline(cmd), wait=wait)
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run a BLAST (Basic Local Alignment Search Tool) command with the specified parameters.

# Arguments
- `out_dir::String`: The output directory where the BLAST results will be stored.
- `fasta::String`: The path to the input FASTA file.
- `blast_db::String`: The path to the BLAST database.
- `blast_command::String`: The BLAST command to be executed (e.g., `blastn`, `blastp`).
- `force::Bool`: If `true`, forces the BLAST command to run even if the output file already exists. Default is `false`.
- `remote::Bool`: If `true`, runs the BLAST command remotely. Default is `false`.
- `wait::Bool`: If `true`, waits for the BLAST command to complete before returning. Default is `true`.

# Returns
- `outfile::String`: The path to the output file containing the BLAST results.

# Description
This function constructs and runs a BLAST command based on the provided parameters. It creates the necessary output directory, constructs the output file name, and determines whether the BLAST command needs to be run based on the existence and size of the output file. The function supports both local and remote execution of the BLAST command.

If `force` is set to `true` or the output file does not exist or is empty, the BLAST command is executed. The function logs the command being run and measures the time taken for execution. The output file path is returned upon completion.
"""
function run_blast(;out_dir, fasta, blast_db, blast_command, force=false, remote=false, wait=true)
    blast_dir = mkpath(joinpath(out_dir, blast_command))
    outfile = "$(blast_dir)/$(basename(fasta)).$(blast_command).$(basename(blast_db)).txt"
    if remote
        outfile = replace(outfile, ".txt" => ".remote.txt")
    end
    
    need_to_run = !isfile(outfile) || (filesize(outfile) == 0)
    
    # default max target seqs = 500, which seemed like too much
    # default evalue is 10, which also seems like too much
    if force || need_to_run
        if remote
            cmd = 
                `
                $(blast_command)
                -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
                -query $(fasta)
                -db $(basename(blast_db))
                -out $(outfile)
                -max_target_seqs 10
                -evalue 0.001
                -remote
                `
        else
            cmd = 
            `
            $(blast_command)
            -num_threads $(get_default_threads())
            -outfmt '7 qseqid qtitle sseqid sacc saccver stitle qlen slen qstart qend sstart send evalue bitscore length pident nident mismatch staxid'
            -query $(fasta)
            -db $(blast_db)
            -out $(outfile)
            -max_target_seqs 10
            -evalue 0.001
            `
        end
#         p = pipeline(cmd, 
#                 stdout="$(blastn_dir)/$(ID).blastn.out",
#                 stderr="$(blastn_dir)/$(ID).blastn.err")
        @info "running cmd $(cmd)"
        @time run(pipeline(cmd), wait=wait)
    end
    return outfile
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Imports results of Diamond (or blast) in outfmt 6 as a DataFrame
# """
# function read_diamond(path::String)
#   diamond_colnames = [ "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
  
#   diamond_coltypes = Dict(
#      1 => String, 
#      2 => String, 
#      3 => Float64,
#      4 => Int,
#      5 => Int,
#      6 => Int,
#      7 => Int,
#      8 => Int,
#      9 => Int,
#      10 => Int,
#      11 => String,
#      12 => String
#   )
#     return DataFrames.DataFrame(uCSV.read(open(path), delim = '\t', header = diamond_colnames, types = diamond_coltypes)...)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Expects output type 7 from BLAST, default output type 6 doesn't have the header comments and won't auto-parse

Parse a BLAST output file into a structured DataFrame.

# Arguments
- `blast_report::AbstractString`: Path to a BLAST output file in format 7 (tabular with comments)

# Returns
- `DataFrame`: Table containing BLAST results with columns matching the header fields.
  Returns empty DataFrame if no hits found.

# Details
- Requires BLAST output format 7 (`-outfmt 7`), which includes header comments
- Handles missing values (encoded as "N/A") automatically
- Infers column types based on BLAST field names
- Supports standard BLAST tabular fields including sequence IDs, scores, alignments and taxonomic information
"""
function parse_blast_report(blast_report)
    # example header line 
    # "# Fields: query id, subject id, subject acc., subject acc.ver, subject title, query length, subject length, q. start, q. end, s. start, s. end, evalue, bit score, score, alignment length, % identity, identical, mismatches, subject tax id"
    header_lines = collect(Iterators.filter(x -> occursin(r"# Fields:", x), eachline(blast_report)))
    if isempty(header_lines)
        @info "not hits found, returning empty table"
        return DataFrames.DataFrame()
    end
    header_line = first(header_lines)
    header = split(last(split(header_line, ": ")), ", ")
    blast_col_types = Dict(
        "query id" => String,
        "query title" => String,
        "subject id" => String,
        "subject gi" => String,
        "subject acc." => String,
        "subject acc.ver" => String,
        "subject title" => String,
        "query length" => Int,
        "subject length" => Int,
        "q. start" => Int,
        "q. end" => Int,
        "s. start" => Int,
        "s. end" => Int,
        "evalue" => Float64,
        "bit score" => Float64,
        "score" => Float64,
        "alignment length" => Int,
        "% identity" => Float64,
        "identical" => Int,
        "mismatches" => Int,
        "subject tax id" => Int,
        "subject sci name" => String,
        "subject com names" => String,
        "subject blast name" => String,
        "subject super kingdom" => String,
        "subject tax ids" => String,
        "subject sci names" => String,
        "subject com names" => String,
        "subject blast names" => String,
        "subject super kingdoms" => String,
        "subject title" => String,
        "subject titles" => String
    )
    data, _ = uCSV.read(
        blast_report,
        delim='\t',
        comment='#',
        skipmalformed=true,
        allowmissing=true,
        encodings=Dict("N/A" => missing),
        types=[blast_col_types[h] for h in header])
    return DataFrames.DataFrame(data, header, makeunique=true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse contig coverage statistics from a Qualimap BAM QC report file.

# Arguments
- `qualimap_report_txt::String`: Path to Qualimap bamqc report text file

# Returns
- `DataFrame`: Coverage statistics with columns:
  - `Contig`: Contig identifier
  - `Length`: Contig length in bases
  - `Mapped bases`: Number of bases mapped to contig
  - `Mean coverage`: Average coverage depth
  - `Standard Deviation`: Standard deviation of coverage
  - `% Mapped bases`: Percentage of total mapped bases on this contig

# Supported Assemblers
Handles output from both SPAdes and MEGAHIT assemblers:
- SPAdes format: NODE_X_length_Y_cov_Z
- MEGAHIT format: kXX_Y 

Parse the contig coverage information from qualimap bamqc text report, which looks like the following:

```
# this is spades
>>>>>>> Coverage per contig

	NODE_1_length_107478_cov_9.051896	107478	21606903	201.0355886786133	60.39424208607496
	NODE_2_length_5444_cov_1.351945	5444	153263	28.152645113886848	5.954250612823136
	NODE_3_length_1062_cov_0.154390	1062	4294	4.043314500941619	1.6655384692688975
	NODE_4_length_776_cov_0.191489	776	3210	4.13659793814433	2.252009588980858

# below is megahit
>>>>>>> Coverage per contig

	k79_175	235	3862	16.43404255319149	8.437436249612457
	k79_89	303	3803	12.551155115511552	5.709975376279777
	k79_262	394	6671	16.931472081218274	7.579217802849293
	k79_90	379	1539	4.060686015831134	1.2929729111266581
	k79_91	211	3749	17.767772511848342	11.899185693011933
	k79_0	2042	90867	44.49902056807052	18.356525483516613
```

To make this more robust, consider reading in the names of the contigs from the assembled fasta
"""
function parse_qualimap_contig_coverage(qualimap_report_txt)
    coverage_line_regex = r"\t.*?\t\d+\t\d+\t[\d\.]+\t[\d\.]+$"
    lines = filter(x -> occursin(coverage_line_regex, x), readlines("$(qualimap_report_txt)"))
    io = IOBuffer(join(map(x -> join(split(x, '\t')[2:end], '\t'), lines), '\n'))
    header = ["Contig", "Length", "Mapped bases", "Mean coverage", "Standard Deviation"]
    types = [String, Int, Int, Float64, Float64]
    data, _ = uCSV.read(io, delim='\t', types=types)
    qualimap_results = DataFrames.DataFrame(data, header)
    qualimap_results[!, "% Mapped bases"] = qualimap_results[!, "Mapped bases"] ./ sum(qualimap_results[!, "Mapped bases"]) .* 100
    return qualimap_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Imports results of fastani

Reads and processes FastANI output results from a tab-delimited file.

# Arguments
- `path::String`: Path to the FastANI output file

# Returns
DataFrame with columns:
- `query`: Original query filepath
- `query_identifier`: Extracted filename without extension
- `reference`: Original reference filepath
- `reference_identifier`: Extracted filename without extension
- `%_identity`: ANI percentage identity
- `fragments_mapped`: Number of fragments mapped
- `total_query_fragments`: Total number of query fragments

# Notes
- Expects tab-delimited input file from FastANI
- Automatically strips .fasta, .fna, or .fa extensions from filenames
- Column order is preserved as listed above
"""
function read_fastani(path::String)
    data, header = uCSV.read(path, delim='\t', typedetectrows=100)
    header = [
        "query",
        "reference",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]
    ani_table = DataFrames.DataFrame(data, header)
    ani_table[!, "query_identifier"] = replace.(basename.(ani_table[!, "query"]), r"\.(fasta|fna|fa)$" => "")
    ani_table[!, "reference_identifier"] = replace.(basename.(ani_table[!, "reference"]), r"\.(fasta|fna|fa)$" => "")
    columns = [
        "query",
        "query_identifier",
        "reference",
        "reference_identifier",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]
    ani_table = ani_table[!, columns]    
    return ani_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a Kraken taxonomic classification report into a structured DataFrame.

# Arguments
- `kraken_report::AbstractString`: Path to a tab-delimited Kraken report file

# Returns
- `DataFrame`: A DataFrame with the following columns:
  - `percentage_of_fragments_at_or_below_taxon`: Percentage of fragments covered
  - `number_of_fragments_at_or_below_taxon`: Count of fragments at/below taxon
  - `number_of_fragments_assigned_directly_to_taxon`: Direct fragment assignments
  - `rank`: Taxonomic rank
  - `ncbi_taxonid`: NCBI taxonomy identifier
  - `scientific_name`: Scientific name (whitespace-trimmed)

# Notes
- Scientific names are automatically stripped of leading/trailing whitespace
- Input file must be tab-delimited
"""
function read_kraken_report(kraken_report)
    kraken_report_header = [
        "percentage_of_fragments_at_or_below_taxon",
        "number_of_fragments_at_or_below_taxon",
        "number_of_fragments_assigned_directly_to_taxon",
        "rank",
        "ncbi_taxonid",
        "scientific_name"
    ]

    data, header = uCSV.read(kraken_report, delim='\t')
    kraken_report_table = DataFrames.DataFrame(data, kraken_report_header)
    kraken_report_table[!, "scientific_name"] = string.(strip.(kraken_report_table[!, "scientific_name"]))
    return kraken_report_table
end

# function diamond_line_to_named_tuple(diamond_line)
#     sline = split(line)
#     values_named_tuple = (
#         qseqid = sline[1],
#         sseqid = sline[2],
#         pident = parse(Float64, sline[3]),
#         length = parse(Int, sline[4]),
#         mismatch = parse(Int, sline[5]),
#         gapopen = parse(Int, sline[6]),
#         qlen = parse(Int, sline[7]),
#         qstart = parse(Int, sline[8]),
#         qend = parse(Int, sline[9]),
#         slen = parse(Int, sline[10]),
#         sstart = parse(Int, sline[11]),
#         send = parse(Int, sline[12]),
#         evalue = parse(Float64, sline[13]),
#         bitscore = parse(Float64, sline[14])
#         )
#     return values_named_tuple
# end

# function read_diamond_alignments_file(diamond_file)
#     column_names_to_types = [
#         "qseqid" => String,
#         "sseqid" => String,
#         "pident" => Float64,
#         "length" => Int,
#         "mismatch" => Int,
#         "gapopen" => Int,
#         "qlen" => Int,
#         "qstart" => Int,
#         "qend" => Int,
#         "slen" => Int,
#         "sstart" => Int,
#         "send" => Int,
#         "evalue" => Float64,
#         "bitscore" => Float64,
#     ]
#     types = Dict(i => t for (i, t) in enumerate(last.(column_names_to_types)))
    
#     data, header = uCSV.read(diamond_file, header=0, delim='\t', types = types)
#     header = first.(column_names_to_types)    
    
#     # data, header = uCSV.read(diamond_file, header=1, delim='\t', types = types)
#     # @assert header == first.(column_names_to_types)
    
#     table = DataFrames.DataFrame(data, header)
#     return table
# end

# function add_header_to_diamond_file(infile, outfile=replace(infile, ".tsv" => ".with-header.tsv"))
#     column_names = [
#         "qseqid",
#         "sseqid",
#         "pident",
#         "length",
#         "mismatch",
#         "gapopen",
#         "qlen",
#         "qstart",
#         "qend",
#         "slen",
#         "sstart",
#         "send",
#         "evalue",
#         "bitscore"
#     ]
#     # dangerous but fast
#     # try
#     #     inserted_text = join(columns_names, '\t') * '\n'
#     #     sed_cmd = "1s/^/$(inserted_text)/"
#     #     full_cmd = `sed -i $sed_cmd $infile`
#     # catch
#     open(outfile, "w") do io
#         println(io, join(column_names, "\t"))
#         for line in eachline(infile)
#             println(io, line)
#         end
#     end
#     return outfile
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function read_mmseqs_easy_search(mmseqs_file; top_hit_only=false)
#     mmseqs_results = DataFrames.DataFrame(uCSV.read(mmseqs_file, header=1, delim='\t')...)
#     if top_hit_only
#         gdf = DataFrames.groupby(mmseqs_results, "query")
#         for g in gdf
#             @assert issorted(g[!, "evalue"])
#         end
#         top_hits = DataFrames.combine(gdf, first)
#         mmseqs_results = top_hits
#     end
#     return mmseqs_results
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Read results from MMSeqs2 easy-search output file into a DataFrame.

# # Arguments
# - `mmseqs_file::String`: Path to the tab-delimited output file from MMSeqs2 easy-search

# # Returns
# - `DataFrame`: Contains search results with columns:
#   - `query`: Query sequence identifier
#   - `target`: Target sequence identifier
#   - `seqIdentity`: Sequence identity (0.0-1.0)
#   - `alnLen`: Alignment length
#   - `mismatch`: Number of mismatches
#   - `gapOpen`: Number of gap openings
#   - `qStart`: Query start position
#   - `qEnd`: Query end position
#   - `tStart`: Target start position 
#   - `tEnd`: Target end position
#   - `evalue`: Expected value
#   - `bits`: Bit score
# """
# function read_mmseqs_easy_search(mmseqs_file)
#     mmseqs_results = CSV.read(mmseqs_file, DataFrames.DataFrame, header=1, delim='\t')
#     return mmseqs_results
# end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Read results from MMSeqs2 easy-search output file (plain or gzipped) into a DataFrame
with optimized memory usage. Automatically detects if the file is gzipped based on
the '.gz' extension.

# Arguments
- `mmseqs_file::String`: Path to the tab-delimited output file from MMSeqs2 easy-search.
                         Can be a plain text file or a gzipped file (ending in .gz).

# Returns
- `DataFrame`: Contains search results with columns:
  - `query::String`: Query sequence identifier (pooled)
  - `target::String`: Target sequence identifier (pooled)
  - `seqIdentity::Float64`: Sequence identity (0.0-1.0)
  - `alnLen::Int`: Alignment length
  - `mismatch::Int`: Number of mismatches
  - `gapOpen::Int`: Number of gap openings
  - `qStart::Int`: Query start position
  - `qEnd::Int`: Query end position
  - `tStart::Int`: Target start position
  - `tEnd::Int`: Target end position
  - `evalue::Float64`: Expected value
  - `bits::Float64`: Bit score

# Remarks
- Ensure the `CodecZlib.jl` package is installed for gzipped file support.
"""
function read_mmseqs_easy_search(mmseqs_file::String)
    # Define the expected column types.
    # These keys MUST match the column headers in your MMSeqs2 output file.
    # Based on the error message, names are likely lowercase and 'seqIdentity' might be 'fident' or 'pident'.
    # We'll use 'fident' as it often represents fractional identity (0.0-1.0).
    col_types = Dict(
        "query" => String,        # Query sequence identifier (pooled)
        "qheader" => String,
        "qlen" => Int,
        "target" => String,       # Target sequence identifier (pooled)
        "tlen" => Int,
        "theader" => String,
        "nident" => Int,
        "fident" => Float64,      
        "pident" => Float64,
        "alnlen" => Int,          # Alignment length (lowercase 'l')
        "mismatch" => Int,        # Number of mismatches
        "gapopen" => Int,         # Number of gap openings (lowercase 'o')
        "qstart" => Int,          # Query start position (lowercase 's')
        "qend" => Int,            # Query end position (lowercase 'e')
        "tstart" => Int,          # Target start position (lowercase 's')
        "tend" => Int,            # Target end position (lowercase 'e')
        "evalue" => Float64,      # E-values
        "bits" => Float64,         # Bit scores
        "taxid" => Int, # Or String, depending on content
        "taxname" => String
    )

    # Common CSV reading options
    csv_options = (
        header=1,             # Assumes the first row contains column names
        delim='\t',           # Tab-delimited file
        types=col_types,      # Apply our defined column types
        pool=true,            # Pool string columns to save memory
        ntasks=get_default_threads(), # Utilize available CPU threads for parsing
    )

    local mmseqs_results::DataFrames.DataFrame

    # Check if the file is gzipped
    if endswith(lowercase(mmseqs_file), ".gz")
        # Open the gzipped file with a decompressor stream
        open(mmseqs_file, "r") do file_stream
            gzip_stream = CodecZlib.GzipDecompressorStream(file_stream)
            mmseqs_results = CSV.read(
                gzip_stream,
                DataFrames.DataFrame;
                csv_options... # Splat the common options
            )
        end
    else
        # Read as a plain text file
        mmseqs_results = CSV.read(
            mmseqs_file,
            DataFrames.DataFrame;
            csv_options... # Splat the common options
        )
    end

    return mmseqs_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
Stream a (plain or gzipped) MMSeqs2 easy-search tab-delimited output file and
return a DataFrame containing exactly one (top) hit per query.

Operation Modes
- Grouped mode (assume_grouped=true):
  All hits for a query are contiguous (block per query, blocks may appear in any order).
  The first row of each block is assumed to be the top hit; subsequent rows in the block
  are skipped (unless needed for validation checks).
- Unsorted mode (assume_grouped=false):
  Queries may interleave; every row is examined and the best row per query retained
  according to ranking rules.

Ranking
- rank_by = :bits   (higher bits better; tie-break by lower evalue if present).
- rank_by = :evalue (lower evalue better; tie-break by higher bits if present).

Validation (only in grouped mode)
- validation = :none disables all checks.
- validation = :fast or :full (currently identical) performs:
  1. Interleaving detection (query reappears after its block ended).
  2. Intra-block ordering check (later hit outranks the first).

Column Typing
- A default schema of known MMSeqs2 columns is provided (can be overridden via schema).
- Unknown columns default to String.
- keep_columns restricts output to a specified ordered subset; otherwise the full header order is kept.
- Numeric parse failures become missing if allow_parse_fail=true (columns will then have a Union element type).
- String columns can be pooled if pool_strings=true.
- When narrow_types=true a post-pass removes Missing from column types when no missings are present.

Arguments
- mmseqs_file::String: Input path (optionally .gz).
- assume_grouped::Bool=true: Whether hits for each query are contiguous.
- validation::Symbol=:fast: One of :none, :fast, :full.
- rank_by::Symbol=:bits: One of :bits or :evalue.
- keep_columns::Union{Nothing,Vector{String}}=nothing: Subset of columns to keep (ordered).
- schema::Union{Nothing,Dict{String,DataType}}=nothing: Override default column types.
- pool_strings::Bool=true: Pool string columns for memory savings.
- allow_parse_fail::Bool=true: If true, failed numeric parses become missing; otherwise errors are thrown.
- narrow_types::Bool=true: Post-pass type narrowing when no missings are present.

Returns
- DataFrames.DataFrame: One row per query (top hit), columns in original or requested order with concrete element types.

Remarks
- In grouped mode memory is O(#queries).
- In unsorted mode memory is also O(#queries) since only best rows are stored.
"""
function top_hits_mmseqs(
    mmseqs_file::String;
    assume_grouped::Bool = true,
    validation::Symbol = :fast,
    rank_by::Symbol = :bits,
    keep_columns::Union{Nothing,Vector{String}} = nothing,
    schema::Union{Nothing,Dict{String,DataType}} = nothing,
    pool_strings::Bool = true,
    allow_parse_fail::Bool = true,
    narrow_types::Bool = true
) :: DataFrames.DataFrame
    if !(rank_by in (:bits, :evalue))
        throw(ArgumentError("rank_by must be :bits or :evalue"))
    end
    if !(validation in (:none, :fast, :full))
        throw(ArgumentError("validation must be :none, :fast, or :full"))
    end

    selected_cols::Vector{String} = String[]
    selected_idx::Vector{Int} = Int[]
    # Allow union types by using Any for values
    col_types = Dict{String,Any}()
    stored_cols = Dict{String,AbstractVector}()

    io = endswith(lowercase(mmseqs_file), ".gz") ?
        CodecZlib.GzipDecompressorStream(open(mmseqs_file, "r")) :
        open(mmseqs_file, "r")

    try
        if eof(io)
            return DataFrames.DataFrame()
        end
        header_line = readline(io)
        if isempty(header_line)
            return DataFrames.DataFrame()
        end

        cols = split(header_line, '\t', keepempty=true)
        col_index = Dict{String,Int}(c => i for (i, c) in enumerate(cols))

        if !haskey(col_index, "query")
            throw(ArgumentError("Missing required column 'query'"))
        end
        if rank_by == :bits && !haskey(col_index, "bits")
            throw(ArgumentError("Missing column 'bits' required for rank_by=:bits"))
        elseif rank_by == :evalue && !haskey(col_index, "evalue")
            throw(ArgumentError("Missing column 'evalue' required for rank_by=:evalue"))
        end

        default_schema = Dict(
            "query"    => String,
            "qheader"  => String,
            "qlen"     => Int,
            "target"   => String,
            "theader"  => String,
            "tlen"     => Int,
            "nident"   => Int,
            "fident"   => Float64,
            "pident"   => Float64,
            "alnlen"   => Int,
            "mismatch" => Int,
            "gapopen"  => Int,
            "qstart"   => Int,
            "qend"     => Int,
            "tstart"   => Int,
            "tend"     => Int,
            "evalue"   => Float64,
            "bits"     => Float64,
            "taxid"    => Int,
            "taxname"  => String
        )
        if schema !== nothing
            for (k,v) in schema
                default_schema[k] = v
            end
        end

        if keep_columns === nothing
            selected_cols = copy(cols)
        else
            missing = setdiff(keep_columns, collect(keys(col_index)))
            if !isempty(missing)
                throw(ArgumentError("Requested keep_columns not found: $(missing)"))
            end
            selected_cols = keep_columns
        end
        selected_idx = map(c -> col_index[c], selected_cols)

        for c in selected_cols
            T = get(default_schema, c, String)
            if T === Int
                col_types[c] = allow_parse_fail ? Union{Int,Missing} : Int
            elseif T === Float64
                col_types[c] = allow_parse_fail ? Union{Float64,Missing} : Float64
            elseif T === String
                col_types[c] = String
            else
                col_types[c] = T
            end
        end

        if assume_grouped
            for c in selected_cols
                T = col_types[c]
                stored_cols[c] = Vector{T}()
            end
        end

        @inline function parse_cell(T, raw::AbstractString)
            if T === String
                return raw  # Will auto-convert SubString -> String on push if needed
            elseif T === Int
                return parse(Int, raw)
            elseif T === Float64
                return parse(Float64, raw)
            elseif T === Union{Int,Missing}
                raw == "" && return missing
                v = tryparse(Int, raw)
                return v === nothing ? missing : v
            elseif T === Union{Float64,Missing}
                raw == "" && return missing
                v = tryparse(Float64, raw)
                return v === nothing ? missing : v
            else
                return raw
            end
        end

        @inline function extract_rank(fields::Vector{SubString{String}})
            if rank_by == :bits
                bits_val = begin
                    b = tryparse(Float64, fields[col_index["bits"]]); b === nothing && (b = -Inf); b
                end
                if haskey(col_index, "evalue")
                    e_val = begin
                        e = tryparse(Float64, fields[col_index["evalue"]]); e === nothing && (e = Inf); e
                    end
                    return (bits_val, e_val)
                else
                    return (bits_val, Inf)
                end
            else
                e_val = begin
                    e = tryparse(Float64, fields[col_index["evalue"]]); e === nothing && (e = Inf); e
                end
                if haskey(col_index, "bits")
                    b_val = begin
                        b = tryparse(Float64, fields[col_index["bits"]]); b === nothing && (b = -Inf); b
                    end
                    return (e_val, -b_val)
                else
                    return (e_val, 0.0)
                end
            end
        end

        @inline function better_rank(a, b)
            if rank_by == :bits
                return (a[1] > b[1]) || (a[1] == b[1] && a[2] < b[2])
            else
                return (a[1] < b[1]) || (a[1] == b[1] && a[2] < b[2])
            end
        end

        current_query::Union{Nothing,String} = nothing
        current_first_rank = nothing
        closed_queries = validation == :none ? nothing : Set{String}()
        interleave_violation = false
        intra_rank_violation = false

        best_hits = Dict{String, Tuple{Any, Vector{Any}}}()

        function append_grouped!(fields)
            @inbounds for (j, c) in enumerate(selected_cols)
                raw = fields[selected_idx[j]]
                T = col_types[c]
                vec = stored_cols[c]
                push!(vec, parse_cell(T, raw))
            end
        end

        function build_row(fields)::Vector{Any}
            row = Vector{Any}(undef, length(selected_cols))
            @inbounds for (j, c) in enumerate(selected_cols)
                raw = fields[selected_idx[j]]
                T = col_types[c]
                row[j] = parse_cell(T, raw)
            end
            return row
        end

        for line in eachline(io)
            isempty(line) && continue
            fields = split(line, '\t', keepempty=true)
            length(fields) < length(cols) && continue
            q = fields[col_index["query"]]

            if assume_grouped
                if current_query === nothing
                    current_query = q
                    current_first_rank = extract_rank(fields)
                    append_grouped!(fields)
                elseif q == current_query
                    if validation != :none
                        rtuple = extract_rank(fields)
                        if better_rank(rtuple, current_first_rank)
                            intra_rank_violation = true
                        end
                    end
                    continue
                else
                    if validation != :none && closed_queries !== nothing
                        if in(q, closed_queries)
                            interleave_violation = true
                        end
                        push!(closed_queries, current_query)
                    end
                    current_query = q
                    current_first_rank = extract_rank(fields)
                    append_grouped!(fields)
                end
            else
                rtuple = extract_rank(fields)
                if haskey(best_hits, q)
                    old_rank, _ = best_hits[q]
                    if better_rank(rtuple, old_rank)
                        best_hits[q] = (rtuple, build_row(fields))
                    end
                else
                    best_hits[q] = (rtuple, build_row(fields))
                end
            end
        end

        if !assume_grouped
            n = length(best_hits)
            for c in selected_cols
                T = col_types[c]
                stored_cols[c] = Vector{T}(undef, n)
            end
            i = 1
            for (_, tup) in best_hits
                row = tup[2]
                @inbounds for (j, c) in enumerate(selected_cols)
                    stored_cols[c][i] = row[j]
                end
                i += 1
            end
        end

        if assume_grouped && validation != :none
            if interleave_violation
                @warn "Validation: Detected interleaving (a query block reappeared)."
            end
            if intra_rank_violation
                @warn "Validation: Detected a later hit within a block that outranks the first."
            end
        end

        df_cols = Vector{AbstractVector}(undef, length(selected_cols))
        @inbounds for (j, c) in enumerate(selected_cols)
            df_cols[j] = stored_cols[c]
        end
        df = DataFrames.DataFrame(df_cols, Symbol.(selected_cols))

        if narrow_types
            narrow_column_types!(df)
        end

        if pool_strings
            for c in names(df)
                if eltype(df[!, c]) <: AbstractString
                    df[!, c] = PooledArrays.PooledArray(df[!, c])
                end
            end
        end

        return df
    finally
        close(io)
    end
end

# Internal utility: narrow Union{T,Missing} columns without missings down to Vector{T}.
function narrow_column_types!(df::DataFrames.DataFrame)
    for c in names(df)
        T = eltype(df[!, c])
        if T isa Union
            nm = Base.nonmissingtype(T)
            if (nm !== Missing) && (nm !== T)
                has_missing = false
                @inbounds for v in df[!, c]
                    if v === missing
                        has_missing = true
                        break
                    end
                end
                if !has_missing
                    df[!, c] = convert(Vector{nm}, df[!, c])
                end
            end
        end
    end
    return df
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Parse MMseqs2 tophit alignment output file into a structured DataFrame.

# # Arguments
# - `tophit_aln::AbstractString`: Path to tab-delimited MMseqs2 alignment output file

# # Returns
# DataFrame with columns:
# - `query`: Query sequence/profile identifier
# - `target`: Target sequence/profile identifier  
# - `percent identity`: Sequence identity percentage
# - `alignment length`: Length of alignment
# - `number of mismatches`: Count of mismatched positions
# - `number of gaps`: Count of gap openings
# - `query start`: Start position in query sequence
# - `query end`: End position in query sequence
# - `target start`: Start position in target sequence
# - `target end`: End position in target sequence
# - `evalue`: E-value of alignment
# - `bit score`: Bit score of alignment
# """
# function parse_mmseqs_tophit_aln(tophit_aln)
#     data, header = uCSV.read(tophit_aln, delim='\t')
#     # (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.
#     # query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
#     header = [
#         "query",
#         "target",
#         "percent identity",
#         "alignment length",
#         "number of mismatches",
#         "number of gaps",
#         "query start",
#         "query end",
#         "target start",
#         "target end",
#         "evalue",
#         "bit score"
#     ]
#     DataFrames.DataFrame(data, header)
# end
