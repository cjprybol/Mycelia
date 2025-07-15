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
#         threads = Sys.CPU_THREADS,
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
#         threads = Sys.CPU_THREADS,
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
function minimap_index(;fasta, mapping_type, mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85), threads=Sys.CPU_THREADS, as_string=false, denominator=DEFAULT_MINIMAP_DENOMINATOR)
    Mycelia.add_bioconda_env("minimap2")
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    index_file = "$(fasta).x" * replace(mapping_type, ":" => "-") * ".I$(index_size).mmi"
    # if lr:hq, deal with : in the name
    # index_file = replace(mapping_type, ":" => "-")
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
        threads=Sys.CPU_THREADS,
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
        threads=Sys.CPU_THREADS,
        denominator=DEFAULT_MINIMAP_DENOMINATOR
    )
    @assert mapping_type in ["map-hifi", "map-ont", "map-pb", "sr", "lr:hq"]
    index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
    temp_sam_outfile = fastq * "." * basename(fasta) * "." * "minimap2.sam"
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(fasta) $(fastq) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = pipeline(map, compress)
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
- `outfile`: Path to compressed output SAM file (*.sam.gz)

# Notes
- Requires minimap2, samtools, and pigz conda environments
- Automatically compresses output using pigz
- Index file must exist at `\$(fasta).x\$(mapping_type).I\$(index_size).mmi`
"""
function minimap_map_paired_end_with_index(;
        forward,
        reverse,
        mem_gb=(Int(Sys.free_memory()) / 1e9),
        threads=Sys.CPU_THREADS,
        outdir = dirname(forward),
        as_string=false,
        denominator=DEFAULT_MINIMAP_DENOMINATOR,
        fasta="",
        index_file = ""
    )
    if isempty(index_file)
        @assert !isempty(fasta) "must supply index file or fasta + mem_gb + denominator values to infer index file"
        index_size = system_mem_to_minimap_index_size(system_mem_gb=mem_gb, denominator=denominator)
        index_file = "$(fasta).x$(mapping_type).I$(index_size).mmi"
    end
    # @show index_file
    @assert isfile(index_file) "$(index_file) not found!!"
    @assert isfile(forward) "$(forward) not found!!"
    @assert isfile(reverse) "$(reverse) not found!!"
    fastq_prefix = find_matching_prefix(basename(forward), basename(reverse))
    temp_sam_outfile = joinpath(outdir, fastq_prefix) * "." * basename(index_file) * "." * "minimap2.sam"
    # outfile = temp_sam_outfile
    outfile = replace(temp_sam_outfile, ".sam" => ".sam.gz")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    Mycelia.add_bioconda_env("pigz")
    if as_string
        cmd =
        """
        $(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile) \\
        && $(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)
        """
    else
        map = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -t $(threads) -x $(mapping_type) -I$(index_size) -a $(index_file) $(forward) $(reverse) --split-prefix=$(temp_sam_outfile).tmp -o $(temp_sam_outfile)`
        compress = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pigz pigz --processes $(threads) $(temp_sam_outfile)`
        cmd = pipeline(map, compress)
    end
    return (;cmd, outfile)
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
#         --threads $(Sys.CPU_THREADS)
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

# Arguments
- `query_fasta::String`: Path to the query FASTA file.
- `target_database::String`: Path to the target database.
- `out_dir::String`: Directory to store the output file. Defaults to the directory of the query FASTA file.
- `outfile::String`: Name of the output file. Defaults to a combination of the query FASTA and target database filenames.
- `format_output::String`: Format of the output. Defaults to a predefined set of fields.
- `threads::Int`: Number of CPU threads to use. Defaults to the number of CPU threads available.
- `force::Bool`: If true, forces the re-generation of the output file even if it already exists. Defaults to false.

# Returns
- `outfile_path::String`: Path to the generated output file.

# Notes
- Adds the `mmseqs2` environment using Bioconda if not already present.
- Removes temporary files created during the process.
"""
function run_mmseqs_easy_search(;
        query_fasta,
        target_database,
        out_dir=dirname(query_fasta),
        outfile=basename(query_fasta) * ".mmseqs_easy_search." * basename(target_database) * ".txt",
        format_output = "query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,taxid,taxname",
        threads = Sys.CPU_THREADS,
        force=false)
    
    add_bioconda_env("mmseqs2")
    outfile_path = joinpath(out_dir, outfile)
    tmp_dir = joinpath(out_dir, "tmp")
    if force || (!force && !isfile(outfile_path))
        cmd = 
        `$(CONDA_RUNNER) run --no-capture-output -n mmseqs2 mmseqs
            easy-search
            $(query_fasta)
            $(target_database)
            $(outfile_path)
            $(tmp_dir)
            --threads $(threads)
            --format-mode 4
            --format-output $(format_output)
            --start-sens 1 -s 7 --sens-steps 7
            --sort-results 1
            --remove-tmp-files 1
        `
        @time run(pipeline(cmd))
    else
        @info "target outfile $(outfile_path) already exists, remove it or set force=true to re-generate"
    end
    # we set remote tmp files = 1 above, but it still doesn't seem to work?
    if isdir(tmp_dir)
        rm(tmp_dir, recursive=true)
    end
    return outfile_path
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
function run_blastn(;outdir=pwd(), fasta, blastdb, threads=min(Sys.CPU_THREADS, 8), task="megablast", force=false, remote=false, wait=true)
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
            -num_threads $(Sys.CPU_THREADS)
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
        ntasks=Sys.CPU_THREADS, # Utilize available CPU threads for parsing
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

Parse MMseqs2 tophit alignment output file into a structured DataFrame.

# Arguments
- `tophit_aln::AbstractString`: Path to tab-delimited MMseqs2 alignment output file

# Returns
DataFrame with columns:
- `query`: Query sequence/profile identifier
- `target`: Target sequence/profile identifier  
- `percent identity`: Sequence identity percentage
- `alignment length`: Length of alignment
- `number of mismatches`: Count of mismatched positions
- `number of gaps`: Count of gap openings
- `query start`: Start position in query sequence
- `query end`: End position in query sequence
- `target start`: Start position in target sequence
- `target end`: End position in target sequence
- `evalue`: E-value of alignment
- `bit score`: Bit score of alignment
"""
function parse_mmseqs_tophit_aln(tophit_aln)
    data, header = uCSV.read(tophit_aln, delim='\t')
    # (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.
    # query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
    header = [
        "query",
        "target",
        "percent identity",
        "alignment length",
        "number of mismatches",
        "number of gaps",
        "query start",
        "query end",
        "target start",
        "target end",
        "evalue",
        "bit score"
    ]
    DataFrames.DataFrame(data, header)
end