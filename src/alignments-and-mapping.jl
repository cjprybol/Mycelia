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

function system_mem_to_minimap_index_size(;system_mem_gb=(Int(Sys.total_memory()) / 1e9 * 0.85), denominator=DEFAULT_MINIMAP_DENOMINATOR)
    value = Int(floor(system_mem_gb/denominator))
    # 4G is the default
    # this value should be larger for larger memory machines, and smaller for smaller ones
    # it seems related to the total size of the sequences stored in memory, rather than the total size of the in-memory database
    return "$(value)G"
end

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