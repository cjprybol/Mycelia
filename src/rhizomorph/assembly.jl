# ============================================================================
# Phase 3: Unified Assembly Pipeline Interface
# ============================================================================

"""
Assembly method enumeration for unified interface.
"""
@enum AssemblyMethod begin
    # Fixed-length graph types (assembly foundation)
    NgramGraph       # N-gram graph assembly (for unicode character analysis)
    KmerGraph        # K-mer graph assembly with DNAKmer/RNAKmer/AAKmer (for FASTA data)
    QualmerGraph     # Quality-aware k-mer graph assembly (for FASTQ data) - PRIMARY METHOD

    # Variable-length graph types (simplified products)
    StringGraph      # String graph assembly (simplified from N-gram graphs)
    BioSequenceGraph # BioSequence graph assembly (simplified from K-mer graphs)
    QualityBioSequenceGraph # Quality-aware BioSequence graph assembly (simplified from Qualmer graphs)

    # Hybrid approaches
    HybridOLC        # Hybrid OLC + qualmer graph approach
    MultiK           # Multi-k assembly with merging
end

function _graph_mode_symbol(graph_mode::GraphMode)
    if graph_mode == SingleStrand
        return :singlestrand
    elseif graph_mode == DoubleStrand
        return :doublestrand
    elseif graph_mode == Canonical
        return :canonical
    end
    return :singlestrand
end

function _qualmer_vertex_score(vertex_data)
    return Float64(Rhizomorph.count_evidence(vertex_data))
end

function _qualmer_edge_score(edge_data)
    return Float64(max(1, Rhizomorph.count_evidence(edge_data)))
end

function _qualmer_vertex_quality_scores(vertex_data)
    dataset_ids = Rhizomorph.get_all_dataset_ids(vertex_data)
    if isempty(dataset_ids)
        return UInt8[]
    end

    joint_quality = Rhizomorph.get_vertex_mean_quality(vertex_data, first(dataset_ids))
    if joint_quality === nothing
        return UInt8[]
    end

    return UInt8.(round.(clamp.(joint_quality, 0.0, 60.0)))
end

"""
Detect sequence type from input records using existing robust functions.
"""
function _detect_sequence_type(reads)
    if isempty(reads)
        return BioSequences.LongDNA{4}  # Default to DNA
    end

    # Get first sequence
    first_seq = if reads[1] isa FASTX.FASTA.Record
        FASTX.FASTA.sequence(reads[1])
    elseif reads[1] isa FASTX.FASTQ.Record
        FASTX.FASTQ.sequence(reads[1])
    elseif reads[1] isa String
        return String
    else
        error("Unsupported input type: $(typeof(reads[1]))")
    end

    # Use existing robust convert_sequence function which handles detection internally
    biosequence = Mycelia.convert_sequence(first_seq)
    return typeof(biosequence)
end

"""
Determine k-mer type from observations and k value.
"""
function _determine_kmer_type(observations, k::Int)
    if isempty(observations)
        return Kmers.DNAKmer{k}  # Default to DNA k-mers
    end

    # Extract first sequence to determine type
    first_seq = if observations[1] isa FASTX.FASTA.Record
        FASTX.FASTA.sequence(observations[1])
    elseif observations[1] isa FASTX.FASTQ.Record
        FASTX.FASTQ.sequence(observations[1])
    else
        error("Unsupported observation type: $(typeof(observations[1]))")
    end

    # Use existing robust convert_sequence function for type detection
    biosequence = Mycelia.convert_sequence(first_seq)

    # Determine appropriate k-mer type based on sequence type
    if biosequence isa BioSequences.LongDNA
        return Kmers.DNAKmer{k}
    elseif biosequence isa BioSequences.LongRNA
        return Kmers.RNAKmer{k}
    elseif biosequence isa BioSequences.LongAA
        return Kmers.AAKmer{k}
    else
        error("Unsupported sequence type for k-mer graph: $(typeof(biosequence))")
    end
end

"""
Enhanced assembly configuration structure with input validation.
"""
struct AssemblyConfig
    # Core parameters - exactly one of these should be specified
    k::Union{Int, Nothing}                                              # k-mer size (Nothing for overlap-based)
    min_overlap::Union{Int, Nothing}                                    # Min overlap (Nothing for k-mer based)

    # Input constraints
    sequence_type::Union{Type{<:BioSequences.BioSequence}, Type{String}}  # Type of input sequences
    graph_mode::GraphMode                                               # SingleStrand or DoubleStrand

    # Assembly parameters
    error_rate::Float64                     # Expected sequencing error rate
    min_coverage::Int                       # Minimum coverage for k-mer inclusion
    use_quality_scores::Bool                # Whether to use FASTQ quality scores
    polish_iterations::Int                  # Number of polishing iterations
    bubble_resolution::Bool                 # Whether to resolve bubble structures
    repeat_resolution::Bool                 # Whether to resolve repeat regions
    verbose::Bool                           # Whether to emit info logs during assembly
    tda::Union{Nothing, Mycelia.TDAConfig}  # Optional TDA configuration for topology-aware metrics/cleaning

    # Additive efficiency modes (all default to today's behavior, so existing
    # assemblies are byte-for-byte unchanged unless a caller opts in).
    # NOTE: with dedup_revcomp on, a reported contig may be the reverse-complement
    # (canonical) orientation of the sequence that was actually assembled.
    dedup_revcomp::Bool                     # Collapse RC-pair contigs to one canonical rep (post-assembly)
    compact_unitigs::Bool                   # Populate simplified_graph via linear-chain compaction
    memory_profile::Symbol                  # build_kmer_graph evidence footprint (:full|:lightweight|:ultralight|...)

    # Optional read-correction front-end (opt-in; default :none preserves today's
    # single-k-from-uncorrected-reads behavior byte-for-byte).
    corrector::Symbol                       # :none (default) or :iterative (mycelia_iterative_assemble)
    skip_solid::Bool                        # When corrector=:iterative, skip already-solid/clean reads

    # Constructor with validation
    function AssemblyConfig(;
            k::Union{Int, Nothing} = nothing,
            min_overlap::Union{Int, Nothing} = nothing,
            sequence_type::Union{Type{<:BioSequences.BioSequence}, Type{String}} = BioSequences.LongDNA{4},
            graph_mode::GraphMode = DoubleStrand,
            error_rate::Float64 = 0.01,
            min_coverage::Int = 3,
            use_quality_scores::Bool = true,
            polish_iterations::Int = 3,
            bubble_resolution::Bool = true,
            repeat_resolution::Bool = true,
            verbose::Bool = false,
            tda::Union{Nothing, Mycelia.TDAConfig} = nothing,
            dedup_revcomp::Bool = false,
            compact_unitigs::Bool = false,
            memory_profile::Symbol = :full,
            corrector::Symbol = :none,
            skip_solid::Bool = false
    )
        # Validation: Must specify exactly one of k or min_overlap
        if k === nothing && min_overlap === nothing
            k = 31  # Default to k-mer mode with k=31
        elseif k !== nothing && min_overlap !== nothing
            error("Cannot specify both k ($(k)) and min_overlap ($(min_overlap)). Choose one approach.")
        end

        # Validation: Check strand compatibility with sequence types.
        # DoubleStrand and Canonical both require a defined reverse complement,
        # so both are rejected for amino acids and general strings.
        if sequence_type <: BioSequences.LongAA && (graph_mode == DoubleStrand || graph_mode == Canonical)
            error("Amino acid sequences can only use SingleStrand mode (reverse complement undefined for proteins)")
        end
        if sequence_type == String && (graph_mode == DoubleStrand || graph_mode == Canonical)
            error("String sequences can only use SingleStrand mode (reverse complement undefined for general strings)")
        end

        # Validation: memory_profile must be one recognized by build_kmer_graph
        _valid_memory_profiles = (:full, :lightweight, :ultralight, :lightweight_quality, :ultralight_quality)
        if !(memory_profile in _valid_memory_profiles)
            error("memory_profile must be one of $(_valid_memory_profiles), got :$(memory_profile)")
        end

        # Validation: corrector front-end must be a recognized mode
        _valid_correctors = (:none, :iterative)
        if !(corrector in _valid_correctors)
            error("corrector must be one of $(_valid_correctors), got :$(corrector)")
        end

        # Validation: Parameter ranges
        if k !== nothing && (k < 1 || k > 64)
            error("k-mer size must be between 1 and 64, got k=$(k)")
        end
        if min_overlap !== nothing && min_overlap < 1
            error("min_overlap must be positive, got min_overlap=$(min_overlap)")
        end
        if !(0.0 <= error_rate <= 1.0)
            error("error_rate must be between 0.0 and 1.0, got error_rate=$(error_rate)")
        end
        if min_coverage < 1
            error("min_coverage must be positive, got min_coverage=$(min_coverage)")
        end

        # The additive efficiency modes (dedup_revcomp / compact_unitigs /
        # memory_profile) are only implemented on the non-quality k-mer path
        # (_assemble_kmer_graph). FASTQ input auto-sets use_quality_scores=true,
        # which dispatches to the qualmer arm and silently ignores these flags.
        # Warn unconditionally so a caller opting in on quality data is not left
        # believing the flag took effect. (Default config sets no flags, so this
        # never fires for existing assemblies.)
        if use_quality_scores &&
           (dedup_revcomp || compact_unitigs || memory_profile != :full)
            @warn "Efficiency modes (dedup_revcomp / compact_unitigs / " *
                  "memory_profile) are only implemented on the non-quality " *
                  "k-mer path and are ignored on the quality/qualmer path. " *
                  "Set use_quality_scores=false to use them."
        end

        new(
            k,
            min_overlap,
            sequence_type,
            graph_mode,
            error_rate,
            min_coverage,
            use_quality_scores,
            polish_iterations,
            bubble_resolution,
            repeat_resolution,
            verbose,
            tda,
            dedup_revcomp,
            compact_unitigs,
            memory_profile,
            corrector,
            skip_solid
        )
    end
end

"""
Assembly result structure containing contigs and metadata.
"""
struct AssemblyResult
    contigs::Vector{String}             # Final assembled contigs
    contig_names::Vector{String}        # Contig identifiers
    graph::Union{Nothing, MetaGraphsNext.MetaGraph}  # Complete assembly graph (optional)
    simplified_graph::Union{Nothing, MetaGraphsNext.MetaGraph}  # Simplified graph with collapsed paths
    paths::Dict{String, Vector}         # Path mappings for GFA P-lines (path_id -> vertex_sequence)
    assembly_stats::Dict{String, Any}   # Assembly statistics and metrics
    fastq_contigs::Vector{FASTX.FASTQ.Record}  # Quality-aware contigs (FASTQ format)
    gfa_compatible::Bool                # Whether graph structure is valid for GFA export

    function AssemblyResult(contigs::Vector{String}, contig_names::Vector{String};
            graph = nothing, simplified_graph = nothing, paths = Dict{String, Vector}(),
            assembly_stats = Dict{String, Any}(), fastq_contigs = FASTX.FASTQ.Record[],
            gfa_compatible = true)
        new(contigs, contig_names, graph, simplified_graph,
            paths, assembly_stats, fastq_contigs, gfa_compatible)
    end
end

"""
    get_fastq_contigs(result::AssemblyResult) -> Vector{FASTX.FASTQ.Record}

Extract quality-aware FASTQ contigs from assembly result if available.
Returns empty vector if no quality information was preserved during assembly.
"""
function get_fastq_contigs(result::AssemblyResult)
    return result.fastq_contigs
end

"""
    has_quality_information(result::AssemblyResult) -> Bool

Check if the assembly result preserves quality information from the original reads.
"""
function has_quality_information(result::AssemblyResult)
    return !isempty(result.fastq_contigs) &&
           get(result.assembly_stats, "quality_preserved", false)
end

"""
    write_fastq_contigs(result::AssemblyResult, output_file::String)

Write quality-aware contigs to a FASTQ file if quality information is available.
"""
function write_fastq_contigs(result::AssemblyResult, output_file::String)
    if !has_quality_information(result)
        error("Assembly result does not contain quality information")
    end

    FASTX.FASTQ.Writer(open(output_file, "w")) do writer
        for record in result.fastq_contigs
            write(writer, record)
        end
    end

    @info "Quality-aware contigs written to $(output_file)"
    return output_file
end

"""
    write_gfa(result::AssemblyResult, output_file::String)

Write assembly result to GFA (Graphical Fragment Assembly) format.
Exports both graph topology and path information leveraging existing infrastructure.
"""
function write_gfa(result::AssemblyResult, output_file::String)
    if isnothing(result.graph)
        error("AssemblyResult contains no graph - cannot write GFA format")
    end

    if !result.gfa_compatible
        @warn "AssemblyResult is marked as not GFA compatible - output may be invalid"
    end

    # Use the simplified graph if available, otherwise the full graph
    graph_to_write = isnothing(result.simplified_graph) ? result.graph :
                     result.simplified_graph

    # Write base GFA structure using existing function
    Rhizomorph.write_gfa_next(graph_to_write, output_file)

    # Append path information if available
    if !isempty(result.paths)
        _append_gfa_paths(output_file, result.paths, result.contigs, result.contig_names)
    end

    @info "Assembly written to GFA format: $(output_file)"
    return output_file
end

"""
    save_assembly(result::AssemblyResult, output_file::String)

Save complete assembly result using robust JLD2 serialization.
"""
function save_assembly(result::AssemblyResult, output_file::String)
    Mycelia.JLD2.save(output_file, "assembly_result", result)
    @info "Assembly saved to $(output_file)"
    return output_file
end

"""
    load_assembly(input_file::String) -> AssemblyResult

Load complete assembly result from JLD2 file.
"""
function load_assembly(input_file::String)
    return Mycelia.JLD2.load(input_file, "assembly_result")
end

"""
    has_graph_structure(result::AssemblyResult) -> Bool

Check if assembly result contains graph structure information.
"""
has_graph_structure(result::AssemblyResult) = !isnothing(result.graph)

"""
    has_simplified_graph(result::AssemblyResult) -> Bool

Check if assembly result contains simplified graph with collapsed paths.
"""
has_simplified_graph(result::AssemblyResult) = !isnothing(result.simplified_graph)

"""
    has_paths(result::AssemblyResult) -> Bool

Check if assembly result contains path mapping information.
"""
has_paths(result::AssemblyResult) = !isempty(result.paths)

"""
    validate_assembly_structure(result::AssemblyResult) -> Dict{String, Any}

Validate the internal consistency of an AssemblyResult structure.
Returns validation report complementing the existing validate_assembly function.
"""
function validate_assembly_structure(result::AssemblyResult)
    report = Dict{String, Any}(
        "valid" => true,
        "issues" => String[],
        "warnings" => String[]
    )

    # Check contig/name consistency
    if length(result.contigs) != length(result.contig_names)
        push!(report["issues"], "Contigs and contig_names have different lengths")
        report["valid"] = false
    end

    # Check graph consistency if present
    if !isnothing(result.graph) && !isnothing(result.simplified_graph)
        if typeof(result.graph) != typeof(result.simplified_graph)
            push!(report["warnings"], "Graph and simplified_graph have different types")
        end
    end

    # Check path consistency
    if !isempty(result.paths) && isnothing(result.graph)
        push!(report["warnings"], "Paths provided but no graph structure available")
    end

    # Check GFA compatibility
    if result.gfa_compatible && isnothing(result.graph)
        push!(report["issues"], "Marked as GFA compatible but no graph structure")
        report["valid"] = false
    end

    return report
end

"""
    _append_gfa_paths(gfa_file::String, paths::Dict{String, Vector}, contigs::Vector{String}, contig_names::Vector{String})

Append path information to an existing GFA file.
Adds GFA P-lines (path lines) that describe walks through the graph corresponding to assembled contigs.
"""
function _append_gfa_paths(gfa_file::String, paths::Dict{String, Vector},
        contigs::Vector{String}, contig_names::Vector{String})
    open(gfa_file, "a") do io
        println(io, "# Path information for assembled contigs")

        for (i, contig_name) in enumerate(contig_names)
            if haskey(paths, contig_name) && i <= length(contigs)
                path_vertices = paths[contig_name]
                if !isempty(path_vertices)
                    # Format: P <path_name> <vertex_list> <overlaps>
                    vertex_list = join(string.(path_vertices) .* "+", ",")
                    overlaps = repeat("*,", length(path_vertices) - 1) * "*"  # Default overlaps
                    println(io, "P\t$(contig_name)\t$(vertex_list)\t$(overlaps)")
                end
            end
        end
    end
end

function _log_info(config::AssemblyConfig, msg, args...)
    if config.verbose
        @info msg args...
    end
    return nothing
end

"""
Auto-configure assembly based on input type and parameters.
"""
function _auto_configure_assembly(
        reads; k = nothing, min_overlap = nothing, graph_mode = nothing, kwargs...)
    # Detect input format and sequence type
    sequence_type = _detect_sequence_type(reads)

    # Determine if quality scores are available
    use_quality_scores = all(r -> r isa FASTX.FASTQ.Record, reads)

    # Auto-detect graph mode if not specified
    if graph_mode === nothing
        graph_mode = if sequence_type <: BioSequences.LongAA || sequence_type == String
            SingleStrand  # AA and strings can only be single strand
        else
            DoubleStrand  # DNA/RNA default to double strand
        end
    end

    # Create config with detected parameters
    return AssemblyConfig(;
        k = k,
        min_overlap = min_overlap,
        sequence_type = sequence_type,
        graph_mode = graph_mode,
        use_quality_scores = use_quality_scores,
        kwargs...
    )
end

"""
    assemble_genome(reads; k=31, kwargs...) -> AssemblyResult
    assemble_genome(reads, config::AssemblyConfig) -> AssemblyResult

Unified genome assembly interface with auto-detection and type-stable dispatch.

# Auto-Detection Convenience Method
```julia
# Auto-detect sequence type and format, use k-mer approach
result = Mycelia.Rhizomorph.assemble_genome(fasta_records; k=25)

# Auto-detect sequence type and format, use overlap approach
result = Mycelia.Rhizomorph.assemble_genome(fasta_records; min_overlap=100)

# Override auto-detected parameters
result = Mycelia.Rhizomorph.assemble_genome(reads; k=31, graph_mode=Mycelia.Rhizomorph.SingleStrand, error_rate=0.005)
```

# Type-Stable Direct Method
```julia
# Explicit configuration for maximum performance
config = Mycelia.Rhizomorph.AssemblyConfig(k=25, sequence_type=BioSequences.LongDNA{4},
                                           graph_mode=Mycelia.Rhizomorph.DoubleStrand, use_quality_scores=true)
result = Mycelia.Rhizomorph.assemble_genome(reads, config)
```

# Arguments
- `reads`: Vector of FASTA/FASTQ records or file paths
- `config`: Assembly configuration (for type-stable version)

# Keyword Arguments (auto-detection version)
- `k`: k-mer size (mutually exclusive with min_overlap)
- `min_overlap`: Minimum overlap length (mutually exclusive with k)
- `graph_mode`: Mycelia.Rhizomorph.SingleStrand or Mycelia.Rhizomorph.DoubleStrand (auto-detected if not specified)
- `error_rate`, `min_coverage`, etc.: Assembly parameters

# Returns
- `AssemblyResult`: Structure containing contigs, names, and assembly metadata

# Details
This interface automatically:
1. **Detects sequence type**: DNA, RNA, AA, or String from input
2. **Chooses assembly method**: k-mer vs overlap-based on parameters
3. **Validates compatibility**: AA/String sequences -> SingleStrand only
4. **Dispatches optimally**: Based on input type and quality scores

**Method Selection Logic:**
- String input + k -> N-gram graph
- String input + min_overlap -> String graph
- BioSequence input + k + quality -> Qualmer graph
- BioSequence input + k -> K-mer graph
- BioSequence input + min_overlap + quality -> Quality BioSequence graph
- BioSequence input + min_overlap -> BioSequence graph
"""
function assemble_genome(reads; kwargs...)
    # Auto-detect configuration and dispatch to type-stable version
    config = _auto_configure_assembly(reads; kwargs...)
    return assemble_genome(reads, config)
end

"""
Type-stable main assembly function that dispatches based on configuration.
"""
function assemble_genome(reads, config::AssemblyConfig)
    # Opt-in read-correction front-end. The default (corrector=:none) falls
    # through to the byte-identical single-k pipeline below; only an explicit
    # corrector=:iterative diverts to the iterative+skip corrector.
    if config.corrector == :iterative
        return _assemble_with_iterative_corrector(reads, config)
    end

    _log_info(config, "Starting unified genome assembly", config.sequence_type,
        config.graph_mode, config.k, config.min_overlap)

    # Phase 1: Load and validate input
    observations = _prepare_observations(reads)
    _log_info(config, "Loaded $(length(observations)) sequence observations")

    # Phase 2: Type-stable dispatch based on config parameters
    result = if config.sequence_type == String
        if config.k !== nothing
            _assemble_ngram_graph(observations, config)  # N-gram
        else
            _assemble_string_graph(observations, config)  # String graph (overlap-based)
        end
    elseif config.sequence_type <: BioSequences.BioSequence
        if config.use_quality_scores
            if config.k !== nothing
                _assemble_qualmer_graph(observations, config)  # Qualmer
            else
                _assemble_quality_biosequence_graph(observations, config)  # Quality overlap-based
            end
        else
            if config.k !== nothing
                _assemble_kmer_graph(observations, config)  # K-mer
            else
                _assemble_biosequence_graph(observations, config)  # BioSequence overlap-based
            end
        end
    else
        throw(ArgumentError("Unsupported sequence type: $(config.sequence_type)"))
    end

    _log_info(config, "Assembly completed: $(length(result.contigs)) contigs generated")
    return result
end

"""
Write assembly input `reads` to a FASTQ file at `path` for the iterative
corrector, which consumes a FASTQ file path. Handles FASTQ records (written
verbatim), FASTA records and file paths (assigned a placeholder max quality,
since the corrector requires per-base quality strings).
"""
function _write_reads_to_fastq(reads, path::String)
    _placeholder_qual(n) = repeat("I", n)  # Phred+33 'I' == Q40
    placeholder_used = Ref(false)
    open(path, "w") do io
        writer = FASTX.FASTQ.Writer(io)
        if reads isa Vector{String}
            for file_path in reads
                if endswith(file_path, ".fastq") || endswith(file_path, ".fq")
                    open(FASTX.FASTQ.Reader, file_path) do reader
                        for record in reader
                            write(writer, record)
                        end
                    end
                else
                    open(FASTX.FASTA.Reader, file_path) do reader
                        for record in reader
                            seq = FASTX.FASTA.sequence(String, record)
                            placeholder_used[] = true
                            write(writer, FASTX.FASTQ.Record(
                                FASTX.FASTA.identifier(record), seq, _placeholder_qual(length(seq))))
                        end
                    end
                end
            end
        else
            for record in reads
                if record isa FASTX.FASTQ.Record
                    write(writer, record)
                elseif record isa FASTX.FASTA.Record
                    seq = FASTX.FASTA.sequence(String, record)
                    placeholder_used[] = true
                    write(writer, FASTX.FASTQ.Record(
                        FASTX.FASTA.identifier(record), seq, _placeholder_qual(length(seq))))
                else
                    error("Unsupported read type for iterative corrector: $(typeof(record))")
                end
            end
        end
        close(writer)
    end
    if placeholder_used[]
        @warn "corrector=:iterative: input lacked per-base quality (FASTA / quality-less); " *
              "assigned placeholder Q40. The corrector's quality model is degenerate (uniform) " *
              "for these reads — its decisions become coverage-driven only."
    end
    return path
end

"""
Route assembly through `Mycelia.mycelia_iterative_assemble` (the iterative +
skip-solid maximum-likelihood corrector), then RE-ASSEMBLE the corrected reads
into a real `AssemblyResult`.

`mycelia_iterative_assemble` is a read CORRECTOR: its `:final_assembly` is the
corrected READS (not contigs). This function reads the corrected FASTQ back and
re-assembles it through the naive path (`assemble_genome(...; corrector=:none)`),
so the returned `AssemblyResult` has real contigs + graph (`gfa_compatible=true`)
— NOT raw corrected reads (that v0 shortcut made QUAST comparisons invalid,
td-zru6). Corrector provenance (`corrector`, `skip_solid`, `k_progression`,
`corrected_read_count`) is stamped onto the re-assembly's `assembly_stats`.

The re-assembly deliberately lets `assemble_genome` auto-detect its graph mode
(DoubleStrand for DNA) rather than inheriting `config.graph_mode`: the corrector
needs `:canonical` for skip-solid, but the naive contig path is invalid under
`:canonical`. Corrected reads are FASTQ, so the re-assembly runs the same
quality-aware (qualmer) path a naive `assemble_genome` on FASTQ reads would —
i.e. it mirrors the naive-on-FASTQ baseline, keeping the comparison apples-to-
apples.
"""
function _assemble_with_iterative_corrector(reads, config::AssemblyConfig)
    _log_info(config, "Routing assembly through iterative corrector (corrector=:iterative)")

    # The corrector is k-mer based; a min_overlap-only config has no k, so the
    # k-progression floors at 13 and min_overlap is not used — surface that rather
    # than silently drop the OLC intent (review Important #3).
    if config.k === nothing
        @warn "corrector=:iterative is k-mer based; min_overlap is ignored and the k-progression floors at 13."
    end
    max_k = config.k === nothing ? 13 : max(config.k, 13)

    input_dir = mktempdir()
    output_dir = mktempdir()
    try
        temp_fastq = joinpath(input_dir, "corrector_input.fastq")
        _write_reads_to_fastq(reads, temp_fastq)

        result_dict = Mycelia.mycelia_iterative_assemble(
            temp_fastq;
            max_k = max_k,
            skip_solid = config.skip_solid,
            graph_mode = _graph_mode_symbol(config.graph_mode),
            # Opt into the tuned fast settings (td-q70n) explicitly on this route:
            # a coarse ~3-rung k-ladder + a low iteration cap. Kept out of the
            # mycelia_iterative_assemble DEFAULTS (which stay 10 / prime-walk) so
            # other callers are unchanged until the accuracy tradeoff is validated.
            n_k_rungs = 3,
            max_iterations_per_k = 2,
            verbose = false,
            enable_checkpointing = false,
            output_dir = output_dir
        )

        # mycelia_iterative_assemble is a read CORRECTOR: its :final_assembly is the
        # corrected READS, not an assembly. Re-assemble the corrected reads through
        # the naive path so assemble_genome returns real contigs + a graph —
        # otherwise downstream (QUAST, GFA) would treat raw corrected reads as an
        # assembly, which is not an apples-to-apples assembly result (td-zru6).
        corrected_fastq = get(result_dict[:metadata], :final_fastq_file, nothing)
        if corrected_fastq === nothing
            error("iterative corrector metadata is missing the :final_fastq_file key; " *
                  "keys=$(collect(keys(result_dict[:metadata])))")
        elseif !isfile(corrected_fastq)
            error("iterative corrector :final_fastq_file points at a nonexistent path: " *
                  "$(corrected_fastq) (output_dir=$(output_dir))")
        end
        # Eager `collect` inside a do-block: materializes ALL corrected reads into
        # memory (load-bearing — the finally-block rm's output_dir, so a lazy reader
        # would hit a deleted file) AND closes the stream (no leaked fd across many
        # assemblies, and closes on a malformed-FASTQ throw too).
        corrected_reads = open(FASTX.FASTQ.Reader, corrected_fastq) do reader
            collect(reader)
        end
        n_corrected = length(corrected_reads)
        # A corrector that silently ate every read (empty/header-only FASTQ) would
        # otherwise flow into assemble_genome([]) → a 0-contig "successful" assembly.
        # Fail loud instead of handing downstream a silently-empty assembly.
        if n_corrected == 0
            error("iterative corrector produced 0 corrected reads (from $(corrected_fastq)); " *
                  "refusing to return an empty assembly")
        end
        _log_info(config, "Corrected $(n_corrected) reads; re-assembling them (corrector=:none)")
        reassembly_k = config.k === nothing ? max_k : config.k
        # Re-assemble with AUTO-DETECTED graph_mode (not config.graph_mode): the
        # corrector needs :canonical for skip_solid, but the naive graph path emits
        # invalid contigs under :canonical, so the re-assembly must use the mode the
        # naive baseline uses (DoubleStrand for DNA), which auto-config selects.
        assembly = assemble_genome(corrected_reads;
            k = reassembly_k, corrector = :none)
        # Stamp the corrector provenance onto the real assembly's stats.
        assembly.assembly_stats["corrector"] = "iterative"
        assembly.assembly_stats["skip_solid"] = config.skip_solid
        assembly.assembly_stats["k_progression"] = result_dict[:k_progression]
        assembly.assembly_stats["corrected_read_count"] = n_corrected
        if isempty(assembly.contigs)
            @warn "corrector=:iterative re-assembled $(n_corrected) corrected reads into 0 " *
                  "contigs — the corrected read set did not assemble."
        end
        _log_info(config, "Re-assembled corrected reads into $(length(assembly.contigs)) contigs")
        return assembly
    finally
        # Prompt cleanup so repeated assemblies in a long-lived process do not leak
        # the input FASTQ + corrector output dirs until process exit (review #2).
        rm(input_dir; recursive = true, force = true)
        rm(output_dir; recursive = true, force = true)
    end
end

"""
    polish_assembly(assembly::AssemblyResult, reads; iterations=3) -> AssemblyResult

Polish assembled contigs using quality-aware error correction.

# Arguments  
- `assembly`: Initial assembly result to polish
- `reads`: Original reads for polishing (FASTQ with quality scores preferred)
- `iterations`: Number of polishing iterations (default: 3)

# Returns
- `AssemblyResult`: Polished assembly with improved accuracy

# Details
Uses Phase 2 enhanced Viterbi algorithms with quality score integration for:
- Error correction based on k-mer graph traversals
- Consensus calling from multiple observations
- Iterative improvement until convergence
"""
function polish_assembly(assembly::AssemblyResult, reads; iterations::Int = 3)
    @info "Starting assembly polishing" iterations

    observations = _prepare_observations(reads)
    polished_contigs = copy(assembly.contigs)

    for iter in 1:iterations
        @info "Polishing iteration $iter/$iterations"

        # Build k-mer graph from current contigs + reads
        graph = Rhizomorph.build_kmer_graph(
            vcat(observations, _contigs_to_records(polished_contigs)),
            31
        )

        # Polish each contig using Viterbi error correction
        for (i, contig) in enumerate(polished_contigs)
            polished_contigs[i] = _polish_contig_viterbi(contig, graph, observations)
        end
    end

    # Update assembly stats
    new_stats = merge(assembly.assembly_stats, Dict(
        "polishing_iterations" => iterations,
        "polished" => true
    ))

    @info "Polishing completed"
    return AssemblyResult(polished_contigs, assembly.contig_names;
        graph = assembly.graph, assembly_stats = new_stats)
end

"""
    validate_assembly(assembly::AssemblyResult; reference=nothing) -> Dict{String, Any}

Validate assembly quality using various metrics and optional reference comparison.

# Arguments
- `assembly`: Assembly result to validate
- `reference`: Optional reference sequence for comparison

# Returns  
- `Dict{String, Any}`: Comprehensive validation metrics

# Details
Computes assembly quality metrics including:
- N50, N90 statistics
- Total assembly length and number of contigs
- Coverage uniformity (if reference provided)
- Structural variant detection (if reference provided)
- Gap analysis and repeat characterization
"""
function validate_assembly(assembly::AssemblyResult; reference = nothing)
    @info "Validating assembly quality"

    contigs = assembly.contigs
    metrics = Dict{String, Any}()

    # Basic assembly statistics
    contig_lengths = [length(contig) for contig in contigs]
    total_length = sum(contig_lengths)
    sort!(contig_lengths, rev = true)

    metrics["num_contigs"] = length(contigs)
    metrics["total_length"] = total_length
    metrics["mean_contig_length"] = Statistics.mean(contig_lengths)
    metrics["max_contig_length"] = maximum(contig_lengths)
    metrics["min_contig_length"] = minimum(contig_lengths)

    # N-statistics
    metrics["N50"] = _calculate_n_statistic(contig_lengths, 0.5)
    metrics["N90"] = _calculate_n_statistic(contig_lengths, 0.9)
    metrics["L50"] = _calculate_l_statistic(contig_lengths, 0.5)
    metrics["L90"] = _calculate_l_statistic(contig_lengths, 0.9)

    # Reference-based validation (if provided)
    if !isnothing(reference)
        ref_metrics = _validate_against_reference(contigs, reference)
        merge!(metrics, ref_metrics)
    end

    @info "Assembly validation completed" metrics["N50"] metrics["num_contigs"]
    return metrics
end

# ============================================================================
# Phase 3: Helper Functions for Assembly Pipeline
# ============================================================================

"""
Prepare observations from various input formats (FASTA/FASTQ records or file paths).
"""
function _prepare_observations(reads)
    if reads isa Vector{String}
        # File paths - load FASTA/FASTQ records
        observations = FASTX.FASTA.Record[]
        for file_path in reads
            if endswith(file_path, ".fastq") || endswith(file_path, ".fq")
                # Load FASTQ and convert to FASTA records
                open(FASTX.FASTQ.Reader, file_path) do reader
                    for record in reader
                        push!(observations,
                            FASTX.FASTA.Record(FASTX.FASTQ.identifier(record),
                                FASTX.FASTQ.sequence(record)))
                    end
                end
            else
                # Load FASTA records
                open(FASTX.FASTA.Reader, file_path) do reader
                    for record in reader
                        push!(observations, record)
                    end
                end
            end
        end
        return observations
    elseif reads isa Vector{<:FASTX.FASTA.Record}
        return reads
    elseif reads isa Vector{<:FASTX.FASTQ.Record}
        # Convert FASTQ to FASTA records
        return [FASTX.FASTA.Record(FASTX.FASTQ.identifier(record), FASTX.FASTQ.sequence(record))
                for record in reads]
    else
        throw(ArgumentError("Unsupported reads format: $(typeof(reads))"))
    end
end

"""
String graph assembly implementation (variable-length simplified from N-gram graphs).
"""
function _assemble_string_graph(observations, config)
    _log_info(config, "Using string graph assembly strategy (variable-length OLC)")

    strings = [String(FASTX.FASTA.sequence(obs)) for obs in observations]
    min_overlap = isnothing(config.min_overlap) ? 1 : config.min_overlap
    string_graph = Rhizomorph.build_string_graph(strings; min_overlap = min_overlap)

    if config.bubble_resolution
        _simplify_string_graph(string_graph)
    end

    contig_paths = Rhizomorph.find_contigs_next(string_graph; min_contig_length = 1)
    contigs = [string(contig.sequence) for contig in contig_paths]
    if isempty(contigs)
        contigs = strings
    end
    contig_names = ["string_contig_$i" for i in 1:length(contigs)]

    stats = Dict{String, Any}(
        "method" => "StringGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "unicode_strings",
        "min_overlap" => min_overlap,
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Mycelia.Dates.now())
    )

    return AssemblyResult(contigs, contig_names; graph = string_graph, assembly_stats = stats)
end

"""
K-mer graph assembly implementation (fixed-length k-mer foundation).
"""
function _assemble_kmer_graph(observations, config)
    _log_info(config, "Using k-mer graph assembly strategy (fixed-length k-mer foundation)")
    mode = _graph_mode_symbol(config.graph_mode)
    # Canonical graph_mode now supports correct contig reconstruction. The
    # undirected canonical graph merges each k-mer with its reverse complement
    # onto one vertex; find_eulerian_paths_next handles undirected graphs
    # (degree-parity feasibility + symmetric Hierholzer) and path_to_sequence /
    # generate_contig_sequence recover each k-mer's orientation from the (k-1)
    # overlap, reverse-complementing where the overlap is on the reverse strand.
    # Canonical reconstruction therefore matches DoubleStrand (verified in
    # test/4_assembly/rhizomorph_efficiency_modes_test.jl Mode 2), so the result
    # is flagged valid and no warning is emitted.
    reconstruction_valid = true
    # Mode 3a (opt-in): memory_profile selects the k-mer graph's evidence footprint
    # (:full default, or :lightweight / :ultralight / *_quality). This is an internal
    # representation change; the assembled contigs are expected to be identical.
    graph = Rhizomorph.build_kmer_graph(
        observations, config.k; mode = mode, memory_profile = config.memory_profile)

    # NOTE: detect_bubbles_next / resolve_repeats_next are intentionally NOT run
    # here. They return analysis structures that this function only logged and
    # never used (they do not simplify the graph), while scaling ~O(V^2) — on a
    # real read graph (errors create thousands of branch vertices) they took tens
    # of minutes on a single phage genome. They remain available as standalone
    # analysis functions; wire them in only once they actually mutate the graph.
    # Until then the config flags below have no effect on the k-mer arm; surface
    # that so a caller setting them does not silently get unmodified behavior.
    if config.bubble_resolution || config.repeat_resolution
        _log_info(config,
            "bubble/repeat resolution not yet implemented for k-mer graphs; flags ignored")
    end

    # Find contigs. find_eulerian_paths_next is a fast path that only succeeds on
    # (near-)balanced graphs — i.e. toy or error-free inputs; it returns no paths
    # the moment any vertex is unbalanced, which is guaranteed on real read graphs
    # (errors create tips/bubbles). On such inputs the unitig fallback below
    # (find_contigs_next via _generate_contigs_probabilistic) is the real contig
    # generator.
    paths = Rhizomorph.find_eulerian_paths_next(graph)

    # Convert paths to sequences
    contigs = String[]
    for path in paths
        if length(path) > 1
            sequence = Rhizomorph.path_to_sequence(path, graph)
            push!(contigs, string(sequence))
        end
    end

    # If no Eulerian paths, use probabilistic walks
    if isempty(contigs)
        _log_info(config, "No Eulerian paths found, using linear contigs")
        contigs = _generate_contigs_probabilistic(graph, config)
    end

    # Mode 1 (opt-in): collapse contigs that are reverse complements of each other
    # to a single canonical representative. Default (dedup_revcomp=false) leaves the
    # RC pairs a DoubleStrand assembly naturally emits, preserving today's behavior.
    if config.dedup_revcomp
        n_before = length(contigs)
        contigs = _dedup_reverse_complements(contigs)
        _log_info(config, "Reverse-complement dedup: $(n_before) -> $(length(contigs)) contigs")
    end

    contig_names = ["kmer_contig_$i" for i in 1:length(contigs)]

    # Mode 3b (opt-in): populate simplified_graph via linear-chain (unitig)
    # compaction. collapse_linear_chains! is a no-op on fixed-length k-mer graphs
    # (collapsing a linear chain produces a sequence longer than k, which would
    # change the vertex label type from Kmer to BioSequence and violate
    # MetaGraphsNext's parametric label_type). The keystone fix is to first run
    # convert_fixed_to_variable(), which relabels each k-mer vertex as a
    # variable-length BioSequence vertex (overlap = k-1) while preserving evidence
    # and edge topology; collapse_linear_chains! can then merge non-branching runs
    # into single unitig vertices, reducing the vertex count. The k-mer `graph`
    # returned in the result is unchanged (n_full is measured against it), and the
    # emitted contigs are untouched (they come from the k-mer traversal above), so
    # the "contigs unchanged" contract still holds.
    simplified = nothing
    unitig_compaction_effective = false
    if config.compact_unitigs
        variable_graph = Rhizomorph.convert_fixed_to_variable(graph)
        simplified = Rhizomorph.collapse_linear_chains!(variable_graph)
        n_full = length(MetaGraphsNext.labels(graph))
        n_simpl = length(MetaGraphsNext.labels(simplified))
        # Effective when the collapsed variable-length graph has strictly fewer
        # vertices than the fixed-length k-mer graph (a linear tiling collapses to
        # a handful of unitigs). Remains false only when the graph is already
        # maximally branchy (no non-branching runs to merge).
        unitig_compaction_effective = n_simpl < n_full
        _log_info(config,
            "Unitig compaction: $(n_full) k-mer vertices -> $(n_simpl) unitig " *
            "vertices (via convert_fixed_to_variable + collapse_linear_chains!)")
        if !unitig_compaction_effective
            @warn "compact_unitigs requested but the graph had no non-branching " *
                  "runs to collapse: simplified_graph has the same vertex count " *
                  "as the k-mer graph."
        end
    end

    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "KmerGraph",
        "graph_type" => "fixed_length",
        "vertex_type" => isempty(MetaGraphsNext.labels(graph)) ? "unknown" :
                         string(typeof(first(MetaGraphsNext.labels(graph)))),
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(observations),
        # Record that the requested bubble/repeat resolution was NOT applied, so
        # callers can detect the ignored flags regardless of the `verbose` setting
        # (the _log_info above is silent when verbose=false).
        "bubble_resolution_requested" => config.bubble_resolution,
        "repeat_resolution_requested" => config.repeat_resolution,
        "graph_cleaning_applied" => false,
        "unitig_compaction_requested" => config.compact_unitigs,
        "unitig_compaction_effective" => unitig_compaction_effective,
        # Always true here: canonical (undirected) traversal is now orientation-
        # aware, so its contigs are a valid reconstruction (matching DoubleStrand).
        "reconstruction_valid" => reconstruction_valid,
        "assembly_date" => string(Mycelia.Dates.now())
    )

    return AssemblyResult(contigs, contig_names;
        graph = graph, simplified_graph = simplified, assembly_stats = stats)
end

"""
Quality-aware k-mer graph assembly implementation (fixed-length qualmer foundation).
"""
function _assemble_qualmer_graph(observations, config)
    _log_info(config, "Using quality-aware k-mer graph assembly strategy (fixed-length qualmer foundation)")

    # Convert observations to FASTQ records for quality processing
    fastq_records = _prepare_fastq_observations(observations)

    # Build qualmer graph using Phase 2 quality-aware algorithms
    mode = _graph_mode_symbol(config.graph_mode)
    # See _assemble_kmer_graph: Canonical graph_mode does not yet support correct
    # contig reconstruction (undirected canonical traversal is not orientation-
    # aware). Warn UNCONDITIONALLY and flag the result; do NOT hard-error.
    reconstruction_valid = config.graph_mode != Canonical
    if config.graph_mode == Canonical
        @warn "Canonical graph_mode does not yet support correct contig " *
              "reconstruction (undirected canonical traversal is not " *
              "orientation-aware); emitted contigs are INVALID. Use " *
              "graph_mode=DoubleStrand for correct assembly."
    end
    graph = Rhizomorph.build_qualmer_graph(fastq_records, config.k; mode = mode)

    # Apply Phase 2 graph algorithms if enabled
    if config.bubble_resolution
        # Note: This would use detect_bubbles_next adapted for qualmer graphs
        _log_info(config, "Bubble resolution enabled for qualmer graphs")
    end

    if config.repeat_resolution
        # Note: This would use resolve_repeats_next adapted for qualmer graphs
        _log_info(config, "Repeat resolution enabled for qualmer graphs")
    end

    # Find contigs by extracting maximal unitigs (non-branching paths), mirroring
    # the k-mer arm. The prior heaviest-path / iterative-Viterbi heuristics found
    # no substantial paths on branchy read graphs and fell to a placeholder walk
    # that emitted ~one short fragment per vertex, collapsing the largest contig
    # far below the genome length. find_contigs_next returns the vertex path
    # (ContigPath.vertices); per-base quality is then propagated by
    # _qualmer_path_to_consensus_fastq below.
    paths = [contig_path.vertices
             for contig_path in
                 Rhizomorph.find_contigs_next(graph; min_contig_length = config.k + 1)]

    # Convert paths to FASTQ records with quality propagation
    contig_records = FASTX.FASTQ.Record[]
    for (i, path) in enumerate(paths)
        if !isempty(path)
            contig_name = "qualmer_contig_$i"
            # Use consensus quality calculation for better accuracy
            fastq_record = _qualmer_path_to_consensus_fastq(path, graph, contig_name)
            if length(FASTX.sequence(fastq_record)) > config.k  # Only keep substantial contigs
                push!(contig_records, fastq_record)
            end
        end
    end

    # If no paths found, use probabilistic walks on qualmer graph
    if isempty(contig_records)
        _log_info(config, "No quality-aware paths found, using probabilistic walks")
        contig_records = _generate_fastq_contigs_from_qualmer_graph(graph, config)
    end

    # Convert FASTQ records to strings for backward compatibility with existing code
    contigs = [String(FASTX.sequence(record)) for record in contig_records]
    contig_names = [String(FASTX.identifier(record)) for record in contig_records]

    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "QualmerGraph",
        "graph_type" => "fixed_length",
        "vertex_type" => "quality_aware_kmers",
        "k" => config.k,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "quality_preserved" => true,  # Mark that quality information is preserved
        # false for Canonical graph_mode (undirected traversal is not orientation-aware).
        "reconstruction_valid" => reconstruction_valid,
        "num_fastq_contigs" => length(contig_records),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Mycelia.Dates.now())
    )

    # Add quality-specific statistics
    if !isempty(MetaGraphsNext.labels(graph))
        qualmer_stats = Rhizomorph.get_qualmer_statistics(graph)
        for (key, value) in qualmer_stats
            stats[string(key)] = value
        end
    end

    return AssemblyResult(contigs, contig_names; graph = graph,
        assembly_stats = stats, fastq_contigs = contig_records)
end

"""
Hybrid OLC assembly (placeholder for future implementation).
"""
function _assemble_hybrid_olc(observations, config)
    _log_info(config, "Using hybrid OLC assembly strategy")

    # For now, fall back to k-mer graph assembly
    # Future implementation would combine overlap-layout-consensus with string graphs
    @warn "Hybrid OLC not fully implemented, using k-mer graph assembly"
    return _assemble_kmer_graph(observations, config)
end

"""
Multi-k assembly (placeholder for future implementation).
"""
function _assemble_multi_k(observations, config)
    _log_info(config, "Using multi-k assembly strategy")

    # For now, fall back to single-k assembly
    # Future implementation would use multiple k-mer sizes and merge results
    @warn "Multi-k assembly not fully implemented, using single-k assembly"
    return _assemble_kmer_graph(observations, config)
end

"""
    _dedup_reverse_complements(contigs::Vector{String}) -> Vector{String}

Collapse contigs that are reverse complements of one another onto a single
canonical representative (the lexicographically-smaller of a sequence and its
reverse complement), preserving first-seen order.

A DoubleStrand assembly emits both orientations of each traversal, so the
contig set contains RC pairs. This is a lossless deduplication for downstream
consumers that treat a sequence and its reverse complement as equivalent: the
genome content is fully retained (every input contig maps to a canonical rep
that is either itself or its RC), while the redundant second orientation is
removed. Sequences that cannot be reverse-complemented as DNA (e.g. amino-acid
or general-string contigs) are passed through unchanged.
"""
function _dedup_reverse_complements(contigs::Vector{String})::Vector{String}
    canonical_reps = String[]
    seen = Set{String}()
    for contig in contigs
        canonical = _canonical_string(contig)
        if !(canonical in seen)
            push!(seen, canonical)
            push!(canonical_reps, canonical)
        end
    end
    return canonical_reps
end

"""
    _canonical_string(seq::AbstractString) -> String

Return `min(seq, reverse_complement(seq))` as an uppercase DNA string. If `seq`
is not valid DNA it is returned unchanged (no reverse complement is defined).
"""
function _canonical_string(seq::AbstractString)::String
    fwd = String(seq)
    rc = try
        String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(fwd)))
    catch e
        # Only a genuine invalid-DNA conversion failure means "not DNA — pass the
        # sequence through unchanged". BioSequences.LongDNA{4} throws an
        # ErrorException on unencodable characters (verified). Let anything else
        # (InterruptException, OutOfMemoryError, MethodError, ...) propagate rather
        # than silently masking a real fault.
        e isa ErrorException || rethrow()
        return fwd
    end
    return min(fwd, rc)
end

"""
Generate contigs using probabilistic walks when Eulerian paths are not available.
"""
function _generate_contigs_probabilistic(graph, config)
    # dedup_revcomp: prefer STRUCTURAL reverse-complement dedup during traversal
    # over the post-hoc whole-contig string dedup applied later. The structural
    # form marks each walked vertex's RC partner visited, so the reverse strand is
    # never independently traversed. This removes RC twins with OFFSET fragment
    # breakpoints that the string-level _dedup_reverse_complements misses (the
    # empirically-observed case where contig-level dedup halves the count but
    # leaves QUAST duplication at ~2.0).
    contig_paths = Rhizomorph.find_contigs_next(
        graph; min_contig_length = 1, rc_aware = config.dedup_revcomp)
    contigs = String[]

    for contig in contig_paths
        if length(contig.sequence) > config.k
            push!(contigs, string(contig.sequence))
        end
    end

    return contigs
end

"""
Polish a single contig using Viterbi error correction.
"""
function _polish_contig_viterbi(contig, graph, observations)
    # Create observation sequence from contig
    contig_kmers = [contig[i:(i + 30)] for i in 1:(length(contig) - 30)]  # Assuming k=31

    # Use Viterbi decoding for error correction
    try
        path = Mycelia.viterbi_decode_next(graph, contig_kmers)
        return string(Rhizomorph.path_to_sequence(path, graph))
    catch e
        @warn "Viterbi polishing failed for contig" error=e
        return contig  # Return original if polishing fails
    end
end

"""
Convert contigs to FASTA records for graph construction.
"""
function _contigs_to_records(contigs)
    records = FASTX.FASTA.Record[]
    for (i, contig) in enumerate(contigs)
        push!(records, FASTX.FASTA.Record("contig_$i", contig))
    end
    return records
end

"""
Calculate N-statistic (N50, N90, etc.) for contig lengths.
"""
function _calculate_n_statistic(sorted_lengths, threshold)
    total_length = sum(sorted_lengths)
    target_length = total_length * threshold

    cumulative_length = 0
    for length in sorted_lengths
        cumulative_length += length
        if cumulative_length >= target_length
            return length
        end
    end

    return 0
end

"""
    _calculate_l_statistic(sorted_lengths, threshold)

Calculate L-statistic (number of contigs needed to reach a given percentage of total assembly length).
For example, L50 is the number of contigs needed to reach 50% of the total assembly length.

# Arguments
- `sorted_lengths`: Vector of contig lengths sorted in descending order
- `threshold`: Fraction of total length (e.g., 0.5 for L50, 0.9 for L90)

# Returns
- `Int`: Number of contigs needed to reach the threshold
"""
function _calculate_l_statistic(sorted_lengths, threshold)
    total_length = sum(sorted_lengths)
    target_length = total_length * threshold

    cumulative_length = 0
    contig_count = 0

    for length in sorted_lengths
        cumulative_length += length
        contig_count += 1
        if cumulative_length >= target_length
            return contig_count
        end
    end

    return length(sorted_lengths)  # Return total number of contigs if threshold not reached
end

"""
Validate assembly against reference sequence (placeholder).
"""
function _validate_against_reference(contigs, reference)
    # Placeholder for reference-based validation
    # Would implement alignment-based metrics, coverage analysis, etc.
    return Dict{String, Any}(
        "reference_provided" => true,
        "reference_length" => length(reference)
    )
end

"""
Simplify string graph by removing unnecessary complexity.
"""
function _simplify_string_graph(graph)
    # Basic graph simplification - remove low-weight edges, merge linear paths
    # This is a placeholder for more sophisticated graph simplification
    return graph
end

# ============================================================================
# Qualmer Graph Assembly Helper Functions
# ============================================================================

"""
Convert observations to FASTQ records for quality processing.
"""
function _prepare_fastq_observations(observations)
    fastq_records = FASTX.FASTQ.Record[]

    for (i, obs) in enumerate(observations)
        # Handle different observation formats
        if obs isa FASTX.FASTQ.Record
            push!(fastq_records, obs)
        elseif obs isa FASTX.FASTA.Record
            # Convert FASTA to FASTQ with default quality (assume high quality)
            seq = FASTX.FASTA.sequence(obs)
            qual = repeat('I', length(seq))  # Quality score 40 (high quality)
            record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(obs), seq, qual)
            push!(fastq_records, record)
        elseif obs isa Tuple && length(obs) >= 1
            # Handle tuple format like (record, index) or (record, other_data)
            record = obs[1]
            if record isa FASTX.FASTQ.Record
                push!(fastq_records, record)
            elseif record isa FASTX.FASTA.Record
                # Convert FASTA to FASTQ with default quality
                seq = FASTX.FASTA.sequence(record)
                qual = repeat('I', length(seq))  # Quality score 40 (high quality)
                fastq_record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(record), seq, qual)
                push!(fastq_records, fastq_record)
            else
                @warn "Unsupported record type in tuple: $(typeof(record))"
            end
        else
            @warn "Unsupported observation type: $(typeof(obs))"
        end
    end

    return fastq_records
end

"""
Convert qualmer path to DNA sequence.
"""
function _qualmer_path_to_sequence(path, graph)
    if isempty(path)
        return ""
    end
    return string(Rhizomorph.path_to_sequence(path, graph))
end

"""
Enhanced qualmer path to FASTQ record conversion with quality propagation.
"""
function _qualmer_path_to_consensus_fastq(path, graph, contig_name::String)
    if isempty(path)
        return FASTX.FASTQ.Record("", "", UInt8[])
    end

    sequence = Rhizomorph.path_to_sequence(path, graph)
    quality_scores = UInt8[]

    for (idx, vertex) in enumerate(path)
        vertex_quality = _qualmer_vertex_quality_scores(graph[vertex])
        kmer_len = length(string(vertex))

        if idx == 1
            if length(vertex_quality) == kmer_len
                append!(quality_scores, vertex_quality)
            else
                append!(quality_scores, fill(UInt8(2), kmer_len))
            end
        else
            if isempty(vertex_quality)
                push!(quality_scores, UInt8(2))
            else
                push!(quality_scores, vertex_quality[end])
            end
        end
    end

    if length(quality_scores) != length(sequence)
        min_len = min(length(quality_scores), length(sequence))
        quality_scores = quality_scores[1:min_len]
        sequence = sequence[1:min_len]
    end

    return FASTX.FASTQ.Record(contig_name, sequence, quality_scores)
end

"""
Generate FASTQ contigs from qualmer graph using probabilistic walks when no paths are found.
"""
function _generate_fastq_contigs_from_qualmer_graph(graph, config)
    contig_records = FASTX.FASTQ.Record[]
    vertices = collect(MetaGraphsNext.labels(graph))

    if isempty(vertices)
        return contig_records
    end

    # Track visited vertices to avoid duplicates
    visited = Set()

    # Start evidence-weighted walks from high-confidence vertices
    min_quality = 1.0
    start_vertices = filter(v -> _qualmer_vertex_score(graph[v]) >= min_quality, vertices)

    # If no high-confidence vertices, use all vertices
    if isempty(start_vertices)
        start_vertices = vertices
    end

    contig_id = 1
    for start_vertex in start_vertices
        if start_vertex in visited
            continue
        end

        # Perform quality-weighted walk from this vertex
        path = _quality_weighted_walk(graph, start_vertex, config.k * 10)  # Reasonable max length

        if length(path) > 1  # Need at least 2 vertices to form a meaningful contig
            # Mark all vertices in path as visited
            union!(visited, path)

            # Convert path to FASTQ record
            contig_name = "qualmer_probabilistic_contig_$(contig_id)"
            fastq_record = _qualmer_path_to_consensus_fastq(path, graph, contig_name)

            # Only keep contigs longer than k
            if length(FASTX.sequence(fastq_record)) > config.k
                push!(contig_records, fastq_record)
                contig_id += 1
            end
        end
    end

    return contig_records
end

"""
Perform quality-weighted walk through qualmer graph.
"""
function _quality_weighted_walk(graph, start_vertex, max_length::Int = 1000)
    path = [start_vertex]
    current = start_vertex
    visited = Set([current])

    while length(path) < max_length
        # Get outgoing neighbors
        neighbors = collect(MetaGraphsNext.outneighbor_labels(graph, current))

        # Filter out visited vertices
        unvisited = filter(n -> !(n in visited), neighbors)

        if isempty(unvisited)
            break
        end

        # Choose next vertex based on evidence scores
        best_neighbor = nothing
        best_score = -1.0

        for neighbor in unvisited
            # Calculate score based on vertex evidence and edge evidence
            vertex_quality = _qualmer_vertex_score(graph[neighbor])

            # Check if edge exists and get its quality weight
            edge_quality = 1.0
            if haskey(graph, current, neighbor)
                edge_data = graph[current, neighbor]
                edge_quality = _qualmer_edge_score(edge_data)
            end

            # Combined score (vertex quality * edge quality)
            score = vertex_quality * edge_quality

            if score > best_score
                best_score = score
                best_neighbor = neighbor
            end
        end

        if best_neighbor === nothing
            break
        end

        push!(path, best_neighbor)
        push!(visited, best_neighbor)
        current = best_neighbor
    end

    return path
end

"""
Generate contigs from qualmer graph using probabilistic walks.
"""
function _generate_contigs_from_qualmer_graph(graph, config)
    contigs = String[]
    vertices = collect(MetaGraphsNext.labels(graph))

    # Start probabilistic walks from high-quality vertices
    for start_vertex in vertices
        vertex_data = graph[start_vertex]

        # Only start from high-quality k-mers
        if _qualmer_vertex_score(vertex_data) > 0.0
            # Simple random walk (placeholder for more sophisticated algorithm)
            path = [start_vertex]
            current = start_vertex

            for _ in 1:100  # Max walk length
                # Find outgoing edges
                outgoing = []
                for edge in MetaGraphsNext.edge_labels(graph)
                    src, dst = edge
                    if src == current && !(dst in path)  # Avoid cycles
                        push!(outgoing, dst)
                    end
                end

                if isempty(outgoing)
                    break
                end

                # Choose randomly weighted by quality
                current = rand(outgoing)
                push!(path, current)
            end

            if length(path) > 1
                sequence = _qualmer_path_to_sequence(path, graph)
                if length(sequence) > config.k
                    push!(contigs, sequence)
                end
            end
        end
    end

    return contigs
end

# ============================================================================
# Additional Assembly Methods for 6-Graph Hierarchy
# ============================================================================

"""
N-gram graph assembly implementation (fixed-length unicode character analysis).
"""
function _assemble_ngram_graph(observations, config)
    _log_info(config, "Using N-gram graph assembly strategy (fixed-length unicode analysis)")

    strings = [String(FASTX.FASTA.sequence(obs)) for obs in observations]
    graph = Rhizomorph.build_ngram_graph(strings, config.k)

    if config.bubble_resolution
        _simplify_ngram_graph(graph)
    end

    paths = Rhizomorph.find_eulerian_paths_next(graph)
    contigs = String[]
    for path in paths
        if length(path) > 1
            push!(contigs, string(Rhizomorph.path_to_sequence(path, graph)))
        end
    end

    if isempty(contigs)
        contig_paths = Rhizomorph.find_contigs_next(graph; min_contig_length = 1)
        contigs = [string(contig.sequence) for contig in contig_paths]
    end

    contig_names = ["ngram_contig_$i" for i in 1:length(contigs)]

    stats = Dict{String, Any}(
        "method" => "NgramGraph",
        "graph_type" => "fixed_length",
        "vertex_type" => "unicode_character_vectors",
        "k" => config.k,
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Mycelia.Dates.now())
    )

    return AssemblyResult(contigs, contig_names; graph = graph, assembly_stats = stats)
end

"""
BioSequence graph assembly implementation (variable-length simplified from k-mer graphs).
"""
function _assemble_biosequence_graph(observations, config)
    _log_info(config, "Using BioSequence graph assembly strategy (variable-length simplified from k-mer)")

    min_overlap = isnothing(config.min_overlap) ? 1 : config.min_overlap
    biosequence_graph = Rhizomorph.build_fasta_graph(observations; min_overlap = min_overlap)

    # Extract sequences from graph vertices
    contigs = String[]
    contig_names = String[]

    for (i, vertex_label) in enumerate(MetaGraphsNext.labels(biosequence_graph))
        sequence = string(vertex_label)
        push!(contigs, sequence)
        push!(contig_names, "biosequence_contig_$i")
    end

    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "BioSequenceGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "biosequences",
        "min_overlap" => min_overlap,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(biosequence_graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(biosequence_graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Mycelia.Dates.now())
    )

    return AssemblyResult(contigs, contig_names; graph = biosequence_graph, assembly_stats = stats)
end

"""
Quality-aware BioSequence graph assembly implementation (variable-length simplified from qualmer graphs).
"""
function _assemble_quality_biosequence_graph(observations, config)
    _log_info(config,
        "Using quality-aware BioSequence graph assembly strategy (variable-length simplified from qualmer)")

    # Convert observations to FASTQ records for quality processing
    fastq_records = _prepare_fastq_observations(observations)

    min_overlap = isnothing(config.min_overlap) ? 1 : config.min_overlap
    quality_biosequence_graph = Rhizomorph.build_fastq_graph(fastq_records; min_overlap = min_overlap)

    # Extract sequences from graph vertices
    contigs = String[]
    contig_names = String[]

    for (i, vertex_label) in enumerate(MetaGraphsNext.labels(quality_biosequence_graph))
        sequence = string(vertex_label)
        push!(contigs, sequence)
        push!(contig_names, "quality_biosequence_contig_$i")
    end

    # Assembly statistics
    stats = Dict{String, Any}(
        "method" => "QualityBioSequenceGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "quality_aware_biosequences",
        "min_overlap" => min_overlap,
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(quality_biosequence_graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(quality_biosequence_graph)),
        "num_input_sequences" => length(observations),
        "assembly_date" => string(Mycelia.Dates.now())
    )

    return AssemblyResult(contigs, contig_names; graph = quality_biosequence_graph, assembly_stats = stats)
end

"""
Simplify N-gram graph by removing unnecessary complexity.
"""
function _simplify_ngram_graph(graph)
    # Basic graph simplification - remove low-weight edges, merge linear paths
    # This is a placeholder for more sophisticated graph simplification
    return graph
end

"""
Convert N-gram graph to variable-length string graph.
"""
function _simplify_ngram_to_string_graph(ngram_graph)
    # Placeholder for N-gram to string graph conversion
    # This would implement path collapsing and simplification
    return ngram_graph
end
