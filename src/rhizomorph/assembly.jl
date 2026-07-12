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
    TokenGraph       # Token graph assembly (SentencePiece/word-token sequences, SingleStrand only)

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
    # (canonical) orientation of the sequence that was actually assembled. Wired into
    # BOTH the k-mer and qualmer arms (td-47di): structural rc_aware traversal in
    # find_contigs_next plus a post-hoc canonical collapse. Defaults ON for the
    # corrector=:iterative route, OFF otherwise (see constructor).
    dedup_revcomp::Bool                     # Collapse RC-pair contigs to one canonical rep
    compact_unitigs::Bool                   # Populate simplified_graph via linear-chain compaction
    memory_profile::Symbol                  # build_kmer_graph evidence footprint (:full|:lightweight|:ultralight|...)

    # Optional read-correction front-end (opt-in; default :none preserves today's
    # single-k-from-uncorrected-reads behavior byte-for-byte).
    corrector::Symbol                       # :none (default) or :iterative (mycelia_iterative_assemble)
    skip_solid::Bool                        # When corrector=:iterative, skip already-solid/clean reads

    # Corrector TIER selector (td-fuo8). Only consulted when corrector=:iterative.
    # :scalable (DEFAULT) — coarse k-ladder + low iteration cap + skip-solid +
    #   hard-read gating + soft-EM edge memory (record-only v1). Built for
    #   real-scale inputs.
    # :exhaustive — maximum-sensitivity EXACT-ML tier: prime-by-prime k-walk,
    #   10 iterations/k, exact UNBOUNDED (typemax) Viterbi beam, no skip, no
    #   hard-window, no soft-EM. This is NOT a reproduction of the prior corrector
    #   default: master's corrector route used the size-aware auto-beam
    #   (beam_width=nothing, bounded on large reads), so forcing an exact unbounded
    #   beam here can OOM above the auto-beam threshold (the td-63qy regime). Use
    #   :exhaustive for small-scale / high-sensitivity inputs; use :scalable at
    #   scale. None of the :scalable-gated behavior perturbs :exhaustive.
    strategy::Symbol

    # Optional pre-tokenized input (opt-in; default `nothing` preserves the
    # read-driven pipeline). When supplied, assemble_genome routes to the token
    # graph and IGNORES `reads` — the tokens ARE the input. SingleStrand only
    # (reverse complement is undefined for general string tokens).
    token_sequences::Union{Nothing, Vector{Vector{String}}}

    # Linear-time defragmentation of the (qualmer) assembly graph BEFORE contig
    # extraction (td-969e): coverage-1 dead-end tip clipping + guarded low-coverage
    # bubble collapse. Three-state sentinel so the DEFAULT can differ by context
    # (see _qualmer_graph_to_assembly): `nothing` = context default (ON for the
    # iterative corrector's re-assembly, OFF for a plain qualmer assembly), `true`
    # / `false` = force. Removes only unambiguous errors; never collapses a
    # data-supported variant (td-h6w9), so it cannot lose real sequence.
    graph_cleanup::Union{Bool, Nothing}

    # Persistent directory for the Stage-1 corrector handoff (hybrid-OLC route,
    # td-ohob). Only consulted on the corrector=:iterative path. `nothing`
    # (default) = ephemeral: Stage-1 moves the corrected FASTQ to a `tempname()`
    # path the native re-assembly deletes after use, preserving the historical
    # no-stray-file behavior. When set, the corrected FASTQ is persisted to the
    # FIXED path `joinpath(output_dir, "corrected.fastq")` and left in place so an
    # external OLC assembler (Stage-2 route (a)) can consume it; give each Stage-1
    # run its own output_dir since the basename is fixed.
    output_dir::Union{String, Nothing}

    # Sequencing-technology ERROR PROFILE driving the corrector's indel-aware decode
    # (td-9q84 / 4a). Only consulted when corrector=:iterative. The profile maps to
    # per-base indel fractions (see `Mycelia.indel_error_profile`); indel-prone
    # technologies (:nanopore, :pacbio) enable the pair-HMM gap moves, while the
    # DEFAULT :illumina (and :ultima) carry ~0 indel fractions and stay on the
    # substitution-only path byte-for-byte. Wiring the profile — not read length or
    # a hardcoded tech — keeps "correct with the profile you'd simulate".
    sequencing_tech::Symbol

    # Stage-2 layout selector (hybrid-OLC route (a), td-yymj). `:native` (default)
    # keeps the current behavior — for corrector=:iterative that is the graph-HMM
    # re-assembly. `:olc` diverts corrected reads to an EXTERNAL overlap-layout-
    # consensus assembler so contiguity is inherited from a tool tuned for it.
    # `:olc` requires corrector=:iterative (an OLC layout with no correction is just
    # the plain external assembler, reachable via the wrapper directly).
    layout::Symbol

    # Which external assembler the `:olc` layout dispatches to. `:auto` (default)
    # picks by read type (sequencing_tech): short-read techs (:illumina, :ultima)
    # → :megahit, long-read techs (:nanopore, :pacbio) → :flye. Explicit tools:
    # short-read :megahit/:metaspades (td-yymj), long-read :flye/:metaflye/:canu/
    # :hifiasm (td-wvto). An explicit tool must match the corrected read type
    # (hifiasm is PacBio-HiFi only; canu also needs an olc_options genome_size).
    olc_tool::Symbol

    # Pass-through options forwarded to the external-assembler wrapper (e.g.
    # (; k_list = "21,33", threads = 8)). Empty by default; keys must match the
    # chosen wrapper's keyword arguments. Route-managed keys (fastq1/fastq2/fastq/
    # outdir/out_dir/executor) are RESERVED and rejected at construction — they are
    # set by the route and cannot be overridden here.
    olc_options::NamedTuple

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
            dedup_revcomp::Union{Bool, Nothing} = nothing,
            compact_unitigs::Bool = false,
            memory_profile::Symbol = :full,
            corrector::Symbol = :none,
            skip_solid::Union{Bool, Nothing} = nothing,
            strategy::Union{Symbol, Nothing} = nothing,
            token_sequences::Union{Nothing, Vector{Vector{String}}} = nothing,
            graph_cleanup::Union{Bool, Nothing} = nothing,
            output_dir::Union{String, Nothing} = nothing,
            sequencing_tech::Symbol = :illumina,
            layout::Symbol = :native,
            olc_tool::Symbol = :auto,
            olc_options::NamedTuple = (;)
    )
        # Sentinel `nothing` defaults let us DETECT whether the caller set
        # strategy/skip_solid explicitly (FIX 4) so we can warn on the silent
        # :scalable default and on a skip_solid the tier will override, without
        # changing the stored field types. Resolve to concrete values here.
        strategy_explicit = strategy !== nothing
        skip_solid_explicit = skip_solid !== nothing
        effective_strategy = strategy === nothing ? :scalable : strategy
        effective_skip_solid = skip_solid === nothing ? false : skip_solid
        # Default decision (td-47di): a DoubleStrand assembly emits BOTH strands of
        # every contig, so on that path canonical (RC-deduped) output — one contig
        # per genomic locus rather than a forward/reverse twin pair — is the sensible
        # default. The iterative corrector always re-assembles DNA on a DoubleStrand
        # graph, so dedup_revcomp defaults ON for corrector=:iterative. It stays OFF
        # for the plain (corrector=:none) k-mer/qualmer route, preserving the existing
        # both-strands behavior that current tests assert. Either default can be
        # overridden by passing an explicit dedup_revcomp=true/false.
        effective_dedup_revcomp = dedup_revcomp === nothing ?
                                  (corrector == :iterative) : dedup_revcomp
        # Validation: Must specify exactly one of k or min_overlap
        if k === nothing && min_overlap === nothing
            k = 31  # Default to k-mer mode with k=31
        elseif k !== nothing && min_overlap !== nothing
            error("Cannot specify both k ($(k)) and min_overlap ($(min_overlap)). Choose one approach.")
        end

        # Validation: Check strand compatibility with sequence types.
        # DoubleStrand and Canonical both require a defined reverse complement,
        # so both are rejected for amino acids and general strings.
        if sequence_type <: BioSequences.LongAA &&
           (graph_mode == DoubleStrand || graph_mode == Canonical)
            error("Amino acid sequences can only use SingleStrand mode (reverse complement undefined for proteins)")
        end
        if sequence_type == String &&
           (graph_mode == DoubleStrand || graph_mode == Canonical)
            error("String sequences can only use SingleStrand mode (reverse complement undefined for general strings)")
        end

        # Validation: token-graph input is SingleStrand only. Tokens are opaque
        # strings (SentencePiece pieces / words) with no defined reverse
        # complement, so DoubleStrand/Canonical are meaningless for them.
        if token_sequences !== nothing && graph_mode != SingleStrand
            error("token_sequences requires SingleStrand mode (reverse complement undefined for string tokens)")
        end

        # Validation: memory_profile must be one recognized by build_kmer_graph
        _valid_memory_profiles = (
            :full, :lightweight, :ultralight, :lightweight_quality, :ultralight_quality)
        if !(memory_profile in _valid_memory_profiles)
            error("memory_profile must be one of $(_valid_memory_profiles), got :$(memory_profile)")
        end

        # Validation: corrector front-end must be a recognized mode
        _valid_correctors = (:none, :iterative)
        if !(corrector in _valid_correctors)
            error("corrector must be one of $(_valid_correctors), got :$(corrector)")
        end

        # Validation: corrector tier (only meaningful when corrector=:iterative,
        # but validate unconditionally so a typo is caught at construction time).
        _valid_strategies = (:scalable, :exhaustive)
        if !(effective_strategy in _valid_strategies)
            error("strategy must be one of $(_valid_strategies), got :$(effective_strategy)")
        end

        # Validation: sequencing_tech must be a recognized error profile (validate
        # unconditionally so a typo is caught at construction, even for
        # corrector=:none where it is not consulted). Derived from the OLC taxonomy
        # so the short/long partition `:auto` resolution depends on stays a single
        # source of truth — adding a tech in one place can't silently mislabel it.
        _seq_tax = _olc_taxonomy()
        _valid_seq_techs = (_seq_tax.short_read_techs..., _seq_tax.long_read_techs...)
        if !(sequencing_tech in _valid_seq_techs)
            error("sequencing_tech must be one of $(_valid_seq_techs), got :$(sequencing_tech)")
        end

        # FIX 4 (silent-default + skip_solid-override provenance). corrector=:iterative
        # now DEFAULTS to strategy=:scalable, which is a materially different engine
        # than the prior corrector route (coarse ladder / low iteration cap /
        # skip-solid / doublestrand graph mode / hard-read gate). Announce the default
        # once so a caller relying on old behavior is not silently switched.
        if corrector == :iterative && !strategy_explicit
            @warn "corrector=:iterative now defaults to strategy=:scalable — a coarse " *
                  "k-ladder + low iteration cap with skip-solid, doublestrand graph mode, " *
                  "and hard-read gating enabled. This is NOT the prior corrector " *
                  "behavior; pass strategy=:exhaustive for the maximum-sensitivity " *
                  "exact-ML engine, or strategy=:scalable to silence this warning." maxlog = 1
        end
        # The tier fully determines skip_solid (see _corrector_strategy_knobs):
        # :scalable ⇒ true, :exhaustive ⇒ false. An explicit config.skip_solid that
        # disagrees is silently overridden by the corrector route; warn so the
        # override is visible.
        if corrector == :iterative && skip_solid_explicit
            _tier_skip_solid = effective_strategy == :scalable
            if effective_skip_solid != _tier_skip_solid
                @warn "strategy=:$(effective_strategy) forces skip_solid=$(_tier_skip_solid); " *
                      "your explicit skip_solid=$(effective_skip_solid) is overridden by the " *
                      "corrector tier. assembly_stats records both requested and effective." maxlog = 1
            end
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

        # dedup_revcomp is now wired into the qualmer arm too (td-47di): the
        # quality-aware find_contigs_next calls in _qualmer_graph_to_assembly pass
        # rc_aware=config.dedup_revcomp, matching the k-mer arm, so it is NO LONGER
        # ignored on the quality path. The remaining efficiency modes
        # (compact_unitigs / memory_profile) are still only implemented on the
        # non-quality k-mer path (_assemble_kmer_graph); FASTQ input auto-sets
        # use_quality_scores=true, which dispatches to the qualmer arm and silently
        # ignores them. Warn so a caller opting in on quality data is not left
        # believing those flags took effect.
        if use_quality_scores &&
           (compact_unitigs || memory_profile != :full)
            @warn "Efficiency modes (compact_unitigs / memory_profile) are only " *
                  "implemented on the non-quality k-mer path and are ignored on " *
                  "the quality/qualmer path. Set use_quality_scores=false to use them."
        end

        # Validation: output_dir, when supplied, must be a non-empty path. An empty
        # string is a footgun — it takes the persist branch (only `nothing` is
        # ephemeral) and `joinpath("", "corrected.fastq")` silently writes into the
        # process cwd, the exact stray-file behavior the ephemeral path avoids.
        if output_dir !== nothing && isempty(output_dir)
            error("output_dir must be a non-empty directory path (got \"\"); " *
                  "pass nothing for the ephemeral tempfile behavior.")
        end
        # output_dir only affects the corrector=:iterative Stage-1 handoff; warn
        # rather than silently ignore it on other paths (mirrors the skip_solid /
        # efficiency-mode warnings above).
        if output_dir !== nothing && corrector != :iterative
            @warn "output_dir is only used by corrector=:iterative (the Stage-1 " *
                  "corrected-FASTQ handoff) and is ignored for corrector=:$(corrector)."
        end

        # Validation: Stage-2 layout selector (hybrid-OLC route (a), td-yymj).
        if !(layout in (:native, :olc))
            error("layout must be :native or :olc, got :$(layout)")
        end
        # Read-type <-> tool taxonomy from the SINGLE source of truth (_olc_taxonomy).
        # sequencing_tech is the read-type proxy: short-read techs run short-read
        # layout assemblers, long-read techs run long-read assemblers (design doc
        # "Read-type mapping"). Both arms are wired (short-read td-yymj, long-read
        # td-wvto).
        _tax = _olc_taxonomy()
        _valid_olc_tools = (:auto, _tax.short_read_tools..., _tax.long_read_tools...)
        if layout == :olc
            # An OLC layout with no correction is just the plain external assembler,
            # reachable via the wrapper directly — the route exists to compose Stage-1
            # correction with external layout, so it requires the corrector.
            if corrector != :iterative
                error("layout=:olc requires corrector=:iterative (the hybrid-OLC " *
                      "route hands Stage-1-corrected reads to an external assembler); " *
                      "got corrector=:$(corrector).")
            end
            if !(olc_tool in _valid_olc_tools)
                error("olc_tool must be one of $(_valid_olc_tools), got :$(olc_tool)")
            end
            is_short = sequencing_tech in _tax.short_read_techs
            if olc_tool == :auto
                # :auto always resolves to a WIRED tool (short-read tech → :megahit,
                # long-read tech → :flye), so no rejection is needed here.
                nothing
            elseif olc_tool in _tax.short_read_tools
                is_short || error("olc_tool=:$(olc_tool) is a short-read assembler but " *
                      "sequencing_tech=:$(sequencing_tech) is a long-read tech; pair a short-read " *
                      "tool with a short-read tech.")
            elseif olc_tool in _tax.long_read_tools
                sequencing_tech in _tax.long_read_techs ||
                    error("olc_tool=:$(olc_tool) is a " *
                          "long-read assembler but sequencing_tech=:$(sequencing_tech) is a short-read " *
                          "tech; pair a long-read tool with :nanopore or :pacbio.")
                # hifiasm assembles PacBio HiFi reads specifically — not Nanopore.
                if olc_tool == :hifiasm && sequencing_tech != :pacbio
                    error("olc_tool=:hifiasm assembles PacBio HiFi reads; " *
                          "sequencing_tech=:$(sequencing_tech) is not HiFi. Use :flye or :canu " *
                          "for :nanopore.")
                end
                # canu requires an estimated genome size (the wrapper has no default).
                if olc_tool == :canu && !haskey(olc_options, :genome_size)
                    error("olc_tool=:canu requires an estimated genome size; pass it via " *
                          "olc_options=(; genome_size=\"5m\").")
                end
            end
            # Reserved-key guard: olc_options is splatted into the wrapper, so a key
            # that collides with a route-managed keyword would redirect the assembler
            # (wrong output dir → broken cleanup + megahit's recursive rm; wrong fastq
            # → bypasses the Stage-1 corrected handoff). Reject those keys up front.
            _reserved_olc_keys = (:fastq1, :fastq2, :fastq, :outdir, :out_dir,
                :executor, :read_type)
            for _k in keys(olc_options)
                _k in _reserved_olc_keys && error("olc_options key :$(_k) is managed by " *
                      "the hybrid-OLC route and cannot be overridden; drop it from olc_options.")
            end
            # metaSPAdes targets PAIRED-END libraries; the hybrid-OLC route emits a
            # SINGLE corrected FASTQ (single-end), which metaSPAdes rejects at runtime.
            # Keep it selectable (paired-end support is a follow-on) but warn loudly —
            # :megahit is the single-end-robust short-read default.
            if olc_tool == :metaspades
                @warn "olc_tool=:metaspades expects paired-end reads, but the hybrid-OLC " *
                      "route feeds a single corrected FASTQ (single-end); metaSPAdes will " *
                      "likely fail at runtime. Use :megahit for the single-end short-read arm."
            end
        elseif olc_tool != :auto || !isempty(olc_options)
            # Discoverability: olc_tool/olc_options are inert unless layout=:olc
            # (mirrors the output_dir / skip_solid / efficiency-mode warnings above).
            @warn "olc_tool / olc_options are only used when layout=:olc and are " *
                  "ignored for layout=:$(layout)."
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
            effective_dedup_revcomp,
            compact_unitigs,
            memory_profile,
            corrector,
            effective_skip_solid,
            effective_strategy,
            token_sequences,
            graph_cleanup,
            output_dir,
            sequencing_tech,
            layout,
            olc_tool,
            olc_options
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
- `corrector`: `:none` (default, single-k from uncorrected reads) or `:iterative`
  (route through the iterative maximum-likelihood read corrector before assembly)
- `sequencing_tech`: error profile driving the corrector's indel-aware decode; only
  consulted when `corrector=:iterative`. Default `:illumina` = substitution-only
  correction (byte-identical to the pre-wiring corrector). Set `:nanopore` or
  `:pacbio` to enable indel-aware pair-HMM correction of long / indel-prone reads;
  `:ultima` is substitution-only. Correcting nanopore/pacbio reads under the default
  `:illumina` leaves their indels uncorrected.
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
        # Stage-2 layout selector (hybrid-OLC route (a), td-yymj). :olc diverts the
        # Stage-1-corrected reads to an external OLC assembler; :native (default)
        # keeps the graph-HMM re-assembly. The constructor guarantees layout=:olc
        # implies corrector=:iterative, so this branch is only reachable here.
        if config.layout == :olc
            return _assemble_hybrid_olc(reads, config)
        end
        return _assemble_with_iterative_corrector(reads, config)
    end

    # Token-graph route. When the caller supplies pre-tokenized sequences, the
    # tokens ARE the input and `reads` is ignored — token assembly is driven by
    # config.token_sequences, not by loaded read observations. SingleStrand only
    # (enforced in the AssemblyConfig constructor).
    if config.token_sequences !== nothing
        _log_info(config, "Routing to token graph assembly (token_sequences supplied)")
        return _assemble_token_graph(config.token_sequences, config)
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
                            write(writer,
                                FASTX.FASTQ.Record(
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
                    write(writer,
                        FASTX.FASTQ.Record(
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
    _corrector_strategy_knobs(strategy::Symbol) -> NamedTuple

Map a corrector `strategy` tier onto the concrete engine knobs threaded into
`Mycelia.mycelia_iterative_assemble` (td-fuo8). Pure + side-effect-free so the
routing can be unit-tested without running the (slow) corrector.

- `:scalable` — coarse LoRMA-style 3-rung k-ladder, a low (2) iteration cap,
  skip-solid volume reduction, Stage 0 cheap k-mer-spectrum correction (td-bjnt,
  linear single-substitution fix before the decode), hard-read gating (Stage 3,
  now narrowed to bubble/repeat vertices only), soft-EM edge memory (Stage 2,
  record-only v1), the size-aware auto-beam (`beam_width=nothing`), and
  `graph_mode=:doublestrand` (td-nt69 — canonical was over-constrained by the skip
  machinery and was the cause of the quality gap; the skip/classification is
  coverage-based and mode-agnostic).
- `:exhaustive` — maximum-sensitivity EXACT-ML tier: prime-by-prime k-walk
  (`n_k_rungs=nothing`), 10 iterations/k, no skip, no hard-window, no soft-EM,
  `graph_mode=nothing` (derive from `config.graph_mode`, byte-identical
  passthrough), and an exact UNBOUNDED (`typemax(Int)`) Viterbi beam. This is NOT a reproduction
  of the prior corrector default: master's corrector route used the size-aware
  auto-beam (`beam_width=nothing`, bounded on large reads), so forcing an exact
  unbounded beam here reintroduces the td-63qy OOM ABOVE the auto-beam threshold.
  Intended for SMALL-SCALE / high-sensitivity inputs; the exact beam can OOM on
  very large reads — use `:scalable` at scale.
"""
function _corrector_strategy_knobs(strategy::Symbol)
    if strategy == :scalable
        return (
            n_k_rungs = 3,
            max_iterations_per_k = 2,
            skip_solid = true,
            hard_window = true,
            # Stage 3c (td-nn6l): decode each hard read WINDOW-BY-WINDOW (only the
            # boundary-constrained hard sub-window(s), <=500 bp) instead of
            # whole-read, bounding a long read's decode to O(window) not
            # O(read_length) — the #375 long-read super-linear term.
            windowed_decode = true,
            soft_em = true,
            cheap_correct = true,  # Stage 0 linear k-mer-spectrum correction (td-bjnt)
            beam_width = nothing,  # size-aware auto-beam (bounded on huge reads)
            # graph_mode=:doublestrand (td-nt69): forcing :canonical was THE cause of
            # the :scalable quality gap. On a controlled 1kb fixture, flipping ONLY
            # canonical→doublestrand (all other knobs held) collapses 197→16 contigs
            # and lifts N50 43→891 — the historical near-complete regime. Canonical
            # was forced only because the skip machinery was (wrongly) wired to
            # require it; the k-mer classification is coverage-based and mode-agnostic
            # (each vertex + its RC are separate doublestrand vertices, still
            # separable by coverage), so skip_solid + hard_window work under
            # :doublestrand too — the naive contig path stays valid.
            graph_mode = :doublestrand,
            # Indel-decode bounds (td-9q84 / 4a), consulted ONLY when the error
            # profile enables indels (nanopore/pacbio). :scalable keeps the pair-HMM
            # tractable at scale: a bounded diagonal band on the net gap + small run
            # caps (the bounded Bellman-Ford relaxation replacing the O(V³)
            # Floyd-Warshall). Ignored on substitution-only profiles (:illumina),
            # whose decode threads no indel params at all.
            band_width = 16,
            deletion_max_run = 3,
            max_insertion_run = 3
        )
    elseif strategy == :exhaustive
        return (
            n_k_rungs = nothing,    # prime-by-prime hyper-sensitive walk
            max_iterations_per_k = 10,
            skip_solid = false,
            hard_window = false,
            windowed_decode = false,
            soft_em = false,
            cheap_correct = false,  # exact-ML tier: no cheap pre-correction
            beam_width = typemax(Int),  # exact ML decode
            # nothing ⇒ derive the corrector graph_mode from config.graph_mode below,
            # keeping :exhaustive a byte-identical passthrough of prior behavior.
            graph_mode = nothing,
            # Indel-decode bounds (td-9q84 / 4a): maximum-sensitivity tier ⇒ UNBOUNDED
            # band (`band_width=nothing`, the exact/oracle setting) + larger run caps.
            # Consulted only when the error profile enables indels.
            band_width = nothing,
            deletion_max_run = 10,
            max_insertion_run = 10
        )
    else
        error("unknown corrector strategy :$(strategy); expected :scalable or :exhaustive")
    end
end

"""
    _run_stage1_correction(reads, config::AssemblyConfig)
        -> (; corrected_reads, corrected_fastq, result_dict, knobs, max_k,
              ephemeral, indel_params)

Run Stage-1 correction (`Mycelia.mycelia_iterative_assemble`, the iterative +
skip-solid maximum-likelihood corrector) and PERSIST the corrected FASTQ to a
caller-owned location so it outlives the corrector's internal temp dirs.

This is the reusable Stage-1 correction half of the corrector route. The native
re-assembly path materializes corrected reads, while disk-backed OLC and Stage-2
callers pass `materialize_corrected_reads = false` to retain only the validated
read count and persistent FASTQ hand-off.
The corrected FASTQ produced by the corrector normally lives in an `mktempdir`
that is deleted here on return; this helper moves it OUT of that doomed directory
to `corrected_fastq` before cleanup, so an external OLC assembler can consume it.

Destination: `config.output_dir === nothing` ⇒ a fresh `tempname()` path the
CALLER owns and must delete (`ephemeral == true`; preserves the native path's
historical no-stray-file behavior). Otherwise the FIXED path
`joinpath(config.output_dir, "corrected.fastq")`, left in place for the Stage-2
handoff (`ephemeral == false`). Because the basename is fixed, a caller reusing
one `output_dir` across runs overwrites the prior corrected set — give each
Stage-1 run its own `output_dir`.

Returns a NamedTuple: optionally materialized `corrected_reads`, the validated
`corrected_read_count`, persistent FASTQ path (`corrected_fastq`), whether the
caller owns cleanup (`ephemeral`), and the corrector's `result_dict`, tier
`knobs`, and `max_k`. Callers should consume the returned `corrected_fastq`, not
`result_dict[:metadata][:final_fastq_file]` (both point at the same persisted
path after this returns). Fails loud (never returns) on a missing/absent
corrected FASTQ or a 0-read correction, guarding both callers.
"""
function _run_stage1_correction(
        reads,
        config::AssemblyConfig;
        materialize_corrected_reads::Bool = true
)
    _log_info(config,
        "Routing assembly through iterative corrector " *
        "(corrector=:iterative, strategy=:$(config.strategy))")

    # Discoverability nudge (PR #408 review, FIX 3): the corrector defaults to the
    # substitution-only :illumina profile, so a user correcting long / indel-prone
    # reads (nanopore/pacbio) without setting sequencing_tech silently gets NO indel
    # correction. Surface it once (@info, not @warn — :illumina is a valid, common
    # choice) so the off-by-default indel path is discoverable.
    if config.sequencing_tech == :illumina
        @info "corrector=:iterative is running with sequencing_tech=:illumina " *
              "(substitution-only correction). For long / indel-prone reads set " *
              "sequencing_tech=:nanopore or :pacbio to enable indel-aware correction." maxlog = 1
    end

    # The corrector is k-mer based; a min_overlap-only config has no k, so the
    # k-progression floors at 13 and min_overlap is not used — surface that rather
    # than silently drop the OLC intent (review Important #3).
    if config.k === nothing
        @warn "corrector=:iterative is k-mer based; min_overlap is ignored and the k-progression floors at 13."
    end
    max_k = config.k === nothing ? 13 : max(config.k, 13)

    # Resolve the persistent corrected-FASTQ destination up front (td-ohob). When
    # config.output_dir is nothing this is a fresh tempname() path the CALLER owns
    # and must clean up (native re-assembly deletes it, preserving no-stray-file
    # behavior); when set, the corrected reads persist at a FIXED basename for an
    # external OLC assembler (hybrid route) to consume and are left in place.
    ephemeral = config.output_dir === nothing
    persistent_fastq = ephemeral ?
                       tempname() * ".fastq" :
                       joinpath(mkpath(config.output_dir), "corrected.fastq")

    input_dir = mktempdir()
    corrector_output_dir = mktempdir()
    # Did THIS call actually move the corrected FASTQ into persistent_fastq? The
    # catch cleanup must delete only a file this call created — otherwise an early
    # failure (before the mv) with a reused config.output_dir would delete a PRIOR
    # run's good corrected.fastq (a data-loss footgun surfaced in review).
    moved_here = false
    try
        temp_fastq = joinpath(input_dir, "corrector_input.fastq")
        _write_reads_to_fastq(reads, temp_fastq)

        # Fork the engine knobs by tier (td-fuo8). :exhaustive is the
        # maximum-sensitivity exact-ML engine (n_k_rungs=nothing / 10 iters /
        # exact UNBOUNDED beam / no skip / no hard-window / no soft-EM) — NOT a
        # byte-for-byte reproduction of the prior corrector default (which used the
        # bounded auto-beam), so it can OOM on very large reads; :scalable opts into
        # the coarse ladder + low iteration cap + skip-solid + hard-read gate +
        # soft-EM. The per-tier knobs live in _corrector_strategy_knobs so the
        # routing is unit-testable without running the corrector.
        knobs = _corrector_strategy_knobs(config.strategy)
        # td-nt69: :scalable now runs the corrector on a :doublestrand graph
        # (`knobs.graph_mode == :doublestrand`). Forcing :canonical was THE cause of
        # the quality gap — the skip machinery was over-constrained to require it,
        # but skip-solid + hard-window classification is coverage-based and
        # mode-agnostic, so both gates stay ACTIVE under :doublestrand and the naive
        # contig path stays valid. :exhaustive threads `knobs.graph_mode === nothing`
        # ⇒ derive from config.graph_mode, a byte-identical passthrough.
        corrector_graph_mode = knobs.graph_mode === nothing ?
                               _graph_mode_symbol(config.graph_mode) : knobs.graph_mode
        # Indel-aware correction wiring (td-9q84 / 4a): map the sequencing-tech error
        # profile to indel fractions and gate on the ABSOLUTE indel rate
        # (base_error_rate × summed fractions). Only indel-prone profiles (:nanopore,
        # :pacbio) build a non-nothing IndelDecodeParams; the DEFAULT :illumina (and
        # :ultima) resolve to `nothing`, so the corrector threads NO indel params and
        # runs the substitution decode byte-for-byte (oracle preservation). The
        # base_error_rate (threaded into ViterbiCorrectionConfig.error_rate to scale
        # the gap masses), gap-open fractions, and extend probabilities come from the
        # profile; the run caps + band from the tier knobs above.
        indel_params = if Mycelia.profile_enables_indels(config.sequencing_tech)
            profile = Mycelia.indel_error_profile(config.sequencing_tech)
            Mycelia.IndelDecodeParams(
                profile.base_error_rate,
                profile.insertion_fraction,
                profile.deletion_fraction,
                profile.insertion_extend_probability,
                profile.deletion_extend_probability,
                knobs.deletion_max_run,
                knobs.max_insertion_run,
                knobs.band_width
            )
        else
            nothing
        end
        result_dict = Mycelia.mycelia_iterative_assemble(
            temp_fastq;
            max_k = max_k,
            skip_solid = knobs.skip_solid,
            graph_mode = corrector_graph_mode,
            n_k_rungs = knobs.n_k_rungs,
            max_iterations_per_k = knobs.max_iterations_per_k,
            hard_window = knobs.hard_window,
            windowed_decode = knobs.windowed_decode,
            soft_em = knobs.soft_em,
            cheap_correct = knobs.cheap_correct,
            beam_width = knobs.beam_width,
            indel_params = indel_params,
            verbose = false,
            enable_checkpointing = false,
            output_dir = corrector_output_dir
        )

        # mycelia_iterative_assemble is a read CORRECTOR: its :final_assembly is the
        # corrected READS, not an assembly. The caller re-assembles them (or hands
        # them to an external assembler) — here we just materialize + persist them.
        corrected_fastq = get(result_dict[:metadata], :final_fastq_file, nothing)
        if corrected_fastq === nothing
            error("iterative corrector metadata is missing the :final_fastq_file key; " *
                  "keys=$(collect(keys(result_dict[:metadata])))")
        elseif !isfile(corrected_fastq)
            error("iterative corrector :final_fastq_file points at a nonexistent path: " *
                  "$(corrected_fastq) (output_dir=$(corrector_output_dir))")
        end
        # Persist the corrected FASTQ OUT of the doomed corrector_output_dir before
        # the finally block deletes it, so the corrected reads outlive corrector
        # cleanup (the hybrid-OLC gap, td-ohob). `mv` (not cp) because
        # corrector_output_dir is removed immediately after — no need to double-store.
        Base.mv(corrected_fastq, persistent_fastq; force = true)
        moved_here = true
        # Keep the returned metadata honest: the corrector set :final_fastq_file to
        # the now-moved original path (inside the deleted temp dir). Repoint it at
        # persistent_fastq so a downstream consumer that reads the metadata key sees
        # a live path, not a dangling one. (The native tail never reads this key, so
        # this does not perturb the byte-identical re-assembly.)
        result_dict[:metadata][:final_fastq_file] = persistent_fastq
        # Eager `collect` inside a do-block: materializes ALL corrected reads into
        # memory AND closes the stream (no leaked fd across many assemblies, and
        # closes on a malformed-FASTQ throw too).
        corrected_reads, n_corrected = open(FASTX.FASTQ.Reader, persistent_fastq) do reader
            if materialize_corrected_reads
                records = collect(reader)
                return records, length(records)
            end
            return nothing, count(_ -> true, reader)
        end
        # A corrector that silently ate every read (empty/header-only FASTQ) would
        # otherwise flow into assemble_genome([]) → a 0-contig "successful" assembly.
        # Fail loud instead of handing downstream a silently-empty assembly.
        if n_corrected == 0
            error("iterative corrector produced 0 corrected reads (from $(persistent_fastq)); " *
                  "refusing to return an empty assembly")
        end
        return (; corrected_reads, corrected_read_count = n_corrected,
            corrected_fastq = persistent_fastq,
            result_dict, knobs, max_k, ephemeral, indel_params)
    catch
        # Ownership of the persistent FASTQ only transfers to the caller on a clean
        # return. Delete ONLY a file this call created (moved_here) — never a
        # pre-existing corrected.fastq a reused output_dir may already hold.
        moved_here && rm(persistent_fastq; force = true)
        rethrow()
    finally
        # Prompt cleanup so repeated assemblies in a long-lived process do not leak
        # the input FASTQ + corrector output dirs until process exit (review #2).
        rm(input_dir; recursive = true, force = true)
        rm(corrector_output_dir; recursive = true, force = true)
    end
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
td-zru6). Corrector provenance (`corrector`, `strategy`, `skip_solid`,
`k_progression`, `corrected_read_count`) is stamped onto the re-assembly's
`assembly_stats`.

The re-assembly deliberately lets `assemble_genome` auto-detect its graph mode
(DoubleStrand for DNA) rather than inheriting `config.graph_mode`. Corrected
reads are FASTQ, so the re-assembly runs the same
quality-aware (qualmer) path a naive `assemble_genome` on FASTQ reads would —
i.e. it mirrors the naive-on-FASTQ baseline, keeping the comparison apples-to-
apples.
"""
function _assemble_with_iterative_corrector(reads, config::AssemblyConfig)
    # Stage 1: materialize + persist the corrected reads (shared with the hybrid
    # OLC route, td-ohob). The helper owns the corrector's temp dirs and returns
    # the corrected reads in memory plus everything the re-assembly tail consumes
    # (result_dict for graph reuse + stat stamps, tier knobs, max_k) with no
    # recomputation — so the native output stays byte-identical.
    stage1 = _run_stage1_correction(reads, config)
    # `indel_params` is a correction-phase value the tail stamps into
    # assembly_stats["indel_moves"] (merged from the sequencing_tech/indel work) —
    # thread it through the helper's return so the extracted tail keeps that stamp.
    (; corrected_reads, result_dict, knobs, max_k, indel_params) = stage1
    n_corrected = length(corrected_reads)
    try
        # mycelia_iterative_assemble is a read CORRECTOR: its :final_assembly is the
        # corrected READS, not an assembly. Re-assemble the corrected reads through
        # the naive path so assemble_genome returns real contigs + a graph —
        # otherwise downstream (QUAST, GFA) would treat raw corrected reads as an
        # assembly, which is not an apples-to-apples assembly result (td-zru6).
        _log_info(config, "Corrected $(n_corrected) reads; re-assembling them (corrector=:none)")
        # Coverage-aware re-assembly k (td-jt7r). Pinning re-assembly to the
        # corrector's k CEILING shatters the de Bruijn graph on high-error long
        # reads (nanopore/pacbio), fragmenting the assembly and MASKING the
        # correction gain. select_reassembly_k measures solid-k-mer connectivity
        # directly and drops k to the largest prime that keeps the corrected-read
        # graph connected; clean / high-coverage (Illumina) reads keep the ceiling
        # UNCHANGED (byte-identical, graph-reuse stays eligible).
        reassembly_ceiling = config.k === nothing ? max_k : config.k
        reassembly_k = select_reassembly_k(corrected_reads, reassembly_ceiling)
        if reassembly_k != reassembly_ceiling
            _log_info(config,
                "Re-assembly k adapted $(reassembly_ceiling) -> $(reassembly_k) " *
                "(coverage-aware connectivity floor) to keep the corrected-read graph connected")
        end
        # Re-assemble with AUTO-DETECTED graph_mode (not config.graph_mode): match
        # the mode the naive baseline uses (DoubleStrand for DNA), which auto-config
        # selects, so the corrected-read re-assembly is apples-to-apples.
        #
        # Final-pass graph reuse (td-04tb). The corrector already built a qualmer
        # graph in its final pass; when that pass converged (0 improvements), the
        # graph is byte-identical to the one a from-scratch re-assembly would build
        # from the corrected reads. Resolve the re-assembly config exactly as
        # `assemble_genome(corrected_reads; k=reassembly_k, corrector=:none)` would
        # (same auto-detected sequence type / graph mode / quality flag), then reuse
        # the corrector's graph ONLY when every parameter matches AND the corrector
        # marked it reusable. Otherwise fall back to a full rebuild. The reuse path
        # calls the SAME downstream contig extraction the rebuild would, so contigs /
        # N50 / sequences are identical — this is a pure-performance change.
        # The re-assembly config is built with corrector=:none, so two additive
        # corrector behaviors must be threaded through explicitly (both default ON
        # for corrector=:iterative; a caller can force either off):
        #   * dedup_revcomp (td-47di): forward the OUTER config's resolved value so
        #     the corrected re-assembly emits the canonical (RC-deduped) contig set
        #     instead of re-inflating with both-strand twins.
        #   * graph_cleanup (td-969e): the corrector's final graph retains error-
        #     induced branch points that fragment the re-assembly; clean it BEFORE
        #     contig extraction. nothing = context default ON for the corrector.
        # Ordering inside _qualmer_graph_to_assembly is cleanup -> extract -> dedup:
        # clean_corrector_graph! runs first on the graph, find_contigs_next then
        # extracts, and the RC-dedup collapses strands last.
        corrector_cleanup = config.graph_cleanup === nothing ? true : config.graph_cleanup
        reassembly_config = _auto_configure_assembly(corrected_reads;
            k = reassembly_k, corrector = :none,
            dedup_revcomp = config.dedup_revcomp,
            graph_cleanup = corrector_cleanup)
        reused_graph = get(result_dict, :final_graph, nothing)
        _rmeta = result_dict[:metadata]
        can_reuse_graph = reused_graph !== nothing &&
                          get(_rmeta, :final_graph_reusable, false) === true &&
                          get(_rmeta, :final_graph_k, nothing) == reassembly_config.k &&
                          get(_rmeta, :final_graph_mode, nothing) ==
                          _graph_mode_symbol(reassembly_config.graph_mode) &&
                          reassembly_config.k !== nothing &&
                          reassembly_config.use_quality_scores &&
                          reassembly_config.sequence_type <: BioSequences.BioSequence
        assembly = if can_reuse_graph
            _log_info(config,
                "Reusing corrector final-pass qualmer graph (k=$(reassembly_config.k), " *
                "mode=$(_graph_mode_symbol(reassembly_config.graph_mode))); " *
                "skipping redundant from-scratch build_qualmer_graph (td-04tb)")
            _qualmer_graph_to_assembly(reused_graph, n_corrected, reassembly_config;
                graph_cleanup = corrector_cleanup)
        else
            _log_info(config,
                "Corrector final-pass graph not reusable " *
                "(reusable=$(get(_rmeta, :final_graph_reusable, false)), " *
                "graph_k=$(get(_rmeta, :final_graph_k, nothing)) vs reassembly_k=$(reassembly_config.k)); " *
                "rebuilding from corrected reads")
            assemble_genome(corrected_reads, reassembly_config)
        end
        # Stamp the corrector provenance onto the real assembly's stats.
        assembly.assembly_stats["corrector"] = "iterative"
        assembly.assembly_stats["strategy"] = String(config.strategy)
        # Indel-aware correction provenance (td-9q84 / 4a): the driving error profile
        # and whether it engaged the pair-HMM gap moves. `indel_moves=false` (the
        # :illumina default) marks the substitution-only oracle path.
        assembly.assembly_stats["sequencing_tech"] = String(config.sequencing_tech)
        assembly.assembly_stats["indel_moves"] = indel_params !== nothing
        # Final-pass graph reuse provenance (td-04tb): true when the re-assembly
        # reused the corrector's already-built final-pass graph (converged run),
        # false when it rebuilt from scratch. Pure telemetry — does not affect the
        # contigs/N50/sequences (identical either way).
        assembly.assembly_stats["reassembly_graph_reused"] = can_reuse_graph
        # Coverage-aware re-assembly k provenance (td-jt7r): the ceiling requested
        # (corrector k) and the k actually used after the connectivity criterion.
        # When they differ the corrected-read graph was adapted DOWN to stay
        # connected (the un-shatter fix for high-error long reads).
        assembly.assembly_stats["reassembly_k_ceiling"] = reassembly_ceiling
        assembly.assembly_stats["reassembly_k"] = reassembly_k
        # Provenance (FIX 4): stamp BOTH the caller's requested skip_solid and the
        # value the tier actually used, so the tier's silent override cannot hide.
        # `skip_solid` is retained as the EFFECTIVE value for back-compat.
        assembly.assembly_stats["skip_solid_requested"] = config.skip_solid
        assembly.assembly_stats["skip_solid_effective"] = knobs.skip_solid
        assembly.assembly_stats["skip_solid"] = knobs.skip_solid
        assembly.assembly_stats["k_progression"] = result_dict[:k_progression]
        assembly.assembly_stats["corrected_read_count"] = n_corrected
        # Scalable-tier telemetry (hard-read skip fraction + honest gate flags).
        _corr_meta = result_dict[:metadata]
        # `hard_window`/`hard_read_gate` = the skip gate (active on :scalable);
        # `windowed_decode` = per-hard-region windowed decode (td-nn6l Stage 3c),
        # now ACTIVE on :scalable ⇒ hard reads decoded window-by-window (bounded),
        # not whole. Kept distinct from the skip gate so the surfaced flag is honest.
        assembly.assembly_stats["hard_window"] = get(_corr_meta, :hard_window, false)
        assembly.assembly_stats["hard_read_gate"] = get(_corr_meta, :hard_read_gate, false)
        assembly.assembly_stats["windowed_decode"] = get(_corr_meta, :windowed_decode, false)
        corrector_errors = get(_corr_meta, :corrector_errors, Dict{Symbol, Int}())
        assembly.assembly_stats["corrector_errors"] = corrector_errors
        assembly.assembly_stats["indel_decodes"] =
            get(corrector_errors, :indel_decodes, 0)
        assembly.assembly_stats["truncated_decodes"] =
            get(corrector_errors, :truncated_decodes, 0)
        assembly.assembly_stats["trace_contract_errors"] =
            get(corrector_errors, :trace_contract_errors, 0)
        assembly.assembly_stats["window_divergences"] =
            get(corrector_errors, :window_divergences, 0)
        # soft-EM v1 is record-only ⇒ "scaffold-v1-record-only" (or false on
        # :exhaustive), never a bare `true` (FIX 1/5).
        assembly.assembly_stats["soft_em"] = get(_corr_meta, :soft_em, false)
        assembly.assembly_stats["skip_fraction"] = get(_corr_meta, :last_skip_fraction, 0.0)
        assembly.assembly_stats["skip_fraction_per_pass"] = get(_corr_meta, :skip_fraction_per_pass, Float64[])
        # Stage 0 cheap correction + graph-decode fraction (td-bjnt). The decode
        # fraction is the critical-path win metric (reads that reached graph
        # Viterbi); Stage 0 + hard-set narrowing drive it toward the true ~5-15%.
        assembly.assembly_stats["cheap_correct"] = get(_corr_meta, :cheap_correct, false)
        assembly.assembly_stats["cheap_corrections_total"] = get(_corr_meta, :cheap_corrections_total, 0)
        assembly.assembly_stats["cheap_corrections_per_pass"] = get(
            _corr_meta, :cheap_corrections_per_pass, Int[])
        assembly.assembly_stats["decode_fraction"] = get(_corr_meta, :decode_fraction_mean, 0.0)
        assembly.assembly_stats["decode_fraction_per_pass"] = get(
            _corr_meta, :decode_fraction_per_pass, Float64[])
        if isempty(assembly.contigs)
            @warn "corrector=:iterative re-assembled $(n_corrected) corrected reads into 0 " *
                  "contigs — the corrected read set did not assemble."
        end
        _log_info(config, "Re-assembled corrected reads into $(length(assembly.contigs)) contigs")
        return assembly
    finally
        # The persisted corrected FASTQ is an ephemeral tempfile ONLY when the
        # helper minted one (config.output_dir unset); then it is ours to delete,
        # preserving the historical no-stray-file behavior of the native
        # re-assembly. A caller-supplied output_dir is the caller's to keep for the
        # hybrid-OLC handoff (td-ohob). We branch on the helper's returned
        # `ephemeral` flag rather than re-deriving the predicate from config, so the
        # ownership decision lives at the single point that made it. The corrector's
        # own temp dirs were already cleaned by _run_stage1_correction.
        if stage1.ephemeral
            rm(stage1.corrected_fastq; force = true)
        end
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

Quality note (graph-as-HMM correction redesign,
`docs/design/2026-07-06-graph-as-hmm-correction-redesign.md`): this helper
normalizes to FASTA records, which DROPS per-base quality. That is intentional
for the naive k-mer / string / n-gram graph builders, which do not consume
quality. The iterative CORRECTOR does NOT use this path — it preserves real
FASTQ quality end-to-end via `_write_reads_to_fastq` (FASTQ written verbatim;
FASTA flagged quality-absent with a warning, never silently Q40) and the per-read
per-base emission wrapping in `try_viterbi_path_improvement`. The remaining gap is
the qualmer GRAPH builder (`_assemble_qualmer_graph`), which currently receives
these quality-stripped observations and then re-injects a Q40 placeholder in
`_prepare_fastq_observations`; giving it real FASTQ quality requires routing raw
reads through the `assemble_genome` dispatch (a separate track) rather than this
FASTA-normalizing helper.
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
        # PRESERVE FASTQ records (do NOT strip quality). FASTQ input always routes to
        # the quality-aware paths (use_quality_scores=true → qualmer / quality-biosequence),
        # which consume per-base quality; stripping to FASTA here forced the downstream
        # _prepare_fastq_observations to re-inject placeholder Q40 (td-tps5). FASTA input
        # (which legitimately lacks quality) still flows through the FASTA branch above.
        return reads
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
    graph = Rhizomorph.build_qualmer_graph(fastq_records, config.k; mode = mode)

    # Contig extraction + AssemblyResult assembly is factored into
    # `_qualmer_graph_to_assembly` so the iterative corrector can REUSE its already-
    # built final-pass qualmer graph (td-04tb) and skip the redundant from-scratch
    # build_qualmer_graph here — the two paths share this exact downstream code, so
    # the reused-graph assembly is byte-identical to a from-scratch re-assembly.
    # A plain qualmer assembly cleans only when the caller explicitly opts in
    # (config.graph_cleanup === true). The corrector's re-assembly path opts in
    # via the reuse call in the :iterative branch and via reassembly_config below.
    return _qualmer_graph_to_assembly(graph, length(observations), config;
        graph_cleanup = config.graph_cleanup === true)
end

"""
    _qualmer_graph_to_assembly(graph, num_input_sequences::Int, config) -> AssemblyResult

Extract contigs from an already-built qualmer `graph` and assemble the
`AssemblyResult` (paths -> consensus FASTQ contigs -> stats). Split out of
`_assemble_qualmer_graph` (td-04tb) so the iterative corrector's re-assembly can
REUSE the corrector's final-pass qualmer graph instead of rebuilding an identical
one from scratch. Given the same `graph`, `num_input_sequences`, and `config`
this is a pure function of the graph structure, so a reused graph yields
byte-identical contigs/N50/sequences.
"""
function _qualmer_graph_to_assembly(graph, num_input_sequences::Int, config;
        graph_cleanup::Bool = false)
    # Canonical graph_mode now supports correct contig reconstruction on the qualmer
    # arm as well. Contigs are reconstructed by the SAME orientation-aware path used
    # by the k-mer arm (_qualmer_path_to_consensus_fastq -> path_to_sequence ->
    # _reconstruct_oriented_kmer_path), and per-base quality is oriented to match
    # (reverse-oriented k-mers take reversed / [1]-indexed quality). Canonical
    # reconstruction therefore matches DoubleStrand for both sequence AND quality,
    # so the result is flagged valid and no warning is emitted (see
    # _assemble_kmer_graph for the sequence-side rationale).
    reconstruction_valid = true

    # Apply Phase 2 graph algorithms if enabled
    if config.bubble_resolution
        # Note: This would use detect_bubbles_next adapted for qualmer graphs
        _log_info(config, "Bubble resolution enabled for qualmer graphs")
    end

    if config.repeat_resolution
        # Note: This would use resolve_repeats_next adapted for qualmer graphs
        _log_info(config, "Repeat resolution enabled for qualmer graphs")
    end

    # Linear-time graph defragmentation BEFORE contig extraction (td-969e). The
    # scalable corrector's final graph keeps error-induced branch points (dead-end
    # tips + coverage-1 bubbles); find_contigs_next breaks a unitig at every branch
    # vertex, so each residual branch point is a contig boundary. clean_corrector_
    # graph! removes ONLY unambiguous errors (coverage-1 dead-end tips + guarded
    # low-coverage bubble collapse) in O(V+E), never collapsing a data-supported
    # variant (the td-h6w9 invariant). Operate on a deepcopy so a REUSED corrector
    # final-pass graph (td-04tb) is not mutated for other consumers; the cleaned
    # copy is what we both extract contigs from AND return in the AssemblyResult.
    cleanup_stats = nothing
    if graph_cleanup
        graph = deepcopy(graph)
        cleanup_stats = Rhizomorph.clean_corrector_graph!(
            graph; k = config.k === nothing ? 21 : config.k)
        _log_info(config,
            "Graph cleanup (td-969e): $(cleanup_stats["graph_cleanup_vertices_before"]) -> " *
            "$(cleanup_stats["graph_cleanup_vertices_after"]) vertices " *
            "($(cleanup_stats["graph_cleanup_tips_removed"]) tips clipped, " *
            "$(cleanup_stats["graph_cleanup_bubbles_collapsed"]) error bubbles collapsed, " *
            "$(get(cleanup_stats, "graph_cleanup_components_pruned", 0)) error components pruned)")
    end

    # Find contigs by extracting maximal unitigs (non-branching paths), mirroring
    # the k-mer arm. The prior heaviest-path / iterative-Viterbi heuristics found
    # no substantial paths on branchy read graphs and fell to a placeholder walk
    # that emitted ~one short fragment per vertex, collapsing the largest contig
    # far below the genome length. find_contigs_next returns the vertex path
    # (ContigPath.vertices); per-base quality is then propagated by
    # _qualmer_path_to_consensus_fastq below.
    #
    # Canonical exception: on the UNDIRECTED canonical qualmer graph, unitig
    # extraction does not yield an overlap-ordered traversal (adjacent canonical
    # labels are stored lexicographically-minimal, so a non-branching-path walk
    # fragments and its reconstruction is invalid). Mirror the k-mer arm and prefer
    # the Eulerian traversal, which IS orientation-reconstructable by
    # _qualmer_path_to_consensus_fastq -> path_to_sequence. Fall back to unitig
    # extraction only when no Eulerian path exists (branchy/error-laden real inputs),
    # preserving the existing behavior for single/doublestrand and hard cases.
    # rc_aware=config.dedup_revcomp (td-47di): on a DoubleStrand qualmer graph the
    # unitig walk otherwise emits BOTH strands of every contig (forward + reverse-
    # complement twin), ~2x-inflating the assembly. Passing rc_aware marks each walked
    # vertex's RC partner visited so the reverse strand is never independently
    # traversed — the same structural dedup the k-mer arm uses via
    # _generate_contigs_probabilistic.
    paths = if config.graph_mode == Canonical
        eulerian_paths = Rhizomorph.find_eulerian_paths_next(graph)
        if isempty(eulerian_paths)
            [contig_path.vertices
             for contig_path in
                 Rhizomorph.find_contigs_next(
                graph; min_contig_length = config.k + 1,
                rc_aware = config.dedup_revcomp)]
        else
            eulerian_paths
        end
    else
        [contig_path.vertices
         for contig_path in
             Rhizomorph.find_contigs_next(
            graph; min_contig_length = config.k + 1,
            rc_aware = config.dedup_revcomp)]
    end

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

    # Belt-and-suspenders RC dedup (td-47di): the structural rc_aware traversal above
    # removes twins that share a breakpoint, but RC twins with OFFSET fragment
    # breakpoints (or the probabilistic-walk fallback, which does not honor rc_aware)
    # can still slip through — the same case the k-mer arm covers with a post-hoc
    # _dedup_reverse_complements. We dedup by canonical key AND emit each survivor in
    # its CANONICAL orientation (min(seq, reverse_complement(seq))), reversing the
    # per-base quality when the reverse complement is canonical so quality stays
    # registered against the emitted sequence. Emitting the canonical orientation
    # (rather than keeping whichever strand was traversed) makes the output
    # ORDER-INDEPENDENT: rc_aware picks one strand per pair based on graph iteration
    # order, so the corrector's REUSED final-pass graph and a from-scratch REBUILD
    # could otherwise emit opposite strands for the same locus and diverge. Canonical
    # emission keeps them byte-identical (reassembly_graph_reuse_test invariant).
    # Gated on config.dedup_revcomp (default ON for the corrector route), so plain
    # both-strands assemblies are unchanged.
    if config.dedup_revcomp && length(contig_records) > 1
        seen = Set{String}()
        deduped = FASTX.FASTQ.Record[]
        for record in contig_records
            canonical_record = _canonical_fastq_record(record)
            canonical_seq = String(FASTX.sequence(canonical_record))
            if !(canonical_seq in seen)
                push!(seen, canonical_seq)
                push!(deduped, canonical_record)
            end
        end
        if length(deduped) < length(contig_records)
            _log_info(config,
                "Reverse-complement dedup (qualmer): $(length(contig_records)) -> " *
                "$(length(deduped)) contig records")
        end
        contig_records = deduped
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
        # Always true: canonical (undirected) qualmer reconstruction is now
        # orientation-aware for both sequence and per-base quality.
        "reconstruction_valid" => reconstruction_valid,
        "num_fastq_contigs" => length(contig_records),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => num_input_sequences,
        "assembly_date" => string(Mycelia.Dates.now())
    )

    # Add quality-specific statistics
    if !isempty(MetaGraphsNext.labels(graph))
        qualmer_stats = Rhizomorph.get_qualmer_statistics(graph)
        for (key, value) in qualmer_stats
            stats[string(key)] = value
        end
    end

    # Graph-cleanup provenance (td-969e). Present only when the defrag pass ran.
    stats["graph_cleanup_applied"] = graph_cleanup
    if cleanup_stats !== nothing
        for (key, value) in cleanup_stats
            stats[key] = value
        end
    end

    return AssemblyResult(contigs, contig_names; graph = graph,
        assembly_stats = stats, fastq_contigs = contig_records)
end

"""
    _assemble_hybrid_olc(reads, config::AssemblyConfig) -> AssemblyResult

Hybrid-OLC route (a), td-yymj: run Stage-1 correction, then hand the corrected
reads to an EXTERNAL overlap-layout-consensus assembler instead of the native
graph-HMM re-assembly. Contiguity is inherited from a tool tuned for it while the
per-base accuracy gain is contributed by Stage 1.

Reuses `_run_stage1_correction` (which persists the corrected FASTQ) and honors its
ephemeral-cleanup contract: when `stage1.ephemeral` (no `config.output_dir`), this
function deletes the corrected FASTQ and the external assembler's temp output dir
after reading the contigs; when the caller supplied `output_dir`, both are left in
place under it. Dispatched only via `assemble_genome` when `config.layout == :olc`
(the constructor guarantees that implies `corrector == :iterative`).

Both arms are wired: short-read (:megahit/:metaspades, td-yymj) and long-read
(:flye/:metaflye/:canu/:hifiasm, td-wvto), routed by sequencing_tech. The
two-read-set hybrid arm (short+long combined) is td-06er.
"""
function _assemble_hybrid_olc(reads, config::AssemblyConfig)
    tool = _resolve_olc_tool(config)
    _log_info(config,
        "Hybrid-OLC route (a): Stage-1 correction -> external assembler :$(tool)")
    stage1 = _run_stage1_correction(reads, config; materialize_corrected_reads = false)
    # The external assembler writes into its own output dir. Co-locate it with the
    # persisted corrected FASTQ when the caller owns output_dir; otherwise a temp
    # dir cleaned alongside the ephemeral corrected FASTQ.
    olc_outdir = config.output_dir === nothing ?
                 mktempdir() : mkpath(joinpath(config.output_dir, "olc_$(tool)"))
    # Stale-output guard: the external wrappers skip re-running when their contigs
    # file already exists, so a reused (non-empty) output_dir would silently return
    # a PRIOR run's assembly, ignoring THESE corrected reads. Warn loudly.
    if config.output_dir !== nothing && !isempty(readdir(olc_outdir))
        @warn "hybrid-OLC: assembler output dir is non-empty; the external tool may " *
              "skip re-running and return a STALE prior assembly — use a fresh " *
              "output_dir per run." olc_outdir
    end
    try
        contigs_fasta = _run_olc_tool(tool, stage1.corrected_fastq, olc_outdir, config)
        if !isfile(contigs_fasta)
            error("external OLC assembler :$(tool) produced no contigs file " *
                  "(expected at $(contigs_fasta))")
        end
        _log_info(config, "External :$(tool) assembly complete; wrapping contigs")
        return _wrap_external_contigs(contigs_fasta, tool, config, stage1)
    finally
        # Ephemeral (output_dir unset): the corrected FASTQ and the assembler's temp
        # output dir are both ours to clean. A caller-supplied output_dir keeps them.
        # Wrap cleanup so a failing rm (e.g. a locked file on an HPC/NFS mount) is
        # WARNed rather than replacing the in-flight assembler exception — Julia's
        # finally-throw discards the original error.
        if stage1.ephemeral
            try
                rm(stage1.corrected_fastq; force = true)
                rm(olc_outdir; recursive = true, force = true)
            catch cleanup_err
                @warn "hybrid-OLC: cleanup of ephemeral artifacts failed" cleanup_err
            end
        end
    end
end

"""
    _olc_taxonomy() -> NamedTuple

Single source of truth for the hybrid-OLC tool taxonomy: which sequencing techs
are short- vs long-read, and which external assemblers serve each class. The
constructor validation (`_valid_seq_techs`, `:auto` resolution) and the per-tool
adapter all read from THIS table so they cannot drift out of sync. (Drift is
exactly what let a stale `:metaflye` branch masquerade as "wired" before td-yymj's
review.) The "taxonomy single-source" test in `hybrid_olc_config_test.jl` asserts
`:auto` resolves to a member of the wired-tool set for every valid tech, and that
`_run_olc_tool` rejects a tool not in it. (A stronger test iterating every wired
tool through `_run_olc_tool` needs an external-wrapper stub — a follow-up.)
"""
function _olc_taxonomy()
    return (
        short_read_techs = (:illumina, :ultima),
        long_read_techs = (:nanopore, :pacbio),
        short_read_tools = (:megahit, :metaspades),
        long_read_tools = (:flye, :metaflye, :canu, :hifiasm)
    )
end

"""
    _resolve_olc_tool(config::AssemblyConfig) -> Symbol

Resolve `config.olc_tool` to a concrete external assembler. An explicit tool is
returned as-is (the constructor already validated it against `sequencing_tech`).
`:auto` picks by read type: a short-read tech → `:megahit` (single-end-robust
short-read layout), a long-read tech → `:flye` (general long-read assembler that
handles both Nanopore and PacBio via its read-type flag).
"""
function _resolve_olc_tool(config::AssemblyConfig)
    config.olc_tool == :auto || return config.olc_tool
    return config.sequencing_tech in _olc_taxonomy().short_read_techs ? :megahit : :flye
end

"""
    _flye_read_type(sequencing_tech::Symbol) -> String

Map the corrector's read tech to a Flye/metaFlye `--read-type` flag. The reads are
Stage-1 corrector output, so the error-corrected (`-corr`) profiles fit. NOTE this
means a native PacBio-HiFi set is assembled as `pacbio-corr`, not `pacbio-hifi` —
`:pacbio` does not distinguish HiFi from CLR (finer tech granularity is a
follow-up). Fails loud on an unexpected tech rather than silently defaulting.
"""
function _flye_read_type(sequencing_tech::Symbol)
    if sequencing_tech == :nanopore
        return "nano-corr"
    elseif sequencing_tech == :pacbio
        return "pacbio-corr"
    else
        error("_flye_read_type: unexpected sequencing_tech :$(sequencing_tech) " *
              "(expected :nanopore or :pacbio).")
    end
end

"""
    _canu_read_type(sequencing_tech::Symbol) -> String

Canu's coarse read-type flag (`nanopore` | `pacbio`). Fails loud on an unexpected
tech rather than silently defaulting to PacBio.
"""
function _canu_read_type(sequencing_tech::Symbol)
    if sequencing_tech == :nanopore
        return "nanopore"
    elseif sequencing_tech == :pacbio
        return "pacbio"
    else
        error("_canu_read_type: unexpected sequencing_tech :$(sequencing_tech) " *
              "(expected :nanopore or :pacbio).")
    end
end

"""
    _run_olc_tool(tool, corrected_fastq, outdir, config) -> contigs_fasta_path::String

Per-tool argument adapter: the external-assembler wrappers have non-uniform
signatures, so map each wired `tool` onto its wrapper call and return the path to
its primary-contigs FASTA. Short-read assemblers (`:megahit`, `:metaspades`) take
`fastq1` (single-end — the corrector emits ONE corrected FASTQ); long-read
assemblers take a single `fastq` plus a read-type flag derived from
`sequencing_tech`. hifiasm exposes contigs via `hifiasm_primary_contigs` rather
than a return field; Flye/metaFlye/Canu return their assembly under `.assembly`.

`config.olc_options` is splatted in FIRST so the route's managed keywords
(`fastq1`/`fastq`/`outdir`/`read_type`) always win over a caller-supplied option —
Julia resolves duplicate keywords rightmost-wins, so managed keys must come last.
(The constructor also rejects reserved keys — incl. `:read_type` — but this
ordering keeps the invariant even if that guard is ever relaxed.)
"""
function _run_olc_tool(tool::Symbol, corrected_fastq::AbstractString,
        outdir::AbstractString, config::AssemblyConfig)::String
    fastq = String(corrected_fastq)
    out = String(outdir)
    if tool == :megahit
        result = Mycelia.run_megahit(; config.olc_options..., fastq1 = fastq, outdir = out)
        return result.contigs
    elseif tool == :metaspades
        result = Mycelia.run_metaspades(; config.olc_options..., fastq1 = fastq, outdir = out)
        return result.contigs
    elseif tool == :flye
        result = Mycelia.run_flye(; config.olc_options..., fastq = fastq, outdir = out,
            read_type = _flye_read_type(config.sequencing_tech))
        return result.assembly
    elseif tool == :metaflye
        result = Mycelia.run_metaflye(; config.olc_options..., fastq = fastq, outdir = out,
            read_type = _flye_read_type(config.sequencing_tech))
        return result.assembly
    elseif tool == :canu
        # genome_size comes from olc_options (required — constructor-validated).
        result = Mycelia.run_canu(; config.olc_options..., fastq = fastq, outdir = out,
            read_type = _canu_read_type(config.sequencing_tech))
        return result.assembly
    elseif tool == :hifiasm
        result = Mycelia.run_hifiasm(; config.olc_options..., fastq = fastq, outdir = out)
        contigs = Mycelia.hifiasm_primary_contigs(result)
        contigs === nothing && error("_run_olc_tool: hifiasm produced no primary " *
              "contigs (.p_ctg.fa) in $(out)")
        return contigs
    else
        _tax = _olc_taxonomy()
        error("_run_olc_tool: :$(tool) is not a wired OLC tool; expected one of " *
              "$((_tax.short_read_tools..., _tax.long_read_tools...)).")
    end
end

"""
    _wrap_external_contigs(contigs_fasta, tool, config, stage1) -> AssemblyResult

Read an external assembler's primary-contigs FASTA and wrap it in an
`AssemblyResult`. An external contig-only assembly has NO Mycelia graph, so
`gfa_compatible` MUST be false (a GFA-compatible result requires a graph, enforced
by `validate_assembly_structure`). Corrector + layout provenance is stamped into
`assembly_stats`. Fails loud on a 0-contig assembly (parity with Stage-1's
0-corrected-reads guard) rather than returning a silently-empty result.
"""
function _wrap_external_contigs(contigs_fasta::AbstractString, tool::Symbol,
        config::AssemblyConfig, stage1)
    records = open(FASTX.FASTA.Reader, contigs_fasta) do reader
        collect(reader)
    end
    # FASTX.sequence(String, r) already yields a String; FASTX.identifier returns a
    # StringView, so wrap it — AssemblyResult requires Vector{String} for both.
    contigs = [FASTX.sequence(String, r) for r in records]
    contig_names = [String(FASTX.identifier(r)) for r in records]
    if isempty(contigs)
        # Fail loud, matching _run_stage1_correction's 0-reads guard: an empty
        # external assembly must not flow downstream as a "successful" 0-contig
        # result. (megahit errors upstream on empty input; metaspades can write an
        # empty contigs.fasta, so guard here for tool-agnostic parity.)
        error("hybrid-OLC: external assembler :$(tool) produced 0 contigs from the " *
              "$(stage1.corrected_read_count) corrected reads (from $(contigs_fasta)); " *
              "refusing to return an empty assembly.")
    end
    stats = Dict{String, Any}(
        "method" => "HybridOLC",
        "corrector" => "iterative",
        "strategy" => String(config.strategy),
        "layout" => "olc",
        "olc_tool" => String(tool),
        "sequencing_tech" => String(config.sequencing_tech),
        "corrected_read_count" => stage1.corrected_read_count,
        "corrected_fastq" => stage1.corrected_fastq,
        "num_contigs" => length(contigs),
        "assembly_date" => string(Mycelia.Dates.now())
    )
    return AssemblyResult(contigs, contig_names;
        assembly_stats = stats, gfa_compatible = false)
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
    _canonical_fastq_record(record::FASTX.FASTQ.Record) -> FASTX.FASTQ.Record

Return `record` re-oriented to its canonical strand: if the reverse complement of
the sequence is lexicographically smaller than the forward sequence, emit a new
record holding the reverse-complemented sequence with its per-base quality REVERSED
(Phred values are strand-invariant; only their order flips under reverse-complement),
otherwise return `record` unchanged. Sequences that are not valid DNA (no defined
reverse complement) are passed through unchanged. Used by the qualmer RC-dedup so a
contig's emitted orientation is a deterministic function of its content, not of the
graph-traversal order — keeping the corrector's reused-graph and from-scratch-rebuild
assemblies byte-identical under dedup.
"""
function _canonical_fastq_record(record::FASTX.FASTQ.Record)::FASTX.FASTQ.Record
    seq = String(FASTX.sequence(record))
    canonical = _canonical_string(seq)
    canonical == seq && return record
    id = String(FASTX.identifier(record))
    quality = reverse(collect(FASTX.quality_scores(record)))
    return FASTX.FASTQ.Record(id, canonical, quality)
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
    # Mirror `_write_reads_to_fastq`: any FASTA (quality-less) input that gets a
    # placeholder Q40 must be SURFACED once per call, never silently substituted.
    # A degenerate uniform-Q40 model makes the corrector's decisions
    # coverage-driven only, so a silent substitution masks that the quality model
    # is inert (review: silent Q40).
    placeholder_used = false

    for (i, obs) in enumerate(observations)
        # Handle different observation formats
        if obs isa FASTX.FASTQ.Record
            push!(fastq_records, obs)
        elseif obs isa FASTX.FASTA.Record
            # Convert FASTA to FASTQ with placeholder quality (Q40).
            seq = FASTX.FASTA.sequence(obs)
            qual = repeat('I', length(seq))  # Phred+33 'I' == Q40 placeholder
            record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(obs), seq, qual)
            push!(fastq_records, record)
            placeholder_used = true
        elseif obs isa Tuple && length(obs) >= 1
            # Handle tuple format like (record, index) or (record, other_data)
            record = obs[1]
            if record isa FASTX.FASTQ.Record
                push!(fastq_records, record)
            elseif record isa FASTX.FASTA.Record
                # Convert FASTA to FASTQ with placeholder quality (Q40).
                seq = FASTX.FASTA.sequence(record)
                qual = repeat('I', length(seq))  # Phred+33 'I' == Q40 placeholder
                fastq_record = FASTX.FASTQ.Record(FASTX.FASTA.identifier(record), seq, qual)
                push!(fastq_records, fastq_record)
                placeholder_used = true
            else
                @warn "Unsupported record type in tuple: $(typeof(record))"
            end
        else
            @warn "Unsupported observation type: $(typeof(obs))"
        end
    end

    if placeholder_used
        @warn "corrector: observation input lacked per-base quality (FASTA / quality-less); " *
              "assigned placeholder Q40. The corrector's quality model is degenerate (uniform) " *
              "for these reads — its decisions become coverage-driven only."
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

    # Per-vertex orientation for the emitted contig. On an undirected (canonical)
    # nucleotide graph, path_to_sequence reverse-complements reverse-oriented
    # k-mers; the per-base quality vector MUST follow the same orientation or the
    # quality string will be mis-registered against the (corrected) sequence. The
    # stored quality vector is aligned to the canonical label, so for a reverse-
    # oriented k-mer the emitted base at position j is label[k+1-j] complemented and
    # its quality is quality[k+1-j] (i.e. reverse the vector; the emitted LAST base
    # takes quality[1], not quality[end]). Directed (single/doublestrand) graphs and
    # non-nucleotide labels keep the plain forward orientation.
    T = eltype(path)
    orientation_aware = !Graphs.is_directed(graph) &&
                        Rhizomorph._label_has_reverse_complement(T) &&
                        Rhizomorph._sequence_type_from_kmer_type(T) <:
                        BioSequences.LongSequence
    forward_flags = orientation_aware ? Rhizomorph._resolve_path_orientations(path, graph) :
                    trues(length(path))

    quality_scores = UInt8[]

    for (idx, vertex) in enumerate(path)
        vertex_quality = _qualmer_vertex_quality_scores(graph[vertex])
        kmer_len = length(string(vertex))
        forward = forward_flags[idx]

        if idx == 1
            if length(vertex_quality) == kmer_len
                # Full first k-mer: reverse the quality vector when the k-mer is
                # emitted reverse-complemented so quality[j] tracks the emitted base.
                append!(quality_scores, forward ? vertex_quality : reverse(vertex_quality))
            else
                append!(quality_scores, fill(UInt8(2), kmer_len))
            end
        else
            # Only the last emitted base of this k-mer is appended. Forward: that is
            # label[end] -> quality[end]. Reverse: it is complement(label[1]) ->
            # quality[1].
            if isempty(vertex_quality)
                push!(quality_scores, UInt8(2))
            else
                push!(quality_scores, forward ? vertex_quality[end] : vertex_quality[1])
            end
        end
    end

    # Sequence and quality are both length (k + (n-1)) by construction, so this
    # guard is defensive and should not fire; it must never silently SHORTEN a
    # correctly-reconstructed contig, so surface a mismatch if one ever occurs.
    if length(quality_scores) != length(sequence)
        @warn "qualmer consensus: quality/sequence length mismatch " *
              "($(length(quality_scores)) vs $(length(sequence))) for $(contig_name); " *
              "clamping to the shorter length."
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
Token graph assembly implementation (variable-length SentencePiece/word-token
analysis).

Assembles pre-tokenized sequences (e.g. SentencePiece pieces or whitespace
tokens) into contigs. It builds a token graph over the String token labels,
finds an Eulerian path when one exists (or, when it does not, the maximal linear
contigs), and reconstructs each contig by joining that path's token labels with
`token_separator`.

Modeled on `_assemble_ngram_graph`, with one deliberate difference: token
vertices are whole opaque strings, not fixed-length character windows, so the
character-overlap reconstruction used by `path_to_sequence` /
`generate_contig_sequence` (append only the last char of each subsequent vertex)
is wrong here. Instead the ordered token labels are joined directly. Tokens are
SingleStrand only — reverse complement is undefined for general string tokens —
so no reverse-complement handling is performed.
"""
function _assemble_token_graph(
        token_sequences::Vector{Vector{String}}, config;
        token_separator::AbstractString = " "
)
    _log_info(config, "Using token graph assembly strategy (variable-length token analysis)")

    graph = Rhizomorph.build_token_graph(token_sequences)

    # Fast path: a single Eulerian traversal covering every edge, if one exists.
    paths = Rhizomorph.find_eulerian_paths_next(graph)
    contigs = String[]
    for path in paths
        if length(path) > 1
            push!(contigs, join(path, token_separator))
        end
    end

    # Fallback: emit maximal linear contigs (token label paths). `min_contig_length`
    # is measured against find_contigs_next's character-overlap string, so use 1
    # to keep short token contigs (single-token paths included).
    if isempty(contigs)
        contig_paths = Rhizomorph.find_contigs_next(graph; min_contig_length = 1)
        for contig in contig_paths
            if !isempty(contig.vertices)
                push!(contigs, join(contig.vertices, token_separator))
            end
        end
    end

    contig_names = ["token_contig_$i" for i in 1:length(contigs)]

    stats = Dict{String, Any}(
        "method" => "TokenGraph",
        "graph_type" => "variable_length",
        "vertex_type" => "string_tokens",
        "token_separator" => String(token_separator),
        "graph_mode" => string(config.graph_mode),
        "num_vertices" => length(MetaGraphsNext.labels(graph)),
        "num_edges" => length(MetaGraphsNext.edge_labels(graph)),
        "num_input_sequences" => length(token_sequences),
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
