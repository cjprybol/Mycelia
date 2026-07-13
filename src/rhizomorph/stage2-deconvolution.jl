"""
Compact transition-probability index extracted from the Stage-1 accurized graph.

The index deliberately excludes read-level evidence and graph topology so the
full corrector graph can be released before the external layout/deconvolution
tools run.
"""
struct TransitionLikelihoodIndex{KmerType <: Kmers.DNAKmer}
    k::Int
    vertex_ids::Dict{KmerType, UInt32}
    log2_probability::Dict{UInt64, Float64}
end

function _stage2_edge_key(source_id::UInt32, destination_id::UInt32)::UInt64
    return (UInt64(source_id) << 32) | UInt64(destination_id)
end

function TransitionLikelihoodIndex(
        k::Int,
        probabilities::AbstractDict{
            <:Tuple{<:AbstractString, <:AbstractString}, <:Real}
)::TransitionLikelihoodIndex
    k > 0 || throw(ArgumentError("k must be positive"))
    KmerType = Kmers.derive_type(Kmers.DNAKmer{k})
    vertex_ids = Dict{KmerType, UInt32}()
    encoded_probabilities = Dict{UInt64, Float64}()
    for ((source, destination), probability) in probabilities
        source_kmer = KmerType(source)
        destination_kmer = KmerType(destination)
        source_id = get!(
            vertex_ids, source_kmer, UInt32(length(vertex_ids) + 1))
        destination_id = get!(
            vertex_ids, destination_kmer, UInt32(length(vertex_ids) + 1))
        encoded_probabilities[_stage2_edge_key(source_id, destination_id)] =
            Float64(probability)
    end
    return TransitionLikelihoodIndex(k, vertex_ids, encoded_probabilities)
end

"""
One ranked Stage-2 haplotype candidate.
"""
struct RankedVariant
    likelihood_rank::Int
    abundance_rank::Int
    role::Symbol
    id::String
    sequence::String
    graph_mean_log2_probability::Float64
    scored_transition_fraction::Float64
    strainy_coverage::Float64
end

"""Fitted soft-mixture abundance attached to one Stage-2 candidate."""
struct ProbabilisticVariantRanking
    candidate_id::String
    likelihood_rank::Int
    abundance_rank::Int
    abundance::Float64
    expected_read_count::Float64
    has_alignment_support::Bool
end

"""Soft abundance inference and the resulting Stage-2 primary agreement."""
struct ProbabilisticVariantRankingResult
    rankings::Vector{ProbabilisticVariantRanking}
    inference::CandidateInferenceResult
    primary_agreement::Bool
    primary_status::Symbol
    candidate_set_sha256::String
    inference_input_sha256::String
end

function _stage2_candidate_set_sha256(
        variants::AbstractVector{RankedVariant}
)::String
    context = Mycelia.SHA.SHA2_256_CTX()
    integer_buffer = Vector{UInt8}(undef, 8)
    tag_buffer = Vector{UInt8}(undef, 1)
    order = sortperm(eachindex(variants); by = index -> variants[index].id)
    _stage2_update_uint64!(context, integer_buffer, UInt64(length(order)))
    for index in order
        variant = variants[index]
        _stage2_update_tag!(context, tag_buffer, 0x20)
        _stage2_update_bytes!(context, integer_buffer, variant.id)
        _stage2_update_uint64!(
            context,
            integer_buffer,
            reinterpret(UInt64, Int64(variant.likelihood_rank)),
        )
        _stage2_update_bytes!(context, integer_buffer, variant.sequence)
    end
    return bytes2hex(Mycelia.SHA.digest!(context))
end

"""
    infer_variant_abundances(variants, alignments, noise_likelihoods; config)

Attach the technology-agnostic soft-EM inference core to an existing Stage-2
candidate set. This bridge deliberately accepts calibrated log-likelihoods
rather than manufacturing them from a best-hit alignment score. Callers must
derive the likelihood records from an explicit read-error model.
"""
function infer_variant_abundances(
        variants::AbstractVector{RankedVariant},
        alignments::AbstractVector{ReadCandidateLikelihood},
        noise_likelihoods::AbstractVector{ReadNoiseLikelihood};
        config::CandidateInferenceConfig = CandidateInferenceConfig()
)::ProbabilisticVariantRankingResult
    isempty(variants) && throw(ArgumentError("variants must not be empty"))
    candidate_ids = [variant.id for variant in variants]
    length(unique(candidate_ids)) == length(candidate_ids) ||
        throw(ArgumentError("variant identifiers must be unique"))
    inference = infer_candidate_abundances(
        alignments,
        noise_likelihoods;
        candidate_ids,
        config,
    )
    inferred_by_id = Dict(candidate.candidate_id => candidate
        for candidate in inference.candidates)
    rankings = ProbabilisticVariantRanking[
        ProbabilisticVariantRanking(
            variant.id,
            variant.likelihood_rank,
            inferred_by_id[variant.id].rank,
            inferred_by_id[variant.id].abundance,
            inferred_by_id[variant.id].expected_read_count,
            inferred_by_id[variant.id].has_alignment_support,
        )
        for variant in variants
    ]
    sort!(rankings; by = ranking -> ranking.abundance_rank)
    likelihood_primary = only(filter(
        ranking -> ranking.likelihood_rank == 1, rankings))
    inferred_primary = inference.primary_candidate_id
    primary_agreement = inferred_primary !== nothing &&
                        likelihood_primary.candidate_id == inferred_primary
    primary_status = if inference.primary_status != :resolved
        inference.primary_status
    elseif primary_agreement
        :resolved
    else
        :rank_disagreement
    end
    return ProbabilisticVariantRankingResult(
        rankings,
        inference,
        primary_agreement,
        primary_status,
        _stage2_candidate_set_sha256(variants),
        inference.input_sha256,
    )
end

"""Write the auditable fitted-abundance ranking table."""
function write_probabilistic_variant_ranking(
        path::AbstractString,
        result::ProbabilisticVariantRankingResult
)::String
    table = DataFrames.DataFrame(
        candidate_id = [ranking.candidate_id for ranking in result.rankings],
        likelihood_rank = [ranking.likelihood_rank for ranking in result.rankings],
        abundance_rank = [ranking.abundance_rank for ranking in result.rankings],
        estimated_abundance = [ranking.abundance for ranking in result.rankings],
        expected_read_count = [ranking.expected_read_count for ranking in result.rankings],
        noise_abundance = fill(
            result.inference.noise_abundance, length(result.rankings)),
        noise_expected_read_count = fill(
            result.inference.noise_expected_read_count, length(result.rankings)),
        has_alignment_support = [
            ranking.has_alignment_support for ranking in result.rankings],
        primary_agreement = fill(result.primary_agreement, length(result.rankings)),
        primary_status = fill(String(result.primary_status), length(result.rankings)),
        inference_converged = fill(
            result.inference.converged, length(result.rankings)),
        inference_iterations = fill(
            result.inference.iterations, length(result.rankings)),
        final_log_likelihood = fill(
            result.inference.final_log_likelihood, length(result.rankings)),
        n_reads = fill(result.inference.n_reads, length(result.rankings)),
        n_alignments = fill(
            result.inference.n_alignments, length(result.rankings)),
        n_active_pairs = fill(
            result.inference.n_active_pairs, length(result.rankings)),
        max_iterations = fill(
            result.inference.config.max_iterations, length(result.rankings)),
        abundance_tolerance = fill(
            result.inference.config.abundance_tolerance, length(result.rankings)),
        primary_tolerance = fill(
            result.inference.config.primary_tolerance, length(result.rankings)),
        retain_full_trace = fill(
            result.inference.config.retain_full_trace, length(result.rankings)),
        trace_sha256 = fill(
            result.inference.trace_sha256, length(result.rankings)),
        candidate_set_sha256 = fill(
            result.candidate_set_sha256, length(result.rankings)),
        inference_input_sha256 = fill(
            result.inference_input_sha256, length(result.rankings)),
    )
    output = String(path)
    mkpath(dirname(output))
    mktemp(dirname(output)) do temporary_path, stream
        close(stream)
        CSV.write(temporary_path, table; delim = '\t')
        mv(temporary_path, output; force = true)
    end
    return output
end

"""
Persistent outputs from the Stage-2 metaFlye + Strainy route.
"""
struct Stage2DeconvolutionResult
    primary_fasta::String
    layout_assembly::String
    layout_gfa::String
    likelihood_ranked_variants_fasta::String
    abundance_ranked_variants_fasta::String
    ranking_tsv::String
    corrected_fastq::String
    variants::Vector{RankedVariant}
    provenance::Dict{String, Any}
end

function _prepare_stage2_inputs(
        reads::R,
        config::AssemblyConfig
)::NamedTuple where {R}
    stage1 = _run_stage1_correction(
        reads, config; materialize_corrected_reads = false)
    metadata = stage1.result_dict[:metadata]
    graph = get(stage1.result_dict, :final_graph, nothing)
    recorded_graph_k = get(metadata, :final_graph_k, stage1.max_k)
    graph_k = recorded_graph_k isa Int && recorded_graph_k > 0 ?
              recorded_graph_k : stage1.max_k
    reusable = get(metadata, :final_graph_reusable, false) === true
    graph_mode = get(metadata, :final_graph_mode, nothing)
    can_reuse_graph = graph !== nothing && reusable &&
                      recorded_graph_k == graph_k && graph_mode == :doublestrand
    rebuild_reasons = String[]
    graph === nothing && push!(rebuild_reasons, "final-graph-absent")
    reusable || push!(rebuild_reasons, "final-graph-not-byte-identical")
    recorded_graph_k == graph_k || push!(rebuild_reasons, "invalid-final-graph-k")
    graph_mode == :doublestrand || push!(rebuild_reasons, "non-doublestrand-final-graph")
    graph_source = can_reuse_graph ? "stage1-final" : "corrected-fastq-rebuild"
    if !can_reuse_graph
        graph = open(FASTX.FASTQ.Reader, stage1.corrected_fastq) do reader
            build_kmer_graph(
                collect(reader),
                graph_k;
                mode = :doublestrand,
                memory_profile = :ultralight_quality,
            )
        end
    end
    index = _transition_likelihood_index(graph, graph_k)
    return (;
        corrected_fastq = stage1.corrected_fastq,
        corrected_read_count = stage1.corrected_read_count,
        graph_source,
        graph_rebuild_reason =
            can_reuse_graph ? nothing : join(rebuild_reasons, ","),
        graph_k,
        index)
end

function _transition_likelihood_index(
        graph::MetaGraphsNext.MetaGraph,
        k::Int
)::TransitionLikelihoodIndex
    base_graph = graph.graph
    n_vertices = Graphs.nv(base_graph)
    n_vertices > 0 || error("cannot index an empty Stage-1 graph")
    n_vertices <= typemax(UInt32) ||
        error("Stage-1 graph has too many vertices for the compact index")
    first_code = first(Graphs.vertices(base_graph))
    first_label = MetaGraphsNext.label_for(graph, first_code)
    KmerType = typeof(first_label)
    KmerType <: Kmers.DNAKmer ||
        error("Stage-1 transition index requires DNA k-mer graph labels")
    length(first_label) == k ||
        error("Stage-1 graph label length differs from graph_k")
    vertex_ids = Dict{KmerType, UInt32}()
    sizehint!(vertex_ids, n_vertices)
    for code in Graphs.vertices(base_graph)
        vertex_ids[MetaGraphsNext.label_for(graph, code)] = UInt32(code)
    end
    outgoing = zeros(Float64, n_vertices)
    probabilities = Dict{UInt64, Float64}()
    sizehint!(probabilities, Graphs.ne(base_graph))
    for edge in Graphs.edges(base_graph)
        source_code = Graphs.src(edge)
        destination_code = Graphs.dst(edge)
        source_label = MetaGraphsNext.label_for(graph, source_code)
        destination_label = MetaGraphsNext.label_for(graph, destination_code)
        weight = _edge_transition_weight(graph[source_label, destination_label])
        isfinite(weight) && weight > 0.0 ||
            error("Stage-1 graph contains a nonpositive transition weight")
        source_id = UInt32(source_code)
        destination_id = UInt32(destination_code)
        outgoing[Int(source_id)] += weight
        probabilities[_stage2_edge_key(source_id, destination_id)] = weight
    end

    for (edge, weight) in probabilities
        source_id = UInt32(edge >> 32)
        total = outgoing[Int(source_id)]
        total > 0.0 || continue
        probabilities[edge] = log2(weight / total)
    end
    return TransitionLikelihoodIndex(k, vertex_ids, probabilities)
end

function _score_variant_sequence(
        sequence::AbstractString,
        index::TransitionLikelihoodIndex
)::NamedTuple{(:mean_log2_probability, :scored_fraction), Tuple{Float64, Float64}}
    return _score_variant_sequence(sequence, index, Val(index.k))
end

function _score_variant_sequence(
        sequence::AbstractString,
        index::TransitionLikelihoodIndex{KmerType},
        ::Val{K}
)::NamedTuple{(:mean_log2_probability, :scored_fraction), Tuple{Float64, Float64}} where {
        KmerType <: Kmers.DNAKmer, K}
    dna_sequence = try
        BioSequences.LongDNA{4}(sequence)
    catch exception
        throw(ArgumentError(
            "variant sequence is not valid DNA: $(sprint(showerror, exception))"))
    end
    n_transitions = length(dna_sequence) - K
    n_transitions > 0 || return (mean_log2_probability = -Inf, scored_fraction = 0.0)

    score = 0.0
    n_scored = 0
    previous_kmer = nothing
    previous_position = 0
    for (kmer, position) in Kmers.UnambiguousDNAMers{K}(dna_sequence)
        if previous_kmer !== nothing && position == previous_position + 1
            source_id = get(index.vertex_ids, previous_kmer, UInt32(0))
            destination_id = get(index.vertex_ids, kmer, UInt32(0))
            edge_score = get(
                index.log2_probability,
                _stage2_edge_key(source_id, destination_id),
                nothing,
            )
            if source_id != 0 && destination_id != 0 && edge_score !== nothing
                score += edge_score
                n_scored += 1
            end
        end
        previous_kmer = kmer
        previous_position = position
    end
    fraction = n_scored / n_transitions
    mean_score = n_scored == n_transitions ? score / n_transitions : -Inf
    return (mean_log2_probability = mean_score, scored_fraction = fraction)
end

function _read_gfa_records(path::String)::Dict{String, NamedTuple}
    isfile(path) || error("GFA file does not exist: $(path)")
    segments = Dict{String, NamedTuple}()
    open(path, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            length(fields) >= 3 || continue
            fields[1] == "S" || continue
            fields[3] == "*" && continue
            coverage = nothing
            for tag in fields[4:end]
                startswith(tag, "dp:") || continue
                parsed = tryparse(Float64, split(tag, ':'; limit = 3)[3])
                parsed === nothing || (coverage = parsed)
            end
            segments[String(fields[2])] = (;
                sequence = String(fields[3]), coverage)
        end
    end
    isempty(segments) && error("GFA contains no sequence-bearing segments: $(path)")
    return segments
end

function _read_gfa_segments(path::String)::Dict{String, String}
    return Dict(id => record.sequence for (id, record) in _read_gfa_records(path))
end

function _find_column(
        table::DataFrames.DataFrame,
        candidates::Vector{String},
        purpose::String
)::String
    lookup = Dict(
        lowercase(String(name)) => String(name)
        for name in DataFrames.names(table)
    )
    for candidate in candidates
        haskey(lookup, candidate) && return lookup[candidate]
    end
    error("Strainy phased-unitig table lacks a $(purpose) column; " *
          "columns=$(join(DataFrames.names(table), ", "))")
end

function _strainy_coverages(path::String)::Dict{String, Float64}
    isfile(path) || error("Strainy phased-unitig table does not exist: $(path)")
    table = CSV.read(path, DataFrames.DataFrame)
    id_column = _find_column(table,
        ["unitig", "unitig_id", "strain_unitig", "id", "name"], "unitig identifier")
    coverage_column = _find_column(table,
        ["coverage", "mean_coverage", "median_coverage", "cov"], "coverage")
    coverages = Dict{String, Float64}()
    for row in DataFrames.eachrow(table)
        id = String(row[id_column])
        coverage = tryparse(Float64, string(row[coverage_column]))
        coverage === nothing && error("non-numeric Strainy coverage for $(id)")
        coverages[id] = coverage
    end
    return coverages
end

"""
Legacy first-slice best-hit support used by the existing production route.

This is retained only until a calibrated read-candidate likelihood builder is
wired to [`infer_variant_abundances`](@ref). It must not be described as
probabilistic abundance inference.
"""
function _read_support_coverages(
        records::Dict{String, NamedTuple},
        fastq::String,
        threads::Int
)::Dict{String, Float64}
    isfile(fastq) || error("support FASTQ does not exist: $(fastq)")
    mktempdir() do directory
        candidates = joinpath(directory, "candidates.fasta")
        open(candidates, "w") do io
            for (id, record) in sort(collect(records); by = first)
                println(io, ">$(id)")
                println(io, record.sequence)
            end
        end
        command = `$(Mycelia.CONDA_RUNNER) run -n strainy minimap2 -x map-ont
            -t $(threads) $(candidates) $(fastq)`
        paf = read(command, String)
        best = Dict{String, Tuple{Int, String, Int}}()
        for line in eachline(IOBuffer(paf))
            fields = split(line, '\t')
            length(fields) >= 12 || continue
            query = String(fields[1])
            target = String(fields[6])
            matches = something(tryparse(Int, fields[10]), 0)
            aligned = something(tryparse(Int, fields[11]), 0)
            current = get(best, query, (-1, "", 0))
            (matches, aligned, target) > (current[1], current[3], current[2]) || continue
            best[query] = (matches, target, aligned)
        end
        aligned_bases = Dict(id => 0.0 for id in keys(records))
        for (_, target, aligned) in values(best)
            aligned_bases[target] = get(aligned_bases, target, 0.0) + aligned
        end
        return Dict(id => aligned_bases[id] / length(record.sequence)
            for (id, record) in records)
    end
end

function _variant_role(rank::Int)::Symbol
    rank == 1 && return :primary
    rank == 2 && return :secondary
    rank == 3 && return :tertiary
    return :alternate
end

function _rank_stage2_variants(
        strain_unitigs_gfa::String,
        phased_unitig_info::String,
        index::TransitionLikelihoodIndex;
        max_variants::Int = 3,
        min_scored_fraction::Float64 = 0.9,
        support_fastq::Union{Nothing, String} = nothing,
        threads::Int = Mycelia.get_default_threads()
)::Vector{RankedVariant}
    max_variants > 0 || throw(ArgumentError("max_variants must be positive"))
    0.0 <= min_scored_fraction <= 1.0 ||
        throw(ArgumentError("min_scored_fraction must be between 0 and 1"))

    segments = _read_gfa_records(strain_unitigs_gfa)
    coverages = _strainy_coverages(phased_unitig_info)
    support_coverages = support_fastq === nothing ? Dict{String, Float64}() :
                        _read_support_coverages(segments, support_fastq, threads)
    candidates = NamedTuple[]
    for (id, record) in segments
        sequence = record.sequence
        coverage = haskey(support_coverages, id) ? support_coverages[id] :
                   something(record.coverage, get(coverages, id, nothing))
        coverage === nothing && continue
        score = _score_variant_sequence(sequence, index)
        score.scored_fraction >= min_scored_fraction || continue
        isfinite(score.mean_log2_probability) || continue
        push!(candidates, (;
            id,
            sequence,
            graph_mean_log2_probability = score.mean_log2_probability,
            scored_transition_fraction = score.scored_fraction,
            strainy_coverage = Float64(coverage)
        ))
    end
    likelihood_order = sort(
        copy(candidates);
        by = candidate -> (
            -candidate.graph_mean_log2_probability,
            -candidate.strainy_coverage,
            -length(candidate.sequence),
            candidate.id,
        ),
    )
    abundance_order = sort(
        copy(candidates);
        by = candidate -> (
            -candidate.strainy_coverage,
            -candidate.graph_mean_log2_probability,
            -length(candidate.sequence),
            candidate.id,
        ),
    )
    length(likelihood_order) >= max_variants ||
        error(
            "Stage-2 produced only $(length(likelihood_order)) " *
            "graph-supported variants; $(max_variants) required")
    likelihood_rank = Dict(candidate.id => rank
        for (rank, candidate) in enumerate(likelihood_order))
    abundance_rank = Dict(candidate.id => rank
        for (rank, candidate) in enumerate(abundance_order))
    selected_ids = union(
        Set(candidate.id for candidate in likelihood_order[1:max_variants]),
        Set(candidate.id for candidate in abundance_order[1:max_variants]))
    selected = filter(candidate -> candidate.id in selected_ids, likelihood_order)
    return [RankedVariant(likelihood_rank[candidate.id], abundance_rank[candidate.id],
                _variant_role(likelihood_rank[candidate.id]), candidate.id,
                candidate.sequence, candidate.graph_mean_log2_probability,
                candidate.scored_transition_fraction, candidate.strainy_coverage)
            for candidate in selected]
end

function _write_ranked_variants(
        variants::Vector{RankedVariant},
        output_dir::String,
        max_variants::Int
)::NamedTuple
    mkpath(output_dir)
    likelihood_primary = only(filter(variant -> variant.likelihood_rank == 1, variants))
    abundance_primary = only(filter(variant -> variant.abundance_rank == 1, variants))
    primary_status = likelihood_primary.id == abundance_primary.id ?
                     :resolved : :rank_disagreement
    primary = joinpath(output_dir, "primary_consensus.fasta")
    likelihood_fasta = joinpath(output_dir, "ranked_variants_likelihood.fasta")
    abundance_fasta = joinpath(output_dir, "ranked_variants_abundance.fasta")
    tsv = joinpath(output_dir, "ranked_variants.tsv")
    open(likelihood_fasta, "w") do io
        for variant in filter(variant -> variant.likelihood_rank <= max_variants,
                sort(variants; by = variant -> variant.likelihood_rank))
            println(io, ">likelihood_rank=$(variant.likelihood_rank)|id=$(variant.id)")
            println(io, variant.sequence)
        end
    end
    open(abundance_fasta, "w") do io
        for variant in filter(variant -> variant.abundance_rank <= max_variants,
                sort(variants; by = variant -> variant.abundance_rank))
            println(io, ">abundance_rank=$(variant.abundance_rank)|id=$(variant.id)")
            println(io, variant.sequence)
        end
    end
    open(primary, "w") do io
        variant = likelihood_primary
        println(
            io,
            ">primary|id=$(variant.id)|status=$(primary_status)|" *
            "abundance_primary_id=$(abundance_primary.id)",
        )
        println(io, variant.sequence)
    end
    table = DataFrames.DataFrame(
        likelihood_rank = [variant.likelihood_rank for variant in variants],
        abundance_rank = [variant.abundance_rank for variant in variants],
        role = [String(variant.role) for variant in variants],
        id = [variant.id for variant in variants],
        length = [length(variant.sequence) for variant in variants],
        graph_mean_log2_probability =
            [variant.graph_mean_log2_probability for variant in variants],
        scored_transition_fraction =
            [variant.scored_transition_fraction for variant in variants],
        strainy_coverage = [variant.strainy_coverage for variant in variants],
        primary_status = fill(String(primary_status), length(variants)),
        likelihood_primary_id = fill(likelihood_primary.id, length(variants)),
        abundance_primary_id = fill(abundance_primary.id, length(variants)),
    )
    CSV.write(tsv, table; delim = '\t')
    return (;
        primary,
        likelihood_fasta,
        abundance_fasta,
        tsv,
        primary_status,
        likelihood_primary_id = likelihood_primary.id,
        abundance_primary_id = abundance_primary.id,
    )
end

"""
    deconvolve_stage2(reads, config::AssemblyConfig; kwargs...)

Run Stage-1 graph accurization once, reuse metaFlye for long-read repeat-graph
layout, and reuse Strainy for strain-haplotype deconvolution. Candidate
haplotypes are ranked by normalized likelihood under the accurized graph. The
production route still uses explicitly labeled legacy best-hit abundance until
a calibrated read-candidate likelihood builder is connected to
[`infer_variant_abundances`](@ref); the soft-EM implementation in this slice is
an inference seam, not yet the production abundance method.
"""
function deconvolve_stage2(
        reads::R,
        config::AssemblyConfig;
        genome_size::Union{Nothing, Integer, AbstractString} = nothing,
        min_overlap::Union{Nothing, Int} = nothing,
        layout_iterations::Int = 0,
        max_variants::Int = 3,
        min_scored_fraction::Float64 = 0.9,
        threads::Int = Mycelia.get_default_threads()
)::Stage2DeconvolutionResult where {R}
    config.corrector == :iterative ||
        throw(ArgumentError("Stage-2 requires corrector=:iterative"))
    config.sequencing_tech == :nanopore ||
        throw(ArgumentError("the first Stage-2 slice requires sequencing_tech=:nanopore"))
    config.output_dir === nothing &&
        throw(ArgumentError("Stage-2 requires a persistent AssemblyConfig.output_dir"))
    layout_iterations >= 0 ||
        throw(ArgumentError("layout_iterations must be nonnegative"))

    output_dir = config.output_dir
    layout_dir = joinpath(output_dir, "metaflye")
    strainy_dir = joinpath(output_dir, "strainy")
    ranked_dir = joinpath(output_dir, "ranked")
    any(isdir, (layout_dir, strainy_dir, ranked_dir)) &&
        error("Stage-2 output directory already contains a prior run: $(output_dir)")

    prepared = _prepare_stage2_inputs(reads, config)
    (; corrected_fastq, graph_k, index) = prepared

    layout = Mycelia.run_metaflye(;
        fastq = corrected_fastq,
        outdir = layout_dir,
        genome_size = genome_size,
        read_type = "nano-corr",
        meta = true,
        min_overlap = min_overlap,
        iterations = layout_iterations,
        keep_haplotypes = true,
        no_alt_contigs = true,
        threads = threads)
    isfile(layout.assembly) || error("metaFlye did not produce assembly.fasta")
    isfile(layout.graph) || error("metaFlye did not produce assembly_graph.gfa")

    phased = Mycelia.run_strainy(;
        gfa_ref = layout.graph,
        fastq = corrected_fastq,
        outdir = strainy_dir,
        read_mode = :nano,
        stage = :e2e,
        min_unitig_coverage = 3,
        unitig_split_length = 0,
        threads = threads)
    variants = _rank_stage2_variants(
        phased.strain_contigs_gfa, phased.phased_unitig_info, index;
        max_variants = max_variants,
        min_scored_fraction = min_scored_fraction,
        support_fastq = corrected_fastq,
        threads = threads)
    ranked = _write_ranked_variants(variants, ranked_dir, max_variants)

    provenance = Dict{String, Any}(
        "stage1_corrector" => "iterative",
        "stage1_strategy" => String(config.strategy),
        "sequencing_tech" => "nanopore",
        "layout_tool" => "metaflye",
        "layout_read_type" => "nano-corr",
        "layout_iterations" => layout_iterations,
        "deconvolution_tool" => "strainy",
        "abundance_method" => "legacy-best-hit-minimap2",
        "probabilistic_abundance_status" => "available-requires-calibrated-likelihoods",
        "graph_k" => graph_k,
        "graph_source" => prepared.graph_source,
        "graph_rebuild_reason" => prepared.graph_rebuild_reason,
        "corrected_read_count" => prepared.corrected_read_count,
        "max_variants" => max_variants,
        "min_scored_fraction" => min_scored_fraction
    )
    provenance["primary_status"] = String(ranked.primary_status)
    provenance["likelihood_primary_id"] = ranked.likelihood_primary_id
    provenance["abundance_primary_id"] = ranked.abundance_primary_id
    return Stage2DeconvolutionResult(
        ranked.primary,
        layout.assembly,
        layout.graph,
        ranked.likelihood_fasta,
        ranked.abundance_fasta,
        ranked.tsv,
        corrected_fastq,
        variants,
        provenance)
end
