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
        log2_probability = Float64(probability)
        isfinite(log2_probability) || throw(ArgumentError(
            "transition log2 probabilities must be finite"))
        log2_probability <= 0.0 || throw(ArgumentError(
            "transition log2 probabilities must be nonpositive"))
        source_kmer = KmerType(source)
        destination_kmer = KmerType(destination)
        source_id = get!(
            vertex_ids, source_kmer, UInt32(length(vertex_ids) + 1))
        destination_id = get!(
            vertex_ids, destination_kmer, UInt32(length(vertex_ids) + 1))
        edge = _stage2_edge_key(source_id, destination_id)
        haskey(encoded_probabilities, edge) && throw(ArgumentError(
            "duplicate transition after DNA k-mer normalization: " *
            "$(source) -> $(destination)"))
        encoded_probabilities[edge] = log2_probability
    end
    return TransitionLikelihoodIndex(k, vertex_ids, encoded_probabilities)
end

"""
One ranked experimental Strainy GFA segment candidate.

This record does not represent a complete path, consensus genome, or inferred
haplotype. Path grouping is a blocked follow-on.
"""
struct RankedSegmentCandidate
    likelihood_rank::Int
    abundance_rank::Int
    role::Symbol
    id::String
    sequence::String
    graph_mean_log2_probability::Float64
    scored_transition_fraction::Float64
    strainy_coverage::Float64
end

struct GfaSegmentRecord
    sequence::String
    coverage::Union{Nothing, Float64}
end

"""Fitted soft-mixture abundance attached to one Stage-2 candidate."""
struct ProbabilisticSegmentCandidateRanking
    candidate_id::String
    likelihood_rank::Int
    abundance_rank::Int
    abundance::Float64
    expected_read_count::Float64
    has_alignment_support::Bool
end

"""Soft abundance inference and the resulting Stage-2 primary agreement."""
struct ProbabilisticSegmentCandidateRankingResult
    rankings::Vector{ProbabilisticSegmentCandidateRanking}
    inference::CandidateInferenceResult
    primary_agreement::Bool
    primary_status::Symbol
    candidate_set_sha256::String
    inference_input_sha256::String
end

function _stage2_candidate_set_sha256(
        candidates::AbstractVector{RankedSegmentCandidate}
)::String
    context = Mycelia.SHA.SHA2_256_CTX()
    integer_buffer = Vector{UInt8}(undef, 8)
    tag_buffer = Vector{UInt8}(undef, 1)
    order = sortperm(eachindex(candidates); by = index -> candidates[index].id)
    _stage2_update_uint64!(context, integer_buffer, UInt64(length(order)))
    for index in order
        candidate = candidates[index]
        _stage2_update_tag!(context, tag_buffer, 0x20)
        _stage2_update_bytes!(context, integer_buffer, candidate.id)
        _stage2_update_uint64!(
            context,
            integer_buffer,
            reinterpret(UInt64, Int64(candidate.likelihood_rank)),
        )
        _stage2_update_bytes!(context, integer_buffer, candidate.sequence)
    end
    return bytes2hex(Mycelia.SHA.digest!(context))
end

"""
    infer_segment_candidate_abundances(
        candidates, alignments, noise_likelihoods; config)

Attach the technology-agnostic soft-EM inference core to an existing Stage-2
candidate set. This bridge deliberately accepts calibrated log-likelihoods
rather than manufacturing them from a best-hit alignment score. Callers must
derive the likelihood records from an explicit read-error model.
"""
function infer_segment_candidate_abundances(
        candidates::AbstractVector{RankedSegmentCandidate},
        alignments::AbstractVector{ReadCandidateLikelihood},
        noise_likelihoods::AbstractVector{ReadNoiseLikelihood};
        config::CandidateInferenceConfig = CandidateInferenceConfig()
)::ProbabilisticSegmentCandidateRankingResult
    isempty(candidates) && throw(
        ArgumentError("segment candidates must not be empty"))
    candidate_ids = [candidate.id for candidate in candidates]
    length(unique(candidate_ids)) == length(candidate_ids) ||
        throw(ArgumentError("segment-candidate identifiers must be unique"))
    inference = infer_candidate_abundances(
        alignments,
        noise_likelihoods;
        candidate_ids,
        config,
    )
    inferred_by_id = Dict(candidate.candidate_id => candidate
        for candidate in inference.candidates)
    rankings = ProbabilisticSegmentCandidateRanking[
        ProbabilisticSegmentCandidateRanking(
            candidate.id,
            candidate.likelihood_rank,
            inferred_by_id[candidate.id].rank,
            inferred_by_id[candidate.id].abundance,
            inferred_by_id[candidate.id].expected_read_count,
            inferred_by_id[candidate.id].has_alignment_support,
        )
        for candidate in candidates
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
    return ProbabilisticSegmentCandidateRankingResult(
        rankings,
        inference,
        primary_agreement,
        primary_status,
        _stage2_candidate_set_sha256(candidates),
        inference.input_sha256,
    )
end

"""Write the auditable fitted-abundance ranking table."""
function write_probabilistic_segment_candidate_ranking(
        path::AbstractString,
        result::ProbabilisticSegmentCandidateRankingResult
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
Persistent experimental segment-candidate outputs from the Stage-2 route.

The result deliberately exposes no complete-path or haplotype claim.
"""
struct Stage2SegmentCandidateResult
    primary_segment_candidate_fasta::String
    layout_assembly::String
    layout_gfa::String
    likelihood_ranked_segment_candidates_fasta::String
    abundance_ranked_segment_candidates_fasta::String
    segment_candidate_ranking_tsv::String
    corrected_fastq::String
    segment_candidates::Vector{RankedSegmentCandidate}
    provenance::Dict{String, Any}
end

"""
Atomically reserve a unique Stage-2 attempt directory under a caller-owned root.

Each production attempt keeps its Stage-1, layout, phasing, and ranking artifacts
in this directory. Failed attempts are intentionally preserved for inspection;
concurrent attempts never share paths.
"""
function _reserve_stage2_attempt_directory(
        output_root::AbstractString
)::NamedTuple
    isempty(strip(output_root)) && error("Stage-2 output root must be nonempty")
    root = abspath(normpath(String(output_root)))
    islink(root) && error("Stage-2 output root must not be a symbolic link")
    mkpath(root)
    islink(root) && error("Stage-2 output root must not be a symbolic link")
    attempt_output_dir = mktempdir(
        root; prefix = "attempt-", cleanup = false)
    return (;
        attempt_id = basename(attempt_output_dir),
        attempt_output_dir,
        output_root = root,
        corrected_fastq = joinpath(attempt_output_dir, "corrected.fastq"),
        layout_dir = joinpath(attempt_output_dir, "metaflye"),
        strainy_dir = joinpath(attempt_output_dir, "strainy"),
        ranked_dir = joinpath(attempt_output_dir, "ranked"),
    )
end

function _stage2_attempt_provenance(
        attempt::NamedTuple
)::Dict{String, Any}
    attempt_id = String(attempt.attempt_id)
    isempty(attempt_id) && error("Stage-2 attempt ID must be nonempty")
    relative_path = relpath(
        String(attempt.attempt_output_dir), String(attempt.output_root))
    relative_path == attempt_id || error(
        "Stage-2 attempt path is not bound to its attempt ID")
    return Dict{String, Any}(
        "stage2_attempt_id" => attempt_id,
        "stage2_attempt_relative_path" => relative_path,
        "stage2_attempt_path_semantics" => "relative-to-config-output-dir",
    )
end

function _release_unreusable_stage2_graph!(
        result_dict::AbstractDict
)::Nothing
    pop!(result_dict, :final_graph, nothing)
    return nothing
end

function _merge_stage2_graph_chunk(
        graph::Union{Nothing, MetaGraphsNext.MetaGraph},
        records::Vector{FASTX.FASTQ.Record},
        k::Int
)::MetaGraphsNext.MetaGraph
    isempty(records) && error("Stage-2 graph chunk must not be empty")
    chunk_graph = build_kmer_graph(
        records,
        k;
        mode = :singlestrand,
        memory_profile = :ultralight_quality,
    )
    if graph === nothing
        return chunk_graph
    end
    _merge_reduced_graphs!(graph, chunk_graph, :ultralight_quality)
    return graph
end

"""
Build the Stage-2 fallback graph from a FASTQ with bounded read materialization.

Only `chunk_size` records are held in addition to the graph. `chunk_observer`
receives each live batch size as a regression seam. Chunk graphs are merged in
input order before one doublestrand conversion, matching the eager
`build_kmer_graph(...; mode=:doublestrand, memory_profile=:ultralight_quality)`
semantics without `collect(reader)`.
"""
function _build_stage2_graph_from_fastq(
        fastq_path::AbstractString,
        k::Int;
        chunk_size::Int = 1024,
        chunk_observer::Function = _ -> nothing,
)::MetaGraphsNext.MetaGraph
    isfile(fastq_path) || error("corrected FASTQ does not exist: $(fastq_path)")
    chunk_size > 0 || throw(ArgumentError("chunk_size must be positive"))
    graph::Union{Nothing, MetaGraphsNext.MetaGraph} = nothing
    chunk = FASTX.FASTQ.Record[]
    open(FASTX.FASTQ.Reader, fastq_path) do reader
        for record in reader
            push!(chunk, record)
            if length(chunk) == chunk_size
                chunk_observer(length(chunk))
                graph = _merge_stage2_graph_chunk(graph, chunk, k)
                empty!(chunk)
            end
        end
    end
    if !isempty(chunk)
        chunk_observer(length(chunk))
        graph = _merge_stage2_graph_chunk(graph, chunk, k)
    end
    graph === nothing && error("corrected FASTQ contains no reads: $(fastq_path)")
    isempty(MetaGraphsNext.labels(graph)) && error(
        "corrected FASTQ contains no k-mers at k=$(k): $(fastq_path)")
    return convert_to_doublestrand(graph)
end

function _prepare_stage2_inputs(
        reads::R,
        config::AssemblyConfig;
        persistent_output_dir::Union{Nothing, AbstractString} = config.output_dir
)::NamedTuple where {R}
    stage1 = _run_stage1_correction(
        reads,
        config;
        materialize_corrected_reads = false,
        persistent_output_dir,
    )
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
        # A retained Stage-1 graph can be the dominant allocation. Remove every
        # strong reference before rebuilding so the streamed fallback never
        # co-resides with an obsolete full graph at scale.
        _release_unreusable_stage2_graph!(stage1.result_dict)
        graph = nothing
        GC.gc()
        graph = _build_stage2_graph_from_fastq(
            stage1.corrected_fastq, graph_k)
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

function _score_segment_candidate_sequence(
        sequence::AbstractString,
        index::TransitionLikelihoodIndex
)::NamedTuple{(:mean_log2_probability, :scored_fraction), Tuple{Float64, Float64}}
    return _score_segment_candidate_sequence(sequence, index, Val(index.k))
end

function _score_segment_candidate_sequence(
        sequence::AbstractString,
        index::TransitionLikelihoodIndex{KmerType},
        ::Val{K}
)::NamedTuple{(:mean_log2_probability, :scored_fraction), Tuple{Float64, Float64}} where {
        KmerType <: Kmers.DNAKmer, K}
    dna_sequence = try
        BioSequences.LongDNA{4}(sequence)
    catch exception
        throw(ArgumentError(
            "segment-candidate sequence is not valid DNA: " *
            sprint(showerror, exception)))
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
    mean_score = n_scored > 0 ? score / n_scored : -Inf
    return (mean_log2_probability = mean_score, scored_fraction = fraction)
end

function _read_gfa_records(path::String)::Dict{String, GfaSegmentRecord}
    isfile(path) || error("GFA file does not exist: $(path)")
    segments = Dict{String, GfaSegmentRecord}()
    identifiers = Set{String}()
    open(path, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            fields = split(line, '\t'; keepempty = true)
            isempty(fields) && continue
            fields[1] == "S" || continue
            length(fields) >= 3 || error(
                "malformed GFA S record at line $(line_number): " *
                "expected at least 3 tab-delimited fields")
            identifier = String(fields[2])
            isempty(identifier) && error(
                "malformed GFA S record at line $(line_number): empty identifier")
            identifier in identifiers && error(
                "duplicate GFA segment identifier at line $(line_number): " *
                identifier)
            push!(identifiers, identifier)
            fields[3] == "*" && continue
            sequence = String(fields[3])
            isempty(sequence) && error(
                "malformed GFA S record at line $(line_number): empty sequence")
            try
                BioSequences.LongDNA{4}(sequence)
            catch exception
                error(
                    "invalid DNA in GFA segment $(identifier) at line " *
                    "$(line_number): $(sprint(showerror, exception))",
                )
            end
            coverage = nothing
            for tag in fields[4:end]
                startswith(tag, "dp:") || continue
                coverage === nothing || error(
                    "duplicate dp coverage tag for GFA segment $(identifier)")
                components = split(tag, ':'; limit = 3, keepempty = true)
                length(components) == 3 && components[2] in ("i", "f") || error(
                    "malformed dp coverage tag for GFA segment $(identifier): " *
                    tag)
                parsed = tryparse(Float64, components[3])
                parsed === nothing && error(
                    "non-numeric dp coverage for GFA segment $(identifier): " *
                    components[3])
                isfinite(parsed) && parsed >= 0.0 || error(
                    "GFA segment coverage must be finite and nonnegative for " *
                    identifier)
                coverage = parsed
            end
            segments[identifier] = GfaSegmentRecord(sequence, coverage)
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
        isempty(id) && error("empty Strainy unitig identifier")
        haskey(coverages, id) && error("duplicate Strainy coverage row for $(id)")
        coverage = tryparse(Float64, string(row[coverage_column]))
        coverage === nothing && error("non-numeric Strainy coverage for $(id)")
        isfinite(coverage) && coverage >= 0.0 || error(
            "Strainy coverage must be finite and nonnegative for $(id)")
        coverages[id] = coverage
    end
    return coverages
end

const _Stage2PafBestHit = NamedTuple{
    (:matches, :target, :aligned), Tuple{Int, String, Int}}

function _parse_stage2_paf_integer(
        value::AbstractString,
        field_name::String,
        line_number::Int
)::Int
    parsed = tryparse(Int, value)
    parsed === nothing && error(
        "malformed PAF row $(line_number): $(field_name) is not an integer")
    return parsed
end

function _read_stage2_paf_best_hits(
        io::IO,
        records::AbstractDict{String, GfaSegmentRecord}
)::Dict{String, _Stage2PafBestHit}
    best = Dict{String, _Stage2PafBestHit}()
    for (line_number, line) in enumerate(eachline(io))
        isempty(strip(line)) && continue
        fields = split(line, '\t'; keepempty = true)
        length(fields) >= 12 || error(
            "malformed PAF row $(line_number): expected at least 12 fields")
        query = String(fields[1])
        target = String(fields[6])
        isempty(query) && error(
            "malformed PAF row $(line_number): empty query identifier")
        isempty(target) && error(
            "malformed PAF row $(line_number): empty target identifier")
        haskey(records, target) || error(
            "PAF row $(line_number) references unknown target: $(target)")

        query_length = _parse_stage2_paf_integer(fields[2], "query length", line_number)
        query_start = _parse_stage2_paf_integer(fields[3], "query start", line_number)
        query_end = _parse_stage2_paf_integer(fields[4], "query end", line_number)
        target_length =
            _parse_stage2_paf_integer(fields[7], "target length", line_number)
        target_start =
            _parse_stage2_paf_integer(fields[8], "target start", line_number)
        target_end = _parse_stage2_paf_integer(fields[9], "target end", line_number)
        matches = _parse_stage2_paf_integer(fields[10], "matches", line_number)
        aligned =
            _parse_stage2_paf_integer(fields[11], "alignment length", line_number)
        mapq = _parse_stage2_paf_integer(fields[12], "mapping quality", line_number)

        query_length > 0 && 0 <= query_start < query_end <= query_length || error(
            "malformed PAF row $(line_number): invalid query coordinates")
        fields[5] in ("+", "-") || error(
            "malformed PAF row $(line_number): invalid strand")
        target_length == length(records[target].sequence) || error(
            "malformed PAF row $(line_number): target length disagrees with " *
            "candidate sequence")
        0 <= target_start < target_end <= target_length || error(
            "malformed PAF row $(line_number): invalid target coordinates")
        matches > 0 && aligned > 0 && matches <= aligned || error(
            "malformed PAF row $(line_number): support must be positive, " *
            "finite, and no greater than alignment length")
        0 <= mapq <= 255 || error(
            "malformed PAF row $(line_number): invalid mapping quality")

        current = get(best, query, (matches = -1, target = "", aligned = 0))
        (matches, aligned, target) >
        (current.matches, current.aligned, current.target) || continue
        best[query] = (; matches, target, aligned)
    end
    return best
end

function _stage2_support_coverages(
        records::AbstractDict{String, GfaSegmentRecord},
        best::AbstractDict{String, _Stage2PafBestHit}
)::Dict{String, Float64}
    aligned_bases = Dict(id => 0.0 for id in keys(records))
    for hit in values(best)
        haskey(records, hit.target) || error(
            "support record references unknown target: $(hit.target)")
        hit.matches > 0 && hit.aligned > 0 || error(
            "mapped support must be finite and positive")
        aligned_bases[hit.target] += hit.aligned
    end
    total_support = sum(values(aligned_bases))
    isfinite(total_support) && total_support > 0.0 || error(
        "minimap2 produced zero total mapped support for segment candidates")
    coverages = Dict(id => aligned_bases[id] / length(record.sequence)
        for (id, record) in records)
    all(coverage -> isfinite(coverage) && coverage >= 0.0,
        values(coverages)) || error(
        "minimap2 produced invalid segment-candidate support")
    return coverages
end

"""
Legacy first-slice best-hit support used by the existing production route.

This is retained only until a calibrated read-candidate likelihood builder is
wired to [`infer_segment_candidate_abundances`](@ref). It must not be described as
probabilistic abundance inference.
"""
function _read_support_coverages(
        records::AbstractDict{String, GfaSegmentRecord},
        fastq::String,
        threads::Int
)::Dict{String, Float64}
    isfile(fastq) || error("support FASTQ does not exist: $(fastq)")
    threads > 0 || throw(ArgumentError("threads must be positive"))
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
        best = open(command, "r") do paf_stream
            _read_stage2_paf_best_hits(paf_stream, records)
        end
        return _stage2_support_coverages(records, best)
    end
end

function _segment_candidate_role(rank::Int)::Symbol
    rank == 1 && return :primary
    rank == 2 && return :secondary
    rank == 3 && return :tertiary
    return :alternate
end

function _rank_stage2_segment_candidates(
        strain_unitigs_gfa::String,
        phased_unitig_info::String,
        index::TransitionLikelihoodIndex;
        max_variants::Union{Nothing, Int} = nothing,
        min_scored_fraction::Float64 = 0.9,
        support_fastq::Union{Nothing, String} = nothing,
        threads::Int = Mycelia.get_default_threads()
)::Vector{RankedSegmentCandidate}
    max_variants === nothing || max_variants > 0 || throw(
        ArgumentError("max_variants diagnostic output cap must be positive"))
    0.0 <= min_scored_fraction <= 1.0 ||
        throw(ArgumentError("min_scored_fraction must be between 0 and 1"))

    segments = _read_gfa_records(strain_unitigs_gfa)
    coverages = _strainy_coverages(phased_unitig_info)
    support_coverages = support_fastq === nothing ? Dict{String, Float64}() :
                        _read_support_coverages(segments, support_fastq, threads)
    candidate_type = NamedTuple{
        (:id, :sequence, :graph_mean_log2_probability,
            :scored_transition_fraction, :strainy_coverage),
        Tuple{String, String, Float64, Float64, Float64},
    }
    candidates = candidate_type[]
    for (id, record) in segments
        sequence = record.sequence
        coverage = haskey(support_coverages, id) ? support_coverages[id] :
                   something(record.coverage, get(coverages, id, nothing))
        coverage === nothing && continue
        score = _score_segment_candidate_sequence(sequence, index)
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
    isempty(likelihood_order) && error(
        "Stage-2 produced no graph-supported experimental segment candidates")
    likelihood_rank = Dict(candidate.id => rank
        for (rank, candidate) in enumerate(likelihood_order))
    abundance_rank = Dict(candidate.id => rank
        for (rank, candidate) in enumerate(abundance_order))
    output_count = max_variants === nothing ? length(likelihood_order) :
                   min(max_variants, length(likelihood_order))
    selected_ids = union(
        Set(candidate.id for candidate in likelihood_order[1:output_count]),
        Set(candidate.id for candidate in abundance_order[1:output_count]))
    selected = filter(candidate -> candidate.id in selected_ids, likelihood_order)
    return [RankedSegmentCandidate(
                likelihood_rank[candidate.id],
                abundance_rank[candidate.id],
                _segment_candidate_role(likelihood_rank[candidate.id]),
                candidate.id,
                candidate.sequence,
                candidate.graph_mean_log2_probability,
                candidate.scored_transition_fraction,
                candidate.strainy_coverage,
            )
            for candidate in selected]
end

function _write_ranked_variants(
        candidates::Vector{RankedSegmentCandidate},
        output_dir::String,
        max_variants::Union{Nothing, Int}
)::NamedTuple
    mkpath(output_dir)
    likelihood_primary = only(filter(
        candidate -> candidate.likelihood_rank == 1, candidates))
    abundance_primary = only(filter(
        candidate -> candidate.abundance_rank == 1, candidates))
    primary_status = likelihood_primary.id == abundance_primary.id ?
                     :resolved : :rank_disagreement
    primary_path = joinpath(output_dir, "primary_segment_candidate.fasta")
    likelihood_fasta =
        joinpath(output_dir, "ranked_segment_candidates_likelihood.fasta")
    abundance_fasta =
        joinpath(output_dir, "ranked_segment_candidates_abundance.fasta")
    tsv = joinpath(output_dir, "segment_candidate_ranking.tsv")
    within_output_cap = max_variants === nothing ? (_ -> true) :
                        (candidate -> candidate.likelihood_rank <= max_variants)
    open(likelihood_fasta, "w") do io
        for candidate in filter(within_output_cap,
                sort(candidates; by = candidate -> candidate.likelihood_rank))
            println(
                io,
                ">experimental_segment_candidate|" *
                "likelihood_rank=$(candidate.likelihood_rank)|id=$(candidate.id)",
            )
            println(io, candidate.sequence)
        end
    end
    within_abundance_cap = max_variants === nothing ? (_ -> true) :
                           (candidate -> candidate.abundance_rank <= max_variants)
    open(abundance_fasta, "w") do io
        for candidate in filter(within_abundance_cap,
                sort(candidates; by = candidate -> candidate.abundance_rank))
            println(
                io,
                ">experimental_segment_candidate|" *
                "abundance_rank=$(candidate.abundance_rank)|id=$(candidate.id)",
            )
            println(io, candidate.sequence)
        end
    end
    primary = nothing
    if primary_status === :resolved
        open(primary_path, "w") do io
            println(
                io,
                ">experimental_primary_segment_candidate|" *
                "id=$(likelihood_primary.id)|status=resolved",
            )
            println(io, likelihood_primary.sequence)
        end
        primary = primary_path
    else
        rm(primary_path; force = true)
    end
    table = DataFrames.DataFrame(
        likelihood_rank = [candidate.likelihood_rank for candidate in candidates],
        abundance_rank = [candidate.abundance_rank for candidate in candidates],
        role = [String(candidate.role) for candidate in candidates],
        id = [candidate.id for candidate in candidates],
        length = [length(candidate.sequence) for candidate in candidates],
        graph_mean_log2_probability =
            [candidate.graph_mean_log2_probability for candidate in candidates],
        scored_transition_fraction =
            [candidate.scored_transition_fraction for candidate in candidates],
        strainy_coverage = [candidate.strainy_coverage for candidate in candidates],
        primary_status = fill(String(primary_status), length(candidates)),
        likelihood_primary_id = fill(likelihood_primary.id, length(candidates)),
        abundance_primary_id = fill(abundance_primary.id, length(candidates)),
    )
    CSV.write(tsv, table; delim = '\t')
    return (;
        primary,
        primary_segment_candidate_fasta = primary,
        likelihood_fasta,
        likelihood_ranked_segment_candidates_fasta = likelihood_fasta,
        abundance_fasta,
        abundance_ranked_segment_candidates_fasta = abundance_fasta,
        tsv,
        segment_candidate_ranking_tsv = tsv,
        primary_status,
        likelihood_primary_id = likelihood_primary.id,
        abundance_primary_id = abundance_primary.id,
    )
end

"""
    deconvolve_stage2(reads, config::AssemblyConfig; kwargs...)

Run Stage-1 graph accurization once, reuse metaFlye for long-read repeat-graph
layout, and rank Strainy sequence-bearing GFA `S` records as experimental
segment candidates under the accurized graph. These candidates are not complete
paths, consensus genomes, or inferred haplotypes. The production route still
uses explicitly labeled legacy best-hit abundance until a calibrated
read-candidate likelihood builder is connected to
[`infer_segment_candidate_abundances`](@ref); the soft-EM implementation in this
slice is
an inference seam, not yet the production abundance method.

`max_variants` is an optional per-ranking diagnostic output cap. `nothing`
retains every graph-supported segment candidate. It does not supply or infer a
truth strain count.
"""
function deconvolve_stage2(
        reads::R,
        config::AssemblyConfig;
        genome_size::Union{Nothing, Integer, AbstractString} = nothing,
        min_overlap::Union{Nothing, Int} = nothing,
        layout_iterations::Int = 0,
        max_variants::Union{Nothing, Int} = nothing,
        min_scored_fraction::Float64 = 0.9,
        threads::Int = Mycelia.get_default_threads()
)::Stage2SegmentCandidateResult where {R}
    config.corrector == :iterative ||
        throw(ArgumentError("Stage-2 requires corrector=:iterative"))
    config.sequencing_tech == :nanopore ||
        throw(ArgumentError("the first Stage-2 slice requires sequencing_tech=:nanopore"))
    config.output_dir === nothing &&
        throw(ArgumentError("Stage-2 requires a persistent AssemblyConfig.output_dir"))
    layout_iterations >= 0 ||
        throw(ArgumentError("layout_iterations must be nonnegative"))
    max_variants === nothing || max_variants > 0 || throw(
        ArgumentError("max_variants diagnostic output cap must be positive"))

    attempt = _reserve_stage2_attempt_directory(config.output_dir)
    output_dir = attempt.attempt_output_dir
    layout_dir = attempt.layout_dir
    strainy_dir = attempt.strainy_dir
    ranked_dir = attempt.ranked_dir

    prepared = _prepare_stage2_inputs(
        reads, config; persistent_output_dir = output_dir)
    (; corrected_fastq, graph_k, index) = prepared
    corrected_fastq == attempt.corrected_fastq || error(
        "Stage-1 corrected FASTQ escaped its reserved Stage-2 attempt")

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
    candidates = _rank_stage2_segment_candidates(
        phased.strain_contigs_gfa, phased.phased_unitig_info, index;
        max_variants = max_variants,
        min_scored_fraction = min_scored_fraction,
        support_fastq = corrected_fastq,
        threads = threads)
    ranked = _write_ranked_variants(candidates, ranked_dir, max_variants)
    ranked.primary_status === :resolved || error(
        "Stage-2 segment-candidate ranks disagree; no consumable primary " *
        "segment candidate was written")

    provenance = _stage2_attempt_provenance(attempt)
    Base.merge!(provenance, Dict{String, Any}(
        "stage1_corrector" => "iterative",
        "stage1_strategy" => String(config.strategy),
        "sequencing_tech" => "nanopore",
        "layout_tool" => "metaflye",
        "layout_read_type" => "nano-corr",
        "layout_iterations" => layout_iterations,
        "segment_candidate_tool" => "strainy",
        "segment_candidate_semantics" =>
            "experimental-gfa-s-records-not-complete-paths-or-haplotypes",
        "path_grouping_status" => "not-implemented",
        "strain_count_inference_status" => "not-implemented",
        "abundance_method" => "legacy-best-hit-minimap2",
        "probabilistic_abundance_status" => "available-requires-calibrated-likelihoods",
        "graph_k" => graph_k,
        "graph_source" => prepared.graph_source,
        "graph_rebuild_reason" => prepared.graph_rebuild_reason,
        "corrected_read_count" => prepared.corrected_read_count,
        "diagnostic_segment_candidate_cap" => max_variants,
        "min_scored_fraction" => min_scored_fraction,
    ))
    provenance["primary_status"] = String(ranked.primary_status)
    provenance["likelihood_primary_id"] = ranked.likelihood_primary_id
    provenance["abundance_primary_id"] = ranked.abundance_primary_id
    return Stage2SegmentCandidateResult(
        something(ranked.primary_segment_candidate_fasta),
        layout.assembly,
        layout.graph,
        ranked.likelihood_ranked_segment_candidates_fasta,
        ranked.abundance_ranked_segment_candidates_fasta,
        ranked.segment_candidate_ranking_tsv,
        corrected_fastq,
        candidates,
        provenance)
end
