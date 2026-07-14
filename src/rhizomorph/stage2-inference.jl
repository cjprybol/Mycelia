"""
One log-likelihood contribution from an alignment of a read to a candidate.

`log_likelihood` uses natural-log units and must already marginalize every
plausible alignment/start state for this read-candidate pair under an explicit
error model. Exactly one record per pair is accepted, preventing duplicated
aligner rows from manufacturing support.
"""
struct ReadCandidateLikelihood
    read_id::String
    candidate_id::String
    log_likelihood::Float64

    function ReadCandidateLikelihood(
            read_id::String,
            candidate_id::String,
            log_likelihood::Float64
    )::ReadCandidateLikelihood
        isempty(read_id) && throw(ArgumentError("read_id must not be empty"))
        isempty(candidate_id) && throw(ArgumentError("candidate_id must not be empty"))
        isfinite(log_likelihood) ||
            throw(ArgumentError("alignment log_likelihood must be finite"))
        return new(read_id, candidate_id, log_likelihood)
    end
end

function ReadCandidateLikelihood(
        read_id::AbstractString,
        candidate_id::AbstractString,
        log_likelihood::Real
)::ReadCandidateLikelihood
    return ReadCandidateLikelihood(
        String(read_id), String(candidate_id), Float64(log_likelihood))
end

"""
Per-read log-likelihood under the explicit noise component.

The noise likelihood may represent an unassembled haplotype, contamination, or
an aligner failure. It uses natural-log units and must be supplied once for each
read considered by [`infer_candidate_abundances`](@ref).
"""
struct ReadNoiseLikelihood
    read_id::String
    log_likelihood::Float64

    function ReadNoiseLikelihood(
            read_id::String,
            log_likelihood::Float64
    )::ReadNoiseLikelihood
        isempty(read_id) && throw(ArgumentError("read_id must not be empty"))
        isfinite(log_likelihood) ||
            throw(ArgumentError("noise log_likelihood must be finite"))
        return new(read_id, log_likelihood)
    end
end

function ReadNoiseLikelihood(
        read_id::AbstractString,
        log_likelihood::Real
)::ReadNoiseLikelihood
    return ReadNoiseLikelihood(String(read_id), Float64(log_likelihood))
end

"""
Configuration for soft-EM candidate-abundance inference.
"""
struct CandidateInferenceConfig
    max_iterations::Int
    abundance_tolerance::Float64
    primary_tolerance::Float64
    retain_full_trace::Bool

    function CandidateInferenceConfig(
            max_iterations::Int,
            abundance_tolerance::Float64,
            primary_tolerance::Float64,
            retain_full_trace::Bool
    )::CandidateInferenceConfig
        max_iterations > 0 ||
            throw(ArgumentError("max_iterations must be positive"))
        isfinite(abundance_tolerance) && abundance_tolerance >= 0.0 ||
            throw(ArgumentError(
                "abundance_tolerance must be finite and nonnegative"))
        isfinite(primary_tolerance) && primary_tolerance >= 0.0 ||
            throw(ArgumentError(
                "primary_tolerance must be finite and nonnegative"))
        return new(
            max_iterations,
            abundance_tolerance,
            primary_tolerance,
            retain_full_trace,
        )
    end
end

function CandidateInferenceConfig(
        max_iterations::Int,
        abundance_tolerance::Float64,
        primary_tolerance::Float64
)::CandidateInferenceConfig
    return CandidateInferenceConfig(
        max_iterations, abundance_tolerance, primary_tolerance, false)
end

function CandidateInferenceConfig(;
        max_iterations::Int = 500,
        abundance_tolerance::Float64 = 1.0e-10,
        primary_tolerance::Float64 = 1.0e-8,
        retain_full_trace::Bool = false
)::CandidateInferenceConfig
    return CandidateInferenceConfig(
        max_iterations,
        abundance_tolerance,
        primary_tolerance,
        retain_full_trace,
    )
end

"""
One point in the soft-EM convergence trace.

When full tracing is enabled, `candidate_abundances` is ordered
lexicographically by candidate identifier. It is empty under the bounded-memory
default. Iteration zero is the deterministic uniform initialization; its
`max_abundance_change` is zero by convention.
"""
struct EMConvergencePoint
    iteration::Int
    log_likelihood::Float64
    max_abundance_change::Float64
    candidate_abundances::Vector{Pair{String, Float64}}
    noise_abundance::Float64
end

"""
Final abundance and soft read count for one candidate.
"""
struct CandidateAbundance
    candidate_id::String
    abundance::Float64
    expected_read_count::Float64
    rank::Int
    has_alignment_support::Bool
    aligned_read_count::Int
end

"""
Result of [`infer_candidate_abundances`](@ref).

Candidates are returned in deterministic rank order: decreasing abundance,
then lexicographic candidate identifier for exact ties. Candidate abundances
and `noise_abundance` sum to one.
"""
struct CandidateInferenceResult
    candidates::Vector{CandidateAbundance}
    noise_abundance::Float64
    noise_expected_read_count::Float64
    trace::Vector{EMConvergencePoint}
    trace_sha256::String
    input_sha256::String
    config::CandidateInferenceConfig
    final_log_likelihood::Float64
    iterations::Int
    converged::Bool
    n_reads::Int
    n_alignments::Int
    n_active_pairs::Int
    primary_candidate_id::Union{Nothing, String}
    primary_status::Symbol
end

struct _Stage2InferenceData
    read_ids::Vector{String}
    candidate_ids::Vector{String}
    active_likelihoods::Vector{Vector{Tuple{Int, Float64}}}
    noise_log_likelihoods::Vector{Float64}
    candidate_active_read_counts::Vector{Int}
    n_alignments::Int
    n_active_pairs::Int
end

function _stage2_update_uint64!(
        context::Mycelia.SHA.SHA2_256_CTX,
        buffer::Vector{UInt8},
        value::UInt64
)::Nothing
    length(buffer) == 8 ||
        throw(DimensionMismatch("binary digest buffer must contain eight bytes"))
    @inbounds for index in eachindex(buffer)
        shift = 64 - 8 * index
        buffer[index] = UInt8((value >> shift) & 0xff)
    end
    Mycelia.SHA.update!(context, buffer)
    return nothing
end

function _stage2_update_bytes!(
        context::Mycelia.SHA.SHA2_256_CTX,
        buffer::Vector{UInt8},
        value::String
)::Nothing
    bytes = codeunits(value)
    _stage2_update_uint64!(context, buffer, UInt64(length(bytes)))
    Mycelia.SHA.update!(context, bytes)
    return nothing
end

function _stage2_update_tag!(
        context::Mycelia.SHA.SHA2_256_CTX,
        buffer::Vector{UInt8},
        tag::UInt8
)::Nothing
    length(buffer) == 1 ||
        throw(DimensionMismatch("binary digest tag buffer must contain one byte"))
    buffer[1] = tag
    Mycelia.SHA.update!(context, buffer)
    return nothing
end

function _stage2_inference_data_sha256(data::_Stage2InferenceData)::String
    context = Mycelia.SHA.SHA2_256_CTX()
    integer_buffer = Vector{UInt8}(undef, 8)
    tag_buffer = Vector{UInt8}(undef, 1)
    _stage2_update_uint64!(
        context, integer_buffer, UInt64(length(data.candidate_ids)))
    for candidate_id in data.candidate_ids
        _stage2_update_tag!(context, tag_buffer, 0x01)
        _stage2_update_bytes!(context, integer_buffer, candidate_id)
    end
    _stage2_update_uint64!(context, integer_buffer, UInt64(length(data.read_ids)))
    for read_index in eachindex(data.read_ids)
        _stage2_update_tag!(context, tag_buffer, 0x02)
        _stage2_update_bytes!(context, integer_buffer, data.read_ids[read_index])
        likelihoods = data.active_likelihoods[read_index]
        _stage2_update_uint64!(
            context, integer_buffer, UInt64(length(likelihoods)))
        for (candidate_index, log_likelihood) in likelihoods
            _stage2_update_tag!(context, tag_buffer, 0x03)
            _stage2_update_uint64!(
                context, integer_buffer, UInt64(candidate_index))
            _stage2_update_uint64!(
                context, integer_buffer, reinterpret(UInt64, log_likelihood))
        end
        _stage2_update_tag!(context, tag_buffer, 0x04)
        _stage2_update_uint64!(
            context,
            integer_buffer,
            reinterpret(UInt64, data.noise_log_likelihoods[read_index]),
        )
    end
    return bytes2hex(Mycelia.SHA.digest!(context))
end

function _stage2_inference_input_sha256(
        alignments::AbstractVector{ReadCandidateLikelihood},
        noise_likelihoods::AbstractVector{ReadNoiseLikelihood};
        candidate_ids::Union{Nothing, AbstractVector{<:AbstractString}} = nothing
)::String
    data = _prepare_inference_data(
        alignments, noise_likelihoods, candidate_ids)
    return _stage2_inference_data_sha256(data)
end

"""
    logsumexp(values::AbstractVector{<:Real}) -> Float64

Compute ``\\log(\\sum_i \\exp(values_i))`` without overflow or avoidable
underflow. An empty vector or a vector containing only `-Inf` returns `-Inf`;
any `Inf` entry returns `Inf`. `NaN` inputs are rejected.

The shifted exponentials are summed in a deterministic sorted order so input
record order cannot perturb a downstream abundance tie.
"""
function logsumexp(values::AbstractVector{<:Real})::Float64
    isempty(values) && return -Inf
    converted = Float64[value for value in values]
    any(isnan, converted) && throw(ArgumentError("logsumexp values must not contain NaN"))
    maximum_value = maximum(converted)
    maximum_value == Inf && return Inf
    maximum_value == -Inf && return -Inf

    sort!(converted; rev = true)
    shifted_sum = 0.0
    for value in converted
        shifted_sum += exp(value - maximum_value)
    end
    return maximum_value + log(shifted_sum)
end

function _prepare_inference_data(
        alignments::AbstractVector{ReadCandidateLikelihood},
        noise_likelihoods::AbstractVector{ReadNoiseLikelihood},
        candidate_universe::Union{Nothing, AbstractVector{<:AbstractString}}
)::_Stage2InferenceData
    isempty(noise_likelihoods) &&
        throw(ArgumentError("at least one noise likelihood is required"))

    noise_by_read = Dict{String, Float64}()
    for noise in noise_likelihoods
        haskey(noise_by_read, noise.read_id) &&
            throw(ArgumentError("duplicate noise likelihood for read $(noise.read_id)"))
        noise_by_read[noise.read_id] = noise.log_likelihood
    end

    observed_candidate_ids = Set{String}()
    for alignment in alignments
        haskey(noise_by_read, alignment.read_id) ||
            throw(ArgumentError(
                "alignment read $(alignment.read_id) has no noise likelihood"))
        push!(observed_candidate_ids, alignment.candidate_id)
    end

    read_ids = sort!(collect(keys(noise_by_read)))
    if candidate_universe === nothing
        candidate_ids = sort!(collect(observed_candidate_ids))
    else
        candidate_ids = String[]
        seen_candidate_ids = Set{String}()
        for candidate_id in candidate_universe
            isempty(candidate_id) &&
                throw(ArgumentError("candidate_ids must not contain an empty identifier"))
            normalized_id = String(candidate_id)
            normalized_id in seen_candidate_ids &&
                throw(ArgumentError("duplicate candidate identifier $(normalized_id)"))
            push!(candidate_ids, normalized_id)
            push!(seen_candidate_ids, normalized_id)
        end
        missing_ids = sort!(collect(setdiff(observed_candidate_ids, seen_candidate_ids)))
        isempty(missing_ids) ||
            throw(ArgumentError(
                "alignment candidates are absent from candidate_ids: " *
                join(missing_ids, ", ")))
        sort!(candidate_ids)
    end
    read_index = Dict(read_id => index for (index, read_id) in enumerate(read_ids))
    candidate_index = Dict(
        candidate_id => index for (index, candidate_id) in enumerate(candidate_ids))
    active_likelihoods = [Tuple{Int, Float64}[] for _ in read_ids]
    active_read_counts = zeros(Int, length(candidate_ids))
    for alignment in alignments
        candidate_position = candidate_index[alignment.candidate_id]
        push!(
            active_likelihoods[read_index[alignment.read_id]],
            (candidate_position, alignment.log_likelihood),
        )
    end
    for (read_position, read_likelihoods) in enumerate(active_likelihoods)
        sort!(read_likelihoods; by = first)
        for support_index in eachindex(read_likelihoods)
            candidate_position = first(read_likelihoods[support_index])
            support_index > 1 && candidate_position ==
                first(read_likelihoods[support_index - 1]) &&
                throw(ArgumentError(
                    "duplicate read-candidate likelihood for " *
                    "$(read_ids[read_position]) / " *
                    "$(candidate_ids[candidate_position]); marginalize " *
                    "alignment states before inference"))
            active_read_counts[candidate_position] += 1
        end
    end
    noise_values = [noise_by_read[read_id] for read_id in read_ids]
    return _Stage2InferenceData(
        read_ids,
        candidate_ids,
        active_likelihoods,
        noise_values,
        active_read_counts,
        length(alignments),
        length(alignments),
    )
end

function _log_abundance(abundance::Float64)::Float64
    abundance > 0.0 && return log(abundance)
    return -Inf
end

function _expectation_counts_and_log_likelihood!(
        counts::Vector{Float64},
        data::_Stage2InferenceData,
        mixture::Vector{Float64},
        log_mixture::Vector{Float64}
)::Float64
    length(counts) == length(mixture) ||
        throw(DimensionMismatch("count and mixture vectors must have equal length"))
    length(log_mixture) == length(mixture) ||
        throw(DimensionMismatch("log-mixture and mixture vectors must have equal length"))
    fill!(counts, 0.0)
    total_log_likelihood = 0.0
    noise_index = length(mixture)
    for read_index in eachindex(data.read_ids)
        noise_term = log_mixture[noise_index] +
                     data.noise_log_likelihoods[read_index]
        maximum_term = noise_term
        for (candidate_index, log_likelihood) in
            data.active_likelihoods[read_index]
            candidate_term = log_mixture[candidate_index] +
                             log_likelihood
            maximum_term = max(maximum_term, candidate_term)
        end
        isfinite(maximum_term) ||
            error("read $(data.read_ids[read_index]) has zero probability " *
                  "under all components")
        shifted_sum = exp(noise_term - maximum_term)
        for (candidate_index, log_likelihood) in
            data.active_likelihoods[read_index]
            shifted_sum += exp(
                log_mixture[candidate_index] +
                log_likelihood - maximum_term,
            )
        end
        normalizer = maximum_term + log(shifted_sum)
        isfinite(normalizer) ||
            error(
                "read $(data.read_ids[read_index]) has zero probability " *
                "under all components")
        total_log_likelihood += normalizer
        counts[noise_index] += exp(noise_term - normalizer)
        for (candidate_index, log_likelihood) in
            data.active_likelihoods[read_index]
            counts[candidate_index] += exp(
                log_mixture[candidate_index] +
                log_likelihood - normalizer,
            )
        end
    end
    return total_log_likelihood
end

function _update_log_mixture!(
        log_mixture::Vector{Float64},
        mixture::Vector{Float64}
)::Nothing
    length(log_mixture) == length(mixture) ||
        throw(DimensionMismatch("log-mixture and mixture vectors must have equal length"))
    for index in eachindex(log_mixture, mixture)
        log_mixture[index] = _log_abundance(mixture[index])
    end
    return nothing
end

function _maximum_abundance_change(
        previous::Vector{Float64},
        current::Vector{Float64}
)::Float64
    length(previous) == length(current) ||
        throw(DimensionMismatch("abundance vectors must have the same length"))
    maximum_change = 0.0
    for index in eachindex(previous, current)
        maximum_change = max(
            maximum_change,
            abs(previous[index] - current[index]),
        )
    end
    return maximum_change
end

function _trace_point(
        iteration::Int,
        data::_Stage2InferenceData,
        mixture::Vector{Float64},
        maximum_change::Float64,
        log_likelihood::Float64,
        retain_full_trace::Bool
)::EMConvergencePoint
    candidate_abundances = retain_full_trace ? [
        data.candidate_ids[index] => mixture[index]
        for index in eachindex(data.candidate_ids)
    ] : Pair{String, Float64}[]
    return EMConvergencePoint(
        iteration,
        log_likelihood,
        maximum_change,
        candidate_abundances,
        mixture[end],
    )
end

function _trace_sha256(trace::Vector{EMConvergencePoint})::String
    context = Mycelia.SHA.SHA2_256_CTX()
    integer_buffer = Vector{UInt8}(undef, 8)
    tag_buffer = Vector{UInt8}(undef, 1)
    _stage2_update_uint64!(context, integer_buffer, UInt64(length(trace)))
    for point in trace
        _stage2_update_tag!(context, tag_buffer, 0x10)
        _stage2_update_uint64!(
            context, integer_buffer, UInt64(point.iteration))
        for value in (
                point.log_likelihood,
                point.max_abundance_change,
                point.noise_abundance,
            )
            _stage2_update_uint64!(
                context, integer_buffer, reinterpret(UInt64, value))
        end
        _stage2_update_uint64!(
            context,
            integer_buffer,
            UInt64(length(point.candidate_abundances)),
        )
        for pair in point.candidate_abundances
            _stage2_update_tag!(context, tag_buffer, 0x11)
            _stage2_update_bytes!(context, integer_buffer, pair.first)
            _stage2_update_uint64!(
                context, integer_buffer, reinterpret(UInt64, pair.second))
        end
    end
    return bytes2hex(Mycelia.SHA.digest!(context))
end

"""
    infer_candidate_abundances(
        alignments, noise_likelihoods; candidate_ids=nothing, config)

Estimate candidate and noise abundances with a soft maximum-likelihood
expectation-maximization mixture. Each read-candidate record must contain the
already-marginalized likelihood over alignment states. The E-step assigns
fractional responsibility to all supported candidates and to the explicit
per-read noise component; the M-step updates mixture abundances from those soft
counts. These are fitted mixture estimates, not posterior probabilities.

# Active-support approximation

The E-step is sparse: a read-candidate pair is active only when one marginalized
`ReadCandidateLikelihood` likelihood is supplied. An omitted pair is treated as
having likelihood zero (`log_likelihood == -Inf`). Thus callers must retain every
materially plausible pair when constructing the input.

Initialization and output ranking are deterministic. Iteration zero and every
subsequent M-step retain scalar convergence diagnostics in `result.trace`.
Set `config.retain_full_trace=true` only for small forensic runs that need the
full candidate-abundance vector at every iteration; the default trace is O(I)
rather than O(I * C).

Pass `candidate_ids` to retain a declared candidate universe, including
candidates with no supplied alignment. Such candidates remain in the ranked
output with zero abundance and zero expected read count. `has_alignment_support`
reports only whether any read-candidate likelihood was supplied; it is not a
biological inclusion decision. Bayesian candidate-count inference remains a
separate later slice.
"""
function infer_candidate_abundances(
        alignments::AbstractVector{ReadCandidateLikelihood},
        noise_likelihoods::AbstractVector{ReadNoiseLikelihood};
        candidate_ids::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
        config::CandidateInferenceConfig = CandidateInferenceConfig()
)::CandidateInferenceResult
    data = _prepare_inference_data(alignments, noise_likelihoods, candidate_ids)
    input_sha256 = _stage2_inference_data_sha256(data)
    n_candidates = length(data.candidate_ids)
    n_components = n_candidates + 1
    mixture = fill(1.0 / n_components, n_components)
    updated = similar(mixture)
    log_mixture = similar(mixture)
    expected_counts = zeros(Float64, n_components)
    mixture_counts = similar(expected_counts)
    _update_log_mixture!(log_mixture, mixture)
    current_log_likelihood = _expectation_counts_and_log_likelihood!(
        expected_counts, data, mixture, log_mixture)
    trace = EMConvergencePoint[_trace_point(
        0,
        data,
        mixture,
        0.0,
        current_log_likelihood,
        config.retain_full_trace,
    )]
    converged = false
    iterations = 0

    for iteration in 1:config.max_iterations
        copyto!(mixture_counts, expected_counts)
        total_count = sum(expected_counts)
        for index in eachindex(updated, expected_counts)
            updated[index] = expected_counts[index] / total_count
        end
        maximum_change = _maximum_abundance_change(mixture, updated)
        mixture, updated = updated, mixture
        _update_log_mixture!(log_mixture, mixture)
        current_log_likelihood = _expectation_counts_and_log_likelihood!(
            expected_counts, data, mixture, log_mixture)
        push!(trace, _trace_point(
            iteration,
            data,
            mixture,
            maximum_change,
            current_log_likelihood,
            config.retain_full_trace,
        ))
        iterations = iteration
        if maximum_change <= config.abundance_tolerance
            converged = true
            break
        end
    end

    rank_order = sortperm(
        eachindex(data.candidate_ids);
        by = index -> (-mixture[index], data.candidate_ids[index]),
    )
    rank_by_index = zeros(Int, n_candidates)
    for (rank, candidate_index) in enumerate(rank_order)
        rank_by_index[candidate_index] = rank
    end
    candidates = [
        CandidateAbundance(
            data.candidate_ids[candidate_index],
            mixture[candidate_index],
            mixture_counts[candidate_index],
            rank_by_index[candidate_index],
            data.candidate_active_read_counts[candidate_index] > 0,
            data.candidate_active_read_counts[candidate_index],
        )
        for candidate_index in rank_order
    ]
    primary_candidate_id = nothing
    primary_status = if !converged
        :nonconverged
    elseif isempty(candidates) || candidates[1].abundance <= config.primary_tolerance
        :no_candidate_support
    elseif mixture[end] + config.primary_tolerance >= candidates[1].abundance
        :noise_dominant
    elseif length(candidates) > 1 &&
           abs(candidates[1].abundance - candidates[2].abundance) <=
           config.primary_tolerance
        :tied
    else
        primary_candidate_id = candidates[1].candidate_id
        :resolved
    end
    return CandidateInferenceResult(
        candidates,
        mixture[end],
        mixture_counts[end],
        trace,
        _trace_sha256(trace),
        input_sha256,
        config,
        current_log_likelihood,
        iterations,
        converged,
        length(data.read_ids),
        data.n_alignments,
        data.n_active_pairs,
        primary_candidate_id,
        primary_status,
    )
end
