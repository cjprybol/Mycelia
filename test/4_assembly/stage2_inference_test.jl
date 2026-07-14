import Mycelia
import Test

if !isdefined(Mycelia.Rhizomorph, :infer_candidate_abundances)
    Base.include(
        Mycelia.Rhizomorph,
        normpath(joinpath(@__DIR__, "../../src/rhizomorph/stage2-inference.jl")),
    )
end

const Stage2Inference = Mycelia.Rhizomorph

function abundance_by_id(
        result::Stage2Inference.CandidateInferenceResult
)::Dict{String, Float64}
    return Dict(
        candidate.candidate_id => candidate.abundance
        for candidate in result.candidates
    )
end

function _test_stage2_inference_error(
        f::Function,
        exception_type::Type{<:Exception},
        expected_message::AbstractString,
)::Nothing
    caught = try
        f()
        nothing
    catch error
        error
    end
    message = caught === nothing ? "" : sprint(showerror, caught)
    Test.@test caught isa exception_type
    Test.@test occursin(expected_message, message)
    return nothing
end

Test.@testset "Stage-2 stable logsumexp" begin
    expected = -1000.0 + log1p(exp(-1.0))
    Test.@test Stage2Inference.logsumexp([-1000.0, -1001.0]) ≈ expected
    Test.@test Stage2Inference.logsumexp([-Inf, -Inf]) == -Inf
    Test.@test Stage2Inference.logsumexp(Float64[]) == -Inf
    Test.@test Stage2Inference.logsumexp([Inf, -1.0]) == Inf
    _test_stage2_inference_error(
        ArgumentError,
        "logsumexp values must not contain NaN",
    ) do
        Stage2Inference.logsumexp([0.0, NaN])
    end
end

Test.@testset "Stage-2 soft-EM abundance and noise inference" begin
    alignments = Stage2Inference.ReadCandidateLikelihood[]
    noise = Stage2Inference.ReadNoiseLikelihood[]

    for read_index in 1:7
        read_id = "major-$(read_index)"
        push!(alignments, Stage2Inference.ReadCandidateLikelihood(read_id, "major", 0.0))
        push!(alignments, Stage2Inference.ReadCandidateLikelihood(read_id, "minor", -30.0))
        push!(noise, Stage2Inference.ReadNoiseLikelihood(read_id, -30.0))
    end
    for read_index in 1:2
        read_id = "minor-$(read_index)"
        push!(alignments, Stage2Inference.ReadCandidateLikelihood(read_id, "major", -30.0))
        push!(alignments, Stage2Inference.ReadCandidateLikelihood(read_id, "minor", 0.0))
        push!(noise, Stage2Inference.ReadNoiseLikelihood(read_id, -30.0))
    end
    push!(alignments, Stage2Inference.ReadCandidateLikelihood("noise-1", "major", -30.0))
    push!(alignments, Stage2Inference.ReadCandidateLikelihood("noise-1", "minor", -30.0))
    push!(noise, Stage2Inference.ReadNoiseLikelihood("noise-1", 0.0))

    config = Stage2Inference.CandidateInferenceConfig(
        max_iterations = 500,
        abundance_tolerance = 1.0e-12,
    )
    result = Stage2Inference.infer_candidate_abundances(
        reverse(alignments),
        reverse(noise);
        config,
    )
    abundances = abundance_by_id(result)

    Test.@test result isa Stage2Inference.CandidateInferenceResult
    Test.@test result.converged
    Test.@test result.n_reads == 10
    Test.@test result.n_alignments == 20
    Test.@test result.n_active_pairs == 20
    Test.@test abundances["major"] ≈ 0.7 atol = 1.0e-10
    Test.@test abundances["minor"] ≈ 0.2 atol = 1.0e-10
    Test.@test result.noise_abundance ≈ 0.1 atol = 1.0e-10
    Test.@test result.noise_expected_read_count ≈ 1.0 atol = 1.0e-9
    Test.@test result.config == config
    Test.@test result.iterations == length(result.trace) - 1
    Test.@test isfinite(result.final_log_likelihood)
    Test.@test length(result.trace_sha256) == 64
    Test.@test all(isxdigit, result.trace_sha256)
    Test.@test length(result.input_sha256) == 64
    Test.@test all(isxdigit, result.input_sha256)
    Test.@test all(isempty(point.candidate_abundances) for point in result.trace)
    Test.@test result.primary_candidate_id == "major"
    Test.@test result.primary_status == :resolved
    Test.@test sum(values(abundances)) + result.noise_abundance ≈ 1.0
    Test.@test [candidate.candidate_id for candidate in result.candidates] ==
               ["major", "minor"]
    Test.@test [candidate.rank for candidate in result.candidates] == [1, 2]
    Test.@test result.candidates[1].expected_read_count ≈ 7.0 atol = 1.0e-9
    Test.@test result.candidates[2].expected_read_count ≈ 2.0 atol = 1.0e-9

    evidence = [point.log_likelihood for point in result.trace]
    Test.@test all(
        evidence[index] + 1.0e-10 >= evidence[index - 1]
        for index in 2:length(evidence)
    )
    Test.@test first(result.trace).iteration == 0
    Test.@test last(result.trace).iteration == length(result.trace) - 1
    Test.@test last(result.trace).max_abundance_change <= config.abundance_tolerance

    forward_result = Stage2Inference.infer_candidate_abundances(
        alignments, noise; config)
    Test.@test forward_result.input_sha256 == result.input_sha256
end

Test.@testset "Stage-2 binary inference provenance" begin
    collision_left = Stage2Inference.EMConvergencePoint[
        Stage2Inference.EMConvergencePoint(
            0,
            0.0,
            0.0,
            ["a=0.0\tb" => 0.0],
            1.0,
        ),
    ]
    collision_right = Stage2Inference.EMConvergencePoint[
        Stage2Inference.EMConvergencePoint(
            0,
            0.0,
            0.0,
            ["a" => 0.0, "b" => 0.0],
            1.0,
        ),
    ]
    Test.@test Stage2Inference._trace_sha256(collision_left) !=
               Stage2Inference._trace_sha256(collision_right)

    alignments = [
        Stage2Inference.ReadCandidateLikelihood("read-2", "beta", -2.0),
        Stage2Inference.ReadCandidateLikelihood("read-1", "alpha", -1.0),
    ]
    noise = [
        Stage2Inference.ReadNoiseLikelihood("read-1", -3.0),
        Stage2Inference.ReadNoiseLikelihood("read-2", -4.0),
    ]
    first_digest = Stage2Inference._stage2_inference_input_sha256(
        alignments, noise)
    second_digest = Stage2Inference._stage2_inference_input_sha256(
        reverse(alignments), reverse(noise))
    Test.@test first_digest == second_digest
    Test.@test first_digest != Stage2Inference._stage2_inference_input_sha256(
        alignments,
        noise;
        candidate_ids = ["alpha", "beta", "unobserved"],
    )

    bulk_alignments = [
        Stage2Inference.ReadCandidateLikelihood(
            "bulk-$(index)", "candidate-$(mod1(index, 3))", -1.0)
        for index in 1:10_000
    ]
    bulk_noise = [
        Stage2Inference.ReadNoiseLikelihood("bulk-$(index)", -10.0)
        for index in 1:10_000
    ]
    bulk_data = Stage2Inference._prepare_inference_data(
        bulk_alignments, bulk_noise, nothing)
    Stage2Inference._stage2_inference_data_sha256(bulk_data)
    digest_allocations = @allocated Stage2Inference._stage2_inference_data_sha256(
        bulk_data)
    Test.@test digest_allocations <= 4_096
end

Test.@testset "Stage-2 consumes marginalized pair likelihoods" begin
    alignments = [
        Stage2Inference.ReadCandidateLikelihood("read-1", "candidate-b", -0.5),
        Stage2Inference.ReadCandidateLikelihood(
            "read-1", "candidate-a", Stage2Inference.logsumexp([-1.0, -1.0])),
    ]
    noise = [Stage2Inference.ReadNoiseLikelihood("read-1", -30.0)]
    result = Stage2Inference.infer_candidate_abundances(
        alignments,
        noise;
        config = Stage2Inference.CandidateInferenceConfig(
            max_iterations = 1000,
            abundance_tolerance = 1.0e-12,
        ),
    )

    Test.@test result.n_alignments == 2
    Test.@test result.n_active_pairs == 2
    Test.@test first(result.candidates).candidate_id == "candidate-a"
    Test.@test first(result.candidates).abundance > 0.999
end

Test.@testset "Stage-2 retained universe, active support, and deterministic ranks" begin
    alignments = [
        Stage2Inference.ReadCandidateLikelihood("read-2", "beta", 0.0),
        Stage2Inference.ReadCandidateLikelihood("read-1", "alpha", 0.0),
    ]
    noise = [
        Stage2Inference.ReadNoiseLikelihood("read-2", -100.0),
        Stage2Inference.ReadNoiseLikelihood("read-1", -100.0),
    ]
    result = Stage2Inference.infer_candidate_abundances(
        alignments,
        noise;
        candidate_ids = ["unobserved", "beta", "alpha"],
    )

    Test.@test result.converged
    Test.@test result.n_active_pairs == 2
    Test.@test [candidate.candidate_id for candidate in result.candidates] ==
               ["alpha", "beta", "unobserved"]
    Test.@test all(candidate -> candidate.abundance ≈ 0.5, result.candidates[1:2])
    Test.@test all(candidate -> candidate.expected_read_count ≈ 1.0,
        result.candidates[1:2])
    Test.@test all(candidate -> candidate.has_alignment_support,
        result.candidates[1:2])
    Test.@test all(candidate -> candidate.aligned_read_count == 1,
        result.candidates[1:2])
    Test.@test result.candidates[3].abundance == 0.0
    Test.@test result.candidates[3].expected_read_count == 0.0
    Test.@test !result.candidates[3].has_alignment_support
    Test.@test result.candidates[3].aligned_read_count == 0
    Test.@test result.noise_abundance < 1.0e-30

    tied_alignments = [
        Stage2Inference.ReadCandidateLikelihood("read-1", "beta", 0.0),
        Stage2Inference.ReadCandidateLikelihood("read-1", "alpha", 0.0),
    ]
    tied_result = Stage2Inference.infer_candidate_abundances(
        tied_alignments,
        [Stage2Inference.ReadNoiseLikelihood("read-1", -100.0)],
    )
    Test.@test [candidate.candidate_id for candidate in tied_result.candidates] ==
               ["alpha", "beta"]
    Test.@test tied_result.candidates[1].abundance ==
               tied_result.candidates[2].abundance
    Test.@test tied_result.primary_candidate_id === nothing
    Test.@test tied_result.primary_status == :tied
end

Test.@testset "Stage-2 inference validation" begin
    _test_stage2_inference_error(
        ArgumentError,
        "read_id must not be empty",
    ) do
        Stage2Inference.ReadCandidateLikelihood("", "candidate", 0.0)
    end
    _test_stage2_inference_error(
        ArgumentError,
        "alignment log_likelihood must be finite",
    ) do
        Stage2Inference.ReadCandidateLikelihood("read", "candidate", Inf)
    end
    _test_stage2_inference_error(
        ArgumentError,
        "noise log_likelihood must be finite",
    ) do
        Stage2Inference.ReadNoiseLikelihood("read", NaN)
    end
    _test_stage2_inference_error(
        ArgumentError,
        "max_iterations must be positive",
    ) do
        Stage2Inference.CandidateInferenceConfig(max_iterations = 0)
    end
    _test_stage2_inference_error(
        ArgumentError,
        "abundance_tolerance must be finite and nonnegative",
    ) do
        Stage2Inference.CandidateInferenceConfig(abundance_tolerance = -1.0)
    end
    _test_stage2_inference_error(
        ArgumentError,
        "primary_tolerance must be finite and nonnegative",
    ) do
        Stage2Inference.CandidateInferenceConfig(primary_tolerance = -1.0)
    end
    _test_stage2_inference_error(
        ArgumentError,
        "primary_tolerance must be finite and nonnegative",
    ) do
        Stage2Inference.CandidateInferenceConfig(500, 1.0e-10, NaN)
    end

    _test_stage2_inference_error(
        ArgumentError,
        "at least one noise likelihood is required",
    ) do
        Stage2Inference.infer_candidate_abundances(
            Stage2Inference.ReadCandidateLikelihood[],
            Stage2Inference.ReadNoiseLikelihood[],
        )
    end
    _test_stage2_inference_error(
        ArgumentError,
        "alignment read read-1 has no noise likelihood",
    ) do
        Stage2Inference.infer_candidate_abundances(
            [Stage2Inference.ReadCandidateLikelihood(
                "read-1", "candidate", 0.0)],
            [Stage2Inference.ReadNoiseLikelihood("read-2", 0.0)],
        )
    end
    _test_stage2_inference_error(
        ArgumentError,
        "duplicate noise likelihood for read read-1",
    ) do
        Stage2Inference.infer_candidate_abundances(
            Stage2Inference.ReadCandidateLikelihood[],
            [
                Stage2Inference.ReadNoiseLikelihood("read-1", 0.0),
                Stage2Inference.ReadNoiseLikelihood("read-1", -1.0),
            ],
        )
    end
    _test_stage2_inference_error(
        ArgumentError,
        "alignment candidates are absent from candidate_ids: candidate",
    ) do
        Stage2Inference.infer_candidate_abundances(
            [Stage2Inference.ReadCandidateLikelihood(
                "read-1", "candidate", 0.0)],
            [Stage2Inference.ReadNoiseLikelihood("read-1", 0.0)];
            candidate_ids = ["different-candidate"],
        )
    end
    _test_stage2_inference_error(
        ArgumentError,
        "duplicate read-candidate likelihood for read-1 / candidate",
    ) do
        Stage2Inference.infer_candidate_abundances(
            [
                Stage2Inference.ReadCandidateLikelihood(
                    "read-1", "candidate", 0.0),
                Stage2Inference.ReadCandidateLikelihood(
                    "read-1", "candidate", -1.0),
            ],
            [Stage2Inference.ReadNoiseLikelihood("read-1", 0.0)],
        )
    end
    _test_stage2_inference_error(
        ArgumentError,
        "duplicate candidate identifier duplicate",
    ) do
        Stage2Inference.infer_candidate_abundances(
            Stage2Inference.ReadCandidateLikelihood[],
            [Stage2Inference.ReadNoiseLikelihood("read-1", 0.0)];
            candidate_ids = ["duplicate", "duplicate"],
        )
    end

    noise_only = Stage2Inference.infer_candidate_abundances(
        Stage2Inference.ReadCandidateLikelihood[],
        [Stage2Inference.ReadNoiseLikelihood("read-1", 0.0)];
        candidate_ids = ["retained"],
    )
    Test.@test length(noise_only.candidates) == 1
    Test.@test noise_only.candidates[1].candidate_id == "retained"
    Test.@test !noise_only.candidates[1].has_alignment_support
    Test.@test noise_only.candidates[1].expected_read_count == 0.0
    Test.@test noise_only.noise_abundance == 1.0
    Test.@test noise_only.converged
    Test.@test noise_only.primary_candidate_id === nothing
    Test.@test noise_only.primary_status == :no_candidate_support

    nonconverged = Stage2Inference.infer_candidate_abundances(
        [Stage2Inference.ReadCandidateLikelihood("read-1", "candidate", 0.0)],
        [Stage2Inference.ReadNoiseLikelihood("read-1", -2.0)];
        config = Stage2Inference.CandidateInferenceConfig(max_iterations = 1),
    )
    Test.@test !nonconverged.converged
    Test.@test nonconverged.primary_candidate_id === nothing
    Test.@test nonconverged.primary_status == :nonconverged
    Test.@test nonconverged.candidates[1].expected_read_count ≈
               nonconverged.candidates[1].abundance * nonconverged.n_reads

    exact_tie = Stage2Inference.infer_candidate_abundances(
        [
            Stage2Inference.ReadCandidateLikelihood("tie", "a", 0.0),
            Stage2Inference.ReadCandidateLikelihood("tie", "b", 0.0),
        ],
        [Stage2Inference.ReadNoiseLikelihood("tie", -30.0)];
        config = Stage2Inference.CandidateInferenceConfig(
            500, 1.0e-10, 1.0e-8),
    )
    Test.@test exact_tie.primary_status == :tied
    Test.@test exact_tie.primary_candidate_id === nothing
    full_trace = Stage2Inference.infer_candidate_abundances(
        [Stage2Inference.ReadCandidateLikelihood("read", "candidate", 0.0)],
        [Stage2Inference.ReadNoiseLikelihood("read", -2.0)];
        config = Stage2Inference.CandidateInferenceConfig(
            max_iterations = 2, retain_full_trace = true),
    )
    Test.@test all(!isempty(point.candidate_abundances)
        for point in full_trace.trace)
end
