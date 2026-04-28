import DataFrames
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_benchmark_harness.jl"))

Test.@testset "Rhizomorph benchmark manifest" begin
    manifest = load_rhizomorph_benchmark_manifest()
    Test.@test manifest["schema_version"] == 1
    Test.@test validate_rhizomorph_benchmark_manifest(manifest)

    datasets = list_rhizomorph_benchmark_datasets()
    Test.@test DataFrames.nrow(datasets) >= 10
    Test.@test issubset(
        Set(["toy_control", "public_isolate", "heterogeneous_candidate"]),
        Set(datasets.category)
    )
    Test.@test "synthetic_isolate_5386" in datasets.id
    Test.@test "phix174" in datasets.id
    Test.@test "zymo_d6300" in datasets.id

    slices = list_rhizomorph_benchmark_slices()
    Test.@test DataFrames.nrow(slices) == 7
    Test.@test Set(slices.id) == Set(["H1", "H2", "H3", "H4", "H5", "H6", "H7"])
    Test.@test all(occursin("benchmarking/rhizomorph_benchmark_harness.jl", entrypoint)
        for entrypoint in slices.entrypoint)
end

Test.@testset "Rhizomorph benchmark dry-run plans" begin
    ci_plan = build_rhizomorph_benchmark_plan(scale = "ci")
    Test.@test !isempty(ci_plan)
    Test.@test all(ci_plan.ci_suitability .== "ci")
    Test.@test "synthetic_isolate_5386" in ci_plan.dataset_id
    Test.@test !("phix174" in ci_plan.dataset_id)
    Test.@test all(.!ci_plan.implemented)

    full_plan = build_rhizomorph_benchmark_plan(scale = "full", hypothesis_ids = ["H2", "H7"])
    Test.@test Set(full_plan.hypothesis_id) == Set(["H2", "H7"])
    Test.@test "phix174" in full_plan.dataset_id
    Test.@test !("zymo_d6300" in full_plan.dataset_id)

    candidate_plan = build_rhizomorph_benchmark_plan(
        scale = "candidate",
        hypothesis_ids = ["H5"],
        dataset_ids = ["zymo_d6300", "atcc_msa_1003"]
    )
    Test.@test Set(candidate_plan.dataset_id) == Set(["zymo_d6300", "atcc_msa_1003"])
    Test.@test all(candidate_plan.hypothesis_id .== "H5")

    Test.@test_throws ErrorException build_rhizomorph_benchmark_plan(scale = "tiny")
    Test.@test_throws ErrorException build_rhizomorph_benchmark_plan(hypothesis_ids = ["H9"])
    Test.@test_throws ErrorException build_rhizomorph_benchmark_plan(dataset_ids = ["unknown_dataset"])
    Test.@test_throws ErrorException run_rhizomorph_benchmark_harness(dry_run = false)
end
