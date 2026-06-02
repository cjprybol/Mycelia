import CSV
import DataFrames
import JSON
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
    Test.@test datasets.accessions[findfirst(==("cccv_viroid"), datasets.id)] == "NC_001462.1"

    slices = list_rhizomorph_benchmark_slices()
    Test.@test DataFrames.nrow(slices) == 7
    Test.@test Set(slices.id) == Set(["H1", "H2", "H3", "H4", "H5", "H6", "H7"])
    Test.@test all(occursin("benchmarking/rhizomorph_benchmark_harness.jl", entrypoint)
        for entrypoint in slices.entrypoint)

    for dataset in manifest["datasets"]
        provenance = dataset["provenance"]
        if provenance["kind"] == "ncbi_refseq"
            Test.@test all(occursin(".", accession) for accession in provenance["accessions"])
        end
    end
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

    invalid_manifest = deepcopy(load_rhizomorph_benchmark_manifest(validate = false))
    invalid_manifest["schema"]["ci_suitability_values"] = ["ci", "full", "tiny"]
    Test.@test_throws ErrorException validate_rhizomorph_benchmark_manifest(invalid_manifest)

    Test.@test_throws ErrorException _flag_values(["--scale", "--slice"], "--scale")
    Test.@test_throws ErrorException main(["--unknown"])
end

Test.@testset "Rhizomorph benchmark public-record artifacts" begin
    mktempdir() do output_dir
        artifacts = write_rhizomorph_benchmark_plan_artifacts(
            output_dir = output_dir,
            scale = "ci",
            run_id = "ci_smoke",
            command_args = ["--plan", "--scale", "ci", "--write-artifacts"],
            generated_at = "2026-06-02T00:00:00Z"
        )

        Test.@test isdir(joinpath(output_dir, "tables"))
        Test.@test isdir(joinpath(output_dir, "plots"))
        Test.@test isdir(joinpath(output_dir, "logs"))
        Test.@test isdir(joinpath(output_dir, "provenance"))
        Test.@test isfile(joinpath(output_dir, "artifact-index.json"))
        Test.@test isfile(joinpath(output_dir, "provenance", "run.provenance.json"))

        plan_csv = joinpath(output_dir, "tables", "rhizomorph_benchmark_plan.csv")
        plan_provenance = joinpath(output_dir, "provenance", "rhizomorph_benchmark_plan.provenance.json")
        Test.@test isfile(plan_csv)
        Test.@test isfile(plan_provenance)

        plan_table = DataFrames.DataFrame(CSV.File(plan_csv))
        required_columns = [
            "benchmark_schema_version",
            "benchmark_run_id",
            "benchmark_scale",
            "benchmark_dataset_ids",
            "benchmark_git_commit",
            "benchmark_dataset_id",
            "benchmark_hypothesis_id"
        ]
        Test.@test issubset(Set(required_columns), Set(names(plan_table)))
        Test.@test all(plan_table.benchmark_run_id .== "ci_smoke")
        Test.@test all(plan_table.benchmark_scale .== "ci")
        Test.@test all(plan_table.benchmark_schema_version .== BENCHMARK_ARTIFACT_SCHEMA_VERSION)
        Test.@test plan_table.benchmark_dataset_id == plan_table.dataset_id
        Test.@test plan_table.benchmark_hypothesis_id == plan_table.hypothesis_id
        Test.@test all(plan_table.ci_suitability .== "ci")

        run_provenance = JSON.parsefile(joinpath(output_dir, "provenance", "run.provenance.json"))
        Test.@test run_provenance["run_id"] == "ci_smoke"
        Test.@test run_provenance["scale"] == "ci"
        Test.@test run_provenance["generated_at"] == "2026-06-02T00:00:00Z"
        Test.@test run_provenance["command_args"] == ["--plan", "--scale", "ci", "--write-artifacts"]
        Test.@test issubset(Set(["synthetic_isolate_5386", "synthetic_metagenome_pair"]), Set(run_provenance["dataset_ids"]))
        Test.@test haskey(run_provenance["tool_versions"], "julia")
        Test.@test haskey(run_provenance["tool_versions"], "Mycelia")

        index = JSON.parsefile(joinpath(output_dir, "artifact-index.json"))
        Test.@test index["directories"]["tables"] == "tables"
        Test.@test haskey(index["tables"], "rhizomorph_benchmark_plan")
        Test.@test index["tables"]["rhizomorph_benchmark_plan"]["table"] ==
                   joinpath("tables", "rhizomorph_benchmark_plan.csv")

        mktempdir() do second_output_dir
            second_artifacts = write_rhizomorph_benchmark_plan_artifacts(
                output_dir = second_output_dir,
                scale = "ci",
                run_id = "ci_smoke",
                command_args = ["--plan", "--scale", "ci", "--write-artifacts"],
                generated_at = "2026-06-02T00:00:00Z"
            )
            Test.@test read(artifacts.index, String) == read(second_artifacts.index, String)
            Test.@test read(artifacts.provenance, String) == read(second_artifacts.provenance, String)
            Test.@test read(plan_csv, String) ==
                       read(joinpath(second_output_dir, "tables", "rhizomorph_benchmark_plan.csv"), String)
        end
    end
end
