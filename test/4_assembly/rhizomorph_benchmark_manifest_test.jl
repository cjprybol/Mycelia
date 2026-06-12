import CSV
import DataFrames
import JSON
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_benchmark_harness.jl"))

function test_error_message(callable, expected_message::AbstractString)
    try
        callable()
        Test.@test false
    catch error_value
        Test.@test occursin(expected_message, sprint(showerror, error_value))
    end
end

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
    Test.@test Set(slices.id[slices.status .== "implemented"]) == Set(["H1", "H2", "H7"])
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
    Test.@test Set(ci_plan.hypothesis_id[ci_plan.implemented]) == Set(["H1", "H2", "H7"])

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

    test_error_message(
        () -> build_rhizomorph_benchmark_plan(scale = "tiny"),
        "Unknown Rhizomorph benchmark scale"
    )
    test_error_message(
        () -> build_rhizomorph_benchmark_plan(hypothesis_ids = ["H9"]),
        "Unknown Rhizomorph benchmark hypothesis id"
    )
    test_error_message(
        () -> build_rhizomorph_benchmark_plan(dataset_ids = ["unknown_dataset"]),
        "Unknown Rhizomorph benchmark dataset id"
    )

    invalid_manifest = deepcopy(load_rhizomorph_benchmark_manifest(validate = false))
    invalid_manifest["schema"]["ci_suitability_values"] = ["ci", "full", "tiny"]
    test_error_message(
        () -> validate_rhizomorph_benchmark_manifest(invalid_manifest),
        "Unsupported ci_suitability_values"
    )

    test_error_message(
        () -> _flag_values(["--scale", "--slice"], "--scale"),
        "Missing value after --scale"
    )
    test_error_message(
        () -> main(["--unknown"]),
        "Unknown flag"
    )
end

Test.@testset "Rhizomorph executable benchmark slices" begin
    mktempdir() do output_dir
        artifacts = run_rhizomorph_benchmark_harness(
            dry_run = false,
            output_dir = output_dir,
            scale = "ci",
            hypothesis_ids = ["H1"],
            dataset_ids = ["synthetic_isolate_5386"],
            run_id = "h1_smoke",
            command_args = ["--execute", "--slice", "H1", "--scale", "ci"],
            generated_at = "2026-06-02T00:00:00Z"
        )

        graph_csv = joinpath(output_dir, "tables", "graph_construction_metrics.csv")
        summary_csv = joinpath(output_dir, "tables", "benchmark_suite_summary.csv")
        Test.@test isfile(artifacts.index)
        Test.@test isfile(graph_csv)
        Test.@test isfile(summary_csv)

        graph_table = DataFrames.DataFrame(CSV.File(graph_csv))
        Test.@test DataFrames.nrow(graph_table) == 4
        Test.@test all(graph_table.hypothesis_id .== "H1")
        Test.@test all(graph_table.dataset_id .== "synthetic_isolate_5386")
        Test.@test all(graph_table.status .== "ok")
        Test.@test all(graph_table.runtime_seconds .>= 0)
        Test.@test all(graph_table.allocated_bytes .>= 0)
        Test.@test all(graph_table.vertices .> 0)

        summary_table = DataFrames.DataFrame(CSV.File(summary_csv))
        h1_summary = summary_table[summary_table.hypothesis_id .== "H1", :]
        Test.@test h1_summary.result_rows[1] == 4
        Test.@test h1_summary.ok_rows[1] == 4
    end

    mktempdir() do output_dir
        run_rhizomorph_benchmark_harness(
            dry_run = false,
            output_dir = output_dir,
            scale = "ci",
            hypothesis_ids = ["H2"],
            dataset_ids = ["synthetic_isolate_5386"],
            run_id = "h2_smoke",
            command_args = ["--execute", "--slice", "H2", "--scale", "ci"],
            generated_at = "2026-06-02T00:00:00Z"
        )

        assembly_csv = joinpath(output_dir, "tables", "assembly_accuracy_metrics.csv")
        summary_csv = joinpath(output_dir, "tables", "benchmark_suite_summary.csv")
        Test.@test isfile(assembly_csv)
        Test.@test isfile(summary_csv)

        assembly_table = DataFrames.DataFrame(CSV.File(assembly_csv))
        Test.@test DataFrames.nrow(assembly_table) == 4
        Test.@test all(assembly_table.hypothesis_id .== "H2")
        Test.@test all(assembly_table.dataset_id .== "synthetic_isolate_5386")
        Test.@test Set(assembly_table.k) == Set([15, 21])
        Test.@test Set(assembly_table.graph_mode) == Set(["singlestrand", "doublestrand"])
        Test.@test all(assembly_table.status .== "ok")
        Test.@test all(assembly_table.runtime_seconds .>= 0)
        Test.@test all(assembly_table.allocated_bytes .>= 0)
        Test.@test length(unique(assembly_table.contigs_path)) == 4
        Test.@test all(isfile.(assembly_table.contigs_path))

        summary_table = DataFrames.DataFrame(CSV.File(summary_csv))
        h2_summary = summary_table[summary_table.hypothesis_id .== "H2", :]
        Test.@test h2_summary.result_rows[1] == 4
        Test.@test h2_summary.ok_rows[1] == 4
    end

    mktempdir() do output_dir
        run_rhizomorph_benchmark_harness(
            dry_run = false,
            output_dir = output_dir,
            scale = "ci",
            hypothesis_ids = ["H7"],
            dataset_ids = ["synthetic_isolate_5386"],
            assemblers = ["MEGAHIT"],
            run_external = false,
            run_id = "h7_skip_smoke",
            command_args = ["--execute", "--slice", "H7", "--assembler", "MEGAHIT"],
            generated_at = "2026-06-02T00:00:00Z"
        )

        comparison_csv = joinpath(output_dir, "tables", "assembler_comparison_metrics.csv")
        comparison_table = DataFrames.DataFrame(CSV.File(comparison_csv))
        Test.@test DataFrames.nrow(comparison_table) == 1
        Test.@test comparison_table.hypothesis_id[1] == "H7"
        Test.@test comparison_table.assembler[1] == "MEGAHIT"
        Test.@test comparison_table.status[1] == "skipped"
        Test.@test comparison_table.run_external[1] == false
        Test.@test occursin("MYCELIA_RUN_EXTERNAL", comparison_table.error[1])
    end
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
