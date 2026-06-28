import CSV
import DataFrames
import JSON
import Test

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

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
    h1_ci_plan = ci_plan[ci_plan.hypothesis_id .== "H1", :]
    Test.@test DataFrames.nrow(h1_ci_plan) == 1
    Test.@test only(h1_ci_plan.dataset_id) == "rhizomorph_graph_unit_fixtures"
    Test.@test only(h1_ci_plan.implemented)
    Test.@test all(ci_plan[ci_plan.hypothesis_id .!= "H1", :].implemented .== false)

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

    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        build_rhizomorph_benchmark_plan(scale = "tiny")
    end
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        build_rhizomorph_benchmark_plan(hypothesis_ids = ["H9"])
    end
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        build_rhizomorph_benchmark_plan(dataset_ids = ["unknown_dataset"])
    end
    test_throws_message(
        ErrorException,
        "pass --slice H1"
    ) do
        run_rhizomorph_benchmark_harness(dry_run = false)
    end

    h1_smoke = run_rhizomorph_benchmark_harness(dry_run = false, hypothesis_ids = ["H1"])
    Test.@test DataFrames.nrow(h1_smoke) == 10
    Test.@test Set(h1_smoke.dataset_id) == Set(["rhizomorph_graph_unit_fixtures"])
    Test.@test Set(h1_smoke.hypothesis_id) == Set(["H1"])
    Test.@test Set(h1_smoke.fixture_id) == Set(["H1-G0", "H1-G1", "H1-G2", "H1-G3", "H1-G4"])
    Test.@test all(
        fixture_id -> Set(h1_smoke.strategy_name[h1_smoke.fixture_id .== fixture_id]) ==
                      Set(["ExhaustiveViterbiObjectiveOracle", "GreedyViterbi"]),
        unique(h1_smoke.fixture_id)
    )
    h1_g1_oracle = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G1") .&
        (h1_smoke.strategy_name .== "ExhaustiveViterbiObjectiveOracle"),
        :
    ])
    h1_g1_greedy = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G1") .& (h1_smoke.strategy_name .== "GreedyViterbi"),
        :
    ])
    Test.@test h1_g1_oracle.path_vertices == "S,B1,B2,T"
    Test.@test h1_g1_greedy.path_vertices == "S,A1,A2,T"
    Test.@test h1_g1_greedy.expected_strategy_path_match
    Test.@test !h1_g1_greedy.exact_path_match
    Test.@test h1_g1_oracle.log_likelihood_gap_dp_minus_greedy > 1.0

    h1_g2_oracle = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G2") .&
        (h1_smoke.strategy_name .== "ExhaustiveViterbiObjectiveOracle"),
        :
    ])
    h1_g2_greedy = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G2") .& (h1_smoke.strategy_name .== "GreedyViterbi"),
        :
    ])
    Test.@test h1_g2_oracle.path_vertices == "S,B1,B2,B3,T"
    Test.@test h1_g2_greedy.path_vertices == "S,A1,A2,A3,T"
    Test.@test !h1_g2_greedy.exact_path_match

    h1_g3_greedy = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G3") .& (h1_smoke.strategy_name .== "GreedyViterbi"),
        :
    ])
    Test.@test h1_g3_greedy.path_vertices == "S,L1,R1,T"
    Test.@test h1_g3_greedy.failure_code == "length_mismatch"
    Test.@test h1_g3_greedy.repeat_copy_number_error == 1

    h1_g4_oracle = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G4") .&
        (h1_smoke.strategy_name .== "ExhaustiveViterbiObjectiveOracle"),
        :
    ])
    h1_g4_greedy = only(h1_smoke[
        (h1_smoke.fixture_id .== "H1-G4") .& (h1_smoke.strategy_name .== "GreedyViterbi"),
        :
    ])
    Test.@test h1_g4_oracle.exact_path_match
    Test.@test h1_g4_greedy.exact_path_match
    Test.@test h1_g4_oracle.ambiguity_margin < 0.002
    test_throws_message(
        ErrorException,
        "only supports dataset rhizomorph_graph_unit_fixtures"
    ) do
        run_rhizomorph_benchmark_harness(
            dry_run = false,
            hypothesis_ids = ["H1"],
            dataset_ids = ["rhizomorph_graph_unit_fixtures", "phix174"]
        )
    end

    invalid_manifest = deepcopy(load_rhizomorph_benchmark_manifest(validate = false))
    invalid_manifest["schema"]["ci_suitability_values"] = ["ci", "full", "tiny"]
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        validate_rhizomorph_benchmark_manifest(invalid_manifest)
    end

    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        _flag_values(["--scale", "--slice"], "--scale")
    end
    test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
        main(["--unknown"])
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

Test.@testset "Rhizomorph H1 Viterbi smoke artifacts" begin
    mktempdir() do output_dir
        artifacts = write_h1_viterbi_dp_greedy_artifacts(
            output_dir = output_dir,
            run_id = "h1_smoke_test",
            command_args = ["--slice", "H1", "--execute", "--write-artifacts"],
            generated_at = "2026-06-27T00:00:00Z"
        )
        index_json = joinpath(output_dir, "artifact-index.json")
        metrics_csv = joinpath(output_dir, "tables", "h1_viterbi_dp_greedy_path_metrics.csv")
        metrics_provenance_json = joinpath(
            output_dir,
            "provenance",
            "h1_viterbi_dp_greedy_path_metrics.provenance.json"
        )
        Test.@test isfile(index_json)
        Test.@test isfile(metrics_csv)
        Test.@test isfile(metrics_provenance_json)
        metrics = DataFrames.DataFrame(CSV.File(metrics_csv))
        Test.@test DataFrames.nrow(metrics) == 10
        Test.@test all(metrics.benchmark_hypothesis_id .== "H1")
        Test.@test all(metrics.benchmark_dataset_id .== "rhizomorph_graph_unit_fixtures")
        Test.@test Set(metrics.fixture_id) == Set(["H1-G0", "H1-G1", "H1-G2", "H1-G3", "H1-G4"])
        Test.@test Set(metrics.hypothesis_id) == Set(["H1"])
        Test.@test Set(metrics.dataset_id) == Set(["rhizomorph_graph_unit_fixtures"])
        Test.@test all(
            fixture_id -> Set(metrics.strategy_name[metrics.fixture_id .== fixture_id]) ==
                          Set(["ExhaustiveViterbiObjectiveOracle", "GreedyViterbi"]),
            unique(metrics.fixture_id)
        )

        h1_g1_oracle = only(metrics[
            (metrics.fixture_id .== "H1-G1") .&
            (metrics.strategy_name .== "ExhaustiveViterbiObjectiveOracle"),
            :
        ])
        h1_g1_greedy = only(metrics[
            (metrics.fixture_id .== "H1-G1") .&
            (metrics.strategy_name .== "GreedyViterbi"),
            :
        ])
        Test.@test h1_g1_oracle.path_vertices == "S,B1,B2,T"
        Test.@test h1_g1_greedy.path_vertices == "S,A1,A2,T"
        Test.@test h1_g1_greedy.expected_strategy_path_match
        Test.@test !h1_g1_greedy.exact_path_match
        Test.@test h1_g1_oracle.log_likelihood_gap_dp_minus_greedy > 1.0

        h1_g3_greedy = only(metrics[
            (metrics.fixture_id .== "H1-G3") .& (metrics.strategy_name .== "GreedyViterbi"),
            :
        ])
        Test.@test h1_g3_greedy.failure_code == "length_mismatch"
        Test.@test h1_g3_greedy.repeat_copy_number_error == 1
        Test.@test all(.!ismissing.(metrics.peak_rss_mib))

        index = JSON.parsefile(index_json)
        Test.@test index["tables"]["h1_viterbi_dp_greedy_path_metrics"]["table"] ==
                   joinpath("tables", "h1_viterbi_dp_greedy_path_metrics.csv")
        Test.@test index["tables"]["h1_viterbi_dp_greedy_path_metrics"]["provenance"] ==
                   joinpath("provenance", "h1_viterbi_dp_greedy_path_metrics.provenance.json")

        run_provenance = JSON.parsefile(artifacts.provenance)
        Test.@test run_provenance["dataset_ids"] == ["rhizomorph_graph_unit_fixtures"]
        Test.@test run_provenance["metadata"]["artifact_kind"] == "h1_viterbi_dp_greedy_path_metrics"

        metrics_provenance = JSON.parsefile(metrics_provenance_json)
        Test.@test metrics_provenance["metadata"]["artifact_kind"] == "h1_viterbi_dp_greedy_path_metrics"
        Test.@test metrics_provenance["dataset_ids"] == ["rhizomorph_graph_unit_fixtures"]

        test_throws_message(ErrorException, COMMON_ERROR_MESSAGE_FRAGMENTS) do
            run_rhizomorph_benchmark_harness(
                dry_run = false,
                hypothesis_ids = ["H1"],
                scale = "tiny"
            )
        end
    end
end
