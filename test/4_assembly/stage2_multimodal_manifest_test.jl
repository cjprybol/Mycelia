import Test
import TOML

const MANIFEST_PATH = joinpath(
    @__DIR__, "..", "..", "benchmarking", "rhizomorph_stage2_multimodal_benchmark.toml")

Test.@testset "Stage-2 frozen multimodal manifest" begin
    manifest = TOML.parsefile(MANIFEST_PATH)

    Test.@test manifest["schema_version"] == 1
    Test.@test manifest["reference_seed"] == 43003
    Test.@test manifest["development_seeds"] == [43003, 43006, 43007, 43008, 43009]
    Test.@test manifest["confirmation_seeds"] == [44001, 44002, 44003]
    Test.@test isempty(intersect(
        manifest["development_seeds"], manifest["confirmation_seeds"]))
    Test.@test manifest["confirmation_status"] == "reserved-unrun"
    Test.@test !manifest["candidate_count_is_known_to_inference"]
    Test.@test manifest["commit_artifacts_only_after_pass"]

    modalities = manifest["modalities"]
    Test.@test Set(keys(modalities)) == Set([
        "long_noisy", "long_accurate", "short_paired", "short_single",
        "hybrid_ont_short", "hybrid_hifi_ont"])
    Test.@test all(modality["required"] for modality in values(modalities))
    Test.@test all(!isempty(modality["comparators"]) for modality in values(modalities))
    Test.@test all(!isempty(modality["input_components"])
        for modality in values(modalities))
    Test.@test all(length(unique(modality["input_components"])) ==
                   length(modality["input_components"])
        for modality in values(modalities))
    Test.@test "myloasm" in modalities["long_accurate"]["comparators"]
    Test.@test "ultima" in modalities["short_single"]["branded_validation"]

    comparators = manifest["comparators"]
    required_comparators = union(
        (Set(modality["comparators"]) for modality in values(modalities))...)
    Test.@test required_comparators == Set(keys(comparators))
    for (comparator_id, comparator) in comparators
        Test.@test comparator["version"] == "UNRESOLVED"
        Test.@test !isempty(comparator["command"])
        Test.@test !isempty(comparator["configuration"])
        Test.@test !isempty(comparator["input_contract"])
        Test.@test comparator["configuration_budget"] == 1
        Test.@test all(modality in keys(modalities)
            for modality in comparator["modalities"])
        Test.@test comparator_id in required_comparators
        consumed = comparator["consumed_input_components"]
        Test.@test Set(keys(consumed)) == Set(comparator["modalities"])
        for modality in comparator["modalities"]
            Test.@test !isempty(consumed[modality])
            Test.@test issubset(
                Set(consumed[modality]), Set(modalities[modality]["input_components"]))
        end
    end
    Test.@test manifest["comparator_contract"]["unresolved_version_blocks_run"]
    Test.@test manifest["comparator_contract"]["input_policy"] ==
               "same-frozen-modality-bundle-declared-consumed-subset"
    Test.@test manifest["comparator_contract"]["subset_consumer_policy"] ==
               "stratify-report-non-gating"
    evaluation = manifest["evaluation"]
    Test.@test evaluation["tool"] == "QUAST"
    Test.@test evaluation["version"] == "UNRESOLVED"
    Test.@test !isempty(evaluation["command"])
    Test.@test evaluation["report_schema"] ==
               "rhizomorph-stage2-metrics-v1"
    Test.@test evaluation["unresolved_version_blocks_run"]
    for (modality_id, modality) in modalities
        Test.@test any(
            Set(comparators[comparator_id]["consumed_input_components"][modality_id]) ==
            Set(modality["input_components"])
            for comparator_id in modality["comparators"]
        )
    end
    Test.@test Set(comparators["metaspades"]["consumed_input_components"][
        "hybrid_ont_short"]) == Set(["ont", "short_r1", "short_r2"])
    Test.@test Set(comparators["metaflye-strainy"][
        "consumed_input_components"]["hybrid_ont_short"]) == Set(["ont"])
    Test.@test Set(comparators["verkko"]["consumed_input_components"][
        "hybrid_hifi_ont"]) == Set(["hifi", "ont"])
    Test.@test Set(comparators["hifiasm-meta"][
        "consumed_input_components"]["hybrid_hifi_ont"]) == Set(["hifi"])

    development_gate = manifest["gates"]["development"]
    Test.@test development_gate["evaluation_unit"] ==
               "each-required-modality-each-development-seed"
    Test.@test development_gate["comparator_id"] == "rhizomorph-native"
    Test.@test development_gate["cell_id_format"] ==
               "<modality>:<seed>:rhizomorph-native"
    Test.@test development_gate["require_exact_grid"]

    release_gate = manifest["gates"]["release"]
    Test.@test release_gate["maximum_extra_misassemblies"] == 0
    Test.@test release_gate["minimum_nga50_ratio"] == 0.99
    Test.@test release_gate["evaluation_unit"] ==
               "each-modality-each-confirmation-seed"
    Test.@test release_gate["baseline"] ==
               "best-compatible-locked-comparator-within-modality-and-input-stratum"
    Test.@test release_gate["required_input_stratum"] == "full"
    Test.@test !release_gate["subset_input_strata_are_gating"]
    Test.@test release_gate["minimum_accuracy_win_seeds_per_modality"] == 2
    Test.@test release_gate["require_three_distinct_truths"]
    Test.@test release_gate["require_exactly_three_supported_haplotypes"]
    Test.@test release_gate["require_no_spurious_or_duplicate_haplotypes"]
    Test.@test release_gate["required_strain_count_mode"] == 3
    Test.@test release_gate["require_90pct_strain_count_set_contains_truth"]
    Test.@test release_gate["require_primary_rank_agreement"]
    Test.@test release_gate["require_every_modality"]
    Test.@test release_gate["abundance_baseline_relation"] ==
               "mae-no-worse-than-best-compatible-comparator"
    Test.@test release_gate["ensemble_policy"] ==
               "noninferior-to-native-or-abstain"

    stretch_gate = manifest["gates"]["stretch"]
    Test.@test stretch_gate["minimum_nga50_ratio"] == 1.0
    Test.@test stretch_gate["maximum_genome_fraction_deficit_points"] == 0.0
    Test.@test stretch_gate["require_lower_error_every_seed"]
end
