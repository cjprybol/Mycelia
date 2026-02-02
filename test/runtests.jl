# Mycelia Test Suite - Tiered Testing System
#
# Test Tiers:
#   1. Core (default): Pure Julia tests - no network, no external tools
#   2. External: Tests requiring Bioconda tools, NCBI downloads, HPC resources
#
# Usage:
#   Core tests (CI/local):      julia --project=. -e "import Pkg; Pkg.test()"
#   Full tests (HPC):           MYCELIA_RUN_ALL=true julia --project=. -e 'import Pkg; Pkg.test()'
#   Full tests (alt):           MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'import Pkg; Pkg.test()'
#   Full tests (portable):      LD_LIBRARY_PATH="" MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'import Pkg; Pkg.update(); Pkg.instantiate(); Pkg.precompile(); Pkg.test()'
#   Benchmarks:                 julia --project=. benchmarking/run_all_benchmarks.jl
#   Tutorials:                  julia --project=. tutorials/run_all_tutorials.jl
#
# External dependencies (skipped by default):
#   - NCBI datasets CLI, SRA tools (prefetch, fasterq_dump)
#   - Bioconda tools: QUAST, BUSCO, CheckM, CoverM, mosdepth, MEGAHIT, SPAdes, etc.
#   - BLAST databases, reference genome downloads
#   - Bandage (visualization), minimap2 (alignment)
#
# Install Julia LTS: curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel lts

const MYCELIA_RUN_ALL = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
const MYCELIA_RUN_EXTERNAL = MYCELIA_RUN_ALL ||
                             lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
const PROJECT_ROOT = dirname(@__DIR__)

const TEST_ARTIFACT_DIRS = [
    joinpath(PROJECT_ROOT, "busco_downloads"),
    joinpath(PROJECT_ROOT, "busco_results"),
    joinpath(PROJECT_ROOT, "test", "busco_downloads"),
    joinpath(PROJECT_ROOT, "test", "busco_results")
]

function cleanup_test_artifacts()
    if any(isdir, TEST_ARTIFACT_DIRS)
        @info "Cleaning test artifacts (busco downloads/results)."
    end
    for path in TEST_ARTIFACT_DIRS
        if isdir(path)
            rm(path; recursive = true, force = true)
        end
    end
    return nothing
end

atexit(cleanup_test_artifacts)

# Files requiring external tools, network downloads, or HPC resources.
# These are skipped by default; enable with MYCELIA_RUN_EXTERNAL=true.
const EXTERNAL_DEPENDENCY_FILES = Set([
    # Dir 1: Data acquisition (NCBI downloads, SRA tools)
    "ncbi_download.jl",
    "simulation_fastq.jl",  # PacBio simulation, NCBI downloads
    "un_parallel_corpus.jl",
    # Dir 2: Preprocessing (SRA tools)
    "preprocessing.jl",
    # Dir 4: Assembly (external assemblers, Bandage, MEGAHIT)
    "bandage_integration.jl",
    "megahit_phix_workflow.jl",
    "assembly_merging.jl",
    "hpc_assembly_wrappers_test.jl",
    # Dir 5: Validation (QUAST, BUSCO, CheckM, CoverM, mosdepth)
    "quast_busco_wrappers_test.jl",
    "checkm_tools.jl",
    "coverm_wrappers.jl",
    "coverm_integration_extended.jl",
    "mosdepth_coverage_qc.jl",
    "coverage_taxonomy_integration.jl",
    # Dir 6: Annotation (Prodigal ORF caller)
    "orf_callers.jl",
    # Dir 7: Comparative (BLAST, minimap2, external datasets)
    "blastdb_integration.jl",
    "gold_standard_genome_comparison.jl",
    "multiomics_alignment.jl",
    "sequence_comparison.jl"
])

# Helper function to include all test files in a directory
function include_all_tests(dir)
    test_count = 0
    skipped_count = 0
    for (root, dirs, files) in walkdir(dir)
        if occursin(".ipynb_checkpoints", root)
            continue
        end
        for file in sort(files)
            if endswith(file, ".jl")
                # Skip external dependency files unless MYCELIA_RUN_EXTERNAL is set
                if !MYCELIA_RUN_EXTERNAL
                    if occursin("third_party_assemblers", file) ||
                       file in EXTERNAL_DEPENDENCY_FILES
                        skipped_count += 1
                        continue
                    end
                end
                include(joinpath(root, file))
                test_count += 1
            end
        end
    end
    if skipped_count > 0 && !MYCELIA_RUN_EXTERNAL
        @info "Skipped $skipped_count test file(s) with external dependencies in $(basename(dir)). Set MYCELIA_RUN_EXTERNAL=true to include."
    end
    return test_count
end

# Aqua.jl quality assurance tests
include("aqua.jl")

# JET.jl static analysis - uncomment to enable (can be slow)
include("jet.jl")

# ExplicitImports.jl - enforce import style (no implicit imports from using)
include("explicit_imports.jl")

include_all_tests(joinpath(@__DIR__, "1_data_acquisition"))
# For debugging individual suites, include explicit files instead of the directory sweep.
# for file in (
#     "1_data_acquisition/ncbi_download.jl",
#     "1_data_acquisition/simulation_fasta.jl",
#     "1_data_acquisition/simulation_fastq.jl",
# )
#     include(joinpath(@__DIR__, file))
# end

include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"))
# For debugging individual suites, include explicit files instead of the directory sweep.
# for file in (
#     "2_preprocessing_qc/alphabets.jl",
#     "2_preprocessing_qc/constants.jl",
#     "2_preprocessing_qc/dimensionality_reduction_and_clustering.jl",
#     "2_preprocessing_qc/diversity_sampling.jl",
#     "2_preprocessing_qc/preprocessing.jl",
#     "2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl",
#     "2_preprocessing_qc/saturation_plot_test.jl",
#     "2_preprocessing_qc/sequence-complexity.jl",
#     "2_preprocessing_qc/sequence_io.jl",
# )
#     include(joinpath(@__DIR__, file))
# end

include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"))
# For debugging individual suites, include explicit files instead of the directory sweep.
# for file in (
#     "3_feature_extraction_kmer/kmer_analysis.jl",
#     "3_feature_extraction_kmer/qualmer-analysis.jl",
# )
#     include(joinpath(@__DIR__, file))
# end

# Stage 4 (assembly): include individual files until stable, then switch back to include_all_tests.
include_all_tests(joinpath(@__DIR__, "4_assembly"))
# for file in [
#         "test/4_assembly/aa_fasta_singlestrand_test.jl",
#         "test/4_assembly/aa_fastq_singlestrand_test.jl",
#         "test/4_assembly/aa_kmer_graph_test.jl",
#         "test/4_assembly/aa_kmer_singlestrand_test.jl",
#         "test/4_assembly/aa_qualmer_graph_test.jl",
#         "test/4_assembly/aa_qualmer_singlestrand_test.jl",
#         "test/4_assembly/amino_acid_fastq_test.jl",
#         "test/4_assembly/assembly_merging.jl",
#         "test/4_assembly/bandage_integration.jl",
#         "test/4_assembly/basic_graph_tests.jl",
#         "test/4_assembly/canonicalization_consistency_test.jl",
#         "test/4_assembly/complete_assembly_workflow_test.jl",
#         "test/4_assembly/comprehensive_correctness_tests.jl",
#         "test/4_assembly/comprehensive_fixes_tests.jl",
#         "test/4_assembly/comprehensive_type_stable_corrected_tests.jl",
#         "test/4_assembly/comprehensive_type_stable_tests.jl",
#         "test/4_assembly/dna_fasta_doublestrand_test.jl",
#         "test/4_assembly/dna_fasta_singlestrand_test.jl",
#         "test/4_assembly/dna_fastq_doublestrand_test.jl",
#         "test/4_assembly/dna_fastq_singlestrand_test.jl",
#         "test/4_assembly/dna_kmer_doublestrand_test.jl",
#         "test/4_assembly/dna_kmer_singlestrand_test.jl",
#         "test/4_assembly/dna_qualmer_doublestrand_test.jl",
#         "test/4_assembly/dna_qualmer_graph_test.jl",
#         "test/4_assembly/dna_qualmer_singlestrand_test.jl",
#         "test/4_assembly/doublestrand_canonicalization_test.jl",
#         "test/4_assembly/end_to_end_assembly_tests.jl", # broken with infinite loop somewhere?
#         "test/4_assembly/end_to_end_graph_tests.jl",
#         "test/4_assembly/evidence_functions_test.jl",
#         "test/4_assembly/evidence_structures_test.jl",
#         "test/4_assembly/gfa_io_next.jl",
#         "test/4_assembly/graph_algorithms_next.jl",
#         "test/4_assembly/graph_conversion_test.jl",
#         "test/4_assembly/graph_query_test.jl",
#         "test/4_assembly/iterative_assembly_tests.jl",
#         "test/4_assembly/kmer_edge_data_test.jl",
#         "test/4_assembly/kmer_vertex_data_test.jl",
#         "test/4_assembly/megahit_phix_workflow.jl",
#         "test/4_assembly/path_finding_test.jl",
#         "test/4_assembly/probabilistic_algorithms_next.jl",
#         "test/4_assembly/quality_functions_test.jl",
#         "test/4_assembly/rhizomorph_bubbles_and_gfa_test.jl",
#         "test/4_assembly/rhizomorph_canonical_path_test.jl",
#         "test/4_assembly/rhizomorph_conversion_errors_test.jl",
#         "test/4_assembly/rhizomorph_doublestrand_files_test.jl",
#         "test/4_assembly/rhizomorph_doublestrand_traversal_test.jl",
#         "test/4_assembly/rhizomorph_kmer_mode_support_test.jl",
#         "test/4_assembly/rhizomorph_qualmer_canonical_traversal_test.jl",
#         "test/4_assembly/rhizomorph_qualmer_rc_evidence_test.jl",
#         "test/4_assembly/rna_fasta_doublestrand_test.jl",
#         "test/4_assembly/rna_fasta_singlestrand_test.jl",
#         "test/4_assembly/rna_fastq_doublestrand_test.jl",
#         "test/4_assembly/rna_fastq_singlestrand_test.jl",
#         "test/4_assembly/rna_kmer_doublestrand_test.jl",
#         "test/4_assembly/rna_kmer_graph_test.jl",
#         "test/4_assembly/rna_kmer_singlestrand_test.jl",
#         "test/4_assembly/rna_qualmer_doublestrand_test.jl",
#         "test/4_assembly/rna_qualmer_graph_test.jl",
#         "test/4_assembly/rna_qualmer_singlestrand_test.jl",
#         "test/4_assembly/sequence_graphs_next.jl",
#         "test/4_assembly/simplification_test.jl",
#         "test/4_assembly/singlestrand_canonicalization_test.jl",
#         "test/4_assembly/six_graph_hierarchy_tests.jl",
#         "test/4_assembly/strand_specific_graph_construction_test.jl",
#         "test/4_assembly/string-graph-helpers.jl",
#         "test/4_assembly/string_graphs.jl",
#         "test/4_assembly/string_ngram_singlestrand_quality_test.jl",
#         "test/4_assembly/string_ngram_singlestrand_test.jl",
#         "test/4_assembly/string_variable_singlestrand_quality_test.jl",
#         "test/4_assembly/string_variable_singlestrand_test.jl",
#         "test/4_assembly/tda_metrics_test.jl",
#         "test/4_assembly/third_party_assemblers.jl",
#         "test/4_assembly/third_party_assemblers_hybrid.jl",
#         "test/4_assembly/third_party_assemblers_legacy.jl",
#         "test/4_assembly/third_party_assemblers_long_read_isolate.jl",
#         "test/4_assembly/third_party_assemblers_long_read_metagenomic.jl",
#         "test/4_assembly/third_party_assemblers_plass_penguin.jl",
#         "test/4_assembly/third_party_assemblers_short_read_isolate.jl",
#         "test/4_assembly/third_party_assemblers_short_read_metagenomic.jl",
#         "test/4_assembly/unicode-graph-assembly.jl",
#         "test/4_assembly/variable_length_reduced_alphabet_test.jl",
#         "test/4_assembly/variable_length_singlestrand_test.jl",
#         "test/4_assembly/variable_length_strand_conversion_test.jl"
#     ]
#     include(joinpath(PROJECT_ROOT, file))
# end

# Stage 5 (validation): focused suites (all other external-heavy validation stays opt-in for now).
include_all_tests(joinpath(@__DIR__, "5_validation"))
# for file in [
#     # "5_validation/checkm_tools.jl",
#     # "5_validation/coverm_integration_extended.jl",
#     # "5_validation/coverm_wrappers.jl",
#     # "5_validation/mosdepth_coverage_qc.jl",
#     "5_validation/quast_busco_wrappers_test.jl",
#     "5_validation/validation.jl",
#     "5_validation/viterbi_polishing_and_error_correction.jl",
# ]
#     include(joinpath(@__DIR__, file))
# end

# Stage 6 (annotation)
include_all_tests(joinpath(@__DIR__, "6_annotation"))
# for file in (
#     # "6_annotation/annotation.jl",
#     # "6_annotation/codon_optimization.jl",
#     # "6_annotation/genetic_code.jl",
#     "6_annotation/genome_features.jl",
#     "6_annotation/orf_callers.jl"
# )
#     include(joinpath(@__DIR__, file))
# end

# Stage 7 (comparative/pangenomics): lightweight, synthetic-only suites by default.
include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"))
# for file in (
#     "7_comparative_pangenomics/blastdb_integration.jl",
#     "7_comparative_pangenomics/comparative_analyses.jl",
#     "7_comparative_pangenomics/distance_metrics.jl",
#     "7_comparative_pangenomics/multiomics_alignment.jl",
#     "7_comparative_pangenomics/pangenome.jl",
#     "7_comparative_pangenomics/pangenome_wrappers.jl",
#     "7_comparative_pangenomics/panproteome.jl",
#     "7_comparative_pangenomics/pantranscriptome.jl",
#     "7_comparative_pangenomics/phylogenetics.jl",
#     "7_comparative_pangenomics/sequence_classification.jl",
#     "7_comparative_pangenomics/sequence_comparison.jl",
#     "7_comparative_pangenomics/test_pcoa_vis.jl"
# )
#     include(joinpath(@__DIR__, file))
# end

if MYCELIA_RUN_EXTERNAL
    include_all_tests(joinpath(@__DIR__, "8_tool_integration"))
else
    @info "Skipping tool integration tests; set MYCELIA_RUN_EXTERNAL=true to enable."
end
# for file in (
#     "8_tool_integration/autocycler.jl",
#     "8_tool_integration/bioconda.jl",
#     "8_tool_integration/binning_tools.jl",
#     "8_tool_integration/classification_tools.jl",
#     "8_tool_integration/foldseek.jl",
#     "8_tool_integration/mash.jl",
#     "8_tool_integration/metabuli_metaphlan_strainphlan.jl",
#     "8_tool_integration/minimap_merge_map_split.jl",
#     "8_tool_integration/pantools.jl",
#     "8_tool_integration/sentencepiece.jl",
#     "8_tool_integration/utility_functions.jl",
#     "8_tool_integration/xam.jl",
# )
#     include(joinpath(@__DIR__, file))
# end

# In-development suites (opt-in; keep list up to date for bulk enabling).
# include_all_tests(joinpath(@__DIR__, "in_development"))
# for file in (
#     "in_development/cross_validation_tests.jl",
#     "in_development/ensemble_assembly.jl",
#     "in_development/hybrid_assembly.jl",
#     "in_development/intelligent_assembly_tests.jl",
#     "in_development/polishing.jl",
#     "in_development/reinforcement_learning_mcts_tests.jl",
#     "in_development/reinforcement_learning_tests.jl",
#     "in_development/strain_resolution.jl",
#     "in_development/test_reinforcement_learning_comparison.jl",
#     "in_development/test_viroid_assembly_workflow.jl",
#     "in_development/viterbi_next.jl",
# )
#     include(joinpath(@__DIR__, file))
# end

# Deprecated suites (opt-in).
# for file in (
#     "deprecated/sequence_graphs.jl",
# )
#     include(joinpath(@__DIR__, file))
# end

# Metadata/fixture helpers (manual).
# include(joinpath(@__DIR__, "metadata/download_comebin_data.jl"))

# Integration harness (manual).
# include(joinpath(@__DIR__, "run_integration_tests.jl"))
