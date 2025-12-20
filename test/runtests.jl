# Mycelia Test Suite - Tiered Testing System
#
# curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel lts
# Usage:
#   Core tests:                 julia --project=. -e "import Pkg; Pkg.test()"
#   Extended tests:             MYCELIA_RUN_ALL=true julia --project=. -e 'import Pkg; Pkg.test()'
#   Extended tests:             MYCELIA_RUN_EXTERNAL=true MYCELIA_RUN_BANDAGE_DOWNLOAD=true MYCELIA_RUN_SENTENCEPIECE_INTEGRATION=true julia --project=. -e 'import Pkg; Pkg.test()'
#   Extended tests - shorthand: MYCELIA_RUN_ALL=true julia --project=. -e 'import Pkg; Pkg.test()'
#   Extended tests - portable:  LD_LIBRARY_PATH="" MYCELIA_RUN_ALL=true julia --project=. -e 'import Pkg; Pkg.update(); Pkg.instantiate(); Pkg.precompile(); Pkg.test()'
#   Benchmarks:                 julia --project=. benchmarking/run_all_benchmarks.jl
#   Tutorials:                  julia --project=. tutorials/run_all_tutorials.jl

const MYCELIA_RUN_ALL = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
const MYCELIA_RUN_EXTERNAL = MYCELIA_RUN_ALL || lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"
const PROJECT_ROOT = dirname(@__DIR__)

# Helper function to include all test files in a directory
function include_all_tests(dir)
    test_count = 0
    for (root, dirs, files) in walkdir(dir)
        if occursin(".ipynb_checkpoints", root)
            continue
        end
        for file in sort(files)
            if endswith(file, ".jl")
                include(joinpath(root, file))
                test_count += 1
            end
        end
    end
    return test_count
end

# # Aqua.jl quality assurance tests
# include("aqua.jl")

# # JET.jl static analysis - uncomment to enable (can be slow)
# include("jet.jl")

# include_all_tests(joinpath(@__DIR__, "1_data_acquisition"))
# # For debugging individual suites, include explicit files instead of the directory sweep.
# # for file in (
# #     "1_data_acquisition/ncbi_download.jl",
# #     "1_data_acquisition/simulation_fasta.jl",
# #     "1_data_acquisition/simulation_fastq.jl",
# # )
# #     include(joinpath(@__DIR__, file))
# # end
# include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"))
# # For debugging individual suites, include explicit files instead of the directory sweep.
# # for file in (
# #     "2_preprocessing_qc/alphabets.jl",
# #     "2_preprocessing_qc/constants.jl",
# #     "2_preprocessing_qc/dimensionality_reduction_and_clustering.jl",
# #     "2_preprocessing_qc/preprocessing.jl",
# #     "2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl",
# #     "2_preprocessing_qc/sequence-complexity.jl",
# #     "2_preprocessing_qc/sequence_io.jl",
# # )
# #     include(joinpath(@__DIR__, file))
# # end
# include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"))
# # For debugging individual suites, include explicit files instead of the directory sweep.
# # for file in (
# #     "3_feature_extraction_kmer/kmer_analysis.jl",
# #     "3_feature_extraction_kmer/qualmer-analysis.jl",
# # )
# #     include(joinpath(@__DIR__, file))
# # end

# Stage 4 (assembly): include individual files until stable, then switch back to include_all_tests.
# include_all_tests(joinpath(@__DIR__, "4_assembly"))
for file in [
        # "test/4_assembly/aa_fasta_singlestrand_test.jl",
        # "test/4_assembly/aa_fastq_singlestrand_test.jl",
        # "test/4_assembly/aa_kmer_graph_test.jl",
        # "test/4_assembly/aa_kmer_singlestrand_test.jl",
        # "test/4_assembly/aa_qualmer_graph_test.jl",
        # "test/4_assembly/aa_qualmer_singlestrand_test.jl",
        # "test/4_assembly/amino_acid_fastq_test.jl",
        # "test/4_assembly/assembly_merging.jl",
        # "test/4_assembly/bandage_integration.jl",
        # "test/4_assembly/basic_graph_tests.jl",
        # "test/4_assembly/canonicalization_consistency_test.jl",
        # "test/4_assembly/complete_assembly_workflow_test.jl",
        # "test/4_assembly/comprehensive_correctness_tests.jl",
        # "test/4_assembly/comprehensive_fixes_tests.jl",
        # "test/4_assembly/comprehensive_type_stable_corrected_tests.jl",
        # "test/4_assembly/comprehensive_type_stable_tests.jl",
        # "test/4_assembly/dna_fasta_doublestrand_test.jl",
        # "test/4_assembly/dna_fasta_singlestrand_test.jl",
        # "test/4_assembly/dna_fastq_doublestrand_test.jl",
        # "test/4_assembly/dna_fastq_singlestrand_test.jl",
        # "test/4_assembly/dna_kmer_doublestrand_test.jl",
        # "test/4_assembly/dna_kmer_singlestrand_test.jl",
        # "test/4_assembly/dna_qualmer_doublestrand_test.jl",
        # "test/4_assembly/dna_qualmer_graph_test.jl",
        # "test/4_assembly/dna_qualmer_singlestrand_test.jl",
        # "test/4_assembly/doublestrand_canonicalization_test.jl",
        "test/4_assembly/end_to_end_assembly_tests.jl",
        "test/4_assembly/end_to_end_graph_tests.jl",
        "test/4_assembly/evidence_functions_test.jl",
        "test/4_assembly/evidence_structures_test.jl",
        "test/4_assembly/gfa_io_next.jl",
        "test/4_assembly/graph_algorithms_next.jl",
        "test/4_assembly/graph_conversion_test.jl",
        "test/4_assembly/graph_query_test.jl",
        "test/4_assembly/iterative_assembly_tests.jl",
        "test/4_assembly/kmer_edge_data_test.jl",
        "test/4_assembly/kmer_vertex_data_test.jl",
        "test/4_assembly/megahit_phix_workflow.jl",
        "test/4_assembly/path_finding_test.jl",
        "test/4_assembly/probabilistic_algorithms_next.jl",
        "test/4_assembly/quality_functions_test.jl",
        "test/4_assembly/rhizomorph_bubbles_and_gfa_test.jl",
        "test/4_assembly/rhizomorph_canonical_path_test.jl",
        "test/4_assembly/rhizomorph_conversion_errors_test.jl",
        "test/4_assembly/rhizomorph_doublestrand_files_test.jl",
        "test/4_assembly/rhizomorph_doublestrand_traversal_test.jl",
        "test/4_assembly/rhizomorph_kmer_mode_support_test.jl",
        "test/4_assembly/rhizomorph_qualmer_canonical_traversal_test.jl",
        "test/4_assembly/rhizomorph_qualmer_rc_evidence_test.jl",
        "test/4_assembly/rna_fasta_doublestrand_test.jl",
        "test/4_assembly/rna_fasta_singlestrand_test.jl",
        "test/4_assembly/rna_fastq_doublestrand_test.jl",
        "test/4_assembly/rna_fastq_singlestrand_test.jl",
        "test/4_assembly/rna_kmer_doublestrand_test.jl",
        "test/4_assembly/rna_kmer_graph_test.jl",
        "test/4_assembly/rna_kmer_singlestrand_test.jl",
        "test/4_assembly/rna_qualmer_doublestrand_test.jl",
        "test/4_assembly/rna_qualmer_graph_test.jl",
        "test/4_assembly/rna_qualmer_singlestrand_test.jl",
        "test/4_assembly/sequence_graphs_next.jl",
        "test/4_assembly/simplification_test.jl",
        "test/4_assembly/singlestrand_canonicalization_test.jl",
        "test/4_assembly/six_graph_hierarchy_tests.jl",
        "test/4_assembly/strand_specific_graph_construction_test.jl",
        "test/4_assembly/string-graph-helpers.jl",
        "test/4_assembly/string_graphs.jl",
        "test/4_assembly/string_ngram_singlestrand_quality_test.jl",
        "test/4_assembly/string_ngram_singlestrand_test.jl",
        "test/4_assembly/string_variable_singlestrand_quality_test.jl",
        "test/4_assembly/string_variable_singlestrand_test.jl",
        "test/4_assembly/tda_metrics_test.jl",
        "test/4_assembly/third_party_assemblers.jl",
        "test/4_assembly/third_party_assemblers_hybrid.jl",
        "test/4_assembly/third_party_assemblers_legacy.jl",
        "test/4_assembly/third_party_assemblers_long_read_isolate.jl",
        "test/4_assembly/third_party_assemblers_long_read_metagenomic.jl",
        "test/4_assembly/third_party_assemblers_plass_penguin.jl",
        "test/4_assembly/third_party_assemblers_short_read_isolate.jl",
        "test/4_assembly/third_party_assemblers_short_read_metagenomic.jl",
        "test/4_assembly/unicode-graph-assembly.jl",
        "test/4_assembly/variable_length_reduced_alphabet_test.jl",
        "test/4_assembly/variable_length_singlestrand_test.jl",
        "test/4_assembly/variable_length_strand_conversion_test.jl"
    ]
    include(joinpath(PROJECT_ROOT, file))
end

# Stage 5 (validation): focused suites (all other external-heavy validation stays opt-in for now).
include_all_tests(joinpath(@__DIR__, "5_validation")) # broken
# for file in (
#     # "5_validation/coverm_integration_extended.jl", # broken
#     "5_validation/coverm_wrappers.jl",
#     "5_validation/mosdepth_coverage_qc.jl",
#     # "5_validation/quast_busco_wrappers_test.jl", # need to add Glob to test dependencies
# )
#     include(file)
# end

# Stage 6 (annotation)
include_all_tests(joinpath(@__DIR__, "6_annotation")) # broken
# for file in (
#     # "6_annotation/annotation.jl", # broken
#     # "6_annotation/codon_optimization.jl", # broken
#     # "6_annotation/genome_features.jl", # broken
# )
#     include(joinpath(@__DIR__, file))
# end

# # Stage 7 (comparative/pangenomics): lightweight, synthetic-only suites by default.
include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"))
# for file in (
#     "7_comparative_pangenomics/distance_metrics.jl",
#     "7_comparative_pangenomics/pangenome.jl",
#     "7_comparative_pangenomics/pangenome_wrappers.jl",
#     "7_comparative_pangenomics/panproteome.jl",
#     "7_comparative_pangenomics/pantranscriptome.jl",
#     "7_comparative_pangenomics/phylogenetics.jl",
#     "7_comparative_pangenomics/sequence_classification.jl",
# )
#     include(file)
# end

# # Network/tooling-heavy suites (NCBI downloads, conda tools) are opt-in.
# if MYCELIA_RUN_EXTERNAL
#     for file in (
#         # "7_comparative_pangenomics/blastdb_integration.jl", # broken
#         "7_comparative_pangenomics/multiomics_alignment.jl",
#         "7_comparative_pangenomics/sequence_comparison.jl",
#     )
#         include(file)
#     end
# end

# metaphlan issue to fix
# Sun Dec 14 17:46:31 2025: [Warning] The number of reads in the sample (315) is below the recommended minimum of 10,000 reads.
# Sun Dec 14 17:46:44 2025: [Warning] Warning: No species were detected.

# metabuli to fix
# [ Info: Skipping Metabuli classify test; set METABULI_DB or METABULI_DB_PATH to a valid database directory

# external tool running issue
# WARNING: redefinition of constant Main.RUN_EXTERNAL. This may fail, cause incorrect answers, or produce other errors.
# [ Info: Skipping integration tests; set MYCELIA_RUN_EXTERNAL=true to run external tools

# sentencepiece isn't working
include_all_tests(joinpath(@__DIR__, "8_tool_integration"))
# for file in (
#     "8_tool_integration/bioconda.jl",
#     "8_tool_integration/binning_tools.jl",
#     "8_tool_integration/classification_tools.jl",
#     "8_tool_integration/metabuli_metaphlan_strainphlan.jl",
#     "8_tool_integration/minimap_merge_map_split.jl",
#     "8_tool_integration/pantools.jl",
#     "8_tool_integration/sentencepiece.jl",
#     "8_tool_integration/utility_functions.jl",
#     "8_tool_integration/xam.jl",
# )
#     include(joinpath(@__DIR__, file))
# end
