# Mycelia Test Suite - Tiered Testing System
#
# Usage:
#   Core tests:                 julia --project=. -e "import Pkg; Pkg.test()"
#   Extended tests:             MYCELIA_RUN_ALL=true julia --project=. -e 'import Pkg; Pkg.test()'
#   Extended tests - portable:  LD_LIBRARY_PATH="" MYCELIA_RUN_ALL=true julia --project=. -e 'import Pkg; Pkg.update(); Pkg.instantiate(); Pkg.precompile(); Pkg.test()'
#   Benchmarks:                 julia --project=. benchmarking/run_all_benchmarks.jl
#   Tutorials:                  julia --project=. tutorials/run_all_tutorials.jl

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

# Aqua.jl quality assurance tests
include("aqua.jl")

# JET.jl static analysis - uncomment to enable (can be slow)
include("jet.jl")

include_all_tests(joinpath(@__DIR__, "1_data_acquisition"))
include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"))
include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"))

# add individual files until they all pass, then move to include_all_tests
# include_all_tests(joinpath(@__DIR__, "4_assembly"))
include("4_assembly/assembly_merging.jl")
include("4_assembly/third_party_assemblers.jl")
include("4_assembly/basic_graph_tests.jl")
include("4_assembly/amino_acid_fastq_test.jl")
include("4_assembly/canonicalization_consistency_test.jl")
include("4_assembly/doublestrand_canonicalization_test.jl")
include("4_assembly/singlestrand_canonicalization_test.jl")
include("4_assembly/comprehensive_fixes_tests.jl")
include("4_assembly/comprehensive_type_stable_corrected_tests.jl")
include("4_assembly/string-graph-helpers.jl")
include("4_assembly/end_to_end_graph_tests.jl") # may not yet be fully complete - review and confirm
include("4_assembly/comprehensive_correctness_tests.jl")
include("4_assembly/bandage_integration.jl")
include("4_assembly/megahit_phix_workflow.jl") # MEGAHIT/Bandage/Qualimap end-to-end workflow on PhiX
# TODO: fix
# include("4_assembly/end_to_end_assembly_tests.jl")

# Rhizomorph graph ecosystem tests (added 2025-12-10)
# Path finding and sequence reconstruction
include("4_assembly/path_finding_test.jl")
include("4_assembly/rhizomorph_doublestrand_traversal_test.jl")
include("4_assembly/rhizomorph_canonical_path_test.jl")
include("4_assembly/rhizomorph_qualmer_canonical_traversal_test.jl")
include("4_assembly/rhizomorph_qualmer_rc_evidence_test.jl")
include("4_assembly/rhizomorph_doublestrand_files_test.jl")
include("4_assembly/rhizomorph_conversion_errors_test.jl")

# K-mer graph construction tests
include("4_assembly/dna_kmer_singlestrand_test.jl")
include("4_assembly/dna_kmer_doublestrand_test.jl")
include("4_assembly/rna_kmer_singlestrand_test.jl")
include("4_assembly/rna_kmer_doublestrand_test.jl")
include("4_assembly/rna_kmer_graph_test.jl")
include("4_assembly/aa_kmer_singlestrand_test.jl")
include("4_assembly/aa_kmer_graph_test.jl")
include("4_assembly/kmer_vertex_data_test.jl")
include("4_assembly/kmer_edge_data_test.jl")

# Targeted Rhizomorph coverage
# - Mode support (directed singlestrand/doublestrand, undirected canonical)
# - Bubble detection + GFA round-trip
include("4_assembly/rhizomorph_kmer_mode_support_test.jl")
include("4_assembly/rhizomorph_bubbles_and_gfa_test.jl")

# # include_all_tests(joinpath(@__DIR__, "5_validation"))
# include_all_tests(joinpath(@__DIR__, "5_validation"))
# Focused validation suites
include("5_validation/coverm_wrappers.jl")
include("5_validation/coverm_integration_extended.jl")
include("5_validation/mosdepth_coverage_qc.jl")
include("5_validation/quast_busco_wrappers_test.jl")

include_all_tests(joinpath(@__DIR__, "6_annotation"))

# include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"))
# Comparative pangenomics: enable lightweight, synthetic-only suites by default.
include("7_comparative_pangenomics/distance_metrics.jl")
include("7_comparative_pangenomics/pangenome.jl")
include("7_comparative_pangenomics/panproteome.jl")
include("7_comparative_pangenomics/pantranscriptome.jl")
include("7_comparative_pangenomics/phylogenetics.jl")
include("7_comparative_pangenomics/sequence_classification.jl")
include("7_comparative_pangenomics/pangenome_wrappers.jl")

include("7_comparative_pangenomics/sequence_comparison.jl")   # sylph/skani, ART reads, NCBI fetch
include("7_comparative_pangenomics/multiomics_alignment.jl")  # badread + minimap2
include("7_comparative_pangenomics/blastdb_integration.jl")   # BLAST DB downloads

include_all_tests(joinpath(@__DIR__, "8_tool_integration"))
