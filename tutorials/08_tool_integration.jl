# # Tutorial 8: Tool Integration and Workflow Management
#
# This tutorial covers integration with external bioinformatics tools,
# workflow management, and creating comprehensive analysis pipelines.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - Integration with external bioinformatics tools
# - Workflow management and pipeline construction
# - HPC job submission and resource management
# - Data management and cloud integration
# - Error handling and quality control in pipelines
# - Reproducible research practices

# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/08_tool_integration.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import Random
import Statistics

Random.seed!(42)

# ## Part 1: External Tool Integration
#
# Modern bioinformatics relies on integration of multiple specialized tools.
# Understanding how to effectively combine tools is crucial for comprehensive analysis.

println("=== Tool Integration Tutorial ===")

println("Common Bioinformatics Tools:")
println("- Assembly: hifiasm, Canu, Flye")
println("- Annotation: Prodigal, Pyrodigal, Prodigal-gv, Augustus, MetaEuk")
println("- Alignment: BWA, minimap2, BLAST")
println("- Phylogenetics: IQ-TREE, RAxML, MrBayes")
println("- Visualization: Circos, IGV, Artemis")
println("- Quality Control: FastQC, Quast, BUSCO")

# ### Wrapper entry points
#
# The wrappers below expose thin Julia entry points for external tools. Use these
# examples as templates and replace paths with real inputs. These snippets are
# shown as documentation only and are not executed in docs builds.
# Autocycler 0.5.2 is a reproducible compatibility pin for bacterial isolates
# where mostly complete alternative long-read assemblies are expected; it is not
# a generic metagenome/eukaryote assembler or an automatically current release.
#
# ```julia
# Mycelia.install_autocycler()
# Mycelia.run_autocycler(
#     long_reads = "reads.fastq",
#     out_dir = "autocycler_out",
#     read_type = "ont_r10",
# )
# Mycelia.run_autocycler_polished(
#     long_reads = "reads.fastq",
#     short_reads_1 = "reads_R1.fastq",
#     short_reads_2 = "reads_R2.fastq",
#     out_dir = "autocycler_polished_out",
#     read_type = "ont_r10",
# )
#
# # High-level corrected multi-input workflows use typed sibling adapters.
# unicycler_result = Mycelia.Rhizomorph.assemble_unicycler_hybrid(
#     "reads_R1.fastq.gz",
#     "reads_R2.fastq.gz",
#     "reads_ont.fastq.gz";
#     config = Mycelia.Rhizomorph.UnicyclerHybridConfig(
#         output_dir = "hybrid_unicycler_out",
#         # Cumulative stable source + correction-output copy budget.
#         input_snapshot_byte_ceiling = 500_000_000_000,
#     ),
# )
# # Unicycler's mutable Conda environment is bound to the realized package
# # inventory for this run (name, version, build, and channel for every package).
# unicycler_toolchain = unicycler_result.assembly_stats["toolchain"]
# unicycler_toolchain["package_inventory_sha256"]
# unicycler_toolchain["packages"]
# autocycler_result = Mycelia.Rhizomorph.assemble_autocycler_polished(
#     "reads_R1.fastq.gz",
#     "reads_R2.fastq.gz",
#     "reads_hifi.fastq.gz";
#     config = Mycelia.Rhizomorph.AutocyclerPolishConfig(
#         long_read_tech = :pacbio_hifi,
#         autocycler_read_type = :pacbio_hifi,
#         output_dir = "hybrid_autocycler_out",
#         input_snapshot_byte_ceiling = 500_000_000_000,
#     ),
# )
# autocycler_result.assembly_stats["toolchain"]
#
# # Both high-level workflows bind the original source content and corrected
# # FASTQ path, size, and content to deterministic SHA-256 provenance.
# read_content = unicycler_result.assembly_stats["read_content_provenance"]
# read_content["source_inputs"]
# read_content["corrected_fastqs"]
# # Ceiling, free-space, write, or hash failures attempt exact identity-bound
# # partial-snapshot cleanup. Cleanup is never silent: a failure propagates when
# # no result can be returned or retains and reports exact evidence with the
# # result. High-level children consume the already-bound corrected snapshots
# # directly, without a second scratch copy. Standalone run_unicycler and
# # run_autocycler_polished calls still expose input_spool_parent and
# # input_spool_byte_ceiling when a direct wrapper needs bounded scratch.
#
# metaMDBG v1.4 accepts one or more HiFi FASTQs, or one or more ONT R10.4-or-
# later FASTQs with an explicit attestation. It is never a HiFi-plus-ONT or
# Illumina-plus-long adapter. Synchronous execution returns :complete only after
# semantic FASTA/GFA validation; collected or dry-run jobs return :planned.
# metamdbg_hifi_result = Mycelia.run_metamdbg(
#     hifi_reads = ["reads_hifi_1.fastq.gz", "reads_hifi_2.fastq.gz"],
# )
# @assert metamdbg_hifi_result.status == :complete
# metamdbg_hifi_result.provenance.package_inventory
# metamdbg_hifi_result.provenance.package_inventory_sha256
# metamdbg_ont_result = Mycelia.run_metamdbg(
#     ont_reads = ["reads_ont_1.fastq.gz", "reads_ont_2.fastq.gz"],
#     ont_r10_4_plus = true,
# )
# @assert metamdbg_ont_result.status == :complete
#
# Mycelia.run_bcalm(["reads_R1.fastq", "reads_R2.fastq"], "bcalm_out"; kmer_size = 31)
# Mycelia.ggcat_build("reads.fastq", "graph.lz4", 31)
# Mycelia.ggcat_query("graph.lz4", "queries.fasta", "ggcat_hits.tsv", 31)
#
# Mycelia.foldseek_easy_search("query_structures/", "target_db/", "foldseek_hits.m8")
# Mycelia.run_pantools(["--help"])
#
# Mycelia.install_prokrustean()
# Mycelia.prokrustean_build_graph("example.ebwt", "prokrustean.bin"; kmin = 21)
# ```

# ## Part 2: Bioconda Integration
#
# Bioconda provides a standardized way to install and manage
# bioinformatics software packages.

println("\n=== Bioconda Integration ===")

# ### Environment Management
#
# Create and manage conda environments for different tools

println("--- Environment Management ---")

# TODO: Implement bioconda environment management
# - Create tool-specific environments
# - Install packages with dependency resolution
# - Manage environment versions
# - Export and reproduce environments

# Example environment specification
environment_spec = Dict(
    "name" => "mycelia-analysis",
    "channels" => ["conda-forge", "bioconda"],
    "dependencies" => [
        "hifiasm",
        "prodigal",
        "prodigal-gv",
        "pyrodigal",
        "augustus",
        "metaeuk",
        "blast",
        "busco",
        "iqtree",
        "circos"
    ]
)

println("Environment Specification:")
for (key, value) in environment_spec
    println("  $key: $value")
end

# ### Tool Installation and Configuration
#
# Install and configure external tools

println("--- Tool Installation ---")

# TODO: Implement automated tool installation
# - Check tool availability
# - Install missing tools
# - Configure tool paths
# - Validate tool functionality

# ## Part 3: Workflow Management
#
# Create robust, reproducible workflows that integrate multiple tools

println("\n=== Workflow Management ===")

# ### Pipeline Architecture
#
# Design modular, scalable analysis pipelines

println("--- Pipeline Architecture ---")

# TODO: Implement workflow architecture
# - Modular task design
# - Dependency management
# - Resource allocation
# - Error handling and recovery

# Example workflow structure
workflow_steps = [
    Dict("name" => "quality_control", "tool" => "fastqc", "input" => "reads.fastq", "output" => "qc_report"),
    Dict("name" => "assembly", "tool" => "hifiasm", "input" => "reads.fastq", "output" => "contigs.fasta"),
    Dict("name" => "annotation", "tool" => "prodigal", "input" => "contigs.fasta", "output" => "genes.gff"),
    Dict("name" => "validation", "tool" => "busco", "input" => "contigs.fasta", "output" => "busco_results")
]

println("Workflow Steps:")
for (i, step) in enumerate(workflow_steps)
    println("  Step $i: $(step["name"]) ($(step["tool"]))")
end

# ### Dependency Management
#
# Handle complex dependencies between analysis steps

println("--- Dependency Management ---")

# TODO: Implement dependency management
# - Build dependency graphs
# - Topological sorting
# - Parallel execution where possible
# - Handle conditional dependencies

# ## Part 4: HPC Integration
#
# Submit jobs to high-performance computing clusters

println("\n=== HPC Integration ===")

# ### SLURM Integration
#
# Submit and manage SLURM jobs

println("--- SLURM Integration ---")

# TODO: Implement SLURM integration
# - Generate SLURM job scripts
# - Submit jobs with appropriate resources
# - Monitor job status
# - Handle job failures and resubmission

# Example SLURM job configuration
slurm_config = Dict(
    "job_name" => "mycelia_analysis",
    "partition" => "compute",
    "time" => "24:00:00",
    "memory" => "64G",
    "cpus" => 16,
    "array" => "1-10",
    "output" => "analysis_%A_%a.out",
    "error" => "analysis_%A_%a.err"
)

println("SLURM Configuration:")
for (key, value) in slurm_config
    println("  $key: $value")
end

# ### Resource Management
#
# Optimize resource usage for different analysis types

println("--- Resource Management ---")

# TODO: Implement resource management
# - Estimate resource requirements
# - Dynamic resource allocation
# - Resource monitoring
# - Cost optimization

# ## Part 5: Cloud Integration
#
# Use cloud platforms for scalable analysis

println("\n=== Cloud Integration ===")

# ### Cloud Storage
#
# Integrate with cloud storage services

println("--- Cloud Storage ---")

# TODO: Implement cloud storage integration
# - Upload/download data to/from cloud
# - Manage cloud storage costs
# - Implement data lifecycle policies
# - Ensure data security and privacy

# ### Cloud Computing
#
# Use cloud computing resources

println("--- Cloud Computing ---")

# TODO: Implement cloud computing integration
# - Launch cloud instances
# - Configure analysis environments
# - Monitor resource usage
# - Optimize costs

# ## Part 6: Database Integration
#
# Integrate with biological databases

println("\n=== Database Integration ===")

# ### NCBI Integration
#
# Download and process data from NCBI

println("--- NCBI Integration ---")

# TODO: Implement comprehensive NCBI integration
# - Programmatic data download
# - Metadata processing
# - Format conversion
# - Batch processing

# ### Custom Database Integration
#
# Work with local and custom databases

println("--- Custom Database Integration ---")

# TODO: Implement custom database integration
# - Database schema design
# - Data import/export
# - Query interfaces
# - Performance optimization

# ## Part 7: Quality Control and Validation
#
# Implement comprehensive quality control throughout workflows

println("\n=== Quality Control ===")

# ### Automated Quality Checks
#
# Implement automated quality control checkpoints

println("--- Automated Quality Checks ---")

# TODO: Implement automated QC
# - Input data validation
# - Intermediate result checking
# - Output quality assessment
# - Automated reporting

# ### Error Handling
#
# Robust error handling and recovery

println("--- Error Handling ---")

# TODO: Implement error handling
# - Comprehensive error detection
# - Automatic error recovery
# - Error logging and reporting
# - User notification systems

# ## Part 8: Reproducibility and Documentation
#
# Ensure analysis reproducibility

println("\n=== Reproducibility ===")

# ### Version Control
#
# Track software versions and parameters

println("--- Version Control ---")

# TODO: Implement version control
# - Track software versions
# - Record analysis parameters
# - Version control for results
# - Reproducible environment creation

# ### Documentation Generation
#
# Automatic documentation of analysis workflows

println("--- Documentation Generation ---")

# TODO: Implement documentation generation
# - Automatic workflow documentation
# - Parameter documentation
# - Result interpretation guides
# - Method citations

# ## Part 9: Performance Optimization
#
# Optimize workflow performance

println("\n=== Performance Optimization ===")

# ### Parallel Processing
#
# Implement parallel processing strategies

println("--- Parallel Processing ---")

# TODO: Implement parallel processing
# - Task parallelization
# - Data parallelization
# - Pipeline parallelization
# - Load balancing

# ### Memory Management
#
# Optimize memory usage

println("--- Memory Management ---")

# TODO: Implement memory optimization
# - Memory usage monitoring
# - Streaming data processing
# - Memory-efficient algorithms
# - Garbage collection optimization

# ## Part 10: User Interface and Visualization
#
# Create user-friendly interfaces

println("\n=== User Interface ===")

# ### Command Line Interface
#
# Comprehensive CLI for workflow management

println("--- Command Line Interface ---")

# TODO: Implement comprehensive CLI
# - Intuitive command structure
# - Interactive configuration
# - Progress reporting
# - Help and documentation

# ### Web Interface
#
# Web-based workflow management

println("--- Web Interface ---")

# TODO: Implement web interface
# - Workflow configuration UI
# - Real-time monitoring
# - Result visualization
# - User management

# ## Part 11: Testing and Validation
#
# Comprehensive testing strategies

println("\n=== Testing and Validation ===")

# ### Unit Testing
#
# Test individual components

println("--- Unit Testing ---")

# TODO: Implement unit testing
# - Test individual functions
# - Mock external dependencies
# - Test edge cases
# - Automated test execution

# ### Integration Testing
#
# Test complete workflows

println("--- Integration Testing ---")

# TODO: Implement integration testing
# - Test complete workflows
# - Test with real data
# - Performance testing
# - Stress testing

# ## Part 12: Deployment and Distribution
#
# Deploy workflows for production use

println("\n=== Deployment ===")

# ### Container Integration
#
# Package workflows in containers

println("--- Container Integration ---")

# TODO: Implement container integration
# - Create Docker containers
# - Singularity integration
# - Container orchestration
# - Container registries

# ### Package Distribution
#
# Distribute workflows as packages

println("--- Package Distribution ---")

# TODO: Implement package distribution
# - Package creation
# - Dependency management
# - Version management
# - Distribution channels

# ## Part 13: Case Studies
#
# Real-world workflow examples

println("\n=== Case Studies ===")

# ### Bacterial Genome Analysis
#
# Complete bacterial genome analysis pipeline

println("--- Bacterial Genome Analysis ---")

println("Bacterial Genome Pipeline:")
println("1. Quality control (FastQC)")
println("2. Assembly (hifiasm)")
println("3. Assembly validation (QUAST, BUSCO)")
println("4. Gene prediction (Prodigal, Pyrodigal)")
println("5. Functional annotation (BLAST, eggNOG)")
println("6. Comparative analysis (ANI, phylogeny)")
println("7. Visualization (Circos, IGV)")

# ### Viral Genome Analysis
#
# Viral genome analysis and classification

println("--- Viral Genome Analysis ---")

println("Viral Genome Pipeline:")
println("1. Host depletion")
println("2. Viral read identification")
println("3. Assembly (SPAdes, Canu)")
println("4. Virus classification (BLAST, vContact)")
println("5. Gene prediction (Prodigal-gv, GeneMarkS)")
println("6. Functional annotation")
println("7. Phylogenetic analysis")

# ### Metagenome Analysis
#
# Metagenomic analysis pipeline

println("--- Metagenome Analysis ---")

println("Metagenome Pipeline:")
println("1. Quality control and trimming")
println("2. Host removal")
println("3. Assembly (MEGAHIT, metaSPAdes)")
println("4. Binning (MetaBAT, CONCOCT)")
println("5. Bin quality assessment (CheckM)")
println("6. Taxonomic classification")
println("7. Functional annotation")
println("8. Comparative analysis")

# ## Part 14: Best Practices
#
# Guidelines for effective tool integration

println("\n=== Best Practices ===")

println("Workflow Design:")
println("- Start with simple, working pipelines")
println("- Use modular design for flexibility")
println("- Implement comprehensive error handling")
println("- Plan for scalability from the beginning")
println()
println("Tool Integration:")
println("- Use standardized file formats")
println("- Validate tool outputs")
println("- Handle tool version differences")
println("- Document tool-specific requirements")
println()
println("Resource Management:")
println("- Monitor resource usage")
println("- Optimize for your computing environment")
println("- Plan for data storage requirements")
println("- Consider cost implications")
println()
println("Reproducibility:")
println("- Version control everything")
println("- Document all parameters")
println("- Use containerization when possible")
println("- Provide example datasets")

# ## Summary
println("\n=== Tool Integration Summary ===")
println("Topics introduced in this survey:")
println("- External-tool wrapper entry points")
println("- Environment and dependency-management concepts")
println("- Workflow, HPC, cloud, and resource-management concepts")
println("- Quality-control, reproducibility, and testing concepts")
println("- Interface, deployment, and distribution concepts")
println()
println("Sections marked TODO are conceptual placeholders, not implemented exercises.")
println("You have reached the end of Tutorial 8; this survey does not certify")
println("completion of every workflow or every Mycelia tutorial.")
println()
println("Continue exploring the Mycelia package for advanced features")
println("and consider contributing to the project!")

nothing
