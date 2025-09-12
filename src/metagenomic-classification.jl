# """
# Metagenomic Classification with MAPQ-Aware Taxonomic Assignment

# This module implements advanced metagenomic classification techniques based on
# research from Mycelia-Dev metagenome analysis, specifically addressing the
# problem of MAPQ score misinterpretation in taxonomic assignment.

# Key Innovation: Use alignment scores instead of filtering MAPQ=0 reads for
# more accurate metagenomic profiling, as MAPQ=0 often indicates taxonomic
# ambiguity rather than poor alignment quality.
# """

# # =============================================================================
# # Core Data Structures
# # =============================================================================

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Represents a taxonomic assignment with confidence scoring.

# # Fields
# - `read_id::String`: Read identifier
# - `taxonomic_id::String`: Assigned taxonomic identifier
# - `taxonomic_level::Symbol`: Level of assignment (:species, :genus, :family, etc.)
# - `alignment_score::Float64`: Raw alignment score (NOT MAPQ)
# - `mapq_score::Int`: MAPQ score for reference (but not used for filtering)
# - `alignment_length::Int`: Length of alignment
# - `num_mismatches::Int`: Number of mismatches in alignment
# - `is_primary::Bool`: Whether this is the primary alignment
# - `secondary_scores::Vector{Float64}`: Scores of secondary alignments

# # Example
# ```julia
# assignment = TaxonomicAssignment(
#     "read_001", "species_123", :species, 85.5, 0, 150, 2, true, [78.2, 72.1]
# )
# ```
# """
# struct TaxonomicAssignment
#     read_id::String
#     taxonomic_id::String
#     taxonomic_level::Symbol  # :species, :genus, :family, :order, :class, :phylum, :kingdom
#     alignment_score::Float64  # Use this for weighting, NOT MAPQ
#     mapq_score::Int          # Store for reference but don't filter on it
#     alignment_length::Int
#     num_mismatches::Int
#     is_primary::Bool
#     secondary_scores::Vector{Float64}  # For ambiguity assessment
    
#     function TaxonomicAssignment(read_id::String, taxonomic_id::String, 
#                                taxonomic_level::Symbol, alignment_score::Float64,
#                                mapq_score::Int, alignment_length::Int, num_mismatches::Int,
#                                is_primary::Bool = true, secondary_scores::Vector{Float64} = Float64[])
#         new(read_id, taxonomic_id, taxonomic_level, alignment_score, mapq_score,
#             alignment_length, num_mismatches, is_primary, secondary_scores)
#     end
# end

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Represents weighted abundance data for metagenomic samples.

# # Fields
# - `taxonomic_assignments::Vector{TaxonomicAssignment}`: All taxonomic assignments
# - `abundance_matrix::Dict{String, Float64}`: Weighted abundance per taxon
# - `total_reads::Int`: Total number of classified reads
# - `unclassified_reads::Int`: Number of unclassified reads
# - `ambiguous_reads::Int`: Number of reads with MAPQ=0 (ambiguous but informative)

# # Example
# ```julia
# profile = MetagenomicProfile(assignments, abundances, 10000, 500, 1200)
# ```
# """
# struct MetagenomicProfile
#     taxonomic_assignments::Vector{TaxonomicAssignment}
#     abundance_matrix::Dict{String, Float64}  # taxonomic_id -> weighted abundance
#     total_reads::Int
#     unclassified_reads::Int
#     ambiguous_reads::Int  # MAPQ=0 reads that were retained and weighted
    
#     function MetagenomicProfile(assignments::Vector{TaxonomicAssignment})
#         abundance_matrix = calculate_weighted_abundances(assignments)
#         total_reads = length(assignments)
#         unclassified_reads = 0  # Would be calculated separately
#         ambiguous_reads = count(a -> a.mapq_score == 0, assignments)
        
#         new(assignments, abundance_matrix, total_reads, unclassified_reads, ambiguous_reads)
#     end
# end

# # =============================================================================
# # MAPQ-Aware Taxonomic Assignment Functions
# # =============================================================================

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Process SAM/BAM alignment records for metagenomic classification using alignment
# scores instead of MAPQ filtering.

# **Key Innovation**: Retains MAPQ=0 alignments (often filtered out) and weights
# them by alignment score, preventing loss of valuable taxonomic information.

# # Arguments
# - `alignment_records::Vector`: SAM/BAM alignment records (from XAM.jl or similar)
# - `taxonomic_database::Dict{String, String}`: Reference -> taxonomic_id mapping  
# - `min_alignment_score::Float64`: Minimum alignment score threshold (default: 50.0)
# - `min_alignment_length::Int`: Minimum alignment length (default: 50)

# # Returns
# - `Vector{TaxonomicAssignment}`: Weighted taxonomic assignments

# # Example
# ```julia
# # Load alignment records (pseudo-code)
# records = load_sam_records("sample.sam")
# tax_db = load_taxonomic_database("ncbi_taxonomy.tsv")

# # Process with MAPQ-aware method
# assignments = process_alignment_records(records, tax_db)

# # Count how many MAPQ=0 reads were retained
# mapq_zero_count = count(a -> a.mapq_score == 0, assignments)
# println("Retained \$mapq_zero_count MAPQ=0 reads for taxonomic profiling")
# ```

# # References
# Based on Mycelia-Dev metagenome README.md insights about Minimap2 MAPQ 
# score interpretation in taxonomic contexts.
# """
# function process_alignment_records(alignment_records::Vector, 
#                                  taxonomic_database::Dict{String, String};
#                                  min_alignment_score::Float64 = 50.0,
#                                  min_alignment_length::Int = 50)::Vector{TaxonomicAssignment}
    
#     assignments = TaxonomicAssignment[]
    
#     for record in alignment_records
#         # Extract key alignment properties (pseudo-code interface)
#         read_id = get_read_id(record)
#         reference_id = get_reference_id(record)
#         alignment_score = get_alignment_score(record)  # Use this, NOT MAPQ
#         mapq_score = get_mapq_score(record)  # Store but don't filter
#         alignment_length = get_alignment_length(record)
#         num_mismatches = get_num_mismatches(record)
#         is_primary = is_primary_alignment(record)
        
#         # Apply quality filters based on alignment properties, NOT MAPQ
#         if alignment_score < min_alignment_score || alignment_length < min_alignment_length
#             continue
#         end
        
#         # Look up taxonomic assignment
#         if haskey(taxonomic_database, reference_id)
#             taxonomic_id = taxonomic_database[reference_id]
            
#             # Create assignment with alignment score weighting
#             assignment = TaxonomicAssignment(
#                 read_id, taxonomic_id, :species,  # Default to species level
#                 alignment_score, mapq_score, alignment_length, 
#                 num_mismatches, is_primary
#             )
            
#             push!(assignments, assignment)
#         end
#     end
    
#     return assignments
# end

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Calculate weighted taxonomic abundances using alignment scores.

# Uses alignment scores (not MAPQ) to weight taxonomic assignments, providing
# more accurate relative abundance estimates by retaining ambiguous but 
# informative alignments.

# # Arguments
# - `assignments::Vector{TaxonomicAssignment}`: Taxonomic assignments with scores

# # Returns
# - `Dict{String, Float64}`: taxonomic_id -> weighted relative abundance

# # Algorithm
# 1. Group assignments by taxonomic ID
# 2. Sum alignment scores for each taxonomic group
# 3. Normalize by total alignment score sum
# 4. MAPQ=0 assignments contribute proportionally to their alignment scores

# # Example
# ```julia
# abundances = calculate_weighted_abundances(assignments)
# # abundances["species_123"] = 0.15  # 15% relative abundance
# ```
# """
# function calculate_weighted_abundances(assignments::Vector{TaxonomicAssignment})::Dict{String, Float64}
#     # Group by taxonomic ID and sum alignment scores
#     taxonomic_scores = Dict{String, Float64}()
    
#     for assignment in assignments
#         if haskey(taxonomic_scores, assignment.taxonomic_id)
#             taxonomic_scores[assignment.taxonomic_id] += assignment.alignment_score
#         else
#             taxonomic_scores[assignment.taxonomic_id] = assignment.alignment_score
#         end
#     end
    
#     # Calculate total score for normalization
#     total_score = sum(values(taxonomic_scores))
    
#     # Normalize to relative abundances
#     abundance_matrix = Dict{String, Float64}()
#     for (taxonomic_id, score) in taxonomic_scores
#         abundance_matrix[taxonomic_id] = score / total_score
#     end
    
#     return abundance_matrix
# end

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Identify ambiguous taxonomic assignments based on alignment score distributions.

# Analyzes the distribution of alignment scores to identify reads that map well
# to multiple taxa, providing insight into taxonomic ambiguity.

# # Arguments
# - `assignments::Vector{TaxonomicAssignment}`: All assignments for analysis
# - `ambiguity_threshold::Float64`: Score ratio threshold for ambiguity (default: 0.9)

# # Returns
# - `Vector{String}`: Read IDs with ambiguous taxonomic assignments

# # Algorithm
# Identifies reads where the ratio of second-best to best alignment score
# exceeds the ambiguity threshold, indicating potential taxonomic ambiguity.

# # Example
# ```julia
# ambiguous_reads = identify_ambiguous_assignments(assignments, 0.85)
# println("Found \$(length(ambiguous_reads)) ambiguous reads")
# ```
# """
# function identify_ambiguous_assignments(assignments::Vector{TaxonomicAssignment}, 
#                                        ambiguity_threshold::Float64 = 0.9)::Vector{String}
#     # Group assignments by read ID
#     read_assignments = Dict{String, Vector{TaxonomicAssignment}}()
    
#     for assignment in assignments
#         if haskey(read_assignments, assignment.read_id)
#             push!(read_assignments[assignment.read_id], assignment)
#         else
#             read_assignments[assignment.read_id] = [assignment]
#         end
#     end
    
#     ambiguous_reads = String[]
    
#     for (read_id, read_assigns) in read_assignments
#         if length(read_assigns) >= 2
#             # Sort by alignment score (descending)
#             sorted_assigns = sort(read_assigns, by = a -> a.alignment_score, rev = true)
            
#             # Calculate ratio of second-best to best score
#             if length(sorted_assigns) >= 2
#                 score_ratio = sorted_assigns[2].alignment_score / sorted_assigns[1].alignment_score
                
#                 if score_ratio >= ambiguity_threshold
#                     push!(ambiguous_reads, read_id)
#                 end
#             end
#         end
#     end
    
#     return ambiguous_reads
# end

# # =============================================================================
# # Multi-Database Partitioned Mapping Strategy
# # =============================================================================

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Represents a partitioned taxonomic database for memory-efficient mapping.

# # Fields
# - `name::String`: Database partition name (e.g., "prokaryotic", "viral")
# - `reference_file::String`: Path to reference sequences
# - `index_file::String`: Path to mapping index
# - `taxonomic_mapping::Dict{String, String}`: Reference -> taxonomic_id mapping
# - `memory_requirement::Int`: Estimated memory requirement in bytes

# # Example
# ```julia
# db_partition = DatabasePartition(
#     "prokaryotic", 
#     "/data/prokaryotic_refs.fasta", 
#     "/data/prokaryotic.idx",
#     prokaryotic_taxonomy,
#     8_000_000_000  # 8GB
# )
# ```
# """
# struct DatabasePartition
#     name::String
#     reference_file::String
#     index_file::String
#     taxonomic_mapping::Dict{String, String}
#     memory_requirement::Int  # Estimated memory in bytes
# end

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Perform sequential metagenomic mapping against partitioned databases.

# Maps reads against multiple taxonomic database partitions sequentially to
# manage memory usage and improve mapping specificity.

# # Arguments
# - `fastq_file::String`: Input reads file
# - `database_partitions::Vector{DatabasePartition}`: Ordered database partitions
# - `output_dir::String`: Output directory for results
# - `threads::Int`: Number of mapping threads (default: 8)
# - `available_memory_gb::Int`: Available system memory in GB (default: 16)

# # Returns
# - `Vector{MetagenomicProfile}`: Taxonomic profiles from each database partition

# # Algorithm
# 1. Map reads against each database partition sequentially
# 2. Extract unmapped reads from each round for next partition
# 3. Combine results with partition-specific weighting
# 4. Manage memory by loading/unloading database indices

# # Example
# ```julia
# partitions = [
#     prokaryotic_partition,
#     viral_partition, 
#     eukaryotic_partition
# ]

# profiles = sequential_partitioned_mapping(
#     "sample.fastq", partitions, "output/", 16, 32
# )
# ```

# # References
# Based on Mycelia-Dev metagenome multi-database partitioning strategy for
# memory optimization and improved taxonomic resolution.
# """
# function sequential_partitioned_mapping(fastq_file::String,
#                                        database_partitions::Vector{DatabasePartition},
#                                        output_dir::String;
#                                        threads::Int = 8,
#                                        available_memory_gb::Int = 16)::Vector{MetagenomicProfile}
    
#     if !isdir(output_dir)
#         mkpath(output_dir)
#     end
    
#     profiles = MetagenomicProfile[]
#     current_reads = fastq_file
    
#     for (i, partition) in enumerate(database_partitions)
#         @info "Processing database partition: $(partition.name)"
        
#         # Check memory requirements
#         required_memory_gb = partition.memory_requirement ÷ 1_000_000_000
#         if required_memory_gb > available_memory_gb
#             @warn "Partition $(partition.name) requires $(required_memory_gb)GB but only $(available_memory_gb)GB available"
#         end
        
#         # Map reads against current partition
#         mapped_sam = joinpath(output_dir, "$(partition.name)_mapped.sam")
#         unmapped_fastq = joinpath(output_dir, "$(partition.name)_unmapped.fastq")
        
#         # Run mapping (pseudo-code - would use actual mapping tool)
#         mapping_cmd = build_mapping_command(
#             current_reads, partition.reference_file, mapped_sam, unmapped_fastq, threads
#         )
#         run(mapping_cmd)
        
#         # Process mapping results
#         alignment_records = load_sam_records(mapped_sam)
#         assignments = process_alignment_records(alignment_records, partition.taxonomic_mapping)
#         profile = MetagenomicProfile(assignments)
        
#         push!(profiles, profile)
        
#         # Use unmapped reads for next partition
#         if i < length(database_partitions) && isfile(unmapped_fastq)
#             current_reads = unmapped_fastq
#         end
        
#         @info "Partition $(partition.name): $(length(assignments)) reads classified"
#     end
    
#     return profiles
# end

# # =============================================================================
# # Strain-Level Clustering with FastANI
# # =============================================================================

# """
#     $(DocStringExtensions.TYPEDSIGNATURES)

# Perform hierarchical strain clustering using FastANI with configurable thresholds.

# # Arguments
# - `genome_files::Vector{String}`: Paths to genome FASTA files
# - `ani_threshold::Float64`: ANI threshold for strain clustering (default: 99.5)
# - `output_dir::String`: Output directory for results
# - `threads::Int`: Number of threads for FastANI (default: 8)

# # Returns
# - `Dict{String, String}`: genome_file -> cluster_id mapping

# # Example
# ```julia
# strain_clusters = fastani_strain_clustering(
#     ["genome1.fasta", "genome2.fasta", "genome3.fasta"],
#     99.5, "clustering_output/", 16
# )
# ```

# # References
# Based on Mycelia-Dev metagenome strain resolution methodology using
# FastANI-based hierarchical clustering.
# """
# function fastani_strain_clustering(genome_files::Vector{String};
#                                   ani_threshold::Float64 = 99.5,
#                                   output_dir::String = "fastani_clustering",
#                                   threads::Int = 8)::Dict{String, String}
    
#     if !isdir(output_dir)
#         mkpath(output_dir)
#     end
    
#     # Create genome list file for FastANI
#     genome_list = joinpath(output_dir, "genome_list.txt")
#     open(genome_list, "w") do io
#         for genome_file in genome_files
#             println(io, genome_file)
#         end
#     end
    
#     # Run all-vs-all FastANI comparison
#     fastani_output = joinpath(output_dir, "fastani_matrix.txt")
#     fastani_cmd = `fastani --ql $genome_list --rl $genome_list -o $fastani_output -t $threads --matrix`
    
#     @info "Running FastANI all-vs-all comparison"
#     run(fastani_cmd)
    
#     # Parse ANI matrix
#     ani_matrix = parse_fastani_matrix(fastani_output, genome_files)
    
#     # Convert to distance matrix
#     distance_matrix = 1.0 .- (ani_matrix ./ 100.0)
    
#     # Perform hierarchical clustering
#     clusters = hierarchical_cluster_ani(distance_matrix, (100 - ani_threshold) / 100.0)
    
#     # Create cluster assignments
#     cluster_mapping = Dict{String, String}()
#     for (i, genome_file) in enumerate(genome_files)
#         cluster_id = "cluster_$(clusters[i])"
#         cluster_mapping[genome_file] = cluster_id
#     end
    
#     return cluster_mapping
# end

# # =============================================================================
# # Helper Functions
# # =============================================================================

# """
# Helper function to parse FastANI matrix output (implementation would depend on actual format)
# """
# function parse_fastani_matrix(fastani_output::String, genome_files::Vector{String})::Matrix{Float64}
#     n = length(genome_files)
#     ani_matrix = zeros(Float64, n, n)
    
#     # Implementation would parse actual FastANI output format
#     # This is a placeholder showing the expected structure
    
#     return ani_matrix
# end

# """
# Helper function for hierarchical clustering (would use actual clustering library)
# """
# function hierarchical_cluster_ani(distance_matrix::Matrix{Float64}, threshold::Float64)::Vector{Int}
#     # Implementation would use actual hierarchical clustering
#     # This is a placeholder showing the expected interface
    
#     n = size(distance_matrix, 1)
#     return collect(1:n)  # Placeholder: each genome in its own cluster
# end

# # Additional helper functions for SAM/BAM processing would be implemented
# # based on the specific alignment record format being used (XAM.jl, etc.)

# """
# Placeholder helper functions for alignment record processing.
# Actual implementation would depend on the alignment format library used.
# """
# function get_read_id(record)::String
#     # Implementation depends on alignment record type
#     return "placeholder_read_id"
# end

# function get_reference_id(record)::String
#     # Implementation depends on alignment record type  
#     return "placeholder_reference"
# end

# function get_alignment_score(record)::Float64
#     # Implementation depends on alignment record type
#     return 0.0
# end

# function get_mapq_score(record)::Int
#     # Implementation depends on alignment record type
#     return 0
# end

# function get_alignment_length(record)::Int
#     # Implementation depends on alignment record type
#     return 0
# end

# function get_num_mismatches(record)::Int
#     # Implementation depends on alignment record type
#     return 0
# end

# function is_primary_alignment(record)::Bool
#     # Implementation depends on alignment record type
#     return true
# end

# function load_sam_records(sam_file::String)::Vector
#     # Implementation would load actual SAM records
#     return []
# end

# function build_mapping_command(reads_file::String, reference_file::String, 
#                               output_sam::String, unmapped_fastq::String, 
#                               threads::Int)
#     # Implementation would build actual mapping command
#     return `echo "placeholder mapping command"`
# end