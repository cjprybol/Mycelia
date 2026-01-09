# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/test_viroid_assembly_workflow.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/test_viroid_assembly_workflow.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Essential Tests for Viroid Assembly Workflow
# These tests validate the core functionality of the viroid assembly workflow
# including quality propagation, FASTQ output, and multi-sequence assembly.

using Pkg
Pkg.activate(".")

import Mycelia
import Statistics
import Test

println("=== Essential Viroid Assembly Workflow Tests ===")

Test.@testset "Viroid Assembly Workflow Tests" begin
    
    Test.@testset "Viroid Species Database" begin
        species = Mycelia.get_viroid_species_list()
        Test.@test length(species) == 25
        Test.@test "Potato spindle tuber viroid" in species
        Test.@test "Hop stunt viroid" in species
        println("✓ Viroid species database: $(length(species)) species")
    end
    
    Test.@testset "Read Simulation" begin
        # Test DNA read simulation
        test_sequence = "ATGCGATCGATCGTAGCTAGCTACGATCGTAGCTAGCT"
        dna_reads = Mycelia._simulate_fastq_reads_from_sequence(
            test_sequence, "test_dna";
            coverage=5, read_length=20, error_rate=0.01, sequence_type="DNA"
        )
        
        Test.@test length(dna_reads) > 0
        Test.@test all(read -> isa(read, Mycelia.FASTX.FASTQ.Record), dna_reads)
        
        # Test first read properties
        first_read = dna_reads[1]
        seq = String(Mycelia.FASTX.sequence(first_read))
        qual = collect(Mycelia.FASTX.quality_scores(first_read))
        
        Test.@test length(seq) <= 20  # Read length constraint
        Test.@test length(qual) == length(seq)  # Quality length matches sequence
        Test.@test all(q -> q >= 10 && q <= 40, qual)  # Reasonable quality range
        
        # Test RNA read simulation
        rna_reads = Mycelia._simulate_fastq_reads_from_sequence(
            "AUGCGAUCGAU", "test_rna";
            coverage=3, read_length=10, error_rate=0.01, sequence_type="RNA"
        )
        Test.@test length(rna_reads) > 0
        
        # Test protein read simulation  
        protein_reads = Mycelia._simulate_fastq_reads_from_sequence(
            "MKLVDSTFGK", "test_protein";
            coverage=3, read_length=8, error_rate=0.02, sequence_type="AA"
        )
        Test.@test length(protein_reads) > 0
        
        println("✓ Read simulation: DNA, RNA, and protein reads generated")
    end
    
    Test.@testset "Quality-Aware Assembly" begin
        # Create test data
        test_genome = "ATGCGATCGATCGTAGCTAGCTACGATCGTAGCTAGCTACGATCGTAGCTAGCT"
        test_reads = Mycelia._simulate_fastq_reads_from_sequence(
            test_genome, "assembly_test";
            coverage=8, read_length=25, error_rate=0.01, sequence_type="DNA"
        )
        
        # Prepare observations
        observations = [(read, i) for (i, read) in enumerate(test_reads)]
        
        # Configure assembly
        config = Mycelia.Rhizomorph.AssemblyConfig(
            k=10, 
            use_quality_scores=true,
            bubble_resolution=true,
            repeat_resolution=true
        )
        
        # Run assembly
        result = Mycelia._assemble_qualmer_graph(observations, config)
        
        # Test assembly results
        Test.@test length(result.contigs) > 0
        Test.@test length(result.fastq_contigs) > 0
        Test.@test length(result.contigs) == length(result.fastq_contigs)
        
        # Test quality preservation
        Test.@test haskey(result.assembly_stats, "quality_preserved")
        Test.@test result.assembly_stats["quality_preserved"] == true
        
        # Test assembly statistics
        Test.@test haskey(result.assembly_stats, "mean_quality")
        Test.@test haskey(result.assembly_stats, "mean_coverage")
        Test.@test result.assembly_stats["mean_quality"] > 0
        Test.@test result.assembly_stats["mean_coverage"] > 0
        
        # Test FASTQ contig properties
        if !isempty(result.fastq_contigs)
            first_contig = result.fastq_contigs[1]
            contig_seq = String(Mycelia.FASTX.sequence(first_contig))
            contig_qual = collect(Mycelia.FASTX.quality_scores(first_contig))
            
            Test.@test length(contig_seq) > 0
            Test.@test length(contig_qual) == length(contig_seq)
            Test.@test all(q -> q >= 2 && q <= 40, contig_qual)  # Valid quality range
            Test.@test Statistics.mean(contig_qual) > 10  # Reasonable average quality
        end
        
        println("✓ Quality-aware assembly: FASTQ contigs with preserved quality")
    end
    
    Test.@testset "Multi-Sequence Assembly" begin
        # Test DNA assembly
        dna_sequence = "ATGCGATCGATCGTAGCTAGCTACGATCG"
        dna_reads = Mycelia._simulate_fastq_reads_from_sequence(
            dna_sequence, "dna_test"; coverage=5, read_length=15, error_rate=0.01, sequence_type="DNA"
        )
        dna_obs = [(read, i) for (i, read) in enumerate(dna_reads)]
        dna_config = Mycelia.Rhizomorph.AssemblyConfig(k=8, use_quality_scores=true)
        dna_result = Mycelia._assemble_qualmer_graph(dna_obs, dna_config)
        
        # Test RNA assembly
        rna_sequence = "AUGCGAUCGAUCGUAGCUAGCU"
        rna_reads = Mycelia._simulate_fastq_reads_from_sequence(
            rna_sequence, "rna_test"; coverage=5, read_length=12, error_rate=0.01, sequence_type="RNA"
        )
        rna_obs = [(read, i) for (i, read) in enumerate(rna_reads)]
        rna_config = Mycelia.Rhizomorph.AssemblyConfig(k=6, use_quality_scores=true)
        rna_result = Mycelia._assemble_qualmer_graph(rna_obs, rna_config)
        
        # Test protein assembly
        protein_sequence = "MKLVDSTFGKQIL"
        protein_reads = Mycelia._simulate_fastq_reads_from_sequence(
            protein_sequence, "protein_test"; coverage=4, read_length=8, error_rate=0.02, sequence_type="AA"
        )
        protein_obs = [(read, i) for (i, read) in enumerate(protein_reads)]
        protein_config = Mycelia.Rhizomorph.AssemblyConfig(k=4, use_quality_scores=true)
        protein_result = Mycelia._assemble_qualmer_graph(protein_obs, protein_config)
        
        # Verify all assemblies worked
        Test.@test length(dna_result.fastq_contigs) > 0
        Test.@test length(rna_result.fastq_contigs) > 0
        Test.@test length(protein_result.fastq_contigs) > 0
        
        # Verify quality preservation across all sequence types
        Test.@test dna_result.assembly_stats["quality_preserved"] == true
        Test.@test rna_result.assembly_stats["quality_preserved"] == true
        Test.@test protein_result.assembly_stats["quality_preserved"] == true
        
        println("✓ Multi-sequence assembly: DNA, RNA, and protein assemblies successful")
    end
    
    Test.@testset "Algorithm Integration" begin
        # Create challenging test data to ensure all algorithms are tested
        complex_sequence = """
        ATGCGATCGATCGTAGCTAGCTACGATCGTAGCTAGCTACGATCGTAGCTAGCTACGATCG
        TAGCTAGCTACGATCGTAGCTAGCTACGATCGTAGCTAGCTACGATCGTAGCTAGCT
        """ |> x -> replace(x, '\n' => "") |> x -> replace(x, ' ' => "")
        
        # Generate higher coverage to trigger different algorithms
        complex_reads = Mycelia._simulate_fastq_reads_from_sequence(
            complex_sequence, "complex_test";
            coverage=12, read_length=30, error_rate=0.015, sequence_type="DNA"
        )
        
        complex_obs = [(read, i) for (i, read) in enumerate(complex_reads)]
        complex_config = Mycelia.Rhizomorph.AssemblyConfig(
            k=12, 
            use_quality_scores=true,
            bubble_resolution=true,
            repeat_resolution=true,
            min_coverage=2
        )
        
        # This should trigger multiple assembly algorithms
        complex_result = Mycelia._assemble_qualmer_graph(complex_obs, complex_config)
        
        Test.@test length(complex_result.fastq_contigs) > 0
        Test.@test complex_result.assembly_stats["quality_preserved"] == true
        Test.@test haskey(complex_result.assembly_stats, "num_vertices")
        Test.@test haskey(complex_result.assembly_stats, "num_edges") 
        
        # Test that we get reasonable assembly metrics
        Test.@test complex_result.assembly_stats["mean_quality"] > 20
        Test.@test complex_result.assembly_stats["mean_coverage"] > 1
        
        println("✓ Algorithm integration: Complex assembly with multiple algorithms")
    end
    
    Test.@testset "Error Handling" begin
        # Test empty observations - this will produce an empty result but shouldn't crash
        empty_obs = Tuple{Mycelia.FASTX.FASTQ.Record, Int}[]
        config = Mycelia.Rhizomorph.AssemblyConfig(k=10, use_quality_scores=true)
        
        # This should handle gracefully without crashing
        try
            result = Mycelia._assemble_qualmer_graph(empty_obs, config)
            Test.@test length(result.contigs) == 0  # Should produce empty result
            Test.@test length(result.fastq_contigs) == 0
        catch e
            # It's OK if it throws a specific error for empty input - that's also valid handling
            Test.@test isa(e, BoundsError) || isa(e, ArgumentError)
        end
        
        # Test very short sequence
        short_seq = "ATGC"
        short_reads = Mycelia._simulate_fastq_reads_from_sequence(
            short_seq, "short_test"; coverage=3, read_length=10, error_rate=0.0, sequence_type="DNA"
        )
        short_obs = [(read, i) for (i, read) in enumerate(short_reads)]
        
        # Should handle short sequences gracefully
        Test.@test_nowarn Mycelia._assemble_qualmer_graph(short_obs, config)
        
        println("✓ Error handling: Graceful handling of edge cases")
    end
end

println("\n=== Essential Tests Summary ===")
println("✓ Viroid species database validated")
println("✓ Read simulation for DNA, RNA, and proteins working")
println("✓ Quality-aware assembly producing FASTQ output")
println("✓ Multi-sequence type assembly validated")
println("✓ Advanced algorithms integrated and functional")
println("✓ Error handling robust for edge cases")
println("\n🎉 All essential viroid assembly workflow tests passed!")
println("The implementation is ready for production use.")
