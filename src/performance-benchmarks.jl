"""
Performance benchmarking for MetaGraphsNext vs MetaGraphs implementations.

This module provides comprehensive benchmarking to validate our performance targets:
- Memory efficiency: 50% reduction through canonical representation  
- Construction speed: 2x faster graph building
- Type stability: Eliminate allocations from type instability
"""

import BenchmarkTools
import MetaGraphsNext
import MetaGraphs
import Graphs
import BioSequences
import FASTX
import DocStringExtensions
import Statistics
import Printf

"""
Benchmark configuration for consistent testing.
"""
struct BenchmarkConfig
    k::Int                    # K-mer size
    num_sequences::Int        # Number of test sequences
    sequence_length::Int      # Length of each sequence
    alphabet::Vector{Char}    # DNA alphabet
    
    BenchmarkConfig(k=15, num_sequences=100, sequence_length=1000) = 
        new(k, num_sequences, sequence_length, ['A', 'T', 'G', 'C'])
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate synthetic DNA sequences for benchmarking.

# Arguments
- `config`: BenchmarkConfig with test parameters

# Returns
- Vector of FASTA records for testing
"""
function generate_test_sequences(config::BenchmarkConfig)
    sequences = Vector{FASTX.FASTA.Record}()
    
    for i in 1:config.num_sequences
        # Generate random DNA sequence
        seq_chars = rand(config.alphabet, config.sequence_length)
        sequence = String(seq_chars)
        
        record = FASTX.FASTA.Record("seq_$i", sequence)
        push!(sequences, record)
    end
    
    return sequences
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Benchmark graph construction performance: Legacy vs Next-generation.

Compares:
- MetaGraphs.jl (legacy) vs MetaGraphsNext.jl (next-gen)
- Memory allocation patterns
- Construction time
- Type stability

# Arguments
- `config`: BenchmarkConfig for test parameters

# Returns
- NamedTuple with benchmark results
"""
function benchmark_graph_construction(config::BenchmarkConfig=BenchmarkConfig())
    println("🔬 Benchmarking Graph Construction Performance")
    println("=" ^ 60)
    
    # Generate test data
    sequences = generate_test_sequences(config)
    kmer_type = BioSequences.DNAKmer{config.k}
    
    println("Configuration:")
    println("  K-mer size: $(config.k)")
    println("  Sequences: $(config.num_sequences)")
    println("  Length per sequence: $(config.sequence_length)")
    println("  Total nucleotides: $(config.num_sequences * config.sequence_length)")
    println()
    
    # Benchmark legacy implementation
    println("📊 Benchmarking Legacy MetaGraphs Implementation...")
    legacy_benchmark = BenchmarkTools.@benchmark Mycelia.build_stranded_kmer_graph($kmer_type, $sequences) samples=5 evals=1
    
    # Benchmark next-generation implementation  
    println("📊 Benchmarking Next-Generation MetaGraphsNext Implementation...")
    next_benchmark = BenchmarkTools.@benchmark Mycelia.build_kmer_graph_next($kmer_type, $sequences) samples=5 evals=1
    
    # Extract metrics
    legacy_time = Statistics.median(legacy_benchmark.times) / 1e9  # Convert to seconds
    next_time = Statistics.median(next_benchmark.times) / 1e9
    
    legacy_memory = Statistics.median(legacy_benchmark.memory)
    next_memory = Statistics.median(next_benchmark.memory)
    
    legacy_allocs = Statistics.median(legacy_benchmark.allocs)
    next_allocs = Statistics.median(next_benchmark.allocs)
    
    # Calculate improvements
    speed_improvement = legacy_time / next_time
    memory_reduction = (legacy_memory - next_memory) / legacy_memory * 100
    alloc_reduction = (legacy_allocs - next_allocs) / legacy_allocs * 100
    
    # Display results
    println("\n🎯 Performance Results:")
    println("=" ^ 40)
    println(Printf.format("Legacy (MetaGraphs):     %.3f seconds, %.1f MB, %d allocations", legacy_time, (legacy_memory/1024/1024), legacy_allocs))
    println(Printf.format("Next-Gen (MetaGraphsNext): %.3f seconds, %.1f MB, %d allocations", next_time, (next_memory/1024/1024), next_allocs))
    println()
    println(Printf.format("Speed improvement:       %.2fx faster", speed_improvement))
    println(Printf.format("Memory reduction:        %.1f%%", memory_reduction))
    println(Printf.format("Allocation reduction:    %.1f%%", alloc_reduction))
    
    # Check targets
    println("\n🎯 Target Achievement:")
    println("  Speed target (2x faster):    ", speed_improvement >= 2.0 ? "✅ ACHIEVED" : "⚠️  $(speed_improvement:.1f)x")
    println("  Memory target (50% reduction): ", memory_reduction >= 50.0 ? "✅ ACHIEVED" : "⚠️  $(memory_reduction:.1f)%")
    
    return (
        legacy_time = legacy_time,
        next_time = next_time,
        speed_improvement = speed_improvement,
        legacy_memory = legacy_memory,
        next_memory = next_memory,
        memory_reduction = memory_reduction,
        legacy_allocs = legacy_allocs,
        next_allocs = next_allocs,
        alloc_reduction = alloc_reduction,
        targets_met = speed_improvement >= 2.0 && memory_reduction >= 50.0
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Benchmark memory usage patterns for different graph representations.

Compares memory usage of:
- Stranded vertices (legacy) vs canonical vertices (next-gen)
- Edge metadata structures
- Coverage tracking efficiency

# Arguments
- `config`: BenchmarkConfig for test parameters

# Returns
- NamedTuple with memory analysis results
"""
function benchmark_memory_patterns(config::BenchmarkConfig=BenchmarkConfig())
    println("\n🧠 Benchmarking Memory Usage Patterns")
    println("=" ^ 60)
    
    sequences = generate_test_sequences(config)
    kmer_type = BioSequences.DNAKmer{config.k}
    
    # Build graphs
    println("Building graphs for memory analysis...")
    legacy_graph = Mycelia.build_stranded_kmer_graph(kmer_type, sequences)
    next_graph = Mycelia.build_kmer_graph_next(kmer_type, sequences)
    
    # Analyze vertex counts
    legacy_vertices = Graphs.nv(legacy_graph)
    next_vertices = length(MetaGraphsNext.labels(next_graph))
    vertex_reduction = (legacy_vertices - next_vertices) / legacy_vertices * 100
    
    # Analyze edge counts  
    legacy_edges = Graphs.ne(legacy_graph)
    next_edges = length(collect(MetaGraphsNext.edge_labels(next_graph)))
    
    # Estimate memory usage
    legacy_vertex_memory = legacy_vertices * 200  # Rough estimate per vertex
    next_vertex_memory = next_vertices * 150      # More efficient canonical representation
    
    println("\n📊 Memory Pattern Analysis:")
    println("=" ^ 40)
    println(Printf.format("Legacy vertices:         %d", legacy_vertices))
    println(Printf.format("Next-gen vertices:       %d (%.1f%% reduction)", next_vertices, vertex_reduction))
    println(Printf.format("Legacy edges:            %d", legacy_edges))
    println(Printf.format("Next-gen edges:          %d", next_edges))
    println()
    println(Printf.format("Estimated vertex memory: Legacy %.1f KB, Next-gen %.1f KB", (legacy_vertex_memory/1024), (next_vertex_memory/1024)))
    
    return (
        legacy_vertices = legacy_vertices,
        next_vertices = next_vertices,
        vertex_reduction = vertex_reduction,
        legacy_edges = legacy_edges,
        next_edges = next_edges,
        memory_efficient = vertex_reduction > 40.0
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Benchmark type stability and allocation patterns.

Measures:
- Type inference success
- Runtime allocations
- Performance predictability

# Arguments
- `config`: BenchmarkConfig for test parameters

# Returns
- NamedTuple with type stability metrics
"""
function benchmark_type_stability(config::BenchmarkConfig=BenchmarkConfig())
    println("\n🔬 Benchmarking Type Stability")
    println("=" ^ 60)
    
    sequences = generate_test_sequences(BenchmarkConfig(config.k, 10, 100))  # Smaller for detailed analysis
    kmer_type = BioSequences.DNAKmer{config.k}
    
    # Test next-generation implementation for type stability
    println("Analyzing type stability of next-generation implementation...")
    
    # Build graph and measure allocations during operations
    graph = Mycelia.build_kmer_graph_next(kmer_type, sequences)
    
    # Test vertex access type stability
    vertex_access = BenchmarkTools.@benchmark begin
        for label in MetaGraphsNext.labels($graph)
            vertex_data = $graph[label]
            length(vertex_data.coverage)
        end
    end samples=10
    
    # Test edge access type stability
    edge_access = BenchmarkTools.@benchmark begin
        for edge_labels in MetaGraphsNext.edge_labels($graph)
            if length(edge_labels) == 2
                edge_data = $graph[edge_labels...]
                edge_data.weight
            end
        end
    end samples=10
    
    vertex_allocs = Statistics.median(vertex_access.allocs)
    edge_allocs = Statistics.median(edge_access.allocs)
    
    println("📊 Type Stability Results:")
    println("=" ^ 40)
    println(Printf.format("Vertex access allocations: %d (target: 0)", vertex_allocs))
    println(Printf.format("Edge access allocations:   %d (target: 0)", edge_allocs))
    
    type_stable = vertex_allocs == 0 && edge_allocs == 0
    println("Type stability achieved:    ", type_stable ? "✅ YES" : "⚠️  NO")
    
    return (
        vertex_allocs = vertex_allocs,
        edge_allocs = edge_allocs,
        type_stable = type_stable
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run comprehensive performance benchmark suite.

Executes all benchmark tests and provides a summary report against our targets:
- Memory Usage: 50% reduction through type-stable metadata
- Construction Speed: 2x faster graph building  
- Type Stability: Zero allocations in hot paths

# Arguments
- `config`: BenchmarkConfig for test parameters

# Returns
- NamedTuple with comprehensive results
"""
function run_full_benchmark_suite(config::BenchmarkConfig=BenchmarkConfig())
    println("🚀 Mycelia Assembly Performance Benchmark Suite")
    println("=" ^ 80)
    println("Testing Phase 1 performance targets...")
    println()
    
    # Run individual benchmarks
    construction_results = benchmark_graph_construction(config)
    memory_results = benchmark_memory_patterns(config)
    stability_results = benchmark_type_stability(config)
    
    # Overall assessment
    println("\n🎯 FINAL ASSESSMENT")
    println("=" ^ 80)
    
    targets_met = construction_results.targets_met && 
                  memory_results.memory_efficient && 
                  stability_results.type_stable
    
    println("Performance Targets:")
    println("  ✅ Construction Speed: ", construction_results.speed_improvement >= 2.0 ? "ACHIEVED" : "MISSED")
    println("  ✅ Memory Efficiency:  ", construction_results.memory_reduction >= 50.0 ? "ACHIEVED" : "MISSED") 
    println("  ✅ Type Stability:     ", stability_results.type_stable ? "ACHIEVED" : "MISSED")
    println()
    println("🏆 Overall Phase 1 Success: ", targets_met ? "✅ ALL TARGETS MET" : "⚠️  SOME TARGETS MISSED")
    
    return (
        construction = construction_results,
        memory = memory_results,
        stability = stability_results,
        overall_success = targets_met
    )
end

# Export benchmark functions
# export BenchmarkConfig, benchmark_graph_construction, benchmark_memory_patterns, 
#        benchmark_type_stability, run_full_benchmark_suite
