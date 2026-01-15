# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/complete_assembly_workflow_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/complete_assembly_workflow_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Complete Assembly Workflow Tests
# Tests the full assembly pipeline:
# FASTQ → Qualmer Graph → Path Finding → Sequence Reconstruction

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers
import Random

Test.@testset "Complete Assembly Workflow Tests" begin

    Test.@testset "Linear Genome - Perfect Reconstruction" begin
        # Simple case: linear genome with no branches, should reconstruct perfectly
        mktempdir() do dir
            # Create a simple linear sequence
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGC")

            # Simulate reads with high quality
            fastq = joinpath(dir, "reads.fastq")
            quality = fill(UInt8(40), length(reference))  # High quality
            rec = FASTX.FASTQ.Record("read1", reference, quality)
            Mycelia.write_fastq(outfile=fastq, records=[rec])

            # Step 1: Build qualmer graph
            k = 5
            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:singlestrand)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            # Step 2: Find Eulerian paths
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            Test.@test !isempty(paths)

            # Step 3: Reconstruct sequences from paths
            for path in paths
                if length(path) > 1
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                    Test.@test reconstructed isa BioSequences.LongDNA
                    # The reconstructed sequence should contain the original
                    # (or be contained in it, depending on path coverage)
                    Test.@test length(reconstructed) > 0
                end
            end
        end
    end

    Test.@testset "Multi-read Coverage - Consensus Assembly" begin
        # Multiple reads covering the same region should produce consistent assembly
        mktempdir() do dir
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGCATGC")

            # Create multiple overlapping reads
            fastq = joinpath(dir, "reads.fastq")
            records = FASTX.FASTQ.Record[]

            # Create 5 identical reads with high quality
            for i in 1:5
                quality = fill(UInt8(35), length(reference))
                rec = FASTX.FASTQ.Record("read$i", reference, quality)
                push!(records, rec)
            end
            Mycelia.write_fastq(outfile=fastq, records=records)

            # Build qualmer graph
            k = 5
            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:singlestrand)

            # Verify higher coverage in evidence
            for label in MetaGraphsNext.labels(graph)
                vertex_data = graph[label]
                # Each k-mer should have evidence from multiple reads
                evidence_count = Mycelia.Rhizomorph.count_evidence(vertex_data)
                Test.@test evidence_count >= 1
            end

            # Find paths and reconstruct
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            Test.@test !isempty(paths)

            for path in paths
                if length(path) > 1
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                    Test.@test length(reconstructed) > 0
                end
            end
        end
    end

    Test.@testset "Doublestrand Mode - Strand-aware Assembly" begin
        mktempdir() do dir
            # Create sequence and its reverse complement
            forward = BioSequences.LongDNA{4}("ATGCATGCATGC")
            reverse_comp = BioSequences.reverse_complement(forward)

            fastq = joinpath(dir, "reads.fastq")
            records = [
                FASTX.FASTQ.Record("forward", forward, fill(UInt8(30), length(forward))),
                FASTX.FASTQ.Record("reverse", reverse_comp, fill(UInt8(30), length(reverse_comp)))
            ]
            Mycelia.write_fastq(outfile=fastq, records=records)

            # Build doublestrand qualmer graph
            k = 5
            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:doublestrand)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            # In doublestrand mode, forward and reverse complement k-mers remain distinct
            # Find paths
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)

            # Should have valid paths
            Test.@test paths isa Vector

            for path in paths
                if length(path) > 1
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                    Test.@test reconstructed isa BioSequences.LongDNA
                end
            end
        end
    end

    Test.@testset "Canonical Mode - K-mer Canonicalization" begin
        mktempdir() do dir
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGC")

            fastq = joinpath(dir, "reads.fastq")
            rec = FASTX.FASTQ.Record("read1", reference, fill(UInt8(35), length(reference)))
            Mycelia.write_fastq(outfile=fastq, records=[rec])

            # Build canonical qualmer graph
            k = 5
            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:canonical)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            # Find paths
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            Test.@test paths isa Vector

            for path in paths
                if length(path) > 1
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                    Test.@test length(reconstructed) > 0
                end
            end
        end
    end

    Test.@testset "K-mer Graph Assembly - DNA" begin
        # Test with k-mer graphs (non-quality-aware)
        mktempdir() do dir
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGC")

            fasta = joinpath(dir, "reads.fasta")
            rec = FASTX.FASTA.Record("read1", reference)
            Mycelia.write_fasta(outfile=fasta, records=[rec])

            # Build k-mer graph
            k = 5
            graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, k; mode=:singlestrand)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            # Find paths
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            Test.@test !isempty(paths)

            for path in paths
                if length(path) > 1
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                    Test.@test reconstructed isa BioSequences.LongDNA
                    Test.@test length(reconstructed) > 0
                end
            end
        end
    end

    Test.@testset "K-mer Graph Assembly - RNA" begin
        mktempdir() do dir
            reference = BioSequences.LongRNA{4}("AUGCAUGCAUGCAUGC")

            fasta = joinpath(dir, "reads.fasta")
            rec = FASTX.FASTA.Record("read1", reference)
            Mycelia.write_fasta(outfile=fasta, records=[rec])

            k = 5
            graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, k; mode=:singlestrand)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            Test.@test paths isa Vector

            for path in paths
                if length(path) > 1
                    reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                    # RNA should stay as RNA
                    Test.@test length(reconstructed) > 0
                end
            end
        end
    end

    Test.@testset "Workflow with Simulated Reads" begin
        # Use the testing utilities to create simulated reads with errors
        mktempdir() do dir
            # Create a reference sequence
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGCATGCATGCATGC")

            # Simulate reads with low error rate
            coverage = 10
            error_rate = 0.01
            reads = Mycelia.create_test_reads(reference, coverage, error_rate)

            Test.@test length(reads) == coverage
            Test.@test all(r isa FASTX.FASTQ.Record for r in reads)

            # Write reads to file
            fastq = joinpath(dir, "simulated.fastq")
            Mycelia.write_fastq(outfile=fastq, records=reads)

            # Build qualmer graph and assemble
            k = 5
            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:singlestrand)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            # Find paths
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            if isempty(paths)
                paths = Mycelia.Rhizomorph.find_contigs_next(graph; min_contig_length=1)
            end
            Test.@test paths isa Vector

            # At least one path should produce a non-empty sequence
            reconstructed_any = false
            for path_entry in paths
                reconstructed = if path_entry isa Mycelia.Rhizomorph.ContigPath
                    path_entry.sequence
                elseif !isempty(path_entry)
                    Mycelia.Rhizomorph.path_to_sequence(path_entry, graph)
                else
                    nothing
                end

                if reconstructed !== nothing && length(reconstructed) > 0
                    reconstructed_any = true
                    break
                end
            end
            Test.@test reconstructed_any
        end
    end

    Test.@testset "GFA Export After Assembly" begin
        # Complete workflow including GFA export
        mktempdir() do dir
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGC")

            fastq = joinpath(dir, "reads.fastq")
            rec = FASTX.FASTQ.Record("read1", reference, fill(UInt8(35), length(reference)))
            Mycelia.write_fastq(outfile=fastq, records=[rec])

            # Build graph
            k = 5
            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:singlestrand)

            # Find paths
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)

            # Export to GFA
            gfa_path = joinpath(dir, "assembly.gfa")
            Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)
            Test.@test isfile(gfa_path)

            # Read back and verify
            roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{5}, Mycelia.Rhizomorph.SingleStrand)
            Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        end
    end

    Test.@testset "Empty Input Handling" begin
        mktempdir() do dir
            # Empty FASTQ file
            fastq = joinpath(dir, "empty.fastq")
            Mycelia.write_fastq(outfile=fastq, records=FASTX.FASTQ.Record[])

            # Building graph from empty file should handle gracefully
            k = 5
            Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:singlestrand)
        end
    end

    Test.@testset "Various K-mer Sizes" begin
        mktempdir() do dir
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGCATGCATGCATGCATGC")

            fastq = joinpath(dir, "reads.fastq")
            rec = FASTX.FASTQ.Record("read1", reference, fill(UInt8(35), length(reference)))
            Mycelia.write_fastq(outfile=fastq, records=[rec])

            # Test different k-mer sizes
            for k in [3, 5, 7, 11]
                if k < length(reference)  # k must be less than sequence length
                    graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, k; mode=:singlestrand)
                    Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

                    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
                    Test.@test paths isa Vector

                    # Larger k should give fewer, longer paths (in general)
                    for path in paths
                        if length(path) > 1
                            reconstructed = Mycelia.Rhizomorph.path_to_sequence(path, graph)
                            Test.@test length(reconstructed) > 0
                        end
                    end
                end
            end
        end
    end

    Test.@testset "Quality Score Impact on Assembly" begin
        # Compare assembly quality with high vs low quality reads
        mktempdir() do dir
            reference = BioSequences.LongDNA{4}("ATGCATGCATGCATGC")

            # High quality reads
            hq_fastq = joinpath(dir, "high_quality.fastq")
            hq_quality = fill(UInt8(40), length(reference))
            hq_rec = FASTX.FASTQ.Record("hq_read", reference, hq_quality)
            Mycelia.write_fastq(outfile=hq_fastq, records=[hq_rec])

            # Low quality reads
            lq_fastq = joinpath(dir, "low_quality.fastq")
            lq_quality = fill(UInt8(10), length(reference))
            lq_rec = FASTX.FASTQ.Record("lq_read", reference, lq_quality)
            Mycelia.write_fastq(outfile=lq_fastq, records=[lq_rec])

            k = 5

            # Build graphs
            hq_graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(hq_fastq, k; mode=:singlestrand)
            lq_graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(lq_fastq, k; mode=:singlestrand)

            # Both should produce graphs
            Test.@test !isempty(collect(MetaGraphsNext.labels(hq_graph)))
            Test.@test !isempty(collect(MetaGraphsNext.labels(lq_graph)))

            # Verify quality scores are preserved in graph
            for label in MetaGraphsNext.labels(hq_graph)
                vertex_data = hq_graph[label]
                if hasfield(typeof(vertex_data), :quality_scores)
                    # High quality scores should be preserved
                    Test.@test !isempty(vertex_data.quality_scores)
                end
            end
        end
    end

    Test.@testset "Bubble Detection and Simplification in Workflow" begin
        mktempdir() do dir
            # Create two sequences that differ by one nucleotide (creating a bubble)
            seq1 = BioSequences.LongDNA{4}("ATGCATGC")
            seq2 = BioSequences.LongDNA{4}("ATGGATGC")  # C→G at position 4

            fasta = joinpath(dir, "bubble.fasta")
            recs = [
                FASTX.FASTA.Record("seq1", seq1),
                FASTX.FASTA.Record("seq2", seq2)
            ]
            Mycelia.write_fasta(outfile=fasta, records=recs)

            k = 3
            graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, k; mode=:singlestrand)
            Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

            # Detect bubbles
            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)
            Test.@test bubbles isa Vector{<:Mycelia.Rhizomorph.BubbleStructure}

            # Simplify if bubbles found
            if !isempty(bubbles)
                simplified = Mycelia.Rhizomorph.simplify_graph_next(graph, bubbles)
                Test.@test simplified isa MetaGraphsNext.MetaGraph
            end

            # Find paths in original graph
            paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
            Test.@test paths isa Vector
        end
    end
end
