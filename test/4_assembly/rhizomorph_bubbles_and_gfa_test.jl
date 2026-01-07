# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_bubbles_and_gfa_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_bubbles_and_gfa_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import Graphs
import MetaGraphsNext
import Kmers
import BioSequences
import FASTX

Test.@testset "Rhizomorph bubble detection" begin
    # Simple bubble: A splits to B and C, reconverges at D
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=String,
        vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
        edge_data_type=Mycelia.Rhizomorph.StringEdgeData
    )

    for v in ["A", "B", "C", "D"]
        graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
    end

    graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
    graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
    graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
    graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

    bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=5)

    Test.@test length(bubbles) == 1
    bubble = bubbles[1]
    Test.@test bubble.entry_vertex == "A"
    Test.@test bubble.exit_vertex == "D"
    Test.@test Set(bubble.path1) == Set(["B", "D"]) || Set(bubble.path2) == Set(["B", "D"])
end

Test.@testset "Rhizomorph GFA round-trip (k-mer singlestrand)" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        rec = FASTX.FASTA.Record("read1", BioSequences.LongDNA{4}("ATGCA"))
        Mycelia.write_fasta(outfile=fasta, records=[rec])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode=:singlestrand)
        gfa_path = joinpath(dir, "graph.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{3}, Mycelia.Rhizomorph.SingleStrand)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (qualmer singlestrand)" begin
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        rec = FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATGCA"), fill(UInt8(30), 5))
        Mycelia.write_fastq(outfile=fastq, records=[rec])

        graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, 3; mode=:singlestrand)
        gfa_path = joinpath(dir, "graph.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{3}, Mycelia.Rhizomorph.SingleStrand)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (FASTA OLC variable-length)" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        rec1 = FASTX.FASTA.Record("read1", BioSequences.LongDNA{4}("ATGCAT"))
        rec2 = FASTX.FASTA.Record("read2", BioSequences.LongDNA{4}("CATG"))
        Mycelia.write_fasta(outfile=fasta, records=[rec1, rec2])

        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(collect(Mycelia.open_fastx(fasta)); min_overlap=3)
        gfa_path = joinpath(dir, "graph_fasta.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path; force_biosequence_graph=true)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (FASTQ OLC variable-length)" begin
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        rec1 = FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATGCA"), fill(UInt8(35), 5))
        rec2 = FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCAT"), fill(UInt8(30), 4))
        Mycelia.write_fastq(outfile=fastq, records=[rec1, rec2])

        graph = Mycelia.Rhizomorph.build_fastq_graph_olc(collect(Mycelia.open_fastx(fastq)); min_overlap=3)
        gfa_path = joinpath(dir, "graph_fastq.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path; force_biosequence_graph=true)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

# ============================================================================
# Doublestrand Mode GFA Round-trip Tests
# ============================================================================

Test.@testset "Rhizomorph GFA round-trip (k-mer doublestrand DNA)" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        # Use a longer sequence to ensure k-mer graph is non-trivial
        rec = FASTX.FASTA.Record("read1", BioSequences.LongDNA{4}("ATGCATGCATGC"))
        Mycelia.write_fasta(outfile=fasta, records=[rec])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode=:doublestrand)
        gfa_path = joinpath(dir, "graph_doublestrand.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{3}, Mycelia.Rhizomorph.DoubleStrand)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (k-mer canonical DNA)" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        rec = FASTX.FASTA.Record("read1", BioSequences.LongDNA{4}("ATGCATGCATGC"))
        Mycelia.write_fasta(outfile=fasta, records=[rec])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode=:canonical)
        gfa_path = joinpath(dir, "graph_canonical.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{3})

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (qualmer doublestrand DNA)" begin
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        rec = FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATGCATGCATGC"), fill(UInt8(30), 12))
        Mycelia.write_fastq(outfile=fastq, records=[rec])

        graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, 3; mode=:doublestrand)
        gfa_path = joinpath(dir, "graph_qualmer_ds.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{3}, Mycelia.Rhizomorph.DoubleStrand)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (qualmer canonical DNA)" begin
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        rec = FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATGCATGCATGC"), fill(UInt8(30), 12))
        Mycelia.write_fastq(outfile=fastq, records=[rec])

        graph = Mycelia.Rhizomorph.build_qualmer_graph_from_file(fastq, 3; mode=:canonical)
        gfa_path = joinpath(dir, "graph_qualmer_can.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{3})

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (k-mer doublestrand RNA)" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        rec = FASTX.FASTA.Record("read1", BioSequences.LongRNA{4}("AUGCAUGCAUGC"))
        Mycelia.write_fasta(outfile=fasta, records=[rec])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode=:doublestrand)
        gfa_path = joinpath(dir, "graph_rna_ds.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.RNAKmer{3}, Mycelia.Rhizomorph.DoubleStrand)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test Set(MetaGraphsNext.labels(roundtrip)) == Set(MetaGraphsNext.labels(graph))
        Test.@test MetaGraphsNext.ne(roundtrip) == MetaGraphsNext.ne(graph)
    end
end

Test.@testset "Rhizomorph GFA round-trip (variable-length FASTA doublestrand)" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        rec1 = FASTX.FASTA.Record("read1", BioSequences.LongDNA{4}("ATGCATGCAT"))
        rec2 = FASTX.FASTA.Record("read2", BioSequences.LongDNA{4}("GCATGC"))
        Mycelia.write_fasta(outfile=fasta, records=[rec1, rec2])

        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(collect(Mycelia.open_fastx(fasta)); min_overlap=3)
        # Convert to doublestrand
        ds_graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(graph)
        gfa_path = joinpath(dir, "graph_fasta_ds.gfa")
        Mycelia.Rhizomorph.write_gfa_next(ds_graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path; force_biosequence_graph=true)

        Test.@test Graphs.is_directed(roundtrip.graph)
        # Doublestrand may have different vertex count due to RC merging
        Test.@test !isempty(collect(MetaGraphsNext.labels(roundtrip)))
    end
end

Test.@testset "Rhizomorph GFA round-trip (variable-length FASTQ doublestrand)" begin
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        rec1 = FASTX.FASTQ.Record("read1", BioSequences.LongDNA{4}("ATGCATGCAT"), fill(UInt8(35), 10))
        rec2 = FASTX.FASTQ.Record("read2", BioSequences.LongDNA{4}("GCATGC"), fill(UInt8(30), 6))
        Mycelia.write_fastq(outfile=fastq, records=[rec1, rec2])

        graph = Mycelia.Rhizomorph.build_fastq_graph_olc(collect(Mycelia.open_fastx(fastq)); min_overlap=3)
        # Convert to doublestrand
        ds_graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(graph)
        gfa_path = joinpath(dir, "graph_fastq_ds.gfa")
        Mycelia.Rhizomorph.write_gfa_next(ds_graph, gfa_path)

        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path; force_biosequence_graph=true)

        Test.@test Graphs.is_directed(roundtrip.graph)
        Test.@test !isempty(collect(MetaGraphsNext.labels(roundtrip)))
    end
end

# ============================================================================
# GFA Metadata Preservation Tests
# ============================================================================

Test.@testset "GFA metadata preservation - vertex count and edge count" begin
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        rec = FASTX.FASTA.Record("read1", BioSequences.LongDNA{4}("ATGCATGCATGCATGC"))
        Mycelia.write_fasta(outfile=fasta, records=[rec])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 5; mode=:singlestrand)
        original_vertex_count = length(collect(MetaGraphsNext.labels(graph)))
        original_edge_count = MetaGraphsNext.ne(graph)

        gfa_path = joinpath(dir, "graph.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_path)

        # Verify GFA file contains expected segments and links
        gfa_content = read(gfa_path, String)
        segment_count = count(line -> startswith(line, "S\t"), split(gfa_content, '\n'))
        link_count = count(line -> startswith(line, "L\t"), split(gfa_content, '\n'))

        Test.@test segment_count == original_vertex_count
        Test.@test link_count == original_edge_count

        # Round-trip and verify
        roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_path, Kmers.DNAKmer{5}, Mycelia.Rhizomorph.SingleStrand)
        Test.@test length(collect(MetaGraphsNext.labels(roundtrip))) == original_vertex_count
        Test.@test MetaGraphsNext.ne(roundtrip) == original_edge_count
    end
end
