## From the Mycelia base directory, run the tests with:
## 
## ```bash
## julia --project=. -e 'include("test/4_assembly/basic_graph_tests.jl")'
## ```
##
## And to turn this file into a jupyter notebook, run:
## ```bash
## julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/basic_graph_tests.jl", "test/4_assembly", execute=false)'
## ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext
import Graphs
import Kmers

Test.@testset "Basic Graph Type Construction Tests" begin
    
    ## Test data
    dna_seq = "ATCGATCGATCGATCG"
    rna_seq = "AUCGAUCGAUCGAUCG"
    protein_seq = "ALAVALINEGLUTAMINE"
    high_quality = "HHHHHHHHHHHHHHHH"
    
    Test.@testset "1. N-gram Graphs" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph([dna_seq], 3; dataset_id="ng_basic")
        Test.@test !isempty(MetaGraphsNext.labels(graph))
        println("✓ N-gram Graph: $(length(MetaGraphsNext.labels(graph))) vertices")
    end
    
    Test.@testset "2. K-mer Graphs" begin
        
        Test.@testset "DNA K-mer Graph" begin
            records = [FASTX.FASTA.Record("test", dna_seq)]
            graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand(records, 5; dataset_id="dna_k")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            println("✓ DNA K-mer Graph: $(length(MetaGraphsNext.labels(graph))) vertices")
        end
        
        Test.@testset "RNA K-mer Graph" begin
            records = [FASTX.FASTA.Record("test", rna_seq)]
            graph = Mycelia.Rhizomorph.build_kmer_graph_doublestrand(records, 4; dataset_id="rna_k")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            println("✓ RNA K-mer Graph: $(length(MetaGraphsNext.labels(graph))) vertices")
        end
        
        Test.@testset "Amino Acid K-mer Graph" begin
            records = [FASTX.FASTA.Record("test", protein_seq)]
            graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="aa_k")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            println("✓ Amino Acid K-mer Graph: $(length(MetaGraphsNext.labels(graph))) vertices")
        end
    end
    
    Test.@testset "3. Qualmer Graphs" begin
        
        Test.@testset "DNA Qualmer Graph" begin
            records = [FASTX.FASTQ.Record("test", dna_seq, high_quality)]
            graph = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand(records, 5; dataset_id="dna_q")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            first_kmer = first(MetaGraphsNext.labels(graph))
            vertex_data = graph[first_kmer]
            Test.@test Mycelia.Rhizomorph.count_evidence(vertex_data) > 0
            println("✓ DNA Qualmer Graph: $(length(MetaGraphsNext.labels(graph))) vertices")
        end
        
        Test.@testset "RNA Qualmer Graph" begin
            records = [FASTX.FASTQ.Record("test", rna_seq, high_quality)]
            graph = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand(records, 4; dataset_id="rna_q")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            first_kmer = first(MetaGraphsNext.labels(graph))
            vertex_data = graph[first_kmer]
            Test.@test Mycelia.Rhizomorph.count_evidence(vertex_data) > 0
            println("✓ RNA Qualmer Graph: $(length(MetaGraphsNext.labels(graph))) vertices")
        end
    end
    
    Test.@testset "4. FASTA Graphs" begin
        Test.@testset "Direct BioSequence Graph" begin
            records = [FASTX.FASTA.Record("test", dna_seq)]
            graph = Mycelia.Rhizomorph.build_fasta_graph(records; dataset_id="fasta_basic")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            
            ## Test vertex data
            sequences = collect(MetaGraphsNext.labels(graph))
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            println("✓ FASTA Graph: $(length(sequences)) BioSequences")
        end
    end
    
    Test.@testset "5. FASTQ Graphs" begin
        Test.@testset "Direct Quality BioSequence Graph" begin
            records = [FASTX.FASTQ.Record("test", dna_seq, high_quality)]
            graph = Mycelia.Rhizomorph.build_fastq_graph(records; dataset_id="fastq_basic")
            Test.@test !isempty(MetaGraphsNext.labels(graph))
            
            ## Test vertex data
            sequences = collect(MetaGraphsNext.labels(graph))
            Test.@test all(seq -> seq isa BioSequences.LongDNA, sequences)
            println("✓ FASTQ Graph: $(length(sequences)) quality BioSequences")
        end
    end

    Test.@testset "GFA Parsing and Structure" begin
        lines = [
            "H\tVN:Z:1.0",
            "S\t1\tACGT\t*\tDP:i:1",
            "S\t2\tCGTA\t*\tDP:i:1",
            "L\t1\t+\t2\t+\t3M",
        ]
        gfa = joinpath(tempdir(), "simple.gfa")
        open(gfa, "w") do io
            for l in lines
                println(io, l)
            end
        end

        g = Mycelia.Rhizomorph.read_gfa_next(gfa, Kmers.DNAKmer{4}, Mycelia.DoubleStrand)
        Test.@test g isa MetaGraphsNext.MetaGraph
        Test.@test length(MetaGraphsNext.labels(g)) == 2
        Test.@test length(collect(MetaGraphsNext.edge_labels(g))) == 1

        ## Re-implementing the structure table logic
        components = Graphs.connected_components(g)
        Test.@test length(components) == 1
        
        component = first(components)
        contig_indices = component
        Test.@test contig_indices == [1, 2]
        all_labels = collect(MetaGraphsNext.labels(g))
        Test.@test all_labels == [Kmers.DNAKmer{4}("ACGT"), Kmers.DNAKmer{4}("CGTA")]
        Test.@test length.(all_labels) == [4, 4]

        g = Mycelia.Rhizomorph.read_gfa_next(gfa, Kmers.DNAKmer{4}, Mycelia.SingleStrand)
        Test.@test g isa MetaGraphsNext.MetaGraph
        Test.@test length(MetaGraphsNext.labels(g)) == 2
        Test.@test length(collect(MetaGraphsNext.edge_labels(g))) == 1

        ## Re-implementing the structure table logic
        components = Graphs.connected_components(g)
        Test.@test length(components) == 1
        
        component = first(components)
        contig_indices = component
        Test.@test contig_indices == [1, 2]
        all_labels = collect(MetaGraphsNext.labels(g))
        Test.@test all_labels == [BioSequences.LongDNA{4}("ACGT"), BioSequences.LongDNA{4}("CGTA")]
        Test.@test length.(all_labels) == [4, 4]
    end
end
