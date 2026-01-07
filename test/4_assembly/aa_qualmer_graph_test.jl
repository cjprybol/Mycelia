# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/aa_qualmer_graph_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/aa_qualmer_graph_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Amino Acid Qualmer Graph Construction Test
# Tests for AA qualmer graph construction.  
# Note: AA sequences with quality scores are rare in practice,
# but we support them for completeness.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "AA Qualmer Graph - Singlestrand" begin

    Test.@testset "AA Qualmer Graph Construction - Basic" begin
        seq = "MKVLW"
        qual = [30, 35, 32, 28, 40]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("protein_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        # Should have 3 unique AA k-mers
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 3
        Test.@test Mycelia.Rhizomorph.edge_count(graph) == 2

        # Verify AA k-mers are present
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("MKV"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("KVL"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.AAKmer{3}("VLW"))
    end

    Test.@testset "AA Qualmer Graph - Quality Evidence" begin
        seq = "MKVLW"
        qual = [30, 35, 32, 28, 40]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("protein_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.AAKmer{3}("MKV")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        Test.@test !isnothing(vertex_data)
        evidence = first(vertex_data.evidence["test"]["protein_001"])

        Test.@test evidence isa Mycelia.Rhizomorph.QualityEvidenceEntry
        Test.@test evidence.position == 1
        Test.@test evidence.strand == Mycelia.Rhizomorph.Forward
        # Phred+33 encoding
        Test.@test evidence.quality_scores == UInt8[63, 68, 65]
    end

    Test.@testset "AA Qualmer Graph - Query Functions" begin
        seq = "MKVLWALL"
        qual = [30 for _ in 1:length(seq)]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("protein_001", seq, qual_str)

        records = [record]
        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(records, 3; dataset_id="test")

        sources = Mycelia.Rhizomorph.get_all_sources(graph)
        sinks = Mycelia.Rhizomorph.get_all_sinks(graph)

        Test.@test Kmers.AAKmer{3}("MKV") in sources
        Test.@test Kmers.AAKmer{3}("ALL") in sinks
    end
end

println("✓ AA qualmer graph tests completed")
