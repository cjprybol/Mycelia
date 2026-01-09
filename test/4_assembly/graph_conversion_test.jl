# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/graph_conversion_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/graph_conversion_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Graph Conversion Test
# Tests for converting between singlestrand and doublestrand graph modes.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "Graph Conversion" begin

    Test.@testset "Singlestrand → Doublestrand → Singlestrand (DNA K-mer)" begin
        # Create original singlestrand graph
        seq = "ATGCAT"
        record = FASTX.FASTA.Record("read_001", seq)

        original = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], 3; dataset_id="test")
        original_vertex_count = Mycelia.Rhizomorph.vertex_count(original)

        # Convert to doublestrand
        doublestrand = Mycelia.Rhizomorph.convert_to_doublestrand(original)
        ds_vertex_count = Mycelia.Rhizomorph.vertex_count(doublestrand)

        # Doublestrand includes reverse-complement vertices
        Test.@test ds_vertex_count >= original_vertex_count

        # Convert back to singlestrand
        recovered = Mycelia.Rhizomorph.convert_to_singlestrand(doublestrand)
        recovered_vertex_count = Mycelia.Rhizomorph.vertex_count(recovered)

        # Recovered should match original vertex count
        Test.@test recovered_vertex_count == original_vertex_count
    end

    Test.@testset "Evidence Preservation Through Conversion (DNA K-mer)" begin
        seq = "ATGC"
        record = FASTX.FASTA.Record("read_001", seq)

        original = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], 3; dataset_id="test")

        # Get original evidence count
        atg_kmer = Kmers.DNAKmer{3}("ATG")
        original_count = Mycelia.Rhizomorph.get_vertex_observation_count(original, atg_kmer)

        # Convert to doublestrand and back
        doublestrand = Mycelia.Rhizomorph.convert_to_doublestrand(original)
        recovered = Mycelia.Rhizomorph.convert_to_singlestrand(doublestrand)

        # Evidence should be preserved
        recovered_count = Mycelia.Rhizomorph.get_vertex_observation_count(recovered, atg_kmer)
        Test.@test recovered_count == original_count
    end

    Test.@testset "Quality Preservation Through Conversion (DNA Qualmer)" begin
        seq = "ATGC"
        qual = [30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        original = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record], 3; dataset_id="test")

        # Get original quality
        atg_kmer = Kmers.DNAKmer{3}("ATG")
        original_vertex = Mycelia.Rhizomorph.get_vertex_data(original, atg_kmer)
        original_has_quality = any(e -> e isa Mycelia.Rhizomorph.QualityEvidenceEntry,
                                      original_vertex.evidence["test"]["read_001"])

        # Convert to doublestrand and back
        doublestrand = Mycelia.Rhizomorph.convert_to_doublestrand(original)
        recovered = Mycelia.Rhizomorph.convert_to_singlestrand(doublestrand)

        # Quality should be preserved
        recovered_vertex = Mycelia.Rhizomorph.get_vertex_data(recovered, atg_kmer)
        recovered_has_quality = any(e -> e isa Mycelia.Rhizomorph.QualityEvidenceEntry,
                                       recovered_vertex.evidence["test"]["read_001"])
        Test.@test recovered_has_quality == original_has_quality
    end

    Test.@testset "Edge Preservation Through Conversion" begin
        seq = "ATGC"
        record = FASTX.FASTA.Record("read_001", seq)

        original = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], 3; dataset_id="test")
        original_edge_count = Mycelia.Rhizomorph.edge_count(original)

        # Convert to doublestrand and back
        doublestrand = Mycelia.Rhizomorph.convert_to_doublestrand(original)
        recovered = Mycelia.Rhizomorph.convert_to_singlestrand(doublestrand)

        # Edge count might differ due to RC edges, but should be related
        recovered_edge_count = Mycelia.Rhizomorph.edge_count(recovered)
        Test.@test recovered_edge_count >= original_edge_count
    end

    Test.@testset "RNA K-mer Conversion" begin
        seq = "AUGCAU"
        record = FASTX.FASTA.Record("read_001", seq)

        original = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], 3; dataset_id="test")
        doublestrand = Mycelia.Rhizomorph.convert_to_doublestrand(original)
        recovered = Mycelia.Rhizomorph.convert_to_singlestrand(doublestrand)

        # Vertex counts should match
        Test.@test Mycelia.Rhizomorph.vertex_count(recovered) == Mycelia.Rhizomorph.vertex_count(original)
    end

    Test.@testset "Palindromic K-mer Handling" begin
        # GCGC is palindromic
        seq = "GCGCGC"
        record = FASTX.FASTA.Record("read_001", seq)

        original = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], 4; dataset_id="test")
        doublestrand = Mycelia.Rhizomorph.convert_to_doublestrand(original)
        recovered = Mycelia.Rhizomorph.convert_to_singlestrand(doublestrand)

        # For palindromes, doublestrand does not add distinct RC vertices
        Test.@test Mycelia.Rhizomorph.vertex_count(recovered) == Mycelia.Rhizomorph.vertex_count(doublestrand)
    end

    Test.@testset "Direct Doublestrand Construction Matches Conversion" begin
        seq = "ATGCAT"
        record = FASTX.FASTA.Record("read_001", seq)

        # Build doublestrand directly
        direct_ds = Mycelia.Rhizomorph.build_kmer_graph_doublestrand([record], 3; dataset_id="test")

        # Build via conversion
        ss = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([record], 3; dataset_id="test")
        converted_ds = Mycelia.Rhizomorph.convert_to_doublestrand(ss)

        # Should have same vertex count
        Test.@test Mycelia.Rhizomorph.vertex_count(direct_ds) == Mycelia.Rhizomorph.vertex_count(converted_ds)

        # Should have same observations for explicit k-mers
        atg_kmer = Kmers.DNAKmer{3}("ATG")
        rc_atg = BioSequences.reverse_complement(atg_kmer)
        direct_count = Mycelia.Rhizomorph.get_vertex_observation_count(direct_ds, atg_kmer)
        converted_count = Mycelia.Rhizomorph.get_vertex_observation_count(converted_ds, atg_kmer)
        Test.@test direct_count == converted_count

        direct_rc_count = Mycelia.Rhizomorph.get_vertex_observation_count(direct_ds, rc_atg)
        converted_rc_count = Mycelia.Rhizomorph.get_vertex_observation_count(converted_ds, rc_atg)
        Test.@test direct_rc_count == converted_rc_count
    end
end

println("✓ Graph conversion tests completed")
