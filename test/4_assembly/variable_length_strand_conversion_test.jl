# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/variable_length_strand_conversion_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/variable_length_strand_conversion_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Variable-length FASTA/FASTQ strand conversion tests

import Test
import Mycelia
import FASTX
import BioSequences

Test.@testset "Variable-length strand conversions" begin
    # FASTA doublestrand/canonical
    fasta_records = [FASTX.FASTA.Record("contig1", "ATGCAG")]
    ss_fasta = Mycelia.Rhizomorph.build_fasta_graph(fasta_records; dataset_id="fasta", min_overlap=3)
    ds_fasta = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(ss_fasta)
    Test.@test Mycelia.Rhizomorph.vertex_count(ds_fasta) >= 2
    canon_fasta = Mycelia.Rhizomorph.convert_variable_length_to_canonical(ss_fasta)
    Test.@test Mycelia.Rhizomorph.vertex_count(canon_fasta) == 1

    # FASTQ doublestrand/canonical
    fastq_records = [FASTX.FASTQ.Record("read1", "ATGCAG", "IIIIII")]
    ss_fastq = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="fastq", min_overlap=3)
    ds_fastq = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(ss_fastq)
    Test.@test Mycelia.Rhizomorph.vertex_count(ds_fastq) >= 2
    canon_fastq = Mycelia.Rhizomorph.convert_variable_length_to_canonical(ss_fastq)
    Test.@test Mycelia.Rhizomorph.vertex_count(canon_fastq) == 1

    # Non-nucleotide should error
    string_records = ["hello world"]
    string_graph = Mycelia.Rhizomorph.build_string_graph(string_records; dataset_id="str", min_overlap=3)
    Test.@test_throws ErrorException Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(string_graph)
    Test.@test_throws ErrorException Mycelia.Rhizomorph.convert_variable_length_to_canonical(string_graph)
end

println("✓ Variable-length strand conversion tests completed")
