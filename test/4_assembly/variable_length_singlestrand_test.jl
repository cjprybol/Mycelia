# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/variable_length_singlestrand_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/variable_length_singlestrand_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Variable-length FASTA/FASTQ singlestrand OLC graph tests

import Test
import Mycelia
import FASTX
import BioSequences

Test.@testset "Variable-length OLC graphs (singlestrand)" begin
    # FASTA OLC graph
    fasta_records = [
        FASTX.FASTA.Record("contig1", "ATGCGT"),
        FASTX.FASTA.Record("contig2", "GCGTAA")
    ]
    fasta_graph = Mycelia.Rhizomorph.build_fasta_graph(fasta_records; dataset_id="olc", min_overlap=3)
    Test.@test Mycelia.Rhizomorph.vertex_count(fasta_graph) == 2
    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(fasta_graph)
    Test.@test !isempty(paths)

    # FASTQ OLC graph (quality-aware)
    fastq_records = [
        FASTX.FASTQ.Record("read1", "ATGCGT", "IIIIII"),
        FASTX.FASTQ.Record("read2", "GCGTAA", "IIIIII")
    ]
    fastq_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="olc_q", min_overlap=3)
    Test.@test Mycelia.Rhizomorph.vertex_count(fastq_graph) == 2
    qpaths = Mycelia.Rhizomorph.find_eulerian_paths_next(fastq_graph)
    Test.@test !isempty(qpaths)
end

println("✓ Variable-length OLC singlestrand tests completed")
