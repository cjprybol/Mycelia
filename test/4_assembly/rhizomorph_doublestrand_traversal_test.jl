# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_doublestrand_traversal_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_doublestrand_traversal_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Doublestrand traversal and reconstruction tests for DNA/RNA k-mer and qualmer graphs

import Test
import Mycelia
import FASTX
import Kmers
import BioSequences

Test.@testset "Rhizomorph doublestrand traversal" begin
    # DNA k-mer doublestrand
    dna_records = [FASTX.FASTA.Record("dna1", "ATGCA")]
    dna_ds = Mycelia.Rhizomorph.build_kmer_graph_doublestrand(dna_records, 3; dataset_id = "dna")
    dna_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(dna_ds)
    Test.@test !isempty(dna_paths)
    recon_ok = false
    for p in dna_paths
        seq = Mycelia.Rhizomorph.path_to_sequence(p, dna_ds)
        if seq !== nothing && seq isa BioSequences.LongDNA{4}
            recon_ok = true
            break
        end
    end
    Test.@test recon_ok

    # RNA k-mer doublestrand
    rna_records = [FASTX.FASTA.Record("rna1", "AUGCA")]
    rna_ds = Mycelia.Rhizomorph.build_kmer_graph_doublestrand(rna_records, 3; dataset_id = "rna")
    rna_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(rna_ds)
    Test.@test !isempty(rna_paths)

    # DNA qualmer doublestrand
    dq_record = FASTX.FASTQ.Record("dq1", "ATGCA", "IIIII")
    dq_ds = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand([dq_record], 3; dataset_id = "dnaq")
    dq_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(dq_ds)
    Test.@test !isempty(dq_paths)
    dq_recon_ok = false
    for p in dq_paths
        seq = Mycelia.Rhizomorph.path_to_sequence(p, dq_ds)
        if seq !== nothing && seq isa BioSequences.LongDNA{4}
            dq_recon_ok = true
            break
        end
    end
    Test.@test dq_recon_ok
end

println("✓ Rhizomorph doublestrand traversal tests completed")
