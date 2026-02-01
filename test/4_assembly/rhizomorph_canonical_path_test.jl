# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_canonical_path_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_canonical_path_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Rhizomorph canonical graph path reconstruction tests

import Test
import Mycelia
import FASTX
import Kmers
import BioSequences

Test.@testset "Rhizomorph canonical path reconstruction" begin
    # DNA k-mer graph -> canonical
    records = [FASTX.FASTA.Record("dna1", "ATGCA")]
    ss_graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id = "dna")
    canon_graph = Mycelia.Rhizomorph.convert_to_canonical(ss_graph)

    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon_graph)
    Test.@test !isempty(paths)

    # Ensure path_to_sequence can reconstruct
    recon_ok = false
    for p in paths
        seq = Mycelia.Rhizomorph.path_to_sequence(p, canon_graph)
        if seq !== nothing && seq isa BioSequences.LongDNA{4}
            recon_ok = true
            break
        end
    end
    Test.@test recon_ok

    # DNA qualmer graph -> canonical
    qual_record = FASTX.FASTQ.Record("read_qual", "ATGCA", "IIIII")
    ss_qgraph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([qual_record], 3; dataset_id = "dnaq")
    canon_qgraph = Mycelia.Rhizomorph.convert_to_canonical(ss_qgraph)
    qpaths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon_qgraph)
    Test.@test !isempty(qpaths)
    q_recon_ok = false
    for p in qpaths
        seq = Mycelia.Rhizomorph.path_to_sequence(p, canon_qgraph)
        if seq !== nothing && seq isa BioSequences.LongDNA{4}
            q_recon_ok = true
            break
        end
    end
    Test.@test q_recon_ok

    # RNA k-mer canonicalization should work similarly (no qual scores)
    rna_records = [FASTX.FASTA.Record("rna1", "AUGCA")]
    ss_rna = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(rna_records, 3; dataset_id = "rna")
    canon_rna = Mycelia.Rhizomorph.convert_to_canonical(ss_rna)
    rna_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon_rna)
    Test.@test !isempty(rna_paths)
end

println("✓ Rhizomorph canonical path tests completed")
