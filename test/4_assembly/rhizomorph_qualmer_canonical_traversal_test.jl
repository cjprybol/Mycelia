# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_qualmer_canonical_traversal_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_qualmer_canonical_traversal_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Canonical traversal and reconstruction for DNA/RNA qualmer graphs

import Test
import Mycelia
import FASTX
import BioSequences
import Kmers

Test.@testset "Rhizomorph qualmer canonical traversal" begin
    # DNA qualmer canonical
    dq_record = FASTX.FASTQ.Record("dq1", "ATGCA", "IIIII")
    ss_q = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([dq_record], 3; dataset_id="dnaq")
    canon_q = Mycelia.Rhizomorph.convert_to_canonical(ss_q)
    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon_q)
    Test.@test !isempty(paths)
    recon_ok = false
    for p in paths
        seq = Mycelia.Rhizomorph.path_to_sequence(p, canon_q)
        if seq !== nothing && seq isa BioSequences.LongDNA{4}
            recon_ok = true
            break
        end
    end
    Test.@test recon_ok

    # RNA qualmer canonical
    rq_record = FASTX.FASTQ.Record("rq1", "AUGCA", "IIIII")
    ss_rq = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([rq_record], 3; dataset_id="rnaq")
    canon_rq = Mycelia.Rhizomorph.convert_to_canonical(ss_rq)
    rpaths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon_rq)
    Test.@test !isempty(rpaths)
end

println("✓ Rhizomorph qualmer canonical traversal tests completed")
