# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_qualmer_rc_evidence_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_qualmer_rc_evidence_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Qualmer doublestrand/canonical tests with mixed datasets and RC evidence

import Test
import Mycelia
import FASTX
import BioSequences
import Kmers

Test.@testset "Qualmer RC evidence handling" begin
    # Two datasets: forward read and its reverse complement
    f_record = FASTX.FASTQ.Record("fwd", "ATGCA", "IIIII")
    rc_seq = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}("ATGCA")))
    rc_record = FASTX.FASTQ.Record("rev", rc_seq, "IIIII")

    ss = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([f_record], 3; dataset_id="fwd_ds")
    Mycelia.Rhizomorph.add_observations_to_graph!(ss, [rc_record], 3; dataset_id="rev_ds")

    ds = Mycelia.Rhizomorph.convert_to_doublestrand(ss)
    canon = Mycelia.Rhizomorph.convert_to_canonical(ss)

    # Check that evidence from both datasets exists and strand flags are preserved in doublestrand
    kmer = Kmers.DNAKmer{3}("ATG")
    rc_kmer = BioSequences.reverse_complement(kmer)
    vdata_fwd = ds[kmer]
    vdata_rc = ds[rc_kmer]
    Test.@test haskey(vdata_fwd.evidence, "fwd_ds")
    Test.@test haskey(vdata_rc.evidence, "rev_ds")

    # Canonical graph should merge evidence
    canon_kmer = Kmers.canonical(kmer)
    cv = canon[canon_kmer]
    Test.@test haskey(cv.evidence, "fwd_ds")
    Test.@test haskey(cv.evidence, "rev_ds")

    # Traversal and reconstruction on canonical should still succeed
    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon)
    Test.@test !isempty(paths)
    recon_ok = false
    for p in paths
        seq = Mycelia.Rhizomorph.path_to_sequence(p, canon)
        if seq !== nothing && seq isa BioSequences.LongDNA{4}
            recon_ok = true
            break
        end
    end
    Test.@test recon_ok
end

println("✓ Qualmer RC evidence handling tests completed")
