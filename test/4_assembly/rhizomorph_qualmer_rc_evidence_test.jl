# Qualmer doublestrand/canonical tests with mixed datasets and RC evidence

import Test
import Mycelia
import FASTX
import BioSequences

Test.@testset "Qualmer RC evidence handling" begin
    # Two datasets: forward read and its reverse complement
    f_record = FASTX.FASTQ.Record("fwd", "ATGCA", "IIIII")
    rc_seq = String(reverse_complement(BioSequences.DNASeq("ATGCA")))
    rc_record = FASTX.FASTQ.Record("rev", rc_seq, "IIIII")

    ss = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([f_record], 3; dataset_id="fwd_ds")
    Mycelia.Rhizomorph.add_observations_to_graph!(ss, [rc_record], 3; dataset_id="rev_ds")

    ds = Mycelia.Rhizomorph.convert_to_doublestrand(ss)
    canon = Mycelia.Rhizomorph.convert_to_canonical(ss)

    # Check that evidence from both datasets exists and strand flags are preserved in doublestrand
    kmer = BioSequences.DNAKmer{3}("ATG")
    rc_kmer = BioSequences.reverse_complement(kmer)
    vdata_fwd = ds[kmer]
    vdata_rc = ds[rc_kmer]
    Test.@test haskey(vdata_fwd.evidence, "fwd_ds")
    Test.@test haskey(vdata_rc.evidence, "rev_ds")

    # Canonical graph should merge evidence
    canon_kmer = BioSequences.canonical(kmer)
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
