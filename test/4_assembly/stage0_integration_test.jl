# Unit tests for Stage 0 → corrector integration primitives (td-1do7).
# The skip mechanism: classify k-mers once, then skip reads with no weak k-mer.

import Test
import Mycelia
import FASTX

rec(id, seq) = FASTX.FASTQ.Record(id, seq, String(fill('I', length(seq))))

Test.@testset "Stage 0 skip primitives" begin
    k = 5
    clean_seq = "ATGCGTACGTACGT"
    recs = [rec("clean$i", clean_seq) for i in 1:8]
    err_seq = "ATGTGTACGTACGT"                       # single substitution at pos 4
    push!(recs, rec("err1", err_seq))
    graph = Mycelia.Rhizomorph.build_qualmer_graph(recs, k;
        dataset_id = "ds", mode = :canonical, memory_profile = :full)

    solid = Mycelia._solid_kmer_set(graph)
    Test.@test !isempty(solid)

    # The clean read (coverage 9 on every k-mer) has no weak k-mer → skippable.
    Test.@test Mycelia._read_is_all_solid(rec("c", clean_seq), k, solid)
    # The error read carries coverage-1 error k-mers → NOT all-solid → corrected.
    Test.@test !Mycelia._read_is_all_solid(rec("e", err_seq), k, solid)
    # A read shorter than k has no k-mers → never skipped on absent evidence.
    Test.@test !Mycelia._read_is_all_solid(rec("s", "ACG"), k, solid)
    # Empty solid set → never skip.
    Test.@test !Mycelia._read_is_all_solid(rec("c", clean_seq), k, Set(eltype(solid)[]))

    # Wiring: improve_read_set_likelihood with skip_solid=true runs and returns
    # every read; the all-solid clean reads pass through unchanged.
    out_skip, _ = Mycelia.improve_read_set_likelihood(recs, graph, k;
        graph_mode = :canonical, skip_solid = true)
    Test.@test length(out_skip) == length(recs)
    for i in 1:8
        Test.@test FASTX.sequence(String, out_skip[i]) == clean_seq   # solid ⇒ untouched
    end
    # skip_solid=false path still works (no regression).
    out_noskip, _ = Mycelia.improve_read_set_likelihood(recs, graph, k;
        graph_mode = :canonical, skip_solid = false)
    Test.@test length(out_noskip) == length(recs)
end
