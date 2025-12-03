# Conversion error handling tests for Rhizomorph graphs

import Test
import Mycelia
import FASTX
import Kmers

Test.@testset "Rhizomorph conversion errors" begin
    # AA k-mer graph should not support doublestrand/canonical
    aa_records = [FASTX.FASTA.Record("aa1", "MKVLW")] 
    aa_graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(aa_records, 3; dataset_id="aa")
    Test.@test_throws Error Mycelia.Rhizomorph.convert_to_doublestrand(aa_graph)
    Test.@test_throws Error Mycelia.Rhizomorph.convert_to_canonical(aa_graph)

    # N-gram (String) graph should not support doublestrand/canonical
    strings = ["hello world"]
    ngram_graph = Mycelia.Rhizomorph.build_ngram_graph(strings, 3; dataset_id="txt")
    Test.@test_throws Error Mycelia.Rhizomorph.convert_to_doublestrand(ngram_graph)
    Test.@test_throws Error Mycelia.Rhizomorph.convert_to_canonical(ngram_graph)
end

println("✓ Rhizomorph conversion error tests completed")
