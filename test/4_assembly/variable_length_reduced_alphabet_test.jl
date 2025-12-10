# Variable-length FASTA/FASTQ graphs with reduced amino acid alphabets (as strings)

import Test
import Mycelia
import FASTX
import BioSequences

Test.@testset "Variable-length graphs with reduced AA alphabets" begin
    # FASTA: reduce AA to MURPHY2 and build string OLC graph via FASTA builder
    full_seq = BioSequences.LongAA("ACDEFGHIKLMNPQRSTVWY")
    reduced = Mycelia.reduce_amino_acid_alphabet(full_seq, :MURPHY2)
    fasta_records = [FASTX.FASTA.Record("aa_reduced", reduced)]
    fasta_graph = Mycelia.Rhizomorph.build_fasta_graph(fasta_records; dataset_id="aa_reduced", min_overlap=3)
    Test.@test Mycelia.Rhizomorph.vertex_count(fasta_graph) == 1
    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(fasta_graph)
    Test.@test !isempty(paths)

    # FASTQ: reduce AA, attach dummy qualities, and build quality-aware OLC graph
    qual_str = repeat("I", length(reduced))
    fastq_records = [FASTX.FASTQ.Record("aa_reduced_q", reduced, qual_str)]
    fastq_graph = Mycelia.Rhizomorph.build_fastq_graph(fastq_records; dataset_id="aa_reduced_q", min_overlap=3)
    Test.@test Mycelia.Rhizomorph.vertex_count(fastq_graph) == 1
    qpaths = Mycelia.Rhizomorph.find_eulerian_paths_next(fastq_graph)
    Test.@test !isempty(qpaths)
end

println("✓ Variable-length reduced alphabet tests completed")
