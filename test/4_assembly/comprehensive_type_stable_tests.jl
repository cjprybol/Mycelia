# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/comprehensive_type_stable_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/comprehensive_type_stable_tests.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "Comprehensive Type-Stable Reconstruction Tests" begin
    ## Test helper function to validate round-trip reconstruction
    function test_round_trip_reconstruction(graph, expected_type)
        Test.@test !isempty(MetaGraphsNext.labels(graph))

        # Find Eulerian paths and test reconstruction
        paths = Mycelia.find_eulerian_paths_next(graph)
        Test.@test !isempty(paths)

        for path_vector in paths
            if !isempty(path_vector)
                # Create WalkStep objects from the path labels
                walk_steps = [
                    Mycelia.WalkStep(vertex_label, Mycelia.Forward, 1.0, Float64(i))
                    for (i, vertex_label) in enumerate(path_vector)
                ]

                # Create GraphPath from WalkStep objects
                graph_path = Mycelia.GraphPath(walk_steps)

                # Test type-stable reconstruction
                reconstructed = Mycelia.path_to_sequence(graph_path, graph)
                Test.@test reconstructed !== nothing

                # Verify the type is preserved
                if expected_type == :dna
                    Test.@test isa(reconstructed, BioSequences.LongDNA)
                elseif expected_type == :rna
                    Test.@test isa(reconstructed, BioSequences.LongRNA)
                elseif expected_type == :aa
                    Test.@test isa(reconstructed, BioSequences.LongAA)
                elseif expected_type == :string
                    Test.@test isa(reconstructed, String)
                end
            end
        end
    end

    # Test sequences for each type
    dna_seq = BioSequences.dna"ATCGATCGATCGATCG"
    rna_seq = BioSequences.rna"AUCGAUCGAUCGAUCG"
    aa_seq = BioSequences.aa"MKTVRQERLKSIVRIL"

    # Quality scores for quality-aware tests
    quality_scores = [30, 35, 40, 30, 35, 40, 30, 35, 40, 30, 35, 40, 30, 35, 40, 30]

    Test.@testset "DNA Graph Types" begin
        # DNA K-mer graphs (SingleStrand + DoubleStrand)
        Test.@testset "DNA K-mer SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_kmer_type = Kmers.DNAKmer{5}
            dna_graph = Mycelia.build_kmer_graph_next(dna_kmer_type, dna_reads; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(dna_graph, dna_kmer_type)
        end

        Test.@testset "DNA K-mer DoubleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_kmer_type = Kmers.DNAKmer{5}
            dna_graph = Mycelia.build_kmer_graph_next(dna_kmer_type, dna_reads; graph_mode=Mycelia.DoubleStrand)
            test_round_trip_reconstruction(dna_graph, dna_kmer_type)
        end

        # DNA BioSequence graphs (SingleStrand + DoubleStrand)
        Test.@testset "DNA BioSequence SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_graph = Mycelia.build_biosequence_graph_next(BioSequences.LongDNA{4}, dna_reads; graph_mode=Mycelia.SingleStrand)
            # Use a dummy type for BioSequence graphs since they don't use k-mers
            test_round_trip_reconstruction(dna_graph, BioSequences.LongDNA{4})
        end

        Test.@testset "DNA BioSequence DoubleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_graph = Mycelia.build_biosequence_graph_next(BioSequences.LongDNA{4}, dna_reads; graph_mode=Mycelia.DoubleStrand)
            test_round_trip_reconstruction(dna_graph, BioSequences.LongDNA{4})
        end

        # DNA Quality-aware graphs (SingleStrand + DoubleStrand)
        Test.@testset "DNA Quality K-mer SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            # Create quality-aware reads
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(dna_reads)
                seq = FASTX.sequence(BioSequences.LongDNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_dna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                dna_qual_graph = Mycelia.build_qualmer_graph_next(Kmers.DNAKmer{5}, quality_reads; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(dna_qual_graph, Kmers.DNAKmer{5})
            end
        end

        Test.@testset "DNA Quality K-mer DoubleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(dna_reads)
                seq = FASTX.sequence(BioSequences.LongDNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_dna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                dna_qual_graph = Mycelia.build_qualmer_graph_next(Kmers.DNAKmer{5}, quality_reads; graph_mode=Mycelia.DoubleStrand)
                test_round_trip_reconstruction(dna_qual_graph, Kmers.DNAKmer{5})
            end
        end

        Test.@testset "DNA Quality BioSequence SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(dna_reads)
                seq = FASTX.sequence(BioSequences.LongDNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_dna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                dna_qual_bio_graph = Mycelia.build_qualmer_biosequence_graph_next(BioSequences.LongDNA{4}, quality_reads; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(dna_qual_bio_graph, BioSequences.LongDNA{4})
            end
        end

        Test.@testset "DNA Quality BioSequence DoubleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(dna_reads)
                seq = FASTX.sequence(BioSequences.LongDNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_dna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                dna_qual_bio_graph = Mycelia.build_qualmer_biosequence_graph_next(BioSequences.LongDNA{4}, quality_reads; graph_mode=Mycelia.DoubleStrand)
                test_round_trip_reconstruction(dna_qual_bio_graph, BioSequences.LongDNA{4})
            end
        end
    end

    Test.@testset "RNA Graph Types" begin
        # RNA K-mer graphs (SingleStrand + DoubleStrand)
        Test.@testset "RNA K-mer SingleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_kmer_type = Kmers.RNAKmer{5}
            rna_graph = Mycelia.build_kmer_graph_next(rna_kmer_type, rna_reads; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(rna_graph, rna_kmer_type)
        end

        Test.@testset "RNA K-mer DoubleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_kmer_type = Kmers.RNAKmer{5}
            rna_graph = Mycelia.build_kmer_graph_next(rna_kmer_type, rna_reads; graph_mode=Mycelia.DoubleStrand)
            test_round_trip_reconstruction(rna_graph, rna_kmer_type)
        end

        # RNA BioSequence graphs (SingleStrand + DoubleStrand)
        Test.@testset "RNA BioSequence SingleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_graph = Mycelia.build_biosequence_graph_next(BioSequences.LongRNA{4}, rna_reads; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(rna_graph, BioSequences.LongRNA{4})
        end

        Test.@testset "RNA BioSequence DoubleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_graph = Mycelia.build_biosequence_graph_next(BioSequences.LongRNA{4}, rna_reads; graph_mode=Mycelia.DoubleStrand)
            test_round_trip_reconstruction(rna_graph, BioSequences.LongRNA{4})
        end

        # RNA Quality-aware graphs (SingleStrand + DoubleStrand)
        Test.@testset "RNA Quality K-mer SingleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(rna_reads)
                seq = FASTX.sequence(BioSequences.LongRNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_rna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                rna_qual_graph = Mycelia.build_qualmer_graph_next(Kmers.RNAKmer{5}, quality_reads; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(rna_qual_graph, Kmers.RNAKmer{5})
            end
        end

        Test.@testset "RNA Quality K-mer DoubleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(rna_reads)
                seq = FASTX.sequence(BioSequences.LongRNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_rna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                rna_qual_graph = Mycelia.build_qualmer_graph_next(Kmers.RNAKmer{5}, quality_reads; graph_mode=Mycelia.DoubleStrand)
                test_round_trip_reconstruction(rna_qual_graph, Kmers.RNAKmer{5})
            end
        end

        Test.@testset "RNA Quality BioSequence SingleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(rna_reads)
                seq = FASTX.sequence(BioSequences.LongRNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_rna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                rna_qual_bio_graph = Mycelia.build_qualmer_biosequence_graph_next(BioSequences.LongRNA{4}, quality_reads; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(rna_qual_bio_graph, BioSequences.LongRNA{4})
            end
        end

        Test.@testset "RNA Quality BioSequence DoubleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(rna_reads)
                seq = FASTX.sequence(BioSequences.LongRNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_rna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                rna_qual_bio_graph = Mycelia.build_qualmer_biosequence_graph_next(BioSequences.LongRNA{4}, quality_reads; graph_mode=Mycelia.DoubleStrand)
                test_round_trip_reconstruction(rna_qual_bio_graph, BioSequences.LongRNA{4})
            end
        end
    end

    Test.@testset "Amino Acid Graph Types" begin
        # Amino Acid K-mer graphs (SingleStrand only)
        Test.@testset "AA K-mer SingleStrand" begin
            aa_reads = [FASTX.FASTA.Record("aa_seq", aa_seq)]
            aa_kmer_type = Kmers.AAKmer{5}
            aa_graph = Mycelia.build_kmer_graph_next(aa_kmer_type, aa_reads; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(aa_graph, aa_kmer_type)
        end

        # Amino Acid BioSequence graphs (SingleStrand only)
        Test.@testset "AA BioSequence SingleStrand" begin
            aa_reads = [FASTX.FASTA.Record("aa_seq", aa_seq)]
            aa_graph = Mycelia.build_biosequence_graph_next(BioSequences.LongAA, aa_reads; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(aa_graph, BioSequences.LongAA)
        end

        # Amino Acid Quality-aware graphs (SingleStrand only)
        Test.@testset "AA Quality K-mer SingleStrand" begin
            aa_reads = [FASTX.FASTA.Record("aa_seq", aa_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(aa_reads)
                seq = FASTX.sequence(BioSequences.LongAA, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_aa_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                aa_qual_graph = Mycelia.build_qualmer_graph_next(Kmers.AAKmer{5}, quality_reads; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(aa_qual_graph, Kmers.AAKmer{5})
            end
        end

        Test.@testset "AA Quality BioSequence SingleStrand" begin
            aa_reads = [FASTX.FASTA.Record("aa_seq", aa_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(aa_reads)
                seq = FASTX.sequence(BioSequences.LongAA, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_aa_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                aa_qual_bio_graph = Mycelia.build_qualmer_biosequence_graph_next(BioSequences.LongAA, quality_reads; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(aa_qual_bio_graph, BioSequences.LongAA)
            end
        end
    end

    Test.@testset "String Graph Types" begin
        # String N-gram graphs (SingleStrand only)
        Test.@testset "String N-gram SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", "ABCDEFGHIJKLMNOP")]
            string_graph = Mycelia.build_string_graph_next(string_reads, ngram_length=5; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(string_graph, String)
        end

        # String BioSequence-like graphs (SingleStrand only) - treating strings as sequences
        Test.@testset "String BioSequence SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", "ABCDEFGHIJKLMNOP")]
            # For strings, we can use build_biosequence_graph_next with String type
            string_seq_graph = Mycelia.build_biosequence_graph_next(String, string_reads; graph_mode=Mycelia.SingleStrand)
            test_round_trip_reconstruction(string_seq_graph, String)
        end

        # String Quality-aware graphs (SingleStrand only)
        Test.@testset "String Quality N-gram SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", "ABCDEFGHIJKLMNOP")]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(string_reads)
                seq_str = FASTX.sequence(String, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq_str)]])
                fastq_record = FASTX.FASTQ.Record("qual_string_$i", seq_str, qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                string_qual_graph = Mycelia.build_qualmer_string_graph_next(quality_reads, ngram_length=5; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(string_qual_graph, String)
            end
        end

        Test.@testset "String Quality BioSequence SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", "ABCDEFGHIJKLMNOP")]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(string_reads)
                seq_str = FASTX.sequence(String, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq_str)]])
                fastq_record = FASTX.FASTQ.Record("qual_string_$i", seq_str, qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                string_qual_bio_graph = Mycelia.build_qualmer_string_graph_next(quality_reads, ngram_length=5; graph_mode=Mycelia.SingleStrand)
                test_round_trip_reconstruction(string_qual_bio_graph, String)
            end
        end
    end
end
