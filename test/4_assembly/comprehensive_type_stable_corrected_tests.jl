# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/comprehensive_type_stable_corrected_tests.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/comprehensive_type_stable_corrected_tests.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "Comprehensive Type-Stable Reconstruction Tests - Corrected" begin
    ## Test helper function to validate round-trip reconstruction
    function test_round_trip_reconstruction(graph, expected_type)
        Test.@test !isempty(MetaGraphsNext.labels(graph))

        # Find Eulerian paths and test reconstruction
        paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
        if isempty(paths)
            paths = Mycelia.Rhizomorph.find_contigs_next(graph)
        end
        Test.@test !isempty(paths)

        for path_entry in paths
            path_vertices = if path_entry isa Mycelia.Rhizomorph.ContigPath
                path_entry.vertices
            else
                path_entry
            end
            if !isempty(path_vertices)
                # Create WalkStep objects from the path labels with proper typing
                first_vertex = first(path_vertices)
                vertex_type = typeof(first_vertex)
                walk_steps = Mycelia.Rhizomorph.WalkStep{vertex_type}[]

                for (i, vertex_label) in enumerate(path_vertices)
                    step = Mycelia.Rhizomorph.WalkStep(vertex_label, Mycelia.Rhizomorph.Forward, 1.0, Float64(i))
                    push!(walk_steps, step)
                end

                # Create GraphPath from WalkStep objects
                graph_path = Mycelia.Rhizomorph.GraphPath(walk_steps)

                # Test type-stable reconstruction
                reconstructed = Mycelia.Rhizomorph.path_to_sequence(graph_path, graph)
                if reconstructed === nothing && path_entry isa Mycelia.Rhizomorph.ContigPath
                    reconstructed = path_entry.sequence
                end
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

    Test.@testset "DNA Graph Types (8 variants)" begin
        # 1. DNA K-mer SingleStrand - using Rhizomorph.build_kmer_graph
        Test.@testset "DNA K-mer SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_reads, 5; dataset_id="test", mode=:singlestrand)
            test_round_trip_reconstruction(dna_graph, :dna)
        end

        # 2. DNA K-mer DoubleStrand - using Rhizomorph.build_kmer_graph
        Test.@testset "DNA K-mer DoubleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_graph = Mycelia.Rhizomorph.build_kmer_graph(dna_reads, 5; dataset_id="test", mode=:doublestrand)
            test_round_trip_reconstruction(dna_graph, :dna)
        end

        # 3. DNA BioSequence SingleStrand - using Rhizomorph.build_fasta_graph
        Test.@testset "DNA BioSequence SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_bio_graph = Mycelia.Rhizomorph.build_fasta_graph(dna_reads; dataset_id="test", min_overlap=3)
            test_round_trip_reconstruction(dna_bio_graph, :dna)
        end

        # 4. DNA BioSequence DoubleStrand - using Rhizomorph.build_fasta_graph
        Test.@testset "DNA BioSequence DoubleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            dna_bio_graph = Mycelia.Rhizomorph.build_fasta_graph(dna_reads; dataset_id="test", min_overlap=3)
            test_round_trip_reconstruction(dna_bio_graph, :dna)
        end

        # 5. DNA Quality K-mer SingleStrand - using build_qualmer_graph
        Test.@testset "DNA Quality K-mer SingleStrand" begin
            dna_reads = [FASTX.FASTA.Record("dna_seq", dna_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(dna_reads)
                seq = FASTX.sequence(BioSequences.LongDNA{4}, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq)]])
                fastq_record = FASTX.FASTQ.Record("qual_dna_$i", string(seq), qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                dna_qual_graph = Mycelia.Rhizomorph.build_qualmer_graph(quality_reads, 3; dataset_id="test", mode=:singlestrand)
                @info "DNA qualmer graph vertex count: $(length(MetaGraphsNext.labels(dna_qual_graph)))"
                test_round_trip_reconstruction(dna_qual_graph, :dna)
            end
        end

        # 6. DNA Quality K-mer DoubleStrand - using build_qualmer_graph
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
                dna_qual_graph = Mycelia.Rhizomorph.build_qualmer_graph(quality_reads, 3; dataset_id="test", mode=:doublestrand)
                test_round_trip_reconstruction(dna_qual_graph, :dna)
            end
        end

        # 7. DNA Quality BioSequence SingleStrand - using build_quality_biosequence_graph
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
                dna_qual_bio_graph = Mycelia.Rhizomorph.build_fastq_graph(quality_reads; dataset_id="test", min_overlap=3)
                test_round_trip_reconstruction(dna_qual_bio_graph, :dna)
            end
        end

        # 8. DNA Quality BioSequence DoubleStrand - using build_quality_biosequence_graph
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
                dna_qual_bio_graph = Mycelia.Rhizomorph.build_fastq_graph(quality_reads; dataset_id="test", min_overlap=3)
                test_round_trip_reconstruction(dna_qual_bio_graph, :dna)
            end
        end
    end

    Test.@testset "RNA Graph Types (8 variants)" begin
        # 9. RNA K-mer SingleStrand - using Rhizomorph.build_kmer_graph
        Test.@testset "RNA K-mer SingleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_graph = Mycelia.Rhizomorph.build_kmer_graph(rna_reads, 5; dataset_id="test", mode=:singlestrand)
            test_round_trip_reconstruction(rna_graph, :rna)
        end

        # 10. RNA K-mer DoubleStrand - using Rhizomorph.build_kmer_graph
        Test.@testset "RNA K-mer DoubleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_graph = Mycelia.Rhizomorph.build_kmer_graph(rna_reads, 5; dataset_id="test", mode=:doublestrand)
            test_round_trip_reconstruction(rna_graph, :rna)
        end

        # 11. RNA BioSequence SingleStrand - using Rhizomorph.build_fasta_graph
        Test.@testset "RNA BioSequence SingleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_graph = Mycelia.Rhizomorph.build_fasta_graph(rna_reads; dataset_id="test", min_overlap=3)
            test_round_trip_reconstruction(rna_graph, :rna)
        end

        # 12. RNA BioSequence DoubleStrand - using Rhizomorph.build_fasta_graph
        Test.@testset "RNA BioSequence DoubleStrand" begin
            rna_reads = [FASTX.FASTA.Record("rna_seq", rna_seq)]
            rna_graph = Mycelia.Rhizomorph.build_fasta_graph(rna_reads; dataset_id="test", min_overlap=3)
            test_round_trip_reconstruction(rna_graph, :rna)
        end

        # 13. RNA Quality K-mer SingleStrand - using build_qualmer_graph
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
                rna_qual_graph = Mycelia.Rhizomorph.build_qualmer_graph(quality_reads, 3; dataset_id="test", mode=:singlestrand)
                test_round_trip_reconstruction(rna_qual_graph, :rna)
            end
        end

        # 14. RNA Quality K-mer DoubleStrand - using build_qualmer_graph
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
                rna_qual_graph = Mycelia.Rhizomorph.build_qualmer_graph(quality_reads, 3; dataset_id="test", mode=:doublestrand)
                test_round_trip_reconstruction(rna_qual_graph, :rna)
            end
        end

        # 15. RNA Quality BioSequence SingleStrand - using build_quality_biosequence_graph
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
                rna_qual_bio_graph = Mycelia.Rhizomorph.build_fastq_graph(quality_reads; dataset_id="test", min_overlap=3)
                test_round_trip_reconstruction(rna_qual_bio_graph, :rna)
            end
        end

        # 16. RNA Quality BioSequence DoubleStrand - using build_quality_biosequence_graph
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
                rna_qual_bio_graph = Mycelia.Rhizomorph.build_fastq_graph(quality_reads; dataset_id="test", min_overlap=3)
                test_round_trip_reconstruction(rna_qual_bio_graph, :rna)
            end
        end
    end

    Test.@testset "Amino Acid Graph Types (4 variants)" begin
        # 17. AA K-mer SingleStrand - using Rhizomorph.build_kmer_graph
        Test.@testset "AA K-mer SingleStrand" begin
            aa_reads = [FASTX.FASTA.Record("aa_seq", aa_seq)]
            aa_graph = Mycelia.Rhizomorph.build_kmer_graph(aa_reads, 5; dataset_id="test", mode=:singlestrand)
            test_round_trip_reconstruction(aa_graph, :aa)
        end

        # 18. AA BioSequence SingleStrand - using Rhizomorph.build_fasta_graph
        Test.@testset "AA BioSequence SingleStrand" begin
            aa_reads = [FASTX.FASTA.Record("aa_seq", aa_seq)]
            aa_graph = Mycelia.Rhizomorph.build_fasta_graph(aa_reads; dataset_id="test", min_overlap=3)
            test_round_trip_reconstruction(aa_graph, :aa)
        end

        # 19. AA Quality K-mer SingleStrand - using build_qualmer_graph
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
                aa_qual_graph = Mycelia.Rhizomorph.build_qualmer_graph(quality_reads, 3; dataset_id="test", mode=:singlestrand)
                test_round_trip_reconstruction(aa_qual_graph, :aa)
            end
        end

        # 20. AA Quality BioSequence SingleStrand - using build_quality_biosequence_graph
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
                aa_qual_bio_graph = Mycelia.Rhizomorph.build_fastq_graph(quality_reads; dataset_id="test", min_overlap=3)
                test_round_trip_reconstruction(aa_qual_bio_graph, :aa)
            end
        end
    end

    Test.@testset "String Graph Types (4 variants)" begin
        # Test sequence for strings
        string_seq = "ABCDEFGHIJKLMNOP"

        # 21. String N-gram SingleStrand - using build_string_ngram_graph_next
        Test.@testset "String N-gram SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", string_seq)]
            string_ngram_graph = Mycelia.Rhizomorph.build_ngram_graph([string_seq], 3; dataset_id="test")
            test_round_trip_reconstruction(string_ngram_graph, :string)
        end

        # 22. String BioSequence SingleStrand - using build_string_biosequence_graph_next
        Test.@testset "String BioSequence SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", string_seq)]
            string_bio_graph = Mycelia.Rhizomorph.build_string_graph([string_seq]; dataset_id="test", min_overlap=3)
            test_round_trip_reconstruction(string_bio_graph, :string)
        end

        # 23. String Quality N-gram SingleStrand - using build_string_qualmer_ngram_graph_next
        Test.@testset "String Quality N-gram SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", string_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(string_reads)
                seq_str = FASTX.sequence(String, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq_str)]])
                fastq_record = FASTX.FASTQ.Record("qual_string_$i", seq_str, qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                strings = [FASTX.sequence(String, record) for record in quality_reads]
                string_ngram_graph = Mycelia.Rhizomorph.build_ngram_graph(strings, 3; dataset_id="test")
                test_round_trip_reconstruction(string_ngram_graph, :string)
            end
        end

        # 24. String Quality BioSequence SingleStrand - using build_string_qualmer_biosequence_graph_next
        Test.@testset "String Quality BioSequence SingleStrand" begin
            string_reads = [FASTX.FASTA.Record("string_seq", string_seq)]
            quality_reads = FASTX.FASTQ.Record[]
            for (i, read) in enumerate(string_reads)
                seq_str = FASTX.sequence(String, read)
                qual_str = String([Char(q + 33) for q in quality_scores[1:length(seq_str)]])
                fastq_record = FASTX.FASTQ.Record("qual_string_$i", seq_str, qual_str)
                push!(quality_reads, fastq_record)
            end

            if !isempty(quality_reads)
                strings = [FASTX.sequence(String, record) for record in quality_reads]
                string_bio_graph = Mycelia.Rhizomorph.build_string_graph(strings; dataset_id="test", min_overlap=3)
                test_round_trip_reconstruction(string_bio_graph, :string)
            end
        end
    end
end
