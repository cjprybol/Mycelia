# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/gfa_io_next.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/gfa_io_next.jl", "test/4_assembly", execute=false)'
# ```

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers

Test.@testset "GFA I/O Next-Generation Tests (Rhizomorph)" begin
    function make_fastq_record(identifier, sequence)
        return FASTX.FASTQ.Record(
            identifier,
            sequence,
            fill(UInt8(35), length(sequence))
        )
    end

    function edge_label_set(graph)
        return Set(Tuple(edge_label) for edge_label in MetaGraphsNext.edge_labels(graph))
    end

    function parse_gfa_segments_and_links(path)
        segment_map = Dict{String, String}()
        links = Tuple{String, String}[]
        for line in eachline(path)
            fields = split(line, '\t')
            isempty(fields) && continue
            if fields[1] == "S" && length(fields) >= 3
                segment_map[fields[2]] = fields[3]
            elseif fields[1] == "L" && length(fields) >= 5
                push!(links, (fields[2], fields[4]))
            end
        end
        return segment_map, links
    end

    function assert_roundtrip(case_name, graph_builder, graph_reader, expected_label_type;
            segment_naming::Symbol = :numeric)
        original_graph = graph_builder()
        original_labels = Set(MetaGraphsNext.labels(original_graph))
        original_edges = edge_label_set(original_graph)

        Test.@test !isempty(original_labels)

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, replace(case_name, r"[^A-Za-z0-9]+" => "_") * ".gfa")
            written_file = Mycelia.Rhizomorph.write_gfa_next(
                original_graph,
                gfa_file;
                segment_naming = segment_naming
            )
            restored_graph = graph_reader(gfa_file)

            Test.@test written_file == gfa_file
            Test.@test isfile(gfa_file)
            Test.@test restored_graph isa MetaGraphsNext.MetaGraph
            Test.@test Graphs.is_directed(restored_graph.graph)
            Test.@test Set(MetaGraphsNext.labels(restored_graph)) == original_labels
            Test.@test edge_label_set(restored_graph) == original_edges
            Test.@test all(label isa expected_label_type for label in MetaGraphsNext.labels(restored_graph))
        end
    end

    dna_fasta_records = [
        FASTX.FASTA.Record("dna1", BioSequences.LongDNA{4}("ATGCAAA")),
        FASTX.FASTA.Record("dna2", BioSequences.LongDNA{4}("CAAATTT")),
    ]
    rna_fasta_records = [
        FASTX.FASTA.Record("rna1", BioSequences.LongRNA{4}("AUGCAAA")),
        FASTX.FASTA.Record("rna2", BioSequences.LongRNA{4}("CAAAUUU")),
    ]
    aa_fasta_records = [
        FASTX.FASTA.Record("aa1", BioSequences.LongAA("QWELFP")),
        FASTX.FASTA.Record("aa2", BioSequences.LongAA("LFPVYH")),
    ]

    dna_fastq_records = [
        make_fastq_record("dnaq1", BioSequences.LongDNA{4}("ATGCAAA")),
        make_fastq_record("dnaq2", BioSequences.LongDNA{4}("CAAATTT")),
    ]
    rna_fastq_records = [
        make_fastq_record("rnaq1", BioSequences.LongRNA{4}("AUGCAAA")),
        make_fastq_record("rnaq2", BioSequences.LongRNA{4}("CAAAUUU")),
    ]
    aa_fastq_records = [
        make_fastq_record("aaq1", BioSequences.LongAA("QWELFP")),
        make_fastq_record("aaq2", BioSequences.LongAA("LFPVYH")),
    ]

    string_ngram_inputs = ["abc123xyz", "123xyz789"]
    string_graph_inputs = ["alpha-123", "123-omega!"]
    quality_string_inputs = ["uvw456rst", "456rst!?"]
    quality_string_graph_inputs = ["left->mid", "mid->right"]

    Test.@testset "GFA Writing" begin
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            dna_fasta_records,
            3;
            dataset_id = "gfa_write",
            mode = :doublestrand
        )

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "test.gfa")
            result_file = Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)

            Test.@test result_file == gfa_file
            Test.@test isfile(gfa_file)

            content = read(gfa_file, String)
            lines = split(strip(content), '\n')

            Test.@test startswith(lines[1], "H\t")
            segment_lines = filter(line -> startswith(line, "S\t"), lines)
            Test.@test !isempty(segment_lines)
            Test.@test all(occursin("DP:f:", line) for line in segment_lines)

            link_lines = filter(line -> startswith(line, "L\t"), lines)
            Test.@test !isempty(link_lines)
            link_fields = split(first(link_lines), '\t')
            Test.@test link_fields[3] == "+"
            Test.@test link_fields[5] == "+"
            Test.@test endswith(link_fields[6], "M")
        end
    end

    Test.@testset "Sequence-Based Segment Naming For Bandage Labels" begin
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            dna_fasta_records,
            3;
            dataset_id = "gfa_sequence_names",
            mode = :singlestrand
        )

        mktempdir() do tmpdir
            numeric_gfa = joinpath(tmpdir, "numeric.gfa")
            sequence_gfa = joinpath(tmpdir, "sequence.gfa")

            Mycelia.Rhizomorph.write_gfa_next(graph, numeric_gfa)
            Mycelia.Rhizomorph.write_gfa_next(
                graph,
                sequence_gfa;
                segment_naming = :sequence
            )

            numeric_segments, numeric_links = parse_gfa_segments_and_links(numeric_gfa)
            sequence_segments, sequence_links = parse_gfa_segments_and_links(sequence_gfa)

            Test.@test !isempty(numeric_segments)
            Test.@test !isempty(sequence_segments)
            Test.@test all(key != value for (key, value) in numeric_segments)
            Test.@test all(key == value for (key, value) in sequence_segments)

            numeric_link_sequences = Set(
                (numeric_segments[src], numeric_segments[dst]) for (src, dst) in numeric_links
            )
            sequence_link_sequences = Set(
                (src, dst) for (src, dst) in sequence_links
            )

            Test.@test numeric_link_sequences == sequence_link_sequences
            Test.@test Set(values(numeric_segments)) == Set(keys(sequence_segments))
        end
    end

    Test.@testset "Sequence-Based Segment Naming Round-Trip" begin
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            dna_fasta_records,
            3;
            dataset_id = "gfa_sequence_roundtrip",
            mode = :singlestrand
        )
        expected_segments = Set(string.(MetaGraphsNext.labels(graph)))

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "sequence_ids_roundtrip.gfa")
            Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file; segment_naming = :sequence)

            sequence_segments, _ = parse_gfa_segments_and_links(gfa_file)
            Test.@test Set(keys(sequence_segments)) == expected_segments
            Test.@test Set(values(sequence_segments)) == expected_segments
        end

        assert_roundtrip(
            "DNA K-mer SingleStrand Sequence IDs",
            () -> Mycelia.Rhizomorph.build_kmer_graph(
                dna_fasta_records,
                3;
                dataset_id = "gfa_sequence_roundtrip",
                mode = :singlestrand
            ),
            path -> Mycelia.Rhizomorph.read_gfa_next(
                path,
                Kmers.DNAKmer{3},
                Mycelia.Rhizomorph.SingleStrand
            ),
            Kmers.DNAKmer{3};
            segment_naming = :sequence
        )
    end

    Test.@testset "GFA Reading" begin
        kmer_gfa = """
        H\tVN:Z:1.0
        S\t1\tATC\tDP:f:2.0
        S\t2\tTCG\tDP:f:1.0
        L\t1\t+\t2\t+\t2M
        """
        string_gfa = """
        H\tVN:Z:1.0
        S\t1\tab!\tDP:f:1.0
        S\t2\tb!?\tDP:f:1.0
        L\t1\t+\t2\t+\t2M
        """

        mktempdir() do tmpdir
            kmer_file = joinpath(tmpdir, "test_input.gfa")
            write(kmer_file, kmer_gfa)
            graph = Mycelia.Rhizomorph.read_gfa_next(
                kmer_file,
                Kmers.DNAKmer{3},
                Mycelia.Rhizomorph.DoubleStrand
            )

            Test.@test graph isa MetaGraphsNext.MetaGraph
            Test.@test Set(MetaGraphsNext.labels(graph)) ==
                       Set([Kmers.DNAKmer{3}("ATC"), Kmers.DNAKmer{3}("TCG")])
            Test.@test graph[Kmers.DNAKmer{3}("ATC"), Kmers.DNAKmer{3}("TCG")] isa
                       Mycelia.Rhizomorph.KmerEdgeData

            string_file = joinpath(tmpdir, "string_input.gfa")
            write(string_file, string_gfa)
            string_graph = Mycelia.Rhizomorph.read_gfa_next(string_file)
            forced_string_graph = Mycelia.Rhizomorph.read_gfa_next(
                string_file;
                force_biosequence_graph = true
            )

            for graph in (string_graph, forced_string_graph)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test Set(MetaGraphsNext.labels(graph)) == Set(["ab!", "b!?"])
                Test.@test graph["ab!", "b!?"] isa Mycelia.Rhizomorph.StringEdgeData
            end
        end
    end

    Test.@testset "GFA Read Then Collapse Linear Chains" begin
        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "chain.gfa")
            graph = Mycelia.Rhizomorph.build_ngram_graph(["ab!def"], 4; dataset_id = "chain")
            Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)

            roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_file)
            Test.@test roundtrip isa MetaGraphsNext.MetaGraph

            collapsed = deepcopy(roundtrip)
            Mycelia.Rhizomorph.collapse_linear_chains!(collapsed)

            labels = Set(string(label) for label in MetaGraphsNext.labels(collapsed))
            Test.@test labels == Set(["ab!def"])
            Test.@test Graphs.nv(collapsed.graph) == 1
            Test.@test Graphs.ne(collapsed.graph) == 0
        end

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "unicode_chain.gfa")
            graph = Mycelia.Rhizomorph.build_ngram_graph(["élan"], 2; dataset_id = "unicode_chain")
            Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)

            roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_file)
            collapsed = deepcopy(roundtrip)
            Mycelia.Rhizomorph.collapse_linear_chains!(collapsed)

            labels = Set(string(label) for label in MetaGraphsNext.labels(collapsed))
            Test.@test labels == Set(["élan"])
            Test.@test Graphs.nv(collapsed.graph) == 1
        end

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "cycle.gfa")
            graph = Mycelia.Rhizomorph.build_ngram_graph(["!x!x"], 2; dataset_id = "cycle")
            Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)

            roundtrip = Mycelia.Rhizomorph.read_gfa_next(gfa_file; force_biosequence_graph = true)
            before_labels = Set(string(label) for label in MetaGraphsNext.labels(roundtrip))
            collapsed = deepcopy(roundtrip)
            Mycelia.Rhizomorph.collapse_linear_chains!(collapsed)
            after_labels = Set(string(label) for label in MetaGraphsNext.labels(collapsed))

            Test.@test before_labels == after_labels
            Test.@test Graphs.nv(collapsed.graph) == Graphs.nv(roundtrip.graph)
        end
    end

    Test.@testset "GFA Round-trip Across 24 Graph Types" begin
        cases = [
            (
                name = "DNA K-mer SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_kmer_graph(
                    dna_fasta_records, 5; dataset_id = "dna_kmer_ss", mode = :singlestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.DNAKmer{5}, Mycelia.Rhizomorph.SingleStrand),
                label_type = Kmers.DNAKmer{5},
            ),
            (
                name = "DNA K-mer DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.build_kmer_graph(
                    dna_fasta_records, 5; dataset_id = "dna_kmer_ds", mode = :doublestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.DNAKmer{5}, Mycelia.Rhizomorph.DoubleStrand),
                label_type = Kmers.DNAKmer{5},
            ),
            (
                name = "DNA BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_fasta_graph(
                    dna_fasta_records; dataset_id = "dna_fasta_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongDNA{4},
            ),
            (
                name = "DNA BioSequence DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(
                    Mycelia.Rhizomorph.build_fasta_graph(
                        dna_fasta_records; dataset_id = "dna_fasta_ds", min_overlap = 3)),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongDNA{4},
            ),
            (
                name = "DNA Quality K-mer SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_qualmer_graph(
                    dna_fastq_records, 5; dataset_id = "dna_qual_ss", mode = :singlestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.DNAKmer{5}, Mycelia.Rhizomorph.SingleStrand),
                label_type = Kmers.DNAKmer{5},
            ),
            (
                name = "DNA Quality K-mer DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.build_qualmer_graph(
                    dna_fastq_records, 5; dataset_id = "dna_qual_ds", mode = :doublestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.DNAKmer{5}, Mycelia.Rhizomorph.DoubleStrand),
                label_type = Kmers.DNAKmer{5},
            ),
            (
                name = "DNA Quality BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_fastq_graph(
                    dna_fastq_records; dataset_id = "dna_fastq_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongDNA{4},
            ),
            (
                name = "DNA Quality BioSequence DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(
                    Mycelia.Rhizomorph.build_fastq_graph(
                        dna_fastq_records; dataset_id = "dna_fastq_ds", min_overlap = 3)),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongDNA{4},
            ),
            (
                name = "RNA K-mer SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_kmer_graph(
                    rna_fasta_records, 5; dataset_id = "rna_kmer_ss", mode = :singlestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.RNAKmer{5}, Mycelia.Rhizomorph.SingleStrand),
                label_type = Kmers.RNAKmer{5},
            ),
            (
                name = "RNA K-mer DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.build_kmer_graph(
                    rna_fasta_records, 5; dataset_id = "rna_kmer_ds", mode = :doublestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.RNAKmer{5}, Mycelia.Rhizomorph.DoubleStrand),
                label_type = Kmers.RNAKmer{5},
            ),
            (
                name = "RNA BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_fasta_graph(
                    rna_fasta_records; dataset_id = "rna_fasta_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongRNA{4},
            ),
            (
                name = "RNA BioSequence DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(
                    Mycelia.Rhizomorph.build_fasta_graph(
                        rna_fasta_records; dataset_id = "rna_fasta_ds", min_overlap = 3)),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongRNA{4},
            ),
            (
                name = "RNA Quality K-mer SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_qualmer_graph(
                    rna_fastq_records, 5; dataset_id = "rna_qual_ss", mode = :singlestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.RNAKmer{5}, Mycelia.Rhizomorph.SingleStrand),
                label_type = Kmers.RNAKmer{5},
            ),
            (
                name = "RNA Quality K-mer DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.build_qualmer_graph(
                    rna_fastq_records, 5; dataset_id = "rna_qual_ds", mode = :doublestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.RNAKmer{5}, Mycelia.Rhizomorph.DoubleStrand),
                label_type = Kmers.RNAKmer{5},
            ),
            (
                name = "RNA Quality BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_fastq_graph(
                    rna_fastq_records; dataset_id = "rna_fastq_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongRNA{4},
            ),
            (
                name = "RNA Quality BioSequence DoubleStrand",
                builder = () -> Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(
                    Mycelia.Rhizomorph.build_fastq_graph(
                        rna_fastq_records; dataset_id = "rna_fastq_ds", min_overlap = 3)),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongRNA{4},
            ),
            (
                name = "AA K-mer SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_kmer_graph(
                    aa_fasta_records, 3; dataset_id = "aa_kmer_ss", mode = :singlestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.AAKmer{3}, Mycelia.Rhizomorph.SingleStrand),
                label_type = Kmers.AAKmer{3},
            ),
            (
                name = "AA BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_fasta_graph(
                    aa_fasta_records; dataset_id = "aa_fasta_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongAA,
            ),
            (
                name = "AA Quality K-mer SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_qualmer_graph(
                    aa_fastq_records, 3; dataset_id = "aa_qual_ss", mode = :singlestrand),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, Kmers.AAKmer{3}, Mycelia.Rhizomorph.SingleStrand),
                label_type = Kmers.AAKmer{3},
            ),
            (
                name = "AA Quality BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_fastq_graph(
                    aa_fastq_records; dataset_id = "aa_fastq_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path; force_biosequence_graph = true),
                label_type = BioSequences.LongAA,
            ),
            (
                name = "String N-gram SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_ngram_graph(
                    string_ngram_inputs, 3; dataset_id = "string_ngram_ss"),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, String, Mycelia.Rhizomorph.SingleStrand),
                label_type = String,
            ),
            (
                name = "String BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_string_graph(
                    string_graph_inputs; dataset_id = "string_graph_ss", min_overlap = 3),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(path, String),
                label_type = String,
            ),
            (
                name = "String Quality N-gram SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_ngram_graph(
                    quality_string_inputs, 3; dataset_id = "string_ngram_quality"),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(
                    path, String, Mycelia.Rhizomorph.SingleStrand),
                label_type = String,
            ),
            (
                name = "String Quality BioSequence SingleStrand",
                builder = () -> Mycelia.Rhizomorph.build_string_graph(
                    quality_string_graph_inputs;
                    dataset_id = "string_graph_quality",
                    min_overlap = 3
                ),
                reader = path -> Mycelia.Rhizomorph.read_gfa_next(path, String),
                label_type = String,
            ),
        ]

        Test.@test length(cases) == 24

        for case in cases
            Test.@testset "$(case.name)" begin
                assert_roundtrip(case.name, case.builder, case.reader, case.label_type)
            end
        end
    end

    Test.@testset "Error Handling" begin
        Test.@test_throws SystemError Mycelia.Rhizomorph.read_gfa_next(
            "nonexistent.gfa", Kmers.DNAKmer{3}, Mycelia.Rhizomorph.DoubleStrand)

        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph();
            label_type = Kmers.DNAKmer{3},
            vertex_data_type = Mycelia.Rhizomorph.KmerVertexData{Kmers.DNAKmer{3}},
            edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
        )

        mktempdir() do tmpdir
            empty_gfa = joinpath(tmpdir, "empty.gfa")
            result = Mycelia.Rhizomorph.write_gfa_next(empty_graph, empty_gfa)
            Test.@test result == empty_gfa
            Test.@test isfile(empty_gfa)

            content = read(empty_gfa, String)
            lines = filter(!isempty, split(strip(content), '\n'))
            Test.@test length(lines) == 1
            Test.@test startswith(lines[1], "H\t")
        end
    end
end
