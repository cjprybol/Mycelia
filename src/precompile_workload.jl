# PrecompileTools workload for faster startup
#
# This file defines a minimal workload that gets precompiled to reduce
# time-to-first-analysis. It covers the core "hello world" flow:
# - Load and parse a small FASTX file
# - Build a minimal k-mer graph
# - Basic graph operations
#
# Note: This workload should be side-effect free (no file I/O, no network)

import PrecompileTools

PrecompileTools.@setup_workload begin
    # Setup code - create minimal test data in memory

    # Small DNA sequence for testing
    test_dna = BioSequences.LongDNA{4}("ATGCATGCATGC")
    test_rna = BioSequences.LongRNA{4}("AUGCAUGCAUGC")
    test_aa = BioSequences.LongAA("MKHLLVGGG")

    # Quality scores
    test_quality = fill(UInt8(30), 12)

    PrecompileTools.@compile_workload begin
        # Core BioSequences operations
        rc = BioSequences.reverse_complement(test_dna)
        local canonical_seq = BioSequences.canonical(test_dna)

        # K-mer operations
        kmer = Kmers.DNAKmer{5}(BioSequences.LongDNA{4}("ATGCA"))
        kmer_rc = Kmers.reverse_complement(kmer)
        kmer_can = Kmers.canonical(kmer)

        # FASTA record creation
        fasta_rec = FASTX.FASTA.Record("test", test_dna)
        seq_from_rec = FASTX.sequence(BioSequences.LongDNA{4}, fasta_rec)

        # FASTQ record creation
        fastq_rec = FASTX.FASTQ.Record("test", test_dna, test_quality)
        seq_from_fastq = FASTX.sequence(BioSequences.LongDNA{4}, fastq_rec)
        qual = FASTX.quality(fastq_rec)

        # Graph creation (minimal)
        graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = String,
            edge_data_type = Int
        )
        graph["A"] = "vertex_A"
        graph["B"] = "vertex_B"
        graph["A", "B"] = 1

        # Graph queries
        labels = collect(MetaGraphsNext.labels(graph))
        nv = MetaGraphsNext.nv(graph)
        ne = MetaGraphsNext.ne(graph)
        has_v = haskey(graph, "A")

        # Alphabet detection
        alphabet = detect_alphabet(string(test_dna))

        # DataFrame creation (common operation)
        df = DataFrames.DataFrame(x = [1, 2, 3], y = ["a", "b", "c"])

        # String operations used in parsing
        s = "ATGC"
        parts = split(s, "")
        joined = join(parts, "")
    end
end
