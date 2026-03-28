import Test
import Mycelia
import FASTX
import BioSequences

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

Test.@testset "Assembly Core Coverage" begin
    Test.@testset "Assembly command helpers" begin
        single_end = Mycelia.megahit_cmd(
            fastq1 = "reads.fastq",
            min_contig_len = 500,
            k_list = "21,41",
            threads = 8
        )
        Test.@test single_end.outdir == "reads_megahit"
        Test.@test single_end.contigs == joinpath("reads_megahit", "final.contigs.fa")
        Test.@test occursin("megahit -r \"reads.fastq\"", single_end.cmd)
        Test.@test !occursin(" -1 ", single_end.cmd)
        Test.@test occursin("--min-contig-len 500", single_end.cmd)
        Test.@test occursin("--k-list 21,41", single_end.cmd)
        Test.@test occursin("-t 8", single_end.cmd)

        paired_end = Mycelia.megahit_cmd(
            fastq1 = "sample_R1.fastq.gz",
            fastq2 = "sample_R2.fastq.gz",
            outdir = "paired_output"
        )
        Test.@test paired_end.outdir == "paired_output"
        Test.@test paired_end.fastg == joinpath("paired_output", "final.contigs.fastg")
        Test.@test paired_end.gfa == joinpath("paired_output", "final.contigs.fastg.gfa")
        Test.@test occursin("megahit -1 \"sample_R1.fastq.gz\" -2 \"sample_R2.fastq.gz\"", paired_end.cmd)
        Test.@test occursin("gfatools view", paired_end.cmd)

        success_result, success_runtime = Mycelia.run_assembler("ok", () -> 42)
        Test.@test success_result == 42
        Test.@test success_runtime isa Real
        Test.@test success_runtime >= 0

        failure_result, failure_runtime = Mycelia.run_assembler("fail", () -> error("boom"))
        Test.@test failure_result === nothing
        Test.@test ismissing(failure_runtime)
    end

    Test.@testset "Iterative assembly helper edge cases" begin
        short_reads = [
            Mycelia.fastq_record(identifier = "short1", sequence = "A", quality_scores = [30]),
            Mycelia.fastq_record(identifier = "short2", sequence = "", quality_scores = UInt8[])
        ]
        Test.@test Mycelia.calculate_sparsity(short_reads, 3) == 1.0
        Test.@test !Mycelia.errors_are_singletons(short_reads, 3)

        high_coverage_reads = [
            Mycelia.fastq_record(identifier = "r1", sequence = "AAAAAA", quality_scores = fill(30, 6)),
            Mycelia.fastq_record(identifier = "r2", sequence = "AAAAAA", quality_scores = fill(30, 6))
        ]
        Test.@test !Mycelia.errors_are_singletons(high_coverage_reads, 2)
        Test.@test Mycelia.find_initial_k(high_coverage_reads; k_range = [5], sparsity_threshold = 0.9) == 5
        Test.@test Mycelia.find_initial_k(high_coverage_reads; k_range = Int[]) == 3
        Test.@test Mycelia._label_k_length(:ATGC) == 4
    end

    Test.@testset "Rhizomorph assembly helper validation" begin
        Test.@test Mycelia.Rhizomorph._detect_sequence_type(FASTX.FASTA.Record[]) <: BioSequences.LongDNA

        rna_reads = [FASTX.FASTA.Record("rna1", "AUGC")]
        aa_reads = [FASTX.FASTA.Record("aa1", "MKWV")]
        Test.@test Mycelia.Rhizomorph._detect_sequence_type(rna_reads) <: BioSequences.LongRNA
        Test.@test Mycelia.Rhizomorph._detect_sequence_type(aa_reads) <: BioSequences.LongAA
        Test.@test_throws ErrorException Mycelia.Rhizomorph._detect_sequence_type([1, 2, 3])

        Test.@test Mycelia.Rhizomorph._determine_kmer_type(FASTX.FASTA.Record[], 4) ===
                   Mycelia.Kmers.DNAKmer{4}
        Test.@test Mycelia.Rhizomorph._determine_kmer_type(rna_reads, 3) ===
                   Mycelia.Kmers.RNAKmer{3}
        Test.@test Mycelia.Rhizomorph._determine_kmer_type(aa_reads, 2) ===
                   Mycelia.Kmers.AAKmer{2}
        Test.@test_throws ErrorException Mycelia.Rhizomorph._determine_kmer_type(["ATGC"], 3)

        default_config = Mycelia.Rhizomorph.AssemblyConfig()
        Test.@test default_config.k == 31
        Test.@test default_config.min_overlap === nothing

        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(k = 21, min_overlap = 10)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(
            sequence_type = BioSequences.LongAA,
            graph_mode = Mycelia.Rhizomorph.DoubleStrand
        )
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(
            sequence_type = String,
            graph_mode = Mycelia.Rhizomorph.DoubleStrand
        )
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(k = 0)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(min_overlap = 0)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(error_rate = 1.5)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(min_coverage = 0)

        mktempdir() do dir
            fastq_path = joinpath(dir, "reads.fastq")
            open(fastq_path, "w") do io
                println(io, "@read1")
                println(io, "ATGC")
                println(io, "+")
                println(io, "IIII")
            end

            prepared = Mycelia.Rhizomorph._prepare_observations([fastq_path])
            Test.@test length(prepared) == 1
            Test.@test prepared[1] isa FASTX.FASTA.Record
            Test.@test FASTX.FASTA.sequence(prepared[1]) == "ATGC"
        end
        Test.@test_throws ArgumentError Mycelia.Rhizomorph._prepare_observations([1, 2, 3])

        fastq_reads = [FASTX.FASTQ.Record("read1", "ATGC", "IIII")]
        auto_config = Mycelia.Rhizomorph._auto_configure_assembly(
            fastq_reads;
            min_overlap = 2,
            graph_mode = Mycelia.Rhizomorph.SingleStrand
        )
        Test.@test auto_config.k === nothing
        Test.@test auto_config.min_overlap == 2
        Test.@test auto_config.use_quality_scores
        Test.@test auto_config.graph_mode == Mycelia.Rhizomorph.SingleStrand
    end

    Test.@testset "Assembly result reporting helpers" begin
        ngram_graph = Mycelia.Rhizomorph.build_ngram_graph(["ATCG"], 2; dataset_id = "ngram")
        fasta_graph = Mycelia.Rhizomorph.build_fasta_graph(
            [FASTX.FASTA.Record("dna1", "ATCG"), FASTX.FASTA.Record("dna2", "TCGA")];
            dataset_id = "fasta",
            min_overlap = 2
        )
        fastq_record = FASTX.FASTQ.Record("contig_1", "ATCG", "IIII")

        qualityless_result = Mycelia.Rhizomorph.AssemblyResult(["ATCG"], ["contig_1"])
        Test.@test_throws ErrorException Mycelia.Rhizomorph.write_fastq_contigs(
            qualityless_result,
            "unused.fastq"
        )

        warning_result = Mycelia.Rhizomorph.AssemblyResult(
            ["ATCG"],
            ["contig_1"];
            graph = ngram_graph,
            simplified_graph = fasta_graph,
            gfa_compatible = false
        )
        warning_report = Mycelia.Rhizomorph.validate_assembly_structure(warning_result)
        Test.@test warning_report["valid"]
        Test.@test "Graph and simplified_graph have different types" in warning_report["warnings"]

        path_only_result = Mycelia.Rhizomorph.AssemblyResult(
            ["ATCG"],
            ["contig_1"];
            paths = Dict("contig_1" => ["AT", "TC"]),
            gfa_compatible = false
        )
        path_report = Mycelia.Rhizomorph.validate_assembly_structure(path_only_result)
        Test.@test path_report["valid"]
        Test.@test "Paths provided but no graph structure available" in path_report["warnings"]

        metrics = Mycelia.Rhizomorph.validate_assembly(
            Mycelia.Rhizomorph.AssemblyResult(["ATCG", "GG"], ["contig_1", "contig_2"]);
            reference = "ATCGAT"
        )
        Test.@test metrics["reference_provided"]
        Test.@test metrics["reference_length"] == 6

        mktempdir() do dir
            gfa_path = joinpath(dir, "assembly.gfa")
            write(gfa_path, "H\tVN:Z:1.0\n")
            paths = Dict{String, Vector}(
                "contig_1" => ["AT", "TC"],
                "contig_2" => String[],
                "ghost" => ["ZZ"]
            )
            Mycelia.Rhizomorph._append_gfa_paths(
                gfa_path,
                paths,
                ["ATCG", "GG"],
                ["contig_1", "contig_2"]
            )
            contents = read(gfa_path, String)
            Test.@test occursin("P\tcontig_1\tAT+,TC+\t*,*", contents)
            Test.@test !occursin("P\tcontig_2", contents)
            Test.@test !occursin("ghost", contents)
        end

        quality_result = Mycelia.Rhizomorph.AssemblyResult(
            ["ATCG"],
            ["contig_1"];
            graph = ngram_graph,
            assembly_stats = Dict("quality_preserved" => true),
            fastq_contigs = [fastq_record]
        )
        Test.@test Mycelia.Rhizomorph.has_quality_information(quality_result)
    end
end
