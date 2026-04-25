import Test
import Mycelia
import FASTX
import Kmers
import Graphs
import MetaGraphsNext

function make_fastq_record(id::AbstractString, sequence::AbstractString, phred_scores::Vector{Int})
    quality_string = String(Char.(phred_scores .+ 33))
    return FASTX.FASTQ.Record(id, sequence, quality_string)
end

function write_fasta(path::AbstractString, entries::Vector{Tuple{String, String}})
    open(path, "w") do io
        for (id, sequence) in entries
            println(io, ">$id")
            println(io, sequence)
        end
    end
    return path
end

function write_fastq(path::AbstractString, records::Vector{FASTX.FASTQ.Record})
    open(FASTX.FASTQ.Writer, path) do writer
        for record in records
            write(writer, record)
        end
    end
    return path
end

Test.@testset "Rhizomorph k-mer and qualmer wrapper coverage" begin
    fasta_records = [FASTX.FASTA.Record("dna_read", "ATGCAT")]
    fastq_records = [make_fastq_record("dna_read", "ATGCAT", fill(30, 6))]

    Test.@testset "build_kmer_graph handles wrapper branches" begin
        lightweight_graph = Mycelia.Rhizomorph.build_kmer_graph(
            fasta_records,
            3;
            lightweight = true
        )
        lightweight_label = first(MetaGraphsNext.labels(lightweight_graph))
        Test.@test lightweight_graph[lightweight_label] isa Mycelia.Rhizomorph.LightweightKmerVertexData

        full_graph = Mycelia.Rhizomorph.build_kmer_graph(
            fasta_records,
            3;
            memory_profile = :ultralight,
            lightweight = false
        )
        full_label = first(MetaGraphsNext.labels(full_graph))
        Test.@test full_graph[full_label] isa Mycelia.Rhizomorph.KmerVertexData

        for mode in (:singlestrand, :doublestrand, :canonical)
            graph = Mycelia.Rhizomorph.build_kmer_graph(
                fasta_records,
                3;
                memory_profile = :ultralight,
                mode = mode
            )
            Test.@test MetaGraphsNext.nv(graph) > 0
        end

        for mode in (:singlestrand, :doublestrand, :canonical)
            graph = Mycelia.Rhizomorph.build_kmer_graph(
                fastq_records,
                3;
                memory_profile = :ultralight_quality,
                mode = mode
            )
            Test.@test MetaGraphsNext.nv(graph) > 0
        end

        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            fasta_records,
            3;
            memory_profile = :not_a_profile
        )
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            fasta_records,
            3;
            mode = :not_a_mode
        )
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            fasta_records,
            3;
            memory_profile = :ultralight,
            mode = :not_a_mode
        )
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph(
            FASTX.FASTQ.Record[],
            3;
            memory_profile = :ultralight_quality
        )
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            fasta_records,
            3;
            memory_profile = :ultralight_quality
        )
    end

    Test.@testset "build_kmer_graph_from_file and _from_files handle edge cases" begin
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_from_file(
            "/tmp/definitely_missing_mycelia_reads.fasta",
            3
        )

        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_from_files(
            String[],
            3
        )

        mktempdir() do dir
            dna_path_1 = write_fasta(
                joinpath(dir, "sample_a.fna"),
                [("a1", "ATGCAT"), ("a2", "ATGCAA")]
            )
            dna_path_2 = write_fasta(
                joinpath(dir, "sample_b.fna"),
                [("b1", "ATGCAT"), ("b2", "ATGCAA")]
            )

            reduced_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path_1, dna_path_2],
                3;
                memory_profile = :ultralight
            )
            reduced_label = first(MetaGraphsNext.labels(reduced_graph))
            Test.@test Set(keys(reduced_graph[reduced_label].dataset_counts)) ==
                       Set(["sample_a", "sample_b"])

            lightweight_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path_1, dna_path_2],
                3;
                lightweight = true
            )
            lightweight_label = first(MetaGraphsNext.labels(lightweight_graph))
            Test.@test lightweight_graph[lightweight_label] isa Mycelia.Rhizomorph.LightweightKmerVertexData

            full_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path_1],
                3;
                lightweight = false
            )
            full_label = first(MetaGraphsNext.labels(full_graph))
            Test.@test full_graph[full_label] isa Mycelia.Rhizomorph.KmerVertexData

            Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path_1],
                3;
                mode = :not_a_mode
            )
            Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path_1],
                3;
                type_hint = :rna
            )

            aa_path_1 = write_fasta(joinpath(dir, "protein_a.faa"), [("p1", "ACDEFG")])
            aa_path_2 = write_fasta(joinpath(dir, "protein_b.faa"), [("p2", "ACDFGG")])
            conflicting_path = write_fasta(joinpath(dir, "sample_rna.frn"), [("r1", "AUGCAU")])
            Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [dna_path_1, conflicting_path],
                3
            )
            Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [aa_path_1, aa_path_2],
                3;
                mode = :canonical,
                memory_profile = :ultralight
            )
        end
    end

    Test.@testset "build_qualmer_graph handles wrapper branches" begin
        for profile in (:ultralight_quality, :lightweight_quality)
            singlestrand = Mycelia.Rhizomorph.build_qualmer_graph(
                fastq_records,
                3;
                memory_profile = profile
            )
            singlestrand_label = first(MetaGraphsNext.labels(singlestrand))
            Test.@test hasproperty(singlestrand[singlestrand_label], :joint_quality)

            doublestrand = Mycelia.Rhizomorph.build_qualmer_graph(
                fastq_records,
                3;
                memory_profile = profile,
                mode = :doublestrand
            )
            Test.@test Graphs.is_directed(doublestrand.graph)

            canonical = Mycelia.Rhizomorph.build_qualmer_graph(
                fastq_records,
                3;
                memory_profile = profile,
                mode = :canonical
            )
            Test.@test MetaGraphsNext.nv(canonical) > 0
        end

        full_doublestrand = Mycelia.Rhizomorph.build_qualmer_graph(
            fastq_records,
            3;
            mode = :doublestrand
        )
        Test.@test Graphs.is_directed(full_doublestrand.graph)

        full_canonical = Mycelia.Rhizomorph.build_qualmer_graph(
            fastq_records,
            3;
            mode = :canonical
        )
        Test.@test MetaGraphsNext.nv(full_canonical) > 0

        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_qualmer_graph(
            fastq_records,
            3;
            memory_profile = :lightweight
        )
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_qualmer_graph(
            fastq_records,
            3;
            mode = :not_a_mode
        )
    end
end
