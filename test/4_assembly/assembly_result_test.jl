import Test
import Mycelia
import FASTX

Test.@testset "Assembly Result Utilities" begin
    contigs = ["ATCG", "GG"]
    contig_names = ["contig_1", "contig_2"]

    graph = Mycelia.Rhizomorph.build_ngram_graph(["ATCG"], 2; dataset_id = "test")
    fastq_record = FASTX.FASTQ.Record("contig_1", "ATCG", "IIII")

    result = Mycelia.Rhizomorph.AssemblyResult(
        contigs,
        contig_names;
        graph = graph,
        assembly_stats = Dict("quality_preserved" => true),
        fastq_contigs = [fastq_record]
    )

    Test.@test Mycelia.Rhizomorph.get_fastq_contigs(result) == [fastq_record]
    Test.@test Mycelia.Rhizomorph.has_quality_information(result)

    Test.@testset "FASTQ contig output" begin
        mktempdir() do dir
            output_path = joinpath(dir, "contigs.fastq")
            written = Mycelia.Rhizomorph.write_fastq_contigs(result, output_path)
            Test.@test isfile(written)
            contents = read(written, String)
            Test.@test occursin("contig_1", contents)
        end
    end

    Test.@testset "Assembly serialization" begin
        mktempdir() do dir
            output_path = joinpath(dir, "assembly.jld2")
            Mycelia.Rhizomorph.save_assembly(result, output_path)
            loaded = Mycelia.Rhizomorph.load_assembly(output_path)
            Test.@test loaded.contigs == result.contigs
            Test.@test loaded.contig_names == result.contig_names
        end
    end

    Test.@testset "GFA export and validation" begin
        result_with_paths = Mycelia.Rhizomorph.AssemblyResult(
            contigs,
            contig_names;
            graph = graph,
            simplified_graph = graph,
            paths = Dict("contig_1" => ["AT", "TC", "CG"]),
            gfa_compatible = true
        )

        Test.@test Mycelia.Rhizomorph.has_graph_structure(result_with_paths)
        Test.@test Mycelia.Rhizomorph.has_simplified_graph(result_with_paths)
        Test.@test Mycelia.Rhizomorph.has_paths(result_with_paths)

        mktempdir() do dir
            gfa_path = joinpath(dir, "assembly.gfa")
            written = Mycelia.Rhizomorph.write_gfa(result_with_paths, gfa_path)
            Test.@test isfile(written)
            contents = read(written, String)
            Test.@test occursin("\tcontig_1\t", contents)
        end

        bad_result = Mycelia.Rhizomorph.AssemblyResult(["A"], ["contig_1", "contig_2"]; gfa_compatible = true)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.write_gfa(bad_result, "unused.gfa")

        report = Mycelia.Rhizomorph.validate_assembly_structure(bad_result)
        Test.@test report["valid"] == false
        Test.@test "Contigs and contig_names have different lengths" in report["issues"]
        Test.@test "Marked as GFA compatible but no graph structure" in report["issues"]
    end

    Test.@testset "Assembly metrics" begin
        metrics = Mycelia.Rhizomorph.validate_assembly(result)
        Test.@test metrics["num_contigs"] == 2
        Test.@test metrics["total_length"] == 6
        Test.@test metrics["N50"] == 4
        Test.@test metrics["L50"] == 1
        Test.@test metrics["N90"] == 2
        Test.@test metrics["L90"] == 2
    end
end
