import Test
import CSV
import DataFrames
import FASTX
import Mycelia

include(joinpath(@__DIR__, "..", "..", "benchmarking", "standard_assembler_fixtures.jl"))

Test.@testset "Standard assembler fixture catalog" begin
    fixture_table = list_standard_assembler_fixtures()

    Test.@test DataFrames.nrow(fixture_table) == 2
    Test.@test Set(fixture_table.name) == Set([
        "synthetic_isolate_5386",
        "synthetic_metagenome_pair"
    ])

    plan = build_standard_assembler_benchmark_plan()
    Test.@test DataFrames.nrow(plan) == 6
    Test.@test Set(plan.assembler) == Set(["Rhizomorph", "MEGAHIT", "metaSPAdes"])
end

Test.@testset "Synthetic isolate fixture materialization" begin
    mktempdir() do dir
        fixture = materialize_standard_assembler_fixture("synthetic_isolate_5386"; outdir = dir)

        Test.@test isfile(fixture.reference_fasta)
        Test.@test isfile(fixture.fastq1)
        Test.@test isfile(fixture.fastq2)
        Test.@test fixture.total_reference_bases == 5_386

        reference_records = collect(Mycelia.open_fastx(fixture.reference_fasta))
        Test.@test length(reference_records) == 1
        Test.@test length(FASTX.sequence(reference_records[1])) == 5_386

        forward_reads = collect(Mycelia.open_fastx(fixture.fastq1))
        reverse_reads = collect(Mycelia.open_fastx(fixture.fastq2))
        Test.@test length(forward_reads) == length(reverse_reads)
        Test.@test !isempty(forward_reads)

        truth_table = CSV.read(fixture.truth_table_path, DataFrames.DataFrame)
        Test.@test DataFrames.nrow(truth_table) == 1
        Test.@test truth_table.sequence_id[1] == "synthetic_isolate_5386"
    end
end

Test.@testset "Synthetic metagenome fixture dry materialization" begin
    mktempdir() do dir
        fixture = materialize_standard_assembler_fixture(
            "synthetic_metagenome_pair";
            outdir = dir,
            emit_reads = false
        )

        Test.@test isfile(fixture.reference_fasta)
        Test.@test fixture.fastq1 === nothing
        Test.@test fixture.fastq2 === nothing
        Test.@test fixture.total_reference_bases == 7_000
        Test.@test DataFrames.nrow(fixture.truth_table) == 2
    end
end
