import Test
import Mycelia
import FASTX

Test.@testset "Testing utilities" begin
    Test.@testset "create_test_fasta_records" begin
        records = Mycelia.create_test_fasta_records(3; length = 12, seed = 7, id_prefix = "fixture")

        Test.@test length(records) == 3
        Test.@test FASTX.identifier.(records) == ["fixture_1", "fixture_2", "fixture_3"]
        Test.@test all(length(FASTX.sequence(record)) == 12 for record in records)

        repeated = Mycelia.create_test_fasta_records(3; length = 12, seed = 7, id_prefix = "fixture")
        Test.@test FASTX.identifier.(repeated) == FASTX.identifier.(records)
        Test.@test FASTX.sequence.(repeated) == FASTX.sequence.(records)

        shifted_seed = Mycelia.create_test_fasta_records(3; length = 12, seed = 8, id_prefix = "fixture")
        Test.@test FASTX.sequence.(shifted_seed) != FASTX.sequence.(records)

        empty_records = Mycelia.create_test_fasta_records(0; length = 12)
        Test.@test isempty(empty_records)

        Test.@test_throws ArgumentError Mycelia.create_test_fasta_records(-1)
        Test.@test_throws ArgumentError Mycelia.create_test_fasta_records(1; length = -1)
    end
end
