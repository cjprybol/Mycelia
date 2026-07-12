import Test
import Mycelia

Test.@testset "taxid harmonization (merged.dmp)" begin
    mktempdir() do dir
        # Real NCBI merged.dmp format: old_tax_id \t | \t new_tax_id \t |
        dmp = joinpath(dir, "merged.dmp")
        write(dmp,
            "12\t|\t74109\t|\n" *
            "30\t|\t40\t|\n" *
            "\t|\t\n" *          # blank/garbage line must be skipped
            "not_an_int\t|\t5\t|\n")
        m = Mycelia.load_merged_taxid_map(dmp)
        Test.@test m[12] == 74109
        Test.@test m[30] == 40
        Test.@test length(m) == 2                        # unparseable lines skipped

        # Retired taxids remap; taxids absent from the map pass through unchanged.
        Test.@test Mycelia.harmonize_taxids([12, 999, 30], m) == [74109, 999, 40]
        Test.@test Mycelia.harmonize_taxids(Int[], m) == Int[]

        # String-path overload loads the map itself.
        Test.@test Mycelia.harmonize_taxids([12, 999], dmp) == [74109, 999]

        # Missing file errors rather than silently returning identity.
        Test.@test_throws ErrorException Mycelia.load_merged_taxid_map(
            joinpath(dir, "nope.dmp"))
    end
end
