import Test
import Mycelia

Test.@testset "taxid harmonization (merged.dmp)" begin
    mktempdir() do dir
        # Real NCBI merged.dmp format: old_tax_id \t | \t new_tax_id \t |
        dmp = joinpath(dir, "merged.dmp")
        write(dmp,
            "12\t|\t74109\t|\n" *
            "30\t|\t40\t|\n" *
            "\t|\t\n" *              # empty fields: caught by tryparse-nothing skip
            "not_an_int\t|\t5\t|\n" *  # non-integer old id: tryparse-nothing skip
            "garbage_no_pipe\n")       # no '|': caught by the length<2 skip
        m = Mycelia.load_merged_taxid_map(dmp)
        Test.@test m[12] == 74109
        Test.@test m[30] == 40
        Test.@test length(m) == 2                        # all 3 junk lines skipped

        # Retired taxids remap; taxids absent from the map pass through unchanged.
        Test.@test Mycelia.harmonize_taxids([12, 999, 30], m) == [74109, 999, 40]
        Test.@test Mycelia.harmonize_taxids(Int[], m) == Int[]
        # Non-Int64 element type: the Integer signature + Int(t) conversion contract.
        Test.@test Mycelia.harmonize_taxids(Int32[12, 999], m) == [74109, 999]
        Test.@test eltype(Mycelia.harmonize_taxids(Int32[12], m)) == Int

        # String-path overload loads the map itself.
        Test.@test Mycelia.harmonize_taxids([12, 999], dmp) == [74109, 999]

        # Missing file errors rather than silently returning identity.
        Test.@test_throws ErrorException Mycelia.load_merged_taxid_map(
            joinpath(dir, "nope.dmp"))

        # I1: a present-but-wrong input (no parseable mappings) must ERROR, not
        # return a silently empty map that degrades harmonize to identity.
        empty_dmp = joinpath(dir, "empty.dmp")
        write(empty_dmp, "garbage_no_pipe\nalso garbage\n")
        Test.@test_throws ErrorException Mycelia.load_merged_taxid_map(empty_dmp)
    end
end
