# Default outbase for merge_and_map_single_end_samples must be unique per
# input set (td-s9ig). A date-only default made every same-day call in one
# CWD share paths, so a later run with different FASTQs silently reused the
# earlier outputs.
#
# Run:
#   julia --project=. -e 'include("test/7_comparative_pangenomics/merge_and_map_outbase_test.jl")'

import Test
import Mycelia

Test.@testset "merge_and_map_single_end_samples default outbase (td-s9ig)" begin
    mktempdir() do dir
        a = [joinpath(dir, "a.fq")]
        b = [joinpath(dir, "b.fq")]
        ab = [joinpath(dir, "a.fq"), joinpath(dir, "b.fq")]
        ba = [joinpath(dir, "b.fq"), joinpath(dir, "a.fq")]
        ref = joinpath(dir, "ref.fa")
        idx = joinpath(dir, "ref.mmi")
        other_ref = joinpath(dir, "other.fa")
        for path in (a[1], b[1], ref, idx, other_ref)
            write(path, "x")
        end

        out_a = Mycelia.merge_and_map_default_outbase(a;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont")
        out_b = Mycelia.merge_and_map_default_outbase(b;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont")
        out_a2 = Mycelia.merge_and_map_default_outbase(a;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont")
        out_ab = Mycelia.merge_and_map_default_outbase(ab;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont")
        out_ba = Mycelia.merge_and_map_default_outbase(ba;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont")
        out_a_other_ref = Mycelia.merge_and_map_default_outbase(a;
            fasta_reference = other_ref,
            minimap_index = idx, mapping_type = "map-ont")

        Test.@test out_a != out_b
        Test.@test out_a == out_a2
        Test.@test out_ab == out_ba
        Test.@test out_a != out_a_other_ref
        # Not the colliding date-only name.
        Test.@test !startswith(basename(out_a), Mycelia.normalized_current_date())
        # An explicit outbase is used as-is (resolution helper, no I/O).
        Test.@test Mycelia.resolve_merge_and_map_outbase("my-out", a;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont") == "my-out"
        Test.@test Mycelia.resolve_merge_and_map_outbase(nothing, a;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont") == out_a

        # An in-place FASTQ rewrite must not reuse the previous cache key.
        write(a[1], "xy")
        out_a_rewritten = Mycelia.merge_and_map_default_outbase(a;
            fasta_reference = ref, minimap_index = idx, mapping_type = "map-ont")
        Test.@test out_a_rewritten != out_a
    end
end
