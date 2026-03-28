import Test
import FASTX
import BioSequences
import Mycelia

Test.@testset "Read Quality Control" begin
    Test.@testset "gc_content_per_read — type safety regression (td-zyf99)" begin
        # Regression guard: the function must use FASTX.sequence(LongDNA{4}, record)
        # so BioSequences symbol comparisons work.  If it reverts to plain
        # FASTX.sequence(record) (returns String), all GC fractions would be 0.0.

        make_record(id, seq) = FASTX.FASTQ.Record(
            id,
            seq,
            collect(fill(UInt8(40), length(seq)))  # Q40 quality scores
        )

        all_gc = make_record("all_gc", "GCGCGCGC")   # 100% GC
        all_at = make_record("all_at", "ATATATATAT")  # 0% GC
        half_gc = make_record("half_gc", "ATATGCGC")   # 50% GC
        empty_r = make_record("empty", "")

        records = [all_gc, all_at, half_gc, empty_r]
        fracs = Mycelia.gc_content_per_read(records)

        Test.@test length(fracs) == 4
        Test.@test fracs[1] ≈ 1.0
        Test.@test fracs[2] ≈ 0.0
        Test.@test fracs[3] ≈ 0.5
        Test.@test fracs[4] == 0.0   # empty record → 0.0, not NaN
    end
end
