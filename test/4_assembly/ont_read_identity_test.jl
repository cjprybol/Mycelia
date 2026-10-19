# Contract tests for the ONT read-identity harness's pure helpers (td-4e19d.29).
#
# These pin two things that fail SILENTLY — no exception, no warning, just a
# number that is wrong in a consistent direction:
#
#   * the read-base accounting behind the k-mer coverage estimate. Identity is
#     measured from alignments, so it exists only for reads that aligned; raw
#     coverage is charged for every read. An estimate built from the mapped
#     subset alone is optimistic by exactly the unaligned base fraction, and
#     nothing about its output looks wrong.
#   * the SAM cache key. It decides whether an existing SAM is reused, and a key
#     that ignores which reads were aligned turns a cache hit into a report
#     about a different dataset.
#
# Run:
#   julia --project=. test/4_assembly/ont_read_identity_test.jl

# Wrapped in a module: runtests.jl includes every test file into one shared
# `Main`, and this script defines top-level consts (COVERAGE, SEED, ORGANISM,
# K_LADDER, ...) and an `arg_value` helper that collide with the other
# benchmarking harnesses already included there.
module OntReadIdentityTest

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "ont_read_identity.jl"))

"""SAM record with the 11 mandatory fields plus any `tags`."""
function sam_record(qname, flag, cigar, seq; tags = String[])
    fields = [qname, string(flag), "ref", "1", "60", cigar, "*", "0", "0",
        seq, "*"]
    append!(fields, tags)
    return join(fields, '\t')
end

Test.@testset "ONT read-identity helpers" begin
    Test.@testset "unaligned reads enter the base denominator, not the identity set" begin
        # The population this script describes is EVERY sequenced read. Badread
        # emits ~1% junk and ~1% random reads by default, and those do not align;
        # they still consumed coverage and still contribute no k-mer matching the
        # reference. Counting them as absent rather than as zero is what made the
        # error-free k-mer coverage optimistic.
        mktempdir() do dir
            sam = joinpath(dir, "test.sam")
            write(sam,
                join(
                    [
                        "@SQ\tSN:ref\tLN:1000",
                        # Mapped and scorable: 10 aligned columns, 1 mismatch.
                        sam_record("mapped", 0, "10M", "A"^10; tags = ["NM:i:1"]),
                        # Unmapped: no identity, but 20 sequenced bases.
                        sam_record("unmapped", 4, "*", "C"^20),
                        # Secondary: excluded entirely, or a chimeric read would
                        # be counted twice.
                        sam_record("mapped", 256, "5M", "A"^5; tags = ["NM:i:0"]),
                        # Primary, mapped, but no NM tag — unscorable. Its bases
                        # were still sequenced.
                        sam_record("no_nm", 0, "10M", "G"^30)
                    ],
                    '\n') * '\n')

            result = identities_from_sam(sam)

            Test.@test result.n_primary == 1        # only `mapped` scored
            Test.@test result.n_unmapped == 1
            Test.@test result.n_skipped == 1        # `no_nm`
            Test.@test result.n_records == 3        # secondary excluded

            # The denominator is every primary record's read length; the
            # numerator only the records that produced an identity.
            Test.@test result.total_read_bases == 10 + 20 + 30
            Test.@test result.scored_read_bases == 10

            # THE assertion this file exists for. 10/60, not 10/10: a helper that
            # summed only mapped bases would report 1.0 here and inflate every
            # downstream error-free coverage figure by 6x.
            Test.@test result.scored_read_bases / result.total_read_bases ≈ 10 / 60

            # Identity itself is unchanged — this is an accounting fix, not a
            # re-definition of e. 9 matches over 10 aligned columns.
            mapped_rows = [r for r in result.rows if r.mapped]
            Test.@test length(mapped_rows) == 1
            Test.@test mapped_rows[1].blast_identity ≈ 0.9
            Test.@test mapped_rows[1].read_length == 10

            # The unmapped read is carried as a row with its length, so
            # per_read_identity.tsv records what it contributed.
            unmapped_rows = [r for r in result.rows if !r.mapped]
            Test.@test length(unmapped_rows) == 1
            Test.@test unmapped_rows[1].read_length == 20
        end
    end

    Test.@testset "a missing SEQ is counted, because it biases the fraction" begin
        # SEQ = "*" makes a read contribute zero to the denominator, which pushes
        # the aligned fraction UP. The caller warns on this rather than silently
        # publishing the optimistic number, so the count has to be tracked.
        mktempdir() do dir
            sam = joinpath(dir, "nosq.sam")
            write(sam,
                sam_record("mapped", 0, "10M", "A"^10; tags = ["NM:i:0"]) * "\n" *
                sam_record("unmapped", 4, "*", "*") * "\n")
            result = identities_from_sam(sam)
            Test.@test result.n_missing_seq == 1
            Test.@test result.total_read_bases == 10
        end
    end

    Test.@testset "the SAM cache tag identifies the read FILE" begin
        # --coverage/--seed do not determine the reads when --reads overrides
        # them, so the alignment cache is keyed to the file itself.
        mktempdir() do dir
            a = joinpath(dir, "a.fq")
            b = joinpath(dir, "b.fq")
            write(a, "@r1\nACGT\n+\nIIII\n")
            write(b, "@r1\nACGT\n+\nIIII\n")

            # Different files are different cache entries even with identical
            # CONTENT: two runs writing distinct paths must not share a SAM.
            Test.@test reads_cache_tag(a) != reads_cache_tag(b)

            # Deterministic for an unchanged file, or the cache never hits and
            # every invocation re-runs minimap2.
            Test.@test reads_cache_tag(a) == reads_cache_tag(a)

            # A REGENERATED file busts the cache. The sweep's one observed
            # failure was a truncated non-empty .fq.gz being reused; a key that
            # ignored size would keep scoring the stale SAM against it.
            before = reads_cache_tag(a)
            sleep(1.1)          # mtime has 1-second granularity on some filesystems
            write(a, "@r1\nACGTACGT\n+\nIIIIIIII\n")
            Test.@test reads_cache_tag(a) != before
        end
    end

    Test.@testset "cache tags are filesystem-safe" begin
        Test.@test sanitize_for_filename("ont reads (run#3).fq.gz") ==
                   "ont_reads__run_3_.fq.gz"
        Test.@test sanitize_for_filename("../../etc/passwd") == ".._.._etc_passwd"
    end
end

end  # module
