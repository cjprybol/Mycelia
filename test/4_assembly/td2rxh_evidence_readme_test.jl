# Wording contract for the td-2rxh evidence README (td-1qib).
#
# The PR #470 round-6 review deferred three minor README fixes. This test
# pins the two that change the text: the run is named as a commit, not
# "Mycelia HEAD", and the e0.05 total-contig-base multiples are stated
# per arm, not as a single 11–22× span. The third finding (internal IDs
# in public paths) is a recorded decision to leave the names, not a
# wording assertion.
#
# Run:
#   julia --project=. -e 'include("test/4_assembly/td2rxh_evidence_readme_test.jl")'

import Test

const README_PATH = joinpath(
    @__DIR__, "..", "..",
    "benchmarking", "results", "td2rxh-evidence", "27118659", "README.md")

Test.@testset "td2rxh evidence README wording (td-1qib)" begin
    text = read(README_PATH, String)
    Test.@test isfile(README_PATH)

    # The run label already embeds the hash; the prose must still name the
    # commit explicitly so "Mycelia HEAD" is not the only description.
    Test.@test occursin("Mycelia commit e984af449ac3", text)
    Test.@test !occursin("which is Mycelia HEAD", text)

    # At 5% error the arm totals are distinct (nanopore ~11–13×, illumina
    # ~19–22×). A single 11–22× span across both arms is the loose wording
    # the review asked to drop.
    Test.@test occursin(r"11–13×|11-13x|11–13x", text)
    Test.@test occursin(r"19–22×|19-22x|19–22x", text)
    Test.@test !occursin("11–22×", text)
    Test.@test !occursin("11-22x", text)
end

println("✓ td2rxh evidence README wording tests completed")
