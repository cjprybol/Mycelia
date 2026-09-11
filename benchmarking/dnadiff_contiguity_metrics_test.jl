# Pure unit test for the dnadiff-derived contiguity metrics — NO Mycelia
# dependency, no external tools, runs in ms.
#
# These metrics fill the two gaps the SOTA-parity benchmark needs beyond the
# reference-free N50 + dnadiff identity block the correction-validation harness
# already reports: (1) a misassembly PROXY (relocations + translocations +
# inversions from the dnadiff .report [Feature Estimates] block) and (2) NGA50
# (N50 of reference-aligned block lengths against the reference genome length,
# parsed from the dnadiff .1coords alignment file). QUAST is the canonical
# source for both, but its bioconda build does not solve on osx-arm64, so the
# harness derives them from MUMmer dnadiff output instead.
#
# Run:
#   julia benchmarking/dnadiff_contiguity_metrics_test.jl

import Test

include(joinpath(@__DIR__, "dnadiff_contiguity_metrics.jl"))

# A real MUMmer 3.23 dnadiff .report [Feature Estimates] block (verified against
# a committed lambda run): all structural counts zero, 2 breakpoints. Columns
# are REF then QRY; we take the REF column.
const REPORT_CLEAN = """
[Feature Estimates]
Breakpoints                        2                    0
Relocations                        0                    0
Translocations                     0                    0
Inversions                         0                    0

Insertions                         2                    0
InsertionSum                      28                    0
"""

# A synthetic report exercising nonzero structural misassemblies.
const REPORT_MISASSEMBLED = """
[Feature Estimates]
Breakpoints                       14                    9
Relocations                        3                    2
Translocations                     1                    0
Inversions                         2                    1

Insertions                        11                    4
"""

Test.@testset "dnadiff contiguity metrics" begin
    Test.@testset "feature-estimates parse (REF column)" begin
        fe = parse_dnadiff_feature_estimates(REPORT_CLEAN)
        Test.@test fe.breakpoints == 2
        Test.@test fe.relocations == 0
        Test.@test fe.translocations == 0
        Test.@test fe.inversions == 0

        fe2 = parse_dnadiff_feature_estimates(REPORT_MISASSEMBLED)
        Test.@test fe2.breakpoints == 14
        Test.@test fe2.relocations == 3
        Test.@test fe2.translocations == 1
        Test.@test fe2.inversions == 2
    end

    Test.@testset "misassembly proxy = reloc + transloc + inversions" begin
        # Clean lambda: 0 structural misassemblies (breakpoints are NOT counted —
        # they are alignment endpoints, not misassembly events).
        Test.@test misassembly_proxy(parse_dnadiff_feature_estimates(REPORT_CLEAN)) == 0
        # Synthetic: 3 + 1 + 2 = 6.
        Test.@test misassembly_proxy(parse_dnadiff_feature_estimates(REPORT_MISASSEMBLED)) ==
                   6
    end

    Test.@testset "feature-estimates missing block fails closed" begin
        # A report with no [Feature Estimates] block must be rejected, not
        # silently treated as zero misassemblies (which would mask a broken run).
        Test.@test_throws Exception parse_dnadiff_feature_estimates("NUCMER\n[Bases]\nTotalBases 100 100\n")
    end

    Test.@testset "nga50 (pure) against reference length" begin
        # Single block spanning the whole reference => NGA50 == reference length.
        Test.@test nga50([300], 300) == 300
        # blocks [100,80,60,40,20], ref=300: threshold=150; cumsum 100<150,
        # 100+80=180>=150 => NGA50 = 80.
        Test.@test nga50([20, 100, 40, 80, 60], 300) == 80
        # Broken assembly: aligned blocks sum below half the reference => the
        # 50% threshold is never reached => NGA50 = 0 (QUAST prints "-" here).
        Test.@test nga50([50, 40], 1000) == 0
        # Empty alignment set => 0.
        Test.@test nga50(Int[], 300) == 0
        # Exactly-half boundary is inclusive: block cumsum reaching exactly 50%.
        Test.@test nga50([50, 50], 100) == 50
    end

    Test.@testset "nga50 input guards" begin
        Test.@test_throws Exception nga50([100], 0)     # non-positive reference length
        Test.@test_throws Exception nga50([100], -5)
        Test.@test_throws Exception nga50([-1, 100], 300)  # negative block length
    end

    Test.@testset "aligned ref-block lengths from .1coords" begin
        # MUMmer 3.23 dnadiff .1coords (show-coords -rclHT): tab-separated,
        # columns S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ TAGR TAGQ.
        # Ref block length is derived from S1/E1 (cols 1-2) as abs(E1-S1)+1,
        # robust to reverse-strand ref coords and to trailing-column drift.
        coords = "1\t1000\t1\t1000\t1000\t1000\t99.9\t5000\t1000\t20.0\t100.0\tref\tq1\n" *
                 "2001\t2500\t500\t1\t500\t500\t99.5\t5000\t500\t10.0\t100.0\tref\tq2\n"
        lens = aligned_ref_block_lengths(coords)
        Test.@test lens == [1000, 500]
        # Chaining into nga50: ref=5000, threshold=2500; blocks [1000,500] sum
        # 1500 < 2500 => 0.
        Test.@test nga50(lens, 5000) == 0
        # Blank / headerless-empty coords => no blocks.
        Test.@test aligned_ref_block_lengths("") == Int[]
    end
end
