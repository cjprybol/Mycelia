import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_stage2_toy_benchmark.jl"))

Test.@testset "Stage-2 dnadiff full-discrepancy metric" begin
    report = """
    [Bases]
    TotalBases 8000 8030
    AlignedBases 8000(100.00%) 8000(99.63%)
    [Alignments]
    AvgIdentity 99.50 99.50
    """
    metrics = parse_dnadiff_discrepancy(report)
    Test.@test metrics.aligned_pct_ref == 100.0
    Test.@test metrics.avg_identity == 99.5
    Test.@test metrics.total_errors_per_100kbp ≈ 875.0

    Test.@test_throws ErrorException parse_dnadiff_discrepancy(
        "TotalBases 8000 8000\nAvgIdentity 100.0")
end
