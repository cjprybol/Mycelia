# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/diversity_sampling.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/diversity_sampling.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import LinearAlgebra
import Clustering
import Distances
import CairoMakie
import Random

Test.@testset "Diversity-Based Sampling" begin

    # --- Setup Synthetic Data ---
    # Create 3 Groups:
    # Group A: 3 points close together (Indices 1, 2, 3)
    # Group B: 2 points far apart (Indices 4, 5)
    # Group C: 1 point isolated (Index 6)

    # Distance Matrix (lower = closer)
    # 1-3 are a tight cluster (dist ~0.1)
    # 4-5 are dist ~0.8
    # 6 is dist ~1.0 from everyone
    dist_matrix = [
        0.0 0.1 0.1 0.8 0.8 1.0; # 1
        0.1 0.0 0.1 0.8 0.8 1.0; # 2
        0.1 0.1 0.0 0.8 0.8 1.0; # 3
        0.8 0.8 0.8 0.0 0.9 1.0; # 4 (Note: 4-5 is 0.9 dist)
        0.8 0.8 0.8 0.9 0.0 1.0; # 5
        1.0 1.0 1.0 1.0 1.0 0.0  # 6
    ]

    group_ids = ["A", "A", "A", "B", "B", "C"]
    weights = [10.0, 5.0, 1.0, 10.0, 10.0, 10.0] # Weights/Abundance

    Test.@testset "Medoid Selection (Distance)" begin
        # For Group A (1,2,3), all have equal sum distances (0.2).
        # Should pick index 1 (first appearance) or logic specific
        medoid = Mycelia.find_group_medoid(dist_matrix, [1, 2, 3], metric_type=:distance)
        Test.@test medoid in [1, 2, 3]

        # For Group B (4,5), symmetric.
        medoid_b = Mycelia.find_group_medoid(dist_matrix, [4, 5], metric_type=:distance)
        Test.@test medoid_b in [4, 5]
    end

    Test.@testset "Stage 1 & Constraints" begin
        # Requesting n_total=3, max_per_group=1
        # Should select exactly one from A, B, and C
        selected = Mycelia.select_diverse_representatives(
            dist_matrix, group_ids;
            n_total=3, max_per_group=1, metric_type=:distance
        )
        Test.@test length(selected) == 3
        selected_groups = [group_ids[i] for i in selected]
        Test.@test sort(selected_groups) == ["A", "B", "C"]
    end

    Test.@testset "Stage 2 Diversity Filling" begin
        # Request n_total=5, max_per_group=2
        # Stage 1 picks 1 medoid from A, B, C (3 total)
        # Stage 1.5 picks 2nd diverse from A (if dist > 0), B (dist 0.9)
        # Then fills remaining
        selected = Mycelia.select_diverse_representatives(
            dist_matrix, group_ids;
            n_total=5, max_per_group=2, metric_type=:distance
        )

        # Must pick all of B (since they are very diverse, 0.9 dist)
        Test.@test 4 in selected && 5 in selected
        # Must pick 6 (only one in C)
        Test.@test 6 in selected
        # Must pick 2 from A
        count_a = count(i -> group_ids[i] == "A", selected)
        Test.@test count_a == 2
    end

    Test.@testset "Similarity Metric (ANI) Support" begin
        # Convert distance to similarity (1.0 - dist)
        sim_matrix = 1.0 .- dist_matrix

        # Logic check:
        # A (1,2,3) are 0.9 similarity (High)
        # B (4,5) are 0.1 similarity (Low)

        # Medoid for B in similarity mode should maximize similarity sum.
        # Since only 2 points, both are equal.

        # Run selection
        selected = Mycelia.select_diverse_representatives(
            sim_matrix, group_ids;
            n_total=3, max_per_group=1, metric_type=:similarity
        )

        Test.@test length(selected) == 3
        Test.@test sort([group_ids[i] for i in selected]) == ["A", "B", "C"]
    end

    Test.@testset "Abundance Tie-Breaking" begin
        # Create a case where distances are identical
        # Group D: 7, 8. Distance 0.
        dist_ties = [0.0 0.0; 0.0 0.0]
        g_ties = ["D", "D"]
        # 7 has weight 100, 8 has weight 1
        w_ties = [100.0, 1.0]

        # Should pick 7 because of higher weight
        # Note: We need to pass indices 1,2 relative to this matrix
        # But the function handles indices relative to the matrix rows
        medoid = Mycelia.find_group_medoid(dist_ties, [1, 2], metric_type=:distance)
        # In a perfect tie of scores (0.0), implementation usually takes first. 
        # But GreedyMaxMin uses weights.

        # Let's test greedy filling specifically
        # Already selected: nothing. n=1.
        selected = Mycelia.greedy_maxmin_diversity(
            dist_ties, Int[], [1, 2], 1, w_ties, g_ties, 1, metric_type=:distance
        )
        Test.@test selected[1] == 1 # Should pick index 1 (Weight 100)
    end
end

Test.@testset "Radial Dendrogram Visualization" begin

    # 1. Create Synthetic Data
    # Generate 50 points in 2D space forming 3 distinct clusters
    Random.seed!(42)
    # Cluster 1: centered at (0,0)
    # Cluster 2: centered at (5,5)
    # Cluster 3: centered at (10,10)
    points = hcat(randn(2, 20), randn(2, 20) .+ 5, randn(2, 10) .+ 10)

    # Compute Euclidean distance matrix
    dist_matrix = Distances.pairwise(Distances.Euclidean(), points)

    # Perform Hierarchical Clustering using Ward's linkage
    hcl = Clustering.hclust(dist_matrix, linkage=:ward)

    # 2. Define "Rings" (Annotation Layers)
    # We map specific indices to visual properties
    rings = [
        (
            indices=collect(1:20),
            label="Cluster A (Negative)",
            color=:red,
            marker=:circle,
            size=8
        ),
        (
            indices=collect(21:40),
            label="Cluster B (Positive)",
            color=:blue,
            marker=:rect,
            size=8
        ),
        (
            indices=collect(41:50),
            label="Cluster C (Outliers)",
            color=:green,
            marker=:utriangle,
            size=10
        )
    ]

    # 3. Generate Plot
    # We call the function fully qualified
    fig = Mycelia.plot_radial_dendrogram(
        hcl;
        rings=rings,
        title="Synthetic Radial Dendrogram Test",
        inner_radius=0.2,
        outer_radius=1.0,
        ring_spacing=0.1,
        line_color=(:black, 0.5),
        line_width=1.5
    )

    # 4. Verify Output
    Test.@test isa(fig, CairoMakie.Figure)

    # 5. Save Artifact (Optional, for visual inspection)
    output_path = "radial_test_output.png"
    CairoMakie.save(output_path, fig)
    Test.@test isfile(output_path)
    rm(output_path; force=true)
end
