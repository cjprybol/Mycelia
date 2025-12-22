import Test
import Mycelia
import CairoMakie

Test.@testset "PCoA Visualization Utils" begin
    
    # 1. Test Ellipse Generation
    # Create a perfect diagonal line of points (highly correlated)
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = [1.0, 2.0, 3.0, 4.0, 5.0]
    
    points, anchors = Mycelia.get_ellipse_points(x, y)
    
    # Check return types
    Test.@test points isa Vector{CairoMakie.Point2f}
    Test.@test anchors isa Dict{Symbol, CairoMakie.Point2f}
    
    # Check anchor keys exist
    Test.@test haskey(anchors, :top)
    Test.@test haskey(anchors, :bottom)
    
    # Check that points were generated (length 100 is hardcoded in function)
    Test.@test length(points) == 100
    
    # 2. Test Plot Generation (Integration)
    # Define minimal series
    series = [Mycelia.PointSeries([1, 2, 3], "Test Group", :blue)]
    
    # Call main function
    fig = Mycelia.plot_generalized_pcoa(
        x, y, series;
        title = "Test Plot",
        output_file = "" # Don't save to disk during test
    )
    
    # Check that a Figure was returned
    Test.@test fig isa CairoMakie.Figure
    
    # 3. Test Empty Data Handling
    # Should handle empty ellipse points gracefully
    empty_x = Float64[]
    pts_empty, _ = Mycelia.get_ellipse_points(empty_x, empty_x)
    Test.@test isempty(pts_empty)
end

# example_pcoa_pseudomonas.jl

import Mycelia
import CairoMakie
import OrderedCollections

# 1. Mock Data Generation (Simulating Pseudomonas analysis)
# Imagine we have 100 strains: 
# - 20 Clinical "Type A" (High Risk)
# - 30 Clinical "Type B" (High Risk)
# - 50 Environmental (Low Risk background)

n_clinical_a = 20
n_clinical_b = 30
n_env = 50
n_total = n_clinical_a + n_clinical_b + n_env

# Generate synthetic PCoA coordinates
# Clinical A: Clustered at (5, 5)
x_a = 5.0 .+ randn(n_clinical_a)
y_a = 5.0 .+ randn(n_clinical_a)

# Clinical B: Clustered at (-2, 8)
x_b = -2.0 .+ randn(n_clinical_b)
y_b = 8.0 .+ randn(n_clinical_b)

# Environmental: Scattered widely around (0, 0)
x_env = randn(n_env) .* 4
y_env = randn(n_env) .* 4

# Combine coordinates
x_all = vcat(x_a, x_b, x_env)
y_all = vcat(y_a, y_b, y_env)

# Generate Indices
indices_a = collect(1:n_clinical_a)
indices_b = collect(n_clinical_a+1 : n_clinical_a+n_clinical_b)
indices_env = collect(n_clinical_a+n_clinical_b+1 : n_total)
indices_clinical = vcat(indices_a, indices_b)

# Map indices to "Species" or Subtypes for shapes
species_map = Dict{Int, String}()
for i in indices_a species_map[i] = "P. aeruginosa Type A" end
for i in indices_b species_map[i] = "P. aeruginosa Type B" end
for i in indices_env species_map[i] = "P. aeruginosa Env" end

# 2. Configuration

# Define Colors
const clinical_red = "#D32F2F"
const env_green = "#388E3C"

# Define Shapes (Marker Map)
marker_map = OrderedCollections.OrderedDict(
    "P. aeruginosa Type A" => :diamond,
    "P. aeruginosa Type B" => :rect,
    "P. aeruginosa Env" => :circle
)

# 3. Create PointSeries (The logical grouping for coloring/layering)
series_list = [
    # Background (Environmental) - Plot first (low z_order)
    Mycelia.PointSeries(
        indices_env,
        "Environmental Reservoirs",
        env_green,
        marker=:circle,
        size=8.0,
        z_order=1
    ),
    # Foreground (Clinical) - Plot on top (high z_order)
    Mycelia.PointSeries(
        indices_clinical,
        "Clinical Outbreak",
        clinical_red,
        marker=:diamond, # will be overridden by marker_map if provided
        size=12.0,
        z_order=2
    )
]

# 4. Define Ellipses (The annotation logic)
ellipses = [
    Mycelia.EllipseConfig(
        indices_a,
        "Type A Cluster",
        clinical_red,
        anchor=:top,
        align=(:center, :bottom),
        offset=(0.0, 10.0)
    ),
    Mycelia.EllipseConfig(
        indices_b,
        "Type B Cluster",
        clinical_red,
        anchor=:right,
        align=(:left, :center),
        offset=(5.0, 0.0)
    )
]

# 5. Execute Visualization
fig = Mycelia.plot_generalized_pcoa(
    x_all, 
    y_all, 
    series_list;
    species_map = species_map,
    marker_map = marker_map,
    ellipses = ellipses,
    title = "Pseudomonas aeruginosa Outbreak Investigation\nPCoA - Pairwise ANI Distance",
    output_file = "pseudomonas_outbreak_pcoa.svg"
)

# Display result
CairoMakie.display(fig)