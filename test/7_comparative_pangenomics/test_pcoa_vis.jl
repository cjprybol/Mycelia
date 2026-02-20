# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/test_pcoa_vis.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/test_pcoa_vis.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

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

    # 2b. Legend marker shapes follow PointSeries markers
    legend_series = [
        Mycelia.PointSeries([1], "Circle Group", :blue; marker = :circle),
        Mycelia.PointSeries([2], "Diamond Group", :red; marker = :diamond)
    ]
    legend_elements = Mycelia.build_pcoa_group_legend_elements(legend_series)
    Test.@test length(legend_elements) == 2
    Test.@test legend_elements[1].marker[] == :circle
    Test.@test legend_elements[2].marker[] == :diamond

    # 3. Test Empty Data Handling
    # Should handle empty ellipse points gracefully
    empty_x = Float64[]
    pts_empty, _ = Mycelia.get_ellipse_points(empty_x, empty_x)
    Test.@test isempty(pts_empty)
end

# example_pcoa_visualization.jl

import Mycelia
import CairoMakie
import OrderedCollections

# 1. Mock Data Generation (Simulating microbial community analysis)
# Imagine we have 100 samples:
# - 20 from Treatment Group A (clustered)
# - 30 from Treatment Group B (clustered)
# - 50 from Control Group (dispersed background)

n_group_a = 20
n_group_b = 30
n_control = 50
n_total = n_group_a + n_group_b + n_control

# Generate synthetic PCoA coordinates
# Group A: Clustered at (5, 5)
x_a = 5.0 .+ randn(n_group_a)
y_a = 5.0 .+ randn(n_group_a)

# Group B: Clustered at (-2, 8)
x_b = -2.0 .+ randn(n_group_b)
y_b = 8.0 .+ randn(n_group_b)

# Control: Scattered widely around (0, 0)
x_control = randn(n_control) .* 4
y_control = randn(n_control) .* 4

# Combine coordinates
x_all = vcat(x_a, x_b, x_control)
y_all = vcat(y_a, y_b, y_control)

# Generate Indices
indices_a = collect(1:n_group_a)
indices_b = collect((n_group_a + 1):(n_group_a + n_group_b))
indices_control = collect((n_group_a + n_group_b + 1):n_total)
indices_treatment = vcat(indices_a, indices_b)

# Map indices to species/subtypes for shapes
species_map = Dict{Int, String}()
for i in indices_a
    species_map[i] = "Species Alpha"
end
for i in indices_b
    species_map[i] = "Species Beta"
end
for i in indices_control
    species_map[i] = "Species Gamma"
end

# 2. Configuration

# Define Colors
const highlight_red = "#D32F2F"
const background_green = "#388E3C"

# Define Shapes (Marker Map)
marker_map = OrderedCollections.OrderedDict(
    "Species Alpha" => :diamond,
    "Species Beta" => :rect,
    "Species Gamma" => :circle
)

# 3. Create PointSeries (The logical grouping for coloring/layering)
series_list = [
    # Background (Control) - Plot first (low z_order)
    Mycelia.PointSeries(
        indices_control,
        "Control Group",
        background_green,
        marker = :circle,
        size = 8.0,
        z_order = 1
    ),
    # Foreground (Treatment) - Plot on top (high z_order)
    Mycelia.PointSeries(
        indices_treatment,
        "Treatment Group",
        highlight_red,
        marker = :diamond, # will be overridden by marker_map if provided
        size = 12.0,
        z_order = 2
    )
]

# 4. Define Ellipses (The annotation logic)
ellipses = [
    Mycelia.EllipseConfig(
        indices_a,
        "Cluster 1",
        highlight_red,
        anchor = :top,
        align = (:center, :bottom),
        offset = (0.0, 10.0)
    ),
    Mycelia.EllipseConfig(
        indices_b,
        "Cluster 2",
        highlight_red,
        anchor = :right,
        align = (:left, :center),
        offset = (5.0, 0.0)
    )
]

example_dir = mktempdir()
output_file = joinpath(example_dir, "pcoa_visualization_example.svg")
try
    # 5. Execute Visualization
    fig = Mycelia.plot_generalized_pcoa(
        x_all,
        y_all,
        series_list;
        species_map = species_map,
        marker_map = marker_map,
        ellipses = ellipses,
        title = "Microbial Community Structure\nPCoA - Distance Matrix",
        output_file = output_file
    )

    # Display result
    CairoMakie.display(fig)
finally
    isfile(output_file) && rm(output_file; force = true)
    rm(example_dir; recursive = true, force = true)
end
