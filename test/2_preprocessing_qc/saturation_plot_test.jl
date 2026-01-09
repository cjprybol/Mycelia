# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/saturation_plot_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/saturation_plot_test.jl", "test/2_preprocessing_qc", execute=false)'
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
import Colors

Test.@testset "Diversity Saturation Visualization" begin

    # 1. Synthetic Data Generation
    # Simulating a saturation curve (diminishing returns)
    n_points = 100
    x_data = collect(1:n_points)
    # Logarithmic growth function
    y_data = [500 * log(i + 1) for i in 1:n_points]

    # 2. Define Groups/Waves
    # 3 waves of discovery
    groups = vcat(
        fill("Wave 1", 30),
        fill("Wave 2", 30),
        fill("Wave 3", 40)
    )

    # 3. Define References
    ref_val = [2500.0]
    ref_label = ["Theoretical Maximum"]

    Test.@testset "Basic Plot (No Grouping)" begin
        fig = Mycelia.plot_diversity_saturation(
            y_data;
            title="Simple Saturation",
            reference_values=ref_val,
            reference_labels=ref_label
        )
        Test.@test isa(fig, CairoMakie.Figure)
        output_path = "test_saturation_basic.png"
        CairoMakie.save(output_path, fig)
        Test.@test isfile(output_path)
        rm(output_path; force=true)
    end

    Test.@testset "Grouped Plot (Discovery Waves)" begin
        fig = Mycelia.plot_diversity_saturation(
            y_data;
            x_values=x_data,
            grouping_values=groups,
            title="Phage Discovery Waves",
            reference_values=ref_val,
            reference_labels=ref_label
        )
        Test.@test isa(fig, CairoMakie.Figure)
        output_path = "test_saturation_grouped.png"
        CairoMakie.save(output_path, fig)
        Test.@test isfile(output_path)
        rm(output_path; force=true)
    end

    Test.@testset "Multiple Series" begin
        y_data_2 = [300 * log(i + 1) for i in 1:n_points] # Lower diversity

        fig = Mycelia.plot_diversity_saturation(
            [y_data, y_data_2];
            labels=["Method A", "Method B"],
            title="Method Comparison",
            reference_values=ref_val,
            reference_labels=ref_label
        )
        Test.@test isa(fig, CairoMakie.Figure)
        output_path = "test_saturation_multi.png"
        CairoMakie.save(output_path, fig)
        Test.@test isfile(output_path)
        rm(output_path; force=true)
    end
end
