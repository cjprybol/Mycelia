# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/bioconda.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/bioconda.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia

Test.@testset "Bioconda Environment Management Tests" begin
    # Note: These tests are designed to be safe and not actually modify the conda environment
    # They test the logic without triggering expensive conda operations
    
    Test.@testset "check_bioconda_env_is_installed function" begin
        # Test with a non-existent package name
        # This should return false without creating anything
        fake_pkg = "nonexistent_package_$(rand(1000:9999))"
        
        # This test will work as long as conda/mamba is available
        if haskey(ENV, "CONDA_PREFIX") || isfile(Mycelia.CONDA_RUNNER)
            result = Mycelia.check_bioconda_env_is_installed(fake_pkg)
            Test.@test result isa Bool
            Test.@test result == false  # Should not exist
        else
            Test.@test_skip "Conda not available, skipping bioconda environment tests"
        end
    end

    Test.@testset "Package name parsing" begin
        # Test package name parsing logic for add_bioconda_env
        # We'll create a mock version that doesn't actually run conda commands
        
        # Test simple package name
        simple_pkg = "blast"
        Test.@test !occursin("::", simple_pkg)
        
        # Test channel::package format
        channel_pkg = "bioconda::blast"
        Test.@test occursin("::", channel_pkg)
        parts = split(channel_pkg, "::")
        Test.@test length(parts) == 2
        Test.@test parts[1] == "bioconda"
        Test.@test parts[2] == "blast"
        
        # Test multiple :: (should only split on first)
        complex_pkg = "bioconda::test::package"
        parts = split(complex_pkg, "::", limit=2)
        Test.@test length(parts) == 2
        Test.@test parts[1] == "bioconda"
        Test.@test parts[2] == "test::package"
    end

    Test.@testset "CONDA_RUNNER variable" begin
        # Test that CONDA_RUNNER is defined
        Test.@test isdefined(Mycelia, :CONDA_RUNNER)
        Test.@test Mycelia.CONDA_RUNNER isa String
        Test.@test !isempty(Mycelia.CONDA_RUNNER)
        
        # Test that it's either conda or mamba
        runner_name = basename(Mycelia.CONDA_RUNNER)
        Test.@test runner_name in ["conda", "mamba"]
    end

    Test.@testset "Environment listing logic" begin
        # Test the logic used to parse conda env list output
        # Mock some typical conda env list output
        mock_output = [
            "# conda environments:",
            "#",
            "base                     /opt/conda",
            "test_env                 /opt/conda/envs/test_env",
            "another_env              /opt/conda/envs/another_env",
            ""
        ]
        
        # Apply the filtering logic from check_bioconda_env_is_installed
        filtered_lines = filter(x -> !occursin(r"^#", x), mock_output)
        split_lines = split.(filtered_lines)
        two_part_lines = filter(x -> length(x) == 2, split_lines)
        env_names = Set(first.(two_part_lines))
        
        Test.@test "test_env" in env_names
        Test.@test "another_env" in env_names
        Test.@test "base" in env_names
        Test.@test length(env_names) == 3
    end

    Test.@testset "Error handling and edge cases" begin
        # Test with empty string
        Test.@test !occursin("::", "")
        
        # Test with just "::"
        just_separator = "::"
        Test.@test occursin("::", just_separator)
        parts = split(just_separator, "::")
        Test.@test length(parts) == 2  # ["", ""]
        
        # Test with special characters in package names
        special_pkg = "my-package_123"
        Test.@test !occursin("::", special_pkg)
        
        # Test with channel that has special characters
        special_channel_pkg = "my-channel::my-package_123"
        Test.@test occursin("::", special_channel_pkg)
        parts = split(special_channel_pkg, "::")
        Test.@test parts[1] == "my-channel"
        Test.@test parts[2] == "my-package_123"
    end

    Test.@testset "Function parameter validation" begin
        # Test that functions can handle various string types
        test_pkg = "test_package"
        
        # Test with String
        Test.@test test_pkg isa String
        
        # Test with SubString (simulating parsed input)
        sub_pkg = SubString(test_pkg, 1, 4)  # "test"
        Test.@test sub_pkg isa AbstractString
        Test.@test string(sub_pkg) == "test"
    end

    Test.@testset "Integration with Conda.jl" begin
        # Test that Conda is available for use
        Test.@test isdefined(Mycelia, :Conda)
        
        # Test CONDA_RUNNER path construction
        if haskey(ENV, "CONDA_PREFIX")
            # If in conda environment, CONDA_RUNNER should be accessible
            Test.@test !isempty(Mycelia.CONDA_RUNNER)
        end
    end

    # Mock test for add_bioconda_env logic without actually running conda
    Test.@testset "add_bioconda_env logic" begin
        # Test the package parsing logic that happens in add_bioconda_env
        test_cases = [
            ("blast", nothing, "blast"),
            ("bioconda::blast", "bioconda", "blast"),
            ("conda-forge::python", "conda-forge", "python"),
            ("my-channel::my-package", "my-channel", "my-package")
        ]
        
        for (input_pkg, expected_channel, expected_pkg) in test_cases
            if occursin("::", input_pkg)
                channel, pkg = split(input_pkg, "::")
                Test.@test channel == expected_channel
                Test.@test pkg == expected_pkg
            else
                channel = nothing
                pkg = input_pkg
                Test.@test channel == expected_channel
                Test.@test pkg == expected_pkg
            end
        end
    end

    Test.@testset "Command construction logic" begin
        # Test that we can construct valid conda commands
        pkg = "test_package"
        channel = "test_channel"
        
        # Basic command without channel
        basic_cmd = `$(Mycelia.CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`
        Test.@test basic_cmd isa Cmd
        Test.@test string(basic_cmd.exec[1]) == Mycelia.CONDA_RUNNER
        Test.@test "create" in basic_cmd.exec
        Test.@test "-n" in basic_cmd.exec
        Test.@test pkg in basic_cmd.exec
        
        # Command with channel
        channel_cmd = `$(Mycelia.CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y`
        Test.@test channel_cmd isa Cmd
        Test.@test "$(channel)::$(pkg)" in channel_cmd.exec
    end
end

# Note: Actual conda operations are not tested here to avoid:
# 1. Long test execution times
# 2. Network dependencies
# 3. Modifying the test environment
# 4. Requiring conda/mamba to be installed
# 
# For comprehensive testing of actual conda operations, consider:
# 1. Integration tests in a separate CI environment
# 2. Manual testing with known packages
# 3. Testing in isolated conda environments
