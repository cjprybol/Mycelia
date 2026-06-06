import Test
import Mycelia

function make_fake_rclone(stdout::AbstractString)
    fake_bin_dir = mktempdir()
    fake_rclone = joinpath(fake_bin_dir, "rclone")

    open(fake_rclone, "w") do io
        println(io, "#!/bin/sh")
        println(io, "printf '%s' '$stdout'")
    end
    Base.Filesystem.chmod(fake_rclone, 0o755)

    return fake_bin_dir
end

function capture_error(f::Function)
    try
        f()
    catch err
        return err
    end
    return nothing
end

Test.@testset "rclone wrapper error handling" begin
    Test.@testset "missing rclone binary" begin
        fake_path = mktempdir()
        err = Base.withenv("PATH" => fake_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
            end
        end

        Test.@test err isa ErrorException
        Test.@test occursin("Failed to run `rclone lsjson`", sprint(showerror, err))
    end

    Test.@testset "malformed lsjson output" begin
        fake_path = make_fake_rclone("{not-json}")
        err = Base.withenv("PATH" => fake_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
            end
        end

        Test.@test err isa ErrorException
        Test.@test occursin("Failed to parse `rclone lsjson` output as JSON", sprint(showerror, err))
    end

    Test.@testset "invalid link parameters" begin
        fake_path = make_fake_rclone("""[{"Path":"example.txt","ID":"abc123"}]""")

        missing_template_err = Base.withenv("PATH" => fake_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links(
                    "fake-remote:path";
                    link_type = :custom,
                    return_df = false
                )
            end
        end

        unknown_link_type_err = Base.withenv("PATH" => fake_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links(
                    "fake-remote:path";
                    link_type = :invalid,
                    return_df = false
                )
            end
        end

        Test.@test missing_template_err isa ArgumentError
        Test.@test occursin(
            "link_template must be provided when link_type == :custom",
            sprint(showerror, missing_template_err)
        )

        Test.@test unknown_link_type_err isa ArgumentError
        Test.@test occursin("Unknown link_type: invalid", sprint(showerror, unknown_link_type_err))
    end
end
