import Test
import Mycelia

function make_fake_rclone(stdout::AbstractString)
    return make_fake_rclone_script([
        "#!/bin/sh",
        "printf '%s' '$stdout'"
    ])
end

function make_fake_rclone_script(lines::Vector{String})
    fake_bin_dir = mktempdir()
    fake_rclone = joinpath(fake_bin_dir, "rclone")

    open(fake_rclone, "w") do io
        for line in lines
            println(io, line)
        end
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
        Test.@test_throws ErrorException Base.withenv("PATH" => fake_path) do
            Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
        end
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
        Test.@test_throws ErrorException Base.withenv("PATH" => fake_path) do
            Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
        end
        err = Base.withenv("PATH" => fake_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
            end
        end

        Test.@test err isa ErrorException
        Test.@test occursin("Failed to parse `rclone lsjson` output as JSON", sprint(showerror, err))
    end

    Test.@testset "unexpected lsjson structure" begin
        non_array_path = make_fake_rclone("""{"Path":"example.txt","ID":"abc123"}""")
        Test.@test_throws ErrorException Base.withenv("PATH" => non_array_path) do
            Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
        end
        non_array_err = Base.withenv("PATH" => non_array_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
            end
        end
        Test.@test occursin(
            "Expected JSON array from `rclone lsjson`",
            sprint(showerror, non_array_err)
        )

        non_dict_path = make_fake_rclone("""["example.txt"]""")
        Test.@test_throws ErrorException Base.withenv("PATH" => non_dict_path) do
            Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
        end
        non_dict_err = Base.withenv("PATH" => non_dict_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links("fake-remote:path"; return_df = false)
            end
        end
        Test.@test occursin(
            "Expected each entry from `rclone lsjson` to be an object/dict",
            sprint(showerror, non_dict_err)
        )
    end

    Test.@testset "invalid link parameters" begin
        fake_path = make_fake_rclone("""[{"Path":"example.txt","ID":"abc123"}]""")

        Test.@test_throws ArgumentError Base.withenv("PATH" => fake_path) do
            Mycelia.list_gdrive_with_links(
                "fake-remote:path";
                link_type = :custom,
                return_df = false
            )
        end
        missing_template_err = Base.withenv("PATH" => fake_path) do
            capture_error() do
                Mycelia.list_gdrive_with_links(
                    "fake-remote:path";
                    link_type = :custom,
                    return_df = false
                )
            end
        end

        Test.@test_throws ArgumentError Base.withenv("PATH" => fake_path) do
            Mycelia.list_gdrive_with_links(
                "fake-remote:path";
                link_type = :invalid,
                return_df = false
            )
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

    Test.@testset "copy list returns false and cleans temp file on rclone failure" begin
        fake_path = make_fake_rclone_script([
            "#!/bin/sh",
            "exit 1"
        ])
        tempfiles_before = Set(filter(name -> occursin(r"^rclone_sources_.*\.txt$", name), readdir(tempdir())))

        success = Base.withenv("PATH" => fake_path) do
            Mycelia.rclone_copy_list(
                source = "fake-remote:path",
                destination = mktempdir(),
                relative_paths = ["file1.txt", "nested/file2.txt"]
            )
        end

        tempfiles_after = Set(filter(name -> occursin(r"^rclone_sources_.*\.txt$", name), readdir(tempdir())))
        Test.@test !success
        Test.@test tempfiles_after == tempfiles_before
    end
end
