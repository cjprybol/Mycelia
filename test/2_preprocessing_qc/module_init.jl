import Test
import Mycelia

Test.@testset "Module __init__ and LD_LIBRARY_PATH clearing" begin
    Test.@testset "__init__ function exists" begin
        Test.@test isdefined(Mycelia, :__init__)
    end

    Test.@testset "_clear_ld_library_path! clears non-empty value" begin
        original = get(ENV, "LD_LIBRARY_PATH", nothing)
        try
            ENV["LD_LIBRARY_PATH"] = "/usr/lib:/usr/local/lib"
            Mycelia._clear_ld_library_path!()
            Test.@test ENV["LD_LIBRARY_PATH"] == ""
        finally
            if original === nothing
                delete!(ENV, "LD_LIBRARY_PATH")
            else
                ENV["LD_LIBRARY_PATH"] = original
            end
        end
    end

    Test.@testset "_clear_ld_library_path! emits warning with previous value" begin
        original = get(ENV, "LD_LIBRARY_PATH", nothing)
        try
            ENV["LD_LIBRARY_PATH"] = "/usr/lib"
            Test.@test_logs (:warn, r"Previous value: /usr/lib") Mycelia._clear_ld_library_path!()
        finally
            if original === nothing
                delete!(ENV, "LD_LIBRARY_PATH")
            else
                ENV["LD_LIBRARY_PATH"] = original
            end
        end
    end

    Test.@testset "_clear_ld_library_path! no-op when key is unset" begin
        original = get(ENV, "LD_LIBRARY_PATH", nothing)
        try
            delete!(ENV, "LD_LIBRARY_PATH")
            Test.@test_logs Mycelia._clear_ld_library_path!()
            Test.@test !haskey(ENV, "LD_LIBRARY_PATH")
        finally
            if original !== nothing
                ENV["LD_LIBRARY_PATH"] = original
            end
        end
    end

    Test.@testset "_clear_ld_library_path! no-op when value is empty" begin
        original = get(ENV, "LD_LIBRARY_PATH", nothing)
        try
            ENV["LD_LIBRARY_PATH"] = ""
            Test.@test_logs Mycelia._clear_ld_library_path!()
            Test.@test ENV["LD_LIBRARY_PATH"] == ""
        finally
            if original === nothing
                delete!(ENV, "LD_LIBRARY_PATH")
            else
                ENV["LD_LIBRARY_PATH"] = original
            end
        end
    end
end
