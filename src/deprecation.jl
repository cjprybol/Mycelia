"""
Deprecation system for Mycelia.jl API migration.

This file provides macros and utilities for managing the transition from a flat namespace
to an organized submodule structure while maintaining backward compatibility.
"""

"""
    @deprecated_toplevel(old_func, new_path)

Create a deprecation shim for functions being moved to submodules.

# Arguments
- `old_func`: The old function name (symbol)
- `new_path`: String describing the new location (e.g., "Mycelia.IO.open_fastx")

# Example
```julia
@deprecated_toplevel open_fastx "Mycelia.IO.open_fastx"
```
"""
macro deprecated_toplevel(old_func, new_path)
    # Parse the new path to extract module and function
    path_parts = split(new_path, ".")
    if length(path_parts) >= 3 && path_parts[1] == "Mycelia"
        submodule = Symbol(path_parts[2])
        func_name = Symbol(path_parts[3])
        quote
            function $(esc(old_func))(args...; kwargs...)
                @warn "$($(string(old_func))) is deprecated. Use $($(new_path)) instead." maxlog=1 _id=$(hash((old_func, new_path)))
                return $(esc(submodule)).$(esc(func_name))(args...; kwargs...)
            end
        end
    else
        error("Invalid new_path format. Expected 'Mycelia.SubModule.function_name'")
    end
end

"""
    @deprecated_internal(old_func, new_path)

Create a deprecation shim for internal functions being renamed with underscore prefix.

# Arguments  
- `old_func`: The old function name (symbol)
- `new_path`: String describing the new location (e.g., "Mycelia.IO._fastx_type")

# Example
```julia
@deprecated_internal fastx_type "Mycelia.IO._fastx_type (note: this function will become private)"
```
"""
macro deprecated_internal(old_func, new_path)
    quote
        function $(esc(old_func))(args...; kwargs...)
            @warn "$($(string(old_func))) is deprecated and will become private. Use $($(new_path)) instead." maxlog=1 _id=$(hash((old_func, new_path)))
            return $(esc(Symbol("_" * string(old_func))))(args...; kwargs...)
        end
    end
end

"""
    @deprecated_moved(old_func, new_module, new_func)

Create a deprecation shim for functions moved to a different submodule.

# Arguments
- `old_func`: The old function name (symbol) 
- `new_module`: The new submodule (symbol)
- `new_func`: The new function name in the submodule (symbol)

# Example
```julia
@deprecated_moved count_kmers Sequences count_kmers
```
"""
macro deprecated_moved(old_func, new_module, new_func)
    quote
        function $(esc(old_func))(args...; kwargs...)
            @warn "$($(string(old_func))) has moved to $($(string(new_module))).$($(string(new_func))). Please update your code." maxlog=1 _id=$(hash((old_func, new_module, new_func)))
            return $(esc(new_module)).$(esc(new_func))(args...; kwargs...)
        end
    end
end