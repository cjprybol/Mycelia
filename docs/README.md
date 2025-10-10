# Mycelia Documentation

This directory contains the source for Mycelia's documentation, built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).

## Building the Documentation

### Quick Build

The simplest way to build the documentation:

```bash
julia docs/verify_build.jl
```

This script will:
1. Check for required files
2. Build the documentation
3. Report the build status
4. Show the path to the generated HTML

### Manual Build

If you prefer to build manually:

```bash
julia docs/make.jl
```

**Important**: Do NOT use `julia --project=docs docs/make.jl` in this environment. The build system uses a temporary environment to avoid network dependency issues.

### Build Output

The generated documentation will be in `docs/build/`:
- `docs/build/index.html` - Main documentation page
- View locally: `file://path/to/Mycelia/docs/build/index.html`

## Documentation Structure

```
docs/
├── make.jl                 # Main build script
├── verify_build.jl         # Build verification script
├── src/                    # Documentation source files
│   ├── index.md           # Home page
│   ├── installation.md    # Installation guide
│   ├── getting-started.md # Getting started guide
│   ├── api-reference.md   # API reference
│   ├── api/               # Detailed API documentation
│   └── generated/         # Auto-generated tutorial pages
└── build/                  # Generated HTML (gitignored)
```

## Special Build Configuration

This documentation uses a special build configuration to work around network restrictions:

1. **Mock Module**: The build creates a minimal mock Mycelia module instead of loading the full package
2. **Temporary Environment**: Uses `mktempdir()` to create a fresh environment
3. **Minimal Dependencies**: Installs only Documenter and Literate
4. **No Code Execution**: Tutorials are processed with `execute = false`

This approach allows documentation to build even when:
- Network access is restricted
- Some package dependencies are unavailable
- The full package environment cannot be instantiated

## Current Documentation Status

See [../planning-docs/DOCS_ACCURACY_REPORT.md](../planning-docs/DOCS_ACCURACY_REPORT.md) for detailed analysis.

**Key Statistics**:
- **656 public functions** in the codebase
- **11 functions** currently documented with proper @ref links
- **19 tutorial pages** generated from .jl files
- **~177 build warnings** (all non-fatal)

### Known Issues

1. **Missing @ref Targets**: Many functions referenced in documentation don't have docstrings yet
2. **Planned Features**: Some documented functions are marked as "(planned)" and not yet implemented
3. **Tutorial Links**: Some internal links between tutorials need updating
4. **No Code Execution**: Tutorial code is not executed during build (intentional for stability)

## Contributing to Documentation

### Adding Function Documentation

To document a function, add a docstring in the source code:

```julia
"""
    my_function(arg1, arg2)

Description of what the function does.

# Arguments
- `arg1::Type`: Description of arg1
- `arg2::Type`: Description of arg2

# Returns
- `ReturnType`: Description of return value

# Examples
```julia
result = my_function("test", 42)
```
"""
function my_function(arg1, arg2)
    # implementation
end
```

### Adding Documentation Pages

1. Create a new `.md` file in `docs/src/`
2. Add it to the `pages` array in `docs/make.jl`
3. Build and verify

### Adding Tutorials

1. Create a new `.jl` file in `tutorials/`
2. Use Literate.jl comment conventions:
   - Single `#` becomes markdown text
   - Double `##` remains as code comment
3. The build will automatically process it

## Troubleshooting

### Build Fails with Package Errors

If you see errors about missing packages or gitlab.com access:
- This is expected in network-restricted environments
- The build system is designed to work around this
- Make sure you're using `julia docs/make.jl` (not with `--project=docs`)

### Missing Pages

If pages don't appear:
- Check that the file is listed in the `pages` array in `make.jl`
- Verify the file path is correct relative to `docs/src/`
- Check for Markdown syntax errors in the file

### Warnings About @ref Links

These warnings indicate:
- A function is referenced but doesn't have a docstring
- The function name might be misspelled
- The function might be marked as "planned" but not implemented

To fix:
- Add docstrings to the referenced functions
- Update documentation to remove references to non-existent functions
- Mark clearly which functions are planned vs implemented

## Resources

- [Documenter.jl Documentation](https://juliadocs.github.io/Documenter.jl/)
- [Literate.jl Documentation](https://fredrikekre.github.io/Literate.jl/)
- [Julia Documentation Guidelines](https://docs.julialang.org/en/v1/manual/documentation/)
