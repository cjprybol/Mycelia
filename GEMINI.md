
# Mycelia Project Overview

This document provides a high-level overview of the Mycelia project, its structure, and development conventions.

## Project Overview

Mycelia is a comprehensive bioinformatics package written in Julia. It focuses on graph-based genome assembly and quality-aware sequence analysis. The project provides a rich set of tools for various bioinformatics workflows, including data acquisition, quality control, assembly, annotation, and comparative genomics. It also integrates with over 30 external bioinformatics tools.

The project is structured as a standard Julia package, with source code in the `src` directory, tests in the `test` directory, and documentation in the `docs` directory. It also includes tutorials, benchmarks, and a command-line interface.

**Key Technologies:**

*   **Language:** Julia
*   **Package Management:** Pkg.jl, Conda.jl
*   **Testing:** Julia's built-in test framework
*   **Documentation:** Documenter.jl, Literate.jl
*   **CI/CD:** GitHub Actions
*   **Containerization:** Docker

## Building and Running

### Development Environment

The recommended development environment is the provided Dev Container, which can be launched using VS Code. This will set up a consistent environment with all the necessary dependencies.

### Running Tests

The project has a tiered testing system:

*   **Core Tests:** These are Julia-only tests that can be run with the following command. They are designed to be fast and are run in CI.

    ```bash
    julia --project=. -e "import Pkg; Pkg.test()"
    ```

*   **Integration Tests:** These tests include integration with external tools and require a Conda environment.

    ```bash
    julia --project=. test/run_integration_tests.jl
    ```

*   **Extended Tests and Tutorials:** The `run_extended_tests.jl` script provides a command-line interface for running tutorials, benchmarks, and generating reports.

    ```bash
    # Run all tutorials
    julia --project=. run_extended_tests.jl tutorials

    # Run all benchmarks (resource intensive)
    julia --project=. run_extended_tests.jl benchmarks
    ```

### Command-Line Interface (CLI)

Mycelia provides a command-line interface for its main functionalities.

```bash
# Display help message
julia --project=. bin/mycelia.jl --help

# Construct a genome graph
julia --project=. bin/mycelia.jl construct --k 31 --out graph.gfa --fastx reads.fasta

# Assemble a genome
julia --project=. bin/mycelia.jl assemble --kmers trusted.kmers --sequences reads.fasta --out assembly.gfa
```

### Building Documentation

The documentation can be built locally using the following command:

```bash
julia --project=docs -e 'include("docs/make.jl")'
```

The generated documentation will be in the `docs/build` directory.

## Development Conventions

*   **Code Style:** The project follows standard Julia coding conventions.
*   **Testing:**
    *   Tests are organized into tiers: `core`, `integration`, and `quick`.
    *   The `TEST_TIER` environment variable can be used to select which tier to run.
    *   New functionality should be accompanied by tests.
*   **Dependencies:**
    *   Julia dependencies are managed by `Pkg.jl`.
    *   External tool dependencies are managed through `Conda.jl`. The `deps/build.jl` script handles the Conda environment setup.
*   **Documentation:**
    *   Tutorials are written as Julia scripts in the `tutorials` directory and are converted to Markdown using `Literate.jl`.
    *   API documentation is generated from docstrings using `Documenter.jl`.
*   **Continuous Integration:** The CI pipeline, defined in `.github/workflows/ci.yml`, runs the core tests and submits code coverage reports to Codecov.
