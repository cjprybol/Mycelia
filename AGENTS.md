# Repository Guidelines

## Overview
Mycelia is a Julia package for bioinformatics and computational biology, providing tools for sequence analysis, genome assembly, annotation, and comparative genomics.

## Project Structure & Module Organization
- `src/` holds the Mycelia module (`Mycelia.jl`) plus domain-specific files and experimental `src/rhizomorph/` work.
- Source files are auto-included from `src/`; add new files with clear names and they will load via `Mycelia.jl`.
- Key domain modules include `fastx.jl`, `assembly.jl`, `annotation.jl`, `alignments-and-mapping.jl`, `kmer-analysis.jl`, `clustering.jl`, `classification.jl`, `quality-control-and-benchmarking.jl`, `plotting-and-visualization.jl`, `utility-functions.jl`, `reference-databases.jl`, `taxonomy-and-trees.jl`, and `variant-analysis.jl`.
- `test/` is tiered by workflow stage; `runtests.jl` pulls in Aqua/JET checks and staged suites. Add new coverage to the closest stage folder.
- Current test stages: `1_data_acquisition/`, `2_preprocessing_qc/`, `3_feature_extraction_kmer/`, `4_assembly/`, `5_validation/`, `6_annotation/`, `7_comparative_pangenomics/`, `8_tool_integration/`. Optional suites live in `in_development/` and `deprecated/`.
- `tutorials/` contains Literate.jl walkthroughs; `benchmarking/` carries heavy scripts; `results/` stores generated artifacts; `docs/` builds the Documenter.jl site; `assembly_test_data/` and `test/metadata/` hold fixtures.

## Build, Test, and Development Commands
- Install deps: `julia --project=. -e "import Pkg; Pkg.instantiate()"`
- Core tests: `julia --project=. -e "import Pkg; Pkg.test()"` (runs Aqua and optional JET).
- Full pipeline tests: `MYCELIA_RUN_ALL=true julia --project=. -e "import Pkg; Pkg.test()"`.
- External tool tests: `MYCELIA_RUN_EXTERNAL=true julia --project=. -e "import Pkg; Pkg.test()"` (no extra flags required).
- Coverage: `julia --project=. --code-coverage=user -e "import Pkg; Pkg.test()"`.
- Static analysis (JET): `julia --project=. test/jet.jl`.
- Docs: `julia --project=docs -e 'include("docs/make.jl")'` (builds into `docs/build`).
- Tutorials sweep: `julia --project=. tutorials/run_all_tutorials.jl` (manual; can be slow).
- Benchmarks: `julia --project=. benchmarking/benchmark_runner.jl small|medium|large` or `sbatch benchmarking/run_all_benchmarks.sh` for SLURM (resource-intensive).

## Coding Style & Naming Conventions
- Follow standard Julia style: 4-space indentation, clear docstrings (`"""signature..."""`), favor pure functions and explicit keyword arguments.
- Names: modules in CamelCase, functions/variables in `snake_case`, constants in `SCREAMING_SNAKE_CASE`; no emojis.
- Imports are strict: never use `using` or import specific functions. Mycelia code should be fully qualified (e.g., `Mycelia.Rhizomorph`). For third-party libraries used by Mycelia, you may either `import` the package directly in the file or reference it via the `Mycelia` namespace (more verbose but acceptable).
- Do not introduce shorthand module aliases as constants (e.g., avoid `R = Mycelia.Rhizomorph`); call modules/functions by their full qualified names to keep tests and sources explicit.
- Prefer `joinpath` for portability, avoid type piracy, and keep external tool calls isolated in helpers under `src/`.

## Sequence Type Guidelines
- Always use BioSequences types: `BioSequences.LongDNA{4}`, `BioSequences.LongRNA{4}`, `BioSequences.LongAA`.
- Extract sequences from FASTQ records with types: `FASTX.sequence(BioSequences.LongDNA{4}, record)`.
- Avoid string conversions in k-mer graphs, qualmer graphs, and assembly algorithms; use `string()` only when interfacing with external tools or final output.

## Architecture & Dependencies
- The main module uses dynamic file inclusion; all `.jl` files in `src/` are automatically included.
- All package dependencies are imported at the top-level in `src/Mycelia.jl` and are available in all source files.
- Key dependencies include BioSequences.jl, FASTX.jl, BioAlignments.jl, Graphs.jl, DataFrames.jl, Makie.jl/Plots.jl, Arrow.jl/JLD2.jl, and XAM.jl.
- Memory estimation/checking utilities live in `utility-functions.jl`; use them for large-scale analyses.

## Tutorial & Benchmark Guidelines
- Do not add generally useful functions inside `tutorials/`, `test/`, or `benchmarking/`. If a helper is broadly useful, add it under `src/` and use it from the tutorial.
- Tutorials should demonstrate core Mycelia functionality; avoid defining new helper functions except true one-offs that only support the tutorial narrative.
- Prefer Rhizomorph graph/assembly functionality when newer variants exist; avoid legacy base-graph/assembly APIs in tutorials unless explicitly requested.
- Literate.jl comment rules: `#` starts markdown, `##` remains a code comment. Use `##` for inline comments inside code blocks.

## Reusable Code Principle
- **Mycelia is the library; notebooks are consumers.** When implementing functionality in external notebooks and/or repositories, always ask: "Is this broadly useful for genomics/metagenomics workflows?"
- If yes, add the function to Mycelia under `src/` with proper documentation and tests, then call `Mycelia.function_name()` from the notebook.
- This applies to: diversity metrics, distance calculations, ordination methods, parsing utilities, visualization helpers, and any other analysis functions that could benefit multiple workflows.
- **Do not duplicate functionality.** Before writing a custom implementation in a notebook, check if Mycelia already has the capability (e.g., `Mycelia.bray_curtis_distance`, `Mycelia.pcoa_from_dist`).
- One-offs specific to a particular dataset or Locus-proprietary logic may remain in notebooks, but general-purpose bioinformatics functions belong in Mycelia.
- When in doubt, add to Mycelia. It's easier to use a library function than to copy-paste between notebooks.

## Testing Guidelines
- Place new tests in the relevant stage directory; name files after the feature (e.g., `assembly_merging.jl`). Use `Test.@testset` with descriptive labels.
- Never disable or skip tests because functionality is broken; fix the implementation first. Use `Test.@test_skip` only for explicitly planned features.
- Use small fixtures from `assembly_test_data/` and deterministic seeds (`StableRNGs`) for reproducibility.
- External tool tests must run when `MYCELIA_RUN_EXTERNAL=true` without requiring extra tool-specific flags; use simulated inputs and default database paths where possible.
- Extended tutorials/benchmarks are opt-in; note runtime or external-tool needs.

## External Tool Integration
- The package integrates with external tools and infrastructure such as Bioconda, SLURM, and Rclone. Keep tool-specific logic isolated in helper modules under `src/`.

## Communication and Documentation Standards
- Be conservative and understated in all statements and claims; avoid overpromising or marketing language.
- Back performance claims with benchmarks and clearly distinguish between completed and planned features.
- Use clear, professional language with plain text markers like "NOTE:" or "WARNING:" and never use emojis.

## Collaboration Notes
- Before large deletions or restructures, explain why the deletion is necessary and what logic is being kept/replaced, then proceed once agreed.
- Do not re-enable deprecated/removed code or revert deprecations without explicit approval.
- Rhizomorph and Mycelia do not export symbols by design; tests and callers must use fully qualified names.

## Commit & Pull Request Guidelines
- Commits: short, present-tense subjects (~72 chars) similar to recent history (e.g., "Refactor path-finding algorithms"). Keep scopes focused.
- PRs: include what/why, tests run (`Pkg.test`, tutorials/benchmarks if applicable), linked issues, and screenshots for plots/notebooks when relevant.
- Note external tool or dataset requirements explicitly; avoid committing generated outputs or large files; keep claims measured and backed by tests/benchmarks.
