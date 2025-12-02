# Repository Guidelines

## Project Structure & Module Organization
- `src/` holds the Mycelia module (`Mycelia.jl`) plus domain-specific files (assembly, QC, taxonomy, graph utilities) and experimental `src/rhizomorph/` work.
- `test/` is tiered by workflow stage; `runtests.jl` pulls in Aqua/JET checks and staged suites (assembly, validation, tool integration). Add new coverage to the closest stage folder.
- `tutorials/` contains small-data walkthroughs (Literate.jl: `#` for markdown, `##` for code comments); `benchmarking/` carries heavy scripts; `results/` stores generated artifacts; `docs/` builds the Documenter.jl site; `assembly_test_data/` and `test/metadata/` hold fixtures.
- Source files are auto-included from `src/`; add new files with clear names and they will load via `Mycelia.jl`.

## Build, Test, and Development Commands
- Install deps: `julia --project=. -e "using Pkg; Pkg.instantiate()"`
- Core tests: `julia --project=. -e "using Pkg; Pkg.test()"` (runs Aqua and optional JET).
- Tutorials sweep: `julia --project=. run_extended_tests.jl tutorials` (manual; small datasets).
- Benchmarks: `julia --project=. run_extended_tests.jl benchmarks` or `--hpc` to submit via SLURM using `benchmarking/run_all_benchmarks.sh` (resource-intensive).
- Docs: `julia --project=docs docs/make.jl` to build the site into `docs/build`.
- Static analysis: `julia --project=. test/jet.jl` for JET; coverage: `julia --project=. --code-coverage=user -e "using Pkg; Pkg.test()"`.

## Coding Style & Naming Conventions
- Follow standard Julia style: 4-space indentation, clear docstrings (`"""signature..."""`), favor pure functions and explicit keyword arguments.
- Names: modules in CamelCase, functions/variables in `snake_case`, constants in `SCREAMING_SNAKE_CASE`; no emojis.
- Imports: dependencies are imported once in `src/Mycelia.jl`; do not `using` or re-import inside leaf files—fully qualify (e.g., `Test.@test`, `Dates.now()`).
- Sequences: use `BioSequences.LongDNA{4}/LongRNA{4}/LongAA`; avoid string conversions in k-mer/qualmer/assembly code—only convert when interfacing with external tools or outputs.
- Prefer `joinpath` for portability, avoid type piracy, and keep external tool calls isolated in helpers under `src/`.

## Testing Guidelines
- Place new tests in the relevant stage directory; name files after the feature (e.g., `assembly_merging.jl`). Use `Test.@testset` with descriptive labels; do not skip/disable tests for broken functionality—fix the code instead.
- Use small fixtures from `assembly_test_data/` and deterministic seeds (`StableRNGs`) for reproducibility.
- Extended tutorials/benchmarks are opt-in; note runtime or external-tool needs.

## Commit & Pull Request Guidelines
- Commits: short, present-tense subjects (~72 chars) similar to recent history (e.g., "Refactor path-finding algorithms"). Keep scopes focused.
- PRs: include what/why, tests run (`Pkg.test`, tutorials/benchmarks if applicable), linked issues, and screenshots for plots/notebooks when relevant.
- Note external tool or dataset requirements explicitly; avoid committing generated outputs or large files; keep claims measured and backed by tests/benchmarks.
