---
name: julia-style
description: Julia coding style and import rules for Mycelia package
---

# Mycelia Julia Coding Style

When working in the Mycelia repository, follow these strict coding rules.

## Import Rules (STRICT)

- **NEVER use `using`** - always use `import`
- **Fully qualify all names** - use `Mycelia.Rhizomorph`, `DataFrames.DataFrame()`, etc.
- **No shorthand aliases** - avoid `R = Mycelia.Rhizomorph`; use full qualified names
- Dependencies are imported at top-level in `src/Mycelia.jl` and available in all source files

Example:
```julia
# WRONG
using DataFrames
df = DataFrame()

# WRONG
import DataFrames: DataFrame
df = DataFrame()

# CORRECT
import DataFrames
df = DataFrames.DataFrame()
```

## Naming Conventions

- **Modules**: CamelCase (e.g., `Rhizomorph`, `KmerAnalysis`)
- **Functions/Variables**: snake_case (e.g., `process_sequences`, `kmer_count`)
- **Constants**: SCREAMING_SNAKE_CASE (e.g., `MAX_KMER_SIZE`)
- **No emojis** in code or comments

## Sequence Types (BioSequences)

Always use proper BioSequences types:
- `BioSequences.LongDNA{4}` for DNA
- `BioSequences.LongRNA{4}` for RNA
- `BioSequences.LongAA` for amino acids

Extract from FASTQ:
```julia
FASTX.sequence(BioSequences.LongDNA{4}, record)
```

Avoid string conversions in k-mer graphs, qualmer graphs, and assembly algorithms. Use `string()` only for external tool interfaces or final output.

## Formatting

- 4-space indentation
- Clear docstrings: `"""signature..."""`
- Prefer pure functions and explicit keyword arguments
- Use `joinpath` for file paths (portability)
- Avoid type piracy

Default formatter: JuliaFormatter.jl with SciMLStyle or BlueStyle, unless repo indicates Runic.jl.

## Architecture Notes

- Main module uses dynamic file inclusion
- All `.jl` files in `src/` are auto-included via `Mycelia.jl`
- Rhizomorph and Mycelia do NOT export symbols by design
- Tests and callers MUST use fully qualified names
