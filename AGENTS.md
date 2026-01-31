# AGENTS.md - Mycelia Agent Navigation Guide

Comprehensive guidance for AI agents working in this repository.

## Quick Reference

| Task | Command |
|------|---------|
| Run tests | `julia --project=. -e "import Pkg; Pkg.test()"` |
| Run all tests | `MYCELIA_RUN_ALL=true julia --project=. -e "import Pkg; Pkg.test()"` |
| Build docs | `julia --project=docs -e 'include("docs/make.jl")'` |
| Static analysis | `julia --project=. test/jet.jl` |
| Install deps | `julia --project=. -e "import Pkg; Pkg.instantiate()"` |

---

## Repository Overview

Mycelia is a Julia package for bioinformatics and computational biology, providing tools for sequence analysis, genome assembly, annotation, and comparative genomics. The core assembly system is **Rhizomorph** - a probabilistic graph-based genome assembler.

### Key Architectural Principles

1. **No `using` statements** - Always `import Package` and use qualified names (`DataFrames.DataFrame()`)
2. **No exports** - All symbols accessed via `Mycelia.function_name()` or `Mycelia.Rhizomorph.function_name()`
3. **Avoid unnecessary type conversions** - Use BioSequences for DNA/RNA/protein; plain strings for NLP/text
4. **Test-first development** - All features require tests before claiming completion

---

## Codebase Map

### Source Files (`src/`) - 53 modules

#### Core Infrastructure
| Module | Purpose | Functions |
|--------|---------|-----------|
| `Mycelia.jl` | Main module, auto-includes all files | 2 |
| `constants.jl` | Shared constants | 1 |
| `alphabets.jl` | Sequence alphabet utilities | 4 |
| `utility-functions.jl` | Shared helpers, memory estimation | 160 |
| `checkpointing.jl` | JLD2-based stage caching for workflows | 5 |

#### Sequence I/O & Processing
| Module | Purpose | Functions |
|--------|---------|-----------|
| `fastx.jl` | FASTA/FASTQ reading, normalization | 85 |
| `xam.jl` | BAM/CRAM processing, coverage | 37 |
| `genome-features.jl` | GFF/GTF parsing | 18 |

#### K-mer & Graph Analysis
| Module | Purpose | Functions |
|--------|---------|-----------|
| `kmer-analysis.jl` | K-mer counting, spectrum analysis | 58 |
| `kmer-saturation-analysis.jl` | Sequencing depth estimation | 8 |
| `qualmer-analysis.jl` | Quality-weighted k-mers | 21 |
| `distance-metrics.jl` | Mash, Jaccard, JS divergence | 21 |
| `coverage-clustering.jl` | Coverage-based binning | 7 |

#### Assembly (Rhizomorph System)
| Module | Purpose | Functions |
|--------|---------|-----------|
| `rhizomorph/rhizomorph.jl` | Submodule entry point | 2 |
| `rhizomorph/core/*.jl` | Graph types, evidence, quality | ~120 |
| `rhizomorph/fixed-length/*.jl` | K-mer, qualmer, n-gram graphs | ~30 |
| `rhizomorph/variable-length/*.jl` | OLC graphs from FASTA/FASTQ | ~30 |
| `rhizomorph/algorithms/*.jl` | Path-finding, simplification, repeats | ~80 |
| `rhizomorph/assembly.jl` | High-level assembly orchestration | 54 |
| `assembly.jl` | External assembler wrappers | 62 |
| `graph-cleanup.jl` | Graph simplification utilities | 22 |
| `iterative-assembly.jl` | Iterative refinement methods | 70 |

#### External Tool Wrappers
| Module | Tool(s) | Status |
|--------|---------|--------|
| `alignments-and-mapping.jl` | BLAST, DIAMOND, MMSeqs2, minimap2 | Tested |
| `annotation.jl` | Prokka, Bakta, eggNOG-mapper | Tested |
| `quality-control-and-benchmarking.jl` | QUAST, BUSCO, fastp | Tested |
| `classification.jl` | Kraken2, MetaPhlAn, Sylph | Tested |
| `binning.jl` | MetaBAT2, CONCOCT, SemiBin | Tested |
| `pangenome-analysis.jl` | PGGB, Cactus | Tested |
| `autocycler.jl`, `bcalm.jl`, `ggcat.jl` | Specialized assemblers | Opt-in |

#### Analysis & Visualization
| Module | Purpose | Functions |
|--------|---------|-----------|
| `clustering.jl` | Sequence clustering | 21 |
| `dimensionality-reduction.jl` | PCA, UMAP, t-SNE | 8 |
| `taxonomy-and-trees.jl` | Taxonomy classification, trees | 52 |
| `plotting-and-visualization.jl` | Makie/Plots visualization | 68 |
| `simulation.jl` | Read simulation | 73 |

#### Infrastructure
| Module | Purpose | Functions |
|--------|---------|-----------|
| `reference-databases.jl` | NCBI/SRA downloads, database setup | 84 |
| `ncbi-datasets-cli.jl` | NCBI datasets API | 44 |
| `bioconda.jl` | Conda environment management | 10 |
| `slurm-sbatch.jl` | SLURM job submission | 8 |
| `rclone.jl` | Cloud storage sync | 9 |

---

## Test Structure

Tests are organized by workflow stage (8 tiers):

```
test/
├── 1_data_acquisition/      # NCBI downloads, read simulation
├── 2_preprocessing_qc/      # Alphabet handling, QC filtering
├── 3_feature_extraction_kmer/  # K-mer analysis, saturation
├── 4_assembly/              # 100 test files - core assembly tests
├── 5_validation/            # QUAST, BUSCO validation
├── 6_annotation/            # Gene prediction tests
├── 7_comparative_pangenomics/  # Comparative analysis
├── 8_tool_integration/      # External tool wrapper tests
├── in_development/          # Experimental tests (opt-in)
└── deprecated/              # Legacy tests
```

### Test Categories

**Julia-Only (Fast, Local)** - Safe for parallel agents:
- `test/2_preprocessing_qc/alphabets.jl`
- `test/2_preprocessing_qc/constants.jl`
- `test/3_feature_extraction_kmer/*.jl`
- `test/4_assembly/path_finding_test.jl`
- `test/4_assembly/graph_cleanup_test.jl`
- `test/4_assembly/rhizomorph_*_test.jl` (most)
- `test/4_assembly/dna_kmer_*_test.jl`

**External Tool Required** - Need `MYCELIA_RUN_EXTERNAL=true`:
- `test/5_validation/*.jl` (QUAST, BUSCO)
- `test/6_annotation/*.jl` (Prokka, Bakta)
- `test/8_tool_integration/*.jl` (all wrapper tests)

**Identifying Julia-Only Tests**:
- No `run_*` wrapper functions (e.g., `run_megahit`, `run_minimap2`)
- No `Base.run` or `Mycelia.run` calls
- Only Julia stdlib + registered packages

---

## Common Agent Tasks

### 1. TODO Comment Cleanup

**Find TODOs:**
```bash
grep -r "TODO\|FIXME\|XXX\|HACK" src/ --include="*.jl"
```

**Current TODOs (17 across 11 files):**
- `sequence-comparison.jl` (4) - Use constants instead of hardcoding
- `kmer-analysis.jl` (2) - Ambiguity handling, sequence type inference
- `annotation.jl` (1) - EGAPx wrapper placeholder
- `bioconda.jl` (2) - Stub module
- `rclone.jl` (1) - Stub module
- `plotting-and-visualization.jl` (1) - Incomplete write function

**Approach:**
1. Simple TODOs: Resolve directly
2. Complex TODOs: Create Beads issue with `workspace:mycelia` label
3. Stale TODOs: Remove if already implemented

### 2. Test Expansion

**Run existing tests:**
```bash
julia --project=. -e "import Pkg; Pkg.test()"
```

**Add new test:**
1. Create file in appropriate stage directory (e.g., `test/4_assembly/`)
2. Use `Test.@testset` with descriptive name
3. Use small fixtures from `assembly_test_data/` or `test/metadata/`
4. Use `StableRNGs` for deterministic seeds

**Test template:**
```julia
import Test
import StableRNGs
import Mycelia

Test.@testset "Feature Name Tests" begin
    rng = StableRNGs.StableRNG(42)

    Test.@testset "specific behavior" begin
        # Test code here
        Test.@test expected == actual
    end
end
```

### 3. Documentation Updates

**Build docs:**
```bash
julia --project=docs -e 'include("docs/make.jl")'
```

**Documentation structure:**
- `docs/src/tutorials/` - Literate.jl tutorials (8 stages)
- `docs/src/api/` - API reference (auto-generated)
- `docs/src/workflow-map.md` - Tool/capability matrix

**Docstring format:**
```julia
"""
    function_name(arg1, arg2; kwarg=default)

Brief description of what the function does.

# Arguments
- `arg1`: Description
- `arg2`: Description

# Returns
- Description of return value

# Example
```julia
result = Mycelia.function_name(x, y)
```
"""
function function_name(arg1, arg2; kwarg=default)
    # implementation
end
```

### 4. Code Formatting

**Check style:** Follow BlueStyle conventions
- 4-space indentation
- `snake_case` for functions/variables
- `CamelCase` for types/modules
- `SCREAMING_SNAKE_CASE` for constants

**Apply Runic formatting (if configured):**
```bash
julia --project=. -e "import Runic; Runic.format(\"src/\")"
```

---

## Planning Documents Index

Located in `planning-docs/`:

| Document | Purpose | Lines |
|----------|---------|-------|
| `TODO.md` | Active execution roadmap, checkbox items | 889 |
| `COMPREHENSIVE_ROADMAP.md` | Strategic 5-phase vision | 486 |
| `DEVELOPMENT_TRIAGE.md` | Code promotion decisions | 84 |
| `COMPREHENSIVE_TESTING_FRAMEWORK.md` | Test organization spec | 1049 |
| `FUNCTION_COVERAGE_AUDIT.md` | Module→doc mapping | 113 |
| `TOOL_WRAPPER_STATUS.md` | External tool coverage | 443 |
| `RHIZOMORPH_SUPPORT_MATRIX.md` | Graph type × alphabet matrix | - |
| `rhizomorph-graph-ecosystem-plan.md` | Deep technical spec | 4101 |
| `HPC_CI_PLAN.md` | SLURM CI/CD strategy | - |
| `TDA_INTEGRATION_PLAN.md` | Topology data analysis | - |

**Key finding from TODO.md**: Test-first approach discovered real bugs (4/42 path-finding failures initially). All 302 assembly tests now pass after fixes.

---

## Rhizomorph Graph System

### Graph Type Hierarchy

```
                        Graph Types
                             │
            ┌────────────────┴────────────────┐
            ▼                                 ▼
      Fixed-Length                     Variable-Length (OLC)
            │                                 │
    ┌───────┼───────┐             ┌───────────┼───────────┐
    ▼       ▼       ▼             ▼           ▼           ▼
  K-mer  Qualmer  N-gram       FASTA       FASTQ       String
                              (DNA/RNA)   (DNA/RNA)   (Unicode)
```

**Fixed-Length**: Unit size determined by k (k-mers) or n (n-grams)
**Variable-Length**: Overlap-Layout-Consensus from full sequences

### Strand Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| SingleStrand | Forward only | Proteins, directed graphs, text/NLP |
| DoubleStrand | Forward + RC separate | Strand-specific RNA-seq |
| Canonical | Forward + RC merged | DNA assembly (default) |

### Support Matrix (Graph × Alphabet × Strand)

| Graph Type | DNA | RNA | AA | String/Unicode |
|------------|-----|-----|-----|----------------|
| **K-mer** | S/D/C | S/D/C | S only | - |
| **Qualmer** | S/D/C | S/D/C | - | - |
| **N-gram** | - | - | - | S only |
| **FASTA OLC** | S/D/C | S/D/C | S only | S only |
| **FASTQ OLC** | S/D/C | S/D/C | - | - |
| **String OLC** | - | - | - | S only |

*S = SingleStrand, D = DoubleStrand, C = Canonical*
*Reverse complement only applies to DNA/RNA; AA and String use SingleStrand only*

### Key API Patterns

```julia
# Build a k-mer graph
graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(sequences, k)
graph = Mycelia.Rhizomorph.build_kmer_graph_canonical(sequences, k)

# Find paths
paths = Mycelia.Rhizomorph.find_eulerian_paths(graph)
sequence = Mycelia.Rhizomorph.path_to_sequence(path, graph)

# Simplify
simplified = Mycelia.Rhizomorph.simplify_graph(graph)
```

---

## Coding Conventions

### Import Style (Strict)
```julia
# CORRECT
import BioSequences
import DataFrames
seq = BioSequences.LongDNA{4}("ACGT")
df = DataFrames.DataFrame(a=[1,2,3])

# WRONG - Never use `using`
using BioSequences
using DataFrames: DataFrame
```

### Sequence Types

**Principle:** Avoid unnecessary conversions. Use the appropriate type for the domain.

```julia
# DNA/RNA sequences - Use BioSequences for efficiency
seq = BioSequences.LongDNA{4}("ACGT")
seq = FASTX.sequence(BioSequences.LongDNA{4}, record)

# N-grams / Natural language - Use plain strings
text = "The quick brown fox"
ngram_graph = Mycelia.Rhizomorph.build_ngram_graph([text], n=3)

# String OLC graphs - Use plain strings for Unicode text
string_graph = Mycelia.Rhizomorph.build_string_graph(sentences)
```

**When to use BioSequences:** DNA, RNA, protein sequences where k-mer operations, reverse complement, or biological alphabet constraints apply.

**When to use strings:** N-gram graphs for NLP, text analysis, or any non-biological sequence data.

### Function Calls
```julia
# CORRECT - Fully qualified
result = Mycelia.function_name(args...)
graph = Mycelia.Rhizomorph.build_kmer_graph_canonical(seqs, k)

# WRONG - No aliases
R = Mycelia.Rhizomorph  # Don't do this
result = R.build_kmer_graph_canonical(seqs, k)
```

---

## Commit Guidelines

**Format:** `<type>: <description>`

Types: `add`, `fix`, `update`, `remove`, `refactor`

**Examples:**
```
add: checkpointing utilities for JLD2-based stage caching
fix: path_to_sequence API mismatch in Rhizomorph
update: improve heatmap visualization layout
refactor: extract validation logic to separate module
```

**Rules:**
- Short present-tense subjects (~72 chars)
- Include what/why in body if needed
- Link related issues
- No co-author lines

---

## Beads Integration

Mycelia tasks are tracked in the central todo repo (`~/workspace/todo/.beads/`).

**Labels:**
- `workspace:mycelia` - All Mycelia tasks
- `compute:local` - Can run on laptop
- `compute:scg|lrc|nersc` - Needs HPC

**Finding Mycelia work:**
```bash
cd ~/workspace/todo
bd ready --json | jq '.[] | select(.title | test("Mycelia|\\[B\\]"))'
```

**Existing issues:**
- `td-29s` - Add comprehensive AGENTS.md (this task)
- `td-mqo` - Systematize work tracking from planning docs
- `td-q5f` - Systematic test coverage expansion
- `td-a1b` - Fix broken/skipped tests
- `td-8p0` - Apply Runic/SciML formatting

---

## Troubleshooting

### Common Issues

**"Function not found"**
- Check you're using fully qualified name: `Mycelia.function_name()`
- Check the function is included in `src/Mycelia.jl`

**"Type mismatch with BioSequences"**
- Ensure using `BioSequences.LongDNA{4}`, not `LongDNA`
- For FASTQ: `FASTX.sequence(BioSequences.LongDNA{4}, record)`

**"Tests pass locally but fail in CI"**
- Check if test needs `MYCELIA_RUN_EXTERNAL=true`
- Check for hardcoded paths
- Use `StableRNGs` for determinism

**"External tool not found"**
- Install via Bioconda: `Mycelia.ensure_conda_environment()`
- Or use system install and update PATH

### Getting Help

- Check `planning-docs/` for architectural context
- Run `julia --project=. -e "?Mycelia.function_name"` for docstrings
- See `docs/src/workflow-map.md` for capability matrix
