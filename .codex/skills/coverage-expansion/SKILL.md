---
name: coverage-expansion
description: Mycelia-specific test coverage expansion workflow
---

# Mycelia Test Coverage Expansion

Systematically expand test coverage for the Mycelia Julia package.

## Prerequisites

1. Check for existing coverage files:
   ```bash
   find . -name "*.cov" -mtime -1 2>/dev/null
   ```
   If recent .cov files exist, use them to identify gaps without re-running.

2. Environment flags - ALWAYS enable these for accurate coverage:
   ```bash
   export MYCELIA_RUN_ALL=true
   export MYCELIA_RUN_EXTERNAL=true
   ```
   Coverage needs non-precompiled code for accurate line hits. Use `--compiled-modules=no`
   (matches `ci/hpc/run_hpc_ci.sh`).

## Running Coverage Analysis

```bash
# Full coverage with all tests enabled (matches HPC CI flags)
MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true \
  julia --project=. --compiled-modules=no \
  -e "import Pkg; Pkg.test(coverage=true)"

# Coverage files will be generated as src/*.jl.cov
```

Optional parity with HPC CI (skips Codecov upload):
```bash
ci/hpc/run_hpc_ci.sh --tests-only --no-codecov
```

## Analyzing Coverage

Coverage files are generated alongside source files (e.g., `src/utility-functions.jl.cov`).

Format: Each line shows execution count, then source line:
```
    1  function foo(x)
    0      if x < 0      # uncovered branch
    -          return -x
    1      end
    1      return x
    -  end
```

- Number = times executed
- `0` = never executed (needs test)
- `-` = non-executable line

## Test Placement

Match source file to test stage:

| Source File | Test Stage |
|-------------|------------|
| `reference-databases.jl` | `test/1_data_acquisition/` |
| `fastx.jl`, `preprocessing.jl` | `test/2_preprocessing_qc/` |
| `kmer-analysis.jl` | `test/3_feature_extraction_kmer/` |
| `src/rhizomorph/` | `test/4_assembly/` |
| `quality-control-and-benchmarking.jl` | `test/5_validation/` |
| `annotation.jl`, `genome-features.jl` | `test/6_annotation/` |
| `taxonomy-and-trees.jl`, `pangenome.jl` | `test/7_comparative_pangenomics/` |
| `bioconda.jl`, `sentencepiece.jl` | `test/8_tool_integration/` |

## Test Template

```julia
import Test
import StableRNGs
import Mycelia

Test.@testset "FunctionName - description" begin
    rng = StableRNGs.StableRNG(42)

    Test.@testset "normal case" begin
        result = Mycelia.function_name(input)
        Test.@test result == expected
    end

    Test.@testset "edge case - empty input" begin
        result = Mycelia.function_name([])
        Test.@test isempty(result)
    end

    Test.@testset "error case" begin
        Test.@test_throws ArgumentError Mycelia.function_name(invalid)
    end
end
```

## Priority Order

1. **Core utilities** - `utility-functions.jl` (used everywhere)
2. **Data I/O** - `fastx.jl`, `reference-databases.jl`
3. **Assembly core** - `src/rhizomorph/` modules
4. **Analysis** - `kmer-analysis.jl`, `taxonomy-and-trees.jl`
5. **External tools** - `bioconda.jl`, `sentencepiece.jl` (require MYCELIA_RUN_EXTERNAL)

## Rules

- **NEVER disable tests** because functionality is broken - fix the implementation
- Use `Test.@test_skip` ONLY for explicitly planned but unimplemented features
- Always use fully qualified names: `Mycelia.function_name()`, not imported functions
- Use `StableRNGs.StableRNG(seed)` for any randomness
- Use small synthetic data, not large external files
- Commit after each file's tests pass

## Commit Format

```
add: tests for utility-functions.jl coverage

- Add tests for memory_estimate() normal and edge cases
- Add tests for safe_mkdir() error handling
- Coverage: 45% -> 72% for utility-functions.jl
```

## Verifying Progress

After adding tests, re-run coverage and compare:

```bash
# Quick check of specific file coverage
MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true \
  julia --project=. --compiled-modules=no -e '
    import Pkg; Pkg.test()
' && grep -c "^    0" src/utility-functions.jl.cov
```
