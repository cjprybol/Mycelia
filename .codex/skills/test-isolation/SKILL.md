---
name: test-isolation
description:
  Rules for test isolation — preventing global state leakage from @eval
  monkeypatching
---

# Test Isolation Rules

Mycelia's test suite runs all files in a single Julia process via
`include_all_tests`. Global state mutations persist across files.

## @eval Monkeypatching Convention

**NEVER** use `@eval Mycelia` to redefine module functions in test files unless
**both** of the following are true:

1. The file uses a `zz_` prefix (e.g., `zz_simulation_wrappers.jl`) so it runs
   last in alphabetical include order
2. The redefined functions are not exercised by any subsequent test file

### Preferred Alternatives (in order)

1. **Fake CONDA_RUNNER scripts** — write a shell stub to a temp dir, override
   `Mycelia.CONDA_RUNNER` path (see `with_fake_conda_runner` pattern)
2. **Dependency injection** — accept a callable parameter that defaults to the
   real implementation
3. **Environment variable gating** — skip the code path when
   `MYCELIA_RUN_EXTERNAL=false`
4. **Save/restore pattern** — if `@eval` is unavoidable:

```julia
Test.@testset "Guarded monkeypatch" begin
    original = Mycelia.some_function
    try
        @eval Mycelia function some_function(args...)
            # stub
        end
        # ... test code ...
    finally
        @eval Mycelia some_function = $original
    end
end
```

### Why This Matters

`@eval Mycelia` permanently replaces methods in the module's method table for
the entire process lifetime. Files sorted alphabetically after the patching file
will silently get the stub instead of the real function. This was found in ~12
of 34 PRs during a comprehensive audit (2026-04-04), causing:

- Tests passing with stubs instead of real implementations (false coverage)
- Integration tests silently disabled when `MYCELIA_RUN_EXTERNAL=true`
- Non-deterministic failures depending on file include order

### Checklist for PR Authors

- These checks are surfaced automatically by the CI lint job and the repo
  pre-commit hook, so violations should be visible before merge.
- [ ] No `@eval Mycelia` in files without `zz_` prefix
- [ ] If `@eval` is used, a save/restore pattern is in place
- [ ] Prefer fake CONDA_RUNNER scripts over method replacement
- [ ] Test file does not redefine core APIs (`simulate_*`, `add_bioconda_env`,
      `concatenate_fastq_files`)
