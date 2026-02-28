# Type Check Audit Protocol — Rhizomorph Module

## Lesson Learned

When a new variant is added to a family of types, any `isa OldType` check in the
codebase is a potential silent bug if the new type can reach that code path but
is not included in the checked union.

Example failure mode:

```julia
# Before adding UltralightQualityEdgeData:
if edge_data isa LightweightEdgeData
    do_reduced_thing(edge_data)   # works for lightweight
end

# After adding UltralightQualityEdgeData:
# The check silently excludes the new type — wrong branch is taken.
# The fix: check isa AllReducedEdgeData instead.
```

Run this checklist every time a new concrete type is added to any of the union
families below.

---

## Union Type Registry

The canonical definitions live in:

```
src/rhizomorph/core/evidence-functions.jl   (lines ~617–1209)
```

### Tier 1 — Leaf Unions (single kind)

| Const                          | Members                                                                                        |
| ------------------------------ | ---------------------------------------------------------------------------------------------- |
| `LightweightVertexData`        | `LightweightKmerVertexData`, `LightweightBioSequenceVertexData`, `LightweightStringVertexData` |
| `LightweightData`              | `LightweightVertexData`, `LightweightEdgeData`                                                 |
| `UltralightVertexData`         | `UltralightKmerVertexData`, `UltralightBioSequenceVertexData`, `UltralightStringVertexData`    |
| `UltralightData`               | `UltralightVertexData`, `UltralightEdgeData`                                                   |
| `UltralightQualityVertexData`  | `UltralightQualityKmerVertexData`, `UltralightQualityBioSequenceVertexData`                    |
| `UltralightQualityData`        | `UltralightQualityVertexData`, `UltralightQualityEdgeData`                                     |
| `LightweightQualityVertexData` | `LightweightQualityKmerVertexData`, `LightweightQualityBioSequenceVertexData`                  |
| `LightweightQualityData`       | `LightweightQualityVertexData`, `LightweightQualityEdgeData`                                   |

### Tier 2 — Super-Unions (aggregate all reduced types)

| Const                      | Members                                                                                                        |
| -------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `AllReducedVertexData`     | `LightweightVertexData`, `UltralightVertexData`, `UltralightQualityVertexData`, `LightweightQualityVertexData` |
| `AllReducedEdgeData`       | `LightweightEdgeData`, `UltralightEdgeData`, `UltralightQualityEdgeData`, `LightweightQualityEdgeData`         |
| `AllReducedData`           | `AllReducedVertexData`, `AllReducedEdgeData`                                                                   |
| `QualityReducedVertexData` | `UltralightQualityVertexData`, `LightweightQualityVertexData`                                                  |
| `QualityReducedEdgeData`   | `UltralightQualityEdgeData`, `LightweightQualityEdgeData`                                                      |
| `QualityReducedData`       | `QualityReducedVertexData`, `QualityReducedEdgeData`                                                           |

**Design rule:** Code that handles any reduced type should check against a Tier
2 super-union, not a Tier 1 leaf union. Leaf union checks in algorithm code are
almost always bugs waiting to happen.

---

## Audit Checklist — Run on Every New Type Addition

### Step 1: Update the union definition

Add the new concrete type to every leaf union it belongs to in
`evidence-functions.jl`. Then verify the Tier 2 super-unions still cover it
(they are defined in terms of Tier 1 unions, so they usually update
automatically).

```julia
# Verify super-union expansion by checking the definition block around line 1192:
grep -n "AllReducedVertexData\|AllReducedEdgeData\|QualityReducedVertexData" \
    src/rhizomorph/core/evidence-functions.jl
```

### Step 2: Grep for isa checks against leaf union names

These are the highest-risk patterns — a leaf-union isa check excludes any newly
added sibling type.

```bash
# Checks that should almost always be AllReducedEdgeData or AllReducedData instead:
grep -rn "isa LightweightEdgeData\|isa UltralightEdgeData" src/ --include="*.jl"
grep -rn "isa LightweightData\|isa UltralightData" src/ --include="*.jl"

# Checks that should almost always be AllReducedVertexData or AllReducedData instead:
grep -rn "isa LightweightVertexData\|isa UltralightVertexData" src/ --include="*.jl"

# Quality-aware leaf checks (may be legitimate if only quality types are intended,
# but confirm the new type is quality-aware or that exclusion is intentional):
grep -rn "isa LightweightQuality\|isa UltralightQuality" src/ --include="*.jl"
grep -rn "isa QualityReducedVertexData\|isa QualityReducedEdgeData" src/ --include="*.jl"
```

For each hit: ask "can the new type reach this code path?" If yes and it is
excluded, either widen the check to a super-union or add an explicit branch for
the new type.

### Step 3: Grep for isa checks against concrete type names

Individual concrete-type isa checks are the most brittle form. They are
appropriate only inside constructors, conversion functions, and factory dispatch
(e.g., `_build_collapsed_vertex`). They are bugs if they appear in general
algorithm code.

```bash
# Vertex concrete types
grep -rn "isa LightweightKmerVertexData\|isa LightweightBioSequenceVertexData\|isa LightweightStringVertexData" \
    src/ --include="*.jl"
grep -rn "isa UltralightKmerVertexData\|isa UltralightBioSequenceVertexData\|isa UltralightStringVertexData" \
    src/ --include="*.jl"
grep -rn "isa LightweightQualityKmerVertexData\|isa LightweightQualityBioSequenceVertexData" \
    src/ --include="*.jl"
grep -rn "isa UltralightQualityKmerVertexData\|isa UltralightQualityBioSequenceVertexData" \
    src/ --include="*.jl"

# Edge concrete types
grep -rn "isa LightweightEdgeData\b\|isa UltralightEdgeData\b" src/ --include="*.jl"
grep -rn "isa LightweightQualityEdgeData\|isa UltralightQualityEdgeData" src/ --include="*.jl"
```

### Step 4: Grep for direct `.evidence` field access

Reduced types do NOT have an `.evidence` field. Direct field access without a
type guard will throw at runtime.

```bash
grep -rn "\.evidence\b" src/ --include="*.jl" | \
    grep -v "evidence-functions.jl\|vertex-data.jl\|edge-data.jl\|evidence-structures.jl"
```

For each hit outside the core definition files: verify the surrounding code
guards against `AllReducedData` before accessing `.evidence`, or delegates via
`get_dataset_evidence` / `collect_evidence_entries` (which handle reduced types
safely).

### Step 5: Grep for `hasfield` checks that may need updating

`hasfield` guards are a sign that code was patched to handle a specific reduced
type. When a new reduced type is added, these guards may need to include the new
type in their condition.

```bash
grep -rn "hasfield.*dataset_observations\|hasfield.*joint_quality\|hasfield.*evidence" \
    src/ --include="*.jl"
```

Known live site (as of 2026-02-28):

```
src/rhizomorph/algorithms/simplification.jl:568
    if hasfield(typeof(vertex_data), :dataset_observations)
```

This guard exists because `UltralightData` lacks `:dataset_observations` but
`LightweightData` has it. If a new reduced type is added that also lacks
`:dataset_observations`, this guard remains correct. If a new reduced type adds
`:dataset_observations`, verify the merge loop body is still correct for it.

### Step 6: Check factory / construction dispatch

Factory functions that return different concrete types based on a mode flag must
include a branch for the new type.

```bash
grep -rn "UltralightEdgeData\|LightweightEdgeData\|UltralightQualityEdgeData\|LightweightQualityEdgeData" \
    src/rhizomorph/fixed-length/kmer-graphs.jl
```

Known live site (as of 2026-02-28):

```
src/rhizomorph/fixed-length/kmer-graphs.jl:507-513
```

Add a branch for any new edge type. Mirror the pattern for vertex factory sites.

### Step 7: Check `_build_collapsed_vertex` and similar reconstruction functions

These functions reconstruct a vertex of the same type after simplification. They
must include a branch for every concrete reduced type.

```bash
grep -n "_build_collapsed_vertex\|_build_collapsed_edge" src/ -r --include="*.jl"
```

Current site:

```
src/rhizomorph/algorithms/simplification.jl:544-554
```

If the new type is a vertex type, add an `elseif` branch that returns an
instance of the new type.

### Step 8: Run the test suite

```bash
julia --project test/runtests.jl
```

Pay particular attention to:

```
test/4_assembly/
test/in_development/
```

If any test constructs a graph with the new type and runs through
simplification, path-finding, or error-correction, it will exercise the isa
dispatch paths.

---

## Quick One-Liner Audit

### Via TYPE-CHECK-AUDIT markers (preferred)

All live `isa` dispatch sites in algorithm code are tagged with a
`# TYPE-CHECK-AUDIT` comment. Get an instant inventory:

```bash
grep -rn "TYPE-CHECK-AUDIT" src/ --include="*.jl"
```

This is line-number-stable — markers travel with the code, so the output is
always accurate regardless of refactoring.

### Via pattern grep (comprehensive)

Run all grep checks at once from the repo root:

```bash
grep -rn \
    "isa LightweightEdgeData\|isa UltralightEdgeData\|isa LightweightVertexData\|isa UltralightVertexData\|isa LightweightData\b\|isa UltralightData\b\|isa LightweightQuality\|isa UltralightQuality\|hasfield.*dataset_observations\|hasfield.*joint_quality\|hasfield.*evidence\b" \
    src/ --include="*.jl"
```

Zero hits in algorithm files (everything under `src/rhizomorph/algorithms/` and
`src/rhizomorph/analysis/`) is the passing condition. Hits in
`src/rhizomorph/core/evidence-functions.jl` and
`src/rhizomorph/core/graph-type-conversions.jl` are expected (that is where the
dispatch overloads live).

---

## Files Most Likely to Need Updates

When adding a new reduced type, check these files in order:

| File                                                   | What to check                                                                                          |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| `src/rhizomorph/core/evidence-functions.jl`            | Add to correct leaf union; verify Tier 2 covers it; add `add_evidence!`, counting, and query overloads |
| `src/rhizomorph/core/vertex-data.jl` or `edge-data.jl` | Add the struct definition                                                                              |
| `src/rhizomorph/fixed-length/kmer-graphs.jl`           | Add factory branch                                                                                     |
| `src/rhizomorph/algorithms/simplification.jl`          | Add `_build_collapsed_vertex` branch; review `hasfield` guard                                          |
| `src/rhizomorph/algorithms/path-finding.jl`            | Verify `isa AllReducedEdgeData` covers the new type                                                    |
| `src/rhizomorph/algorithms/error-correction.jl`        | Same as path-finding                                                                                   |
| `src/rhizomorph/core/graph-type-conversions.jl`        | Verify conversion logic handles new type                                                               |
| `src/rhizomorph/core/graph-construction.jl`            | Verify construction dispatch                                                                           |

---

## Adding a New Concrete Type — Minimal Checklist

1. `[ ]` Define struct in `vertex-data.jl` or `edge-data.jl`
2. `[ ]` Add to leaf union in `evidence-functions.jl`
3. `[ ]` Verify Tier 2 super-unions still cover it (grep for
   `AllReducedVertexData` / `AllReducedEdgeData`)
4. `[ ]` Add `add_evidence!` overload in `evidence-functions.jl`
5. `[ ]` Add counting / query overloads (`count_total_observations`,
   `count_evidence_entries`, `get_dataset_evidence`, etc.)
6. `[ ]` Add factory branch in `kmer-graphs.jl`
7. `[ ]` Add `_build_collapsed_vertex` / `_build_collapsed_edge` branch in
   `simplification.jl`
8. `[ ]` Run the full grep audit (Step 1–7 above)
9. `[ ]` Run `julia --project test/runtests.jl` — all tests green
