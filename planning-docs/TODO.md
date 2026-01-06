# Mycelia Rhizomorph Implementation TODO

## 🚨 CRITICAL FINDING

**Test-First Approach Validates Concerns**: Initial testing reveals that claiming work "complete" without tests was premature. First comprehensive test run discovered:

1. **API Mismatch**: Rhizomorph functions exist but weren't accessible (no exports, wrong API assumptions)
2. **Real Bugs Found**: 4/42 path finding tests fail, revealing:
   - Bubble/branching graph algorithm issues
   - DoubleStrand mode k-mer connectivity broken
   - Incorrect error handling
3. **Documentation Gap**: Planning docs claimed Phase 1-3 "complete" but ~50% was actually untested

**Conclusion**: The rigorous test-first approach was absolutely necessary. Without it, we would have moved forward with broken implementations.

### 2025-12-10 Bug Fixes Completed

The following issues have been **FIXED** and all 302 assembly tests now pass:

1. **path_to_sequence API mismatch** - FIXED in `src/rhizomorph/algorithms/path-finding.jl`
   - Added overload for `path_to_sequence(path::Vector{T}, graph)` to handle output of `find_eulerian_paths_next`

2. **DoubleStrand evidence merging bug** - FIXED in `src/rhizomorph/core/graph-construction.jl:744-790`
   - `convert_to_doublestrand()` now merges evidence when RC vertex already exists instead of overwriting

3. **Test file fixes**:
   - `rhizomorph_doublestrand_files_test.jl` - Fixed `BioSequences.DNAKmer` to `Kmers.DNAKmer`
   - `rhizomorph_conversion_errors_test.jl` - Fixed `@test_throws Error` to `@test_throws ErrorException`
   - `rhizomorph_qualmer_rc_evidence_test.jl` - Fixed BioSequences/Kmers type usage
   - `dna_kmer_doublestrand_test.jl` and `rna_kmer_doublestrand_test.jl` - Fixed tests expecting canonical behavior to use `build_kmer_graph_canonical`

---

### Recent Rhizomorph migration progress
- Ported legacy suites to Rhizomorph APIs: `basic_graph_tests`, `sequence_graphs_next`, `graph_algorithms_next` (bubble helpers, DFS fallback, tip thresholds), `end_to_end_graph_tests` (exact k-mer content/evidence positions, qualmer edge/vertex quality checks, variable-length overlap evidence, doublestrand/canonical conversion and GFA vertex-count round-trip), `gfa_io_next`, singlestrand canonicalization and orientation suites (`dna|rna|aa_kmer_singlestrand_test`, `doublestrand_canonicalization_test`, `canonicalization_consistency_test`) now build graphs via `Mycelia.Rhizomorph` (no shims).
- `rna_kmer_singlestrand_test` now asserts strand-aware evidence positions and edge evidence using `Mycelia.Rhizomorph` helpers.
- Remaining legacy-heavy suites still calling `build_kmer_graph_next`/legacy conversions: `end_to_end_assembly_tests.jl`, `six_graph_hierarchy_tests.jl`, `comprehensive_*`, `probabilistic_algorithms_next.jl`, and older string/assembly end-to-end suites.
- Implementation gaps blocking full migration: repeat detection/contig/coverage helpers live only in legacy `sequence-graphs-next.jl`; simplification does not delete edges when resolving bubbles; no JLD2 round-trip coverage for Rhizomorph graphs; probabilistic walk utilities are still legacy-only.

## Executive Summary

The Rhizomorph graph ecosystem has **substantial implementation (~7,000 lines) and excellent documentation (94% coverage)**, but has **critical testing gaps (~40% coverage)**. We cannot consider work "complete" without comprehensive tests. This document tracks the work needed to reach production-ready status.

### Current State

| Component | Status | Coverage |
|-----------|--------|----------|
| **Implementation** | 85% Complete | ~7,000 lines of Rhizomorph code |
| **Documentation** | 94% Complete | All user-facing functions documented |
| **Testing** | 40% Complete | ❌ CRITICAL GAP |

### Quality Standard

A feature is ✅ **COMPLETE** only if ALL THREE criteria are met:
1. ✅ Implementation exists with working logic (not stub)
2. ✅ Docstring present with examples
3. ✅ Comprehensive tests pass (including edge cases)

---

## Critical Findings

### What Works Well ✅

**Data Structures (100% tested)**
- All 4 evidence entry types (EvidenceEntry, QualityEvidenceEntry, EdgeEvidenceEntry, EdgeQualityEvidenceEntry)
- All 6 vertex data types (Kmer, Qualmer, BioSequence, QualityBioSequence, String, QualityString)
- All 6 edge data types
- Nested Dict evidence structure with O(1) queries
- Quality score functions (Phred math, aggregation)
- Graph query functions (traversal, statistics)

**Implementation Quality**
- Type-stable design with generic parameters
- Excellent documentation with examples
- Clean modular architecture
- MetaGraphsNext integration

### Critical Gaps ❌ (Updated 2025-12-10)

**Core Algorithms**
- ✅ Path finding (find_eulerian_paths_next) - 43 tests passing (path_finding_test.jl)
- ✅ Sequence reconstruction (path_to_sequence) - Tests passing (doublestrand/canonical traversal tests)
- 🔸 Simplification (detect_bubbles_next) - Smoke coverage in `graph_algorithms_next`; edge removal still unimplemented
- 🔸 I/O roundtrip (GFA export/import) - Covered in `gfa_io_next` and `end_to_end_graph_tests`; JLD2 still missing

**Graph Construction Testing (Updated 2025-12-10)**
- ✅ K-mer graphs (DNA/RNA/AA) - 232 tests passing across singlestrand/doublestrand/canonical
- ✅ Qualmer graphs - Doublestrand and canonical traversal tests passing
- ✅ Multi-read, multi-dataset scenarios - RC evidence handling test passing
- 🔸 Variable-length OLC graphs - Covered in `sequence_graphs_next`/`end_to_end_graph_tests`; traversal coverage still light
- 🔸 N-gram graphs - Covered in `sequence_graphs_next`/`end_to_end_graph_tests`; traversal coverage still light
- ❌ Edge cases and error handling - Minimal coverage

**Missing Algorithm Implementations**
- ✅ remove_tips! exists and is smoke-tested; full simplification still pending
- ❌ collapse_linear_chains() - Not found
- ✅ Strand conversion implemented for fixed-length and variable-length graphs (convert_* in core/graph-construction.jl and variable-length/strand-conversions.jl)
- ❌ Error correction - Not found (algorithms/error-correction.jl planned)
- ❌ Graph type conversions - Not found (core/graph-type-conversions.jl planned)

**Old Code to Deprecate (~7,300 lines)**
- src/graph-core.jl (15 lines)
- src/kmer-graphs.jl (280 lines)
- src/sequence-graphs-next.jl (3,017 lines)
- src/fasta-graphs.jl (465 lines)
- src/fastq-graphs.jl (779 lines)
- src/string-graphs.jl (816 lines)
- src/qualmer-analysis.jl (635 lines)
- src/qualmer-graphs.jl (641 lines)

---

## 🔍 OLD PLANNING DOCUMENT VERIFICATION (2025-01-25)

**See**: `planning-docs/VERIFICATION_FINDINGS.md` for full details

**Key Discoveries from Systematic Verification**:

### ✅ Good News - More Code Than Expected
- Most claimed implementations DO exist with substantial code (900-1,200+ lines each)
- intelligent-assembly.jl (964 lines) EXISTS
- iterative-assembly.jl (1,218 lines) EXISTS and is ACCESSIBLE
- 4 reinforcement learning implementations (1,253+ lines each) EXIST
- cross-validation.jl EXISTS
- Simulation, quality control, and algorithm files EXIST

### 🚨 Critical Issues Discovered
1. **Many implementations are COMMENTED OUT** in src/Mycelia.jl:
   - Line 120: intelligent-assembly.jl (DISABLED despite existing)
   - Line 129: cross-validation.jl (DISABLED)
   - Lines 133-137: ALL 4 reinforcement learning files (DISABLED)

2. **Zero Tests for Most Claimed Features**:
   - No tests for intelligent-assembly.jl (claimed "complete")
   - No tests for iterative-assembly.jl (accessible but untested!)
   - No tests for cross-validation.jl (claimed "89/89 tests passing" - FALSE!)

3. **Tool Wrapper Reality (verified in TOOL_WRAPPER_STATUS.md, 2025-01-25)**:
   - ✅ 13 wrappers implemented **and tested**: megahit, metaspades, skesa, spades, velvet, flye, metaflye, canu, hifiasm, metamdbg, minimap2, diamond, mmseqs
   - ✅ 9 wrappers implemented but **untested**: QUAST, BUSCO, HyLight, STRONG, Strainy, apollo, homopolish, unicycler, metavelvet
- ✅ QUAST/BUSCO wrappers implemented **with opt-in tests** (extended env flags) and default outdir behavior; CI-safe guards in harness; HyLight/STRONG/Strainy smoke tests enabled with resource-aware skips.
- ⚠️ Still lightly validated: apollo, homopolish, unicycler, metavelvet
- ⚠️ hifiasm-meta wrapper exists but is commented out
- ✅ NCBI FCS-GX contamination screening wrapper added (bioconda `ncbi-fcs-gx`); supports GX DB download via manifest URL and clean-genome output.
- ⚠️ Implemented but coverage gaps remain: classification (sourmash, metaphlan, metabuli, mosdepth, sylph), binning/post-binning (Taxometer, TaxVAMB, VAMB, MetaBAT2, COMEBin, dRep, MAGmax, etc.), variant calling (GATK, Freebayes, Clair3, BCFtools), pangenome (PGGB, Cactus, vg toolkit). Track tests/docs/tutorials/benchmarks below.

### 📊 Verification Summary
- **Implementations**: Higher quality than expected (many are substantial, not placeholders)
- **Accessibility**: Lower than claimed (~40% commented out and inaccessible)
- **Testing**: Even accessible code often has zero tests
- **False Claims**: Tool integrations and "89/89 tests" demonstrably false

### 🎯 Immediate Actions Added to Roadmap
1. **Test Accessible Implementations** (before continuing with Rhizomorph):
   - [ ] Create tests for iterative-assembly.jl (1,218 lines, accessible but untested)
   - [ ] Verify viterbi-next.jl has tests
   - [ ] Test simulation functions
   - [ ] Test quality control functions

### 📚 UN Corpus Graph Tutorial + Token/Annotation Graphs (Planned)
- [ ] Draft UN Parallel Corpus tutorial (Literate) comparing character n-gram graphs vs token graphs (SentencePiece).
- [ ] Add token-graph builder for adjacency graphs from token streams (supports String, AA, DNA, RNA).
- [ ] Extend token graphs for DNA/RNA strand modes: singlestrand, doublestrand, canonical (match Rhizomorph conventions).
- [ ] Add annotation-path graph builder for ordered protein annotations (e.g., UniRef chains).
- [ ] Provide graph statistics/topology comparisons for n-gram vs token graphs (leveraging `get_graph_statistics`/TDA).

2. **Test and Uncomment High-Value Code**:
   - [ ] Create tests for intelligent-assembly.jl (964 lines)
   - [ ] If tests pass, uncomment in src/Mycelia.jl:120
   - [ ] Create tests for cross-validation.jl
   - [ ] If tests pass, uncomment in src/Mycelia.jl:129

3. **Documentation Cleanup**:
   - [x] Remove false tool integration claims from old planning docs (updated Sylph/Skani status)
   - [ ] Correct "89/89 tests passing" claim
   - [ ] Add warnings about commented-out code
   - [ ] Keep function coverage audit current via `planning-docs/generate_function_coverage.jl` → `planning-docs/FUNCTION_COVERAGE_AUDIT.md`

4. **Tool Wrapper Validation Updates (done)**:
   - QUAST/BUSCO wrappers covered by opt-in extended tests (simulated + phiX NCBI download) and default to outdirs derived from input FASTA stem.
   - hifiasm-meta un-commented and tested with resource-aware skips.
   - HyLight/STRONG/Strainy smoke tests enabled (small synthetic inputs; skip on OOM/killed).

### Active Verification Notes (moved from archive for visibility)
- `src/development/intelligent-assembly.jl` (~964 lines) exists but is commented out in `src/Mycelia.jl`; no tests exist.
- `src/iterative-assembly.jl` is included in `src/Mycelia.jl` but has no dedicated tests.
- `src/development/cross-validation.jl` is commented out; the "89/89 tests passing" claim is false because no tests exist.
- Four reinforcement learning implementations under `src/development/` are all commented out; associated tests live in `test/in_development/` and are not part of the main suite.
- Tool wrapper status: wrappers exist across classification/binning/variant calling/pangenome, but coverage gaps remain (tests/docs/tutorials/benchmarks). Sylph and Skani are now implemented and tested.

**Conclusion**: Old planning docs overstated completion. Code quality is good, but accessibility and testing are the gaps.

---

## New Rhizomorph Migration Actions (current)
- Avoid shims; port remaining legacy graph tests directly to `Mycelia.Rhizomorph` and update GraphMode usage so old code paths can be retired cleanly.
- Port legacy graph tests (comprehensive/end-to-end/probabilistic/assembly hierarchy) to Rhizomorph APIs and add missing Rhizomorph coverage for variable-length/n-gram traversal and GFA/JLD2 round-trips.
- Remove legacy graph includes from `src/Mycelia.jl` once Rhizomorph tests pass; then delete deprecated graph files (`graph-core.jl`, `kmer-graphs.jl`, `sequence-graphs-next.jl`, `string-graphs.jl`, `qualmer-analysis.jl`, `qualmer-graphs.jl`, `fasta-graphs.jl`, `fastq-graphs.jl`).
- Implement the planned conversion layer (`core/graph-type-conversions.jl`, `algorithms/strand-conversions.jl`) plus simplification tip-removal/linear-chain collapse and expose non-strand-specific merge helpers as public APIs; add focused tests.

## Approach: Rigorous Test-First Validation

**Philosophy**: Test everything, fix failures, then update docs based on tested reality.

**User Requirements**:
- ✅ 100% critical path coverage before marking anything complete
- ✅ Test-first approach (write tests → fix failures → document)
- ✅ Remove old code only after new code proven

**Timeline**: 4 weeks to production-ready status

---

## Today's Priority Actions (Third-Party + Rhizomorph Coverage)

**Third-party assemblers (end-to-end with QC + benchmarks)**
- Re-enable and validate all shipped wrappers: un-comment `run_hifiasm_meta` in `src/assembly.jl` and corresponding tests; make sure HyLight/STRONG/Strainy tests run instead of being commented out.
- Wire wrappers into a repeatable QC pipeline: after each assembler run, auto-run `run_quast` and `run_busco` (where applicable) and capture runtime/memory so benchmarks report both accuracy and efficiency.
- Stabilize the PhiX comparison harness (`benchmarking/phix174_assembler_comparison.jl`): fix macOS failures for SPAdes/SKESA/metaSPAdes, address metaMDBG errors on small genomes, and add a Linux CI-friendly smoke dataset.
- Expand coverage: PLASS/PenguiN wrappers exist; add docs/tutorial/benchmark coverage and surface classification/binning gaps (metaphlan/metabuli, mosdepth, Taxometer/TaxVAMB/VAMB/MetaBAT2/etc.) so the workflow is complete from reads → assembly → QC → binning.

**Rhizomorph graph types (single, double, canonical across 6 graph variants)**
- Verify Singlestrand/Doublestrand/Canonical construction for all BioSequence graph variants (k-mer, qualmer, FASTA, FASTQ) across DNA/RNA/AA; add tests that walk paths and canonicalize to catch strand bugs.
- Patch correctness gaps already visible in code: `build_kmer_graph_from_files` ignores `mode` and always builds singlestrand; N-gram/string graphs document why RC/canonical modes are N/A or add the conversions if required.
- Ensure path-finding/simplification functions operate on canonical graphs and have coverage for all alphabets (DNA/RNA/AA) and quality-aware variants.
- Update documentation to explicitly list the 3×6×alphabet matrix and current support status so we can checkpoint progress if interrupted.
- Make strand-mode interconversion rules explicit (Single → Double → Canonical and reverse), including evidence handling and directed vs undirected storage; see updated section in `planning-docs/rhizomorph-graph-ecosystem-plan.md`.
- Add reduced amino acid alphabet coverage to graph-creation tests (AA graphs and Unicode/string graphs) to validate integration beyond preprocessing.
- Clarify current status of variable-length graphs (FASTA/FASTQ): singlestrand only; doublestrand/canonical conversion still pending.
- Rhizomorph 100% plan (remaining):
  - Fixed-length k-mer/qualmer: add doublestrand traversal/reconstruction tests for DNA/RNA; add canonical traversal tests for DNA/RNA qualmers; AA/string already error on doublestrand/canonical.
- Variable-length FASTA/FASTQ: implement doublestrand/canonical converters for DNA/RNA OR add explicit errors+tests if deferring; currently singlestrand only.
  - N-gram/string: keep singlestrand-only; ensure doc/tests cover non-applicable conversions (errors).
  - Algorithms: add quality-aware traversal edge cases (mixed datasets, RC evidence) and verify path_to_sequence on canonical for DNA/RNA k-mer/qualmer.
- Matrix: add/refresh 3×6×alphabet support matrix (supported vs not-applicable vs pending).
- Implemented: variable-length strand conversions for DNA/RNA (convert_variable_length_to_doublestrand / convert_variable_length_to_canonical) + tests; update matrix/docs accordingly.
- Added support matrix: planning-docs/RHIZOMORPH_SUPPORT_MATRIX.md (✅/🚫/⏳ by graph type/alphabet/strand mode).
- Next up (recommended):
  - Add quality-aware traversal edge cases (mixed datasets, RC evidence) for doublestrand/canonical qualmer graphs.
  - Add perf/scale benchmarks for Rhizomorph builders/traversal (optional; after correctness).
  - Keep BUSCO enabled for benchmarks; retain `--skip-busco` for CI.
- Implement doublestrand/canonical support for variable-length FASTA/FASTQ OLC graphs (or add explicit converters) and add traversal tests once available.

---

## Verification Policy (migrated from OLD_DOCS_VERIFICATION)

- Treat any historical “✅ COMPLETE” claims as unverified until code is located, tests are written, and those tests pass.
- Check accessibility as well as existence: several files exist but are commented out in `src/Mycelia.jl` (intelligent-assembly, cross-validation, RL).
- Prefer reality-checked status in this TODO over older planning docs; retire older claims once they are captured here.

---

## Phase 1: Critical Algorithm Testing (Week 1) ✅ MOSTLY COMPLETE

**Priority**: HIGHEST - Core assembly functionality cannot be trusted without tests

**STATUS UPDATE (2025-12-10)**: Path finding and sequence reconstruction now have comprehensive test coverage. All 302 assembly tests pass. Key bugs fixed: path_to_sequence API mismatch, DoubleStrand evidence merging.

### 1.1 Path Finding Tests ✅ COMPREHENSIVE (43 tests passing)
**File**: `test/4_assembly/path_finding_test.jl` ✅ Complete

**Current coverage (43 assertions):**
- [x] Simple linear DNA k=3
- [x] Two overlapping sequences (k=4)
- [x] Basic cycle smoke test (k=3)
- [x] SNP bubble expecting 0 Eulerian paths
- [x] Disconnected components
- [x] RNA k-mer path smoke test (k=3)
- [x] AA k-mer path smoke test (k=3)
- [x] Large k DNA (k=31)
- [x] DoubleStrand mode with correct evidence merging
- [x] Empty input throws `ArgumentError`

### 1.2 Sequence Reconstruction Tests ✅ COVERED
**Files**: Multiple Rhizomorph test files now cover path_to_sequence

- [x] `rhizomorph_doublestrand_traversal_test.jl` - DNA/RNA k-mer and qualmer reconstruction
- [x] `rhizomorph_canonical_path_test.jl` - Canonical graph reconstruction
- [x] `rhizomorph_qualmer_canonical_traversal_test.jl` - Quality-aware reconstruction

**Remaining to enhance:**
- [ ] Degree validation and multiple valid path enumeration
- [ ] Error handling cases beyond empty input (no-path, invalid graph)

### 1.2b Additional Sequence Reconstruction Tests (from original plan)
**Status**: Partially covered by existing tests, remaining items optional

- [x] Test k-mer graph reconstruction (single sequence) - covered in path_finding_test.jl
- [x] Test k-mer graph reconstruction (multiple sequences with overlap) - covered
- [x] Test qualmer graph reconstruction with quality preservation - rhizomorph_qualmer_canonical_traversal_test.jl
- [ ] Test variable-length FASTA graph reconstruction
- [ ] Test variable-length FASTQ graph reconstruction with quality
- [ ] Test string graph reconstruction
- [x] Test type stability (output type matches input BioSequence type) - verified in tests
- [x] Test strand orientation handling (Forward vs Reverse) - rhizomorph_qualmer_rc_evidence_test.jl
- [x] Test reverse complement scenarios (DoubleStrand mode) - rhizomorph_doublestrand_traversal_test.jl
- [x] Verify reconstructed sequence matches original input - covered
- [x] Verify length correctness (k-mer overlap handling) - path connectivity tests
- [ ] Test error cases (invalid paths, disconnected graphs)

### 1.3 Simplification Algorithm Tests
**File**: `test/4_assembly/simplification_test.jl` (CREATE)

- [ ] Test bubble detection on simple SNP bubble
- [ ] Test bubble detection on multiple parallel paths
- [ ] Test bubble structure correctness (source, paths, target)
- [ ] Test path support calculation
- [ ] Test bubble complexity scoring
- [ ] Test simplify_graph_next with support threshold
- [ ] Test bubble resolution (keep high-support path)
- [ ] Test graph modification correctness after simplification
- [ ] Test with quality-aware graphs (Qualmer)
- [ ] Test edge cases (no bubbles, nested bubbles)
- [ ] Verify edge removal (or document MetaGraphsNext limitation)
- [ ] Test error handling

**Exit Criteria**: All Phase 1 tests pass 100% before moving to Phase 2

---

## Phase 2: Graph Construction Testing (Week 2) ✅ MOSTLY COMPLETE

**Priority**: HIGH - Validate all graph builders work correctly

**STATUS UPDATE (2025-12-10)**: 232 k-mer tests passing across all graph types. Tests use correct APIs.

### 2.1 K-mer Graph Comprehensive Tests ✅ COMPLETE
**Files** (all updated and passing):
- `test/4_assembly/dna_kmer_singlestrand_test.jl` ✅ 26 tests passing
- `test/4_assembly/dna_kmer_doublestrand_test.jl` ✅ 9 tests passing (uses canonical builder for canonical tests)
- `test/4_assembly/rna_kmer_singlestrand_test.jl` ✅ 26 tests passing
- `test/4_assembly/rna_kmer_doublestrand_test.jl` ✅ 7 tests passing (uses canonical builder for canonical tests)
- `test/4_assembly/aa_kmer_singlestrand_test.jl` ✅ 12 tests passing
- `test/4_assembly/aa_kmer_graph_test.jl` ✅ 37 tests passing
- `test/4_assembly/kmer_vertex_data_test.jl` ✅ 46 tests passing
- `test/4_assembly/kmer_edge_data_test.jl` ✅ 35 tests passing

**DNA K-mer Tests**
- [ ] Single read, simple sequence (k=3, k=31, k=101)
- [ ] Multiple reads with overlaps
- [ ] Multiple datasets
- [ ] Same observation contributing multiple k-mers
- [ ] Ambiguous base handling (skip N bases)
- [ ] Empty input error handling
- [ ] Invalid k-mer size error handling
- [ ] Evidence structure validation (nested Dict)
- [ ] Coverage tracking correctness
- [ ] Strand tracking (SingleStrand: observed only, DoubleStrand: both orientations)
- [ ] Edge construction validation (k-1 overlap)
- [ ] Palindromic k-mers in DoubleStrand mode
- [ ] Verify NO canonical conversion in SingleStrand mode
- [ ] Path reconstruction roundtrip

**RNA K-mer Tests**
- [ ] All above tests with RNA sequences
- [ ] U vs T handling

**AA K-mer Tests**
- [ ] All above tests with amino acid sequences
- [ ] 20 amino acid alphabet handling

### 2.2 Qualmer Graph Comprehensive Tests
**Files**:
- `test/4_assembly/dna_qualmer_singlestrand_test.jl` (UPDATE)
- `test/4_assembly/dna_qualmer_doublestrand_test.jl` (UPDATE)
- `test/4_assembly/rna_qualmer_singlestrand_test.jl` (UPDATE)
- `test/4_assembly/aa_qualmer_singlestrand_test.jl` (UPDATE)

- [ ] Quality score extraction from FASTQ records
- [ ] Quality score storage (Phred 0-60, UInt8 type)
- [ ] Multiple observations with different quality scores
- [ ] Quality aggregation (independence assumption, additive in Phred space)
- [ ] Quality scores exceeding 60 (joint observations)
- [ ] Quality preservation through graph construction
- [ ] Quality-aware path selection
- [ ] Mean/min quality calculation
- [ ] Quality filtering
- [ ] All k-mer tests above + quality tracking

### 2.3 Variable-Length OLC Graph Tests
**Files**:
- `test/4_assembly/fasta_graph_test.jl` (CREATE or UPDATE)
- `test/4_assembly/fastq_graph_test.jl` (CREATE or UPDATE)
- `test/4_assembly/string_graph_test.jl` (CREATE or UPDATE)

**FASTA OLC Tests**
- [ ] Two sequences with overlap (odd length requirement)
- [ ] Multiple sequences with various overlaps
- [ ] No overlap scenarios
- [ ] Contained sequences
- [ ] DNA, RNA, AA sequence types
- [ ] Overlap length calculation correctness
- [ ] Edge evidence tracking
- [ ] Variable-length vertex data structure
- [ ] Path reconstruction

**FASTQ OLC Tests**
- [ ] All FASTA tests + quality preservation
- [ ] Quality score aggregation at overlaps
- [ ] Quality-aware overlap scoring

**String Graph Tests**
- [ ] Unicode string overlap detection
- [ ] Variable-length string vertices
- [ ] All OLC tests with string data

### 2.4 N-gram Graph Tests
**File**: `test/4_assembly/ngram_graph_test.jl` (CREATE or UPDATE)

- [ ] Fixed-length n-gram extraction
- [ ] N-1 overlap edge construction
- [ ] String vertex types
- [ ] Evidence tracking
- [ ] Path reconstruction
- [ ] Unicode handling

**Exit Criteria**: All Phase 2 tests pass 100% before moving to Phase 3

---

## Phase 3: Integration & Validation (Week 3) 🟡 PLANNED

**Priority**: HIGH - System integration and validation

### 3.1 I/O Roundtrip Tests
**File**: `test/4_assembly/gfa_io_roundtrip_test.jl` (CREATE or UPDATE)

- [ ] GFA export from k-mer graph
- [ ] GFA import back to k-mer graph
- [ ] Topology preservation (vertices, edges)
- [ ] Metadata preservation (evidence, quality scores)
- [ ] Lossless roundtrip validation
- [ ] Test all 6 graph types (Kmer, Qualmer, BioSequence, QualityBioSequence, String, QualityString)
- [ ] Test both strand modes (SingleStrand, DoubleStrand)
- [ ] Error handling (malformed GFA)

### 3.2 End-to-End Assembly Tests
**File**: `test/4_assembly/end_to_end_assembly_test.jl` (CREATE)

- [ ] Complete workflow: FASTQ → Qualmer graph → Path finding → Sequence reconstruction
- [ ] Multi-dataset assembly
- [ ] Quality-aware assembly with filtering
- [ ] Simplification pipeline
- [ ] Assembly validation metrics
- [ ] Real-world test case (small viral genome)

### 3.3 Pathological Test Suite
**File**: `test/4_assembly/pathological_test_suite.jl` (CREATE)

**Graph Types** (24 combinations):
- 3 sequence types: DNA, RNA, AA
- 2 element types: Fixed (k-mer/n-gram), Variable (OLC)
- 2 quality modes: Non-quality, Quality-aware
- 2 strand modes: SingleStrand, DoubleStrand

**Pathological Cases** (6 standard cases):
1. **Linear path** - Simple sequence, no branches
2. **Bifurcation (bubble)** - SNP creating parallel paths
3. **Tandem repeat** - Cycle in graph
4. **Hairpin/inverted repeat** - Reverse complement testing
5. **Disconnected components** - Multiple unrelated sequences
6. **Quality edge case** - Same k-mer with vastly different quality scores

**Test Matrix**:
- [ ] DNA k-mer, SingleStrand, Non-quality × 6 pathological cases
- [ ] DNA k-mer, DoubleStrand, Non-quality × 6 pathological cases
- [ ] DNA qualmer, SingleStrand, Quality × 6 pathological cases
- [ ] DNA qualmer, DoubleStrand, Quality × 6 pathological cases
- [ ] RNA k-mer, SingleStrand, Non-quality × 6 pathological cases
- [ ] RNA k-mer, DoubleStrand, Non-quality × 6 pathological cases
- [ ] RNA qualmer, SingleStrand, Quality × 6 pathological cases
- [ ] RNA qualmer, DoubleStrand, Quality × 6 pathological cases
- [ ] AA k-mer, SingleStrand, Non-quality × 6 pathological cases
- [ ] AA qualmer, SingleStrand, Quality × 6 pathological cases
- [ ] String n-gram × 6 pathological cases
- [ ] Variable-length graphs × 6 pathological cases

**Total**: 144 test cases (24 graph types × 6 pathological cases)

### 3.4 Performance Benchmarking
**File**: `test/4_assembly/performance_benchmarks.jl` (CREATE)

- [ ] Graph construction time vs sequence count
- [ ] Graph construction memory usage
- [ ] Path finding performance
- [ ] Reconstruction performance
- [ ] Scalability tests (small, medium, large datasets)
- [ ] Comparison: new Rhizomorph vs old implementation
- [ ] Baseline metrics for regression testing

**Exit Criteria**: All Phase 3 tests pass 100% before moving to Phase 4

---

## Phase 3.5: Third-Party Assembler Benchmarking 🟡 IN PROGRESS

**Priority**: MEDIUM - Establish baseline before adding Rhizomorph to comparisons

### 3.5.1 PhiX174 Benchmarking Framework
**File**: `benchmarking/phix174_assembler_comparison.jl` ✅ COMPLETE

**Status**: Core framework implemented, bugs identified and partially fixed

**Completed**:
- [x] PhiX174 reference download using `Mycelia.download_genome_by_accession()`
- [x] Read simulation (ART Illumina 100x, Badread PacBio 50x)
- [x] Platform detection and filtering (macOS vs Linux)
- [x] Comparison framework for 8 assemblers:
  - Short-read: SPAdes, SKESA, MEGAHIT (Linux-only), metaSPAdes
  - Long-read: Flye, Canu, metaFlye, metaMDBG
- [x] Metrics calculation (N50, contig count, total length, runtime)
- [x] FastANI integration for genome comparison
- [x] Plots generation (N50, runtime, identity)
- [x] CSV export

**Issues Fixed**:
- [x] Wrapper return values (Flye, Canu, metaFlye) - changed `assembly=` to `contigs=`
- [x] metaMDBG function signature - added `hifi_reads` parameter handling
- [x] Platform checks - MEGAHIT filtered for Linux-only
- [x] Refactored to use existing Mycelia functions instead of duplicating code

**Known Issues** (deferred to next session):
- [ ] Short-read assemblers failing on macOS (SPAdes, SKESA, metaSPAdes)
- [ ] metaMDBG errors on PhiX174 (likely genome size too small)
- [ ] Some assemblers may have platform-specific issues

**Next Steps**:
- [ ] Debug remaining assembler failures
- [ ] Validate all assemblers complete successfully
- [ ] Expand to larger genomes (bacterial, archaeal)
- [ ] Add Rhizomorph to comparison suite

---

## Phase 4: Cleanup & Documentation (Week 4) 🟡 PLANNED

**Priority**: MEDIUM - Only after all tests pass

### 4.1 Deprecate Old Code
**Files to Remove** (~7,300 lines):
- [ ] Add deprecation warnings to old functions
- [ ] Update Mycelia.jl to stop including old files
- [ ] Remove src/graph-core.jl
- [ ] Remove src/kmer-graphs.jl
- [ ] Remove src/sequence-graphs-next.jl
- [ ] Remove src/fasta-graphs.jl
- [ ] Remove src/fastq-graphs.jl
- [ ] Remove src/string-graphs.jl
- [ ] Remove src/qualmer-analysis.jl
- [ ] Remove src/qualmer-graphs.jl
- [ ] Verify no imports or dependencies remain

### 4.2 Update Planning Documents
- [ ] Create ACTIVE_ROADMAP.md with only remaining work
- [ ] Update COMPREHENSIVE_ROADMAP.md with reality-checked status
- [ ] Create TESTING_STATUS.md documenting test coverage
- [ ] Archive or consolidate rhizomorph-graph-ecosystem-plan.md (4,055 lines)
- [ ] Keep comprehensive-testing-framework.md as reference
- [ ] Remove or archive TEST-TEMPLATE.md

### 4.3 Create Migration Guide
**File**: `planning-docs/MIGRATION_GUIDE.md` (CREATE)

- [ ] Old API → New API function mapping
- [ ] Breaking changes documentation
- [ ] Example code updates
- [ ] Deprecation timeline
- [ ] FAQ for common issues

### 4.4 Documentation Updates
- [ ] "Getting Started" tutorial (5 minutes to first assembly)
- [ ] Update all examples to use new Rhizomorph API
- [ ] API documentation verification
- [ ] Tutorial notebook updates

**Exit Criteria**: All old code removed, docs updated, 100% test passage maintained

---

## Future Work (Post-Production) 🔵 BACKLOG

### Missing Algorithm Implementations
- [ ] Implement remove_tips() in simplification.jl
- [ ] Implement collapse_linear_chains() in simplification.jl
- [ ] Implement algorithms/strand-conversions.jl
- [ ] Implement algorithms/error-correction.jl
- [ ] Implement core/graph-type-conversions.jl

### Advanced Features (from planning documents)
- [ ] Reinforcement learning framework
- [ ] Intelligent self-optimizing assembly
- [ ] Iterative maximum likelihood assembly
- [ ] Hybrid assembly methods (OLC + de Bruijn)
- [ ] Cross-validation pipeline
- [ ] Real-time assembly dashboard
- [ ] Interactive visualization tools

### Performance Optimization
- [ ] Multi-threading for graph construction
- [ ] Memory-efficient streaming for large datasets
- [ ] On-disk graph representations
- [ ] Parallel processing for algorithms

### Extended Validation
- [ ] Biological dataset validation (viral, bacterial, eukaryotic)
- [ ] Comparison with gold-standard assemblers
- [ ] Publication-quality benchmarking
- [ ] Manuscript preparation

---

## Success Metrics (Non-Negotiable)

- ✅ 100% of critical path functions have comprehensive tests
- ✅ 100% test passage (zero failures)
- ✅ All 24 graph type combinations tested against pathological suite
- ✅ All test files use new Rhizomorph API (not old build_kmer_graph_next)
- ✅ No placeholders or incomplete implementations in production code
- ✅ Old code removed with clear migration path
- ✅ Planning documents reflect actual tested state
- ✅ >90% code coverage overall

---

## Progress Tracking

**Phase 1**: ✅✅✅✅✅✅✅⬜⬜⬜ 7/10 complete (path finding, sequence reconstruction, evidence merging - all working)
**Phase 2**: ✅✅✅✅✅⬜⬜⬜⬜⬜ 5/10 complete (k-mer graphs DNA/RNA/AA, qualmer graphs)
**Phase 3**: ⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜ 0/10 complete
**Phase 4**: ⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜ 0/10 complete

**Overall**: ~30% → Target: 100%

**2025-12-10 Update**: 302 assembly tests now passing. Key bugs fixed.

---

## Verification Backlog (migrated from archived docs)

- **Commented-out but implemented (needs tests before enabling)**:
  - `src/development/intelligent-assembly.jl` (~964 lines) – uncomment in `src/Mycelia.jl` after tests.
  - `src/development/cross-validation.jl` – “89/89 tests passing” claim was false; no tests exist.
  - Four RL implementations under `src/development/` – all commented out; tests live only in `test/in_development/`.
- **Accessible but untested**:
  - `src/iterative-assembly.jl`
  - `src/viterbi-next.jl` (confirm dedicated tests exist)
  - `src/simulation.jl` functions: `simulate_pacbio_reads`, `simulate_nanopore_reads`, `simulate_illumina_paired_reads`
  - Quality/QC/visualization functions: `plot_per_base_quality`, `analyze_fastq_quality`, `calculate_gc_content`, `plot_embeddings`, `plot_optimal_cluster_assessment_results`, `plot_taxa_abundances`, `visualize_many_timeseries`
- **Algorithm verification (exists; test status unknown)**:
  - `src/coverage-clustering.jl` (k-medoids coverage clustering)
  - `src/kmer-saturation-analysis.jl` (saturation curves and thresholds)
  - `src/graph-cleanup.jl` (statistical cleanup)
  - `src/development/genomic-graph-algorithms.jl`
  - `src/development/pangenome-core-genome.jl`
- **Tutorial/docs checks**:
  - Run `tutorials/00_assembly_in_5_minutes.jl` end-to-end.
  - Inspect `src/development/intelligent-assembly.jl` for the k=5 loop issue noted historically.
  - Verify any “probabilistic-assembly-hub” doc page if referenced elsewhere.
- **Tool wrappers**:
  - Authoritative inventory: `planning-docs/TOOL_WRAPPER_STATUS.md` (13 tested, 9 untested, 1 commented out). Coverage gaps remain in classification/binning/post-binning/variant calling/pangenome; add tests/docs/tutorials/benchmarks and decide on enabling hifiasm-meta.

### Tool Wrapper Coverage TODOs (implemented)

**Classification/Profiling**:
- **sourmash**
  - [ ] Tests: expand `test/8_tool_integration/classification_tools.jl` with opt-in end-to-end assertions (sketch/search/gather outputs).
  - [ ] Docs: add usage snippet + output description in `docs/src/metagenomic-workflow.md`.
  - [ ] Tutorials: add a small profiling example (tool integration or taxonomy tutorial).
  - [ ] Benchmarks: add a small profiling benchmark under `benchmarking/`.
- **metaphlan**
  - [ ] Tests: stabilize opt-in runs in `test/8_tool_integration/classification_tools.jl` and `test/8_tool_integration/metabuli_metaphlan_strainphlan.jl` (document required DB env).
  - [ ] Docs: add usage snippet + output example in `docs/src/metagenomic-workflow.md`.
  - [ ] Tutorials: add a short profiling example (tool integration or taxonomy tutorial).
  - [ ] Benchmarks: add profiling benchmark under `benchmarking/`.
- **metabuli**
  - [ ] Tests: stabilize opt-in runs in `test/8_tool_integration/classification_tools.jl` and `test/8_tool_integration/metabuli_metaphlan_strainphlan.jl` (document required DB env).
  - [ ] Docs: add usage snippet + output example in `docs/src/metagenomic-workflow.md`.
  - [ ] Tutorials: add a short profiling example (tool integration or taxonomy tutorial).
  - [ ] Benchmarks: add profiling benchmark under `benchmarking/`.
- **mosdepth**
  - [ ] Tests: wire `test/5_validation/mosdepth_coverage_qc.jl` into opt-in external runs and add assertions for parsed metrics.
  - [ ] Docs: add a coverage-QC usage snippet in `docs/src/metagenomic-workflow.md`.
  - [ ] Tutorials: expand `tutorials/05_assembly_validation.jl` to explain mosdepth outputs.
  - [ ] Benchmarks: add a coverage-QC benchmark under `benchmarking/`.
- **sylph**
  - [ ] Tests: add opt-in integration assertions in `test/7_comparative_pangenomics/sequence_comparison.jl`.
  - [ ] Docs: add usage snippet in `docs/src/metagenomic-workflow.md`.
  - [ ] Tutorials: add a small profiling example (tool integration or taxonomy tutorial).
  - [ ] Benchmarks: add profiling benchmark under `benchmarking/`.

**Assembly (protein/nucleotide)**:
- **PLASS**
  - [ ] Tests: add stronger output assertions in `test/4_assembly/third_party_assemblers_plass_penguin.jl` and a smaller fixture for CI.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md`.
  - [ ] Tutorials: add a short PLASS example in tool integration or assembly tutorials.
  - [ ] Benchmarks: keep `benchmarking/phix174_assembler_comparison.jl` and add a lighter benchmark path.
- **PenguiN**
  - [ ] Tests: add stronger output assertions in `test/4_assembly/third_party_assemblers_plass_penguin.jl` and a smaller fixture for CI.
  - [ ] Docs: fix `docs/src/metagenomic-workflow.md` wording (wrapper exists) and add usage details.
  - [ ] Tutorials: add a short PenguiN example in tool integration or assembly tutorials.
  - [ ] Benchmarks: keep `benchmarking/phix174_assembler_comparison.jl` and add a lighter benchmark path.

**Binning**:
- **Taxometer**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering Taxometer.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.
- **TaxVAMB**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering TaxVAMB.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.
- **VAMB**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering VAMB.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.
- **MetaBAT2**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering MetaBAT2.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.
- **MetaCoAG**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering MetaCoAG.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.
- **GenomeFace**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering GenomeFace.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.
- **COMEBin**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: add wrapper details to `docs/src/metagenomic-workflow.md` and workflow map.
  - [ ] Tutorials: add a binning tutorial covering COMEBin.
  - [ ] Benchmarks: add a binning benchmark under `benchmarking/`.

**Post-Binning**:
- **dRep**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: expand `docs/src/metagenomic-workflow.md` with dereplication outputs.
  - [ ] Tutorials: add a post-binning example covering dRep.
  - [ ] Benchmarks: add a post-binning benchmark under `benchmarking/`.
- **MAGmax**
  - [ ] Tests: add opt-in external runs to `test/8_tool_integration/binning_tools.jl`.
  - [ ] Docs: expand `docs/src/metagenomic-workflow.md` with MAGmax outputs.
  - [ ] Tutorials: add a post-binning example covering MAGmax.
  - [ ] Benchmarks: add a post-binning benchmark under `benchmarking/`.

**Variant Calling**:
- **GATK**
  - [ ] Tests: add basic wrapper tests in `test/7_comparative_pangenomics/` or a new variant-calling suite.
  - [ ] Docs: update `docs/src/workflow-map.md` to mark wrapper availability and add usage details.
  - [ ] Tutorials: add a variant-calling tutorial (normalize VCF -> evaluate accuracy).
  - [ ] Benchmarks: add a variant-calling benchmark under `benchmarking/`.
- **Freebayes**
  - [ ] Tests: add basic wrapper tests in a variant-calling suite.
  - [ ] Docs: add usage details to workflow map or API workflows.
  - [ ] Tutorials: include in variant-calling tutorial alongside GATK.
  - [ ] Benchmarks: add a variant-calling benchmark under `benchmarking/`.
- **Clair3**
  - [ ] Tests: add basic wrapper tests in a variant-calling suite.
  - [ ] Docs: add usage details to workflow map or API workflows.
  - [ ] Tutorials: include in variant-calling tutorial alongside GATK.
  - [ ] Benchmarks: add a variant-calling benchmark under `benchmarking/`.
- **BCFtools**
  - [ ] Tests: add basic wrapper tests in a variant-calling suite.
  - [ ] Docs: add usage details to workflow map or API workflows.
  - [ ] Tutorials: include in variant-calling tutorial (normalization + consensus).
  - [ ] Benchmarks: add a variant-calling benchmark under `benchmarking/`.

**Pangenome**:
- **PGGB**
  - [ ] Tests: add opt-in end-to-end runs in `test/7_comparative_pangenomics/pangenome_wrappers.jl`.
  - [ ] Docs: add usage details to `docs/src/workflow-map.md` and API workflows.
  - [ ] Tutorials: add a pangenome tutorial showcasing PGGB.
  - [ ] Benchmarks: add a pangenome benchmark under `benchmarking/`.
- **Cactus**
  - [ ] Tests: add opt-in end-to-end runs in `test/7_comparative_pangenomics/pangenome_wrappers.jl`.
  - [ ] Docs: add usage details to `docs/src/workflow-map.md` and API workflows.
  - [ ] Tutorials: add a pangenome tutorial showcasing Cactus.
  - [ ] Benchmarks: add a pangenome benchmark under `benchmarking/`.
- **vg toolkit**
  - [ ] Tests: add opt-in end-to-end runs in `test/7_comparative_pangenomics/pangenome_wrappers.jl` beyond vg deconstruct validation.
  - [ ] Docs: add wrapper details to `docs/src/workflow-map.md` (beyond `docs/src/related-projects.md`).
  - [ ] Tutorials: add a pangenome tutorial showcasing vg conversions and indexing.
  - [ ] Benchmarks: add a pangenome benchmark under `benchmarking/`.

### Planned Tool Integrations
- **Graph/Topology**:
   - [ ] [`prokrustean` (de Bruijn graph topology)](https://github.com/KoslickiLab/prokrustean) — wrap or re-implement to explore k-length graph structure across target k values.
- **Classification/Profiling**:
   - [ ] [`mOTUs` (short/long-read profiling)](https://github.com/motu-tool/mOTUs) — supports long reads via `motus prep_long -i input_long_reads.fastq -o output_converted.fastq`.
- **Assembly**: 
   - [ ] [`shasta` (eukaryotic nanopore)](https://github.com/paoloshasta/shasta)
   - [ ] [`verkko` (T2T/diploid)](https://github.com/marbl/verkko)
   - [ ] [`myloasm` (long-read)](https://github.com/bluenote-1577/myloasm)
   - [ ] [`floria` (strain-aware)](https://github.com/bluenote-1577/floria)
- **Alignment/Variant Calling**:
   - [ ] [`strobealign` (fast mapping)](https://github.com/ksahlin/strobealign)
   - [ ] [`rtg-tools` (vcfeval)](https://github.com/RealTimeGenomics/rtg-tools)
   - [ ] [`gatk`](https://github.com/broadgsa/gatk/)
   - [ ] [`freebayes`](https://github.com/ekg/freebayes)
- **Genomics/Binning**:
   - [ ] [`galah` (genome selection)](https://github.com/wwood/galah)
   - [ ] [`dashing2` (sketching)](https://github.com/dnbaker/dashing2)
- **NCBI**:
   - [ ] [`ncbi-genome-download` (robust fetch)](https://github.com/kblin/ncbi-genome-download)
   - [ ] [`egapx` (eukaryotic genome annotation pipeline)](https://github.com/ncbi/egapx)
     - Protein datasets suitable for most vertebrates, arthropods, echinoderms, cnidaria, and some plants:
       - Chordata - Mammalia, Sauropsida, Actinopterygii (ray-finned fishes), other Vertebrates
       - Insecta - Hymenoptera, Diptera, Lepidoptera, Coleoptera, Hemiptera
       - Arthropoda - Arachnida, other Arthropoda
       - Echinodermata
       - Cnidaria
       - Monocots - Liliopsida
       - Eudicots - Asterids, Rosids, Fabids, Caryophyllales
     - Warning: Fungi, protists, and nematodes are out-of-scope for EGAPx; use a different annotation method.

---

## Notes

- This TODO represents work discovered through comprehensive code analysis
- Testing gaps were found in ~50% of the implementation
- Cannot trust "complete" claims without tests
- Test-first approach will reveal hidden issues
- Conservative estimate: 4 weeks to production-ready
- Update this document as work progresses
