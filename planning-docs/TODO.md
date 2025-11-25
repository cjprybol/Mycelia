# Mycelia Rhizomorph Implementation TODO

**Created**: 2025-01-25
**Status**: In Progress
**Branch**: rhizomorph-graph-suite

---

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

### Critical Gaps ❌

**Untested Core Algorithms (~50% of implementation)**
- ❌ Path finding (find_eulerian_paths_next) - NO TESTS FOUND
- ❌ Sequence reconstruction (path_to_sequence) - NO TESTS FOUND
- ❌ Simplification (detect_bubbles_next) - NO TESTS FOUND
- ❌ I/O roundtrip (GFA export/import) - UNCLEAR STATUS

**Incomplete Graph Construction Testing**
- ❌ Qualmer graphs - Quality functions tested, full construction unclear
- ❌ RNA/AA graphs - Old API tests only, comprehensive scenarios missing
- ❌ Variable-length OLC graphs - Implementation exists, tests unknown
- ❌ N-gram graphs - Implementation exists, tests unknown
- ❌ Multi-read, multi-dataset scenarios - Limited coverage
- ❌ Edge cases and error handling - Minimal coverage

**Missing Algorithm Implementations**
- ❌ remove_tips() - Not found
- ❌ collapse_linear_chains() - Not found
- ❌ Strand conversion algorithms - Not found (algorithms/strand-conversions.jl planned)
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

## Approach: Rigorous Test-First Validation

**Philosophy**: Test everything, fix failures, then update docs based on tested reality.

**User Requirements**:
- ✅ 100% critical path coverage before marking anything complete
- ✅ Test-first approach (write tests → fix failures → document)
- ✅ Remove old code only after new code proven

**Timeline**: 4 weeks to production-ready status

---

## Phase 1: Critical Algorithm Testing (Week 1) 🔴 IN PROGRESS

**Priority**: HIGHEST - Core assembly functionality cannot be trusted without tests

### 1.1 Path Finding Tests
**File**: `test/4_assembly/path_finding_test.jl` (CREATE)

- [ ] Test Eulerian path detection on simple linear graph
- [ ] Test multiple valid paths (branching structures)
- [ ] Test disconnected components (multiple separate subgraphs)
- [ ] Test circular paths (cycles)
- [ ] Test error cases (no valid Eulerian paths)
- [ ] Test degree validation (even in-degree/out-degree requirements)
- [ ] Test with DNA k-mer graphs (k=3, k=31)
- [ ] Test with RNA k-mer graphs
- [ ] Test with AA k-mer graphs
- [ ] Verify path_vector contains correct vertex labels
- [ ] Verify path_vector order preserves graph connectivity

### 1.2 Sequence Reconstruction Tests
**File**: `test/4_assembly/sequence_reconstruction_test.jl` (CREATE)

- [ ] Test k-mer graph reconstruction (single sequence)
- [ ] Test k-mer graph reconstruction (multiple sequences with overlap)
- [ ] Test qualmer graph reconstruction with quality preservation
- [ ] Test variable-length FASTA graph reconstruction
- [ ] Test variable-length FASTQ graph reconstruction with quality
- [ ] Test string graph reconstruction
- [ ] Test type stability (output type matches input BioSequence type)
- [ ] Test strand orientation handling (Forward vs Reverse)
- [ ] Test reverse complement scenarios (DoubleStrand mode)
- [ ] Verify reconstructed sequence matches original input
- [ ] Verify length correctness (k-mer overlap handling)
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

## Phase 2: Graph Construction Testing (Week 2) 🟡 PLANNED

**Priority**: HIGH - Validate all graph builders work correctly

### 2.1 K-mer Graph Comprehensive Tests
**Files**:
- `test/4_assembly/dna_kmer_singlestrand_test.jl` (UPDATE - use new API)
- `test/4_assembly/dna_kmer_doublestrand_test.jl` (UPDATE - use new API)
- `test/4_assembly/rna_kmer_singlestrand_test.jl` (UPDATE - use new API)
- `test/4_assembly/rna_kmer_doublestrand_test.jl` (UPDATE - use new API)
- `test/4_assembly/aa_kmer_singlestrand_test.jl` (UPDATE - use new API)

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

**Phase 1**: ⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜ 0/10 complete
**Phase 2**: ⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜ 0/10 complete
**Phase 3**: ⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜ 0/10 complete
**Phase 4**: ⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜ 0/10 complete

**Overall**: 0% → Target: 100%

---

## Notes

- This TODO represents work discovered through comprehensive code analysis
- Testing gaps were found in ~50% of the implementation
- Cannot trust "complete" claims without tests
- Test-first approach will reveal hidden issues
- Conservative estimate: 4 weeks to production-ready
- Update this document as work progresses
