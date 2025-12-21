# Tool Wrapper Implementation Status (Archived)

**Date**: 2025-01-25 (historical snapshot)
**Current source of truth**: `planning-docs/TODO.md` (assembler/QC status) and `planning-docs/TODO.md` “Today's Priority Actions” for next steps.
**Purpose**: Historical verification of external tool integrations in Mycelia. Retained for provenance; do not rely on counts here.

---

## Executive Summary

**Mycelia has extensive assembler coverage** with 13+ assemblers wrapped and tested, plus quality control tools. Coverage gaps remain for tests/docs/tutorials/benchmarks in classification, binning, variant calling, and pangenome tooling.

### Status Overview (Updated 2025-12-10)

| Category | Implemented | Tested | Missing from Workflow |
|----------|-------------|--------|----------------------|
| **Assemblers (Short)** | 5/5 | 4/5 | 0 |
| **Assemblers (Long)** | 6+/6 | 6/6 | 0 |
| **Strain-Aware** | 3/3 | **3/3 opt-in** | 0 |
| **Quality Control** | 2/2 | **2/2 opt-in** | 0 |
| **Classification** | 5/5 | 5/5 (opt-in/validation) | 0 |
| **Binning** | 7/7 | 7/7 (validation-only) | 0 |
| **Post-Binning** | 2/2 | 2/2 (validation-only) | 0 |
| **Sequence Encoding** | 1/1 | 1/1 | 0 |

**Coverage Audit Update**:
- Classification tools (sourmash, metaphlan, metabuli) implemented in `src/classification.jl`; sylph in `src/sequence-comparison.jl`.
- Binning/post-binning tools implemented in `src/binning.jl` with validation-only tests in `test/8_tool_integration/binning_tools.jl`.
- mosdepth confirmed in `src/xam.jl` with opt-in tests in `test/5_validation/mosdepth_coverage_qc.jl`.
- SentencePiece tokenization implemented in `src/sentencepiece.jl` with tests - supports DNA/RNA/AA/text encoding for ML applications.

---

## Detailed Verification by Workflow Stage

### 1. Pre-Assembly: Containment & Coverage Assessment

**Desired Tools** (from user workflow):
- sourmash - assess containment
- sylph - assess relative coverage

**Status**:
- ✅ sourmash: **IMPLEMENTED** in `src/classification.jl` (2025-12-10)
  - Functions: `run_sourmash_sketch()`, `run_sourmash_search()`, `run_sourmash_gather()`
- ✅ sylph: Implemented and tested (`run_sylph_profile` in `src/sequence-comparison.jl`, integration test in `test/7_comparative_pangenomics/sequence_comparison.jl`)

**Verdict**: ✅ **COMPLETE** - Both sourmash and sylph implemented

---

### 2. Read Classification

**Desired Tools** (from user workflow):
- minimap (DNA alignment) + mosdepth (coverage analysis)
- metaphlan (marker gene classification)
- metabuli (DNA + AA classification)
- skani (reference ANI)

**Status**:
- ✅ minimap: **FOUND** in `src/alignments-and-mapping.jl`
  - Functions: `minimap_index()`, `minimap_map()`, `minimap_map_with_index()`, `minimap_map_paired_end_with_index()`
  - **VERIFIED**: Multiple functions, well-integrated
- ✅ mosdepth: **FOUND** in `src/xam.jl:1001`
  - Functions: `run_mosdepth()`, `parse_mosdepth_distribution()`, `parse_mosdepth_summary()`, `parse_mosdepth_thresholds()`, `summarize_mosdepth_qc()`
  - **VERIFIED**: Comprehensive coverage analysis with QC metrics
- ✅ metaphlan: **IMPLEMENTED** in `src/classification.jl` (2025-12-10)
  - Functions: `run_metaphlan()`, `parse_metaphlan_profile()`
- ✅ metabuli: **IMPLEMENTED** in `src/classification.jl` (2025-12-10)
  - Functions: `run_metabuli_classify()`, `run_metabuli_build_db()`, `parse_metabuli_report()`, `parse_metabuli_classifications()`
- ✅ skani: Implemented and tested (`skani_triangle`/`skani_dist` in `src/sequence-comparison.jl`, test in `test/7_comparative_pangenomics/sequence_comparison.jl`)

**Verdict**: ✅ **COMPLETE** - All read classification tools implemented (minimap, mosdepth, metaphlan, metabuli, skani)

---

### 3. Assembly

#### 3a. Short Read Assemblers

**Desired Tools**:
- megahit
- metaspades
- PLASS/penguin
- SKESA

**Status**:
- ✅ megahit: **FOUND** & **TESTED** in `src/assembly.jl:24`
  - Test: `test/4_assembly/third_party_assemblers.jl:255-289`
  - Function: `run_megahit()`
- ✅ metaspades: **FOUND** & **TESTED** in `src/assembly.jl:69`
  - Test: `test/4_assembly/third_party_assemblers.jl:291-328`
  - Function: `run_metaspades()`
- ✅ PLASS/penguin: **FOUND** in `src/assembly.jl` (`run_plass_assemble`, `run_penguin_nuclassemble`, `run_penguin_guided_nuclassemble`)
  - Test: `test/4_assembly/third_party_assemblers_plass_penguin.jl`
- ✅ SKESA: **FOUND** & **TESTED** in `src/assembly.jl:144`
  - Test: `test/4_assembly/third_party_assemblers.jl:83-117`
  - Function: `run_skesa()`

**Bonus Assemblers Also Wrapped**:
- ✅ SPAdes: `run_spades()` (tested)
- ✅ Velvet: `run_velvet()` (tested)
- ⚠️ MetaVelvet: `run_metavelvet()` (implemented but commented out in tests)

**Verdict**: ✅ **COMPLETE** - 4/4 from workflow implemented (PLASS/PenguiN tests are opt-in/resource-aware)

#### 3b. Long Read Assemblers

**Desired Tools**:
- metaflye
- hifiasm-meta

**Status**:
- ✅ metaflye: **FOUND** & **TESTED** in `src/assembly.jl:268`
  - Test: `test/4_assembly/third_party_assemblers.jl:699-732`
  - Function: `run_metaflye()`
- ⚠️ hifiasm-meta: **WRAPPER EXISTS** but **COMMENTED OUT** in `src/assembly.jl:443`
  - Function: `run_hifiasm_meta()` (line 443 is commented)
  - Test: Commented out in third_party_assemblers.jl:637-696

**Bonus Long Read Assemblers Also Wrapped**:
- ✅ flye: `run_flye()` (tested)
- ✅ canu: `run_canu()` (tested)
- ✅ hifiasm: `run_hifiasm()` (tested)
- ✅ metamdbg: `run_metamdbg()` (tested for both HiFi and ONT)
- ✅ unicycler: `run_unicycler()` (hybrid assembler)

**Verdict**: ⚠️ **MOSTLY COMPLETE** - metaflye works, hifiasm-meta exists but disabled

---

### 4. Binning

**Desired Tools** (from user workflow):
- Taxometer + TaxVAMB
- VAMB
- MetaBAT2
- MetaCoAG
- GenomeFace
- COMEBin

**Status**:
- ✅ Taxometer: **FOUND** in `src/binning.jl`
  - Function: `run_taxometer`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`
- ✅ TaxVAMB: **FOUND** in `src/binning.jl`
  - Function: `run_taxvamb`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`
- ✅ VAMB: **FOUND** in `src/binning.jl`
  - Function: `run_vamb`
  - Tests: ✅ Parser/validation in `test/8_tool_integration/binning_tools.jl`
- ✅ MetaBAT2: **FOUND** in `src/binning.jl`
  - Function: `run_metabat2`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`
- ✅ MetaCoAG: **FOUND** in `src/binning.jl`
  - Function: `run_metacoag`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`
- ✅ GenomeFace: **FOUND** in `src/binning.jl`
  - Function: `run_genomeface`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`
- ✅ COMEBin: **FOUND** in `src/binning.jl`
  - Function: `run_comebin`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`

**Verdict**: ✅ **IMPLEMENTED** - 7/7 binning tools wrapped (tests are validation-only; external runs opt-in).

---

### 5. Post-Binning Merging

**Desired Tools** (from user workflow):
- dRep (dereplication)
- MAGmax (merging)

**Status**:
- ✅ dRep: **FOUND** in `src/binning.jl`
  - Function: `run_drep_dereplicate`
  - Tests: ✅ Validation/parsing in `test/8_tool_integration/binning_tools.jl`
- ✅ MAGmax: **FOUND** in `src/binning.jl`
  - Function: `run_magmax_merge`
  - Tests: ✅ Validation in `test/8_tool_integration/binning_tools.jl`

**Verdict**: ✅ **IMPLEMENTED** - 2/2 post-binning tools wrapped (tests are validation-only; external runs opt-in).

---

### 6. Quality Control & Validation

**Tools**:
- QUAST (assembly quality assessment)
- BUSCO (gene completeness)

**Status**:
- ✅ QUAST: **FOUND** in `src/quality-control-and-benchmarking.jl:1038`
  - Functions: `run_quast(assembly_files::Vector{String})`, `run_quast(assembly_file::String)`
  - Tests: ❓ NOT VERIFIED
- ✅ BUSCO: **FOUND** in `src/quality-control-and-benchmarking.jl:1165`
  - Functions: `run_busco(assembly_files::Vector{String})`, `run_busco(assembly_file::String)`
  - Tests: ❓ NOT VERIFIED

**Verdict**: ✅ **IMPLEMENTED** but tests needed

---

### 7. Strain-Aware Tools
  - Tests: ✅ Opt-in extended smoke in `test/5_validation/quast_busco_wrappers_test.jl` (simulated + phiX download) with default outdir derivation
- ✅ BUSCO: **FOUND** in `src/quality-control-and-benchmarking.jl:1165`
  - Functions: `run_busco(assembly_files::Vector{String})`, `run_busco(assembly_file::String)`; supports auto-lineage and dataset predownload helper
  - Tests: ✅ Opt-in extended smoke in `test/5_validation/quast_busco_wrappers_test.jl` (auto-lineage, default outdir)

**Verdict**: ✅ **IMPLEMENTED/TESTED (opt-in extended)**

---

### 7. Sequence Encoding/Tokenization

**Tools**:
- SentencePiece (subword tokenization for ML/NLP applications)

**Status**:
- ✅ SentencePiece: **IMPLEMENTED & TESTED** in `src/sentencepiece.jl` (2025-12-10)
  - Functions: `train_sentencepiece_model()`, `train_sentencepiece_model_from_sequences()`, `train_sentencepiece_model_from_fasta()`, `train_sentencepiece_model_from_fastq()`, `encode_sentencepiece()`, `decode_sentencepiece()`, `load_sentencepiece_model()`, `get_sentencepiece_vocab()`, `sentencepiece_vocab_size()`, `is_valid_sentencepiece_model()`
  - Auto-installs via pip in conda environment on first use
  - Supports DNA, RNA, AA sequences, ASCII, and Unicode text
  - Supports subword regularization/sampling for neural network training
  - Tests: `test/8_tool_integration/sentencepiece.jl`

**Verdict**: ✅ **COMPLETE** - SentencePiece implemented with comprehensive wrapper

---

### 8. Strain-Aware Tools

**Tools** (from old planning docs):
- HyLight
- STRONG
- Strainy

**Status**:
- ✅ HyLight: **FOUND** in `src/assembly.jl:2286`
  - Function: `run_hylight()`
  - Tests: Commented out in third_party_assemblers.jl
- ✅ STRONG: **FOUND** in `src/assembly.jl:2323`
  - Function: `run_strong()`
  - Tests: Commented out in third_party_assemblers.jl
- ✅ Strainy: **FOUND** in `src/assembly.jl:2767`
  - Function: `run_strainy()`
  - Tests: Commented out in third_party_assemblers.jl

**Verdict**: ✅ **IMPLEMENTED** but tests are commented out

---

## Original "False Claims" Verification - CORRECTED

### Previously Claimed False, Actually TRUE:

**Assembly Validation Tools**:
- ✅ QUAST: **EXISTS** (falsely claimed missing) - `src/quality-control-and-benchmarking.jl:1038`
- ✅ BUSCO: **EXISTS** (falsely claimed missing) - `src/quality-control-and-benchmarking.jl:1165`
- ❌ MUMmer: Still NOT FOUND (claim was correct)

**Strain-Aware Tools**:
- ✅ HyLight: **EXISTS** (falsely claimed missing) - `src/assembly.jl:2286`
- ✅ STRONG: **EXISTS** (falsely claimed missing) - `src/assembly.jl:2323`
- ✅ Strainy: **EXISTS** (falsely claimed missing) - `src/assembly.jl:2767`

**Long-Read Assemblers**:
- ✅ metaFlye: **EXISTS** & **TESTED** (falsely claimed missing) - `src/assembly.jl:268`
- ⚠️ hifiasm-meta: **EXISTS** but commented out (partially correct) - `src/assembly.jl:443`
- ✅ SKESA: **EXISTS** & **TESTED** (falsely claimed missing) - `src/assembly.jl:144`
- ❌ IDBA-UD: Test exists but commented out, wrapper unclear

### Still Actually Missing (Claims Were Correct):

**Variant Calling Tools**:
- ❌ GATK: NOT FOUND (claim was correct)
- ❌ Freebayes: NOT FOUND (claim was correct)
- ❌ Clair3: NOT FOUND (claim was correct)
- ❌ BCFtools: NOT FOUND (claim was correct)

**Pangenome Tools**:
- ❌ PGGB: NOT FOUND (claim was correct)
- ❌ Cactus: NOT FOUND (claim was correct)
- ❌ vg toolkit: NOT FOUND (claim was correct)

---

## Comprehensive Tool Inventory

### ✅ IMPLEMENTED & TESTED (13 tools)

**Short Read Assemblers**:
1. megahit - `src/assembly.jl:24` + `test/4_assembly/third_party_assemblers.jl:255`
2. metaspades - `src/assembly.jl:69` + `test/4_assembly/third_party_assemblers.jl:291`
3. skesa - `src/assembly.jl:144` + `test/4_assembly/third_party_assemblers.jl:83`
4. spades - `src/assembly.jl:107` + `test/4_assembly/third_party_assemblers.jl:44`
5. velvet - `src/assembly.jl:2486` + `test/4_assembly/third_party_assemblers.jl:195`

**Long Read Assemblers**:
6. flye - `src/assembly.jl:228` + `test/4_assembly/third_party_assemblers.jl:413`
7. metaflye - `src/assembly.jl:268` + `test/4_assembly/third_party_assemblers.jl:699`
8. canu - `src/assembly.jl:316` + `test/4_assembly/third_party_assemblers.jl:475`
9. hifiasm - `src/assembly.jl:355` + `test/4_assembly/third_party_assemblers.jl:530`
10. metamdbg - `src/assembly.jl:2669` + `test/4_assembly/third_party_assemblers.jl:734,772`

**Mapping**:
11. minimap2 - `src/alignments-and-mapping.jl:497` (multiple functions)
12. diamond - `src/alignments-and-mapping.jl:22`
13. mmseqs - `src/alignments-and-mapping.jl:203`

### ✅ IMPLEMENTED, UNTESTED (9 tools)

1. quast - `src/quality-control-and-benchmarking.jl:1038`
2. busco - `src/quality-control-and-benchmarking.jl:1165`
3. hylight - `src/assembly.jl:2286`
4. strong - `src/assembly.jl:2323`
5. strainy - `src/assembly.jl:2767`
6. apollo - `src/assembly.jl:2189`
7. homopolish - `src/assembly.jl:2236`
8. unicycler - `src/assembly.jl:443`
9. metavelvet - `src/assembly.jl:2540`

### ✅ IMPLEMENTED & TESTED - Sequence Encoding (1 tool)

1. sentencepiece - `src/sentencepiece.jl` + `test/8_tool_integration/sentencepiece.jl`

### ⚠️ COMMENTED OUT (1 tool)

1. hifiasm-meta - `src/assembly.jl:443` (function exists but commented)

### Implemented (coverage gaps remain)

**Classification & Profiling**:
- sourmash — Tests: `test/8_tool_integration/classification_tools.jl` (opt-in external). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: none. Benchmarks: none. TODO: add tutorial + benchmark + usage snippet.
- metaphlan — Tests: `test/8_tool_integration/classification_tools.jl`, `test/8_tool_integration/metabuli_metaphlan_strainphlan.jl` (opt-in external). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: none. Benchmarks: none. TODO: stabilize opt-in runs; add tutorial + benchmark + output examples.
- metabuli — Tests: `test/8_tool_integration/classification_tools.jl`, `test/8_tool_integration/metabuli_metaphlan_strainphlan.jl` (opt-in external). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: none. Benchmarks: none. TODO: stabilize opt-in runs; add tutorial + benchmark + output examples.
- mosdepth — Tests: `test/5_validation/mosdepth_coverage_qc.jl` (opt-in external). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: `tutorials/05_assembly_validation.jl`. Benchmarks: none. TODO: add benchmark + ensure opt-in test wiring.
- sylph — Tests: `test/7_comparative_pangenomics/sequence_comparison.jl` (opt-in external). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: none. Benchmarks: none. TODO: add tutorial + benchmark.

**Assembly (protein/nucleotide)**:
- PLASS — Tests: `test/4_assembly/third_party_assemblers_plass_penguin.jl` (skips on resource constraints). Docs: none. Tutorials: none. Benchmarks: `benchmarking/phix174_assembler_comparison.jl`. TODO: add docs + tutorial; add lighter benchmark case.
- PenguiN — Tests: `test/4_assembly/third_party_assemblers_plass_penguin.jl` (skips on resource constraints). Docs: `docs/src/metagenomic-workflow.md` (currently claims wrapper missing; fix wording). Tutorials: none. Benchmarks: `benchmarking/phix174_assembler_comparison.jl`. TODO: fix doc text; add tutorial.

**Binning**:
- Taxometer — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.
- TaxVAMB — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.
- VAMB — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.
- MetaBAT2 — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.
- MetaCoAG — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.
- GenomeFace — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.
- COMEBin — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: none. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + docs/tutorial/benchmark.

**Post-Binning**:
- dRep — Tests: `test/8_tool_integration/binning_tools.jl` (parser + validation). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + tutorial + benchmark.
- MAGmax — Tests: `test/8_tool_integration/binning_tools.jl` (validation-only). Docs: `docs/src/metagenomic-workflow.md`. Tutorials: none. Benchmarks: none. TODO: add opt-in external test + tutorial + benchmark.

**Variant Calling**:
- GATK — Tests: none. Docs: `docs/src/workflow-map.md` (planned). Tutorials: none. Benchmarks: none. TODO: add tests + docs + tutorial + benchmark.
- Freebayes — Tests: none. Docs: none. Tutorials: none. Benchmarks: none. TODO: add tests + docs + tutorial + benchmark.
- Clair3 — Tests: none. Docs: none. Tutorials: none. Benchmarks: none. TODO: add tests + docs + tutorial + benchmark.
- BCFtools — Tests: none. Docs: none. Tutorials: none. Benchmarks: none. TODO: add tests + docs + tutorial + benchmark.

**Pangenome**:
- PGGB — Tests: `test/7_comparative_pangenomics/pangenome_wrappers.jl` (validation-only). Docs: `docs/src/workflow-map.md`, `docs/src/api/function-coverage.md`. Tutorials: none. Benchmarks: none. TODO: add end-to-end test + tutorial + benchmark.
- Cactus — Tests: `test/7_comparative_pangenomics/pangenome_wrappers.jl` (validation-only). Docs: `docs/src/workflow-map.md`, `docs/src/api/function-coverage.md`. Tutorials: none. Benchmarks: none. TODO: add end-to-end test + tutorial + benchmark.
- vg toolkit — Tests: `test/7_comparative_pangenomics/pangenome_wrappers.jl` (vg deconstruct validation). Docs: `docs/src/related-projects.md`. Tutorials: none. Benchmarks: none. TODO: add docs in workflow map + tutorial + benchmark + broader tests.

### Planned Integrations (not yet wrapped)
- EGAPx (NCBI eukaryotic genome annotation pipeline): planned wrapper; readme datasets cover vertebrates, arthropods (insecta/arachnida), echinoderms, cnidaria, monocots (Liliopsida), eudicots (Asterids/Rosids/Fabids/Caryophyllales). Warning: fungi, protists, and nematodes are out-of-scope.

---

## Corrected Summary

**My Initial Assessment Was Too Harsh**:
- Claimed "ALL tool integration claims are false" → **INCORRECT**
- Actually **22 tools are wrapped**, 13 with comprehensive tests
- QUAST, BUSCO, strain-aware tools (HyLight, STRONG, Strainy) DO exist and now have opt-in smoke tests

**What's Actually Missing**:
- No longer missing implementations for classification, binning, variant calling, or pangenome wrappers listed above.
- Missing coverage remains in tests, docs, tutorials, and benchmarks for most tools; track these gaps in `planning-docs/TODO.md`.

**What Exists But Needs Tests**:
- ✅ QUAST, BUSCO: implemented, opt-in extended tests present (sim + phiX download)
- ✅ Strain-aware tools: HyLight, STRONG, Strainy implemented with resource-aware opt-in smoke tests

---

## Recommendations

### Immediate (This Session)

1. **Update TODO.md** with coverage gaps (tests/docs/tutorials/benchmarks) for implemented wrappers.
2. **Fix documentation mismatches** (e.g., PenguiN marked missing in `docs/src/metagenomic-workflow.md`).

### Short-Term

3. **Expand tests beyond validation-only**:
   - Add opt-in external runs for binning/post-binning tools.
   - Add end-to-end runs for pggb/cactus/vg with small fixtures.
   - Add basic test coverage for variant calling wrappers (gatk/freebayes/clair3/bcftools).

4. **Document and teach**:
   - Add tutorial coverage for classification, binning, variant calling, and pangenome workflows.
   - Add benchmark harnesses (classification profiling, binning, pangenome, variant calling).

### Medium-Term

5. **Raise coverage quality**:
   - Convert validation-only tests into opt-in integration tests with real tool runs.
   - Add doc/API references for all wrappers in the workflow map and API workflow pages.

---

## Files Referenced

**Implementation**:
- `src/assembly.jl` - PLASS/PenguiN wrappers
- `src/binning.jl` - binning + post-binning wrappers
- `src/classification.jl` - sourmash/metaphlan/metabuli wrappers
- `src/metagenomic-classification.jl` - additional metaphlan/metabuli helpers
- `src/pangenome-analysis.jl` - pggb/cactus/vg wrappers
- `src/sequence-comparison.jl` - sylph profiling wrapper
- `src/variant-analysis.jl` - gatk/freebayes/clair3/bcftools wrappers
- `src/xam.jl` - mosdepth wrapper

**Tests**:
- `test/4_assembly/third_party_assemblers_plass_penguin.jl`
- `test/5_validation/mosdepth_coverage_qc.jl`
- `test/7_comparative_pangenomics/pangenome_wrappers.jl`
- `test/7_comparative_pangenomics/sequence_comparison.jl`
- `test/8_tool_integration/binning_tools.jl`
- `test/8_tool_integration/classification_tools.jl`
- `test/8_tool_integration/metabuli_metaphlan_strainphlan.jl`

**Docs/Tutorials/Benchmarks**:
- `docs/src/metagenomic-workflow.md`
- `docs/src/workflow-map.md`
- `docs/src/api/function-coverage.md`
- `docs/src/related-projects.md`
- `tutorials/05_assembly_validation.jl`
- `benchmarking/phix174_assembler_comparison.jl`
