# Tool Wrapper Implementation Status (Archived)

**Date**: 2025-01-25 (historical snapshot)
**Current source of truth**: `planning-docs/TODO.md` (assembler/QC status) and `planning-docs/TODO.md` “Today's Priority Actions” for next steps.
**Purpose**: Historical verification of external tool integrations in Mycelia. Retained for provenance; do not rely on counts here.

---

## Executive Summary

**Mycelia has extensive assembler coverage** with 13+ assemblers wrapped and tested, plus quality control tools. However, **binning and classification tools from the metagenomics workflow are missing**.

### Status Overview (Updated 2025-12-10)

| Category | Implemented | Tested | Missing from Workflow |
|----------|-------------|--------|----------------------|
| **Assemblers (Short)** | 5/5 | 4/5 | 0 |
| **Assemblers (Long)** | 6+/6 | 6/6 | 0 |
| **Strain-Aware** | 3/3 | **3/3 opt-in** | 0 |
| **Quality Control** | 2/2 | **2/2 opt-in** | 0 |
| **Classification** | 5/5 | 1/5 | 0 |
| **Binning** | 0/7 | 0/7 | 7 |
| **Post-Binning** | 0/2 | 0/2 | 2 |
| **Sequence Encoding** | 1/1 | 1/1 | 0 |

**2025-12-10 Update**:
- Classification tools (sourmash, metaphlan, metabuli) implemented in `src/classification.jl`
- mosdepth confirmed in `src/xam.jl` with test in `test/5_validation/mosdepth_coverage_qc.jl`
- SentencePiece tokenization implemented in `src/sentencepiece.jl` with tests - supports DNA/RNA/AA/text encoding for ML applications

**2025-12-10 Update**:
- Classification tools (sourmash, metaphlan, metabuli) implemented in `src/classification.jl`
- mosdepth confirmed in `src/xam.jl` with test in `test/5_validation/mosdepth_coverage_qc.jl`

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
- ❌ PLASS/penguin: NOT FOUND
- ✅ SKESA: **FOUND** & **TESTED** in `src/assembly.jl:144`
  - Test: `test/4_assembly/third_party_assemblers.jl:83-117`
  - Function: `run_skesa()`

**Bonus Assemblers Also Wrapped**:
- ✅ SPAdes: `run_spades()` (tested)
- ✅ Velvet: `run_velvet()` (tested)
- ⚠️ MetaVelvet: `run_metavelvet()` (implemented but commented out in tests)

**Verdict**: ⚠️ **MOSTLY COMPLETE** - 3/4 from workflow, missing PLASS/penguin

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

### ❌ MISSING FROM WORKFLOW (15+ tools)

**Classification & Profiling**:
1. sourmash (commented in bioconda.jl)
2. metaphlan
3. metabuli
4. mosdepth

**Assembly**:
6. PLASS/penguin

**Binning** (ALL MISSING):
7. Taxometer
8. TaxVAMB
9. VAMB
10. MetaBAT2
11. MetaCoAG
12. GenomeFace
13. COMEBin

**Post-Binning** (ALL MISSING):
14. dRep
15. MAGmax

---

## Corrected Summary

**My Initial Assessment Was Too Harsh**:
- Claimed "ALL tool integration claims are false" → **INCORRECT**
- Actually **22 tools are wrapped**, 13 with comprehensive tests
- QUAST, BUSCO, strain-aware tools (HyLight, STRONG, Strainy) DO exist and now have opt-in smoke tests

**What's Actually Missing**:
- ❌ ALL binning tools (7 tools): VAMB, MetaBAT2, MetaCoAG, etc.
- ❌ Post-binning tools (2 tools): dRep, MAGmax
- ❌ Most classification tools (4 tools): sourmash, sylph, metaphlan, metabuli
- ❌ PLASS/penguin assembler
- ❌ Variant calling tools: GATK, Freebayes, Clair3, BCFtools
- ❌ Pangenome tools: PGGB, Cactus, vg toolkit

**What Exists But Needs Tests**:
- ✅ QUAST, BUSCO: implemented, opt-in extended tests present (sim + phiX download)
- ✅ Strain-aware tools: HyLight, STRONG, Strainy implemented with resource-aware opt-in smoke tests

---

## Recommendations

### Immediate (This Session)

1. **Update VERIFICATION_FINDINGS.md** with corrected tool status
2. **Update TODO.md** to reflect actual gaps:
   - Remove false claims about assemblers/QC tools
   - Add binning tools as missing priority
   - Add classification tools as missing priority

### Short-Term

3. **Enable & Test Existing Tools**:
   - Uncomment hifiasm-meta and test
   - Create tests for QUAST and BUSCO
   - Uncomment/fix tests for HyLight, STRONG, Strainy

4. **Implement Missing High-Priority Tools**:
   - mosdepth (coverage analysis)
   - metaphlan or metabuli (classification)

### Medium-Term

5. **Implement Binning Pipeline**:
   - Priority: VAMB (most popular)
   - Secondary: MetaBAT2, COMEBin
   - Post-binning: dRep (essential for MAG dereplication)

6. **Implement Classification Tools**:
   - sourmash (containment)
   - sylph (relative coverage)

---

## Files Referenced

**Implementation**:
- `src/assembly.jl` - 13 assemblers + strain-aware tools
- `src/alignments-and-mapping.jl` - minimap, diamond, mmseqs
- `src/quality-control-and-benchmarking.jl` - QUAST, BUSCO
- `src/metagenomic-classification.jl` - (entirely commented out)
- `src/bioconda.jl` - sourmash reference (commented)

**Tests**:
- `test/4_assembly/third_party_assemblers.jl` - Comprehensive assembler tests (1299 lines)
- Tests for binning/classification: NOT FOUND
