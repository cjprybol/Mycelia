# Tool Wrapper Implementation Status (Archived)

**Date**: 2025-01-25 (historical snapshot)
**Current source of truth**: `planning-docs/TODO.md` (assembler/QC status) and `planning-docs/TODO.md` “Today's Priority Actions” for next steps.
**Purpose**: Historical verification of external tool integrations in Mycelia. Retained for provenance; do not rely on counts here.

---

## Executive Summary

**Mycelia has extensive assembler coverage** with 13+ assemblers wrapped and tested, plus quality control tools. However, **binning and classification tools from the metagenomics workflow are missing**.

### Status Overview

| Category | Implemented | Tested | Missing from Workflow |
|----------|-------------|--------|----------------------|
| **Assemblers (Short)** | 5/5 | 4/5 | 0 |
| **Assemblers (Long)** | 6+/6 | 6/6 | 0 |
| **Strain-Aware** | 3/3 | 0/3 | 0 |
| **Quality Control** | 2/2 | 0/2 | 0 |
| **Classification** | 0/4 | 0/4 | 4 |
| **Binning** | 0/7 | 0/7 | 7 |
| **Post-Binning** | 0/2 | 0/2 | 2 |

---

## Detailed Verification by Workflow Stage

### 1. Pre-Assembly: Containment & Coverage Assessment

**Desired Tools** (from user workflow):
- sourmash - assess containment
- sylph - assess relative coverage

**Status**:
- ❌ sourmash: Commented out in `src/bioconda.jl:299`, NOT implemented
- ❌ sylph: NOT FOUND

**Verdict**: ❌ **MISSING** - Need wrappers for both tools

---

### 2. Read Classification

**Desired Tools** (from user workflow):
- minimap (DNA alignment) + mosdepth (coverage analysis)
- metaphlan (marker gene classification)
- metabuli (DNA + AA classification)

**Status**:
- ✅ minimap: **FOUND** in `src/alignments-and-mapping.jl`
  - Functions: `minimap_index()`, `minimap_map()`, `minimap_map_with_index()`, `minimap_map_paired_end_with_index()`
  - **VERIFIED**: Multiple functions, well-integrated
- ❌ mosdepth: NOT FOUND
- ❌ metaphlan: NOT FOUND (metagenomic-classification.jl is entirely commented out)
- ❌ metabuli: NOT FOUND

**Verdict**: ⚠️ **PARTIAL** - minimap works, but missing mosdepth, metaphlan, metabuli

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
- ❌ Taxometer: NOT FOUND
- ❌ TaxVAMB: NOT FOUND
- ❌ VAMB: NOT FOUND
- ❌ MetaBAT2: NOT FOUND
- ❌ MetaCoAG: NOT FOUND
- ❌ GenomeFace: NOT FOUND
- ❌ COMEBin: NOT FOUND

**Verdict**: ❌ **MISSING** - Zero binning tools wrapped (0/7)

---

### 5. Post-Binning Merging

**Desired Tools** (from user workflow):
- dRep (dereplication)
- MAGmax (merging)

**Status**:
- ❌ dRep: NOT FOUND
- ❌ MAGmax: NOT FOUND

**Verdict**: ❌ **MISSING** - Zero post-binning tools wrapped (0/2)

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

### ⚠️ COMMENTED OUT (1 tool)

1. hifiasm-meta - `src/assembly.jl:443` (function exists but commented)

### ❌ MISSING FROM WORKFLOW (15+ tools)

**Classification & Profiling**:
1. sourmash (commented in bioconda.jl)
2. sylph
3. metaphlan
4. metabuli
5. mosdepth

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
- QUAST, BUSCO, strain-aware tools (HyLight, STRONG, Strainy) DO exist

**What's Actually Missing**:
- ❌ ALL binning tools (7 tools): VAMB, MetaBAT2, MetaCoAG, etc.
- ❌ Post-binning tools (2 tools): dRep, MAGmax
- ❌ Most classification tools (4 tools): sourmash, sylph, metaphlan, metabuli
- ❌ PLASS/penguin assembler
- ❌ Variant calling tools: GATK, Freebayes, Clair3, BCFtools
- ❌ Pangenome tools: PGGB, Cactus, vg toolkit

**What Exists But Needs Tests**:
- ⚠️ QUAST, BUSCO (implemented, need test verification)
- ⚠️ Strain-aware tools: HyLight, STRONG, Strainy (tests commented out)

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
