# Implementation Plan: Reduced Amino Acid Alphabets (Archived)

## Status: ✅ COMPLETE — authoritative status now tracked in `planning-docs/TODO.md`
**Last Updated**: 2026-XX-XX (synced post-rebase)
**Branch**: `main`
**Authoritative implementation**: `src/constants.jl` (`REDUCED_ALPHABET_*` maps) and `reduce_amino_acid_alphabet()`
**Tests**: `test/2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl` covers all shipped schemes; add graph-creation coverage (Unicode + AA graphs) if gaps are found.
**Note**: This document is retained for literature references only; do not use it as a status tracker.

---

## Overview

Adding 30 new reduced amino acid alphabet schemes from validated literature sources. These schemes have been extensively studied and validated for protein structure analysis, fold recognition, sequence alignment, and contact interaction preservation.

> **Current status (2026)**: All schemes in this document are implemented in `src/constants.jl` and exposed via `reduce_amino_acid_alphabet`. The detailed lists below are historical reference only (not work items). If new schemes are added, update `planning-docs/TODO.md` and the tests, not this archived plan.

**Sources:**
- Murphy et al. (2000): 12 alphabets (sizes 2, 3, 4, 5, 8, 10, 12, 15)
- Prlić et al. (2000): 2 alphabets (sizes 12, 17)
- Solis (2015): 18 alphabets (sizes 2-19, comprehensive hierarchical set)
- Wang & Wang: 1 alphabet (size 5)
- Melo & Marti-Renom: 1 alphabet (size 5)

---

## Character Encoding Strategy (DECIDED)

### Principle
**Each reduced alphabet character must be a member of the group it represents, OR use non-amino-acid symbols.**

### Rule
Don't use amino acid letter 'X' to represent a group unless 'X' is a member of that group.

### Examples
- ✅ Group {K,R} represented by 'K' - OK (K is in the group)
- ✅ Group {D,E} represented by '-' - OK (non-amino-acid symbol)
- ❌ Group {A,G,S,T} represented by 'H' - WRONG (H not in group)

### Output Format
Return `String` for all schemes (consistent, allows '+'/'-' chars where needed like in CHEMICAL6)

---

## Verified Citations

1. **Murphy et al. (2000)** - "Simplified amino acid alphabets for protein fold recognition and implications for folding" *Protein Engineering Design & Selection* 13(3):149-152
   - Schemes: MURPHY2, MURPHY3, MU4, ML5, MURPHY8, MURPHY10, MURPHY12, MURPHY15

2. **Prlić et al. (2000)** - "Structure-derived substitution matrices for alignment of distantly related sequences" *Protein Engineering*
   - Schemes: SDM12_PRLIC, HSDM17

3. **Solis (2015)** - "Amino acid alphabet reduction preserves fold information contained in contact interactions in proteins" *Proteins* 83(12):2198-2216. doi: 10.1002/prot.24936. PMID: 26407535
   - Schemes: SOLIS2, SOLIS3, SOLIS4, SOLIS5, SOLIS6, SOLIS7, SOLIS8, SOLIS9, SOLIS10, SOLIS11, SOLIS12, SOLIS13, SOLIS14, SOLIS15, SOLIS16, SOLIS17, SOLIS18, SOLIS19

4. **Wang & Wang (1999)** - "A computational approach to simplifying the protein folding alphabet" *Nature Structural & Molecular Biology* 6:1033-1038
   - Schemes: WW5

5. **Melo & Marti-Renom (2006)** - "Accuracy of sequence alignment and fold assessment using reduced amino acid alphabets" *Proteins* 63(4):986-995. doi: 10.1002/prot.20881. PMID: 16506243
   - Schemes: MM5

---

## Schemes to Implement (30 Total)

### Initially Planned (9 schemes - Status: Constants Added)

### 2-Letter (1 scheme)
**MURPHY2** - Murphy HP model (BLOSUM50-derived)
- 'I': L,V,I,M,C,A,G,S,T,P,F,Y,W (Hydrophobic/Small)
- 'E': E,D,N,Q,K,R,H (Polar/Charged)
- Status: ✅ **CONSTANT ADDED**

### 3-Letter (1 scheme)
**MURPHY3** - Murphy 3-class (BLOSUM50-derived)
- 'L': L,A,S,G,V,T,I,P,M,C
- 'E': E,K,R,D,N,Q,H (Charged/Polar)
- 'F': F,Y,W (Aromatic)
- Status: ✅ **CONSTANT ADDED**

### 4-Letter (1 scheme)
**MU4** - Murphy 4-class (minimal practical fold assessment)
- 'S': A,G,P,S,T (Small/Flexible)
- 'L': C,I,L,M,V (Aliphatic)
- 'E': D,E,H,K,N,Q,R (Charged/Polar)
- 'F': F,Y,W (Aromatic)
- Status: ✅ **CONSTANT ADDED**

### 5-Letter (3 schemes)
**ML5** - Murphy/BLOSUM 5-class (optimal balance, ~90% MI retained)
- 'L': L,V,I,M,C (Aliphatic)
- 'A': A,S,G,T,P (Small)
- 'F': F,Y,W (Aromatic)
- 'E': E,D,N,Q (Acidic/Polar)
- 'K': K,R,H (Basic)
- Status: ✅ **CONSTANT ADDED**

**WW5** - Wang & Wang 5-class
- 'I': C,F,I,L,M,V,W,Y (Hydrophobic)
- 'A': A,H,T
- 'D': D,E (Acidic)
- 'G': G,P
- 'K': K,N,Q,R,S (Polar/Basic)
- Status: ✅ **CONSTANT ADDED**

**MM5** - Melo & Marti-Renom 5-class
- 'A': A,G (Tiny)
- 'C': C (Cysteine)
- 'D': D,E,K,N,P,Q,R,S,T (Large mixed)
- 'I': F,I,L,M,V,W,Y (Hydrophobic)
- 'H': H (Histidine)
- Status: ✅ **CONSTANT ADDED**

### 12-Letter (2 schemes)
**MURPHY12** - Murphy 12-class (98.1% MI retained)
- 'L': L,V,I,M (Aliphatic)
- 'C': C | 'A': A | 'G': G | 'P': P | 'W': W | 'H': H
- 'S': S,T (Hydroxyl)
- 'F': F,Y (Aromatic)
- 'E': E,Q (Glu/Gln)
- 'D': D,N (Asp/Asn)
- 'K': K,R (Basic)
- Status: ✅ **CONSTANT ADDED**

**SDM12_PRLIC** - Prlić structure-derived (best AUC for fold recognition)
- 'A': A | 'D': D | 'N': N | 'C': C | 'W': W | 'H': H | 'G': G | 'P': P
- 'K': K,E,R (Mixed charge cluster)
- 'T': T,S,Q (Hydroxyl/amide)
- 'Y': Y,F (Aromatic)
- 'L': L,I,V,M (Aliphatic)
- Status: ✅ **CONSTANT ADDED**

### 17-Letter (1 scheme)
**HSDM17** - Prlić structure-derived (best mean pooled precision)
- 'A': A | 'D': D | 'R': R | 'N': N | 'T': T | 'S': S | 'Q': Q
- 'Y': Y | 'F': F | 'M': M | 'C': C | 'W': W | 'H': H | 'G': G | 'P': P
- 'K': K,E (Lys/Glu cluster)
- 'L': L,I,V (Aliphatic)
- Status: ✅ **CONSTANT ADDED**

---

### Additional Schemes to Implement (21 schemes - Status: PENDING)

#### Murphy et al. (2000) - 3 Additional Sizes

**MURPHY8** - Murphy 8-class
- 'L': L,V,I,M,C (Aliphatic)
- 'A': A,G (Small)
- 'S': S,T (Hydroxyl)
- 'P': P (Proline)
- 'F': F,Y,W (Aromatic)
- 'E': E,D,N,Q (Acidic/Polar)
- 'K': K,R (Basic)
- 'H': H (Histidine)
- Status: ⏳ **TO BE ADDED**

**MURPHY10** - Murphy 10-class
- 'L': L,V,I,M (Aliphatic)
- 'C': C (Cysteine)
- 'A': A (Alanine)
- 'G': G (Glycine)
- 'S': S,T (Hydroxyl)
- 'P': P (Proline)
- 'F': F,Y,W (Aromatic)
- 'E': E,D,N,Q (Acidic/Polar)
- 'K': K,R (Basic)
- 'H': H (Histidine)
- Status: ⏳ **TO BE ADDED**

**MURPHY15** - Murphy 15-class (high information retention)
- 'L': L,V,I,M (Aliphatic)
- 'C': C | 'A': A | 'G': G | 'S': S | 'T': T | 'P': P
- 'F': F,Y (Aromatic)
- 'W': W (Tryptophan)
- 'E': E | 'D': D | 'N': N | 'Q': Q
- 'K': K,R (Basic)
- 'H': H (Histidine)
- Status: ⏳ **TO BE ADDED**

#### Solis (2015) - Complete Hierarchical Set (18 Alphabets)

**SOLIS2** - Solis 2-class (minimal reduction)
- 'C': C,F,I,L,M,V,W,Y (Hydrophobic)
- 'A': A,D,E,G,H,K,N,P,Q,R,S,T (Hydrophilic)
- Status: ⏳ **TO BE ADDED**

**SOLIS3** - Solis 3-class
- 'C': C,F,I,L,M,V,W,Y (Hydrophobic)
- 'D': D,E,G,K,N,Q,S (Polar)
- 'A': A,H,P,R,T
- Status: ⏳ **TO BE ADDED**

**SOLIS4** - Solis 4-class
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E,G,K,N,Q,S (Polar)
- 'A': A,H,P,R,T
- Status: ⏳ **TO BE ADDED**

**SOLIS5** - Solis 5-class
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E,G,K,N,S (Charged/Polar)
- 'A': A,P,Q,T
- 'H': H,R (Positive)
- Status: ⏳ **TO BE ADDED**

**SOLIS6** - Solis 6-class
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G,K,N,Q,S (Mixed polar)
- 'A': A,P,T
- 'H': H,R (Basic)
- Status: ⏳ **TO BE ADDED**

**SOLIS7** - Solis 7-class
- 'F': F,W,Y (Aromatic)
- 'C': C,I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'K': K (Lysine)
- 'G': G,N,P,Q,S
- 'A': A,T
- 'H': H,R (Basic)
- Status: ⏳ **TO BE ADDED**

**SOLIS8** - Solis 8-class
- 'F': F,W,Y (Aromatic)
- 'I': I,L,M,V (Aliphatic)
- 'C': C (Cysteine)
- 'D': D,E (Acidic)
- 'K': K (Lysine)
- 'G': G,N,P,Q,S
- 'A': A,T
- 'H': H,R (Basic)
- Status: ⏳ **TO BE ADDED**

**SOLIS9** - Solis 9-class
- 'F': F,W,Y (Aromatic)
- 'I': I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G,N,Q,S
- 'P': P,T
- 'A': A (Alanine)
- 'H': H,R (Basic)
- 'C': C (Cysteine)
- 'K': K (Lysine)
- Status: ⏳ **TO BE ADDED**

**SOLIS10** - Solis 10-class
- 'W': W,Y (Large aromatic)
- 'F': F (Phenylalanine)
- 'I': I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G,N,Q,S
- 'P': P,T
- 'A': A (Alanine)
- 'H': H,R (Basic)
- 'C': C (Cysteine)
- 'K': K (Lysine)
- Status: ⏳ **TO BE ADDED**

**SOLIS11** - Solis 11-class
- 'W': W,Y (Large aromatic)
- 'I': I,L,M,V (Aliphatic)
- 'D': D,E (Acidic)
- 'G': G (Glycine)
- 'N': N,P,Q,S
- 'T': T (Threonine)
- 'H': H,R (Basic)
- 'F': F (Phenylalanine)
- 'A': A (Alanine)
- 'C': C (Cysteine)
- 'K': K (Lysine)
- Status: ⏳ **TO BE ADDED**

**SOLIS12** - Solis 12-class
- 'W': W,Y (Large aromatic)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'N': N,P,Q,S
- 'H': H,R (Basic)
- 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

**SOLIS13** - Solis 13-class
- 'W': W,Y (Large aromatic)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'P': P (Proline)
- 'N': N,Q,S
- 'H': H,R (Basic)
- 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

**SOLIS14** - Solis 14-class
- 'W': W | 'Y': Y (Separate aromatics)
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'P': P (Proline)
- 'N': N,Q,S
- 'H': H,R (Basic)
- 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

**SOLIS15** - Solis 15-class
- 'I': I,L (Aliphatic branched)
- 'M': M,V (Aliphatic sulfur/branched)
- 'D': D,E (Acidic)
- 'N': N,Q,S
- 'H': H | 'R': R (Separate basic)
- 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

**SOLIS16** - Solis 16-class
- 'I': I,L (Aliphatic branched)
- 'M': M | 'V': V (Separate sulfur and branched)
- 'D': D,E (Acidic)
- 'N': N,Q,S
- 'H': H | 'R': R | 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

**SOLIS17** - Solis 17-class
- 'I': I | 'L': L | 'M': M | 'V': V (All aliphatic separate)
- 'D': D,E (Acidic)
- 'N': N,Q,S
- 'H': H | 'R': R | 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

**SOLIS18** - Solis 18-class
- 'D': D,E (Acidic - only remaining pair)
- 'N': N (Asparagine)
- 'Q': Q,S (Amide/hydroxyl)
- All others separate: I, L, M, V, H, R, W, Y, P, G, T, F, A, C, K
- Status: ⏳ **TO BE ADDED**

**SOLIS19** - Solis 19-class (nearly complete alphabet)
- 'D': D | 'E': E (Separate acidic)
- 'Q': Q,S (Amide/hydroxyl - final pair)
- 'N': N | 'I': I | 'L': L | 'M': M | 'V': V | 'H': H | 'R': R | 'W': W | 'Y': Y | 'P': P | 'G': G | 'T': T | 'F': F | 'A': A | 'C': C | 'K': K
- Status: ⏳ **TO BE ADDED**

---

## SDM12 Naming Conflict Resolution

### Issue
Current `SDM12` in codebase:
- Groups: D,E together | N,Q together | K,R together | S,T together
- This is actually closest to Murphy12 but with different pairings!

True Murphy12 (from report):
- Groups: E,Q together | D,N together | K,R together | S,T together

True SDM12 (Prlić et al.):
- Groups: K,E,R together | T,S,Q together | D separate | N separate

### Resolution
1. **Keep current SDM12 as-is** - it's a valid legacy scheme
2. **Add MURPHY12** - exact from report
3. **Add SDM12_PRLIC** - true Prlić scheme with `:SDM12_PRLIC` symbol
4. Users can choose which 12-letter scheme fits their needs

---

## Implementation Progress

### ✅ COMPLETED

#### Phase 1a: Add Individual Constants to `src/constants.jl` (9 new)
- ✅ REDUCED_ALPHABET_MURPHY2
- ✅ REDUCED_ALPHABET_MURPHY3
- ✅ REDUCED_ALPHABET_MU4
- ✅ REDUCED_ALPHABET_ML5
- ✅ REDUCED_ALPHABET_WW5
- ✅ REDUCED_ALPHABET_MM5
- ✅ REDUCED_ALPHABET_MURPHY12
- ✅ REDUCED_ALPHABET_SDM12_PRLIC
- ✅ REDUCED_ALPHABET_HSDM17

All constants added with full docstrings and citations.

### 🔄 IN PROGRESS

#### Phase 1b: Update Dictionaries in `src/constants.jl`
- ⏳ Update `REDUCED_ALPHABETS` dictionary to include all 9 new schemes
- ⏳ Update `REDUCED_ALPHABET_INFO` dictionary with metadata for each

### ⏸️ PENDING

#### Phase 2: Update Function Documentation (`src/amino-acid-analysis.jl`)
- Update `reduce_amino_acid_alphabet()` docstring with all new schemes
- Update `list_reduced_alphabets()` expected count (6 → 15)
- Update `get_reduced_alphabet_info()` docstring

#### Phase 3: Add Comprehensive Tests (`test/2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl`)
For each new scheme (9 total):
1. Test all 20 amino acids map correctly
2. Verify reduced characters are group members (no ambiguity)
3. Test with standard sequence "ACDEFGHIKLMNPQRSTVWY"
4. Verify metadata consistency (classes count matches unique chars)
5. Update test count expectations (6 → 15 schemes)

#### Phase 4: Update Tutorial (`tutorials/11_reduced_amino_acid_alphabets.jl`)
1. Add examples for Murphy schemes (MURPHY2, MURPHY3, MU4, ML5, MURPHY12)
2. Add examples for Wang & Wang (WW5) and Melo & Marti-Renom (MM5)
3. Add examples for Prlić schemes (SDM12_PRLIC, HSDM17)
4. Update comparison table with all 15 schemes
5. Add section on choosing alphabets based on:
   - Information retention (2-letter: 75-80%, 5-letter: 90%, 12-letter: 98%)
   - Application (fold recognition, alignment, ML models)
   - Trade-offs (simplicity vs. information loss)

#### Phase 5: Final Verification
1. Run all tests: `julia --project=. -e "import Pkg; Pkg.test()"`
2. Verify documentation builds
3. Check that `list_reduced_alphabets()` returns 15 schemes
4. Verify all metadata is complete and consistent

---

## Testing Strategy

For each of the 9 new schemes, verify:
1. ✅ All 20 standard amino acids are mapped
2. ✅ Reduced characters are members of groups they represent (no ambiguity)
3. ✅ Metadata `:classes` count matches unique output characters
4. ✅ Metadata `:groups` covers all 20 amino acids
5. ✅ Full alphabet test: `"ACDEFGHIKLMNPQRSTVWY"` produces expected output

---

## Expected Final State

**Total schemes**: 36 (6 existing + 30 new)

### Existing (6):
- HP2, HYDROPATHY3, GBMR4, CHEMICAL5, CHEMICAL6, SDM12 (legacy)

### New - Initially Planned (9):
- MURPHY2, MURPHY3, MU4, ML5, WW5, MM5, MURPHY12, SDM12_PRLIC, HSDM17

### New - Additional Murphy (3):
- MURPHY8, MURPHY10, MURPHY15

### New - Solis Hierarchical Set (18):
- SOLIS2, SOLIS3, SOLIS4, SOLIS5, SOLIS6, SOLIS7, SOLIS8, SOLIS9, SOLIS10, SOLIS11, SOLIS12, SOLIS13, SOLIS14, SOLIS15, SOLIS16, SOLIS17, SOLIS18, SOLIS19

---

## Next Steps When Resuming

1. ✅ **DONE**: Add initial 9 constants to `src/constants.jl`
2. **NEXT**: Add remaining 21 constants to `src/constants.jl` (MURPHY8, MURPHY10, MURPHY15, SOLIS2-SOLIS19)
3. **THEN**: Update `REDUCED_ALPHABETS` dictionary with all 30 new schemes (line ~574)
4. **THEN**: Update `REDUCED_ALPHABET_INFO` dictionary with metadata for all 30 schemes (line ~596)
5. **THEN**: Update function docstrings in `src/amino-acid-analysis.jl`
6. **THEN**: Add tests for all 30 new schemes
7. **THEN**: Update tutorial with examples from Murphy, Solis, WW, and MM alphabets
8. **FINALLY**: Run tests and verify all 36 schemes work correctly

---

## Files Modified So Far

1. ✅ `src/constants.jl` - Added initial 9 new alphabet constants (MURPHY2, MURPHY3, MU4, ML5, WW5, MM5, MURPHY12, SDM12_PRLIC, HSDM17)
2. ✅ `planning-docs/IMPLEMENTATION_PLAN_REDUCED_ALPHABETS.md` - Updated to include all 30 new alphabets

## Files Still to Modify

1. ⏳ `src/constants.jl` - Add remaining 21 alphabet constants (MURPHY8, MURPHY10, MURPHY15, SOLIS2-SOLIS19)
2. ⏳ `src/constants.jl` - Update REDUCED_ALPHABETS dict with all 30 schemes (~line 574)
3. ⏳ `src/constants.jl` - Update REDUCED_ALPHABET_INFO dict with metadata for all 30 schemes (~line 596)
4. ⏳ `src/amino-acid-analysis.jl` - Update docstrings to reflect 36 total schemes
5. ⏳ `test/2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl` - Add tests for all 30 new schemes
6. ⏳ `tutorials/11_reduced_amino_acid_alphabets.jl` - Add examples for all new alphabet families

---

## Notes

- Character encoding follows strict "member of group" rule to avoid ambiguity
- All schemes return String (not LongAA) for flexibility and consistency
- Citations verified: Murphy et al. (2000), Prlić et al. (2000), Solis (2015), Wang & Wang (1999), Melo & Marti-Renom (2006)
- Solis (2015) provides a complete hierarchical set from size 2 to 19, offering researchers fine-grained control over alphabet reduction
- Murphy et al. (2000) provides well-validated schemes at multiple information retention levels

---

## Quick Reference: Expected Outputs

Test sequence: `"ACDEFGHIKLMNPQRSTVWY"`

| Scheme | Expected Output | Classes |
|--------|----------------|---------|
| MURPHY2 | `"IIEEIIEIEIIEIEEIIIII"` | 2 |
| MURPHY3 | `"LLEEFLELEELEEELLLLFF"` | 3 |
| MU4 | `"SEEEEFSLEEELFEEESLLF"` | 4 |
| ML5 | `"LEEEEFALEEELKEEEALFF"` | 5 |
| WW5 | `"AIDDDIGIAIKIDDKIKIDI"` | 5 |
| MM5 | `"ADDDDDIADDIDIIDDAAIH"` | 5 |
| MURPHY12 | `"ADDDDGLKLLSPKKKSSLFY"` | 12 |
| SDM12_PRLIC | `"ADDDDGLKLLTNKKKTSLFY"` | 12 |
| HSDM17 | `"ADDDDGLKLLSQKKRTSLYF"` | 17 |

---

**End of implementation plan. All alphabets have been incorporated above.**
