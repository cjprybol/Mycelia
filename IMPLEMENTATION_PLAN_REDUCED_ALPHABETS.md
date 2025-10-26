# Implementation Plan: Reduced Amino Acid Alphabets

## Status: IN PROGRESS
**Last Updated**: 2025-10-25
**Branch**: `reduced-amino-acid-alphabets`

---

## Overview

Adding 9 new reduced amino acid alphabet schemes to match the comprehensive report on validated reduced alphabets from the literature. These schemes have been extensively studied and validated for protein structure analysis, fold recognition, and sequence alignment.

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
   - Schemes: MURPHY2, MURPHY3, MU4, ML5, MURPHY12

2. **Prlić et al. (2000)** - "Structure-derived substitution matrices for alignment of distantly related sequences" *Protein Engineering*
   - Schemes: SDM12_PRLIC, HSDM17

3. **Wang & Wang** - Citation needed for WW5 scheme

4. **Melo & Marti-Renom** - Citation needed for MM5 scheme

5. **Solis (2019)** - "Reduced alphabet of prebiotic amino acids optimally encodes the conformational space of diverse extant protein folds" *BMC Ecology and Evolution*
   - PREBIOTIC10 (**DEFERRED** - no explicit substitution mappings available in paper)

---

## Schemes to Implement (9 Total)

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

**Total schemes**: 15 (6 existing + 9 new)

### Existing (6):
- HP2, HYDROPATHY3, GBMR4, CHEMICAL5, CHEMICAL6, SDM12 (legacy)

### New (9):
- MURPHY2, MURPHY3, MU4, ML5, WW5, MM5, MURPHY12, SDM12_PRLIC, HSDM17

### Deferred:
- PREBIOTIC10 (awaiting explicit substitution mappings)

---

## Next Steps When Resuming

1. ✅ **DONE**: Add all 9 constants to `src/constants.jl`
2. **NEXT**: Update `REDUCED_ALPHABETS` dictionary (line ~574)
3. **THEN**: Update `REDUCED_ALPHABET_INFO` dictionary (line ~596)
4. **THEN**: Update function docstrings in `src/amino-acid-analysis.jl`
5. **THEN**: Add tests for all 9 new schemes
6. **THEN**: Update tutorial
7. **FINALLY**: Run tests and verify

---

## Files Modified So Far

1. ✅ `src/constants.jl` - Added 9 new alphabet constants (lines 255-572)

## Files Still to Modify

1. ⏳ `src/constants.jl` - Update REDUCED_ALPHABETS dict (~line 574)
2. ⏳ `src/constants.jl` - Update REDUCED_ALPHABET_INFO dict (~line 596)
3. ⏳ `src/amino-acid-analysis.jl` - Update docstrings
4. ⏳ `test/2_preprocessing_qc/reduced_amino_acid_alphabets_test.jl` - Add tests
5. ⏳ `tutorials/11_reduced_amino_acid_alphabets.jl` - Add examples

---

## Notes

- Character encoding follows strict "member of group" rule to avoid ambiguity
- All schemes return String (not LongAA) for flexibility and consistency
- Citations verified via WebSearch (Murphy 2000, Prlić 2000, Solis 2019)
- WW5 and MM5 citations still needed (marked in docstrings)
- PREBIOTIC10 deferred because Solis (2019) doesn't provide explicit substitutions

---

## Quick Reference: Expected Outputs

Test sequence: `"ACDEFGHIKLMNPQRSTVWY"`

| Scheme | Expected Output | Classes |
|--------|----------------|---------|
| MURPHY2 | `"IEEEEEIEEEIIEEEIIEI"` | 2 |
| MURPHY3 | `"LEEEEEFLEEELEEELLLF"` | 3 |
| MU4 | `"SEEEEFSLEEELFEEESLLF"` | 4 |
| ML5 | `"LEEEEFALEEELKEEEALFF"` | 5 |
| WW5 | `"AIDDDIGIAIKIDDKIKIDI"` | 5 |
| MM5 | `"ADDDDDIADDIDIIDDAAIH"` | 5 |
| MURPHY12 | `"ADDDDGLKLLSPKKKSSLFY"` | 12 |
| SDM12_PRLIC | `"ADDDDGLKLLTNKKKTSLFY"` | 12 |
| HSDM17 | `"ADDDDGLKLLSQKKRTSLYF"` | 17 |

---

End of implementation plan. Resume from "Next Steps When Resuming" section.
