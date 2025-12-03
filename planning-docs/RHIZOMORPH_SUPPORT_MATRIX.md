# Rhizomorph Support Matrix (single/double/canonical × graph type × alphabet)

Legend: ✅ supported, 🚫 not applicable, ⏳ pending/partial (documented)

| Graph Type | Alphabet | Singlestrand | Doublestrand | Canonical | Notes |
|------------|----------|--------------|--------------|-----------|-------|
| K-mer | DNA | ✅ | ✅ | ✅ | core k-mer graph-construction, traversal/tests (singlestrand, doublestrand, canonical) |
| K-mer | RNA | ✅ | ✅ | ✅ | traversal/tests in place |
| K-mer | AA  | ✅ | 🚫 | 🚫 | no reverse complement; errors tested |
| Qualmer | DNA | ✅ | ✅ | ✅ | doublestrand/canonical traversal tests added |
| Qualmer | RNA | ✅ | ✅ | ✅ | canonical traversal test added |
| Qualmer | AA  | 🚫 | 🚫 | 🚫 | not applicable |
| N-gram | String | ✅ | 🚫 | 🚫 | non-RC; errors tested |
| Variable-length OLC (FASTA) | DNA/RNA | ✅ | ✅ | ✅ | singlestrand build; conversions implemented in variable-length/strand-conversions.jl |
| Variable-length OLC (FASTA) | AA/String | ✅ | 🚫 | 🚫 | RC not defined; errors enforced |
| Variable-length OLC (FASTQ) | DNA/RNA | ✅ | ✅ | ✅ | quality-aware; conversions implemented |
| Variable-length OLC (FASTQ) | AA/String | 🚫 | 🚫 | 🚫 | not applicable |
| Variable-length OLC (String) | Unicode | ✅ | 🚫 | 🚫 | non-RC; errors enforced |
| Reduced AA alphabets | AA/String inputs | ✅ | 🚫 | 🚫 | covered in k-mer/ngram and OLC tests |

