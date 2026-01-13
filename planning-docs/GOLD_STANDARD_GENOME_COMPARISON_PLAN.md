# Gold-Standard Genome Comparison Plan

## Goal
Add an alignment-first whole-genome comparison toolkit that reports ANI/AAI with coverage,
handles rearrangements, and is robust to circular breakpoint differences.

## Current Baseline
- `src/sequence-comparison.jl` already provides BLAST-fragment ANI (`calculate_gold_standard_ani`)
  and BLASTP best-hit AAI (`calculate_gold_standard_aai`).
- `src/alignments-and-mapping.jl` already wraps MUMmer (`run_nucmer`, `run_dnadiff`,
  `parse_dnadiff_report`) and BLAST utilities.

## Implementation Steps
1. **Planning + wiring**
   - Add new high-level functions (ANIm via `dnadiff`, ANIb directional, AAI RBH).
   - Establish deterministic per-pair output directories and caching.
2. **FASTA prep + circular handling**
   - Normalize FASTA (uppercase, whitespace-free; optional record sorting).
   - Implement canonical rotation for single-contig circular genomes
     (lexicographically minimal rotation, plus reverse-complement check).
3. **Gold-standard metrics**
   - ANIm: parse `dnadiff` report for 1-to-1 and M-to-M identity + aligned fractions.
   - ANIb: fragment BLAST both directions; compute ANI + aligned-fragment fraction.
   - AAI: reciprocal best hits (BLASTP or DIAMOND) + orthologous fraction.
4. **Integration wrapper**
   - `compare_genomes_gold` orchestrator with method selection and consolidated output.
5. **Tests + docs**
   - Unit tests for circular canonicalization and report parsing.
   - Optional external-tool tests gated by `MYCELIA_RUN_EXTERNAL=true`.

## Open Questions
- Circular detection heuristic for `circular=:auto` (header keyword vs. explicit flag).
- Should ANIb default fragment size align to 1000 nt, or keep current 1020?
- Preferred default AAI tool (BLASTP vs DIAMOND) for current CI environment.
