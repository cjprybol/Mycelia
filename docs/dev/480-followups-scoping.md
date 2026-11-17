# Issue #480 — item 3 scoping (I6 concordance)

Item 2 (CairoMakie) is implemented in this PR. Item 3 is not. It changes what
`ref_kmer_fraction` / `ref_walkability` mean, and the two halves of the
limitation are different sizes of change. Item 1 (codon-frame input) and item 4
(`@info` line, draft PR #528) are out of scope here.

## What the code does today

`gqc_reference_concordance` (`src/gqc.jl`) always RC-folds via
`_gqc_canonical_kmer` and never reads GFA `L` lines.

- Graph k-mers, reference k-mers, and consecutive reference pairs are all
  stored as the lexicographic min of the k-mer and its reverse complement.
- `ref_kmer_fraction` is the fraction of reference k-mer positions whose
  canonical k-mer occurs in any segment sequence.
- `ref_path_fraction` and `ref_walkability` are the same number: the fraction
  of consecutive reference pairs where both canonical k-mers occur in the node
  set. `ref_breaks` is the count of pairs that fail that test.
- There is no `canonical` argument. The only call site,
  `run_gqc` → `gqc_reference_concordance(g, ref_seqs; k=k, circular=circ)`,
  does not pass `effective_canonical`.
- Repeat nulls do the opposite under `:single`. `gqc_repeat_counts` folds only
  when `canonical=true`, and `run_gqc` passes `effective_canonical` from
  `_gqc_resolve_strand`. `:single` defaults to `canonical=false` (an explicit
  `canonical=true` is honored). `:canonical` is always folded. `:double`
  defaults to folded.

So under `strand_mode=:single`, a k-mer and its RC are distinct for the repeat
null and identical for concordance. An RC-only node counts as "present" in
`ref_kmer_fraction` and can hide a strand-specific absence in
`ref_absent_kmers`.

Walkability is node-set co-occurrence, not a path. `GFALink` orientations and
overlaps are already parsed and used for topology, unitig spans, GraphML, and
`edges.tsv`. `_gqc_link_overlaps` is only the k-mer-node guard. Concordance
does not consult them. The walk loop is also linear only: `circular=true`
wraps k-mer fraction, but not the consecutive-pair loop.

Documented in `docs/gqc-v0.1-spec.md` (Concordance limitation), the manifest
field `concordance_note`, and a comment on `gqc_reference_concordance`. The
HPC note prefers `strand_mode=:canonical` for k-mer-node graphs.
`test/5_validation/gqc.jl` only checks that circular concordance does not
throw. Nothing asserts RC folding or link-awareness.

`benchmarking/indel_toy_overcorrection_diagnostic.jl` has columns named
`raw_mean_ref_kmer_fraction`. Those are alignment stats, not this function.

## Decision 1 — RC folding

**A. Leave concordance always folded.** Status quo. The `:single` null and
the concordance score stay inconsistent. No packet numbers change.

**B. Fold concordance with `effective_canonical`.** Add `canonical::Bool=true`
to `gqc_reference_concordance` (default keeps today's direct-call behavior)
and pass `effective_canonical` from `run_gqc`. `:single` packets change
`ref_kmer_fraction`, `ref_absent_kmers`, `ref_path_fraction`,
`ref_walkability`, and `ref_breaks`. Column names and `GQC_SUMMARY_COLUMNS`
order stay put, but the meaning of those cells changes for non-canonical
runs. Spec, manifest note, and a test (graph contains only the RC of a
reference k-mer; fraction is 0 when `canonical=false` and 1 when
`canonical=true`) go with it.

**C. Keep the current columns and append strand-aware ones.** Frozen columns
stay comparable with existing packets. New names need a decision too, and
summary.csv grows.

**Recommendation:** B, if the follow-up is "null and concordance use the same
k-mer identity." The flag already exists; this is one keyword, one call site,
the spec paragraph, and two assertions. Choose C instead if v0.1 packet
numbers must stay comparable under the same column names. Do not silently
change the default of the public `gqc_reference_concordance` call (tests and
any external caller pass only `k` and `circular`).

## Decision 2 — walkability vs GFA `L` lines

This is not a drop-in edit. "Both canonical k-mers occur as sequence content"
and "an L line joins the nodes that contain them" disagree on unitig graphs.

- On the synthetic fixture, segments are ~63 bp and L overlaps are `62M`.
  Consecutive reference k-mers usually sit inside one segment. There is no L
  line between them. Treating walkability as "an L line exists" would mark
  those intra-segment steps as breaks.
- On a k-mer-node graph, each node is a k-mer and L lines (typically overlap
  `k-1`, with `+`/`-` orientation) are the real adjacencies.
- Orientation matters under `:single` and `:double`. `:canonical` collapses
  RC, so a link on the opposite orient may or may not count.
- Circular references would need an explicit last-to-first rule. Today the
  walk does not wrap.

**A. Keep the content-overlap definition.** Leave `ref_walkability` as
documented. No number changes.

**B. Redefine `ref_walkability` / `ref_path_fraction` / `ref_breaks`.** Needs
a written rule before code: intra-segment steps count as walkable; an L line
counts only when its orientations match the strand mode; whether a circular
reference closes the last pair. Existing column values change for any graph
that is not a complete k-mer-node path.

**C. Keep the current columns and add link-aware ones** (for example
`ref_link_walkability` and `ref_link_breaks`, appended to
`GQC_SUMMARY_COLUMNS`). Still blocked on the rule in B.

**Recommendation:** A for the existing columns. If a path score is wanted,
C, after the walk rule above is chosen. Do not implement B or C in a drive-by
change to `gqc_reference_concordance`.

## Affected files (when implemented)

- `src/gqc.jl` — `gqc_reference_concordance`, the `run_gqc` call, manifest
  `concordance_note`. `GQC_SUMMARY_COLUMNS` only if new columns are appended.
- `docs/gqc-v0.1-spec.md` — Concordance limitation (I6).
- `docs/gqc-v0.1-hpc-run.md` — the `:canonical` preference comment.
- `test/5_validation/gqc.jl` — RC present-only-as-reverse-complement, and, if
  a link score is added, one k-mer-node GFA and one unitig GFA.

## Effort

RC alignment (decision 1, option B) is local: the function, the call site, the
spec sentence, and the two assertions. No new columns.

Link-aware walkability is a new helper plus tests on both graph shapes,
orientation, and circular wrap. It waits on the walk rule. Schema edits only
if the new columns in option C are chosen.
