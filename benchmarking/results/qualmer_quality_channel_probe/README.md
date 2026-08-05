# Is the qualmer decoder arm's quality channel inert?

Findings for bead td-4e19d.2. Every claim below rests on a run executed in this
worktree or on code read and cited by file and line. Nothing is inferred from a
function name — that caution is load-bearing here, because several functions in
this area are named for quality and contain no quality term, and at least one is
named for evidence counts and does consume Phred.

## Short answer

| Question                                                             | Answer                                                                              |
| -------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| Is the quality channel inert in the DEFAULT greedy qualmer assembly? | **Yes**, proven in 18/18 probe cells across three chemistries                       |
| Are the Track-A `qualmer` and `kmer` arms the same experiment?       | **Yes** on every assembly metric; they differ only in runtime and memory            |
| Is the Viterbi/DP machinery also quality-blind?                      | **No.** It demonstrably consumes Phred, and Phred changes what it decides           |
| Is the Stage-1 corrector quality-blind end-to-end?                   | **Not determined by this work.** It was never reached — see the gating result below |
| Is the greedy behaviour a defect or a design choice?                 | **A design choice**, with a documented rationale in PR #335                         |
| What changed on this branch?                                         | An **opt-in** `traversal_weighting = :quality`; the default is untouched            |

## 1. The observation reproduces

`~/workspace/rhizomorph-paper/benchmarking/track-a-pilot/track_a_pilot_results_20260724.tsv`
contains 132 rows = 66 paired runs (illumina 24, pacbio 24, ont 18). Re-checking
every pair independently: **zero disagreements** on `n_reads`, `n_contigs`,
`NGA50`, `misassemblies`, `genome_fraction`, `duplication_ratio`,
`largest_contig`, `status`.

Only cost differs, and consistently — the qualmer arm is roughly 1.4-1.9x slower
(median wall seconds, never pooled across chemistry):

| chemistry | qualmer  | kmer     |
| --------- | -------- | -------- |
| illumina  | 194.04 s | 109.86 s |
| pacbio    | 307.91 s | 165.79 s |
| ont       | 424.78 s | 312.36 s |

Note this table is a _description of the committed pilot_, not independent
evidence about the code. Two arms agreeing could equally mean a harness bug.
Section 2 settles that.

## 2. The channel is inert — positive control

`benchmarking/qualmer_quality_channel_probe.jl` holds the read **sequences**
fixed and varies **only** the per-base quality vector, then compares assemblies
by SHA-256 over the strand-canonicalised, sorted contig set. Any difference is
then attributable to quality alone.

Four conditions per cell:

- `oracle` — quality that is _correct by construction_: bases the simulator
  corrupted get a low Q, uncorrupted bases a high Q.
- `constant40` / `constant02` — every base Q40 / Q2, ignoring the truth.
- `antioracle` — the oracle inverted: corrupted bases get HIGH Q, clean bases
  LOW Q.

`oracle` is strictly _more_ informative than any real instrument's Phred
estimate and `antioracle` is maximally adversarial, so this is decisive in the
positive direction: an assembler that consumes quality at all cannot return the
same contigs under all four.

**Result: 18/18 cells byte-identical across all four conditions** — 3
chemistries x coverage {5, 10, 30} x 2 seeds, corrupted-base fractions
0.17%-4.6%. Low coverage and high error rate, where quality should matter most,
are included and behave no differently. The k-mer arm (same reads, quality
stripped to FASTA) also reproduces the identical contig set in every cell.

Data: `qualmer_quality_channel_probe_verdicts.tsv`,
`qualmer_quality_channel_probe_runs.tsv`.

### Why this is not a broken comparator

Every cell also runs a **sensitivity control**: a single base in a single read
is mutated, quality left at oracle. That digest must change, and it **did in
18/18 cells**. Without this the "identical" verdicts would be unfalsifiable — a
comparator that always returns "same" would produce exactly the table above.

### Why the outputs are identical, in code

The mechanism is stronger than a weak scoring function; the traversal has **no
scoring hook at all**:

- Contigs come from `find_eulerian_paths_next`
  (`src/rhizomorph/algorithms/path-finding.jl:2055`) with a fallback to
  `find_contigs_next` (`src/rhizomorph/algorithms/contigs.jl:65`), selected in
  `_qualmer_graph_to_assembly`. Both are purely topological —
  `grep -c quality src/rhizomorph/algorithms/contigs.jl` returns **0** — and
  `find_contigs_next` is the same routine the k-mer arm uses. Identical k-mer
  sets therefore give identical contigs.
- The only quality read in the whole greedy path is
  `_qualmer_vertex_quality_scores` at `src/rhizomorph/assembly.jl:4057`, inside
  contig **emission**, building the output FASTQ quality string. It cannot
  influence which path was taken because the path already exists by then. The
  Track-A harness then writes contigs as FASTA, discarding even that.
- `_qualmer_vertex_score` / `_qualmer_edge_score`
  (`src/rhizomorph/assembly.jl:43,47`) are named for quality but return
  `count_evidence`, a pure count of evidence entries
  (`src/rhizomorph/core/evidence-functions.jl:393`). They are also **only
  reachable from the fallback walk** (`src/rhizomorph/assembly.jl:2461`), which
  runs only when the primary extraction returns nothing — so on the common path
  they are not even the operative code. The original brief for this bead traced
  the root cause to these two functions; that trace is correct about them being
  count-based but understates the problem, because fixing them alone would not
  have changed a single benchmark row.

## 3. Scope: which paths DO consume Phred

This is the part that matters most for manuscript accuracy, and the answer is
not uniform. `benchmarking/qualmer_corrector_quality_sensitivity.jl` measures it
in three stages, because the stages disagree.

### Stage A — the machinery reads Phred (unit level)

| function                                   | Q40     | Q2       | consumes quality |
| ------------------------------------------ | ------- | -------- | ---------------- |
| `default_viterbi_emission_logp` (match)    | -0.001  | -9.968   | yes              |
| `default_viterbi_emission_logp` (mismatch) | -10.310 | -10.531  | yes              |
| `calculate_read_likelihood`                | -0.474  | -391.439 | yes              |

Mechanism: `_viterbi_position_error_rate` (`src/viterbi-next.jl:847`) converts
each position's Phred to `10^(-Q/10)`, and `default_viterbi_emission_logp`
(`src/viterbi-next.jl:10`) scores every emitted position with it. The
corrector's accept/reject likelihood accumulates Phred-derived log-probabilities
in `_rhizomorph_joint_probability` (`src/iterative-assembly.jl:6190`). This is
real, not nominal.

### Stage B — Phred changes what the decoder DECIDES

Running the read decoder (`find_optimal_sequence_path`) over a whole read set
under each condition and counting corrected reads that differ from the `oracle`
condition:

| chemistry | constant40 | constant02 | antioracle | n_reads |
| --------- | ---------- | ---------- | ---------- | ------- |
| illumina  | 0          | 1          | 0          | 360     |
| pacbio    | 0          | 24         | 0          | 27      |
| ont       | 0          | 18         | 18         | 18      |

The effect scales with error rate exactly as it should — negligible on illumina
(0.18% error), near-total on ONT (4.4% error), and ONT is the only chemistry
where even the _adversarial_ `antioracle` flips every read. **The Viterbi/DP
path is genuinely quality-consuming.** It is not the same finding as the greedy
arm, and the manuscript must not describe them with one sentence.

### Stage C — end-to-end, and the trap in reading it

Running full `assemble_genome(corrector = :iterative)` under all four conditions
gave **byte-identical assemblies in all three chemistries**. Taken at face value
that contradicts Stage B.

It does not, and the corrector diagnostics say why:

| chemistry | assemblies identical | quality reached graph | decoder ever ran | skip_fraction per pass     |
| --------- | -------------------- | --------------------- | ---------------- | -------------------------- |
| illumina  | yes                  | yes                   | **no**           | [1.0, 0.942, 0.983, 0.986] |
| pacbio    | yes                  | yes                   | **no**           | [1.0, 1.0, 1.0, 1.0]       |
| ont       | yes                  | yes                   | **no**           | [1.0, 1.0, 1.0, 1.0]       |

`decode_fraction_per_pass` was `[0.0, 0.0, 0.0, 0.0]` in every condition: the
`skip_solid` gate — auto-enabled, not requested (`skip_solid_requested = false`,
`skip_solid_effective = true`) — skipped essentially every read on every pass,
so the quality-consuming decoder **never executed once**. The only corrections
applied were 120 from the count-based Stage-0 cheap-correct. Meanwhile
`mean_joint_quality` _did_ differ between conditions (34.33 vs 3.95 on ONT),
proving quality reached the graph and simply never reached a decision.

So Stage C is a **gating** result, not a quality-blindness result. Reported
without the diagnostics it would have been a false and much more damaging claim
than the one this bead started from.

### Quality API with no production callers

Zero non-test, non-docstring call sites anywhere in `src/`: `filter_by_quality`
(`quality-functions.jl:271`), `get_vertex_min_quality` (`:218`),
`mean_path_quality` (`:86`), `combine_phred_scores` (`:43`). The quality
machinery is larger than the set of paths that use it.

### One more count-based default worth knowing

`weighted_graph_from_rhizomorph` defaults to `edge_weight = count_evidence`
(`src/rhizomorph/algorithms/path-finding.jl:327`). The Phred-derived alternative
`edge_quality_weight` (`:270`) exists and is wired up in exactly one place —
`polish_fastq(weighting = :quality)` in
`src/viterbi-polishing-and-error-correction.jl`. Any Viterbi decode over a
default-weighted graph therefore has a quality-aware **emission** model and a
count-based **transition** model. That split is worth stating precisely if the
manuscript credits "Phred-derived edge weights".

## 4. Defect or design choice?

**Design choice, with contemporaneous evidence.** Commit `9ac3c9426` (PR #335,
2026-06-29) moved the qualmer arm onto unitig extraction and says so directly:

> qualmer arm's primary path-finder found no substantial paths on branchy graphs
> and fell to a placeholder walk emitting ~one fragment per vertex, so the
> largest contig collapsed far below genome length. Now extracts maximal unitigs
> via find_contigs_next (**the same routine the k-mer fallback already used**)
> and propagates per-base quality via `_qualmer_path_to_consensus_fastq`.

That is an explicit, deliberate decision to adopt the k-mer arm's traversal in
order to fix a real degenerate-assembly failure, keeping quality for emission
only. The same PR's review notes knowingly retain `_quality_weighted_walk` and
`_qualmer_vertex_score` as the live fallback. Related awareness appears
elsewhere: `src/rhizomorph/assembly.jl` warns that placeholder Q40 makes the
corrector's "decisions coverage-driven only", and a review comment there
explicitly worries that a silent substitution "masks that the quality model is
inert".

**The code is behaving as designed. What is defective is the claim built on top
of it**, in two places:

1. The Track-A pilot presents `qualmer` and `kmer` as two decoder arms. On
   assembly metrics they are one arm measured twice, plus overhead. Any
   manuscript text reading a quality effect into that comparison is unsupported.
2. Crediting Phred-derived edge weights as the framework's core construction is
   accurate for the Viterbi/corrector machinery (Stage A/B) and **inaccurate for
   the default assembly path**, which is what every committed Track-A row was
   produced by.

The existing guard `scripts/track_a_pilot_cv.py` in the manuscript repo raises
if the two arms ever _disagree_. It tests the wrong direction and cannot fire
while the channel is dead; it is not evidence of anything and was not relied on
here.

## 5. What changed on this branch

Per explicit maintainer decision, both the documentation and an opt-in fix
landed.

**Default behaviour is unchanged.** `traversal_weighting = :evidence` is an
exact no-op on every path, so no existing benchmark result moves. Verified: 344
assertions pass across eight existing suites — `reassembly_graph_reuse_test`
(20, and it asserts byte-identity), `rhizomorph_qualmer_coverage_prefilter_test`
(89), `rhizomorph_efficiency_modes_test` (57), `dna_qualmer_graph_test` (44),
`rhizomorph_qualmer_canonical_traversal_test` (3),
`rhizomorph_qualmer_rc_evidence_test` (6),
`rhizomorph_read_assembly_recovery_test` (6), and
`assembly_core_helper_branches_test` (119), which exercises the one-argument
`_qualmer_vertex_score` / `_qualmer_edge_score` methods left deliberately
untouched.

**New opt-in:**
`assemble_genome(reads; traversal_weighting = :quality, traversal_min_quality = 20.0)`.
Because `find_contigs_next` / `find_eulerian_paths_next` expose no scoring hook,
the only way quality can influence them is by changing the vertex set they
traverse, so the option prunes qualmer vertices whose **weakest position** falls
below the Phred floor, and additionally makes the fallback walk's vertex and
edge scores Phred-derived.

The weakest-position statistic is deliberate. A sequencing error corrupts one
base, but the k-mer containing it spans `k` positions, so at k=11 a Q2 error
base among ten Q40 neighbours _averages_ to ~Q36 — above any sane floor. A
mean-based filter would silently pass every error k-mer in the graph; this was
caught during implementation because the first version pruned nothing.

**Test:** `test/4_assembly/qualmer_quality_channel_test.jl` pins both halves.
Its characterization block asserts the default IS quality-invariant (and that
the k-mer arm reaches the same contigs); its opt-in block asserts identical
reads with different quality produce **different** assemblies — an assertion
that fails against the default path, as the characterization block in the same
file demonstrates.

## 6. What I could not determine

- **Whether the corrector is quality-blind end-to-end.** Stage C never ran the
  decoder (`skip_fraction` ~1.0 on every pass), so the end-to-end invariance is
  uninformative about quality. Settling it requires re-running with the
  solid-read gate disabled, at coverage/error settings where reads are not all
  classified solid. Stage B answers the narrower question — the decoder itself
  is quality-sensitive — and that is the only corrector claim made here.
- **Whether `skip_fraction ~ 1.0` is correct behaviour** at these settings, or
  whether the gate is over-eager and quietly bypassing the graph-HMM core on
  realistic inputs. That is a separate and potentially significant question;
  this work only observed it.
- **Real-instrument quality.** No Conda root exists on this host, so ART / pbsim
  / badread were unavailable and reads come from a self-contained Julia
  simulator (deterministic via StableRNGs, chemistry-stratified, never pooled).
  This is a strength for the _positive control_ — the oracle condition requires
  knowing which bases were corrupted, which an external binary cannot tell us,
  and an oracle signal is strictly stronger than real Phred. It does mean error
  _topology_ is idealised. The committed 66-pair Track-A table, which did use
  the real simulators, agrees.
- **Magnitude of the opt-in's effect on assembly quality.** The option is proven
  to make quality load-bearing; whether Q20 is a good floor, and whether it
  improves NGA50 against a reference, is unmeasured. It is opt-in and
  unbenchmarked, and should be treated as such.
- **Other graph families.** Only the fixed-length qualmer path was probed. The
  quality-aware OLC / FASTQ-graph builders in `src/rhizomorph/variable-length/`
  were not tested.

## 7. Reproducing

```bash
julia --project=. benchmarking/qualmer_quality_channel_probe.jl            # 18-cell grid
julia --project=. benchmarking/qualmer_quality_channel_probe.jl --smoke    # 1 cell
julia --project=. benchmarking/qualmer_corrector_quality_sensitivity.jl    # 3-stage scoping
julia --project=. test/4_assembly/qualmer_quality_channel_test.jl          # fast CI guard
```

Reference is the committed
`benchmarking/fixtures/viterbi_accuracy/phix174_nc001422.fasta` (5,386 bp); k =
31 in the probes, k = 11 in the test. No network or Conda required.
