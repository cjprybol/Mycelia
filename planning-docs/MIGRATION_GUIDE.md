# Rhizomorph migration guide

Old Mycelia graph calls → current Rhizomorph entry points.

This guide is for external callers moving off the pre-Rhizomorph
MetaGraphsNext builders. It records breaking changes, updated examples,
deprecation notes, and a short FAQ. Signatures were checked against the
loaded sources on current `master`:

- `src/rhizomorph/` (builders, GFA, paths, simplification, greedy correction)
- `src/viterbi-next.jl` and `src/viterbi-polishing-and-error-correction.jl`
- `src/sequence-graphs.jl` (legacy stranded DNA k-mer graph)
- `src/gqc.jl` (`read_gfa_records` only)
- `planning-docs/RHIZOMORPH_SUPPORT_MATRIX.md`
- `planning-docs/rhizomorph-graph-ecosystem-plan.md`
- `docs/src/workflow-map.md` and the docstrings on the functions below

Mycelia does not `export` symbols. Call them as `Mycelia.fn(...)` or
`Mycelia.Rhizomorph.fn(...)`.

The ecosystem plan still sketches removed names (`build_kmer_graph_next`,
`fragment_to_fixed_length`, `remove_quality_scores`,
`collapse_to_variable_length`) and still says error correction is missing.
Those sentences are a design snapshot. Use this guide and the source
docstrings for names that exist today. Alphabet × strand support is the
table in `RHIZOMORPH_SUPPORT_MATRIX.md`, not a second matrix here.

## What changed

The MetaGraphsNext builders used to live in the `Mycelia` module
(`src/kmer-graphs.jl`, `src/qualmer-graphs.jl`, `src/sequence-graphs-next.jl`,
`src/fasta-graphs.jl`, `src/fastq-graphs.jl`, `src/string-graphs.jl`,
`src/graph-core.jl`). Those files are gone. The same jobs now go through
`Mycelia.Rhizomorph`.

There used to be two `GraphMode` enums:

| Old name | Where it lived | Values |
| --- | --- | --- |
| `Mycelia.GraphMode` / `Mycelia.SingleStrand` / `Mycelia.DoubleStrand` | `src/graph-core.jl` | `SingleStrand`, `DoubleStrand` |
| `Mycelia.Rhizomorph.GraphMode` | `src/rhizomorph/core/enums.jl` | same two values at the time of that split |

`Mycelia.GraphMode` is gone. The only enum is
`Mycelia.Rhizomorph.GraphMode`, and it now has a third value, `Canonical`.

User-facing k-mer and qualmer builders do **not** take that enum. They take
`mode::Symbol` (`:singlestrand`, `:doublestrand`, `:canonical`).
`AssemblyConfig`, `assemble_genome`, and `read_gfa_next` still take the enum.

## Breaking changes

1. **Module and argument order.** `Mycelia.build_kmer_graph_next(kmer_type, records; graph_mode=...)` is not defined. The replacement is `Mycelia.Rhizomorph.build_kmer_graph(records, k; mode=...)`. `k` is an `Int`. The alphabet is inferred from the records (`type_hint` is optional).

2. **Default strand mode flipped, and the old `DoubleStrand` label does not mean the new one.** `build_kmer_graph_next` defaulted to `graph_mode=DoubleStrand`. That branch called `Mycelia.count_canonical_kmers` and still built a **directed** graph whose labels were canonical k-mers. Anything else used `Mycelia.count_kmers` (k-mers as observed). `build_kmer_graph` defaults to `mode=:singlestrand` (observed orientation, directed). The new `:doublestrand` keeps forward and reverse-complement k-mers as separate directed vertices. The new `:canonical` merges each pair onto one vertex and yields an **undirected** graph. Old `DoubleStrand` is closer to canonical *labels* than to the new doublestrand graph, and it is not a drop-in for either.

3. **`assemble_genome` uses a different default than the low-level builder.** For DNA and RNA, `assemble_genome` defaults `graph_mode` to `Mycelia.Rhizomorph.DoubleStrand` (both strands). Amino-acid and `String` inputs default to `SingleStrand`. Passing `mode=:singlestrand` into `assemble_genome` is not a recognized keyword.

4. **Empty k-mer input.** `build_kmer_graph_next` warned and returned an empty graph when no k-mers were counted. `build_kmer_graph_singlestrand` and `build_kmer_graph_canonical` throw `ArgumentError` when `records` is empty.

5. **Qualmer input.** `build_qualmer_graph` requires `Vector{FASTX.FASTQ.Record}`. FASTA records are an error. The old `build_qualmer_graph_next(kmer_type, fastq_records; graph_mode=DoubleStrand)` type argument and enum are gone.

6. **Variable-length overlap defaults.** `build_biosequence_graph` and `build_quality_biosequence_graph` defaulted to `min_overlap=100` and `graph_mode=DoubleStrand`. `build_fasta_graph`, `build_fastq_graph`, and `build_string_graph` default to `min_overlap=3` and build a singlestrand graph. Overlap search only keeps odd lengths. `_normalize_min_overlap` rounds an even value up by 1, so the old default `100` becomes `101`. There is no `graph_mode` keyword on these builders. DNA/RNA doublestrand or canonical form is a second call: `convert_variable_length_to_doublestrand` or `convert_variable_length_to_canonical`. `build_fastq_graph` has no `min_quality` filter.

7. **No k-mer graph → FASTA-graph converter.** `kmer_graph_to_biosequence_graph` and `qualmer_graph_to_quality_biosequence_graph` were removed. `collapse_linear_chains!` only collapses graphs whose vertex data already have a `sequence` field (variable-length and string vertices). It does not turn a k-mer graph into an OLC graph.

8. **GFA moved, and the header string changed.** `Mycelia.write_gfa_next` / `Mycelia.read_gfa_next` are gone. Graph I/O is `Mycelia.Rhizomorph.write_gfa_next` and `Mycelia.Rhizomorph.read_gfa_next`. New GFA headers use `MY:Z:Mycelia-Rhizomorph` (the old k-mer writer wrote `MY:Z:Mycelia-Next`). `write_gfa_next` gained `segment_naming` (`:numeric` or `:sequence`). `Mycelia.Rhizomorph.write_gfa` is a different method: it takes an `AssemblyResult`, not a raw graph.

9. **`read_gfa_next` does not apply link orientation.** The k-mer method accepts `graph_mode` (default `DoubleStrand`) but does not branch on it, and it does not rewrite segment k-mers from link `+`/`-`. Each link becomes an edge between the segment sequences as stored. Do not pass `Canonical` expecting a canonical graph. Build that with `build_kmer_graph(...; mode=:canonical)` or `convert_to_canonical`, and export it with `write_gfa_canonical`.

10. **Path and simplification names keep the `_next` suffix.** `find_eulerian_paths` and `simplify_graph` are not defined. The functions are `find_eulerian_paths_next` and `simplify_graph_next`. `simplify_graph_next` requires the bubble vector from `detect_bubbles_next`; it is not a one-argument “simplify this graph” call.

11. **Correction is not a Rhizomorph graph method named `correct_errors_at_k`.** That function sits in `src/development/intelligent-assembly.jl`, which `src/Mycelia.jl` does not include. The loaded correction entry points are below.

## Graph construction

`records` below are `FASTX` records (`Mycelia.open_fastx` yields them).
File helpers open plain or `.gz` / `.bz2` / `.xz` paths themselves.

### Fixed-length builders

| Removed (`Mycelia.`) | Current |
| --- | --- |
| `build_kmer_graph_next(kmer_type, records; graph_mode=DoubleStrand)` | `Mycelia.Rhizomorph.build_kmer_graph(records, k; mode=:singlestrand, dataset_id="dataset_01", type_hint=nothing, ambiguous_action=:dna, memory_profile=:full)` |
| mode-specific construction inside that function | `build_kmer_graph_singlestrand`, `build_kmer_graph_doublestrand`, `build_kmer_graph_canonical` (same record / `k` arguments; no `mode` symbol) |
| `build_qualmer_graph_next(kmer_type, fastq_records; graph_mode=DoubleStrand)` | `Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, k; mode=:singlestrand, memory_profile=:full, min_count=1)` |
| — | `build_qualmer_graph_singlestrand`, `build_qualmer_graph_doublestrand`, `build_qualmer_graph_canonical` |
| `build_string_graph_next(fasta_records, ngram_length; graph_mode=SingleStrand)`, `strings_to_ngram_graph`, `string_to_ngram_graph`, `build_string_ngram_graph_next` | `Mycelia.Rhizomorph.build_ngram_graph(strings, n; dataset_id="dataset_01", memory_profile=:full)` where `strings::Vector{String}` |
| — | `build_ngram_graph_singlestrand(strings, n)`. N-grams are singlestrand only. `:doublestrand` / `:canonical` do not apply. |

`memory_profile` for k-mer graphs:

| Value | What is kept |
| --- | --- |
| `:full` (default) | Evidence, including positions and strands |
| `:lightweight` | Counts and observation ids, no positions |
| `:ultralight` | Counts only |
| `:ultralight_quality` | Counts plus aggregated quality (FASTQ only) |
| `:lightweight_quality` | Counts, observation ids, and aggregated quality (FASTQ only) |

Qualmer graphs accept `:full`, `:ultralight_quality`, and `:lightweight_quality` only. N-gram graphs accept `:full`, `:lightweight`, and `:ultralight`.

`type_hint` is `:DNA`, `:RNA`, or `:AA`. `ambiguous_action` is `:dna`, `:rna`, `:aa`, or `:error`. Reverse complement is defined for DNA and RNA only. N-gram graphs are singlestrand. Amino-acid k-mer graphs are singlestrand: `mode=:doublestrand` calls `convert_to_doublestrand`, which errors with "Cannot create doublestrand graph for non-nucleotide k-mer labels". `mode=:canonical` calls `convert_to_canonical`, which errors with "Cannot create canonical graph for non-nucleotide k-mer labels". The support matrix marks qualmer × amino acid as not applicable. The singlestrand qualmer core still dispatches `Kmers.AAKmer` when inference returns `:AA`; doublestrand and canonical still fail in those converters. The `build_qualmer_graph` docstring says DNA and RNA only.

File entry points, all under `Mycelia.Rhizomorph`:

- `build_kmer_graph_from_file(path, k; ...)` and `build_kmer_graph_from_files(paths, k; ...)`
- `build_qualmer_graph_from_file` / `build_qualmer_graph_from_files`
- `build_ngram_graph_from_file` / `build_ngram_graph_from_files`

`.fna`, `.frn`, and `.faa` supply `type_hint` when the caller does not. `.fa` and `.fasta` stay ambiguous.

Add more reads to an existing k-mer or qualmer graph with
`Mycelia.Rhizomorph.add_observations_to_graph!(graph, records, k; dataset_id=...)`.
The graph must already contain vertices. An empty graph raises; build it
with a `build_*` function first.

Strand changes after construction:

- Fixed-length: `convert_to_doublestrand`, `convert_to_canonical`
- FASTA/FASTQ OLC, DNA/RNA labels only: `convert_variable_length_to_doublestrand`, `convert_variable_length_to_canonical`

The variable-length converters throw for amino acids and generic strings.

### Variable-length (OLC) builders

| Removed (`Mycelia.`) | Current |
| --- | --- |
| `build_biosequence_graph(fasta_records; graph_mode=DoubleStrand, min_overlap=100)` and `build_biosequence_graph_next(...)` | `Mycelia.Rhizomorph.build_fasta_graph(records; dataset_id="dataset_01", min_overlap=3, type_hint=nothing, ambiguous_action=:dna)` |
| `build_quality_biosequence_graph(fastq_records; graph_mode=DoubleStrand, min_overlap=100, min_quality=0x14)` and `build_qualmer_biosequence_graph_next(...)` | `Mycelia.Rhizomorph.build_fastq_graph(records; dataset_id="dataset_01", min_overlap=3)` |
| string OLC helpers that took FASTA records | `Mycelia.Rhizomorph.build_string_graph(strings; dataset_id="dataset_01", min_overlap=3)` with `strings::Vector{String}` |
| — | `build_token_graph(token_sequences; dataset_id="dataset_01")` where `token_sequences::Vector{Vector{String}}` |

Matching `*_from_file` / `*_from_files` helpers exist for FASTA, FASTQ, and
string graphs. Lower-level names `build_fasta_graph_olc`,
`build_fastq_graph_olc`, and `build_string_graph_olc` are the same
constructions; the short names above are the ones to call.

### High-level assembly

`Mycelia.Rhizomorph.assemble_genome(reads; kwargs...)` returns an
`AssemblyResult`. `reads` are FASTA/FASTQ records. Pass `k` or
`min_overlap`, not both. If neither is set, `k` becomes 31.

```julia
result = Mycelia.Rhizomorph.assemble_genome(records; k = 31,
    graph_mode = Mycelia.Rhizomorph.SingleStrand)

result = Mycelia.Rhizomorph.assemble_genome(records; min_overlap = 100)

config = Mycelia.Rhizomorph.AssemblyConfig(;
    k = 25,
    sequence_type = BioSequences.LongDNA{4},
    graph_mode = Mycelia.Rhizomorph.DoubleStrand,
    use_quality_scores = true,
)
result = Mycelia.Rhizomorph.assemble_genome(records, config)
```

`graph_mode` is the enum (`SingleStrand`, `DoubleStrand`, `Canonical`), not
a `Symbol`. `corrector` is `:none` (default) or `:iterative`. There is no
`method=` keyword. `AssemblyConfig` rejects `DoubleStrand` and `Canonical`
for amino acids and for `String`.

Other loaded assembly entry points, for hybrid layouts rather than the
graph builders above: `assemble_hybrid`, `assemble_unicycler_hybrid`,
`assemble_autocycler_polished`.

## Error correction

Three loaded surfaces. They are not interchangeable.

### Legacy stranded DNA graph (still loaded)

`src/sequence-graphs.jl` remains for the DNA Viterbi oracle. It returns a
`MetaGraphs.MetaDiGraph`, not a Rhizomorph `MetaGraph`.

```julia
graph = Mycelia.build_stranded_kmer_graph(Kmers.DNAKmer{31}, records)
corrected = Mycelia.viterbi_maximum_likelihood_traversals(graph;
    error_rate = 0.01,
    verbosity = "dataset")
```

`kmer_type` must be a concrete DNA k-mer type. `verbosity` is `"debug"`,
`"reads"`, or `"dataset"`. `error_rate` must be `< 0.5`. The default is
`1 / (k + 1)`. The return value is a vector of `FASTX.FASTA.Record`.

Related helpers on that same graph: `sequence_to_stranded_path`,
`Mycelia.path_to_sequence(kmers, path)` (DNA only), and `edge_probability`.
`Mycelia.path_to_sequence` is **not** the Rhizomorph method; that one takes
`(path, graph)`.

### Viterbi seam for both graph types

`Mycelia.correct_observations` is defined on `Mycelia`, not under
`Rhizomorph`. A legacy `MetaDiGraph` still uses
`viterbi_maximum_likelihood_traversals`. A Rhizomorph `MetaGraph` uses the
weighted-graph decoder.

```julia
config = Mycelia.ViterbiCorrectionConfig(; error_rate = 0.01)
result = Mycelia.correct_observations(graph; config = config)
```

`ViterbiCorrectionConfig` fields that callers actually set include
`error_rate` (default `0.01`), `verbosity` (default `"dataset"`),
`strand_mode` (default `:auto`), `alphabet` (default `:auto`), and
`start_strand` (default `Mycelia.Rhizomorph.Forward`). The return value is a
`Mycelia.ViterbiCorrectionResult`.

To correct reads and then assemble, pass the corrector into assembly
instead of calling the decoder yourself:

```julia
result = Mycelia.Rhizomorph.assemble_genome(records;
    k = 31,
    corrector = :iterative)
```

`corrector=:none` is the default and does not run that front-end.

### Greedy k-mer correction

`Mycelia.Rhizomorph.correct_sequence_greedy` is a single-substitution greedy
walk. It is not the Viterbi corrector. One method takes count dictionaries.
The other reads counts off a Rhizomorph graph by stringifying labels:

```julia
corrected, positions = Mycelia.Rhizomorph.correct_sequence_greedy(
    "ACGTA", 3, graph)
```

The return value is `(corrected_sequence::String, corrected_kmer_indexes::Vector{Int})`.

## GFA

| Job | Call |
| --- | --- |
| Write any Rhizomorph `MetaGraph` to GFA 1.0 | `Mycelia.Rhizomorph.write_gfa_next(graph, path; segment_naming=:numeric)` |
| Write an undirected canonical nucleotide k-mer graph | `Mycelia.Rhizomorph.write_gfa_canonical(graph, path; segment_naming=:numeric)` |
| Read GFA as a k-mer graph | `Mycelia.Rhizomorph.read_gfa_next(path, Kmers.DNAKmer{31}, Mycelia.Rhizomorph.DoubleStrand)` |
| Read GFA and infer fixed- vs variable-length | `Mycelia.Rhizomorph.read_gfa_next(path, Mycelia.Rhizomorph.DoubleStrand; force_biosequence_graph=false)` |
| Write an `AssemblyResult` (graph, then path lines) | `Mycelia.Rhizomorph.write_gfa(result, path)` |
| Parse GFA 1.0 into segment/link structs | `Mycelia.read_gfa_records(path)` → `Mycelia.GFAGraph` |

`segment_naming=:numeric` uses segment ids `1, 2, 3, ...`.
`:sequence` uses the sequence text as the id (Bandage `--names`). Any other
value throws `ArgumentError`.

`write_gfa_canonical` throws if the graph is directed, if a label is not a
`Kmers.DNAKmer` or `Kmers.RNAKmer`, or if an edge has no `(k - 1)` overlap in
any orientation. Directed graphs, including singlestrand and doublestrand
k-mer graphs, use `write_gfa_next`.

`read_gfa_records` does not build a `MetaGraph`. It keeps header lines,
segments, and links (including `LN:i`, `DP:f`, `RC:i`, `KC:i` when present).

Removed writers with no same-named replacement: `write_biosequence_gfa`,
`write_quality_biosequence_gfa`. Variable-length graphs go through
`write_gfa_next`.

`algorithms/io.jl` mentions JLD2 in its file header. The functions in that
file are GFA only. Tests that round-trip a graph call `JLD2.save` /
`JLD2.load` directly. The wrapped helpers are for an `AssemblyResult`:

```julia
Mycelia.Rhizomorph.save_assembly(result, "assembly.jld2")
loaded = Mycelia.Rhizomorph.load_assembly("assembly.jld2")
```

## Paths, simplification, and generation

These are the current names. Unsuffixed `find_eulerian_paths` and
`simplify_graph` are not defined.

```julia
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31; mode = :singlestrand)

paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
sequence = Mycelia.Rhizomorph.path_to_sequence(first(paths), graph)

bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph;
    min_bubble_length = 2, max_bubble_length = 100)
simplified = Mycelia.Rhizomorph.simplify_graph_next(graph, bubbles)

Mycelia.Rhizomorph.remove_tips!(graph; min_support = 1)

weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
generated = Mycelia.Rhizomorph.generate_sequences(weighted, 10;
    walk_length = 200, seed = 42)

contigs = Mycelia.Rhizomorph.find_contigs_next(graph; min_contig_length = 500)
```

`find_eulerian_paths_next` returns a `Vector` of label paths (empty when no
Eulerian path exists). `path_to_sequence` also accepts a `GraphPath`.
`generate_sequences` accepts either a weighted graph or an evidence graph; an
evidence graph is converted inside the call. Each item is a named tuple
`(sequence, path, walk_probability, length)`.

`collapse_linear_chains!(graph)` is in-place and only for vertex data that
already store a `sequence`.

Query helpers on a `MetaGraph` live in `src/rhizomorph/core/graph-query.jl`
and match the support-matrix names: `vertex_count`, `edge_count`,
`has_vertex`, `has_edge`, `get_vertex_data`, `get_edge_data`,
`get_graph_statistics`, `assemble_path_sequence` (k-mer paths only). On an
undirected canonical graph, in-neighbors and out-neighbors are the same set,
and source/sink queries do not mark the ends of a path.

Metrics that exist today: `compute_betti_numbers`,
`compute_eigenvector_centrality`, `graph_transition_entropy`,
`jensen_shannon_divergence`, `evaluate_generation_quality`.

## Updated examples

### K-mer graph, observed strand, then GFA

```julia
import Mycelia

records = collect(Mycelia.open_fastx("reads.fastq.gz"))
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31;
    dataset_id = "sample_a",
    mode = :singlestrand)

Mycelia.Rhizomorph.write_gfa_next(graph, "sample_a.gfa")
Mycelia.Rhizomorph.write_gfa_next(graph, "sample_a.names.gfa";
    segment_naming = :sequence)
```

Old call this replaces:

```julia
# removed
graph = Mycelia.build_kmer_graph_next(Kmers.DNAKmer{31}, records;
    graph_mode = Mycelia.SingleStrand)
Mycelia.write_gfa_next(graph, "sample_a.gfa")
```

### Canonical DNA graph

```julia
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31; mode = :canonical)
# equivalent:
# graph = Mycelia.Rhizomorph.build_kmer_graph_canonical(records, 31)

Mycelia.Rhizomorph.write_gfa_canonical(graph, "canonical.gfa")
```

This is the usual replacement when old code passed
`graph_mode=Mycelia.DoubleStrand` in order to count canonical k-mers. The
new graph is undirected. Read it back with `read_gfa_next` only if a
directed import is acceptable; that reader does not rebuild canonical mode.

### Qualmer graph and greedy correction

```julia
import FASTX

graph = Mycelia.Rhizomorph.build_qualmer_graph(fastq_records, 21;
    mode = :doublestrand,
    memory_profile = :ultralight_quality,
    min_count = 2)

corrected, positions = Mycelia.Rhizomorph.correct_sequence_greedy(
    string(FASTX.sequence(fastq_records[1])), 21, graph)
```

### Viterbi correction on a Rhizomorph graph

```julia
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31; mode = :doublestrand)
config = Mycelia.ViterbiCorrectionConfig(; error_rate = 0.01)
corrected = Mycelia.correct_observations(graph; config = config)
```

### FASTA overlap graph

```julia
contigs = collect(Mycelia.open_fastx("contigs.fna"))
graph = Mycelia.Rhizomorph.build_fasta_graph(contigs;
    min_overlap = 101,
    type_hint = :DNA)
graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(graph)
Mycelia.Rhizomorph.write_gfa_next(graph, "contigs.gfa")
```

`101` is what the old default `100` becomes after even values are rounded up.
The new default of 3 draws many more edges.

### N-gram graph

```julia
graph = Mycelia.Rhizomorph.build_ngram_graph(
    ["ACGTACGT", "CGTACGTA"], 3; memory_profile = :lightweight)
```

### Read a GFA file

```julia
import Kmers
import Mycelia

graph = Mycelia.Rhizomorph.read_gfa_next(
    "assembly.gfa",
    Kmers.DNAKmer{31},
    Mycelia.Rhizomorph.SingleStrand)

detected = Mycelia.Rhizomorph.read_gfa_next("assembly.gfa")
records = Mycelia.read_gfa_records("assembly.gfa")
```

## Deprecation notes

No Rhizomorph builder uses `Base.@deprecate`. Removed Mycelia builders are
absent: calling them is an `UndefVarError`, not a warning with a forwarding
method.

Still accepted, and marked deprecated in the docstring, is the `lightweight`
keyword on `build_kmer_graph`, `build_kmer_graph_from_file`,
`build_kmer_graph_from_files`, and `build_ngram_graph`:

| Keyword | Effect |
| --- | --- |
| `lightweight=true` | `memory_profile=:lightweight` |
| `lightweight=false` | `memory_profile=:full` |

The keyword overrides `memory_profile` when both are passed. It does not emit
a deprecation warning. New code should pass `memory_profile` only.

`Mycelia.build_stranded_kmer_graph` is a compatibility surface for the DNA
Viterbi oracle, not a deprecated alias of `build_kmer_graph`. Its docstring
limits it to that `MetaDiGraph` property contract. New graph work should use
`Mycelia.Rhizomorph.build_kmer_graph`.

`read_gfa_next`'s `graph_mode` argument is part of the signature and defaults
to `DoubleStrand`. The current k-mer method does not use the value. Treat it
as unused input, not as a strand-mode switch.

## FAQ

**Where did `build_kmer_graph_next` go?**

It was removed with `src/kmer-graphs.jl`. Use
`Mycelia.Rhizomorph.build_kmer_graph(records, k; mode=:singlestrand)`. Pass
the integer `k`, not `Kmers.DNAKmer{k}`.

**I used `graph_mode=DoubleStrand`. Which `mode` is that?**

Not `:doublestrand`. The old flag stored canonical k-mer labels on a directed
graph. Read the breaking-changes section before picking `:canonical` or
`:doublestrand`. `assemble_genome`'s default `DoubleStrand` is the new
both-strands enum, not the old canonical counter.

**Why is `Mycelia.SingleStrand` undefined?**

That enum lived in `src/graph-core.jl`, which was removed. Use
`Mycelia.Rhizomorph.SingleStrand`, `DoubleStrand`, or `Canonical` anywhere an
enum is required. K-mer and qualmer builders want the symbols `:singlestrand`,
`:doublestrand`, and `:canonical` instead.

**Which correction function should I call?**

- Existing DNA `MetaDiGraph` from `build_stranded_kmer_graph`: `Mycelia.viterbi_maximum_likelihood_traversals`, or `Mycelia.correct_observations` (it detects that graph).
- Rhizomorph `MetaGraph`: `Mycelia.correct_observations`.
- Correct, then assemble: `assemble_genome(...; corrector=:iterative)`.
- One sequence, greedy, string k-mers: `Mycelia.Rhizomorph.correct_sequence_greedy`.

**Is `write_gfa` the graph exporter?**

`Mycelia.Rhizomorph.write_gfa` writes an `AssemblyResult`. A graph uses
`write_gfa_next`. An undirected canonical nucleotide graph uses
`write_gfa_canonical`. A parsed segment/link table uses
`Mycelia.read_gfa_records`, which does not return a `MetaGraph`.

**Does `lightweight=true` still work?**

Yes on the four builders listed under Deprecation notes. It is a silent
remap onto `memory_profile`. Prefer `memory_profile=:lightweight`.

**Can I simplify a graph with one call?**

`simplify_graph` is not defined. Call `detect_bubbles_next`, then
`simplify_graph_next(graph, bubbles)`. Tips are `remove_tips!`. Linear
unitigs on variable-length graphs are `collapse_linear_chains!`.

**Do planning-doc examples such as `method=:qualmer_graph` or `assemble_probabilistic` work?**

No. They are not defined on the loaded modules. Quality-aware assembly is
`build_qualmer_graph` or `assemble_genome` on FASTQ records
(`use_quality_scores` follows the record type). `correct_errors_at_k` is not
included by `src/Mycelia.jl`.

**Where is the support matrix?**

`planning-docs/RHIZOMORPH_SUPPORT_MATRIX.md`. DNA and RNA k-mer and qualmer
graphs support singlestrand, doublestrand, and canonical. Amino-acid k-mer
graphs and n-gram / string graphs are singlestrand only. FASTA and FASTQ OLC
graphs support doublestrand and canonical for DNA and RNA via the converters
above.
