# Rhizomorph Support Matrix (single/double/canonical × graph type × alphabet)

> For detailed implementation guidance, code examples, and architectural
> rationale, see
> [rhizomorph-graph-ecosystem-plan.md](rhizomorph-graph-ecosystem-plan.md).

Legend: ✅ supported, 🚫 not applicable, ⏳ pending/partial (documented).
Testing gaps: variable-length and n-gram traversal coverage is still light;
simplification edge-removal behavior needs broader tests; end-to-end suites
assert evidence positions/quality, doublestrand/canonical conversions, and
GFA/JLD2 round-trips.

| Graph Type                   | Alphabet         | Singlestrand | Doublestrand | Canonical | Notes                                                                                  |
| ---------------------------- | ---------------- | ------------ | ------------ | --------- | -------------------------------------------------------------------------------------- |
| K-mer                        | DNA              | ✅           | ✅           | ✅        | core k-mer graph-construction, traversal/tests (singlestrand, doublestrand, canonical) |
| K-mer                        | RNA              | ✅           | ✅           | ✅        | traversal/tests in place                                                               |
| K-mer                        | AA               | ✅           | 🚫           | 🚫        | no reverse complement; errors tested                                                   |
| Qualmer                      | DNA              | ✅           | ✅           | ✅        | doublestrand/canonical traversal tests added                                           |
| Qualmer                      | RNA              | ✅           | ✅           | ✅        | canonical traversal test added                                                         |
| Qualmer                      | AA               | 🚫           | 🚫           | 🚫        | not applicable                                                                         |
| N-gram                       | String           | ✅           | 🚫           | 🚫        | non-RC; errors tested                                                                  |
| Variable-length OLC (FASTA)  | DNA/RNA          | ✅           | ✅           | ✅        | singlestrand build; conversions implemented in variable-length/strand-conversions.jl   |
| Variable-length OLC (FASTA)  | AA/String        | ✅           | 🚫           | 🚫        | RC not defined; errors enforced                                                        |
| Variable-length OLC (FASTQ)  | DNA/RNA          | ✅           | ✅           | ✅        | quality-aware; conversions implemented                                                 |
| Variable-length OLC (FASTQ)  | AA/String        | 🚫           | 🚫           | 🚫        | not applicable                                                                         |
| Variable-length OLC (String) | Unicode          | ✅           | 🚫           | 🚫        | non-RC; errors enforced                                                                |
| Reduced AA alphabets         | AA/String inputs | ✅           | 🚫           | 🚫        | covered in k-mer/ngram and OLC tests                                                   |

## Graph Query Functions

Source: `src/rhizomorph/core/graph-query.jl`. All functions operate on
`MetaGraphsNext.MetaGraph` instances.

### Basic Properties and Data Access

| Function                                           | Applicable Alphabets | Applicable Modes | Notes                                                                                                        |
| -------------------------------------------------- | -------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------ |
| `vertex_count(graph)`                              | All                  | All              | Wrapper around `MetaGraphsNext.labels`                                                                       |
| `edge_count(graph)`                                | All                  | All              | Wrapper around `MetaGraphsNext.ne`                                                                           |
| `has_vertex(graph, kmer)`                          | All                  | All              | Label existence check                                                                                        |
| `has_edge(graph, src, dst)`                        | All                  | All              | Edge existence check                                                                                         |
| `get_vertex_data(graph, kmer)`                     | All                  | All              | Returns vertex data or `nothing`                                                                             |
| `get_edge_data(graph, src, dst)`                   | All                  | All              | Returns edge data or `nothing`                                                                               |
| `get_vertex_observation_count(graph, kmer)`        | All                  | All              | Total observations across datasets                                                                           |
| `get_edge_observation_count(graph, src, dst)`      | All                  | All              | Total observations across datasets                                                                           |
| `get_all_vertices(graph)`                          | All                  | All              | Collect all vertex labels                                                                                    |
| `filter_vertices_by_observation_count(graph, min)` | All                  | All              | Vertices with `>= min` observations                                                                          |
| `get_graph_statistics(graph)`                      | All                  | All (directed)   | Returns Dict with vertex/edge/source/sink/branch/join counts. Uses in/out degree so assumes directed graphs. |

### Traversal and Topology

All functions work on canonical (undirected) graphs. For undirected graphs,
Graphs.jl treats `outneighbors == inneighbors == neighbors`, so the in/out
distinction collapses to symmetric degree. Sources and sinks will always be
empty on canonical graphs with edges.

| Function                              | Applicable Alphabets | Applicable Modes | Notes                                                                  |
| ------------------------------------- | -------------------- | ---------------- | ---------------------------------------------------------------------- |
| `get_outgoing_neighbors(graph, kmer)` | All                  | All              | On undirected: returns all neighbors                                   |
| `get_incoming_neighbors(graph, kmer)` | All                  | All              | On undirected: returns all neighbors (same as outgoing)                |
| `get_outgoing_degree(graph, kmer)`    | All                  | All              | On undirected: equals vertex degree                                    |
| `get_incoming_degree(graph, kmer)`    | All                  | All              | On undirected: equals vertex degree                                    |
| `is_linear_path(graph, kmer)`         | All                  | All              | On undirected: degree 1 (pendant vertex)                               |
| `is_source_vertex(graph, kmer)`       | All                  | All              | On undirected: degree 0 (isolated vertex)                              |
| `is_sink_vertex(graph, kmer)`         | All                  | All              | On undirected: degree 0 (isolated vertex)                              |
| `is_valid_path(graph, path)`          | All                  | All              | Checks consecutive edge existence                                      |
| `get_all_sources(graph)`              | All                  | All              | On undirected: only isolated vertices (no edges)                       |
| `get_all_sinks(graph)`                | All                  | All              | On undirected: only isolated vertices (no edges)                       |
| `get_graph_statistics(graph)`         | All                  | All              | On undirected: source/sink/branch/join counts reflect symmetric degree |

### Unitig Assembly

| Function                              | Applicable Alphabets      | Applicable Modes | Notes                                                       |
| ------------------------------------- | ------------------------- | ---------------- | ----------------------------------------------------------- |
| `extend_unitig_forward(graph, kmer)`  | All                       | All              | On undirected: follows degree-1 chains                      |
| `extend_unitig_backward(graph, kmer)` | All                       | All              | On undirected: equivalent to forward (symmetric neighbors)  |
| `get_maximal_unitig(graph, kmer)`     | All                       | All              | On undirected: finds maximal degree-1 chain containing kmer |
| `assemble_path_sequence(path)`        | DNA, RNA, AA (Kmer types) | All              | Reconstructs sequence from kmer path                        |

### Reverse Complement Analysis

| Function                   | Applicable Alphabets | Applicable Modes                           | Notes                                                                                                                                                                                                                            |
| -------------------------- | -------------------- | ------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_self_rc_closed(graph)` | DNA, RNA             | Singlestrand, Doublestrand (directed only) | Detects self-complementary edge sets — returns `true` when every edge's RC counterpart also exists. Doublestrand graphs are always closed by construction. Returns `false` for undirected (canonical) graphs or edgeless graphs. |
