"""
Rhizomorph Graph Ecosystem

A comprehensive graph construction framework for bioinformatics that properly
handles all graph types as directed, strand-aware graphs where vertices and
edges are added only when observed.

# Core Principles
1. **Record exactly what we observe** - no inference, no canonicalization during construction
2. **Strand-specific by default** - preserve orientation as observed
3. **Evidence-based** - efficient tracking of observations through nested Dict structure
4. **Type-stable** - generic parameterization for performance

# Module Organization
- `core/`: Fundamental data structures and operations
  - `enums.jl`: StrandOrientation, GraphMode
  - `evidence-structures.jl`: Evidence entry types
  - `vertex-data.jl`: Vertex data structures for all graph types
  - `edge-data.jl`: Edge data structures for all graph types
  - `evidence-functions.jl`: Evidence manipulation and query functions
  - `graph-construction.jl`: Graph construction algorithms
  - `quality-functions.jl`: Quality score utilities
  - `graph-query.jl`: Graph query helpers
- `fixed-length/`: Fixed-length graph builders (k-mers, qualmers, n-grams)
- `variable-length/`: Variable-length graph builders (FASTA, FASTQ, strings)
- `algorithms/`: Graph algorithms
  - `path-finding.jl`: Eulerian paths and sequence reconstruction
  - `io.jl`: GFA and serialization I/O
  - `simplification.jl`: Bubble detection and graph simplification
"""
module Rhizomorph

# Import necessary packages for Rhizomorph functionality
import BioSequences
import DataStructures
import FASTX
import Graphs
import Kmers
import MetaGraphsNext
import Statistics

# Load core components
include("core/enums.jl")
include("core/evidence-structures.jl")
include("core/vertex-data.jl")
include("core/edge-data.jl")
include("core/evidence-functions.jl")
include("core/quality-functions.jl")
include("core/graph-query.jl")
include("core/graph-construction.jl")
include("core/graph-type-conversions.jl")

# Load algorithm modules
include("algorithms/path-finding.jl")
include("algorithms/io.jl")
include("algorithms/simplification.jl")
include("algorithms/repeats.jl")
include("algorithms/contigs.jl")

# Load graph builders
include("fixed-length/kmer-graphs.jl")
include("fixed-length/qualmer-graphs.jl")
include("fixed-length/ngram-graphs.jl")
include("variable-length/string-graphs.jl")
include("variable-length/fasta-graphs.jl")
include("variable-length/fastq-graphs.jl")
include("variable-length/strand-conversions.jl")

end  # module Rhizomorph
