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
"""
module Rhizomorph

# Import necessary packages for Rhizomorph functionality
import BioSequences
import FASTX
import Graphs
import Kmers
import MetaGraphsNext

# Load core components
include("core/enums.jl")
include("core/evidence-structures.jl")
include("core/vertex-data.jl")
include("core/edge-data.jl")
include("core/evidence-functions.jl")
include("core/graph-construction.jl")
include("core/graph-query.jl")
include("core/quality-functions.jl")

end  # module Rhizomorph
