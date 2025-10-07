"""
Core enumerations for the Rhizomorph graph ecosystem.

These enums define fundamental concepts used throughout the graph construction
and manipulation pipeline.
"""

"""
Strand orientation for sequence observations and transitions.

# Values
- `Forward`: Sequence as observed (5' to 3' direction)
- `Reverse`: Reverse complement of sequence (3' to 5' direction)

# Details
This enum tracks the biological orientation of sequences, k-mers, and transitions.
- In strand-specific construction: Forward and Reverse evidence are kept separate
- In non-strand-specific construction: Forward and Reverse evidence are merged

The enum is backed by Bool for memory efficiency:
- Forward = true
- Reverse = false
"""
@enum StrandOrientation Forward=true Reverse=false

"""
Graph mode for handling strand representation.

# Values
- `SingleStrand`: Graph contains one strand orientation
- `DoubleStrand`: Graph contains both forward and reverse-complement strands

# Details
This enum defines the structural representation of the graph:

**SingleStrand**: Only one strand orientation is present
- May reflect biological reality (ssRNA viruses, mRNA)
- May be a filtering choice for double-stranded molecules
- Simpler graph structure

**DoubleStrand**: Both RC orientations present with bidirectional relationships
- Both strands represented with explicit relationships
- May reflect biological reality (dsDNA)
- Created by replicating vertices/edges in RC orientation
- Reversible: can collapse to single-stranded by filtering

# Important
This is **orthogonal** to strand specificity (whether evidence is merged across RC pairs).
Strand representation (SingleStrand/DoubleStrand) is about graph structure.
Strand specificity is about evidence tracking methodology.
"""
@enum GraphMode SingleStrand DoubleStrand
