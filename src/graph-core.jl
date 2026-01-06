"""
Core graph enums shared across fixed-length and variable-length graph modules.
"""

"""
Strand orientation for k-mer observations and transitions.
- `Forward`: k-mer as observed (5' to 3')
- `Reverse`: reverse complement of k-mer (3' to 5')
"""
@enum StrandOrientation Forward=true Reverse=false

"""
Graph mode for handling strand information.
- `SingleStrand`: Sequences are single-stranded (RNA, amino acids, or directional DNA)
- `DoubleStrand`: Sequences are double-stranded DNA/RNA with canonical representation
"""
@enum GraphMode SingleStrand DoubleStrand

"""
    edge_data_weight(edge_data)

Return the `weight` field from edge data types that store explicit weights.
"""
function edge_data_weight(edge_data)
    return edge_data.weight
end
