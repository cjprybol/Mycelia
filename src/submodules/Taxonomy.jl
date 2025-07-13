"""
Taxonomy Submodule for Mycelia.jl

Handles taxonomic classification and phylogenetic analysis.
"""
module Taxonomy

import DataFrames
import CSV
import Graphs
import MetaGraphs

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Taxonomic classification
export name2taxid, taxids2lca, classify_sequences
export list_superkingdoms, list_species, get_taxonomy_lineage

# Phylogenetic analysis
export build_phylogeny, calculate_phylogenetic_distance
export parse_newick, write_newick

# Tree operations
export root_tree, prune_tree, compare_trees

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _setup_taxonkit_taxonomy, _validate_taxonomy_database
export _parse_taxonomy_output

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from taxonomy-and-trees.jl to here

end # module Taxonomy