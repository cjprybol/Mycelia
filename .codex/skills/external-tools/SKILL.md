---
name: external-tools
description: External tool integration and citation requirements for Mycelia
---

# Mycelia External Tool Integration

Rules for integrating external bioinformatics tools.

## Tool Integration Guidelines

- Keep external tool calls isolated in helpers under `src/`
- Integrates with: Bioconda, SLURM, Rclone
- Tool-specific logic should be in dedicated helper modules

## Citation Requirements

When using external tools or databases, ensure proper citations:

### Required Citations

Always cite:
- Third-party bioinformatics tools (e.g., BWA, BLAST, etc.)
- Databases used (NCBI, UniProt, etc.)
- Algorithms from papers
- Data sources with DOIs

### Citation Format

Include in:
- Code comments with DOI/paper reference
- README or documentation
- Manuscript methods section

Example:
```julia
# Uses minimap2 for alignment
# Citation: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
# DOI: 10.1093/bioinformatics/bty191
```

## Testing External Tools

- Tests must run with `MYCELIA_RUN_EXTERNAL=true` without extra flags
- Use simulated inputs where possible
- Use default database paths
- Document any tool version requirements

## Reproducibility

Track for each external tool:
- Version pinning
- Input/output specifications
- Provenance information
- Default parameters used

## Adding New Tool Wrappers

1. Create isolated helper in `src/`
2. Document tool version requirements
3. Add citation in code comments
4. Create test in `test/8_tool_integration/`
5. Use `MYCELIA_RUN_EXTERNAL` flag for tests
