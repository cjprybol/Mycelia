# Gene Annotation Workflows

This page documents gene annotation workflows and functions in Mycelia.

## Gene Prediction

Gene prediction identifies protein-coding sequences and other genomic features.

### ORF Callers

- `Mycelia.run_prodigal`: Prokaryotic ORF calling (single genomes or metagenomes)
- `Mycelia.run_pyrodigal`: Fast prokaryotic ORF calling for metagenomic contigs
- `Mycelia.run_prodigal_gv`: Viral ORF calling for genomes or metaviromes
- `Mycelia.run_augustus`: Eukaryotic ab initio gene prediction (species model required)
- `Mycelia.run_metaeuk`: Eukaryotic metagenomic gene prediction using a reference database

### Taxonomy-Aware Translation Tables

When contigs are classified against NCBI, use the taxonomy graph to select the correct
translation table before calling ORFs in prokaryotic or viral sequences.

```julia
ncbi_taxonomy = Mycelia.load_ncbi_taxonomy()
tax_id = 562
table_id = Mycelia.get_ncbi_genetic_code(ncbi_taxonomy, tax_id)

orf_calls = Mycelia.run_pyrodigal(
    fasta_file=contig_fasta,
    out_dir="orf_calls",
    translation_table=table_id
)
```

Use `type=:mitochondrial` if you need the mitochondrial translation table.

## Functional Annotation

Functional annotation assigns biological functions to predicted genes.

### Functions

- Homology-based annotation
- Domain identification
- Pathway assignment

## Related Documentation

- See [Tutorial 6: Gene Annotation](../../generated/tutorials/06_gene_annotation.md) for practical examples
- See [API Reference](../../api-reference.md) for detailed function documentation
