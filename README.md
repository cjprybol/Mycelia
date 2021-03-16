# [Eisenia](https://en.wikipedia.org/wiki/Eisenia_fetida), A pan-meta-omics graph framework

## Powered by JuliaGraphs, BioJulia, and Neo4J

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cameronprybol.gitlab.io/Eisenia.jl/dev) -->
<!-- [![Build Status](https://github.com/cjprybol/Eisenia.jl/badges/master/pipeline.svg)](https://github.com/cjprybol/Eisenia.jl/pipelines) -->
<!-- [![Coverage](https://github.com/cjprybol/Eisenia.jl/badges/master/coverage.svg)](https://github.com/cjprybol/Eisenia.jl/commits/master) -->
<!-- [![Build Status](https://ci.appveyor.com/api/projects/status/github/cjprybol/Eisenia.jl?svg=true)](https://ci.appveyor.com/project/cjprybol/Eisenia-jl) -->
<!-- [![Build Status](https://cloud.drone.io/api/badges/cjprybol/Eisenia.jl/status.svg)](https://cloud.drone.io/cjprybol/Eisenia.jl) -->
<!-- [![Coverage](https://codecov.io/gh/cjprybol/Eisenia.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cjprybol/Eisenia.jl) -->
<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

## Building a pan-meta-genome graph (multi-dataset)

![pan-meta-genome-graph](pan-meta-genome-graph.svg)

```
CREATE (dnamer:DNAmer {label: "DNAmer"})
CREATE (read:Read {label: "Read"})
CREATE (aamer:AAmer {label: "AAmer"})
CREATE (dataset:Dataset {label: "Dataset"})
CREATE (environment:Environment {label: "Environment"})
CREATE (genome:Genome {label: "Genome"})
CREATE (species:Species {label: "Species"})
CREATE (genus:Genus {label: "Genus"})
CREATE (family:Family {label: "Family"})
CREATE (order:Order {label: "Order"})
CREATE (subclass:Subclass {label: "Subclass"})
CREATE (class:Class {label: "Class"})
CREATE (subphylum:Subphylum {label: "Subphylum"})
CREATE (phylum:Phylum {label: "Phylum"})
CREATE (subkingdom:Subkingdom {label: "Subkingdom"})
CREATE (kingdom:Kingdom {label: "Kingdom"})
CREATE (clade:Clade {label: "Clade"})
CREATE (superkingdom:Superkingdom {label: "Superkingdom"})
CREATE (root:Root {label: "Root"})
CREATE (annotation:Annotation {label: "Annotation"})
CREATE (entity:Entity {label: "Entity"})
CREATE (pathway:Pathway {label: "Pathway"})

CREATE (dnamer)-[dnamer_annotation:CONTAINED_IN]->(annotation)
CREATE (aamer)-[aamer_annotation:CONTAINED_IN]->(annotation)
CREATE (genome)-[genome_annotation:CONTAINED_IN]->(annotation)
CREATE (superkingdom)-[superkingdom_root:PARENT]->(root)
CREATE (clade)-[clade_superkingdom:PARENT]->(superkingdom)
CREATE (kingdom)-[kingdom_clade:PARENT]->(clade)
CREATE (subkingdom)-[subkingdom_kingdom:PARENT]->(kingdom)
CREATE (phylum)-[phylum_subkingdom:PARENT]->(subkingdom)
CREATE (subphylum)-[subphylum_phylum:PARENT]->(phylum)
CREATE (class)-[class_subphylum:PARENT]->(subphylum)
CREATE (subclass)-[subclass_class:PARENT]->(class)
CREATE (order)-[order_subclass:PARENT]->(subclass)
CREATE (family)-[family_order:PARENT]->(order)
CREATE (genus)-[genus_family:PARENT]->(family)
CREATE (species)-[species_genus:PARENT]->(genus)
CREATE (entity)-[entity_species:ISA]->(species)
CREATE (dnamer)-[dnamer_genome:CONTAINED_IN]->(genome)
CREATE (dataset)-[dataset_environment:SOURCED_FROM]->(environment)
CREATE (read)-[read_dataset:CONTAINED_IN]->(dataset)
CREATE (aamer)-[aamer_read:CONTAINED_IN]->(read)
CREATE (dnamer)-[dnamer_read:CONTAINED_IN]->(read)
CREATE (dnamer)-[dnamer_aamer:TRANSLATES_TO]->(aamer)
CREATE (aamer)-[aamer_genome:CONTAINED_IN]->(genome)
CREATE (genome)-[genome_entity:GENOME_OF]->(entity)
CREATE (annotation)-[annotation_pathway:CONTAINED_IN]->(pathway)


RETURN dnamer,
read,
aamer,
dataset,
environment,
genome,
species,
genus,
family,
order,
subclass,
class,
subphylum,
phylum,
subkingdom,
kingdom,
clade,
superkingdom,
root,
annotation,
entity,
pathway,
annotation_pathway
```

## Building a metagenome graph (single dataset)

```
CREATE (dnamer:DNAmer {label: "DNAmer"})
CREATE (aamer:AAmer {label: "AAmer"})
CREATE (read:Read {label: "Read"})
CREATE (genome:Genome {label: "Genome"})
CREATE (annotation:Annotation {label: "Annotation"})
CREATE (pathway:Pathway {label: "Pathway"})

CREATE (dnamer)-[dnamer_annotation:CONTAINED_IN]->(annotation)
CREATE (aamer)-[aamer_annotation:CONTAINED_IN]->(annotation)
CREATE (genome)-[genome_annotation:CONTAINED_IN]->(annotation)
CREATE (annotation)-[annotation_pathway:CONTAINED_IN]->(pathway)
CREATE (dnamer)-[dnamer_genome:CONTAINED_IN]->(genome)
CREATE (read)-[read_dataset:CONTAINED_IN]->(dataset)
CREATE (aamer)-[aamer_read:CONTAINED_IN]->(read)
CREATE (dnamer)-[dnamer_read:CONTAINED_IN]->(read)
CREATE (dnamer)-[dnamer_aamer:TRANSLATES_TO]->(aamer)
CREATE (aamer)-[aamer_genome:CONTAINED_IN]->(genome)
CREATE (genome)-[genome_entity:GENOME_OF]->(entity)

RETURN dnamer,
read,
aamer,
dataset,
genome,
annotation,
pathway,
annotation_pathway
```

## Example Queries

English
```
CYPHER
```

- Give me the genome of entity x
- Give me the pangenome of species x
- Give me the pangenome of order x
- Give me all kmers in species x

- Finding spacers:
- Give me all kmers in species x matching pattern y
- Give me all datasets containing species x
- Give me all datasets containing pathway x
- Give me all reads supporting path/variant x



## Database

- Kmer tables
  - DNA+RNA: all primes from 7 <= x <= 63 (7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)
    - necessitates using BigKmers in design
  - AA: 3, 5, 7 (this means we need to include all 9mers, 15mers, and 21mers of DNA so that we have equivalencies between AAmers and DNAmers



## Creating probabilistic assemblies
1. Build weighted [de-bruijn graphs](https://en.wikipedia.org/wiki/De_Bruijn_graph) with observed data
2. Use the weighted de-bruijn graph as a [hidden markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) to [error correct observations](https://en.wikipedia.org/wiki/Viterbi_algorithm)
3. Use the error-corrected observations to build a new, more accurate weighted de-bruijn graph
4. Repeat 1-3 until convergence
5. Return maximum likelihood assembly

Consider for graph cleaning:
  - [x] Solve most likely path (slow but robust) vs
  - [ ] resample most likely paths (fast, less theoretically sound)
    - fairly certain that this is what the Flye assembler does


## Implementation requirements

- kmer spectra -> signal detection
- kmer size with signal -> initial assembly graph
    - bandage visualization
- initial assembly graph + reads -> iterative assembly
    - bandage visualization
- final assembly
    - bandage visualization
- annotation
    - find all start codons
    - find all stop codons
    - report in fasta format all coding regions
    - try and verify against prodigal to make sure that we aren’t missing any

- graph -> path -> sequence -> amino acid sequence

## API

- Need to allow users to add individual datasets to database
- does a sequence exist?
  - does a sequence exist allowing n mismatches
  - relative likelihood and frequency of most similar paths
    - this is just the viterbi traversal
- reads supporting traversal, and the datasets they came from
- connected components (species) containing path
- phylogenetic hierarchy of a component
  - nodes <-> kmers
- which datasets are a given kmer present in?
- export and import with GFA
- condense and expand between simplified and kmer graphs
- visualization 
  - with Bandage (works automatically with GFA)
  - cytoscape (needs to be built out)
