# Mycelia

biological knowledge graphs for "omics" powered by [JuliaGraphs](https://github.com/JuliaGraphs), [BioJulia](https://github.com/BioJulia), and [Neo4J](https://neo4j.com/)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://research.cjp.garden/Mycelia/)
<!-- [![Build Status](https://github.com/cjprybol/Mycelia.jl/badges/master/pipeline.svg)](https://github.com/cjprybol/Mycelia.jl/pipelines) -->
<!-- [![Coverage](https://github.com/cjprybol/Mycelia.jl/badges/master/coverage.svg)](https://github.com/cjprybol/Mycelia.jl/commits/master) -->
<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

## Install
```julia
pkg> add git@github.com:cjprybol/Mycelia.git
```

## Related Software:
- [PanTools](https://www.bioinformatics.nl/pangenomics/manual/) ([paper](https://pubmed.ncbi.nlm.nih.gov/27587666/))
- [BioJulia/GenomeGraphs.jl](https://github.com/BioJulia/GenomeGraphs.jl)
- [Pangenome Graphs Review](https://doi.org/10.1146/annurev-genom-120219-080406)

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
- condense kmer graphs to simplified graphs
- visualization 
  - with Bandage (works automatically with GFA)
  - with Neo4J visualization library

