# [Eisenia](https://en.wikipedia.org/wiki/Eisenia_fetida), A meta-pan-omics graph framework

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cameronprybol.gitlab.io/Eisenia.jl/dev) -->
<!-- [![Build Status](https://github.com/cjprybol/Eisenia.jl/badges/master/pipeline.svg)](https://github.com/cjprybol/Eisenia.jl/pipelines) -->
<!-- [![Coverage](https://github.com/cjprybol/Eisenia.jl/badges/master/coverage.svg)](https://github.com/cjprybol/Eisenia.jl/commits/master) -->
<!-- [![Build Status](https://ci.appveyor.com/api/projects/status/github/cjprybol/Eisenia.jl?svg=true)](https://ci.appveyor.com/project/cjprybol/Eisenia-jl) -->
<!-- [![Build Status](https://cloud.drone.io/api/badges/cjprybol/Eisenia.jl/status.svg)](https://cloud.drone.io/cjprybol/Eisenia.jl) -->
<!-- [![Coverage](https://codecov.io/gh/cjprybol/Eisenia.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cjprybol/Eisenia.jl) -->
<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

## Creating probabilistic assemblies
1. Build weighted [de-bruijn graphs](https://en.wikipedia.org/wiki/De_Bruijn_graph) with observed data
2. Use the weighted de-bruijn graph as a [hidden markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) to [error correct observations](https://en.wikipedia.org/wiki/Viterbi_algorithm)
3. Use the error-corrected observations to build a new, more accurate weighted de-bruijn graph
4. Repeat 1-3 until convergence
5. Return maximum likelihood assembly

## Database

### Kmer Nodes
- Kmer tables
  - DNA+RNA: all primes from 7 <= x <= 63 (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)
    - necessitates using BigKmers in design
  - AA: 3, 5, 7

| hash | sequence |
|------|----------|

### Edges

- edge tables with evidence

| hash | orientation | hash | orientation | sequence ID | sequence index | reference_or_observation |
|------|-------------|------|-------------|-------------|----------------|--------------------------|

### Reference Sequences

- Reference sequences

| hash | id | description | sequence |
|------|----|-------------|----------|

- Reference sequence paths

| id | k-length | path |
|----|----------|------|

- Kmer to reference paths

| hash (kmer) | hash (reference path) | index (reference path) |
|-------------|-----------------------|------------------------|

- reference sequence kmer spectra

| id | kmer | count |
|----|------|-------|

- annotations
  - these will be re-interpreted into the graph as needed to interconvert between linear and graph contexts

| reference sequence ID | start | stop | strand | other stuff in GFF files? | is_protein_coding |
|-----------------------|-------|------|--------|---------------------------|-------------------|

### Datasets

- kmer to dataset

| kmer | dataset | sequence id | index |
|------|---------|-------------|-------|

- dataset kmer spectra

| id | kmer | count |
|----|------|-------|

- Dataset metadata

| joint hash of fastq files (ID) | id | description | lat | long | country | state | city | zip code | source type |
|--------------------------------|----|-------------|-----|------|---------|-------|------|----------|-------------|

### External data

- Dataset by dataset analysis
  - HDF5
    - raw fastq
    - kmer spectra primes 7 - 63
    - correction stage k=7
    - correction stage k=...
    - correction stage k=61

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
