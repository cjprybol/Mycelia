# Mycelia

:warning: **pre-alpha!** This software is in the process of being brought to production-readiness

Biological knowledge graphs for "omics" powered by [JuliaGraphs](https://github.com/JuliaGraphs), [BioJulia](https://github.com/BioJulia), and [Neo4J](https://neo4j.com/). Goals of the project include:
- dynamic creation of arbitrary pangenomes with metadata from NCBI and other sources
- characterization of pangenomes
- using pangenome graphs to inform bioengineering of organisms and biological communities (bioremediation and clinical microbiome therapies)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://research.cjp.garden/Mycelia/)
[![Build Status](https://github.com/cjprybol/Mycelia/badges/master/pipeline.svg)](https://github.com/cjprybol/Mycelia.jl/pipelines)
[![Coverage](https://github.com/cjprybol/Mycelia/badges/master/coverage.svg)](https://github.com/cjprybol/Mycelia.jl/commits/master)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

## Install for local development

1. clone into `/your/path/to/Mycelia` 

```julia
pkg> dev /your/path/to/Mycelia
```

## Related Software
- [PanTools](https://www.bioinformatics.nl/pangenomics/manual/) ([paper](https://pubmed.ncbi.nlm.nih.gov/27587666/))
- [BioJulia/GenomeGraphs.jl](https://github.com/BioJulia/GenomeGraphs.jl)
- [JCVI PanGenomePipeline](https://github.com/JCVenterInstitute/PanGenomePipeline)
- [Pangenome Graphs Review](https://doi.org/10.1146/annurev-genom-120219-080406)
- [Roary](https://github.com/sanger-pathogens/Roary) - deprecated :(
- [vg](https://github.com/vgteam/vg)
- [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN)
- [optimized dynamic genome graph implementation](https://github.com/pangenome/odgi)
