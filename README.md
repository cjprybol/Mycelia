# Mycelia

Multiomics knowledge graphs for [characterization](), [discovery](), and [biological design]()

The goal is to enable passing a description of a genome's fiction and then have it generate the sequence. Like [DALLE-E2](https://openai.com/dall-e-2/), but for genomes instead of pictures.

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://research.cjp.garden/Mycelia/)
[![Build Status](https://github.com/cjprybol/Mycelia/badges/master/pipeline.svg)](https://github.com/cjprybol/Mycelia.jl/pipelines)
[![Coverage](https://github.com/cjprybol/Mycelia/badges/master/coverage.svg)](https://github.com/cjprybol/Mycelia.jl/commits/master)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

## Development and contributing

1. branch/fork this repository and develop in github codespaces
  - I like to work in the [jupyter-datascience](https://github.com/cjprybol/Mycelia/tree/master/.devcontainer/jupyter-datascience) devcontainer
1. clone into `/your/path/to/Mycelia` and load the julia package into your Julia Project.toml
  - `pkg> dev /your/path/to/Mycelia`

Related Software
- pangenomes
  - [PanTools](https://www.bioinformatics.nl/pangenomics/manual/) ([paper](https://pubmed.ncbi.nlm.nih.gov/27587666/))
  - [BioJulia/GenomeGraphs.jl](https://github.com/BioJulia/GenomeGraphs.jl)
  - [JCVI PanGenomePipeline](https://github.com/JCVenterInstitute/PanGenomePipeline)
  - [Roary](https://github.com/sanger-pathogens/Roary) - deprecated :(
  - [vg](https://github.com/vgteam/vg)
  - [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN)
  - [optimized dynamic genome graph implementation](https://github.com/pangenome/odgi)
  - [Other softwae mentioned within this pangenome graphs review](https://doi.org/10.1146/annurev-genom-120219-080406)
- genome machine learning
  - https://genoml.com/
  - https://github.com/davek44/Basset
  - https://github.com/calico/basenji
  - https://github.com/kundajelab/dragonn
  - https://github.com/kundajelab/deeplift

Related platforms:
- [TeselaGen](https://teselagen.com/)
- [CELLO](https://github.com/CIDARLAB/cello)
- [Lattice Automation](https://www.latticeautomation.com/)
- [Asimov](https://www.asimov.com/)
- [Infobiotics Workbench](https://github.com/Infobiotics/ibw)
- [iBioSim](https://async.ece.utah.edu/tools/ibiosim/)
- [TinkerCell](http://www.tinkercell.com/)

Semantic models & ontologies:
- [Synthetic Biology Open Language](https://sbolstandard.org/)
- [Systems Biology Markup Language](https://sbml.org/)
- [Systems Biology Graphical Notation](https://sbgn.github.io/)
- [BioOntologies](https://bioportal.bioontology.org/ontologies)
- [Open Biological and Biomedical Ontology (Foundry](https://obofoundry.org/)

Model repositories:
- [BioModels](https://www.ebi.ac.uk/biomodels/)
- [GLAMM](https://glamm.lbl.gov/)

communities
- [AI 4 SynBio](https://www.ai4synbio.org/)
