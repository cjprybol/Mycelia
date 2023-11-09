# Mycelia

! In development !

Multiomic analysis and data integration for biological characterization, prediction, and  design.

Mostly just a wrapper around existing best-in-class bioinformatics tools with just enough Julia code to glue them together and enable downsteam analytics.

Designed for (academic and goverment) HPC and (GCP & AWS) cloud systems.

## Install

Install Julia (if not already installed)
```bash
# replacing "release" with "lts" is a good alternate for deployments to production
curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel release
# add Julia to your $PATH by reloading the appropriate source file, as instructed by the installer
source ~/.bashrc
```

I have had trouble getting the visualization libraries Plots.jl and Makie.jl (and associated packages) to load correctly on Stanford's SCG3 due to the complexity of the native LD_LIBRARY_PATH

I imagine other research supercomputer users may have similar issues, although I don't have these issues on cloud vendors like GCP or AWS

To enable Julia to install all of it's own necessary dependencies independent of the system, I needed to reset the LD_LIBRARY_PATH variable prior to launching Julia

This can be done easily when launching Julia from the command line by 
```bash
export LD_LIBRARY_PATH="" && julia
```

[And can be done for Julia jupyter kernels by setting the `env` key => value pair in the appropriate kernel.json file](https://stackoverflow.com/a/53595397)


Clone the repo directly
```
cd /path/where/you/want/the/repo
# for production usage
git clone https://github.com/cjprybol/Mycelia.git
# for development
git clone git@github.com:cjprybol/Mycelia.git
```

Or as Julia package
```
import Pkg
# for production usage
Pkg.add(url="https://github.com/cjprybol/Mycelia.git")
# for development
Pkg.develop(url="git@github.com:cjprybol/Mycelia.git")
```

## TODO: re-enable the command line API

See the docs [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://research.cjp.garden/Mycelia/)


<!-- [![Build Status](https://github.com/cjprybol/Mycelia/badges/master/pipeline.svg)](https://github.com/cjprybol/Mycelia.jl/pipelines) -->
<!-- [![Coverage](https://github.com/cjprybol/Mycelia/badges/master/coverage.svg)](https://github.com/cjprybol/Mycelia.jl/commits/master) -->
<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

## Related & Relevant Projects & Companies

- https://github.com/merenlab/anvio
- https://github.com/SciGraph/SciGraph
- https://github.com/geneontology/obographs
- https://github.com/biolink/ontobio
- https://github.com/biolink/kgx
- https://github.com/monarch-initiative/biolink-api
- https://github.com/monarch-initiative/embiggen
- https://github.com/AnacletoLAB/grape
- https://github.com/linkml/linkml

- pangenomes
  - [PanTools](https://www.bioinformatics.nl/pangenomics/manual/) ([paper](https://pubmed.ncbi.nlm.nih.gov/27587666/))
  - [BioJulia/GenomeGraphs.jl](https://github.com/BioJulia/GenomeGraphs.jl)
  - [JCVI PanGenomePipeline](https://github.com/JCVenterInstitute/PanGenomePipeline)
  - [Roary](https://github.com/sanger-pathogens/Roary) - deprecated :(
  - [vg](https://github.com/vgteam/vg)
  - [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN)
  - [optimized dynamic genome graph implementation](https://github.com/pangenome/odgi)
  - [metagraph](https://github.com/ratschlab/metagraph)
  - [Other softwae mentioned within this pangenome graphs review](https://doi.org/10.1146/annurev-genom-120219-080406)
- genome machine learning
  - https://genoml.com/
  - https://github.com/davek44/Basset
  - https://github.com/calico/basenji
  - https://github.com/kundajelab/dragonn
  - https://github.com/kundajelab/deeplift

Related platforms:
- genomes and organisms:
  - [TeselaGen](https://teselagen.com/)
  - [CELLO](https://github.com/CIDARLAB/cello)
  - [Lattice Automation](https://www.latticeautomation.com/)
  - [Asimov](https://www.asimov.com/)
  - [Infobiotics Workbench](https://github.com/Infobiotics/ibw)
  - [iBioSim](https://async.ece.utah.edu/tools/ibiosim/)
  - [TinkerCell](http://www.tinkercell.com/)
  - [Cyrus Bio](https://cyrusbio.com/)
  - [Gingko](https://www.ginkgobioworks.com/)
- Proteins:
  - [Big Hat Biosciences](https://www.bighatbio.com/)
  - [OpenFold](https://openfold.io/)
  - [Rosetta](https://www.rosettacommons.org/software)
  - [Arzeda](https://www.arzeda.com/)
  - [DeepChain](https://deepchain.bio/)

Semantic models & ontologies:
- [Synthetic Biology Open Language](https://sbolstandard.org/)
- [Systems Biology Markup Language](https://sbml.org/)
- [Systems Biology Graphical Notation](https://sbgn.github.io/)
- [BioOntologies](https://bioportal.bioontology.org/ontologies)
- [Open Biological and Biomedical Ontology](https://obofoundry.org/)

Model repositories:
- [BioModels](https://www.ebi.ac.uk/biomodels/)
- [GLAMM](https://glamm.lbl.gov/)

communities
- [AI 4 SynBio](https://www.ai4synbio.org/)
