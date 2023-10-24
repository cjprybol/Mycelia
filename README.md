# Mycelia

setup on SCG3
```bash
julia
```
then
```julia
ENV["LD_LIBRARY_PATH"] = ""
import Pkg
Pkg.add("IJulia")
Pkg.build("IJulia")
exit()
```

! In development !

Multiomic analysis and data integration for [characterization](), [discovery](), [prediction](), and [biological design]().

## Install

Install Julia (if not already installed)
```bash
curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel release
# add Julia to your $PATH by reloading the appropriate source file, as instructed by the installer
source ~/.bashrc
```

Install Mycelia & precompile to validate installation
```
export LD_LIBRARY_PATH="" && julia -e 'import Pkg; Pkg.add("IJulia"); Pkg.develop(url="https://github.com/cjprybol/Mycelia.git"); Pkg.precompile(); import Mycelia'
```

Add Mycelia to your system path
```
export PATH="$HOME/.julia/"
```

See the docs [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://research.cjp.garden/Mycelia/)


<!-- [![Build Status](https://github.com/cjprybol/Mycelia/badges/master/pipeline.svg)](https://github.com/cjprybol/Mycelia.jl/pipelines) -->
<!-- [![Coverage](https://github.com/cjprybol/Mycelia/badges/master/coverage.svg)](https://github.com/cjprybol/Mycelia.jl/commits/master) -->
<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

## Development and contributing

1. branch/fork this repository and develop in github codespaces
1. clone into `/your/path/to/Mycelia` and load the julia package into your Julia Project.toml
  - `pkg> dev /your/path/to/Mycelia`

## Related & Relevant Projects & Companies

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
