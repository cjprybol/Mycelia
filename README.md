# Mycelia

![](banner-logo.jpg)

Multiomic analysis and data integration for biological characterization, prediction, and  design.

Mycelia wraps existing best-in-class bioinformatics tools via Conda where existing solutions are available, and extends and integrates those tools with code written in Julia.

Designed for linux-based HPC and cloud systems.

## Install

[Install Julia (if not already installed)](https://github.com/JuliaLang/juliaup)

I have had trouble getting the visualization libraries Plots.jl and Makie.jl (and associated packages) to load correctly on HPC due to the complexities of the default LD_LIBRARY_PATH

I imagine other research supercomputer users may have similar issues, although I don't have these issues on cloud vendors like GCP or AWS

To enable Julia to install all of it's own necessary dependencies independent of the system, I reset the LD_LIBRARY_PATH variable prior to launching Julia !!

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

documentation in prep

https://doi.org/10.1101/2024.05.29.596541
