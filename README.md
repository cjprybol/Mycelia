# [Eisenia](https://en.wikipedia.org/wiki/Eisenia_fetida)

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cameronprybol.gitlab.io/Eisenia.jl/dev) -->
<!-- [![Build Status](https://github.com/cjprybol/Eisenia.jl/badges/master/pipeline.svg)](https://github.com/cjprybol/Eisenia.jl/pipelines) -->
<!-- [![Coverage](https://github.com/cjprybol/Eisenia.jl/badges/master/coverage.svg)](https://github.com/cjprybol/Eisenia.jl/commits/master) -->
<!-- [![Build Status](https://ci.appveyor.com/api/projects/status/github/cjprybol/Eisenia.jl?svg=true)](https://ci.appveyor.com/project/cjprybol/Eisenia-jl) -->
<!-- [![Build Status](https://cloud.drone.io/api/badges/cjprybol/Eisenia.jl/status.svg)](https://cloud.drone.io/cjprybol/Eisenia.jl) -->
<!-- [![Coverage](https://codecov.io/gh/cjprybol/Eisenia.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cjprybol/Eisenia.jl) -->
<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

Probablistic Pan-(gen/transcript/prote)omics

- Build weighted [de-bruijn graphs](https://en.wikipedia.org/wiki/De_Bruijn_graph) with observed data
- Use the weighted de-bruijn graph as a [hidden markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) to [error correct observations](https://en.wikipedia.org/wiki/Viterbi_algorithm)
- Use the error-corrected observations to build a new, more accurate weighted de-bruijn graph
- Repeat until convergence
- Return maximum likelihood assembly