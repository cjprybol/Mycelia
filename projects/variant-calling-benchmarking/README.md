# variant-calling-benchmarking

The goal of this work is:
1. get my probabilistic assembler up and running again
2. Compare the accuracy of graph-based variant calling directly from assemblies against short and long read variant callers

real genomes:
- [x] smallest refseq viral, bacterial, and archaeal genomes
- 20240121.identify-small-reference-genomes.ipynb

Simulated chromosomes:
- [x] 10kb, 100kb, 1Mb

producing variants
- [x] substitutions, insertion/deletions, inversions, translocations, and duplications <= 10% of chromosome length
- https://biojulia.dev/GeneticVariation.jl/stable/man/io/vcf-bcf.html
- https://github.com/rasmushenningsson/VariantCallFormat.jl
- https://github.com/OpenMendel/VCFTools.jl

simulating reads
- 10, 100, 1000x
- paired-end short reads
    - https://www.niehs.nih.gov/research/resources/software/biostatistics/art w/ default error rates
- single-end long reads
    - https://github.com/rrwick/Badread w/ 10% error rates
    - https://github.com/bcgsc/NanoSim w/ 1% error rates

assembly
- short read only assembly
    - Mycelia
    - MegaHIT
    - spades
- long read only assembly
    - Mycelia
    - flye
    - raven https://github.com/lbcb-sci/raven
- hybrid-assembly
    - Mycelia
    - https://github.com/rrwick/Polypolish

graph-based variant calling
- vg
- pggb

read mapping
- minimap2

variant calling
- https://github.com/google/deepvariant
- bcftools https://samtools.github.io/bcftools/howtos/variant-calling.html
- Freebayes

short-only variant calling
- gatk

long-only variant calling
- sniffles
- https://github.com/WGLab/NanoCaller
- Clair3

results, analysis, discussion
- analyze accuracy of Mycelia against all of the other approaches individually
- assess Mycelia against hybrid and ensemble approaches

references
- https://www.nejm.org/doi/10.1056/NEJMc2112090
- https://www.frontiersin.org/articles/10.3389/fgene.2022.887644/full
- https://ieeexplore.ieee.org/document/18626

Considered but not tested:
- https://github.com/kishwarshafin/pepper