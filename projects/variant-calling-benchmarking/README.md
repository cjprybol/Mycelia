# variant-calling-benchmarking

The goal of this work is:
1. get my probabilistic assembler up and running again
2. Compare the accuracy of graph-based variant calling directly from assemblies against short and long read variant callers

## Data Generation
- [x] real genomes:
    - smallest refseq viral, bacterial, and archaeal genomes
    - 20240121.identify-small-reference-genomes.ipynb
- [x] Simulated chromosomes:
    - 10kb, 100kb, 1Mb
- [x] producing variants
    - substitutions, insertion/deletions, inversions, translocations, and duplications <= 10% of chromosome length
    - https://biojulia.dev/GeneticVariation.jl/stable/man/io/vcf-bcf.html
    - https://github.com/rasmushenningsson/VariantCallFormat.jl
    - https://github.com/OpenMendel/VCFTools.jl
- [x] simulating reads
    - 10, 100, 1000x
    - paired-end short reads
        - https://www.niehs.nih.gov/research/resources/software/biostatistics/art
    - single-end long reads
        - https://github.com/bcgsc/NanoSim
            - couldn't get to run
            - nanosim-h fork didn't seem to have recent error profiles
            - skipped
        - https://github.com/rrwick/Badread
- [x] qc-filtering reads
    - filtlong for long reads
    - trimgalore for short reads

From here, reads take two paths:
1. assembly-graph variant calling
2. mapping-based variant calling

## assembly
- [ ] short read only assembly
    - Mycelia
    - MegaHIT
    - spades --isolate
- [ ] long read only assembly
    - Mycelia
    - metaflye
    - raven https://github.com/lbcb-sci/raven
- [ ] hybrid-assembly
    - Mycelia
    - https://github.com/rrwick/Polypolish
    - hybrid metaspades
    
## assembly-graph variant calling
- [ ] graph-based variant calling
    - vg
    - pggb
    
## Mapping
- [x] read mapping
    - minimap2

## mapping-based variant-calling
- [ ] long and short variant calling
    - https://github.com/google/deepvariant
    - bcftools https://samtools.github.io/bcftools/howtos/variant-calling.html
    - Freebayes
- [ ] short-read variant calling
    - gatk
- [ ] long-read variant calling
    - nanovar
    - sniffles
    - https://github.com/WGLab/NanoCaller
    - Clair3

## results, analysis, discussion
- [ ] analyze accuracy of Mycelia against all of the other approaches individually
- [ ] assess Mycelia against hybrid and ensemble approaches

references
- https://www.nejm.org/doi/10.1056/NEJMc2112090
- https://www.frontiersin.org/articles/10.3389/fgene.2022.887644/full
- https://ieeexplore.ieee.org/document/18626

Considered but not tested:
- https://github.com/kishwarshafin/pepper