# References and Citations

This page tracks citation guidance for Mycelia and the third-party tools it wraps.
If a workflow uses external tools or datasets, cite both Mycelia and the tools.

## Citing Mycelia

Mycelia does not yet ship a formal paper reference. Until then:
- Cite the GitHub repository and the specific release tag or commit used.
- If you mint a DOI for a release, add it here and cite the DOI in publications.

## Third-Party Tools (cite when used)

Below is a non-exhaustive list of external tools wrapped by Mycelia. Use the
recommended citation from each tool's documentation or manuscript.

- Assemblers: MEGAHIT, metaSPAdes, SPAdes, Flye/metaFlye, hifiasm, Canu, SKESA, Unicycler
- Assembly pipelines: Autocycler (https://github.com/rrwick/Autocycler)
- Annotation: Prodigal, Pyrodigal, Prodigal-gv, Augustus, MetaEuk, BLAST+, MMSeqs2,
  TransTerm, tRNAscan-SE, MLST
- QC and validation: fastp, Filtlong, Trim Galore, QUAST, BUSCO, CheckM/CheckM2, MUMmer
- Alignment and mapping: minimap2, samtools, Qualimap, BWA
- Variant calling and evaluation: GATK, FreeBayes, Clair3, BCFtools, RTG vcfeval
- Comparative genomics: Mash, PGGB, cactus, odgi, PanTools (https://github.com/BioinformaticsLab/pantools)
- Graph construction and k-mer utilities: BCALM (https://github.com/GATB/bcalm),
  GGCAT (https://github.com/algbio/ggcat), Prokrustean (https://github.com/KoslickiLab/prokrustean)
- Read simulation: ART, PBSim, DeepSimulator, Badread
- Structure search and clustering: Foldseek (https://github.com/steineggerlab/foldseek)
- Taxonomy and profiling: MetaPhlAn, Metabuli, Sylph, Kraken-style tools (planned)
- Binning: MetaBAT2, VAMB, TaxVAMB, Taxometer, MetaCoAG, GenomeFace, COMEBin, dRep, MAGmax
- Data access: NCBI datasets, Entrez, SRA Toolkit
- Infrastructure: SLURM, rclone, Bioconda/conda

## Reduced Amino Acid Alphabet Citations

These references correspond to the reduced alphabet schemes used in `src/constants.jl`.

- Wang J., Wang W. (1999). A computational approach to simplifying the protein folding alphabet.
  *Nature Structural Biology*, 6(11), 1033-1038. https://doi.org/10.1038/14918 (WW5)
- Melo F., Marti-Renom M.A. (2006). Accuracy of sequence alignment and fold assessment using reduced amino acid alphabets.
  *Proteins*, 63(4), 986-995. https://doi.org/10.1002/prot.20976 (MM5)
