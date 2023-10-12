# Viral Exposome

## Figures
- species, genus, family, order, class, phylum, kingdom
    - Kraken
    - mmseqs
        - UniRef50
        - UniRef90
        - UniRef100
    - Blast NT

- which groups have the most divergence
- sample by sample clusters
- kmer saturation diversity of proteins and dna
- diveristy gain over known databases for protein and dna

Run me to backup (run the second version if we hit API throttling limits)
```bash
rclone copy --progress $HOME/workspace/Mycelia/projects/viral-pangenome-discovery/data google_drive:Projects/viral-pangenome-discovery/data
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/Mycelia/projects/viral-pangenome-discovery/data google_drive:Projects/viral-pangenome-discovery/data
```
copy from Google Drive back to local
```bash
rclone copy --progress --exclude=*.{fastq.gz,fq.gz,bam} google_drive:Projects/viral-pangenome-discovery/data $HOME/workspace/Mycelia/projects/viral-pangenome-discovery/data
rclone copy --progress google_drive:Projects/viral-pangenome-discovery/data $HOME/workspace/Mycelia/projects/viral-pangenome-discovery/data
```

copy Globus-JGI reference DB into cloud for shuffling back and forth to non-globus machines
note - mmseqs db is already in the drive!
```bash
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/JGI google_drive:Projects/reference-databases
```

```bash
rclone copy --progress google_drive:Projects/reference-databases/mmseqs.tar.gz $HOME/workspace/mmseqs_alt
```

## Classification
All of the samples have been classified using 3 algorithms, one of which used 3 databases, for a total of 5 independent classifications for each assembled sequence:
- DNA
    - complete:
        - blastn:
            - ref_viruses_rep_genomes w/ dc-megablast algorithm
            - ref_viruses_rep_genomes w/ blastn algorithm
    - in progress:
        - mmseqs:
            - ICTV
- protein
    - complete:
        - geNomad
        - virsorter2
        - mmseqs:
            UniRef100
            UniRef90
            UniRef50

## Quantification of Novelty
All contigs that have been identified as viral by at least 4 of the 7 above classification approaches are then assessed for genomic and proteomic novelty by comparing against:

- DNA
    - complete
        - NCBI complete viral database
        - nt_viruses w/ megablast algorithm
    - in progress
        - kmers
        - JGI's IMG/VR
            - data acquired by Globus
            - complete variant
            - high confidence variant

- Protein
    - complete
    - in progress
        - UniRef50

Multiples paths to discovery:
- classify first
    - classify reads to be of viral origin
        - assemble by sample and then add to pangenome
        - I didn't trust the kraken results on these, so dropped the short-read based classification
            - could still be great for long read classification, since long reads are effectively assembled contigs
- assemble first
    - assemble by dataset & classify assembled contigs
^ compare against eachother and then to literature for best practices
- read classification for short reads using Kraken had very low concordance with the assembled contigs

- Annotate proteins of viral contigs to search for novel proteins
    - UniRef100

Done (performed on LCFTA & Saturn clusters, current code may not run on HPC):
- Library QC
- Library Assembly
- Blast-based DNA classification

Acknowledgements:
- Mike
- Mingming
- Aeron
- NERSC
- SCG (SCGPM)
- NERSC
- Lawrencium (Lawrence Berkeley National Lab)