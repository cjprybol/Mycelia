```
# -name "*.sam"
find . -type f -exec du -h {} + | sort -rh | head -n 100
```

Acquire
- [x] ICTV reference viral
- [x] refseq viral
- [x] human genome
- [x] IMG/VR viral
    - copied using GLOBUS from JGI to NERSC, then rclone to Stanford Google Drive, Then rclone Google Drive to SCG3
    - rclone copy --progress $HOME/workspace/JGI stanford_viral_exposome:JGI
    - rclone copy --progress exposome:JGI $HOME/workspace/JGI
- [x] nt viral

map reads to each
- bwa-mem2
    - [x] ICTV reference viral
    - [x] refseq viral
    - [x] human genome
Got stuck with bwa-mem2 for larger fastas, unable to build bwa-mem2 indices w/ 1Tb of memory

Trying to subset the fastas with sourmash
- not working, limited to single core and running out of time (24 hours) on Slurm
- did work for IMG/VR k=11,13,17,23,31,41,53 but no nt_viral and none of the other k sizes
- [ ] IMG/VR sourmash on each short read sample

Switched to minimap2 which was more resource efficient
- minimap2
    - [x] nt viral
    - [x] IMG/VR viral
    - copied using GLOBUS from JGI to NERSC, then rclone to Stanford Google Drive, Then rclone Google Drive to SCG3
    - rclone copy --progress $HOME/workspace/JGI stanford_viral_exposome:JGI
    - rclone copy --progress exposome:JGI $HOME/workspace/JGI
    
- unable to build blast database with IMG/VR records due to IDs being >= 50 characters long, so did not bother blasting

- [ ] do variant calling of mapped reads against references
- [ ] do binned co-assembly w/ megabhit + graph-based variant-calling against reference w/ PGGB

associate taxa/variant abundance with other taxa/host/location/season/

# Viral Exposome

Andy's data location
/oak/stanford/scg/lab_mpsnyder/microbiome_brooks/
/oak/stanford/scg/lab_mpsnyder/microbiome_brooks/ultimagen/pool_2/kraken_reports

results back and forth
```bash
rclone copy --progress --verbose $HOME/workspace/Mycelia/projects/viral-exposome/data/results snyder-virome:viral-exposome/data/results
rclone copy --progress --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA snyder-virome:viral-exposome/data/SRA
# rclone copy --progress --verbose $HOME/workspace/Mycelia/projects/viral-exposome/data/exposome_data snyder-virome:viral-exposome/data/exposome_data
# rclone copy --progress --verbose google_drive:Projects/viral-exposome/data/results $HOME/workspace/Mycelia/projects/viral-exposome/data/results
```

Run me to restore
```bash
rclone copy --progress snyder-virome:viral-exposome/data $HOME/workspace/Mycelia/projects/viral-exposome/data
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress snyder-virome:viral-exposome/data $HOME/workspace/Mycelia/projects/viral-exposome/data
```

Run me to backup (run the second version if we hit API throttling limits)
```bash
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/Mycelia/projects/viral-exposome/data snyder-virome:viral-exposome/data

rclone copy --progress --include=*{.gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/blastdb snyder-virome:blastdb
rclone copy --progress --include=*{.mmi} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/blastdb snyder-virome:blastdb
```


copy from Google Drive back to local
```bash
# --exclude=*.{fq.gz,bam,fastq.gz}
# rclone copy --progress --exclude=*.{bam,fastq.gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 google_drive:Projects/viral-exposome/data $HOME/workspace/Mycelia/projects/viral-exposome/data
rclone copy --progress --include=*.{fq.gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 google_drive:Projects/viral-exposome/data/SRA $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA

rclone copy --progress --include=*.{fq.gz} --verbose google_drive:Projects/viral-exposome/data/SRA $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA

rclone copy --progress google_drive:Projects/viral-exposome/data $HOME/workspace/Mycelia/projects/viral-exposome/data

# --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1
rclone copy --progress --include=*{final.contigs.fastg.gfa.fna.blastn.nt.megablast.txt} --verbose  google_drive:Projects/viral-exposome/data $HOME/workspace/Mycelia/projects/viral-exposome/data
```








rclone copy --progress --include=*.{gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA snyder-virome:viral-exposome/data/SRA

rclone copy --progress --include=*.{sam.gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA snyder-virome:viral-exposome/data/SRA


```bash
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA google_drive:Projects/viral-exposome/data/SRA

# workspace/Mycelia/projects/viral-exposome/data/SRA
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
