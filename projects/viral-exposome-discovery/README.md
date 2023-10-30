# Viral Exposome

Andy's data location
/oak/stanford/scg/lab_mpsnyder/microbiome_brooks/
/oak/stanford/scg/lab_mpsnyder/microbiome_brooks/ultimagen/pool_2/kraken_reports

Mingming's data
/oak/stanford/projects/genomics/data/r84085_20231013_184512-GSSC-Snyder-MT-00000
/oak/stanford/projects/genomics/data/r64283e_20231026_225118-GSSC-Snyder-MT-00000
/labs/mpsnyder/share/exposome_data

conda create -n mmseqs2 -c bioconda mmseqs2

conda run --live-stream --no-capture-output -n mmseqs2 mmseqs databases NT $HOME/workspace/mmseqs/NT $HOME/workspace/mmseqs/tmp

mmseqs databases UniProtKB/Swiss-Prot $HOME/workspace/mmseqs $HOME/workspace/mmseqs/tmp

season - color scale
person - shape
location - annotation

pip install nbconvert
jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace 2.analyze-read-classifications.ipynb

jupyter nbconvert --clear-output 2.analyze-read-classifications.ipynb

 my_notebook.ipynb

kraken samples to rerun
SRR6399724
SRR6399725
SRR6399726
SRR6399727
SRR6399728
SRR6399729
SRR6399730
SRR6399731
SRR6399732
SRR6399773
SRR6399810
SRR6399900
SRR6399901
SRR6399902
SRR6399903
SRR6399905
SRR6399906
SRR6399907
SRR6399908
SRR6399909
SRR6399910
SRR6399911
SRR6399912
SRR6399913
SRR6399914
SRR6399915
SRR6399916
SRR6399917
SRR6399918
SRR6399919
SRR6399920
SRR6399921
SRR6399922
SRR6399923
SRR6399924
SRR6399925
SRR6399926
SRR6399927
SRR6399928
SRR6399929
SRR6399930
SRR6399931
SRR6399932
SRR6399933
SRR6399934
SRR6399935
SRR6399936
SRR6399937
SRR6399938
SRR6399939
SRR6399940
SRR6399941
SRR6399942
SRR6399943
SRR6399944
SRR6399945
SRR6399946
SRR6399947
SRR6399948
SRR6399949
SRR6399950
SRR6399951
SRR6399952
SRR6399953
SRR6399954
SRR6399955
SRR6399956
SRR6399957
SRR6399958
SRR6399959
SRR6399960
SRR6399961
SRR6399962
SRR6399963
SRR6399964
SRR6399965
SRR6399966
SRR6399967
SRR6399968
SRR6399969
SRR6399970
SRR6399971
SRR6399972
SRR6399973
SRR6399974
SRR6399975
SRR6399976
SRR6399977
SRR6399978
SRR6399979
SRR6399980
SRR6399981
SRR6399982
SRR6399983
SRR6399984
SRR6399985
SRR6399986
SRR6399987
SRR6399988
SRR6399989
SRR6399990
SRR6399991
SRR6399992
SRR6399993
SRR6399994
SRR6399995
SRR6399996
SRR6399997
SRR6399998
SRR6399999
SRR6400000
SRR6400001
SRR6400002
SRR6400003
SRR6400004
SRR6400005
SRR6400006
SRR6400007
SRR6400008
SRR6400009
SRR6400010
SRR6400011
SRR6400012
SRR6400013
SRR6400014
SRR6400015
SRR6400016
SRR6400017
SRR6400018
SRR6400019
SRR6400021
SRR6400022
SRR6400023
SRR6400024
SRR7365458
SRR7365460
SRR7365461
SRR7365462
SRR7365463
SRR7365464
SRR7365465
SRR7365466
SRR7365467
SRR7365468
SRR7365469
SRR7365470
SRR7365471
SRR7365473
SRR7365474
SRR7365475
SRR7365476
SRR7365477
SRR7365478
SRR7365479
SRR7365480
SRR7365481
SRR7365482
SRR7365483
SRR7365484
SRR7365485


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

results back and forth
```bash
rclone copy --progress --verbose $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data/results google_drive:Projects/viral-exposome-discovery/data/results

rclone copy --progress --verbose google_drive:Projects/viral-exposome-discovery/data/results $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data/results

```


Run me to backup (run the second version if we hit API throttling limits)
```bash
rclone copy --progress $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data google_drive:Projects/viral-exposome-discovery/data
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data google_drive:Projects/viral-exposome-discovery/data
```
copy from Google Drive back to local
```bash
# --exclude=*.{fq.gz,bam,fastq.gz}
# rclone copy --progress --exclude=*.{bam,fastq.gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 google_drive:Projects/viral-exposome-discovery/data $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data
rclone copy --progress --include=*.{fq.gz} --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 google_drive:Projects/viral-exposome-discovery/data/SRA $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data/SRA

rclone copy --progress google_drive:Projects/viral-exposome-discovery/data $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data

# --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1
rclone copy --progress --include=*{final.contigs.fastg.gfa.fna.blastn.nt.megablast.txt} --verbose  google_drive:Projects/viral-exposome-discovery/data $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data


```


```bash
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/Mycelia/projects/viral-exposome-discovery/data/SRA google_drive:Projects/viral-exposome-discovery/data/SRA

# workspace/Mycelia/projects/viral-exposome-discovery/data/SRA
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
