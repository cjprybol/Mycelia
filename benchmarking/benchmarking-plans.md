## Databases

- plasmids
    - https://ccb-microbe.cs.uni-saarland.de/plsdb2025/ (72K)
        - https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download
    - https://plasmid.deepomics.org/database/ (852K)
        - https://plasmid.deepomics.org/download 
    - https://www.ncbi.nlm.nih.gov/genome/plasmids/
- environmental viruses
    - IMG/VR5
- prakaryotes
    - GTDB
- proteins
    - UniProt

- NCBI blast database
```julia
julia> show(sort(Mycelia.list_blastdbs(), "SIZE (GB)"), allrows=true, allcols=true)
38×4 DataFrame
 Row │ NAME                     DESCRIPTION                        SIZE (GB)  LAST_UPDATED 
     │ String                   String                             Float64    Dates.Date   
─────┼─────────────────────────────────────────────────────────────────────────────────────
   1 │ ref_viroids_rep_genomes  Refseq viroids reference genomes      0.0001  2025-12-15
   2 │ 18S_fungal_sequences     18S ribosomal RNA sequences (SSU…     0.0026  2025-12-16
   3 │ LSU_prokaryote_rRNA      Large subunit ribosomal nucleic …     0.0041  2022-12-05
   4 │ LSU_eukaryote_rRNA       Large subunit ribosomal nucleic …     0.0053  2022-12-05
   5 │ SSU_eukaryote_rRNA       Small subunit ribosomal nucleic …     0.0063  2022-12-05
   6 │ 28S_fungal_sequences     28S ribosomal RNA sequences (LSU…     0.0066  2025-12-16
   7 │ ITS_RefSeq_Fungi         Internal transcribed spacer regi…     0.0087  2025-12-12
   8 │ pdbnt                    PDB nucleotide database               0.0142  2025-12-15
   9 │ 16S_ribosomal_RNA        16S ribosomal RNA (Bacteria and …     0.0182  2025-12-16
  10 │ ITS_eukaryote_sequences  ITS eukaryote BLAST                   0.0344  2025-12-17
  11 │ refseq_select_rna        RefSeq Select RNA sequences           0.0681  2025-12-15
  12 │ ref_viruses_rep_genomes  Refseq viruses reference genomes      0.144   2025-12-14
  13 │ mito                     NCBI Genomic Mitochondrial Refer…     0.154   2025-12-14
  14 │ pdbaa                    PDB protein database                  0.2396  2025-12-15
  15 │ taxdb                    Taxonomy database                     0.2788  2025-12-15
  16 │ env_nt                   environmental samples                 0.3298  2025-10-29
  17 │ swissprot                Non-redundant UniProtKB/SwissPro…     0.3639  2025-12-17
  18 │ landmark                 Landmark database for SmartBLAST      0.4693  2025-08-06
  19 │ nt_others                Artificial and other seqs nt          0.7596  2025-12-02
  20 │ human_genome             Homo sapiens GRCh38.p13 [GCF_000…     1.1572  2021-06-02
  21 │ pataa                    Protein sequences derived from t…     2.1798  2025-12-15
  22 │ mouse_genome             Mus musculus GRCm39 [GCF_0000016…     3.6543  2021-06-02
  23 │ env_nr                   Proteins from WGS metagenomic pr…     4.4633  2025-12-15
  24 │ tsa_nr                   Transcriptome Shotgun Assembly (…     6.0917  2025-12-15
  25 │ tsa_nt                   Transcriptome Shotgun Assembly (…     6.3489  2025-09-20
  26 │ patnt                    Nucleotide sequences derived fro…    17.5441  2025-12-13
  27 │ ref_prok_rep_genomes     Refseq prokaryote reference geno…    25.078   2025-12-14
  28 │ refseq_select_prot       RefSeq Select proteins               44.0842  2025-12-14
  29 │ nt_viruses               Viruses nt                           67.7264  2025-12-15
  30 │ refseq_rna               NCBI Transcript Reference Sequen…    68.7504  2025-12-17
  31 │ Betacoronavirus          Betacoronavirus                      71.0289  2025-06-26
  32 │ nt_prok                  Prokaryota (bacteria and archaea…    89.1318  2025-12-15
  33 │ refseq_protein           NCBI Protein Reference Sequences    243.559   2025-12-16
  34 │ core_nt                  Core nucleotide BLAST database      277.672   2025-12-08
  35 │ ref_euk_rep_genomes      RefSeq Eukaryotic Representative…   482.557   2025-04-12
  36 │ nr                       All non-redundant GenBank CDS tr…   644.21    2025-12-15
  37 │ nt_euk                   Eukaryota nt                        688.86    2025-12-01
  38 │ nt                       Nucleotide collection (nt)          869.168   2025-12-08
```

Simulated sequences
- Depth:
    - low (10 x)
    - medium (100 x)
    - high (1,000 x)
- Diversity:
    - Low (isolate)
    - medium (defined community of 10)
    - high (random community of 100)
- Balance
    - equal representation
    - random representation
    - log-scaled random representation

## Target datasets

### Cami
- https://cami-challenge.org/reference-databases/
- https://cami-challenge.org/datasets/
- datasets with known targets:
    - https://cami-challenge.org/datasets/Toy%20Test%20Dataset%20Low%20Complexity/
    - https://cami-challenge.org/datasets/Toy%20Test%20Dataset%20Medium%20Complexity/
    - https://cami-challenge.org/datasets/Toy%20Test%20Dataset%20High%20Complexity/

### Multi-vendor reference
- [https://www.ncbi.nlm.nih.gov/bioproject/1092695](ONT for ATCC (MSA-1002, MSA-1004, MSA-1005) and Zymobiomics (D6300)) 

### ATCC

- [10 Strain Even Mix Genomic Material](https://www.atcc.org/products/msa-1000)
- [10 Strain Staggered Mix Genomic Material](https://www.atcc.org/products/msa-1001)
- [20 Strain Even Mix Genomic Material](https://www.atcc.org/products/msa-1002)
- [20 Strain Staggered Mix Genomic Material](https://www.atcc.org/products/msa-1003)
- [Oral Microbiome Genomic Mix](https://www.atcc.org/products/msa-1004)
- [Skin Microbiome Genomic Mix](https://www.atcc.org/products/msa-1005)
- [Gut Microbiome Genomic Mix](https://www.atcc.org/products/msa-1006)
- [Vaginal Microbiome Genomic Mix](https://www.atcc.org/products/msa-1007)
- [Virome Nucleic Acid Mix](https://www.atcc.org/products/msa-1008)
- [Mycobiome Genomic DNA Mix](https://www.atcc.org/products/msa-1010)
- [10 Strain Even Mix Whole Cell Material](https://www.atcc.org/products/msa-2003)
- [Oral Microbiome Whole Cell Mix](https://www.atcc.org/products/msa-2004)
- [Skin Microbiome Whole Cell Mix](https://www.atcc.org/products/msa-2005)
- [Gut Microbiome Whole cell Mix](https://www.atcc.org/products/msa-2006)
- [Vaginal Microbiome Whole Cell Mix](https://www.atcc.org/products/msa-2007)
- [Virome Virus Mix](https://www.atcc.org/products/msa-2008)
- [Mycobiome Whole Cell Mix](https://www.atcc.org/products/msa-2010)
- [ABRF-MGRG 10 Strain Staggered Mix Genomic Material](https://www.atcc.org/products/msa-3002)
- [ABRF-MGRG 6 Strain Even Mix Genomic Material](https://www.atcc.org/products/msa-3000)
- [ABRF-MGRG 10 Strain Even Mix Genomic Material](https://www.atcc.org/products/msa-3001)
- [Metagenomic Control Material for Pathogen Detection](https://www.atcc.org/products/msa-4000)

### ZymoBiomics

#### [HMW DNA standard](https://www.zymoresearch.com/products/zymobiomics-hmw-dna-standard)
    - [nanopore](https://www.ncbi.nlm.nih.gov/bioproject/779673)
#### [ZymoBIOMICS Gut Microbiome Standard](https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard)
    - https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip
    - https://www.ncbi.nlm.nih.gov/bioproject/1215043 - illumina
    - https://www.ncbi.nlm.nih.gov/bioproject/1215044 - promethion
    - https://www.ncbi.nlm.nih.gov/bioproject/804004 - ONT
    - https://www.ncbi.nlm.nih.gov/bioproject/992364 - synthetic short reads 
    - https://www.ncbi.nlm.nih.gov/bioproject/680590 - pacbio sequel II hifi

#### [ZymoBiomics Microbial Community DNA standard](https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-dna-standard)
- https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
- https://www.ncbi.nlm.nih.gov/bioproject/1271121 - illumina
- https://www.ncbi.nlm.nih.gov/bioproject/1192878 - ONT
- https://www.ncbi.nlm.nih.gov/bioproject/564477 - promethion and gridion
- https://www.ncbi.nlm.nih.gov/bioproject/934869 - nanopore
- https://www.ncbi.nlm.nih.gov/bioproject/1110296 - pacbio hifi
- https://www.ncbi.nlm.nih.gov/bioproject/1084203 - varying extraction methods
- https://www.ncbi.nlm.nih.gov/bioproject/934863 - testing various extraction methods
- https://www.ncbi.nlm.nih.gov/bioproject/699918 - doesn’t say
- https://www.ncbi.nlm.nih.gov/bioproject/648136 - doesn’t say
- https://www.ncbi.nlm.nih.gov/bioproject/601657 - doesn’t say
- https://www.ncbi.nlm.nih.gov/bioproject/1380770 - doesn’t say

### NIST

#### [RM 8376 Microbial Pathogen DNA Standards for Detection and Identification](https://www.nist.gov/programs-projects/rm-8376-microbial-pathogen-dna-standards-detection-and-identification)
- https://tsapps.nist.gov/srmext/certificates/8376.pdf
- https://www.ncbi.nlm.nih.gov/bioproject/605254
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1074773/

#### Genome in a Bottle

https://github.com/genome-in-a-bottle/giab_latest_release
https://github.com/genome-in-a-bottle/giab_data_indexes
https://www.nist.gov/programs-projects/genome-bottle
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
https://github.com/marbl/HG002
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA200694

### MBARC-26

- [https://www.ncbi.nlm.nih.gov/sra/?term=SRX1836715](pacbio)
- [https://www.ncbi.nlm.nih.gov/sra/?term=SRX1836716](Illumina HiSeq)

### Human microbiome project

https://github.com/awslabs/open-data-registry/blob/main/datasets/human-microbiome-project.yaml
https://portal.hmpdacc.org/
^ no longer available?