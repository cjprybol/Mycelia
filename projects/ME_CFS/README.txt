14 Initial samples to try

ME/CFS patients - blinded to the nature of the 14 samples

Sample NGS QC reports
0. NGS QC trim & filter <= 16 hours @ 1core x 4Gb
1. Kraken analysis with PlusPFP-16
if success, Kraken analysis with PlusPFP


1. Assemble - ?? hours @ 8core
2. Annotate
3. Classify
4. Integrate and Analyze

what pathways are present
what organisms are present
what proteins are present

Run me to backup (run the second version if we hit API throttling limits)
```bash
rclone copy --progress $HOME/workspace/Mycelia/projects/ME_CFS/data ME_CFS:data
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress ... ...
```
copy from Google Drive back to local
```bash
rclone copy --progress ME_CFS:data $HOME/workspace/Mycelia/projects/ME_CFS/data
rclone copy --progress --exclude=*.{fastq.gz,fq.gz,bam}
```