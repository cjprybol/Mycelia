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


organize kraken reports and upload to Google Drive
```
mkdir -p $HOME/workspace/Mycelia/projects/ME_CFS/results/kraken-reports
cp $HOME/workspace/Mycelia/projects/ME_CFS/data/samples/*/*_trimgalore/*unmapped_kraken/*.html $HOME/workspace/Mycelia/projects/ME_CFS/results/kraken-reports/.
rclone copy --verbose --progress $HOME/workspace/Mycelia/projects/ME_CFS/results ME_CFS:results
```

Run me to backup (run the second version if we hit API throttling limits)
```bash
rclone copy --verbose --progress --retries 5 --drive-chunk-size 1G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/Mycelia/projects/ME_CFS/data ME_CFS:data
```

copy from Google Drive back to local
```bash
rclone copy --progress --retries 5 ME_CFS:data $HOME/workspace/Mycelia/projects/ME_CFS/data
```

helpful filters
```bash
--include "*kraken*"
--exclude=*.{fastq.gz,fq.gz,bam}

# e.g.
rclone copy --verbose --progress --retries 5 --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --include "*kraken*" $HOME/workspace/Mycelia/projects/ME_CFS/data ME_CFS:data
```

rclone copy --verbose --progress --retries 5 --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --include "*kraken*" $HOME/workspace/Mycelia/projects/ME_CFS/data google_drive:Projects/ME_CFS/data


rclone copy --verbose --progress --retries 5 --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 $HOME/workspace/Mycelia/projects/ME_CFS/results google_drive:Projects/ME_CFS/results


rclone copy --verbose --progress --retries 5 --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --include "*kraken*" $HOME/workspace/Mycelia/projects/ME_CFS/results google_drive:Projects/ME_CFS/results

rclone copy --verbose --progress --retries 5 --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --include "*kraken*" $HOME/workspace/Mycelia/projects/ME_CFS/data google_drive:Projects/ME_CFS/data

rclone copy --verbose --progress --retries 5 --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --include "*kraken*" google_drive:Projects/ME_CFS/data $HOME/workspace/Mycelia/projects/ME_CFS/data