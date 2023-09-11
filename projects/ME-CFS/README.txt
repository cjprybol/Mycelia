14 Initial samples to try

ME/CFS patients - blinded to the nature of the 14 samples

Sample NGS QC reports
1. Assemble
2. Annotate
3. Classify
4. Integrate


what pathways are present
what organisms are present
what proteins are present


Run me to backup (run the second version if we hit API throttling limits)
```bash
rclone copy --progress $HOME/workspace/Mycelia/projects/ME-CFS/data google_drive:"Projects/BR-WGS data-Aug2023/"
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress ... ...
```
copy from Google Drive back to local
```bash
rclone copy --progress google_drive:"Projects/BR-WGS data-Aug2023/" $HOME/workspace/Mycelia/projects/ME-CFS/data
rclone copy --progress --exclude=*.{fastq.gz,fq.gz,bam}
```

copy Globus-JGI reference DB into cloud for shuffling back and forth to non-globus machines
note - mmseqs db is already in the drive!
```bash
rclone copy --verbose --drive-chunk-size 2G --drive-upload-cutoff 1T --tpslimit 1 --progress $HOME/workspace/JGI google_drive:Projects/reference-databases
```
