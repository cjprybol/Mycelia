New data - `/labs/mpsnyder/amir/ME-CFS-VCF/FASTQs/Shared-ME-CFS/`

find /labs/mpsnyder/amir/ME-CFS-VCF/FASTQs/Shared-ME-CFS/ -type f -name "*.fq.gz" -exec ls -l {} \; | awk '{sum += $5} END {if (sum >= 1099511627776) {size = sum / 1099511627776; unit = "TB"} else if (sum >= 1073741824) {size = sum / 1073741824; unit = "GB"} else if (sum >= 1048576) {size = sum / 1048576; unit = "MB"} else if (sum >= 1024) {size = sum / 1024; unit = "KB"} else {size = sum; unit = "bytes"}; printf "%.2f %s\n", size, unit}'
30.24 TB

find $HOME/workspace/Mycelia/projects/viral-exposome/data/SRA/ -type f -name "*.fq.gz" -exec ls -l {} \; | awk '{sum += $5} END {if (sum >= 1099511627776) {size = sum / 1099511627776; unit = "TB"} else if (sum >= 1073741824) {size = sum / 1073741824; unit = "GB"} else if (sum >= 1048576) {size = sum / 1048576; unit = "MB"} else if (sum >= 1024) {size = sum / 1024; unit = "KB"} else {size = sum; unit = "bytes"}; printf "%.2f %s\n", size, unit}'
4.5 TB

All charges incurred in past 18 months


get_compute_charges -a mpsnyder -y 2023

[-m month] [-v] [-p]


$1.50 / CPU-day for batch partition usage ($0.0625 / CPU-hour)

SRR6399459 4.33 Gb
Assembly 
Mapping
Blast classification


SRR7365485 9Gb




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