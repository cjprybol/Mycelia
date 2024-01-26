Acquire
- [x] refseq viral
- [x] nt viral
- [x] IMG/VR viral
    - copied using GLOBUS from JGI to NERSC, then rclone to Stanford Google Drive, Then rclone Google Drive to SCG3
    - rclone copy --progress $HOME/workspace/JGI stanford_viral_exposome:JGI
    - rclone copy --progress exposome:JGI $HOME/workspace/JGI
- [x] ICTV reference viral

map reads to each

do variant calling

do binned co-assembly

associate taxa/variant abundance with other taxa/host/location/season/