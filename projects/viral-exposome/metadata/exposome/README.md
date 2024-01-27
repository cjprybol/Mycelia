- BioProject: https://www.ncbi.nlm.nih.gov/bioproject/421162
- SRA sequencing files and the source of the files in this folder
    - https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=421162
- BioSamples for metadata
    - https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&from_uid=421162

```bash
yq --input-format "xml" --output-format "json" biosample_result.xml
```