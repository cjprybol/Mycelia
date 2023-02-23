# sars-cov2-pangenome-analysis


https://registry.opendata.aws/ncbi-covid-19/

`aws s3 ls --no-sign-request s3://sra-pub-sars-cov2/`

https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/

https://registry.opendata.aws/ncbi-sra/

expression
https://trace.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=CORONAVIRIDAE&o=organism_s%3Aa%253Bacc_s%3Bacc_s%3Aa


mount and link storage

This didn't work because we can't download the entire thing to the baby little codespaces disk size

Need to download iteratively

mkdir -p data/ncbi-cli-download
cd data/ncbi-cli-download
datasets download virus genome taxon SARS-CoV-2

conda install -c conda-forge parallel

cut -d, -f1 sequences.csv | tail -n+2 | less