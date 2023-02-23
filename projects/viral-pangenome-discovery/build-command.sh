# jupyter-repo2docker https://github.com/cjprybol/sars-cov2-pangenome-analysis.git \
#     --image-name sas-cov2-pangenome-analysis-$(date +"%Y%m%d%H%M%S")

# should make image-name dependend on the directory being run in
# so any repo we build will be named after the repo itself

DIR=$(basename $PWD)
TIMESTAMP=$(date +"%Y%m%d%H%M%S")
IMAGE_NAME=$(DIR)-$(date +"%Y%m%d%H%M%S")

jupyter-repo2docker \
    --no-run \
    --push \
    --image-name $IMAGE_NAME \
    --target-repo-dir workspace \
    .
