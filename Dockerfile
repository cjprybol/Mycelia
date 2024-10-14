# See here for image contents: https://hub.docker.com/r/jupyter/datascience-notebook/
FROM quay.io/jupyter/datascience-notebook:latest

USER root

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    gzip

RUN pip3 --disable-pip-version-check --no-cache-dir install --upgrade \
    papermill \
    jupyter-ai \
    langchain-anthropic \
    langchain-google-genai \
    langchain-ollama \
    langchain-openai

RUN curl -fsSL https://ollama.com/install.sh | sudo sh
RUN ollama start & sleep 5 \
    && ollama pull llama3.2 \
    && ollama pull mxbai-embed-large

# 3B parameters, 4Gb RAM
RUN jupyter lab --AiExtension.default_language_model=ollama:llama3.2
RUN jupyter lab --AiExtension.default_embeddings_model=ollama:mxbai-embed-large

RUN mamba install \
    -c conda-forge \
    -c bioconda \
    samtools \
    htslib

# install Rclone
RUN curl https://rclone.org/install.sh | bash