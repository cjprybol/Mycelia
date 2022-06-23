# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.236.0/containers/codespaces-linux/.devcontainer/base.Dockerfile

# FROM mcr.microsoft.com/vscode/devcontainers/universal:2-focal
# FROM mcr.microsoft.com/vscode/devcontainers/universal:2
FROM mcr.microsoft.com/vscode/devcontainers/universal:2.0.2

# originally I had been using this variable to get back to the right working directory
# WORKDIR $CODESPACE_VSCODE_FOLDER
# I've lost this variable somehow, but not sure what happened
# recording the $PWD as a variable also isn't working so just record the repo that we're working in
ENV REPO="Mycelia"
ENV REPO_DIR=/workspaces/$REPO

# ** [Optional] Uncomment this section to install additional packages. **
USER root

# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>

# python
# pip freeze > requirements.txt
# COPY requirements.txt .
# pip install requirements.txt

# Pick your poison
# RUN conda install -c conda-forge mamba
# ENV CONDA_CHOICE=mamba
ENV CONDA_CHOICE=conda

# conda
# write out current installations
# conda list --explicit > spec-file.txt
# COPY conda-environment.txt .
# create new environment
# RUN conda create --name myenv --file conda-environment.txt
# install environment into default environment
# RUN conda install --file spec-file.txt

RUN $CONDA_CHOICE install -c conda-forge -c bioconda snakemake

# install R & R kernel
RUN $CONDA_CHOICE install -y -c conda-forge r-base
ENV PATH=/opt/conda/bin:$PATH

# install additional R packages here
# COPY R-requirements.R .

# install precompiled binaries
WORKDIR /installations

ENV JULIA_VERSION_FULL=1.7.3
ENV JULIA_VERSION=1.7

# install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/$JULIA_VERSION/julia-$JULIA_VERSION_FULL-linux-x86_64.tar.gz && \
    tar -xzf julia-$JULIA_VERSION_FULL-linux-x86_64.tar.gz && \
    rm julia-$JULIA_VERSION_FULL-linux-x86_64.tar.gz
ENV PATH=/installations/julia-$JULIA_VERSION_FULL/bin:$PATH

# install julia kernel
RUN julia -e 'import Pkg; Pkg.add("IJulia"); Pkg.build("IJulia")'

# COPY Project.toml $HOME/.julia/environments/v1.6/
# COPY Manifest.toml $HOME/.julia/environments/v1.6/
COPY Project.toml .
RUN julia -e 'import Pkg; Pkg.instantiate()'

# install bash kernel
# RUN pip install bash_kernel

# copy templates into this directory for them to automatically pop up
# /usr/local/share/jupyter/labextensions/jupyterlab_templates

# adding templates
# https://github.com/jpmorganchase/jupyterlab_templates#adding-templates

# Install jupyter templates
RUN pip install jupyter_contrib_nbextensions && \
    jupyter contrib nbextension install && \
    pip install jupyterlab_templates && \
    jupyter labextension install jupyterlab_templates && \
    jupyter serverextension enable --py jupyterlab_templates

# install papermill, which will be our script driver
RUN python3 -m pip install papermill

WORKDIR $REPO_DIR
USER codespace

