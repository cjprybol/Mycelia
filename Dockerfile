# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.236.0/containers/codespaces-linux/.devcontainer/base.Dockerfile

FROM mcr.microsoft.com/vscode/devcontainers/universal:2-focal

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

RUN conda install -c conda-forge mamba
RUN mamba install -c conda-forge -c bioconda snakemake

# conda
# write out current installations
# conda list --explicit > spec-file.txt
# COPY conda-environment.txt .
# create new environment
# RUN conda create --name myenv --file conda-environment.txt
# install environment into default environment
# RUN conda install --file spec-file.txt

# install R & R kernel
# with conda
# RUN mamba install -y -c conda-forge r-base
# RUN conda install -y -c r r-essentials
# ENV PATH=/opt/conda/bin:$PATH

# with linux package manager
# update indices
RUN apt update -qq -y
# install two helper packages we need
RUN apt install -y --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt install -y --no-install-recommends r-base
RUN Rscript -e 'install.packages("IRkernel",repos = "http://cran.us.r-project.org");IRkernel::installspec()'
# install additional R packages here
# COPY R-requirements.R .

# install docker https://docs.docker.com/engine/install/ubuntu/
RUN apt update -qq -y
RUN apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg
RUN echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
RUN apt-get update -qq -y
RUN apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin


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

# Install jupyter templates
RUN pip install jupyter_contrib_nbextensions
RUN jupyter contrib nbextension install
RUN pip install jupyterlab_templates
RUN jupyter labextension install jupyterlab_templates
RUN jupyter serverextension enable --py jupyterlab_templates

# copy templates into this directory for them to automatically pop up
# /usr/local/share/jupyter/labextensions/jupyterlab_templates

# adding templates
# https://github.com/jpmorganchase/jupyterlab_templates#adding-templates

# install papermill, which will be our script driver
RUN python3 -m pip install papermill

WORKDIR $REPO_DIR
USER codespace