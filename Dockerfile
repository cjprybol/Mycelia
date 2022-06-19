# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.236.0/containers/codespaces-linux/.devcontainer/base.Dockerfile

FROM mcr.microsoft.com/vscode/devcontainers/universal:2-focal

ENV HOME /home/codespace/

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
RUN mamba install -y -c conda-forge r-base
# RUN conda install -y -c r r-essentials
ENV PATH=/opt/conda/bin:$PATH
RUN Rscript -e 'install.packages("IRkernel",repos = "http://cran.us.r-project.org");IRkernel::installspec()'
# install additional R packages here
# COPY R-requirements.R .

# install precompiled binaries
WORKDIR /installations

# install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.6-linux-x86_64.tar.gz && \
    tar -xzf julia-1.6.6-linux-x86_64.tar.gz && \
    rm julia-1.6.6-linux-x86_64.tar.gz
ENV PATH=/installations/julia-1.6.6/bin:$PATH

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
#  https://github.com/jpmorganchase/jupyterlab_templates#adding-templates

# install papermill, which will be our script driver

WORKDIR $CODESPACE_VSCODE_FOLDER
USER codespace
