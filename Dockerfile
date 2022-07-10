# See here for image contents: https://hub.docker.com/r/jupyter/datascience-notebook/

FROM jupyter/datascience-notebook

# We want to run common-debian.sh from here:
# https://github.com/microsoft/vscode-dev-containers/tree/main/script-library#development-container-scripts
# But that script assumes that the main non-root user (in this case jovyan)
# is in a group with the same name (in this case jovyan).  So we must first make that so.
COPY .devcontainer/library-scripts/common-debian.sh /tmp/library-scripts/
USER root

# [Optional] Uncomment this section to install additional OS packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libfuse2 \
    fuse \
    graphviz

# RUN modprobe fuse

RUN apt-get update \
 && groupadd jovyan \
 && usermod -g jovyan -a -G users jovyan \
 && bash /tmp/library-scripts/common-debian.sh \
 && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts

# [Optional] If your pip requirements rarely change, uncomment this section to add them to the image.
# COPY requirements.txt /tmp/pip-tmp/
# RUN pip3 --disable-pip-version-check --no-cache-dir install -r /tmp/pip-tmp/requirements.txt \
#    && rm -rf /tmp/pip-tmp

RUN pip3 --disable-pip-version-check --no-cache-dir install \
    papermill \
    cookiecutter \
    && rm -rf /tmp/pip-tmp

RUN pip3 --disable-pip-version-check --no-cache-dir install jupyter_contrib_nbextensions && \
    jupyter contrib nbextension install && \
    pip3 --disable-pip-version-check --no-cache-dir install jupyterlab_templates && \
    jupyter labextension install jupyterlab_templates && \
    jupyter serverextension enable --py jupyterlab_templates && \
    rm -rf /tmp/pip-tmp

# get python system site packages folder
# python -c "import site; print(site.getsitepackages()[0])"
# delete original sample here
# copy in good sample
RUN export JUPYTER_TEMPLATE_DIR=$(python -c "import site; print(site.getsitepackages()[0])")/jupyterlab_templates/templates/jupyterlab_templates && \
    echo $JUPYTER_TEMPLATE_DIR && \
    rm $JUPYTER_TEMPLATE_DIR/Sample.ipynb && \
    wget https://raw.githubusercontent.com/cjprybol/structured-science/master/YYYY-MM-DD-descriptive-title.ipynb && \
    mv YYYY-MM-DD-descriptive-title.ipynb $JUPYTER_TEMPLATE_DIR/YYYY-MM-DD-descriptive-title.ipynb

# RUN conda install -c conda-forge mamba
# ENV CONDA_CHOICE=mamba
ENV CONDA_CHOICE=conda
RUN $CONDA_CHOICE install \
    -c conda-forge \
    -c bioconda \
    snakemake \
    taxonkit \
    samtools

RUN wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && \
    tar -zxvf taxdump.tar.gz && \
    mkdir -p $HOME/.taxonkit && \
    cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit && \
    rm *.dmp taxdump.tar.gz

# install Rclone
RUN curl https://rclone.org/install.sh | sudo bash

# COPY Project.toml .
# RUN julia -e 'import Pkg; Pkg.instantiate()'
RUN julia -e 'import Pkg; Pkg.develop(path="/workspaces/Mycelia"); import Mycelia'

USER jovyan

RUN mkdir -p /home/jovyan/.config/rclone/
COPY rclone.conf /home/jovyan/.config/rclone/rclone.conf
RUN mkdir -p /home/jovyan/rclone-mounts
