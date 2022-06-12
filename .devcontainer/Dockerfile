# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.236.0/containers/codespaces-linux/.devcontainer/base.Dockerfile

FROM mcr.microsoft.com/vscode/devcontainers/universal:2-focal

ENV HOME /home/codespace/

# ** [Optional] Uncomment this section to install additional packages. **
USER root

# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>

WORKDIR /installations

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.6-linux-x86_64.tar.gz && \
    tar -xzf julia-1.6.6-linux-x86_64.tar.gz && \
    rm julia-1.6.6-linux-x86_64.tar.gz
ENV PATH=/installations/julia-1.6.6/bin:$PATH

# # COPY Project.toml $HOME/.julia/environments/v1.6/
# # COPY Manifest.toml $HOME/.julia/environments/v1.6/
# # RUN julia -e 'import Pkg; Pkg.instantiate()'

WORKDIR $CODESPACE_VSCODE_FOLDER
USER codespace
