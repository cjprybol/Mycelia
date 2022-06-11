# get list of versions from here
# https://github.com/microsoft/vscode-dev-containers/tree/main/containers/codespaces-linux/history
# or run
# devcontainer-info content-url from a base environment
FROM mcr.microsoft.com/vscode/devcontainers/universal:2.0.2

RUN cd $HOME && \
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.6-linux-x86_64.tar.gz && \
    tar -xzf julia-1.6.6-linux-x86_64.tar.gz && \
    rm julia-1.6.6-linux-x86_64.tar.gz

ENV PATH="$HOME/julia-1.6.6/bin:$PATH"

RUN cd $CODESPACE_VSCODE_FOLDER