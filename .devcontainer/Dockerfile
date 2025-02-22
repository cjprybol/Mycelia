FROM mcr.microsoft.com/devcontainers/universal:2-linux

# Install Python and pip
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends python3-pip

# Install Jupyter tools
RUN pip3 install jupyter nbconvert
