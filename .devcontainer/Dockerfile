FROM mcr.microsoft.com/devcontainers/base:debian-11

# Install prerequisites
RUN apt-get update && \
    apt-get install -y wget ca-certificates curl

# Install juliaup and use it to install Julia LTS
RUN curl -fsSL https://install.julialang.org | bash -s -- --yes && \
    /root/.juliaup/bin/juliaup add lts && \
    ln -s /root/.juliaup/bin/julia /usr/local/bin/julia

# Install Rust (via rustup)
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:/root/.juliaup/bin:${PATH}"

# Install Python 3 and pip
RUN apt-get update && \
    apt-get install -y python3 python3-pip

# Install Node.js (LTS) and npm
RUN curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - && \
    apt-get install -y nodejs

# Clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*