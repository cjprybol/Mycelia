#!/bin/bash
# Example: Only run in remote environments
if [ "$CLAUDE_CODE_REMOTE" != "true" ]; then
  exit 0
fi

curl -fsSL https://install.julialang.org | sh -s -- --default-channel lts --yes && \
    ln -s $HOME/.juliaup/bin/julia /usr/local/bin/julia