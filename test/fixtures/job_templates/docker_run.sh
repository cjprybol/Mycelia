#!/usr/bin/env bash
# TEMPLATE_ID=docker/run.sh
set -euo pipefail

docker run --rm --volume '/tmp/data:/workspace:ro' --workdir '/workspace' --env 'OMP_NUM_THREADS=4' 'ghcr.io/acme/example:1.0.0' bash -lc 'python main.py --input /workspace/input.tsv'
