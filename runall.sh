#!/bin/bash

# Load the container path from config.yaml
DATA_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('config.yaml'))['data_dir'])")

docker run --rm \
  -v $DATA_DIR:/Volumes/Alpha/Research/data \
  -v $(pwd):/project \
  -w /project \
  outflows \
  snakemake --cores 1