#!/bin/bash

declare -a configs=(
    "GCF_007210705.1 Tcal_SD_v2.1"
)

for entry in "${configs[@]}"; do
    set -- $entry
    id=$1
    assembly=$2
    echo ">>> Running Snakemake for $id ($assembly)"
    snakemake --cores 8 --use-conda --rerun-incomplete --config id=$id assembly=$assembly
done
