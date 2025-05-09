#!/bin/bash

declare -a configs=(
    "GCF_007210705.1 Tcal_SD_v2.1"
    "GCF_000764305.2 Hazt_2.0.2"
    "GCF_016086655.4 UVic_Lsal_1.4"
    "GCF_032884065.1 ASM3288406v1"
    "GCF_005281655.2 TAMU_Nfulva_1.1"
    "GCF_000190715.1 v1.0"
)

for entry in "${configs[@]}"; do
    set -- $entry
    id=$1
    assembly=$2
    echo ">>> Running Snakemake for $id ($assembly)"
    snakemake --cores 8 --use-conda --rerun-incomplete --config id=$id assembly=$assembly
done