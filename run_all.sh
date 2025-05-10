#!/bin/bash

declare -a configs=(
    "GCF_007210705.1 Tcal_SD_v2.1" # Tigriopus californicus
    "GCF_000764305.2 Hazt_2.0.2" # Hyalella azteca v1
    "GCA_000764305.1 Hazt_1.0" # Hyalella azteca v0
    "GCA_016086655.1 UVic_Lsal_1.0" # Lepeophtheirus salmonis
    "GCF_016086655.4 UVic_Lsal_1.4" # Lepeophtheirus salmonis
    "GCF_032884065.1 ASM3288406v1" # Artemia
    "GCA_005281655.1 TAMU_Nfulva_1.0" # Nylanderia fulva v0
    "GCF_005281655.2 TAMU_Nfulva_1.1" # Nylanderia fulva v1
    "GCF_000190715.1 v1.0" # Dictyostelium purpureum
    "GCF_000836235.1 Pxut_1.0" # Papilio xuntus
    "GCA_001017525.1 ASM101752v1" # Teleopsis dalmanni
    "GCA_024586455.1 APGP_CSIRO_Bneo_wtdbg2-racon-allhic-juicebox.fasta" # Bactrocera neohumeralis
)

for entry in "${configs[@]}"; do
    set -- $entry
    id=$1
    assembly=$2
    echo ">>> Running Snakemake for $id ($assembly)"
    snakemake --cores 8 --use-conda --rerun-incomplete --config id=$id assembly=$assembly
done