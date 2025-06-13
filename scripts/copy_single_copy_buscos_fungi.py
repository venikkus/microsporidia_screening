import os
import glob
import shutil
import pathlib

sample_id = snakemake.wildcards.id

with open(snakemake.input.best_bins) as f:
    full_bin_names = [line.strip() for line in f if line.strip()]

copied = 0
for full_bin in full_bin_names:
    bin_id = full_bin.split(".")[-1]
    run_dir = f"run_{full_bin}"
    src_dir = f"busco_output/{sample_id}/{run_dir}/run_fungi_odb10/busco_sequences/single_copy_busco_sequences"
    dst_dir = f"data_for_tree/busco_{full_bin}/run_fungi_odb10/busco_sequences/single_copy_busco_sequences"

    if not os.path.isdir(src_dir):
        print(f"WARNING: source directory missing: {src_dir}")
        continue

    os.makedirs(dst_dir, exist_ok=True)

    for faa_file in glob.glob(os.path.join(src_dir, "*.faa")):
        shutil.copy(faa_file, dst_dir)
        copied += 1

assert copied > 0, f"No .faa files copied for {sample_id}"

pathlib.Path(snakemake.output[0]).touch()