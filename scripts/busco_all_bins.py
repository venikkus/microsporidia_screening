import glob
import os

fasta_files = glob.glob(f"maxbin_output/{snakemake.wildcards.id}/*.fasta")
for fasta in fasta_files:
    binname = os.path.basename(fasta).replace(".fasta", "")
    outdir = f"busco_output/{snakemake.wildcards.id}/"
    os.makedirs(outdir, exist_ok=True)
    os.system(f"""
        busco -i {fasta} -o run_{binname} \
            -l {snakemake.params.lineage} -m genome -c {snakemake.threads} \
            --offline --out_path {outdir}
    """)
    faa_dir = os.path.join(outdir, f"run_{binname}", "busco_sequences", "single_copy_busco_sequences")
    merged_faa = f"{outdir}/sc_{binname}_merged.faa"
    os.system(f"cat {faa_dir}/*.faa > {merged_faa} || touch {merged_faa}")
os.system(f"touch {snakemake.output[0]}")
