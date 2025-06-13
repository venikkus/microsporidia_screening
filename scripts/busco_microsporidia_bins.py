import glob
import os

fasta_files = glob.glob(f"maxbin_output/{snakemake.wildcards.id}/*.fasta")
for fasta in fasta_files:
    binname = os.path.basename(fasta).replace(".fasta", "")
    outdir = f"busco_output/{snakemake.wildcards.id}/{binname}_micro"
    os.makedirs(outdir, exist_ok=True)
    os.system(f"""
        busco -i {fasta} -o {binname}_micro \
            -l {snakemake.params.lineage} -m genome -c {snakemake.threads} \
            --offline -f --out_path {outdir}
    """)
    run_name = f"run_{binname}_micro"
    faa_dir = os.path.join(outdir, run_name, "busco_sequences", "single_copy_busco_sequences")
    merged_faa = f"busco_output/{snakemake.wildcards.id}/sc_{binname}_micro_merged.faa"
    os.system(f"cat {faa_dir}/*.faa > {merged_faa} || touch {merged_faa}")
os.system(f"touch {snakemake.output[0]}")
