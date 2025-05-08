import glob
import os

SAMPLE_ID = config["id"]
ASSEMBLY_NAME = config.get("assembly", f"{SAMPLE_ID}_genomic")

rule all:
    input:
        f"input_reads/{SAMPLE_ID}.fna",
        f"alignment/{SAMPLE_ID}_abundance.tsv",
        f"maxbin_output/{SAMPLE_ID}/{SAMPLE_ID}.summary",
        f"busco_output/{SAMPLE_ID}/done",
        f"busco_output/{SAMPLE_ID}/done_microsporidia",
        f"busco_output/{SAMPLE_ID}/best_fungi_bin.txt"

def gcf_url(wildcards):
    gcf_id = wildcards.id
    assembly_name = config.get("assembly", f"{gcf_id}_genomic")
    prefix = gcf_id.replace("GCF_", "").split(".")[0]
    parts = [prefix[i:i+3] for i in range(0, 9, 3)]
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{'/'.join(parts)}/"
        f"{gcf_id}_{assembly_name}/{gcf_id}_{assembly_name}_genomic.fna.gz"
    )

rule download_genome:
    output:
        fasta = "input_reads/{id}.fna",
        gz = "input_reads/{id}.fna.gz"
    params:
        url = gcf_url
    shell:
        """
        wget -c {params.url} -O {output.gz}
        gunzip -c {output.gz} > {output.fasta}
        """

rule simulate_abundance:
    input:
        "input_reads/{id}.fna"
    output:
        abundance="alignment/{id}_abundance.tsv",
        coverage="maxbin_output/{id}_coverage.tsv"
    run:
        from Bio import SeqIO
        lines = []
        for rec in SeqIO.parse(input[0], "fasta"):
            lines.append(f"{rec.id}\t{len(rec.seq)}\t300\n")
        with open(output.abundance, "w") as fa, open(output.coverage, "w") as fc:
            fa.writelines(lines)
            fc.writelines(lines)

rule maxbin:
    input:
        contigs = "input_reads/{id}.fna",
        abund = "alignment/{id}_abundance.tsv"
    output:
        summary = "maxbin_output/{id}/{id}.summary"
    threads: 20
    shell:
        """
        mkdir -p maxbin_output/{wildcards.id}
        run_MaxBin.pl -contig {input.contigs} -abund {input.abund} -out maxbin_output/{wildcards.id}/{wildcards.id}
        """

rule busco_all_bins:
    input:
        summary = "maxbin_output/{id}/{id}.summary"
    output:
        touch("busco_output/{id}/done")
    params:
        lineage = "busco/fungi_odb10"
    threads: 20
    run:
        fasta_files = glob.glob(f"maxbin_output/{wildcards.id}/*.fasta")
        for fasta in fasta_files:
            binname = os.path.basename(fasta).replace(".fasta", "")
            outdir = f"busco_output/{wildcards.id}/"
            shell(f"""
                busco -i {fasta} -o run_{binname} \
                    -l {params.lineage} -m genome -c {threads} \
                    --offline -f --out_path {outdir}
            """)
            faa_dir = os.path.join(outdir, f"run_{binname}_fungi", "busco_sequences", "single_copy_busco_sequences")
            merged_faa = f"{outdir}/sc_{binname}_merged.faa"
            shell(f"cat {faa_dir}/*.faa > {merged_faa} || touch {merged_faa}")
        shell(f"touch {output[0]}")

rule busco_microsporidia_bins:
    input:
        summary = "maxbin_output/{id}/{id}.summary",
        fungi = "busco_output/{id}/done"
    output:
        touch("busco_output/{id}/done_microsporidia")
    params:
        lineage = "busco/microsporidia_odb10"
    threads: 20
    run:
        fasta_files = glob.glob(f"maxbin_output/{wildcards.id}/*.fasta")
        for fasta in fasta_files:
            binname = os.path.basename(fasta).replace(".fasta", "")
            outdir = f"busco_output/{wildcards.id}/{binname}_micro"
            os.makedirs(outdir, exist_ok=True)
            shell(f"""
                busco -i {fasta} -o {binname}_micro \
                    -l {params.lineage} -m genome -c {threads} \
                    --offline -f --out_path {outdir}
            """)
            run_name = f"run_{binname}_micro"
            faa_dir = os.path.join(outdir, run_name, "busco_sequences", "single_copy_busco_sequences")
            merged_faa = f"busco_output/{wildcards.id}/sc_{binname}_micro_merged.faa"
            shell(f"cat {faa_dir}/*.faa > {merged_faa} || touch {merged_faa}")
        shell(f"touch {output[0]}")

rule best_bin_fungi:
    input:
        done = "busco_output/{id}/done",
        summaries = lambda wildcards: glob.glob(f"busco_output/{wildcards.id}/run_*/short_summary.*.txt")
    output:
        "busco_output/{id}/best_fungi_bin.txt"
    run:
        import re
        print("Found summary files:")
        for path in input.summaries:
            print(" -", path)
        best_bin = None
        best_s_value = -1.0
        for path in input.summaries:
            match = re.search(r"run_(.+?)/short_summary", path)
            if not match:
                continue
            binname = match.group(1)
            with open(path) as f:
                for line in f:
                    if "C:" in line and "S:" in line:
                        s_match = re.search(r"S:(\d+\.\d+)%", line)
                        if s_match:
                            s_value = float(s_match.group(1))
                            print(f"Bin {binname}: S = {s_value}")
                            if s_value > best_s_value:
                                best_s_value = s_value
                                best_bin = binname
                        break
        assert best_bin is not None, "No valid BUSCO summaries with S: found"
        with open(output[0], "w") as out:
            out.write(best_bin + "\n")

rule copy_single_copy_buscos_fungi:
    input:
        best_bin = "busco_output/{id}/best_fungi_bin.txt"
    output:
        directory("data_for_tree/busco_{id}/run_fungi_odb10/busco_sequences/single_copy_busco_sequences")
    run:
        with open(input.best_bin) as f:
            binname = f.read().strip()
        src_dir = f"busco_output/{wildcards.id}/run_{binname}/run_fungi_odb10/busco_sequences/single_copy_busco_sequences"
        shell(f"""
            mkdir -p {output[0]}
            cp {src_dir}/*.faa {output[0]}/
        """)