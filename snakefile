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
        f"busco_output/{SAMPLE_ID}/best_fungi_bin.txt",
        f"busco_output/{SAMPLE_ID}/copied_all_fungi_bins.ok"

def gcf_url(wildcards):
    acc = wildcards.id
    assembly_name = config.get("assembly", f"{acc}_genomic")
    prefix = acc.replace("GCA_", "").replace("GCF_", "").split(".")[0]
    parts = [prefix[i:i+3] for i in range(0, 9, 3)]
    root_dir = acc.split("_")[0]  # GCA или GCF

    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{root_dir}/{'/'.join(parts)}/"
        f"{acc}_{assembly_name}/{acc}_{assembly_name}_genomic.fna.gz"
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
    threads: 4
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
    threads: 4
    run:
        fasta_files = glob.glob(f"maxbin_output/{wildcards.id}/*.fasta")
        for fasta in fasta_files:
            binname = os.path.basename(fasta).replace(".fasta", "")
            outdir = f"busco_output/{wildcards.id}/"
            shell(f"""
                busco -i {fasta} -o run_{binname} \
                    -l {params.lineage} -m genome -c {threads} \
                    --offline --out_path {outdir}
            """)
            faa_dir = os.path.join(outdir, f"run_{binname}", "busco_sequences", "single_copy_busco_sequences")
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
    threads: 4
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
        all_bins = []

        for path in input.summaries:
            print(" -", path)
            match = re.search(r"run_(.+?)/short_summary", path)
            if not match:
                continue
            binname = match.group(1)
            with open(path) as f:
                for line in f:
                    if "C:" in line:
                        c_match = re.search(r"C:(\d+\.\d+)%", line)
                        if c_match:
                            c_value = float(c_match.group(1))
                            print(f"Bin {binname}: C = {c_value}")
                            all_bins.append((binname, c_value))
                        break

        # filter by C > 10%
        passing = [(b, c) for b, c in all_bins if c > 10.0]

        if passing:
            selected = passing
        else:
            print("Warning: No bins found with C > 10%, selecting best-scoring bin.")
            if all_bins:
                selected = [max(all_bins, key=lambda x: x[1])]
            else:
                raise ValueError("No valid BUSCO summary files found at all.")

        with open(output[0], "w") as out:
            for binname, _ in selected:
                out.write(binname + "\n")

rule copy_single_copy_buscos_fungi:
    input:
        best_bins = "busco_output/{id}/best_fungi_bin.txt"
    output:
        touch("busco_output/{id}/copied_all_fungi_bins.ok")
    run:
        import os, glob, shutil, pathlib

        sample_id = wildcards.id

        with open(input.best_bins) as f:
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

        pathlib.Path(output[0]).touch()