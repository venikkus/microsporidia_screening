import glob
import os
import scripts.utils as utils

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

rule download_genome:
    output:
        fasta = "input_reads/{id}.fna",
        gz = "input_reads/{id}.fna.gz"
    params:
        url = lambda wildcards: utils.gcf_url(wildcards, config)
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
    script:
        "scripts/simulate_abundance.py"

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
    script:
        "scripts/busco_all_bins.py"

rule busco_microsporidia_bins:
    input:
        summary = "maxbin_output/{id}/{id}.summary",
        fungi = "busco_output/{id}/done"
    output:
        touch("busco_output/{id}/done_microsporidia")
    params:
        lineage = "busco/microsporidia_odb10"
    threads: 4
    script:
        "scripts/busco_microsporidia_bins.py"

rule best_bin_fungi:
    input:
        done = "busco_output/{id}/done",
        summaries = lambda wildcards: glob.glob(f"busco_output/{wildcards.id}/run_*/short_summary.*.txt")
    output:
        "busco_output/{id}/best_fungi_bin.txt"
    script:
        "scripts/best_bin_fungi.py"

rule copy_single_copy_buscos_fungi:
    input:
        best_bins = "busco_output/{id}/best_fungi_bin.txt"
    output:
        touch("busco_output/{id}/copied_all_fungi_bins.ok")
    script:
        "scripts/copy_single_copy_buscos_fungi.py"