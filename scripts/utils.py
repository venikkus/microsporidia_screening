def gcf_url(wildcards, config):
    acc = wildcards.id
    assembly_name = config.get("assembly", f"{acc}_genomic")
    prefix = acc.replace("GCA_", "").replace("GCF_", "").split(".")[0]
    parts = [prefix[i:i+3] for i in range(0, 9, 3)]
    root_dir = acc.split("_")[0]  # GCA или GCF

    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{root_dir}/{'/'.join(parts)}/"
        f"{acc}_{assembly_name}/{acc}_{assembly_name}_genomic.fna.gz"
    )