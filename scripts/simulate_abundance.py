from Bio import SeqIO

lines = []
for rec in SeqIO.parse(snakemake.input[0], "fasta"):
    lines.append(f"{rec.id}\t{len(rec.seq)}\t300\n")
with open(snakemake.output.abundance, "w") as fa, open(snakemake.output.coverage, "w") as fc:
    fa.writelines(lines)
    fc.writelines(lines)