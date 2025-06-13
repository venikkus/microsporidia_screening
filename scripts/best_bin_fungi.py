import re

print("Found summary files:")
all_bins = []

for path in snakemake.input.summaries:
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

with open(snakemake.output[0], "w") as out:
    for binname, _ in selected:
        out.write(binname + "\n")