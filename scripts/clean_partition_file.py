with open(input[0]) as fin, open(output[0], "w") as fout:
    inside_sets = False
    for line in fin:
        stripped = line.strip()
        if stripped.lower().startswith("begin sets"):
            inside_sets = True
            fout.write(line)
            continue
        elif stripped.lower().startswith("end"):
            inside_sets = False
            fout.write(line)
            continue

        if inside_sets and "SUPERMATRIX.phylip:" in line:
            # rm "SUPERMATRIX.phylip:"
            cleaned_line = line.replace("SUPERMATRIX.phylip:", "").strip()
            fout.write("    " + cleaned_line + "\n")
        else:
            fout.write(line)
