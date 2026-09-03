#!/usr/bin/env python3
import argparse, csv
from collections import Counter

p = argparse.ArgumentParser()
p.add_argument("--fam", required=True); p.add_argument("--bim", required=True)
p.add_argument("--metadata", required=True); p.add_argument("--update-sex", required=True)
p.add_argument("--duplicate-samples", required=True); p.add_argument("--duplicate-variants", required=True)
p.add_argument("--duplicate-positions", required=True)
a = p.parse_args()

fam = [x.split() for x in open(a.fam) if x.strip()]
ids = [(x[0], x[1]) for x in fam]
counts = Counter(ids)
with open(a.duplicate_samples, "w") as out:
    out.write("FID\tIID\tCOUNT\n")
    for (fid, iid), n in sorted(counts.items()):
        if n > 1: out.write(f"{fid}\t{iid}\t{n}\n")

bim = [x.split() for x in open(a.bim) if x.strip()]
id_counts = Counter(x[1] for x in bim)
pos_counts = Counter((x[0].lower().removeprefix("chr"), x[3]) for x in bim)
with open(a.duplicate_variants, "w") as out:
    out.write("ID\tCOUNT\n")
    for vid, n in sorted(id_counts.items()):
        if n > 1: out.write(f"{vid}\t{n}\n")
with open(a.duplicate_positions, "w") as out:
    out.write("CHROM\tPOS\tCOUNT\n")
    for (chrom, pos), n in sorted(pos_counts.items()):
        if n > 1: out.write(f"{chrom}\t{pos}\t{n}\n")

known = set(ids)
with open(a.metadata, newline="") as src, open(a.update_sex, "w") as out:
    rows = csv.DictReader(src, delimiter="\t")
    required = {"FID", "IID", "SEX"}
    if not required.issubset(rows.fieldnames or []): raise SystemExit("Metadata requires FID, IID, SEX")
    out.write("#FID\tIID\tSEX\n")
    for r in rows:
        key = (r["FID"], r["IID"])
        if key in known: out.write(f'{key[0]}\t{key[1]}\t{r["SEX"]}\n')

