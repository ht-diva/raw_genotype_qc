#!/usr/bin/env python3
import argparse, csv
p=argparse.ArgumentParser(); p.add_argument("--smiss",required=True); p.add_argument("--threshold",type=float,required=True); p.add_argument("--keep",required=True); p.add_argument("--plot",required=True); a=p.parse_args()
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
with open(a.smiss) as f:
    r=csv.DictReader(f,delimiter="\t",skipinitialspace=True); rows=list(r)
key=lambda row,name: row.get(name) or row.get("#"+name)
vals=[float(x["F_MISS"]) for x in rows]
with open(a.keep,"w") as out:
    for x in rows:
        if float(x["F_MISS"]) <= a.threshold: out.write(f'{key(x,"FID")}\t{key(x,"IID")}\n')
plt.hist(vals,bins=50); plt.axvline(a.threshold,color="red",ls="--"); plt.xlabel("Sample missingness"); plt.ylabel("Samples"); plt.tight_layout(); plt.savefig(a.plot,dpi=160)

