#!/usr/bin/env python3
import argparse, csv, statistics
p=argparse.ArgumentParser(); p.add_argument("--het",required=True); p.add_argument("--sd",type=float,required=True); p.add_argument("--remove",required=True); a=p.parse_args()
with open(a.het) as f: rows=list(csv.DictReader(f,delimiter="\t",skipinitialspace=True))
def ratio(r):
    n=float(r["N(NM)"]); return (n-float(r["O(HOM)"]))/n if n else float("nan")
vals=[ratio(r) for r in rows]; valid=[x for x in vals if x==x]; mu=statistics.mean(valid); sd=statistics.stdev(valid) if len(valid)>1 else 0
with open(a.remove,"w") as out:
    for r,x in zip(rows,vals):
        if x==x and abs(x-mu)>a.sd*sd: out.write(f'{r.get("#FID",r.get("FID"))}\t{r["IID"]}\n')

