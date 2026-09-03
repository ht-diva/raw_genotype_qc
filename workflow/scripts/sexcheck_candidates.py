#!/usr/bin/env python3
import argparse, csv
p=argparse.ArgumentParser(); p.add_argument("--sexcheck",required=True); p.add_argument("--remove",required=True); a=p.parse_args()
with open(a.sexcheck) as f: rows=list(csv.DictReader(f,delimiter="\t",skipinitialspace=True))
with open(a.remove,"w") as out:
    for r in rows:
        if r.get("STATUS","").upper() not in {"OK","PASS"}: out.write(f'{r.get("#FID",r.get("FID"))}\t{r["IID"]}\n')

